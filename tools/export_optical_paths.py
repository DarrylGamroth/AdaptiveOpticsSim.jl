#!/usr/bin/env python3
# Derived from SPIDERSSim.jl b479871 tools/export_physical15.py.
"""Compile spiders-sim's frozen fitted model for the native Julia runtime.

Run with spectre's Python environment. The compiler imports the actual reference
implementation; neither JAX nor Python is needed to execute the exported model.
FITS matrices retain Python (row=y, column=x) axes. FITSIO.jl readers must transpose
each matrix. Complex matrices use a real primary HDU and an imaginary second HDU.
"""

import argparse
import copy
import dataclasses
import hashlib
import pathlib
import subprocess
import sys
from unittest.mock import patch

import numpy as np
from import_compiled_optical_model import model_document, write_model


# Keep these values in parity with spiders_sim.sim.PRESETS. Resolution is pure
# Python so argument tests need not import the reference's JAX runtime.
PRESET_OPTIONS = {
    "fast": dict(grid=512, band=None, ghosts=False, mtf=False),
    "medium": dict(grid=512, band=(None, 0.015, 3), ghosts=True, mtf=True),
    "full": dict(grid=1024, band=(None, 0.015, 3), ghosts=True, mtf=True),
}
FTS_PHANTOMS = frozenset(("stack0.724", "stack1.956", "stack2.707", "stack3.082",
                          "stack3.618", "stack4.368", "stack6.325"))


def positive_int(text):
    value = int(text)
    if value <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return value


def positive_float(text):
    value = float(text)
    if not np.isfinite(value) or value <= 0:
        raise argparse.ArgumentTypeError("must be finite and positive")
    return value


def argument_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    workspace = pathlib.Path(__file__).resolve().parents[2]
    parser.add_argument("--spectre", type=pathlib.Path, default=workspace / "spectre")
    parser.add_argument("--spiders-sim", type=pathlib.Path, default=workspace / "spiders-sim")
    parser.add_argument("--model", type=pathlib.Path, default=workspace / "spiders-sim/models/physical15")
    parser.add_argument("--spectre-revision", help="Required pinned source revision when supplied")
    parser.add_argument("--spiders-sim-revision", help="Required pinned source revision when supplied")
    parser.add_argument("--artifact-manifest-sha256", help="Required hash of frozen manifest.json when supplied")
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--manifest-name", default="model.yaml")
    parser.add_argument("--preset", choices=tuple(PRESET_OPTIONS), default="fast")
    parser.add_argument("--grid", type=positive_int)
    parser.add_argument("--configuration", choices=("reference", "scc-bench", "scc-sky"), default="reference",
                        help="Fitted reference or scoped SCC candidate; LLOWFS is unsupported")
    parser.add_argument("--mode", choices=("sky", "source"),
                        help="reference entrance mode (default sky); SCC configuration fixes the mode")
    parser.add_argument("--dm-base", choices=("fitted", "flat"), default="fitted")
    parser.add_argument("--reference", action="store_true", help="Write deterministic reference products")
    ghosts = parser.add_mutually_exclusive_group()
    ghosts.add_argument("--ghosts", dest="ghosts", action="store_true")
    ghosts.add_argument("--no-ghosts", dest="ghosts", action="store_false")
    parser.set_defaults(ghosts=None)
    spectrum = parser.add_mutually_exclusive_group()
    spectrum.add_argument("--monochromatic", action="store_true")
    spectrum.add_argument("--band-width", type=positive_float, metavar="FRACTION",
                          help="Fractional full band width, strictly less than 2")
    parser.add_argument("--band-centre-m", type=positive_float)
    parser.add_argument("--anchors", type=positive_int)
    parser.add_argument("--channels", type=positive_int)
    mtf = parser.add_mutually_exclusive_group()
    mtf.add_argument("--mtf", dest="mtf", action="store_true", help="Enable fitted deterministic pixel response")
    mtf.add_argument("--no-mtf", dest="mtf", action="store_false", help="Disable fitted deterministic pixel response")
    parser.set_defaults(mtf=None)
    parser.add_argument("--fixed-ghost", choices=("fitted", "off"),
                        help="Fixed empirical ghost field: default off for scc-sky, fitted otherwise")
    parser.add_argument("--apodizer-gamma", choices=("fitted", "nominal"), default="fitted",
                        help="Retain fitted apodizer coverage exponent or set its adjustment to zero")
    return parser


def resolve_options(args, parser):
    """Resolve independent selections without importing JAX or reading artifacts."""
    selected = dict(PRESET_OPTIONS[args.preset])
    configuration_mode = {"scc-bench": "source", "scc-sky": "sky"}.get(args.configuration)
    if configuration_mode is not None and args.mode not in (None, configuration_mode):
        parser.error(f"--configuration {args.configuration} requires --mode {configuration_mode}")
    args.mode = configuration_mode or args.mode or "sky"
    if args.grid is not None:
        selected["grid"] = args.grid
    if args.ghosts is not None:
        selected["ghosts"] = args.ghosts
    if args.mtf is not None:
        selected["mtf"] = args.mtf
    if args.band_width is not None and args.band_width >= 2:
        parser.error("--band-width must be strictly less than 2 (both wavelength edges must be positive)")
    if args.monochromatic:
        if any(value is not None for value in (args.band_centre_m, args.anchors, args.channels)):
            parser.error("--monochromatic cannot be combined with --band-centre-m, --anchors or --channels")
        selected["band"] = None
    elif args.band_width is not None:
        selected["band"] = (args.band_centre_m, args.band_width, args.anchors or 3)
    elif selected["band"] is not None:
        centre, width, anchors = selected["band"]
        selected["band"] = (args.band_centre_m or centre, width, args.anchors or anchors)
    elif any(value is not None for value in (args.band_centre_m, args.anchors, args.channels)):
        parser.error("spectral sampling options on a monochromatic preset require --band-width")
    if selected["band"] is not None:
        centre, width, _ = selected["band"]
        lower, upper = 1 - width / 2, 1 + width / 2
        if lower == upper:
            parser.error("--band-width is too small to resolve distinct wavelength edges")
        if centre is not None:
            edges = (centre * lower, centre * upper)
            if not all(np.isfinite(edge) and edge > 0 and np.isfinite(1 / edge) for edge in edges):
                parser.error("--band-centre-m and --band-width must give finite positive wavelengths and wavenumbers")
    selected["channels"] = args.channels or 8
    selected["fixed_ghost"] = args.fixed_ghost or ("off" if args.configuration == "scc-sky" else "fitted")
    # band=None does not override medium/full upstream; restore all requested
    # settings explicitly when a fast construction base is necessary.
    selected["construction_preset"] = "fast" if selected["band"] is None else args.preset
    args.resolved = selected
    return args


def configured_manifest(manifest, configuration, fixed_ghost):
    """Copy and filter fitted additions, retaining every physical envelope path."""
    result = copy.deepcopy(manifest)
    original = manifest["phantoms"]
    keep = [i for i, phantom in enumerate(original)
            if not (configuration != "reference" and phantom["label"] in FTS_PHANTOMS)
            and not (fixed_ghost == "off" and phantom["source"] is None)]
    result["phantoms"] = [result["phantoms"][i] for i in keep]
    for key in ("phantom_opd_um", "phantom_defocus_mm"):
        result["fast"][key] = np.asarray(manifest["fast"][key])[keep].copy()
    removed = [p["label"] for i, p in enumerate(original) if i not in keep]
    return result, removed


def configured_model(model, manifest, configuration):
    """Remove the FTS-fitted L2B corrections while retaining the SCC lens."""
    if configuration == "reference":
        return model, manifest, [], []
    screens = {name: value for name, value in model.screens.items() if name != "L2B"}
    slots = [slot for slot in model.adjust if slot.split(":", 1)[-1].split(".")[0] == "L2B"]
    adjust = {slot: value for slot, value in model.adjust.items() if slot not in slots}
    selected = dataclasses.replace(model, screens=screens, adjust=adjust)
    selected_manifest = dict(manifest)
    dropped = [name for name in model.screens if name == "L2B"]
    selected_manifest["screens"] = [name for name in manifest["screens"] if name not in dropped]
    selected_manifest["slots"] = [slot for slot in manifest["slots"] if slot not in slots]
    return selected, selected_manifest, dropped, slots


def configure_simulation(sim, args, tier_factory):
    """Replace the cold fitted owner before arrays or oracle products are made."""
    options = args.resolved
    manifest, removed = configured_manifest(sim.manifest, args.configuration, options["fixed_ghost"])
    if removed:
        # Preserve the exact source-generated channels and their top-hat width.
        dsigma = None
        if options["band"] is not None:
            centre, width, _ = options["band"]
            centre = float(np.mean(sim.model.band)) if centre is None else centre
            edges = (centre * (1 - width / 2), centre * (1 + width / 2))
            sig = np.linspace(1 / edges[1], 1 / edges[0], options["channels"] + 1)
            dsigma = float(sig[1] - sig[0])
        sim.tier = tier_factory(manifest, sim.fwd, sim.wavelengths, sim.channels, dsigma,
                                ghosts=options["ghosts"])
        sim.manifest, sim.labels = manifest, sim.tier.labels
    if args.apodizer_gamma == "nominal":
        slot = "gamma:Apodizer"
        if slot not in sim.adjust:
            raise ValueError("The selected model has no fitted gamma:Apodizer adjustment")
        sim.adjust = {**sim.adjust, slot: np.zeros_like(sim.adjust[slot])}
    return removed


class Arrays:
    def __init__(self, output):
        self.output = output
        self.output.mkdir(parents=True, exist_ok=True)
        self.known = {}
        self.sha256 = {}

    def save(self, array):
        from astropy.io import fits
        value = np.asarray(array)
        value = np.ascontiguousarray(value, dtype=np.complex128 if np.iscomplexobj(value) else np.float64)
        if value.ndim != 2:
            raise ValueError(f"Export arrays must be matrices, got {value.shape}")
        if not np.isfinite(value).all():
            raise ValueError(f"Export array contains non-finite values: {value.shape}")
        key = hashlib.sha256(str(value.shape).encode() + value.dtype.str.encode() + value.tobytes()).hexdigest()
        if key not in self.known:
            name = f"a_{key[:24]}.fits"
            header = fits.Header()
            header["GENERATE"] = ("AOS export_optical_paths.py", "Generated numerical optical array")
            header["SOURCE"] = ("model.yaml provenance", "Frozen source and compiler identity")
            hdus = [fits.PrimaryHDU(value.real, header=header)]
            if np.iscomplexobj(value):
                hdus.append(fits.ImageHDU(value.imag, name="IMAG"))
            destination = self.output / name
            if destination.exists():
                with fits.open(destination) as existing:
                    actual = existing[0].data
                    if len(existing) == 2:
                        actual = actual + 1j * existing[1].data
                    if not np.array_equal(actual, value):
                        raise ValueError(f"Existing content-addressed array differs: {name}")
            else:
                fits.HDUList(hdus).writeto(destination)
            self.known[key] = "arrays/" + name
            self.sha256["arrays/" + name] = hashlib.sha256((self.output / name).read_bytes()).hexdigest()
        return self.known[key]


def coordinates(n, dx=1.0):
    return (np.arange(n) - n // 2) * dx


def r2(n, dx=1.0):
    x = coordinates(n, dx)
    return x[:, None] ** 2 + x[None, :] ** 2


def fourier_kernel(n, dx_in, dx_out):
    frequency = np.fft.fftshift(np.fft.fftfreq(n, dx_in))
    return np.exp(2j * np.pi * np.outer(coordinates(n, dx_out), frequency))


def interpolation_kernel(n, source_size, pitch, dx):
    """The exact periodic n-point DFT interpolant of a centered zero-padded map.

    Dirichlet evaluation avoids materializing an n-by-n Fourier matrix product.
    The phase accounts for the unpaired -Nyquist bin of an even n-point FFT.
    """
    delta = coordinates(n, dx / pitch)[:, None] - coordinates(source_size)[None, :]
    delta = (delta + n / 2) % n - n / 2
    phase = np.exp(-1j * np.pi * delta / n) if n % 2 == 0 else 1.0
    return phase * np.sinc(delta) / np.sinc(delta / n)


class Compiler:
    def __init__(self, sim, arrays):
        self.sim, self.arrays = sim, arrays
        self.n = sim.knobs["grid"]
        self.base = sim.screens()

    def mult(self, out, value):
        value = np.asarray(value)
        if value.ndim == 0:
            value = np.full((self.n, self.n), value)
        if np.all(value == 1):
            return
        # Consecutive fixed pointwise operations can be fused without changing order.
        if out and out[-1]["kind"] == "mult":
            prior = out.pop()
            value = prior["_value"] * value
        out.append(dict(kind="mult", _value=value))

    def lower(self, skeleton, parameters, anchor):
        from spectre.tape import _apply, _merged, screen_phasor
        from spectre.propagators import edge_window
        import jax.numpy as jnp
        out = []
        for (kind, static, slot), stacked in zip(skeleton, parameters):
            p = {key: np.asarray(value)[anchor] for key, value in stacked.items()}
            p = _merged(kind, static, slot, p, self.sim.adjust)
            if kind == "mark":
                continue
            if kind == "normalize":
                out.append(dict(kind="normalize"))
            elif kind in ("mult", "scale", "tilt", "lens", "aperture", "opd", "defocus", "pattern", "reflect_pattern", "pow"):
                f, a = _apply(kind, static, p, jnp.ones((self.n, self.n), dtype=complex), 1.0, {})
                self.mult(out, np.asarray(f * a))
            elif kind in ("tf", "ft"):
                op = dict(kind=kind, dx=float(p["dx"]), dz=float(p["dz"]), lam=float(p["lam"]))
                edge = static[2]
                if edge != 0 and (kind == "ft" or static[0]):
                    op["edge"] = self.arrays.save(edge_window(self.n, edge))
                if kind == "ft":
                    op.update(from_waist=bool(static[0]), forward=bool(static[1]),
                              gouy_real=0.0, gouy_imag=-1.0 if static[1] else 1.0)
                out.append(op)
            elif kind == "fshift":
                dx, lam, curv = float(p["dx"]), float(p["lam"]), float(p["curv"])
                dy, dx_offset = np.asarray(p["offset"], float)
                f = np.fft.fftshift(np.fft.fftfreq(self.n, dx))
                ramp = np.exp(-2j * np.pi * (dy * f[:, None] + dx_offset * f[None, :]))
                out.append(dict(kind="spectral_mult", array=self.arrays.save(ramp)))
                x = coordinates(self.n, dx)
                sphere = np.exp(1j * np.pi * curv * (dy**2 + dx_offset**2 - 2 * (dy*x[:, None]+dx_offset*x[None, :])) / lam)
                self.mult(out, sphere)
            elif kind == "regrid":
                di, do, lam = float(p["dx_in"]), float(p["dx_out"]), float(p["lam"])
                f = coordinates(self.n) / (self.n * di)
                f2 = (f[:, None]**2 + f[None, :]**2) * (2 * do / 0.9)**2
                roll = np.exp(-float(p["roll"]) * f2**20)
                rr = r2(self.n, do)
                post = np.exp(-(rr / (0.5*self.n*di*0.9)**2)**20).astype(complex)
                curv = float(p["dcurv"])
                if curv != 0:
                    post *= np.exp(1j*np.pi*rr*curv/lam)
                    resolved = lam / (2*do*abs(curv)) * 0.9
                    post *= np.exp(-(rr / resolved**2)**20)
                out.append(dict(kind="resample", kernel=self.arrays.save(fourier_kernel(self.n,di,do)),
                                roll=self.arrays.save(roll), post=self.arrays.save(post)))
            elif kind == "screen":
                name, pitch = static
                if name not in ("DM", "Subaru_Pupil"):
                    raise ValueError(f"Unexpected remaining live screen {name}")
                base = self.base[name]
                self.mult(out, screen_phasor(base, pitch, p["dx"], p["lam"], self.n))
                source_size = np.asarray(base["opd"]).shape[0]
                out.append(dict(kind="live_screen", name=name, pitch=float(pitch), dx=float(p["dx"]),
                                lam=float(p["lam"]), size=source_size,
                                kernel=self.arrays.save(interpolation_kernel(self.n,source_size,pitch,float(p["dx"])))))
            else:
                raise ValueError(f"Unsupported spectre operation: {kind}")
        for op in out:
            if op["kind"] == "mult":
                op["array"] = self.arrays.save(op.pop("_value"))
        return out


def revision(path):
    return subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"], text=True).strip()


def verify_sources(args):
    for root, expected in ((args.spectre, args.spectre_revision), (args.spiders_sim, args.spiders_sim_revision)):
        if expected is not None and revision(root) != expected:
            raise ValueError(f"Pinned source revision differs: {root}")
        if expected is not None:
            subprocess.run(["git", "-C", str(root), "diff", "--exit-code", "HEAD", "--", "src", "pyproject.toml"],
                           check=True, stdout=subprocess.DEVNULL)
    if args.artifact_manifest_sha256 is not None:
        actual = hashlib.sha256((args.model / "manifest.json").read_bytes()).hexdigest()
        if actual != args.artifact_manifest_sha256:
            raise ValueError("Frozen fitted artifact manifest differs from its pinned hash")


def export(args):
    verify_sources(args)
    sys.path[:0] = [str(args.spectre / "src"), str(args.spiders_sim / "src")]
    from spiders_sim import Detector, Spiders
    from spiders_sim import artifact
    from spiders_sim.fitted import FastTier
    from spectre.fit import lowphase_basis, pixel_mtf
    import jax.numpy as jnp
    options = args.resolved
    source_loader = artifact.load
    dropped_screens, dropped_slots = [], []

    def selected_loader(path):
        nonlocal dropped_screens, dropped_slots
        model, manifest = source_loader(path)
        model, manifest, dropped_screens, dropped_slots = configured_model(model, manifest, args.configuration)
        return model, manifest

    # Spiders has no model-object constructor. Scope the cold loader substitution
    # to this one constructor; no source files or artifact files are rewritten.
    with patch.object(artifact, "load", selected_loader):
        sim = Spiders(args.model, preset=options["construction_preset"], grid=options["grid"], mode=args.mode,
                      band=options["band"], ghosts=options["ghosts"], channels=options["channels"],
                      dm_base=args.dm_base, detector=Detector() if options["mtf"] else None)
    removed = configure_simulation(sim, args, FastTier)
    args.output.mkdir(parents=True, exist_ok=True)
    arrays = Arrays(args.output / "arrays")
    compiler = Compiler(sim, arrays)
    n, a, ne = 2*sim.model.window+1, len(sim.wavelengths), sim.tier.n_env
    root = dict(schema_version=1, grid=sim.knobs["grid"], image_size=n, crop_start=sim.fwd.crop[0],
                mode=args.mode, configuration=args.configuration, dm_base=args.dm_base,
                preset=args.preset, ghosts=bool(sim.knobs["ghosts"]), lod=float(sim.lod),
                dx_model=float(sim.tier.dx_model), dx_output=float(sim.tier.dx), centre=float(sim.tier.centre),
                star_peak=float(sim.star_peak), wavelengths=sim.wavelengths.tolist(), channels=sim.channels.tolist(),
                labels=list(sim.labels), n_env=ne, array_orientation="yx", complex_fits="primary_real_extension_imag")
    propagations = []
    branch_indices = {name: i+1 for i,name in enumerate(sim.fwd.branches)}
    for k, tape in enumerate(sim.fwd.tapes):
        for j in range(a):
            prefix = compiler.lower(tape.ops,tape.params,j)
            tails = {name: compiler.lower(tape.tail_ops[branch_indices[name]],tape.tail_params[branch_indices[name]],j)
                     for name in ("closed","only")}
            counts = {name: sum((1 if op["forward"] else -1) for op in prefix+ops if op["kind"]=="ft")
                      for name,ops in tails.items()}
            entry = np.asarray(tape.entrance["field"])[j]*np.asarray(tape.entrance["amplitude"])[j]
            propagations.append(dict(path_index=k,anchor_index=j,label=sim.fwd.labels[k],entrance=arrays.save(entry),
                                     prefix=prefix,**tails,signed_ft_count=counts))
    basis = np.asarray(lowphase_basis(n,n,sim.lod))
    lp = np.asarray(sim.tier.lowphase)
    main_phase = np.exp(1j*np.einsum("m,myx->yx",lp[0],basis))
    fast = dict(resample_kernel=arrays.save(fourier_kernel(n,sim.tier.dx_model,sim.tier.dx)),
                basis=arrays.save(sim.tier.basis),
                visibility=[arrays.save(v) for v in np.asarray(sim.tier.vis)],
                coefficient=arrays.save(np.asarray(sim.tier.gain)[None,:]*np.asarray(sim.tier.phasor)),
                coefficient_anchor=arrays.save(sim.tier.gain[None,:]*np.exp(2j*np.pi*sim.tier.opd_anchor/sim.wavelengths[:,None])),
                pinhole_gain_real=sim.pinhole_gain.real.tolist(), pinhole_gain_imag=sim.pinhole_gain.imag.tolist(),
                main_phase=arrays.save(main_phase), phantom_pad=32, pixel_pad=16,
                apply_pixel_mtf=sim.detector is not None and sim.detector.mtf)
    size = n+32
    freq = np.fft.fftshift(np.fft.fftfreq(size))
    fast["pixel_mtf"] = arrays.save(pixel_mtf(jnp.asarray(sim.manifest["pixel"]),freq[:,None],freq[None,:]))
    y,x = coordinates(n)[:,None]/sim.lod, coordinates(n)[None,:]/sim.lod
    frequency = np.fft.fftshift(np.fft.fftfreq(n))
    fast_paths=[]
    copies=np.asarray(sim.tier.copies)
    for k,label in enumerate(sim.labels):
        dy,dx=np.asarray(sim.tier.shift)[k] if k else (0.,0.)
        ramp=np.exp(-2j*np.pi*(dy*frequency[:,None]+dx*frequency[None,:]))
        own=lp[k] if k else np.zeros(5)
        own=own+copies[k]*lp[0]
        phase=np.einsum("m,myx->yx",own,basis)
        ty,tx=np.asarray(sim.tier.tilt)[k] if k else (0.,0.)
        phases=[arrays.save(np.exp(1j*phase+2j*np.pi*(ty*y+tx*x)*sim.tier.centre/lam)) for lam in sim.wavelengths]
        path=dict(label=label,path_index=k,source_index=k,main_copy=bool(copies[k]),
                  shift_transfer=arrays.save(ramp),phase_multiplier=phases)
        if k>=ne:
            ph=sim.tier.phantoms[k-ne]
            if ph["source"] is None:
                path.update(source_index=-1,ghost_field_anchor=[arrays.save(v) for v in np.asarray(sim.tier.ghost)])
            else:
                source=sim.tier.source[k-ne]
                freq=np.fft.fftshift(np.fft.fftfreq(n+64,sim.tier.dx))
                rho=freq[:,None]**2+freq[None,:]**2
                transfers=[arrays.save(np.exp(-1j*np.pi*lam*sim.tier.defocus[k-ne]*rho)) for lam in sim.wavelengths]
                path.update(source_index=source,defocus_transfer=transfers)
        fast_paths.append(path)
    provenance=dict(spectre_revision=revision(args.spectre),spiders_sim_revision=revision(args.spiders_sim),
                    artifact_spectre_revision=sim.manifest["spectre_commit"],model_path=str(args.model.resolve()),
                    compiler_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
                    source_fit=sim.manifest["source_fit"],fit_loss=float(sim.manifest["loss"]),
                    artifact_sha256={str(f.relative_to(args.model)):hashlib.sha256(f.read_bytes()).hexdigest()
                                     for f in sorted(args.model.rglob("*")) if f.is_file()},
                    geometry=dict(angle=float(sim.geometry.angle),flip=sim.geometry.flip,
                                  offset=list(sim.geometry.offset),pitch=float(sim.geometry.pitch)))
    warnings = []
    if args.configuration != "reference":
        warnings.extend([
            "SCC candidate configuration has not been validated against SCC measurements.",
            "Near-detector stack phantoms are excluded: spectre/ROADMAP.md identifies them as FTS-camera-only additions.",
            "Physical lenses are retained. The fitted L2B screen and its nuisance slots are excluded following spectre/ROADMAP.md's SCC prescription; this approximation is not experimentally validated.",
            "Shared fitted screens and geometry remain bench-derived and are not independently qualified on sky.",
        ])
        if options["fixed_ghost"] == "off":
            warnings.append("The fixed empirical clear-pupil ghost field is excluded; a physical upstream SCC ghost remains unresolved.")
        else:
            warnings.append("The fitted fixed ghost field is an unvalidated bench clear-pupil empirical contribution; it does not respond to atmosphere.")
    if args.apodizer_gamma == "nominal":
        warnings.append("Apodizer coverage exponent adjustment is zero; the supplied apodizer map and cosine microdot pattern are retained.")
    provenance.update(configuration=args.configuration,
                      branch_scope={"reference": "fitted reference in selected entrance mode",
                                    "scc-bench": "SCC bench candidate with internal source",
                                    "scc-sky": "SCC sky candidate with Subaru entrance"}[args.configuration],
                      qualification_status="upstream fitted reference" if args.configuration == "reference"
                                           else "unvalidated SCC candidate",
                      requested_preset=args.preset, construction_preset=options["construction_preset"],
                      spectral_mode="monochromatic" if options["band"] is None else "finite band in wavenumber",
                      intensity_synthesis="channel mean of coherent path pairs with top-hat wavenumber visibility",
                      band_width=0.0 if options["band"] is None else options["band"][1],
                      band_centre_m=float(np.mean(sim.model.band)) if options["band"] is None or options["band"][0] is None
                                    else options["band"][0],
                      anchors=len(sim.wavelengths), channels=len(sim.channels), grid=options["grid"],
                      physical_ghosts=options["ghosts"], fitted_pixel_mtf=options["mtf"],
                      fixed_ghost_selection=options["fixed_ghost"],
                      fixed_ghost_active=options["ghosts"] and any(p["source"] is None for p in sim.tier.phantoms),
                      apodizer_gamma=args.apodizer_gamma, removed_fitted_paths=removed,
                      removed_fitted_screens=dropped_screens, removed_fitted_slots=dropped_slots,
                      actual_path_labels=list(sim.labels), retained_fitted_screens=list(sim.model.screens),
                      qualifications=warnings)
    live=dict(dm_size=sim._dm_size,dm_pitch=float(sim.pitch["DM"]), entrance_size=0, entrance_pitch=0.0)
    if "Subaru_Pupil" in compiler.base:
        live.update(entrance_size=int(compiler.base["Subaru_Pupil"]["opd"].shape[0]),
                    entrance_pitch=float(sim.pitch["Subaru_Pupil"]))
    root.update(live=live,fast=fast,provenance=provenance,array_sha256=arrays.sha256,
                propagations=propagations,fast_paths=fast_paths)
    if args.reference:
        root["references"]=write_references(sim,args.output)
    metadata = root
    document = model_document(metadata, args.output, include_references=args.reference)
    document["provenance"]["compiler_source"] = "SPIDERSSim b479871 tools/export_physical15.py; generic AOS YAML export"
    document["provenance"]["array_generation_policy"] = "Frozen content-addressed FITS bytes are reused after exact numerical verification; new arrays carry GENERATE and SOURCE cards"
    write_model(document, args.output / args.manifest_name)
    print(f"Exported {len(propagations)} path/anchor propagations, {len(arrays.known)} deduplicated arrays to {args.output}")


def write_references(sim,output):
    """Deterministic command and OPD fixtures, plus actual Spiders.image products."""
    from astropy.io import fits
    valid=np.flatnonzero(sim.dm.valid)
    commands={"zero":np.zeros(468)}
    for label,index in zip(("poke_first","poke_middle","poke_last"),(valid[0],valid[len(valid)//2],valid[-1])):
        c=np.zeros(468); c[index]=0.01; commands[label]=c
    c=0.01*(np.sin(2*np.pi*(3*sim.dm.cols+2*sim.dm.rows)/24)+0.37*np.cos(2*np.pi*(sim.dm.cols-4*sim.dm.rows)/24))
    c[~sim.dm.valid]=0
    commands["asymmetric_fourier"]=c
    if sim.mode == "sky":
        commands["asymmetric_entrance"] = c.copy()
    references=[]
    for label,command in commands.items():
        print(f"Reference {label}",flush=True)
        opd=np.asarray(sim.dm_opd(command))
        inputs={}
        if label == "asymmetric_entrance":
            m=sim.screens()["Subaru_Pupil"]["opd"].shape[0]
            y,x=coordinates(m)[:,None]/m,coordinates(m)[None,:]/m
            inputs=dict(atmosphere=2e-8*np.sin(2*np.pi*(2*x+3*y)),
                        atmosphere_amplitude=0.02*(1+0.4*np.cos(2*np.pi*(x-2*y))))
        image=sim.image(dm=command,**inputs)
        products={"command":command[:,None],"dm_opd":opd,"fringed":image.fringed,
                  "unfringed":image.unfringed,"reference":image.reference}
        if inputs:
            products.update(entrance_opd=inputs["atmosphere"],entrance_attenuation=inputs["atmosphere_amplitude"])
        entry=dict(label=label,command_unit="ALPAO",dm_opd_unit="m",geometry_source="provenance.geometry")
        for name,value in products.items():
            filename=f"oracle_{label}_{name}.fits"
            header = fits.Header()
            header["GENERATE"] = "AOS export_optical_paths.py"
            header["SOURCE"] = "cold compiler reference scenario"
            fits.PrimaryHDU(value, header=header).writeto(output/filename,overwrite=True)
            entry[name]=filename
            entry[name+"_sha256"]=hashlib.sha256((output/filename).read_bytes()).hexdigest()
        for branch,fields in image.fields.items():
            for k,path in enumerate(fields):
                for j,field in enumerate(path):
                    fits.HDUList([fits.PrimaryHDU(field.real,header=header),fits.ImageHDU(field.imag,name="IMAG")]).writeto(
                        output/f"oracle_{label}_{branch}_{k}_{j}.fits",overwrite=True)
        references.append(entry)
    return references


def main():
    parser = argument_parser()
    export(resolve_options(parser.parse_args(), parser))


if __name__=="__main__":
    main()

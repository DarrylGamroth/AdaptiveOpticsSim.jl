#!/usr/bin/env python3
"""Generate a deterministic OOPAO reference bundle for AdaptiveOpticsSim.jl tests.

The script can either use an existing OOPAO checkout or clone a pinned checkout
from GitHub into a writable temporary directory before import. The output
format matches `test/reference_harness.jl`.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

DEFAULT_OOPAO_REPO = "https://github.com/cheritier/OOPAO.git"
DEFAULT_OOPAO_REF = "085d5e50ace0d20fe13cc2da20129d5400166973"

Atmosphere = None
BioEdge = None
Detector = None
DeformableMirror = None
GainSensingCamera = None
LiFT = None
Pyramid = None
ShackHartmann = None
Source = None
Telescope = None
strehlMeter = None


def run_git(*args: str, cwd: Path | None = None) -> str:
    proc = subprocess.run(
        ["git", *args],
        cwd=str(cwd) if cwd is not None else None,
        check=True,
        capture_output=True,
        text=True,
    )
    return proc.stdout.strip()


def clone_oopao_checkout(repo_url: str, repo_ref: str, workspace: Path) -> tuple[Path, str]:
    checkout = workspace / "OOPAO"
    run_git("clone", "--filter=blob:none", "--no-checkout", repo_url, str(checkout))
    run_git("fetch", "--depth", "1", "origin", repo_ref, cwd=checkout)
    run_git("checkout", "--detach", "FETCH_HEAD", cwd=checkout)
    commit = run_git("rev-parse", "HEAD", cwd=checkout)
    return checkout, commit


def prepare_oopao_checkout(args: argparse.Namespace) -> tuple[Path, str, tempfile.TemporaryDirectory[str] | None]:
    if args.oopao_path is not None:
        checkout = args.oopao_path.resolve()
        commit = run_git("rev-parse", "HEAD", cwd=checkout)
        return checkout, commit, None
    workspace = tempfile.TemporaryDirectory(prefix="oopao-checkout-")
    checkout, commit = clone_oopao_checkout(args.oopao_repo, args.oopao_ref, Path(workspace.name))
    return checkout, commit, workspace


def bootstrap_oopao(checkout: Path) -> None:
    global Atmosphere, BioEdge, Detector, GainSensingCamera, LiFT
    global Pyramid, ShackHartmann, Source, Telescope, strehlMeter, DeformableMirror

    sys.path.insert(0, str(checkout))
    from OOPAO.Atmosphere import Atmosphere as _Atmosphere
    from OOPAO.BioEdge import BioEdge as _BioEdge
    from OOPAO.DeformableMirror import DeformableMirror as _DeformableMirror
    from OOPAO.Detector import Detector as _Detector
    from OOPAO.GainSensingCamera import GainSensingCamera as _GainSensingCamera
    from OOPAO.LiFT import LiFT as _LiFT
    from OOPAO.Pyramid import Pyramid as _Pyramid
    from OOPAO.ShackHartmann import ShackHartmann as _ShackHartmann
    from OOPAO.Source import Source as _Source
    from OOPAO.Telescope import Telescope as _Telescope
    from OOPAO.tools.tools import strehlMeter as _strehlMeter

    Atmosphere = _Atmosphere
    BioEdge = _BioEdge
    DeformableMirror = _DeformableMirror
    Detector = _Detector
    GainSensingCamera = _GainSensingCamera
    LiFT = _LiFT
    Pyramid = _Pyramid
    ShackHartmann = _ShackHartmann
    Source = _Source
    Telescope = _Telescope
    strehlMeter = _strehlMeter


def make_telescope(*, resolution: int, diameter: float, sampling_time: float, central_obstruction: float, fov_arcsec: float = 0.0) -> Telescope:
    return Telescope(
        resolution=resolution,
        diameter=diameter,
        samplingTime=sampling_time,
        centralObstruction=central_obstruction,
        display_optical_path=False,
        fov=fov_arcsec,
    )


def make_source(*, band: str, magnitude: float, coordinates: tuple[float, float] = (0.0, 0.0)) -> Source:
    return Source(optBand=band, magnitude=magnitude, coordinates=list(coordinates))


def make_atmosphere(tel: Telescope) -> Atmosphere:
    atm = Atmosphere(
        telescope=tel,
        r0=0.15,
        L0=25,
        fractionalR0=[0.45, 0.1, 0.1, 0.25, 0.1],
        windSpeed=[10, 12, 11, 15, 20],
        windDirection=[0, 72, 144, 216, 288],
        altitude=[0, 1000, 5000, 10000, 12000],
    )
    atm.initializeAtmosphere(tel)
    atm.generateNewPhaseScreen(seed=10)
    return atm


def apply_ramp_opd(tel: Telescope, *, scale_x: float, scale_y: float, bias: float = 0.0) -> None:
    x = np.arange(tel.resolution, dtype=float)[:, None]
    y = np.arange(tel.resolution, dtype=float)[None, :]
    tel.OPD_no_pupil = bias + scale_x * x + scale_y * y


def write_array(path: Path, data: np.ndarray) -> None:
    np.savetxt(path, np.asarray(data, dtype=np.float64).reshape(-1), fmt="%.18e")


def write_manifest(path: Path, cases: dict[str, dict], metadata: dict[str, str] | None = None) -> None:
    lines = ["version = 1", ""]
    if metadata:
        lines.append("[metadata.oopao]")
        for key, value in metadata.items():
            lines.append(f'{key} = "{value}"')
        lines.append("")
    for case_id, case in cases.items():
        lines.extend(
            [
                f"[cases.{case_id}]",
                f'kind = "{case["kind"]}"',
                f'data = "{case["data"]}"',
                f"shape = [{', '.join(str(x) for x in case['shape'])}]",
            ]
        )
        if "storage_order" in case:
            lines.append(f'storage_order = "{case["storage_order"]}"')
        lines.extend(
            [
                f"atol = {case['atol']:.16e}",
                f"rtol = {case['rtol']:.16e}",
                "",
            ]
        )
        for section in (
            "telescope",
            "source",
            "science_source",
            "opd",
            "wfs",
            "basis",
            "controllable_optic",
            "detector",
            "compute",
            "compare",
        ):
            if section not in case:
                continue
            lines.append(f"[cases.{case_id}.{section}]")
            write_toml_mapping(lines, f"cases.{case_id}.{section}", case[section])
    path.write_text("\n".join(lines), encoding="utf-8")


def write_toml_mapping(lines, prefix, mapping):
    nested_tables = []
    array_tables = []
    for key, value in mapping.items():
        if isinstance(value, dict):
            nested_tables.append((key, value))
        elif isinstance(value, (list, tuple)) and all(isinstance(v, dict) for v in value):
            array_tables.append((key, value))
        else:
            lines.append(f"{key} = {toml_literal(value)}")
    lines.append("")
    for key, value in nested_tables:
        lines.append(f"[{prefix}.{key}]")
        write_toml_mapping(lines, f"{prefix}.{key}", value)
    for key, values in array_tables:
        for value in values:
            lines.append(f"[[{prefix}.{key}]]")
            write_toml_mapping(lines, f"{prefix}.{key}", value)


def toml_literal(value):
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return f'"{value}"'
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.16e}"
    if isinstance(value, (list, tuple)):
        return "[" + ", ".join(toml_literal(v) for v in value) + "]"
    raise TypeError(f"unsupported TOML value: {value!r}")


def cartesian_polynomial_basis(tel: Telescope, n_modes: int = 4) -> np.ndarray:
    x = np.linspace(-1.0, 1.0, tel.resolution, endpoint=False)
    y = np.linspace(-1.0, 1.0, tel.resolution, endpoint=False)
    xx, yy = np.meshgrid(x, y, indexing="ij")
    pupil = np.asarray(tel.pupil, dtype=float)
    basis = np.zeros((tel.resolution, tel.resolution, n_modes), dtype=np.float64)
    if n_modes >= 1:
        basis[:, :, 0] = pupil * xx
    if n_modes >= 2:
        basis[:, :, 1] = pupil * yy
    if n_modes >= 3:
        basis[:, :, 2] = pupil * xx * yy
    if n_modes >= 4:
        basis[:, :, 3] = pupil * (xx**2 - yy**2)
    return basis


def combine_modes(basis: np.ndarray, coeffs: np.ndarray) -> np.ndarray:
    out = np.zeros(basis.shape[:2], dtype=np.float64)
    for idx in range(min(basis.shape[2], coeffs.size)):
        out += float(coeffs[idx]) * basis[:, :, idx]
    return out


def pupil_rms_nm(opd: np.ndarray, pupil: np.ndarray) -> float:
    values = np.asarray(opd, dtype=np.float64)[np.asarray(pupil) > 0]
    return float(1e9 * np.sqrt(np.mean(values**2)))


def make_reference_wfs(kind: str, tel: Telescope):
    if kind == "shack_hartmann_slopes":
        return ShackHartmann(
            nSubap=4,
            telescope=tel,
            lightRatio=0.0,
            is_geometric=False,
            pixel_scale=0.06,
            n_pixel_per_subaperture=8,
        )
    if kind == "pyramid_slopes":
        return Pyramid(
            nSubap=4,
            telescope=tel,
            lightRatio=0.0,
            modulation=1.0,
            n_pix_separation=4,
            n_pix_edge=0,
            binning=2,
            psfCentering=True,
            postProcessing="slopesMaps",
        )
    if kind == "bioedge_slopes":
        return BioEdge(
            nSubap=4,
            telescope=tel,
            modulation=1.0,
            grey_width=0.0,
            grey_length=False,
            lightRatio=0.0,
            n_pix_separation=4,
            n_pix_edge=2,
            binning=2,
            postProcessing="slopesMaps",
        )
    raise ValueError(f"unsupported WFS kind {kind!r}")


def measure_signal(kind: str, wfs, tel: Telescope, src: Source) -> np.ndarray:
    tel * wfs
    if kind == "shack_hartmann_slopes":
        return np.asarray(wfs.signal, dtype=np.float64).copy()
    return np.asarray(wfs.signal, dtype=np.float64).reshape(-1).copy()


def reference_interaction_matrix(kind: str, wfs, tel: Telescope, src: Source, basis: np.ndarray, amplitude: float) -> np.ndarray:
    n_modes = basis.shape[2]
    mat = np.zeros((measure_signal(kind, wfs, tel, src).size, n_modes), dtype=np.float64)
    for idx in range(n_modes):
        tel.resetOPD()
        tel.OPD_no_pupil = amplitude * basis[:, :, idx]
        mat[:, idx] = measure_signal(kind, wfs, tel, src)
    tel.resetOPD()
    return mat


def psf_case(root: Path) -> dict:
    tel = make_telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.2)
    src = make_source(band="I", magnitude=0.0)
    src ** tel
    tel.computePSF(zeroPaddingFactor=2)
    data = np.asarray(tel.PSF, dtype=np.float64)
    rel = "psf_baseline.txt"
    write_array(root / rel, data)
    return {
        "kind": "psf",
        "data": rel,
        "shape": list(data.shape),
        "storage_order": "C",
        "atol": 1e-8,
        "rtol": 1e-8,
        "telescope": {
            "resolution": 32,
            "diameter": 8.0,
            "sampling_time": 1e-3,
            "central_obstruction": 0.2,
        },
        "source": {
            "kind": "ngs",
            "band": "I",
            "magnitude": 0.0,
        },
        "compute": {
            "zero_padding": 2,
        },
    }


def sh_geometric_case(root: Path, *, case_id: str, scale_x: float, scale_y: float) -> dict:
    tel = make_telescope(resolution=24, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = make_source(band="I", magnitude=0.0)
    src ** tel
    wfs = ShackHartmann(nSubap=4, telescope=tel, lightRatio=0.0, is_geometric=True)
    src ** tel
    apply_ramp_opd(tel, scale_x=scale_x, scale_y=scale_y)
    tel * wfs
    data = np.asarray(wfs.signal_2D, dtype=np.float64).reshape(-1)
    rel = f"{case_id}.txt"
    write_array(root / rel, data)
    return {
        "kind": "shack_hartmann_slopes",
        "data": rel,
        "shape": [int(data.size)],
        "atol": 1e-12,
        "rtol": 1e-12,
        "telescope": {
            "resolution": 24,
            "diameter": 8.0,
            "sampling_time": 1e-3,
            "central_obstruction": 0.0,
        },
        "source": {
            "kind": "ngs",
            "band": "I",
            "magnitude": 0.0,
        },
        "opd": {
            "kind": "ramp",
            "scale_x": scale_x,
            "scale_y": scale_y,
            "bias": 0.0,
        },
        "wfs": {
            "n_subap": 4,
            "threshold": 0.0,
            "mode": "geometric",
        },
        "compare": {
            "swap_halves": True,
            "scale": 7.594936708860756e6,
        },
    }


def sh_diffractive_case(root: Path) -> dict:
    tel = make_telescope(resolution=24, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = make_source(band="I", magnitude=0.0)
    src ** tel
    wfs = ShackHartmann(
        nSubap=4,
        telescope=tel,
        lightRatio=0.0,
        is_geometric=False,
        pixel_scale=0.06,
        n_pixel_per_subaperture=8,
    )
    src ** tel
    apply_ramp_opd(tel, scale_x=5e-9, scale_y=-2e-9)
    tel * wfs
    data = np.asarray(wfs.signal_2D, dtype=np.float64).reshape(-1)
    rel = "shack_hartmann_diffractive_ramp.txt"
    write_array(root / rel, data)
    return {
        "kind": "shack_hartmann_slopes",
        "data": rel,
        "shape": [int(data.size)],
        "atol": 2e-4,
        "rtol": 1e-2,
        "telescope": {
            "resolution": 24,
            "diameter": 8.0,
            "sampling_time": 1e-3,
            "central_obstruction": 0.0,
        },
        "source": {
            "kind": "ngs",
            "band": "I",
            "magnitude": 0.0,
        },
        "opd": {
            "kind": "ramp",
            "scale_x": 5e-9,
            "scale_y": -2e-9,
            "bias": 0.0,
        },
        "wfs": {
            "n_subap": 4,
            "threshold": 0.0,
            "mode": "diffractive",
            "pixel_scale": 0.06,
            "n_pix_subap": 8,
        },
    }


def pyramid_case(root: Path) -> dict:
    tel = make_telescope(resolution=24, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = make_source(band="I", magnitude=0.0)
    src ** tel
    wfs = Pyramid(
        nSubap=4,
        telescope=tel,
        modulation=1.0,
        lightRatio=0.0,
        postProcessing="slopesMaps",
        psfCentering=True,
        n_pix_separation=4,
        binning=2,
    )
    src ** tel
    apply_ramp_opd(tel, scale_x=5e-9, scale_y=-2e-9)
    tel * wfs
    data = np.asarray(wfs.signal_2D, dtype=np.float64).reshape(-1)
    rel = "pyramid_diffractive_ramp.txt"
    write_array(root / rel, data)
    return {
        "kind": "pyramid_slopes",
        "data": rel,
        "shape": [int(data.size)],
        "atol": 2e-3,
        "rtol": 5e-2,
        "telescope": {
            "resolution": 24,
            "diameter": 8.0,
            "sampling_time": 1e-3,
            "central_obstruction": 0.0,
        },
        "source": {
            "kind": "ngs",
            "band": "I",
            "magnitude": 0.0,
        },
        "opd": {
            "kind": "ramp",
            "scale_x": 5e-9,
            "scale_y": -2e-9,
            "bias": 0.0,
        },
        "wfs": {
            "n_subap": 4,
            "threshold": 0.0,
            "mode": "diffractive",
            "modulation": 1.0,
            "modulation_points": int(wfs.nTheta),
            "diffraction_padding": 2,
            "psf_centering": True,
            "n_pix_separation": 4,
            "binning": 2,
        },
    }


def bioedge_case(root: Path) -> dict:
    tel = make_telescope(resolution=24, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = make_source(band="I", magnitude=0.0)
    src ** tel
    wfs = BioEdge(
        nSubap=4,
        telescope=tel,
        modulation=1.0,
        lightRatio=0.0,
        postProcessing="slopesMaps",
        n_pix_separation=4,
        binning=2,
    )
    src ** tel
    apply_ramp_opd(tel, scale_x=5e-9, scale_y=-2e-9)
    tel * wfs
    data = np.asarray(wfs.signal_2D, dtype=np.float64).reshape(-1)
    rel = "bioedge_diffractive_ramp.txt"
    write_array(root / rel, data)
    return {
        "kind": "bioedge_slopes",
        "data": rel,
        "shape": [int(data.size)],
        "atol": 2e-3,
        "rtol": 5e-2,
        "telescope": {
            "resolution": 24,
            "diameter": 8.0,
            "sampling_time": 1e-3,
            "central_obstruction": 0.0,
        },
        "source": {
            "kind": "ngs",
            "band": "I",
            "magnitude": 0.0,
        },
        "opd": {
            "kind": "ramp",
            "scale_x": 5e-9,
            "scale_y": -2e-9,
            "bias": 0.0,
        },
        "wfs": {
            "n_subap": 4,
            "threshold": 0.0,
            "mode": "diffractive",
            "modulation": 1.0,
            "modulation_points": int(wfs.nTheta),
            "diffraction_padding": 2,
            "psf_centering": True,
            "n_pix_separation": 4,
            "binning": 2,
        },
    }


def modal_tiptilt_case(root: Path, *, case_id: str, kind: str, mode_index: int, amplitude: float) -> dict:
    tel = make_telescope(resolution=24, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = make_source(band="I", magnitude=0.0)
    src ** tel
    wfs = make_reference_wfs(kind, tel)
    src ** tel
    basis = cartesian_polynomial_basis(tel, n_modes=2)
    tel.resetOPD()
    tel.OPD_no_pupil = amplitude * basis[:, :, mode_index - 1]
    data = measure_signal(kind, wfs, tel, src)
    rel = f"{case_id}.txt"
    write_array(root / rel, data)
    return {
        "kind": kind,
        "data": rel,
        "shape": [int(data.size)],
        "atol": 4e-4 if kind == "shack_hartmann_slopes" else 2e-3,
        "rtol": 5e-2 if kind == "shack_hartmann_slopes" else 5e-2,
        "telescope": {
            "resolution": 24,
            "diameter": 8.0,
            "sampling_time": 1e-3,
            "central_obstruction": 0.0,
        },
        "source": {
            "kind": "ngs",
            "band": "I",
            "magnitude": 0.0,
        },
        "basis": {
            "kind": "cartesian_polynomials",
            "n_modes": 2,
        },
        "opd": {
            "kind": "basis_mode",
            "mode_index": mode_index,
            "amplitude": amplitude,
        },
        "controllable_optic": {
            "kind": "modal",
            "basis": "cartesian_tilt",
            "label": "tiptilt",
            "scale": 1.0,
            "command": [amplitude, 0.0] if mode_index == 1 else [0.0, amplitude],
        },
        "wfs": {
            "n_subap": 4,
            "threshold": 0.0,
            "mode": "diffractive",
            **(
                {
                    "pixel_scale": 0.06,
                    "n_pix_subap": 8,
                }
                if kind == "shack_hartmann_slopes"
                else {
                    "modulation": 1.0,
                    "modulation_points": 8,
                    "diffraction_padding": 2,
                    "psf_centering": True,
                    "n_pix_separation": 4,
                    "binning": 2,
                    **({"n_pix_edge": 2} if kind == "bioedge_slopes" else {}),
                }
            ),
        },
    }


def julia_dm_mechanical_coupling(n_act: int, influence_width: float) -> float:
    pitch_norm = 2.0 / (n_act - 1)
    return float(np.exp(-(pitch_norm ** 2) / (2.0 * influence_width ** 2)))


def julia_dm_coordinates(tel: Telescope, n_act: int) -> np.ndarray:
    coords = np.linspace(-tel.D / 2.0, tel.D / 2.0, n_act, dtype=np.float64)
    pairs = [(x, y) for x in coords for y in coords]
    return np.asarray(pairs, dtype=np.float64)


def composite_tiptilt_dm_case(root: Path, *, case_id: str, kind: str, tip_amplitude: float,
    dm_command: np.ndarray, n_act: int = 4, influence_width: float = 0.3) -> dict:
    tel = make_telescope(resolution=24, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = make_source(band="I", magnitude=0.0)
    src ** tel
    wfs = make_reference_wfs(kind, tel)
    src ** tel
    basis = cartesian_polynomial_basis(tel, n_modes=2)
    dm = DeformableMirror(
        telescope=tel,
        nSubap=n_act - 1,
        coordinates=julia_dm_coordinates(tel, n_act),
        pitch=tel.D / (n_act - 1),
        mechCoupling=julia_dm_mechanical_coupling(n_act, influence_width),
    )
    dm.coefs = np.asarray(dm_command, dtype=np.float64)
    combined_opd = tip_amplitude * basis[:, :, 0] + np.asarray(dm.OPD, dtype=np.float64)
    tel.resetOPD()
    tel.OPD_no_pupil = combined_opd
    data = measure_signal(kind, wfs, tel, src)
    rel = f"{case_id}.txt"
    write_array(root / rel, data)
    if kind == "shack_hartmann_slopes":
        atol = 4e-4
        rtol = 4e-2
    elif kind == "pyramid_slopes":
        atol = 4e-4
        rtol = 3e-2
    else:
        atol = 1e-3
        rtol = 4e-2
    return {
        "kind": kind,
        "data": rel,
        "shape": [int(data.size)],
        "atol": atol,
        "rtol": rtol,
        "telescope": {
            "resolution": 24,
            "diameter": 8.0,
            "sampling_time": 1e-3,
            "central_obstruction": 0.0,
        },
        "source": {
            "kind": "ngs",
            "band": "I",
            "magnitude": 0.0,
        },
        "controllable_optic": {
            "kind": "composite",
            "command": np.concatenate((np.array([tip_amplitude, 0.0]), np.asarray(dm_command, dtype=np.float64))).tolist(),
            "components": [
                {
                    "kind": "modal",
                    "basis": "cartesian_tilt",
                    "label": "tiptilt",
                    "scale": 1.0,
                },
                {
                    "kind": "dm",
                    "label": "dm",
                    "n_act": n_act,
                    "influence_width": influence_width,
                },
            ],
        },
        "wfs": {
            "n_subap": 4,
            "threshold": 0.0,
            "mode": "diffractive",
            **(
                {
                    "pixel_scale": 0.06,
                    "n_pix_subap": 8,
                }
                if kind == "shack_hartmann_slopes"
                else {
                    "modulation": 1.0,
                    "modulation_points": 8,
                    "diffraction_padding": 2,
                    "psf_centering": True,
                    "n_pix_separation": 4,
                    "binning": 2,
                    **({"n_pix_edge": 2} if kind == "bioedge_slopes" else {}),
                }
            ),
        },
    }


def gsc_case(root: Path) -> dict:
    tel = make_telescope(resolution=24, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = make_source(band="R", magnitude=8.0)
    src ** tel
    wfs = Pyramid(
        nSubap=4,
        telescope=tel,
        modulation=3.0,
        lightRatio=0.0,
        postProcessing="slopesMaps",
        psfCentering=True,
        n_pix_separation=2,
        n_pix_edge=1,
        binning=1,
    )
    basis = cartesian_polynomial_basis(tel, 4)
    gsc = GainSensingCamera(wfs.mask, basis)
    src.OPD_no_pupil = np.zeros((tel.resolution, tel.resolution), dtype=np.float64)
    tel * wfs
    wfs * wfs.focal_plane_camera
    wfs.focal_plane_camera * gsc
    amplitude = 5e-8
    mode_index = 1
    src.OPD_no_pupil = amplitude * basis[:, :, mode_index - 1]
    tel * wfs
    wfs * wfs.focal_plane_camera
    wfs.focal_plane_camera * gsc
    data = np.asarray(gsc.og, dtype=np.float64).reshape(-1)
    rel = "gain_sensing_camera_optical_gains.txt"
    write_array(root / rel, data)
    return {
        "kind": "gsc_optical_gains",
        "data": rel,
        "shape": [int(data.size)],
        "atol": 5e-10,
        "rtol": 5e-10,
        "telescope": {
            "resolution": 24,
            "diameter": 8.0,
            "sampling_time": 1e-3,
            "central_obstruction": 0.0,
        },
        "source": {
            "kind": "ngs",
            "band": "R",
            "magnitude": 8.0,
        },
        "opd": {
            "kind": "basis_mode",
            "mode_index": mode_index,
            "amplitude": amplitude,
        },
        "wfs": {
            "n_subap": 4,
            "threshold": 0.0,
            "mode": "diffractive",
            "modulation": 3.0,
            "modulation_points": int(wfs.nTheta),
            "diffraction_padding": 2,
            "psf_centering": True,
            "n_pix_separation": 2,
            "n_pix_edge": 1,
            "binning": 1,
        },
        "basis": {
            "kind": "cartesian_polynomials",
            "n_modes": 4,
        },
    }


def gsc_branch_step_cases(root: Path) -> dict[str, dict]:
    tel = make_telescope(
        resolution=16,
        diameter=8.0,
        sampling_time=1e-3,
        central_obstruction=0.0,
        fov_arcsec=1.0,
    )
    ngs = make_source(band="R", magnitude=8.0)
    ngs ** tel
    basis = cartesian_polynomial_basis(tel, 4)
    wfs = Pyramid(
        nSubap=4,
        telescope=tel,
        lightRatio=0.5,
        modulation=3.0,
        binning=1,
        n_pix_separation=2,
        n_pix_edge=1,
        postProcessing="slopesMaps_incidence_flux",
    )
    gsc = GainSensingCamera(mask=wfs.mask, basis=basis)
    ngs ** tel * wfs
    wfs.focal_plane_camera.resolution = wfs.nRes
    wfs * wfs.focal_plane_camera
    wfs.focal_plane_camera * gsc

    H = reference_interaction_matrix("pyramid_slopes", wfs, tel, ngs, basis, 1e-9)
    recon = np.linalg.pinv(H)
    gain = 0.2
    og_floor = 0.05
    branch_iteration = 4

    atm = make_atmosphere(tel)
    forcing_ngs = np.zeros((tel.resolution, tel.resolution, 6), dtype=np.float64)
    for idx in range(forcing_ngs.shape[2]):
        atm.update()
        atm * ngs * tel
        forcing_ngs[:, :, idx] = np.asarray(tel.OPD, dtype=np.float64).copy()

    control_coeffs = np.zeros(basis.shape[2], dtype=np.float64)
    delayed_signal = np.zeros(H.shape[0], dtype=np.float64)
    for idx in range(branch_iteration - 1):
        correction_opd = combine_modes(basis, control_coeffs)
        residual_ngs = forcing_ngs[:, :, idx] - correction_opd
        ngs ** tel
        tel.OPD_no_pupil = residual_ngs
        tel * wfs
        signal = np.asarray(wfs.signal, dtype=np.float64).reshape(-1)
        wfs * wfs.focal_plane_camera
        wfs.focal_plane_camera * gsc
        og_safe = np.maximum(np.abs(np.asarray(gsc.og, dtype=np.float64)), og_floor)
        control_coeffs += gain * ((recon @ delayed_signal) / og_safe)
        delayed_signal[:] = signal

    correction_opd = combine_modes(basis, control_coeffs)
    residual_ngs = forcing_ngs[:, :, branch_iteration - 1] - correction_opd
    ngs ** tel
    tel.OPD_no_pupil = residual_ngs
    tel * wfs
    signal = np.asarray(wfs.signal, dtype=np.float64).reshape(-1)
    wfs * wfs.focal_plane_camera
    frame = np.asarray(wfs.focal_plane_camera.frame, dtype=np.float64).copy()
    wfs.focal_plane_camera * gsc
    optical_gains = np.asarray(gsc.og, dtype=np.float64).reshape(-1)

    residual_rel = "gsc_branch_step_residual_opd.txt"
    frame_rel = "gsc_branch_step_modulation_frame.txt"
    og_rel = "gsc_branch_step_optical_gains.txt"
    signal_rel = "gsc_branch_step_signal.txt"
    write_array(root / residual_rel, residual_ngs)
    write_array(root / frame_rel, frame)
    write_array(root / og_rel, optical_gains)
    write_array(root / signal_rel, signal)

    common = {
        "telescope": {
            "resolution": 16,
            "diameter": 8.0,
            "sampling_time": 1e-3,
            "central_obstruction": 0.0,
            "fov_arcsec": 1.0,
        },
        "source": {
            "kind": "ngs",
            "band": "R",
            "magnitude": 8.0,
        },
        "basis": {
            "kind": "cartesian_polynomials",
            "n_modes": 4,
        },
        "wfs": {
            "kind": "pyramid_slopes",
            "n_subap": 4,
            "mode": "diffractive",
            "threshold": 0.5,
            "light_ratio": 0.5,
            "normalization": "incidence_flux",
            "modulation": 3.0,
            "n_pix_separation": 2,
            "n_pix_edge": 1,
            "binning": 1,
        },
        "compute": {
            "calibration_amplitude": 1e-9,
            "branch_iteration": branch_iteration,
            "control_coefficients": control_coeffs.tolist(),
            "residual_opd_data": residual_rel,
            "residual_shape": list(residual_ngs.shape),
            "residual_storage_order": "C",
        },
    }

    return {
        "gsc_branch_step_modulation_frame": {
            "kind": "gsc_modulation_frame",
            "data": frame_rel,
            "shape": list(frame.shape),
            "storage_order": "C",
            "atol": 5e-8,
            "rtol": 1e-6,
            **common,
        },
        "gsc_branch_step_optical_gains": {
            "kind": "gsc_optical_gains",
            "data": og_rel,
            "shape": [int(optical_gains.size)],
            "atol": 5e-9,
            "rtol": 5e-7,
            **common,
        },
        "gsc_branch_step_signal": {
            "kind": "pyramid_slopes",
            "data": signal_rel,
            "shape": [int(signal.size)],
            "atol": 5e-9,
            "rtol": 5e-7,
            **common,
        },
    }


def transfer_function_case(root: Path) -> dict:
    fs = 300.0
    n_freq = 512
    loop_gains = np.array([0.2, 0.6], dtype=np.float64)
    freq = np.linspace(fs / n_freq, fs / 2, n_freq - 1, dtype=np.float64)
    Ti = 1.0 / fs
    Tau = Ti / 2
    Tdm = Ti / 2
    S = 1j * 2.0 * np.pi * freq
    H_WFS = np.exp(-Ti / 2 * S)
    H_RTC = np.exp(-Tau * S)
    H_DM = np.exp(-Tdm * S)
    H_DAC = (1.0 - np.exp(-S * Ti)) / (S * Ti)
    out = np.zeros((freq.size, loop_gains.size), dtype=np.float64)
    for idx, gain in enumerate(loop_gains):
        CC = gain / (1.0 - np.exp(-S * Ti))
        H_OL = H_WFS * H_RTC * H_DAC * H_DM * CC
        H_ER = 1.0 / (1.0 + H_OL)
        out[:, idx] = 20.0 * np.log10(np.abs(H_ER))
    rel = "transfer_function_rejection.txt"
    write_array(root / rel, out)
    return {
        "kind": "transfer_function_rejection",
        "data": rel,
        "shape": list(out.shape),
        "storage_order": "C",
        "atol": 1e-12,
        "rtol": 1e-12,
        "compute": {
            "fs": fs,
            "n_freq": n_freq,
            "loop_gains": loop_gains.tolist(),
        },
    }


def lift_case(root: Path) -> dict:
    tel = make_telescope(resolution=24, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = make_source(band="I", magnitude=8.0)
    src ** tel
    basis = cartesian_polynomial_basis(tel, 4)
    det = Detector(photonNoise=False, readoutNoise=1.0, QE=1.0, psf_sampling=2, binning=1)
    diversity = np.zeros((tel.resolution, tel.resolution), dtype=np.float64)
    lift = LiFT(tel, basis, det, diversity_OPD=diversity, iterations=3, img_resolution=48, numerical=False)
    coeffs = np.zeros(4, dtype=np.float64)
    mode_ids = [0, 1, 2, 3]
    H = lift.generateLIFTinteractionMatrices(coeffs, mode_ids, flux_norm=1.0)
    stack = np.stack([H[:, idx].reshape((48, 48), order="C") for idx in range(H.shape[1])], axis=2)
    rel = "lift_interaction_matrix.txt"
    write_array(root / rel, stack)
    return {
        "kind": "lift_interaction_matrix",
        "data": rel,
        "shape": list(stack.shape),
        "storage_order": "C",
        "atol": 1e-4,
        "rtol": 1e-10,
        "telescope": {
            "resolution": 24,
            "diameter": 8.0,
            "sampling_time": 1e-3,
            "central_obstruction": 0.0,
        },
        "source": {
            "kind": "ngs",
            "band": "I",
            "magnitude": 8.0,
        },
        "basis": {
            "kind": "cartesian_polynomials",
            "n_modes": 4,
        },
        "detector": {
            "noise": "readout",
            "readout_sigma": 1.0,
            "integration_time": 1.0,
            "qe": 1.0,
            "psf_sampling": 2,
            "binning": 1,
        },
        "compute": {
            "img_resolution": 48,
            "mode_ids": [1, 2, 3, 4],
            "coefficients": coeffs.tolist(),
            "iterations": 3,
            "numerical": False,
            "flux_norm": 1.0,
        },
    }


def closed_loop_trace_case(root: Path, *, case_id: str, kind: str) -> dict:
    tel = make_telescope(resolution=24, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = make_source(band="I", magnitude=0.0)
    src ** tel
    basis = cartesian_polynomial_basis(tel, 4)
    wfs = make_reference_wfs(kind, tel)
    forcing_coeffs = np.array(
        [
            [2.0e-8, -1.0e-8, 0.5e-8, -0.25e-8],
            [1.5e-8, 0.5e-8, -0.75e-8, 0.5e-8],
            [-1.0e-8, 1.25e-8, 0.5e-8, -0.5e-8],
            [0.5e-8, -0.75e-8, 1.0e-8, 0.25e-8],
            [0.25e-8, 0.5e-8, -1.25e-8, 0.75e-8],
            [-0.75e-8, -0.25e-8, 0.75e-8, -1.0e-8],
        ],
        dtype=np.float64,
    )
    if kind == "bioedge_slopes":
        forcing_coeffs = forcing_coeffs[:4, :]
    gain = 0.2 if kind == "bioedge_slopes" else 0.4
    frame_delay = 2
    calibration_amplitude = 1e-9
    zero_padding = 2
    H = reference_interaction_matrix(kind, wfs, tel, src, basis, calibration_amplitude)
    recon = np.linalg.pinv(H)
    control_coeffs = np.zeros(basis.shape[2], dtype=np.float64)
    delayed_signal = np.zeros(H.shape[0], dtype=np.float64)
    tel.resetOPD()
    tel.computePSF(zeroPaddingFactor=zero_padding)
    psf_ref = np.asarray(tel.PSF, dtype=np.float64).copy()
    trace = np.zeros((forcing_coeffs.shape[0], 5), dtype=np.float64)

    for idx, coeffs in enumerate(forcing_coeffs):
        forcing_opd = combine_modes(basis, coeffs)
        correction_opd = combine_modes(basis, control_coeffs)
        residual_opd = forcing_opd - correction_opd
        tel.resetOPD()
        tel.OPD_no_pupil = residual_opd
        trace[idx, 0] = pupil_rms_nm(forcing_opd, tel.pupil)
        trace[idx, 1] = pupil_rms_nm(residual_opd, tel.pupil)
        tel.computePSF(zeroPaddingFactor=zero_padding)
        trace[idx, 2] = strehlMeter(
            np.asarray(tel.PSF, dtype=np.float64),
            tel,
            PSF_ref=psf_ref,
            zeroPaddingFactor=zero_padding,
            display=False,
        )
        signal = measure_signal(kind, wfs, tel, src)
        trace[idx, 3] = float(np.linalg.norm(signal))
        if frame_delay == 1:
            delayed_signal[:] = signal
        control_coeffs += gain * (recon @ delayed_signal)
        trace[idx, 4] = float(np.linalg.norm(control_coeffs))
        if frame_delay == 2:
            delayed_signal[:] = signal

    rel = f"{case_id}.txt"
    write_array(root / rel, trace)
    wfs_section = {"kind": kind, "n_subap": 4, "mode": "diffractive", "threshold": 0.0}
    if kind == "shack_hartmann_slopes":
        wfs_section["pixel_scale"] = 0.06
        wfs_section["n_pix_subap"] = 8
    elif kind in ("pyramid_slopes", "bioedge_slopes"):
        wfs_section["modulation"] = 1.0
        wfs_section["n_pix_separation"] = 4
        wfs_section["n_pix_edge"] = 0 if kind == "pyramid_slopes" else 2
        if kind == "bioedge_slopes":
            wfs_section["grey_width"] = 0.0
            wfs_section["grey_length"] = False
            wfs_section["binning"] = 2
    return {
        "kind": "closed_loop_trace",
        "data": rel,
        "shape": list(trace.shape),
        "storage_order": "C",
        "atol": 2e-4,
        "rtol": 1e-3,
        "telescope": {
            "resolution": 24,
            "diameter": 8.0,
            "sampling_time": 1e-3,
            "central_obstruction": 0.0,
        },
        "source": {
            "kind": "ngs",
            "band": "I",
            "magnitude": 0.0,
        },
        "basis": {
            "kind": "cartesian_polynomials",
            "n_modes": int(basis.shape[2]),
        },
        "wfs": wfs_section,
        "compute": {
            "gain": gain,
            "frame_delay": frame_delay,
            "calibration_amplitude": calibration_amplitude,
            "psf_zero_padding": zero_padding,
            "forcing_coefficients": forcing_coeffs.tolist(),
        },
    }


def gsc_closed_loop_trace_case(root: Path) -> dict:
    tel = make_telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0, fov_arcsec=1.0)
    src = make_source(band="R", magnitude=8.0)
    src ** tel
    basis = cartesian_polynomial_basis(tel, 4)
    wfs = Pyramid(
        nSubap=4,
        telescope=tel,
        lightRatio=0.5,
        modulation=3.0,
        binning=1,
        n_pix_separation=2,
        n_pix_edge=1,
        postProcessing="slopesMaps_incidence_flux",
    )
    gsc = GainSensingCamera(mask=wfs.mask, basis=basis)
    ngs = src
    ngs ** tel * wfs
    wfs.focal_plane_camera.resolution = wfs.nRes
    wfs * wfs.focal_plane_camera
    wfs.focal_plane_camera * gsc
    H = reference_interaction_matrix("pyramid_slopes", wfs, tel, src, basis, 1e-9)
    recon = np.linalg.pinv(H)
    forcing_coeffs = np.array(
        [
            [2.0e-8, -1.0e-8, 0.5e-8, -0.25e-8],
            [1.5e-8, 0.5e-8, -0.75e-8, 0.5e-8],
            [-1.0e-8, 1.25e-8, 0.5e-8, -0.5e-8],
            [0.5e-8, -0.75e-8, 1.0e-8, 0.25e-8],
        ],
        dtype=np.float64,
    )
    gain = 0.2
    frame_delay = 2
    og_floor = 0.05
    zero_padding = 2
    control_coeffs = np.zeros(basis.shape[2], dtype=np.float64)
    delayed_signal = np.zeros(H.shape[0], dtype=np.float64)
    tel.resetOPD()
    tel.computePSF(zeroPaddingFactor=zero_padding)
    psf_ref = np.asarray(tel.PSF, dtype=np.float64).copy()
    trace = np.zeros((forcing_coeffs.shape[0], 6), dtype=np.float64)

    for idx, coeffs in enumerate(forcing_coeffs):
        forcing_opd = combine_modes(basis, coeffs)
        correction_opd = combine_modes(basis, control_coeffs)
        residual_opd = forcing_opd - correction_opd
        tel.resetOPD()
        tel.OPD_no_pupil = residual_opd
        trace[idx, 0] = pupil_rms_nm(forcing_opd, tel.pupil)
        trace[idx, 1] = pupil_rms_nm(residual_opd, tel.pupil)
        tel.computePSF(zeroPaddingFactor=zero_padding)
        trace[idx, 2] = strehlMeter(np.asarray(tel.PSF, dtype=np.float64), tel, PSF_ref=psf_ref, zeroPaddingFactor=zero_padding, display=False)
        signal = measure_signal("pyramid_slopes", wfs, tel, src)
        trace[idx, 3] = float(np.linalg.norm(signal))
        wfs * wfs.focal_plane_camera
        wfs.focal_plane_camera * gsc
        og_safe = np.maximum(np.abs(np.asarray(gsc.og, dtype=np.float64)), og_floor)
        trace[idx, 4] = float(np.mean(og_safe))
        if frame_delay == 1:
            delayed_signal[:] = signal
        control_coeffs += gain * ((recon @ delayed_signal) / og_safe)
        trace[idx, 5] = float(np.linalg.norm(control_coeffs))
        if frame_delay == 2:
            delayed_signal[:] = signal

    rel = "gsc_closed_loop_trace.txt"
    write_array(root / rel, trace)
    return {
        "kind": "gsc_closed_loop_trace",
        "data": rel,
        "shape": list(trace.shape),
        "storage_order": "C",
        "atol": 2e-4,
        "rtol": 1e-3,
        "telescope": {
            "resolution": 16,
            "diameter": 8.0,
            "sampling_time": 1e-3,
            "central_obstruction": 0.0,
            "fov_arcsec": 1.0,
        },
        "source": {
            "kind": "ngs",
            "band": "R",
            "magnitude": 8.0,
        },
        "basis": {
            "kind": "cartesian_polynomials",
            "n_modes": 4,
        },
        "wfs": {
            "kind": "pyramid_slopes",
            "n_subap": 4,
            "mode": "diffractive",
            "threshold": 0.5,
            "light_ratio": 0.5,
            "normalization": "incidence_flux",
            "modulation": 3.0,
            "n_pix_separation": 2,
            "n_pix_edge": 1,
            "binning": 1,
        },
        "compute": {
            "gain": gain,
            "frame_delay": frame_delay,
            "calibration_amplitude": 1e-9,
            "psf_zero_padding": zero_padding,
            "og_floor": og_floor,
            "forcing_coefficients": forcing_coeffs.tolist(),
        },
    }


def gsc_atmosphere_replay_trace_case(root: Path, *, case_id: str = "gsc_atmosphere_replay_trace", n_iter: int = 6) -> dict:
    tel = make_telescope(
        resolution=16,
        diameter=8.0,
        sampling_time=1e-3,
        central_obstruction=0.0,
        fov_arcsec=1.0,
    )
    ngs = make_source(band="R", magnitude=8.0)
    sci = make_source(band="K", magnitude=8.0, coordinates=(0.5, 0.0))
    ngs ** tel
    basis = cartesian_polynomial_basis(tel, 4)
    wfs = Pyramid(
        nSubap=4,
        telescope=tel,
        lightRatio=0.5,
        modulation=3.0,
        binning=1,
        n_pix_separation=2,
        n_pix_edge=1,
        postProcessing="slopesMaps_incidence_flux",
    )
    gsc = GainSensingCamera(mask=wfs.mask, basis=basis)
    ngs ** tel * wfs
    wfs.focal_plane_camera.resolution = wfs.nRes
    wfs * wfs.focal_plane_camera
    wfs.focal_plane_camera * gsc

    H = reference_interaction_matrix("pyramid_slopes", wfs, tel, ngs, basis, 1e-9)
    recon = np.linalg.pinv(H)
    gain = 0.2
    frame_delay = 2
    og_floor = 0.05
    zero_padding = 2
    atm = make_atmosphere(tel)
    forcing_ngs = np.zeros((tel.resolution, tel.resolution, n_iter), dtype=np.float64)
    forcing_src = np.zeros_like(forcing_ngs)
    for idx in range(n_iter):
        atm.update()
        atm * ngs * tel
        forcing_ngs[:, :, idx] = np.asarray(tel.OPD, dtype=np.float64).copy()
        atm * sci * tel
        forcing_src[:, :, idx] = np.asarray(tel.OPD, dtype=np.float64).copy()

    tel.resetOPD()
    ngs ** tel
    tel.computePSF(zeroPaddingFactor=zero_padding)
    ngs_psf_ref = np.asarray(tel.PSF, dtype=np.float64).copy()
    tel.resetOPD()
    sci ** tel
    tel.computePSF(zeroPaddingFactor=zero_padding)
    src_psf_ref = np.asarray(tel.PSF, dtype=np.float64).copy()

    control_coeffs = np.zeros(basis.shape[2], dtype=np.float64)
    delayed_signal = np.zeros(H.shape[0], dtype=np.float64)
    trace = np.zeros((n_iter, 7), dtype=np.float64)

    for idx in range(n_iter):
        correction_opd = combine_modes(basis, control_coeffs)
        residual_ngs = forcing_ngs[:, :, idx] - correction_opd
        residual_src = forcing_src[:, :, idx] - correction_opd

        ngs ** tel
        tel.OPD_no_pupil = residual_ngs
        tel * wfs
        tel.computePSF(zeroPaddingFactor=zero_padding)
        trace[idx, 0] = pupil_rms_nm(forcing_ngs[:, :, idx], tel.pupil)
        trace[idx, 1] = pupil_rms_nm(residual_ngs, tel.pupil)
        trace[idx, 3] = strehlMeter(
            np.asarray(tel.PSF, dtype=np.float64),
            tel,
            PSF_ref=ngs_psf_ref,
            zeroPaddingFactor=zero_padding,
            display=False,
        )
        signal = np.asarray(wfs.signal, dtype=np.float64).reshape(-1)
        trace[idx, 5] = float(np.linalg.norm(signal))
        wfs * wfs.focal_plane_camera
        wfs.focal_plane_camera * gsc
        og_safe = np.maximum(np.abs(np.asarray(gsc.og, dtype=np.float64)), og_floor)

        sci ** tel
        tel.OPD_no_pupil = residual_src
        tel.computePSF(zeroPaddingFactor=zero_padding)
        trace[idx, 2] = pupil_rms_nm(residual_src, tel.pupil)
        trace[idx, 4] = strehlMeter(
            np.asarray(tel.PSF, dtype=np.float64),
            tel,
            PSF_ref=src_psf_ref,
            zeroPaddingFactor=zero_padding,
            display=False,
        )

        if frame_delay == 1:
            delayed_signal[:] = signal
        control_coeffs += gain * ((recon @ delayed_signal) / og_safe)
        trace[idx, 6] = float(np.mean(og_safe))
        if frame_delay == 2:
            delayed_signal[:] = signal

    trace_rel = f"{case_id}.txt"
    ngs_rel = f"{case_id}_ngs_opd.txt"
    src_rel = f"{case_id}_src_opd.txt"
    write_array(root / trace_rel, trace)
    write_array(root / ngs_rel, forcing_ngs)
    write_array(root / src_rel, forcing_src)
    return {
        "kind": "gsc_atmosphere_replay_trace",
        "data": trace_rel,
        "shape": list(trace.shape),
        "storage_order": "C",
        "atol": 5e-4,
        "rtol": 5e-3,
        "telescope": {
            "resolution": 16,
            "diameter": 8.0,
            "sampling_time": 1e-3,
            "central_obstruction": 0.0,
            "fov_arcsec": 1.0,
        },
        "source": {
            "kind": "ngs",
            "band": "R",
            "magnitude": 8.0,
        },
        "science_source": {
            "kind": "ngs",
            "band": "K",
            "magnitude": 8.0,
            "coordinates": [0.5, 0.0],
        },
        "basis": {
            "kind": "cartesian_polynomials",
            "n_modes": 4,
        },
        "wfs": {
            "kind": "pyramid_slopes",
            "n_subap": 4,
            "mode": "diffractive",
            "threshold": 0.5,
            "light_ratio": 0.5,
            "normalization": "incidence_flux",
            "modulation": 3.0,
            "modulation_points": int(wfs.nTheta),
            "diffraction_padding": 2,
            "psf_centering": True,
            "n_pix_separation": 2,
            "n_pix_edge": 1,
            "binning": 1,
        },
        "compute": {
            "gain": gain,
            "frame_delay": frame_delay,
            "calibration_amplitude": 1e-9,
            "psf_zero_padding": zero_padding,
            "og_floor": og_floor,
            "forcing_ngs_data": ngs_rel,
            "forcing_src_data": src_rel,
            "forcing_shape": list(forcing_ngs.shape),
            "forcing_storage_order": "C",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path, help="Bundle output directory")
    parser.add_argument(
        "--oopao-path",
        type=Path,
        default=None,
        help="use an existing writable OOPAO checkout instead of cloning",
    )
    parser.add_argument(
        "--oopao-repo",
        default=DEFAULT_OOPAO_REPO,
        help=f"OOPAO git repo URL to clone (default: {DEFAULT_OOPAO_REPO})",
    )
    parser.add_argument(
        "--oopao-ref",
        default=DEFAULT_OOPAO_REF,
        help=f"OOPAO git ref/commit to checkout when cloning (default: {DEFAULT_OOPAO_REF})",
    )
    args = parser.parse_args()

    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    checkout, commit, workspace = prepare_oopao_checkout(args)
    try:
        bootstrap_oopao(checkout)
        composite_dm_command = np.array(
            [
                0.0,
                1.5e-9,
                -0.5e-9,
                0.0,
                1.0e-9,
                -2.0e-9,
                0.75e-9,
                -0.5e-9,
                -0.75e-9,
                1.25e-9,
                -1.5e-9,
                0.5e-9,
                0.0,
                -0.5e-9,
                1.0e-9,
                0.0,
            ],
            dtype=np.float64,
        )
        cases = {
            "psf_baseline": psf_case(root),
            "shack_hartmann_geometric_ramp_xy": sh_geometric_case(
                root,
                case_id="shack_hartmann_geometric_ramp_xy",
                scale_x=5e-9,
                scale_y=-2e-9,
            ),
            "shack_hartmann_geometric_ramp_y": sh_geometric_case(
                root,
                case_id="shack_hartmann_geometric_ramp_y",
                scale_x=0.0,
                scale_y=4e-9,
            ),
            "shack_hartmann_diffractive_ramp": sh_diffractive_case(root),
            "shack_hartmann_diffractive_tip_mode": modal_tiptilt_case(
                root,
                case_id="shack_hartmann_diffractive_tip_mode",
                kind="shack_hartmann_slopes",
                mode_index=1,
                amplitude=5e-9,
            ),
            "shack_hartmann_diffractive_tilt_mode": modal_tiptilt_case(
                root,
                case_id="shack_hartmann_diffractive_tilt_mode",
                kind="shack_hartmann_slopes",
                mode_index=2,
                amplitude=5e-9,
            ),
            "shack_hartmann_diffractive_tiptilt_dm": composite_tiptilt_dm_case(
                root,
                case_id="shack_hartmann_diffractive_tiptilt_dm",
                kind="shack_hartmann_slopes",
                tip_amplitude=5e-9,
                dm_command=composite_dm_command,
            ),
            "pyramid_diffractive_ramp": pyramid_case(root),
            "pyramid_diffractive_tip_mode": modal_tiptilt_case(
                root,
                case_id="pyramid_diffractive_tip_mode",
                kind="pyramid_slopes",
                mode_index=1,
                amplitude=5e-9,
            ),
            "pyramid_diffractive_tiptilt_dm": composite_tiptilt_dm_case(
                root,
                case_id="pyramid_diffractive_tiptilt_dm",
                kind="pyramid_slopes",
                tip_amplitude=5e-9,
                dm_command=composite_dm_command,
            ),
            "bioedge_diffractive_ramp": bioedge_case(root),
            "bioedge_diffractive_tip_mode": modal_tiptilt_case(
                root,
                case_id="bioedge_diffractive_tip_mode",
                kind="bioedge_slopes",
                mode_index=1,
                amplitude=5e-9,
            ),
            "bioedge_diffractive_tiptilt_dm": composite_tiptilt_dm_case(
                root,
                case_id="bioedge_diffractive_tiptilt_dm",
                kind="bioedge_slopes",
                tip_amplitude=5e-9,
                dm_command=composite_dm_command,
            ),
            "gain_sensing_camera_optical_gains": gsc_case(root),
            "transfer_function_rejection": transfer_function_case(root),
            "lift_interaction_matrix": lift_case(root),
            "closed_loop_shack_hartmann_trace": closed_loop_trace_case(root, case_id="closed_loop_shack_hartmann_trace", kind="shack_hartmann_slopes"),
            "closed_loop_pyramid_trace": closed_loop_trace_case(root, case_id="closed_loop_pyramid_trace", kind="pyramid_slopes"),
            "closed_loop_bioedge_trace": closed_loop_trace_case(root, case_id="closed_loop_bioedge_trace", kind="bioedge_slopes"),
            "gsc_closed_loop_trace": gsc_closed_loop_trace_case(root),
            "gsc_atmosphere_replay_trace_bounded": gsc_atmosphere_replay_trace_case(
                root,
                case_id="gsc_atmosphere_replay_trace_bounded",
                n_iter=2,
            ),
        }
        cases.update(gsc_branch_step_cases(root))
        metadata = {
            "repo_url": args.oopao_repo if args.oopao_path is None else str(checkout),
            "requested_ref": args.oopao_ref if args.oopao_path is None else "local-checkout",
            "resolved_commit": commit,
        }
        write_manifest(root / "manifest.toml", cases, metadata=metadata)
    finally:
        if workspace is not None:
            workspace.cleanup()


if __name__ == "__main__":
    main()

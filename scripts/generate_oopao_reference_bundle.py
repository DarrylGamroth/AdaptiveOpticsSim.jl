#!/usr/bin/env python3
"""Generate a deterministic OOPAO reference bundle for AdaptiveOptics.jl tests.

Run this inside an environment where OOPAO and its Python dependencies are
available. The output format matches `test/reference_harness.jl`.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from OOPAO.BioEdge import BioEdge
from OOPAO.Detector import Detector
from OOPAO.GainSensingCamera import GainSensingCamera
from OOPAO.LiFT import LiFT
from OOPAO.Pyramid import Pyramid
from OOPAO.ShackHartmann import ShackHartmann
from OOPAO.Source import Source
from OOPAO.Telescope import Telescope


def make_telescope(*, resolution: int, diameter: float, sampling_time: float, central_obstruction: float) -> Telescope:
    return Telescope(
        resolution=resolution,
        diameter=diameter,
        samplingTime=sampling_time,
        centralObstruction=central_obstruction,
        display_optical_path=False,
        fov=0,
    )


def make_source(*, band: str, magnitude: float) -> Source:
    return Source(optBand=band, magnitude=magnitude, coordinates=[0, 0])


def apply_ramp_opd(tel: Telescope, *, scale_x: float, scale_y: float, bias: float = 0.0) -> None:
    x = np.arange(tel.resolution, dtype=float)[:, None]
    y = np.arange(tel.resolution, dtype=float)[None, :]
    tel.OPD_no_pupil = bias + scale_x * x + scale_y * y


def write_array(path: Path, data: np.ndarray) -> None:
    np.savetxt(path, np.asarray(data, dtype=np.float64).reshape(-1), fmt="%.18e")


def write_manifest(path: Path, cases: dict[str, dict]) -> None:
    lines = ["version = 1", ""]
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
        for section in ("telescope", "source", "opd", "wfs", "basis", "detector", "compute", "compare"):
            if section not in case:
                continue
            lines.append(f"[cases.{case_id}.{section}]")
            for key, value in case[section].items():
                lines.append(f"{key} = {toml_literal(value)}")
            lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path, help="Bundle output directory")
    args = parser.parse_args()

    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
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
        "pyramid_diffractive_ramp": pyramid_case(root),
        "bioedge_diffractive_ramp": bioedge_case(root),
        "gain_sensing_camera_optical_gains": gsc_case(root),
        "transfer_function_rejection": transfer_function_case(root),
        "lift_interaction_matrix": lift_case(root),
    }
    write_manifest(root / "manifest.toml", cases)


if __name__ == "__main__":
    main()

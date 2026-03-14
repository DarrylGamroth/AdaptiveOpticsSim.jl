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
                f"atol = {case['atol']:.16e}",
                f"rtol = {case['rtol']:.16e}",
                "",
            ]
        )
        for section in ("telescope", "source", "opd", "wfs", "compute", "compare"):
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
            "modulation_points": 1,
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
            "diffraction_padding": 2,
            "psf_centering": True,
            "n_pix_separation": 4,
            "binning": 2,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path, help="Bundle output directory")
    args = parser.parse_args()

    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    cases = {
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
    }
    write_manifest(root / "manifest.toml", cases)


if __name__ == "__main__":
    main()

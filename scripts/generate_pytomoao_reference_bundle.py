#!/usr/bin/env python3
"""Generate deterministic pyTomoAO reference data for AdaptiveOpticsSim.jl."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import yaml


def write_array(path: Path, data: np.ndarray) -> None:
    np.savetxt(path, np.asarray(data, dtype=np.float64).reshape(-1), fmt="%.18e")


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


def merge_cases(path: Path, new_cases: dict[str, dict]) -> None:
    if path.exists():
        import tomllib

        raw = tomllib.loads(path.read_text(encoding="utf-8"))
        cases = dict(raw.get("cases", {}))
    else:
        cases = {}
    cases.update(new_cases)

    def emit_table(lines: list[str], prefix: list[str], mapping: dict) -> None:
        scalar_items = [(k, v) for k, v in mapping.items() if not isinstance(v, dict)]
        nested_items = [(k, v) for k, v in mapping.items() if isinstance(v, dict)]
        if prefix:
            lines.append(f"[{'.'.join(prefix)}]")
        for key, value in scalar_items:
            lines.append(f"{key} = {toml_literal(value)}")
        if prefix or scalar_items:
            lines.append("")
        for key, value in nested_items:
            emit_table(lines, prefix + [key], value)

    lines = ["version = 1", ""]
    for case_id in sorted(cases):
        emit_table(lines, ["cases", case_id], cases[case_id])
    path.write_text("\n".join(lines), encoding="utf-8")


def compact_config() -> dict:
    return {
        "lgs_wfs_parameters": {
            "D": 8.0,
            "nLenslet": 2,
            "nPx": 4,
            "fieldStopSize": 2.0,
            "nLGS": 1,
            "validLLMap": [[1, 1], [1, 1]],
            "wfsLensletsRotation": [0.0],
            "wfsLensletsOffset": [[0.0], [0.0]],
        },
        "dm_parameters": {
            "dmHeights": [0.0],
            "dmPitch": [0.5],
            "dmCrossCoupling": 0.2,
            "nActuators": [2],
            "validActuators": [[1, 1], [1, 1]],
        },
        "atmosphere_parameters": {
            "nLayer": 1,
            "zenithAngleInDeg": 0.0,
            "altitude": [0.0],
            "L0": 25.0,
            "r0": 0.2,
            "fractionnalR0": [1.0],
            "wavelength": 5.0e-7,
            "windDirection": [0.0],
            "windSpeed": [10.0],
        },
        "lgs_asterism": {
            "radiusAst": 7.6,
            "LGSwavelength": 5.89e-7,
            "baseLGSHeight": 90000.0,
            "nLGS": 1,
        },
        "tomography_parameters": {
            "fovOptimization": 0.0,
            "nFitSrc": 1,
        },
    }


def write_temp_config(path: Path, config: dict) -> None:
    with path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)


def common_sections() -> dict:
    cfg = compact_config()
    return {
        "atmosphere": {
            "zenith_angle_deg": 0.0,
            "altitude_km": [0.0],
            "L0": 25.0,
            "r0_zenith": 0.2,
            "fractional_r0": [1.0],
            "wavelength": 5.0e-7,
            "wind_direction_deg": [0.0],
            "wind_speed": [10.0],
        },
        "asterism": {
            "radius_arcsec": 7.6,
            "wavelength": 5.89e-7,
            "base_height_m": 90000.0,
            "n_lgs": 1,
        },
        "wfs": {
            "diameter": 8.0,
            "n_lenslet": 2,
            "n_px": 4,
            "field_stop_size_arcsec": 2.0,
            "valid_lenslet_map": [[1, 1], [1, 1]],
            "lenslet_rotation_rad": [0.0],
            "lenslet_offset": [[0.0], [0.0]],
        },
        "tomography": {
            "n_fit_src": 1,
            "fov_optimization_arcsec": 0.0,
            "fit_src_height_m": "Inf",
        },
        "dm": {
            "heights_m": [0.0],
            "pitch_m": [0.5],
            "cross_coupling": 0.2,
            "n_actuators": [2],
            "valid_actuators": [[1, 1], [1, 1]],
        },
    }


def sections_from_pytomoao_config(cfg: dict) -> dict:
    atm = cfg["atmosphere_parameters"]
    asterism = cfg["lgs_asterism"]
    wfs = cfg["lgs_wfs_parameters"]
    tomo = cfg["tomography_parameters"]
    dm = cfg["dm_parameters"]
    return {
        "atmosphere": {
            "zenith_angle_deg": float(atm["zenithAngleInDeg"]),
            "altitude_km": [float(x) for x in atm["altitude"]],
            "L0": float(atm["L0"]),
            "r0_zenith": float(atm["r0"]),
            "fractional_r0": [float(x) for x in atm["fractionnalR0"]],
            "wavelength": float(atm["wavelength"]),
            "wind_direction_deg": [float(x) for x in atm["windDirection"]],
            "wind_speed": [float(x) for x in atm["windSpeed"]],
        },
        "asterism": {
            "radius_arcsec": float(asterism["radiusAst"]),
            "wavelength": float(asterism["LGSwavelength"]),
            "base_height_m": float(asterism["baseLGSHeight"]),
            "n_lgs": int(asterism["nLGS"]),
        },
        "wfs": {
            "diameter": float(wfs["D"]),
            "n_lenslet": int(wfs["nLenslet"]),
            "n_px": int(wfs["nPx"]),
            "field_stop_size_arcsec": float(wfs["fieldStopSize"]),
            "valid_lenslet_map": wfs["validLLMap"],
            "lenslet_rotation_rad": [float(x) for x in wfs.get("wfsLensletsRotation", [0.0] * int(wfs["nLGS"]))],
            "lenslet_offset": wfs.get("wfsLensletsOffset", [[0.0] * int(wfs["nLGS"]), [0.0] * int(wfs["nLGS"])]),
        },
        "tomography": {
            "n_fit_src": int(tomo["nFitSrc"]),
            "fov_optimization_arcsec": float(tomo["fovOptimization"]),
            "fit_src_height_m": "Inf",
        },
        "dm": {
            "heights_m": [float(x) for x in dm["dmHeights"]],
            "pitch_m": [float(x) for x in dm["dmPitch"]],
            "cross_coupling": float(dm["dmCrossCoupling"]),
            "n_actuators": [int(x) for x in dm["nActuators"]],
            "valid_actuators": dm["validActuators"],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pytomoao-root",
        type=Path,
        default=Path("../pyTomoAO"),
        help="path to the pyTomoAO checkout",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("test/reference_data"),
        help="bundle output directory",
    )
    args = parser.parse_args()

    sys.path.insert(0, str(args.pytomoao_root.resolve()))
    from pyTomoAO.tomographicReconstructor import tomographicReconstructor

    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    config_path = output_root / "pytomoao_compact_config.yaml"
    write_temp_config(config_path, compact_config())

    reconstructor = tomographicReconstructor(str(config_path))
    reconstructor.build_reconstructor()

    slopes_model = np.array([0.1, -0.2, 0.05, 0.15, -0.1, 0.25, -0.05, 0.2], dtype=np.float64)
    model_wavefront = np.asarray(reconstructor.reconstruct_wavefront(slopes_model), dtype=np.float64)

    imat = np.array(
        [
            [1.0, 0.0, 0.2, 0.0],
            [0.0, 1.0, 0.0, 0.2],
            [0.5, 0.1, 0.0, 0.0],
            [0.0, 0.5, 0.1, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [0.25, 0.0, 0.25, 0.0],
            [0.0, 0.25, 0.0, 0.25],
        ],
        dtype=np.float64,
    )
    reconstructor_im = tomographicReconstructor(str(config_path))
    reconstructor_im.build_reconstructor(IM=imat)
    slopes_im = np.array([0.2, -0.1, 0.3, -0.2, 0.15, -0.05, 0.1, -0.15], dtype=np.float64)
    im_wavefront = np.asarray(reconstructor_im.reconstruct_wavefront(slopes_im), dtype=np.float64)

    benchmark_cfg_path = args.pytomoao_root.resolve() / "tests" / "tomography_config_kapa.yaml"
    benchmark_cfg = yaml.safe_load(benchmark_cfg_path.read_text(encoding="utf-8"))
    benchmark_sections = sections_from_pytomoao_config(benchmark_cfg)
    benchmark_reconstructor = tomographicReconstructor(str(benchmark_cfg_path))
    benchmark_reconstructor.build_reconstructor()
    n_valid = int(benchmark_reconstructor.lgsWfsParams.nValidSubap)
    benchmark_slopes = np.zeros(n_valid * 2, dtype=np.float64)
    benchmark_slopes[: n_valid - 1] = 4.0
    benchmark_slopes[n_valid:] = -4.0
    benchmark_wavefront_slopes = np.tile(benchmark_slopes, benchmark_reconstructor.nLGS)
    benchmark_wavefront = np.asarray(
        benchmark_reconstructor.reconstruct_wavefront(benchmark_wavefront_slopes),
        dtype=np.float64,
    )
    benchmark_reconstructor.assemble_reconstructor_and_fitting(
        nChannels=1,
        slopesOrder="simu",
        scalingFactor=1.5e7,
    )
    benchmark_dm_commands = np.asarray(
        benchmark_reconstructor.FR @ benchmark_slopes,
        dtype=np.float64,
    )
    benchmark_reconstructor.mask_DM_actuators(174)
    benchmark_dm_commands_masked = np.asarray(
        benchmark_reconstructor.FR @ benchmark_slopes,
        dtype=np.float64,
    )

    sections = common_sections()
    cases = {
        "tomography_model_gamma": {
            "kind": "tomography_model_gamma",
            "data": "tomography_model_gamma.txt",
            "shape": list(reconstructor.Gamma.toarray().shape),
            "storage_order": "C",
            "atol": 1e-12,
            "rtol": 1e-12,
            **sections,
        },
        "tomography_model_cxx": {
            "kind": "tomography_model_cxx",
            "data": "tomography_model_cxx.txt",
            "shape": list(reconstructor.Cxx.shape),
            "storage_order": "C",
            "atol": 5e-10,
            "rtol": 5e-8,
            **sections,
        },
        "tomography_model_cox": {
            "kind": "tomography_model_cox",
            "data": "tomography_model_cox.txt",
            "shape": list(reconstructor.Cox.shape),
            "storage_order": "C",
            "atol": 5e-10,
            "rtol": 5e-8,
            **sections,
        },
        "tomography_model_cnz": {
            "kind": "tomography_model_cnz",
            "data": "tomography_model_cnz.txt",
            "shape": list(reconstructor.CnZ.shape),
            "storage_order": "C",
            "atol": 1e-12,
            "rtol": 1e-12,
            **sections,
        },
        "tomography_model_reconstructor": {
            "kind": "tomography_model_reconstructor",
            "data": "tomography_model_reconstructor.txt",
            "shape": list(reconstructor.reconstructor.shape),
            "storage_order": "C",
            "atol": 1e-9,
            "rtol": 1e-7,
            **sections,
        },
        "tomography_model_wavefront": {
            "kind": "tomography_model_wavefront",
            "data": "tomography_model_wavefront.txt",
            "shape": list(model_wavefront.shape),
            "storage_order": "C",
            "atol": 1e-9,
            "rtol": 1e-7,
            **sections,
            "compute": {
                "slopes": slopes_model.tolist(),
            },
        },
        "tomography_im_reconstructor": {
            "kind": "tomography_im_reconstructor",
            "data": "tomography_im_reconstructor.txt",
            "shape": list(reconstructor_im.reconstructor.shape),
            "storage_order": "C",
            "atol": 1e-9,
            "rtol": 1e-7,
            **sections,
            "compute": {
                "interaction_matrix": imat.tolist(),
            },
        },
        "tomography_im_wavefront": {
            "kind": "tomography_im_wavefront",
            "data": "tomography_im_wavefront.txt",
            "shape": list(im_wavefront.shape),
            "storage_order": "C",
            "atol": 1e-9,
            "rtol": 1e-7,
            **sections,
            "compute": {
                "interaction_matrix": imat.tolist(),
                "slopes": slopes_im.tolist(),
            },
        },
        "tomography_kapa_model_wavefront": {
            "kind": "tomography_model_wavefront",
            "data": "tomography_kapa_model_wavefront.txt",
            "shape": list(benchmark_wavefront.shape),
            "storage_order": "C",
            "atol": 1e-9,
            "rtol": 1e-7,
            **benchmark_sections,
            "compute": {
                "slopes_generator": "pytomoao_tiptilt",
                "positive_amplitude": 4.0,
                "negative_amplitude": -4.0,
                "n_channels": int(benchmark_reconstructor.nLGS),
            },
        },
        "tomography_kapa_model_dm_commands": {
            "kind": "tomography_model_dm_commands",
            "data": "tomography_kapa_model_dm_commands.txt",
            "shape": list(benchmark_dm_commands.shape),
            "atol": 5e-7,
            "rtol": 5e-6,
            **benchmark_sections,
            "compute": {
                "slopes_generator": "pytomoao_tiptilt",
                "positive_amplitude": 4.0,
                "negative_amplitude": -4.0,
                "n_channels": 1,
                "assembly": {
                    "n_channels": 1,
                    "slope_order": "simu",
                    "scaling_factor": 1.5e7,
                },
            },
        },
        "tomography_kapa_model_dm_commands_masked": {
            "kind": "tomography_model_dm_commands",
            "data": "tomography_kapa_model_dm_commands_masked.txt",
            "shape": list(benchmark_dm_commands_masked.shape),
            "atol": 5e-7,
            "rtol": 5e-6,
            **benchmark_sections,
            "compute": {
                "slopes_generator": "pytomoao_tiptilt",
                "positive_amplitude": 4.0,
                "negative_amplitude": -4.0,
                "n_channels": 1,
                "assembly": {
                    "n_channels": 1,
                    "slope_order": "simu",
                    "scaling_factor": 1.5e7,
                    "mask_actuators": [175],
                },
            },
        },
    }

    write_array(output_root / "tomography_model_gamma.txt", reconstructor.Gamma.toarray())
    write_array(output_root / "tomography_model_cxx.txt", reconstructor.Cxx)
    write_array(output_root / "tomography_model_cox.txt", reconstructor.Cox)
    write_array(output_root / "tomography_model_cnz.txt", reconstructor.CnZ)
    write_array(output_root / "tomography_model_reconstructor.txt", reconstructor.reconstructor)
    write_array(output_root / "tomography_model_wavefront.txt", model_wavefront)
    write_array(output_root / "tomography_im_reconstructor.txt", reconstructor_im.reconstructor)
    write_array(output_root / "tomography_im_wavefront.txt", im_wavefront)
    write_array(output_root / "tomography_kapa_model_wavefront.txt", benchmark_wavefront)
    write_array(output_root / "tomography_kapa_model_dm_commands.txt", benchmark_dm_commands)
    write_array(output_root / "tomography_kapa_model_dm_commands_masked.txt", benchmark_dm_commands_masked)
    merge_cases(output_root / "manifest.toml", cases)


if __name__ == "__main__":
    main()

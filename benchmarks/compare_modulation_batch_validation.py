"""Check archived fixed-sampling integration parity; no GPU work.

Usage: python compare_modulation_batch_validation.py EVIDENCE OUTPUT.json
Requires NumPy. Bounds concern implementation parity, not instrument fidelity.
"""
import hashlib
import json
from pathlib import Path
import sys
import tomllib

import numpy as np


def load_array(directory, record):
    path = directory / record["file"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == record["sha256"]
    dtype = {"Float32": "<f4", "Float64": "<f8"}[record["eltype"]]
    return np.fromfile(path, dtype=dtype).reshape(record["shape"], order="F").astype(np.float64)


def difference(actual, reference):
    delta = actual - reference
    return {"relative_l2": float(np.linalg.norm(delta) / np.linalg.norm(reference)),
            "maximum_absolute": float(np.max(np.abs(delta))),
            "identical": bool(np.array_equal(actual, reference))}


def compare(evidence):
    result = {"scope": "Same discretization, seeds, calibration and controller; finite recorded workloads"}
    before, after = evidence / "copper-before", evidence / "copper-after"
    old = tomllib.loads((before / "validation.toml").read_text())
    new = tomllib.loads((after / "validation.toml").read_text())
    assert old["complete"] and new["complete"]
    old, new = old["cases"]["n240-m96"], new["cases"]["n240-m96"]
    copper = {}
    for name in ("rates", "slopes", "interaction", "atmosphere_rates", "atmosphere_frames"):
        a, b = load_array(before, old[name]), load_array(after, new[name])
        copper[name] = difference(b, a)
        if name == "rates":
            copper["maximum_probe_flux_relative_error"] = float(np.max(np.abs(
                b.sum(axis=(0, 1)) / a.sum(axis=(0, 1)) - 1)))
        if name == "interaction":
            relative = np.linalg.norm(b - a, axis=0) / np.linalg.norm(a, axis=0)
            copper["maximum_column_relative_error"] = float(relative.max())
            copper["maximum_error_column_1based"] = int(relative.argmax() + 1)
    assert copper["rates"]["relative_l2"] < 2e-6
    assert copper["maximum_probe_flux_relative_error"] < 2e-6
    assert copper["interaction"]["relative_l2"] < 2e-5
    assert copper["maximum_column_relative_error"] < 5e-4
    old_loop = tomllib.loads((before / "closed-loop.toml").read_text())
    new_loop = tomllib.loads((after / "closed-loop.toml").read_text())
    assert old_loop["all_finite"] and new_loop["all_finite"]
    assert old_loop["frozen_calibration"] == new_loop["frozen_calibration"]
    loop = {}
    for name in ("commands", "slopes", "diagnostics"):
        a, b = load_array(before, old_loop[name]), load_array(after, new_loop[name])
        loop[name] = difference(b, a)
        if name == "diagnostics":
            loop["maximum_residual_wfe_difference_m"] = float(np.max(np.abs(b[0] - a[0])))
            loop["maximum_strehl_difference"] = float(np.max(np.abs(b[1] - a[1])))
    assert loop["commands"]["relative_l2"] < 1e-4
    assert loop["maximum_residual_wfe_difference_m"] < 5e-11
    assert loop["maximum_strehl_difference"] < 1e-4
    loop["tail_residual_before_m"] = old_loop["tail_residual_rms_m"]
    loop["tail_residual_after_m"] = new_loop["tail_residual_rms_m"]
    copper["closed_loop"] = loop
    for name in ("calibration_captured", "atmosphere_captured"):
        copper[name] = {"before_median_ms": old[name]["median_ms"],
                        "after_median_ms": new[name]["median_ms"]}
    result["copper_240_96"] = copper
    before, after = evidence / "gpi-before", evidence / "gpi-after"
    for directory in (before, after):
        assert tomllib.loads((directory / "validation.toml").read_text())["complete"]
    names = sorted(p.relative_to(before) for p in before.rglob("*")
                   if p.suffix in (".f32", ".f64"))
    assert names
    assert names == sorted(p.relative_to(after) for p in after.rglob("*")
                           if p.suffix in (".f32", ".f64"))
    for name in names:
        assert (before / name).read_bytes() == (after / name).read_bytes(), str(name)
    result["gpi_240_24"] = {"identical_saved_arrays": len(names),
                             "files": list(map(str, names)),
                             "scope": "All saved coefficient histories, calibration arrays and final frames; unsaved intermediate images excluded"}
    result["passed"] = True
    return result


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("EVIDENCE OUTPUT.json")
    report = compare(Path(sys.argv[1]))
    Path(sys.argv[2]).write_text(json.dumps(report, indent=2) + "\n")
    print("Fixed-sampling integration parity passed")

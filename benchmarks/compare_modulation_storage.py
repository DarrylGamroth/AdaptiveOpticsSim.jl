"""Compare fixed-model CUDA captures produced by validate_modulation_storage_cuda.jl.

Usage: python3 compare_modulation_storage.py BASELINE_DIR CHANGED_DIR OUTPUT_JSON
Requires NumPy. Each input array is checked against its archived SHA-256.
"""
import hashlib
import json
import pathlib
import sys
import tomllib

import numpy as np


def read_array(directory, record):
    path = directory / record["file"]
    raw = path.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == record["sha256"], path
    dtype = {"Float32": "<f4", "Float64": "<f8"}[record["eltype"]]
    result = np.frombuffer(raw, dtype=dtype).reshape(record["shape"], order="F")
    assert np.isfinite(result).all(), path
    return result


def compare(baseline, changed):
    directories = (pathlib.Path(baseline), pathlib.Path(changed))
    records = [tomllib.loads((d / "validation.toml").read_text()) for d in directories]
    assert all(r["complete"] for r in records)
    assert records[0]["harness_sha256"] == records[1]["harness_sha256"]
    assert records[0]["helper_sha256"] == records[1]["helper_sha256"]
    assert records[0]["cases"].keys() == records[1]["cases"].keys()
    result = {"harness_sha256": records[0]["harness_sha256"], "cases": {}}
    for case in records[0]["cases"]:
        before, after = (r["cases"][case] for r in records)
        assert before["probes"] == after["probes"]
        assert before["divisor"] == after["divisor"]
        comparisons = {}
        for name in ("rates", "slopes", "interaction", "atmosphere_rates",
                     "atmosphere_frames", "stream_rates", "stream_frames"):
            a = read_array(directories[0], before[name])
            b = read_array(directories[1], after[name])
            assert a.shape == b.shape
            delta = b.astype(np.float64) - a.astype(np.float64)
            reference_norm = np.linalg.norm(a.astype(np.float64))
            comparisons[name] = {
                "shape": list(a.shape), "identical": bool(np.array_equal(a, b)),
                "max_abs": float(np.max(np.abs(delta))),
                "relative_l2": float(np.linalg.norm(delta) / reference_norm)
                if reference_norm else float(np.linalg.norm(delta)),
            }
            # The shifted-mask rate kernels and discrete physical model did not
            # change. Require exact replay, including every calibration column.
            assert comparisons[name]["identical"], (case, name, comparisons[name])
        timings = {}
        for name in ("calibration_captured", "atmosphere_captured", "atmosphere_stream"):
            timings[name] = {"baseline_median_ms": before[name]["median_ms"],
                             "changed_median_ms": after[name]["median_ms"],
                             "baseline_p95_ms": before[name]["p95_ms"],
                             "changed_p95_ms": after[name]["p95_ms"],
                             "baseline_maximum_ms": before[name]["maximum_ms"],
                             "changed_maximum_ms": after[name]["maximum_ms"]}
            if "captured" in name:
                assert all(x == 0 for x in before[name]["host_bytes"])
                assert all(x == 0 for x in after[name]["host_bytes"])
        result["cases"][case] = {"arrays": comparisons, "timings": timings,
                                 "baseline_modulation": before["modulation"],
                                 "changed_modulation": after["modulation"]}
    return result


if __name__ == "__main__":
    assert len(sys.argv) == 4, __doc__
    result = compare(sys.argv[1], sys.argv[2])
    pathlib.Path(sys.argv[3]).write_text(json.dumps(result, indent=2) + "\n")
    print("Exact parity:", ", ".join(result["cases"]))

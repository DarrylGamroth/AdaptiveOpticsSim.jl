"""Summarize repeated baseline/selected AOS batch timings (stdlib only).

Usage: python summarize_modulation_batches.py OUTPUT.json TIMING.toml ...
Uses the repeated randomized baseline trials, excluding initial measurements.
"""
import json
from pathlib import Path
from statistics import median
import sys
import tomllib


def summarize(paths):
    rows = {}
    sample_count = allocation_observations = 0
    for path in paths:
        data = tomllib.loads(path.read_text())
        assert data["complete"]
        for case in data["cases"]:
            key = (case["name"], case["strategy"])
            row = rows.setdefault(key, {
                "name": case["name"], "strategy": case["strategy"],
                "fft_side": case["baseline"]["padded_side"],
                "points": case["configuration"]["points"],
                "baseline_batch": case["baseline"]["batch"],
                "selected_batch": case["selected"]["batch"],
                "baseline_stack_bytes": case["baseline"]["field_stack_bytes"],
                "selected_stack_bytes": case["selected"]["field_stack_bytes"],
                "processes": [], "maximum_output_relative_l2": 0.,
            })
            trials = case["trials"]
            baseline = [t for t in trials if t["policy"] == "baseline"]
            selected = [t for t in trials if t["policy"] == "selected"]
            assert len(baseline) == len(selected) == data["rounds"]
            old, new = median(t["median_ms"] for t in baseline), median(t["median_ms"] for t in selected)
            row["processes"].append({"file": path.name,
                "baseline_round_medians_ms": [t["median_ms"] for t in baseline],
                "selected_round_medians_ms": [t["median_ms"] for t in selected],
                "baseline_median_ms": old, "selected_median_ms": new,
                "reduction_percent": 100 * (1 - new / old),
                "baseline_p95_ms": [t["p95_ms"] for t in baseline],
                "selected_p95_ms": [t["p95_ms"] for t in selected],
                "baseline_maximum_ms": max(t["maximum_ms"] for t in baseline),
                "selected_maximum_ms": max(t["maximum_ms"] for t in selected)})
            for trial in [case["baseline"], case["selected"], *trials]:
                assert all(n == 0 for n in trial["host_bytes"])
                allocation_observations += len(trial["host_bytes"])
                sample_count += len(trial["samples_ms"])
                row["maximum_output_relative_l2"] = max(row["maximum_output_relative_l2"], trial["relative_l2"])
    for row in rows.values():
        row["baseline_median_ms"] = median(p["baseline_median_ms"] for p in row["processes"])
        row["selected_median_ms"] = median(p["selected_median_ms"] for p in row["processes"])
        row["reduction_percent"] = 100 * (1 - row["selected_median_ms"] / row["baseline_median_ms"])
    return {"scope": "Unpaced completed captured AOS photon-rate node; two policies randomized within each round; initial timings excluded from speedup",
            "process_count": len(paths), "case_strategy_count": len(rows),
            "timed_completions_including_initial": sample_count,
            "zero_allocation_observations": allocation_observations,
            "cases": list(rows.values())}


if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise SystemExit("OUTPUT.json TIMING.toml ...")
    report = summarize([Path(p) for p in sys.argv[2:]])
    Path(sys.argv[1]).write_text(json.dumps(report, indent=2) + "\n")

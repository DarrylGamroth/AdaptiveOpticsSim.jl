"""Verify saved AOS/OOPAO comparisons, separating batching and model errors.

Usage: python compare.py REFERENCE_DIRECTORY OUTPUT.json
NumPy required; no GPU or simulator is launched.
"""
import json
from pathlib import Path
import sys
import tomllib

import numpy as np


def relative(actual, reference):
    return float(np.linalg.norm(actual - reference) / np.linalg.norm(reference))


def compare(root):
    reference_summary = json.loads((root / "summary.json").read_text())
    result = {"oopao_commit": reference_summary["oopao_commit"], "cases": {}}
    for name, configuration in reference_summary.items():
        if not isinstance(configuration, dict) or "resolution" not in configuration:
            continue
        directory = root / name
        comparisons = {}
        for flat_path in sorted(directory.glob("*_flat_frame.txt")):
            prefix = flat_path.name.removesuffix("_flat_frame.txt")
            if not prefix.startswith(("aos", "cuda_f32")):
                continue
            metrics = {}
            for label in ("flat_frame", "ramp_frame", "tip_frame_derivative"):
                reference = np.loadtxt(directory / (label + ".txt"))
                actual = np.loadtxt(directory / (prefix + "_" + label + ".txt"))
                metrics[label] = relative(actual, reference)
            parameters = directory / (prefix + "_parameters.toml")
            metadata = tomllib.loads(parameters.read_text()) if parameters.exists() else {}
            shifted = "shifted_mask" in prefix
            single_precision = "f32" in prefix
            if not shifted:
                image_bound, response_bound = (2e-5, 2e-4) if single_precision else (1e-6, 1e-6)
                assert metrics["flat_frame"] < image_bound, (name, prefix, metrics)
                assert metrics["ramp_frame"] < image_bound, (name, prefix, metrics)
                assert metrics["tip_frame_derivative"] < response_bound, (name, prefix, metrics)
            comparisons[prefix] = {"relative_l2_vs_oopao": metrics,
                "parameters": metadata, "physical_equivalence_asserted": not shifted}
        result["cases"][name] = comparisons
    assert len(result["cases"]) == 3
    tuned = root / "n120_p30_m3"
    batching = {}
    for strategy in ("pupil_tilt", "shifted_mask"):
        metrics = {}
        for policy, expected in (("baseline", 32), ("selected", 8)):
            prefix = "cuda_f32_" + strategy + "_" + policy
            metadata = tomllib.loads((tuned / (prefix + "_parameters.toml")).read_text())
            assert metadata["actual_batch"] == expected
        for label in ("flat_frame", "ramp_frame", "tip_frame_derivative"):
            before = np.loadtxt(tuned / ("cuda_f32_" + strategy + "_baseline_" + label + ".txt"))
            after = np.loadtxt(tuned / ("cuda_f32_" + strategy + "_selected_" + label + ".txt"))
            metrics[label] = relative(after, before)
            assert metrics[label] < (2e-4 if label == "tip_frame_derivative" else 2e-6)
        batching[strategy] = metrics
    result["tuned_float32_batching_relative_l2"] = batching
    result["pupil_tilt_comparison_passed"] = True
    result["batching_comparison_passed"] = True
    result["shifted_mask_physical_equivalence_qualified"] = False
    result["scope"] = "Noiseless photon images and central +/-1nm tip response; no external closed-loop or detector-statistics qualification"
    return result


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("REFERENCE_DIRECTORY OUTPUT.json")
    report = compare(Path(sys.argv[1]))
    Path(sys.argv[2]).write_text(json.dumps(report, indent=2) + "\n")
    print("Batching and physical pupil-tilt comparisons passed; shifted-mask equivalence remains unqualified")

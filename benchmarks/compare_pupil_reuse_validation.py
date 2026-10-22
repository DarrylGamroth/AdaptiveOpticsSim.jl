"""Compare retained physical-model implementation changes at fixed sampling.

Usage: python compare_pupil_reuse_validation.py OLD_EVIDENCE NEW_EVIDENCE OUTPUT.json
OLD_EVIDENCE is the completed physical-model campaign (not retired shifted-mask data).
Requires NumPy; reads and verifies recorded array hashes.
"""
import json
from pathlib import Path
import sys
import tomllib
import numpy as np
from compare_modulation_batch_validation import load_array, difference


def metadata(directory):
    return tomllib.loads((directory / "validation.toml").read_text())


def compare(old_root, new_root):
    result = {"scope": "Fixed physical sampling; implementation parity, not instrument fidelity"}
    old_dir, new_dir = old_root / "copper-after", new_root / "copper-after"
    old, new = metadata(old_dir), metadata(new_dir)
    assert old["complete"] and new["complete"]
    old, new = old["cases"]["n240-m96"], new["cases"]["n240-m96"]
    copper = {}
    for name in ("rates", "slopes", "interaction", "atmosphere_rates", "atmosphere_frames"):
        a, b = load_array(old_dir, old[name]), load_array(new_dir, new[name])
        copper[name] = difference(b, a)
        if name == "rates":
            copper["maximum_probe_flux_relative_error"] = float(np.max(np.abs(b.sum((0, 1)) / a.sum((0, 1)) - 1)))
        if name == "interaction":
            columns = np.linalg.norm(b-a, axis=0) / np.linalg.norm(a, axis=0)
            copper["maximum_column_relative_error"] = float(columns.max())
    # Same gates as the earlier FFT implementation campaign.
    assert copper["rates"]["relative_l2"] < 2e-6
    assert copper["maximum_probe_flux_relative_error"] < 2e-6
    assert copper["interaction"]["relative_l2"] < 2e-5
    assert copper["maximum_column_relative_error"] < 5e-4
    a = tomllib.loads((old_dir/"closed-loop.toml").read_text())
    b = tomllib.loads((new_dir/"closed-loop.toml").read_text())
    assert a["all_finite"] and b["all_finite"]
    assert a["frozen_calibration"] == b["frozen_calibration"]
    loop = {name: difference(load_array(new_dir,b[name]),load_array(old_dir,a[name]))
            for name in ("commands","slopes","diagnostics")}
    x, y = load_array(old_dir,a["diagnostics"]), load_array(new_dir,b["diagnostics"])
    loop["maximum_wfe_difference_m"] = float(np.max(np.abs(y[0]-x[0])))
    loop["maximum_strehl_difference"] = float(np.max(np.abs(y[1]-x[1])))
    assert loop["commands"]["relative_l2"] < 1e-4
    assert loop["maximum_wfe_difference_m"] < 5e-11
    assert loop["maximum_strehl_difference"] < 1e-4
    copper["frozen_calibration_loop"] = loop
    result["copper"] = copper

    old_dir, new_dir = old_root/"gpi-after", new_root/"gpi-after"
    a, b = metadata(old_dir), metadata(new_dir)
    assert a["complete"] and b["complete"]
    gram = load_array(old_dir,a["basis_gram"])
    assert np.array_equal(gram,load_array(new_dir,b["basis_gram"]))
    gpi = {}
    for case in ("default","fine","transfer"):
        ac, bc = metadata(old_dir/case), metadata(new_dir/case)
        assert ac["pupil_resolution"] == bc["pupil_resolution"]
        assert ac["modulation_points"] == bc["modulation_points"]
        rec = {}
        if case != "transfer":
            aa, bb = load_array(old_dir/case,ac["calibration"]["interaction_per_m"]), load_array(new_dir/case,bc["calibration"]["interaction_per_m"])
            rec["interaction"] = difference(bb,aa)
            rec["maximum_column_relative_error"] = float(np.max(np.linalg.norm(bb-aa,axis=0)/np.linalg.norm(aa,axis=0)))
            assert rec["interaction"]["relative_l2"] < 2e-5
            assert rec["maximum_column_relative_error"] < 5e-4
        runs = {}
        for run_name, old_run in ac["runs"].items():
            new_run = bc["runs"][run_name]
            assert old_run["all_finite"] and new_run["all_finite"]
            assert new_run["bounded_unsaturated_observed"] and new_run["saturation_samples"] == 0
            if "host_allocation_bytes" in new_run:
                assert all(x == 0 for x in new_run["host_allocation_bytes"])
            else:
                assert run_name == "static"  # Allocation sampling is scoped to moving runs.
            aa = load_array(old_dir/case/run_name,old_run["controller_coefficients_used"])
            bb = load_array(new_dir/case/run_name,new_run["controller_coefficients_used"])
            delta = bb-aa
            wavefront = np.sqrt(np.maximum(0,np.einsum('if,ij,jf->f',delta,gram,delta)))
            wfe_delta = np.abs(np.array(old_run["common_fine_residual_rms_m"])-new_run["common_fine_residual_rms_m"])
            rr = {"coefficients": difference(bb,aa),
                  "maximum_residual_wavefront_difference_m":float(wavefront.max()),
                  "tail_residual_wavefront_difference_rms_m":float(np.sqrt(np.mean(wavefront[-64:]**2))),
                  "maximum_wfe_difference_m":float(wfe_delta.max()),
                  "tail_wfe_before_m":old_run["tail_residual_rms_m"],
                  "tail_wfe_after_m":new_run["tail_residual_rms_m"]}
            rr["provisional_0p05_nm_vector_check_passed"] = bool(wavefront.max() < 5e-11)
            if run_name == "moving_noisy":
                # Integer stochastic detection is not a continuous coupling of rates.
                # Preserve and report the provisional vector check; it is not a
                # statistical-equivalence criterion for one paired noisy trajectory.
                rr["comparison_scope"] = "One seeded noisy realization; statistical equivalence not established"
                old_quiet = load_array(old_dir/case/"moving_noiseless",
                    ac["runs"]["moving_noiseless"]["controller_coefficients_used"])
                noise = aa-old_quiet
                noise_power = np.einsum('if,ij,jf->f',noise,gram,noise)
                rr["wavefront_difference_rms_m"] = float(np.sqrt(np.mean(wavefront**2)))
                rr["realized_noise_effect_rms_m"] = float(np.sqrt(np.mean(noise_power)))
            else:
                assert wavefront.max() < 5e-11
            assert wfe_delta.max() < 5e-11
            frame_a = load_array(old_dir/case/run_name,old_run["last_frame"])
            frame_b = load_array(new_dir/case/run_name,new_run["last_frame"])
            rr["last_frame"] = difference(frame_b,frame_a)
            runs[run_name] = rr
        rec["runs"] = runs
        gpi[case] = rec
    result["gpi"] = gpi
    result["passed"] = True
    return result


if __name__ == "__main__":
    if len(sys.argv)!=4:
        raise SystemExit(__doc__)
    report=compare(Path(sys.argv[1]),Path(sys.argv[2]))
    Path(sys.argv[3]).write_text(json.dumps(report,indent=2)+"\n")
    print("Physical pupil reuse integration parity passed")

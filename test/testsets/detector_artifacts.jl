@testset "Detector validation artifacts" begin
    detector_artifact_path = normpath(joinpath(@__DIR__, "..", "..", "benchmarks", "results",
        "detectors", "2026-07-12-detector-mkid-validation.toml"))
    @test isfile(detector_artifact_path)
    detector_artifact = TOML.parsefile(detector_artifact_path)
    @test all(values(detector_artifact["interpretation"]))
    @test issubset(Set(["apd", "spad_array", "mkid_array"]),
        Set(detector_artifact["scope"]["families"]))
    detector_manifest = TOML.parsefile(joinpath(dirname(detector_artifact_path),
        "manifest.toml"))
    detector_entries = Dict(entry["id"] => entry for entry in
        detector_manifest["artifacts"])
    @test detector_entries["DET-VAL-2026-07-12"]["status"] == "active"
    @test detector_entries["DET-VAL-2026-04-23"]["status"] == "superseded"
    @test detector_entries["DET-VAL-2026-04-23"]["superseded_by"] ==
        "DET-VAL-2026-07-12"

    low_fidelity_entry =
        detector_entries["DET-HIL-LOW-FIDELITY-2026-07-30"]
    @test low_fidelity_entry["status"] == "active"
    low_fidelity_path = joinpath(dirname(detector_artifact_path),
        low_fidelity_entry["path"])
    @test isfile(low_fidelity_path)
    low_fidelity = TOML.parsefile(low_fidelity_path)
    @test low_fidelity["schema_version"] == 2
    @test low_fidelity["benchmark"] == "detector_hil_latency"
    @test low_fidelity["all_gates_passed"]
    @test low_fidelity["contract"]["load_model"] ==
        "warmed serial self-paced service time"
    @test low_fidelity["contract"]["selected_card_ids"] == ["DET-HIL-00"]
    @test length(low_fidelity["cards"]) == 1
    card = only(low_fidelity["cards"])
    @test card["id"] == "DET-HIL-00"
    @test card["steady_alloc_bytes"] == 0
    @test card["regression"]["allocation_gate_passed"]
    @test all(card["runs"]) do run
        run["completed_operations"] == run["samples"] &&
            run["failed_operations"] == 0 &&
            run["histogram_recorded_bins"] > 0 &&
            !isempty(run["histogram_base64"])
    end

    ccd_entry = detector_entries["DET-CCD-QUAL-2026-07-30"]
    @test ccd_entry["status"] == "active"
    ccd_path = joinpath(dirname(detector_artifact_path),
        ccd_entry["path"])
    @test isfile(ccd_path)
    ccd = TOML.parsefile(ccd_path)
    @test ccd["schema_version"] == 1
    @test ccd["family"] == "conventional_ccd_single_read"
    @test ccd["all_gates_passed"]
    @test !ccd["environment"]["source_dirty"]
    qualification = ccd["qualification"]
    @test qualification["samples_per_case"] == 16_384
    @test qualification["sigma_limit"] == 6.0
    @test Set(case["id"] for case in qualification["moment_cases"]) ==
        Set(("shot", "dark", "clock_induced_charge", "read_noise",
            "combined"))
    @test all(case -> case["mean_passed"] && case["variance_passed"],
        qualification["moment_cases"])
    deterministic = qualification["deterministic"]
    @test deterministic["cic_exposure_invariant"]
    @test deterministic["single_read_read_time_rejected"]
    @test deterministic["allocation_gate_passed"]
    @test deterministic["steady_alloc_bytes"] == 0

    emccd_entry = detector_entries["DET-EMCCD-QUAL-2026-07-30"]
    @test emccd_entry["status"] == "active"
    emccd_path = joinpath(dirname(detector_artifact_path),
        emccd_entry["path"])
    @test isfile(emccd_path)
    emccd = TOML.parsefile(emccd_path)
    @test emccd["schema_version"] == 1
    @test emccd["family"] == "emccd"
    @test emccd["all_gates_passed"]
    @test !emccd["environment"]["source_dirty"]
    qualification = emccd["qualification"]
    @test qualification["samples_per_case"] == 16_384
    @test qualification["sigma_limit"] == 6.0
    moment_cases = qualification["multiplication_moment_cases"]
    @test Set(case["id"] for case in moment_cases) == Set((
        "conditional_gamma_exponential",
        "conditional_gamma_erlang",
        "conditional_gamma_fractional_shape",
        "clipped_gaussian_moderate_charge",
        "cic_before_conditional_gamma",
    ))
    @test all(case -> case["mean_passed"] && case["variance_passed"],
        moment_cases)
    @test only(filter(case ->
        case["id"] == "clipped_gaussian_moderate_charge",
        moment_cases))["all_nonnegative"]
    @test all(case -> case["passed"],
        qualification["multiplication_cdf_cases"])
    photon_counting_cases = qualification["photon_counting_cases"]
    @test Set(case["incident_mean"] for case in photon_counting_cases) ==
        Set((0.05, 0.5, 2.0))
    @test all(case -> case["mean_passed"] &&
        case["variance_passed"] && case["all_binary"] &&
        case["coincidence_limited"], photon_counting_cases)
    emccd_deterministic = qualification["deterministic"]
    @test all(emccd_deterministic[key] for key in (
        "cic_exposure_invariant",
        "input_full_well_passed",
        "register_full_well_passed",
        "conventional_output_passed",
        "gain_range_rejected",
        "accelerator_conditional_gamma_rejected",
        "frame_transfer_preserves_optical_output",
        "timing_contract_passed",
        "snr_contract_passed",
        "allocation_gate_passed",
    ))
    @test emccd_deterministic["steady_alloc_bytes"] == 0

    cmos_entry = detector_entries["DET-CMOS-QUAL-2026-07-30"]
    @test cmos_entry["status"] == "active"
    cmos_path = joinpath(dirname(detector_artifact_path),
        cmos_entry["path"])
    @test isfile(cmos_path)
    cmos = TOML.parsefile(cmos_path)
    @test cmos["schema_version"] == 1
    @test cmos["family"] == "parameterized_cmos"
    @test cmos["all_gates_passed"]
    @test !cmos["environment"]["source_dirty"]
    cmos_qualification = cmos["qualification"]
    @test cmos_qualification["samples_per_case"] == 16_384
    @test cmos_qualification["sigma_limit"] == 6.0
    covariance_cases = cmos_qualification["spatial_noise_cases"]
    @test Set(case["id"] for case in covariance_cases) == Set((
        "column_common",
        "row_common",
        "pixel_independent",
        "combined",
    ))
    @test all(case -> all(values(case["gates"])), covariance_cases)
    cmos_deterministic = cmos_qualification["deterministic"]
    @test all(cmos_deterministic[key] for key in (
        "global_shutter_passed",
        "rolling_exposure_passed",
        "global_reset_passed",
        "window_preserves_full_frame_timing",
        "configured_mtf_preserved",
        "readout_pipeline_passed",
        "deterministic_replay_passed",
        "allocation_gate_passed",
    ))
    @test cmos_deterministic["steady_alloc_bytes"] == 0
end

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

    hgcdte_entry = detector_entries["DET-HGCDTE-QUAL-2026-07-30"]
    @test hgcdte_entry["status"] == "active"
    hgcdte_path = joinpath(dirname(detector_artifact_path),
        hgcdte_entry["path"])
    @test isfile(hgcdte_path)
    hgcdte = TOML.parsefile(hgcdte_path)
    @test hgcdte["schema_version"] == 1
    @test hgcdte["family"] == "conventional_hgcdte"
    @test hgcdte["all_gates_passed"]
    @test !hgcdte["environment"]["source_dirty"]
    hgcdte_qualification = hgcdte["qualification"]
    @test hgcdte_qualification["samples_per_case"] == 16_384
    @test hgcdte_qualification["sigma_limit"] == 6.0
    moment_cases = hgcdte_qualification["moment_cases"]
    @test Set(case["id"] for case in moment_cases) == Set((
        "single_read_noise",
        "averaged_ndr4_noise",
        "correlated_double_sampling_noise",
        "fowler8_noise",
        "up_the_ramp16_noise",
        "dark_current",
        "readout_glow",
    ))
    @test all(case -> case["mean_passed"] && case["variance_passed"],
        moment_cases)
    hgcdte_deterministic = hgcdte_qualification["deterministic"]
    @test all(hgcdte_deterministic[key] for key in (
        "architecture_separated",
        "single_ndr_cds_fowler_passed",
        "direct_synthesized_ramp_passed",
        "irregular_scheduled_ramp_passed",
        "window_preserves_full_frame_timing",
        "configured_mtf_preserved",
        "configured_ipc_passed",
        "reference_correction_passed",
        "saturation_passed",
        "nonlinearity_passed",
        "persistence_passed",
        "deterministic_replay_passed",
        "allocation_gate_passed",
    ))
    @test hgcdte_deterministic["steady_alloc_bytes"] == 0

    avalanche_entry =
        detector_entries["DET-HGCDTE-AVALANCHE-QUAL-2026-07-31"]
    @test avalanche_entry["status"] == "active"
    avalanche_path = joinpath(
        dirname(detector_artifact_path), avalanche_entry["path"])
    @test isfile(avalanche_path)
    avalanche = TOML.parsefile(avalanche_path)
    @test avalanche["schema_version"] == 1
    @test avalanche["family"] ==
        "hgcdte_linear_avalanche_photodiode_array"
    @test avalanche["all_gates_passed"]
    @test !avalanche["environment"]["source_dirty"]
    @test avalanche["model"]["excess_noise_factor_definition"] ==
        "F = E[M^2] / E[M]^2"
    avalanche_qualification = avalanche["qualification"]
    @test avalanche_qualification["samples_per_case"] == 16_384
    @test avalanche_qualification["sigma_limit"] == 6.0
    avalanche_moments =
        avalanche_qualification["multiplication_moment_cases"]
    @test Set(case["id"] for case in avalanche_moments) == Set((
        "conditional_gamma_single_carrier",
        "conditional_gamma_multiple_carriers",
        "conditional_gamma_low_excess_noise",
        "clipped_gaussian_moderate_charge",
    ))
    @test all(case ->
        case["mean_passed"] &&
        case["variance_passed"] &&
        case["all_nonnegative"] &&
        case["qualified_regime_passed"],
        avalanche_moments)
    avalanche_deterministic =
        avalanche_qualification["deterministic"]
    @test all(avalanche_deterministic[key] for key in (
        "architecture_separated",
        "exact_gain_passed",
        "input_referred_saturation_passed",
        "generated_charge_ordering_passed",
        "read_noise_and_conversion_gain_ordering_passed",
        "single_ndr_cds_fowler_ramp_passed",
        "window_preserves_full_frame_timing",
        "configured_mtf_preserved",
        "configured_ipc_passed",
        "scheduled_retains_prior_multiplication",
        "gamma_accelerator_rejected",
        "deterministic_replay_passed",
        "allocation_gate_passed",
    ))
    @test avalanche_deterministic[
        "approximate_steady_alloc_bytes"] == 0
    @test avalanche_deterministic[
        "gamma_steady_alloc_bytes"] == 0
    @test avalanche_deterministic[
        "scheduled_steady_alloc_bytes"] == 0

    skipper_entry = detector_entries["DET-SKIPPER-QUAL-2026-07-31"]
    @test skipper_entry["status"] == "active"
    skipper_path = joinpath(
        dirname(detector_artifact_path), skipper_entry["path"])
    @test isfile(skipper_path)
    skipper = TOML.parsefile(skipper_path)
    @test skipper["schema_version"] == 1
    @test skipper["family"] == "skipper_ccd_independent_read"
    @test skipper["all_gates_passed"]
    @test !skipper["environment"]["source_dirty"]
    skipper_qualification = skipper["qualification"]
    @test skipper_qualification["samples_per_case"] == 16_384
    @test skipper_qualification["sigma_limit"] == 6.0
    skipper_cases = skipper_qualification["moment_cases"]
    @test Set(case["n_samples"] for case in skipper_cases) ==
        Set((1, 4, 16, 64))
    @test all(case ->
        case["mean_passed"] && case["variance_passed"],
        skipper_cases)
    skipper_deterministic = skipper_qualification["deterministic"]
    @test all(skipper_deterministic[key] for key in (
        "exact_mean_and_gain_passed",
        "retained_charge_packet_passed",
        "input_referred_full_well_passed",
        "fixed_frame_storage_passed",
        "sample_count_passed",
        "timing_passed",
        "batched_capture_rejected",
        "batched_input_unmodified",
        "deterministic_replay_passed",
        "allocation_gate_passed",
    ))
    @test skipper_deterministic["steady_alloc_bytes"] == 0

    ingaas_entry = detector_entries["DET-INGAAS-QUAL-2026-07-31"]
    @test ingaas_entry["status"] == "active"
    ingaas_path = joinpath(
        dirname(detector_artifact_path), ingaas_entry["path"])
    @test isfile(ingaas_path)
    ingaas = TOML.parsefile(ingaas_path)
    @test ingaas["schema_version"] == 1
    @test ingaas["family"] == "ingaas_frame_detector"
    @test ingaas["all_gates_passed"]
    @test !ingaas["environment"]["source_dirty"]
    @test length(ingaas["environment"]["source_revision"]) == 40
    @test ingaas["model"]["default_response"] == "none"
    @test ingaas["model"]["coefficient_domain"] ==
        "dimensionless per completed frame"
    ingaas_qualification = ingaas["qualification"]
    @test ingaas_qualification["samples_per_case"] == 16_384
    @test ingaas_qualification["sigma_limit"] == 6.0
    ingaas_poisson = ingaas_qualification["poisson_case"]
    @test ingaas_poisson["id"] == "combined_glow_and_dark_current"
    @test ingaas_poisson["mean_passed"]
    @test ingaas_poisson["variance_passed"]
    ingaas_deterministic = ingaas_qualification["deterministic"]
    @test all(ingaas_deterministic[key] for key in (
        "default_response_is_null",
        "default_response_preserves_impulse",
        "default_mtf_is_not_claimed",
        "pipeline_order_passed",
        "persistence_recurrence_passed",
        "persistence_gain_independence_passed",
        "deterministic_replay_passed",
        "batched_persistence_rejected",
        "batched_input_unmodified",
        "allocation_gate_passed",
    ))
    @test ingaas_deterministic["steady_alloc_bytes"] == 0

    linear_apd_entry =
        detector_entries["DET-LINEAR-APD-QUAL-2026-07-31"]
    @test linear_apd_entry["status"] == "active"
    linear_apd_path = joinpath(
        dirname(detector_artifact_path), linear_apd_entry["path"])
    @test isfile(linear_apd_path)
    linear_apd = TOML.parsefile(linear_apd_path)
    @test linear_apd["schema_version"] == 1
    @test linear_apd["family"] == "linear_mode_apd_channels"
    @test linear_apd["all_gates_passed"]
    @test !linear_apd["environment"]["source_dirty"]
    @test linear_apd["model"]["operating_regime"] == "linear mode"
    @test linear_apd["model"]["excess_noise_factor_definition"] ==
        "F = E[M^2] / E[M]^2"
    linear_apd_qualification = linear_apd["qualification"]
    @test linear_apd_qualification["samples_per_case"] == 16_384
    @test linear_apd_qualification["sigma_limit"] == 6.0
    linear_apd_moments = linear_apd_qualification["moment_cases"]
    @test Set(case["id"] for case in linear_apd_moments) == Set((
        "multiplied_shot", "multiplication_only", "read_only"))
    @test all(case -> case["mean_passed"] && case["variance_passed"],
        linear_apd_moments)
    linear_apd_deterministic = linear_apd_qualification["deterministic"]
    @test all(linear_apd_deterministic[key] for key in (
        "ambiguous_generic_apd_removed",
        "single_element_vector_storage",
        "exact_signal_order_passed",
        "channel_bank_vector_storage",
        "topology_metadata_passed",
        "matrix_input_rejected",
        "deterministic_replay_passed",
        "allocation_gate_passed",
    ))
    @test linear_apd_deterministic["steady_alloc_bytes"] == 0

    spad_entry = detector_entries["DET-SPAD-QUAL-2026-08-01"]
    @test spad_entry["status"] == "active"
    spad_path = joinpath(
        dirname(detector_artifact_path), spad_entry["path"])
    @test isfile(spad_path)
    spad = TOML.parsefile(spad_path)
    @test spad["schema_version"] == 1
    @test spad["family"] == "spad_accumulated_count_array"
    @test spad["all_gates_passed"]
    @test !spad["environment"]["source_dirty"]
    @test spad["environment"]["source_revision"] ==
        "75e2fdab1d502e0e096f2751d677db5137fedea8"
    @test spad["model"]["statistical_scope"] ==
        "Poisson draw from adjusted mean; not the exact dead-time or afterpulse count distribution"
    spad_qualification = spad["qualification"]
    @test spad_qualification["samples_per_moment_case"] == 16_384
    @test spad_qualification["sigma_limit"] == 6.0
    spad_curves = spad_qualification["dead_time_curves"]
    @test Set(curve["id"] for curve in spad_curves) ==
        Set(("no_dead_time", "nonparalyzable", "paralyzable"))
    @test all(curve -> curve["all_points_passed"], spad_curves)
    @test all(curve ->
        Set(point["lambda_tau"] for point in curve["points"]) ==
            Set((0.0, 1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0)),
        spad_curves)
    spad_poisson_cases =
        spad_qualification["poisson_surrogate_cases"]
    @test Set(case["id"] for case in spad_poisson_cases) == Set((
        "photon_only", "dark_only", "photon_and_dark",
        "dead_time_adjusted_surrogate"))
    @test all(case -> case["mean_passed"] && case["variance_passed"],
        spad_poisson_cases)
    spad_deterministic = spad_qualification["deterministic"]
    @test all(spad_deterministic[key] for key in (
        "exact_radiometry_and_gate_passed",
        "gated_dark_live_time_passed",
        "temperature_dependent_dark_count_passed",
        "first_order_afterpulse_mean_passed",
        "ordered_live_time_pipeline_passed",
        "afterpulse_metadata_passed",
        "redistribution_center_passed",
        "redistribution_conserves_counts",
        "integer_rounding_and_saturation_passed",
        "fixed_shape_mismatch_rejected",
        "fixed_shape_storage_preserved",
        "invalid_input_rejected",
        "detector_mtf_not_applicable",
        "deterministic_replay_passed",
        "allocation_gate_passed",
    ))
    @test spad_deterministic["steady_alloc_bytes"] == 0
end

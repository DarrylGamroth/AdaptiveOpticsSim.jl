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
end

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

end

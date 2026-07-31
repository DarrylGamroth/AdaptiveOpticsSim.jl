using Aqua
using SHA

canonical_lf_sha256(text::AbstractString) = bytes2hex(SHA.sha256(
    codeunits(replace(text, "\r\n" => "\n")),
))

_rng_family(::Xoshiro) = :xoshiro
_rng_family(::MersenneTwister) = :mersenne_twister
_rng_family(::AbstractRNG) = :other

@testset "RNG policy helpers" begin
    @test _rng_family(runtime_rng(1)) == :xoshiro
    @test _rng_family(deterministic_reference_rng(1)) == :mersenne_twister
    @test rand(runtime_rng(42), UInt64) == rand(runtime_rng(42), UInt64)
    @test rand(deterministic_reference_rng(42), UInt64) ==
        rand(deterministic_reference_rng(42), UInt64)
    @test coverage_runner_flag("true")
    @test coverage_runner_flag(" YES ")
    @test !coverage_runner_flag("false")
    @test coverage_instrumented() ==
        (coverage_runner_flag(get(ENV, "ADAPTIVEOPTICS_TEST_COVERAGE",
            "false")) || Base.JLOptions().code_coverage != 0)
end

@test AdaptiveOpticsSim.PROJECT_STATUS == :in_development

@testset "Gate 7 benchmark contract" begin
    root = normpath(joinpath(@__DIR__, "..", ".."))
    contract_path =
        joinpath(root, "benchmarks", "contracts", "gate7_single_gpu.toml")
    harness_path =
        joinpath(root, "benchmarks", "benchmark_gate7_single_gpu.jl")
    support_path =
        joinpath(root, "benchmarks", "support", "gate7_single_gpu.jl")
    @test isfile(harness_path)
    @test isfile(support_path)

    contract = TOML.parsefile(contract_path)
    @test contract["schema_version"] == 1
    @test contract["name"] == "gate7_single_gpu"
    @test contract["samples_per_run"] >=
        contract["minimum_samples_for_p95"]
    @test contract["runs"] >= 3
    @test contract["warmup_operations"] > 0
    @test contract["batched_relative_p95_factor"] >= 1
    @test canonical_lf_sha256("gate7\r\nartifact\r\n") ==
        canonical_lf_sha256("gate7\nartifact\n")

    workload = contract["workload"]
    @test workload["path_count"] == 2
    @test workload["numeric_type"] == "Float32"
    @test workload["witness_phase_ns"] > workload["sample_period_ns"]

    independent = contract["submission_proxy"]["independent"]
    batched = contract["submission_proxy"]["batched"]
    @test independent["top_level_path_submissions"] == 2
    @test independent["device_owner_submissions"] == 0
    @test batched["top_level_path_submissions"] == 0
    @test batched["device_owner_submissions"] == 1
    @test batched["atmosphere_direction_render_calls"] <
        independent["atmosphere_direction_render_calls"]
    @test batched["wfs_formation_calls"] ==
        independent["wfs_formation_calls"]

    gpu_boundaries = [
        "independent_device_ready",
        "batched_device_ready",
        "batched_host_ready",
        "transfer_only",
    ]
    placements = contract["placements"]
    @test placements["local_cpu"]["boundaries"] ==
        ["independent_device_ready"]
    @test placements["local_amdgpu"]["boundaries"] == gpu_boundaries
    @test placements["wsl_cuda"]["boundaries"] == gpu_boundaries
    @test contract["contract"]["coordinated_omission_correction"] ===
        false
    @test any(
        exclusion -> occursin("HIL latency", exclusion),
        contract["scope_exclusions"],
    )

    artifact_root = joinpath(root, "benchmarks", "results", "gate7")
    manifest = TOML.parsefile(joinpath(artifact_root, "manifest.toml"))
    closure = manifest["closure"]
    artifacts = manifest["artifacts"]
    @test closure["status"] == "passed"
    @test closure["requirements"] == ["HIL-GPU-001"]
    @test length(artifacts) == 3
    @test Set(artifact["backend"] for artifact in artifacts) ==
        Set(("cpu", "amdgpu", "cuda"))

    for entry in artifacts
        artifact_path = joinpath(artifact_root, entry["path"])
        @test isfile(artifact_path)
        @test canonical_lf_sha256(read(artifact_path, String)) ==
            entry["sha256"]
        artifact = TOML.parsefile(artifact_path)
        @test artifact["all_gates_passed"]
        @test artifact["p95_supported"]
        @test artifact["characterized_source_revision"] ==
            closure["characterized_source_revision"]
        @test artifact["environment"]["git_dirty"] === false
        @test artifact["environment"]["backend"] == entry["backend"]
        @test artifact["correctness"]["passed"]
        @test artifact["submission_proxy_evidence"]["passed"]
        @test artifact["relative_comparison"]["passed"]
        @test all(artifact["boundaries"]) do boundary
            boundary["passed"] &&
                length(boundary["runs"]) ==
                    artifact["configured_runs"] &&
                all(boundary["runs"]) do run
                    run["samples"] ==
                        artifact["configured_samples_per_run"] &&
                        !isempty(run["histogram_base64"])
                end
        end
    end
end

@testset "Selective test registry" begin
    @test resolve_test_suites(String[]) === TEST_SUITE_SPECS
    @test resolve_test_suites(["all"]) === TEST_SUITE_SPECS
    @test resolve_test_suites(("all", "all")) === TEST_SUITE_SPECS
    @test isnothing(resolve_test_suites(["--list"]))
    @test isnothing(resolve_test_suites(("--list", "--list")))
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["plant-topology"])) == ("plant-topology",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["plant-topology-growth"])) == ("plant-topology-growth",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["plant-command-schemas"])) == ("plant-command-schemas",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["plant-command-admission"])) == ("plant-command-admission",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["plant-command-application"])) == ("plant-command-application",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["plant-device-batching"])) == ("plant-device-batching",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["plant-device-model-matrix"])) ==
        ("plant-device-model-matrix",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["plant-command-composition"])) == ("plant-command-composition",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["plant-controller-routing"])) == ("plant-controller-routing",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["plant-autonomous-optics"])) == ("plant-autonomous-optics",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["plant-placed-optics"])) == ("plant-placed-optics",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["plant-conjugate-geometry"])) ==
        ("plant-conjugate-geometry",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["plant-mcao-moao"])) == ("plant-mcao-moao",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["plant-sampled-aberrations"])) ==
        ("plant-sampled-aberrations",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["plant-gate5-closure"])) ==
        ("plant-gate5-closure",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["plant-time"])) == ("plant-time",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["detectors"])) == DETECTOR_TEST_SUITE_NAMES
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["wfs-common"])) == ("wfs-common",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["wfs-shack-hartmann"])) == ("wfs-shack-hartmann",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["wfs-pyramid-bioedge"])) == ("wfs-pyramid-bioedge",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["wfs-zernike-curvature"])) ==
        ("wfs-zernike-curvature",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["wfs-lift"])) == ("wfs-lift",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["calibration"])) == ("calibration-workflows",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["control"])) ==
        ("control-primitives", "control-reconstruction")
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["sensors"])) == (
            DETECTOR_TEST_SUITE_NAMES...,
            "wfs-common",
            "wfs-shack-hartmann",
            "wfs-pyramid-bioedge",
            "wfs-zernike-curvature",
            "wfs-lift",
        )
    @test Tuple(spec.name for spec in resolve_test_suites(["gate4"])) == (
        "plant-command-schemas",
        "plant-command-admission",
        "plant-command-application",
        "plant-event-composition",
        "plant-command-composition",
        "plant-controller-routing",
        "plant-reduced-order",
        "plant-autonomous-optics",
    )
    @test Tuple(spec.name for spec in resolve_test_suites(["gate5"])) ==
        (
            "plant-placed-optics",
            "plant-conjugate-geometry",
            "plant-mcao-moao",
            "plant-sampled-aberrations",
            "plant-gate5-closure",
        )
    @test Tuple(spec.name for spec in resolve_test_suites(["gate7"])) ==
        (
            "direct-imaging-batch",
            "atmosphere-direction-batch",
            "plant-device-batching",
            "plant-device-model-matrix",
            "backend-smoke",
        )
    @test Tuple(spec.name for spec in resolve_test_suites(["plant"])) == (
        "plant-topology",
        "plant-topology-growth",
        "plant-command-schemas",
        "plant-command-admission",
        "plant-command-application",
        "plant-time",
        "plant-scheduler",
        "plant-triggers",
        "plant-detector-transitions",
        "plant-event-composition",
        "plant-device-batching",
        "plant-device-model-matrix",
        "plant-cpu-execution",
        "plant-command-composition",
        "plant-controller-routing",
        "plant-reduced-order",
        "plant-autonomous-optics",
        "plant-placed-optics",
        "plant-conjugate-geometry",
        "plant-mcao-moao",
        "plant-sampled-aberrations",
        "plant-gate5-closure",
        "plant-preparation",
        "plant-providers",
        "plant-rng",
        "plant-illumination",
    )
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["plant", "plant-topology", "plant"])) == (
        "plant-topology",
        "plant-topology-growth",
        "plant-command-schemas",
        "plant-command-admission",
        "plant-command-application",
        "plant-time",
        "plant-scheduler",
        "plant-triggers",
        "plant-detector-transitions",
        "plant-event-composition",
        "plant-device-batching",
        "plant-device-model-matrix",
        "plant-cpu-execution",
        "plant-command-composition",
        "plant-controller-routing",
        "plant-reduced-order",
        "plant-autonomous-optics",
        "plant-placed-optics",
        "plant-conjugate-geometry",
        "plant-mcao-moao",
        "plant-sampled-aberrations",
        "plant-gate5-closure",
        "plant-preparation",
        "plant-providers",
        "plant-rng",
        "plant-illumination",
    )
    @test_throws ArgumentError resolve_test_suites(["unknown"])
    @test_throws ArgumentError resolve_test_suites(["all", "quality"])
    @test_throws ArgumentError resolve_test_suites(["--list", "quality"])
    @test isnothing(validate_test_suite_registry())

    duplicate_name_specs = (
        TEST_SUITE_SPECS...,
        TestSuiteSpec("quality", "testsets/quality.jl"),
    )
    @test_throws ArgumentError validate_test_suite_registry(
        duplicate_name_specs, ())
    duplicate_path_specs = (
        TEST_SUITE_SPECS...,
        TestSuiteSpec("duplicate-path", "testsets/quality.jl"),
    )
    @test_throws ArgumentError validate_test_suite_registry(
        duplicate_path_specs, ())
    missing_path_specs = (
        TestSuiteSpec("missing-path", "testsets/not_present.jl"),
    )
    @test_throws ArgumentError validate_test_suite_registry(
        missing_path_specs, ())
    minimal_specs = (
        TestSuiteSpec("quality-only", "testsets/quality.jl"),
    )
    @test isnothing(validate_test_suite_registry(minimal_specs, ()))
    incomplete_specs = Tuple(
        spec for spec in TEST_SUITE_SPECS
        if spec.name != "quality")
    @test_throws ArgumentError validate_test_suite_registry(
        incomplete_specs, (); require_complete=true)

    duplicated_ci_shards = (
        "ci-first" => test_suite_names(),
        "ci-second" => ("quality",),
    )
    @test_throws ArgumentError validate_test_suite_registry(
        TEST_SUITE_SPECS,
        duplicated_ci_shards;
        ci_shard_specs=duplicated_ci_shards,
    )
    incomplete_ci_shards = (
        "ci-incomplete" => Base.front(test_suite_names()),
    )
    @test_throws ArgumentError validate_test_suite_registry(
        TEST_SUITE_SPECS,
        incomplete_ci_shards;
        ci_shard_specs=incomplete_ci_shards,
    )
    ci_members = Tuple(
        member
        for (_, members) in TEST_CI_SHARD_SPECS
        for member in members)
    @test length(ci_members) == length(test_suite_names())
    @test Set(ci_members) == Set(test_suite_names())

    listing = IOBuffer()
    @test isnothing(print_test_suite_help(listing))
    listing_text = String(take!(listing))
    @test occursin("plant-topology", listing_text)
    @test occursin("plant-topology-growth", listing_text)
    @test occursin("plant-command-schemas", listing_text)
    @test occursin("plant-command-admission", listing_text)
    @test occursin("plant-command-application", listing_text)
    @test occursin("plant-device-batching", listing_text)
    @test occursin("plant-autonomous-optics", listing_text)
    @test occursin("plant-placed-optics", listing_text)
    @test occursin("plant-conjugate-geometry", listing_text)
    @test occursin("plant-mcao-moao", listing_text)
    @test occursin("plant-sampled-aberrations", listing_text)
    @test occursin("plant-gate5-closure", listing_text)
    @test occursin("plant-cpu-execution", listing_text)
    @test occursin("plant-time", listing_text)
    @test occursin("plant =", listing_text)
    @test occursin("sensors =", listing_text)
    @test occursin("calibration =", listing_text)
    @test occursin("ci-foundations =", listing_text)
    @test occursin("ci-sensors-control =", listing_text)
    @test occursin("ci-plant-runtime =", listing_text)
    @test occursin("ci-plant-optics =", listing_text)
    @test occursin("gate4 =", listing_text)
    @test occursin("gate5 =", listing_text)
    @test occursin("gate6 =", listing_text)
    @test occursin("gate7 =", listing_text)

    actual_testsets = sort!(filter(
        path -> endswith(path, ".jl"),
        readdir(@__DIR__; join=true),
    ))
    @test registered_testset_paths() == normpath.(actual_testsets)
    @test registered_test_fixture_paths() == sort!(normpath.([
        joinpath(dirname(@__DIR__), "detector_test_fixtures.jl"),
        joinpath(dirname(@__DIR__), "ka_cpu_style_fixture.jl"),
        joinpath(dirname(@__DIR__), "plant_device_batching_fixtures.jl"),
        joinpath(
            dirname(@__DIR__),
            "plant_device_model_matrix_fixtures.jl",
        ),
        joinpath(dirname(@__DIR__), "plant_test_fixtures.jl"),
        joinpath(dirname(@__DIR__), "wfs_stage_contract_fixtures.jl"),
    ]))
    fixture_users = Tuple(spec.name for spec in TEST_SUITE_SPECS
        if !isempty(spec.fixtures))
    @test fixture_users == (
        "ka-cpu",
        "direct-imaging-batch",
        "atmosphere-direction-batch",
        "plant-device-batching",
        "plant-device-model-matrix",
        "plant-sampled-aberrations",
        DETECTOR_TEST_SUITE_NAMES...,
        "wfs-common",
        "wfs-shack-hartmann",
        "wfs-pyramid-bioedge",
        "wfs-zernike-curvature",
        "plant-preparation",
        "plant-providers",
        "plant-rng",
        "backend-smoke",
    )
end

@testset "Aqua" begin
    Aqua.test_all(
        AdaptiveOpticsSim;
        undocumented_names=false,
    )
end

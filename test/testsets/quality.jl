using Aqua

_rng_family(::SplitMix64RNG) = :splitmix64
_rng_family(::MersenneTwister) = :mersenne_twister
_rng_family(::AbstractRNG) = :other

@testset "RNG policy helpers" begin
    @test _rng_family(runtime_rng(1)) == :splitmix64
    @test _rng_family(deterministic_reference_rng(1)) == :mersenne_twister
    @test runtime_rng(1) isa SplitMix64RNG
    @test rand(runtime_rng(1), UInt64) == 0x910a2dec89025cc1
    @test rand(runtime_rng(42), UInt64) == rand(runtime_rng(42), UInt64)
    splitmix = runtime_rng(42)
    splitmix_copy = copy(splitmix)
    @test splitmix_copy == splitmix
    @test rand(splitmix_copy, UInt64) == rand(splitmix, UInt64)
    Random.seed!(splitmix, 42)
    @test splitmix == runtime_rng(42)
    @test randn(runtime_rng(42)) == randn(runtime_rng(42))
    @test rand(runtime_rng(42), 1:10) == rand(runtime_rng(42), 1:10)
    @test_throws InvalidConfiguration SplitMix64RNG(true)
    @test_throws InvalidConfiguration SplitMix64RNG(-1)
    @test_throws InvalidConfiguration SplitMix64RNG(
        big(typemax(UInt64)) + 1,
    )
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

@testset "Implementation storage policy" begin
    repository = pkgdir(AdaptiveOpticsSim)
    implementation_paths = String[]
    for directory in ("src", "ext")
        for (root, _, files) in walkdir(joinpath(repository, directory))
            for file in files
                endswith(file, ".jl") || continue
                push!(implementation_paths, joinpath(root, file))
            end
        end
    end
    sort!(implementation_paths)

    forbidden = (
        "direct Memory storage" => r"\bMemory\s*(?:\{|\()",
        "Any collection construction" =>
            r"\b(?:Array|Vector|Matrix|Dict|Set|Ref|FixedSizeArray|FixedSizeVector)\s*\{[^}\n]*\bAny\b[^\n]*\}\s*\(",
        "Any array literal" => r"\bAny\s*\[",
        "Any field or collection-field storage" =>
            r"(?m)^\s*[A-Za-z_]\w*\s*::\s*(?:Any|(?:Array|Vector|Matrix|Dict|Set|Ref|FixedSizeArray|FixedSizeVector)\s*\{[^\n}]*\bAny\b[^\n]*\})\s*$",
    )
    violations = String[]
    for path in implementation_paths
        source = read(path, String)
        relative = relpath(path, repository)
        for (label, pattern) in forbidden
            occursin(pattern, source) || continue
            push!(violations, "$relative: $label")
        end
    end
    @test isempty(violations)
end

@testset "Selective test registry" begin
    @test resolve_test_suites(String[]) === TEST_SUITE_SPECS
    @test resolve_test_suites(["all"]) === TEST_SUITE_SPECS
    @test resolve_test_suites(("all", "all")) === TEST_SUITE_SPECS
    @test isnothing(resolve_test_suites(["--list"]))
    @test isnothing(resolve_test_suites(("--list", "--list")))
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["detectors"])) == DETECTOR_TEST_SUITE_NAMES
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["wfs-acquisition-ownership"])) ==
        ("wfs-acquisition-ownership",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["wfs-common"])) == ("wfs-common",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["wfs-shack-hartmann"])) == ("wfs-shack-hartmann",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["wfs-pyramid-bi-o-edge"])) == ("wfs-pyramid-bi-o-edge",)
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
            "wfs-acquisition-ownership",
            "wfs-common",
            "wfs-shack-hartmann",
            "wfs-pyramid-bi-o-edge",
            "wfs-zernike-curvature",
            "wfs-lift",
        )
    @test_throws ArgumentError resolve_test_suites(["plant"])
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
    @test occursin("algorithm-graphs", listing_text)
    @test occursin("sensors =", listing_text)
    @test occursin("calibration =", listing_text)
    @test occursin("ci-foundations =", listing_text)
    @test occursin("ci-sensors-control =", listing_text)

    actual_testsets = sort!(filter(
        path -> endswith(path, ".jl"),
        readdir(@__DIR__; join=true),
    ))
    @test registered_testset_paths() == normpath.(actual_testsets)
    @test registered_test_fixture_paths() == sort!(normpath.([
        joinpath(
            dirname(@__DIR__),
            "backend_execution_context_conformance.jl",
        ),
        joinpath(
            dirname(@__DIR__),
            "calibration_interface_conformance.jl",
        ),
        joinpath(dirname(@__DIR__), "detector_test_fixtures.jl"),
        joinpath(dirname(@__DIR__), "ka_cpu_style_fixture.jl"),
        joinpath(
            dirname(@__DIR__),
            "wfs_four_pupil_interface_conformance.jl",
        ),
        joinpath(dirname(@__DIR__), "wfs_stage_contract_fixtures.jl"),
    ]))
    fixture_users = Tuple(spec.name for spec in TEST_SUITE_SPECS
        if !isempty(spec.fixtures))
    @test fixture_users == (
        "ka-cpu",
        "prepared-execution-interfaces",
        "direct-imaging-batch",
        "atmosphere-direction-batch",
        DETECTOR_TEST_SUITE_NAMES...,
        "wfs-common",
        "wfs-shack-hartmann",
        "wfs-pyramid-bi-o-edge",
        "wfs-zernike-curvature",
        "backend-smoke",
    )
end

@testset "Aqua" begin
    Aqua.test_all(
        AdaptiveOpticsSim;
        undocumented_names=false,
    )
end

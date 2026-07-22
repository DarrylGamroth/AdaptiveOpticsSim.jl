using Aqua

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

@testset "Selective test registry" begin
    @test resolve_test_suites(String[]) === TEST_SUITE_SPECS
    @test resolve_test_suites(["all"]) === TEST_SUITE_SPECS
    @test resolve_test_suites(("all", "all")) === TEST_SUITE_SPECS
    @test isnothing(resolve_test_suites(["--list"]))
    @test isnothing(resolve_test_suites(("--list", "--list")))
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["plant-topology"])) == ("plant-topology",)
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["plant-time"])) == ("plant-time",)
    @test Tuple(spec.name for spec in resolve_test_suites(["plant"])) == (
        "plant-topology",
        "plant-time",
        "plant-scheduler",
        "plant-triggers",
        "plant-detector-transitions",
        "plant-preparation",
        "plant-providers",
        "plant-rng",
        "plant-illumination",
    )
    @test Tuple(spec.name for spec in resolve_test_suites(
        ["plant", "plant-topology", "plant"])) == (
        "plant-topology",
        "plant-time",
        "plant-scheduler",
        "plant-triggers",
        "plant-detector-transitions",
        "plant-preparation",
        "plant-providers",
        "plant-rng",
        "plant-illumination",
    )
    @test_throws ArgumentError resolve_test_suites(["unknown"])
    @test_throws ArgumentError resolve_test_suites(["all", "quality"])
    @test_throws ArgumentError resolve_test_suites(["--list", "quality"])

    listing = IOBuffer()
    @test isnothing(print_test_suite_help(listing))
    listing_text = String(take!(listing))
    @test occursin("plant-topology", listing_text)
    @test occursin("plant-time", listing_text)
    @test occursin("plant =", listing_text)

    actual_testsets = sort!(filter(
        path -> endswith(path, ".jl"),
        readdir(@__DIR__; join=true),
    ))
    @test registered_testset_paths() == normpath.(actual_testsets)
    @test registered_test_fixture_paths() == sort!(normpath.([
        joinpath(dirname(@__DIR__), "ka_cpu_style_fixture.jl"),
        joinpath(dirname(@__DIR__), "wfs_stage_contract_fixtures.jl"),
    ]))
    fixture_users = Tuple(spec.name for spec in TEST_SUITE_SPECS
        if !isempty(spec.fixtures))
    @test fixture_users == (
        "ka-cpu",
        "detectors-wfs",
        "plant-preparation",
        "plant-providers",
        "backend-smoke",
    )
end

@testset "Aqua" begin
    Aqua.test_all(
        AdaptiveOpticsSim;
        undocumented_names=false,
    )
end

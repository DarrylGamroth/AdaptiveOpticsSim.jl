#!/usr/bin/env julia

using Test

include(joinpath(@__DIR__, "impact_planner.jl"))

@testset "local validation impact planner" begin
    @test Set(CLOSURE_SELECTORS) == Set(first.(TEST_CI_SHARD_SPECS))

    docs = plan_validation(("docs/user-guide.md",))
    @test docs.selectors == Set(("quality",))
    @test docs.run_cpu
    @test !docs.run_cuda
    @test !docs.run_amdgpu
    @test isempty(docs.manual_gates)

    glossary = plan_validation(("docs/glossary.md",))
    @test glossary.selectors == Set(("quality",))

    cmos = plan_validation(("src/detectors/cmos.jl",))
    @test "detector-cmos" in cmos.selectors
    @test "detector-shared" in cmos.selectors
    @test "detector-thermal" in cmos.selectors
    @test cmos.accelerator_scope == DetectorAccelerator
    @test cmos.run_cuda
    @test cmos.run_amdgpu
    @test !cmos.run_aarch64

    ccd = plan_validation(("src/detectors/ccd.jl",))
    @test "detector-ccd" in ccd.selectors
    @test "detector-skipper" in ccd.selectors

    detector_common = plan_validation(("src/detectors/pipeline.jl",))
    @test "detectors" in detector_common.selectors
    @test detector_common.accelerator_scope == DetectorAccelerator

    shack_hartmann = plan_validation(("src/wfs/shack_hartmann.jl",))
    @test shack_hartmann.selectors == Set((
        "wfs-acquisition-ownership", "wfs-common", "wfs-shack-hartmann"))
    @test shack_hartmann.accelerator_scope == FullAccelerator

    atmosphere = plan_validation(("src/atmosphere/multilayer.jl",))
    @test atmosphere.selectors == Set(("ci-foundations", "ci-sensors-control"))
    @test atmosphere.accelerator_scope == FullAccelerator

    registered = plan_validation(("test/testsets/detector_cmos.jl",))
    @test registered.selectors == Set(("detector-cmos",))
    @test registered.accelerator_scope == DetectorAccelerator

    shared_fixture = plan_validation(("test/detector_test_fixtures.jl",))
    @test issubset(Set(DETECTOR_TEST_SUITE_NAMES), shared_fixture.selectors)

    cuda_extension = plan_validation(("ext/AdaptiveOpticsSimCUDAExt.jl",))
    @test cuda_extension.run_cuda
    @test !cuda_extension.run_amdgpu
    @test cuda_extension.accelerator_scope == FullAccelerator

    apple_extension =
        plan_validation(("ext/AdaptiveOpticsSimAppleAccelerateExt.jl",))
    @test apple_extension.manual_gates == Set(("appleaccelerate",))

    workflow = plan_validation((".github/workflows/cpu-validation.yml",))
    @test workflow.selectors == Set(("quality",))
    @test isempty(workflow.manual_gates)
    @test !workflow.run_cuda

    planner = plan_validation(("test/ci/impact_planner.jl",))
    @test planner.selectors == Set(("quality",))
    @test !planner.run_cuda
    @test !planner.run_amdgpu

    project = plan_validation(("Project.toml",))
    @test Set(CLOSURE_SELECTORS) == project.selectors
    @test project.run_aarch64
    @test project.run_grouped_cpu
    @test project.run_scheduler
    @test project.accelerator_scope == FullAccelerator
    @test project.manual_gates == Set((
        "macos-selector", "windows-selector", "appleaccelerate",
        "metal", "coverage"))

    unknown = plan_validation(("new_runtime/component.jl",))
    @test unknown.selectors == project.selectors
    @test unknown.manual_gates == project.manual_gates
    @test_throws ArgumentError plan_validation(("../outside.jl",))

    cmos_commands = Dict(local_commands(cmos))
    @test occursin("detector-cmos", cmos_commands["local CPU"])
    @test occursin("runtests_amdgpu_detectors.jl",
        cmos_commands["local AMDGPU"])
    @test occursin("ssh wsl", cmos_commands["WSL CUDA"])
    @test !haskey(cmos_commands, "Raspberry Pi aarch64")

    project_commands = Dict(local_commands(project))
    @test occursin("ssh raspberrypi",
        project_commands["Raspberry Pi aarch64"])
    @test occursin("runtests_cuda.jl", project_commands["WSL CUDA"])
    @test length(manual_commands(project)) == 5
    @test occursin("-f gate=macos-selector",
        Dict(manual_commands(project))["macos-selector"])
    @test occursin("--ref \"\$(git branch --show-current)\"",
        Dict(manual_commands(project))["macos-selector"])

    @test_throws ErrorException parse_cli(String[])
    @test_throws ErrorException parse_cli(["--all", "--path", "src/core/config.jl"])
    @test parse_cli(["--all"]).select_all
    @test parse_cli(["--path", "src/core/config.jl"]).paths ==
        ["src/core/config.jl"]

    mktempdir() do repository
        run(`git -C $repository init -q`)
        run(`git -C $repository config user.email validation@example.invalid`)
        run(`git -C $repository config user.name "Validation Test"`)
        mkpath(joinpath(repository, "src"))
        write(joinpath(repository, "src", "old.jl"), "old\n")
        run(`git -C $repository add src/old.jl`)
        run(`git -C $repository commit -qm baseline`)
        base_sha = readchomp(`git -C $repository rev-parse HEAD`)
        run(`git -C $repository mv src/old.jl src/new.jl`)
        run(`git -C $repository commit -qm rename`)
        paths = cd(repository) do
            changed_paths(base_sha)
        end
        @test Set(paths) == Set(("src/old.jl", "src/new.jl"))
    end

    @test isnothing(validate_impact_policy(joinpath(@__DIR__, "..", "..")))
end

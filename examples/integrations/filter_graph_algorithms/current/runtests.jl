using AdaptiveOpticsCalibration
using FilterGraphAlgorithms
using JuliaFilterGraph
using LinearAlgebra
using Test

const AOS_FGA_SOURCE_DIRECTION_STYLE = Val(:current)
include("../source_direction_keywords.jl")
include("../s1_lockstep.jl")
using .AOSFGALockstep

const EXPECTED_RESIDUAL_NORMS = Float32[
    0.15002562, 0.09018742, 0.05413042, 0.032470427, 0.019473683,
    0.011678018, 0.007003085, 0.0041994397, 0.002518242,
]
const EXPECTED_DEMANDED_COMMANDS = Float32[
    -1.1987396e-8, -1.9193582e-8, -2.3518728e-8, -2.6113192e-8,
    -2.7669184e-8, -2.8602285e-8, -2.9161848e-8, -2.9497393e-8,
    -2.9698606e-8,
]

@testset "registered AOC/FGA/JFG current S1 composition" begin
    @test pkgversion(AdaptiveOpticsCalibration) == v"0.17.0"
    @test pkgversion(FilterGraphAlgorithms) == v"0.5.5"
    @test pkgversion(JuliaFilterGraph) == v"0.2.3"

    prepared = prepare_s1_lockstep()
    @test @inferred(step_lockstep!(prepared)) == UInt64(1)
    @test @allocated(step_lockstep!(prepared)) == 0

    reset_s1_lockstep!(prepared)
    for index in eachindex(EXPECTED_RESIDUAL_NORMS)
        @test step_lockstep!(prepared) == UInt64(index)
        @test norm(prepared.outputs.slopes) ≈
            EXPECTED_RESIDUAL_NORMS[index] rtol=2f-5 atol=0f0
        @test prepared.outputs.demanded[1] ≈
            EXPECTED_DEMANDED_COMMANDS[index] rtol=1f-6
        @test prepared.adopted_command[1] == prepared.outputs.demanded[1]
        @test prepared.applied_command[1] ≈
            (isone(index) ? 0f0 : EXPECTED_DEMANDED_COMMANDS[index - 1]) rtol=1f-6
    end
end

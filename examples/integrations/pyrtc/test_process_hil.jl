using Test

include(joinpath(@__DIR__, "pyrtc_process_hil.jl"))
using .PyRTCProcessHIL

@testset "AOS and pyRTC close loops across native shared memory" begin
    for wavefront_sensor in (:shack_hartmann, :pyramid)
        result = PyRTCProcessHIL.run_validation(; wavefront_sensor)
        @test result.wavefront_sensor == wavefront_sensor
        @test result.interaction_rank == 25
        @test isfinite(result.interaction_condition)
        @test result.initial_residual > 0
        @test result.final_residual >= 0
        @test result.convergence_ratio < 1.0f-3
        @test result.command_error <= 5.0f-11
    end
end

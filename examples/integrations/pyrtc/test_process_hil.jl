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
        @test 0 < result.mean_open_loop_on_axis_strehl < 1
        @test 0.5 < result.mean_closed_loop_on_axis_strehl <= 1
        @test result.improvement > 10
    end
end

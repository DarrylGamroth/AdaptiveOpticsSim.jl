using Test

include(joinpath(@__DIR__, "pyrtc_process_hil.jl"))
using .PyRTCProcessHIL

@testset "AOS and pyRTC close loops across native shared memory" begin
    pyrtc_root = get(ENV, "PYRTC_ROOT", "")
    isdir(joinpath(pyrtc_root, "pyRTC")) || error(
        "set PYRTC_ROOT to a pyRTC checkout before running process HIL tests",
    )

    for wavefront_sensor in (:shack_hartmann, :pyramid)
        result = PyRTCProcessHIL.run_validation(
            pyrtc_root;
            wavefront_sensor,
        )
        @test result.wavefront_sensor == wavefront_sensor
        @test result.interaction_rank == 25
        @test isfinite(result.interaction_condition)
        @test result.initial_residual > 0
        @test result.final_residual >= 0
        @test result.convergence_ratio < 1.0f-3
        @test result.command_error <= 5.0f-11
    end
end

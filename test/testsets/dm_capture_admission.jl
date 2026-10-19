struct UnqualifiedDMCaptureActuators <: AdaptiveOpticsSim.Optics.AbstractDMActuatorModel end

@testset "Dense DM capture admission" begin
    graphs = AdaptiveOpticsSim.AlgorithmGraphs
    optics = AdaptiveOpticsSim.Optics
    backends = AdaptiveOpticsSim.Backends
    # This exercises capability dispatch only; actual device capture is
    # qualified separately on CUDA/HIP with measured influence matrices.
    style = backends.AcceleratorStyle(KernelAbstractions.CPU())
    modes = zeros(Float32, 4, 3)
    measured = optics.MeasuredInfluenceFunctions{Float32,typeof(modes),Nothing}(modes, nothing)
    dense = optics.DenseInfluenceMatrix{Float32,typeof(modes)}(modes)
    stages = (optics.LinearStaticActuators(), optics.ClippedActuators(-0.5f0, 0.5f0),
        optics.ActuatorHealthMap(Float32[1, 0, 1]))
    response = optics.CompositeDMActuatorModel(stages...)
    for model in (measured, dense), stage in (stages..., response)
        @test graphs._dm_surface_capture_capability(model, modes, stage, style) isa
            graphs.GraphNodeCaptureSafe
        @test graphs._dm_surface_capture_capability(model, modes, stage,
            backends.ScalarCPUStyle()) isa graphs.GraphNodeCaptureUnsupported
    end
    @test !graphs._static_dm_actuators_capture_safe(UnqualifiedDMCaptureActuators())
    @test !graphs._static_dm_actuators_capture_safe(
        optics.CompositeDMActuatorModel(stages[1], UnqualifiedDMCaptureActuators()))
end

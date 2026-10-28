# Optional Proper target. Keep the core-only construction check before loading Proper.
using AdaptiveOpticsSim.AlgorithmGraphs
using Test
include("../external_pupil_contract.jl")
@testset "Proper is optional" begin
    @test_throws AlgorithmGraphError proper_optical_train_node_factory()(:missing,
        (resolution=8, diameter_m=0.01, wavelength_um=1.55, steps=()), (;))
end
include("../test_proper_optical_train.jl")
include("runtime_opd.jl")
include("external_pupil.jl")

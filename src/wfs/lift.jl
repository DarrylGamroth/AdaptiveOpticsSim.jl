#
# LiFT physical forward model
#
# AOS owns focal propagation and observation formation. The iterative inverse
# estimator belongs to AdaptiveOpticsCalibration.

include("lift/kernels.jl")
include("lift/contracts.jl")
include("lift/forward.jl")
include("lift/modal_basis.jl")
include("lift/aoc_adapter.jl")

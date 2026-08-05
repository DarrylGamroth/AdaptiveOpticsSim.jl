#
# LiFT phase retrieval
#
# LiFT fits modal coefficients by matching a separately prepared focal-plane
# forward model to a caller-owned observation. Detector acquisition remains a
# separate concern.

include("lift/kernels.jl")
include("lift/contracts.jl")
include("lift/estimation.jl")
include("lift/forward.jl")

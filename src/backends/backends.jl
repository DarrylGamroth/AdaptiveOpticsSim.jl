"""
    Backends

Canonical owner of array-backend selectors, compute-device identity, and
backend extension protocols. Scientific and optical kernels remain with the
domains that define their physical meaning.
"""
module Backends

using AbstractFFTs
import FFTW
using FixedSizeArrays: FixedSizeVector
using KernelAbstractions
using LinearAlgebra
using Random

import ..AdaptiveOpticsSim:
    AdaptiveOpticsSimError,
    InvalidConfiguration,
    normal01,
    poisson_sample,
    runtime_rng,
    splitmix64,
    uniform01

include("arrays.jl")
include("reductions.jl")
include("random.jl")
include("api.jl")

end # module Backends

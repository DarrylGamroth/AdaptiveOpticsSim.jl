"""
    Ensembles

Canonical owner of execution policies and coarse-grained orchestration for
independent offline simulations, sweeps, and calibration work.
"""
module Ensembles

using LinearAlgebra

import ..AdaptiveOpticsSim:
    InvalidConfiguration,
    UnsupportedAlgorithm

import ..Backends:
    set_fft_provider_threads!

include("execution_policies.jl")
include("simulation_ensemble.jl")
include("api.jl")

end # module Ensembles

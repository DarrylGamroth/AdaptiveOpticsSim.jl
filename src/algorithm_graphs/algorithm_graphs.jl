"""
    AlgorithmGraphs

Portable, deterministic composition of transport-neutral algorithm
declarations. Calculon owns the scientific algorithm interface; this module
owns only static graph topology, exact buffer bindings, delayed links, and
model-time execution policy.
"""
module AlgorithmGraphs

import ..AdaptiveOpticsSim: AdaptiveOpticsSimError

include("definitions.jl")
include("preparation.jl")
include("execution.jl")
include("api.jl")

end # module AlgorithmGraphs

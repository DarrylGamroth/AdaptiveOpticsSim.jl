"""
    AlgorithmGraphs

Portable, deterministic composition of transport-neutral algorithm
declarations. Calculon owns the scientific algorithm interface; this module
owns only static graph topology, exact buffer bindings, delayed links, and
model-time execution policy.
"""
module AlgorithmGraphs

using FixedSizeArrays: FixedSizeVector, FixedSizeVectorDefault

import ..AdaptiveOpticsSim: AdaptiveOpticsSimError
import ..Plant:
    PeriodicSchedule,
    PlantDuration,
    PlantTimestamp,
    schedule_period,
    schedule_timestamp

include("definitions.jl")
include("preparation.jl")
include("execution.jl")
include("model_time.jl")
include("api.jl")

end # module AlgorithmGraphs

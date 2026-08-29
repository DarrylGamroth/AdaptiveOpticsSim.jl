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
import ..Backends:
    AbstractComputeDevice,
    HostComputeDevice,
    _prepare_device_execution_context,
    _synchronize_prepared_device_execution_context!,
    _with_prepared_device_execution_context,
    allocate_device_array,
    array_backend_type,
    compute_device,
    compute_device_backend
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

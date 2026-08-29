"""
    AlgorithmGraphs

Portable, deterministic composition of AOS-native algorithm nodes. Scientific
modules retain their domain APIs; this module owns only the graph-node adapter
protocol, static topology, exact buffer bindings, delayed links, and model-time
execution policy.
"""
module AlgorithmGraphs

using FixedSizeArrays: FixedSizeVector, FixedSizeVectorDefault
using TOML

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
import ..Calibration:
    ModalOPDExpansionPlan,
    combine_basis!
import ..Control:
    DiscreteIntegratorPlan,
    DiscreteIntegratorState,
    DiscreteIntegratorWorkspace,
    reset_controller!,
    update!

include("definitions.jl")
include("preparation.jl")
include("native_nodes.jl")
include("graph_files.jl")
include("execution.jl")
include("model_time.jl")
include("api.jl")

end # module AlgorithmGraphs

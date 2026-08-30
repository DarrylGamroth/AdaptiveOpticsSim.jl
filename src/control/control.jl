"""
    Control

Canonical owner of slopes-to-command reconstructors, discrete controller
models, controller composition, and preallocated control-path execution.
"""
module Control

using LinearAlgebra

import ..AdaptiveOpticsSim:
    DimensionMismatchError,
    InvalidConfiguration

import ..Backends:
    AbstractArrayBackend,
    CPUBackend,
    _resolve_array_backend,
    compute_device,
    require_same_backend

import ..Calibration:
    BuildBackend,
    InteractionMatrix,
    InversePolicy,
    condition_number,
    default_modal_inverse_policy,
    default_runtime_calibration_build_backend,
    effective_rank,
    inverse_factorization,
    inverse_operator,
    inverse_policy,
    materialize_build,
    materialize_runtime_build_result,
    singular_values,
    truncation_count

include("controller.jl")
include("reconstructors.jl")
include("api.jl")

end # module Control

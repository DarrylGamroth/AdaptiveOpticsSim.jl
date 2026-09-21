"""
    Control

Canonical owner of slopes-to-command reconstructors, discrete controller
models, controller composition, and preallocated control-path execution.
"""
module Control

import AdaptiveOpticsCalibration
using LinearAlgebra

import ..AdaptiveOpticsSim:
    DimensionMismatchError,
    InvalidConfiguration

import ..Backends:
    AbstractArrayBackend,
    CPUBackend,
    _resolve_array_backend,
    compute_device,
    host_array,
    require_same_backend

import ..Calibration:
    BuildBackend,
    InteractionMatrix,
    _default_svd_inverse_method,
    _prepare_svd_reconstructor,
    calibration_method,
    condition_number,
    default_runtime_calibration_build_backend,
    effective_rank,
    materialize_build,
    materialize_runtime_build_result,
    singular_values,
    truncation_count

include("controller.jl")
include("closed_loop_correction.jl")
include("deformable_mirror_routing.jl")
include("reconstructors.jl")
include("api.jl")

end # module Control

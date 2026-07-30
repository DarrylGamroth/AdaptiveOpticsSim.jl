"""
    Tomography

Canonical owner of guide-star geometry, atmospheric tomography models,
tomographic reconstruction, fitting, and DM-command projection.
"""
module Tomography

using KernelAbstractions
using LinearAlgebra
using SparseArrays
using SpecialFunctions

import ..AdaptiveOpticsSim:
    DimensionMismatchError,
    InvalidConfiguration,
    UnsupportedAlgorithm,
    _scaled_kv56_cpu,
    _scaled_kv56_scalar

import ..Backends:
    AcceleratorStyle,
    ScalarCPUStyle,
    backend_matmul_transpose_right,
    backend_symmetric_product,
    execution_style,
    launch_kernel_async!

import ..Optics:
    DeformableMirror

import ..Detectors:
    Detector,
    readout_noise

import ..WavefrontSensors:
    n_valid_subapertures

import ..Calibration:
    BuildBackend,
    CPUBuildBackend,
    GPUArrayBuildBackend,
    InteractionMatrix,
    NativeBuildBackend,
    _backend_array,
    default_build_backend,
    materialize_build,
    prepare_build_matrix

include("parameters.jl")
include("geometry.jl")
include("fitting.jl")
include("reconstructors.jl")
include("api.jl")

end # module Tomography

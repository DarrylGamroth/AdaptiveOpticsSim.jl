"""
    Detectors

Canonical owner of detector response, conventional frame and area-sensor
models, prepared acquisition, shutter/readout timing, and detector readout
products.
"""
module Detectors

using KernelAbstractions
using LinearAlgebra
using Random
using Statistics

import ..AdaptiveOpticsSim:
    AbstractDetector,
    DimensionMismatchError,
    InvalidConfiguration,
    UnsupportedAlgorithm,
    bin2d!,
    config_value,
    runtime_rng

import ..Backends:
    AbstractArrayBackend,
    AbstractComputeDevice,
    AcceleratorStyle,
    CPUBackend,
    DirectReductionPlan,
    ExecutionStyle,
    HostMirrorReductionPlan,
    ScalarCPUStyle,
    _AbstractPreparedDeviceExecutionContext,
    _HostPreparedDeviceExecutionContext,
    _clamp_array!,
    _poisson_noise!,
    _prepare_device_execution_context,
    _prepared_device_execution_compute_device,
    _rand_uniform_backend!,
    _randn_backend!,
    _resolve_array_backend,
    _resolve_backend_selector,
    _synchronize_prepared_device_execution_context!,
    _throw_compute_device_error,
    _with_prepared_device_execution_context,
    _write_integer_output!,
    allocate_array,
    array_backend_selector,
    backend,
    backend_fill,
    backend_matmul,
    backend_matmul_transpose_right,
    backend_maximum_value,
    backend_rand,
    backend_randn,
    backend_sum_value,
    backend_zeros,
    begin_kernel_phase,
    clamp_array!,
    compute_device,
    execute_fft_plan!,
    execution_style,
    finish_kernel_phase!,
    host_array,
    launch_kernel!,
    launch_kernel_async!,
    masked_sum2d,
    packed_valid_pair_mean,
    plan_fft_backend!,
    plan_ifft_backend!,
    poisson_noise!,
    queue_kernel!,
    rand_uniform_backend!,
    randn_backend!,
    randn_backend_async!,
    reduction_execution_plan,
    reduction_host_view,
    reduction_parent_source,
    require_same_backend,
    synchronize_backend!,
    write_integer_output!

import ..Optics:
    AbstractCombinationPolicy,
    AbstractOpticalNormalization,
    AbstractOpticalPlaneKind,
    AbstractSource,
    AbstractSpatialMeasure,
    AbstractSpectralCoordinate,
    CellIntegratedMeasure,
    DetectorPlane,
    DimensionlessNormalization,
    FocalPlane,
    IncoherentIntensityAddition,
    IntensityMap,
    IntegratedSpectralChannel,
    MonochromaticChannel,
    OpticalPlaneMetadata,
    PhotonRateNormalization,
    SpatialDensityMeasure,
    SpectralSource,
    backend,
    compute_device,
    intensity_values,
    plane_metadata,
    spectral_bundle,
    validate_plane_storage,
    wavelength

include("detector.jl")
include("api.jl")

end # module Detectors

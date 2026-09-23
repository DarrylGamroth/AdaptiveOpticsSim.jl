module AdaptiveOpticsSimAMDGPUExt

import AdaptiveOpticsSim
import AdaptiveOpticsSim: Backends, Calibration, Tomography,
    WavefrontSensors
using AMDGPU
using AbstractFFTs
using KernelAbstractions
using LinearAlgebra
using Random

#
# AMDGPU backend extension
#
# This extension supplies backend-native dense linear algebra and FFT plumbing
# for the maintained ROCArray execution paths. The main mathematical surfaces
# implemented here are:
#
# - stable Hermitian right division used by tomography/calibration
# - normal-equation solves and SVD fallback for LiFT
#
# The key rule is that the algorithms match the core implementation, while the
# execution is specialized to rocBLAS / rocSOLVER / rocFFT where that improves
# performance or avoids host fallback.
#
Backends.gpu_backend_loaded(::Type{Backends.AMDGPUBackendTag}) = true
Backends.gpu_backend_array_type(::Type{Backends.AMDGPUBackendTag}) = AMDGPU.ROCArray
Backends.gpu_backend_name(::Type{Backends.AMDGPUBackendTag}) = :amdgpu
Backends.gpu_backend_name(::Type{<:AMDGPU.ROCArray}) = :amdgpu
Backends.array_backend_selector(::Type{<:AMDGPU.ROCArray}) = Backends.AMDGPUBackend()
Backends.disable_scalar_backend!(::Type{Backends.AMDGPUBackendTag}) = AMDGPU.allowscalar(false)
Backends.backend_rand(::Type{Backends.AMDGPUBackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = AMDGPU.rand(T, dims...)
Backends.backend_randn(::Type{Backends.AMDGPUBackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = AMDGPU.randn(T, dims...)
Backends.backend_zeros(::Type{Backends.AMDGPUBackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = AMDGPU.zeros(T, dims...)
Backends.backend_fill(::Type{Backends.AMDGPUBackendTag}, value, dims::Vararg{Int}) = AMDGPU.fill(value, dims...)
Backends.compute_device_identifier(array::AMDGPU.ROCArray) =
    AMDGPU.device_id(AMDGPU.device(array))

function Backends._prepare_graph_rng(
    device::Backends.AcceleratorComputeDevice{Backends.AMDGPUBackend},
    seed::UInt64,
)
    return Backends._prepare_counter_rng(device, seed)
end

function Backends.compute_device_availability(
    device::Backends.AcceleratorComputeDevice{Backends.AMDGPUBackend,I},
) where {I<:Integer}
    identifier = try
        Int(Backends.compute_device_identifier(device))
    catch
        return Backends.ComputeDeviceUnavailable(:invalid_device_identifier)
    end
    identifier >= 1 || return Backends.ComputeDeviceUnavailable(
        :invalid_device_identifier)
    AMDGPU.functional() || return Backends.ComputeDeviceUnavailable(
        :backend_runtime_unavailable)
    try
        AMDGPU.HIPDevice(identifier)
    catch
        return Backends.ComputeDeviceUnavailable(:device_unavailable)
    end
    return Backends.ComputeDeviceAvailable()
end

Backends.compute_device_availability(
    ::Backends.AcceleratorComputeDevice{Backends.AMDGPUBackend},
) = Backends.ComputeDeviceUnavailable(:invalid_device_identifier)

@noinline function _require_amdgpu_device(
    device::Backends.AcceleratorComputeDevice{Backends.AMDGPUBackend,I},
) where {I<:Integer}
    availability = Backends.compute_device_availability(device)
    Backends.compute_device_is_available(availability) ||
        Backends._throw_compute_device_error(
            :select,
            Backends.compute_device_unavailable_reason(availability),
            device,
            "AMDGPU cannot address the requested device identifier",
        )
    return AMDGPU.HIPDevice(Int(Backends.compute_device_identifier(device)))
end

@noinline function _require_amdgpu_device(
    device::Backends.AcceleratorComputeDevice{Backends.AMDGPUBackend},
)
    Backends._throw_compute_device_error(
        :select,
        :invalid_device_identifier,
        device,
        "AMDGPU device identifiers must be positive integer identifiers",
    )
end

function Backends._with_compute_device(
    f::F,
    device::Backends.AcceleratorComputeDevice{Backends.AMDGPUBackend},
) where {F}
    return AMDGPU.device!(f, _require_amdgpu_device(device))
end

struct AMDGPUPreparedDeviceExecutionContext{
    D<:Backends.AcceleratorComputeDevice{Backends.AMDGPUBackend},
} <: Backends._AbstractPreparedDeviceExecutionContext
    device::AMDGPU.HIPDevice
    stream::AMDGPU.HIPStream
    compute_device::D
end

struct AMDGPUPreparedDeviceGraph{Graph,Executable}
    graph::Graph
    executable::Executable
end

function Backends._prepare_device_execution_context(
    storage::AMDGPU.ROCArray,
)
    device = AMDGPU.device(storage)
    stream = AMDGPU.device!(device) do
        AMDGPU.HIPStream()
    end
    return AMDGPUPreparedDeviceExecutionContext(
        device,
        stream,
        Backends.compute_device(storage),
    )
end

function Backends._prepare_device_execution_context(
    device::Backends.AcceleratorComputeDevice{Backends.AMDGPUBackend},
)
    runtime_device = _require_amdgpu_device(device)
    stream = AMDGPU.device!(runtime_device) do
        AMDGPU.HIPStream()
    end
    return AMDGPUPreparedDeviceExecutionContext(
        runtime_device,
        stream,
        device,
    )
end

@inline Backends._prepared_device_execution_compute_device(
    context::AMDGPUPreparedDeviceExecutionContext,
) = context.compute_device

@inline function Backends._with_prepared_device_execution_context(
    f::F,
    context::AMDGPUPreparedDeviceExecutionContext,
) where {F}
    old_device = AMDGPU.device()
    switched_device = old_device != context.device
    switched_device && AMDGPU.device!(context.device)
    old_stream = AMDGPU.stream()
    AMDGPU.stream!(context.stream)
    try
        return f()
    finally
        AMDGPU.stream!(old_stream)
        switched_device && AMDGPU.device!(old_device)
    end
end

@inline function Backends._synchronize_prepared_device_execution_context!(
    context::AMDGPUPreparedDeviceExecutionContext,
)
    AMDGPU.synchronize(context.stream)
    return nothing
end

# AMDGPU's blocking wait avoids event/task allocation but is incompatible with
# active HostCalls. AlgorithmGraphs selects this method only for native captured
# execution, whose admission contract excludes host callbacks.
@inline function Backends._synchronize_prepared_device_execution_context_blocking!(
    context::AMDGPUPreparedDeviceExecutionContext,
)
    AMDGPU.synchronize(context.stream; blocking=true)
    return nothing
end

@inline function Backends._prepare_device_execution_event(
    context::AMDGPUPreparedDeviceExecutionContext,
)
    return AMDGPU.HIP.HIPEvent(
        context.stream;
        do_record=false,
        timing=false,
    )
end

@inline function Backends._record_prepared_device_execution_event!(
    event::AMDGPU.HIP.HIPEvent,
    ::AMDGPUPreparedDeviceExecutionContext,
)
    AMDGPU.HIP.record(event)
    return nothing
end

@inline function Backends._wait_prepared_device_execution_event!(
    event::AMDGPU.HIP.HIPEvent,
    context::AMDGPUPreparedDeviceExecutionContext,
)
    AMDGPU.HIP.hipStreamWaitEvent(
        context.stream,
        event,
        AMDGPU.HIP.hipEventWaitDefault,
    )
    return nothing
end

function Backends._capture_prepared_device_graph(
    f::F,
    context::AMDGPUPreparedDeviceExecutionContext,
) where {F}
    graph = try
        AMDGPU.HIP.capture(f)
    catch
        # A rejected capture leaves a sticky HIP error on some ROCm releases.
        # Consume it so cold preparation failure does not poison later cleanup
        # or replace the original capture exception in the caller's boundary.
        AMDGPU.HIP.clear_last_error()
        rethrow()
    end
    executable = AMDGPU.HIP.instantiate(graph)
    return AMDGPUPreparedDeviceGraph(graph, executable)
end

@inline function Backends._launch_prepared_device_graph!(
    captured::AMDGPUPreparedDeviceGraph,
    context::AMDGPUPreparedDeviceExecutionContext,
)
    AMDGPU.HIP.launch(captured.executable, context.stream)
    return nothing
end

function Backends.execute_fft_plan!(buffer::AMDGPU.ROCArray, plan::AMDGPU.rocFFT.ROCFFTPlan)
    plan * buffer
    AMDGPU.synchronize()
    return buffer
end
function Backends.execute_fft_plan!(buffer::AMDGPU.ROCArray, plan::AbstractFFTs.ScaledPlan)
    plan * buffer
    AMDGPU.synchronize()
    return buffer
end
function Backends.enqueue_fft_plan!(
    buffer::AMDGPU.ROCArray,
    plan::AMDGPU.rocFFT.ROCFFTPlan,
)
    plan * buffer
    return buffer
end
function Backends.enqueue_fft_plan!(
    buffer::AMDGPU.ROCArray,
    plan::AbstractFFTs.ScaledPlan,
)
    plan * buffer
    return buffer
end
Calibration.default_build_backend(::AMDGPU.ROCArray) =
    Calibration.GPUArrayBuildBackend(Backends.AMDGPUBackendTag)
Calibration.prepare_build_matrix(
    ::Calibration.GPUArrayBuildBackend{Backends.AMDGPUBackendTag},
    A::AbstractMatrix,
) = Matrix(A)
WavefrontSensors.grouped_accumulation_strategy(
    ::Type{<:Backends.AcceleratorStyle{<:AMDGPU.ROCBackend}},
    ::Type{<:WavefrontSensors.PyramidWFS},
) = WavefrontSensors.GroupedStaged2DStrategy()
WavefrontSensors.grouped_accumulation_strategy(
    ::Type{<:Backends.AcceleratorStyle{<:AMDGPU.ROCBackend}},
    ::Type{<:WavefrontSensors.BiOEdgeWFS},
) = WavefrontSensors.GroupedStaged2DStrategy()
AdaptiveOpticsSim.Detectors.detector_execution_strategy(
    ::Type{<:Backends.AcceleratorStyle{<:AMDGPU.ROCBackend}},
    ::Type{<:AdaptiveOpticsSim.Detectors.Detector},
) = AdaptiveOpticsSim.Detectors.DetectorHostMirrorStrategy()
AdaptiveOpticsSim.Detectors._detector_value_strategy(
    strategy::AdaptiveOpticsSim.Detectors.DetectorHostMirrorStrategy,
    ::Backends.AcceleratorStyle{<:AMDGPU.ROCBackend},
) = strategy
AdaptiveOpticsSim.Detectors.can_apply_device_readout_correction(
    ::Backends.AcceleratorStyle{<:AMDGPU.ROCBackend},
    ::AdaptiveOpticsSim.Detectors.FrameReadoutCorrectionModel,
) = false
AdaptiveOpticsSim.Detectors.counting_output_execution_strategy(
    ::Type{<:Backends.AcceleratorStyle{<:AMDGPU.ROCBackend}},
    ::Type{<:AdaptiveOpticsSim.Detectors.AbstractCountingDetector},
    ::Type{<:AMDGPU.ROCArray{T,2}},
) where {T<:Integer} = AdaptiveOpticsSim.Detectors.DetectorHostMirrorStrategy()
Backends.reduction_execution_strategy(
    ::Backends.AcceleratorStyle{<:AMDGPU.ROCBackend},
    ::AMDGPU.ROCArray,
) = Backends.HostMirrorReductionStrategy()
Backends.randn_backend_async!(::Backends.AcceleratorStyle, rng::AbstractRNG, out::AMDGPU.ROCArray) = (Random.randn!(rng, out); out)
Backends._randn_backend!(::Backends.AcceleratorStyle, rng::AbstractRNG, out::AMDGPU.ROCArray) = (Random.randn!(rng, out); out)
function Backends.randn_backend_async!(
    style::Backends.AcceleratorStyle{<:AMDGPU.ROCBackend},
    rng::Backends._PreparedCounterRNG,
    out::AMDGPU.ROCArray{T},
) where {T<:AbstractFloat}
    return Backends._counter_randn_backend_async!(style, rng, out)
end
function Backends._randn_backend!(
    style::Backends.AcceleratorStyle{<:AMDGPU.ROCBackend},
    rng::Backends._PreparedCounterRNG,
    out::AMDGPU.ROCArray{T},
) where {T<:AbstractFloat}
    Backends._counter_randn_backend_async!(style, rng, out)
    AMDGPU.synchronize()
    return out
end
@inline function _prepared_randn_frame_noise!(
    rng::Backends._PreparedCounterRNG,
    out::AMDGPU.ROCArray{T,N},
) where {T<:AbstractFloat,N}
    Backends.randn_backend_async!(Backends.execution_style(out), rng, out)
    return out
end
function AdaptiveOpticsSim.Detectors._randn_frame_noise!(
    ::AdaptiveOpticsSim.Detectors.DetectorHostMirrorStrategy,
    det::AdaptiveOpticsSim.Detectors.Detector,
    rng::Backends._PreparedCounterRNG,
    out::AMDGPU.ROCArray{T,2},
) where {T<:AbstractFloat}
    return _prepared_randn_frame_noise!(rng, out)
end
function AdaptiveOpticsSim.Detectors._randn_frame_noise!(
    ::AdaptiveOpticsSim.Detectors.DetectorHostMirrorStrategy,
    det::AdaptiveOpticsSim.Detectors.Detector,
    rng::Backends._PreparedCounterRNG,
    out::AMDGPU.ROCArray{T,3},
) where {T<:AbstractFloat}
    return _prepared_randn_frame_noise!(rng, out)
end
function AdaptiveOpticsSim.Detectors._randn_frame_noise!(
    ::AdaptiveOpticsSim.Detectors.DetectorHostMirrorStrategy,
    det::AdaptiveOpticsSim.Detectors.Detector,
    rng::AbstractRNG,
    out::AMDGPU.ROCArray{T,2},
) where {T<:AbstractFloat}
    Backends.randn_backend!(rng, out)
    return out
end
function AdaptiveOpticsSim.Detectors._randn_frame_noise!(
    ::AdaptiveOpticsSim.Detectors.DetectorHostMirrorStrategy,
    det::AdaptiveOpticsSim.Detectors.Detector,
    rng::AbstractRNG,
    cube::AMDGPU.ROCArray{T,3},
) where {T<:AbstractFloat}
    Backends.randn_backend!(rng, cube)
    return cube
end
function AdaptiveOpticsSim.Detectors._poisson_noise_frame!(
    ::AdaptiveOpticsSim.Detectors.DetectorHostMirrorStrategy,
    det::AdaptiveOpticsSim.Detectors.Detector,
    rng::AbstractRNG,
    img::AMDGPU.ROCArray{T,2},
) where {T<:AbstractFloat}
    host = AdaptiveOpticsSim.Detectors.detector_host_frame!(det, img)
    Backends._poisson_noise!(Backends.ScalarCPUStyle(), rng, host)
    copyto!(img, host)
    return img
end
@inline function _prepared_poisson_noise_frame!(
    rng::Backends._PreparedCounterRNG,
    img::AMDGPU.ROCArray{T,N},
) where {T<:AbstractFloat,N}
    Backends.poisson_noise_async!(Backends.execution_style(img), rng, img)
    return img
end
function AdaptiveOpticsSim.Detectors._poisson_noise_frame!(
    ::AdaptiveOpticsSim.Detectors.DetectorHostMirrorStrategy,
    det::AdaptiveOpticsSim.Detectors.Detector,
    rng::Backends._PreparedCounterRNG,
    img::AMDGPU.ROCArray{T,2},
) where {T<:AbstractFloat}
    return _prepared_poisson_noise_frame!(rng, img)
end
function AdaptiveOpticsSim.Detectors._poisson_noise_frame!(
    ::AdaptiveOpticsSim.Detectors.DetectorHostMirrorStrategy,
    det::AdaptiveOpticsSim.Detectors.Detector,
    rng::AbstractRNG,
    cube::AMDGPU.ROCArray{T,3},
) where {T<:AbstractFloat}
    host = AdaptiveOpticsSim.Detectors.detector_host_cube!(det, cube)
    Backends._poisson_noise!(Backends.ScalarCPUStyle(), rng, host)
    copyto!(cube, host)
    return cube
end
function AdaptiveOpticsSim.Detectors._poisson_noise_frame!(
    ::AdaptiveOpticsSim.Detectors.DetectorHostMirrorStrategy,
    det::AdaptiveOpticsSim.Detectors.Detector,
    rng::Backends._PreparedCounterRNG,
    img::AMDGPU.ROCArray{T,3},
) where {T<:AbstractFloat}
    return _prepared_poisson_noise_frame!(rng, img)
end
function AdaptiveOpticsSim.Atmospheres.randn_phase_noise!(rng::AbstractRNG, out::AMDGPU.ROCArray{T,2}, host::Matrix{T}) where {T<:AbstractFloat}
    if size(host) != size(out)
        host = Matrix{T}(undef, size(out)...)
    end
    randn!(rng, host)
    copyto!(out, host)
    return host
end
function AdaptiveOpticsSim.Atmospheres.randn_phase_noise!(
    rng::Backends._PreparedCounterRNG,
    out::AMDGPU.ROCArray{T,2},
    host::Matrix{T},
) where {T<:AbstractFloat}
    Backends.randn_backend_async!(Backends.execution_style(out), rng, out)
    return host
end
function AdaptiveOpticsSim.Atmospheres._fill_phase_psd!(
    ::Backends.AcceleratorStyle{<:AMDGPU.ROCBackend},
    psd::AMDGPU.ROCArray{T,2},
    freqs::AMDGPU.ROCArray{T,1},
    coeff::T,
    inv_L0_sq::T,
    exponent::T,
    inv_fm_sq::T,
    n::Int,
) where {T<:AbstractFloat}
    host_psd = Matrix{T}(undef, size(psd))
    AdaptiveOpticsSim.Atmospheres._fill_phase_psd!(Backends.ScalarCPUStyle(), host_psd,
        Array(freqs), coeff, inv_L0_sq, exponent, inv_fm_sq, n)
    copyto!(psd, host_psd)
    return psd
end

# AMDGPU 2.7/GPUCompiler currently fails IR validation for the variable-trip
# KernelAbstractions slope kernels on gfx1030. Keep these backend-specific
# fallbacks explicit so other accelerator backends retain the device kernels.
function WavefrontSensors._geometric_slopes!(
    ::Backends.AcceleratorStyle{<:AMDGPU.ROCBackend},
    slopes::AMDGPU.ROCArray{T,1},
    opd::AMDGPU.ROCArray{T,2},
    valid_mask::AMDGPU.ROCArray{Bool,2},
    sub::Int,
    n_sub::Int,
    offset::Int,
) where {T<:AbstractFloat}
    host_slopes = Vector{T}(undef, length(slopes))
    WavefrontSensors._geometric_slopes!(Backends.ScalarCPUStyle(),
        host_slopes, Array(opd), Array(valid_mask), sub, n_sub, offset)
    copyto!(slopes, host_slopes)
    return slopes
end

function WavefrontSensors._edge_geometric_slopes!(
    ::Backends.AcceleratorStyle{<:AMDGPU.ROCBackend},
    slopes::AMDGPU.ROCArray{T,1},
    opd::AMDGPU.ROCArray{T,2},
    valid_mask::AMDGPU.ROCArray{Bool,2},
    edge_mask::AMDGPU.ROCArray{Bool,2},
    sub::Int,
    n_sub::Int,
    offset::Int,
) where {T<:AbstractFloat}
    host_slopes = Vector{T}(undef, length(slopes))
    WavefrontSensors._edge_geometric_slopes!(Backends.ScalarCPUStyle(),
        host_slopes, Array(opd), Array(valid_mask), Array(edge_mask), sub,
        n_sub, offset)
    copyto!(slopes, host_slopes)
    return slopes
end

Backends.backend_matmul(A::AMDGPU.ROCArray{T,2}, B::AMDGPU.ROCArray{T,2}) where {T<:AbstractFloat} =
    AMDGPU.rocBLAS.gemm('N', 'N', A, B)
Backends.backend_matmul_transpose_right(A::AMDGPU.ROCArray{T,2}, B::AMDGPU.ROCArray{T,2}) where {T<:AbstractFloat} =
    AMDGPU.rocBLAS.gemm('N', 'T', A, B)

function dense_copy_to_roc(A::AbstractMatrix{T}) where {T<:AbstractFloat}
    out = AMDGPU.ROCArray{T}(undef, size(A)...)
    copyto!(out, A)
    return out
end

function dense_host_matrix(A::SubArray{T,2,<:AMDGPU.ROCArray}) where {T<:AbstractFloat}
    host_parent = Array(parent(A))
    return Matrix(@view host_parent[parentindices(A)...])
end

dense_copy_to_roc(A::SubArray{T,2,<:AMDGPU.ROCArray}) where {T<:AbstractFloat} =
    AMDGPU.ROCArray(dense_host_matrix(A))

copy_dense_to_roc!(dest::AMDGPU.ROCArray, src::AMDGPU.ROCArray) = copyto!(dest, src)
copy_dense_to_roc!(dest::AMDGPU.ROCArray, src::SubArray{T,2,<:AMDGPU.ROCArray}) where {T<:AbstractFloat} =
    copyto!(dest, dense_host_matrix(src))
copy_dense_to_roc!(dest::AMDGPU.ROCArray, src::AbstractMatrix) = copyto!(dest, src)

function roc_svd(A::AMDGPU.ROCArray{T,2}) where {T<:AbstractFloat}
    F = copy(A)
    U, S, Vt = AMDGPU.rocSOLVER.gesvd!('S', 'S', F)
    return (; U, S, Vt, s_host=Calibration.singular_values_host(S))
end

end

# Array-backend, compute-device, FFT, and kernel-execution implementation.
abstract type ExecutionStyle end

abstract type AbstractArrayBackend end
abstract type GPUBackendTag end
abstract type GPUPrecisionPolicy end

struct CPUBackend <: AbstractArrayBackend end
struct CUDABackend <: AbstractArrayBackend end
struct MetalBackend <: AbstractArrayBackend end
struct AMDGPUBackend <: AbstractArrayBackend end

struct CUDABackendTag <: GPUBackendTag end
struct MetalBackendTag <: GPUBackendTag end
struct AMDGPUBackendTag <: GPUBackendTag end

struct UnifiedGPUPrecision{T<:AbstractFloat} <: GPUPrecisionPolicy end
struct SplitGPUPrecision{RT<:AbstractFloat,BT<:AbstractFloat} <: GPUPrecisionPolicy end

UnifiedGPUPrecision(::Type{T}) where {T<:AbstractFloat} = UnifiedGPUPrecision{T}()
SplitGPUPrecision(::Type{RT}, ::Type{BT}) where {RT<:AbstractFloat,BT<:AbstractFloat} = SplitGPUPrecision{RT,BT}()

struct ScalarCPUStyle <: ExecutionStyle end

struct AcceleratorStyle{B<:KernelAbstractions.Backend} <: ExecutionStyle
    backend::B
end

struct KernelLaunchPhase{S<:ExecutionStyle}
    style::S
end

#
# Prepared single-owner execution contexts
#
# These contexts are intentionally internal. They let one prepared Plant
# device-batch owner retain an exact backend stream without making stream
# objects, task policy, or accelerator packages part of core's public API.
# Optional backend extensions provide accelerator-specific implementations.
#
abstract type _AbstractPreparedDeviceExecutionContext end

execution_style(A::AbstractArray) = execution_style(KernelAbstractions.get_backend(A))
execution_style(::KernelAbstractions.CPU) = ScalarCPUStyle()
execution_style(backend::KernelAbstractions.Backend) = AcceleratorStyle(backend)
@inline use_host_build_algebra(::ScalarCPUStyle) = true
@inline use_host_build_algebra(::AcceleratorStyle) = false
host_array(A::AbstractArray) = host_array(execution_style(A), A)
host_array(::ScalarCPUStyle, A::AbstractArray) = A
host_array(::AcceleratorStyle, A::AbstractArray) = Array(A)

gpu_backend_loaded(::Type{<:GPUBackendTag}) = false

gpu_backend_array_type(::Type{<:GPUBackendTag}) = nothing

gpu_backend_name(::Type) = nothing

array_backend_type(::CPUBackend) = Array

function _require_gpu_array_backend(::Type{B}, label::AbstractString) where {B<:GPUBackendTag}
    backend = gpu_backend_array_type(B)
    isnothing(backend) && throw(InvalidConfiguration("$(label) is unavailable because the corresponding GPU backend extension is not loaded"))
    return backend
end

array_backend_type(::CUDABackend) = _require_gpu_array_backend(CUDABackendTag, "CUDABackend()")
array_backend_type(::MetalBackend) = _require_gpu_array_backend(MetalBackendTag, "MetalBackend()")
array_backend_type(::AMDGPUBackend) = _require_gpu_array_backend(AMDGPUBackendTag, "AMDGPUBackend()")
resolve_array_backend(backend::AbstractArrayBackend) = array_backend_type(backend)

array_backend_selector(::Type{<:Array}) = CPUBackend()
@inline array_backend_selector(
    ::Type{<:SubArray{T,N,P}}) where {T,N,P<:AbstractArray} =
    array_backend_selector(P)
@inline array_backend_selector(
    ::Type{<:Base.ReshapedArray{T,N,P}}) where {T,N,P<:AbstractArray} =
    array_backend_selector(P)
@inline array_backend_selector(
    ::Type{<:Base.ReinterpretArray{T,N,S,P}}) where {
        T,N,S,P<:AbstractArray,
    } = array_backend_selector(P)
@inline array_backend_selector(
    ::Type{<:PermutedDimsArray{T,N,Perm,IPerm,P}}) where {
        T,N,Perm,IPerm,P<:AbstractArray,
    } = array_backend_selector(P)
@inline array_backend_selector(
    ::Type{<:Transpose{T,P}}) where {T,P<:AbstractArray} =
    array_backend_selector(P)
@inline array_backend_selector(
    ::Type{<:Adjoint{T,P}}) where {T,P<:AbstractArray} =
    array_backend_selector(P)
array_backend_selector(::Type{<:AbstractArray}) =
    throw(InvalidConfiguration("no semantic backend selector is registered for the array storage type"))

# Internal normalization helpers used by constructor chaining after a backend has
# already been resolved to a semantic selector or array container type. Public
# API should prefer the semantic selector surface and call `resolve_array_backend`
# directly.
_resolve_backend_selector(backend::AbstractArrayBackend) = backend
_resolve_backend_selector(backend::Type{<:AbstractArray}) = array_backend_selector(backend)
_resolve_array_backend(backend::AbstractArrayBackend) = array_backend_type(backend)
_resolve_array_backend(backend::Type{<:AbstractArray}) = backend

backend(backend::AbstractArrayBackend) = backend
backend(A::AbstractArray) = array_backend_selector(typeof(A))
@inline backend(A::SubArray) = backend(parent(A))
@inline backend(A::Base.ReshapedArray) = backend(parent(A))
@inline backend(A::Base.ReinterpretArray) = backend(parent(A))
@inline backend(A::PermutedDimsArray) = backend(parent(A))
@inline backend(A::Transpose) = backend(parent(A))
@inline backend(A::Adjoint) = backend(parent(A))
backend_type(x) = typeof(backend(x))

"""
Concrete compute-device identity for array storage.

This is distinct from `AbstractArrayBackend`, which identifies a software
backend family rather than one host or accelerator device.
"""
abstract type AbstractComputeDevice end

"""Failure to inspect, select, or allocate on one exact compute device."""
struct ComputeDeviceError{D<:AbstractComputeDevice} <:
    AdaptiveOpticsSimError
    operation::Symbol
    reason::Symbol
    device::D
    msg::String
end

"""Cold runtime-availability result for one exact compute device."""
abstract type AbstractComputeDeviceAvailability end

"""The backend runtime can currently address the exact compute device."""
struct ComputeDeviceAvailable <: AbstractComputeDeviceAvailability end

"""The exact compute device is unavailable for the structured `reason`."""
struct ComputeDeviceUnavailable <: AbstractComputeDeviceAvailability
    reason::Symbol
end

"""The host compute device used by ordinary CPU array storage."""
struct HostComputeDevice <: AbstractComputeDevice end

struct _HostPreparedDeviceExecutionContext <:
    _AbstractPreparedDeviceExecutionContext
    device::HostComputeDevice
end

@inline _require_accelerator_backend(backend::AbstractArrayBackend) = backend
function _require_accelerator_backend(::CPUBackend)
    throw(InvalidConfiguration(
        "an accelerator compute device cannot use the CPU backend family"))
end

"""
One concrete accelerator, qualified by its semantic array-backend family.

`identifier` is supplied by the owning backend extension and is meaningful
within that backend runtime. It must distinguish simultaneously addressable
devices in one process.
"""
struct AcceleratorComputeDevice{
    B<:AbstractArrayBackend,I,
} <: AbstractComputeDevice
    backend::B
    identifier::I

    @inline function AcceleratorComputeDevice(
        backend::B,
        identifier::I,
    ) where {B<:AbstractArrayBackend,I}
        _require_accelerator_backend(backend)
        isnothing(identifier) && throw(InvalidConfiguration(
            "an accelerator compute device requires a concrete identifier"))
        isbitstype(I) || throw(InvalidConfiguration(
            "an accelerator compute-device identifier must be an isbits value"))
        return new{B,I}(backend, identifier)
    end
end

@inline compute_device(storage::AbstractArray) =
    compute_device(execution_style(storage), storage)
@inline compute_device(storage::SubArray) =
    compute_device(parent(storage))
@inline compute_device(storage::Base.ReshapedArray) =
    compute_device(parent(storage))
@inline compute_device(storage::Base.ReinterpretArray) =
    compute_device(parent(storage))
@inline compute_device(storage::PermutedDimsArray) =
    compute_device(parent(storage))
@inline compute_device(storage::Transpose) =
    compute_device(parent(storage))
@inline compute_device(storage::Adjoint) =
    compute_device(parent(storage))
@inline compute_device_identifier(::AbstractArray) = nothing
@inline compute_device_identifier(storage::SubArray) =
    compute_device_identifier(parent(storage))
@inline compute_device_identifier(storage::Base.ReshapedArray) =
    compute_device_identifier(parent(storage))
@inline compute_device_identifier(storage::Base.ReinterpretArray) =
    compute_device_identifier(parent(storage))
@inline compute_device_identifier(storage::PermutedDimsArray) =
    compute_device_identifier(parent(storage))
@inline compute_device_identifier(storage::Transpose) =
    compute_device_identifier(parent(storage))
@inline compute_device_identifier(storage::Adjoint) =
    compute_device_identifier(parent(storage))
@inline compute_device(::ScalarCPUStyle, ::AbstractArray) = HostComputeDevice()
@inline compute_device(::AcceleratorStyle, storage::AbstractArray) =
    AcceleratorComputeDevice(backend(storage),
        compute_device_identifier(storage))
@inline compute_device_backend(::HostComputeDevice) = CPUBackend()
@inline compute_device_backend(device::AcceleratorComputeDevice) =
    device.backend
@inline compute_device_identifier(::HostComputeDevice) = nothing
@inline compute_device_identifier(device::AcceleratorComputeDevice) =
    device.identifier

@inline compute_device_availability(::HostComputeDevice) =
    ComputeDeviceAvailable()
@inline compute_device_availability(::AcceleratorComputeDevice) =
    ComputeDeviceUnavailable(:exact_device_selection_unavailable)
function compute_device_is_available(::ComputeDeviceAvailable)
    return true
end

function compute_device_is_available(::ComputeDeviceUnavailable)
    return false
end
@inline compute_device_unavailable_reason(::ComputeDeviceAvailable) = nothing
@inline compute_device_unavailable_reason(
    availability::ComputeDeviceUnavailable,
) = availability.reason

@noinline function _throw_compute_device_error(
    operation::Symbol,
    reason::Symbol,
    device::D,
    detail::AbstractString,
) where {D<:AbstractComputeDevice}
    throw(ComputeDeviceError(
        operation,
        reason,
        device,
        "cannot $(operation) exact compute device $(device): $(detail)",
    ))
end

@inline _with_compute_device(f::F, ::HostComputeDevice) where {F} = f()

function _with_compute_device(
    ::F,
    device::AcceleratorComputeDevice,
) where {F}
    _throw_compute_device_error(
        :select,
        :exact_device_selection_unavailable,
        device,
        "the owning backend extension does not provide exact device selection",
    )
end

@inline function _prepare_device_execution_context(storage::AbstractArray)
    return _prepare_device_execution_context(execution_style(storage),
        storage)
end

@inline _prepare_device_execution_context(device::HostComputeDevice) =
    _HostPreparedDeviceExecutionContext(device)

function _prepare_device_execution_context(
    device::AcceleratorComputeDevice,
)
    _throw_compute_device_error(
        :prepare_context,
        :exact_device_selection_unavailable,
        device,
        "the owning backend extension does not provide an exact-device context",
    )
end

@inline function _prepare_device_execution_context(
    ::ScalarCPUStyle,
    storage::AbstractArray,
)
    return _HostPreparedDeviceExecutionContext(compute_device(storage))
end

function _prepare_device_execution_context(
    ::AcceleratorStyle,
    storage::AbstractArray,
)
    throw(InvalidConfiguration(
        "the accelerator storage type $(typeof(storage)) does not provide " *
        "a prepared device execution context"))
end

@inline _prepared_device_execution_compute_device(
    context::_HostPreparedDeviceExecutionContext,
) = context.device

@inline function _with_prepared_device_execution_context(
    f::F,
    ::_HostPreparedDeviceExecutionContext,
) where {F}
    return f()
end

@inline _synchronize_prepared_device_execution_context!(
    ::_HostPreparedDeviceExecutionContext,
) = nothing

#
# Prepared point-to-point array transfers
#
# A transfer implementation is deliberately opt-in for each supported source
# and destination storage pair.  Core provides no generic `copyto!` fallback:
# an accelerator extension (or a test backend) must retain any stream/event
# state it needs and must restore the caller's device context around launch,
# completion observation, and synchronous completion. Submission returns
# `_PreparedArrayTransferSubmitted`
# only after ownership has passed to the backend, or
# `_PreparedArrayTransferSubmissionFailed` only when no transfer can still
# access either slot. One nonblocking completion observation returns pending,
# completed, or failed. The synchronous completion seam blocks inside the
# backend and returns only completed or failed. Either failed result guarantees
# quiescence.
# Throwing or returning any other value leaves slot ownership uncertain and is
# a fail-stop condition for that prepared handoff.
#
@enum _PreparedArrayTransferSubmissionStatus::UInt8 begin
    _PreparedArrayTransferSubmitted = 0x01
    _PreparedArrayTransferSubmissionFailed = 0x02
end

@enum _PreparedArrayTransferCompletionStatus::UInt8 begin
    _PreparedArrayTransferPending = 0x01
    _PreparedArrayTransferCompleted = 0x02
    _PreparedArrayTransferCompletionFailed = 0x03
end

function _prepare_array_transfer(
    source::AbstractArray,
    destination::AbstractArray,
    source_context::_AbstractPreparedDeviceExecutionContext,
    destination_context::_AbstractPreparedDeviceExecutionContext,
)
    source_target = compute_device(source)
    destination_target = compute_device(destination)
    _throw_compute_device_error(
        :prepare_transfer,
        :unsupported_transfer_pair,
        destination_target,
        "no explicit array-transfer implementation is registered from " *
        "$(typeof(source)) on $source_target to $(typeof(destination)) " *
        "on $destination_target",
    )
end

function _submit_prepared_array_transfer!(
    transfer,
    destination::AbstractArray,
    source::AbstractArray,
)
    throw(InvalidConfiguration(
        "prepared array-transfer state $(typeof(transfer)) does not provide " *
        "an explicit submission implementation for $(typeof(source)) to " *
        "$(typeof(destination))"))
end

function _observe_prepared_array_transfer_completion!(transfer)
    throw(InvalidConfiguration(
        "prepared array-transfer state $(typeof(transfer)) does not provide " *
        "a nonblocking completion observation"))
end

function _complete_prepared_array_transfer!(transfer)
    throw(InvalidConfiguration(
        "prepared array-transfer state $(typeof(transfer)) does not provide " *
        "a synchronous completion implementation"))
end

@inline function same_backend(x, y)
    return backend_type(x) === backend_type(y)
end

@inline function same_backend(x, y, z, rest...)
    return same_backend(x, y) && same_backend(x, z) && all(arg -> same_backend(x, arg), rest)
end

@inline _first_composed_backend(::Tuple{}) = nothing

@inline function _first_composed_backend(xs::Tuple)
    value = first(xs)
    return isnothing(value) ? _first_composed_backend(Base.tail(xs)) :
        backend(value)
end

@inline _require_backend_matches!(ref, ::Tuple{}) = ref

@inline function _require_backend_matches!(ref, xs::Tuple)
    value = first(xs)
    if !isnothing(value)
        candidate = backend(value)
        typeof(candidate) === typeof(ref) || throw(InvalidConfiguration(
            "all composed objects in a simulation path must share the same backend; got $(typeof(ref)) and $(typeof(candidate))"))
    end
    return _require_backend_matches!(ref, Base.tail(xs))
end

@inline function require_same_backend(xs...)
    values = xs
    ref = _first_composed_backend(values)
    isnothing(ref) && return CPUBackend()
    return _require_backend_matches!(ref, values)
end

gpu_runtime_type(::UnifiedGPUPrecision{T}) where {T<:AbstractFloat} = T
gpu_build_type(::UnifiedGPUPrecision{T}) where {T<:AbstractFloat} = T
gpu_runtime_type(::SplitGPUPrecision{RT,BT}) where {RT<:AbstractFloat,BT<:AbstractFloat} = RT
gpu_build_type(::SplitGPUPrecision{RT,BT}) where {RT<:AbstractFloat,BT<:AbstractFloat} = BT

default_gpu_precision_policy(::Type{<:GPUBackendTag}) = SplitGPUPrecision(Float32, Float32)
high_accuracy_gpu_precision_policy(::Type{<:GPUBackendTag}) = SplitGPUPrecision(Float32, Float64)

disable_scalar_backend!(::Type{<:GPUBackendTag}) = nothing
backend_rand(::Type{<:GPUBackendTag}, ::Type, dims::Vararg{Int}) = error("GPU backend random generation is not available")
backend_randn(::Type{<:GPUBackendTag}, ::Type, dims::Vararg{Int}) = error("GPU backend normal random generation is not available")
backend_zeros(::Type{<:GPUBackendTag}, ::Type, dims::Vararg{Int}) = error("GPU backend zeros allocation is not available")
backend_fill(::Type{<:GPUBackendTag}, value, dims::Vararg{Int}) = error("GPU backend fill allocation is not available")

function available_gpu_backends()
    out = DataType[]
    gpu_backend_loaded(CUDABackendTag) && push!(out, CUDABackendTag)
    gpu_backend_loaded(MetalBackendTag) && push!(out, MetalBackendTag)
    gpu_backend_loaded(AMDGPUBackendTag) && push!(out, AMDGPUBackendTag)
    return Tuple(out)
end

allocate_array(backend, ::Type{T}, dims::Vararg{Int,N}) where {T,N} =
    _resolve_array_backend(backend){T}(undef, dims...)

"""
Allocate uninitialized array storage on one exact compute device.

This is a preparation-time backend extension seam. The owning extension must
select the requested device, restore the caller's previous device context, and
return storage whose reported `compute_device` exactly matches `device`.
"""
function allocate_device_array(
    device::AbstractComputeDevice,
    ::Type{T},
    dims::Vararg{Int,N},
) where {T,N}
    return _with_compute_device(device) do
        storage = allocate_array(compute_device_backend(device), T, dims...)
        actual = compute_device(storage)
        actual == device || _throw_compute_device_error(
            :allocate,
            :wrong_device,
            device,
            "backend allocation returned storage on $(actual)",
        )
        return storage
    end
end

# Keep the default CPU FFT provider explicit. Optional provider extensions may
# add more-specific methods for ordinary Arrays; accelerator arrays continue to
# use their backend's AbstractFFTs dispatch.
@inline function _plan_fftw_fft!(buffer::StridedArray{T,N}, dims) where {
    T<:Union{ComplexF32,ComplexF64},N}
    return invoke(plan_fft!, Tuple{StridedArray{T,N},Any}, buffer, dims)
end

@inline function _plan_fftw_bfft!(buffer::StridedArray{T,N}, dims) where {
    T<:Union{ComplexF32,ComplexF64},N}
    return invoke(plan_bfft!, Tuple{StridedArray{T,N},Any}, buffer, dims)
end

@inline _fft_region_length(buffer, dim::Integer) = size(buffer, dim)
@inline _fft_region_length(buffer, dims) = prod(size(buffer, dim) for dim in dims)

function plan_fft_backend!(buffer::Array{T,N}) where {
    T<:Union{ComplexF32,ComplexF64},N}
    return plan_fft_backend!(buffer, ntuple(identity, N))
end

function plan_fft_backend!(buffer::Array{T,N}, dims) where {
    T<:Union{ComplexF32,ComplexF64},N}
    return _plan_fftw_fft!(buffer, dims)
end

function plan_ifft_backend!(buffer::Array{T,N}) where {
    T<:Union{ComplexF32,ComplexF64},N}
    return plan_ifft_backend!(buffer, ntuple(identity, N))
end

function plan_ifft_backend!(buffer::Array{T,N}, dims) where {
    T<:Union{ComplexF32,ComplexF64},N}
    scale = one(real(T)) / _fft_region_length(buffer, dims)
    return scale * _plan_fftw_bfft!(buffer, dims)
end

plan_fft_backend!(buffer) = plan_fft!(buffer)
plan_fft_backend!(buffer, dims) = plan_fft!(buffer, dims)
plan_ifft_backend!(buffer) = plan_ifft!(buffer)
plan_ifft_backend!(buffer, dims) = plan_ifft!(buffer, dims)
set_fft_provider_threads!(n::Integer) = FFTW.set_num_threads(n)
execute_fft_plan!(buffer, plan) = (mul!(buffer, plan, buffer); buffer)
backend_matmul(A::AbstractMatrix, B::AbstractMatrix) = A * B
backend_matmul_transpose_right(A::AbstractMatrix, B::AbstractMatrix) = A * transpose(B)

function backend_symmetric_product(A::AbstractMatrix, B::AbstractMatrix)
    tmp = backend_matmul(A, B)
    return backend_matmul_transpose_right(tmp, A)
end

@inline synchronize_backend!(::ScalarCPUStyle) = nothing
@inline synchronize_backend!(style::AcceleratorStyle) = KernelAbstractions.synchronize(style.backend)

@inline begin_kernel_phase(style::ExecutionStyle) = KernelLaunchPhase(style)
@inline finish_kernel_phase!(phase::KernelLaunchPhase) = synchronize_backend!(phase.style)

@inline preferred_workgroupsize(::AcceleratorStyle, ndrange) = nothing
@inline preferred_workgroupsize(::AcceleratorStyle, ndrange::Integer) = min(ndrange, 256)
@inline preferred_workgroupsize(::AcceleratorStyle, ndrange::Tuple{Int,Int}) =
    (min(ndrange[1], 16), min(ndrange[2], 16))
@inline preferred_workgroupsize(::AcceleratorStyle, ndrange::Tuple{Int,Int,Int}) =
    (min(ndrange[1], 8), min(ndrange[2], 8), min(ndrange[3], 4))

@inline function launch_kernel_async!(style::AcceleratorStyle, kernel, args...; ndrange)
    workgroupsize = preferred_workgroupsize(style, ndrange)
    if isnothing(workgroupsize)
        kernel(style.backend)(args...; ndrange=ndrange)
    else
        kernel(style.backend)(args...; ndrange=ndrange, workgroupsize=workgroupsize)
    end
    return nothing
end

@inline function queue_kernel!(phase::KernelLaunchPhase{<:AcceleratorStyle}, kernel, args...; ndrange)
    launch_kernel_async!(phase.style, kernel, args...; ndrange=ndrange)
    return nothing
end

@inline function launch_kernel!(style::AcceleratorStyle, kernel, args...; ndrange)
    launch_kernel_async!(style, kernel, args...; ndrange=ndrange)
    synchronize_backend!(style)
    return nothing
end

@kernel function clamp_array_kernel!(array, lower, upper, n::Int)
    i = @index(Global, Linear)
    if i <= n
        @inbounds array[i] = clamp(array[i], lower, upper)
    end
end

@kernel function integer_output_kernel!(output, input, lower, upper, n::Int, m::Int)
    i, j = @index(Global, NTuple)
    if i <= n && j <= m
        value = round(@inbounds input[i, j])
        @inbounds output[i, j] = clamp(value, lower, upper)
    end
end

function clamp_array!(array::AbstractArray, lower, upper)
    return _clamp_array!(execution_style(array), array, lower, upper)
end

function _clamp_array!(::ScalarCPUStyle, array::AbstractArray, lower, upper)
    clamp!(array, lower, upper)
    return array
end

function _clamp_array!(style::AcceleratorStyle, array::AbstractArray, lower, upper)
    isempty(array) && return array
    launch_kernel!(style, clamp_array_kernel!, array, lower, upper, length(array);
        ndrange=length(array))
    return array
end

function write_integer_output!(output::AbstractMatrix, input::AbstractMatrix)
    return _write_integer_output!(execution_style(output), output, input)
end

function _write_integer_output!(::ScalarCPUStyle, output::AbstractMatrix,
    input::AbstractMatrix)
    output_type = eltype(output)
    lower = typemin(output_type)
    upper = typemax(output_type)
    @. output = output_type(clamp(round(input), lower, upper))
    return output
end

function _write_integer_output!(style::AcceleratorStyle, output::AbstractMatrix,
    input::AbstractMatrix)
    output_type = eltype(output)
    n, m = size(output)
    launch_kernel!(style, integer_output_kernel!, output, input, typemin(output_type),
        typemax(output_type), n, m; ndrange=(n, m))
    return output
end

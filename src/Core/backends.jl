import FFTW

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

execution_style(A::AbstractArray) = execution_style(KernelAbstractions.get_backend(A))
execution_style(::KernelAbstractions.CPU) = ScalarCPUStyle()
execution_style(backend::KernelAbstractions.Backend) = AcceleratorStyle(backend)

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

array_backend_selector(::Type{Array}) = CPUBackend()
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
backend_type(x) = typeof(backend(x))

@inline function same_backend(x, y)
    return backend_type(x) === backend_type(y)
end

@inline function same_backend(x, y, z, rest...)
    return same_backend(x, y) && same_backend(x, z) && all(arg -> same_backend(x, arg), rest)
end

function require_same_backend(xs...)
    filtered = Tuple(x for x in xs if !isnothing(x))
    isempty(filtered) && return CPUBackend()
    ref = backend(filtered[1])
    refT = typeof(ref)
    for x in filtered[2:end]
        bx = backend(x)
        typeof(bx) === refT || throw(InvalidConfiguration("all composed objects in a simulation path must share the same backend; got $(refT) and $(typeof(bx))"))
    end
    return ref
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

allocate_array(backend, ::Type{T}, dims::Vararg{Int,N}) where {T,N} = _resolve_array_backend(backend){T}(undef, dims...)

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

@inline function launch_kernel_async!(style::AcceleratorStyle, kernel, args...; ndrange)
    kernel(style.backend)(args...; ndrange=ndrange)
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

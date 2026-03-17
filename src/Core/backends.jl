abstract type ExecutionStyle end

abstract type GPUBackendTag end

struct CUDABackendTag <: GPUBackendTag end
struct MetalBackendTag <: GPUBackendTag end
struct AMDGPUBackendTag <: GPUBackendTag end

struct ScalarCPUStyle <: ExecutionStyle end

struct AcceleratorStyle{B<:KernelAbstractions.Backend} <: ExecutionStyle
    backend::B
end

execution_style(A::AbstractArray) = execution_style(KernelAbstractions.get_backend(A))
execution_style(::KernelAbstractions.CPU) = ScalarCPUStyle()
execution_style(backend::KernelAbstractions.Backend) = AcceleratorStyle(backend)

gpu_backend_loaded(::Type{<:GPUBackendTag}) = false

gpu_backend_array_type(::Type{<:GPUBackendTag}) = nothing

gpu_backend_name(::Type) = nothing

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

allocate_array(backend, ::Type{T}, dims::Vararg{Int,N}) where {T,N} = backend{T}(undef, dims...)

plan_fft_backend!(buffer) = plan_fft!(buffer)
plan_ifft_backend!(buffer) = plan_ifft!(buffer)

@inline synchronize_backend!(::ScalarCPUStyle) = nothing
@inline synchronize_backend!(style::AcceleratorStyle) = KernelAbstractions.synchronize(style.backend)

@inline function launch_kernel!(style::AcceleratorStyle, kernel, args...; ndrange)
    kernel(style.backend)(args...; ndrange=ndrange)
    synchronize_backend!(style)
    return nothing
end

abstract type ExecutionStyle end

struct ScalarCPUStyle <: ExecutionStyle end

struct AcceleratorStyle{B<:KernelAbstractions.Backend} <: ExecutionStyle
    backend::B
end

execution_style(A::AbstractArray) = execution_style(KernelAbstractions.get_backend(A))
execution_style(::KernelAbstractions.CPU) = ScalarCPUStyle()
execution_style(backend::KernelAbstractions.Backend) = AcceleratorStyle(backend)

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

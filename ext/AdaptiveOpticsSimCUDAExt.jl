module AdaptiveOpticsSimCUDAExt

using AdaptiveOpticsSim
using CUDA

AdaptiveOpticsSim.gpu_backend_loaded(::Type{AdaptiveOpticsSim.CUDABackendTag}) = true
AdaptiveOpticsSim.gpu_backend_array_type(::Type{AdaptiveOpticsSim.CUDABackendTag}) = CUDA.CuArray
AdaptiveOpticsSim.gpu_backend_name(::Type{AdaptiveOpticsSim.CUDABackendTag}) = :cuda
AdaptiveOpticsSim.gpu_backend_name(::Type{<:CUDA.CuArray}) = :cuda
AdaptiveOpticsSim.disable_scalar_backend!(::Type{AdaptiveOpticsSim.CUDABackendTag}) = CUDA.allowscalar(false)
AdaptiveOpticsSim.backend_rand(::Type{AdaptiveOpticsSim.CUDABackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = CUDA.rand(T, dims...)
AdaptiveOpticsSim.backend_randn(::Type{AdaptiveOpticsSim.CUDABackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = CUDA.randn(T, dims...)
AdaptiveOpticsSim.backend_zeros(::Type{AdaptiveOpticsSim.CUDABackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = CUDA.zeros(T, dims...)
AdaptiveOpticsSim.backend_fill(::Type{AdaptiveOpticsSim.CUDABackendTag}, value, dims::Vararg{Int}) = CUDA.fill(value, dims...)

end

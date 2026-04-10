module AdaptiveOpticsSimCUDAExt

using AdaptiveOpticsSim
using CUDA
using LinearAlgebra

AdaptiveOpticsSim.gpu_backend_loaded(::Type{AdaptiveOpticsSim.CUDABackendTag}) = true
AdaptiveOpticsSim.gpu_backend_array_type(::Type{AdaptiveOpticsSim.CUDABackendTag}) = CUDA.CuArray
AdaptiveOpticsSim.gpu_backend_name(::Type{AdaptiveOpticsSim.CUDABackendTag}) = :cuda
AdaptiveOpticsSim.gpu_backend_name(::Type{<:CUDA.CuArray}) = :cuda
AdaptiveOpticsSim.array_backend_selector(::Type{<:CUDA.CuArray}) = AdaptiveOpticsSim.CUDABackend()
AdaptiveOpticsSim.disable_scalar_backend!(::Type{AdaptiveOpticsSim.CUDABackendTag}) = CUDA.allowscalar(false)
AdaptiveOpticsSim.backend_rand(::Type{AdaptiveOpticsSim.CUDABackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = CUDA.rand(T, dims...)
AdaptiveOpticsSim.backend_randn(::Type{AdaptiveOpticsSim.CUDABackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = CUDA.randn(T, dims...)
AdaptiveOpticsSim.backend_zeros(::Type{AdaptiveOpticsSim.CUDABackendTag}, ::Type{T}, dims::Vararg{Int}) where {T} = CUDA.zeros(T, dims...)
AdaptiveOpticsSim.backend_fill(::Type{AdaptiveOpticsSim.CUDABackendTag}, value, dims::Vararg{Int}) = CUDA.fill(value, dims...)

function AdaptiveOpticsSim.masked_sum2d_accelerator(
    ::AdaptiveOpticsSim.AcceleratorStyle{<:CUDA.CUDAKernels.CUDABackend},
    values_parent::CUDA.CuArray{T,2},
    valid_mask::CUDA.CuArray{Bool,2},
    scalar_buffer::AbstractVector{T},
    scalar_host::AbstractVector{T},
    host_parent::AbstractMatrix{T},
    row_offset::Int,
    col_offset::Int,
    n_rows::Int,
    n_cols::Int,
) where {T<:AbstractFloat}
    rows = (row_offset + 1):(row_offset + n_rows)
    cols = (col_offset + 1):(col_offset + n_cols)
    values_view = @view values_parent[rows, cols]
    total = LinearAlgebra.dot(vec(values_view), vec(T.(valid_mask)))
    return total, host_parent
end

end

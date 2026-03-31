@kernel function masked_sum2d_kernel!(out, values, valid_mask, n_rows::Int, n_cols::Int)
    i, j = @index(Global, NTuple)
    if i <= n_rows && j <= n_cols && @inbounds(valid_mask[i, j])
        T = eltype(out)
        @atomic out[1] += T(@inbounds values[i, j])
    end
end

@inline reduction_parent_source(A::AbstractArray) = A
@inline reduction_parent_source(A::SubArray) = parent(A)

@inline reduction_host_view(host_parent::AbstractArray, A::AbstractArray) = host_parent
@inline reduction_host_view(host_parent::AbstractArray, A::SubArray) = @view(host_parent[parentindices(A)...])

@inline reduction_backend_name(A::AbstractArray) = _reduction_backend_name(A)
@inline _reduction_backend_name(A::AbstractArray) = something(gpu_backend_name(typeof(A)), _reduction_backend_name_parent(A))
@inline _reduction_backend_name_parent(A::AbstractArray) = nothing
@inline _reduction_backend_name_parent(A::SubArray) = gpu_backend_name(typeof(parent(A)))

@inline function ensure_reduction_host_parent(host_parent::AbstractMatrix{T}, src::AbstractMatrix{T}) where {T}
    return size(host_parent) == size(src) ? host_parent : Matrix{T}(undef, size(src)...)
end

@inline backend_maximum_value(A::AbstractArray{T}) where {T<:AbstractFloat} = backend_maximum_value(execution_style(A), A)
@inline backend_maximum_value(::ScalarCPUStyle, A::AbstractArray{T}) where {T<:AbstractFloat} = maximum(A)
@inline function backend_maximum_value(::AcceleratorStyle, A::AbstractArray{T}) where {T<:AbstractFloat}
    if reduction_backend_name(A) === :amdgpu
        return maximum(Array(A))
    end
    return maximum(A)
end

@inline function masked_sum2d(values::AbstractMatrix{T}, valid_mask::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    return masked_sum2d(ScalarCPUStyle(), values, valid_mask)
end

function masked_sum2d(::ScalarCPUStyle, values::AbstractMatrix{T}, valid_mask::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    summed = zero(T)
    @inbounds for j in axes(valid_mask, 2), i in axes(valid_mask, 1)
        if valid_mask[i, j]
            summed += values[i, j]
        end
    end
    return summed
end

function masked_sum2d(style::AcceleratorStyle, values::AbstractMatrix{T}, valid_mask::AbstractMatrix{Bool},
    valid_mask_host::AbstractMatrix{Bool}, scalar_buffer::AbstractVector{T}, scalar_host::AbstractVector{T},
    host_parent::AbstractMatrix{T}) where {T<:AbstractFloat}
    if reduction_backend_name(values) === :amdgpu
        src = reduction_parent_source(values)
        refreshed_parent = ensure_reduction_host_parent(host_parent, src)
        copyto!(refreshed_parent, src)
        host_values = reduction_host_view(refreshed_parent, values)
        return masked_sum2d(ScalarCPUStyle(), host_values, valid_mask_host), refreshed_parent
    end
    fill!(scalar_buffer, zero(T))
    n_rows, n_cols = size(valid_mask)
    launch_kernel!(style, masked_sum2d_kernel!, scalar_buffer, values, valid_mask, n_rows, n_cols; ndrange=(n_rows, n_cols))
    copyto!(scalar_host, scalar_buffer)
    return scalar_host[1], host_parent
end

function packed_valid_pair_mean(signal::AbstractVector{T}, valid_mask::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    return packed_valid_pair_mean(execution_style(signal), signal, valid_mask)
end

function packed_valid_pair_mean(::ScalarCPUStyle, signal::AbstractVector{T}, valid_mask::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    n_sub = size(valid_mask, 1)
    offset = n_sub * n_sub
    acc = zero(T)
    count = 0
    idx = 1
    @inbounds for i in 1:n_sub, j in 1:n_sub
        if valid_mask[i, j]
            acc += signal[idx]
            acc += signal[idx + offset]
            count += 2
        end
        idx += 1
    end
    return count == 0 ? one(T) : acc / count
end

function packed_valid_pair_mean(::AcceleratorStyle, signal::AbstractVector{T}, valid_mask::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    if reduction_backend_name(signal) === :amdgpu
        return packed_valid_pair_mean(ScalarCPUStyle(), Array(signal), Array(valid_mask))
    end
    count = sum(valid_mask)
    return count == 0 ? one(T) : sum(signal) / (T(2) * T(count))
end

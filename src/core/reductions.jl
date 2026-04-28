@kernel function masked_sum2d_kernel!(out, values_parent, valid_mask, row_offset::Int, col_offset::Int, n_rows::Int, n_cols::Int)
    i, j = @index(Global, NTuple)
    if i <= n_rows && j <= n_cols && @inbounds(valid_mask[i, j])
        T = eltype(out)
        @atomic out[1] += T(@inbounds values_parent[i + row_offset, j + col_offset])
    end
end

abstract type AbstractReductionExecutionPlan end

struct DirectReductionPlan <: AbstractReductionExecutionPlan end
struct HostMirrorReductionPlan <: AbstractReductionExecutionPlan end

@inline reduction_parent_source(A::AbstractArray) = A
@inline reduction_parent_source(A::SubArray) = parent(A)

@inline reduction_axis_offset(A::AbstractArray, dim::Int) = 0
@inline reduction_axis_offset(A::SubArray, dim::Int) = reduction_index_offset(parentindices(A)[dim])
@inline reduction_index_offset(r::AbstractUnitRange{<:Integer}) = first(r) - 1
@inline reduction_index_offset(::Base.Slice) = 0

@inline reduction_host_view(host_parent::AbstractArray, A::AbstractArray) = host_parent
@inline reduction_host_view(host_parent::AbstractArray, A::SubArray) = @view(host_parent[parentindices(A)...])

@inline reduction_execution_plan(A::AbstractArray) = reduction_execution_plan(execution_style(A), reduction_parent_source(A))
@inline reduction_execution_plan(::ScalarCPUStyle, ::AbstractArray) = DirectReductionPlan()
@inline reduction_execution_plan(::AcceleratorStyle, ::AbstractArray) = DirectReductionPlan()
@inline reduction_execution_plan(::AcceleratorStyle{<:KernelAbstractions.CPU}, ::Array) = HostMirrorReductionPlan()

@inline function ensure_reduction_host_parent(host_parent::AbstractMatrix{T}, src::AbstractMatrix{T}) where {T}
    return size(host_parent) == size(src) ? host_parent : Matrix{T}(undef, size(src)...)
end

@inline backend_maximum_value(A::AbstractArray{T}) where {T<:AbstractFloat} = backend_maximum_value(execution_style(A), A)
@inline backend_maximum_value(::ScalarCPUStyle, A::AbstractArray{T}) where {T<:AbstractFloat} = maximum(A)
@inline backend_maximum_value(style::AcceleratorStyle, A::AbstractArray{T}) where {T<:AbstractFloat} =
    _backend_maximum_value(reduction_execution_plan(style, reduction_parent_source(A)), A)
@inline _backend_maximum_value(::DirectReductionPlan, A::AbstractArray{T}) where {T<:AbstractFloat} = maximum(A)
@inline _backend_maximum_value(::HostMirrorReductionPlan, A::AbstractArray{T}) where {T<:AbstractFloat} = maximum(Array(A))

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
    plan = reduction_execution_plan(style, reduction_parent_source(values))
    return masked_sum2d(plan, style, values, valid_mask, valid_mask_host, scalar_buffer, scalar_host, host_parent)
end

function masked_sum2d(::HostMirrorReductionPlan, style::AcceleratorStyle, values::AbstractMatrix{T}, valid_mask::AbstractMatrix{Bool},
    valid_mask_host::AbstractMatrix{Bool}, scalar_buffer::AbstractVector{T}, scalar_host::AbstractVector{T},
    host_parent::AbstractMatrix{T}) where {T<:AbstractFloat}
    src = reduction_parent_source(values)
    refreshed_parent = ensure_reduction_host_parent(host_parent, src)
    copyto!(refreshed_parent, src)
    host_values = reduction_host_view(refreshed_parent, values)
    return masked_sum2d(ScalarCPUStyle(), host_values, valid_mask_host), refreshed_parent
end

function masked_sum2d(::DirectReductionPlan, style::AcceleratorStyle, values::AbstractMatrix{T}, valid_mask::AbstractMatrix{Bool},
    valid_mask_host::AbstractMatrix{Bool}, scalar_buffer::AbstractVector{T}, scalar_host::AbstractVector{T},
    host_parent::AbstractMatrix{T}) where {T<:AbstractFloat}
    n_rows, n_cols = size(valid_mask)
    values_parent = reduction_parent_source(values)
    row_offset = reduction_axis_offset(values, 1)
    col_offset = reduction_axis_offset(values, 2)
    return masked_sum2d_accelerator(
        style,
        values_parent,
        valid_mask,
        scalar_buffer,
        scalar_host,
        host_parent,
        row_offset,
        col_offset,
        n_rows,
        n_cols,
    )
end

function masked_sum2d_accelerator(style::AcceleratorStyle, values_parent::AbstractMatrix{T}, valid_mask::AbstractMatrix{Bool},
    scalar_buffer::AbstractVector{T}, scalar_host::AbstractVector{T}, host_parent::AbstractMatrix{T},
    row_offset::Int, col_offset::Int, n_rows::Int, n_cols::Int) where {T<:AbstractFloat}
    fill!(scalar_buffer, zero(T))
    launch_kernel!(
        style,
        masked_sum2d_kernel!,
        scalar_buffer,
        values_parent,
        valid_mask,
        row_offset,
        col_offset,
        n_rows,
        n_cols;
        ndrange=(n_rows, n_cols),
    )
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
    return packed_valid_pair_mean(reduction_execution_plan(signal), signal, valid_mask)
end

function packed_valid_pair_mean(::HostMirrorReductionPlan, signal::AbstractVector{T}, valid_mask::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    return packed_valid_pair_mean(ScalarCPUStyle(), Array(signal), Array(valid_mask))
end

function packed_valid_pair_mean(::DirectReductionPlan, signal::AbstractVector{T}, valid_mask::AbstractMatrix{Bool}) where {T<:AbstractFloat}
    count = sum(valid_mask)
    return count == 0 ? one(T) : sum(signal) / (T(2) * T(count))
end

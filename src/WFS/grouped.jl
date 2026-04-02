abstract type AbstractGroupedAccumulationPlan end

struct GroupedStackReducePlan <: AbstractGroupedAccumulationPlan end
struct GroupedDirectAccumulatePlan <: AbstractGroupedAccumulationPlan end
struct GroupedStaged2DPlan <: AbstractGroupedAccumulationPlan end

@kernel function reduce_grouped_stack_kernel!(out, stack, count::Int, n1::Int, n2::Int)
    i, j = @index(Global, NTuple)
    if i <= n1 && j <= n2
        acc = zero(eltype(out))
        @inbounds for idx in 1:count
            acc += stack[i, j, idx]
        end
        @inbounds out[i, j] = acc
    end
end

@kernel function reduce_grouped_blocks_kernel!(out, stack, block_size::Int, n_blocks::Int, n1::Int, n2::Int)
    i, j, k = @index(Global, NTuple)
    if i <= n1 && j <= n2 && k <= block_size
        acc = zero(eltype(out))
        @inbounds for block in 1:n_blocks
            acc += stack[i, j, k + (block - 1) * block_size]
        end
        @inbounds out[i, j, k] = acc
    end
end

@inline grouped_stack_view(stack::AbstractArray{T,3}, count::Int) where {T} = @view stack[:, :, 1:count]

@inline grouped_accumulation_plan(style::S, wfs::W) where {S<:ExecutionStyle,W<:AbstractWFS} =
    grouped_accumulation_plan(S, W)

@inline grouped_accumulation_plan(::Type{<:ExecutionStyle}, ::Type{<:AbstractWFS}) = GroupedStackReducePlan()

function reduce_grouped_stack!(::ScalarCPUStyle, out::AbstractMatrix, stack::AbstractArray{T,3}, count::Int) where {T}
    fill!(out, zero(eltype(out)))
    @inbounds for idx in 1:count
        out .+= @view(stack[:, :, idx])
    end
    return out
end

function reduce_grouped_stack!(style::AcceleratorStyle, out::AbstractMatrix, stack::AbstractArray{T,3}, count::Int) where {T}
    launch_kernel!(style, reduce_grouped_stack_kernel!, out, stack, count, size(out, 1), size(out, 2); ndrange=size(out))
    return out
end

function reduce_grouped_blocks!(style::AcceleratorStyle, out::AbstractArray{T,3}, stack::AbstractArray{T,3}, block_size::Int, n_blocks::Int) where {T}
    launch_kernel!(style, reduce_grouped_blocks_kernel!, out, stack, block_size, n_blocks,
        size(out, 1), size(out, 2); ndrange=size(out))
    return out
end

@inline function accumulate_grouped_sources!(::GroupedStackReducePlan, ::ScalarCPUStyle,
    out::AbstractMatrix, stack::AbstractArray{T,3}, sources, render!, args...) where {T}
    count = length(sources)
    @inbounds for idx in eachindex(sources)
        render!(@view(stack[:, :, idx]), args..., sources[idx])
    end
    return reduce_grouped_stack!(ScalarCPUStyle(), out, stack, count)
end

@inline function accumulate_grouped_sources!(::GroupedStackReducePlan, style::AcceleratorStyle,
    out::AbstractMatrix, stack::AbstractArray{T,3}, sources, render!, args...) where {T}
    count = length(sources)
    @inbounds for idx in eachindex(sources)
        render!(@view(stack[:, :, idx]), args..., sources[idx])
    end
    return reduce_grouped_stack!(style, out, stack, count)
end

@inline function accumulate_grouped_sources!(style::S, wfs::W, out::AbstractMatrix,
    stack::AbstractArray{T,3}, sources, render!, args...) where {S<:ExecutionStyle,W<:AbstractWFS,T}
    plan = grouped_accumulation_plan(style, wfs)
    return accumulate_grouped_sources!(plan, style, out, stack, sources, render!, args...)
end

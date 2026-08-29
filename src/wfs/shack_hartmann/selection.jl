struct _OwnedShackHartmannSlopeSelectionPlan end

const _OWNED_SHACK_HARTMANN_SLOPE_SELECTION_PLAN =
    _OwnedShackHartmannSlopeSelectionPlan()

"""
    ShackHartmannSlopeSelectionPlan(n_lenslets, lenslet_order)

Prepare a run-immutable mapping from the full square Shack–Hartmann lenslet
grid to one ordered compact slope vector. `lenslet_order` contains unique,
one-based linear lenslet indices in the package's Julia column-major grid
order. Construction snapshots the order on its existing compute device.
"""
struct ShackHartmannSlopeSelectionPlan{
    I<:Integer,
    A<:AbstractVector{I},
}
    n_lenslets::Int
    lenslet_order::A

    function ShackHartmannSlopeSelectionPlan(
        ::_OwnedShackHartmannSlopeSelectionPlan,
        n_lenslets::Int,
        lenslet_order::A,
    ) where {I<:Integer,A<:AbstractVector{I}}
        return new{I,A}(n_lenslets, lenslet_order)
    end
end

function ShackHartmannSlopeSelectionPlan(
    n_lenslets::Integer,
    lenslet_order::AbstractVector{I},
) where {I<:Integer}
    Base.require_one_based_indexing(lenslet_order)
    I === Bool && throw(InvalidConfiguration(
        "Shack-Hartmann lenslet order must contain integer indices, not Bool",
    ))
    lenslets = Int(n_lenslets)
    lenslets > 0 || throw(InvalidConfiguration(
        "Shack-Hartmann slope-selection lenslet count must be positive",
    ))
    isempty(lenslet_order) && throw(InvalidConfiguration(
        "Shack-Hartmann lenslet order must not be empty",
    ))
    full_count = Base.checked_mul(lenslets, lenslets)
    length(lenslet_order) <= full_count || throw(DimensionMismatchError(
        "Shack-Hartmann lenslet order cannot exceed the full lenslet grid",
    ))

    host_order = Array(lenslet_order)
    for (position, lenslet) in pairs(host_order)
        one(I) <= lenslet <= full_count || throw(InvalidConfiguration(
            "Shack-Hartmann lenslet-order entry $position is outside 1:$full_count",
        ))
    end
    allunique(host_order) || throw(InvalidConfiguration(
        "Shack-Hartmann lenslet order must not contain duplicates",
    ))
    return ShackHartmannSlopeSelectionPlan(
        _OWNED_SHACK_HARTMANN_SLOPE_SELECTION_PLAN,
        lenslets,
        copy(lenslet_order),
    )
end

@inline function selected_lenslet_count(
    plan::ShackHartmannSlopeSelectionPlan,
)
    return length(plan.lenslet_order)
end

@inline function full_lenslet_count(
    plan::ShackHartmannSlopeSelectionPlan,
)
    return plan.n_lenslets * plan.n_lenslets
end

@kernel function _select_shack_hartmann_slopes_kernel!(
    selected_slopes,
    full_slopes,
    lenslet_order,
    full_count::Int,
    selected_count::Int,
)
    selected_index = @index(Global, Linear)
    if selected_index <= selected_count
        lenslet_index = Int(@inbounds lenslet_order[selected_index])
        output_index = 2 * selected_index - 1
        @inbounds begin
            selected_slopes[output_index] = full_slopes[lenslet_index]
            selected_slopes[output_index + 1] =
                full_slopes[full_count + lenslet_index]
        end
    end
end

@inline function _select_shack_hartmann_slopes!(
    ::ScalarCPUStyle,
    selected_slopes::AbstractVector,
    full_slopes::AbstractVector,
    lenslet_order::AbstractVector,
    full_count::Int,
    selected_count::Int,
)
    @inbounds @simd for selected_index in 1:selected_count
        lenslet_index = Int(lenslet_order[selected_index])
        output_index = 2 * selected_index - 1
        selected_slopes[output_index] = full_slopes[lenslet_index]
        selected_slopes[output_index + 1] =
            full_slopes[full_count + lenslet_index]
    end
    return selected_slopes
end

@inline function _select_shack_hartmann_slopes!(
    style::AcceleratorStyle,
    selected_slopes::AbstractVector,
    full_slopes::AbstractVector,
    lenslet_order::AbstractVector,
    full_count::Int,
    selected_count::Int,
)
    launch_kernel!(
        style,
        _select_shack_hartmann_slopes_kernel!,
        selected_slopes,
        full_slopes,
        lenslet_order,
        full_count,
        selected_count;
        ndrange=selected_count,
    )
    return selected_slopes
end

"""
    select_shack_hartmann_slopes!(selected_slopes, plan, full_slopes)

Gather a full `[axis 1; axis 2]` Shack–Hartmann measurement into a compact
pair-interleaved `[axis 1, axis 2]` vector in the plan's explicit lenslet
order. Input, output, and prepared order must occupy one compute device and
must not alias. Successful warmed CPU calls allocate no Julia heap storage.
"""
function select_shack_hartmann_slopes!(
    selected_slopes::AbstractVector{T},
    plan::ShackHartmannSlopeSelectionPlan,
    full_slopes::AbstractVector{T},
) where {T<:AbstractFloat}
    Base.require_one_based_indexing(selected_slopes, full_slopes)
    full_count = full_lenslet_count(plan)
    selected_count = selected_lenslet_count(plan)
    length(full_slopes) == 2 * full_count || throw(DimensionMismatchError(
        "full Shack-Hartmann slopes must contain two complete lenslet blocks",
    ))
    length(selected_slopes) == 2 * selected_count || throw(
        DimensionMismatchError(
            "selected Shack-Hartmann slopes must contain one pair per selected lenslet",
        ),
    )
    _wfs_storage_mightalias(selected_slopes, full_slopes) && throw(
        InvalidConfiguration(
            "selected and full Shack-Hartmann slope storage must not alias",
        ),
    )
    device = compute_device(selected_slopes)
    device == compute_device(full_slopes) &&
        device == compute_device(plan.lenslet_order) || throw(
        InvalidConfiguration(
            "selected slopes, full slopes, and lenslet order must occupy one compute device",
        ),
    )
    return _select_shack_hartmann_slopes!(
        execution_style(selected_slopes),
        selected_slopes,
        full_slopes,
        plan.lenslet_order,
        full_count,
        selected_count,
    )
end

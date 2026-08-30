struct _OwnedClosedLoopCorrectionStorage end

const _OWNED_CLOSED_LOOP_CORRECTION_STORAGE =
    _OwnedClosedLoopCorrectionStorage()

"""
    ClosedLoopCorrectionPlan(extent, gain, pole, anti_windup_gain)

Run-immutable coefficients for discrete leaky integration with delayed
same-coordinate constraint feedback. For preceding correction `u`, feedback
`f = u - realized_correction`, and current reconstructed residual error `x`,
updates after the first are

```text
state = u - anti_windup_gain * f
correction = pole * state + gain * x
```

The first update starts from zero state and ignores its placeholder feedback.
This plan does not include controller/DM coordinate transforms, hidden-mode
removal, command constraints, or physical DM response.
"""
struct ClosedLoopCorrectionPlan{T<:AbstractFloat}
    extent::Int
    gain::T
    pole::T
    anti_windup_gain::T

    function ClosedLoopCorrectionPlan(
        extent::Integer,
        gain::T,
        pole::T,
        anti_windup_gain::T,
    ) where {T<:AbstractFloat}
        extent > 0 || throw(InvalidConfiguration(
            "closed-loop correction extent must be positive",
        ))
        isfinite(gain) || throw(InvalidConfiguration(
            "closed-loop correction gain must be finite",
        ))
        isfinite(pole) && zero(T) <= pole <= one(T) || throw(
            InvalidConfiguration(
                "closed-loop correction pole must be finite and lie in [0, 1]",
            ),
        )
        isfinite(anti_windup_gain) || throw(InvalidConfiguration(
            "closed-loop correction anti-windup gain must be finite",
        ))
        return new{T}(Int(extent), gain, pole, anti_windup_gain)
    end
end

"""Persistent mathematical state for one closed-loop correction run."""
mutable struct ClosedLoopCorrectionState{
    T<:AbstractFloat,
    V<:AbstractVector{T},
}
    integrator_state::V
    last_correction::V
    has_correction::Bool

    function ClosedLoopCorrectionState(
        ::_OwnedClosedLoopCorrectionStorage,
        integrator_state::V,
        last_correction::V,
    ) where {T<:AbstractFloat,V<:AbstractVector{T}}
        return new{T,V}(integrator_state, last_correction, false)
    end
end

"""
Replaceable scratch for one closed-loop correction execution owner.

Recreating this workspace cannot change the controller trajectory.
"""
struct ClosedLoopCorrectionWorkspace{
    T<:AbstractFloat,
    V<:AbstractVector{T},
}
    next_state::V
    next_correction::V

    function ClosedLoopCorrectionWorkspace(
        ::_OwnedClosedLoopCorrectionStorage,
        next_state::V,
        next_correction::V,
    ) where {T<:AbstractFloat,V<:AbstractVector{T}}
        return new{T,V}(next_state, next_correction)
    end
end

function _closed_loop_storage(
    plan::ClosedLoopCorrectionPlan{T},
    reference::AbstractVector{T},
) where {T<:AbstractFloat}
    Base.require_one_based_indexing(reference)
    length(reference) == plan.extent || throw(DimensionMismatchError(
        "closed-loop correction reference length must match plan extent",
    ))
    first = similar(reference, T, plan.extent)
    second = similar(reference, T, plan.extent)
    fill!(first, zero(T))
    fill!(second, zero(T))
    return first, second
end

function ClosedLoopCorrectionState(
    plan::ClosedLoopCorrectionPlan{T},
    reference::AbstractVector{T},
) where {T<:AbstractFloat}
    integrator_state, last_correction = _closed_loop_storage(plan, reference)
    return ClosedLoopCorrectionState(
        _OWNED_CLOSED_LOOP_CORRECTION_STORAGE,
        integrator_state,
        last_correction,
    )
end

function ClosedLoopCorrectionWorkspace(
    plan::ClosedLoopCorrectionPlan{T},
    reference::AbstractVector{T},
) where {T<:AbstractFloat}
    next_state, next_correction = _closed_loop_storage(plan, reference)
    return ClosedLoopCorrectionWorkspace(
        _OWNED_CLOSED_LOOP_CORRECTION_STORAGE,
        next_state,
        next_correction,
    )
end

@inline function _closed_loop_mightalias(
    value::AbstractArray,
    others::Tuple,
)
    isempty(others) && return false
    return Base.mightalias(value, first(others)) ||
           _closed_loop_mightalias(value, Base.tail(others))
end

@inline _closed_loop_mightalias(::AbstractArray, ::Tuple{}) = false

function _validate_closed_loop_correction(
    correction::AbstractVector,
    state::ClosedLoopCorrectionState,
    workspace::ClosedLoopCorrectionWorkspace,
    plan::ClosedLoopCorrectionPlan,
    residual_error::AbstractVector,
    constraint_feedback::AbstractVector,
)
    Base.require_one_based_indexing(
        correction,
        state.integrator_state,
        state.last_correction,
        workspace.next_state,
        workspace.next_correction,
        residual_error,
        constraint_feedback,
    )
    expected_axes = (Base.OneTo(plan.extent),)
    all(==(expected_axes), (
        axes(correction),
        axes(state.integrator_state),
        axes(state.last_correction),
        axes(workspace.next_state),
        axes(workspace.next_correction),
        axes(residual_error),
        axes(constraint_feedback),
    )) || throw(DimensionMismatchError(
        "closed-loop correction storage axes must match the plan extent",
    ))
    mutable_storage = (
        state.integrator_state,
        state.last_correction,
        workspace.next_state,
        workspace.next_correction,
    )
    _closed_loop_mightalias(correction, (
        residual_error,
        constraint_feedback,
        mutable_storage...,
    )) && throw(InvalidConfiguration(
        "closed-loop correction output must not alias input, state, or workspace storage",
    ))
    _closed_loop_mightalias(residual_error, mutable_storage) && throw(
        InvalidConfiguration(
            "closed-loop residual input must not alias state or workspace storage",
        ),
    )
    _closed_loop_mightalias(constraint_feedback, mutable_storage) && throw(
        InvalidConfiguration(
            "closed-loop constraint feedback must not alias state or workspace storage",
        ),
    )
    device = compute_device(correction)
    all(==(device), (
        compute_device(state.integrator_state),
        compute_device(state.last_correction),
        compute_device(workspace.next_state),
        compute_device(workspace.next_correction),
        compute_device(residual_error),
        compute_device(constraint_feedback),
    )) || throw(InvalidConfiguration(
        "closed-loop correction storage must occupy one compute device",
    ))
    return nothing
end

"""
    apply_closed_loop_correction!(correction, state, workspace, plan,
                                  residual_error, constraint_feedback)

Advance one discrete closed-loop correction sample and write the complete
caller-owned correction. `constraint_feedback` describes the preceding
demanded-minus-realized correction in the same coordinates; its value is
ignored on the first update, when no preceding correction exists. State
commits only after the complete next correction has been formed. Successful
warmed CPU calls allocate no Julia heap storage.
"""
function apply_closed_loop_correction!(
    correction::AbstractVector{T},
    state::ClosedLoopCorrectionState{T},
    workspace::ClosedLoopCorrectionWorkspace{T},
    plan::ClosedLoopCorrectionPlan{T},
    residual_error::AbstractVector{T},
    constraint_feedback::AbstractVector{T},
) where {T<:AbstractFloat}
    _validate_closed_loop_correction(
        correction,
        state,
        workspace,
        plan,
        residual_error,
        constraint_feedback,
    )
    if state.has_correction
        @. workspace.next_state =
            state.last_correction -
            plan.anti_windup_gain * constraint_feedback
    else
        copyto!(workspace.next_state, state.integrator_state)
    end
    @. workspace.next_correction =
        plan.pole * workspace.next_state + plan.gain * residual_error

    copyto!(correction, workspace.next_correction)
    copyto!(state.integrator_state, workspace.next_state)
    copyto!(state.last_correction, workspace.next_correction)
    state.has_correction = true
    return correction
end

"""Reset persistent correction state and replaceable scratch to zero."""
function reset_closed_loop_correction!(
    state::ClosedLoopCorrectionState,
    workspace::ClosedLoopCorrectionWorkspace,
)
    fill!(state.integrator_state, zero(eltype(state.integrator_state)))
    fill!(state.last_correction, zero(eltype(state.last_correction)))
    fill!(workspace.next_state, zero(eltype(workspace.next_state)))
    fill!(workspace.next_correction, zero(eltype(workspace.next_correction)))
    state.has_correction = false
    return state
end

struct _OwnedDMRoutingStorage end

const _OWNED_DM_ROUTING_STORAGE = _OwnedDMRoutingStorage()

function _snapshot_projection_matrix(
    matrix::AbstractMatrix{T},
    role::AbstractString,
) where {T<:AbstractFloat}
    Base.require_one_based_indexing(matrix)
    size(matrix, 1) > 0 || throw(InvalidConfiguration(
        "$role must contain at least one output coordinate",
    ))
    size(matrix, 2) > 0 || throw(InvalidConfiguration(
        "$role must contain at least one input coordinate",
    ))
    all(isfinite, host_array(matrix)) || throw(InvalidConfiguration(
        "$role coefficients must be finite",
    ))
    return copy(matrix)
end

function _validate_projection_application(
    output::AbstractVector,
    matrix::AbstractMatrix,
    input::AbstractVector,
    role::AbstractString,
)
    Base.require_one_based_indexing(output, matrix, input)
    length(output) == size(matrix, 1) || throw(DimensionMismatchError(
        "$role output length must match projection row count",
    ))
    length(input) == size(matrix, 2) || throw(DimensionMismatchError(
        "$role input length must match projection column count",
    ))
    Base.mightalias(output, input) && throw(InvalidConfiguration(
        "$role output must not alias its input",
    ))
    Base.mightalias(output, matrix) && throw(InvalidConfiguration(
        "$role output must not alias prepared projection storage",
    ))
    device = compute_device(output)
    device == compute_device(input) && device == compute_device(matrix) ||
        throw(InvalidConfiguration(
            "$role output, input, and projection must occupy one compute device",
        ))
    return nothing
end

@inline function _apply_projection!(
    output::AbstractVector,
    matrix::AbstractMatrix,
    input::AbstractVector,
    role::AbstractString,
)
    _validate_projection_application(output, matrix, input, role)
    mul!(output, matrix, input)
    return output
end

"""
    ControllerToVDMPlan(controller_to_vdm)

Run-immutable projection from concatenated controller-command coordinates to
one active virtual-deformable-mirror (VDM) coordinate vector. Construction
snapshots the already calibrated matrix on its existing compute device.
"""
struct ControllerToVDMPlan{
    T<:AbstractFloat,
    M<:AbstractMatrix{T},
}
    controller_to_vdm::M

    function ControllerToVDMPlan(
        ::_OwnedDMRoutingStorage,
        controller_to_vdm::M,
    ) where {T<:AbstractFloat,M<:AbstractMatrix{T}}
        return new{T,M}(controller_to_vdm)
    end
end

function ControllerToVDMPlan(
    controller_to_vdm::AbstractMatrix{T},
) where {T<:AbstractFloat}
    matrix = _snapshot_projection_matrix(
        controller_to_vdm,
        "controller-to-VDM projection",
    )
    return ControllerToVDMPlan(_OWNED_DM_ROUTING_STORAGE, matrix)
end

"""
    project_controller_to_vdm!(vdm_command, plan, controller_command)

Project one complete controller command into active VDM coordinates. All
arrays must occupy one compute device. A successful warmed CPU call allocates
no Julia heap storage.
"""
function project_controller_to_vdm!(
    vdm_command::AbstractVector{T},
    plan::ControllerToVDMPlan{T},
    controller_command::AbstractVector{T},
) where {T<:AbstractFloat}
    return _apply_projection!(
        vdm_command,
        plan.controller_to_vdm,
        controller_command,
        "controller-to-VDM projection",
    )
end

"""
    VDMToPDMPlan(active_to_full, vdm_to_pdm)

Run-immutable two-stage projection from active VDM coordinates through the
full VDM space into concatenated physical-deformable-mirror (PDM) actuator
coordinates. Both calibrated matrices are snapshotted during construction.
"""
struct VDMToPDMPlan{
    T<:AbstractFloat,
    A<:AbstractMatrix{T},
    P<:AbstractMatrix{T},
}
    active_to_full::A
    vdm_to_pdm::P

    function VDMToPDMPlan(
        ::_OwnedDMRoutingStorage,
        active_to_full::A,
        vdm_to_pdm::P,
    ) where {
        T<:AbstractFloat,
        A<:AbstractMatrix{T},
        P<:AbstractMatrix{T},
    }
        return new{T,A,P}(active_to_full, vdm_to_pdm)
    end
end

function VDMToPDMPlan(
    active_to_full::AbstractMatrix{T},
    vdm_to_pdm::AbstractMatrix{T},
) where {T<:AbstractFloat}
    active_to_full_snapshot = _snapshot_projection_matrix(
        active_to_full,
        "active-to-full VDM projection",
    )
    vdm_to_pdm_snapshot = _snapshot_projection_matrix(
        vdm_to_pdm,
        "VDM-to-PDM projection",
    )
    size(vdm_to_pdm_snapshot, 2) == size(active_to_full_snapshot, 1) ||
        throw(DimensionMismatchError(
            "VDM-to-PDM projection column count must match full VDM extent",
        ))
    compute_device(active_to_full_snapshot) ==
        compute_device(vdm_to_pdm_snapshot) || throw(InvalidConfiguration(
        "VDM-to-PDM projection matrices must occupy one compute device",
    ))
    return VDMToPDMPlan(
        _OWNED_DM_ROUTING_STORAGE,
        active_to_full_snapshot,
        vdm_to_pdm_snapshot,
    )
end

"""Replaceable full-VDM scratch for one VDM-to-PDM operation owner."""
struct VDMToPDMWorkspace{T<:AbstractFloat,V<:AbstractVector{T}}
    full_vdm_command::V
end

function VDMToPDMWorkspace(
    plan::VDMToPDMPlan{T},
    reference::AbstractVector{T},
) where {T<:AbstractFloat}
    Base.require_one_based_indexing(reference)
    compute_device(reference) == compute_device(plan.active_to_full) ||
        throw(InvalidConfiguration(
            "VDM-to-PDM workspace reference and plan must occupy one compute device",
        ))
    full_vdm_command = similar(reference, T, size(plan.active_to_full, 1))
    fill!(full_vdm_command, zero(T))
    return VDMToPDMWorkspace(full_vdm_command)
end

"""
    project_vdm_to_pdm!(requested_pdm_command, workspace, plan, vdm_command)

Project a complete active-VDM command through the full VDM space into physical
actuator coordinates. The caller retains the output; `workspace` owns only the
replaceable intermediate full-VDM vector.
"""
function project_vdm_to_pdm!(
    requested_pdm_command::AbstractVector{T},
    workspace::VDMToPDMWorkspace{T},
    plan::VDMToPDMPlan{T},
    vdm_command::AbstractVector{T},
) where {T<:AbstractFloat}
    _validate_projection_application(
        workspace.full_vdm_command,
        plan.active_to_full,
        vdm_command,
        "active-to-full VDM projection",
    )
    _validate_projection_application(
        requested_pdm_command,
        plan.vdm_to_pdm,
        workspace.full_vdm_command,
        "VDM-to-PDM projection",
    )
    Base.mightalias(requested_pdm_command, vdm_command) && throw(
        InvalidConfiguration(
            "requested PDM command must not alias active VDM command storage",
        ),
    )
    mul!(
        workspace.full_vdm_command,
        plan.active_to_full,
        vdm_command,
    )
    mul!(
        requested_pdm_command,
        plan.vdm_to_pdm,
        workspace.full_vdm_command,
    )
    return requested_pdm_command
end

"""
    PDMActuatorRangePlan(lower_limits, upper_limits)

Run-immutable per-actuator physical-DM command interval. Construction
snapshots the lower and upper limit vectors on their existing compute device.
This plan performs range conditioning only; it does not imply slew,
quantization, power, inter-actuator stroke, or stroke-optimization behavior.
"""
struct PDMActuatorRangePlan{
    T<:AbstractFloat,
    L<:AbstractVector{T},
    U<:AbstractVector{T},
}
    lower_limits::L
    upper_limits::U

    function PDMActuatorRangePlan(
        ::_OwnedDMRoutingStorage,
        lower_limits::L,
        upper_limits::U,
    ) where {
        T<:AbstractFloat,
        L<:AbstractVector{T},
        U<:AbstractVector{T},
    }
        return new{T,L,U}(lower_limits, upper_limits)
    end
end


function PDMActuatorRangePlan(
    lower_limits::AbstractVector{T},
    upper_limits::AbstractVector{T},
) where {T<:AbstractFloat}
    Base.require_one_based_indexing(lower_limits, upper_limits)
    isempty(lower_limits) && throw(InvalidConfiguration(
        "PDM actuator range must contain at least one actuator",
    ))
    axes(lower_limits) == axes(upper_limits) || throw(DimensionMismatchError(
        "PDM lower and upper actuator limits must have identical axes",
    ))
    compute_device(lower_limits) == compute_device(upper_limits) || throw(
        InvalidConfiguration(
            "PDM lower and upper actuator limits must occupy one compute device",
        ),
    )
    lower_host = host_array(lower_limits)
    upper_host = host_array(upper_limits)
    all(isfinite, lower_host) && all(isfinite, upper_host) || throw(
        InvalidConfiguration("PDM actuator limits must be finite"),
    )
    all(lower_host .<= upper_host) || throw(InvalidConfiguration(
        "each PDM lower actuator limit must not exceed its upper limit",
    ))
    return PDMActuatorRangePlan(
        _OWNED_DM_ROUTING_STORAGE,
        copy(lower_limits),
        copy(upper_limits),
    )
end

"""
    apply_pdm_actuator_range!(demanded, constraint_feedback, plan, requested)

Clamp one complete requested physical-DM command to the prepared per-actuator
range and publish both the demanded command and
`constraint_feedback = requested - demanded`. The two outputs and input must
be disjoint and occupy the plan's compute device.
"""
function apply_pdm_actuator_range!(
    demanded::AbstractVector{T},
    constraint_feedback::AbstractVector{T},
    plan::PDMActuatorRangePlan{T},
    requested::AbstractVector{T},
) where {T<:AbstractFloat}
    Base.require_one_based_indexing(
        demanded,
        constraint_feedback,
        requested,
    )
    expected_axes = axes(plan.lower_limits)
    all(==(expected_axes), (
        axes(plan.upper_limits),
        axes(demanded),
        axes(constraint_feedback),
        axes(requested),
    )) || throw(DimensionMismatchError(
        "PDM command and actuator-limit axes must match",
    ))
    Base.mightalias(demanded, constraint_feedback) && throw(
        InvalidConfiguration("PDM command outputs must not alias"),
    )
    Base.mightalias(demanded, requested) && throw(InvalidConfiguration(
        "demanded PDM command must not alias requested command storage",
    ))
    Base.mightalias(constraint_feedback, requested) && throw(
        InvalidConfiguration(
            "PDM constraint feedback must not alias requested command storage",
        ),
    )
    aliases_limits =
        Base.mightalias(demanded, plan.lower_limits) ||
        Base.mightalias(demanded, plan.upper_limits) ||
        Base.mightalias(constraint_feedback, plan.lower_limits) ||
        Base.mightalias(constraint_feedback, plan.upper_limits)
    aliases_limits && throw(InvalidConfiguration(
            "PDM command outputs must not alias prepared actuator limits",
        ))
    device = compute_device(demanded)
    all(==(device), (
        compute_device(constraint_feedback),
        compute_device(requested),
        compute_device(plan.lower_limits),
        compute_device(plan.upper_limits),
    )) || throw(InvalidConfiguration(
        "PDM command storage and actuator limits must occupy one compute device",
    ))
    @. demanded = clamp(requested, plan.lower_limits, plan.upper_limits)
    @. constraint_feedback = requested - demanded
    return demanded
end

"""
    PDMFeedbackToVDMPlan(pdm_to_vdm, full_to_active)

Run-immutable two-stage reverse projection from concatenated PDM constraint
feedback through full VDM coordinates into active VDM coordinates.
"""
struct PDMFeedbackToVDMPlan{
    T<:AbstractFloat,
    P<:AbstractMatrix{T},
    A<:AbstractMatrix{T},
}
    pdm_to_vdm::P
    full_to_active::A

    function PDMFeedbackToVDMPlan(
        ::_OwnedDMRoutingStorage,
        pdm_to_vdm::P,
        full_to_active::A,
    ) where {
        T<:AbstractFloat,
        P<:AbstractMatrix{T},
        A<:AbstractMatrix{T},
    }
        return new{T,P,A}(pdm_to_vdm, full_to_active)
    end
end

function PDMFeedbackToVDMPlan(
    pdm_to_vdm::AbstractMatrix{T},
    full_to_active::AbstractMatrix{T},
) where {T<:AbstractFloat}
    pdm_to_vdm_snapshot = _snapshot_projection_matrix(
        pdm_to_vdm,
        "PDM-to-VDM feedback projection",
    )
    full_to_active_snapshot = _snapshot_projection_matrix(
        full_to_active,
        "full-to-active VDM feedback projection",
    )
    size(full_to_active_snapshot, 2) == size(pdm_to_vdm_snapshot, 1) ||
        throw(DimensionMismatchError(
            "full-to-active VDM projection column count must match full VDM extent",
        ))
    compute_device(pdm_to_vdm_snapshot) ==
        compute_device(full_to_active_snapshot) || throw(InvalidConfiguration(
        "PDM feedback projection matrices must occupy one compute device",
    ))
    return PDMFeedbackToVDMPlan(
        _OWNED_DM_ROUTING_STORAGE,
        pdm_to_vdm_snapshot,
        full_to_active_snapshot,
    )
end

"""Replaceable full-VDM scratch for one PDM-feedback operation owner."""
struct PDMFeedbackToVDMWorkspace{T<:AbstractFloat,V<:AbstractVector{T}}
    full_vdm_feedback::V
end

function PDMFeedbackToVDMWorkspace(
    plan::PDMFeedbackToVDMPlan{T},
    reference::AbstractVector{T},
) where {T<:AbstractFloat}
    Base.require_one_based_indexing(reference)
    compute_device(reference) == compute_device(plan.pdm_to_vdm) || throw(
        InvalidConfiguration(
            "PDM-feedback workspace reference and plan must occupy one compute device",
        ),
    )
    full_vdm_feedback = similar(reference, T, size(plan.pdm_to_vdm, 1))
    fill!(full_vdm_feedback, zero(T))
    return PDMFeedbackToVDMWorkspace(full_vdm_feedback)
end

"""
    project_pdm_feedback_to_vdm!(vdm_constraint_feedback, workspace, plan,
                                 pdm_constraint_feedback)

Project physical-actuator constraint feedback into active VDM coordinates
through the prepared full-VDM intermediate space.
"""
function project_pdm_feedback_to_vdm!(
    vdm_constraint_feedback::AbstractVector{T},
    workspace::PDMFeedbackToVDMWorkspace{T},
    plan::PDMFeedbackToVDMPlan{T},
    pdm_constraint_feedback::AbstractVector{T},
) where {T<:AbstractFloat}
    _validate_projection_application(
        workspace.full_vdm_feedback,
        plan.pdm_to_vdm,
        pdm_constraint_feedback,
        "PDM-to-VDM feedback projection",
    )
    _validate_projection_application(
        vdm_constraint_feedback,
        plan.full_to_active,
        workspace.full_vdm_feedback,
        "full-to-active VDM feedback projection",
    )
    Base.mightalias(vdm_constraint_feedback, pdm_constraint_feedback) &&
        throw(InvalidConfiguration(
            "active VDM feedback must not alias PDM feedback storage",
        ))
    mul!(
        workspace.full_vdm_feedback,
        plan.pdm_to_vdm,
        pdm_constraint_feedback,
    )
    mul!(
        vdm_constraint_feedback,
        plan.full_to_active,
        workspace.full_vdm_feedback,
    )
    return vdm_constraint_feedback
end

"""
    VDMFeedbackToControllerPlan(vdm_to_controller)

Run-immutable projection from active-VDM constraint feedback into
concatenated controller coordinates.
"""
struct VDMFeedbackToControllerPlan{
    T<:AbstractFloat,
    M<:AbstractMatrix{T},
}
    vdm_to_controller::M

    function VDMFeedbackToControllerPlan(
        ::_OwnedDMRoutingStorage,
        vdm_to_controller::M,
    ) where {T<:AbstractFloat,M<:AbstractMatrix{T}}
        return new{T,M}(vdm_to_controller)
    end
end

function VDMFeedbackToControllerPlan(
    vdm_to_controller::AbstractMatrix{T},
) where {T<:AbstractFloat}
    matrix = _snapshot_projection_matrix(
        vdm_to_controller,
        "VDM-to-controller feedback projection",
    )
    return VDMFeedbackToControllerPlan(_OWNED_DM_ROUTING_STORAGE, matrix)
end

"""
    project_vdm_feedback_to_controller!(controller_constraint_feedback, plan,
                                        vdm_constraint_feedback)

Project one complete active-VDM constraint-feedback vector back into
controller coordinates.
"""
function project_vdm_feedback_to_controller!(
    controller_constraint_feedback::AbstractVector{T},
    plan::VDMFeedbackToControllerPlan{T},
    vdm_constraint_feedback::AbstractVector{T},
) where {T<:AbstractFloat}
    return _apply_projection!(
        controller_constraint_feedback,
        plan.vdm_to_controller,
        vdm_constraint_feedback,
        "VDM-to-controller feedback projection",
    )
end

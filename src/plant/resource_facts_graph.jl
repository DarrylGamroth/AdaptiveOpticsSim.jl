#
# Exact structural resource facts for one prepared event-loop graph
#
# This collector assigns every explicitly owned numerical allocation to one
# stable logical owner. Fixed host registries, scheduler/control state, RNG
# objects, and transport-independent bookkeeping remain outside the numerical
# structural-storage contract documented by `StructuralResourceReport`.
#

@inline function _append_structural_graph_fact!(
    facts::Vector{AbstractStructuralResourceFact},
    fact::AbstractStructuralResourceFact,
    target::AbstractComputeDevice,
)
    _require_fact_target(fact, target)
    if structural_resource_known(fact)
        push!(facts, fact)
    elseif structural_resource_unknown_reason(fact) == :owner_not_on_device
        # The logical owner remains part of every selected-device report even
        # when all of its numerical storage belongs to another explicit
        # partition.  Known zero keeps the owner graph complete without
        # weakening fail-closed handling for unsupported storage.
        push!(facts, KnownStructuralResourceFact(
            structural_resource_owner_id(fact),
            target,
            UInt64(0),
            UInt64(0),
            resource_estimate_method(fact),
        ))
    else
        push!(facts, fact)
    end
    return facts
end

function _append_structural_graph_fact!(
    ::Vector{AbstractStructuralResourceFact},
    fact,
    ::AbstractComputeDevice,
)
    _structural_resource_error(:event_loop, :invalid_fact,
        "event-loop structural reporter returned $(typeof(fact)); expected " *
        "an AbstractStructuralResourceFact")
end

function _event_loop_structural_telescope(
    prepared::PreparedPlantEventLoop,
)
    isempty(prepared.path_groups) && _structural_resource_error(
        :event_loop, :missing_path,
        "an exact event-loop resource report requires at least one path")
    telescope = first(prepared.path_groups).path.telescope
    @inbounds for group in prepared.path_groups
        group.path.telescope === telescope || _structural_resource_error(
            :event_loop, :inconsistent_telescope,
            "prepared event-loop paths do not retain one authoritative " *
            "telescope")
    end
    return telescope
end

@inline function _command_endpoint_structural_resource_fact(
    prepared::_PreparedPlantEventCommandEndpoint,
    state::CommandEndpointState,
    application::CommandApplicationState,
    workspace::CommandDispositionWorkspace,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    return _combine_structural_owner_facts(id, target, (
        structural_resource_fact(prepared, id, target),
        structural_resource_fact(state, id, target),
        structural_resource_fact(application, id, target),
        structural_resource_fact(workspace, id, target),
    ))
end

@inline function _controllable_optic_structural_resource_fact(
    prepared::PreparedControllableOptic,
    state,
    workspace,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    return _combine_structural_owner_facts(id, target, (
        structural_resource_fact(prepared, id, target),
        structural_resource_fact(state, id, target),
        structural_resource_fact(workspace, id, target),
    ))
end

function _require_structural_report_query_target(
    prepared::PreparedPlantEventLoop,
    target::AbstractComputeDevice,
)
    primary = prepared.target
    (target == primary || target == HostComputeDevice()) ||
        _structural_resource_error(:event_loop, :wrong_device,
            "event-loop resource facts may be queried only for the exact " *
            "plant target $primary or its explicit host partition; got " *
            "$target")
    return primary
end

function _append_command_endpoint_structural_facts!(
    facts::Vector{AbstractStructuralResourceFact},
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    target::AbstractComputeDevice,
)
    @inbounds for index in eachindex(prepared.command_endpoints)
        endpoint = prepared.command_endpoints[index]
        id = StructuralResourceOwnerID(
            :command_endpoint, endpoint.id.name)
        fact = _command_endpoint_structural_resource_fact(
            endpoint,
            state.command_endpoints[index],
            state.command_applications[index],
            workspace.command_endpoints[index],
            id,
            target,
        )
        _append_structural_graph_fact!(facts, fact, target)
    end
    return facts
end

function _append_controllable_optic_structural_facts!(
    facts::Vector{AbstractStructuralResourceFact},
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    target::AbstractComputeDevice,
)
    @inbounds for index in eachindex(prepared.optics)
        optic = prepared.optics[index]
        definition = _prepared_controllable_optic_definition(
            prepared.optics, optic)
        id = StructuralResourceOwnerID(
            :controllable_optic,
            controllable_optic_id(definition).name,
        )
        fact = _controllable_optic_structural_resource_fact(
            optic,
            state.controllable_optics[index],
            workspace.controllable_optics[index],
            id,
            target,
        )
        _append_structural_graph_fact!(facts, fact, target)
    end
    return facts
end

function _append_sampled_aberration_structural_facts!(
    facts::Vector{AbstractStructuralResourceFact},
    prepared::PreparedPlantEventLoop,
    target::AbstractComputeDevice,
)
    @inbounds for aberration in prepared.sampled_aberrations
        id = StructuralResourceOwnerID(
            :sampled_aberration, sampled_aberration_id(aberration).name)
        _append_structural_graph_fact!(facts,
            structural_resource_fact(aberration, id, target), target)
    end
    return facts
end

function _append_path_structural_facts!(
    facts::Vector{AbstractStructuralResourceFact},
    prepared::PreparedPlantEventLoop,
    target::AbstractComputeDevice,
)
    @inbounds for group in prepared.path_groups
        id = StructuralResourceOwnerID(:path, group.id.name)
        _append_structural_graph_fact!(facts,
            structural_resource_fact(group.path, id, target), target)
    end
    return facts
end

function _append_acquisition_structural_facts!(
    facts::Vector{AbstractStructuralResourceFact},
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    target::AbstractComputeDevice,
)
    @inbounds for index in eachindex(prepared.acquisitions)
        acquisition = prepared.acquisitions[index]
        id = StructuralResourceOwnerID(:acquisition, acquisition.id.name)
        fact = structural_resource_fact(
            acquisition, state.acquisitions[index], id, target)
        _append_structural_graph_fact!(facts, fact, target)
    end
    return facts
end

function _device_batch_structural_owner_id(
    prepared::PreparedPlantEventLoop,
    owner::PreparedDevicePathBatchOwner,
)
    isempty(owner.group_slots) && _structural_resource_error(
        :device_path_batch, :empty_owner,
        "a prepared device path-batch owner must contain at least one path")
    first_slot = Int(first(owner.group_slots))
    checkbounds(Bool, prepared.path_groups, first_slot) ||
        _structural_resource_error(:device_path_batch, :invalid_group_slot,
            "device path-batch owner references path slot $first_slot " *
            "outside the prepared event loop")
    first_id = @inbounds prepared.path_groups[first_slot].id
    return StructuralResourceOwnerID(
        _device_batch_structural_owner_category(owner.implementation),
        first_id.name,
    )
end

@inline _device_batch_structural_owner_category(
    ::_PreparedDirectImagingDevicePathBatch) = :direct_batch_workspace
@inline _device_batch_structural_owner_category(
    ::_PreparedWFSDevicePathBatch) = :wfs_batch_workspace

function _append_device_batch_structural_facts!(
    facts::Vector{AbstractStructuralResourceFact},
    prepared::PreparedPlantEventLoop,
    target::AbstractComputeDevice,
)
    @inbounds for owner in prepared.device_path_batch_owners
        id = _device_batch_structural_owner_id(prepared, owner)
        fact = structural_resource_fact(
            owner, prepared.atmosphere, id, target)
        _append_structural_graph_fact!(facts, fact, target)
    end
    return facts
end

"""
    require_exact_structural_resource_facts(
        prepared, state, workspace, target; opaque_reserves=())

Validate one prepared plant event loop and return its exact numerical
structural-storage report for `target`. The report derives stable owner IDs
from declared plant identities and uses one non-overlapping owner for the
telescope, atmosphere, each command endpoint, controllable optic, sampled
aberration, path, acquisition, and native device batch.

For an accelerator-prepared loop, `target` may be either that exact accelerator
or `HostComputeDevice()` to inspect explicit host mirrors and staging. Fixed
control-plane registries, scheduling state, RNG objects, and opaque provider
storage are not structural array facts. A declared logical owner with no array
storage in the selected partition remains visible as a known zero-byte fact;
unsupported ownership still fails closed. Supply provider storage explicitly
as `OpaqueResourceReserve` values.
"""
function require_exact_structural_resource_facts(
    prepared::PreparedPlantEventLoop,
    state::PlantEventLoopState,
    workspace::PlantEventLoopWorkspace,
    target::AbstractComputeDevice;
    opaque_reserves=(),
)
    primary = _require_structural_report_query_target(prepared, target)
    _require_exact_prepared_event_loop_target(prepared, primary)
    _require_exact_plant_event_loop_state_target(prepared, state, primary)
    _require_exact_plant_event_loop_workspace_target(
        prepared, workspace, primary)

    facts = AbstractStructuralResourceFact[]
    telescope = _event_loop_structural_telescope(prepared)
    _append_structural_graph_fact!(facts,
        structural_resource_fact(telescope,
            StructuralResourceOwnerID(:telescope, :primary), target),
        target)
    _append_structural_graph_fact!(facts,
        structural_resource_fact(prepared.atmosphere,
            StructuralResourceOwnerID(:atmosphere, :primary), target),
        target)
    _append_sampled_aberration_structural_facts!(facts, prepared, target)
    _append_command_endpoint_structural_facts!(
        facts, prepared, state, workspace, target)
    _append_controllable_optic_structural_facts!(
        facts, prepared, state, workspace, target)
    _append_path_structural_facts!(facts, prepared, target)
    _append_acquisition_structural_facts!(facts, prepared, state, target)
    _append_device_batch_structural_facts!(facts, prepared, target)
    return aggregate_structural_resource_facts(facts, target;
        opaque_reserves)
end

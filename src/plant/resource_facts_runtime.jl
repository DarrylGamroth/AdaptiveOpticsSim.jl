#
# Structural facts for command and controllable-optic runtime owners
#
# Each method lists owned storage explicitly. Reshaped views and arrays borrowed
# from a prepared plan are intentionally absent from the corresponding state or
# workspace fact.
#

@inline _structural_array_tuple_bytes(::AbstractComputeDevice, ::Tuple{}) =
    UInt64(0)

@inline function _structural_array_tuple_bytes(
    target::AbstractComputeDevice, arrays::Tuple)
    head_bytes = structural_array_bytes(first(arrays), target)
    tail_bytes = _structural_array_tuple_bytes(target, Base.tail(arrays))
    return _checked_resource_add(head_bytes, tail_bytes, :array_storage)
end

@inline _optional_structural_array_bytes(
    ::Nothing, ::AbstractComputeDevice) = UInt64(0)
@inline _optional_structural_array_bytes(array::AbstractArray,
    target::AbstractComputeDevice) = structural_array_bytes(array, target)

@inline function _structural_memory_element_bytes(memory::Memory,
    count::Integer)
    0 <= count <= length(memory) || _structural_resource_error(
        :array_storage, :invalid_element_count,
        "structural Memory element count must lie within its capacity")
    checked_count = _checked_resource_byte_count(count, :array_storage)
    element_bytes = _checked_resource_byte_count(
        Base.elsize(typeof(memory)), :array_storage)
    return _checked_resource_multiply(
        checked_count, element_bytes, :array_storage)
end

# A prepared native DM owns its immutable operator/factors. The runtime state
# borrows those arrays and owns only its surface and command buffers.
@inline _dm_topology_metadata_is_exact(::NamedTuple{(),Tuple{}}) = true
@inline _dm_topology_metadata_is_exact(metadata) = false

@inline _dm_influence_storage_is_exact(
    ::Union{GaussianInfluenceWidth,GaussianMechanicalCoupling},
    ::Any,
) = true

@inline function _dm_influence_storage_is_exact(
    model::DenseInfluenceMatrix,
    modes::AbstractMatrix,
)
    return model.modes === modes
end

@inline function _dm_influence_storage_is_exact(
    model::MeasuredInfluenceFunctions,
    modes::AbstractMatrix,
)
    return model.modes === modes &&
        _dm_topology_metadata_is_exact(model.metadata)
end

@inline _dm_influence_storage_is_exact(::AbstractDMInfluenceModel, ::Any) =
    false

function _dm_topology_host_bytes(topology::ActuatorGridTopology,
    target::HostComputeDevice)
    _dm_topology_metadata_is_exact(topology.metadata) || return nothing
    return _structural_array_tuple_bytes(target, (
        topology.coords,
        topology.active_coords,
        topology.valid_actuators,
        topology.active_indices,
    ))
end

function _dm_topology_host_bytes(topology::SampledActuatorTopology,
    target::HostComputeDevice)
    _dm_topology_metadata_is_exact(topology.metadata) || return nothing
    return _structural_array_tuple_bytes(target, (
        topology.coords,
        topology.active_coords,
        topology.valid_actuators,
        topology.active_indices,
    ))
end

@inline _dm_topology_host_bytes(::AbstractDMTopology,
    ::AbstractComputeDevice) = UInt64(0)

@inline function _dm_backend_mode_bytes(operator::GaussianInfluenceOperator,
    target::AbstractComputeDevice)
    compute_device(operator.pupil_backend) == target || return UInt64(0)
    return _structural_array_tuple_bytes(target,
        (operator.pupil_backend, operator.coordinates_backend))
end

@inline function _dm_backend_mode_bytes(modes::AbstractMatrix,
    target::AbstractComputeDevice)
    compute_device(modes) == target || return UInt64(0)
    return structural_array_bytes(modes, target)
end

@inline function _dm_host_mode_bytes(operator::GaussianInfluenceOperator,
    target::HostComputeDevice)
    compute_device(operator.pupil_backend) == target && return UInt64(0)
    return _structural_array_tuple_bytes(target,
        (operator.pupil_host, operator.coordinates_host))
end

@inline _dm_host_mode_bytes(::AbstractMatrix, ::AbstractComputeDevice) =
    UInt64(0)

@inline _dm_actuator_model_bytes(::LinearStaticActuators,
    ::AbstractComputeDevice) = UInt64(0)
@inline _dm_actuator_model_bytes(::ClippedActuators,
    ::AbstractComputeDevice) = UInt64(0)

@inline function _dm_actuator_model_bytes(model::ActuatorHealthMap,
    target::AbstractComputeDevice)
    compute_device(model.gains) == target || return UInt64(0)
    return structural_array_bytes(model.gains, target)
end

@inline _dm_actuator_stage_bytes(::AbstractComputeDevice, ::Tuple{}) =
    UInt64(0)

@inline function _dm_actuator_stage_bytes(target::AbstractComputeDevice,
    stages::Tuple)
    head_bytes = _dm_actuator_model_bytes(first(stages), target)
    tail_bytes = _dm_actuator_stage_bytes(target, Base.tail(stages))
    return _checked_resource_add(head_bytes, tail_bytes, :resident_bytes)
end

@inline function _dm_actuator_model_bytes(model::CompositeDMActuatorModel,
    target::AbstractComputeDevice)
    return _dm_actuator_stage_bytes(target, model.stages)
end

@inline function _prepared_dm_device_bytes(
    prepared::_PreparedPlantDeformableMirror,
    target::AbstractComputeDevice)
    bytes = _dm_backend_mode_bytes(prepared.modes, target)
    bytes = _checked_resource_add(bytes,
        _optional_structural_array_bytes(prepared.separable_x, target),
        :resident_bytes)
    bytes = _checked_resource_add(bytes,
        _optional_structural_array_bytes(prepared.separable_y_t, target),
        :resident_bytes)
    return _checked_resource_add(bytes,
        _dm_actuator_model_bytes(prepared.params.actuator_model, target),
        :resident_bytes)
end

function structural_resource_fact(
    prepared::_PreparedPlantDeformableMirror,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    _dm_influence_storage_is_exact(
        prepared.params.influence_model, prepared.modes) ||
        return UnknownStructuralResourceFact(
            id, target, :unsupported_influence_storage)
    device = prepared.surface_metadata.device
    is_device_target = target == device
    is_host_target = target == HostComputeDevice()
    (is_device_target || is_host_target) || return UnknownStructuralResourceFact(
        id, target, :owner_not_on_device)

    resident = is_device_target ?
        _prepared_dm_device_bytes(prepared, target) : UInt64(0)
    if is_host_target
        topology_bytes = _dm_topology_host_bytes(
            prepared.params.topology, HostComputeDevice())
        isnothing(topology_bytes) && return UnknownStructuralResourceFact(
            id, target, :unsupported_topology_metadata)
        resident = _checked_resource_add(resident, topology_bytes,
            :resident_bytes)
        resident = _checked_resource_add(resident,
            _dm_host_mode_bytes(prepared.modes, HostComputeDevice()),
            :resident_bytes)
    end
    return KnownStructuralResourceFact(id, target, resident, 0)
end

@inline function structural_resource_fact(optic::PreparedControllableOptic,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    return structural_resource_fact(
        controllable_optic_implementation(optic), id, target)
end

@inline function _dm_runtime_owned_bytes(dm::DeformableMirror,
    target::AbstractComputeDevice)
    compute_device(dm.state.opd) == target || return nothing
    bytes = _structural_array_tuple_bytes(target, (
        dm.state.opd,
        dm.state.coefs,
        dm.state.actuator_coefs,
    ))
    return _checked_resource_add(bytes,
        _optional_structural_array_bytes(dm.state.separable_tmp, target),
        :array_storage)
end

function structural_resource_fact(state::_PlantDeformableMirrorState,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    bytes = _dm_runtime_owned_bytes(state.active, target)
    isnothing(bytes) && return UnknownStructuralResourceFact(
        id, target, :owner_not_on_device)
    return KnownStructuralResourceFact(id, target, bytes, 0)
end

function structural_resource_fact(workspace::_PlantDeformableMirrorWorkspace,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    bytes = _dm_runtime_owned_bytes(workspace.staged, target)
    isnothing(bytes) && return UnknownStructuralResourceFact(
        id, target, :owner_not_on_device)
    return KnownStructuralResourceFact(id, target, 0, bytes)
end

# Prepared command bindings own detached initial/safe arrays. Scalar values are
# inline and therefore contribute zero structural array bytes.
@inline _command_seed_bytes(::Nothing, ::AbstractComputeDevice) = UInt64(0)
@inline function _command_seed_bytes(seed::AbstractArray,
    target::AbstractComputeDevice)
    compute_device(seed) == target || return UInt64(0)
    return structural_array_bytes(seed, target)
end
@inline _command_seed_bytes(seed, ::AbstractComputeDevice) = UInt64(0)

@inline _command_seed_target(seed::AbstractArray,
    target::AbstractComputeDevice) = compute_device(seed) == target
@inline _command_seed_target(seed, target::AbstractComputeDevice) =
    target == HostComputeDevice()

function structural_resource_fact(binding::_PreparedPlantCommandEndpoint,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    _command_seed_target(binding.initial_command, target) ||
        return UnknownStructuralResourceFact(
            id, target, :owner_not_on_device)
    resident = _command_seed_bytes(binding.initial_command, target)
    resident = _checked_resource_add(resident,
        _command_seed_bytes(binding.safe_command, target), :resident_bytes)
    return KnownStructuralResourceFact(id, target, resident, 0)
end

@inline function structural_resource_fact(
    endpoint::_PreparedPlantEventCommandEndpoint,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    return structural_resource_fact(endpoint.binding, id, target)
end

@inline function _command_state_host_control_bytes(
    state::CommandEndpointState, target::HostComputeDevice)
    return _structural_array_tuple_bytes(target,
        (state.slots, state.calendar, state.accepted_sequences))
end

@inline _command_state_host_control_bytes(::CommandEndpointState,
    ::AbstractComputeDevice) = UInt64(0)

function _command_payload_bytes(payloads::_ScalarCommandPayloadSlots,
    target::HostComputeDevice)
    return (
        present=true,
        resident=structural_array_bytes(payloads.values, target),
        workspace=structural_array_bytes(payloads.staging, target),
    )
end

@inline _command_payload_bytes(::_ScalarCommandPayloadSlots,
    ::AbstractComputeDevice) =
    (present=false, resident=UInt64(0), workspace=UInt64(0))

function _command_payload_bytes(payloads::_ArrayCommandPayloadSlots,
    target::AbstractComputeDevice)
    values = payloads.values
    array_device = compute_device(first(values))
    is_array_target = target == array_device
    is_host_target = target == HostComputeDevice()
    (is_array_target || is_host_target) || return (
        present=false, resident=UInt64(0), workspace=UInt64(0))

    capacity = length(values) - 1
    resident = is_host_target ?
        _structural_memory_element_bytes(values, capacity) : UInt64(0)
    workspace = is_host_target ?
        _structural_memory_element_bytes(values, 1) : UInt64(0)
    if is_array_target
        for index in 1:capacity
            resident = _checked_resource_add(resident,
                structural_array_bytes(values[index], target),
                :resident_bytes)
        end
        workspace = _checked_resource_add(workspace,
            structural_array_bytes(values[payloads.staging_slot], target),
            :workspace_bytes)
    end
    return (; present=true, resident, workspace)
end

function structural_resource_fact(state::CommandEndpointState,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    payload = _command_payload_bytes(state.payloads, target)
    host_bytes = _command_state_host_control_bytes(state, target)
    (payload.present || !iszero(host_bytes)) ||
        return UnknownStructuralResourceFact(
            id, target, :owner_not_on_device)
    resident = _checked_resource_add(
        payload.resident, host_bytes, :resident_bytes)
    return KnownStructuralResourceFact(
        id, target, resident, payload.workspace)
end

function _command_application_bytes(
    values::_ScalarCommandApplicationValues,
    target::HostComputeDevice)
    return (present=true, resident=UInt64(0), workspace=UInt64(0))
end

@inline _command_application_bytes(::_ScalarCommandApplicationValues,
    ::AbstractComputeDevice) =
    (present=false, resident=UInt64(0), workspace=UInt64(0))

function _command_application_bytes(values::_ArrayCommandApplicationValues,
    target::AbstractComputeDevice)
    compute_device(values.effective) == target || return (
        present=false, resident=UInt64(0), workspace=UInt64(0))
    resident = structural_array_bytes(values.effective, target)
    resident = _checked_resource_add(resident,
        _optional_structural_array_bytes(values.safe, target),
        :resident_bytes)
    workspace = structural_array_bytes(values.staging, target)
    return (; present=true, resident, workspace)
end

function structural_resource_fact(state::CommandApplicationState,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    bytes = _command_application_bytes(state.values, target)
    bytes.present || return UnknownStructuralResourceFact(
        id, target, :owner_not_on_device)
    return KnownStructuralResourceFact(
        id, target, bytes.resident, bytes.workspace)
end

@inline function _target_local_command_replica_bytes(
    ::_TargetLocalScalarEffectiveCommandValues,
    ::AbstractComputeDevice,
)
    return (present=true, resident=UInt64(0), workspace=UInt64(0))
end

function _target_local_command_replica_bytes(
    values::_TargetLocalArrayEffectiveCommandValues,
    target::AbstractComputeDevice,
)
    compute_device(values.active) == target || return (
        present=false, resident=UInt64(0), workspace=UInt64(0))
    compute_device(values.staging) == target || return (
        present=false, resident=UInt64(0), workspace=UInt64(0))
    return (
        present=true,
        resident=structural_array_bytes(values.active, target),
        workspace=structural_array_bytes(values.staging, target),
    )
end

function structural_resource_fact(
    state::TargetLocalCommandEndpointState,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    bytes = _target_local_command_replica_bytes(state.values, target)
    bytes.present || return UnknownStructuralResourceFact(
        id, target, :owner_not_on_device)
    return KnownStructuralResourceFact(
        id, target, bytes.resident, bytes.workspace)
end

@inline function structural_resource_fact(
    prepared::PreparedTargetLocalControllableOptic,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    return structural_resource_fact(prepared.implementation, id, target)
end

@inline function structural_resource_fact(
    state::TargetLocalControllableOpticState,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    return structural_resource_fact(state.physical, id, target)
end

@inline function structural_resource_fact(
    workspace::TargetLocalControllableOpticWorkspace,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    return structural_resource_fact(workspace.physical, id, target)
end

@inline function structural_resource_fact(
    ::PreparedCircularPyramidModulator,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    return KnownStructuralResourceFact(id, target, 0, 0)
end

@inline function structural_resource_fact(
    ::CircularPyramidModulatorState,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    return KnownStructuralResourceFact(id, target, 0, 0)
end

@inline function structural_resource_fact(
    ::CircularPyramidModulatorWorkspace,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    return KnownStructuralResourceFact(id, target, 0, 0)
end

function structural_resource_fact(workspace::CommandDispositionWorkspace,
    id::StructuralResourceOwnerID, target::HostComputeDevice)
    return KnownStructuralResourceFact(id, target, 0,
        structural_array_bytes(workspace.dispositions, target))
end

function structural_resource_fact(::CommandDispositionWorkspace,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    return UnknownStructuralResourceFact(id, target, :owner_not_on_device)
end

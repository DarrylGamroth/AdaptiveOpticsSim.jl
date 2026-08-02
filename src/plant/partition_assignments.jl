#
# Caller-resolved target-local partition assignments
#
# This cold layer records a complete path-group-to-target assignment without
# selecting resources, allocating storage, or preparing executable replicas.
# The assignment binds a private definition identity token supplied by the
# enclosing plant-definition layer; later partition preparation uses that
# token together with the stable topology snapshot to distinguish a stale
# assignment from one resolved for a different topology.
#

import ..Backends: compute_device_availability, compute_device_is_available,
    compute_device_unavailable_reason
import ..Atmospheres: InfiniteMultiLayerAtmosphereDefinition,
    KolmogorovAtmosphereDefinition, MultiLayerAtmosphereDefinition

struct _PartitionAcquisitionTopology
    id::AcquisitionID
    path::OpticalPathID
end

Base.:(==)(left::_PartitionAcquisitionTopology,
    right::_PartitionAcquisitionTopology) =
    left.id == right.id && left.path == right.path
Base.isequal(left::_PartitionAcquisitionTopology,
    right::_PartitionAcquisitionTopology) =
    isequal(left.id, right.id) && isequal(left.path, right.path)

struct _PartitionOpticTopology{E<:Tuple}
    id::ControllableOpticID
    endpoints::E
end

Base.:(==)(left::_PartitionOpticTopology,
    right::_PartitionOpticTopology) =
    left.id == right.id && left.endpoints == right.endpoints
Base.isequal(left::_PartitionOpticTopology,
    right::_PartitionOpticTopology) =
    isequal(left.id, right.id) && isequal(left.endpoints, right.endpoints)

struct _PlantPartitionTopologySnapshot{
    P<:Tuple,A<:Tuple,O<:Tuple,S<:Tuple,L<:Tuple,
}
    paths::P
    acquisitions::A
    optics::O
    sampled_aberrations::S
    atmosphere_layers::L
end

Base.:(==)(left::_PlantPartitionTopologySnapshot,
    right::_PlantPartitionTopologySnapshot) =
    left.paths == right.paths && left.acquisitions == right.acquisitions &&
    left.optics == right.optics &&
    left.sampled_aberrations == right.sampled_aberrations &&
    left.atmosphere_layers == right.atmosphere_layers
Base.isequal(left::_PlantPartitionTopologySnapshot,
    right::_PlantPartitionTopologySnapshot) =
    isequal(left.paths, right.paths) &&
    isequal(left.acquisitions, right.acquisitions) &&
    isequal(left.optics, right.optics) &&
    isequal(left.sampled_aberrations, right.sampled_aberrations) &&
    isequal(left.atmosphere_layers, right.atmosphere_layers)

struct _ResolvedPartitionPathTarget
    id::OpticalPathID
    target_ordinal::UInt8
end

struct ResolvedPlantPartitionAssignment{D<:PlantDefinition,I,T<:Tuple,P<:Tuple,S}
    definition::D
    definition_identity::I
    topology::S
    targets::T
    atmosphere_authority_target_ordinal::UInt8
    paths::P
end

struct _PartitionAssignmentInput
    id::OpticalPathID
    target::AbstractComputeDevice
end

@inline _partition_assignment_error(reason::Symbol, message::AbstractString) =
    throw(PlantPreparationError(:partition_assignment, reason,
        String(message)))

@inline _partition_path_id(id::OpticalPathID) = id
@inline _partition_path_id(name::Symbol) = OpticalPathID(name)

function _partition_path_id(value)
    _partition_assignment_error(:invalid_path_id,
        "partition path identity must be an OpticalPathID or Symbol; got " *
        "$(typeof(value))")
end

@inline _require_partition_target(target::HostComputeDevice) = target
@inline _require_partition_target(target::AcceleratorComputeDevice) = target

function _require_partition_target(target::AbstractComputeDevice)
    _partition_assignment_error(:unsupported_target,
        "partition target $(typeof(target)) is not a supported host or " *
        "accelerator compute device")
end

function _require_partition_target(value)
    _partition_assignment_error(:invalid_target,
        "partition target must be an AbstractComputeDevice; got " *
        "$(typeof(value))")
end

function _partition_assignment_input(pair::Pair)
    return _PartitionAssignmentInput(_partition_path_id(first(pair)),
        _require_partition_target(last(pair)))
end

function _partition_assignment_input(value)
    _partition_assignment_error(:invalid_path_assignment,
        "each path assignment must be a Pair of path identity and exact " *
        "compute device; got $(typeof(value))")
end

@inline _partition_name_key(id::OpticalPathID) = String(id.name)
@inline _partition_name_key(id::AcquisitionID) = String(id.name)
@inline _partition_name_key(id::ControllableOpticID) = String(id.name)
@inline _partition_name_key(id::SampledAberrationID) = String(id.name)
@inline _partition_name_key(id::CommandEndpointID) = String(id.name)
@inline _partition_name_key(id::AtmosphereLayerID) = String(id.name)

function _canonical_partition_path_ids(definition::PlantDefinition)
    ids = OpticalPathID[path_id(path) for path in path_definitions(definition)]
    sort!(ids; by=_partition_name_key)
    return Tuple(ids)
end

function _canonical_partition_acquisitions(definition::PlantDefinition)
    acquisitions = _PartitionAcquisitionTopology[
        _PartitionAcquisitionTopology(acquisition_id(acquisition),
            acquisition_path_id(acquisition)) for acquisition in
        acquisition_definitions(definition)
    ]
    sort!(acquisitions; by=entry -> _partition_name_key(entry.id))
    return Tuple(acquisitions)
end

function _canonical_partition_optics(definition::PlantDefinition)
    optics = _PartitionOpticTopology[]
    for optic in controllable_optic_definitions(definition)
        endpoints = collect(command_endpoint_ids(optic))
        sort!(endpoints; by=_partition_name_key)
        push!(optics, _PartitionOpticTopology(controllable_optic_id(optic),
            Tuple(endpoints)))
    end
    sort!(optics; by=entry -> _partition_name_key(entry.id))
    return Tuple(optics)
end

function _canonical_partition_sampled_aberrations(definition::PlantDefinition)
    ids = SampledAberrationID[
        sampled_aberration_id(aberration) for aberration in
        sampled_aberration_definitions(definition)
    ]
    sort!(ids; by=_partition_name_key)
    return Tuple(ids)
end

@inline _canonical_partition_atmosphere_layers(
    ::KolmogorovAtmosphereDefinition) = ()

function _canonical_partition_atmosphere_layers(
    definition::AbstractTimedAtmosphereDefinition,
)
    _partition_assignment_error(:unsupported_atmosphere_definition,
        "timed-atmosphere definition $(typeof(definition)) does not expose " *
        "canonical stable layer identities for partition preparation")
end

function _canonical_partition_atmosphere_layers(
    definition::Union{MultiLayerAtmosphereDefinition,
        InfiniteMultiLayerAtmosphereDefinition},
)
    ids = AtmosphereLayerID[layer.id for layer in definition.layers]
    sort!(ids; by=_partition_name_key)
    return Tuple(ids)
end

function _plant_partition_topology_snapshot(definition::PlantDefinition)
    return _PlantPartitionTopologySnapshot(
        _canonical_partition_path_ids(definition),
        _canonical_partition_acquisitions(definition),
        _canonical_partition_optics(definition),
        _canonical_partition_sampled_aberrations(definition),
        _canonical_partition_atmosphere_layers(
            atmosphere_definition(definition)),
    )
end

function _require_complete_partition_path_assignments(
    inputs::AbstractVector{_PartitionAssignmentInput},
    topology::_PlantPartitionTopologySnapshot,
)
    expected = topology.paths
    seen = Set{OpticalPathID}()
    for input in inputs
        input.id in expected || _partition_assignment_error(:unknown_path,
            "partition assignment references unknown optical path $(input.id)")
        input.id in seen && _partition_assignment_error(:duplicate_path,
            "partition assignment names optical path $(input.id) more than once")
        push!(seen, input.id)
    end
    for id in expected
        id in seen || _partition_assignment_error(:missing_path,
            "partition assignment omits optical path $id")
    end
    return nothing
end

function _register_partition_target(host_target,
    accelerator_target, target::HostComputeDevice)
    return target, accelerator_target
end

function _register_partition_target(host_target,
    accelerator_target, target::AcceleratorComputeDevice)
    isnothing(accelerator_target) && return host_target, target
    accelerator_target == target && return host_target, accelerator_target
    _partition_assignment_error(:multiple_accelerators,
        "Gate 9A accepts at most one exact accelerator target; got " *
        "$accelerator_target and $target")
end

function _canonical_partition_targets(authority::AbstractComputeDevice,
    inputs::AbstractVector{_PartitionAssignmentInput})
    host_target = nothing
    accelerator_target = nothing
    host_target, accelerator_target = _register_partition_target(
        host_target, accelerator_target, authority)
    for input in inputs
        host_target, accelerator_target = _register_partition_target(
            host_target, accelerator_target, input.target)
    end
    if isnothing(host_target)
        return (accelerator_target,)
    elseif isnothing(accelerator_target)
        return (host_target,)
    end
    return (host_target, accelerator_target)
end

function _require_partition_target_available(target::AbstractComputeDevice)
    availability = compute_device_availability(target)
    compute_device_is_available(availability) && return target
    reason = compute_device_unavailable_reason(availability)
    _partition_assignment_error(:unavailable_target,
        "exact partition target $target is unavailable" *
        (isnothing(reason) ? "" : " ($reason)"))
end

function _require_partition_targets_available(targets::Tuple)
    for target in targets
        _require_partition_target_available(target)
    end
    return targets
end

function _partition_target_ordinal(targets::Tuple,
    target::AbstractComputeDevice)
    for index in eachindex(targets)
        targets[index] == target && return UInt8(index)
    end
    _partition_assignment_error(:unknown_target,
        "partition target $target is not part of this resolved assignment")
end

@inline function _partition_target_at(targets::Tuple, ordinal::UInt8)
    index = Int(ordinal)
    1 <= index <= length(targets) || _partition_assignment_error(
        :invalid_target_ordinal,
        "partition target ordinal $ordinal is outside the resolved targets")
    return @inbounds targets[index]
end

function _assigned_partition_path_target(
    inputs::AbstractVector{_PartitionAssignmentInput},
    id::OpticalPathID,
)
    for input in inputs
        input.id == id && return input.target
    end
    _partition_assignment_error(:missing_path,
        "partition assignment omits optical path $id")
end

function _canonical_resolved_path_targets(
    inputs::AbstractVector{_PartitionAssignmentInput},
    topology::_PlantPartitionTopologySnapshot,
    targets::Tuple,
)
    return Tuple(_ResolvedPartitionPathTarget(id,
        _partition_target_ordinal(targets,
            _assigned_partition_path_target(inputs, id))) for id in
        topology.paths)
end

function _resolve_plant_partition_assignment(definition::PlantDefinition,
    authority::AbstractComputeDevice, assignments)
    resolved_authority = _require_partition_target(authority)
    inputs = _PartitionAssignmentInput[]
    for assignment in assignments
        push!(inputs, _partition_assignment_input(assignment))
    end
    topology = _plant_partition_topology_snapshot(definition)
    _require_complete_partition_path_assignments(inputs, topology)
    targets = _canonical_partition_targets(resolved_authority, inputs)
    _require_partition_targets_available(targets)
    paths = _canonical_resolved_path_targets(inputs, topology, targets)
    authority_ordinal = _partition_target_ordinal(targets,
        resolved_authority)
    return ResolvedPlantPartitionAssignment(definition,
        _plant_definition_identity(definition), topology, targets,
        authority_ordinal, paths)
end

function resolve_plant_partition_assignment(definition::PlantDefinition,
    authority::AbstractComputeDevice, assignments::Pair...)
    return _resolve_plant_partition_assignment(definition, authority,
        assignments)
end

function resolve_plant_partition_assignment(definition::PlantDefinition,
    authority::AbstractComputeDevice, assignments::Union{Tuple,AbstractVector})
    return _resolve_plant_partition_assignment(definition, authority,
        assignments)
end

function resolve_plant_partition_assignment(definition::PlantDefinition,
    authority::AbstractComputeDevice, assignments)
    _partition_assignment_error(:invalid_assignment_container,
        "path assignments must be Pairs, a Tuple of Pairs, or an " *
        "AbstractVector of Pairs; got $(typeof(assignments))")
end

@inline partition_targets(assignment::ResolvedPlantPartitionAssignment) =
    getfield(assignment, :targets)

@inline atmosphere_authority_target(
    assignment::ResolvedPlantPartitionAssignment) = _partition_target_at(
    partition_targets(assignment),
    getfield(assignment, :atmosphere_authority_target_ordinal))

function partition_target(assignment::ResolvedPlantPartitionAssignment,
    id::OpticalPathID)
    for entry in getfield(assignment, :paths)
        entry.id == id && return _partition_target_at(
            partition_targets(assignment), entry.target_ordinal)
    end
    _partition_assignment_error(:unknown_path,
        "resolved partition assignment has no optical path $id")
end

@inline partition_target(assignment::ResolvedPlantPartitionAssignment,
    id::Symbol) = partition_target(assignment, OpticalPathID(id))

function partition_target(assignment::ResolvedPlantPartitionAssignment, id)
    _partition_assignment_error(:invalid_path_id,
        "partition path identity must be an OpticalPathID or Symbol; got " *
        "$(typeof(id))")
end

function assigned_path_ids(assignment::ResolvedPlantPartitionAssignment,
    target::AbstractComputeDevice)
    ordinal = _partition_target_ordinal(partition_targets(assignment), target)
    return Tuple(entry.id for entry in getfield(assignment, :paths) if
        entry.target_ordinal == ordinal)
end

function assigned_path_ids(assignment::ResolvedPlantPartitionAssignment,
    target)
    _require_partition_target(target)
    _partition_assignment_error(:unknown_target,
        "partition target $target is not part of this resolved assignment")
end

function _require_current_partition_assignment_definition(
    assignment::ResolvedPlantPartitionAssignment,
    definition::PlantDefinition,
)
    topology = _plant_partition_topology_snapshot(definition)
    topology == getfield(assignment, :topology) ||
        _partition_assignment_error(:foreign_assignment,
            "resolved partition assignment targets a different plant topology")
    _plant_definition_identity(definition) ===
        getfield(assignment, :definition_identity) ||
        _partition_assignment_error(:stale_assignment,
            "resolved partition assignment was created for an earlier " *
            "plant-definition identity")
    return assignment
end

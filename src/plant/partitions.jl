#
# Preparation-only target-local plant partitions
#
# One PreparedAtmosphereAuthority owns the sole timed atmosphere, mutable
# timeline, and atmosphere RNG state. Target partitions retain only an
# immutable binding to that authority plus co-located path/acquisition
# resources. This layer may make declared static-data copies while preparing an
# exact target, but it does not publish runtime atmosphere products, hand off
# arrays between running partitions, or expose a mixed-domain executor.
#

mutable struct _AtmosphereAuthorityBindingToken end
const _ATMOSPHERE_AUTHORITY_BINDING_TOKEN =
    _AtmosphereAuthorityBindingToken()

"""Immutable binding to the one prepared atmosphere/epoch authority."""
struct AtmosphereAuthorityBinding{D<:AbstractComputeDevice,I}
    target::D
    identity::I

    function AtmosphereAuthorityBinding(
        token::_AtmosphereAuthorityBindingToken,
        target::D,
        identity::I,
    ) where {D<:AbstractComputeDevice,I}
        token === _ATMOSPHERE_AUTHORITY_BINDING_TOKEN || throw(
            ArgumentError("invalid internal atmosphere-authority token"))
        return new{D,I}(target, identity)
    end
end

mutable struct _PreparedAtmosphereAuthorityToken end
const _PREPARED_ATMOSPHERE_AUTHORITY_TOKEN =
    _PreparedAtmosphereAuthorityToken()

"""The sole prepared owner of timed-atmosphere evolution for a partition set."""
struct PreparedAtmosphereAuthority{D,X,T,A,R,B}
    definition_identity::_PlantDefinitionIdentity
    target::D
    context::X
    telescope::T
    atmosphere::A
    rngs::R
    binding::B

    function PreparedAtmosphereAuthority(
        token::_PreparedAtmosphereAuthorityToken,
        definition_identity::_PlantDefinitionIdentity,
        target::D,
        context::X,
        telescope::T,
        atmosphere::A,
        rngs::R,
        binding::B,
    ) where {D,X,T,A,R,B}
        token === _PREPARED_ATMOSPHERE_AUTHORITY_TOKEN || throw(
            ArgumentError("invalid internal prepared-atmosphere token"))
        return new{D,X,T,A,R,B}(
            definition_identity,
            target,
            context,
            telescope,
            atmosphere,
            rngs,
            binding,
        )
    end
end

struct _PreparedTargetPartitionRNGs
    paths::Memory{PreparedOwnerRNGs}
    acquisitions::Memory{PreparedOwnerRNGs}
end

mutable struct _PreparedTargetPartitionToken end
const _PREPARED_TARGET_PARTITION_TOKEN = _PreparedTargetPartitionToken()

"""
One non-executable exact-target partition. It owns target-local static data and
mutable workspaces for complete path/acquisition groups and role-neutral
target-local controllable optics. It owns no atmosphere timeline, atmosphere
RNG, command admission, authority placement, schedule, or transfer state.
"""
struct PreparedTargetPartition{D,X,T,B,S,P,A,R,O}
    target::D
    context::X
    telescope::T
    authority_binding::B
    sampled_aberrations::S
    paths::P
    acquisitions::A
    rngs::R
    controllable_optics::O
    resource_report::StructuralResourceReport
    controllable_optic_ids::_FixedPlantRegistry{ControllableOpticID}
    command_endpoint_ids::_FixedPlantRegistry{CommandEndpointID}

    function PreparedTargetPartition(
        token::_PreparedTargetPartitionToken,
        target::D,
        context::X,
        telescope::T,
        authority_binding::B,
        sampled_aberrations::S,
        paths::P,
        acquisitions::A,
        rngs::R,
        controllable_optics::O,
        resource_report::StructuralResourceReport,
        controllable_optic_ids::_FixedPlantRegistry{ControllableOpticID},
        command_endpoint_ids::_FixedPlantRegistry{CommandEndpointID},
    ) where {D,X,T,B,S,P,A,R,O}
        token === _PREPARED_TARGET_PARTITION_TOKEN || throw(
            ArgumentError("invalid internal target-partition token"))
        return new{D,X,T,B,S,P,A,R,O}(
            target,
            context,
            telescope,
            authority_binding,
            sampled_aberrations,
            paths,
            acquisitions,
            rngs,
            controllable_optics,
            resource_report,
            controllable_optic_ids,
            command_endpoint_ids,
        )
    end
end

mutable struct _PreparedPlantPartitionsToken end
const _PREPARED_PLANT_PARTITIONS_TOKEN = _PreparedPlantPartitionsToken()

"""Prepared target-local partitions for one exact resolved assignment."""
struct PreparedPlantPartitions{D,A,U,C,P,V}
    definition::D
    assignment::A
    authority::U
    command_authority_identity::C
    partitions::P
    run_seed::UInt64
    rng_derivation_version::V

    function PreparedPlantPartitions(
        token::_PreparedPlantPartitionsToken,
        definition::D,
        assignment::A,
        authority::U,
        command_authority_identity::C,
        partitions::P,
        run_seed::UInt64,
        rng_derivation_version::V,
    ) where {D,A,U,C,P,V}
        token === _PREPARED_PLANT_PARTITIONS_TOKEN || throw(
            ArgumentError("invalid internal prepared-partitions token"))
        return new{D,A,U,C,P,V}(
            definition,
            assignment,
            authority,
            command_authority_identity,
            partitions,
            run_seed,
            rng_derivation_version,
        )
    end
end

struct _UnboundTargetPartition{D,X,T,S,P,A,O}
    target::D
    context::X
    telescope::T
    sampled_aberrations::S
    paths::P
    acquisitions::A
    controllable_optics::O
    controllable_optic_ids::_FixedPlantRegistry{ControllableOpticID}
    command_endpoint_ids::_FixedPlantRegistry{CommandEndpointID}
end

@inline atmosphere_authority_target(authority::PreparedAtmosphereAuthority) =
    authority.target
@inline prepared_atmosphere_authority(
    partitions::PreparedPlantPartitions) = partitions.authority
@inline atmosphere_authority_target(partitions::PreparedPlantPartitions) =
    atmosphere_authority_target(prepared_atmosphere_authority(partitions))
@inline atmosphere_authority_identity(
    binding::AtmosphereAuthorityBinding) = binding.identity
@inline atmosphere_authority_identity(
    authority::PreparedAtmosphereAuthority) =
    atmosphere_authority_identity(authority.binding)
@inline atmosphere_authority_identity(partitions::PreparedPlantPartitions) =
    atmosphere_authority_identity(prepared_atmosphere_authority(partitions))
@inline atmosphere_authority_binding(
    authority::PreparedAtmosphereAuthority) = authority.binding
@inline atmosphere_authority_binding(partitions::PreparedPlantPartitions) =
    atmosphere_authority_binding(prepared_atmosphere_authority(partitions))
@inline atmosphere_authority_binding(partition::PreparedTargetPartition) =
    partition.authority_binding
@inline prepared_telescope(authority::PreparedAtmosphereAuthority) =
    authority.telescope
@inline prepared_atmosphere(authority::PreparedAtmosphereAuthority) =
    authority.atmosphere
@inline prepared_telescope(partition::PreparedTargetPartition) =
    partition.telescope
@inline compute_device(authority::PreparedAtmosphereAuthority) =
    authority.target
@inline compute_device(partition::PreparedTargetPartition) = partition.target
@inline plant_definition(partitions::PreparedPlantPartitions) =
    partitions.definition
@inline prepared_partitions(partitions::PreparedPlantPartitions) =
    partitions.partitions
@inline prepared_paths(partition::PreparedTargetPartition) = partition.paths
@inline prepared_acquisitions(partition::PreparedTargetPartition) =
    partition.acquisitions
@inline target_local_controllable_optic_owners(
    partition::PreparedTargetPartition) = partition.controllable_optics
@inline prepared_sampled_aberrations(partition::PreparedTargetPartition) =
    partition.sampled_aberrations
@inline partition_resource_report(partition::PreparedTargetPartition) =
    partition.resource_report
@inline partition_controllable_optic_ids(
    partition::PreparedTargetPartition) = partition.controllable_optic_ids
@inline partition_command_endpoint_ids(partition::PreparedTargetPartition) =
    partition.command_endpoint_ids
@inline command_authority_identity(partitions::PreparedPlantPartitions) =
    partitions.command_authority_identity

function prepared_partition(partitions::PreparedPlantPartitions,
    target::AbstractComputeDevice)
    for partition in partitions.partitions
        partition.target == target && return partition
    end
    throw(PlantPreparationError(:partition, :unknown_target,
        "prepared plant has no target-local partition for $target"))
end

function prepared_partition(partitions::PreparedPlantPartitions, target)
    throw(PlantPreparationError(:partition, :invalid_target,
        "prepared partition target must be an AbstractComputeDevice; got " *
        "$(typeof(target))"))
end

@inline partition_resource_report(partitions::PreparedPlantPartitions,
    target::AbstractComputeDevice) =
    partition_resource_report(prepared_partition(partitions, target))

function _prepare_atmosphere_authority_resources(
    definition::PlantDefinition,
    target::AbstractComputeDevice,
)
    context = _prepare_device_execution_context(target)
    _prepared_device_execution_compute_device(context) == target || throw(
        PlantPreparationError(:partition, :execution_context_target,
            "atmosphere-authority context does not match its exact target"))
    return _with_completed_prepared_device_execution_context(context) do
        telescope = validate_telescope_target(
            prepare_telescope(telescope_definition(definition), target),
            target,
        )
        atmosphere = validate_timed_atmosphere_target(
            prepare_timed_atmosphere(
                atmosphere_definition(definition), telescope, target),
            target,
        )
        return context, telescope, atmosphere
    end
end

function _partition_telescope_resources(
    target::AbstractComputeDevice,
    authority_target::AbstractComputeDevice,
    authority_context,
    authority_telescope,
    definition::PlantDefinition,
)
    if target == authority_target
        return authority_context, authority_telescope
    end
    context = _prepare_device_execution_context(target)
    _prepared_device_execution_compute_device(context) == target || throw(
        PlantPreparationError(:partition, :execution_context_target,
            "target-local partition context does not match its exact target"))
    telescope = _with_completed_prepared_device_execution_context(context) do
        validate_telescope_target(
            prepare_telescope(telescope_definition(definition), target),
            target,
        )
    end
    return context, telescope
end

function _canonical_partition_path_definitions(
    definition::PlantDefinition,
    assignment::ResolvedPlantPartitionAssignment,
    target::AbstractComputeDevice,
)
    ids = assigned_path_ids(assignment, target)
    definitions = OpticalPathDefinition[
        path_definition(definition, id) for id in ids]
    return _partition_registry(definitions, OpticalPathDefinition)
end

function _prepare_target_local_paths(
    definitions::AbstractVector, telescope, context)
    paths = Memory{PreparedTargetLocalPathResources}(
        undef, length(definitions))
    @inbounds for index in eachindex(definitions)
        paths[index] = prepare_target_local_path_resources(
            definitions[index], telescope, context)
    end
    return paths
end

function _target_local_path_for_acquisition(
    definition::AcquisitionDefinition,
    paths::AbstractVector,
)
    id = acquisition_path_id(definition)
    for path in paths
        path_id(path.definition) == id && return path
    end
    throw(PlantPreparationError(:partition, :split_group,
        "target-local acquisition $(acquisition_id(definition)) is not " *
        "co-located with optical path $id"))
end

function _canonical_partition_acquisition_definitions(
    definition::PlantDefinition,
    path_ids::AbstractVector,
)
    values = AcquisitionDefinition[]
    for acquisition in acquisition_definitions(definition)
        acquisition_path_id(acquisition) in path_ids || continue
        push!(values, acquisition)
    end
    sort!(values; by=acquisition ->
        String(acquisition_id(acquisition).name))
    return _partition_registry(values, AcquisitionDefinition)
end

function _prepare_target_local_acquisitions(
    definitions::AbstractVector, paths)
    acquisitions = Memory{PreparedTargetLocalAcquisitionResources}(
        undef, length(definitions))
    @inbounds for index in eachindex(definitions)
        definition = definitions[index]
        path = _target_local_path_for_acquisition(definition, paths)
        acquisitions[index] =
            prepare_target_local_acquisition_resources(definition, path)
    end
    return acquisitions
end

function _partition_visible(
    visibility::AbstractPathVisibility,
    path_ids::AbstractVector,
)
    for id in path_ids
        _visible_on_path(visibility, id) && return true
    end
    return false
end

function _partition_sampled_aberration_definitions(
    definition::PlantDefinition,
    path_ids::AbstractVector,
)
    values = SampledAberrationDefinition[]
    for aberration in _canonical_sampled_aberration_definitions(definition)
        _partition_visible(sampled_aberration_visibility(aberration),
            path_ids) || continue
        push!(values, aberration)
    end
    return _partition_registry(values, SampledAberrationDefinition)
end

function _prepare_partition_sampled_aberrations(
    definitions::AbstractVector,
    target::AbstractComputeDevice,
)
    prepared = Memory{PreparedSampledAberration}(undef, length(definitions))
    @inbounds for index in eachindex(definitions)
        prepared[index] = _prepare_sampled_aberration(
            definitions[index], target)
    end
    return prepared
end

function _partition_controllable_optic_definitions(
    definition::PlantDefinition,
    path_ids::AbstractVector,
)
    values = ControllableOpticDefinition[]
    for optic in controllable_optic_definitions(definition)
        _partition_visible(controllable_optic_visibility(optic), path_ids) ||
            continue
        push!(values, optic)
    end
    sort!(values; by=optic -> String(controllable_optic_id(optic).name))
    return _partition_registry(values, ControllableOpticDefinition)
end

function _partition_controllable_optic_identities(optics::AbstractVector)
    optic_ids = ControllableOpticID[
        controllable_optic_id(optic) for optic in optics]
    endpoint_values = CommandEndpointID[]
    for optic in optics
        append!(endpoint_values, command_endpoint_ids(optic))
    end
    sort!(endpoint_values; by=id -> String(id.name))
    return (
        _partition_registry(optic_ids, ControllableOpticID),
        _partition_registry(endpoint_values, CommandEndpointID),
    )
end

function _prepare_unbound_target_partition(
    definition::PlantDefinition,
    assignment::ResolvedPlantPartitionAssignment,
    target::AbstractComputeDevice,
    authority_target::AbstractComputeDevice,
    authority_context,
    authority_telescope,
    command_configurations,
    command_authority_identity::CommandAuthorityIdentity,
)
    context, telescope = _partition_telescope_resources(
        target,
        authority_target,
        authority_context,
        authority_telescope,
        definition,
    )
    path_definitions = _canonical_partition_path_definitions(
        definition, assignment, target)
    path_ids = _partition_registry(
        OpticalPathID[path_id(path) for path in path_definitions],
        OpticalPathID,
    )
    acquisition_definitions =
        _canonical_partition_acquisition_definitions(definition, path_ids)
    sampled_definitions = _partition_sampled_aberration_definitions(
        definition, path_ids)
    optic_definitions = _partition_controllable_optic_definitions(
        definition, path_ids)
    optic_ids, endpoint_ids =
        _partition_controllable_optic_identities(optic_definitions)
    return _with_completed_prepared_device_execution_context(context) do
        sampled_aberrations = _prepare_partition_sampled_aberrations(
            sampled_definitions, target)
        paths = _prepare_target_local_paths(
            path_definitions, telescope, context)
        acquisitions = _prepare_target_local_acquisitions(
            acquisition_definitions, paths)
        controllable_optics =
            _prepare_target_local_controllable_optic_owners(
                optic_definitions,
                telescope,
                atmosphere_definition(definition),
                target,
                command_configurations,
                command_authority_identity,
            )
        return _UnboundTargetPartition(
            target,
            context,
            telescope,
            sampled_aberrations,
            paths,
            acquisitions,
            controllable_optics,
            optic_ids,
            endpoint_ids,
        )
    end
end

function _target_local_path_rng_owner_roles(path)
    additional = _require_additional_rng_roles(
        additional_path_rng_owner_roles(path.execution), :path)
    return (:provider, additional...)
end

function _target_local_path_rng_owner_bindings(paths::AbstractVector)
    groups = Memory{Tuple}(undef, length(paths))
    @inbounds for index in eachindex(paths)
        path = paths[index]
        groups[index] = _rng_owner_bindings(path, :path,
            path_id(path.definition).name,
            _target_local_path_rng_owner_roles(path))
    end
    return groups
end

function _target_local_acquisition_rng_owner_bindings(
    acquisitions::AbstractVector)
    groups = Memory{Tuple}(undef, length(acquisitions))
    @inbounds for index in eachindex(acquisitions)
        acquisition = acquisitions[index]
        groups[index] = _rng_owner_bindings(acquisition, :acquisition,
            acquisition_id(acquisition.definition).name,
            _acquisition_rng_owner_roles(
                acquisition.provider.implementation))
    end
    return groups
end

function _partition_rng_binding_groups(partition::_UnboundTargetPartition)
    return (
        paths=_target_local_path_rng_owner_bindings(partition.paths),
        acquisitions=_target_local_acquisition_rng_owner_bindings(
            partition.acquisitions),
    )
end

function _append_partition_rng_bindings!(destination, groups)
    _append_rng_binding_groups!(destination, groups.paths)
    _append_rng_binding_groups!(destination, groups.acquisitions)
    return destination
end

function _prepare_partition_rngs(groups, run_seed::UInt64,
    version::RNGDerivationVersion)
    return _PreparedTargetPartitionRNGs(
        _prepare_owner_rng_groups(groups.paths, run_seed, version),
        _prepare_owner_rng_groups(groups.acquisitions, run_seed, version),
    )
end

function _partition_resource_report(
    partition::_UnboundTargetPartition,
    authority_target::AbstractComputeDevice,
    atmosphere,
)
    target = partition.target
    facts = AbstractStructuralResourceFact[]
    push!(facts, structural_resource_fact(
        partition.telescope,
        StructuralResourceOwnerID(:telescope, :primary),
        target,
    ))
    if target == authority_target
        push!(facts, structural_resource_fact(
            atmosphere,
            StructuralResourceOwnerID(:atmosphere, :primary),
            target,
        ))
    end
    for aberration in partition.sampled_aberrations
        id = StructuralResourceOwnerID(
            :sampled_aberration,
            sampled_aberration_id(aberration).name,
        )
        push!(facts, structural_resource_fact(aberration, id, target))
    end
    for path in partition.paths
        id = StructuralResourceOwnerID(:path, path_id(path.definition).name)
        push!(facts, structural_resource_fact(path, id, target))
    end
    for acquisition in partition.acquisitions
        id = StructuralResourceOwnerID(
            :acquisition,
            acquisition_id(acquisition.definition).name,
        )
        push!(facts, structural_resource_fact(acquisition, id, target))
    end
    for owner in partition.controllable_optics
        optic_id = controllable_optic_id(owner)
        id = StructuralResourceOwnerID(
            :controllable_optic_replica, optic_id.name)
        push!(facts, _combine_structural_owner_facts(id, target, (
            structural_resource_fact(owner.prepared, id, target),
            structural_resource_fact(owner.state, id, target),
            structural_resource_fact(owner.workspace, id, target),
        )))
        @inbounds for index in eachindex(owner.prepared.endpoints)
            endpoint = owner.prepared.endpoints[index]
            endpoint_state = owner.state.endpoints[index]
            endpoint_id = StructuralResourceOwnerID(
                :effective_command_replica,
                command_endpoint_id(endpoint).name,
            )
            push!(facts, structural_resource_fact(
                endpoint_state, endpoint_id, target))
        end
    end
    return aggregate_structural_resource_facts(facts, target)
end

function _bind_target_partition(
    partition::_UnboundTargetPartition,
    authority_binding::AtmosphereAuthorityBinding,
    rngs::_PreparedTargetPartitionRNGs,
    report::StructuralResourceReport,
)
    compute_device(report) == partition.target || throw(
        PlantPreparationError(:partition, :resource_report_target,
            "target-local resource report does not match its partition"))
    return PreparedTargetPartition(
        _PREPARED_TARGET_PARTITION_TOKEN,
        partition.target,
        partition.context,
        partition.telescope,
        authority_binding,
        partition.sampled_aberrations,
        partition.paths,
        partition.acquisitions,
        rngs,
        partition.controllable_optics,
        report,
        partition.controllable_optic_ids,
        partition.command_endpoint_ids,
    )
end

function _partition_rng_replay_metadata(
    prepared::PreparedPlantPartitions,
)
    owners = NamedTuple[]
    for stream in prepared.authority.rngs.streams
        push!(owners, _rng_stream_metadata(stream))
    end
    for partition in prepared.partitions
        for group in partition.rngs.paths
            for stream in _rng_group_streams(group)
                push!(owners, _rng_stream_metadata(stream))
            end
        end
        for group in partition.rngs.acquisitions
            for stream in _rng_group_streams(group)
                push!(owners, _rng_stream_metadata(stream))
            end
        end
    end
    sort!(owners; by=_rng_metadata_order)
    owned = Memory{NamedTuple}(undef, length(owners))
    copyto!(owned, owners)
    return (
        run_seed=prepared.run_seed,
        derivation_version=prepared.rng_derivation_version.value,
        derivation_algorithm=_RNG_DERIVATION_ALGORITHM,
        stream_algorithm=_RNG_STREAM_ALGORITHM,
        owners=owned,
    )
end

@inline rng_replay_metadata(prepared::PreparedPlantPartitions) =
    _partition_rng_replay_metadata(prepared)

"""
    prepare_plant_partitions(definition, assignment; run_seed,
        rng_derivation_version, command_endpoints=())

Prepare one cold plant declaration into complete exact-target resource
partitions. Preparation validates the resolved assignment before allocating,
constructs exactly one timed atmosphere/timeline/RNG authority, and may copy
declared static data and configured initial effective commands into its
target-local owners. It creates one placement-neutral command-authority
identity but no command authority, publication source, or admission state. It
does not publish runtime products or move arrays between prepared partitions.
"""
function prepare_plant_partitions(
    definition::PlantDefinition,
    assignment::ResolvedPlantPartitionAssignment;
    run_seed,
    rng_derivation_version=_DEFAULT_RNG_DERIVATION_VERSION,
    command_endpoints=(),
)
    _require_current_partition_assignment_definition(assignment, definition)
    _require_partition_targets_available(partition_targets(assignment))
    seed = _prepare_run_seed(run_seed)
    version = _prepare_rng_derivation_version(rng_derivation_version)
    command_configurations = _partition_command_endpoint_configurations(
        definition, command_endpoints)
    command_identity = CommandAuthorityIdentity(
        _COMMAND_AUTHORITY_IDENTITY_TOKEN)
    authority_target = atmosphere_authority_target(assignment)
    authority_context, authority_telescope, atmosphere =
        _prepare_atmosphere_authority_resources(definition, authority_target)
    targets = partition_targets(assignment)
    unbound = Memory{_UnboundTargetPartition}(undef, length(targets))
    @inbounds for index in eachindex(targets)
        unbound[index] = _prepare_unbound_target_partition(
            definition,
            assignment,
            targets[index],
            authority_target,
            authority_context,
            authority_telescope,
            command_configurations,
            command_identity,
        )
    end
    atmosphere_bindings = _atmosphere_rng_owner_bindings(atmosphere)
    partition_binding_groups = Memory{NamedTuple}(undef, length(unbound))
    @inbounds for index in eachindex(unbound)
        partition_binding_groups[index] =
            _partition_rng_binding_groups(unbound[index])
    end
    all_bindings = RNGOwnerBinding[]
    sizehint!(all_bindings, length(atmosphere_bindings) +
        sum(groups -> length(groups.paths) + length(groups.acquisitions),
            partition_binding_groups; init=0))
    append!(all_bindings, atmosphere_bindings)
    for groups in partition_binding_groups
        _append_partition_rng_bindings!(all_bindings, groups)
    end
    _require_unique_rng_owner_identities(all_bindings)
    _require_unique_rng_stream_seeds(all_bindings, seed, version)
    atmosphere_rngs = _prepare_atmosphere_rngs(
        atmosphere, atmosphere_bindings, seed, version)
    binding = AtmosphereAuthorityBinding(
        _ATMOSPHERE_AUTHORITY_BINDING_TOKEN,
        authority_target, atmosphere_identity(atmosphere))
    authority = PreparedAtmosphereAuthority(
        _PREPARED_ATMOSPHERE_AUTHORITY_TOKEN,
        _plant_definition_identity(definition),
        authority_target,
        authority_context,
        authority_telescope,
        atmosphere,
        atmosphere_rngs,
        binding,
    )
    partitions = Memory{PreparedTargetPartition}(undef, length(unbound))
    @inbounds for index in eachindex(unbound)
        resources = unbound[index]
        groups = partition_binding_groups[index]
        rngs = _prepare_partition_rngs(groups, seed, version)
        report = _partition_resource_report(
            resources, authority_target, atmosphere)
        partitions[index] =
            _bind_target_partition(resources, binding, rngs, report)
    end
    return PreparedPlantPartitions(
        _PREPARED_PLANT_PARTITIONS_TOKEN,
        definition,
        assignment,
        authority,
        command_identity,
        partitions,
        seed,
        version,
    )
end

@inline function prepare_plant_partitions(
    assignment::ResolvedPlantPartitionAssignment;
    kwargs...,
)
    return prepare_plant_partitions(
        plant_definition(assignment), assignment; kwargs...)
end

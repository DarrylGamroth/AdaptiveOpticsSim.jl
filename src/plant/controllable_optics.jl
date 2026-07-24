#
# Prepared controllable optics and endpoint execution configuration
#
# A cold ControllableOpticDefinition describes physical identity, optical
# placement/path visibility, and semantic command schemas. This layer
# separately binds run capacities, initial/safe values, backend storage,
# model-specific physical preparation, and canonical path/group lookup.
# Mutable optic, endpoint, and workspace owners are constructed by the event
# loop.
#

"""
    CommandEndpointConfiguration(endpoint, initial_command;
        capacity, sequence_window=capacity, safe_command=nothing,
        backend=CPUBackend())

Run-configuration input for one declared command endpoint. It supplies bounded
calendar/history capacity, copied initial and optional safe effective commands,
and payload-storage backend without adding cadence, transport, atomicity, or
optical grouping. Stable endpoint ordinals are derived canonically during
plant preparation rather than supplied by declaration order.
"""
struct CommandEndpointConfiguration{I,S,B<:AbstractArrayBackend}
    endpoint::CommandEndpointID
    capacity::Int
    sequence_window::Int
    initial_command::I
    safe_command::S
    backend::B
end

function CommandEndpointConfiguration(endpoint, initial_command;
    capacity,
    sequence_window=capacity,
    safe_command=nothing,
    backend::AbstractArrayBackend=CPUBackend())
    return CommandEndpointConfiguration(
        _as_command_endpoint_id(endpoint),
        _checked_command_endpoint_capacity(capacity),
        _checked_command_sequence_window(sequence_window),
        initial_command,
        safe_command,
        backend,
    )
end

@inline command_endpoint_id(configuration::CommandEndpointConfiguration) =
    configuration.endpoint
@inline command_endpoint_capacity(
    configuration::CommandEndpointConfiguration) = configuration.capacity
@inline command_sequence_window(
    configuration::CommandEndpointConfiguration) =
    configuration.sequence_window
@inline initial_effective_command(
    configuration::CommandEndpointConfiguration) =
    configuration.initial_command
@inline safe_effective_command(
    configuration::CommandEndpointConfiguration) =
    configuration.safe_command
@inline backend(configuration::CommandEndpointConfiguration) =
    configuration.backend

@inline _require_command_endpoint_configuration(
    configuration::CommandEndpointConfiguration) = configuration

function _require_command_endpoint_configuration(configuration)
    throw(PlantPreparationError(:command_endpoint, :invalid_configuration,
        "command endpoint configurations must contain " *
        "CommandEndpointConfiguration values; got $(typeof(configuration))"))
end

function _command_endpoint_configuration_tuple(configurations::Tuple)
    foreach(_require_command_endpoint_configuration, configurations)
    return configurations
end

function _command_endpoint_configuration_tuple(configurations::NamedTuple)
    values_tuple = values(configurations)
    foreach(_require_command_endpoint_configuration, values_tuple)
    for (name, configuration) in pairs(configurations)
        name == command_endpoint_id(configuration).name || throw(
            PlantPreparationError(:command_endpoint, :identity_mismatch,
                "named command-endpoint configuration key $(repr(name)) " *
                "does not match $(command_endpoint_id(configuration))"))
    end
    return values_tuple
end

function _command_endpoint_configuration_tuple(
    configurations::AbstractVector)
    values_tuple = Tuple(configurations)
    foreach(_require_command_endpoint_configuration, values_tuple)
    return values_tuple
end

function _command_endpoint_configuration_tuple(configurations)
    throw(PlantPreparationError(:command_endpoint, :invalid_configuration,
        "command endpoint configurations must be a Tuple, NamedTuple, or " *
        "AbstractVector; got $(typeof(configurations))"))
end

function _sorted_command_endpoint_configurations(configurations)
    values = collect(_command_endpoint_configuration_tuple(configurations))
    sort!(values; by=configuration ->
        String(command_endpoint_id(configuration).name))
    @inbounds for index in 2:length(values)
        command_endpoint_id(values[index - 1]) ==
            command_endpoint_id(values[index]) && throw(
            PlantPreparationError(:command_endpoint,
                :duplicate_configuration,
                "command endpoint $(command_endpoint_id(values[index])) " *
                "has more than one run configuration"))
    end
    return values
end

struct _PreparedPlantCommandEndpoint{E<:PreparedCommandEndpoint,I,S}
    endpoint::E
    optic_slot::UInt32
    initial_command::I
    safe_command::S
end

@inline command_endpoint_id(binding::_PreparedPlantCommandEndpoint) =
    command_endpoint_id(binding.endpoint)
@inline command_endpoint_capacity(binding::_PreparedPlantCommandEndpoint) =
    command_endpoint_capacity(binding.endpoint)
@inline command_sequence_window(binding::_PreparedPlantCommandEndpoint) =
    command_sequence_window(binding.endpoint)
@inline command_endpoint_ordinal(binding::_PreparedPlantCommandEndpoint) =
    command_endpoint_ordinal(binding.endpoint)
@inline initial_effective_command(binding::_PreparedPlantCommandEndpoint) =
    binding.initial_command
@inline safe_effective_command(binding::_PreparedPlantCommandEndpoint) =
    binding.safe_command

struct _PreparedControllableOpticToken end
const _PREPARED_CONTROLLABLE_OPTIC_TOKEN =
    _PreparedControllableOpticToken()

"""
Model-specific physical preparation for one declared controllable optic.

`implementation` is run-immutable preparation data. `endpoint_slots` names the
independent prepared endpoint owners of this physical device; it does not pack
their commands or synchronize their timing.
"""
struct PreparedControllableOptic{D<:ControllableOpticDefinition,P,S<:Tuple}
    definition::D
    implementation::P
    endpoint_slots::S

    function PreparedControllableOptic(
        ::_PreparedControllableOpticToken,
        definition::D,
        implementation::P,
        endpoint_slots::S,
    ) where {D<:ControllableOpticDefinition,P,S<:Tuple}
        return new{D,P,S}(definition, implementation, endpoint_slots)
    end
end

@inline controllable_optic_implementation(
    optic::PreparedControllableOptic) = optic.implementation
@inline controllable_optic_placement(optic::PreparedControllableOptic) =
    controllable_optic_placement(optic.definition)
@inline controllable_optic_visibility(optic::PreparedControllableOptic) =
    controllable_optic_visibility(optic.definition)

"""
    prepare_controllable_optic(model, definition, telescope, atmosphere)

Extension seam that prepares immutable model-specific physical-optic execution
data. Mutable physical state and scratch are constructed separately through
`prepare_controllable_optic_state` and
`prepare_controllable_optic_workspace`.
"""
function prepare_controllable_optic(model,
    ::ControllableOpticDefinition,
    ::AbstractTelescope,
    ::AbstractAtmosphere)
    throw(PlantPreparationError(:controllable_optic, :unsupported_model,
        "controllable-optic model $(typeof(model)) does not implement " *
        "prepare_controllable_optic"))
end

"""
Construct the single-writer physical state for one prepared optic.

Every array in `initial_commands` is a fresh state-owned copy rather than
prepared-plan or caller storage. Implementations may retain and mutate those
arrays as part of the returned physical state.
"""
function prepare_controllable_optic_state(implementation,
    ::ControllableOpticDefinition, endpoint_ids::Tuple,
    initial_commands::Tuple)
    throw(PlantPreparationError(:controllable_optic, :unsupported_state,
        "prepared controllable-optic implementation " *
        "$(typeof(implementation)) does not implement " *
        "prepare_controllable_optic_state"))
end

"""Construct caller-owned scratch for one prepared optic."""
function prepare_controllable_optic_workspace(implementation)
    throw(PlantPreparationError(:controllable_optic,
        :unsupported_workspace,
        "prepared controllable-optic implementation " *
        "$(typeof(implementation)) does not implement " *
        "prepare_controllable_optic_workspace"))
end

"""
Stage one validated effective endpoint command without changing the visible
physical surface. `timestamp` is the plant time at which the command becomes
physically effective. A successful stage must make
`commit_controllable_optic_command!` an infallible bounded publication step.
"""
function stage_controllable_optic_command!(implementation, state, workspace,
    endpoint::CommandEndpointID, effective_command,
    timestamp::PlantTimestamp)
    throw(PlantCommandError(:physical_application,
        :unsupported_optic_application,
        "prepared controllable-optic implementation " *
        "$(typeof(implementation)) does not implement staged command " *
        "application for $endpoint at $timestamp"))
end

"""
Publish one previously staged physical command at `timestamp`.
Implementations must keep this bounded and nonthrowing after successful
staging so explicit multi-optic transactions cannot become partially visible.
"""
function commit_controllable_optic_command!(implementation, state, workspace,
    endpoint::CommandEndpointID, timestamp::PlantTimestamp)
    throw(PlantCommandError(:physical_application,
        :unsupported_optic_commit,
        "prepared controllable-optic implementation " *
        "$(typeof(implementation)) does not implement command commit for " *
        "$endpoint at $timestamp"))
end

"""
Apply the currently visible physical surface to one already materialized
optical-path input. Composition calls only the path's prepared visible
`PupilSurfaceExecutionRole` bindings, grouped by placement and ordered by
canonical optic identity inside each group. Path-local autonomous devices use
their separately prepared exact coupling. Atmospheric-conjugate transforms
are prepared by a later geometry layer.
"""
function apply_controllable_optic_surface!(input, implementation, state)
    throw(PlantPreparationError(:controllable_optic,
        :unsupported_surface_application,
        "prepared controllable-optic implementation " *
        "$(typeof(implementation)) does not apply to path input " *
        "$(typeof(input))"))
end

function _copy_prepared_effective_command(
    endpoint::PreparedCommandEndpoint{<:PlantCommandSchema{T,0}},
    value, label::AbstractString) where {T}
    return _validate_effective_seed(command_schema(endpoint), value, label)
end

function _copy_prepared_effective_command(
    endpoint::PreparedCommandEndpoint{<:PlantCommandSchema{T,N}},
    value, label::AbstractString) where {T,N}
    validated = _validate_effective_seed(command_schema(endpoint), value,
        label)
    copied = allocate_array(backend(endpoint), T,
        command_dimensions(command_schema(endpoint))...)
    copyto!(copied, validated)
    return copied
end

function _copy_prepared_safe_command(endpoint::PreparedCommandEndpoint,
    ::Nothing)
    policy = command_silence_policy(command_schema(endpoint))
    _require_safe_command_configuration(policy, nothing)
    return nothing
end

function _copy_prepared_safe_command(endpoint::PreparedCommandEndpoint,
    value)
    policy = command_silence_policy(command_schema(endpoint))
    _require_safe_command_configuration(policy, value)
    return _copy_prepared_effective_command(endpoint, value, "safe command")
end

function _command_endpoint_configuration(configurations,
    id::CommandEndpointID)
    @inbounds for configuration in configurations
        command_endpoint_id(configuration) == id && return configuration
    end
    throw(PlantPreparationError(:command_endpoint, :missing_configuration,
        "declared command endpoint $id has no run configuration"))
end

function _declared_command_endpoint_count(definition::PlantDefinition)
    count = 0
    for optic in controllable_optic_definitions(definition)
        count += length(command_schemas(optic))
    end
    return count
end

function _canonical_command_endpoint_declarations(
    definition::PlantDefinition)
    declarations = Tuple{ControllableOpticID,PlantCommandSchema}[]
    for optic in controllable_optic_definitions(definition)
        owner = controllable_optic_id(optic)
        for schema in command_schemas(optic)
            push!(declarations, (owner, schema))
        end
    end
    sort!(declarations; by=declaration ->
        String(command_endpoint_id(declaration[2]).name))
    return declarations
end

function _canonical_controllable_optic_definitions(
    definition::PlantDefinition)
    definitions = collect(controllable_optic_definitions(definition))
    sort!(definitions; by=optic -> String(controllable_optic_id(optic).name))
    return definitions
end

function _controllable_optic_slot(definitions,
    id::ControllableOpticID)
    @inbounds for index in eachindex(definitions)
        controllable_optic_id(definitions[index]) == id && return index
    end
    throw(PlantPreparationError(:controllable_optic, :unknown_id,
        "command endpoint references unknown controllable optic $id"))
end

function _prepare_plant_command_endpoints(definition::PlantDefinition,
    configurations, optic_definitions)
    declared = _canonical_command_endpoint_declarations(definition)
    length(configurations) == length(declared) || throw(
        PlantPreparationError(:command_endpoint,
            :configuration_count,
            "plant declares $(length(declared)) command endpoints but " *
            "$(length(configurations)) run configurations were supplied"))
    endpoints = Any[]
    sizehint!(endpoints, length(declared))
    for (ordinal, (owner, schema)) in enumerate(declared)
        configuration = _command_endpoint_configuration(configurations,
            command_endpoint_id(schema))
        endpoint = prepare_command_endpoint(schema;
            capacity=configuration.capacity,
            sequence_window=configuration.sequence_window,
            ordinal,
            backend=configuration.backend)
        initial = _copy_prepared_effective_command(endpoint,
            configuration.initial_command, "initial effective command")
        safe = _copy_prepared_safe_command(endpoint,
            configuration.safe_command)
        optic_slot = _controllable_optic_slot(optic_definitions, owner)
        push!(endpoints, _PreparedPlantCommandEndpoint(endpoint,
            UInt32(optic_slot), initial, safe))
    end
    return Tuple(endpoints)
end

function _prepared_command_endpoint_slot(endpoints,
    id::CommandEndpointID)
    @inbounds for index in eachindex(endpoints)
        command_endpoint_id(endpoints[index]) == id && return UInt32(index)
    end
    throw(PlantPreparationError(:command_endpoint, :unknown_id,
        "prepared plant has no command endpoint $id"))
end

function _prepare_controllable_optics(definition::PlantDefinition,
    optic_definitions, endpoints)
    telescope = plant_telescope(definition)
    atmosphere = plant_atmosphere(definition)
    optics = Any[]
    sizehint!(optics, length(optic_definitions))
    for optic_definition in optic_definitions
        implementation = prepare_controllable_optic(
            controllable_optic_model(optic_definition), optic_definition,
            telescope, atmosphere)
        implementation === nothing && throw(PlantPreparationError(
            :controllable_optic, :invalid_preparation,
            "prepare_controllable_optic returned nothing for " *
            "$(controllable_optic_id(optic_definition))"))
        ismutabletype(typeof(implementation)) && throw(
            PlantPreparationError(:controllable_optic,
                :mutable_preparation,
                "prepare_controllable_optic must return immutable " *
                "preparation data; mutable physical state is constructed " *
                "separately"))
        slots = map(command_schemas(optic_definition)) do schema
            _prepared_command_endpoint_slot(endpoints,
                command_endpoint_id(schema))
        end
        push!(optics, PreparedControllableOptic(
            _PREPARED_CONTROLLABLE_OPTIC_TOKEN, optic_definition,
            implementation, slots))
    end
    return Tuple(optics)
end

struct _PreparedControllableOpticPathBindingsToken end
const _PREPARED_CONTROLLABLE_OPTIC_PATH_BINDINGS_TOKEN =
    _PreparedControllableOpticPathBindingsToken()

"""
One contiguous co-placed group in a prepared path-to-optic binding table.

`representative_optic_slot` identifies an optic whose immutable placement
describes the group. It is a run-local prepared slot, not a physical identity.
"""
struct PreparedControllableOpticPlaneGroup
    path_slot::UInt32
    representative_optic_slot::UInt32
    first_binding::Int
    binding_count::Int

    function PreparedControllableOpticPlaneGroup(
        ::_PreparedControllableOpticPathBindingsToken,
        path_slot::UInt32,
        representative_optic_slot::UInt32,
        first_binding::Int,
        binding_count::Int,
    )
        return new(path_slot, representative_optic_slot, first_binding,
            binding_count)
    end
end

"""
Canonical bounded path-to-optic bindings prepared from stable identities.

Paths are indexed canonically by `OpticalPathID`. Within each path, visible
optics are grouped by optical placement and ordered by stable
`ControllableOpticID` inside a group. The offsets provide direct bounded
lookup; repeated execution performs no visibility-set lookup or all-optic
scan.
"""
struct PreparedControllableOpticPathBindings
    path_ids::Memory{OpticalPathID}
    path_slots::Memory{UInt32}
    binding_offsets::Memory{Int}
    group_offsets::Memory{Int}
    optic_slots::Memory{UInt32}
    plane_groups::Memory{PreparedControllableOpticPlaneGroup}

    function PreparedControllableOpticPathBindings(
        ::_PreparedControllableOpticPathBindingsToken,
        path_ids::Memory{OpticalPathID},
        path_slots::Memory{UInt32},
        binding_offsets::Memory{Int},
        group_offsets::Memory{Int},
        optic_slots::Memory{UInt32},
        plane_groups::Memory{PreparedControllableOpticPlaneGroup},
    )
        return new(path_ids, path_slots, binding_offsets, group_offsets,
            optic_slots, plane_groups)
    end
end

function _optic_binding_memory(values::Vector{T}) where {T}
    result = Memory{T}(undef, length(values))
    copyto!(result, values)
    return result
end

function _canonical_prepared_path_slots(paths::Tuple)
    slots = collect(eachindex(paths))
    sort!(slots; by=slot ->
        String(path_id(paths[slot].definition).name))
    return slots
end

@inline function _optic_binding_slot_isless(left::UInt32, right::UInt32,
    optics::Tuple)
    left_optic = optics[Int(left)]
    right_optic = optics[Int(right)]
    left_placement = controllable_optic_placement(left_optic)
    right_placement = controllable_optic_placement(right_optic)
    if _same_optic_placement(left_placement, right_placement)
        return String(controllable_optic_id(left_optic.definition).name) <
            String(controllable_optic_id(right_optic.definition).name)
    end
    return _optic_placement_isless(left_placement, right_placement)
end

function _visible_prepared_optic_slots(
    optics::Tuple, path::OpticalPathID)
    slots = UInt32[]
    sizehint!(slots, length(optics))
    @inbounds for slot in eachindex(optics)
        optic = optics[slot]
        _optic_visible_on_path(controllable_optic_visibility(optic), path) ||
            continue
        slot <= typemax(UInt32) || throw(PlantPreparationError(
            :controllable_optic, :capacity,
            "prepared controllable-optic count exceeds UInt32 capacity"))
        push!(slots, UInt32(slot))
    end
    sort!(slots; lt=(left, right) ->
        _optic_binding_slot_isless(left, right, optics))
    return slots
end

function _append_prepared_optic_plane_groups!(
    groups::Vector{PreparedControllableOpticPlaneGroup},
    optic_slots::Vector{UInt32},
    visible_slots::Vector{UInt32},
    path_slot::UInt32,
    optics::Tuple,
)
    isempty(visible_slots) && return nothing
    first_binding = length(optic_slots) + 1
    first_slot = first(visible_slots)
    previous_placement =
        controllable_optic_placement(optics[Int(first_slot)])
    group_first = first_binding
    group_representative = first_slot
    @inbounds for slot in visible_slots
        placement = controllable_optic_placement(optics[Int(slot)])
        if !_same_optic_placement(previous_placement, placement)
            push!(groups, PreparedControllableOpticPlaneGroup(
                _PREPARED_CONTROLLABLE_OPTIC_PATH_BINDINGS_TOKEN,
                path_slot, group_representative, group_first,
                length(optic_slots) - group_first + 1))
            group_first = length(optic_slots) + 1
            group_representative = slot
            previous_placement = placement
        end
        push!(optic_slots, slot)
    end
    push!(groups, PreparedControllableOpticPlaneGroup(
        _PREPARED_CONTROLLABLE_OPTIC_PATH_BINDINGS_TOKEN,
        path_slot, group_representative, group_first,
        length(optic_slots) - group_first + 1))
    return nothing
end

function _prepare_controllable_optic_path_bindings(
    optics::Tuple, paths::Tuple)
    length(paths) <= typemax(UInt32) || throw(PlantPreparationError(
        :path, :capacity, "prepared path count exceeds UInt32 capacity"))
    canonical_path_slots = _canonical_prepared_path_slots(paths)
    path_ids = OpticalPathID[]
    path_slots = UInt32[]
    binding_offsets = Int[1]
    group_offsets = Int[1]
    optic_slots = UInt32[]
    groups = PreparedControllableOpticPlaneGroup[]
    sizehint!(path_ids, length(paths))
    sizehint!(path_slots, length(paths))
    sizehint!(binding_offsets, length(paths) + 1)
    sizehint!(group_offsets, length(paths) + 1)
    sizehint!(optic_slots, length(paths) * length(optics))
    sizehint!(groups, length(paths) * length(optics))
    @inbounds for path_slot_value in canonical_path_slots
        path_slot = UInt32(path_slot_value)
        id = path_id(paths[path_slot_value].definition)
        push!(path_ids, id)
        push!(path_slots, path_slot)
        visible_slots = _visible_prepared_optic_slots(optics, id)
        _append_prepared_optic_plane_groups!(groups, optic_slots,
            visible_slots, path_slot, optics)
        push!(binding_offsets, length(optic_slots) + 1)
        push!(group_offsets, length(groups) + 1)
    end
    return PreparedControllableOpticPathBindings(
        _PREPARED_CONTROLLABLE_OPTIC_PATH_BINDINGS_TOKEN,
        _optic_binding_memory(path_ids),
        _optic_binding_memory(path_slots),
        _optic_binding_memory(binding_offsets),
        _optic_binding_memory(group_offsets),
        _optic_binding_memory(optic_slots),
        _optic_binding_memory(groups),
    )
end

@inline prepared_controllable_optic_path_count(
    bindings::PreparedControllableOpticPathBindings) =
    length(bindings.path_ids)
@inline prepared_controllable_optic_binding_count(
    bindings::PreparedControllableOpticPathBindings) =
    length(bindings.optic_slots)
@inline prepared_controllable_optic_plane_group_count(
    bindings::PreparedControllableOpticPathBindings) =
    length(bindings.plane_groups)
@inline prepared_controllable_optic_path_id(
    bindings::PreparedControllableOpticPathBindings, ordinal::Integer) =
    bindings.path_ids[ordinal]

function _prepared_controllable_optic_path_ordinal(
    bindings::PreparedControllableOpticPathBindings,
    path::OpticalPathID,
)
    @inbounds for ordinal in eachindex(bindings.path_ids)
        bindings.path_ids[ordinal] == path && return ordinal
    end
    throw(PlantPreparationError(:path, :unknown_id,
        "prepared controllable-optic bindings have no optical path $path"))
end

@inline function prepared_controllable_optic_path_slot(
    bindings::PreparedControllableOpticPathBindings, path)
    ordinal = _prepared_controllable_optic_path_ordinal(bindings,
        _as_optical_path_id(path))
    return Int(@inbounds bindings.path_slots[ordinal])
end

@inline function prepared_controllable_optic_binding_range(
    bindings::PreparedControllableOpticPathBindings, path)
    ordinal = _prepared_controllable_optic_path_ordinal(bindings,
        _as_optical_path_id(path))
    first_binding = @inbounds bindings.binding_offsets[ordinal]
    last_binding = @inbounds bindings.binding_offsets[ordinal + 1] - 1
    return first_binding:last_binding
end

@inline function prepared_controllable_optic_plane_group_range(
    bindings::PreparedControllableOpticPathBindings, path)
    ordinal = _prepared_controllable_optic_path_ordinal(bindings,
        _as_optical_path_id(path))
    first_group = @inbounds bindings.group_offsets[ordinal]
    last_group = @inbounds bindings.group_offsets[ordinal + 1] - 1
    return first_group:last_group
end

@inline prepared_controllable_optic_slot(
    bindings::PreparedControllableOpticPathBindings, binding::Integer) =
    Int(bindings.optic_slots[binding])

@inline prepared_controllable_optic_plane_group(
    bindings::PreparedControllableOpticPathBindings, group::Integer) =
    bindings.plane_groups[group]

@inline prepared_controllable_optic_plane_group_path_slot(
    group::PreparedControllableOpticPlaneGroup) = Int(group.path_slot)
@inline prepared_controllable_optic_plane_group_representative_slot(
    group::PreparedControllableOpticPlaneGroup) =
    Int(group.representative_optic_slot)
@inline prepared_controllable_optic_plane_group_binding_range(
    group::PreparedControllableOpticPlaneGroup) =
    group.first_binding:(group.first_binding + group.binding_count - 1)

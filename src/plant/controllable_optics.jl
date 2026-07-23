#
# Prepared controllable optics and endpoint execution configuration
#
# A cold ControllableOpticDefinition describes physical identity and semantic
# command schemas. This layer separately binds run capacities, initial/safe
# values, backend storage, and model-specific physical preparation. Mutable
# optic, endpoint, and workspace owners are constructed by the event loop.
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
physical surface. A successful stage must make
`commit_controllable_optic_command!` an infallible bounded publication step.
"""
function stage_controllable_optic_command!(implementation, state, workspace,
    endpoint::CommandEndpointID, effective_command)
    throw(PlantCommandError(:physical_application,
        :unsupported_optic_application,
        "prepared controllable-optic implementation " *
        "$(typeof(implementation)) does not implement staged command " *
        "application for $endpoint"))
end

"""
Publish one previously staged physical command. Implementations must keep this
bounded and nonthrowing after successful staging so explicit multi-optic
transactions cannot become partially visible.
"""
function commit_controllable_optic_command!(implementation, state, workspace,
    endpoint::CommandEndpointID)
    throw(PlantCommandError(:physical_application,
        :unsupported_optic_commit,
        "prepared controllable-optic implementation " *
        "$(typeof(implementation)) does not implement command commit for " *
        "$endpoint"))
end

"""
Apply the currently visible physical surface to one already materialized
optical-path input. Gate 4 composition calls every co-conjugated optic in
canonical identity order; placement and path visibility are later contracts.
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

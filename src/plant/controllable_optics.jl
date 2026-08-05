#
# Prepared controllable optics and endpoint execution configuration
#
# A cold ControllableOpticDefinition describes physical identity, optical
# placement/path visibility, and semantic command schemas. This layer
# separately binds run capacities and initial/safe values,
# model-specific physical preparation, and canonical path/group lookup.
# Mutable optic, endpoint, and workspace owners are constructed by the event
# loop.
#

"""
    CommandEndpointConfiguration(endpoint, initial_command;
        capacity, sequence_window=capacity, safe_command=nothing)

Run-configuration input for one declared command endpoint. It supplies bounded
calendar/history capacity and host-resident initial and optional safe effective
commands without adding cadence, transport, atomicity, storage placement, or
optical grouping. Cold array inputs are defensively copied before exact-target
storage is derived; device-resident configuration arrays are rejected. Scalar
command values remain host-resident. Stable endpoint ordinals are derived
canonically during plant preparation rather than supplied by declaration
order.
"""
struct CommandEndpointConfiguration{I,S}
    endpoint::CommandEndpointID
    capacity::Int
    sequence_window::Int
    initial_command::I
    safe_command::S
end

function CommandEndpointConfiguration(endpoint, initial_command;
    capacity,
    sequence_window=capacity,
    safe_command=nothing)
    return CommandEndpointConfiguration(
        _as_command_endpoint_id(endpoint),
        _checked_command_endpoint_capacity(capacity),
        _checked_command_sequence_window(sequence_window),
        initial_command,
        safe_command,
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

`definition_slot` binds the cold declaration retained by the prepared
registry. `implementation` is run-immutable preparation data.
`endpoint_slots` names the independent prepared endpoint owners of this
physical device; it does not pack their commands or synchronize their timing.
"""
struct PreparedControllableOptic{
    P,
    S<:FixedSizeVector{UInt32},
}
    definition_slot::UInt32
    implementation::P
    endpoint_slots::S

    function PreparedControllableOptic(
        ::_PreparedControllableOpticToken,
        definition_slot::UInt32,
        implementation::P,
        endpoint_slots::S,
    ) where {P,S<:FixedSizeVector{UInt32}}
        iszero(definition_slot) && throw(PlantPreparationError(
            :controllable_optic,
            :definition_slot,
            "prepared controllable-optic definition slot must be positive",
        ))
        return new{P,S}(definition_slot, implementation, endpoint_slots)
    end
end

"""Compact location of one controllable optic in concrete family storage."""
struct _PreparedControllableOpticSlot
    family_slot::UInt32
    member_slot::UInt32
end

"""Fixed-capacity prepared optics with one exact concrete element type."""
struct _PreparedControllableOpticFamily{
    O<:PreparedControllableOptic,
    V<:FixedSizeVector{O},
}
    values::V
end


"""Fixed-capacity persistent states for one exact optic family."""
struct _ControllableOpticStateFamily{T,V<:FixedSizeVector{T}}
    values::V
end


"""Fixed-capacity replaceable workspaces for one exact optic family."""
struct _ControllableOpticWorkspaceFamily{T,V<:FixedSizeVector{T}}
    values::V
end

"""
Definition-order view over concrete prepared controllable-optic families.

The family tuple type depends on the bounded set of exact execution families,
while `slots` carries runtime optic cardinality. Armed homogeneous storage is
therefore concrete without specializing the registry type on the number of
optics in a family.
"""
struct _PreparedControllableOpticRegistry{
    D<:_FixedPlantRegistry{ControllableOpticDefinition},
    G<:Tuple,
    S<:FixedSizeVector{_PreparedControllableOpticSlot},
}
    definitions::D
    groups::G
    slots::S
end

"""Persistent physical states grouped exactly like prepared optic families."""
struct _ControllableOpticStateRegistry{
    G<:Tuple,
    S<:FixedSizeVector{_PreparedControllableOpticSlot},
}
    groups::G
    slots::S
end


"""Replaceable physical workspaces grouped like prepared optic families."""
struct _ControllableOpticWorkspaceRegistry{
    G<:Tuple,
    S<:FixedSizeVector{_PreparedControllableOpticSlot},
}
    groups::G
    slots::S
end

Base.size(registry::_PreparedControllableOpticRegistry) =
    (length(registry.slots),)
Base.axes(registry::_PreparedControllableOpticRegistry) =
    axes(registry.slots)
Base.length(registry::_PreparedControllableOpticRegistry) =
    length(registry.slots)
Base.keys(registry::_PreparedControllableOpticRegistry) =
    eachindex(registry.slots)
Base.eachindex(registry::_PreparedControllableOpticRegistry) =
    eachindex(registry.slots)
Base.firstindex(registry::_PreparedControllableOpticRegistry) =
    firstindex(registry.slots)
Base.lastindex(registry::_PreparedControllableOpticRegistry) =
    lastindex(registry.slots)

Base.size(registry::_ControllableOpticStateRegistry) =
    (length(registry.slots),)
Base.axes(registry::_ControllableOpticStateRegistry) =
    axes(registry.slots)
Base.length(registry::_ControllableOpticStateRegistry) =
    length(registry.slots)

Base.size(registry::_ControllableOpticWorkspaceRegistry) =
    (length(registry.slots),)
Base.axes(registry::_ControllableOpticWorkspaceRegistry) =
    axes(registry.slots)
Base.length(registry::_ControllableOpticWorkspaceRegistry) =
    length(registry.slots)

Base.keys(registry::Union{
    _ControllableOpticStateRegistry,
    _ControllableOpticWorkspaceRegistry,
}) = eachindex(registry.slots)
Base.eachindex(registry::Union{
    _ControllableOpticStateRegistry,
    _ControllableOpticWorkspaceRegistry,
}) = eachindex(registry.slots)
Base.firstindex(registry::Union{
    _ControllableOpticStateRegistry,
    _ControllableOpticWorkspaceRegistry,
}) = firstindex(registry.slots)
Base.lastindex(registry::Union{
    _ControllableOpticStateRegistry,
    _ControllableOpticWorkspaceRegistry,
}) = lastindex(registry.slots)

@noinline function _prepared_controllable_optic_slot_error(
    family::Int,
    member::Int,
)
    throw(BoundsError((family_slot=family, member_slot=member)))
end

@inline function _controllable_optic_runtime_value(
    registry::Union{
        _ControllableOpticStateRegistry,
        _ControllableOpticWorkspaceRegistry,
    },
    index::Int,
)
    checkbounds(registry.slots, index)
    slot = @inbounds registry.slots[index]
    return _prepared_controllable_optic_family_value(
        registry.groups,
        Int(slot.family_slot),
        Int(slot.member_slot),
    )
end

@inline Base.getindex(
    registry::_ControllableOpticStateRegistry,
    index::Int,
) = _controllable_optic_runtime_value(registry, index)

@inline Base.getindex(
    registry::_ControllableOpticWorkspaceRegistry,
    index::Int,
) = _controllable_optic_runtime_value(registry, index)

@inline function Base.iterate(
    registry::Union{
        _ControllableOpticStateRegistry,
        _ControllableOpticWorkspaceRegistry,
    },
    state::Int=1,
)
    state > length(registry) && return nothing
    return (@inbounds registry[state], state + 1)
end

@inline function _prepared_controllable_optic_family_value(
    groups::Tuple{},
    family::Int,
    member::Int,
)
    return _prepared_controllable_optic_slot_error(family, member)
end

@inline function _prepared_controllable_optic_family_value(
    groups::Tuple,
    family::Int,
    member::Int,
)
    family == 1 && return @inbounds groups[1].values[member]
    return _prepared_controllable_optic_family_value(
        Base.tail(groups), family - 1, member)
end

@inline function Base.getindex(
    registry::_PreparedControllableOpticRegistry,
    index::Int,
)
    checkbounds(registry.slots, index)
    slot = @inbounds registry.slots[index]
    return _prepared_controllable_optic_family_value(
        registry.groups,
        Int(slot.family_slot),
        Int(slot.member_slot),
    )
end

@inline function _prepared_controllable_optic_definition(
    registry::_PreparedControllableOpticRegistry,
    optic::PreparedControllableOptic,
)
    slot = Int(optic.definition_slot)
    checkbounds(registry.definitions, slot)
    return @inbounds registry.definitions[slot]
end

@inline function _prepared_controllable_optic_definition(
    registry::_PreparedControllableOpticRegistry,
    index::Int,
)
    return _prepared_controllable_optic_definition(
        registry, registry[index])
end

@inline function Base.iterate(
    registry::_PreparedControllableOpticRegistry,
    state::Int=1,
)
    state > length(registry) && return nothing
    return (@inbounds registry[state], state + 1)
end

struct _PreparedControllableOpticFamilyType{O<:PreparedControllableOptic} end

@inline _prepared_controllable_optic_family_matches(
    ::_PreparedControllableOpticFamilyType{O},
    ::Type{O},
) where {O<:PreparedControllableOptic} = true

@inline _prepared_controllable_optic_family_matches(
    ::_PreparedControllableOpticFamilyType,
    ::Type{<:PreparedControllableOptic},
) = false

@inline function _has_prepared_controllable_optic_family(
    ::Tuple{},
    ::Type{<:PreparedControllableOptic},
)
    return false
end

@inline function _has_prepared_controllable_optic_family(
    families::Tuple,
    ::Type{O},
) where {O<:PreparedControllableOptic}
    _prepared_controllable_optic_family_matches(families[1], O) &&
        return true
    return _has_prepared_controllable_optic_family(Base.tail(families), O)
end

function _append_prepared_controllable_optic_family(
    families::Tuple,
    optic::O,
) where {O<:PreparedControllableOptic}
    _has_prepared_controllable_optic_family(families, O) && return families
    return (families..., _PreparedControllableOpticFamilyType{O}())
end

function _prepared_controllable_optic_family_types(optics)
    Base.@nospecialize optics
    families = ()
    @inbounds for optic in optics
        families = _append_prepared_controllable_optic_family(
            families, optic)
    end
    return families
end

function _prepare_controllable_optic_family(
    ::_PreparedControllableOpticFamilyType{O},
    optics,
) where {O<:PreparedControllableOptic}
    Base.@nospecialize optics
    count = 0
    @inbounds for optic in optics
        typeof(optic) === O && (count += 1)
    end
    values = Vector{O}(undef, count)
    next = 1
    @inbounds for optic in optics
        typeof(optic) === O || continue
        values[next] = optic
        next += 1
    end
    fixed = FixedSizeVectorDefault{O}(values)
    return _PreparedControllableOpticFamily{O,typeof(fixed)}(fixed)
end

@inline _prepare_controllable_optic_families(::Tuple{}, optics) = ()

function _prepare_controllable_optic_families(families::Tuple, optics)
    Base.@nospecialize optics
    first = _prepare_controllable_optic_family(families[1], optics)
    rest = _prepare_controllable_optic_families(Base.tail(families), optics)
    return (first, rest...)
end

@inline function _prepared_controllable_optic_family_index(
    ::Tuple{},
    ::Type{<:PreparedControllableOptic},
    ::Int=1,
)
    return 0
end


@inline function _prepared_controllable_optic_family_index(
    families::Tuple,
    ::Type{O},
    index::Int=1,
) where {O<:PreparedControllableOptic}
    _prepared_controllable_optic_family_matches(families[1], O) &&
        return index
    return _prepared_controllable_optic_family_index(
        Base.tail(families), O, index + 1)
end

function _prepare_controllable_optic_slots(families::Tuple, optics)
    Base.@nospecialize optics
    length(families) <= typemax(UInt32) || throw(PlantPreparationError(
        :controllable_optic,
        :family_capacity,
        "prepared controllable-optic family count exceeds UInt32 capacity",
    ))
    counts = zeros(UInt32, length(families))
    slots = Vector{_PreparedControllableOpticSlot}(undef, length(optics))
    @inbounds for index in eachindex(optics)
        family = _prepared_controllable_optic_family_index(
            families, typeof(optics[index]))
        iszero(family) && throw(PlantPreparationError(
            :controllable_optic,
            :missing_family,
            "prepared controllable optic has no concrete family",
        ))
        counts[family] == typemax(UInt32) && throw(PlantPreparationError(
            :controllable_optic,
            :family_capacity,
            "prepared controllable-optic family exceeds UInt32 capacity",
        ))
        counts[family] += UInt32(1)
        slots[index] = _PreparedControllableOpticSlot(
            UInt32(family), counts[family])
    end
    return FixedSizeVectorDefault{_PreparedControllableOpticSlot}(slots)
end

function _prepare_controllable_optic_registry(
    definitions::_FixedPlantRegistry{ControllableOpticDefinition},
    optics,
)
    Base.@nospecialize optics
    family_types = _prepared_controllable_optic_family_types(optics)
    groups = _prepare_controllable_optic_families(family_types, optics)
    slots = _prepare_controllable_optic_slots(family_types, optics)
    return _PreparedControllableOpticRegistry(definitions, groups, slots)
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
    prepare_target_local_controllable_optic(
        model, definition, telescope, atmosphere_definition, target)

Fail-closed extension seam for immutable physical-optic preparation in a
target-local partition.  Unlike `prepare_controllable_optic`, this seam does
not receive or create a mutable timed-atmosphere owner.  Models may use the
cold atmosphere definition for static conjugate geometry; every numerical
array they return must occupy `target`.
"""
function prepare_target_local_controllable_optic(
    model,
    ::ControllableOpticDefinition,
    ::AbstractTelescope,
    ::AbstractTimedAtmosphereDefinition,
    ::AbstractComputeDevice,
)
    throw(PlantPreparationError(
        :controllable_optic,
        :unsupported_target_local_model,
        "controllable-optic model $(typeof(model)) does not implement " *
        "prepare_target_local_controllable_optic",
    ))
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
optical-path input through its immutable path-local coupling. Composition calls
only the path's prepared visible `PupilSurfaceExecutionRole` bindings, groups
compatible contiguous couplings, and preserves canonical optic identity order.
Path-local autonomous devices use their separately prepared exact coupling.

The four-argument extension seam is defined with the geometric coupling
implementation after path preparation types are available.
"""

function _copy_prepared_effective_command(
    endpoint::PreparedCommandEndpoint{<:PlantCommandSchema{T,0}},
    value, label::AbstractString, ::AbstractComputeDevice) where {T}
    return _validate_effective_seed(command_schema(endpoint), value, label)
end

function _copy_prepared_effective_command(
    endpoint::PreparedCommandEndpoint{<:PlantCommandSchema{T,N}},
    value, label::AbstractString,
    target::AbstractComputeDevice,
) where {T,N}
    validated = _validate_effective_seed(command_schema(endpoint), value,
        label)
    copied = allocate_device_array(target, T,
        command_dimensions(command_schema(endpoint))...)
    copyto!(copied, validated)
    return copied
end

function _copy_prepared_safe_command(endpoint::PreparedCommandEndpoint,
    ::Nothing, ::AbstractComputeDevice)
    policy = command_silence_policy(command_schema(endpoint))
    _require_safe_command_configuration(policy, nothing)
    return nothing
end

function _copy_prepared_safe_command(endpoint::PreparedCommandEndpoint,
    value, target::AbstractComputeDevice)
    policy = command_silence_policy(command_schema(endpoint))
    _require_safe_command_configuration(policy, value)
    return _copy_prepared_effective_command(
        endpoint, value, "safe command", target)
end

@inline _command_endpoint_backend(
    ::PlantCommandSchema{<:Any,0},
    ::AbstractComputeDevice,
) = CPUBackend()

@inline _command_endpoint_backend(
    ::PlantCommandSchema,
    target::AbstractComputeDevice,
) = compute_device_backend(target)

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
    configurations, optic_definitions, target::AbstractComputeDevice)
    declared = _canonical_command_endpoint_declarations(definition)
    length(configurations) == length(declared) || throw(
        PlantPreparationError(:command_endpoint,
            :configuration_count,
            "plant declares $(length(declared)) command endpoints but " *
            "$(length(configurations)) run configurations were supplied"))
    endpoints = Memory{_PreparedPlantCommandEndpoint}(
        undef, length(declared))
    for (ordinal, (owner, schema)) in enumerate(declared)
        configuration = _command_endpoint_configuration(configurations,
            command_endpoint_id(schema))
        endpoint = prepare_command_endpoint(schema;
            capacity=configuration.capacity,
            sequence_window=configuration.sequence_window,
            ordinal,
            backend=_command_endpoint_backend(schema, target))
        initial = _copy_prepared_effective_command(endpoint,
            configuration.initial_command, "initial effective command",
            target)
        safe = _copy_prepared_safe_command(endpoint,
            configuration.safe_command, target)
        optic_slot = _controllable_optic_slot(optic_definitions, owner)
        endpoints[ordinal] = _PreparedPlantCommandEndpoint(
            endpoint, UInt32(optic_slot), initial, safe)
    end
    return endpoints
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
    optic_definitions, endpoints, telescope, atmosphere)
    declared_definitions = controllable_optic_definitions(definition)
    length(declared_definitions) <= typemax(UInt32) || throw(
        PlantPreparationError(
            :controllable_optic,
            :capacity,
            "prepared controllable-optic count exceeds UInt32 capacity",
        ))
    optics = ()
    for index in eachindex(optic_definitions)
        optic_definition = optic_definitions[index]
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
        endpoint_slots = UInt32[
            _prepared_command_endpoint_slot(endpoints,
                command_endpoint_id(schema))
            for schema in command_schemas(optic_definition)
        ]
        fixed_endpoint_slots =
            FixedSizeVectorDefault{UInt32}(endpoint_slots)
        definition_slot = _controllable_optic_slot(
            declared_definitions,
            controllable_optic_id(optic_definition),
        )
        optic = PreparedControllableOptic(
            _PREPARED_CONTROLLABLE_OPTIC_TOKEN,
            UInt32(definition_slot),
            implementation,
            fixed_endpoint_slots,
        )
        optics = (optics..., optic)
    end
    return _prepare_controllable_optic_registry(
        declared_definitions, optics)
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

function _canonical_prepared_path_slots(paths::AbstractVector)
    slots = collect(eachindex(paths))
    sort!(slots; by=slot ->
        String(path_id(paths[slot].definition).name))
    return slots
end

@inline function _optic_binding_slot_isless(left::UInt32, right::UInt32,
    optics)
    left_definition =
        _prepared_controllable_optic_definition(optics, Int(left))
    right_definition =
        _prepared_controllable_optic_definition(optics, Int(right))
    left_placement = controllable_optic_placement(left_definition)
    right_placement = controllable_optic_placement(right_definition)
    if _same_optical_placement(left_placement, right_placement)
        return String(controllable_optic_id(left_definition).name) <
            String(controllable_optic_id(right_definition).name)
    end
    return _optical_placement_isless(left_placement, right_placement)
end

function _visible_prepared_optic_slots(
    optics, path::OpticalPathID)
    slots = UInt32[]
    sizehint!(slots, length(optics))
    @inbounds for slot in eachindex(optics)
        definition = _prepared_controllable_optic_definition(optics, slot)
        _visible_on_path(controllable_optic_visibility(definition), path) ||
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
    optics,
)
    isempty(visible_slots) && return nothing
    first_binding = length(optic_slots) + 1
    first_slot = first(visible_slots)
    previous_definition =
        _prepared_controllable_optic_definition(optics, Int(first_slot))
    previous_placement = controllable_optic_placement(previous_definition)
    group_first = first_binding
    group_representative = first_slot
    @inbounds for slot in visible_slots
        definition =
            _prepared_controllable_optic_definition(optics, Int(slot))
        placement = controllable_optic_placement(definition)
        if !_same_optical_placement(previous_placement, placement)
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
    optics, paths::AbstractVector)
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

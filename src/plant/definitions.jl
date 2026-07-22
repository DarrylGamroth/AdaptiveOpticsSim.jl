#
# Immutable plant topology declarations
#
# These values describe stable physical identity, topology, and cold semantic
# command contracts only. Backend-bound preparation, mutable state/workspace
# ownership, execution, and timing belong to later layers.
#

@inline _require_declared_topology_value(value, ::Symbol, ::Symbol,
    ::AbstractString) = value

function _require_declared_topology_value(::Nothing, component::Symbol,
    reason::Symbol, label::AbstractString)
    throw(PlantDefinitionError(component, reason,
        "$component $label must be declared"))
end

"""
    ColdPlantModelDefinition()

Accepted result of `plant_model_definition_style` for a configuration-only
controllable-optic, optical-path, or acquisition model definition.
"""
struct ColdPlantModelDefinition end
struct _UnsupportedPlantModelDefinition end

"""
    plant_model_definition_style(::Type{T})

Trait for values stored as the optical or acquisition model in a plant
definition. The default rejects the type. A cold, immutable model-definition
type opts in by returning `ColdPlantModelDefinition()`:

```julia
AdaptiveOpticsSim.plant_model_definition_style(::Type{MyModelDefinition}) =
    ColdPlantModelDefinition()
```

Opting in asserts that instances contain configuration only: no prepared
workspace, mutable simulation or detector state, schedule, RNG stream, queue,
transport, or HIL descriptor. Preparation may use the declaration to construct
separately owned plans, state, and workspaces.
"""
plant_model_definition_style(::Type) = _UnsupportedPlantModelDefinition()

@inline function _require_cold_plant_model_definition(model,
    component::Symbol, label::AbstractString)
    style = plant_model_definition_style(typeof(model))
    return _require_cold_plant_model_definition(style, model, component,
        label)
end

@inline _require_cold_plant_model_definition(::ColdPlantModelDefinition,
    model, ::Symbol, ::AbstractString) = model

function _require_cold_plant_model_definition(
    ::_UnsupportedPlantModelDefinition, model, component::Symbol,
    label::AbstractString)
    throw(PlantDefinitionError(component, :unsupported_model_definition,
        "$component $label type $(typeof(model)) has not declared the cold " *
        "plant-model-definition contract"))
end

function _require_cold_plant_model_definition(style, model,
    component::Symbol, label::AbstractString)
    throw(PlantDefinitionError(component, :invalid_model_definition_style,
        "plant_model_definition_style($(typeof(model))) must return " *
        "ColdPlantModelDefinition(); got $(typeof(style))"))
end

@inline _require_path_source(::AbstractSource) = nothing

function _require_path_source(value)
    throw(PlantDefinitionError(:path, :invalid_source,
        "optical-path source must implement AbstractSource; got " *
        "$(typeof(value))"))
end

"""
    OpticalPathDefinition(id, source, optical_model)
    OpticalPathDefinition(id, source; optical_model)

Immutable declaration of one reusable optical path. `id` is explicit physical
identity; tuple or named-tuple position is never identity. `source` and
`optical_model` are cold model declarations. They own no prepared propagation
workspace or mutable acquisition state. The concrete optical-model type must
opt in through `plant_model_definition_style`.
"""
struct OpticalPathDefinition{S,M}
    id::OpticalPathID
    source::S
    optical_model::M

    function OpticalPathDefinition(id::OpticalPathID, source::S,
        optical_model::M) where {S,M}
        _require_declared_topology_value(source, :path, :missing_source,
            "source")
        _require_path_source(source)
        _require_declared_topology_value(optical_model, :path,
            :missing_model, "optical model")
        _require_cold_plant_model_definition(optical_model, :path,
            "optical model")
        return new{S,M}(id, source, optical_model)
    end
end

function OpticalPathDefinition(id, source, optical_model)
    return OpticalPathDefinition(_as_optical_path_id(id), source,
        optical_model)
end

OpticalPathDefinition(id, source; optical_model) =
    OpticalPathDefinition(id, source, optical_model)

"""
    AcquisitionDefinition(id, path, acquisition_model)
    AcquisitionDefinition(id, path; acquisition_model)

Immutable declaration of one independently invocable acquisition. `path`
references an `OpticalPathID`; acquisition preparation, detector/WFS mutable
state, timing, RNG, and publication ownership are deliberately absent. The
concrete acquisition-model type must opt in through
`plant_model_definition_style`.
"""
struct AcquisitionDefinition{M}
    id::AcquisitionID
    path::OpticalPathID
    acquisition_model::M

    function AcquisitionDefinition(id::AcquisitionID, path::OpticalPathID,
        acquisition_model::M) where {M}
        _require_declared_topology_value(acquisition_model, :acquisition,
            :missing_model, "model")
        _require_cold_plant_model_definition(acquisition_model, :acquisition,
            "model")
        return new{M}(id, path, acquisition_model)
    end
end

function AcquisitionDefinition(id, path, acquisition_model)
    return AcquisitionDefinition(_as_acquisition_id(id),
        _as_optical_path_id(path), acquisition_model)
end

AcquisitionDefinition(id, path; acquisition_model) =
    AcquisitionDefinition(id, path, acquisition_model)

@inline _require_command_schema(::PlantCommandSchema) = nothing

function _require_command_schema(value)
    throw(PlantDefinitionError(:command_schema, :invalid_definition,
        "controllable-optic command schemas must contain " *
        "PlantCommandSchema values; got $(typeof(value))"))
end

function _require_named_command_schema_identity(key::Symbol,
    schema::PlantCommandSchema)
    key == command_endpoint_id(schema).name || throw(PlantDefinitionError(
        :command_schema, :identity_mismatch,
        "named command-schema key $(repr(key)) does not match endpoint " *
        "$(command_endpoint_id(schema))"))
    return nothing
end

function _require_unique_command_schemas(schemas::Tuple)
    seen = Set{CommandEndpointID}()
    seen_schemas = Set{PlantCommandSchemaID}()
    for schema in schemas
        endpoint = command_endpoint_id(schema)
        endpoint in seen && throw(PlantDefinitionError(:command_endpoint,
            :duplicate_id,
            "duplicate command-endpoint identity $endpoint in one " *
            "controllable-optic declaration"))
        push!(seen, endpoint)
        id = command_schema_id(schema)
        id in seen_schemas && throw(PlantDefinitionError(:command_schema,
            :duplicate_id,
            "duplicate plant-command schema identity $id in one " *
            "controllable-optic declaration"))
        push!(seen_schemas, id)
    end
    return nothing
end

function _normalize_command_schemas(schemas::Tuple)
    isempty(schemas) && throw(PlantDefinitionError(:controllable_optic,
        :missing_schema,
        "a controllable optic must declare at least one command schema"))
    foreach(_require_command_schema, schemas)
    _require_unique_command_schemas(schemas)
    return schemas
end

function _normalize_command_schemas(schemas::NamedTuple)
    normalized = values(schemas)
    isempty(normalized) && throw(PlantDefinitionError(:controllable_optic,
        :missing_schema,
        "a controllable optic must declare at least one command schema"))
    foreach(_require_command_schema, normalized)
    foreach(_require_named_command_schema_identity, keys(schemas), normalized)
    _require_unique_command_schemas(normalized)
    return normalized
end

function _normalize_command_schemas(schemas)
    throw(PlantDefinitionError(:command_schema, :invalid_container,
        "controllable-optic command schemas must be a Tuple or NamedTuple; " *
        "got $(typeof(schemas))"))
end

"""
    ControllableOpticDefinition(id, optic_model, command_schemas)
    ControllableOpticDefinition(id, optic_model; command_schemas)

Immutable declaration of one physical controllable optic and the stable
semantic schemas of its independently timed or latched command endpoints. The
optic model is cold configuration and must opt in through
`plant_model_definition_style`. Mutable command state, timing, placement, path
visibility, and prepared optical grouping are deliberately absent from this
topology record.
"""
struct ControllableOpticDefinition{M,S<:Tuple}
    id::ControllableOpticID
    optic_model::M
    command_schemas::S

    function ControllableOpticDefinition(id::ControllableOpticID,
        optic_model::M, command_schemas::Tuple) where {M}
        _require_declared_topology_value(optic_model, :controllable_optic,
            :missing_model, "model")
        _require_cold_plant_model_definition(optic_model,
            :controllable_optic, "model")
        schemas = _normalize_command_schemas(command_schemas)
        return new{M,typeof(schemas)}(id, optic_model, schemas)
    end
end

function ControllableOpticDefinition(id, optic_model, command_schemas)
    schemas = _normalize_command_schemas(command_schemas)
    return ControllableOpticDefinition(_as_controllable_optic_id(id),
        optic_model, schemas)
end

ControllableOpticDefinition(id, optic_model; command_schemas) =
    ControllableOpticDefinition(id, optic_model, command_schemas)

@inline path_id(path::OpticalPathDefinition) = path.id
@inline acquisition_id(acquisition::AcquisitionDefinition) = acquisition.id
@inline controllable_optic_id(optic::ControllableOpticDefinition) = optic.id
@inline acquisition_path_id(acquisition::AcquisitionDefinition) =
    acquisition.path
@inline path_source(path::OpticalPathDefinition) = path.source
@inline path_model(path::OpticalPathDefinition) = path.optical_model
@inline acquisition_model(acquisition::AcquisitionDefinition) =
    acquisition.acquisition_model
@inline controllable_optic_model(optic::ControllableOpticDefinition) =
    optic.optic_model
@inline command_schemas(optic::ControllableOpticDefinition) =
    optic.command_schemas
@inline command_endpoint_ids(optic::ControllableOpticDefinition) =
    map(command_endpoint_id, optic.command_schemas)

function command_schema(optic::ControllableOpticDefinition, id)
    resolved = _as_command_endpoint_id(id)
    for schema in optic.command_schemas
        command_endpoint_id(schema) == resolved && return schema
    end
    throw(PlantDefinitionError(:command_endpoint, :unknown_id,
        "controllable optic $(optic.id) has no command endpoint $resolved"))
end

@inline _require_path_definition(::OpticalPathDefinition) = nothing

function _require_path_definition(value)
    throw(PlantDefinitionError(:path, :invalid_definition,
        "plant paths must contain OpticalPathDefinition values; got " *
        "$(typeof(value))"))
end

@inline _require_acquisition_definition(::AcquisitionDefinition) = nothing

function _require_acquisition_definition(value)
    throw(PlantDefinitionError(:acquisition, :invalid_definition,
        "plant acquisitions must contain AcquisitionDefinition values; got " *
        "$(typeof(value))"))
end

@inline _require_controllable_optic_definition(
    ::ControllableOpticDefinition) = nothing

function _require_controllable_optic_definition(value)
    throw(PlantDefinitionError(:controllable_optic, :invalid_definition,
        "plant controllable optics must contain ControllableOpticDefinition " *
        "values; got $(typeof(value))"))
end

function _require_named_path_identity(key::Symbol,
    path::OpticalPathDefinition)
    key == path.id.name || throw(PlantDefinitionError(:path,
        :identity_mismatch,
        "named path key $(repr(key)) does not match $(path.id)"))
    return nothing
end

function _require_named_acquisition_identity(key::Symbol,
    acquisition::AcquisitionDefinition)
    key == acquisition.id.name || throw(PlantDefinitionError(:acquisition,
        :identity_mismatch,
        "named acquisition key $(repr(key)) does not match " *
        "$(acquisition.id)"))
    return nothing
end

function _require_named_controllable_optic_identity(key::Symbol,
    optic::ControllableOpticDefinition)
    key == optic.id.name || throw(PlantDefinitionError(:controllable_optic,
        :identity_mismatch,
        "named controllable-optic key $(repr(key)) does not match " *
        "$(optic.id)"))
    return nothing
end

function _normalize_path_definitions(paths::Tuple)
    foreach(_require_path_definition, paths)
    return paths
end

function _normalize_path_definitions(paths::NamedTuple)
    foreach(_require_path_definition, values(paths))
    foreach(_require_named_path_identity, keys(paths), values(paths))
    return values(paths)
end

function _normalize_path_definitions(paths)
    throw(PlantDefinitionError(:path, :invalid_container,
        "plant paths must be a Tuple or NamedTuple; got $(typeof(paths))"))
end

function _normalize_acquisition_definitions(acquisitions::Tuple)
    foreach(_require_acquisition_definition, acquisitions)
    return acquisitions
end

function _normalize_acquisition_definitions(acquisitions::NamedTuple)
    foreach(_require_acquisition_definition, values(acquisitions))
    foreach(_require_named_acquisition_identity, keys(acquisitions),
        values(acquisitions))
    return values(acquisitions)
end

function _normalize_acquisition_definitions(acquisitions)
    throw(PlantDefinitionError(:acquisition, :invalid_container,
        "plant acquisitions must be a Tuple or NamedTuple; got " *
        "$(typeof(acquisitions))"))
end

function _normalize_controllable_optic_definitions(optics::Tuple)
    foreach(_require_controllable_optic_definition, optics)
    return optics
end

function _normalize_controllable_optic_definitions(optics::NamedTuple)
    foreach(_require_controllable_optic_definition, values(optics))
    foreach(_require_named_controllable_optic_identity, keys(optics),
        values(optics))
    return values(optics)
end

function _normalize_controllable_optic_definitions(optics)
    throw(PlantDefinitionError(:controllable_optic, :invalid_container,
        "plant controllable optics must be a Tuple or NamedTuple; got " *
        "$(typeof(optics))"))
end

@inline _require_plant_telescope(::AbstractTelescope) = nothing

function _require_plant_telescope(value)
    throw(PlantDefinitionError(:plant, :invalid_telescope,
        "plant telescope must implement AbstractTelescope; got " *
        "$(typeof(value))"))
end

@inline _require_plant_atmosphere(::AbstractAtmosphere) = nothing

function _require_plant_atmosphere(value)
    throw(PlantDefinitionError(:plant, :invalid_atmosphere,
        "plant atmosphere must implement AbstractAtmosphere; got " *
        "$(typeof(value))"))
end

function _require_unique_path_ids(paths::Tuple)
    seen = Set{OpticalPathID}()
    for path in paths
        id = path_id(path)
        id in seen && throw(PlantDefinitionError(:path, :duplicate_id,
            "duplicate optical-path identity $id"))
        push!(seen, id)
    end
    return nothing
end

function _require_unique_acquisition_ids(acquisitions::Tuple)
    seen = Set{AcquisitionID}()
    for acquisition in acquisitions
        id = acquisition_id(acquisition)
        id in seen && throw(PlantDefinitionError(:acquisition, :duplicate_id,
            "duplicate acquisition identity $id"))
        push!(seen, id)
    end
    return nothing
end

function _require_unique_controllable_optic_ids(optics::Tuple)
    seen = Set{ControllableOpticID}()
    for optic in optics
        id = controllable_optic_id(optic)
        id in seen && throw(PlantDefinitionError(:controllable_optic,
            :duplicate_id, "duplicate controllable-optic identity $id"))
        push!(seen, id)
    end
    return nothing
end

function _require_unique_command_endpoint_owners(optics::Tuple)
    seen = Set{CommandEndpointID}()
    for optic in optics
        for schema in command_schemas(optic)
            endpoint = command_endpoint_id(schema)
            endpoint in seen && throw(PlantDefinitionError(
                :command_endpoint, :duplicate_owner,
                "command endpoint $endpoint is owned by more than one " *
                "controllable optic"))
            push!(seen, endpoint)
        end
    end
    return nothing
end

function _require_unique_plant_command_schema_ids(optics::Tuple)
    seen = Set{PlantCommandSchemaID}()
    for optic in optics
        for schema in command_schemas(optic)
            id = command_schema_id(schema)
            id in seen && throw(PlantDefinitionError(:command_schema,
                :duplicate_id,
                "duplicate plant-command schema identity $id"))
            push!(seen, id)
        end
    end
    return nothing
end

function _contains_path_id(paths::Tuple, id::OpticalPathID)
    for path in paths
        path_id(path) == id && return true
    end
    return false
end

function _require_acquisition_paths(paths::Tuple, acquisitions::Tuple)
    for acquisition in acquisitions
        id = acquisition_path_id(acquisition)
        _contains_path_id(paths, id) || throw(PlantDefinitionError(
            :acquisition, :unknown_path,
            "acquisition $(acquisition.id) references unknown path $id"))
    end
    return nothing
end

"""
    PlantDefinition(telescope, atmosphere, controllable_optics, paths,
        acquisitions)
    PlantDefinition(; telescope, atmosphere, controllable_optics=(),
        paths=(), acquisitions=())

Immutable declared topology for one telescope and atmosphere, reusable optical
paths, independent acquisitions, and independently identified controllable
optics with versioned semantic command schemas. Tuples and named tuples are
accepted as cold organization only; every component carries its own stable
identity. This value is not prepared execution state and owns no mutable
command state, schedule, queue, transport, RNG stream, or HIL descriptor.
"""
struct PlantDefinition{T,A,O<:Tuple,P<:Tuple,Q<:Tuple}
    telescope::T
    atmosphere::A
    controllable_optics::O
    paths::P
    acquisitions::Q

    function PlantDefinition(telescope::T, atmosphere::A,
        controllable_optics::O, paths::P,
        acquisitions::Q) where {T,A,O<:Tuple,P<:Tuple,Q<:Tuple}
        _require_plant_telescope(telescope)
        _require_plant_atmosphere(atmosphere)
        foreach(_require_controllable_optic_definition,
            controllable_optics)
        foreach(_require_path_definition, paths)
        foreach(_require_acquisition_definition, acquisitions)
        _require_unique_controllable_optic_ids(controllable_optics)
        _require_unique_command_endpoint_owners(controllable_optics)
        _require_unique_plant_command_schema_ids(controllable_optics)
        _require_unique_path_ids(paths)
        _require_unique_acquisition_ids(acquisitions)
        _require_acquisition_paths(paths, acquisitions)
        return new{T,A,O,P,Q}(telescope, atmosphere, controllable_optics,
            paths, acquisitions)
    end
end

function PlantDefinition(; telescope, atmosphere, controllable_optics=(),
    paths=(), acquisitions=())
    normalized_optics =
        _normalize_controllable_optic_definitions(controllable_optics)
    normalized_paths = _normalize_path_definitions(paths)
    normalized_acquisitions = _normalize_acquisition_definitions(acquisitions)
    return PlantDefinition(telescope, atmosphere, normalized_optics,
        normalized_paths, normalized_acquisitions)
end

function PlantDefinition(telescope, atmosphere, controllable_optics, paths,
    acquisitions)
    normalized_optics =
        _normalize_controllable_optic_definitions(controllable_optics)
    normalized_paths = _normalize_path_definitions(paths)
    normalized_acquisitions = _normalize_acquisition_definitions(acquisitions)
    return PlantDefinition(telescope, atmosphere, normalized_optics,
        normalized_paths, normalized_acquisitions)
end

@inline plant_telescope(plant::PlantDefinition) = plant.telescope
@inline plant_atmosphere(plant::PlantDefinition) = plant.atmosphere
@inline controllable_optic_definitions(plant::PlantDefinition) =
    plant.controllable_optics
@inline path_definitions(plant::PlantDefinition) = plant.paths
@inline acquisition_definitions(plant::PlantDefinition) = plant.acquisitions

function path_definition(plant::PlantDefinition, id)
    resolved = _as_optical_path_id(id)
    for path in plant.paths
        path_id(path) == resolved && return path
    end
    throw(PlantDefinitionError(:path, :unknown_id,
        "plant has no optical path $resolved"))
end

function acquisition_definition(plant::PlantDefinition, id)
    resolved = _as_acquisition_id(id)
    for acquisition in plant.acquisitions
        acquisition_id(acquisition) == resolved && return acquisition
    end
    throw(PlantDefinitionError(:acquisition, :unknown_id,
        "plant has no acquisition $resolved"))
end

function controllable_optic_definition(plant::PlantDefinition, id)
    resolved = _as_controllable_optic_id(id)
    for optic in plant.controllable_optics
        controllable_optic_id(optic) == resolved && return optic
    end
    throw(PlantDefinitionError(:controllable_optic, :unknown_id,
        "plant has no controllable optic $resolved"))
end

function command_endpoint_owner(plant::PlantDefinition, id)
    resolved = _as_command_endpoint_id(id)
    for optic in plant.controllable_optics
        for schema in command_schemas(optic)
            command_endpoint_id(schema) == resolved && return optic
        end
    end
    throw(PlantDefinitionError(:command_endpoint, :unknown_id,
        "plant has no command endpoint $resolved"))
end


function command_schema(plant::PlantDefinition, id)
    resolved = _as_command_endpoint_id(id)
    for optic in plant.controllable_optics
        for schema in command_schemas(optic)
            command_endpoint_id(schema) == resolved && return schema
        end
    end
    throw(PlantDefinitionError(:command_endpoint, :unknown_id,
        "plant has no command endpoint $resolved"))
end

@inline _as_plant_command_schema_lookup_id(id::PlantCommandSchemaID) = id

function _as_plant_command_schema_lookup_id(name::Symbol)
    isempty(String(name)) && throw(PlantDefinitionError(:command_schema,
        :empty_id, "plant-command schema identity must not be empty"))
    return PlantCommandSchemaID(name)
end

function _as_plant_command_schema_lookup_id(value)
    throw(PlantDefinitionError(:command_schema, :invalid_id,
        "plant-command schema identity must be a Symbol or " *
        "PlantCommandSchemaID; got $(typeof(value))"))
end

function plant_command_schema(plant::PlantDefinition, id)
    resolved = _as_plant_command_schema_lookup_id(id)
    for optic in plant.controllable_optics
        for schema in command_schemas(optic)
            command_schema_id(schema) == resolved && return schema
        end
    end
    throw(PlantDefinitionError(:command_schema, :unknown_id,
        "plant has no plant-command schema $resolved"))
end

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
AdaptiveOpticsSim.Plant.plant_model_definition_style(
    ::Type{MyModelDefinition},
) = AdaptiveOpticsSim.Plant.ColdPlantModelDefinition()
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
    ControllableOpticDefinition(id, optic_model, command_schemas;
        placement, visibility)
    ControllableOpticDefinition(id, optic_model;
        command_schemas, placement, visibility)

Immutable declaration of one physical controllable optic and the stable
semantic schemas of its independently timed or latched command endpoints. The
optic model is cold configuration and must opt in through
`plant_model_definition_style`. Placement and path visibility are required
optical-topology declarations. Mutable command state, timing, and prepared
optical grouping are deliberately absent from this topology record.
"""
struct ControllableOpticDefinition{
    M,
    S<:Tuple,
    P<:AbstractOpticalPlacement,
    V<:AbstractPathVisibility,
}
    id::ControllableOpticID
    optic_model::M
    command_schemas::S
    placement::P
    visibility::V

    function ControllableOpticDefinition(
        id::ControllableOpticID,
        optic_model::M,
        command_schemas::Tuple,
        placement::P,
        visibility::V,
    ) where {
        M,
        P<:AbstractOpticalPlacement,
        V<:AbstractPathVisibility,
    }
        _require_declared_topology_value(optic_model, :controllable_optic,
            :missing_model, "model")
        _require_cold_plant_model_definition(optic_model,
            :controllable_optic, "model")
        schemas = _normalize_command_schemas(command_schemas)
        _require_optical_placement(placement)
        _require_path_visibility(visibility)
        return new{M,typeof(schemas),P,V}(id, optic_model, schemas,
            placement, visibility)
    end
end

function ControllableOpticDefinition(id, optic_model, command_schemas;
    placement, visibility)
    schemas = _normalize_command_schemas(command_schemas)
    return ControllableOpticDefinition(_as_controllable_optic_id(id),
        optic_model, schemas,
        _require_optical_placement(placement),
        _require_path_visibility(visibility))
end

ControllableOpticDefinition(id, optic_model;
    command_schemas, placement, visibility) =
    ControllableOpticDefinition(id, optic_model, command_schemas;
        placement, visibility)

@inline _require_sampled_aberration_surface(
    surface::Union{NCPA,OPDMap}) = surface

function _require_sampled_aberration_surface(surface)
    throw(PlantDefinitionError(:sampled_aberration, :invalid_surface,
        "sampled aberrations must contain a native NCPA or OPDMap; got " *
        "$(typeof(surface))"))
end

@inline _require_sampled_aberration_application(
    application::Union{DMAdditive,DMReplace}) = application

function _require_sampled_aberration_application(application)
    throw(PlantDefinitionError(:sampled_aberration,
        :invalid_application,
        "sampled-aberration application must be DMAdditive or DMReplace; " *
        "got $(typeof(application))"))
end

@inline _require_sampled_aberration_placement(
    placement::Union{
        PupilPlanePlacement,
        AtmosphericConjugatePlacement,
    }) = placement

function _require_sampled_aberration_placement(placement)
    _require_optical_placement(placement, :sampled_aberration)
    throw(PlantDefinitionError(:sampled_aberration,
        :unsupported_placement,
        "native sampled OPD aberrations require PupilPlanePlacement or " *
        "AtmosphericConjugatePlacement; got $(typeof(placement))"))
end

@inline _require_sampled_aberration_registration(::Nothing) = nothing

function _require_sampled_aberration_registration(registration)
    throw(PlantDefinitionError(:sampled_aberration,
        :invalid_pupil_relay_registration,
        "sampled-aberration pupil-relay registration must be " *
        "PupilRelayRegistration or nothing; got $(typeof(registration))"))
end

function _require_sampled_aberration_metadata(
    metadata::OpticalPlaneMetadata,
    surface::Union{NCPA,OPDMap},
)
    opd = surface_opd(surface)
    typeof(metadata.kind) === PupilPlane || throw(PlantDefinitionError(
        :sampled_aberration, :invalid_surface_plane,
        "sampled-aberration metadata must describe PupilPlane; got " *
        "$(typeof(metadata.kind))"))
    typeof(metadata.coordinate_domain) === MetricCoordinates || throw(
        PlantDefinitionError(:sampled_aberration,
            :invalid_surface_coordinates,
            "sampled-aberration metadata must use MetricCoordinates; got " *
            "$(typeof(metadata.coordinate_domain))"))
    metadata.numeric_type <: AbstractFloat || throw(PlantDefinitionError(
        :sampled_aberration, :invalid_surface_numeric_type,
        "sampled-aberration metadata must use real floating-point samples; " *
        "got $(metadata.numeric_type)"))
    axes = metadata.orientation.axes
    (axes == (:x, :y) || axes == (:y, :x)) || throw(
        PlantDefinitionError(:sampled_aberration,
            :unsupported_axis_orientation,
            "sampled-aberration axes must be (:x, :y) or (:y, :x); got " *
            repr(axes)))
    size(opd) == metadata.dimensions || throw(PlantDefinitionError(
        :sampled_aberration, :surface_dimensions,
        "sampled-aberration OPD dimensions $(size(opd)) do not match " *
        "declared dimensions $(metadata.dimensions)"))
    eltype(opd) === metadata.numeric_type || throw(PlantDefinitionError(
        :sampled_aberration, :surface_numeric_type,
        "sampled-aberration OPD numeric type $(eltype(opd)) does not match " *
        "declared type $(metadata.numeric_type)"))
    typeof(backend(opd)) === typeof(metadata.backend) || throw(
        PlantDefinitionError(:sampled_aberration, :surface_backend,
            "sampled-aberration OPD and metadata use different array " *
            "backends"))
    compute_device(opd) == metadata.device || throw(PlantDefinitionError(
        :sampled_aberration, :surface_device,
        "sampled-aberration OPD and metadata occupy different compute " *
        "devices"))
    return metadata
end

function _require_sampled_aberration_metadata(metadata, surface)
    _require_sampled_aberration_surface(surface)
    throw(PlantDefinitionError(:sampled_aberration,
        :invalid_surface_metadata,
        "sampled-aberration metadata must be OpticalPlaneMetadata; got " *
        "$(typeof(metadata))"))
end

"""
    SampledAberrationDefinition(id, surface, metadata;
        placement, visibility, application, registration=nothing)

Immutable declaration of one native sampled pupil-OPD aberration. `surface`
must be `NCPA` or `OPDMap`; `metadata` gives its exact metric pupil-plane
sampling, array backend, and compute device. Placement, path visibility, and
additive or replacement application are explicit optical semantics.

The declaration retains caller storage. `prepare_plant` makes a defensive
backend-local OPD copy and prepares each visible path's finite-support
pupil-footprint coupling before repeated execution.
"""
struct SampledAberrationDefinition{
    S<:Union{NCPA,OPDMap},
    M<:OpticalPlaneMetadata,
    P<:AbstractOpticalPlacement,
    V<:AbstractPathVisibility,
    A<:DMApplyMode,
    R,
}
    id::SampledAberrationID
    surface::S
    metadata::M
    placement::P
    visibility::V
    application::A
    registration::R

    function SampledAberrationDefinition(
        id::SampledAberrationID,
        surface::S,
        metadata::M,
        placement::P,
        visibility::V,
        application::A,
        registration::R,
    ) where {
        S<:Union{NCPA,OPDMap},
        M<:OpticalPlaneMetadata,
        P<:AbstractOpticalPlacement,
        V<:AbstractPathVisibility,
        A<:DMApplyMode,
        R,
    }
        _require_sampled_aberration_surface(surface)
        _require_sampled_aberration_metadata(metadata, surface)
        _require_sampled_aberration_placement(placement)
        _require_path_visibility(visibility, :sampled_aberration)
        _require_sampled_aberration_application(application)
        _require_sampled_aberration_registration(registration)
        return new{S,M,P,V,A,R}(id, surface, metadata, placement,
            visibility, application, registration)
    end
end

function SampledAberrationDefinition(id, surface, metadata;
    placement, visibility, application, registration=nothing)
    validated_surface = _require_sampled_aberration_surface(surface)
    validated_metadata =
        _require_sampled_aberration_metadata(metadata, validated_surface)
    return SampledAberrationDefinition(
        _as_sampled_aberration_id(id),
        validated_surface,
        validated_metadata,
        _require_sampled_aberration_placement(placement),
        _require_path_visibility(visibility, :sampled_aberration),
        _require_sampled_aberration_application(application),
        _require_sampled_aberration_registration(registration),
    )
end

@inline path_id(path::OpticalPathDefinition) = path.id
@inline acquisition_id(acquisition::AcquisitionDefinition) = acquisition.id
@inline controllable_optic_id(optic::ControllableOpticDefinition) = optic.id
@inline sampled_aberration_id(
    aberration::SampledAberrationDefinition) = aberration.id
@inline acquisition_path_id(acquisition::AcquisitionDefinition) =
    acquisition.path
@inline path_source(path::OpticalPathDefinition) = path.source
@inline path_model(path::OpticalPathDefinition) = path.optical_model
@inline acquisition_model(acquisition::AcquisitionDefinition) =
    acquisition.acquisition_model
@inline controllable_optic_model(optic::ControllableOpticDefinition) =
    optic.optic_model
@inline controllable_optic_placement(optic::ControllableOpticDefinition) =
    optic.placement
@inline controllable_optic_visibility(optic::ControllableOpticDefinition) =
    optic.visibility
@inline sampled_aberration_surface(
    aberration::SampledAberrationDefinition) = aberration.surface
@inline sampled_aberration_metadata(
    aberration::SampledAberrationDefinition) = aberration.metadata
@inline sampled_aberration_placement(
    aberration::SampledAberrationDefinition) = aberration.placement
@inline sampled_aberration_visibility(
    aberration::SampledAberrationDefinition) = aberration.visibility
@inline sampled_aberration_application(
    aberration::SampledAberrationDefinition) = aberration.application
@inline sampled_aberration_registration(
    aberration::SampledAberrationDefinition) = aberration.registration
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

@inline _require_sampled_aberration_definition(
    ::SampledAberrationDefinition) = nothing

function _require_sampled_aberration_definition(value)
    throw(PlantDefinitionError(:sampled_aberration, :invalid_definition,
        "plant sampled aberrations must contain " *
        "SampledAberrationDefinition values; got $(typeof(value))"))
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

function _require_named_sampled_aberration_identity(
    key::Symbol,
    aberration::SampledAberrationDefinition,
)
    key == aberration.id.name || throw(PlantDefinitionError(
        :sampled_aberration, :identity_mismatch,
        "named sampled-aberration key $(repr(key)) does not match " *
        "$(aberration.id)"))
    return nothing
end

function _definition_memory(values, ::Type{T}, validator) where {T}
    memory = Memory{T}(undef, length(values))
    @inbounds for (index, value) in enumerate(values)
        validator(value)
        memory[index] = value
    end
    return memory
end

function _normalize_path_definitions(
    paths::Union{Tuple,AbstractVector},
)
    return _definition_memory(
        paths, OpticalPathDefinition, _require_path_definition)
end

function _normalize_path_definitions(paths::NamedTuple)
    normalized = _definition_memory(
        values(paths), OpticalPathDefinition, _require_path_definition)
    foreach(_require_named_path_identity, keys(paths), normalized)
    return normalized
end

function _normalize_path_definitions(paths)
    throw(PlantDefinitionError(:path, :invalid_container,
        "plant paths must be a Tuple, NamedTuple, or AbstractVector; got " *
        "$(typeof(paths))"))
end

function _normalize_acquisition_definitions(
    acquisitions::Union{Tuple,AbstractVector},
)
    return _definition_memory(acquisitions, AcquisitionDefinition,
        _require_acquisition_definition)
end

function _normalize_acquisition_definitions(acquisitions::NamedTuple)
    normalized = _definition_memory(values(acquisitions),
        AcquisitionDefinition, _require_acquisition_definition)
    foreach(_require_named_acquisition_identity, keys(acquisitions),
        normalized)
    return normalized
end

function _normalize_acquisition_definitions(acquisitions)
    throw(PlantDefinitionError(:acquisition, :invalid_container,
        "plant acquisitions must be a Tuple, NamedTuple, or AbstractVector; " *
        "got $(typeof(acquisitions))"))
end

function _normalize_controllable_optic_definitions(
    optics::Union{Tuple,AbstractVector},
)
    return _definition_memory(optics, ControllableOpticDefinition,
        _require_controllable_optic_definition)
end

function _normalize_controllable_optic_definitions(optics::NamedTuple)
    normalized = _definition_memory(values(optics),
        ControllableOpticDefinition,
        _require_controllable_optic_definition)
    foreach(_require_named_controllable_optic_identity, keys(optics),
        normalized)
    return normalized
end

function _normalize_controllable_optic_definitions(optics)
    throw(PlantDefinitionError(:controllable_optic, :invalid_container,
        "plant controllable optics must be a Tuple, NamedTuple, or " *
        "AbstractVector; got $(typeof(optics))"))
end

function _normalize_sampled_aberration_definitions(
    aberrations::Union{Tuple,AbstractVector},
)
    return _definition_memory(aberrations, SampledAberrationDefinition,
        _require_sampled_aberration_definition)
end

function _normalize_sampled_aberration_definitions(
    aberrations::NamedTuple)
    normalized = _definition_memory(values(aberrations),
        SampledAberrationDefinition,
        _require_sampled_aberration_definition)
    foreach(_require_named_sampled_aberration_identity, keys(aberrations),
        normalized)
    return normalized
end

function _normalize_sampled_aberration_definitions(aberrations)
    throw(PlantDefinitionError(:sampled_aberration, :invalid_container,
        "plant sampled aberrations must be a Tuple, NamedTuple, or " *
        "AbstractVector; got $(typeof(aberrations))"))
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

function _require_unique_path_ids(paths::AbstractVector)
    seen = Set{OpticalPathID}()
    for path in paths
        id = path_id(path)
        id in seen && throw(PlantDefinitionError(:path, :duplicate_id,
            "duplicate optical-path identity $id"))
        push!(seen, id)
    end
    return nothing
end

function _require_unique_acquisition_ids(acquisitions::AbstractVector)
    seen = Set{AcquisitionID}()
    for acquisition in acquisitions
        id = acquisition_id(acquisition)
        id in seen && throw(PlantDefinitionError(:acquisition, :duplicate_id,
            "duplicate acquisition identity $id"))
        push!(seen, id)
    end
    return nothing
end

function _require_unique_controllable_optic_ids(optics::AbstractVector)
    seen = Set{ControllableOpticID}()
    for optic in optics
        id = controllable_optic_id(optic)
        id in seen && throw(PlantDefinitionError(:controllable_optic,
            :duplicate_id, "duplicate controllable-optic identity $id"))
        push!(seen, id)
    end
    return nothing
end

function _require_unique_sampled_aberration_ids(
    aberrations::AbstractVector,
)
    seen = Set{SampledAberrationID}()
    for aberration in aberrations
        id = sampled_aberration_id(aberration)
        id in seen && throw(PlantDefinitionError(:sampled_aberration,
            :duplicate_id, "duplicate sampled-aberration identity $id"))
        push!(seen, id)
    end
    return nothing
end

function _require_unique_command_endpoint_owners(optics::AbstractVector)
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

function _require_unique_plant_command_schema_ids(optics::AbstractVector)
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

function _contains_path_id(paths::AbstractVector, id::OpticalPathID)
    for path in paths
        path_id(path) == id && return true
    end
    return false
end

function _require_acquisition_paths(
    paths::AbstractVector,
    acquisitions::AbstractVector,
)
    for acquisition in acquisitions
        id = acquisition_path_id(acquisition)
        _contains_path_id(paths, id) || throw(PlantDefinitionError(
            :acquisition, :unknown_path,
            "acquisition $(acquisition.id) references unknown path $id"))
    end
    return nothing
end

function _require_controllable_optic_visibility_paths(
    paths::AbstractVector, optics::AbstractVector)
    for optic in optics
        _require_controllable_optic_visibility_paths(paths, optic,
            controllable_optic_visibility(optic))
    end
    return nothing
end

@inline _require_controllable_optic_visibility_paths(
    ::AbstractVector,
    ::ControllableOpticDefinition,
    ::AllPathVisibility,
) = nothing

function _require_controllable_optic_visibility_paths(
    paths::AbstractVector,
    optic::ControllableOpticDefinition,
    visibility::SelectedPathVisibility,
)
    for id in selected_path_ids(visibility)
        _contains_path_id(paths, id) || throw(PlantDefinitionError(
            :controllable_optic, :unknown_visible_path,
            "controllable optic $(controllable_optic_id(optic)) " *
            "references unknown visible path $id"))
    end
    return nothing
end

function _require_sampled_aberration_visibility_paths(
    paths::AbstractVector, aberrations::AbstractVector)
    for aberration in aberrations
        _require_sampled_aberration_visibility_paths(paths, aberration,
            sampled_aberration_visibility(aberration))
    end
    return nothing
end

@inline _require_sampled_aberration_visibility_paths(
    ::AbstractVector,
    ::SampledAberrationDefinition,
    ::AllPathVisibility,
) = nothing

function _require_sampled_aberration_visibility_paths(
    paths::AbstractVector,
    aberration::SampledAberrationDefinition,
    visibility::SelectedPathVisibility,
)
    for id in selected_path_ids(visibility)
        _contains_path_id(paths, id) || throw(PlantDefinitionError(
            :sampled_aberration, :unknown_visible_path,
            "sampled aberration $(sampled_aberration_id(aberration)) " *
            "references unknown visible path $id"))
    end
    return nothing
end

"""
    PlantDefinition(telescope, atmosphere, controllable_optics,
        sampled_aberrations, paths, acquisitions)
    PlantDefinition(; telescope, atmosphere, controllable_optics=(),
        sampled_aberrations=(), paths=(), acquisitions=())

Immutable declared topology for one telescope and atmosphere, reusable optical
paths, independent acquisitions, and independently identified controllable
optics with versioned semantic command schemas. Native sampled aberrations are
separate immutable optical declarations rather than controllable devices.
Tuples, named tuples, and vectors are accepted as cold organization and copied
into fixed-size homogeneous registries; every component carries its own stable
identity. This value is not prepared
execution state and owns no mutable command state, schedule, queue, transport,
RNG stream, or HIL descriptor.
"""
struct PlantDefinition{T,A}
    telescope::T
    atmosphere::A
    controllable_optics::Memory{ControllableOpticDefinition}
    sampled_aberrations::Memory{SampledAberrationDefinition}
    paths::Memory{OpticalPathDefinition}
    acquisitions::Memory{AcquisitionDefinition}

    function PlantDefinition(telescope::T, atmosphere::A,
        controllable_optics::Memory{ControllableOpticDefinition},
        sampled_aberrations::Memory{SampledAberrationDefinition},
        paths::Memory{OpticalPathDefinition},
        acquisitions::Memory{AcquisitionDefinition},
    ) where {T,A}
        _require_plant_telescope(telescope)
        _require_plant_atmosphere(atmosphere)
        foreach(_require_controllable_optic_definition,
            controllable_optics)
        foreach(_require_sampled_aberration_definition,
            sampled_aberrations)
        foreach(_require_path_definition, paths)
        foreach(_require_acquisition_definition, acquisitions)
        _require_unique_controllable_optic_ids(controllable_optics)
        _require_unique_sampled_aberration_ids(sampled_aberrations)
        _require_unique_command_endpoint_owners(controllable_optics)
        _require_unique_plant_command_schema_ids(controllable_optics)
        _require_unique_path_ids(paths)
        _require_unique_acquisition_ids(acquisitions)
        _require_controllable_optic_visibility_paths(paths,
            controllable_optics)
        _require_sampled_aberration_visibility_paths(paths,
            sampled_aberrations)
        _require_acquisition_paths(paths, acquisitions)
        return new{T,A}(telescope, atmosphere,
            controllable_optics, sampled_aberrations, paths, acquisitions)
    end
end

function PlantDefinition(; telescope, atmosphere, controllable_optics=(),
    sampled_aberrations=(), paths=(), acquisitions=())
    normalized_optics =
        _normalize_controllable_optic_definitions(controllable_optics)
    normalized_aberrations =
        _normalize_sampled_aberration_definitions(sampled_aberrations)
    normalized_paths = _normalize_path_definitions(paths)
    normalized_acquisitions = _normalize_acquisition_definitions(acquisitions)
    return PlantDefinition(telescope, atmosphere, normalized_optics,
        normalized_aberrations, normalized_paths, normalized_acquisitions)
end

function PlantDefinition(telescope, atmosphere, controllable_optics,
    sampled_aberrations, paths, acquisitions)
    normalized_optics =
        _normalize_controllable_optic_definitions(controllable_optics)
    normalized_aberrations =
        _normalize_sampled_aberration_definitions(sampled_aberrations)
    normalized_paths = _normalize_path_definitions(paths)
    normalized_acquisitions = _normalize_acquisition_definitions(acquisitions)
    return PlantDefinition(telescope, atmosphere, normalized_optics,
        normalized_aberrations, normalized_paths, normalized_acquisitions)
end

function PlantDefinition(telescope, atmosphere, controllable_optics, paths,
    acquisitions)
    return PlantDefinition(telescope, atmosphere, controllable_optics, (),
        paths, acquisitions)
end

@inline plant_telescope(plant::PlantDefinition) = plant.telescope
@inline plant_atmosphere(plant::PlantDefinition) = plant.atmosphere
@inline controllable_optic_definitions(plant::PlantDefinition) =
    plant.controllable_optics
@inline sampled_aberration_definitions(plant::PlantDefinition) =
    plant.sampled_aberrations
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

function sampled_aberration_definition(plant::PlantDefinition, id)
    resolved = _as_sampled_aberration_id(id)
    for aberration in plant.sampled_aberrations
        sampled_aberration_id(aberration) == resolved && return aberration
    end
    throw(PlantDefinitionError(:sampled_aberration, :unknown_id,
        "plant has no sampled aberration $resolved"))
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

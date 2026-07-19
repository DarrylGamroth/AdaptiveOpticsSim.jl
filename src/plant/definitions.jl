#
# Immutable plant topology declarations
#
# These values describe stable physical identity and topology only. Backend-
# bound preparation, mutable state/workspace ownership, execution, and timing
# belong to later layers.
#

@inline function _require_component_name(name::Symbol, component::Symbol)
    isempty(String(name)) && throw(PlantDefinitionError(component, :empty_id,
        "$component identity must not be empty"))
    return name
end

"""Stable declared identity of one reusable optical path."""
struct OpticalPathID
    name::Symbol

    function OpticalPathID(name::Symbol)
        return new(_require_component_name(name, :path))
    end
end

"""Stable declared identity of one independently invocable acquisition."""
struct AcquisitionID
    name::Symbol

    function AcquisitionID(name::Symbol)
        return new(_require_component_name(name, :acquisition))
    end
end

Base.:(==)(left::OpticalPathID, right::OpticalPathID) =
    left.name == right.name
Base.:(==)(left::AcquisitionID, right::AcquisitionID) =
    left.name == right.name
Base.isequal(left::OpticalPathID, right::OpticalPathID) =
    isequal(left.name, right.name)
Base.isequal(left::AcquisitionID, right::AcquisitionID) =
    isequal(left.name, right.name)
Base.hash(id::OpticalPathID, seed::UInt) =
    hash(id.name, hash(OpticalPathID, seed))
Base.hash(id::AcquisitionID, seed::UInt) =
    hash(id.name, hash(AcquisitionID, seed))

function Base.show(io::IO, id::Union{OpticalPathID,AcquisitionID})
    print(io, nameof(typeof(id)), "(", repr(id.name), ")")
end

@inline _as_optical_path_id(id::OpticalPathID) = id
@inline _as_optical_path_id(name::Symbol) = OpticalPathID(name)

function _as_optical_path_id(value)
    throw(PlantDefinitionError(:path, :invalid_id,
        "optical-path identity must be a Symbol or OpticalPathID; got " *
        "$(typeof(value))"))
end

@inline _as_acquisition_id(id::AcquisitionID) = id
@inline _as_acquisition_id(name::Symbol) = AcquisitionID(name)

function _as_acquisition_id(value)
    throw(PlantDefinitionError(:acquisition, :invalid_id,
        "acquisition identity must be a Symbol or AcquisitionID; got " *
        "$(typeof(value))"))
end

@inline _require_declared_topology_value(value, ::Symbol, ::Symbol,
    ::AbstractString) = value

function _require_declared_topology_value(::Nothing, component::Symbol,
    reason::Symbol, label::AbstractString)
    throw(PlantDefinitionError(component, reason,
        "$component $label must be declared"))
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
workspace or mutable acquisition state.
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
state, timing, RNG, and publication ownership are deliberately absent.
"""
struct AcquisitionDefinition{M}
    id::AcquisitionID
    path::OpticalPathID
    acquisition_model::M

    function AcquisitionDefinition(id::AcquisitionID, path::OpticalPathID,
        acquisition_model::M) where {M}
        _require_declared_topology_value(acquisition_model, :acquisition,
            :missing_model, "model")
        return new{M}(id, path, acquisition_model)
    end
end

function AcquisitionDefinition(id, path, acquisition_model)
    return AcquisitionDefinition(_as_acquisition_id(id),
        _as_optical_path_id(path), acquisition_model)
end

AcquisitionDefinition(id, path; acquisition_model) =
    AcquisitionDefinition(id, path, acquisition_model)

@inline path_id(path::OpticalPathDefinition) = path.id
@inline acquisition_id(acquisition::AcquisitionDefinition) = acquisition.id
@inline acquisition_path_id(acquisition::AcquisitionDefinition) =
    acquisition.path
@inline path_source(path::OpticalPathDefinition) = path.source
@inline path_model(path::OpticalPathDefinition) = path.optical_model
@inline acquisition_model(acquisition::AcquisitionDefinition) =
    acquisition.acquisition_model

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
    PlantDefinition(telescope, atmosphere, paths, acquisitions)
    PlantDefinition(; telescope, atmosphere, paths=(), acquisitions=())

Immutable declared topology for one telescope and atmosphere, reusable optical
paths, and independent acquisitions. Tuples and named tuples are accepted as
cold organization only; every path and acquisition carries its own stable
identity. This value is not prepared execution state and owns no schedule,
queue, transport, RNG stream, or HIL descriptor.
"""
struct PlantDefinition{T,A,P<:Tuple,Q<:Tuple}
    telescope::T
    atmosphere::A
    paths::P
    acquisitions::Q

    function PlantDefinition(telescope::T, atmosphere::A, paths::P,
        acquisitions::Q) where {T,A,P<:Tuple,Q<:Tuple}
        _require_plant_telescope(telescope)
        _require_plant_atmosphere(atmosphere)
        foreach(_require_path_definition, paths)
        foreach(_require_acquisition_definition, acquisitions)
        _require_unique_path_ids(paths)
        _require_unique_acquisition_ids(acquisitions)
        _require_acquisition_paths(paths, acquisitions)
        return new{T,A,P,Q}(telescope, atmosphere, paths, acquisitions)
    end
end

function PlantDefinition(; telescope, atmosphere, paths=(), acquisitions=())
    normalized_paths = _normalize_path_definitions(paths)
    normalized_acquisitions = _normalize_acquisition_definitions(acquisitions)
    return PlantDefinition(telescope, atmosphere, normalized_paths,
        normalized_acquisitions)
end

function PlantDefinition(telescope, atmosphere, paths, acquisitions)
    normalized_paths = _normalize_path_definitions(paths)
    normalized_acquisitions = _normalize_acquisition_definitions(acquisitions)
    return PlantDefinition(telescope, atmosphere, normalized_paths,
        normalized_acquisitions)
end

@inline plant_telescope(plant::PlantDefinition) = plant.telescope
@inline plant_atmosphere(plant::PlantDefinition) = plant.atmosphere
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

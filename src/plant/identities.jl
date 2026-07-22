#
# Stable plant identities
#
# These compact value types are declaration identities. Hashing supports cold
# lookup; a hash value never defines or serializes an identity.
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

"""Stable declared identity of one physical controllable optic."""
struct ControllableOpticID
    name::Symbol

    function ControllableOpticID(name::Symbol)
        return new(_require_component_name(name, :controllable_optic))
    end
end

"""Stable declared identity of one independently timed command endpoint."""
struct CommandEndpointID
    name::Symbol

    function CommandEndpointID(name::Symbol)
        return new(_require_component_name(name, :command_endpoint))
    end
end

const _PlantTopologyID = Union{
    OpticalPathID,
    AcquisitionID,
    ControllableOpticID,
    CommandEndpointID,
}

Base.:(==)(left::T, right::T) where {T<:_PlantTopologyID} =
    left.name == right.name
Base.isequal(left::T, right::T) where {T<:_PlantTopologyID} =
    isequal(left.name, right.name)
Base.hash(id::T, seed::UInt) where {T<:_PlantTopologyID} =
    hash(id.name, hash(T, seed))

function Base.show(io::IO, id::_PlantTopologyID)
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

@inline _as_controllable_optic_id(id::ControllableOpticID) = id
@inline _as_controllable_optic_id(name::Symbol) = ControllableOpticID(name)

function _as_controllable_optic_id(value)
    throw(PlantDefinitionError(:controllable_optic, :invalid_id,
        "controllable-optic identity must be a Symbol or " *
        "ControllableOpticID; got $(typeof(value))"))
end

@inline _as_command_endpoint_id(id::CommandEndpointID) = id
@inline _as_command_endpoint_id(name::Symbol) = CommandEndpointID(name)

function _as_command_endpoint_id(value)
    throw(PlantDefinitionError(:command_endpoint, :invalid_id,
        "command-endpoint identity must be a Symbol or CommandEndpointID; " *
        "got $(typeof(value))"))
end

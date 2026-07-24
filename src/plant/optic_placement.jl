#
# Controllable-optic optical placement and path visibility
#
# These immutable declarations describe where an optic acts and which optical
# paths contain it. They do not describe command timing, controller routing,
# transport, or execution-resource placement.
#

"""Optical location declared for one physical controllable optic."""
abstract type AbstractControllableOpticPlacement end

"""A controllable optic located in a pupil plane."""
struct PupilPlanePlacement <: AbstractControllableOpticPlacement end

"""
    AtmosphericConjugatePlacement(altitude_m)

A geometric controllable-optic plane conjugate to a finite, nonnegative
altitude above the entrance pupil, expressed in metres.
"""
struct AtmosphericConjugatePlacement{T<:AbstractFloat} <:
       AbstractControllableOpticPlacement
    altitude_m::T

    function AtmosphericConjugatePlacement(altitude_m::T) where {
        T<:AbstractFloat,
    }
        isfinite(altitude_m) || throw(PlantDefinitionError(
            :controllable_optic, :invalid_conjugate_altitude,
            "atmospheric conjugate altitude must be finite; got " *
            repr(altitude_m)))
        altitude_m >= zero(T) || throw(PlantDefinitionError(
            :controllable_optic, :invalid_conjugate_altitude,
            "atmospheric conjugate altitude must be nonnegative; got " *
            repr(altitude_m)))
        normalized = iszero(altitude_m) ? zero(T) : altitude_m
        return new{T}(normalized)
    end
end

AtmosphericConjugatePlacement(altitude_m::Real) =
    AtmosphericConjugatePlacement(float(altitude_m))

Base.:(==)(
    left::AtmosphericConjugatePlacement,
    right::AtmosphericConjugatePlacement,
) = left.altitude_m == right.altitude_m
Base.isequal(
    left::AtmosphericConjugatePlacement,
    right::AtmosphericConjugatePlacement,
) = isequal(left.altitude_m, right.altitude_m)
Base.hash(placement::AtmosphericConjugatePlacement, seed::UInt) =
    hash(placement.altitude_m,
        hash(AtmosphericConjugatePlacement, seed))

"""A controllable optic located in a physical focus or focal-mask plane."""
struct FocalPlanePlacement <: AbstractControllableOpticPlacement end

"""Path-selection declaration for one physical controllable optic."""
abstract type AbstractControllableOpticVisibility end

"""The optic is present in every optical path declared by its plant."""
struct AllPathVisibility <: AbstractControllableOpticVisibility end

"""
    SelectedPathVisibility(paths...)
    SelectedPathVisibility(paths)

The optic is present only in the nonempty, duplicate-free set of explicitly
identified optical paths. Input containers are normalized to an immutable
tuple of `OpticalPathID` values in stable identity order.
"""
struct SelectedPathVisibility{P<:Tuple} <:
       AbstractControllableOpticVisibility
    paths::P

    function SelectedPathVisibility(paths::Tuple)
        converted = map(_as_optical_path_id, paths)
        isempty(converted) && throw(PlantDefinitionError(
            :controllable_optic, :empty_path_visibility,
            "selected-path visibility must name at least one optical path"))
        normalized = _canonical_selected_path_ids(converted)
        @inbounds for index in 2:length(normalized)
            normalized[index - 1] == normalized[index] && throw(
                PlantDefinitionError(
                    :controllable_optic, :duplicate_visible_path,
                    "selected-path visibility names " *
                    "$(normalized[index]) more than once"))
        end
        return new{typeof(normalized)}(normalized)
    end
end

function _canonical_selected_path_ids(
    paths::NTuple{N,OpticalPathID}) where {N}
    ordered = collect(paths)
    sort!(ordered; by=id -> String(id.name))
    return ntuple(index -> @inbounds(ordered[index]), Val(N))
end

SelectedPathVisibility(paths::AbstractVector) =
    SelectedPathVisibility(Tuple(paths))
SelectedPathVisibility(paths...) = SelectedPathVisibility(paths)

@inline conjugate_altitude_m(
    placement::AtmosphericConjugatePlacement) = placement.altitude_m
@inline selected_path_ids(visibility::SelectedPathVisibility) =
    visibility.paths

@inline _require_controllable_optic_placement(
    placement::Union{
        PupilPlanePlacement,
        AtmosphericConjugatePlacement,
        FocalPlanePlacement,
    },
) = placement

function _require_controllable_optic_placement(placement)
    throw(PlantDefinitionError(:controllable_optic, :invalid_placement,
        "controllable-optic placement must be PupilPlanePlacement, " *
        "AtmosphericConjugatePlacement, or FocalPlanePlacement; got " *
        "$(typeof(placement))"))
end

@inline _require_controllable_optic_visibility(
    visibility::Union{AllPathVisibility,SelectedPathVisibility}) =
    visibility

function _require_controllable_optic_visibility(visibility)
    throw(PlantDefinitionError(:controllable_optic, :invalid_path_visibility,
        "controllable-optic visibility must be AllPathVisibility or " *
        "SelectedPathVisibility; got $(typeof(visibility))"))
end

@inline _optic_visible_on_path(
    ::AllPathVisibility, ::OpticalPathID) = true

function _optic_visible_on_path(
    visibility::SelectedPathVisibility, path::OpticalPathID)
    @inbounds for selected in visibility.paths
        selected == path && return true
    end
    return false
end

@inline _same_optic_placement(
    ::PupilPlanePlacement, ::PupilPlanePlacement) = true
@inline _same_optic_placement(
    ::FocalPlanePlacement, ::FocalPlanePlacement) = true
@inline _same_optic_placement(
    left::AtmosphericConjugatePlacement,
    right::AtmosphericConjugatePlacement,
) = left.altitude_m == right.altitude_m
@inline _same_optic_placement(
    ::AbstractControllableOpticPlacement,
    ::AbstractControllableOpticPlacement,
) = false

@inline _optic_placement_rank(::PupilPlanePlacement) = UInt8(0)
@inline _optic_placement_rank(::AtmosphericConjugatePlacement) = UInt8(1)
@inline _optic_placement_rank(::FocalPlanePlacement) = UInt8(2)

@inline function _optic_placement_isless(
    left::AbstractControllableOpticPlacement,
    right::AbstractControllableOpticPlacement,
)
    left_rank = _optic_placement_rank(left)
    right_rank = _optic_placement_rank(right)
    left_rank == right_rank || return left_rank < right_rank
    return _optic_placement_isless_same_rank(left, right)
end

@inline _optic_placement_isless_same_rank(
    ::PupilPlanePlacement, ::PupilPlanePlacement) = false
@inline _optic_placement_isless_same_rank(
    ::FocalPlanePlacement, ::FocalPlanePlacement) = false
@inline _optic_placement_isless_same_rank(
    left::AtmosphericConjugatePlacement,
    right::AtmosphericConjugatePlacement,
) = left.altitude_m < right.altitude_m

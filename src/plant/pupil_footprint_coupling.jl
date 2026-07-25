#
# Prepared geometric coupling of sampled OPD surfaces to path-local pupil
# functions.
#
# The cold PupilRelayRegistration describes path/relay geometry. Prepared
# couplings compose it with source direction, finite-height cone geometry, and
# both sampled-plane coordinate contracts. Device-internal actuator
# misregistration has already acted while forming the sampled surface and is
# deliberately absent here.
#

"""Path-local coupling prepared for one sampled pupil-plane OPD surface."""
abstract type AbstractPupilSurfacePathCoupling end

"""Internal marker for an optic executed through another exact path coupling."""
struct _NoPupilSurfacePathCoupling <:
       AbstractPupilSurfacePathCoupling end

struct _PupilRelayRegistrationToken end
const _PUPIL_RELAY_REGISTRATION_TOKEN =
    _PupilRelayRegistrationToken()

"""
    PupilRelayRegistration(; magnification=(1, 1), rotation_deg=0,
        parity=(1, 1), decenter_m=(0, 0), T=Float64)

Cold affine registration from geometric source-footprint coordinates to a
sampled pupil-surface grid. Positive per-axis pupil magnification and
parity are applied first, followed by counterclockwise rotation and then a
metric decenter.

This is path/relay geometry. `Misregistration` remains the separate
device-internal actuator-to-surface convention.
"""
struct PupilRelayRegistration{T<:AbstractFloat}
    magnification::NTuple{2,T}
    rotation_rad::T
    parity::NTuple{2,Int8}
    decenter_m::NTuple{2,T}
    transform::NTuple{4,T}

    function PupilRelayRegistration(
        ::_PupilRelayRegistrationToken,
        magnification::NTuple{2,T},
        rotation_rad::T,
        parity::NTuple{2,Int8},
        decenter_m::NTuple{2,T},
        transform::NTuple{4,T},
    ) where {T<:AbstractFloat}
        return new{T}(
            magnification,
            rotation_rad,
            parity,
            decenter_m,
            transform,
        )
    end
end

@inline _pupil_relay_pair(value::Real, ::Type{T}) where {T<:AbstractFloat} =
    (T(value), T(value))

@inline function _pupil_relay_pair(
    values::NTuple{2,<:Real}, ::Type{T}) where {T<:AbstractFloat}
    return (T(values[1]), T(values[2]))
end

function _pupil_relay_pair(values, ::Type{T}) where {T<:AbstractFloat}
    throw(PlantDefinitionError(:pupil_relay_registration,
        :invalid_pupil_relay_registration,
        "pupil-relay pair values must be a real scalar or a two-value " *
        "real tuple; got $(typeof(values))"))
end

@inline function _pupil_relay_parity(
    values::NTuple{2,<:Integer})
    all(value -> value == -1 || value == 1, values) || throw(
        PlantDefinitionError(:pupil_relay_registration,
            :invalid_pupil_relay_registration,
            "pupil-relay parity values must be -1 or 1; got " *
            repr(values)))
    converted = (Int8(values[1]), Int8(values[2]))
    return converted
end

function _pupil_relay_parity(values)
    throw(PlantDefinitionError(:pupil_relay_registration,
        :invalid_pupil_relay_registration,
        "pupil-relay parity must be a two-value integer tuple; got " *
        "$(typeof(values))"))
end

function PupilRelayRegistration(;
    magnification=(1.0, 1.0),
    rotation_deg::Real=0.0,
    parity=(1, 1),
    decenter_m=(0.0, 0.0),
    T::Type{<:AbstractFloat}=Float64,
)
    isconcretetype(T) || throw(PlantDefinitionError(
        :pupil_relay_registration,
        :invalid_pupil_relay_registration,
        "pupil-relay numeric type must be concrete; got $(T)"))
    magnification_t = _pupil_relay_pair(magnification, T)
    all(value -> isfinite(value) && value > zero(T), magnification_t) ||
        throw(PlantDefinitionError(:pupil_relay_registration,
            :invalid_pupil_relay_registration,
            "pupil-relay magnification must be finite and positive; got " *
            repr(magnification)))
    rotation_rad = T(deg2rad(rotation_deg))
    isfinite(rotation_rad) || throw(PlantDefinitionError(
        :pupil_relay_registration, :invalid_pupil_relay_registration,
        "pupil-relay rotation must be finite; got " *
        repr(rotation_deg)))
    parity_t = _pupil_relay_parity(parity)
    decenter_t = _pupil_relay_pair(decenter_m, T)
    all(isfinite, decenter_t) || throw(PlantDefinitionError(
        :pupil_relay_registration, :invalid_pupil_relay_registration,
        "pupil-relay decenter must be finite; got " *
        repr(decenter_m)))
    sine, cosine = sincos(rotation_rad)
    scale_x = T(parity_t[1]) * magnification_t[1]
    scale_y = T(parity_t[2]) * magnification_t[2]
    transform = (
        cosine * scale_x,
        -sine * scale_y,
        sine * scale_x,
        cosine * scale_y,
    )
    return PupilRelayRegistration(
        _PUPIL_RELAY_REGISTRATION_TOKEN,
        magnification_t,
        rotation_rad,
        parity_t,
        decenter_t,
        transform,
    )
end

@inline pupil_relay_magnification(
    registration::PupilRelayRegistration) = registration.magnification
@inline pupil_relay_rotation_rad(
    registration::PupilRelayRegistration) = registration.rotation_rad
@inline pupil_relay_rotation_deg(
    registration::PupilRelayRegistration) =
    rad2deg(registration.rotation_rad)
@inline pupil_relay_parity(
    registration::PupilRelayRegistration) = registration.parity
@inline pupil_relay_decenter_m(
    registration::PupilRelayRegistration) = registration.decenter_m

@inline _convert_pupil_relay_registration(
    registration::PupilRelayRegistration{T},
    ::Type{T},
) where {T<:AbstractFloat} = registration

@inline function _convert_pupil_relay_registration(
    registration::PupilRelayRegistration,
    ::Type{T},
) where {T<:AbstractFloat}
    magnification = ntuple(
        index -> T(registration.magnification[index]), Val(2))
    all(value -> isfinite(value) && value > zero(T), magnification) ||
        throw(PlantPreparationError(:pupil_relay_registration,
            :invalid_pupil_relay_registration,
            "pupil-relay magnification cannot be represented as finite " *
            "positive $(T) values"))
    rotation_rad = T(registration.rotation_rad)
    isfinite(rotation_rad) || throw(PlantPreparationError(
        :pupil_relay_registration,
        :invalid_pupil_relay_registration,
        "pupil-relay rotation cannot be represented finitely as $(T)"))
    decenter_m = ntuple(
        index -> T(registration.decenter_m[index]), Val(2))
    all(isfinite, decenter_m) || throw(PlantPreparationError(
        :pupil_relay_registration,
        :invalid_pupil_relay_registration,
        "pupil-relay decenter cannot be represented finitely as $(T)"))
    sine, cosine = sincos(rotation_rad)
    scale_x = T(registration.parity[1]) * magnification[1]
    scale_y = T(registration.parity[2]) * magnification[2]
    transform = (
        cosine * scale_x,
        -sine * scale_y,
        sine * scale_x,
        cosine * scale_y,
    )
    return PupilRelayRegistration(
        _PUPIL_RELAY_REGISTRATION_TOKEN,
        magnification,
        rotation_rad,
        registration.parity,
        decenter_m,
        transform,
    )
end

@inline _resolve_pupil_relay_registration(
    ::Nothing, ::Type{T}) where {T<:AbstractFloat} =
    PupilRelayRegistration(T=T)

@inline _resolve_pupil_relay_registration(
    registration::PupilRelayRegistration,
    ::Type{T},
) where {T<:AbstractFloat} =
    _convert_pupil_relay_registration(registration, T)

function _resolve_pupil_relay_registration(
    registration, ::Type{<:AbstractFloat})
    throw(PlantPreparationError(:pupil_relay_registration,
        :invalid_pupil_relay_registration,
        "pupil-relay registration must be PupilRelayRegistration or " *
        "nothing; got $(typeof(registration))"))
end

@inline function _is_identity_pupil_relay_registration(
    registration::PupilRelayRegistration{T}) where {T}
    return registration.magnification == (one(T), one(T)) &&
        iszero(registration.rotation_rad) &&
        registration.parity == (Int8(1), Int8(1)) &&
        registration.decenter_m == (zero(T), zero(T))
end

"""
Default path-local coupling for a pupil-plane model that applies its own
already same-grid surface semantics.
"""
struct _PreparedPupilFootprintCouplingToken end
const _PREPARED_PUPIL_FOOTPRINT_COUPLING_TOKEN =
    _PreparedPupilFootprintCouplingToken()

struct PreparedDirectPupilSurfaceCoupling{P<:PupilFunction} <:
       AbstractPupilSurfacePathCoupling
    destination::P

    function PreparedDirectPupilSurfaceCoupling(
        ::_PreparedPupilFootprintCouplingToken,
        destination::P,
    ) where {P<:PupilFunction}
        return new{P}(destination)
    end
end

"""
Same-grid sampled-surface fast path bound to one exact pupil destination and
one immutable surface-grid contract.
"""
struct PreparedIdentityPupilFootprintCoupling{
    P<:PupilFunction,
    M<:OpticalPlaneMetadata,
} <:
       AbstractPupilSurfacePathCoupling
    destination::P
    surface_metadata::M

    function PreparedIdentityPupilFootprintCoupling(
        ::_PreparedPupilFootprintCouplingToken,
        destination::P,
        surface_metadata::M,
    ) where {P<:PupilFunction,M<:OpticalPlaneMetadata}
        return new{P,M}(destination, surface_metadata)
    end
end

"""
Finite-support affine mapping from one exact pupil destination's sample indices
to a sampled pupil-surface grid.
"""
struct PreparedPupilFootprintCoupling{
    T<:AbstractFloat,
    P<:PupilFunction,
    M<:OpticalPlaneMetadata,
} <: AbstractPupilSurfacePathCoupling
    destination::P
    surface_metadata::M
    index_transform::NTuple{4,T}
    index_offset::NTuple{2,T}

    function PreparedPupilFootprintCoupling(
        ::_PreparedPupilFootprintCouplingToken,
        destination::P,
        surface_metadata::M,
        index_transform::NTuple{4,T},
        index_offset::NTuple{2,T},
    ) where {
        T<:AbstractFloat,
        P<:PupilFunction,
        M<:OpticalPlaneMetadata,
    }
        return new{T,P,M}(
            destination,
            surface_metadata,
            index_transform,
            index_offset,
        )
    end
end

@inline function _require_cartesian_plane_orientation(
    orientation, label::AbstractString)
    axes = orientation.axes
    (axes == (:x, :y) || axes == (:y, :x)) || throw(
        PlantPreparationError(:pupil_surface_coupling,
            :unsupported_axis_orientation,
            "$label axis orientation must be (:x, :y) or (:y, :x); " *
            "got $(axes)"))
    return orientation
end

function _require_metric_pupil_surface_metadata(
    metadata::OpticalPlaneMetadata, label::AbstractString)
    typeof(metadata.kind) === PupilPlane || throw(
        PlantPreparationError(:pupil_surface_coupling,
            :invalid_surface_plane,
            "$label must be declared on PupilPlane; got " *
            "$(typeof(metadata.kind))"))
    typeof(metadata.coordinate_domain) === MetricCoordinates || throw(
        PlantPreparationError(:pupil_surface_coupling,
            :invalid_surface_coordinates,
            "$label must use MetricCoordinates; got " *
            "$(typeof(metadata.coordinate_domain))"))
    metadata.numeric_type <: AbstractFloat || throw(
        PlantPreparationError(:pupil_surface_coupling,
            :invalid_surface_numeric_type,
            "$label must use real floating-point samples; got " *
            "$(metadata.numeric_type)"))
    _require_cartesian_plane_orientation(metadata.orientation, label)
    return metadata
end

@inline function _require_sampled_surface_domain(
    destination::PupilFunction,
    surface_metadata::OpticalPlaneMetadata,
)
    destination_metadata = _require_metric_pupil_surface_metadata(
        destination.metadata, "path pupil destination")
    _require_metric_pupil_surface_metadata(
        surface_metadata, "sampled pupil surface")
    surface_metadata.numeric_type === eltype(destination.opd) || throw(
        PlantPreparationError(:pupil_surface_coupling,
            :surface_numeric_type,
            "sampled pupil-surface numeric type " *
            "$(surface_metadata.numeric_type) does not match path OPD " *
            "$(eltype(destination.opd))"))
    typeof(surface_metadata.backend) ===
        typeof(destination_metadata.backend) || throw(
        PlantPreparationError(:pupil_surface_coupling, :surface_backend,
            "sampled pupil surface and path pupil use different " *
            "array backends"))
    surface_metadata.device == destination_metadata.device || throw(
        PlantPreparationError(:pupil_surface_coupling, :surface_device,
            "sampled pupil surface and path pupil occupy different " *
            "physical devices"))
    return destination_metadata
end

function _require_sampled_surface_domain(
    destination,
    ::OpticalPlaneMetadata,
)
    throw(PlantPreparationError(:pupil_surface_coupling,
        :unsupported_path_input,
        "sampled pupil-footprint coupling requires a PupilFunction path " *
        "input; got $(typeof(destination))"))
end

@inline function _same_pupil_sampling_grid(
    left::OpticalPlaneMetadata,
    right::OpticalPlaneMetadata,
)
    return left.dimensions == right.dimensions &&
        left.sampling == right.sampling &&
        left.origin == right.origin &&
        left.centering == right.centering &&
        left.orientation == right.orientation &&
        left.numeric_type === right.numeric_type &&
        typeof(left.backend) === typeof(right.backend) &&
        left.device == right.device
end

@inline function _pupil_plane_footprint_geometry(
    ::PupilPlanePlacement,
    ::AbstractSource,
    ::Type{T},
) where {T<:AbstractFloat}
    return zero(T), zero(T), one(T)
end

function _pupil_plane_footprint_geometry(
    placement::AtmosphericConjugatePlacement,
    source::AbstractSource,
    ::Type{T},
) where {T<:AbstractFloat}
    return _pupil_plane_footprint_geometry(
        source_composition_style(source),
        placement,
        source,
        T,
    )
end

@inline function _pupil_plane_footprint_geometry(
    ::LeafSourceComposition,
    placement::AtmosphericConjugatePlacement,
    source::AbstractSource,
    ::Type{T},
) where {T<:AbstractFloat}
    return _single_pupil_footprint_geometry(placement, source, T)
end

@inline function _pupil_plane_footprint_geometry(
    ::ExpandedSourceComposition,
    placement::AtmosphericConjugatePlacement,
    source::SpectralSource,
    ::Type{T},
) where {T<:AbstractFloat}
    return _single_pupil_footprint_geometry(placement, source, T)
end

function _pupil_plane_footprint_geometry(
    ::ExpandedSourceComposition,
    ::AtmosphericConjugatePlacement,
    source::AbstractSource,
    ::Type{<:AbstractFloat},
)
    throw(PlantPreparationError(:pupil_surface_coupling,
        :unsupported_source_geometry,
        "atmospheric-conjugate pupil-footprint coupling requires one " *
        "source direction; $(typeof(source)) expands into multiple " *
        "directions that must be prepared as separate paths"))
end

function _pupil_plane_footprint_geometry(
    style::AbstractSourceCompositionStyle,
    ::AtmosphericConjugatePlacement,
    source::AbstractSource,
    ::Type{<:AbstractFloat},
)
    throw(PlantPreparationError(:pupil_surface_coupling,
        :unsupported_source_geometry,
        "atmospheric-conjugate pupil-footprint coupling does not support " *
        "source composition style $(typeof(style)) for $(typeof(source))"))
end

function _single_pupil_footprint_geometry(
    placement::AtmosphericConjugatePlacement,
    source::AbstractSource,
    ::Type{T},
) where {T<:AbstractFloat}
    altitude_m = T(conjugate_altitude_m(placement))
    isfinite(altitude_m) || throw(PlantPreparationError(
        :pupil_surface_coupling, :invalid_conjugate_altitude,
        "optical conjugate altitude cannot be represented " *
        "finitely as $(T)"))
    coordinates_arcsec = coordinates_xy_arcsec(source)
    coordinates_t = (
        T(coordinates_arcsec[1]),
        T(coordinates_arcsec[2]),
    )
    all(isfinite, coordinates_t) || throw(PlantPreparationError(
        :pupil_surface_coupling, :invalid_source_geometry,
        "source angular coordinates must be finite in $(T); got " *
        repr(coordinates_arcsec)))
    source_height = source_height_m(source)
    (isfinite(source_height) || source_height == Inf) || throw(
        PlantPreparationError(:pupil_surface_coupling,
            :invalid_source_geometry,
            "source height must be finite or positive infinity; got " *
            repr(source_height)))
    height_m = T(source_height)
    (isfinite(height_m) || source_height == Inf) || throw(
        PlantPreparationError(:pupil_surface_coupling,
            :invalid_source_geometry,
            "finite source height cannot be represented finitely as $(T); " *
            "got $(source_height)"))
    height_m > zero(T) || throw(PlantPreparationError(
        :pupil_surface_coupling, :invalid_source_geometry,
        "source height must be positive; got " *
        repr(source_height)))
    footprint_scale = one(T)
    if isfinite(height_m)
        height_m > altitude_m || throw(PlantPreparationError(
            :pupil_surface_coupling, :source_below_conjugate,
            "finite source height $(height_m) m must exceed optical " *
            "conjugate altitude $(altitude_m) m"))
        footprint_scale = (height_m - altitude_m) / height_m
    end
    arcsec_to_rad = T(pi / (180 * 3600))
    shift_x_m = coordinates_t[1] * arcsec_to_rad * altitude_m
    shift_y_m = coordinates_t[2] * arcsec_to_rad * altitude_m
    return shift_x_m, shift_y_m, footprint_scale
end

function _pupil_plane_footprint_geometry(
    placement::AbstractOpticalPlacement,
    ::AbstractSource,
    ::Type{<:AbstractFloat},
)
    throw(PlantPreparationError(:pupil_surface_coupling,
        :unsupported_geometric_placement,
        "sampled pupil-footprint coupling does not support " *
        "$(typeof(placement))"))
end

@inline function _logical_plane_coordinates(
    metadata::OpticalPlaneMetadata,
    first_index,
    second_index,
)
    first_oriented = muladd(
        first_index - one(first_index),
        metadata.sampling[1],
        metadata.origin[1],
    )
    second_oriented = muladd(
        second_index - one(second_index),
        metadata.sampling[2],
        metadata.origin[2],
    )
    first_coordinate = metadata.orientation.signs[1] * first_oriented
    second_coordinate = metadata.orientation.signs[2] * second_oriented
    return _logical_plane_coordinates(
        Val(metadata.orientation.axes),
        first_coordinate,
        second_coordinate,
    )
end

@inline _logical_plane_coordinates(
    ::Val{(:x, :y)}, first_coordinate, second_coordinate) =
    (first_coordinate, second_coordinate)

@inline _logical_plane_coordinates(
    ::Val{(:y, :x)}, first_coordinate, second_coordinate) =
    (second_coordinate, first_coordinate)

@inline function _sampled_plane_indices(
    metadata::OpticalPlaneMetadata,
    x_coordinate,
    y_coordinate,
)
    return _sampled_plane_indices(
        Val(metadata.orientation.axes),
        metadata,
        x_coordinate,
        y_coordinate,
    )
end

@inline function _sampled_plane_indices(
    ::Val{(:x, :y)},
    metadata::OpticalPlaneMetadata,
    x_coordinate,
    y_coordinate,
)
    first = one(x_coordinate) +
        (metadata.orientation.signs[1] * x_coordinate -
         metadata.origin[1]) / metadata.sampling[1]
    second = one(y_coordinate) +
        (metadata.orientation.signs[2] * y_coordinate -
         metadata.origin[2]) / metadata.sampling[2]
    return first, second
end

@inline function _sampled_plane_indices(
    ::Val{(:y, :x)},
    metadata::OpticalPlaneMetadata,
    x_coordinate,
    y_coordinate,
)
    first = one(y_coordinate) +
        (metadata.orientation.signs[1] * y_coordinate -
         metadata.origin[1]) / metadata.sampling[1]
    second = one(x_coordinate) +
        (metadata.orientation.signs[2] * x_coordinate -
         metadata.origin[2]) / metadata.sampling[2]
    return first, second
end

@inline function _pupil_footprint_sample_indices(
    destination_metadata::OpticalPlaneMetadata,
    surface_metadata::OpticalPlaneMetadata,
    registration::PupilRelayRegistration,
    shift_x_m,
    shift_y_m,
    footprint_scale,
    first_index,
    second_index,
)
    destination_x, destination_y = _logical_plane_coordinates(
        destination_metadata, first_index, second_index)
    footprint_x = muladd(
        footprint_scale, destination_x, shift_x_m)
    footprint_y = muladd(
        footprint_scale, destination_y, shift_y_m)
    m11, m12, m21, m22 = registration.transform
    surface_x = muladd(
        m12, footprint_y,
        muladd(m11, footprint_x, registration.decenter_m[1]),
    )
    surface_y = muladd(
        m22, footprint_y,
        muladd(m21, footprint_x, registration.decenter_m[2]),
    )
    return _sampled_plane_indices(
        surface_metadata, surface_x, surface_y)
end

function _prepare_pupil_footprint_index_mapping(
    destination_metadata::OpticalPlaneMetadata,
    surface_metadata::OpticalPlaneMetadata,
    registration::PupilRelayRegistration{T},
    shift_x_m::T,
    shift_y_m::T,
    footprint_scale::T,
) where {T<:AbstractFloat}
    offset = _pupil_footprint_sample_indices(
        destination_metadata,
        surface_metadata,
        registration,
        shift_x_m,
        shift_y_m,
        footprint_scale,
        zero(T),
        zero(T),
    )
    first_basis = _pupil_footprint_sample_indices(
        destination_metadata,
        surface_metadata,
        registration,
        shift_x_m,
        shift_y_m,
        footprint_scale,
        one(T),
        zero(T),
    )
    second_basis = _pupil_footprint_sample_indices(
        destination_metadata,
        surface_metadata,
        registration,
        shift_x_m,
        shift_y_m,
        footprint_scale,
        zero(T),
        one(T),
    )
    transform = (
        first_basis[1] - offset[1],
        second_basis[1] - offset[1],
        first_basis[2] - offset[2],
        second_basis[2] - offset[2],
    )
    all(isfinite, transform) && all(isfinite, offset) || throw(
        PlantPreparationError(:pupil_surface_coupling,
            :invalid_pupil_footprint_transform,
            "prepared pupil-footprint sample mapping is not finite"))
    determinant =
        transform[1] * transform[4] - transform[2] * transform[3]
    isfinite(determinant) && !iszero(determinant) || throw(
        PlantPreparationError(:pupil_surface_coupling,
            :singular_pupil_footprint_transform,
            "prepared pupil-footprint sample mapping is singular"))
    return transform, offset
end

"""
    prepare_sampled_pupil_footprint_coupling(
        surface_metadata, surface, path, placement;
        registration=nothing)

Prepare one exact path-local coupling for a sampled OPD surface. The surface
grid and path pupil must use compatible metric pupil-plane metadata, numeric
type, array backend, and physical device. `registration=nothing` selects the
identity pupil-relay registration in the path OPD numeric type.
"""
function prepare_sampled_pupil_footprint_coupling(
    surface_metadata::OpticalPlaneMetadata,
    surface::AbstractMatrix,
    path,
    placement::AbstractOpticalPlacement;
    registration=nothing,
)
    validate_plane_storage(
        surface_metadata, surface; label="sampled pupil surface")
    destination = path_input(path)
    destination_metadata =
        _require_sampled_surface_domain(destination, surface_metadata)
    T = eltype(destination.opd)
    resolved_registration =
        _resolve_pupil_relay_registration(registration, T)
    shift_x_m, shift_y_m, footprint_scale =
        _pupil_plane_footprint_geometry(placement, path.source, T)
    if iszero(shift_x_m) && iszero(shift_y_m) &&
       footprint_scale == one(T) &&
       _is_identity_pupil_relay_registration(resolved_registration) &&
       _same_pupil_sampling_grid(destination_metadata, surface_metadata)
        return PreparedIdentityPupilFootprintCoupling(
            _PREPARED_PUPIL_FOOTPRINT_COUPLING_TOKEN,
            destination, surface_metadata)
    end
    transform, offset = _prepare_pupil_footprint_index_mapping(
        destination_metadata,
        surface_metadata,
        resolved_registration,
        shift_x_m,
        shift_y_m,
        footprint_scale,
    )
    return PreparedPupilFootprintCoupling(
        _PREPARED_PUPIL_FOOTPRINT_COUPLING_TOKEN,
        destination,
        surface_metadata,
        transform,
        offset,
    )
end

function _require_sampled_pupil_surface_binding(
    destination::PupilFunction,
    surface::AbstractMatrix,
    coupling::Union{
        PreparedIdentityPupilFootprintCoupling,
        PreparedPupilFootprintCoupling,
    },
)
    destination === coupling.destination || throw(PlantPreparationError(
        :pupil_surface_coupling, :foreign_pupil_footprint_destination,
        "pupil-footprint coupling belongs to another path input"))
    metadata = coupling.surface_metadata
    size(surface) == metadata.dimensions || throw(
        PlantPreparationError(:pupil_surface_coupling,
            :surface_dimensions,
            "sampled pupil-surface dimensions $(size(surface)) " *
            "do not match prepared dimensions $(metadata.dimensions)"))
    eltype(surface) === metadata.numeric_type || throw(
        PlantPreparationError(:pupil_surface_coupling,
            :surface_numeric_type,
            "sampled pupil-surface numeric type $(eltype(surface)) " *
            "does not match prepared type $(metadata.numeric_type)"))
    typeof(backend(surface)) === typeof(metadata.backend) || throw(
        PlantPreparationError(:pupil_surface_coupling,
            :surface_backend,
            "sampled pupil-surface array backend changed after " *
            "coupling preparation"))
    plane_device(surface) == metadata.device || throw(
        PlantPreparationError(:pupil_surface_coupling,
            :surface_device,
            "sampled pupil-surface physical device changed after " *
            "coupling preparation"))
    return nothing
end

"""
    apply_sampled_pupil_surface!(destination, surface, coupling, application)

Apply one sampled OPD surface through its exact prepared pupil-footprint
coupling. `DMAdditive()` adds the sampled surface to the existing path OPD.
`DMReplace()` replaces the path OPD, writing zero outside a transformed
surface's finite support.
"""
function apply_sampled_pupil_surface!(
    destination::PupilFunction,
    surface::AbstractMatrix,
    coupling::PreparedIdentityPupilFootprintCoupling,
    ::DMAdditive,
)
    _require_sampled_pupil_surface_binding(destination, surface, coupling)
    @. destination.opd += surface
    return destination
end

function apply_sampled_pupil_surface!(
    destination::PupilFunction,
    surface::AbstractMatrix,
    coupling::PreparedIdentityPupilFootprintCoupling,
    ::DMReplace,
)
    _require_sampled_pupil_surface_binding(destination, surface, coupling)
    copyto!(destination.opd, surface)
    return destination
end

function apply_sampled_pupil_surface!(
    destination::PupilFunction,
    surface::AbstractMatrix,
    coupling::PreparedPupilFootprintCoupling,
    application::Union{DMAdditive,DMReplace},
)
    _require_sampled_pupil_surface_binding(destination, surface, coupling)
    _apply_sampled_pupil_surface!(
        execution_style(destination.opd),
        destination.opd,
        surface,
        coupling.index_transform,
        coupling.index_offset,
        application,
    )
    return destination
end

@inline apply_sampled_pupil_surface!(
    destination::PupilFunction,
    surface::AbstractMatrix,
    coupling::Union{
        PreparedIdentityPupilFootprintCoupling,
        PreparedPupilFootprintCoupling,
    },
) = apply_sampled_pupil_surface!(
    destination, surface, coupling, DMAdditive())

@inline _combine_sampled_surface_value(
    current, sampled, ::DMAdditive) = current + sampled
@inline _combine_sampled_surface_value(
    current, sampled, ::DMReplace) = sampled
@inline _outside_sampled_surface_value(current, ::DMAdditive) = current
@inline _outside_sampled_surface_value(current, ::DMReplace) = zero(current)

function _apply_sampled_pupil_surface!(
    ::ScalarCPUStyle,
    destination::AbstractMatrix{T},
    surface::AbstractMatrix{T},
    transform::NTuple{4,T},
    offset::NTuple{2,T},
    application::Union{DMAdditive,DMReplace},
) where {T<:AbstractFloat}
    Base.require_one_based_indexing(destination, surface)
    source_first_size, source_second_size = size(surface)
    m11, m12, m21, m22 = transform
    offset_first, offset_second = offset
    @inbounds for second in axes(destination, 2)
        second_t = T(second)
        for first in axes(destination, 1)
            first_t = T(first)
            source_first = muladd(
                m12, second_t,
                muladd(m11, first_t, offset_first),
            )
            source_second = muladd(
                m22, second_t,
                muladd(m21, first_t, offset_second),
            )
            if one(T) <= source_first <= T(source_first_size) &&
               one(T) <= source_second <= T(source_second_size)
                first0 = floor(Int, source_first)
                second0 = floor(Int, source_second)
                first1 = min(first0 + 1, source_first_size)
                second1 = min(second0 + 1, source_second_size)
                first_fraction = source_first - T(first0)
                second_fraction = source_second - T(second0)
                first_weight0 = one(T) - first_fraction
                second_weight0 = one(T) - second_fraction
                value0 = muladd(
                    first_fraction,
                    surface[first1, second0],
                    first_weight0 * surface[first0, second0],
                )
                value1 = muladd(
                    first_fraction,
                    surface[first1, second1],
                    first_weight0 * surface[first0, second1],
                )
                sampled = muladd(
                    second_fraction,
                    value1,
                    second_weight0 * value0,
                )
                destination[first, second] =
                    _combine_sampled_surface_value(
                        destination[first, second], sampled, application)
            else
                destination[first, second] =
                    _outside_sampled_surface_value(
                        destination[first, second], application)
            end
        end
    end
    return destination
end

@kernel function _sampled_pupil_surface_kernel!(
    destination,
    surface,
    m11,
    m12,
    m21,
    m22,
    offset_first,
    offset_second,
    source_first_size::Int,
    source_second_size::Int,
    application,
)
    first, second = @index(Global, NTuple)
    if first <= size(destination, 1) && second <= size(destination, 2)
        T = eltype(destination)
        first_t = T(first)
        second_t = T(second)
        source_first = muladd(
            m12, second_t,
            muladd(m11, first_t, offset_first),
        )
        source_second = muladd(
            m22, second_t,
            muladd(m21, first_t, offset_second),
        )
        if one(T) <= source_first <= T(source_first_size) &&
           one(T) <= source_second <= T(source_second_size)
            first0 = unsafe_trunc(Int, floor(source_first))
            second0 = unsafe_trunc(Int, floor(source_second))
            first1 = min(first0 + 1, source_first_size)
            second1 = min(second0 + 1, source_second_size)
            first_fraction = source_first - T(first0)
            second_fraction = source_second - T(second0)
            first_weight0 = one(T) - first_fraction
            second_weight0 = one(T) - second_fraction
            @inbounds begin
                value0 = muladd(
                    first_fraction,
                    surface[first1, second0],
                    first_weight0 * surface[first0, second0],
                )
                value1 = muladd(
                    first_fraction,
                    surface[first1, second1],
                    first_weight0 * surface[first0, second1],
                )
                sampled = muladd(
                    second_fraction,
                    value1,
                    second_weight0 * value0,
                )
                destination[first, second] =
                    _combine_sampled_surface_value(
                        destination[first, second], sampled, application)
            end
        else
            @inbounds destination[first, second] =
                _outside_sampled_surface_value(
                    destination[first, second], application)
        end
    end
end

function _apply_sampled_pupil_surface!(
    style::AcceleratorStyle,
    destination::AbstractMatrix{T},
    surface::AbstractMatrix{T},
    transform::NTuple{4,T},
    offset::NTuple{2,T},
    application::Union{DMAdditive,DMReplace},
) where {T<:AbstractFloat}
    m11, m12, m21, m22 = transform
    launch_kernel!(
        style,
        _sampled_pupil_surface_kernel!,
        destination,
        surface,
        m11,
        m12,
        m21,
        m22,
        offset[1],
        offset[2],
        size(surface, 1),
        size(surface, 2),
        application;
        ndrange=size(destination),
    )
    return destination
end

"""
    prepare_controllable_optic_path_coupling(
        implementation, definition, path)

Qualified model-extension seam for immutable path-local optical coupling. The
default supports a model-owned same-grid pupil-plane application. Models with
sampled atmospheric-conjugate surfaces prepare an explicit coupling through
`prepare_sampled_pupil_footprint_coupling`.
"""
function prepare_controllable_optic_path_coupling(
    implementation,
    definition::ControllableOpticDefinition,
    path,
)
    return _prepare_controllable_optic_path_coupling(
        controllable_optic_execution_role(implementation),
        controllable_optic_placement(definition),
        implementation,
        definition,
        path,
    )
end

@inline function _prepare_controllable_optic_path_coupling(
    ::PupilSurfaceExecutionRole,
    ::PupilPlanePlacement,
    implementation,
    ::ControllableOpticDefinition,
    path,
)
    destination = _require_direct_pupil_surface_destination(
        path_input(path), implementation)
    return PreparedDirectPupilSurfaceCoupling(
        _PREPARED_PUPIL_FOOTPRINT_COUPLING_TOKEN,
        destination,
    )
end

@inline _require_direct_pupil_surface_destination(
    destination::PupilFunction, implementation) = destination

function _require_direct_pupil_surface_destination(
    destination, implementation)
    throw(PlantPreparationError(
        :controllable_optic, :unsupported_path_input,
        "pupil-surface controllable optic $(typeof(implementation)) " *
        "requires a PupilFunction path input; got $(typeof(destination))"))
end

function _prepare_controllable_optic_path_coupling(
    ::PupilSurfaceExecutionRole,
    placement::AtmosphericConjugatePlacement,
    implementation,
    definition::ControllableOpticDefinition,
    path,
)
    throw(PlantPreparationError(:controllable_optic,
        :unsupported_conjugate_geometry,
        "controllable optic $(controllable_optic_id(definition)) with " *
        "prepared model $(typeof(implementation)) does not implement " *
        "atmospheric-conjugate path coupling at " *
        "$(conjugate_altitude_m(placement)) m for " *
        "$(path_id(path.definition))"))
end

function _prepare_controllable_optic_path_coupling(
    ::PupilSurfaceExecutionRole,
    placement::FocalPlanePlacement,
    implementation,
    definition::ControllableOpticDefinition,
    path,
)
    throw(PlantPreparationError(:controllable_optic,
        :invalid_optic_placement,
        "pupil-surface controllable optic " *
        "$(controllable_optic_id(definition)) cannot use " *
        "$(typeof(placement)) on $(path_id(path.definition)); prepared " *
        "model is $(typeof(implementation))"))
end

function _prepare_controllable_optic_path_coupling(
    role::AbstractControllableOpticExecutionRole,
    placement::AbstractOpticalPlacement,
    implementation,
    definition::ControllableOpticDefinition,
    path,
)
    throw(PlantPreparationError(:controllable_optic,
        :unsupported_optic_execution_role,
        "controllable optic $(controllable_optic_id(definition)) has " *
        "unsupported execution role $(typeof(role)) at " *
        "$(typeof(placement)) on $(path_id(path.definition)); prepared " *
        "model is $(typeof(implementation))"))
end

function apply_controllable_optic_surface!(
    input,
    implementation,
    state,
    coupling::AbstractPupilSurfacePathCoupling,
)
    throw(PlantPreparationError(:controllable_optic,
        :unsupported_surface_application,
        "prepared controllable-optic implementation " *
        "$(typeof(implementation)) does not apply through " *
        "$(typeof(coupling)) to path input $(typeof(input))"))
end

@inline function _same_pupil_footprint_coupling(
    left::PreparedDirectPupilSurfaceCoupling,
    right::PreparedDirectPupilSurfaceCoupling,
)
    return left.destination === right.destination
end

@inline function _same_pupil_footprint_coupling(
    left::PreparedIdentityPupilFootprintCoupling,
    right::PreparedIdentityPupilFootprintCoupling,
)
    return left.destination === right.destination &&
        _same_pupil_sampling_grid(
            left.surface_metadata, right.surface_metadata)
end

@inline function _same_pupil_footprint_coupling(
    left::PreparedPupilFootprintCoupling,
    right::PreparedPupilFootprintCoupling,
)
    return left.destination === right.destination &&
        left.index_transform == right.index_transform &&
        left.index_offset == right.index_offset &&
        _same_pupil_sampling_grid(
            left.surface_metadata, right.surface_metadata)
end

@inline _same_pupil_footprint_coupling(
    ::AbstractPupilSurfacePathCoupling,
    ::AbstractPupilSurfacePathCoupling,
) = false

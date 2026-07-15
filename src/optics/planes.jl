abstract type AbstractOpticalProduct end

abstract type AbstractOpticalPlaneKind end
struct PupilPlane <: AbstractOpticalPlaneKind end
struct FocalPlane <: AbstractOpticalPlaneKind end
struct PropagationPlane <: AbstractOpticalPlaneKind end
struct DetectorPlane <: AbstractOpticalPlaneKind end

@enum PlaneCentering begin
    SampleCentered
    InterpixelCentered
end

struct PlaneAxisOrientation
    axes::NTuple{2,Symbol}
    signs::NTuple{2,Int8}

    function PlaneAxisOrientation(axes::NTuple{2,Symbol}=(:x, :y),
        signs::NTuple{2,<:Integer}=(1, 1))
        axes[1] != axes[2] || throw(InvalidConfiguration(
            "optical-plane axes must be distinct"))
        all(sign -> sign == -1 || sign == 1, signs) ||
            throw(InvalidConfiguration(
                "optical-plane axis signs must be -1 or 1"))
        return new(axes, (Int8(signs[1]), Int8(signs[2])))
    end
end

abstract type AbstractSpectralCoordinate end
struct AchromaticChannel <: AbstractSpectralCoordinate end

struct MonochromaticChannel{T<:AbstractFloat} <: AbstractSpectralCoordinate
    wavelength_m::T

    function MonochromaticChannel{T}(wavelength_m::T) where {T<:AbstractFloat}
        wavelength_m > zero(T) || throw(InvalidConfiguration(
            "optical-plane wavelength must be positive"))
        return new{T}(wavelength_m)
    end
end

MonochromaticChannel(wavelength_m::T) where {T<:AbstractFloat} =
    MonochromaticChannel{T}(wavelength_m)

# These markers reserve the radiometric contract fields that Pre-HIL 5/12
# finalizes. They prevent an implicit default from being mistaken for a
# physical normalization, spatial measure, or coherence policy.
struct UnspecifiedNormalization end
struct UnspecifiedSpatialMeasure end
struct UnspecifiedCoherence end

abstract type AbstractPlaneDevice end
struct HostPlaneDevice <: AbstractPlaneDevice end

struct AcceleratorPlaneDevice{B<:KernelAbstractions.Backend,I} <:
    AbstractPlaneDevice
    backend::B
    identifier::I
end

@inline plane_device(storage::AbstractArray) =
    plane_device(execution_style(storage), storage)
@inline physical_device_identifier(::AbstractArray) = nothing
@inline plane_device(::ScalarCPUStyle, ::AbstractArray) = HostPlaneDevice()
@inline plane_device(style::AcceleratorStyle, storage::AbstractArray) =
    AcceleratorPlaneDevice(style.backend,
        physical_device_identifier(storage))

struct OpticalPlaneMetadata{
    T<:AbstractFloat,
    E,
    K<:AbstractOpticalPlaneKind,
    S<:AbstractSpectralCoordinate,
    N,
    M,
    C,
    B<:AbstractArrayBackend,
    D<:AbstractPlaneDevice,
}
    kind::K
    dimensions::NTuple{2,Int}
    sampling::NTuple{2,T}
    origin::NTuple{2,T}
    centering::NTuple{2,PlaneCentering}
    orientation::PlaneAxisOrientation
    spectral::S
    numeric_type::Type{E}
    normalization::N
    spatial_measure::M
    coherence::C
    backend::B
    device::D
end

@inline axis_centering(n::Int) =
    isodd(n) ? SampleCentered : InterpixelCentered

@inline function centered_grid_origin(dimensions::NTuple{2,Int},
    sampling::NTuple{2,T}) where {T<:AbstractFloat}
    return (
        -T(dimensions[1] - 1) * sampling[1] / T(2),
        -T(dimensions[2] - 1) * sampling[2] / T(2),
    )
end

function OpticalPlaneMetadata(kind::AbstractOpticalPlaneKind,
    storage::AbstractMatrix{E};
    sampling::NTuple{2,T},
    origin::NTuple{2,T}=centered_grid_origin(size(storage), sampling),
    centering::NTuple{2,PlaneCentering}=(
        axis_centering(size(storage, 1)),
        axis_centering(size(storage, 2)),
    ),
    orientation::PlaneAxisOrientation=PlaneAxisOrientation(),
    spectral::AbstractSpectralCoordinate=AchromaticChannel(),
    normalization=UnspecifiedNormalization(),
    spatial_measure=UnspecifiedSpatialMeasure(),
    coherence=UnspecifiedCoherence(),
    device::AbstractPlaneDevice=plane_device(storage),
) where {T<:AbstractFloat,E}
    selector = backend(storage)
    return OpticalPlaneMetadata{
        T,E,typeof(kind),typeof(spectral),typeof(normalization),
        typeof(spatial_measure),typeof(coherence),typeof(selector),
        typeof(device),
    }(
        kind,
        size(storage),
        sampling,
        origin,
        centering,
        orientation,
        spectral,
        E,
        normalization,
        spatial_measure,
        coherence,
        selector,
        device,
    )
end

@inline plane_metadata(product::AbstractOpticalProduct) = product.metadata

function validate_plane_storage(metadata::OpticalPlaneMetadata,
    storage::AbstractMatrix; label::AbstractString="optical plane")
    size(storage) == metadata.dimensions || throw(DimensionMismatchError(
        "$label storage dimensions $(size(storage)) do not match declared dimensions $(metadata.dimensions)"))
    eltype(storage) === metadata.numeric_type || throw(InvalidConfiguration(
        "$label storage element type $(eltype(storage)) does not match declared numeric type $(metadata.numeric_type)"))
    typeof(backend(storage)) === typeof(metadata.backend) ||
        throw(InvalidConfiguration(
            "$label storage backend does not match declared backend"))
    plane_device(storage) == metadata.device || throw(InvalidConfiguration(
        "$label storage device does not match declared device"))
    return metadata
end

function require_same_plane_grid(a::OpticalPlaneMetadata,
    b::OpticalPlaneMetadata; label::AbstractString="optical planes",
    require_kind::Bool=true, require_spectral::Bool=true,
    require_numeric_type::Bool=true)
    (!require_kind || typeof(a.kind) === typeof(b.kind)) ||
        throw(InvalidConfiguration("$label have incompatible plane kinds"))
    a.dimensions == b.dimensions || throw(DimensionMismatchError(
        "$label have incompatible dimensions $(a.dimensions) and $(b.dimensions)"))
    a.sampling == b.sampling || throw(InvalidConfiguration(
        "$label have incompatible physical sampling"))
    a.origin == b.origin || throw(InvalidConfiguration(
        "$label have incompatible origins"))
    a.centering == b.centering || throw(InvalidConfiguration(
        "$label have incompatible centering"))
    a.orientation == b.orientation || throw(InvalidConfiguration(
        "$label have incompatible axis orientation"))
    (!require_spectral || a.spectral == b.spectral) ||
        throw(InvalidConfiguration(
            "$label have incompatible spectral coordinates"))
    (!require_numeric_type || a.numeric_type === b.numeric_type) ||
        throw(InvalidConfiguration("$label have incompatible numeric types"))
    typeof(a.backend) === typeof(b.backend) || throw(InvalidConfiguration(
        "$label have incompatible backends"))
    a.device == b.device || throw(InvalidConfiguration(
        "$label are on different physical devices"))
    return nothing
end

function require_centered_plane_geometry(metadata::OpticalPlaneMetadata;
    label::AbstractString="optical plane")
    expected_centering = (
        axis_centering(metadata.dimensions[1]),
        axis_centering(metadata.dimensions[2]),
    )
    metadata.centering == expected_centering ||
        throw(InvalidConfiguration(
            "$label centering is incompatible with its dimensions"))
    expected_origin = centered_grid_origin(metadata.dimensions,
        metadata.sampling)
    metadata.origin == expected_origin || throw(InvalidConfiguration(
        "$label origin is incompatible with centered-grid propagation"))
    return metadata
end

struct PupilWavefront{
    M<:OpticalPlaneMetadata,
    A<:AbstractMatrix,
    O<:AbstractMatrix,
    B<:AbstractArrayBackend,
} <: AbstractOpticalProduct
    metadata::M
    amplitude::A
    opd::O
end

@inline backend(::PupilWavefront{<:Any,<:Any,<:Any,B}) where {B} = B()
@inline pupil_amplitude(wavefront::PupilWavefront) = wavefront.amplitude
@inline opd_map(wavefront::PupilWavefront) = wavefront.opd

function PupilWavefront(tel::Telescope;
    T::Type{<:AbstractFloat}=eltype(opd_map(tel)),
    backend::AbstractArrayBackend=backend(tel))
    selector = require_same_backend(tel, _resolve_backend_selector(backend))
    amplitude = similar(pupil_reflectivity(tel), T,
        tel.params.resolution, tel.params.resolution)
    opd = similar(opd_map(tel), T, tel.params.resolution,
        tel.params.resolution)
    reflectivity = pupil_reflectivity(tel)
    @. amplitude = sqrt(reflectivity)
    fill!(opd, zero(T))
    sampling = (T(tel.aperture.sampling_m[1]),
        T(tel.aperture.sampling_m[2]))
    origin = (T(tel.aperture.origin_m[1]), T(tel.aperture.origin_m[2]))
    metadata = OpticalPlaneMetadata(PupilPlane(), opd;
        sampling=sampling, origin=origin)
    validate_plane_storage(metadata, amplitude; label="pupil amplitude")
    return PupilWavefront{
        typeof(metadata),typeof(amplitude),typeof(opd),typeof(selector),
    }(metadata, amplitude, opd)
end

function reset_opd!(wavefront::PupilWavefront)
    fill!(wavefront.opd, zero(eltype(wavefront.opd)))
    return wavefront
end

function apply_opd!(wavefront::PupilWavefront, opd::AbstractMatrix)
    size(opd) == size(wavefront.opd) || throw(DimensionMismatchError(
        "OPD size does not match pupil-wavefront dimensions"))
    require_same_backend(wavefront, opd)
    copyto!(wavefront.opd, opd)
    return wavefront
end

function _validate_surface_application(wavefront::PupilWavefront,
    surface)
    opd = surface_opd(surface)
    size(opd) == wavefront.metadata.dimensions ||
        throw(DimensionMismatchError(
            "optical-surface OPD dimensions do not match PupilWavefront"))
    require_same_backend(wavefront, opd)
    plane_device(opd) == wavefront.metadata.device ||
        throw(InvalidConfiguration(
            "optical surface and PupilWavefront occupy different physical devices"))
    return opd
end

"""
    apply_surface!(wavefront, surface, mode)

Apply an already formed optical-surface OPD to an explicit caller-owned pupil
wavefront. Controllable surfaces must first be formed with `update_surface!`;
this operation never mutates telescope path state.
"""
function apply_surface!(wavefront::PupilWavefront, surface, ::DMAdditive)
    surface_values = _validate_surface_application(wavefront, surface)
    @. wavefront.opd += surface_values
    return wavefront
end

function apply_surface!(wavefront::PupilWavefront, surface, ::DMReplace)
    surface_values = _validate_surface_application(wavefront, surface)
    copyto!(wavefront.opd, surface_values)
    return wavefront
end

struct ElectricField{
    M<:OpticalPlaneMetadata,
    A<:AbstractMatrix,
    B<:AbstractArrayBackend,
} <: AbstractOpticalProduct
    metadata::M
    values::A
end

@inline backend(::ElectricField{<:Any,<:Any,B}) where {B} = B()
@inline field_values(field::ElectricField) = field.values

function ElectricField(metadata::OpticalPlaneMetadata,
    values::AbstractMatrix{<:Complex})
    validate_plane_storage(metadata, values; label="electric field")
    metadata.numeric_type <: Complex || throw(InvalidConfiguration(
        "ElectricField storage must use a complex numeric type"))
    selector = backend(values)
    return ElectricField{typeof(metadata),typeof(values),typeof(selector)}(
        metadata, values)
end

struct IrradiancePlane{
    M<:OpticalPlaneMetadata,
    A<:AbstractMatrix,
    B<:AbstractArrayBackend,
} <: AbstractOpticalProduct
    metadata::M
    values::A
end

@inline backend(::IrradiancePlane{<:Any,<:Any,B}) where {B} = B()
@inline irradiance_values(plane::IrradiancePlane) = plane.values

function IrradiancePlane(metadata::OpticalPlaneMetadata,
    values::AbstractMatrix{<:AbstractFloat})
    validate_plane_storage(metadata, values; label="irradiance plane")
    selector = backend(values)
    return IrradiancePlane{typeof(metadata),typeof(values),typeof(selector)}(
        metadata, values)
end

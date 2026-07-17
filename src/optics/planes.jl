abstract type AbstractOpticalProduct end

abstract type AbstractOpticalPlaneKind end
struct PupilPlane <: AbstractOpticalPlaneKind end
struct FocalPlane <: AbstractOpticalPlaneKind end
struct IntermediatePlane <: AbstractOpticalPlaneKind end
struct DetectorPlane <: AbstractOpticalPlaneKind end

abstract type AbstractPlaneCoordinateDomain end

"""Metric plane coordinates whose sampling and origin are expressed in metres."""
struct MetricCoordinates <: AbstractPlaneCoordinateDomain end

"""Angular plane coordinates whose sampling and origin are expressed in radians."""
struct AngularCoordinates <: AbstractPlaneCoordinateDomain end

@enum PlaneCentering begin
    SampleCentered
    InterSampleCentered
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
struct UnspecifiedSpectralCoordinate <: AbstractSpectralCoordinate end

"""The product is explicitly independent of wavelength."""
struct AchromaticSpectralCoordinate <: AbstractSpectralCoordinate end

struct MonochromaticChannel{T<:AbstractFloat} <: AbstractSpectralCoordinate
    wavelength_m::T

    function MonochromaticChannel{T}(wavelength_m::T) where {T<:AbstractFloat}
        isfinite(wavelength_m) && wavelength_m > zero(T) ||
            throw(InvalidConfiguration(
                "optical-plane wavelength must be finite and positive"))
        return new{T}(wavelength_m)
    end
end

MonochromaticChannel(wavelength_m::T) where {T<:AbstractFloat} =
    MonochromaticChannel{T}(wavelength_m)

"""
    IntegratedSpectralChannel(identifier)

Declare that a product has already been integrated over an application-defined
spectral channel. Products are compatible only when their identifiers match.
The identifier names a passband contract owned by the application; it does not
make a wavelength-dependent detector response safe to apply after integration.
"""
struct IntegratedSpectralChannel <: AbstractSpectralCoordinate
    identifier::Symbol

    function IntegratedSpectralChannel(identifier::Symbol)
        isempty(String(identifier)) && throw(InvalidConfiguration(
            "integrated spectral-channel identifier must not be empty"))
        return new(identifier)
    end
end

abstract type AbstractOpticalNormalization end

"""Radiometry has not been declared and cannot enter physical acquisition."""
struct UnspecifiedNormalization <: AbstractOpticalNormalization end

"""Values, or `abs2` of field values, carry physical photon-arrival rate."""
struct PhotonRateNormalization <: AbstractOpticalNormalization end

"""Values carry explicitly dimensionless relative optical power."""
struct DimensionlessNormalization <: AbstractOpticalNormalization end

abstract type AbstractSpatialMeasure end

"""The spatial interpretation of each sample has not been declared."""
struct UnspecifiedSpatialMeasure <: AbstractSpatialMeasure end

"""Samples are point values and do not include a represented cell measure."""
struct PointSampledMeasure <: AbstractSpatialMeasure end

"""Each value is a density per unit area of its declared coordinate domain."""
struct SpatialDensityMeasure <: AbstractSpatialMeasure end

"""Each value is integrated over its represented spatial cell."""
struct CellIntegratedMeasure <: AbstractSpatialMeasure end

abstract type AbstractCombinationPolicy end

"""The product's coherent or incoherent combination policy is undeclared."""
struct UnspecifiedCoherence <: AbstractCombinationPolicy end

"""Complex fields may be combined coherently when all field contracts match."""
struct CoherentFieldCombination <: AbstractCombinationPolicy end

"""Intensity products may be accumulated elementwise after compatibility checks."""
struct IncoherentIntensityAddition <: AbstractCombinationPolicy end

"""The product must remain separate unless an explicit mapping is prepared."""
struct NonCombinableProduct <: AbstractCombinationPolicy end

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
@inline physical_device_identifier(storage::SubArray) =
    physical_device_identifier(parent(storage))
@inline physical_device_identifier(storage::Base.ReshapedArray) =
    physical_device_identifier(parent(storage))
@inline physical_device_identifier(storage::Base.ReinterpretArray) =
    physical_device_identifier(parent(storage))
@inline physical_device_identifier(storage::PermutedDimsArray) =
    physical_device_identifier(parent(storage))
@inline physical_device_identifier(storage::Transpose) =
    physical_device_identifier(parent(storage))
@inline physical_device_identifier(storage::Adjoint) =
    physical_device_identifier(parent(storage))
@inline plane_device(::ScalarCPUStyle, ::AbstractArray) = HostPlaneDevice()
@inline plane_device(style::AcceleratorStyle, storage::AbstractArray) =
    AcceleratorPlaneDevice(style.backend,
        physical_device_identifier(storage))

struct OpticalPlaneMetadata{
    T<:AbstractFloat,
    E,
    K<:AbstractOpticalPlaneKind,
    Q<:AbstractPlaneCoordinateDomain,
    S<:AbstractSpectralCoordinate,
    N<:AbstractOpticalNormalization,
    M<:AbstractSpatialMeasure,
    C<:AbstractCombinationPolicy,
    B<:AbstractArrayBackend,
    D<:AbstractPlaneDevice,
}
    kind::K
    coordinate_domain::Q
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
    isodd(n) ? SampleCentered : InterSampleCentered

@inline function centered_grid_origin(dimensions::NTuple{2,Int},
    sampling::NTuple{2,T}) where {T<:AbstractFloat}
    return (
        -T(dimensions[1] - 1) * sampling[1] / T(2),
        -T(dimensions[2] - 1) * sampling[2] / T(2),
    )
end

function OpticalPlaneMetadata(kind::AbstractOpticalPlaneKind,
    storage::AbstractMatrix{E};
    coordinate_domain::AbstractPlaneCoordinateDomain,
    sampling::NTuple{2,T},
    origin::NTuple{2,T}=centered_grid_origin(size(storage), sampling),
    centering::NTuple{2,PlaneCentering}=(
        axis_centering(size(storage, 1)),
        axis_centering(size(storage, 2)),
    ),
    orientation::PlaneAxisOrientation=PlaneAxisOrientation(),
    spectral::AbstractSpectralCoordinate=UnspecifiedSpectralCoordinate(),
    normalization::AbstractOpticalNormalization=UnspecifiedNormalization(),
    spatial_measure::AbstractSpatialMeasure=UnspecifiedSpatialMeasure(),
    coherence::AbstractCombinationPolicy=UnspecifiedCoherence(),
    device::AbstractPlaneDevice=plane_device(storage),
) where {T<:AbstractFloat,E}
    all(value -> isfinite(value) && value > zero(T), sampling) ||
        throw(InvalidConfiguration(
            "optical-plane sampling must be finite and positive"))
    all(isfinite, origin) || throw(InvalidConfiguration(
        "optical-plane origin must be finite"))
    selector = backend(storage)
    return OpticalPlaneMetadata{
        T,E,typeof(kind),typeof(coordinate_domain),typeof(spectral),typeof(normalization),
        typeof(spatial_measure),typeof(coherence),typeof(selector),
        typeof(device),
    }(
        kind,
        coordinate_domain,
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
    typeof(a.coordinate_domain) === typeof(b.coordinate_domain) ||
        throw(InvalidConfiguration(
            "$label have incompatible coordinate domains"))
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

function require_compatible_radiometry(a::OpticalPlaneMetadata,
    b::OpticalPlaneMetadata; label::AbstractString="optical products")
    typeof(a.normalization) === typeof(b.normalization) ||
        throw(InvalidConfiguration(
            "$label have incompatible radiometric normalizations"))
    typeof(a.spatial_measure) === typeof(b.spatial_measure) ||
        throw(InvalidConfiguration(
            "$label have incompatible spatial measures"))
    return nothing
end

function require_incoherent_addition(a::OpticalPlaneMetadata,
    b::OpticalPlaneMetadata; label::AbstractString="intensity products")
    require_same_plane_grid(a, b; label=label)
    require_compatible_radiometry(a, b; label=label)
    _require_incoherent_policy(a.coherence, label)
    _require_incoherent_policy(b.coherence, label)
    return nothing
end

@inline _require_incoherent_policy(::IncoherentIntensityAddition,
    ::AbstractString) = nothing

function _require_incoherent_policy(::AbstractCombinationPolicy,
    label::AbstractString)
    throw(InvalidConfiguration(
        "$label must declare incoherent intensity addition"))
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

@inline require_metric_coordinates(metadata::OpticalPlaneMetadata;
    label::AbstractString="optical plane") =
    require_metric_coordinates(metadata.coordinate_domain, label)

@inline require_metric_coordinates(::MetricCoordinates, ::AbstractString) =
    nothing

function require_metric_coordinates(::AbstractPlaneCoordinateDomain,
    label::AbstractString)
    throw(InvalidConfiguration("$label must use metric coordinates"))
end

"""
    PupilFunction

Caller-owned pupil-plane amplitude and optical-path-difference representation.
It is converted to a wavelength-specific scalar `ElectricField` by a prepared
field-formation plan.
"""
struct PupilFunction{
    M<:OpticalPlaneMetadata,
    A<:AbstractMatrix,
    O<:AbstractMatrix,
    B<:AbstractArrayBackend,
} <: AbstractOpticalProduct
    metadata::M
    amplitude::A
    opd::O
end

@inline backend(::PupilFunction{<:Any,<:Any,<:Any,B}) where {B} = B()
@inline pupil_amplitude(pupil::PupilFunction) = pupil.amplitude
@inline opd_map(pupil::PupilFunction) = pupil.opd

function PupilFunction(tel::Telescope;
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
        coordinate_domain=MetricCoordinates(), sampling=sampling, origin=origin,
        spectral=AchromaticSpectralCoordinate(),
        normalization=DimensionlessNormalization(),
        spatial_measure=PointSampledMeasure(),
        coherence=CoherentFieldCombination())
    validate_plane_storage(metadata, amplitude; label="pupil amplitude")
    return PupilFunction{
        typeof(metadata),typeof(amplitude),typeof(opd),typeof(selector),
    }(metadata, amplitude, opd)
end

function reset_opd!(pupil::PupilFunction)
    fill!(pupil.opd, zero(eltype(pupil.opd)))
    return pupil
end

function apply_opd!(pupil::PupilFunction, opd::AbstractMatrix)
    size(opd) == size(pupil.opd) || throw(DimensionMismatchError(
        "OPD size does not match pupil-function dimensions"))
    require_same_backend(pupil, opd)
    copyto!(pupil.opd, opd)
    return pupil
end

function _validate_surface_application(pupil::PupilFunction,
    surface)
    opd = surface_opd(surface)
    size(opd) == pupil.metadata.dimensions ||
        throw(DimensionMismatchError(
            "optical-surface OPD dimensions do not match PupilFunction"))
    require_same_backend(pupil, opd)
    plane_device(opd) == pupil.metadata.device ||
        throw(InvalidConfiguration(
            "optical surface and PupilFunction occupy different physical devices"))
    return opd
end

"""
    apply_surface!(pupil, surface, mode)

Apply an already formed optical-surface OPD to an explicit caller-owned pupil
function. Controllable surfaces must first be formed with `update_surface!`;
this operation never mutates telescope path state.
"""
function apply_surface!(pupil::PupilFunction, surface, ::DMAdditive)
    surface_values = _validate_surface_application(pupil, surface)
    @. pupil.opd += surface_values
    return pupil
end

function apply_surface!(pupil::PupilFunction, surface, ::DMReplace)
    surface_values = _validate_surface_application(pupil, surface)
    copyto!(pupil.opd, surface_values)
    return pupil
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

"""
    IntensityMap

Caller-owned real-valued samples of computational optical intensity (`abs2` of
a scalar `ElectricField`). The plane location, coordinate domain, and physical
normalization are properties of `metadata`; this type alone does not claim SI
irradiance units.
"""
struct IntensityMap{
    M<:OpticalPlaneMetadata,
    A<:AbstractMatrix,
    B<:AbstractArrayBackend,
} <: AbstractOpticalProduct
    metadata::M
    values::A
end

@inline backend(::IntensityMap{<:Any,<:Any,B}) where {B} = B()
@inline intensity_values(map::IntensityMap) = map.values

function IntensityMap(metadata::OpticalPlaneMetadata,
    values::AbstractMatrix{<:AbstractFloat})
    validate_plane_storage(metadata, values; label="intensity map")
    selector = backend(values)
    return IntensityMap{typeof(metadata),typeof(values),typeof(selector)}(
        metadata, values)
end


"""
    OpticalProductBundle(products...)

Preserve prepared optical products that cannot be combined on one physical
grid. A bundle performs no implicit resampling or radiometric conversion.
"""
struct OpticalProductBundle{P<:Union{Tuple,AbstractVector}}
    products::P
end

OpticalProductBundle(products::Vararg{AbstractOpticalProduct,N}) where {N} =
    OpticalProductBundle{typeof(products)}(products)

function OpticalProductBundle(
    products::AbstractVector{<:AbstractOpticalProduct})
    owned = copy(products)
    return OpticalProductBundle{typeof(owned)}(owned)
end

Base.length(bundle::OpticalProductBundle) = length(bundle.products)
Base.getindex(bundle::OpticalProductBundle, index::Int) = bundle.products[index]
Base.iterate(bundle::OpticalProductBundle, state...) =
    iterate(bundle.products, state...)

struct _PreparedIncoherentSumToken end
const _PREPARED_INCOHERENT_SUM_TOKEN = _PreparedIncoherentSumToken()

struct PreparedIncoherentInputs{M,V}
    metadata::M
    values::V
end

struct PreparedIncoherentSum{O<:OpticalPlaneMetadata,I}
    output_metadata::O
    inputs::I

    function PreparedIncoherentSum(::_PreparedIncoherentSumToken,
        output_metadata::O, inputs::I) where {
        O<:OpticalPlaneMetadata,I,
    }
        return new{O,I}(output_metadata, inputs)
    end
end

function prepare_incoherent_sum(output::IntensityMap,
    inputs::Vararg{IntensityMap,N}) where {N}
    N > 0 || throw(InvalidConfiguration(
        "incoherent accumulation requires at least one input"))
    _require_declared_intensity_spectral(output.metadata.spectral,
        "incoherent intensity output")
    _require_nonaliasing_intensity_inputs(output.values, inputs)
    metadata = ntuple(index -> begin
        input = inputs[index]
        _require_declared_intensity_spectral(input.metadata.spectral,
            "incoherent intensity input")
        require_incoherent_addition(output.metadata, input.metadata;
            label="incoherent intensity accumulation")
        input.metadata
    end, N)
    values = ntuple(index -> inputs[index].values, N)
    return PreparedIncoherentSum(_PREPARED_INCOHERENT_SUM_TOKEN,
        output.metadata, PreparedIncoherentInputs(metadata, values))
end

function prepare_incoherent_sum(output::IntensityMap,
    inputs::AbstractVector{<:IntensityMap})
    isempty(inputs) && throw(InvalidConfiguration(
        "incoherent accumulation requires at least one input"))
    _require_declared_intensity_spectral(output.metadata.spectral,
        "incoherent intensity output")
    _require_nonaliasing_intensity_inputs(output.values, inputs)

    first_input = first(inputs)
    first_metadata = first_input.metadata
    first_values = first_input.values
    metadata = Vector{typeof(first_metadata)}(undef, length(inputs))
    values = Vector{typeof(first_values)}(undef, length(inputs))
    @inbounds for index in eachindex(inputs)
        input = inputs[index]
        _require_homogeneous_prepared_input(typeof(first_metadata),
            typeof(input.metadata))
        _require_homogeneous_prepared_values(typeof(first_values),
            typeof(input.values))
        _require_declared_intensity_spectral(input.metadata.spectral,
            "incoherent intensity input")
        require_incoherent_addition(output.metadata, input.metadata;
            label="incoherent intensity accumulation")
        metadata[index] = input.metadata
        values[index] = input.values
    end
    return PreparedIncoherentSum(_PREPARED_INCOHERENT_SUM_TOKEN,
        output.metadata, PreparedIncoherentInputs(metadata, values))
end

@inline function _require_homogeneous_prepared_input(
    ::Type{R}, ::Type{I}) where {
    R<:OpticalPlaneMetadata,I<:OpticalPlaneMetadata,
}
    R === I || throw(InvalidConfiguration(
        "vector incoherent inputs must have one concrete metadata type; got " *
        "$(R) and $(I)"))
    return nothing
end

@inline function _require_homogeneous_prepared_values(
    ::Type{R}, ::Type{I}) where {R<:AbstractArray,I<:AbstractArray}
    R === I || throw(InvalidConfiguration(
        "vector incoherent inputs must have one concrete storage type; got " *
        "$(R) and $(I)"))
    return nothing
end

@inline _require_declared_intensity_spectral(::MonochromaticChannel,
    ::AbstractString) = nothing

@inline _require_declared_intensity_spectral(::AchromaticSpectralCoordinate,
    ::AbstractString) = nothing

@inline _require_declared_intensity_spectral(::IntegratedSpectralChannel,
    ::AbstractString) = nothing

function _require_declared_intensity_spectral(::AbstractSpectralCoordinate,
    label::AbstractString)
    throw(InvalidConfiguration(
        "$label requires an explicit achromatic, monochromatic, or integrated spectral coordinate"))
end

function accumulate_intensity!(output::IntensityMap,
    inputs::Tuple, plan::PreparedIncoherentSum)
    output.metadata === plan.output_metadata || throw(InvalidConfiguration(
        "intensity output does not match its prepared accumulation plan"))
    length(inputs) == length(plan.inputs.metadata) ||
        throw(DimensionMismatchError(
            "intensity input count does not match its prepared accumulation plan"))
    _require_prepared_intensity_inputs(output.values, inputs,
        plan.inputs.metadata, plan.inputs.values)
    fill!(output.values, zero(eltype(output.values)))
    _accumulate_prepared_intensity!(output.values, inputs)
    return output
end


function accumulate_intensity!(output::IntensityMap,
    inputs::AbstractVector{<:IntensityMap}, plan::PreparedIncoherentSum)
    output.metadata === plan.output_metadata || throw(InvalidConfiguration(
        "intensity output does not match its prepared accumulation plan"))
    length(inputs) == length(plan.inputs.metadata) ||
        throw(DimensionMismatchError(
            "intensity input count does not match its prepared accumulation plan"))
    _require_prepared_intensity_inputs(output.values, inputs,
        plan.inputs.metadata, plan.inputs.values)
    fill!(output.values, zero(eltype(output.values)))
    _accumulate_prepared_intensity!(output.values, inputs)
    return output
end

@inline _require_nonaliasing_intensity_inputs(::AbstractArray, ::Tuple{}) =
    nothing

@inline function _require_nonaliasing_intensity_inputs(output::AbstractArray,
    inputs::Tuple)
    input = first(inputs)
    Base.mightalias(output, input.values) && throw(InvalidConfiguration(
        "incoherent intensity output must not alias an input"))
    return _require_nonaliasing_intensity_inputs(output, Base.tail(inputs))
end


function _require_nonaliasing_intensity_inputs(output::AbstractArray,
    inputs::AbstractVector{<:IntensityMap})
    @inbounds for input in inputs
        Base.mightalias(output, input.values) && throw(InvalidConfiguration(
            "incoherent intensity output must not alias an input"))
    end
    return nothing
end

@inline _require_prepared_intensity_inputs(::AbstractArray, ::Tuple{},
    ::Tuple{}, ::Tuple{}) = nothing

@inline function _require_prepared_intensity_inputs(output::AbstractArray,
    inputs::Tuple, metadata::Tuple, values::Tuple)
    input = first(inputs)
    input.metadata === first(metadata) || throw(InvalidConfiguration(
        "intensity input does not match its prepared accumulation plan"))
    input.values === first(values) || throw(InvalidConfiguration(
        "intensity input storage does not match its prepared accumulation plan"))
    Base.mightalias(output, input.values) && throw(InvalidConfiguration(
        "incoherent intensity output must not alias an input"))
    return _require_prepared_intensity_inputs(output, Base.tail(inputs),
        Base.tail(metadata), Base.tail(values))
end


function _require_prepared_intensity_inputs(output::AbstractArray,
    inputs::AbstractVector{<:IntensityMap}, metadata::AbstractVector,
    values::AbstractVector)
    @inbounds for index in eachindex(inputs, metadata, values)
        input = inputs[index]
        input.metadata === metadata[index] || throw(InvalidConfiguration(
            "intensity input does not match its prepared accumulation plan"))
        input.values === values[index] || throw(InvalidConfiguration(
            "intensity input storage does not match its prepared accumulation plan"))
        Base.mightalias(output, input.values) && throw(InvalidConfiguration(
            "incoherent intensity output must not alias an input"))
    end
    return nothing
end

@inline _accumulate_prepared_intensity!(output, ::Tuple{}) = output

@inline function _accumulate_prepared_intensity!(output, inputs::Tuple)
    output .+= first(inputs).values
    return _accumulate_prepared_intensity!(output, Base.tail(inputs))
end


function _accumulate_prepared_intensity!(output,
    inputs::AbstractVector{<:IntensityMap})
    @inbounds for input in inputs
        output .+= input.values
    end
    return output
end

"""
    DetectorAcquisitionPlan

Prepared contract between one frame detector and one immutable `IntensityMap`
description. Preparation validates the optical plane, radiometry, storage, and
device contract and sizes the spatial acquisition work buffers. The plan is
bound to that detector's mutable state and prepared frame storage, so it cannot
be reused by an otherwise identical detector or after moving/replacing the
detector's backend/device storage. It is also bound to the exact validated
`map.values` storage; a different array cannot reuse the plan even when wrapped
with the same metadata. Repeated `capture!` calls consume that storage without
rebuilding metadata or resizing those buffers.
"""
struct _DetectorAcquisitionPlanToken end
const _DETECTOR_ACQUISITION_PLAN_TOKEN = _DetectorAcquisitionPlanToken()

struct DetectorAcquisitionPlan{
    P<:DetectorParams,
    S<:DetectorState,
    A<:AbstractMatrix,
    I<:AbstractMatrix,
    B<:AbstractArrayBackend,
    M<:OpticalPlaneMetadata,
    T<:AbstractFloat,
}
    detector_params::P
    detector_state::S
    detector_frame::A
    detector_backend::B
    input_metadata::M
    input_values::I
    rate_scale::T
    quantum_efficiency::T

    function DetectorAcquisitionPlan(::_DetectorAcquisitionPlanToken,
        detector_params::P, detector_state::S, detector_frame::A,
        detector_backend::B, input_metadata::M, input_values::I, rate_scale::T,
        quantum_efficiency::T) where {
        P<:DetectorParams,S<:DetectorState,A<:AbstractMatrix,
        I<:AbstractMatrix,B<:AbstractArrayBackend,M<:OpticalPlaneMetadata,
        T<:AbstractFloat,
    }
        return new{P,S,A,I,B,M,T}(detector_params, detector_state,
            detector_frame, detector_backend, input_metadata, input_values,
            rate_scale, quantum_efficiency)
    end
end

@inline _require_detector_acquisition_plane(::FocalPlane) = nothing
@inline _require_detector_acquisition_plane(::DetectorPlane) = nothing

function _require_detector_acquisition_plane(::AbstractOpticalPlaneKind)
    throw(InvalidConfiguration(
        "detector acquisition requires a focal-plane or detector-plane intensity map"))
end

@inline function _normalization_rate_scale(::PhotonRateNormalization,
    ::Nothing, ::Type{T}) where {T<:AbstractFloat}
    return one(T)
end

function _normalization_rate_scale(::PhotonRateNormalization,
    normalized_to_photon_rate::Real, ::Type{T}) where {T<:AbstractFloat}
    throw(InvalidConfiguration(
        "normalized_to_photon_rate is only valid for dimensionless intensity maps"))
end

function _normalization_rate_scale(::DimensionlessNormalization,
    ::Nothing, ::Type{T}) where {T<:AbstractFloat}
    throw(InvalidConfiguration(
        "dimensionless intensity requires an explicit normalized_to_photon_rate scale"))
end

function _normalization_rate_scale(::DimensionlessNormalization,
    normalized_to_photon_rate::Real, ::Type{T}) where {T<:AbstractFloat}
    scale = T(normalized_to_photon_rate)
    isfinite(scale) && scale > zero(T) || throw(InvalidConfiguration(
        "normalized_to_photon_rate must be finite and > 0"))
    return scale
end

function _normalization_rate_scale(::AbstractOpticalNormalization,
    normalized_to_photon_rate::Union{Nothing,Real},
    ::Type{T}) where {T<:AbstractFloat}
    throw(InvalidConfiguration(
        "detector acquisition requires photon-rate or explicitly scaled dimensionless normalization"))
end

@inline _spatial_rate_scale(::CellIntegratedMeasure,
    metadata::OpticalPlaneMetadata, ::Type{T}) where {T<:AbstractFloat} = one(T)

function _spatial_rate_scale(::SpatialDensityMeasure,
    metadata::OpticalPlaneMetadata, ::Type{T}) where {T<:AbstractFloat}
    dx = T(metadata.sampling[1])
    dy = T(metadata.sampling[2])
    isfinite(dx) && dx > zero(T) && isfinite(dy) && dy > zero(T) ||
        throw(InvalidConfiguration(
            "spatial-density acquisition requires finite positive plane sampling"))
    cell_measure = dx * dy
    isfinite(cell_measure) && cell_measure > zero(T) ||
        throw(InvalidConfiguration(
            "spatial-density plane sampling has an unrepresentable cell measure"))
    return cell_measure
end

function _spatial_rate_scale(::AbstractSpatialMeasure,
    metadata::OpticalPlaneMetadata, ::Type{T}) where {T<:AbstractFloat}
    throw(InvalidConfiguration(
        "detector acquisition requires cell-integrated or spatial-density intensity samples"))
end

@inline _require_detector_incoherent(::IncoherentIntensityAddition) = nothing

function _require_detector_incoherent(::AbstractCombinationPolicy)
    throw(InvalidConfiguration(
        "detector intensity acquisition requires incoherent-intensity semantics"))
end

@inline function _prepared_quantum_efficiency(det::Detector,
    channel::MonochromaticChannel, ::Type{T}) where {T<:AbstractFloat}
    return T(qe_at(det.params.quantum_efficiency_model,
        channel.wavelength_m))
end

@inline function _prepared_quantum_efficiency(det::Detector,
    ::IntegratedSpectralChannel, ::Type{T}) where {T<:AbstractFloat}
    return _integrated_spectral_qe(det.params.quantum_efficiency_model, T)
end

@inline _integrated_spectral_qe(model::ScalarQuantumEfficiency,
    ::Type{T}) where {T<:AbstractFloat} = T(model.value)

function _integrated_spectral_qe(::SampledQuantumEfficiency,
    ::Type{T}) where {T<:AbstractFloat}
    throw(InvalidConfiguration(
        "wavelength-dependent detector QE cannot be applied after passband integration"))
end

function _prepared_quantum_efficiency(det::Detector,
    ::AbstractSpectralCoordinate, ::Type{T}) where {T<:AbstractFloat}
    throw(InvalidConfiguration(
        "detector acquisition requires a declared monochromatic or integrated spectral channel"))
end

@inline _require_prepared_response_sampling(::NullFrameResponse, ::Int) = nothing

function _require_prepared_response_sampling(::AbstractFrameResponse,
    psf_sampling::Int)
    psf_sampling == 1 || throw(InvalidConfiguration(
        "a non-null detector response currently requires psf_sampling == 1; " *
        "prepare an explicit optical-grid mapping before detector acquisition"))
    return nothing
end


function _require_finite_nonnegative_intensity(values::AbstractMatrix)
    isempty(values) && throw(InvalidConfiguration(
        "detector intensity input must not be empty"))
    minimum_value = minimum(values)
    maximum_value = maximum(values)
    isfinite(minimum_value) && isfinite(maximum_value) &&
        minimum_value >= zero(minimum_value) || throw(InvalidConfiguration(
            "detector intensity input values must be finite and nonnegative"))
    return nothing
end

"""
    prepare_detector_acquisition(detector, map;
        normalized_to_photon_rate=nothing)

Validate and prepare acquisition of an `IntensityMap`. Photon-rate maps cannot
be rescaled. A dimensionless calibration/test map must supply an explicit
`normalized_to_photon_rate` factor. Spatial-density samples are converted to a
cell rate using the declared plane sampling; cell-integrated samples are used
directly. Preparation also validates that every current sample is finite and
nonnegative. Repeated prepared captures trust subsequent producer writes to
the same storage; producers must preserve that value contract, and replacement
storage requires a new plan.
"""
function prepare_detector_acquisition(det::Detector, map::IntensityMap;
    normalized_to_photon_rate::Union{Nothing,Real}=nothing)
    require_whole_capture_idle(det)
    metadata = validate_plane_storage(map.metadata, map.values;
        label="detector intensity input")
    _require_detector_acquisition_plane(metadata.kind)
    _require_detector_incoherent(metadata.coherence)
    _require_finite_nonnegative_intensity(map.values)

    T = eltype(det.state.frame)
    metadata.numeric_type === T || throw(InvalidConfiguration(
        "detector and intensity map must use the same numeric type"))
    typeof(backend(det)) === typeof(metadata.backend) ||
        throw(InvalidConfiguration(
            "detector and intensity map must use the same array backend"))
    plane_device(det.state.frame) == metadata.device ||
        throw(InvalidConfiguration(
            "detector and intensity map must occupy the same physical device"))

    normalization_scale = _normalization_rate_scale(metadata.normalization,
        normalized_to_photon_rate, T)
    spatial_scale = _spatial_rate_scale(metadata.spatial_measure, metadata, T)
    rate_scale = normalization_scale * spatial_scale
    isfinite(rate_scale) && rate_scale > zero(T) || throw(InvalidConfiguration(
        "prepared detector photon-rate scaling is not representable"))
    _require_prepared_response_sampling(det.params.response_model,
        det.params.psf_sampling)
    quantum_efficiency = _prepared_quantum_efficiency(det,
        metadata.spectral, T)
    prepare_detector_buffers!(det, metadata.dimensions)
    return DetectorAcquisitionPlan(_DETECTOR_ACQUISITION_PLAN_TOKEN,
        det.params, det.state, det.state.frame, backend(det), metadata,
        map.values, rate_scale, quantum_efficiency)
end

@inline function _require_prepared_acquisition(det::Detector,
    map::IntensityMap, plan::DetectorAcquisitionPlan)
    det.params === plan.detector_params || throw(InvalidConfiguration(
        "detector does not match its prepared acquisition plan"))
    det.state === plan.detector_state || throw(InvalidConfiguration(
        "detector does not match its prepared acquisition plan"))
    typeof(backend(det)) === typeof(plan.detector_backend) ||
        throw(InvalidConfiguration(
            "detector backend does not match its prepared acquisition plan"))
    det.state.frame === plan.detector_frame || throw(InvalidConfiguration(
        "detector device storage does not match its prepared acquisition plan"))
    map.metadata === plan.input_metadata || throw(InvalidConfiguration(
        "intensity map does not match its prepared acquisition plan"))
    map.values === plan.input_values || throw(InvalidConfiguration(
        "intensity storage does not match its prepared acquisition plan"))
    return nothing
end

function capture!(det::Detector, map::IntensityMap,
    plan::DetectorAcquisitionPlan, rng::AbstractRNG)
    _require_prepared_acquisition(det, map, plan)
    return capture_with_quantum_efficiency!(det, map.values,
        plan.quantum_efficiency * plan.rate_scale, rng)
end

function capture!(det::Detector, map::IntensityMap,
    plan::DetectorAcquisitionPlan; rng::AbstractRNG=Random.default_rng(),
    sample_time::Union{Nothing,Real}=nothing)
    sample_time === nothing && return capture!(det, map, plan, rng)
    _require_prepared_acquisition(det, map, plan)
    return capture_incremental!(det, map.values, rng, sample_time,
        plan.quantum_efficiency * plan.rate_scale)
end

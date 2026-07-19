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
    # Preparation is allowed to synchronize. Validating through a host view
    # avoids backend-specific reduction compilation in this fallible setup
    # path and keeps repeated prepared acquisition device resident.
    minimum_value, maximum_value = extrema(host_array(values))
    isfinite(minimum_value) && isfinite(maximum_value) &&
        minimum_value >= zero(minimum_value) || throw(InvalidConfiguration(
            "detector intensity input values must be finite and nonnegative"))
    return nothing
end

@inline _require_detector_defect_shape(::NullDetectorDefectModel,
    ::Tuple{Int,Int}) = nothing

function _require_detector_defect_shape(model::PixelResponseNonuniformity,
    frame_shape::Tuple{Int,Int})
    size(model.gain_map) == frame_shape || throw(DimensionMismatchError(
        "PixelResponseNonuniformity gain_map size must match detector frame size"))
    return nothing
end

function _require_detector_defect_shape(model::DarkSignalNonuniformity,
    frame_shape::Tuple{Int,Int})
    size(model.dark_map) == frame_shape || throw(DimensionMismatchError(
        "DarkSignalNonuniformity dark_map size must match detector frame size"))
    return nothing
end

function _require_detector_defect_shape(model::BadPixelMask,
    frame_shape::Tuple{Int,Int})
    size(model.mask) == frame_shape || throw(DimensionMismatchError(
        "BadPixelMask mask size must match detector frame size"))
    return nothing
end

@inline _require_detector_defect_shapes(::Tuple{},
    ::Tuple{Int,Int}) = nothing

@inline function _require_detector_defect_shapes(models::Tuple,
    frame_shape::Tuple{Int,Int})
    _require_detector_defect_shape(first(models), frame_shape)
    _require_detector_defect_shapes(Base.tail(models), frame_shape)
    return nothing
end

@inline function _require_detector_defect_shape(
    model::CompositeDetectorDefectModel, frame_shape::Tuple{Int,Int})
    return _require_detector_defect_shapes(model.stages, frame_shape)
end

@inline _require_background_flux_shape(::NoBackground,
    ::Tuple{Int,Int}) = nothing
@inline _require_background_flux_shape(::ScalarBackground,
    ::Tuple{Int,Int}) = nothing

function _require_background_flux_shape(background::BackgroundFrame,
    frame_shape::Tuple{Int,Int})
    size(background.map) == frame_shape || throw(DimensionMismatchError(
        "background_flux size must match detector frame size"))
    return nothing
end

@inline _require_background_map_shape(::NoBackground,
    ::Tuple{Int,Int}) = nothing
@inline _require_background_map_shape(::ScalarBackground,
    ::Tuple{Int,Int}) = nothing

function _require_background_map_shape(background::BackgroundFrame,
    frame_shape::Tuple{Int,Int})
    size(background.map) == frame_shape || throw(DimensionMismatchError(
        "background_map size must match detector frame size"))
    return nothing
end

@inline _require_cmos_readout_shape(::NullCMOSReadNoise,
    ::Tuple{Int,Int}) = nothing

function _require_cmos_readout_shape(model::CMOSReadNoiseMap,
    frame_shape::Tuple{Int,Int})
    size(model.sigma) == frame_shape || throw(DimensionMismatchError(
        "CMOSReadNoiseMap sigma size must match detector frame size"))
    return nothing
end

@inline _require_cmos_output_shape(::NullCMOSOutputModel,
    ::Tuple{Int,Int}) = nothing

function _require_cmos_output_shape(model::StaticCMOSOutputPattern,
    frame_shape::Tuple{Int,Int})
    length(model.gains) * model.output_cols >= frame_shape[2] ||
        throw(DimensionMismatchError(
            "StaticCMOSOutputPattern does not cover detector columns"))
    return nothing
end

@inline _require_sensor_frame_shape(::FrameSensorType,
    ::Tuple{Int,Int}) = nothing

@inline function _require_sensor_frame_shape(sensor::CMOSSensor,
    frame_shape::Tuple{Int,Int})
    _require_cmos_readout_shape(sensor.readout_noise_model, frame_shape)
    _require_cmos_output_shape(sensor.output_model, frame_shape)
    return nothing
end

function _require_detector_output_configuration(det::Detector)
    output = det.state.output_buffer
    output_host = det.state.output_buffer_host
    requires_output = det.params.output_type !== nothing ||
        det.params.readout_window !== nothing
    if requires_output
        output === nothing && throw(InvalidConfiguration(
            "Detector output type or readout window requires allocated output storage"))
        output_host === nothing && throw(InvalidConfiguration(
            "Detector output storage requires a prepared host output buffer"))
        return nothing
    end
    output === nothing || throw(InvalidConfiguration(
        "Detector output storage requires an output type or readout window"))
    output_host === nothing || throw(InvalidConfiguration(
        "Detector host output storage requires detector output storage"))
    return nothing
end

function _require_prepared_detector_storage(det::Detector,
    frame_shape::Tuple{Int,Int}, output_shape::Tuple{Int,Int})
    size(det.state.frame) == frame_shape || throw(DimensionMismatchError(
        "detector frame storage size must match prepared frame size"))
    output = det.state.output_buffer
    output === nothing && return nothing
    size(output) == output_shape || throw(DimensionMismatchError(
        "detector output storage size must match prepared readout size"))
    output_host = det.state.output_buffer_host
    output_host === nothing && throw(InvalidConfiguration(
        "Detector output storage requires a prepared host output buffer"))
    size(output_host) == output_shape || throw(DimensionMismatchError(
        "detector host output storage size must match prepared readout size"))
    return nothing
end

@inline _require_sensor_sampling_configuration(::FrameSensorType,
    ::Tuple{Int,Int}, ::Union{Nothing,FrameWindow}, ::Real,
    ::Type{<:AbstractFloat}) = nothing

@inline function _require_sensor_sampling_configuration(
    sensor::HgCdTeAvalancheArraySensor, frame_shape::Tuple{Int,Int},
    window::Union{Nothing,FrameWindow}, exposure_time::Real,
    ::Type{T}) where {T<:AbstractFloat}
    return _require_hgcdte_sampling_configuration(sensor.sampling_mode,
        sensor, frame_shape, window, exposure_time, T)
end

@inline _require_hgcdte_sampling_configuration(::FrameSamplingMode,
    ::HgCdTeAvalancheArraySensor, ::Tuple{Int,Int},
    ::Union{Nothing,FrameWindow}, ::Real, ::Type{<:AbstractFloat}) = nothing

@inline function _require_hgcdte_sampling_configuration(
    mode::UpTheRampSampling, sensor::HgCdTeAvalancheArraySensor,
    frame_shape::Tuple{Int,Int}, window::Union{Nothing,FrameWindow},
    exposure_time::Real, ::Type{T}) where {T<:AbstractFloat}
    validate_up_the_ramp_schedule(sensor, frame_shape, window, mode,
        exposure_time, T)
    return nothing
end

@inline function _require_detector_configuration(det::Detector,
    frame_shape::Tuple{Int,Int})
    _require_detector_defect_shape(det.params.defect_model, frame_shape)
    _require_background_flux_shape(det.background_flux, frame_shape)
    _require_background_map_shape(det.background_map, frame_shape)
    _require_sensor_frame_shape(det.params.sensor, frame_shape)
    _require_detector_output_configuration(det)
    _require_sensor_sampling_configuration(det.params.sensor, frame_shape,
        det.params.readout_window, det.params.integration_time,
        eltype(det.state.frame))
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
    frame_shape = detector_frame_shape(det, metadata.dimensions)
    output_shape = detector_output_shape(det, metadata.dimensions)
    _require_detector_configuration(det, frame_shape)
    prepare_detector_buffers!(det, metadata.dimensions)
    _require_prepared_detector_storage(det, frame_shape, output_shape)
    prepare_frame_readout_state!(det.params.sensor, det)
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

@inline function _require_prepared_whole_acquisition(det::Detector,
    map::IntensityMap, plan::DetectorAcquisitionPlan)
    _require_prepared_acquisition(det, map, plan)
    require_whole_capture_idle(det)
    return nothing
end

function capture!(det::Detector, map::IntensityMap,
    plan::DetectorAcquisitionPlan, rng::AbstractRNG)
    _require_prepared_whole_acquisition(det, map, plan)
    return capture_with_quantum_efficiency!(det, map.values,
        plan.quantum_efficiency * plan.rate_scale, rng)
end

function capture!(det::Detector, map::IntensityMap,
    plan::DetectorAcquisitionPlan; rng::AbstractRNG=Random.default_rng(),
    integration_duration::Union{Nothing,Real}=nothing)
    integration_duration === nothing && return capture!(det, map, plan, rng)
    _require_prepared_acquisition(det, map, plan)
    return capture_incremental!(det, map.values, rng, integration_duration,
        plan.quantum_efficiency * plan.rate_scale)
end

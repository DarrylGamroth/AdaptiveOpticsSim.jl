"""
    HgCdTeReadout(; glow_rate=0, read_time=0,
        sampling_mode=SingleRead(), T=Float64)

Immutable readout configuration shared by conventional HgCdTe arrays and
HgCdTe linear-avalanche arrays. `read_time` is the duration of one full sensor
read. A `FrameWindow` selects output products only and does not shorten this
physical read duration.
"""
struct HgCdTeReadout{
    T<:AbstractFloat,
    M<:FrameSamplingMode,
}
    glow_rate::T
    read_time::T
    sampling_mode::M
end

function HgCdTeReadout(;
    glow_rate::Real=0.0,
    read_time::Real=0.0,
    sampling_mode::FrameSamplingMode=SingleRead(),
    T::Type{<:AbstractFloat}=Float64,
)
    glow = T(glow_rate)
    read_duration = T(read_time)
    isfinite(glow) && glow >= zero(T) || throw(InvalidConfiguration(
        "HgCdTeReadout glow_rate must be finite and >= 0"))
    isfinite(read_duration) && read_duration >= zero(T) ||
        throw(InvalidConfiguration(
            "HgCdTeReadout read_time must be finite and >= 0"))
    mode = validate_hgcdte_sampling_mode(sampling_mode)
    return HgCdTeReadout{T,typeof(mode)}(glow, read_duration, mode)
end

"""
    HgCdTeSensor(; kwargs...)

Conventional HgCdTe area sensor with nondestructive-read support and no
avalanche multiplication. Detector response, IPC, nonlinearity, persistence,
and readout correction remain explicit configuration rather than implicit
technology or camera profiles.
"""
struct HgCdTeSensor{
    R<:HgCdTeReadout,
    P<:AbstractPersistenceModel,
} <: HgCdTeSensorType
    readout::R
    persistence_model::P
end

function HgCdTeSensor(readout::HgCdTeReadout{T};
    persistence_model::AbstractPersistenceModel=NullPersistence(),
) where {T<:AbstractFloat}
    persistence = prepare_hgcdte_persistence_model(persistence_model, T)
    return HgCdTeSensor{typeof(readout),typeof(persistence)}(
        readout, persistence)
end

function HgCdTeSensor(;
    glow_rate::Real=0.0,
    read_time::Real=0.0,
    sampling_mode::FrameSamplingMode=SingleRead(),
    persistence_model::AbstractPersistenceModel=NullPersistence(),
    T::Type{<:AbstractFloat}=Float64,
)
    return HgCdTeSensor(HgCdTeReadout(
        glow_rate=glow_rate,
        read_time=read_time,
        sampling_mode=sampling_mode,
        T=T,
    ); persistence_model=persistence_model)
end

@inline function prepare_hgcdte_persistence_model(
    model::AbstractPersistenceModel, ::Type{T}) where {T<:AbstractFloat}
    return validate_persistence_model(convert_persistence_model(model, T))
end

validate_hgcdte_sampling_mode(mode::SingleRead) =
    validate_frame_sampling_mode(mode)
validate_hgcdte_sampling_mode(mode::AveragedNonDestructiveReads) =
    validate_frame_sampling_mode(mode)
validate_hgcdte_sampling_mode(mode::CorrelatedDoubleSampling) =
    validate_frame_sampling_mode(mode)
validate_hgcdte_sampling_mode(mode::FowlerSampling) =
    validate_frame_sampling_mode(mode)
validate_hgcdte_sampling_mode(mode::UpTheRampSampling) =
    validate_frame_sampling_mode(mode)
validate_hgcdte_sampling_mode(::FrameSamplingMode) =
    throw(InvalidConfiguration(
        "HgCdTe sampling_mode must be SingleRead, " *
        "AveragedNonDestructiveReads, CorrelatedDoubleSampling, " *
        "FowlerSampling, or UpTheRampSampling"))

@inline hgcdte_readout(sensor::HgCdTeSensor) = sensor.readout

detector_sensor_symbol(::HgCdTeSensor) = :hgcdte
supports_sensor_glow(::HgCdTeSensorType) = true
supports_nondestructive_reads(::HgCdTeSensorType) = true
supports_up_the_ramp(::HgCdTeSensorType) = true
supports_reference_read_subtraction(::HgCdTeSensorType) = true
supports_readout_correction(::HgCdTeSensorType) = true
supports_read_cube(::HgCdTeSensorType) = true
supports_multi_read_readout_products(::HgCdTeSensorType) = true
supports_detector_nonlinearity(::HgCdTeSensorType) = true
supports_detector_persistence(::HgCdTeSensorType) = true

default_response_model(::HgCdTeSensorType;
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=CPUBackend()) = NullFrameResponse()

@inline configured_glow_rate(sensor::HgCdTeSensorType,
    ::Type{T}) where {T<:AbstractFloat} =
    T(hgcdte_readout(sensor).glow_rate)
@inline multi_read_sampling_mode(sensor::HgCdTeSensorType) =
    hgcdte_readout(sensor).sampling_mode
@inline persistence_model(sensor::HgCdTeSensorType) =
    sensor.persistence_model

sampling_read_time(sensor::HgCdTeSensorType,
    ::Type{T}) where {T<:AbstractFloat} =
    T(hgcdte_readout(sensor).read_time)

sampling_read_time(sensor::HgCdTeSensorType, ::Tuple{Int,Int},
    ::Union{Nothing,FrameWindow}, ::Type{T}) where {T<:AbstractFloat} =
    sampling_read_time(sensor, T)

function sampling_wallclock_time(sensor::HgCdTeSensorType,
    integration_time, ::Type{T}) where {T<:AbstractFloat}
    mode = multi_read_sampling_mode(sensor)
    read_count = hgcdte_wallclock_read_count(mode)
    read_count === nothing && return nothing
    return T(integration_time) + T(read_count) * sampling_read_time(sensor, T)
end

@inline hgcdte_wallclock_read_count(mode::FrameSamplingMode) =
    frame_sampling_reads(mode)
@inline hgcdte_wallclock_read_count(::UpTheRampSampling) = 1

sampling_wallclock_time(sensor::HgCdTeSensorType, integration_time,
    ::Tuple{Int,Int}, ::Union{Nothing,FrameWindow},
    ::Type{T}) where {T<:AbstractFloat} =
    sampling_wallclock_time(sensor, integration_time, T)

effective_readout_sigma(sensor::HgCdTeSensorType, sigma) =
    effective_readout_sigma(multi_read_sampling_mode(sensor), sigma)

function effective_dark_current_time(sensor::HgCdTeSensorType, exposure_time)
    mode = multi_read_sampling_mode(sensor)
    return effective_hgcdte_generated_charge_time(
        mode, hgcdte_readout(sensor).read_time, exposure_time)
end

function effective_hgcdte_generated_charge_time(
    mode::FrameSamplingMode, read_time, exposure_time)
    reads = frame_sampling_reads(mode)
    reads === nothing && return exposure_time
    return exposure_time + reads * read_time
end

@inline effective_hgcdte_generated_charge_time(
    ::UpTheRampSampling, _read_time, exposure_time) = exposure_time

effective_sensor_glow_time(sensor::HgCdTeSensorType, exposure_time) =
    effective_dark_current_time(sensor, exposure_time)

function _apply_hgcdte_glow!(sensor::HgCdTeSensorType, det::Detector,
    rng::AbstractRNG, exposure_time::Real)
    rate = effective_glow_rate(det) *
        effective_sensor_glow_time(sensor, exposure_time)
    return add_poisson_rate!(det.products.frame, det, rng, rate)
end

apply_sensor_statistics!(sensor::HgCdTeSensor, det::Detector,
    rng::AbstractRNG, exposure_time::Real) =
    _apply_hgcdte_glow!(sensor, det, rng, exposure_time)

function apply_incremental_sensor_statistics!(sensor::HgCdTeSensorType,
    det::Detector, rng::AbstractRNG, exposure_time::Real)
    rate = effective_glow_rate(det) * exposure_time
    return add_poisson_rate!(det.products.frame, det, rng, rate)
end

function apply_hgcdte_persistence!(::NullPersistence, det::Detector,
    ::Real)
    return det.products.frame
end

function apply_hgcdte_persistence!(model::ExponentialPersistence,
    det::Detector, ::Real)
    ensure_latent_buffer!(det)
    det.products.frame .+= det.state.latent_buffer
    return det.products.frame
end

function update_hgcdte_persistence!(::NullPersistence, det::Detector,
    ::Real)
    return det.products.frame
end

function update_hgcdte_persistence!(model::ExponentialPersistence,
    det::Detector, ::Real)
    ensure_latent_buffer!(det)
    det.state.latent_buffer .= model.decay .* det.state.latent_buffer .+
        model.coupling .* det.products.frame
    return det.state.latent_buffer
end

apply_sensor_persistence!(sensor::HgCdTeSensorType, det::Detector,
    exposure_time::Real) =
    apply_hgcdte_persistence!(persistence_model(sensor), det, exposure_time)

update_sensor_persistence!(sensor::HgCdTeSensorType, det::Detector,
    exposure_time::Real) =
    update_hgcdte_persistence!(persistence_model(sensor), det, exposure_time)

function apply_post_readout_gain!(::HgCdTeSensorType, det::Detector)
    det.products.frame .*= det.params.gain
    return det.products.frame
end

function _batched_hgcdte_glow!(sensor::HgCdTeSensorType, det::Detector,
    cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG,
    exposure_time::Real)
    rate = effective_glow_rate(det) *
        effective_sensor_glow_time(sensor, exposure_time)
    rate <= zero(rate) && return cube
    fill!(scratch, rate)
    poisson_noise_frame!(det, rng, scratch)
    cube .+= scratch
    return cube
end

_batched_sensor_statistics!(sensor::HgCdTeSensor, det::Detector,
    cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG,
    exposure_time::Real) =
    _batched_hgcdte_glow!(sensor, det, cube, scratch, rng, exposure_time)

_batched_post_readout_gain!(::HgCdTeSensorType, det::Detector,
    cube::AbstractArray) = (cube .*= det.params.gain; cube)

@inline _require_batched_sensor_compat(sensor::HgCdTeSensorType) =
    require_batched_hgcdte_sampling(multi_read_sampling_mode(sensor))
@inline require_batched_hgcdte_sampling(::FrameSamplingMode) = nothing
require_batched_hgcdte_sampling(::UpTheRampSampling) =
    throw(InvalidConfiguration(
        "batched detector capture does not retain up-the-ramp read products; " *
        "use repeated capture! calls"))

function _hgcdte_readout_products_type(::UpTheRampSampling, frame::A,
    ::Type{T}) where {T<:AbstractFloat,A<:AbstractMatrix{T}}
    cube_type = typeof(similar(frame, size(frame, 1), size(frame, 2), 1))
    return Union{
        NoFrameReadoutProducts,
        UpTheRampReadoutProducts{A,cube_type,FixedSizeVectorDefault{T}},
    }
end

function _hgcdte_readout_products_type(::FrameSamplingMode, frame::A,
    ::Type{T}) where {T<:AbstractFloat,A<:AbstractMatrix{T}}
    cube_type = typeof(similar(frame, size(frame, 1), size(frame, 2), 1))
    return Union{
        NoFrameReadoutProducts,
        MultiReadFrameReadoutProducts{A,Nothing,Nothing},
        MultiReadFrameReadoutProducts{A,cube_type,Nothing},
        MultiReadFrameReadoutProducts{A,cube_type,FixedSizeVectorDefault{T}},
    }
end

function detector_readout_products_type(sensor::HgCdTeSensorType, frame::A,
    ::Type{T}) where {T<:AbstractFloat,A<:AbstractMatrix{T}}
    return _hgcdte_readout_products_type(
        multi_read_sampling_mode(sensor), frame, T)
end

function _hgcdte_readout_workspace_type(::UpTheRampSampling, frame::A,
    ::Type{T}) where {T<:AbstractFloat,A<:AbstractMatrix{T}}
    cube_type = typeof(similar(frame, size(frame, 1), size(frame, 2), 1))
    return Union{
        NoFrameReadoutWorkspace,
        UpTheRampReadoutWorkspace{A,cube_type},
    }
end

function _hgcdte_readout_workspace_type(::FrameSamplingMode, frame::A,
    ::Type{T}) where {T<:AbstractFloat,A<:AbstractMatrix{T}}
    cube_type = typeof(similar(frame, size(frame, 1), size(frame, 2), 1))
    return Union{
        NoFrameReadoutWorkspace,
        MultiReadFrameReadoutWorkspace{A,cube_type},
    }
end

function detector_readout_workspace_type(sensor::HgCdTeSensorType, frame::A,
    ::Type{T}) where {T<:AbstractFloat,A<:AbstractMatrix{T}}
    return _hgcdte_readout_workspace_type(
        multi_read_sampling_mode(sensor), frame, T)
end

@inline function prepare_frame_readout_state!(
    sensor::HgCdTeSensorType, det::Detector)
    return prepare_hgcdte_readout_state!(
        multi_read_sampling_mode(sensor), sensor, det)
end

@inline function prepare_hgcdte_readout_state!(::FrameSamplingMode,
    sensor::HgCdTeSensorType, det::Detector)
    _ensure_multi_read_products!(sensor, det)
    return det
end

@inline function prepare_hgcdte_readout_state!(mode::UpTheRampSampling,
    ::HgCdTeSensorType, det::Detector)
    ensure_up_the_ramp_products!(det, mode.n_reads)
    return det
end

function finalize_readout_products!(sensor::HgCdTeSensorType,
    det::Detector, rng::AbstractRNG, exposure_time::Real)
    return finalize_hgcdte_readout_products!(
        multi_read_sampling_mode(sensor), sensor, det, rng, exposure_time)
end

finalize_hgcdte_readout_products!(::FrameSamplingMode,
    sensor::HgCdTeSensorType, det::Detector, rng::AbstractRNG,
    ::Real) = finalize_multi_read_readout_products!(sensor, det, rng)

finalize_hgcdte_readout_products!(mode::UpTheRampSampling,
    sensor::HgCdTeSensorType, det::Detector, rng::AbstractRNG,
    exposure_time::Real) =
    finalize_up_the_ramp_readout_products!(
        mode, sensor, det, rng, exposure_time)

function _finalize_capture!(::HgCdTeSensorType, det::Detector,
    rng::AbstractRNG, exposure_time::Real)
    return finalize_hgcdte_capture!(
        det, rng, exposure_time, exposure_time)
end

function _finalize_incremental_capture!(::HgCdTeSensorType,
    det::Detector, rng::AbstractRNG, exposure_time::Real)
    return finalize_hgcdte_capture!(
        det, rng, exposure_time, zero(exposure_time))
end

function finalize_hgcdte_capture!(det::Detector, rng::AbstractRNG,
    exposure_time::Real, charge_exposure_time::Real)
    finalize_charge_generation!(det, rng, charge_exposure_time)
    finalize_charge_transport!(det, rng)
    # Multi-read products generate each raw read, including read noise,
    # conversion gain, and per-read correction. Running the generic
    # electronics stage as well would apply those effects twice.
    finalize_readout_products!(det.params.sensor, det, rng, exposure_time)
    apply_quantization!(det)
    subtract_background_map!(det.background_map, det)
    update_sensor_persistence!(det.params.sensor, det, exposure_time)
    return det.products.frame
end

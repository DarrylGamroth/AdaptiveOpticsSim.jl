"""
    HgCdTeReadout(; glow_rate=0, read_duration=0,
        sampling_mode=SingleRead(), T=Float64)

Immutable readout configuration shared by conventional HgCdTe arrays and
HgCdTe linear-avalanche arrays. `read_duration` is the duration of one full sensor
read. A `FrameWindow` selects output products only and does not shorten this
physical read duration.
"""
struct HgCdTeReadout{
    T<:AbstractFloat,
    M<:FrameSamplingMode,
}
    glow_rate::T
    read_duration::T
    sampling_mode::M
end

function HgCdTeReadout(;
    glow_rate::Real=0.0,
    read_duration::Real=0.0,
    sampling_mode::FrameSamplingMode=SingleRead(),
    T::Type{<:AbstractFloat}=Float64,
)
    glow = T(glow_rate)
    read_duration = T(read_duration)
    isfinite(glow) && glow >= zero(T) || throw(InvalidConfiguration(
        "HgCdTeReadout glow_rate must be finite and >= 0"))
    isfinite(read_duration) && read_duration >= zero(T) ||
        throw(InvalidConfiguration(
            "HgCdTeReadout read_duration must be finite and >= 0"))
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
} <: AbstractHgCdTeSensor
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
    read_duration::Real=0.0,
    sampling_mode::FrameSamplingMode=SingleRead(),
    persistence_model::AbstractPersistenceModel=NullPersistence(),
    T::Type{<:AbstractFloat}=Float64,
)
    return HgCdTeSensor(HgCdTeReadout(
        glow_rate=glow_rate,
        read_duration=read_duration,
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
supports_sensor_glow(::AbstractHgCdTeSensor) = true
supports_nondestructive_reads(::AbstractHgCdTeSensor) = true
supports_up_the_ramp(::AbstractHgCdTeSensor) = true
supports_reference_read_subtraction(::AbstractHgCdTeSensor) = true
supports_readout_correction(::AbstractHgCdTeSensor) = true
supports_read_cube(::AbstractHgCdTeSensor) = true
supports_multi_read_readout_products(::AbstractHgCdTeSensor) = true
supports_detector_nonlinearity(::AbstractHgCdTeSensor) = true
supports_detector_persistence(::AbstractHgCdTeSensor) = true

default_response_model(::AbstractHgCdTeSensor;
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=CPUBackend()) = NullFrameResponse()

@inline configured_glow_rate(sensor::AbstractHgCdTeSensor,
    ::Type{T}) where {T<:AbstractFloat} =
    T(hgcdte_readout(sensor).glow_rate)
@inline multi_read_sampling_mode(sensor::AbstractHgCdTeSensor) =
    hgcdte_readout(sensor).sampling_mode
@inline persistence_model(sensor::AbstractHgCdTeSensor) =
    sensor.persistence_model

sampling_read_duration(sensor::AbstractHgCdTeSensor,
    ::Type{T}) where {T<:AbstractFloat} =
    T(hgcdte_readout(sensor).read_duration)

sampling_read_duration(sensor::AbstractHgCdTeSensor, ::Tuple{Int,Int},
    ::Union{Nothing,FrameWindow}, ::Type{T}) where {T<:AbstractFloat} =
    sampling_read_duration(sensor, T)

function sampling_acquisition_duration(sensor::AbstractHgCdTeSensor,
    exposure_duration, ::Type{T}) where {T<:AbstractFloat}
    mode = multi_read_sampling_mode(sensor)
    read_count = hgcdte_acquisition_read_count(mode)
    read_count === nothing && return nothing
    return T(exposure_duration) + T(read_count) * sampling_read_duration(sensor, T)
end

@inline hgcdte_acquisition_read_count(mode::FrameSamplingMode) =
    frame_sampling_reads(mode)
@inline hgcdte_acquisition_read_count(::UpTheRampSampling) = 1

sampling_acquisition_duration(sensor::AbstractHgCdTeSensor, exposure_duration,
    ::Tuple{Int,Int}, ::Union{Nothing,FrameWindow},
    ::Type{T}) where {T<:AbstractFloat} =
    sampling_acquisition_duration(sensor, exposure_duration, T)

effective_readout_sigma(sensor::AbstractHgCdTeSensor, sigma) =
    effective_readout_sigma(multi_read_sampling_mode(sensor), sigma)

function effective_dark_current_duration(sensor::AbstractHgCdTeSensor, exposure_duration)
    mode = multi_read_sampling_mode(sensor)
    return effective_hgcdte_generated_charge_duration(
        mode, hgcdte_readout(sensor).read_duration, exposure_duration)
end

function effective_hgcdte_generated_charge_duration(
    mode::FrameSamplingMode, read_duration, exposure_duration)
    reads = frame_sampling_reads(mode)
    reads === nothing && return exposure_duration
    return exposure_duration + reads * read_duration
end

@inline effective_hgcdte_generated_charge_duration(
    ::UpTheRampSampling, _read_duration, exposure_duration) = exposure_duration

effective_sensor_glow_duration(sensor::AbstractHgCdTeSensor, exposure_duration) =
    effective_dark_current_duration(sensor, exposure_duration)

function _apply_hgcdte_glow!(sensor::AbstractHgCdTeSensor, det::Detector,
    rng::AbstractRNG, exposure_duration::Real)
    rate = effective_glow_rate(det) *
        effective_sensor_glow_duration(sensor, exposure_duration)
    return add_poisson_rate!(det.products.frame, det, rng, rate)
end

apply_sensor_statistics!(sensor::HgCdTeSensor, det::Detector,
    rng::AbstractRNG, exposure_duration::Real) =
    _apply_hgcdte_glow!(sensor, det, rng, exposure_duration)

function apply_incremental_sensor_statistics!(sensor::AbstractHgCdTeSensor,
    det::Detector, rng::AbstractRNG, exposure_duration::Real)
    rate = effective_glow_rate(det) * exposure_duration
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

apply_sensor_persistence!(sensor::AbstractHgCdTeSensor, det::Detector,
    exposure_duration::Real) =
    apply_hgcdte_persistence!(persistence_model(sensor), det, exposure_duration)

update_sensor_persistence!(sensor::AbstractHgCdTeSensor, det::Detector,
    exposure_duration::Real) =
    update_hgcdte_persistence!(persistence_model(sensor), det, exposure_duration)

function apply_post_readout_gain!(::AbstractHgCdTeSensor, det::Detector)
    det.products.frame .*= det.params.gain
    return det.products.frame
end

function _batched_hgcdte_glow!(sensor::AbstractHgCdTeSensor, det::Detector,
    cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG,
    exposure_duration::Real)
    rate = effective_glow_rate(det) *
        effective_sensor_glow_duration(sensor, exposure_duration)
    rate <= zero(rate) && return cube
    fill!(scratch, rate)
    poisson_noise_frame!(det, rng, scratch)
    cube .+= scratch
    return cube
end

_batched_sensor_statistics!(sensor::HgCdTeSensor, det::Detector,
    cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG,
    exposure_duration::Real) =
    _batched_hgcdte_glow!(sensor, det, cube, scratch, rng, exposure_duration)

_batched_post_readout_gain!(::AbstractHgCdTeSensor, det::Detector,
    cube::AbstractArray) = (cube .*= det.params.gain; cube)

@inline _require_batched_sensor_compat(sensor::AbstractHgCdTeSensor) =
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

function detector_readout_products_type(sensor::AbstractHgCdTeSensor, frame::A,
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

function detector_readout_workspace_type(sensor::AbstractHgCdTeSensor, frame::A,
    ::Type{T}) where {T<:AbstractFloat,A<:AbstractMatrix{T}}
    return _hgcdte_readout_workspace_type(
        multi_read_sampling_mode(sensor), frame, T)
end

@inline function prepare_frame_readout_state!(
    sensor::AbstractHgCdTeSensor, det::Detector)
    return prepare_hgcdte_readout_state!(
        multi_read_sampling_mode(sensor), sensor, det)
end

@inline function prepare_hgcdte_readout_state!(::FrameSamplingMode,
    sensor::AbstractHgCdTeSensor, det::Detector)
    _ensure_multi_read_products!(sensor, det)
    return det
end

@inline function prepare_hgcdte_readout_state!(mode::UpTheRampSampling,
    ::AbstractHgCdTeSensor, det::Detector)
    ensure_up_the_ramp_products!(det, mode.n_reads)
    return det
end

function finalize_readout_products!(sensor::AbstractHgCdTeSensor,
    det::Detector, rng::AbstractRNG, exposure_duration::Real)
    return finalize_hgcdte_readout_products!(
        multi_read_sampling_mode(sensor), sensor, det, rng, exposure_duration)
end

finalize_hgcdte_readout_products!(::FrameSamplingMode,
    sensor::AbstractHgCdTeSensor, det::Detector, rng::AbstractRNG,
    ::Real) = finalize_multi_read_readout_products!(sensor, det, rng)

finalize_hgcdte_readout_products!(mode::UpTheRampSampling,
    sensor::AbstractHgCdTeSensor, det::Detector, rng::AbstractRNG,
    exposure_duration::Real) =
    finalize_up_the_ramp_readout_products!(
        mode, sensor, det, rng, exposure_duration)

function _finalize_capture!(::AbstractHgCdTeSensor, det::Detector,
    rng::AbstractRNG, exposure_duration::Real)
    return finalize_hgcdte_capture!(
        det, rng, exposure_duration, exposure_duration)
end

function _finalize_incremental_capture!(::AbstractHgCdTeSensor,
    det::Detector, rng::AbstractRNG, exposure_duration::Real)
    return finalize_hgcdte_capture!(
        det, rng, exposure_duration, zero(exposure_duration))
end

function finalize_hgcdte_capture!(det::Detector, rng::AbstractRNG,
    exposure_duration::Real, charge_exposure_duration::Real)
    finalize_charge_generation!(det, rng, charge_exposure_duration)
    finalize_charge_transport!(det, rng)
    # Multi-read products generate each raw read, including read noise,
    # conversion gain, and per-read correction. Running the generic
    # electronics stage as well would apply those effects twice.
    finalize_readout_products!(det.params.sensor, det, rng, exposure_duration)
    apply_quantization!(det)
    subtract_background_map!(det.background_map, det)
    update_sensor_persistence!(det.params.sensor, det, exposure_duration)
    return det.products.frame
end

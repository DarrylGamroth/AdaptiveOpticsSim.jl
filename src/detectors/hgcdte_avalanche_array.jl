struct HgCdTeAvalancheArraySensor{T<:AbstractFloat,M<:FrameSamplingMode} <: HgCdTeAvalancheArraySensorType
    avalanche_gain::T
    excess_noise_factor::T
    glow_rate::T
    read_time::T
    sampling_mode::M
end

function HgCdTeAvalancheArraySensor(; avalanche_gain::Real=1.0, excess_noise_factor::Real=1.0,
    glow_rate::Real=0.0, read_time::Real=0.0, sampling_mode::FrameSamplingMode=SingleRead(),
    T::Type{<:AbstractFloat}=Float64)
    avalanche_gain >= 1 || throw(InvalidConfiguration("HgCdTeAvalancheArraySensor avalanche_gain must be >= 1"))
    excess_noise_factor >= 1 || throw(InvalidConfiguration("HgCdTeAvalancheArraySensor excess_noise_factor must be >= 1"))
    glow_rate >= 0 || throw(InvalidConfiguration("HgCdTeAvalancheArraySensor glow_rate must be >= 0"))
    read_time >= 0 || throw(InvalidConfiguration("HgCdTeAvalancheArraySensor read_time must be >= 0"))
    mode = validate_hgcdte_sampling_mode(sampling_mode)
    return HgCdTeAvalancheArraySensor{T,typeof(mode)}(
        T(avalanche_gain), T(excess_noise_factor), T(glow_rate), T(read_time), mode)
end

validate_hgcdte_sampling_mode(mode::SingleRead) = validate_frame_sampling_mode(mode)
validate_hgcdte_sampling_mode(mode::AveragedNonDestructiveReads) =
    validate_frame_sampling_mode(mode)
validate_hgcdte_sampling_mode(mode::CorrelatedDoubleSampling) =
    validate_frame_sampling_mode(mode)
validate_hgcdte_sampling_mode(mode::FowlerSampling) = validate_frame_sampling_mode(mode)
validate_hgcdte_sampling_mode(mode::UpTheRampSampling) = validate_frame_sampling_mode(mode)
validate_hgcdte_sampling_mode(::FrameSamplingMode) = throw(InvalidConfiguration(
    "HgCdTeAvalancheArraySensor sampling_mode must be SingleRead, " *
    "AveragedNonDestructiveReads, CorrelatedDoubleSampling, FowlerSampling, " *
    "or UpTheRampSampling"))

detector_sensor_symbol(::HgCdTeAvalancheArraySensor) = :hgcdte_avalanche_array
supports_sensor_glow(::HgCdTeAvalancheArraySensor) = true
supports_nondestructive_reads(::HgCdTeAvalancheArraySensorType) = true
supports_up_the_ramp(::HgCdTeAvalancheArraySensorType) = true
supports_reference_read_subtraction(::HgCdTeAvalancheArraySensorType) = true
supports_readout_correction(::HgCdTeAvalancheArraySensorType) = true
supports_read_cube(::HgCdTeAvalancheArraySensorType) = true
supports_multi_read_readout_products(::HgCdTeAvalancheArraySensorType) = true
default_response_model(::HgCdTeAvalancheArraySensor; T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=CPUBackend()) =
    NullFrameResponse()
configured_glow_rate(sensor::HgCdTeAvalancheArraySensor, ::Type{T}) where {T<:AbstractFloat} = T(sensor.glow_rate)
multi_read_sampling_mode(sensor::HgCdTeAvalancheArraySensor) = sensor.sampling_mode

sampling_read_time(sensor::HgCdTeAvalancheArraySensor, ::Type{T}) where {T<:AbstractFloat} = T(sensor.read_time)

function _window_row_fraction(frame_size::Tuple{Int,Int}, window::Union{Nothing,FrameWindow}, ::Type{T}) where {T<:AbstractFloat}
    window === nothing && return one(T)
    total_rows = frame_size[1]
    active_rows = length(window.rows)
    return T(active_rows / total_rows)
end

function sampling_read_time(sensor::HgCdTeAvalancheArraySensor, frame_size::Tuple{Int,Int}, window::Union{Nothing,FrameWindow}, ::Type{T}) where {T<:AbstractFloat}
    return T(sensor.read_time) * _window_row_fraction(frame_size, window, T)
end

function sampling_wallclock_time(sensor::HgCdTeAvalancheArraySensor, integration_time, ::Type{T}) where {T<:AbstractFloat}
    reads = frame_sampling_reads(sensor)
    reads === nothing && return nothing
    return T(integration_time) + T(reads) * T(sensor.read_time)
end

function sampling_wallclock_time(
    sensor::HgCdTeAvalancheArraySensor{<:AbstractFloat,<:UpTheRampSampling},
    integration_time, ::Type{T}) where {T<:AbstractFloat}
    return T(integration_time) + T(sensor.read_time)
end

function sampling_wallclock_time(sensor::HgCdTeAvalancheArraySensor, integration_time, frame_size::Tuple{Int,Int},
    window::Union{Nothing,FrameWindow}, ::Type{T}) where {T<:AbstractFloat}
    reads = frame_sampling_reads(sensor)
    reads === nothing && return T(integration_time)
    read_dt = sampling_read_time(sensor, frame_size, window, T)
    return T(integration_time) + T(reads) * read_dt
end

function sampling_wallclock_time(
    sensor::HgCdTeAvalancheArraySensor{<:AbstractFloat,<:UpTheRampSampling},
    integration_time, frame_size::Tuple{Int,Int},
    window::Union{Nothing,FrameWindow}, ::Type{T}) where {T<:AbstractFloat}
    return T(integration_time) + sampling_read_time(sensor, frame_size, window, T)
end

effective_readout_sigma(sensor::HgCdTeAvalancheArraySensorType, sigma) = effective_readout_sigma(multi_read_sampling_mode(sensor), sigma)

function effective_dark_current_time(sensor::HgCdTeAvalancheArraySensorType, exposure_time)
    reads = frame_sampling_reads(sensor)
    reads === nothing && return exposure_time
    return exposure_time + reads * sensor.read_time
end


effective_dark_current_time(
    ::HgCdTeAvalancheArraySensor{<:AbstractFloat,<:UpTheRampSampling},
    exposure_time) = exposure_time

effective_sensor_glow_time(sensor::HgCdTeAvalancheArraySensorType, exposure_time) =
    effective_dark_current_time(sensor, exposure_time)

function sensor_saturation_limit(sensor::HgCdTeAvalancheArraySensor, det::Detector)
    full_well = det.params.full_well
    full_well === nothing && return nothing
    return full_well / sensor.avalanche_gain
end

function apply_sensor_statistics!(sensor::HgCdTeAvalancheArraySensor,
    det::Detector, rng::AbstractRNG, exposure_time::Real)
    rate = effective_glow_rate(det) *
        effective_sensor_glow_time(sensor, exposure_time)
    add_poisson_rate!(det.state.frame, det, rng, rate)
    return apply_avalanche_excess_noise!(sensor.excess_noise_factor, det, rng)
end

function apply_incremental_sensor_statistics!(
    sensor::HgCdTeAvalancheArraySensor, det::Detector, rng::AbstractRNG,
    exposure_time::Real)
    rate = effective_glow_rate(det) * exposure_time
    return add_poisson_rate!(det.state.frame, det, rng, rate)
end

function apply_pre_readout_gain!(sensor::HgCdTeAvalancheArraySensor, det::Detector, rng::AbstractRNG)
    det.state.frame .*= sensor.avalanche_gain
    return det.state.frame
end

function apply_post_readout_gain!(::HgCdTeAvalancheArraySensor, det::Detector)
    det.state.frame .*= det.params.gain
    return det.state.frame
end

_batched_pre_readout_gain!(sensor::HgCdTeAvalancheArraySensor, det::Detector,
    cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG) =
    (cube .*= sensor.avalanche_gain; cube)

function _batched_sensor_statistics!(sensor::HgCdTeAvalancheArraySensor,
    det::Detector, cube::AbstractArray, scratch::AbstractArray,
    rng::AbstractRNG, exposure_time::Real)
    rate = effective_glow_rate(det) *
        effective_sensor_glow_time(sensor, exposure_time)
    if rate > zero(rate)
        fill!(scratch, rate)
        poisson_noise_frame!(det, rng, scratch)
        cube .+= scratch
    end
    return _batched_avalanche_excess_noise!(sensor.excess_noise_factor, cube, scratch, rng)
end

_batched_post_readout_gain!(::HgCdTeAvalancheArraySensor, det::Detector, cube::AbstractArray) = (cube .*= det.params.gain; cube)

_require_batched_sensor_compat(
    ::HgCdTeAvalancheArraySensor{<:AbstractFloat,<:UpTheRampSampling}) =
    throw(InvalidConfiguration(
        "batched detector capture does not retain up-the-ramp read products; " *
        "use repeated capture! calls"))

function detector_readout_products_type(
    ::HgCdTeAvalancheArraySensor{<:AbstractFloat,<:UpTheRampSampling},
    frame::A, ::Type{T}) where {T<:AbstractFloat,A<:AbstractMatrix{T}}
    cube_type = typeof(similar(frame, size(frame, 1), size(frame, 2), 1))
    return Union{NoFrameReadoutProducts,UpTheRampReadoutProducts{A,cube_type,Vector{T}}}
end

function finalize_readout_products!(sensor::HgCdTeAvalancheArraySensor,
    det::Detector, rng::AbstractRNG, exposure_time::Real)
    return finalize_hgcdte_readout_products!(sensor.sampling_mode, sensor, det,
        rng, exposure_time)
end

finalize_hgcdte_readout_products!(::FrameSamplingMode,
    sensor::HgCdTeAvalancheArraySensor, det::Detector, rng::AbstractRNG,
    exposure_time::Real) = finalize_multi_read_readout_products!(sensor, det, rng)

finalize_hgcdte_readout_products!(mode::UpTheRampSampling,
    sensor::HgCdTeAvalancheArraySensor, det::Detector, rng::AbstractRNG,
    exposure_time::Real) =
    finalize_up_the_ramp_readout_products!(mode, sensor, det, rng, exposure_time)

function _finalize_capture!(::HgCdTeAvalancheArraySensorType, det::Detector,
    rng::AbstractRNG, exposure_time::Real)
    return finalize_hgcdte_capture!(det, rng, exposure_time, exposure_time)
end

function _finalize_incremental_capture!(::HgCdTeAvalancheArraySensorType,
    det::Detector, rng::AbstractRNG, exposure_time::Real)
    return finalize_hgcdte_capture!(det, rng, exposure_time,
        zero(exposure_time))
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
    return det.state.frame
end

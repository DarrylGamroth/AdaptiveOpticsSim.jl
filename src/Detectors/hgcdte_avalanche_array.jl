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
    validate_frame_sampling_mode(sampling_mode)
    return HgCdTeAvalancheArraySensor{T,typeof(sampling_mode)}(
        T(avalanche_gain), T(excess_noise_factor), T(glow_rate), T(read_time), sampling_mode)
end

detector_sensor_symbol(::HgCdTeAvalancheArraySensor) = :hgcdte_avalanche_array
supports_sensor_glow(::HgCdTeAvalancheArraySensor) = true
supports_nondestructive_reads(::HgCdTeAvalancheArraySensorType) = true
supports_reference_read_subtraction(::HgCdTeAvalancheArraySensorType) = true
supports_readout_correction(::HgCdTeAvalancheArraySensorType) = true
supports_read_cube(::HgCdTeAvalancheArraySensorType) = true
supports_multi_read_readout_products(::HgCdTeAvalancheArraySensorType) = true
default_response_model(::HgCdTeAvalancheArraySensor; T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=CPUBackend()) =
    SampledFrameResponse([0.0 0.01 0.0; 0.01 0.96 0.01; 0.0 0.01 0.0]; T=T, backend=backend)
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

function sampling_wallclock_time(sensor::HgCdTeAvalancheArraySensor, integration_time, frame_size::Tuple{Int,Int},
    window::Union{Nothing,FrameWindow}, ::Type{T}) where {T<:AbstractFloat}
    reads = frame_sampling_reads(sensor)
    reads === nothing && return T(integration_time)
    read_dt = sampling_read_time(sensor, frame_size, window, T)
    return T(integration_time) + T(reads) * read_dt
end

effective_readout_sigma(sensor::HgCdTeAvalancheArraySensorType, sigma) = effective_readout_sigma(multi_read_sampling_mode(sensor), sigma)

function effective_dark_current_time(sensor::HgCdTeAvalancheArraySensorType, exposure_time)
    reads = frame_sampling_reads(sensor)
    reads === nothing && return exposure_time
    return exposure_time + reads * sensor.read_time
end

effective_sensor_glow_time(sensor::HgCdTeAvalancheArraySensorType, exposure_time) =
    effective_dark_current_time(sensor, exposure_time)

function sensor_saturation_limit(sensor::HgCdTeAvalancheArraySensor, det::Detector)
    full_well = det.params.full_well
    full_well === nothing && return nothing
    return full_well / sensor.avalanche_gain
end

function apply_sensor_statistics!(sensor::HgCdTeAvalancheArraySensor, det::Detector, rng::AbstractRNG)
    rate = effective_glow_rate(det) * effective_sensor_glow_time(sensor, det.params.integration_time)
    add_poisson_rate!(det.state.frame, det, rng, rate)
    return apply_avalanche_excess_noise!(sensor.excess_noise_factor, det, rng)
end

function apply_pre_readout_gain!(sensor::HgCdTeAvalancheArraySensor, det::Detector, rng::AbstractRNG)
    det.state.frame .*= sensor.avalanche_gain
    return det.state.frame
end

function apply_post_readout_gain!(::HgCdTeAvalancheArraySensor, det::Detector)
    det.state.frame .*= det.params.gain
    return det.state.frame
end

_batched_pre_readout_gain!(sensor::HgCdTeAvalancheArraySensor, det::Detector, cube::AbstractArray, rng::AbstractRNG) = (cube .*= sensor.avalanche_gain; cube)

function _batched_sensor_statistics!(sensor::HgCdTeAvalancheArraySensor, det::Detector, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG)
    rate = effective_glow_rate(det) * effective_sensor_glow_time(sensor, det.params.integration_time)
    if rate > zero(rate)
        fill!(scratch, rate)
        poisson_noise_frame!(det, rng, scratch)
        cube .+= scratch
    end
    return _batched_avalanche_excess_noise!(sensor.excess_noise_factor, cube, scratch, rng)
end

_batched_post_readout_gain!(::HgCdTeAvalancheArraySensor, det::Detector, cube::AbstractArray) = (cube .*= det.params.gain; cube)

finalize_readout_products!(sensor::HgCdTeAvalancheArraySensor, det::Detector, rng::AbstractRNG, exposure_time::Real) =
    finalize_multi_read_readout_products!(sensor, det, rng)

function finalize_capture!(det::Detector{N,<:DetectorParams{T,<:HgCdTeAvalancheArraySensor,<:Any,<:Any,<:Any,<:Any,<:Any},S,BF,BM},
    rng::AbstractRNG, exposure_time::Real) where {N,T,S,BF,BM}
    apply_dark_current!(det, rng, exposure_time)
    apply_saturation!(det)
    apply_sensor_statistics!(det.params.sensor, det, rng)
    apply_pre_readout_gain!(det.params.sensor, det, rng)
    finalize_readout_products!(det.params.sensor, det, rng, exposure_time)
    apply_quantization!(det)
    subtract_background_map!(det.background_map, det)
    return det.state.frame
end

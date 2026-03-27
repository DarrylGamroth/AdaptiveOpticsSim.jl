struct SingleRead <: FrameSamplingMode end

struct AveragedNonDestructiveReads <: FrameSamplingMode
    n_reads::Int
end

struct CorrelatedDoubleSampling <: FrameSamplingMode end

struct FowlerSampling <: FrameSamplingMode
    n_pairs::Int
end

struct SAPHIRASensor{T<:AbstractFloat,M<:FrameSamplingMode} <: HgCdTeAvalancheArraySensorType
    avalanche_gain::T
    excess_noise_factor::T
    glow_rate::T
    read_time::T
    sampling_mode::M
end

function SAPHIRASensor(; avalanche_gain::Real=1.0, excess_noise_factor::Real=1.0,
    glow_rate::Real=0.0, read_time::Real=0.0, sampling_mode::FrameSamplingMode=SingleRead(),
    T::Type{<:AbstractFloat}=Float64)
    avalanche_gain >= 1 || throw(InvalidConfiguration("SAPHIRASensor avalanche_gain must be >= 1"))
    excess_noise_factor >= 1 || throw(InvalidConfiguration("SAPHIRASensor excess_noise_factor must be >= 1"))
    glow_rate >= 0 || throw(InvalidConfiguration("SAPHIRASensor glow_rate must be >= 0"))
    read_time >= 0 || throw(InvalidConfiguration("SAPHIRASensor read_time must be >= 0"))
    validate_frame_sampling_mode(sampling_mode)
    return SAPHIRASensor{T,typeof(sampling_mode)}(
        T(avalanche_gain), T(excess_noise_factor), T(glow_rate), T(read_time), sampling_mode)
end

detector_sensor_symbol(::SAPHIRASensor) = :saphira
supports_sensor_glow(::SAPHIRASensor) = true
supports_nondestructive_reads(::HgCdTeAvalancheArraySensorType) = true
supports_reference_read_subtraction(::HgCdTeAvalancheArraySensorType) = true

frame_sampling_symbol(::SingleRead) = :single_read
frame_sampling_symbol(::AveragedNonDestructiveReads) = :averaged_non_destructive_reads
frame_sampling_symbol(::CorrelatedDoubleSampling) = :correlated_double_sampling
frame_sampling_symbol(::FowlerSampling) = :fowler_sampling
frame_sampling_symbol(sensor::SAPHIRASensor) = frame_sampling_symbol(sensor.sampling_mode)

frame_sampling_reads(sensor::SAPHIRASensor) = frame_sampling_reads(sensor.sampling_mode)
frame_sampling_reads(::SingleRead) = nothing
frame_sampling_reads(mode::AveragedNonDestructiveReads) = mode.n_reads
frame_sampling_reads(::CorrelatedDoubleSampling) = 2
frame_sampling_reads(mode::FowlerSampling) = 2 * mode.n_pairs

frame_sampling_reference_reads(sensor::SAPHIRASensor) = frame_sampling_reference_reads(sensor.sampling_mode)
frame_sampling_reference_reads(::SingleRead) = 0
frame_sampling_reference_reads(::AveragedNonDestructiveReads) = 0
frame_sampling_reference_reads(::CorrelatedDoubleSampling) = 1
frame_sampling_reference_reads(mode::FowlerSampling) = mode.n_pairs

frame_sampling_signal_reads(sensor::SAPHIRASensor) = frame_sampling_signal_reads(sensor.sampling_mode)
frame_sampling_signal_reads(::SingleRead) = 1
frame_sampling_signal_reads(mode::AveragedNonDestructiveReads) = mode.n_reads
frame_sampling_signal_reads(::CorrelatedDoubleSampling) = 1
frame_sampling_signal_reads(mode::FowlerSampling) = mode.n_pairs

sampling_read_time(sensor::SAPHIRASensor, ::Type{T}) where {T<:AbstractFloat} = T(sensor.read_time)

function sampling_wallclock_time(sensor::SAPHIRASensor, integration_time, ::Type{T}) where {T<:AbstractFloat}
    reads = frame_sampling_reads(sensor)
    reads === nothing && return nothing
    return T(integration_time) + T(reads) * T(sensor.read_time)
end

effective_readout_sigma(sensor::HgCdTeAvalancheArraySensorType, sigma) = effective_readout_sigma(sensor.sampling_mode, sigma)
effective_readout_sigma(::SingleRead, sigma) = sigma
effective_readout_sigma(mode::AveragedNonDestructiveReads, sigma) = sigma / sqrt(mode.n_reads)
effective_readout_sigma(::CorrelatedDoubleSampling, sigma) = sigma * sqrt(2)
effective_readout_sigma(mode::FowlerSampling, sigma) = sigma * sqrt(2 / mode.n_pairs)

function effective_dark_current_time(sensor::HgCdTeAvalancheArraySensorType, exposure_time)
    reads = frame_sampling_reads(sensor)
    reads === nothing && return exposure_time
    return exposure_time + reads * sensor.read_time
end

effective_sensor_glow_time(sensor::HgCdTeAvalancheArraySensorType, exposure_time) =
    effective_dark_current_time(sensor, exposure_time)

validate_frame_sampling_mode(::SingleRead) = SingleRead()

function validate_frame_sampling_mode(mode::AveragedNonDestructiveReads)
    mode.n_reads >= 1 || throw(InvalidConfiguration("AveragedNonDestructiveReads n_reads must be >= 1"))
    return mode
end

validate_frame_sampling_mode(::CorrelatedDoubleSampling) = CorrelatedDoubleSampling()

function validate_frame_sampling_mode(mode::FowlerSampling)
    mode.n_pairs >= 1 || throw(InvalidConfiguration("FowlerSampling n_pairs must be >= 1"))
    return mode
end

function sensor_saturation_limit(sensor::SAPHIRASensor, det::Detector)
    full_well = det.params.full_well
    full_well === nothing && return nothing
    return full_well / sensor.avalanche_gain
end

function apply_sensor_statistics!(sensor::SAPHIRASensor, det::Detector, rng::AbstractRNG)
    rate = sensor.glow_rate * effective_sensor_glow_time(sensor, det.params.integration_time)
    if rate > zero(rate)
        fill!(det.state.noise_buffer, rate)
        poisson_noise!(rng, det.state.noise_buffer)
        det.state.frame .+= det.state.noise_buffer
    end
    return apply_avalanche_excess_noise!(sensor.excess_noise_factor, det, rng)
end

function apply_pre_readout_gain!(sensor::SAPHIRASensor, det::Detector)
    det.state.frame .*= sensor.avalanche_gain
    return det.state.frame
end

function apply_post_readout_gain!(::SAPHIRASensor, det::Detector)
    det.state.frame .*= det.params.gain
    return det.state.frame
end

_batched_pre_readout_gain!(sensor::SAPHIRASensor, det::Detector, cube::AbstractArray) = (cube .*= sensor.avalanche_gain; cube)

function _batched_sensor_statistics!(sensor::SAPHIRASensor, det::Detector, cube::AbstractArray, scratch::AbstractArray, rng::AbstractRNG)
    rate = sensor.glow_rate * effective_sensor_glow_time(sensor, det.params.integration_time)
    if rate > zero(rate)
        fill!(scratch, rate)
        poisson_noise!(rng, scratch)
        cube .+= scratch
    end
    return _batched_avalanche_excess_noise!(sensor.excess_noise_factor, cube, scratch, rng)
end

_batched_post_readout_gain!(::SAPHIRASensor, det::Detector, cube::AbstractArray) = (cube .*= det.params.gain; cube)

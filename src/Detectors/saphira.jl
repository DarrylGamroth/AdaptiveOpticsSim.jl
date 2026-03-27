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
supports_readout_correction(::HgCdTeAvalancheArraySensorType) = true
supports_read_cube(::HgCdTeAvalancheArraySensorType) = true

frame_sampling_symbol(::SingleRead) = :single_read
frame_sampling_symbol(::AveragedNonDestructiveReads) = :averaged_non_destructive_reads
frame_sampling_symbol(::CorrelatedDoubleSampling) = :correlated_double_sampling
frame_sampling_symbol(::FowlerSampling) = :fowler_sampling
frame_sampling_symbol(sensor::SAPHIRASensor) = frame_sampling_symbol(sensor.sampling_mode)

frame_sampling_reads(sensor::SAPHIRASensor) = frame_sampling_reads(sensor.sampling_mode)
frame_sampling_reads(::SingleRead) = 1
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

function _window_row_fraction(frame_size::Tuple{Int,Int}, window::Union{Nothing,FrameWindow}, ::Type{T}) where {T<:AbstractFloat}
    window === nothing && return one(T)
    total_rows = frame_size[1]
    active_rows = length(window.rows)
    return T(active_rows / total_rows)
end

function sampling_read_time(sensor::SAPHIRASensor, frame_size::Tuple{Int,Int}, window::Union{Nothing,FrameWindow}, ::Type{T}) where {T<:AbstractFloat}
    return T(sensor.read_time) * _window_row_fraction(frame_size, window, T)
end

function sampling_wallclock_time(sensor::SAPHIRASensor, integration_time, ::Type{T}) where {T<:AbstractFloat}
    reads = frame_sampling_reads(sensor)
    reads === nothing && return nothing
    return T(integration_time) + T(reads) * T(sensor.read_time)
end

function sampling_wallclock_time(sensor::SAPHIRASensor, integration_time, frame_size::Tuple{Int,Int},
    window::Union{Nothing,FrameWindow}, ::Type{T}) where {T<:AbstractFloat}
    reads = frame_sampling_reads(sensor)
    reads === nothing && return T(integration_time)
    read_dt = sampling_read_time(sensor, frame_size, window, T)
    return T(integration_time) + T(reads) * read_dt
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

function _sample_frame_read!(sensor::SAPHIRASensor, det::Detector, target::AbstractMatrix, baseline::AbstractMatrix,
    sigma, rng::AbstractRNG)
    copyto!(target, baseline)
    if sigma > zero(sigma)
        randn_backend!(rng, det.state.noise_buffer)
        target .+= sigma .* det.state.noise_buffer
    end
    target .*= det.params.gain
    apply_readout_correction!(det.params.correction_model, target)
    return target
end

function _raw_sampling_sigma(det::Detector{<:NoiseReadout})
    return det.noise.sigma
end

function _raw_sampling_sigma(det::Detector{<:NoisePhotonReadout})
    return det.noise.sigma
end

_raw_sampling_sigma(det::Detector) = zero(eltype(det.state.frame))

function _sampling_average_sigma(sensor::SAPHIRASensor, reads::Int, sigma)
    reads <= 0 && return zero(sigma)
    return sigma / sqrt(reads)
end

function _sampling_average_cube(cube::AbstractArray{T,3}) where {T}
    n_reads = size(cube, 3)
    out = similar(cube, size(cube, 1), size(cube, 2))
    fill!(out, zero(T))
    for read_idx in 1:n_reads
        @views out .+= cube[:, :, read_idx]
    end
    out ./= T(n_reads)
    return out
end

function _sampling_reference_cube(sensor::SAPHIRASensor, det::Detector, sigma, rng::AbstractRNG)
    n_ref = frame_sampling_reference_reads(sensor)
    n_ref <= 0 && return nothing
    cube = similar(det.state.frame, size(det.state.frame)..., n_ref)
    baseline = similar(det.state.frame, size(det.state.frame)...)
    fill!(baseline, zero(eltype(baseline)))
    for read_idx in 1:n_ref
        _sample_frame_read!(sensor, det, baseline, baseline, sigma, rng)
        @views copyto!(cube[:, :, read_idx], baseline)
    end
    return cube
end

function _sampling_signal_cube(sensor::SAPHIRASensor, det::Detector, sigma, rng::AbstractRNG)
    n_sig = frame_sampling_signal_reads(sensor)
    cube = similar(det.state.frame, size(det.state.frame)..., n_sig)
    baseline = similar(det.state.frame, size(det.state.frame)...)
    for read_idx in 1:n_sig
        _sample_frame_read!(sensor, det, baseline, det.state.frame, sigma, rng)
        @views copyto!(cube[:, :, read_idx], baseline)
    end
    return cube
end

function _sampling_read_cube(sensor::SAPHIRASensor, reference_cube::Union{Nothing,AbstractArray{T,3}},
    signal_cube::AbstractArray{T,3}) where {T}
    n_reads = frame_sampling_reads(sensor)
    n_reads <= 1 && return nothing
    cube = similar(signal_cube, size(signal_cube, 1), size(signal_cube, 2), n_reads)
    offset = 0
    if !isnothing(reference_cube)
        for read_idx in axes(reference_cube, 3)
            @views copyto!(cube[:, :, read_idx], reference_cube[:, :, read_idx])
        end
        offset = size(reference_cube, 3)
    end
    for read_idx in axes(signal_cube, 3)
        @views copyto!(cube[:, :, offset + read_idx], signal_cube[:, :, read_idx])
    end
    return cube
end

function _sampling_read_times(sensor::SAPHIRASensor, det::Detector, n_reads::Int)
    n_reads <= 0 && return nothing
    T = eltype(det.state.frame)
    read_dt = sampling_read_time(sensor, size(det.state.frame), det.params.readout_window, T)
    times = Vector{T}(undef, n_reads)
    for read_idx in 1:n_reads
        times[read_idx] = T(read_idx) * read_dt
    end
    return times
end

function finalize_readout_products!(sensor::SAPHIRASensor, det::Detector, rng::AbstractRNG, exposure_time::Real)
    sigma = _raw_sampling_sigma(det)
    reference_cube_full = _sampling_reference_cube(sensor, det, sigma, rng)
    signal_cube_full = _sampling_signal_cube(sensor, det, sigma, rng)
    cube_full = _sampling_read_cube(sensor, reference_cube_full, signal_cube_full)
    reference_full = isnothing(reference_cube_full) ? nothing : _sampling_average_cube(reference_cube_full)
    signal_full = _sampling_average_cube(signal_cube_full)
    combined_full = isnothing(reference_full) ? signal_full : signal_full .- reference_full
    reference = isnothing(reference_full) ? nothing : _copy_windowed_frame(reference_full, det)
    signal = _copy_windowed_frame(signal_full, det)
    combined = _copy_windowed_frame(combined_full, det)
    reference_cube = isnothing(reference_cube_full) ? nothing : _copy_windowed_cube(reference_cube_full, det)
    signal_cube = _copy_windowed_cube(signal_cube_full, det)
    cube = isnothing(cube_full) ? nothing : _copy_windowed_cube(cube_full, det)
    read_times = isnothing(cube_full) ? nothing : _sampling_read_times(sensor, det, size(cube_full, 3))
    det.state.frame .= combined_full
    det.state.readout_products = HgCdTeReadoutProducts(reference, signal, combined, reference_cube, signal_cube, cube, read_times)
    return det.state.readout_products
end

function finalize_capture!(det::Detector{N,<:DetectorParams{T,<:SAPHIRASensor},S,BF,BM},
    rng::AbstractRNG, exposure_time::Real) where {N,T,S,BF,BM}
    apply_dark_current!(det, rng, exposure_time)
    apply_saturation!(det)
    apply_sensor_statistics!(det.params.sensor, det, rng)
    apply_pre_readout_gain!(det.params.sensor, det)
    finalize_readout_products!(det.params.sensor, det, rng, exposure_time)
    apply_quantization!(det)
    subtract_background_map!(det.background_map, det)
    return det.state.frame
end

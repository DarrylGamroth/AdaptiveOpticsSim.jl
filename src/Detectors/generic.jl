readout_ready(det::Detector) = det.state.readout_ready
output_frame(det::Detector) = det.state.output_buffer === nothing ? det.state.frame : det.state.output_buffer
readout_ready(det::APDDetector) = true
channel_output(det::APDDetector) = det.state.output_buffer === nothing ? det.state.channels : det.state.output_buffer
output_frame(det::APDDetector) = channel_output(det)
reset_integration!(det::APDDetector) = det

detector_noise_symbol(::NoiseNone) = :none
detector_noise_symbol(::NoisePhoton) = :photon
detector_noise_symbol(::NoiseReadout) = :readout
detector_noise_symbol(::NoisePhotonReadout) = :photon_readout

detector_sensor_symbol(::CCDSensor) = :ccd
detector_sensor_symbol(::CMOSSensor) = :cmos
detector_sensor_symbol(::EMCCDSensor) = :emccd
detector_sensor_symbol(::InGaAsSensor) = :ingaas
detector_sensor_symbol(::SAPHIRASensor) = :saphira
detector_sensor_symbol(::APDSensor) = :apd

counting_dead_time_symbol(::NoDeadTime) = :none
counting_dead_time_symbol(::NonParalyzableDeadTime) = :nonparalyzable

frame_response_symbol(::NullFrameResponse) = :none
frame_response_symbol(::SeparableGaussianPixelResponse) = :separable_gaussian

frame_sampling_symbol(::SingleRead) = :single_read
frame_sampling_symbol(::AveragedNonDestructiveReads) = :averaged_non_destructive_reads
frame_sampling_symbol(::CorrelatedDoubleSampling) = :correlated_double_sampling
frame_sampling_symbol(::FowlerSampling) = :fowler_sampling
frame_sampling_symbol(::FrameSensorType) = :single_read
frame_sampling_symbol(sensor::SAPHIRASensor) = frame_sampling_symbol(sensor.sampling_mode)

frame_sampling_reads(::FrameSensorType) = nothing
frame_sampling_reads(sensor::SAPHIRASensor) = frame_sampling_reads(sensor.sampling_mode)
frame_sampling_reads(::SingleRead) = nothing
frame_sampling_reads(mode::AveragedNonDestructiveReads) = mode.n_reads
frame_sampling_reads(::CorrelatedDoubleSampling) = 2
frame_sampling_reads(mode::FowlerSampling) = 2 * mode.n_pairs

frame_sampling_reference_reads(::FrameSensorType) = nothing
frame_sampling_reference_reads(sensor::SAPHIRASensor) = frame_sampling_reference_reads(sensor.sampling_mode)
frame_sampling_reference_reads(::SingleRead) = 0
frame_sampling_reference_reads(::AveragedNonDestructiveReads) = 0
frame_sampling_reference_reads(::CorrelatedDoubleSampling) = 1
frame_sampling_reference_reads(mode::FowlerSampling) = mode.n_pairs

frame_sampling_signal_reads(::FrameSensorType) = nothing
frame_sampling_signal_reads(sensor::SAPHIRASensor) = frame_sampling_signal_reads(sensor.sampling_mode)
frame_sampling_signal_reads(::SingleRead) = 1
frame_sampling_signal_reads(mode::AveragedNonDestructiveReads) = mode.n_reads
frame_sampling_signal_reads(::CorrelatedDoubleSampling) = 1
frame_sampling_signal_reads(mode::FowlerSampling) = mode.n_pairs

sampling_read_time(::FrameSensorType, ::Type{T}) where {T<:AbstractFloat} = nothing
sampling_read_time(sensor::SAPHIRASensor, ::Type{T}) where {T<:AbstractFloat} = T(sensor.read_time)
sampling_wallclock_time(::FrameSensorType, integration_time, ::Type{T}) where {T<:AbstractFloat} = nothing

function sampling_wallclock_time(sensor::SAPHIRASensor, integration_time, ::Type{T}) where {T<:AbstractFloat}
    reads = frame_sampling_reads(sensor)
    reads === nothing && return nothing
    return T(integration_time) + T(reads) * T(sensor.read_time)
end

supports_detector_mtf(::AbstractFrameDetector) = false
supports_detector_mtf(det::Detector) = supports_detector_mtf(det.params.response_model)
supports_detector_mtf(::NullFrameResponse) = false
supports_detector_mtf(::SeparableGaussianPixelResponse) = true

supports_clock_induced_charge(::FrameSensorType) = false
supports_clock_induced_charge(::CCDSensor) = true

supports_column_readout_noise(::FrameSensorType) = false
supports_column_readout_noise(::CMOSSensor) = true

supports_avalanche_gain(::FrameSensorType) = false
supports_avalanche_gain(::AvalancheFrameSensorType) = true

supports_sensor_glow(::FrameSensorType) = false
supports_sensor_glow(::InGaAsSensor) = true
supports_sensor_glow(::SAPHIRASensor) = true

supports_nondestructive_reads(::FrameSensorType) = false
supports_nondestructive_reads(::HgCdTeAvalancheArraySensorType) = true

supports_reference_read_subtraction(::FrameSensorType) = false
supports_reference_read_subtraction(::HgCdTeAvalancheArraySensorType) = true

supports_counting_noise(::AbstractCountingDetector) = false
supports_counting_noise(::APDDetector) = true
supports_dead_time(::AbstractCountingDetector) = false
supports_dead_time(det::APDDetector) = supports_dead_time(det.params.dead_time_model)
supports_dead_time(::NoDeadTime) = false
supports_dead_time(::NonParalyzableDeadTime) = true
supports_channel_gain_map(::AbstractCountingDetector) = false
supports_channel_gain_map(det::APDDetector) = !isnothing(det.channel_gain_map)

sensor_is_frame_based(::SensorType) = false
sensor_is_frame_based(::FrameSensorType) = true
sensor_is_counting(::SensorType) = false
sensor_is_counting(::CountingSensorType) = true

detector_readout_sigma(::NoiseNone, ::FrameSensorType, ::Type{T}) where {T<:AbstractFloat} = nothing
detector_readout_sigma(::NoisePhoton, ::FrameSensorType, ::Type{T}) where {T<:AbstractFloat} = nothing
detector_readout_sigma(noise::NoiseReadout, sensor::FrameSensorType, ::Type{T}) where {T<:AbstractFloat} =
    T(effective_readout_sigma(sensor, T(noise.sigma)))
detector_readout_sigma(noise::NoisePhotonReadout, sensor::FrameSensorType, ::Type{T}) where {T<:AbstractFloat} =
    T(effective_readout_sigma(sensor, T(noise.sigma)))

function detector_export_metadata(det::Detector; T::Type{<:AbstractFloat}=eltype(det.state.frame))
    output = output_frame(det)
    full_well = det.params.full_well === nothing ? nothing : T(det.params.full_well)
    row_window = det.params.readout_window === nothing ? nothing :
        (first(det.params.readout_window.rows), last(det.params.readout_window.rows))
    col_window = det.params.readout_window === nothing ? nothing :
        (first(det.params.readout_window.cols), last(det.params.readout_window.cols))
    return DetectorExportMetadata{T}(
        T(det.params.integration_time),
        T(det.params.qe),
        det.params.psf_sampling,
        det.params.binning,
        T(det.params.gain),
        T(det.params.dark_current),
        det.params.bits,
        full_well,
        detector_sensor_symbol(det.params.sensor),
        detector_noise_symbol(det.noise),
        detector_readout_sigma(det.noise, det.params.sensor, T),
        det.params.output_precision,
        size(det.state.frame),
        size(output),
        frame_response_symbol(det.params.response_model),
        frame_response_width(det.params.response_model, T),
        row_window,
        col_window,
        frame_sampling_symbol(det.params.sensor),
        frame_sampling_reads(det.params.sensor),
        frame_sampling_reference_reads(det.params.sensor),
        frame_sampling_signal_reads(det.params.sensor),
        sampling_read_time(det.params.sensor, T),
        sampling_wallclock_time(det.params.sensor, det.params.integration_time, T),
    )
end

function detector_export_metadata(det::APDDetector; T::Type{<:AbstractFloat}=eltype(det.state.channels))
    output = channel_output(det)
    return CountingDetectorExportMetadata{T}(
        T(det.params.integration_time),
        T(det.params.qe),
        T(det.params.gain),
        T(det.params.dark_count_rate),
        counting_dead_time_symbol(det.params.dead_time_model),
        counting_dead_time_value(det.params.dead_time_model, T),
        detector_sensor_symbol(det.params.sensor),
        detector_noise_symbol(det.noise),
        det.params.output_precision,
        CountingReadoutMetadata(det.params.layout, size(output), length(output)),
    )
end

function reset_integration!(det::Detector)
    fill!(det.state.accum_buffer, zero(eltype(det.state.accum_buffer)))
    det.state.integrated_time = zero(det.state.integrated_time)
    det.state.readout_ready = true
    return det
end

normalize_noise(noise::NoiseModel) = noise

function normalize_noise(noises::Tuple{Vararg{NoiseModel}})
    isempty(noises) && return NoiseNone()
    has_photon = false
    has_readout = false
    sigma = 0.0
    sigma_set = false

    for noise in noises
        has_photon, has_readout, sigma, sigma_set = accumulate_noise_components(
            noise, has_photon, has_readout, sigma, sigma_set)
    end

    if has_photon && has_readout
        return NoisePhotonReadout(sigma)
    elseif has_photon
        return NoisePhoton()
    elseif has_readout
        return NoiseReadout(sigma)
    end
    return NoiseNone()
end

accumulate_noise_components(::NoiseNone, has_photon, has_readout, sigma, sigma_set) =
    (has_photon, has_readout, sigma, sigma_set)
accumulate_noise_components(::NoisePhoton, has_photon, has_readout, sigma, sigma_set) =
    (true, has_readout, sigma, sigma_set)
function accumulate_noise_components(noise::NoiseReadout, has_photon, has_readout, sigma, sigma_set)
    sigma_merged, sigma_is_set = merge_readout_sigma(sigma, sigma_set, noise.sigma)
    return has_photon, true, sigma_merged, sigma_is_set
end
function accumulate_noise_components(noise::NoisePhotonReadout, has_photon, has_readout, sigma, sigma_set)
    sigma_merged, sigma_is_set = merge_readout_sigma(sigma, sigma_set, noise.sigma)
    return true, true, sigma_merged, sigma_is_set
end

function merge_readout_sigma(current::Float64, set::Bool, new_sigma::Real)
    σ = float(new_sigma)
    if set && σ != current
        throw(InvalidConfiguration("conflicting readout noise values in noise tuple"))
    end
    return σ, true
end

convert_noise(noise::NoiseNone, ::Type{T}) where {T<:AbstractFloat} = NoiseNone()
convert_noise(noise::NoisePhoton, ::Type{T}) where {T<:AbstractFloat} = NoisePhoton()
convert_noise(noise::NoiseReadout, ::Type{T}) where {T<:AbstractFloat} = NoiseReadout{T}(T(noise.sigma))
convert_noise(noise::NoisePhotonReadout, ::Type{T}) where {T<:AbstractFloat} = NoisePhotonReadout{T}(T(noise.sigma))

validate_noise(noise::NoiseNone) = noise
validate_noise(noise::NoisePhoton) = noise

function validate_noise(noise::NoiseReadout)
    if noise.sigma <= 0
        throw(InvalidConfiguration("readout noise must be > 0 for NoiseReadout"))
    end
    return noise
end

function validate_noise(noise::NoisePhotonReadout)
    if noise.sigma <= 0
        throw(InvalidConfiguration("readout noise must be > 0 for NoisePhotonReadout"))
    end
    return noise
end

background_model(::Nothing; T::Type{<:AbstractFloat}, backend) = NoBackground()
background_model(level::Real; T::Type{<:AbstractFloat}, backend) = ScalarBackground{T}(T(level))

function background_model(map::AbstractMatrix; T::Type{<:AbstractFloat}, backend)
    background = backend{T}(undef, size(map)...)
    copyto!(background, T.(map))
    return BackgroundFrame{T, typeof(background)}(background)
end

counting_dead_time_value(::NoDeadTime, ::Type{T}) where {T<:AbstractFloat} = nothing
counting_dead_time_value(model::NonParalyzableDeadTime, ::Type{T}) where {T<:AbstractFloat} = T(model.dead_time)
frame_response_width(::NullFrameResponse, ::Type{T}) where {T<:AbstractFloat} = nothing
frame_response_width(model::SeparableGaussianPixelResponse, ::Type{T}) where {T<:AbstractFloat} = T(model.response_width_px)

effective_readout_sigma(::FrameSensorType, sigma) = sigma
effective_readout_sigma(sensor::HgCdTeAvalancheArraySensorType, sigma) = effective_readout_sigma(sensor.sampling_mode, sigma)
effective_readout_sigma(::SingleRead, sigma) = sigma
effective_readout_sigma(mode::AveragedNonDestructiveReads, sigma) = sigma / sqrt(mode.n_reads)
effective_readout_sigma(::CorrelatedDoubleSampling, sigma) = sigma * sqrt(2)
effective_readout_sigma(mode::FowlerSampling, sigma) = sigma * sqrt(2 / mode.n_pairs)

effective_dark_current_time(::FrameSensorType, exposure_time) = exposure_time
function effective_dark_current_time(sensor::HgCdTeAvalancheArraySensorType, exposure_time)
    reads = frame_sampling_reads(sensor)
    reads === nothing && return exposure_time
    return exposure_time + reads * sensor.read_time
end

effective_sensor_glow_time(::FrameSensorType, exposure_time) = exposure_time
effective_sensor_glow_time(sensor::HgCdTeAvalancheArraySensorType, exposure_time) =
    effective_dark_current_time(sensor, exposure_time)

convert_dead_time_model(::NoDeadTime, ::Type{T}) where {T<:AbstractFloat} = NoDeadTime()
convert_dead_time_model(model::NonParalyzableDeadTime, ::Type{T}) where {T<:AbstractFloat} =
    NonParalyzableDeadTime{T}(T(model.dead_time))
convert_frame_response_model(::NullFrameResponse, ::Type{T}, backend) where {T<:AbstractFloat} = NullFrameResponse()

function convert_frame_response_model(model::SeparableGaussianPixelResponse, ::Type{T}, backend) where {T<:AbstractFloat}
    kernel = backend{T}(undef, length(model.kernel))
    copyto!(kernel, T.(Array(model.kernel)))
    return SeparableGaussianPixelResponse{T,typeof(kernel)}(T(model.response_width_px), kernel)
end

validate_dead_time_model(model::NoDeadTime) = model

function validate_dead_time_model(model::NonParalyzableDeadTime)
    model.dead_time >= 0 || throw(InvalidConfiguration("NonParalyzableDeadTime dead_time must be >= 0"))
    return model
end

validate_frame_response_model(model::NullFrameResponse) = model

function validate_frame_response_model(model::SeparableGaussianPixelResponse)
    model.response_width_px > 0 || throw(InvalidConfiguration("SeparableGaussianPixelResponse response_width_px must be > 0"))
    length(model.kernel) > 0 || throw(InvalidConfiguration("SeparableGaussianPixelResponse kernel must not be empty"))
    isodd(length(model.kernel)) || throw(InvalidConfiguration("SeparableGaussianPixelResponse kernel length must be odd"))
    return model
end

validate_readout_window(::Nothing) = nothing

function validate_readout_window(window::FrameWindow)
    isempty(window.rows) && throw(InvalidConfiguration("FrameWindow rows must not be empty"))
    isempty(window.cols) && throw(InvalidConfiguration("FrameWindow cols must not be empty"))
    first(window.rows) >= 1 || throw(InvalidConfiguration("FrameWindow rows must start at >= 1"))
    first(window.cols) >= 1 || throw(InvalidConfiguration("FrameWindow cols must start at >= 1"))
    return window
end

function validate_readout_window(window::FrameWindow, n_out::Int, m_out::Int)
    validate_readout_window(window)
    last(window.rows) <= n_out || throw(DimensionMismatchError("FrameWindow rows must lie within detector output rows"))
    last(window.cols) <= m_out || throw(DimensionMismatchError("FrameWindow cols must lie within detector output cols"))
    return window
end

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

validate_frame_detector_sensor(::FrameSensorType) = nothing

function validate_frame_detector_sensor(sensor::CountingSensorType)
    throw(InvalidConfiguration("Detector currently models frame-based sensors only; use a sensor-side readout model for $(detector_sensor_symbol(sensor))"))
end

function resolve_output_precision(bits::Union{Nothing,Int}, output_precision::Union{Nothing,DataType})
    output_precision !== nothing && return output_precision
    bits === 8 && return UInt8
    bits === 16 && return UInt16
    bits === 32 && return UInt32
    bits === 64 && return UInt64
    return nothing
end

function _build_detector(noise::NoiseModel; integration_time::Real, qe::Real,
    psf_sampling::Int, binning::Int, gain::Real, dark_current::Real,
    bits::Union{Nothing,Int}, full_well::Union{Nothing,Real}, sensor::SensorType,
    response_model::FrameResponseModel, readout_window::Union{Nothing,FrameWindow},
    output_precision::Union{Nothing,DataType}, background_flux, background_map,
    T::Type{<:AbstractFloat}, backend)
    validate_frame_detector_sensor(sensor)
    full_well_t = full_well === nothing ? nothing : T(full_well)
    flux_model = background_model(background_flux; T=T, backend=backend)
    map_model = background_model(background_map; T=T, backend=backend)
    response = validate_frame_response_model(convert_frame_response_model(response_model, T, backend))
    output_precision_t = resolve_output_precision(bits, output_precision)
    window = validate_readout_window(readout_window)
    params = DetectorParams{T, typeof(sensor), typeof(response)}(
        T(integration_time),
        T(qe),
        psf_sampling,
        binning,
        T(gain),
        T(dark_current),
        bits,
        full_well_t,
        sensor,
        response,
        window,
        output_precision_t,
    )
    frame = backend{T}(undef, 1, 1)
    response_buffer = backend{T}(undef, 1, 1)
    bin_buffer = backend{T}(undef, 1, 1)
    noise_buffer = backend{T}(undef, 1, 1)
    accum_buffer = backend{T}(undef, 1, 1)
    output_buffer = if output_precision_t === nothing
        window === nothing ? nothing : backend{T}(undef, 1, 1)
    else
        backend{output_precision_t}(undef, 1, 1)
    end
    fill!(frame, zero(T))
    fill!(response_buffer, zero(T))
    fill!(bin_buffer, zero(T))
    fill!(noise_buffer, zero(T))
    fill!(accum_buffer, zero(T))
    output_buffer === nothing || fill!(output_buffer, zero(eltype(output_buffer)))
    state = DetectorState{T, typeof(frame), typeof(output_buffer)}(
        frame,
        response_buffer,
        bin_buffer,
        noise_buffer,
        accum_buffer,
        output_buffer,
        zero(T),
        true,
    )
    return Detector{typeof(noise), typeof(params), typeof(state), typeof(flux_model), typeof(map_model)}(
        noise,
        params,
        state,
        flux_model,
        map_model,
    )
end

validate_apd_noise(noise::NoiseNone) = noise
validate_apd_noise(noise::NoisePhoton) = noise
validate_apd_noise(noise::NoiseReadout) =
    throw(InvalidConfiguration("APDDetector does not support additive readout noise; use NoiseNone or NoisePhoton"))
validate_apd_noise(noise::NoisePhotonReadout) =
    throw(InvalidConfiguration("APDDetector does not support frame-style readout noise; use NoisePhoton"))

function _build_apd_detector(noise::NoiseModel; integration_time::Real, qe::Real, gain::Real,
    dark_count_rate::Real, dead_time_model::CountingDeadTimeModel, sensor::APDSensor, output_precision::Union{Nothing,DataType},
    layout::Symbol, channel_gain_map,
    T::Type{<:AbstractFloat}, backend)
    gain >= 0 || throw(InvalidConfiguration("APDDetector gain must be >= 0"))
    dark_count_rate >= 0 || throw(InvalidConfiguration("APDDetector dark_count_rate must be >= 0"))
    converted = convert_noise(noise, T)
    validated = validate_apd_noise(converted)
    dead_time = validate_dead_time_model(convert_dead_time_model(dead_time_model, T))
    gain_map = channel_gain_map === nothing ? nothing : begin
        g = backend{T}(undef, size(channel_gain_map)...)
        copyto!(g, T.(channel_gain_map))
        g
    end
    params = APDDetectorParams{T,typeof(sensor),typeof(dead_time)}(
        T(integration_time),
        T(qe),
        T(gain),
        T(dark_count_rate),
        dead_time,
        sensor,
        output_precision,
        layout,
    )
    channels = backend{T}(undef, 1, 1)
    noise_buffer = backend{T}(undef, 1, 1)
    output_buffer = output_precision === nothing ? nothing : backend{output_precision}(undef, 1, 1)
    fill!(channels, zero(T))
    fill!(noise_buffer, zero(T))
    output_buffer === nothing || fill!(output_buffer, zero(eltype(output_buffer)))
    state = APDDetectorState{T,typeof(channels),typeof(output_buffer)}(channels, noise_buffer, output_buffer)
    return APDDetector{typeof(validated),typeof(params),typeof(state),typeof(gain_map)}(
        validated, params, state, gain_map)
end

function Detector(; integration_time::Real=1.0, qe::Real=1.0,
    psf_sampling::Int=1, binning::Int=1, noise=NoisePhoton(),
    gain::Real=1.0, dark_current::Real=0.0, bits::Union{Nothing,Int}=nothing,
    full_well::Union{Nothing,Real}=nothing, sensor::SensorType=CCDSensor(),
    response_model::FrameResponseModel=NullFrameResponse(),
    readout_window::Union{Nothing,FrameWindow}=nothing,
    output_precision::Union{Nothing,DataType}=nothing, background_flux=nothing, background_map=nothing,
    T::Type{<:AbstractFloat}=Float64, backend=Array)
    normalized = normalize_noise(noise)
    return Detector(normalized; integration_time=integration_time, qe=qe,
        psf_sampling=psf_sampling, binning=binning, gain=gain,
        dark_current=dark_current, bits=bits, full_well=full_well,
        sensor=sensor, response_model=response_model, readout_window=readout_window, output_precision=output_precision,
        background_flux=background_flux, background_map=background_map,
        T=T, backend=backend)
end

function Detector(noise::NoiseModel; integration_time::Real=1.0, qe::Real=1.0,
    psf_sampling::Int=1, binning::Int=1, gain::Real=1.0, dark_current::Real=0.0,
    bits::Union{Nothing,Int}=nothing, full_well::Union{Nothing,Real}=nothing,
    sensor::SensorType=CCDSensor(), response_model::FrameResponseModel=NullFrameResponse(),
    readout_window::Union{Nothing,FrameWindow}=nothing,
    output_precision::Union{Nothing,DataType}=nothing,
    background_flux=nothing, background_map=nothing,
    T::Type{<:AbstractFloat}=Float64, backend=Array)
    converted = convert_noise(noise, T)
    validated = validate_noise(converted)
    return _build_detector(validated; integration_time=integration_time, qe=qe,
        psf_sampling=psf_sampling, binning=binning, gain=gain,
        dark_current=dark_current, bits=bits, full_well=full_well,
        sensor=sensor, response_model=response_model, readout_window=readout_window, output_precision=output_precision,
        background_flux=background_flux, background_map=background_map,
        T=T, backend=backend)
end

function APDDetector(; integration_time::Real=1.0, qe::Real=1.0, noise::NoiseModel=NoisePhoton(),
    gain::Real=1.0, dark_count_rate::Real=0.0, output_precision::Union{Nothing,DataType}=nothing,
    layout::Symbol=:channels, channel_gain_map=nothing, dead_time_model::CountingDeadTimeModel=NoDeadTime(),
    T::Type{<:AbstractFloat}=Float64, backend=Array)
    return _build_apd_detector(noise; integration_time=integration_time, qe=qe, gain=gain,
        dark_count_rate=dark_count_rate, dead_time_model=dead_time_model, sensor=APDSensor(), output_precision=output_precision,
        layout=layout, channel_gain_map=channel_gain_map, T=T, backend=backend)
end

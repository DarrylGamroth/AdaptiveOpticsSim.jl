readout_ready(det::Detector) = det.state.readout_ready
output_frame(det::Detector) = det.state.output_buffer === nothing ? det.state.frame : det.state.output_buffer

detector_noise_symbol(::NoiseNone) = :none
detector_noise_symbol(::NoisePhoton) = :photon
detector_noise_symbol(::NoiseReadout) = :readout
detector_noise_symbol(::NoisePhotonReadout) = :photon_readout

detector_sensor_symbol(sensor::SensorType) =
    throw(InvalidConfiguration("missing detector_sensor_symbol overload for $(typeof(sensor))"))

frame_response_symbol(::NullFrameResponse) = :none
frame_response_symbol(::SeparableGaussianPixelResponse) = :separable_gaussian

frame_sampling_symbol(::FrameSensorType) = :single_read
frame_sampling_reads(::FrameSensorType) = nothing
frame_sampling_reference_reads(::FrameSensorType) = nothing
frame_sampling_signal_reads(::FrameSensorType) = nothing
sampling_read_time(::FrameSensorType, ::Type{T}) where {T<:AbstractFloat} = nothing
sampling_wallclock_time(::FrameSensorType, integration_time, ::Type{T}) where {T<:AbstractFloat} = nothing

supports_detector_mtf(::AbstractFrameDetector) = false
supports_detector_mtf(det::Detector) = supports_detector_mtf(det.params.response_model)
supports_detector_mtf(::NullFrameResponse) = false
supports_detector_mtf(::SeparableGaussianPixelResponse) = true

supports_clock_induced_charge(::FrameSensorType) = false
supports_column_readout_noise(::FrameSensorType) = false
supports_avalanche_gain(::FrameSensorType) = false
supports_avalanche_gain(::AvalancheFrameSensorType) = true
supports_sensor_glow(::FrameSensorType) = false
supports_nondestructive_reads(::FrameSensorType) = false
supports_reference_read_subtraction(::FrameSensorType) = false

supports_counting_noise(::AbstractCountingDetector) = false
supports_dead_time(::AbstractCountingDetector) = false
supports_channel_gain_map(::AbstractCountingDetector) = false

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

frame_response_width(::NullFrameResponse, ::Type{T}) where {T<:AbstractFloat} = nothing
frame_response_width(model::SeparableGaussianPixelResponse, ::Type{T}) where {T<:AbstractFloat} = T(model.response_width_px)
effective_readout_sigma(::FrameSensorType, sigma) = sigma
effective_dark_current_time(::FrameSensorType, exposure_time) = exposure_time
effective_sensor_glow_time(::FrameSensorType, exposure_time) = exposure_time

convert_frame_response_model(::NullFrameResponse, ::Type{T}, backend) where {T<:AbstractFloat} = NullFrameResponse()

function convert_frame_response_model(model::SeparableGaussianPixelResponse, ::Type{T}, backend) where {T<:AbstractFloat}
    kernel = backend{T}(undef, length(model.kernel))
    copyto!(kernel, T.(Array(model.kernel)))
    return SeparableGaussianPixelResponse{T,typeof(kernel)}(T(model.response_width_px), kernel)
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

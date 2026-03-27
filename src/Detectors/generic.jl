readout_ready(det::Detector) = det.state.readout_ready
output_frame(det::Detector) = det.state.output_buffer === nothing ? det.state.frame : det.state.output_buffer
readout_products(det::Detector) = det.state.readout_products
detector_reference_frame(det::Detector) = detector_reference_frame(det.state.readout_products)
detector_signal_frame(det::Detector) = detector_signal_frame(det.state.readout_products)
detector_read_cube(det::Detector) = detector_read_cube(det.state.readout_products)

detector_reference_frame(::NoFrameReadoutProducts) = nothing
detector_signal_frame(::NoFrameReadoutProducts) = nothing
detector_read_cube(::NoFrameReadoutProducts) = nothing
detector_reference_frame(products::SampledFrameReadoutProducts) = products.reference_frame
detector_signal_frame(products::SampledFrameReadoutProducts) = products.signal_frame
detector_read_cube(products::SampledFrameReadoutProducts) = products.read_cube

detector_noise_symbol(::NoiseNone) = :none
detector_noise_symbol(::NoisePhoton) = :photon
detector_noise_symbol(::NoiseReadout) = :readout
detector_noise_symbol(::NoisePhotonReadout) = :photon_readout

detector_sensor_symbol(sensor::SensorType) =
    throw(InvalidConfiguration("missing detector_sensor_symbol overload for $(typeof(sensor))"))

response_family(::NullFrameResponse) = :none
response_family(::GaussianPixelResponse) = :gaussian
response_family(::SampledFrameResponse) = :sampled
response_family(::RectangularPixelAperture) = :rectangular_aperture
response_family(::SeparablePixelMTF) = :separable_mtf

frame_response_symbol(model::AbstractFrameResponse) = response_family(model)

response_application_domain(::AbstractFrameResponse) = :image

is_shift_invariant(::AbstractFrameResponse) = true
supports_frequency_domain_application(::AbstractFrameResponse) = false
supports_frequency_domain_application(::AbstractFrameMTF) = true
supports_separable_application(::AbstractFrameResponse) = false
supports_separable_application(::NullFrameResponse) = true
supports_separable_application(::GaussianPixelResponse) = true
supports_separable_application(::RectangularPixelAperture) = true
supports_separable_application(::SeparablePixelMTF) = true
supports_subpixel_geometry(::AbstractFrameResponse) = false
supports_subpixel_geometry(::AbstractFrameMTF) = true

response_support(::NullFrameResponse) = nothing, nothing
response_support(model::GaussianPixelResponse) = length(model.kernel), length(model.kernel)
response_support(model::SampledFrameResponse) = size(model.kernel)
response_support(model::RectangularPixelAperture) = length(model.kernel_y), length(model.kernel_x)
response_support(model::SeparablePixelMTF) = length(model.kernel_y), length(model.kernel_x)

response_width_px(::NullFrameResponse, ::Type{T}) where {T<:AbstractFloat} = nothing
response_width_px(model::GaussianPixelResponse, ::Type{T}) where {T<:AbstractFloat} = T(model.response_width_px)
response_width_px(::SampledFrameResponse, ::Type{T}) where {T<:AbstractFloat} = nothing
response_width_px(::AbstractFrameMTF, ::Type{T}) where {T<:AbstractFloat} = nothing

response_pitch_x_px(::AbstractFrameResponse, ::Type{T}) where {T<:AbstractFloat} = nothing
response_pitch_y_px(::AbstractFrameResponse, ::Type{T}) where {T<:AbstractFloat} = nothing
response_fill_factor_x(::AbstractFrameResponse, ::Type{T}) where {T<:AbstractFloat} = nothing
response_fill_factor_y(::AbstractFrameResponse, ::Type{T}) where {T<:AbstractFloat} = nothing
response_aperture_shape(::AbstractFrameResponse) = nothing
response_aperture_shape(::SampledFrameResponse) = :sampled

response_pitch_x_px(model::RectangularPixelAperture, ::Type{T}) where {T<:AbstractFloat} = T(model.pitch_x_px)
response_pitch_y_px(model::RectangularPixelAperture, ::Type{T}) where {T<:AbstractFloat} = T(model.pitch_y_px)
response_fill_factor_x(model::RectangularPixelAperture, ::Type{T}) where {T<:AbstractFloat} = T(model.fill_factor_x)
response_fill_factor_y(model::RectangularPixelAperture, ::Type{T}) where {T<:AbstractFloat} = T(model.fill_factor_y)
response_aperture_shape(::RectangularPixelAperture) = :rectangular

response_pitch_x_px(model::SeparablePixelMTF, ::Type{T}) where {T<:AbstractFloat} = T(model.pitch_x_px)
response_pitch_y_px(model::SeparablePixelMTF, ::Type{T}) where {T<:AbstractFloat} = T(model.pitch_y_px)
response_fill_factor_x(model::SeparablePixelMTF, ::Type{T}) where {T<:AbstractFloat} = T(model.fill_factor_x)
response_fill_factor_y(model::SeparablePixelMTF, ::Type{T}) where {T<:AbstractFloat} = T(model.fill_factor_y)
response_aperture_shape(::SeparablePixelMTF) = :rectangular

frame_sampling_symbol(::FrameSensorType) = :single_read
frame_sampling_reads(::FrameSensorType) = 1
frame_sampling_reference_reads(::FrameSensorType) = nothing
frame_sampling_signal_reads(::FrameSensorType) = nothing
sampling_read_time(::FrameSensorType, ::Type{T}) where {T<:AbstractFloat} = nothing
sampling_read_time(sensor::FrameSensorType, frame_size::Tuple{Int,Int}, window::Union{Nothing,FrameWindow}, ::Type{T}) where {T<:AbstractFloat} =
    sampling_read_time(sensor, T)
sampling_wallclock_time(::FrameSensorType, integration_time, ::Type{T}) where {T<:AbstractFloat} = T(integration_time)
sampling_wallclock_time(sensor::FrameSensorType, integration_time, frame_size::Tuple{Int,Int},
    window::Union{Nothing,FrameWindow}, ::Type{T}) where {T<:AbstractFloat} =
    sampling_wallclock_time(sensor, integration_time, T)

supports_detector_mtf(::AbstractFrameDetector) = false
supports_detector_mtf(det::Detector) = supports_detector_mtf(det.params.response_model)
supports_detector_mtf(::AbstractFrameResponse) = false
supports_detector_mtf(::GaussianPixelResponse) = true
supports_detector_mtf(::SampledFrameResponse) = true
supports_detector_mtf(::AbstractFrameMTF) = true

supports_clock_induced_charge(::FrameSensorType) = false
supports_column_readout_noise(::FrameSensorType) = false
supports_avalanche_gain(::FrameSensorType) = false
supports_avalanche_gain(::AvalancheFrameSensorType) = true
supports_sensor_glow(::FrameSensorType) = false
supports_nondestructive_reads(::FrameSensorType) = false
supports_reference_read_subtraction(::FrameSensorType) = false
supports_readout_correction(::FrameSensorType) = false
supports_read_cube(::FrameSensorType) = false

readout_correction_symbol(::NullFrameReadoutCorrection) = :none
readout_correction_symbol(::ReferencePixelCommonModeCorrection) = :reference_pixel_common_mode
correction_edge_rows(::FrameReadoutCorrectionModel) = nothing
correction_edge_cols(::FrameReadoutCorrectionModel) = nothing
correction_edge_rows(model::ReferencePixelCommonModeCorrection) = model.edge_rows
correction_edge_cols(model::ReferencePixelCommonModeCorrection) = model.edge_cols

supports_counting_noise(::AbstractCountingDetector) = false
supports_dead_time(::AbstractCountingDetector) = false
supports_channel_gain_map(::AbstractCountingDetector) = false

supports_detector_response(::SensorType, ::AbstractDetectorResponse) = false
supports_detector_response(::FrameSensorType, ::AbstractFrameResponse) = true

function default_response_model(::FrameSensorType; T::Type{<:AbstractFloat}=Float64, backend=Array)
    return NullFrameResponse()
end

function default_response_model(::CMOSSensor; T::Type{<:AbstractFloat}=Float64, backend=Array)
    return GaussianPixelResponse(response_width_px=0.35, T=T, backend=backend)
end

function default_response_model(::InGaAsSensor; T::Type{<:AbstractFloat}=Float64, backend=Array)
    return GaussianPixelResponse(response_width_px=0.4, T=T, backend=backend)
end

function default_response_model(::SAPHIRASensor; T::Type{<:AbstractFloat}=Float64, backend=Array)
    return SampledFrameResponse([0.0 0.01 0.0; 0.01 0.96 0.01; 0.0 0.01 0.0]; T=T, backend=backend)
end

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
    support_rows, support_cols = response_support(det.params.response_model)
    products = readout_products(det)
    read_cube = detector_read_cube(products)
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
        response_family(det.params.response_model),
        response_width_px(det.params.response_model, T),
        response_application_domain(det.params.response_model),
        supports_separable_application(det.params.response_model),
        is_shift_invariant(det.params.response_model),
        support_rows,
        support_cols,
        response_pitch_x_px(det.params.response_model, T),
        response_pitch_y_px(det.params.response_model, T),
        response_fill_factor_x(det.params.response_model, T),
        response_fill_factor_y(det.params.response_model, T),
        response_aperture_shape(det.params.response_model),
        row_window,
        col_window,
        frame_sampling_symbol(det.params.sensor),
        frame_sampling_reads(det.params.sensor),
        frame_sampling_reference_reads(det.params.sensor),
        frame_sampling_signal_reads(det.params.sensor),
        sampling_read_time(det.params.sensor, size(det.state.frame), det.params.readout_window, T),
        sampling_wallclock_time(det.params.sensor, det.params.integration_time, size(det.state.frame), det.params.readout_window, T),
        readout_correction_symbol(det.params.correction_model),
        correction_edge_rows(det.params.correction_model),
        correction_edge_cols(det.params.correction_model),
        !isnothing(detector_reference_frame(products)),
        !isnothing(detector_signal_frame(products)),
        !isnothing(read_cube),
        isnothing(read_cube) ? nothing : size(read_cube, 3),
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

effective_readout_sigma(::FrameSensorType, sigma) = sigma
effective_dark_current_time(::FrameSensorType, exposure_time) = exposure_time
effective_sensor_glow_time(::FrameSensorType, exposure_time) = exposure_time

convert_frame_response_model(::NullFrameResponse, ::Type{T}, backend) where {T<:AbstractFloat} = NullFrameResponse()

function convert_frame_response_model(model::GaussianPixelResponse, ::Type{T}, backend) where {T<:AbstractFloat}
    kernel = backend{T}(undef, length(model.kernel))
    copyto!(kernel, T.(Array(model.kernel)))
    return GaussianPixelResponse{T,typeof(kernel)}(T(model.response_width_px), kernel)
end

function convert_frame_response_model(model::SampledFrameResponse, ::Type{T}, backend) where {T<:AbstractFloat}
    kernel = backend{T}(undef, size(model.kernel)...)
    copyto!(kernel, T.(Array(model.kernel)))
    return SampledFrameResponse{T,typeof(kernel)}(kernel)
end

function convert_frame_response_model(model::RectangularPixelAperture, ::Type{T}, backend) where {T<:AbstractFloat}
    kernel_x = backend{T}(undef, length(model.kernel_x))
    kernel_y = backend{T}(undef, length(model.kernel_y))
    copyto!(kernel_x, T.(Array(model.kernel_x)))
    copyto!(kernel_y, T.(Array(model.kernel_y)))
    return RectangularPixelAperture{T,typeof(kernel_x),typeof(kernel_y)}(
        T(model.pitch_x_px), T(model.pitch_y_px), T(model.fill_factor_x), T(model.fill_factor_y), kernel_x, kernel_y)
end

function convert_frame_response_model(model::SeparablePixelMTF, ::Type{T}, backend) where {T<:AbstractFloat}
    kernel_x = backend{T}(undef, length(model.kernel_x))
    kernel_y = backend{T}(undef, length(model.kernel_y))
    copyto!(kernel_x, T.(Array(model.kernel_x)))
    copyto!(kernel_y, T.(Array(model.kernel_y)))
    return SeparablePixelMTF{T,typeof(kernel_x),typeof(kernel_y)}(
        T(model.pitch_x_px), T(model.pitch_y_px), T(model.fill_factor_x), T(model.fill_factor_y), kernel_x, kernel_y)
end

validate_frame_response_model(model::NullFrameResponse) = model

function validate_frame_response_model(model::GaussianPixelResponse)
    model.response_width_px > 0 || throw(InvalidConfiguration("GaussianPixelResponse response_width_px must be > 0"))
    length(model.kernel) > 0 || throw(InvalidConfiguration("GaussianPixelResponse kernel must not be empty"))
    isodd(length(model.kernel)) || throw(InvalidConfiguration("GaussianPixelResponse kernel length must be odd"))
    return model
end

function validate_frame_response_model(model::SampledFrameResponse)
    all(size(model.kernel) .> 0) || throw(InvalidConfiguration("SampledFrameResponse kernel must not be empty"))
    isodd(size(model.kernel, 1)) || throw(InvalidConfiguration("SampledFrameResponse kernel row count must be odd"))
    isodd(size(model.kernel, 2)) || throw(InvalidConfiguration("SampledFrameResponse kernel column count must be odd"))
    sum(model.kernel) > zero(eltype(model.kernel)) ||
        throw(InvalidConfiguration("SampledFrameResponse kernel must have positive sum"))
    return model
end

function validate_frame_response_model(model::RectangularPixelAperture)
    model.pitch_x_px > 0 || throw(InvalidConfiguration("RectangularPixelAperture pitch_x_px must be > 0"))
    model.pitch_y_px > 0 || throw(InvalidConfiguration("RectangularPixelAperture pitch_y_px must be > 0"))
    zero(model.fill_factor_x) < model.fill_factor_x <= one(model.fill_factor_x) ||
        throw(InvalidConfiguration("RectangularPixelAperture fill_factor_x must lie in (0, 1]"))
    zero(model.fill_factor_y) < model.fill_factor_y <= one(model.fill_factor_y) ||
        throw(InvalidConfiguration("RectangularPixelAperture fill_factor_y must lie in (0, 1]"))
    isodd(length(model.kernel_x)) || throw(InvalidConfiguration("RectangularPixelAperture kernel_x length must be odd"))
    isodd(length(model.kernel_y)) || throw(InvalidConfiguration("RectangularPixelAperture kernel_y length must be odd"))
    return model
end

function validate_frame_response_model(model::SeparablePixelMTF)
    model.pitch_x_px > 0 || throw(InvalidConfiguration("SeparablePixelMTF pitch_x_px must be > 0"))
    model.pitch_y_px > 0 || throw(InvalidConfiguration("SeparablePixelMTF pitch_y_px must be > 0"))
    zero(model.fill_factor_x) < model.fill_factor_x <= one(model.fill_factor_x) ||
        throw(InvalidConfiguration("SeparablePixelMTF fill_factor_x must lie in (0, 1]"))
    zero(model.fill_factor_y) < model.fill_factor_y <= one(model.fill_factor_y) ||
        throw(InvalidConfiguration("SeparablePixelMTF fill_factor_y must lie in (0, 1]"))
    isodd(length(model.kernel_x)) || throw(InvalidConfiguration("SeparablePixelMTF kernel_x length must be odd"))
    isodd(length(model.kernel_y)) || throw(InvalidConfiguration("SeparablePixelMTF kernel_y length must be odd"))
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

function validate_detector_response(sensor::SensorType, response_model::AbstractDetectorResponse)
    supports_detector_response(sensor, response_model) && return response_model
    throw(InvalidConfiguration("response model $(typeof(response_model)) is not supported for sensor $(typeof(sensor))"))
end

validate_readout_correction_model(::NullFrameReadoutCorrection) = NullFrameReadoutCorrection()

function validate_readout_correction_model(model::ReferencePixelCommonModeCorrection)
    return ReferencePixelCommonModeCorrection(model.edge_rows, model.edge_cols)
end

resolve_correction_model(sensor::SensorType, ::Nothing) = NullFrameReadoutCorrection()
resolve_correction_model(sensor::SensorType, correction_model::FrameReadoutCorrectionModel) = correction_model

function validate_readout_correction(sensor::SensorType, correction_model::FrameReadoutCorrectionModel)
    supports_readout_correction(sensor) && return correction_model
    correction_model isa NullFrameReadoutCorrection && return correction_model
    throw(InvalidConfiguration("readout correction $(typeof(correction_model)) is not supported for sensor $(typeof(sensor))"))
end

resolve_response_model(sensor::SensorType, ::Nothing; T::Type{<:AbstractFloat}, backend) =
    default_response_model(sensor; T=T, backend=backend)
resolve_response_model(sensor::SensorType, response_model::AbstractFrameResponse; T::Type{<:AbstractFloat}, backend) =
    response_model

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
    response_model::Union{Nothing,FrameResponseModel}, correction_model::Union{Nothing,FrameReadoutCorrectionModel},
    readout_window::Union{Nothing,FrameWindow},
    output_precision::Union{Nothing,DataType}, background_flux, background_map,
    T::Type{<:AbstractFloat}, backend)
    validate_frame_detector_sensor(sensor)
    full_well_t = full_well === nothing ? nothing : T(full_well)
    flux_model = background_model(background_flux; T=T, backend=backend)
    map_model = background_model(background_map; T=T, backend=backend)
    resolved_response = resolve_response_model(sensor, response_model; T=T, backend=backend)
    response = validate_detector_response(sensor,
        validate_frame_response_model(convert_frame_response_model(resolved_response, T, backend)))
    resolved_correction = resolve_correction_model(sensor, correction_model)
    correction = validate_readout_correction(sensor, validate_readout_correction_model(resolved_correction))
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
        correction,
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
    state = DetectorState{T, typeof(frame), typeof(output_buffer), FrameReadoutProducts}(
        frame,
        response_buffer,
        bin_buffer,
        noise_buffer,
        accum_buffer,
        output_buffer,
        NoFrameReadoutProducts(),
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
    response_model::Union{Nothing,FrameResponseModel}=nothing,
    correction_model::Union{Nothing,FrameReadoutCorrectionModel}=nothing,
    readout_window::Union{Nothing,FrameWindow}=nothing,
    output_precision::Union{Nothing,DataType}=nothing, background_flux=nothing, background_map=nothing,
    T::Type{<:AbstractFloat}=Float64, backend=Array)
    normalized = normalize_noise(noise)
    return Detector(normalized; integration_time=integration_time, qe=qe,
        psf_sampling=psf_sampling, binning=binning, gain=gain,
        dark_current=dark_current, bits=bits, full_well=full_well,
        sensor=sensor, response_model=response_model, correction_model=correction_model,
        readout_window=readout_window, output_precision=output_precision,
        background_flux=background_flux, background_map=background_map,
        T=T, backend=backend)
end

function Detector(noise::NoiseModel; integration_time::Real=1.0, qe::Real=1.0,
    psf_sampling::Int=1, binning::Int=1, gain::Real=1.0, dark_current::Real=0.0,
    bits::Union{Nothing,Int}=nothing, full_well::Union{Nothing,Real}=nothing,
    sensor::SensorType=CCDSensor(), response_model::Union{Nothing,FrameResponseModel}=nothing,
    correction_model::Union{Nothing,FrameReadoutCorrectionModel}=nothing,
    readout_window::Union{Nothing,FrameWindow}=nothing,
    output_precision::Union{Nothing,DataType}=nothing,
    background_flux=nothing, background_map=nothing,
    T::Type{<:AbstractFloat}=Float64, backend=Array)
    converted = convert_noise(noise, T)
    validated = validate_noise(converted)
    return _build_detector(validated; integration_time=integration_time, qe=qe,
        psf_sampling=psf_sampling, binning=binning, gain=gain,
        dark_current=dark_current, bits=bits, full_well=full_well,
        sensor=sensor, response_model=response_model, correction_model=correction_model,
        readout_window=readout_window, output_precision=output_precision,
        background_flux=background_flux, background_map=background_map,
        T=T, backend=backend)
end

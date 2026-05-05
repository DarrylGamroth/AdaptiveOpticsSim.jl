readout_ready(det::Detector) = det.state.readout_ready
output_frame(det::Detector) = det.state.output_buffer === nothing ? det.state.frame : det.state.output_buffer
detector_output_type(det::Detector) = det.params.output_type
readout_products(det::Detector) = det.state.readout_products
detector_reference_frame(det::Detector) = detector_reference_frame(det.state.readout_products)
detector_signal_frame(det::Detector) = detector_signal_frame(det.state.readout_products)
detector_combined_frame(det::Detector) = detector_combined_frame(det.state.readout_products)
detector_reference_cube(det::Detector) = detector_reference_cube(det.state.readout_products)
detector_signal_cube(det::Detector) = detector_signal_cube(det.state.readout_products)
detector_read_cube(det::Detector) = detector_read_cube(det.state.readout_products)
detector_read_times(det::Detector) = detector_read_times(det.state.readout_products)

@inline detector_output_value(::Nothing, value) = value
@inline detector_output_value(::Type{T}, value) where {T<:Integer} =
    T(clamp(round(value), typemin(T), typemax(T)))
@inline detector_output_value(::Type{T}, value) where {T} = T(value)

detector_reference_frame(::FrameReadoutProducts) = nothing
detector_signal_frame(::FrameReadoutProducts) = nothing
detector_combined_frame(::FrameReadoutProducts) = nothing
detector_reference_cube(::FrameReadoutProducts) = nothing
detector_signal_cube(::FrameReadoutProducts) = nothing
detector_read_cube(::FrameReadoutProducts) = nothing
detector_read_times(::FrameReadoutProducts) = nothing
detector_reference_frame(products::SampledFrameReadoutProducts) = products.reference_frame
detector_signal_frame(products::SampledFrameReadoutProducts) = products.signal_frame
detector_combined_frame(::SampledFrameReadoutProducts) = nothing
detector_reference_cube(::SampledFrameReadoutProducts) = nothing
detector_signal_cube(::SampledFrameReadoutProducts) = nothing
detector_read_cube(products::SampledFrameReadoutProducts) = products.read_cube
detector_read_times(::SampledFrameReadoutProducts) = nothing
detector_reference_frame(products::MultiReadFrameReadoutProducts) = products.reference_frame
detector_signal_frame(products::MultiReadFrameReadoutProducts) = products.signal_frame
detector_combined_frame(products::MultiReadFrameReadoutProducts) = products.combined_frame
detector_reference_cube(products::MultiReadFrameReadoutProducts) = products.reference_cube
detector_signal_cube(products::MultiReadFrameReadoutProducts) = products.signal_cube
detector_read_cube(products::MultiReadFrameReadoutProducts) = products.read_cube
detector_read_times(products::MultiReadFrameReadoutProducts) = products.read_times

thermal_model(det::Detector) = det.params.thermal_model
thermal_state(det::Detector) = det.state.thermal_state

detector_temperature(det::Detector, ::Type{T}=eltype(det.state.frame)) where {T<:AbstractFloat} =
    detector_temperature_K(det.params.thermal_model, det.state.thermal_state, T)

advance_thermal!(det::Detector, dt) = (advance_thermal!(det.params.thermal_model, det.state.thermal_state, dt); det)

supports_detector_mtf(::AbstractFrameDetector) = false
supports_detector_mtf(det::Detector) = supports_detector_mtf(det.params.response_model)

supports_detector_thermal_model(det::Detector) = !is_null_thermal_model(det.params.thermal_model)

supports_temperature_dependent_dark_current(det::Detector) =
    !is_null_temperature_law(active_dark_current_law(det.params.sensor, det.params.thermal_model))
supports_temperature_dependent_glow(det::Detector) =
    !is_null_temperature_law(active_glow_rate_law(det.params.sensor, det.params.thermal_model))
supports_temperature_dependent_persistence(::Detector) = false
supports_temperature_dependent_dark_counts(det::AbstractCountingDetector) =
    !is_null_temperature_law(active_dark_count_law(det, thermal_model(det)))

quantum_efficiency_model(det::Detector) = det.params.quantum_efficiency_model

ScalarQuantumEfficiency(qe::Real; T::Type{<:AbstractFloat}=Float64) =
    validate_quantum_efficiency_model(ScalarQuantumEfficiency{T}(T(qe)))

function SampledQuantumEfficiency(wavelengths::AbstractVector, values::AbstractVector;
    out_of_band::Real=0.0, T::Type{<:AbstractFloat}=Float64)
    length(wavelengths) == length(values) ||
        throw(DimensionMismatchError("SampledQuantumEfficiency wavelengths and values must have the same length"))
    return validate_quantum_efficiency_model(
        SampledQuantumEfficiency{T,Vector{T}}(T.(wavelengths), T.(values), T(out_of_band)))
end

resolve_quantum_efficiency_model(qe::Real, ::Type{T}) where {T<:AbstractFloat} =
    ScalarQuantumEfficiency{T}(T(qe))
resolve_quantum_efficiency_model(qe::AbstractQuantumEfficiencyModel, ::Type{T}) where {T<:AbstractFloat} =
    convert_quantum_efficiency_model(qe, T)

convert_quantum_efficiency_model(model::ScalarQuantumEfficiency, ::Type{T}) where {T<:AbstractFloat} =
    ScalarQuantumEfficiency{T}(T(model.value))
convert_quantum_efficiency_model(model::SampledQuantumEfficiency, ::Type{T}) where {T<:AbstractFloat} =
    SampledQuantumEfficiency{T,Vector{T}}(T.(model.wavelengths), T.(model.values), T(model.out_of_band))

function validate_quantum_efficiency_model(model::ScalarQuantumEfficiency)
    zero(model.value) <= model.value <= one(model.value) ||
        throw(InvalidConfiguration("ScalarQuantumEfficiency value must lie in [0, 1]"))
    return model
end

function validate_quantum_efficiency_model(model::SampledQuantumEfficiency)
    length(model.wavelengths) >= 2 ||
        throw(InvalidConfiguration("SampledQuantumEfficiency requires at least two samples"))
    length(model.wavelengths) == length(model.values) ||
        throw(DimensionMismatchError("SampledQuantumEfficiency wavelengths and values must have the same length"))
    zero(model.out_of_band) <= model.out_of_band <= one(model.out_of_band) ||
        throw(InvalidConfiguration("SampledQuantumEfficiency out_of_band must lie in [0, 1]"))
    @inbounds for i in eachindex(model.wavelengths)
        model.wavelengths[i] > zero(model.wavelengths[i]) ||
            throw(InvalidConfiguration("SampledQuantumEfficiency wavelengths must be > 0"))
        zero(model.values[i]) <= model.values[i] <= one(model.values[i]) ||
            throw(InvalidConfiguration("SampledQuantumEfficiency values must lie in [0, 1]"))
        if i > firstindex(model.wavelengths)
            model.wavelengths[i] > model.wavelengths[i - 1] ||
                throw(InvalidConfiguration("SampledQuantumEfficiency wavelengths must be strictly increasing"))
        end
    end
    return model
end

qe_at(model::ScalarQuantumEfficiency, wavelength) = model.value

function qe_at(model::SampledQuantumEfficiency{T}, wavelength) where {T<:AbstractFloat}
    λ = T(wavelength)
    λ < first(model.wavelengths) && return model.out_of_band
    λ > last(model.wavelengths) && return model.out_of_band
    @inbounds for hi in eachindex(model.wavelengths)
        λ == model.wavelengths[hi] && return model.values[hi]
        if λ < model.wavelengths[hi]
            lo = hi - 1
            λ0 = model.wavelengths[lo]
            λ1 = model.wavelengths[hi]
            v0 = model.values[lo]
            v1 = model.values[hi]
            return v0 + (v1 - v0) * (λ - λ0) / (λ1 - λ0)
        end
    end
    return last(model.values)
end

effective_qe(model::AbstractQuantumEfficiencyModel, src::AbstractSource, ::Type{T}=Float64) where {T<:AbstractFloat} =
    T(qe_at(model, wavelength(src)))

function effective_qe(model::AbstractQuantumEfficiencyModel, src::SpectralSource, ::Type{T}=Float64) where {T<:AbstractFloat}
    total = zero(T)
    @inbounds for sample in spectral_bundle(src)
        total += T(sample.weight) * T(qe_at(model, sample.wavelength))
    end
    return total
end

effective_qe(det::Detector, src::AbstractSource, ::Type{T}=eltype(det.state.frame)) where {T<:AbstractFloat} =
    effective_qe(det.params.quantum_efficiency_model, src, T)

reference_qe(model::ScalarQuantumEfficiency, ::Type{T}) where {T<:AbstractFloat} = T(model.value)
reference_qe(model::SampledQuantumEfficiency, ::Type{T}) where {T<:AbstractFloat} = T(maximum(model.values))

configured_glow_rate(::FrameSensorType, ::Type{T}) where {T<:AbstractFloat} = zero(T)
configured_cic_rate(::FrameSensorType, ::Type{T}) where {T<:AbstractFloat} = zero(T)

dark_current_law(::SensorType) = NullTemperatureLaw()
glow_rate_law(::SensorType) = NullTemperatureLaw()
cic_rate_law(::SensorType) = NullTemperatureLaw()

active_dark_current_law(sensor::SensorType, ::NullDetectorThermalModel) = dark_current_law(sensor)
active_glow_rate_law(sensor::SensorType, ::NullDetectorThermalModel) = glow_rate_law(sensor)
active_cic_rate_law(sensor::SensorType, ::NullDetectorThermalModel) = cic_rate_law(sensor)

function active_dark_current_law(sensor::SensorType, model::FixedTemperature)
    return is_null_temperature_law(model.dark_current_law) ? dark_current_law(sensor) : model.dark_current_law
end

function active_dark_current_law(sensor::SensorType, model::FirstOrderThermalModel)
    return is_null_temperature_law(model.dark_current_law) ? dark_current_law(sensor) : model.dark_current_law
end

function active_glow_rate_law(sensor::SensorType, model::FixedTemperature)
    return is_null_temperature_law(model.glow_rate_law) ? glow_rate_law(sensor) : model.glow_rate_law
end

function active_glow_rate_law(sensor::SensorType, model::FirstOrderThermalModel)
    return is_null_temperature_law(model.glow_rate_law) ? glow_rate_law(sensor) : model.glow_rate_law
end

function active_cic_rate_law(sensor::SensorType, model::FixedTemperature)
    return is_null_temperature_law(model.cic_rate_law) ? cic_rate_law(sensor) : model.cic_rate_law
end

function active_cic_rate_law(sensor::SensorType, model::FirstOrderThermalModel)
    return is_null_temperature_law(model.cic_rate_law) ? cic_rate_law(sensor) : model.cic_rate_law
end

effective_dark_current(det::Detector, ::Type{T}=eltype(det.state.frame)) where {T<:AbstractFloat} =
    T(evaluate_temperature_law(active_dark_current_law(det.params.sensor, det.params.thermal_model), T(det.params.dark_current), detector_temperature(det, T)))
effective_glow_rate(det::Detector, ::Type{T}=eltype(det.state.frame)) where {T<:AbstractFloat} =
    T(evaluate_temperature_law(active_glow_rate_law(det.params.sensor, det.params.thermal_model),
        configured_glow_rate(det.params.sensor, T), detector_temperature(det, T)))
effective_cic_rate(det::Detector, ::Type{T}=eltype(det.state.frame)) where {T<:AbstractFloat} =
    T(evaluate_temperature_law(active_cic_rate_law(det.params.sensor, det.params.thermal_model),
        configured_cic_rate(det.params.sensor, T), detector_temperature(det, T)))
effective_persistence_model(det::Detector) = persistence_model(det.params.sensor)

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
    reference_cube = detector_reference_cube(products)
    signal_cube = detector_signal_cube(products)
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
        det.params.output_type,
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
        detector_defect_symbol(det.params.defect_model),
        has_prnu(det.params.defect_model),
        has_dsnu(det.params.defect_model),
        has_bad_pixels(det.params.defect_model),
        row_window,
        col_window,
        timing_model_symbol(det.params.timing_model),
        line_time(det.params.timing_model, T),
        thermal_model_symbol(det.params.thermal_model),
        detector_temperature(det, T),
        ambient_temperature_K(det.params.thermal_model, T),
        cooling_setpoint_K(det.params.thermal_model, T),
        thermal_time_constant_s(det.params.thermal_model, T),
        temperature_law_symbol(active_dark_current_law(det.params.sensor, det.params.thermal_model)),
        temperature_law_symbol(active_glow_rate_law(det.params.sensor, det.params.thermal_model)),
        temperature_law_symbol(active_cic_rate_law(det.params.sensor, det.params.thermal_model)),
        frame_sampling_symbol(det.params.sensor),
        frame_sampling_reads(det.params.sensor),
        frame_sampling_reference_reads(det.params.sensor),
        frame_sampling_signal_reads(det.params.sensor),
        sampling_read_time(det.params.sensor, size(det.state.frame), det.params.readout_window, T),
        sampling_wallclock_time(det.params.sensor, det.params.integration_time, size(det.state.frame), det.params.readout_window, T),
        readout_correction_symbol(det.params.correction_model),
        correction_edge_rows(det.params.correction_model),
        correction_edge_cols(det.params.correction_model),
        correction_group_rows(det.params.correction_model),
        correction_group_cols(det.params.correction_model),
        correction_stage_count(det.params.correction_model),
        nonlinearity_symbol(det.params.nonlinearity_model),
        persistence_symbol(persistence_model(det.params.sensor)),
        !isnothing(detector_reference_frame(products)),
        !isnothing(detector_signal_frame(products)),
        !isnothing(detector_combined_frame(products)),
        !isnothing(reference_cube),
        !isnothing(signal_cube),
        !isnothing(read_cube),
        isnothing(reference_cube) ? nothing : size(reference_cube, 3),
        isnothing(signal_cube) ? nothing : size(signal_cube, 3),
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
    storage = _resolve_array_backend(backend)
    background = storage{T}(undef, size(map)...)
    copyto!(background, T.(map))
    return BackgroundFrame{T, typeof(background)}(background)
end

effective_readout_sigma(::FrameSensorType, sigma) = sigma
effective_dark_current_time(::FrameSensorType, exposure_time) = exposure_time
effective_sensor_glow_time(::FrameSensorType, exposure_time) = exposure_time
persistence_model(::FrameSensorType) = NullPersistence()

line_time(::AbstractFrameTimingModel, ::Type{T}) where {T<:AbstractFloat} = nothing
line_time(model::RollingShutter, ::Type{T}) where {T<:AbstractFloat} = T(model.line_time)

has_prnu(::AbstractDetectorDefectModel) = false
has_prnu(::PixelResponseNonuniformity) = true
has_prnu(model::CompositeDetectorDefectModel) = any(has_prnu, model.stages)
has_dsnu(::AbstractDetectorDefectModel) = false
has_dsnu(::DarkSignalNonuniformity) = true
has_dsnu(model::CompositeDetectorDefectModel) = any(has_dsnu, model.stages)
has_bad_pixels(::AbstractDetectorDefectModel) = false
has_bad_pixels(::BadPixelMask) = true
has_bad_pixels(model::CompositeDetectorDefectModel) = any(has_bad_pixels, model.stages)

counting_gate_duty_cycle(::AbstractCountingGateModel, ::Type{T}) where {T<:AbstractFloat} = nothing
counting_gate_duty_cycle(model::DutyCycleGate, ::Type{T}) where {T<:AbstractFloat} = T(model.duty_cycle)
afterpulse_probability(::AbstractCountingCorrelationModel, ::Type{T}) where {T<:AbstractFloat} = nothing
afterpulse_probability(model::AfterpulsingModel, ::Type{T}) where {T<:AbstractFloat} = T(model.probability)
afterpulse_probability(model::CompositeCountingCorrelation, ::Type{T}) where {T<:AbstractFloat} =
    _max_or_nothing((afterpulse_probability(stage, T) for stage in model.stages))
crosstalk_value(::AbstractCountingCorrelationModel, ::Type{T}) where {T<:AbstractFloat} = nothing
crosstalk_value(model::ChannelCrosstalkModel, ::Type{T}) where {T<:AbstractFloat} = T(model.coupling)
crosstalk_value(model::CompositeCountingCorrelation, ::Type{T}) where {T<:AbstractFloat} =
    _max_or_nothing((crosstalk_value(stage, T) for stage in model.stages))

function _max_or_nothing(values)
    found = false
    max_value = nothing
    for value in values
        isnothing(value) && continue
        if !found || value > max_value
            max_value = value
            found = true
        end
    end
    return found ? max_value : nothing
end

convert_frame_response_model(::NullFrameResponse, ::Type{T}, backend) where {T<:AbstractFloat} = NullFrameResponse()

function convert_frame_response_model(model::GaussianPixelResponse, ::Type{T}, backend) where {T<:AbstractFloat}
    array_backend = _resolve_array_backend(backend)
    kernel = array_backend{T}(undef, length(model.kernel))
    copyto!(kernel, T.(Array(model.kernel)))
    return GaussianPixelResponse{T,typeof(kernel)}(T(model.response_width_px), kernel)
end

function convert_frame_response_model(model::SampledFrameResponse, ::Type{T}, backend) where {T<:AbstractFloat}
    array_backend = _resolve_array_backend(backend)
    kernel = array_backend{T}(undef, size(model.kernel)...)
    copyto!(kernel, T.(Array(model.kernel)))
    return SampledFrameResponse{T,typeof(kernel)}(kernel)
end

function convert_frame_response_model(model::RectangularPixelAperture, ::Type{T}, backend) where {T<:AbstractFloat}
    array_backend = _resolve_array_backend(backend)
    kernel_x = array_backend{T}(undef, length(model.kernel_x))
    kernel_y = array_backend{T}(undef, length(model.kernel_y))
    copyto!(kernel_x, T.(Array(model.kernel_x)))
    copyto!(kernel_y, T.(Array(model.kernel_y)))
    return RectangularPixelAperture{T,typeof(kernel_x),typeof(kernel_y)}(
        T(model.pitch_x_px), T(model.pitch_y_px), T(model.fill_factor_x), T(model.fill_factor_y), kernel_x, kernel_y)
end

function convert_frame_response_model(model::SeparablePixelMTF, ::Type{T}, backend) where {T<:AbstractFloat}
    array_backend = _resolve_array_backend(backend)
    kernel_x = array_backend{T}(undef, length(model.kernel_x))
    kernel_y = array_backend{T}(undef, length(model.kernel_y))
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
    kernel_sum = _frame_response_kernel_sum(model.kernel)
    kernel_sum > zero(eltype(model.kernel)) ||
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

convert_detector_defect_model(::NullDetectorDefectModel, ::Type{T}, backend) where {T<:AbstractFloat} = NullDetectorDefectModel()

function convert_detector_defect_model(model::PixelResponseNonuniformity, ::Type{T}, backend) where {T<:AbstractFloat}
    gain_map = _to_backend_matrix(T.(Array(model.gain_map)), backend)
    return PixelResponseNonuniformity{T,typeof(gain_map)}(gain_map)
end

function convert_detector_defect_model(model::DarkSignalNonuniformity, ::Type{T}, backend) where {T<:AbstractFloat}
    dark_map = _to_backend_matrix(T.(Array(model.dark_map)), backend)
    return DarkSignalNonuniformity{T,typeof(dark_map)}(dark_map)
end

function convert_detector_defect_model(model::BadPixelMask, ::Type{T}, backend) where {T<:AbstractFloat}
    mask = _to_backend_bool_matrix(Array(model.mask), backend)
    return BadPixelMask{T,typeof(mask)}(mask, T(model.throughput))
end

function convert_detector_defect_model(model::CompositeDetectorDefectModel, ::Type{T}, backend) where {T<:AbstractFloat}
    return CompositeDetectorDefectModel(tuple((convert_detector_defect_model(stage, T, backend) for stage in model.stages)...))
end

validate_detector_defect_model(::NullDetectorDefectModel) = NullDetectorDefectModel()

function validate_detector_defect_model(model::PixelResponseNonuniformity)
    isempty(model.gain_map) && throw(InvalidConfiguration("PixelResponseNonuniformity gain_map must not be empty"))
    minimum(model.gain_map) >= zero(eltype(model.gain_map)) ||
        throw(InvalidConfiguration("PixelResponseNonuniformity gain_map must be >= 0"))
    return model
end

function validate_detector_defect_model(model::DarkSignalNonuniformity)
    isempty(model.dark_map) && throw(InvalidConfiguration("DarkSignalNonuniformity dark_map must not be empty"))
    minimum(model.dark_map) >= zero(eltype(model.dark_map)) ||
        throw(InvalidConfiguration("DarkSignalNonuniformity dark_map must be >= 0"))
    return model
end

function validate_detector_defect_model(model::BadPixelMask)
    isempty(model.mask) && throw(InvalidConfiguration("BadPixelMask mask must not be empty"))
    zero(model.throughput) <= model.throughput <= one(model.throughput) ||
        throw(InvalidConfiguration("BadPixelMask throughput must lie in [0, 1]"))
    return model
end

function validate_detector_defect_model(model::CompositeDetectorDefectModel)
    return CompositeDetectorDefectModel(tuple((validate_detector_defect_model(stage) for stage in model.stages)...))
end

convert_frame_timing_model(::GlobalShutter, ::Type{T}) where {T<:AbstractFloat} = GlobalShutter()
convert_frame_timing_model(model::RollingShutter, ::Type{T}) where {T<:AbstractFloat} =
    RollingShutter{T}(T(model.line_time), model.row_group_size; exposure_mode=model.exposure_mode)

validate_frame_timing_model(::GlobalShutter) = GlobalShutter()

function validate_frame_timing_model(model::RollingShutter)
    model.line_time >= zero(model.line_time) || throw(InvalidConfiguration("RollingShutter line_time must be >= 0"))
    model.row_group_size > 0 || throw(InvalidConfiguration("RollingShutter row_group_size must be > 0"))
    return model
end

convert_frame_nonlinearity_model(::NullFrameNonlinearity, ::Type{T}) where {T<:AbstractFloat} = NullFrameNonlinearity()
convert_frame_nonlinearity_model(model::SaturatingFrameNonlinearity, ::Type{T}) where {T<:AbstractFloat} =
    SaturatingFrameNonlinearity{T}(T(model.coefficient))

validate_frame_nonlinearity_model(::NullFrameNonlinearity) = NullFrameNonlinearity()

function validate_frame_nonlinearity_model(model::SaturatingFrameNonlinearity)
    model.coefficient >= zero(model.coefficient) ||
        throw(InvalidConfiguration("SaturatingFrameNonlinearity coefficient must be >= 0"))
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

resolve_detector_defect_model(sensor::FrameSensorType, ::Nothing; T::Type{<:AbstractFloat}, backend) =
    default_detector_defect_model(sensor; T=T, backend=backend)
resolve_detector_defect_model(sensor::FrameSensorType, model::AbstractDetectorDefectModel; T::Type{<:AbstractFloat}, backend) = model

validate_detector_defect(sensor::FrameSensorType, model::AbstractDetectorDefectModel) = model

resolve_frame_timing_model(sensor::FrameSensorType, ::Nothing; T::Type{<:AbstractFloat}) =
    default_frame_timing_model(sensor; T=T)
resolve_frame_timing_model(sensor::FrameSensorType, model::AbstractFrameTimingModel; T::Type{<:AbstractFloat}) = model

function validate_frame_timing(sensor::FrameSensorType, model::AbstractFrameTimingModel)
    supports_shutter_timing(sensor) && return model
    is_global_shutter(model) && return model
    throw(InvalidConfiguration("frame timing $(typeof(model)) is not supported for sensor $(typeof(sensor))"))
end

resolve_frame_nonlinearity_model(sensor::FrameSensorType, ::Nothing; T::Type{<:AbstractFloat}) =
    default_frame_nonlinearity_model(sensor; T=T)
resolve_frame_nonlinearity_model(sensor::FrameSensorType, model::AbstractFrameNonlinearityModel; T::Type{<:AbstractFloat}) = model

function validate_frame_nonlinearity(sensor::FrameSensorType, model::AbstractFrameNonlinearityModel)
    supports_detector_nonlinearity(sensor) && return model
    is_null_frame_nonlinearity(model) && return model
    throw(InvalidConfiguration("frame nonlinearity $(typeof(model)) is not supported for sensor $(typeof(sensor))"))
end

validate_readout_correction_model(::NullFrameReadoutCorrection) = NullFrameReadoutCorrection()

function validate_readout_correction_model(model::ReferencePixelCommonModeCorrection)
    return ReferencePixelCommonModeCorrection(model.edge_rows, model.edge_cols)
end

function validate_readout_correction_model(model::ReferenceRowCommonModeCorrection)
    return ReferenceRowCommonModeCorrection(model.edge_cols)
end

function validate_readout_correction_model(model::ReferenceColumnCommonModeCorrection)
    return ReferenceColumnCommonModeCorrection(model.edge_rows)
end

function validate_readout_correction_model(model::ReferenceOutputCommonModeCorrection)
    return ReferenceOutputCommonModeCorrection(model.output_cols; edge_rows=model.edge_rows, edge_cols=model.edge_cols)
end

function validate_readout_correction_model(model::CompositeFrameReadoutCorrection)
    stages = map(validate_readout_correction_model, model.stages)
    return CompositeFrameReadoutCorrection(tuple(stages...))
end

resolve_correction_model(sensor::SensorType, ::Nothing) = NullFrameReadoutCorrection()
resolve_correction_model(sensor::SensorType, correction_model::FrameReadoutCorrectionModel) = correction_model

function validate_readout_correction(sensor::SensorType, correction_model::FrameReadoutCorrectionModel)
    supports_readout_correction(sensor) && return correction_model
    is_null_readout_correction(correction_model) && return correction_model
    throw(InvalidConfiguration("readout correction $(typeof(correction_model)) is not supported for sensor $(typeof(sensor))"))
end

resolve_response_model(sensor::SensorType, ::Nothing; T::Type{<:AbstractFloat}, backend) =
    default_response_model(sensor; T=T, backend=backend)
resolve_response_model(sensor::SensorType, response_model::AbstractFrameResponse; T::Type{<:AbstractFloat}, backend) =
    response_model

function resolve_output_type(bits::Union{Nothing,Int}, output_type::Union{Nothing,DataType})
    output_type !== nothing && return output_type
    bits === 8 && return UInt8
    bits === 16 && return UInt16
    bits === 32 && return UInt32
    bits === 64 && return UInt64
    return nothing
end

@inline function detector_readout_products_type(sensor::FrameSensorType, frame::A, ::Type{T}) where {T<:AbstractFloat,A<:AbstractMatrix{T}}
    supports_multi_read_readout_products(sensor) || return NoFrameReadoutProducts
    cube_type = typeof(similar(frame, size(frame, 1), size(frame, 2), 1))
    return Union{
        NoFrameReadoutProducts,
        MultiReadFrameReadoutProducts{A,Nothing,Nothing},
        MultiReadFrameReadoutProducts{A,cube_type,Nothing},
        MultiReadFrameReadoutProducts{A,cube_type,Vector{T}},
    }
end

@inline initial_readout_products(::FrameSensorType, frame::AbstractMatrix, ::Type{T}) where {T<:AbstractFloat} = NoFrameReadoutProducts()

function _build_detector(noise::NoiseModel; integration_time::Real, qe::Union{Real,AbstractQuantumEfficiencyModel},
    psf_sampling::Int, binning::Int, gain::Real, dark_current::Real,
    bits::Union{Nothing,Int}, full_well::Union{Nothing,Real}, sensor::SensorType,
    response_model::Union{Nothing,FrameResponseModel}, defect_model::Union{Nothing,AbstractDetectorDefectModel},
    timing_model::Union{Nothing,AbstractFrameTimingModel}, correction_model::Union{Nothing,FrameReadoutCorrectionModel},
    nonlinearity_model::Union{Nothing,AbstractFrameNonlinearityModel},
    thermal_model::Union{Nothing,AbstractDetectorThermalModel},
    readout_window::Union{Nothing,FrameWindow},
    output_type::Union{Nothing,DataType}, background_flux, background_map,
    T::Type{<:AbstractFloat}, backend)
    validate_frame_detector_sensor(sensor)
    selector = _resolve_backend_selector(backend)
    array_backend = _resolve_array_backend(backend)
    full_well_t = full_well === nothing ? nothing : T(full_well)
    qe_model = validate_quantum_efficiency_model(resolve_quantum_efficiency_model(qe, T))
    qe_scalar = reference_qe(qe_model, T)
    flux_model = background_model(background_flux; T=T, backend=selector)
    map_model = background_model(background_map; T=T, backend=selector)
    resolved_response = resolve_response_model(sensor, response_model; T=T, backend=selector)
    response = validate_detector_response(sensor,
        validate_frame_response_model(convert_frame_response_model(resolved_response, T, backend)))
    resolved_defect = resolve_detector_defect_model(sensor, defect_model; T=T, backend=selector)
    defects = validate_detector_defect(sensor,
        validate_detector_defect_model(convert_detector_defect_model(resolved_defect, T, backend)))
    resolved_timing = resolve_frame_timing_model(sensor, timing_model; T=T)
    timing = validate_frame_timing(sensor, validate_frame_timing_model(convert_frame_timing_model(resolved_timing, T)))
    resolved_correction = resolve_correction_model(sensor, correction_model)
    correction = validate_readout_correction(sensor, validate_readout_correction_model(resolved_correction))
    resolved_nonlinearity = resolve_frame_nonlinearity_model(sensor, nonlinearity_model; T=T)
    nonlinearity = validate_frame_nonlinearity(sensor,
        validate_frame_nonlinearity_model(convert_frame_nonlinearity_model(resolved_nonlinearity, T)))
    resolved_thermal = resolve_thermal_model(sensor, thermal_model; T=T)
    thermal = validate_thermal_model(convert_thermal_model(resolved_thermal, T))
    output_type_t = resolve_output_type(bits, output_type)
    window = validate_readout_window(readout_window)
    params = DetectorParams{T, typeof(sensor), typeof(qe_model), typeof(response), typeof(defects), typeof(timing), typeof(correction), typeof(nonlinearity), typeof(thermal)}(
        T(integration_time),
        qe_scalar,
        psf_sampling,
        binning,
        T(gain),
        T(dark_current),
        bits,
        full_well_t,
        sensor,
        qe_model,
        response,
        defects,
        timing,
        correction,
        nonlinearity,
        thermal,
        window,
        output_type_t,
    )
    frame = array_backend{T}(undef, 1, 1)
    response_buffer = array_backend{T}(undef, 1, 1)
    bin_buffer = array_backend{T}(undef, 1, 1)
    noise_buffer = array_backend{T}(undef, 1, 1)
    noise_buffer_host = Matrix{T}(undef, 1, 1)
    accum_buffer = array_backend{T}(undef, 1, 1)
    output_buffer = if output_type_t === nothing
        window === nothing ? nothing : array_backend{T}(undef, 1, 1)
    else
        array_backend{output_type_t}(undef, 1, 1)
    end
    fill!(frame, zero(T))
    fill!(response_buffer, zero(T))
    fill!(bin_buffer, zero(T))
    fill!(noise_buffer, zero(T))
    fill!(accum_buffer, zero(T))
    latent_buffer = array_backend{T}(undef, 1, 1)
    fill!(latent_buffer, zero(T))
    output_buffer === nothing || fill!(output_buffer, zero(eltype(output_buffer)))
    thermal_state = thermal_state_from_model(thermal, T)
    readout_products_type = detector_readout_products_type(sensor, frame, T)
    readout_products = initial_readout_products(sensor, frame, T)
    state = DetectorState{T, typeof(frame), typeof(output_buffer), readout_products_type, typeof(thermal_state)}(
        frame,
        response_buffer,
        bin_buffer,
        noise_buffer,
        noise_buffer_host,
        accum_buffer,
        latent_buffer,
        output_buffer,
        readout_products,
        thermal_state,
        zero(T),
        true,
    )
    return Detector{typeof(noise), typeof(params), typeof(state), typeof(flux_model), typeof(map_model), typeof(selector)}(
        noise,
        params,
        state,
        flux_model,
        map_model,
    )
end

function Detector(; integration_time::Real=1.0, qe::Union{Real,AbstractQuantumEfficiencyModel}=1.0,
    psf_sampling::Int=1, binning::Int=1, noise=NoisePhoton(),
    gain::Real=1.0, dark_current::Real=0.0, bits::Union{Nothing,Int}=nothing,
    full_well::Union{Nothing,Real}=nothing, sensor::SensorType=CCDSensor(),
    response_model::Union{Nothing,FrameResponseModel}=nothing,
    defect_model::Union{Nothing,AbstractDetectorDefectModel}=nothing,
    timing_model::Union{Nothing,AbstractFrameTimingModel}=nothing,
    correction_model::Union{Nothing,FrameReadoutCorrectionModel}=nothing,
    nonlinearity_model::Union{Nothing,AbstractFrameNonlinearityModel}=nothing,
    thermal_model::Union{Nothing,AbstractDetectorThermalModel}=nothing,
    readout_window::Union{Nothing,FrameWindow}=nothing,
    output_type::Union{Nothing,DataType}=nothing, background_flux=nothing, background_map=nothing,
    T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=CPUBackend())
    normalized = normalize_noise(noise)
    return Detector(normalized; integration_time=integration_time, qe=qe,
        psf_sampling=psf_sampling, binning=binning, gain=gain,
        dark_current=dark_current, bits=bits, full_well=full_well,
        sensor=sensor, response_model=response_model, defect_model=defect_model, timing_model=timing_model,
        correction_model=correction_model, nonlinearity_model=nonlinearity_model, thermal_model=thermal_model,
        readout_window=readout_window, output_type=output_type,
        background_flux=background_flux, background_map=background_map,
        T=T, backend=backend)
end

function Detector(noise::NoiseModel; integration_time::Real=1.0, qe::Union{Real,AbstractQuantumEfficiencyModel}=1.0,
    psf_sampling::Int=1, binning::Int=1, gain::Real=1.0, dark_current::Real=0.0,
    bits::Union{Nothing,Int}=nothing, full_well::Union{Nothing,Real}=nothing,
    sensor::SensorType=CCDSensor(), response_model::Union{Nothing,FrameResponseModel}=nothing,
    defect_model::Union{Nothing,AbstractDetectorDefectModel}=nothing,
    timing_model::Union{Nothing,AbstractFrameTimingModel}=nothing,
    correction_model::Union{Nothing,FrameReadoutCorrectionModel}=nothing,
    nonlinearity_model::Union{Nothing,AbstractFrameNonlinearityModel}=nothing,
    thermal_model::Union{Nothing,AbstractDetectorThermalModel}=nothing,
    readout_window::Union{Nothing,FrameWindow}=nothing,
    output_type::Union{Nothing,DataType}=nothing,
    background_flux=nothing, background_map=nothing,
    T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=CPUBackend())
    converted = convert_noise(noise, T)
    validated = validate_noise(converted)
    return _build_detector(validated; integration_time=integration_time, qe=qe,
        psf_sampling=psf_sampling, binning=binning, gain=gain,
        dark_current=dark_current, bits=bits, full_well=full_well,
        sensor=sensor, response_model=response_model, defect_model=defect_model,
        timing_model=timing_model, correction_model=correction_model, nonlinearity_model=nonlinearity_model,
        thermal_model=thermal_model,
        readout_window=readout_window, output_type=output_type,
        background_flux=background_flux, background_map=background_map,
        T=T, backend=backend)
end

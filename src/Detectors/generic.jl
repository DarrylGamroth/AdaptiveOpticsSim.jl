readout_ready(det::Detector) = det.state.readout_ready
output_frame(det::Detector) = det.state.output_buffer === nothing ? det.state.frame : det.state.output_buffer
readout_products(det::Detector) = det.state.readout_products
detector_reference_frame(det::Detector) = detector_reference_frame(det.state.readout_products)
detector_signal_frame(det::Detector) = detector_signal_frame(det.state.readout_products)
detector_combined_frame(det::Detector) = detector_combined_frame(det.state.readout_products)
detector_reference_cube(det::Detector) = detector_reference_cube(det.state.readout_products)
detector_signal_cube(det::Detector) = detector_signal_cube(det.state.readout_products)
detector_read_cube(det::Detector) = detector_read_cube(det.state.readout_products)
detector_read_times(det::Detector) = detector_read_times(det.state.readout_products)

detector_reference_frame(::NoFrameReadoutProducts) = nothing
detector_signal_frame(::NoFrameReadoutProducts) = nothing
detector_combined_frame(::NoFrameReadoutProducts) = nothing
detector_reference_cube(::NoFrameReadoutProducts) = nothing
detector_signal_cube(::NoFrameReadoutProducts) = nothing
detector_read_cube(::NoFrameReadoutProducts) = nothing
detector_read_times(::NoFrameReadoutProducts) = nothing
detector_reference_frame(products::SampledFrameReadoutProducts) = products.reference_frame
detector_signal_frame(products::SampledFrameReadoutProducts) = products.signal_frame
detector_combined_frame(::SampledFrameReadoutProducts) = nothing
detector_reference_cube(::SampledFrameReadoutProducts) = nothing
detector_signal_cube(::SampledFrameReadoutProducts) = nothing
detector_read_cube(products::SampledFrameReadoutProducts) = products.read_cube
detector_read_times(::SampledFrameReadoutProducts) = nothing
detector_reference_frame(products::HgCdTeReadoutProducts) = products.reference_frame
detector_signal_frame(products::HgCdTeReadoutProducts) = products.signal_frame
detector_combined_frame(products::HgCdTeReadoutProducts) = products.combined_frame
detector_reference_cube(products::HgCdTeReadoutProducts) = products.reference_cube
detector_signal_cube(products::HgCdTeReadoutProducts) = products.signal_cube
detector_read_cube(products::HgCdTeReadoutProducts) = products.read_cube
detector_read_times(products::HgCdTeReadoutProducts) = products.read_times

detector_noise_symbol(::NoiseNone) = :none
detector_noise_symbol(::NoisePhoton) = :photon
detector_noise_symbol(::NoiseReadout) = :readout
detector_noise_symbol(::NoisePhotonReadout) = :photon_readout

detector_sensor_symbol(sensor::SensorType) =
    throw(InvalidConfiguration("missing detector_sensor_symbol overload for $(typeof(sensor))"))

detector_defect_symbol(::NullDetectorDefectModel) = :none
detector_defect_symbol(::PixelResponseNonuniformity) = :prnu
detector_defect_symbol(::DarkSignalNonuniformity) = :dsnu
detector_defect_symbol(::BadPixelMask) = :bad_pixel_mask
detector_defect_symbol(::CompositeDetectorDefectModel) = :composite

timing_model_symbol(::GlobalShutter) = :global_shutter
timing_model_symbol(::RollingShutter) = :rolling_shutter
is_global_shutter(::AbstractFrameTimingModel) = false
is_global_shutter(::GlobalShutter) = true

nonlinearity_symbol(::NullFrameNonlinearity) = :none
nonlinearity_symbol(::SaturatingFrameNonlinearity) = :saturating
is_null_frame_nonlinearity(::AbstractFrameNonlinearityModel) = false
is_null_frame_nonlinearity(::NullFrameNonlinearity) = true

persistence_symbol(::NullPersistence) = :none
persistence_symbol(::ExponentialPersistence) = :exponential
is_null_persistence(::AbstractPersistenceModel) = false
is_null_persistence(::NullPersistence) = true

thermal_model_symbol(::NullDetectorThermalModel) = :none
thermal_model_symbol(::FixedTemperature) = :fixed_temperature
thermal_model_symbol(::FirstOrderThermalModel) = :first_order
is_null_thermal_model(::AbstractDetectorThermalModel) = false
is_null_thermal_model(::NullDetectorThermalModel) = true

temperature_law_symbol(::NullTemperatureLaw) = :none
temperature_law_symbol(::ArrheniusRateLaw) = :arrhenius
temperature_law_symbol(::LinearTemperatureLaw) = :linear
temperature_law_symbol(::ExponentialTemperatureLaw) = :exponential
is_null_temperature_law(::AbstractTemperatureLaw) = false
is_null_temperature_law(::NullTemperatureLaw) = true

thermal_model(det::Detector) = det.params.thermal_model
thermal_state(det::Detector) = det.state.thermal_state

ambient_temperature_K(::AbstractDetectorThermalModel, ::Type{T}) where {T<:AbstractFloat} = nothing
cooling_setpoint_K(::AbstractDetectorThermalModel, ::Type{T}) where {T<:AbstractFloat} = nothing
thermal_time_constant_s(::AbstractDetectorThermalModel, ::Type{T}) where {T<:AbstractFloat} = nothing
cooling_setpoint_K(model::FixedTemperature, ::Type{T}) where {T<:AbstractFloat} = T(model.temperature_K)
ambient_temperature_K(model::FirstOrderThermalModel, ::Type{T}) where {T<:AbstractFloat} = T(model.ambient_temperature_K)
cooling_setpoint_K(model::FirstOrderThermalModel, ::Type{T}) where {T<:AbstractFloat} = T(model.setpoint_temperature_K)
thermal_time_constant_s(model::FirstOrderThermalModel, ::Type{T}) where {T<:AbstractFloat} = T(model.time_constant_s)
detector_temperature_K(::NullDetectorThermalModel, ::AbstractDetectorThermalState, ::Type{T}) where {T<:AbstractFloat} = nothing
detector_temperature_K(model::FixedTemperature, ::AbstractDetectorThermalState, ::Type{T}) where {T<:AbstractFloat} = T(model.temperature_K)
detector_temperature_K(model::AbstractDetectorThermalModel, state::DetectorThermalState, ::Type{T}) where {T<:AbstractFloat} = T(state.temperature_K)
detector_temperature(det::Detector, ::Type{T}=eltype(det.state.frame)) where {T<:AbstractFloat} =
    detector_temperature_K(det.params.thermal_model, det.state.thermal_state, T)

advance_thermal!(::NullDetectorThermalModel, ::AbstractDetectorThermalState, dt) = nothing
advance_thermal!(::FixedTemperature, ::AbstractDetectorThermalState, dt) = nothing
function advance_thermal!(model::FirstOrderThermalModel, state::DetectorThermalState, dt)
    dt >= 0 || throw(InvalidConfiguration("thermal evolution step dt must be >= 0"))
    iszero(dt) && return nothing
    equilibrium = clamp(model.setpoint_temperature_K, model.min_temperature_K, model.max_temperature_K)
    state.temperature_K = clamp(equilibrium + (state.temperature_K - equilibrium) * exp(-dt / model.time_constant_s),
        model.min_temperature_K, model.max_temperature_K)
    return nothing
end
advance_thermal!(det::Detector, dt) = (advance_thermal!(det.params.thermal_model, det.state.thermal_state, dt); det)

evaluate_temperature_law(::NullTemperatureLaw, base_value, detector_temperature) = base_value

function evaluate_temperature_law(law::ArrheniusRateLaw{T}, base_value, detector_temperature) where {T<:AbstractFloat}
    isnothing(detector_temperature) && return base_value
    detector_temperature > zero(detector_temperature) ||
        throw(InvalidConfiguration("detector temperature must be > 0 K"))
    scale = exp(law.activation_temperature_K * (inv(law.reference_temperature_K) - inv(detector_temperature)))
    return base_value * scale
end

function evaluate_temperature_law(law::LinearTemperatureLaw, base_value, detector_temperature)
    isnothing(detector_temperature) && return base_value
    scaled = base_value * (one(base_value) + law.slope_per_K * (detector_temperature - law.reference_temperature_K))
    return max(zero(base_value), scaled)
end

function evaluate_temperature_law(law::ExponentialTemperatureLaw, base_value, detector_temperature)
    isnothing(detector_temperature) && return base_value
    return base_value * exp(law.exponent_per_K * (detector_temperature - law.reference_temperature_K))
end

counting_gate_symbol(::NullCountingGate) = :none
counting_gate_symbol(::DutyCycleGate) = :duty_cycle

counting_correlation_symbol(::NullCountingCorrelation) = :none
counting_correlation_symbol(::AfterpulsingModel) = :afterpulsing
counting_correlation_symbol(::ChannelCrosstalkModel) = :channel_crosstalk
counting_correlation_symbol(::CompositeCountingCorrelation) = :composite

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
supports_batched_response_application(::ScalarCPUStyle, ::AbstractFrameResponse) = true
supports_batched_response_application(::AcceleratorStyle, ::AbstractFrameResponse) = false
supports_batched_response_application(::ScalarCPUStyle, ::NullFrameResponse) = true
supports_batched_response_application(::AcceleratorStyle, ::NullFrameResponse) = true
supports_batched_response_application(::ScalarCPUStyle, ::GaussianPixelResponse) = true
supports_batched_response_application(::AcceleratorStyle, ::GaussianPixelResponse) = true
supports_batched_response_application(::ScalarCPUStyle, ::SampledFrameResponse) = true
supports_batched_response_application(::AcceleratorStyle, ::SampledFrameResponse) = true
supports_batched_response_application(::ScalarCPUStyle, ::RectangularPixelAperture) = true
supports_batched_response_application(::AcceleratorStyle, ::RectangularPixelAperture) = true
supports_batched_response_application(::ScalarCPUStyle, ::SeparablePixelMTF) = true
supports_batched_response_application(::AcceleratorStyle, ::SeparablePixelMTF) = true
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
supports_batched_readout_correction(::FrameReadoutCorrectionModel) = false
supports_batched_readout_correction(::NullFrameReadoutCorrection) = true
supports_batched_readout_correction(::ReferencePixelCommonModeCorrection) = true
supports_batched_readout_correction(::ReferenceRowCommonModeCorrection) = true
supports_batched_readout_correction(::ReferenceColumnCommonModeCorrection) = true
supports_batched_readout_correction(::ReferenceOutputCommonModeCorrection) = true
supports_batched_readout_correction(model::CompositeFrameReadoutCorrection) =
    all(supports_batched_readout_correction, model.stages)

supports_clock_induced_charge(::FrameSensorType) = false
supports_column_readout_noise(::FrameSensorType) = false
supports_avalanche_gain(::FrameSensorType) = false
supports_avalanche_gain(::AvalancheFrameSensorType) = true
supports_sensor_glow(::FrameSensorType) = false
supports_detector_defect_maps(::FrameSensorType) = false
supports_detector_persistence(::FrameSensorType) = false
supports_detector_nonlinearity(::FrameSensorType) = false
supports_shutter_timing(::FrameSensorType) = false
supports_nondestructive_reads(::FrameSensorType) = false
supports_reference_read_subtraction(::FrameSensorType) = false
supports_readout_correction(::FrameSensorType) = false
supports_read_cube(::FrameSensorType) = false
supports_detector_thermal_model(::AbstractFrameDetector) = false
supports_detector_thermal_model(::AbstractCountingDetector) = false
supports_dynamic_thermal_state(::AbstractDetectorThermalModel) = false
supports_dynamic_thermal_state(::FirstOrderThermalModel) = true
supports_detector_thermal_model(det::Detector) = !is_null_thermal_model(det.params.thermal_model)

supports_temperature_dependent_dark_current(det::Detector) =
    !is_null_temperature_law(active_dark_current_law(det.params.sensor, det.params.thermal_model))
supports_temperature_dependent_glow(det::Detector) =
    !is_null_temperature_law(active_glow_rate_law(det.params.sensor, det.params.thermal_model))
supports_temperature_dependent_persistence(::Detector) = false
supports_temperature_dependent_dark_counts(::AbstractCountingDetector) = false

configured_glow_rate(::FrameSensorType, ::Type{T}) where {T<:AbstractFloat} = zero(T)
configured_cic_rate(::FrameSensorType, ::Type{T}) where {T<:AbstractFloat} = zero(T)

readout_correction_symbol(::NullFrameReadoutCorrection) = :none
readout_correction_symbol(::ReferencePixelCommonModeCorrection) = :reference_pixel_common_mode
readout_correction_symbol(::ReferenceRowCommonModeCorrection) = :reference_row_common_mode
readout_correction_symbol(::ReferenceColumnCommonModeCorrection) = :reference_column_common_mode
readout_correction_symbol(::ReferenceOutputCommonModeCorrection) = :reference_output_common_mode
readout_correction_symbol(::CompositeFrameReadoutCorrection) = :composite
is_null_readout_correction(::FrameReadoutCorrectionModel) = false
is_null_readout_correction(::NullFrameReadoutCorrection) = true
correction_edge_rows(::FrameReadoutCorrectionModel) = nothing
correction_edge_cols(::FrameReadoutCorrectionModel) = nothing
correction_edge_rows(model::ReferencePixelCommonModeCorrection) = model.edge_rows
correction_edge_cols(model::ReferencePixelCommonModeCorrection) = model.edge_cols
correction_edge_rows(model::ReferenceColumnCommonModeCorrection) = model.edge_rows
correction_edge_cols(model::ReferenceRowCommonModeCorrection) = model.edge_cols
correction_edge_rows(model::ReferenceOutputCommonModeCorrection) = model.edge_rows
correction_edge_cols(model::ReferenceOutputCommonModeCorrection) = model.edge_cols
correction_group_rows(::FrameReadoutCorrectionModel) = nothing
correction_group_cols(::FrameReadoutCorrectionModel) = nothing
correction_group_cols(model::ReferenceOutputCommonModeCorrection) = model.output_cols
correction_stage_count(::FrameReadoutCorrectionModel) = 1
correction_stage_count(model::CompositeFrameReadoutCorrection) = length(model.stages)

supports_counting_noise(::AbstractCountingDetector) = false
supports_dead_time(::AbstractCountingDetector) = false
supports_channel_gain_map(::AbstractCountingDetector) = false
supports_counting_gating(::AbstractCountingDetector) = false
supports_afterpulsing(::AbstractCountingDetector) = false
supports_channel_crosstalk(::AbstractCountingDetector) = false
supports_paralyzable_dead_time(::AbstractCountingDetector) = false

supports_detector_response(::SensorType, ::AbstractDetectorResponse) = false
supports_detector_response(::FrameSensorType, ::AbstractFrameResponse) = true

function default_response_model(::FrameSensorType; T::Type{<:AbstractFloat}=Float64, backend=CPUBackend())
    return NullFrameResponse()
end

function default_response_model(::CMOSSensor; T::Type{<:AbstractFloat}=Float64, backend=CPUBackend())
    return GaussianPixelResponse(response_width_px=0.35, T=T, backend=backend)
end

function default_response_model(::InGaAsSensor; T::Type{<:AbstractFloat}=Float64, backend=CPUBackend())
    return GaussianPixelResponse(response_width_px=0.4, T=T, backend=backend)
end

function default_response_model(::HgCdTeAvalancheArraySensor; T::Type{<:AbstractFloat}=Float64, backend=CPUBackend())
    return SampledFrameResponse([0.0 0.01 0.0; 0.01 0.96 0.01; 0.0 0.01 0.0]; T=T, backend=backend)
end

default_detector_defect_model(::FrameSensorType; T::Type{<:AbstractFloat}=Float64, backend=CPUBackend()) = NullDetectorDefectModel()
default_frame_timing_model(::FrameSensorType; T::Type{<:AbstractFloat}=Float64) = GlobalShutter()
default_frame_nonlinearity_model(::FrameSensorType; T::Type{<:AbstractFloat}=Float64) = NullFrameNonlinearity()
default_thermal_model(::SensorType; T::Type{<:AbstractFloat}=Float64) = NullDetectorThermalModel()

convert_temperature_law(::NullTemperatureLaw, ::Type{T}) where {T<:AbstractFloat} = NullTemperatureLaw()
convert_temperature_law(model::ArrheniusRateLaw, ::Type{T}) where {T<:AbstractFloat} =
    ArrheniusRateLaw{T}(T(model.reference_temperature_K), T(model.activation_temperature_K))
convert_temperature_law(model::LinearTemperatureLaw, ::Type{T}) where {T<:AbstractFloat} =
    LinearTemperatureLaw{T}(T(model.reference_temperature_K), T(model.slope_per_K))
convert_temperature_law(model::ExponentialTemperatureLaw, ::Type{T}) where {T<:AbstractFloat} =
    ExponentialTemperatureLaw{T}(T(model.reference_temperature_K), T(model.exponent_per_K))

validate_temperature_law(::NullTemperatureLaw) = NullTemperatureLaw()

function validate_temperature_law(model::ArrheniusRateLaw)
    model.reference_temperature_K > zero(model.reference_temperature_K) ||
        throw(InvalidConfiguration("ArrheniusRateLaw reference_temperature_K must be > 0"))
    return model
end

function validate_temperature_law(model::LinearTemperatureLaw)
    model.reference_temperature_K > zero(model.reference_temperature_K) ||
        throw(InvalidConfiguration("LinearTemperatureLaw reference_temperature_K must be > 0"))
    return model
end

function validate_temperature_law(model::ExponentialTemperatureLaw)
    model.reference_temperature_K > zero(model.reference_temperature_K) ||
        throw(InvalidConfiguration("ExponentialTemperatureLaw reference_temperature_K must be > 0"))
    return model
end

convert_thermal_model(::NullDetectorThermalModel, ::Type{T}) where {T<:AbstractFloat} = NullDetectorThermalModel()

function convert_thermal_model(model::FixedTemperature, ::Type{T}) where {T<:AbstractFloat}
    return FixedTemperature(
        temperature_K=T(model.temperature_K),
        dark_current_law=validate_temperature_law(convert_temperature_law(model.dark_current_law, T)),
        glow_rate_law=validate_temperature_law(convert_temperature_law(model.glow_rate_law, T)),
        dark_count_law=validate_temperature_law(convert_temperature_law(model.dark_count_law, T)),
        cic_rate_law=validate_temperature_law(convert_temperature_law(model.cic_rate_law, T)),
        T=T,
    )
end

function convert_thermal_model(model::FirstOrderThermalModel, ::Type{T}) where {T<:AbstractFloat}
    return FirstOrderThermalModel(
        ambient_temperature_K=T(model.ambient_temperature_K),
        setpoint_temperature_K=T(model.setpoint_temperature_K),
        initial_temperature_K=T(model.initial_temperature_K),
        time_constant_s=T(model.time_constant_s),
        min_temperature_K=T(model.min_temperature_K),
        max_temperature_K=T(model.max_temperature_K),
        dark_current_law=validate_temperature_law(convert_temperature_law(model.dark_current_law, T)),
        glow_rate_law=validate_temperature_law(convert_temperature_law(model.glow_rate_law, T)),
        dark_count_law=validate_temperature_law(convert_temperature_law(model.dark_count_law, T)),
        cic_rate_law=validate_temperature_law(convert_temperature_law(model.cic_rate_law, T)),
        T=T,
    )
end

validate_thermal_model(::NullDetectorThermalModel) = NullDetectorThermalModel()

function validate_thermal_model(model::FixedTemperature)
    model.temperature_K > zero(model.temperature_K) ||
        throw(InvalidConfiguration("FixedTemperature temperature_K must be > 0"))
    validate_temperature_law(model.dark_current_law)
    validate_temperature_law(model.glow_rate_law)
    validate_temperature_law(model.dark_count_law)
    validate_temperature_law(model.cic_rate_law)
    return model
end

function validate_thermal_model(model::FirstOrderThermalModel)
    model.ambient_temperature_K > zero(model.ambient_temperature_K) ||
        throw(InvalidConfiguration("FirstOrderThermalModel ambient_temperature_K must be > 0"))
    model.setpoint_temperature_K > zero(model.setpoint_temperature_K) ||
        throw(InvalidConfiguration("FirstOrderThermalModel setpoint_temperature_K must be > 0"))
    model.initial_temperature_K > zero(model.initial_temperature_K) ||
        throw(InvalidConfiguration("FirstOrderThermalModel initial_temperature_K must be > 0"))
    model.time_constant_s > zero(model.time_constant_s) ||
        throw(InvalidConfiguration("FirstOrderThermalModel time_constant_s must be > 0"))
    model.min_temperature_K > zero(model.min_temperature_K) ||
        throw(InvalidConfiguration("FirstOrderThermalModel min_temperature_K must be > 0"))
    model.max_temperature_K > zero(model.max_temperature_K) ||
        throw(InvalidConfiguration("FirstOrderThermalModel max_temperature_K must be > 0"))
    model.min_temperature_K <= model.initial_temperature_K <= model.max_temperature_K ||
        throw(InvalidConfiguration("FirstOrderThermalModel initial_temperature_K must lie within [min_temperature_K, max_temperature_K]"))
    model.min_temperature_K <= model.setpoint_temperature_K <= model.max_temperature_K ||
        throw(InvalidConfiguration("FirstOrderThermalModel setpoint_temperature_K must lie within [min_temperature_K, max_temperature_K]"))
    model.min_temperature_K <= model.ambient_temperature_K <= model.max_temperature_K ||
        throw(InvalidConfiguration("FirstOrderThermalModel ambient_temperature_K must lie within [min_temperature_K, max_temperature_K]"))
    validate_temperature_law(model.dark_current_law)
    validate_temperature_law(model.glow_rate_law)
    validate_temperature_law(model.dark_count_law)
    validate_temperature_law(model.cic_rate_law)
    return model
end

resolve_thermal_model(sensor::SensorType, ::Nothing; T::Type{<:AbstractFloat}) =
    default_thermal_model(sensor; T=T)
resolve_thermal_model(sensor::SensorType, model::AbstractDetectorThermalModel; T::Type{<:AbstractFloat}) = model

thermal_state_from_model(::NullDetectorThermalModel, ::Type{T}) where {T<:AbstractFloat} = NoThermalState()
thermal_state_from_model(::FixedTemperature, ::Type{T}) where {T<:AbstractFloat} = NoThermalState()
thermal_state_from_model(model::FirstOrderThermalModel, ::Type{T}) where {T<:AbstractFloat} =
    DetectorThermalState{T}(T(model.initial_temperature_K))

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
    background = backend{T}(undef, size(map)...)
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
    kernel_sum = execution_style(model.kernel) isa ScalarCPUStyle ? sum(model.kernel) : sum(Array(model.kernel))
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
convert_frame_timing_model(model::RollingShutter, ::Type{T}) where {T<:AbstractFloat} = RollingShutter{T}(T(model.line_time))

validate_frame_timing_model(::GlobalShutter) = GlobalShutter()

function validate_frame_timing_model(model::RollingShutter)
    model.line_time >= zero(model.line_time) || throw(InvalidConfiguration("RollingShutter line_time must be >= 0"))
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
    response_model::Union{Nothing,FrameResponseModel}, defect_model::Union{Nothing,AbstractDetectorDefectModel},
    timing_model::Union{Nothing,AbstractFrameTimingModel}, correction_model::Union{Nothing,FrameReadoutCorrectionModel},
    nonlinearity_model::Union{Nothing,AbstractFrameNonlinearityModel},
    thermal_model::Union{Nothing,AbstractDetectorThermalModel},
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
    resolved_defect = resolve_detector_defect_model(sensor, defect_model; T=T, backend=backend)
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
    output_precision_t = resolve_output_precision(bits, output_precision)
    window = validate_readout_window(readout_window)
    params = DetectorParams{T, typeof(sensor), typeof(response), typeof(thermal)}(
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
        defects,
        timing,
        correction,
        nonlinearity,
        thermal,
        window,
        output_precision_t,
    )
    frame = backend{T}(undef, 1, 1)
    response_buffer = backend{T}(undef, 1, 1)
    bin_buffer = backend{T}(undef, 1, 1)
    noise_buffer = backend{T}(undef, 1, 1)
    noise_buffer_host = Matrix{T}(undef, 1, 1)
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
    latent_buffer = backend{T}(undef, 1, 1)
    fill!(latent_buffer, zero(T))
    output_buffer === nothing || fill!(output_buffer, zero(eltype(output_buffer)))
    thermal_state = thermal_state_from_model(thermal, T)
    state = DetectorState{T, typeof(frame), typeof(output_buffer), FrameReadoutProducts, typeof(thermal_state)}(
        frame,
        response_buffer,
        bin_buffer,
        noise_buffer,
        noise_buffer_host,
        accum_buffer,
        latent_buffer,
        output_buffer,
        NoFrameReadoutProducts(),
        thermal_state,
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
    defect_model::Union{Nothing,AbstractDetectorDefectModel}=nothing,
    timing_model::Union{Nothing,AbstractFrameTimingModel}=nothing,
    correction_model::Union{Nothing,FrameReadoutCorrectionModel}=nothing,
    nonlinearity_model::Union{Nothing,AbstractFrameNonlinearityModel}=nothing,
    thermal_model::Union{Nothing,AbstractDetectorThermalModel}=nothing,
    readout_window::Union{Nothing,FrameWindow}=nothing,
    output_precision::Union{Nothing,DataType}=nothing, background_flux=nothing, background_map=nothing,
    T::Type{<:AbstractFloat}=Float64, backend=CPUBackend())
    backend = _resolve_array_backend(backend)
    normalized = normalize_noise(noise)
    return Detector(normalized; integration_time=integration_time, qe=qe,
        psf_sampling=psf_sampling, binning=binning, gain=gain,
        dark_current=dark_current, bits=bits, full_well=full_well,
        sensor=sensor, response_model=response_model, defect_model=defect_model, timing_model=timing_model,
        correction_model=correction_model, nonlinearity_model=nonlinearity_model, thermal_model=thermal_model,
        readout_window=readout_window, output_precision=output_precision,
        background_flux=background_flux, background_map=background_map,
        T=T, backend=backend)
end

function Detector(noise::NoiseModel; integration_time::Real=1.0, qe::Real=1.0,
    psf_sampling::Int=1, binning::Int=1, gain::Real=1.0, dark_current::Real=0.0,
    bits::Union{Nothing,Int}=nothing, full_well::Union{Nothing,Real}=nothing,
    sensor::SensorType=CCDSensor(), response_model::Union{Nothing,FrameResponseModel}=nothing,
    defect_model::Union{Nothing,AbstractDetectorDefectModel}=nothing,
    timing_model::Union{Nothing,AbstractFrameTimingModel}=nothing,
    correction_model::Union{Nothing,FrameReadoutCorrectionModel}=nothing,
    nonlinearity_model::Union{Nothing,AbstractFrameNonlinearityModel}=nothing,
    thermal_model::Union{Nothing,AbstractDetectorThermalModel}=nothing,
    readout_window::Union{Nothing,FrameWindow}=nothing,
    output_precision::Union{Nothing,DataType}=nothing,
    background_flux=nothing, background_map=nothing,
    T::Type{<:AbstractFloat}=Float64, backend=CPUBackend())
    backend = _resolve_array_backend(backend)
    converted = convert_noise(noise, T)
    validated = validate_noise(converted)
    return _build_detector(validated; integration_time=integration_time, qe=qe,
        psf_sampling=psf_sampling, binning=binning, gain=gain,
        dark_current=dark_current, bits=bits, full_well=full_well,
        sensor=sensor, response_model=response_model, defect_model=defect_model,
        timing_model=timing_model, correction_model=correction_model, nonlinearity_model=nonlinearity_model,
        thermal_model=thermal_model,
        readout_window=readout_window, output_precision=output_precision,
        background_flux=background_flux, background_map=background_map,
        T=T, backend=backend)
end

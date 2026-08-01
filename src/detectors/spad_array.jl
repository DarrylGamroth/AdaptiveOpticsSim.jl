"""
    SPADArraySensor(; active_area_detection_efficiency=0.5,
        dark_count_rate=0, fill_factor=1, dead_time_model=NoDeadTime(),
        mean_response_model=NullCountingMeanResponse(), ...)

Product-neutral parameters for an accumulated-count SPAD array. Active-area
detection efficiency and fill factor are independent scalar radiometric terms;
neither creates a wavelength response, pixel aperture, or detector MTF.
"""
struct SPADArraySensor{
    T<:AbstractFloat,
    D<:CountingDeadTimeModel,
    C<:CountingMeanResponseModel,
} <: SPADArraySensorType
    active_area_detection_efficiency::T
    dark_count_rate::T
    fill_factor::T
    dead_time_model::D
    mean_response_model::C
end

function SPADArraySensor(; active_area_detection_efficiency::Real=0.5, dark_count_rate::Real=0.0, fill_factor::Real=1.0,
    dead_time_model::CountingDeadTimeModel=NoDeadTime(),
    mean_response_model::CountingMeanResponseModel=NullCountingMeanResponse(),
    T::Type{<:AbstractFloat}=Float64)
    typed_efficiency = T(active_area_detection_efficiency)
    typed_dark_count_rate = T(dark_count_rate)
    typed_fill_factor = T(fill_factor)
    isfinite(typed_efficiency) && zero(T) <= typed_efficiency <= one(T) ||
        throw(InvalidConfiguration("SPADArraySensor active_area_detection_efficiency must be finite and lie in [0, 1]"))
    isfinite(typed_dark_count_rate) && typed_dark_count_rate >= zero(T) ||
        throw(InvalidConfiguration("SPADArraySensor dark_count_rate must be finite and >= 0"))
    isfinite(typed_fill_factor) && zero(T) < typed_fill_factor <= one(T) ||
        throw(InvalidConfiguration("SPADArraySensor fill_factor must be finite and lie in (0, 1]"))
    dead_time = validate_dead_time_model(convert_dead_time_model(dead_time_model, T))
    mean_response = validate_mean_response_model(convert_mean_response_model(mean_response_model, T))
    return SPADArraySensor{T,typeof(dead_time),typeof(mean_response)}(
        typed_efficiency,
        typed_dark_count_rate,
        typed_fill_factor,
        dead_time,
        mean_response,
    )
end

struct SPADArrayDetectorParams{
    T<:AbstractFloat,
    S<:SPADArraySensorType,
    G<:AbstractCountingGateModel,
    TM<:AbstractDetectorThermalModel,
}
    integration_time::T
    dimensions::Tuple{Int,Int}
    gate_model::G
    thermal_model::TM
    sensor::S
    output_type::Union{Nothing,DataType}
end

mutable struct SPADArrayDetectorState{
    T<:AbstractFloat,
    A<:AbstractMatrix{T},
    H<:AbstractMatrix{T},
    O,
    OH,
    TS<:AbstractDetectorThermalState,
}
    counts::A
    noise_buffer::A
    host_buffer::H
    output_buffer::O
    output_buffer_host::OH
    thermal_state::TS
end

struct SPADArrayDetector{
    N<:NoiseModel,
    P<:SPADArrayDetectorParams,
    S<:SPADArrayDetectorState,
    B<:AbstractArrayBackend,
} <: AbstractCountingDetector
    noise::N
    params::P
    state::S
end

@inline backend(::SPADArrayDetector{<:Any,<:Any,<:Any,B}) where {B} = B()

counting_sensor(det::SPADArrayDetector) = det.params.sensor
counting_gate_model(det::SPADArrayDetector) = det.params.gate_model
counting_dead_time_model(det::SPADArrayDetector) = det.params.sensor.dead_time_model
counting_mean_response_model(det::SPADArrayDetector) = det.params.sensor.mean_response_model
counting_integration_time(det::SPADArrayDetector) = det.params.integration_time
counting_layout(::SPADArrayDetector) = :pixel_counts
counting_output_type(det::SPADArrayDetector) = det.params.output_type
counting_array(det::SPADArrayDetector) = det.state.counts
counting_noise_buffer(det::SPADArrayDetector) = det.state.noise_buffer
counting_host_buffer(det::SPADArrayDetector) = det.state.host_buffer
counting_output_buffer(det::SPADArrayDetector) = det.state.output_buffer
counting_output_host_buffer(det::SPADArrayDetector) = det.state.output_buffer_host
set_counting_array!(det::SPADArrayDetector, values) = (det.state.counts = values; det)
set_counting_noise_buffer!(det::SPADArrayDetector, values) = (det.state.noise_buffer = values; det)
set_counting_host_buffer!(det::SPADArrayDetector, values) = (det.state.host_buffer = values; det)
set_counting_output_buffer!(det::SPADArrayDetector, values) = (det.state.output_buffer = values; det)
set_counting_output_host_buffer!(det::SPADArrayDetector, values) =
    (det.state.output_buffer_host = values; det)
counting_detection_efficiency(det::SPADArrayDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} = T(det.params.sensor.active_area_detection_efficiency)
counting_fill_factor(det::SPADArrayDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} = T(det.params.sensor.fill_factor)
counting_reported_fill_factor(det::SPADArrayDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} = T(det.params.sensor.fill_factor)
counting_dark_count_rate(det::SPADArrayDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} = T(det.params.sensor.dark_count_rate)

detector_sensor_symbol(::SPADArraySensor) = :spad_array

dark_count_law(::SPADArrayDetector) = NullTemperatureLaw()

function _build_spad_array_detector(dimensions::Tuple{Int,Int}, noise::NoiseModel;
    integration_time::Real,
    gate_model::AbstractCountingGateModel,
    thermal_model::AbstractDetectorThermalModel,
    sensor::SPADArraySensorType,
    output_type::Union{Nothing,DataType},
    T::Type{<:AbstractFloat},
    backend)
    all(>(0), dimensions) || throw(InvalidConfiguration(
        "SPADArrayDetector dimensions must contain two positive values"))
    typed_integration_time = T(integration_time)
    isfinite(typed_integration_time) && typed_integration_time > zero(T) ||
        throw(InvalidConfiguration("SPADArrayDetector integration_time must be finite and > 0"))
    converted = convert_noise(noise, T)
    validated = validate_counting_noise(converted)
    gate = validate_gate_model(convert_gate_model(gate_model, T))
    thermal = validate_thermal_model(convert_thermal_model(thermal_model, T))
    typed_sensor = convert_spad_sensor(sensor, T)
    params = SPADArrayDetectorParams{T,typeof(typed_sensor),typeof(gate),typeof(thermal)}(
        typed_integration_time,
        dimensions,
        gate,
        thermal,
        typed_sensor,
        output_type,
    )
    counts = backend{T}(undef, dimensions...)
    noise_buffer = backend{T}(undef, dimensions...)
    host_buffer = Matrix{T}(undef, dimensions...)
    output_buffer = output_type === nothing ? nothing : backend{output_type}(undef, dimensions...)
    output_buffer_host = output_type === nothing ? nothing : Matrix{output_type}(undef, dimensions...)
    fill!(counts, zero(T))
    fill!(noise_buffer, zero(T))
    fill!(host_buffer, zero(T))
    output_buffer === nothing || fill!(output_buffer, zero(eltype(output_buffer)))
    output_buffer_host === nothing || fill!(output_buffer_host,
        zero(eltype(output_buffer_host)))
    thermal_state = thermal_state_from_model(thermal, T)
    state = SPADArrayDetectorState{T,typeof(counts),typeof(host_buffer),
        typeof(output_buffer),typeof(output_buffer_host),typeof(thermal_state)}(
        counts, noise_buffer, host_buffer, output_buffer, output_buffer_host,
        thermal_state)
    selector = _resolve_backend_selector(backend)
    return SPADArrayDetector{typeof(validated),typeof(params),typeof(state),typeof(selector)}(
        validated, params, state)
end

"""
    SPADArrayDetector((rows, columns); sensor=SPADArraySensor(), ...)

Construct a fixed-shape Geiger-mode accumulated-count area detector. Input
matrices contain already cell-integrated photon-arrival rates. The model emits
an expected-count or sampled-count image, not photon events or timestamps.
"""
function SPADArrayDetector(dimensions::Tuple{Int,Int};
    integration_time::Real=1.0, noise::NoiseModel=NoisePhoton(),
    sensor::SPADArraySensor=SPADArraySensor(),
    output_type::Union{Nothing,DataType}=nothing,
    gate_model::AbstractCountingGateModel=NullCountingGate(),
    thermal_model::AbstractDetectorThermalModel=NullDetectorThermalModel(),
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=CPUBackend())
    backend = _resolve_array_backend(backend)
    return _build_spad_array_detector(dimensions, noise; integration_time=integration_time,
        gate_model=gate_model, thermal_model=thermal_model, sensor=sensor,
        output_type=output_type, T=T, backend=backend)
end

function ensure_buffers!(det::SPADArrayDetector, dimensions::Tuple{Int,Int})
    dimensions == det.params.dimensions || throw(DimensionMismatchError(
        "SPADArrayDetector input dimensions must match its fixed detector dimensions"))
    return det
end

convert_spad_sensor(sensor::SPADArraySensor{T}, ::Type{T}) where {T<:AbstractFloat} = sensor

function convert_spad_sensor(sensor::SPADArraySensorType, ::Type{T}) where {T<:AbstractFloat}
    return SPADArraySensor(
        active_area_detection_efficiency=sensor.active_area_detection_efficiency,
        dark_count_rate=sensor.dark_count_rate,
        fill_factor=sensor.fill_factor,
        dead_time_model=sensor.dead_time_model,
        mean_response_model=sensor.mean_response_model,
        T=T,
    )
end

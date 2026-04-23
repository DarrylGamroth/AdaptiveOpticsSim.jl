struct SPADArraySensor{
    T<:AbstractFloat,
    D<:CountingDeadTimeModel,
    C<:AbstractCountingCorrelationModel,
} <: SPADArraySensorType
    pde::T
    dark_count_rate::T
    fill_factor::T
    dead_time_model::D
    correlation_model::C
end

function SPADArraySensor(; pde::Real=0.5, dark_count_rate::Real=0.0, fill_factor::Real=1.0,
    dead_time_model::CountingDeadTimeModel=NoDeadTime(),
    correlation_model::AbstractCountingCorrelationModel=NullCountingCorrelation(),
    T::Type{<:AbstractFloat}=Float64)
    zero(pde) <= pde <= one(pde) ||
        throw(InvalidConfiguration("SPADArraySensor pde must lie in [0, 1]"))
    dark_count_rate >= 0 ||
        throw(InvalidConfiguration("SPADArraySensor dark_count_rate must be >= 0"))
    zero(fill_factor) < fill_factor <= one(fill_factor) ||
        throw(InvalidConfiguration("SPADArraySensor fill_factor must lie in (0, 1]"))
    dead_time = validate_dead_time_model(convert_dead_time_model(dead_time_model, T))
    correlation = validate_correlation_model(convert_correlation_model(correlation_model, T))
    return SPADArraySensor{T,typeof(dead_time),typeof(correlation)}(
        T(pde),
        T(dark_count_rate),
        T(fill_factor),
        dead_time,
        correlation,
    )
end

struct SPADArrayDetectorParams{
    T<:AbstractFloat,
    S<:SPADArraySensorType,
    G<:AbstractCountingGateModel,
    TM<:AbstractDetectorThermalModel,
}
    integration_time::T
    gate_model::G
    thermal_model::TM
    sensor::S
    output_precision::Union{Nothing,DataType}
    layout::Symbol
end

mutable struct SPADArrayDetectorState{
    T<:AbstractFloat,
    A<:AbstractMatrix{T},
    O,
    TS<:AbstractDetectorThermalState,
}
    counts::A
    noise_buffer::A
    output_buffer::O
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
counting_correlation_model(det::SPADArrayDetector) = det.params.sensor.correlation_model
counting_integration_time(det::SPADArrayDetector) = det.params.integration_time
counting_layout(det::SPADArrayDetector) = det.params.layout
counting_output_precision(det::SPADArrayDetector) = det.params.output_precision
counting_array(det::SPADArrayDetector) = det.state.counts
counting_noise_buffer(det::SPADArrayDetector) = det.state.noise_buffer
counting_output_buffer(det::SPADArrayDetector) = det.state.output_buffer
set_counting_array!(det::SPADArrayDetector, values) = (det.state.counts = values; det)
set_counting_noise_buffer!(det::SPADArrayDetector, values) = (det.state.noise_buffer = values; det)
set_counting_output_buffer!(det::SPADArrayDetector, values) = (det.state.output_buffer = values; det)
counting_qe(det::SPADArrayDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} = T(det.params.sensor.pde)
counting_fill_factor(det::SPADArrayDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} = T(det.params.sensor.fill_factor)
counting_reported_fill_factor(det::SPADArrayDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} = T(det.params.sensor.fill_factor)
counting_dark_count_rate(det::SPADArrayDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} = T(det.params.sensor.dark_count_rate)

detector_sensor_symbol(::SPADArraySensor) = :spad_array

dark_count_law(::SPADArrayDetector) = NullTemperatureLaw()

function _build_spad_array_detector(noise::NoiseModel; integration_time::Real,
    gate_model::AbstractCountingGateModel,
    thermal_model::AbstractDetectorThermalModel,
    sensor::SPADArraySensorType,
    output_precision::Union{Nothing,DataType},
    layout::Symbol,
    T::Type{<:AbstractFloat},
    backend)
    integration_time > 0 || throw(InvalidConfiguration("SPADArrayDetector integration_time must be > 0"))
    converted = convert_noise(noise, T)
    validated = validate_counting_noise(converted)
    gate = validate_gate_model(convert_gate_model(gate_model, T))
    thermal = validate_thermal_model(convert_thermal_model(thermal_model, T))
    typed_sensor = convert_spad_sensor(sensor, T)
    params = SPADArrayDetectorParams{T,typeof(typed_sensor),typeof(gate),typeof(thermal)}(
        T(integration_time),
        gate,
        thermal,
        typed_sensor,
        output_precision,
        layout,
    )
    counts = backend{T}(undef, 1, 1)
    noise_buffer = backend{T}(undef, 1, 1)
    output_buffer = output_precision === nothing ? nothing : backend{output_precision}(undef, 1, 1)
    fill!(counts, zero(T))
    fill!(noise_buffer, zero(T))
    output_buffer === nothing || fill!(output_buffer, zero(eltype(output_buffer)))
    thermal_state = thermal_state_from_model(thermal, T)
    state = SPADArrayDetectorState{T,typeof(counts),typeof(output_buffer),typeof(thermal_state)}(
        counts, noise_buffer, output_buffer, thermal_state)
    selector = _resolve_backend_selector(backend)
    return SPADArrayDetector{typeof(validated),typeof(params),typeof(state),typeof(selector)}(
        validated, params, state)
end

function SPADArrayDetector(; integration_time::Real=1.0, noise::NoiseModel=NoisePhoton(),
    sensor::SPADArraySensor=SPADArraySensor(),
    output_precision::Union{Nothing,DataType}=nothing,
    layout::Symbol=:pixel_counts,
    gate_model::AbstractCountingGateModel=NullCountingGate(),
    thermal_model::AbstractDetectorThermalModel=NullDetectorThermalModel(),
    T::Type{<:AbstractFloat}=Float64,
    backend::AbstractArrayBackend=CPUBackend())
    backend = _resolve_array_backend(backend)
    return _build_spad_array_detector(noise; integration_time=integration_time,
        gate_model=gate_model, thermal_model=thermal_model, sensor=sensor,
        output_precision=output_precision, layout=layout, T=T, backend=backend)
end

convert_spad_sensor(sensor::SPADArraySensor{T}, ::Type{T}) where {T<:AbstractFloat} = sensor

function convert_spad_sensor(sensor::SPADArraySensorType, ::Type{T}) where {T<:AbstractFloat}
    return SPADArraySensor(
        pde=sensor.pde,
        dark_count_rate=sensor.dark_count_rate,
        fill_factor=sensor.fill_factor,
        dead_time_model=sensor.dead_time_model,
        correlation_model=sensor.correlation_model,
        T=T,
    )
end

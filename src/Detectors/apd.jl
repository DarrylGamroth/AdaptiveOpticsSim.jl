struct APDSensor <: CountingSensorType end

struct APDDetectorParams{T<:AbstractFloat,S<:APDSensor,D<:CountingDeadTimeModel,G<:AbstractCountingGateModel,C<:AbstractCountingCorrelationModel,TM<:AbstractDetectorThermalModel}
    integration_time::T
    qe::T
    gain::T
    dark_count_rate::T
    dead_time_model::D
    gate_model::G
    correlation_model::C
    thermal_model::TM
    sensor::S
    output_precision::Union{Nothing,DataType}
    layout::Symbol
end

mutable struct APDDetectorState{T<:AbstractFloat,A<:AbstractMatrix{T},O,TS<:AbstractDetectorThermalState}
    channels::A
    noise_buffer::A
    output_buffer::O
    thermal_state::TS
end

struct APDDetector{N<:NoiseModel,P<:APDDetectorParams,S<:APDDetectorState,GM,B<:AbstractArrayBackend} <: AbstractCountingDetector
    noise::N
    params::P
    state::S
    channel_gain_map::GM
end

@inline backend(::APDDetector{<:Any,<:Any,<:Any,<:Any,B}) where {B} = B()

counting_sensor(det::APDDetector) = det.params.sensor
counting_gate_model(det::APDDetector) = det.params.gate_model
counting_dead_time_model(det::APDDetector) = det.params.dead_time_model
counting_correlation_model(det::APDDetector) = det.params.correlation_model
counting_integration_time(det::APDDetector) = det.params.integration_time
counting_layout(det::APDDetector) = det.params.layout
counting_output_precision(det::APDDetector) = det.params.output_precision
counting_array(det::APDDetector) = det.state.channels
counting_noise_buffer(det::APDDetector) = det.state.noise_buffer
counting_output_buffer(det::APDDetector) = det.state.output_buffer
set_counting_array!(det::APDDetector, values) = (det.state.channels = values; det)
set_counting_noise_buffer!(det::APDDetector, values) = (det.state.noise_buffer = values; det)
set_counting_output_buffer!(det::APDDetector, values) = (det.state.output_buffer = values; det)
counting_qe(det::APDDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} = T(det.params.qe)
counting_post_gain(det::APDDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} = T(det.params.gain)
counting_dark_count_rate(det::APDDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} = T(det.params.dark_count_rate)
counting_channel_gain_map(det::APDDetector) = det.channel_gain_map

detector_sensor_symbol(::APDSensor) = :apd
dark_count_law(::APDDetector) = NullTemperatureLaw()

function _build_apd_detector(noise::NoiseModel; integration_time::Real, qe::Real, gain::Real,
    dark_count_rate::Real, dead_time_model::CountingDeadTimeModel, gate_model::AbstractCountingGateModel,
    correlation_model::AbstractCountingCorrelationModel, thermal_model::AbstractDetectorThermalModel,
    sensor::APDSensor, output_precision::Union{Nothing,DataType},
    layout::Symbol, channel_gain_map,
    T::Type{<:AbstractFloat}, backend)
    gain >= 0 || throw(InvalidConfiguration("APDDetector gain must be >= 0"))
    dark_count_rate >= 0 || throw(InvalidConfiguration("APDDetector dark_count_rate must be >= 0"))
    converted = convert_noise(noise, T)
    validated = validate_counting_noise(converted)
    dead_time = validate_dead_time_model(convert_dead_time_model(dead_time_model, T))
    gate = validate_gate_model(convert_gate_model(gate_model, T))
    correlation = validate_correlation_model(convert_correlation_model(correlation_model, T))
    thermal = validate_thermal_model(convert_thermal_model(thermal_model, T))
    gain_map = channel_gain_map === nothing ? nothing : begin
        g = backend{T}(undef, size(channel_gain_map)...)
        copyto!(g, T.(channel_gain_map))
        g
    end
    params = APDDetectorParams{T,typeof(sensor),typeof(dead_time),typeof(gate),typeof(correlation),typeof(thermal)}(
        T(integration_time),
        T(qe),
        T(gain),
        T(dark_count_rate),
        dead_time,
        gate,
        correlation,
        thermal,
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
    thermal_state = thermal_state_from_model(thermal, T)
    state = APDDetectorState{T,typeof(channels),typeof(output_buffer),typeof(thermal_state)}(
        channels, noise_buffer, output_buffer, thermal_state)
    selector = _resolve_backend_selector(backend)
    return APDDetector{typeof(validated),typeof(params),typeof(state),typeof(gain_map),typeof(selector)}(
        validated, params, state, gain_map)
end

function APDDetector(; integration_time::Real=1.0, qe::Real=1.0, noise::NoiseModel=NoisePhoton(),
    gain::Real=1.0, dark_count_rate::Real=0.0, output_precision::Union{Nothing,DataType}=nothing,
    layout::Symbol=:channels, channel_gain_map=nothing, dead_time_model::CountingDeadTimeModel=NoDeadTime(),
    gate_model::AbstractCountingGateModel=NullCountingGate(),
    correlation_model::AbstractCountingCorrelationModel=NullCountingCorrelation(),
    thermal_model::AbstractDetectorThermalModel=NullDetectorThermalModel(),
    T::Type{<:AbstractFloat}=Float64, backend::AbstractArrayBackend=CPUBackend())
    backend = _resolve_array_backend(backend)
    return _build_apd_detector(noise; integration_time=integration_time, qe=qe, gain=gain,
        dark_count_rate=dark_count_rate, dead_time_model=dead_time_model, gate_model=gate_model,
        correlation_model=correlation_model, thermal_model=thermal_model, sensor=APDSensor(), output_precision=output_precision,
        layout=layout, channel_gain_map=channel_gain_map, T=T, backend=backend)
end

function apply_counting_channel_gain_map!(det::APDDetector)
    isnothing(det.channel_gain_map) && return det.state.channels
    size(det.channel_gain_map) == size(det.state.channels) ||
        throw(DimensionMismatchError("APDDetector channel_gain_map size must match counting-channel size"))
    det.state.channels .*= det.channel_gain_map
    return det.state.channels
end

function apply_post_counting_gain!(det::APDDetector)
    det.state.channels .*= det.params.gain
    return det.state.channels
end

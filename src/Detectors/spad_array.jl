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

readout_ready(det::SPADArrayDetector) = true
thermal_model(det::SPADArrayDetector) = det.params.thermal_model
thermal_state(det::SPADArrayDetector) = det.state.thermal_state
detector_temperature(det::SPADArrayDetector, ::Type{T}=eltype(det.state.counts)) where {T<:AbstractFloat} =
    detector_temperature_K(det.params.thermal_model, det.state.thermal_state, T)
advance_thermal!(det::SPADArrayDetector, dt) = (advance_thermal!(det.params.thermal_model, det.state.thermal_state, dt); det)
channel_output(det::SPADArrayDetector) = det.state.output_buffer === nothing ? det.state.counts : det.state.output_buffer
output_frame(det::SPADArrayDetector) = channel_output(det)
reset_integration!(det::SPADArrayDetector) = det

detector_sensor_symbol(::SPADArraySensor) = :spad_array
supports_counting_noise(::SPADArrayDetector) = true
supports_dead_time(det::SPADArrayDetector) = supports_dead_time(det.params.sensor.dead_time_model)
supports_channel_gain_map(::SPADArrayDetector) = false
supports_counting_gating(det::SPADArrayDetector) = !is_null_counting_gate(det.params.gate_model)
supports_afterpulsing(det::SPADArrayDetector) = _supports_afterpulsing(det.params.sensor.correlation_model)
supports_channel_crosstalk(det::SPADArrayDetector) = _supports_channel_crosstalk(det.params.sensor.correlation_model)
supports_paralyzable_dead_time(det::SPADArrayDetector) = is_paralyzable_dead_time(det.params.sensor.dead_time_model)
supports_detector_thermal_model(det::SPADArrayDetector) = !is_null_thermal_model(det.params.thermal_model)
supports_temperature_dependent_dark_counts(det::SPADArrayDetector) =
    !is_null_temperature_law(active_dark_count_law(det, det.params.thermal_model))

dark_count_law(::SPADArrayDetector) = NullTemperatureLaw()
active_dark_count_law(det::SPADArrayDetector, ::NullDetectorThermalModel) = dark_count_law(det)

function active_dark_count_law(det::SPADArrayDetector, model::FixedTemperature)
    return is_null_temperature_law(model.dark_count_law) ? dark_count_law(det) : model.dark_count_law
end

function active_dark_count_law(det::SPADArrayDetector, model::FirstOrderThermalModel)
    return is_null_temperature_law(model.dark_count_law) ? dark_count_law(det) : model.dark_count_law
end

effective_dark_count_rate(det::SPADArrayDetector, ::Type{T}=eltype(det.state.counts)) where {T<:AbstractFloat} =
    T(evaluate_temperature_law(active_dark_count_law(det, det.params.thermal_model),
        T(det.params.sensor.dark_count_rate), detector_temperature(det, T)))

validate_spad_noise(noise::NoiseNone) = noise
validate_spad_noise(noise::NoisePhoton) = noise
validate_spad_noise(::NoiseReadout) =
    throw(InvalidConfiguration("SPADArrayDetector does not support additive readout noise; use NoiseNone or NoisePhoton"))
validate_spad_noise(::NoisePhotonReadout) =
    throw(InvalidConfiguration("SPADArrayDetector does not support frame-style readout noise; use NoiseNone or NoisePhoton"))

function detector_export_metadata(det::SPADArrayDetector; T::Type{<:AbstractFloat}=eltype(det.state.counts))
    output = output_frame(det)
    sensor = det.params.sensor
    return CountingDetectorExportMetadata{T}(
        T(det.params.integration_time),
        T(sensor.pde),
        T(sensor.fill_factor),
        one(T),
        T(sensor.dark_count_rate),
        counting_dead_time_symbol(sensor.dead_time_model),
        counting_dead_time_value(sensor.dead_time_model, T),
        counting_gate_symbol(det.params.gate_model),
        counting_gate_duty_cycle(det.params.gate_model, T),
        counting_correlation_symbol(sensor.correlation_model),
        afterpulse_probability(sensor.correlation_model, T),
        crosstalk_value(sensor.correlation_model, T),
        thermal_model_symbol(det.params.thermal_model),
        detector_temperature(det, T),
        ambient_temperature_K(det.params.thermal_model, T),
        cooling_setpoint_K(det.params.thermal_model, T),
        thermal_time_constant_s(det.params.thermal_model, T),
        temperature_law_symbol(active_dark_count_law(det, det.params.thermal_model)),
        detector_sensor_symbol(sensor),
        detector_noise_symbol(det.noise),
        det.params.output_precision,
        CountingReadoutMetadata(det.params.layout, size(output), length(output)),
    )
end

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
    validated = validate_spad_noise(converted)
    gate = validate_gate_model(convert_gate_model(gate_model, T))
    thermal = validate_thermal_model(convert_thermal_model(thermal_model, T))
    typed_sensor = sensor isa SPADArraySensor{T} ? sensor : SPADArraySensor(
        pde=sensor.pde,
        dark_count_rate=sensor.dark_count_rate,
        fill_factor=sensor.fill_factor,
        dead_time_model=sensor.dead_time_model,
        correlation_model=sensor.correlation_model,
        T=T,
    )
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

function ensure_buffers!(det::SPADArrayDetector, dims::Tuple{Int,Int})
    if size(det.state.counts) != dims
        det.state.counts = similar(det.state.counts, dims...)
    end
    if size(det.state.noise_buffer) != dims
        det.state.noise_buffer = similar(det.state.noise_buffer, dims...)
    end
    if det.state.output_buffer !== nothing && size(det.state.output_buffer) != dims
        det.state.output_buffer = similar(det.state.output_buffer, dims...)
        fill!(det.state.output_buffer, zero(eltype(det.state.output_buffer)))
    end
    return det
end

counting_exposure_time(det::SPADArrayDetector) = effective_gate_time(det.params.gate_model, det.params.integration_time)

function apply_dark_counts!(det::SPADArrayDetector, exposure_time::Real)
    dark = effective_dark_count_rate(det) * exposure_time
    dark <= 0 && return det.state.counts
    det.state.counts .+= dark
    return det.state.counts
end

apply_dead_time!(det::SPADArrayDetector) = apply_dead_time!(det.params.sensor.dead_time_model, det)
apply_dead_time!(::NoDeadTime, det::SPADArrayDetector) = det.state.counts

function apply_dead_time!(model::NonParalyzableDeadTime, det::SPADArrayDetector)
    exposure_time = counting_exposure_time(det)
    exposure_time > zero(exposure_time) || return det.state.counts
    scale = model.dead_time / exposure_time
    scale <= zero(scale) && return det.state.counts
    @. det.state.counts = det.state.counts / (1 + det.state.counts * scale)
    return det.state.counts
end

function apply_dead_time!(model::ParalyzableDeadTime, det::SPADArrayDetector)
    exposure_time = counting_exposure_time(det)
    exposure_time > zero(exposure_time) || return det.state.counts
    scale = model.dead_time / exposure_time
    scale <= zero(scale) && return det.state.counts
    @. det.state.counts = det.state.counts * exp(-det.state.counts * scale)
    return det.state.counts
end

apply_counting_noise!(det::SPADArrayDetector{NoiseNone}, rng::AbstractRNG) = det.state.counts

function apply_counting_noise!(det::SPADArrayDetector{NoisePhoton}, rng::AbstractRNG)
    poisson_noise!(rng, det.state.counts)
    return det.state.counts
end

apply_counting_correlation!(::NullCountingCorrelation, det::SPADArrayDetector, rng::AbstractRNG) = det.state.counts

function apply_counting_correlation!(model::AfterpulsingModel, det::SPADArrayDetector, rng::AbstractRNG)
    p = model.probability
    p <= zero(p) && return det.state.counts
    det.state.counts .*= (one(p) + p)
    return det.state.counts
end

function apply_counting_correlation!(model::ChannelCrosstalkModel, det::SPADArrayDetector, rng::AbstractRNG)
    coupling = model.coupling
    coupling <= zero(coupling) && return det.state.counts
    copyto!(det.state.noise_buffer, det.state.counts)
    fill!(det.state.counts, zero(eltype(det.state.counts)))
    n, m = size(det.state.counts)
    @inbounds for i in 1:n, j in 1:m
        center = det.state.noise_buffer[i, j]
        bleed = coupling * center
        keep = center - bleed
        det.state.counts[i, j] += keep
        neighbors = 0
        i > 1 && (neighbors += 1)
        i < n && (neighbors += 1)
        j > 1 && (neighbors += 1)
        j < m && (neighbors += 1)
        neighbors == 0 && continue
        share = bleed / neighbors
        i > 1 && (det.state.counts[i - 1, j] += share)
        i < n && (det.state.counts[i + 1, j] += share)
        j > 1 && (det.state.counts[i, j - 1] += share)
        j < m && (det.state.counts[i, j + 1] += share)
    end
    return det.state.counts
end

function apply_counting_correlation!(model::CompositeCountingCorrelation, det::SPADArrayDetector, rng::AbstractRNG)
    foreach(stage -> apply_counting_correlation!(stage, det, rng), model.stages)
    return det.state.counts
end

function write_output!(det::SPADArrayDetector)
    output = det.state.output_buffer
    output === nothing && return det.state.counts
    output_eltype = eltype(output)
    if output_eltype <: Integer
        lo = typemin(output_eltype)
        hi = typemax(output_eltype)
        @. output = output_eltype(clamp(round(det.state.counts), lo, hi))
    else
        copyto!(output, det.state.counts)
    end
    return output
end

function capture!(det::SPADArrayDetector, channels::AbstractMatrix{T}; rng::AbstractRNG=Random.default_rng()) where {T<:AbstractFloat}
    ensure_buffers!(det, size(channels))
    copyto!(det.state.counts, channels)
    exposure_time = counting_exposure_time(det)
    det.state.counts .*= det.params.sensor.pde * det.params.sensor.fill_factor * exposure_time
    apply_dark_counts!(det, exposure_time)
    apply_counting_noise!(det, rng)
    apply_dead_time!(det)
    apply_counting_correlation!(det.params.sensor.correlation_model, det, rng)
    advance_thermal!(det, det.params.integration_time)
    return write_output!(det)
end

capture!(det::SPADArrayDetector, channels::AbstractMatrix{T}, rng::AbstractRNG) where {T<:AbstractFloat} =
    capture!(det, channels; rng=rng)

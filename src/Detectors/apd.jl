struct APDSensor <: CountingSensorType end
struct NoDeadTime <: CountingDeadTimeModel end

struct NonParalyzableDeadTime{T<:AbstractFloat} <: CountingDeadTimeModel
    dead_time::T
end

struct ParalyzableDeadTime{T<:AbstractFloat} <: CountingDeadTimeModel
    dead_time::T
end

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

struct APDDetector{N<:NoiseModel,P<:APDDetectorParams,S<:APDDetectorState,GM} <: AbstractCountingDetector
    noise::N
    params::P
    state::S
    channel_gain_map::GM
end

NonParalyzableDeadTime(dead_time::Real) = NonParalyzableDeadTime{Float64}(float(dead_time))
ParalyzableDeadTime(dead_time::Real) = ParalyzableDeadTime{Float64}(float(dead_time))

readout_ready(det::APDDetector) = true
thermal_model(det::APDDetector) = det.params.thermal_model
thermal_state(det::APDDetector) = det.state.thermal_state
detector_temperature(det::APDDetector, ::Type{T}=eltype(det.state.channels)) where {T<:AbstractFloat} =
    detector_temperature_K(det.params.thermal_model, det.state.thermal_state, T)
advance_thermal!(det::APDDetector, dt) = (advance_thermal!(det.params.thermal_model, det.state.thermal_state, dt); det)
channel_output(det::APDDetector) = det.state.output_buffer === nothing ? det.state.channels : det.state.output_buffer
output_frame(det::APDDetector) = channel_output(det)
reset_integration!(det::APDDetector) = det

detector_sensor_symbol(::APDSensor) = :apd
counting_dead_time_symbol(::NoDeadTime) = :none
counting_dead_time_symbol(::NonParalyzableDeadTime) = :nonparalyzable
counting_dead_time_symbol(::ParalyzableDeadTime) = :paralyzable
supports_counting_noise(::APDDetector) = true
supports_dead_time(det::APDDetector) = supports_dead_time(det.params.dead_time_model)
supports_dead_time(::NoDeadTime) = false
supports_dead_time(::NonParalyzableDeadTime) = true
supports_dead_time(::ParalyzableDeadTime) = true
supports_channel_gain_map(det::APDDetector) = !isnothing(det.channel_gain_map)
is_null_counting_gate(::AbstractCountingGateModel) = false
is_null_counting_gate(::NullCountingGate) = true
is_paralyzable_dead_time(::CountingDeadTimeModel) = false
is_paralyzable_dead_time(::ParalyzableDeadTime) = true
supports_counting_gating(det::APDDetector) = !is_null_counting_gate(det.params.gate_model)
supports_afterpulsing(det::APDDetector) = _supports_afterpulsing(det.params.correlation_model)
supports_channel_crosstalk(det::APDDetector) = _supports_channel_crosstalk(det.params.correlation_model)
supports_paralyzable_dead_time(det::APDDetector) = is_paralyzable_dead_time(det.params.dead_time_model)
supports_detector_thermal_model(det::APDDetector) = !is_null_thermal_model(det.params.thermal_model)
supports_temperature_dependent_dark_counts(det::APDDetector) =
    !is_null_temperature_law(active_dark_count_law(det, det.params.thermal_model))

_supports_afterpulsing(::AbstractCountingCorrelationModel) = false
_supports_afterpulsing(::AfterpulsingModel) = true
_supports_afterpulsing(model::CompositeCountingCorrelation) = any(_supports_afterpulsing, model.stages)
_supports_channel_crosstalk(::AbstractCountingCorrelationModel) = false
_supports_channel_crosstalk(::ChannelCrosstalkModel) = true
_supports_channel_crosstalk(model::CompositeCountingCorrelation) = any(_supports_channel_crosstalk, model.stages)

counting_dead_time_value(::NoDeadTime, ::Type{T}) where {T<:AbstractFloat} = nothing
counting_dead_time_value(model::NonParalyzableDeadTime, ::Type{T}) where {T<:AbstractFloat} = T(model.dead_time)
counting_dead_time_value(model::ParalyzableDeadTime, ::Type{T}) where {T<:AbstractFloat} = T(model.dead_time)
convert_dead_time_model(::NoDeadTime, ::Type{T}) where {T<:AbstractFloat} = NoDeadTime()
convert_dead_time_model(model::NonParalyzableDeadTime, ::Type{T}) where {T<:AbstractFloat} =
    NonParalyzableDeadTime{T}(T(model.dead_time))
convert_dead_time_model(model::ParalyzableDeadTime, ::Type{T}) where {T<:AbstractFloat} =
    ParalyzableDeadTime{T}(T(model.dead_time))
validate_dead_time_model(model::NoDeadTime) = model

function validate_dead_time_model(model::NonParalyzableDeadTime)
    model.dead_time >= 0 || throw(InvalidConfiguration("NonParalyzableDeadTime dead_time must be >= 0"))
    return model
end

function validate_dead_time_model(model::ParalyzableDeadTime)
    model.dead_time >= 0 || throw(InvalidConfiguration("ParalyzableDeadTime dead_time must be >= 0"))
    return model
end

convert_gate_model(::NullCountingGate, ::Type{T}) where {T<:AbstractFloat} = NullCountingGate()
convert_gate_model(model::DutyCycleGate, ::Type{T}) where {T<:AbstractFloat} = DutyCycleGate{T}(T(model.duty_cycle))
validate_gate_model(::NullCountingGate) = NullCountingGate()

function validate_gate_model(model::DutyCycleGate)
    zero(model.duty_cycle) < model.duty_cycle <= one(model.duty_cycle) ||
        throw(InvalidConfiguration("DutyCycleGate duty_cycle must lie in (0, 1]"))
    return model
end

convert_correlation_model(::NullCountingCorrelation, ::Type{T}) where {T<:AbstractFloat} = NullCountingCorrelation()
convert_correlation_model(model::AfterpulsingModel, ::Type{T}) where {T<:AbstractFloat} = AfterpulsingModel{T}(T(model.probability))
convert_correlation_model(model::ChannelCrosstalkModel, ::Type{T}) where {T<:AbstractFloat} =
    ChannelCrosstalkModel{T}(T(model.coupling))
convert_correlation_model(model::CompositeCountingCorrelation, ::Type{T}) where {T<:AbstractFloat} =
    CompositeCountingCorrelation(tuple((convert_correlation_model(stage, T) for stage in model.stages)...))

dark_count_law(::APDDetector) = NullTemperatureLaw()
active_dark_count_law(det::APDDetector, ::NullDetectorThermalModel) = dark_count_law(det)

function active_dark_count_law(det::APDDetector, model::FixedTemperature)
    return is_null_temperature_law(model.dark_count_law) ? dark_count_law(det) : model.dark_count_law
end

effective_dark_count_rate(det::APDDetector, ::Type{T}=eltype(det.state.channels)) where {T<:AbstractFloat} =
    T(evaluate_temperature_law(active_dark_count_law(det, det.params.thermal_model),
        T(det.params.dark_count_rate), detector_temperature(det, T)))

validate_correlation_model(::NullCountingCorrelation) = NullCountingCorrelation()

function validate_correlation_model(model::AfterpulsingModel)
    zero(model.probability) <= model.probability <= one(model.probability) ||
        throw(InvalidConfiguration("AfterpulsingModel probability must lie in [0, 1]"))
    return model
end

function validate_correlation_model(model::ChannelCrosstalkModel)
    zero(model.coupling) <= model.coupling <= one(model.coupling) ||
        throw(InvalidConfiguration("ChannelCrosstalkModel coupling must lie in [0, 1]"))
    return model
end

function validate_correlation_model(model::CompositeCountingCorrelation)
    return CompositeCountingCorrelation(tuple((validate_correlation_model(stage) for stage in model.stages)...))
end

validate_apd_noise(noise::NoiseNone) = noise
validate_apd_noise(noise::NoisePhoton) = noise
validate_apd_noise(noise::NoiseReadout) =
    throw(InvalidConfiguration("APDDetector does not support additive readout noise; use NoiseNone or NoisePhoton"))
validate_apd_noise(noise::NoisePhotonReadout) =
    throw(InvalidConfiguration("APDDetector does not support frame-style readout noise; use NoisePhoton"))

function detector_export_metadata(det::APDDetector; T::Type{<:AbstractFloat}=eltype(det.state.channels))
    output = channel_output(det)
    return CountingDetectorExportMetadata{T}(
        T(det.params.integration_time),
        T(det.params.qe),
        T(det.params.gain),
        T(det.params.dark_count_rate),
        counting_dead_time_symbol(det.params.dead_time_model),
        counting_dead_time_value(det.params.dead_time_model, T),
        counting_gate_symbol(det.params.gate_model),
        counting_gate_duty_cycle(det.params.gate_model, T),
        counting_correlation_symbol(det.params.correlation_model),
        afterpulse_probability(det.params.correlation_model, T),
        crosstalk_value(det.params.correlation_model, T),
        thermal_model_symbol(det.params.thermal_model),
        detector_temperature(det, T),
        ambient_temperature_K(det.params.thermal_model, T),
        cooling_setpoint_K(det.params.thermal_model, T),
        thermal_time_constant_s(det.params.thermal_model, T),
        temperature_law_symbol(active_dark_count_law(det, det.params.thermal_model)),
        detector_sensor_symbol(det.params.sensor),
        detector_noise_symbol(det.noise),
        det.params.output_precision,
        CountingReadoutMetadata(det.params.layout, size(output), length(output)),
    )
end

function _build_apd_detector(noise::NoiseModel; integration_time::Real, qe::Real, gain::Real,
    dark_count_rate::Real, dead_time_model::CountingDeadTimeModel, gate_model::AbstractCountingGateModel,
    correlation_model::AbstractCountingCorrelationModel, thermal_model::AbstractDetectorThermalModel,
    sensor::APDSensor, output_precision::Union{Nothing,DataType},
    layout::Symbol, channel_gain_map,
    T::Type{<:AbstractFloat}, backend)
    gain >= 0 || throw(InvalidConfiguration("APDDetector gain must be >= 0"))
    dark_count_rate >= 0 || throw(InvalidConfiguration("APDDetector dark_count_rate must be >= 0"))
    converted = convert_noise(noise, T)
    validated = validate_apd_noise(converted)
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
    return APDDetector{typeof(validated),typeof(params),typeof(state),typeof(gain_map)}(
        validated, params, state, gain_map)
end

function APDDetector(; integration_time::Real=1.0, qe::Real=1.0, noise::NoiseModel=NoisePhoton(),
    gain::Real=1.0, dark_count_rate::Real=0.0, output_precision::Union{Nothing,DataType}=nothing,
    layout::Symbol=:channels, channel_gain_map=nothing, dead_time_model::CountingDeadTimeModel=NoDeadTime(),
    gate_model::AbstractCountingGateModel=NullCountingGate(),
    correlation_model::AbstractCountingCorrelationModel=NullCountingCorrelation(),
    thermal_model::AbstractDetectorThermalModel=NullDetectorThermalModel(),
    T::Type{<:AbstractFloat}=Float64, backend=Array)
    return _build_apd_detector(noise; integration_time=integration_time, qe=qe, gain=gain,
        dark_count_rate=dark_count_rate, dead_time_model=dead_time_model, gate_model=gate_model,
        correlation_model=correlation_model, thermal_model=thermal_model, sensor=APDSensor(), output_precision=output_precision,
        layout=layout, channel_gain_map=channel_gain_map, T=T, backend=backend)
end

function ensure_buffers!(det::APDDetector, dims::Tuple{Int,Int})
    if size(det.state.channels) != dims
        det.state.channels = similar(det.state.channels, dims...)
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

effective_gate_time(::NullCountingGate, exposure_time) = exposure_time
effective_gate_time(model::DutyCycleGate, exposure_time) = exposure_time * model.duty_cycle
counting_exposure_time(det::APDDetector) = effective_gate_time(det.params.gate_model, det.params.integration_time)

function apply_gain_map!(det::APDDetector)
    isnothing(det.channel_gain_map) && return det.state.channels
    size(det.channel_gain_map) == size(det.state.channels) ||
        throw(DimensionMismatchError("APDDetector channel_gain_map size must match counting-channel size"))
    det.state.channels .*= det.channel_gain_map
    return det.state.channels
end

function apply_dark_counts!(det::APDDetector, exposure_time::Real)
    dark = effective_dark_count_rate(det) * exposure_time
    dark <= 0 && return det.state.channels
    det.state.channels .+= dark
    return det.state.channels
end

apply_dead_time!(det::APDDetector) = apply_dead_time!(det.params.dead_time_model, det)
apply_dead_time!(::NoDeadTime, det::APDDetector) = det.state.channels

function apply_dead_time!(model::NonParalyzableDeadTime, det::APDDetector)
    exposure_time = counting_exposure_time(det)
    exposure_time > zero(exposure_time) || return det.state.channels
    scale = model.dead_time / exposure_time
    scale <= zero(scale) && return det.state.channels
    @. det.state.channels = det.state.channels / (1 + det.state.channels * scale)
    return det.state.channels
end

function apply_dead_time!(model::ParalyzableDeadTime, det::APDDetector)
    exposure_time = counting_exposure_time(det)
    exposure_time > zero(exposure_time) || return det.state.channels
    scale = model.dead_time / exposure_time
    scale <= zero(scale) && return det.state.channels
    @. det.state.channels = det.state.channels * exp(-det.state.channels * scale)
    return det.state.channels
end

apply_counting_noise!(det::APDDetector{NoiseNone}, rng::AbstractRNG) = det.state.channels

function apply_counting_noise!(det::APDDetector{NoisePhoton}, rng::AbstractRNG)
    poisson_noise!(rng, det.state.channels)
    return det.state.channels
end

function apply_gain!(det::APDDetector)
    det.state.channels .*= det.params.gain
    return det.state.channels
end

apply_counting_correlation!(::NullCountingCorrelation, det::APDDetector, rng::AbstractRNG) = det.state.channels

function apply_counting_correlation!(model::AfterpulsingModel, det::APDDetector, rng::AbstractRNG)
    p = model.probability
    p <= zero(p) && return det.state.channels
    det.state.channels .*= (one(p) + p)
    return det.state.channels
end

function apply_counting_correlation!(model::ChannelCrosstalkModel, det::APDDetector, rng::AbstractRNG)
    coupling = model.coupling
    coupling <= zero(coupling) && return det.state.channels
    copyto!(det.state.noise_buffer, det.state.channels)
    fill!(det.state.channels, zero(eltype(det.state.channels)))
    n, m = size(det.state.channels)
    @inbounds for i in 1:n, j in 1:m
        center = det.state.noise_buffer[i, j]
        bleed = coupling * center
        keep = center - bleed
        det.state.channels[i, j] += keep
        neighbors = 0
        i > 1 && (neighbors += 1)
        i < n && (neighbors += 1)
        j > 1 && (neighbors += 1)
        j < m && (neighbors += 1)
        neighbors == 0 && continue
        share = bleed / neighbors
        i > 1 && (det.state.channels[i - 1, j] += share)
        i < n && (det.state.channels[i + 1, j] += share)
        j > 1 && (det.state.channels[i, j - 1] += share)
        j < m && (det.state.channels[i, j + 1] += share)
    end
    return det.state.channels
end

function apply_counting_correlation!(model::CompositeCountingCorrelation, det::APDDetector, rng::AbstractRNG)
    foreach(stage -> apply_counting_correlation!(stage, det, rng), model.stages)
    return det.state.channels
end

function write_output!(det::APDDetector)
    output = det.state.output_buffer
    output === nothing && return det.state.channels
    output_eltype = eltype(output)
    if output_eltype <: Integer
        lo = typemin(output_eltype)
        hi = typemax(output_eltype)
        @. output = output_eltype(clamp(round(det.state.channels), lo, hi))
    else
        copyto!(output, det.state.channels)
    end
    return output
end

function capture!(det::APDDetector, channels::AbstractMatrix{T}; rng::AbstractRNG=Random.default_rng()) where {T<:AbstractFloat}
    ensure_buffers!(det, size(channels))
    copyto!(det.state.channels, channels)
    exposure_time = effective_gate_time(det.params.gate_model, det.params.integration_time)
    det.state.channels .*= det.params.qe * exposure_time
    apply_gain_map!(det)
    apply_dark_counts!(det, exposure_time)
    apply_counting_noise!(det, rng)
    apply_dead_time!(det)
    apply_counting_correlation!(det.params.correlation_model, det, rng)
    apply_gain!(det)
    return write_output!(det)
end

capture!(det::APDDetector, channels::AbstractMatrix{T}, rng::AbstractRNG) where {T<:AbstractFloat} =
    capture!(det, channels; rng=rng)

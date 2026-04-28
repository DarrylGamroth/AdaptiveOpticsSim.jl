struct NoDeadTime <: CountingDeadTimeModel end

struct NonParalyzableDeadTime{T<:AbstractFloat} <: CountingDeadTimeModel
    dead_time::T
end

struct ParalyzableDeadTime{T<:AbstractFloat} <: CountingDeadTimeModel
    dead_time::T
end

NonParalyzableDeadTime(dead_time::Real) = NonParalyzableDeadTime{Float64}(float(dead_time))
ParalyzableDeadTime(dead_time::Real) = ParalyzableDeadTime{Float64}(float(dead_time))

readout_ready(det::AbstractCountingDetector) = true
thermal_model(det::AbstractCountingDetector) = det.params.thermal_model
thermal_state(det::AbstractCountingDetector) = det.state.thermal_state
reset_integration!(det::AbstractCountingDetector) = det

detector_temperature(det::AbstractCountingDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} =
    detector_temperature_K(thermal_model(det), thermal_state(det), T)

advance_thermal!(det::AbstractCountingDetector, dt) =
    (advance_thermal!(thermal_model(det), thermal_state(det), dt); det)

channel_output(det::AbstractCountingDetector) =
    counting_output_buffer(det) === nothing ? counting_array(det) : counting_output_buffer(det)
output_frame(det::AbstractCountingDetector) = channel_output(det)

counting_sensor(det::AbstractCountingDetector) =
    throw(InvalidConfiguration("missing counting_sensor overload for $(typeof(det))"))
counting_gate_model(det::AbstractCountingDetector) =
    throw(InvalidConfiguration("missing counting_gate_model overload for $(typeof(det))"))
counting_dead_time_model(det::AbstractCountingDetector) =
    throw(InvalidConfiguration("missing counting_dead_time_model overload for $(typeof(det))"))
counting_correlation_model(det::AbstractCountingDetector) =
    throw(InvalidConfiguration("missing counting_correlation_model overload for $(typeof(det))"))
counting_integration_time(det::AbstractCountingDetector) =
    throw(InvalidConfiguration("missing counting_integration_time overload for $(typeof(det))"))
counting_layout(det::AbstractCountingDetector) =
    throw(InvalidConfiguration("missing counting_layout overload for $(typeof(det))"))
counting_output_type(det::AbstractCountingDetector) =
    throw(InvalidConfiguration("missing counting_output_type overload for $(typeof(det))"))
detector_output_type(det::AbstractCountingDetector) = counting_output_type(det)
counting_array(det::AbstractCountingDetector) =
    throw(InvalidConfiguration("missing counting_array overload for $(typeof(det))"))
counting_noise_buffer(det::AbstractCountingDetector) =
    throw(InvalidConfiguration("missing counting_noise_buffer overload for $(typeof(det))"))
counting_output_buffer(det::AbstractCountingDetector) =
    throw(InvalidConfiguration("missing counting_output_buffer overload for $(typeof(det))"))
set_counting_array!(det::AbstractCountingDetector, values) =
    throw(InvalidConfiguration("missing set_counting_array! overload for $(typeof(det))"))
set_counting_noise_buffer!(det::AbstractCountingDetector, values) =
    throw(InvalidConfiguration("missing set_counting_noise_buffer! overload for $(typeof(det))"))
set_counting_output_buffer!(det::AbstractCountingDetector, values) =
    throw(InvalidConfiguration("missing set_counting_output_buffer! overload for $(typeof(det))"))

counting_qe(det::AbstractCountingDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} = one(T)
counting_fill_factor(det::AbstractCountingDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} = one(T)
counting_reported_fill_factor(det::AbstractCountingDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} = nothing
counting_post_gain(det::AbstractCountingDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} = one(T)
counting_dark_count_rate(det::AbstractCountingDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} = zero(T)
counting_channel_gain_map(det::AbstractCountingDetector) = nothing

supports_dead_time(::NoDeadTime) = false
supports_dead_time(::NonParalyzableDeadTime) = true
supports_dead_time(::ParalyzableDeadTime) = true
_supports_afterpulsing(::AbstractCountingCorrelationModel) = false
_supports_afterpulsing(::AfterpulsingModel) = true
_supports_afterpulsing(model::CompositeCountingCorrelation) = any(_supports_afterpulsing, model.stages)
_supports_channel_crosstalk(::AbstractCountingCorrelationModel) = false
_supports_channel_crosstalk(::ChannelCrosstalkModel) = true
_supports_channel_crosstalk(model::CompositeCountingCorrelation) = any(_supports_channel_crosstalk, model.stages)

is_null_counting_gate(::AbstractCountingGateModel) = false
is_null_counting_gate(::NullCountingGate) = true
is_paralyzable_dead_time(::CountingDeadTimeModel) = false
is_paralyzable_dead_time(::ParalyzableDeadTime) = true

counting_dead_time_symbol(::NoDeadTime) = :none
counting_dead_time_symbol(::NonParalyzableDeadTime) = :nonparalyzable
counting_dead_time_symbol(::ParalyzableDeadTime) = :paralyzable
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

validate_counting_noise(noise::NoiseNone) = noise
validate_counting_noise(noise::NoisePhoton) = noise
validate_counting_noise(::NoiseReadout) =
    throw(InvalidConfiguration("counting detectors do not support additive readout noise; use NoiseNone or NoisePhoton"))
validate_counting_noise(::NoisePhotonReadout) =
    throw(InvalidConfiguration("counting detectors do not support frame-style readout noise; use NoisePhoton"))

dark_count_law(::AbstractCountingDetector) = NullTemperatureLaw()
active_dark_count_law(det::AbstractCountingDetector, ::NullDetectorThermalModel) = dark_count_law(det)

function active_dark_count_law(det::AbstractCountingDetector, model::FixedTemperature)
    return is_null_temperature_law(model.dark_count_law) ? dark_count_law(det) : model.dark_count_law
end

function active_dark_count_law(det::AbstractCountingDetector, model::FirstOrderThermalModel)
    return is_null_temperature_law(model.dark_count_law) ? dark_count_law(det) : model.dark_count_law
end

effective_dark_count_rate(det::AbstractCountingDetector, ::Type{T}=eltype(counting_array(det))) where {T<:AbstractFloat} =
    T(evaluate_temperature_law(active_dark_count_law(det, thermal_model(det)),
        counting_dark_count_rate(det, T), detector_temperature(det, T)))

function detector_export_metadata(det::AbstractCountingDetector; T::Type{<:AbstractFloat}=eltype(counting_array(det)))
    output = output_frame(det)
    return CountingDetectorExportMetadata{T}(
        T(counting_integration_time(det)),
        counting_qe(det, T),
        counting_reported_fill_factor(det, T),
        counting_post_gain(det, T),
        counting_dark_count_rate(det, T),
        counting_dead_time_symbol(counting_dead_time_model(det)),
        counting_dead_time_value(counting_dead_time_model(det), T),
        counting_gate_symbol(counting_gate_model(det)),
        counting_gate_duty_cycle(counting_gate_model(det), T),
        counting_correlation_symbol(counting_correlation_model(det)),
        afterpulse_probability(counting_correlation_model(det), T),
        crosstalk_value(counting_correlation_model(det), T),
        thermal_model_symbol(thermal_model(det)),
        detector_temperature(det, T),
        ambient_temperature_K(thermal_model(det), T),
        cooling_setpoint_K(thermal_model(det), T),
        thermal_time_constant_s(thermal_model(det), T),
        temperature_law_symbol(active_dark_count_law(det, thermal_model(det))),
        detector_sensor_symbol(counting_sensor(det)),
        detector_noise_symbol(det.noise),
        counting_output_type(det),
        CountingReadoutMetadata(counting_layout(det), size(output), length(output)),
    )
end

function ensure_buffers!(det::AbstractCountingDetector, dims::Tuple{Int,Int})
    if size(counting_array(det)) != dims
        set_counting_array!(det, similar(counting_array(det), dims...))
    end
    if size(counting_noise_buffer(det)) != dims
        set_counting_noise_buffer!(det, similar(counting_noise_buffer(det), dims...))
    end
    output = counting_output_buffer(det)
    if output !== nothing && size(output) != dims
        new_output = similar(output, dims...)
        fill!(new_output, zero(eltype(new_output)))
        set_counting_output_buffer!(det, new_output)
    end
    return det
end

effective_gate_time(::NullCountingGate, exposure_time) = exposure_time
effective_gate_time(model::DutyCycleGate, exposure_time) = exposure_time * model.duty_cycle
counting_exposure_time(det::AbstractCountingDetector) = effective_gate_time(counting_gate_model(det), counting_integration_time(det))

function seed_counting_input!(det::AbstractCountingDetector, input::AbstractMatrix)
    copyto!(counting_array(det), input)
    counting_array(det) .*= counting_qe(det) * counting_fill_factor(det) * counting_exposure_time(det)
    return counting_array(det)
end

apply_counting_channel_gain_map!(det::AbstractCountingDetector) = counting_array(det)

function apply_dark_counts!(det::AbstractCountingDetector, exposure_time::Real)
    dark = effective_dark_count_rate(det) * exposure_time
    dark <= 0 && return counting_array(det)
    counting_array(det) .+= dark
    return counting_array(det)
end

apply_dead_time!(det::AbstractCountingDetector) = apply_dead_time!(counting_dead_time_model(det), det)
apply_dead_time!(::NoDeadTime, det::AbstractCountingDetector) = counting_array(det)

function apply_dead_time!(model::NonParalyzableDeadTime, det::AbstractCountingDetector)
    exposure_time = counting_exposure_time(det)
    exposure_time > zero(exposure_time) || return counting_array(det)
    scale = model.dead_time / exposure_time
    scale <= zero(scale) && return counting_array(det)
    counts = counting_array(det)
    @. counts = counts / (1 + counts * scale)
    return counts
end

function apply_dead_time!(model::ParalyzableDeadTime, det::AbstractCountingDetector)
    exposure_time = counting_exposure_time(det)
    exposure_time > zero(exposure_time) || return counting_array(det)
    scale = model.dead_time / exposure_time
    scale <= zero(scale) && return counting_array(det)
    counts = counting_array(det)
    @. counts = counts * exp(-counts * scale)
    return counts
end

apply_counting_noise!(det::AbstractCountingDetector, rng::AbstractRNG) = apply_counting_noise!(det.noise, det, rng)
apply_counting_noise!(::NoiseNone, det::AbstractCountingDetector, rng::AbstractRNG) = counting_array(det)

function apply_counting_noise!(::NoisePhoton, det::AbstractCountingDetector, rng::AbstractRNG)
    poisson_noise!(rng, counting_array(det))
    return counting_array(det)
end

apply_counting_correlation!(det::AbstractCountingDetector, rng::AbstractRNG) =
    apply_counting_correlation!(counting_correlation_model(det), det, rng)
apply_counting_correlation!(::NullCountingCorrelation, det::AbstractCountingDetector, rng::AbstractRNG) = counting_array(det)

function apply_counting_correlation!(model::AfterpulsingModel, det::AbstractCountingDetector, rng::AbstractRNG)
    p = model.probability
    p <= zero(p) && return counting_array(det)
    counting_array(det) .*= (one(p) + p)
    return counting_array(det)
end

function apply_counting_correlation!(model::ChannelCrosstalkModel, det::AbstractCountingDetector, rng::AbstractRNG)
    coupling = model.coupling
    coupling <= zero(coupling) && return counting_array(det)
    copyto!(counting_noise_buffer(det), counting_array(det))
    fill!(counting_array(det), zero(eltype(counting_array(det))))
    n, m = size(counting_array(det))
    @inbounds for i in 1:n, j in 1:m
        center = counting_noise_buffer(det)[i, j]
        bleed = coupling * center
        keep = center - bleed
        counting_array(det)[i, j] += keep
        neighbors = 0
        i > 1 && (neighbors += 1)
        i < n && (neighbors += 1)
        j > 1 && (neighbors += 1)
        j < m && (neighbors += 1)
        neighbors == 0 && continue
        share = bleed / neighbors
        i > 1 && (counting_array(det)[i - 1, j] += share)
        i < n && (counting_array(det)[i + 1, j] += share)
        j > 1 && (counting_array(det)[i, j - 1] += share)
        j < m && (counting_array(det)[i, j + 1] += share)
    end
    return counting_array(det)
end

function apply_counting_correlation!(model::CompositeCountingCorrelation, det::AbstractCountingDetector, rng::AbstractRNG)
    foreach(stage -> apply_counting_correlation!(stage, det, rng), model.stages)
    return counting_array(det)
end

apply_post_counting_gain!(det::AbstractCountingDetector) = counting_array(det)

function write_output!(det::AbstractCountingDetector)
    output = counting_output_buffer(det)
    output === nothing && return counting_array(det)
    output_eltype = eltype(output)
    if output_eltype <: Integer
        lo = typemin(output_eltype)
        hi = typemax(output_eltype)
        @. output = output_eltype(clamp(round(counting_array(det)), lo, hi))
    else
        copyto!(output, counting_array(det))
    end
    return output
end

function capture!(det::AbstractCountingDetector, channels::AbstractMatrix{T}; rng::AbstractRNG=Random.default_rng()) where {T<:AbstractFloat}
    ensure_buffers!(det, size(channels))
    exposure_time = counting_exposure_time(det)
    seed_counting_input!(det, channels)
    apply_counting_channel_gain_map!(det)
    apply_dark_counts!(det, exposure_time)
    apply_counting_noise!(det, rng)
    apply_dead_time!(det)
    apply_counting_correlation!(det, rng)
    apply_post_counting_gain!(det)
    advance_thermal!(det, counting_integration_time(det))
    return write_output!(det)
end

capture!(det::AbstractCountingDetector, channels::AbstractMatrix{T}, rng::AbstractRNG) where {T<:AbstractFloat} =
    capture!(det, channels; rng=rng)

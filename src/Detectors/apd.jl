struct APDSensor <: CountingSensorType end
struct NoDeadTime <: CountingDeadTimeModel end

struct NonParalyzableDeadTime{T<:AbstractFloat} <: CountingDeadTimeModel
    dead_time::T
end

struct APDDetectorParams{T<:AbstractFloat,S<:APDSensor,D<:CountingDeadTimeModel}
    integration_time::T
    qe::T
    gain::T
    dark_count_rate::T
    dead_time_model::D
    sensor::S
    output_precision::Union{Nothing,DataType}
    layout::Symbol
end

mutable struct APDDetectorState{T<:AbstractFloat,A<:AbstractMatrix{T},O}
    channels::A
    noise_buffer::A
    output_buffer::O
end

struct APDDetector{N<:NoiseModel,P<:APDDetectorParams,S<:APDDetectorState,GM} <: AbstractCountingDetector
    noise::N
    params::P
    state::S
    channel_gain_map::GM
end

NonParalyzableDeadTime(dead_time::Real) = NonParalyzableDeadTime{Float64}(float(dead_time))

readout_ready(det::APDDetector) = true
channel_output(det::APDDetector) = det.state.output_buffer === nothing ? det.state.channels : det.state.output_buffer
output_frame(det::APDDetector) = channel_output(det)
reset_integration!(det::APDDetector) = det

detector_sensor_symbol(::APDSensor) = :apd
counting_dead_time_symbol(::NoDeadTime) = :none
counting_dead_time_symbol(::NonParalyzableDeadTime) = :nonparalyzable
supports_counting_noise(::APDDetector) = true
supports_dead_time(det::APDDetector) = supports_dead_time(det.params.dead_time_model)
supports_dead_time(::NoDeadTime) = false
supports_dead_time(::NonParalyzableDeadTime) = true
supports_channel_gain_map(det::APDDetector) = !isnothing(det.channel_gain_map)

counting_dead_time_value(::NoDeadTime, ::Type{T}) where {T<:AbstractFloat} = nothing
counting_dead_time_value(model::NonParalyzableDeadTime, ::Type{T}) where {T<:AbstractFloat} = T(model.dead_time)
convert_dead_time_model(::NoDeadTime, ::Type{T}) where {T<:AbstractFloat} = NoDeadTime()
convert_dead_time_model(model::NonParalyzableDeadTime, ::Type{T}) where {T<:AbstractFloat} =
    NonParalyzableDeadTime{T}(T(model.dead_time))
validate_dead_time_model(model::NoDeadTime) = model

function validate_dead_time_model(model::NonParalyzableDeadTime)
    model.dead_time >= 0 || throw(InvalidConfiguration("NonParalyzableDeadTime dead_time must be >= 0"))
    return model
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
        detector_sensor_symbol(det.params.sensor),
        detector_noise_symbol(det.noise),
        det.params.output_precision,
        CountingReadoutMetadata(det.params.layout, size(output), length(output)),
    )
end

function _build_apd_detector(noise::NoiseModel; integration_time::Real, qe::Real, gain::Real,
    dark_count_rate::Real, dead_time_model::CountingDeadTimeModel, sensor::APDSensor, output_precision::Union{Nothing,DataType},
    layout::Symbol, channel_gain_map,
    T::Type{<:AbstractFloat}, backend)
    gain >= 0 || throw(InvalidConfiguration("APDDetector gain must be >= 0"))
    dark_count_rate >= 0 || throw(InvalidConfiguration("APDDetector dark_count_rate must be >= 0"))
    converted = convert_noise(noise, T)
    validated = validate_apd_noise(converted)
    dead_time = validate_dead_time_model(convert_dead_time_model(dead_time_model, T))
    gain_map = channel_gain_map === nothing ? nothing : begin
        g = backend{T}(undef, size(channel_gain_map)...)
        copyto!(g, T.(channel_gain_map))
        g
    end
    params = APDDetectorParams{T,typeof(sensor),typeof(dead_time)}(
        T(integration_time),
        T(qe),
        T(gain),
        T(dark_count_rate),
        dead_time,
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
    state = APDDetectorState{T,typeof(channels),typeof(output_buffer)}(channels, noise_buffer, output_buffer)
    return APDDetector{typeof(validated),typeof(params),typeof(state),typeof(gain_map)}(
        validated, params, state, gain_map)
end

function APDDetector(; integration_time::Real=1.0, qe::Real=1.0, noise::NoiseModel=NoisePhoton(),
    gain::Real=1.0, dark_count_rate::Real=0.0, output_precision::Union{Nothing,DataType}=nothing,
    layout::Symbol=:channels, channel_gain_map=nothing, dead_time_model::CountingDeadTimeModel=NoDeadTime(),
    T::Type{<:AbstractFloat}=Float64, backend=Array)
    return _build_apd_detector(noise; integration_time=integration_time, qe=qe, gain=gain,
        dark_count_rate=dark_count_rate, dead_time_model=dead_time_model, sensor=APDSensor(), output_precision=output_precision,
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

function apply_gain_map!(det::APDDetector)
    isnothing(det.channel_gain_map) && return det.state.channels
    size(det.channel_gain_map) == size(det.state.channels) ||
        throw(DimensionMismatchError("APDDetector channel_gain_map size must match counting-channel size"))
    det.state.channels .*= det.channel_gain_map
    return det.state.channels
end

function apply_dark_counts!(det::APDDetector, exposure_time::Real)
    dark = det.params.dark_count_rate * exposure_time
    dark <= 0 && return det.state.channels
    det.state.channels .+= dark
    return det.state.channels
end

apply_dead_time!(det::APDDetector) = apply_dead_time!(det.params.dead_time_model, det)
apply_dead_time!(::NoDeadTime, det::APDDetector) = det.state.channels

function apply_dead_time!(model::NonParalyzableDeadTime, det::APDDetector)
    exposure_time = det.params.integration_time
    exposure_time > zero(exposure_time) || return det.state.channels
    scale = model.dead_time / exposure_time
    scale <= zero(scale) && return det.state.channels
    @. det.state.channels = det.state.channels / (1 + det.state.channels * scale)
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
    det.state.channels .*= det.params.qe * det.params.integration_time
    apply_gain_map!(det)
    apply_dark_counts!(det, det.params.integration_time)
    apply_counting_noise!(det, rng)
    apply_dead_time!(det)
    apply_gain!(det)
    return write_output!(det)
end

capture!(det::APDDetector, channels::AbstractMatrix{T}, rng::AbstractRNG) where {T<:AbstractFloat} =
    capture!(det, channels; rng=rng)

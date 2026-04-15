using AdaptiveOpticsSim
using Statistics: mean

include(joinpath(dirname(@__DIR__), "benchmarks", "support", "revolt_like_hil_common.jl"))

const _backend_arg = isempty(ARGS) ? "cuda" : lowercase(ARGS[1])
const _default_config_dir = joinpath(dirname(@__DIR__), "benchmarks", "assets", "revolt_like")
const _config_dir_arg = length(ARGS) >= 2 ? ARGS[2] : get(ENV, "REVOLT_CONFIG_DIR", _default_config_dir)
const _sensor_arg = length(ARGS) >= 3 ? lowercase(ARGS[3]) : "cmos"
const _response_arg = length(ARGS) >= 4 ? lowercase(ARGS[4]) : "default"
const _thermal_arg = length(ARGS) >= 5 ? lowercase(ARGS[5]) : "none"
const _steps_arg = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 100_000
const _warmup_arg = length(ARGS) >= 7 ? parse(Int, ARGS[7]) : 1_000

if _backend_arg == "cuda"
    import CUDA
end

if _backend_arg == "amdgpu"
    import AMDGPU
end

function _resolve_sensor(name::AbstractString)
    lowered = lowercase(name)
    lowered == "ccd" && return CCDSensor(), :ccd
    lowered == "cmos" && return CMOSSensor(), :cmos
    lowered == "emccd" && return EMCCDSensor(), :emccd
    lowered == "ingaas" && return InGaAsSensor(), :ingaas
    lowered == "hgcdte" && return HgCdTeAvalancheArraySensor(), :hgcdte_avalanche_array
    error("unsupported sensor '$name'; use ccd, cmos, emccd, ingaas, or hgcdte")
end

function _resolve_response_model(name::AbstractString)
    lowered = lowercase(name)
    lowered == "default" && return nothing, :default
    lowered == "null" && return NullFrameResponse(), :null
    error("unsupported response mode '$name'; use default or null")
end

function _resolve_thermal_model(name::AbstractString, sensor_label::Symbol, T::Type{<:AbstractFloat})
    lowered = lowercase(name)
    lowered == "none" && return nothing, :none
    dark_law = ArrheniusRateLaw(293.0, 4500.0)
    glow_law = ArrheniusRateLaw(293.0, 3000.0)
    cic_law = LinearTemperatureLaw(293.0, 0.01)
    if lowered == "fixed120"
        return FixedTemperature(
            temperature_K=120.0,
            dark_current_law=dark_law,
            glow_rate_law=(sensor_label in (:ingaas, :hgcdte_avalanche_array) ? glow_law : NullTemperatureLaw()),
            cic_rate_law=(sensor_label == :emccd ? cic_law : NullTemperatureLaw()),
            T=T), :fixed120
    elseif lowered == "dynamic120"
        return FirstOrderThermalModel(
            ambient_temperature_K=293.0,
            setpoint_temperature_K=120.0,
            initial_temperature_K=293.0,
            time_constant_s=0.05,
            min_temperature_K=100.0,
            max_temperature_K=300.0,
            dark_current_law=dark_law,
            glow_rate_law=(sensor_label in (:ingaas, :hgcdte_avalanche_array) ? glow_law : NullTemperatureLaw()),
            cic_rate_law=(sensor_label == :emccd ? cic_law : NullTemperatureLaw()),
            T=T), :dynamic120
    end
    error("unsupported thermal mode '$name'; use none, fixed120, or dynamic120")
end

_thermalized_sensor(sensor::CCDSensor, ::Bool, ::Type{T}) where {T<:AbstractFloat} = sensor
_thermalized_sensor(sensor::CMOSSensor, ::Bool, ::Type{T}) where {T<:AbstractFloat} = sensor
_thermalized_sensor(sensor::EMCCDSensor, enabled::Bool, ::Type{T}) where {T<:AbstractFloat} =
    enabled ? EMCCDSensor(cic_rate=T(0.02), excess_noise_factor=sensor.excess_noise_factor,
        em_gain_model=sensor.em_gain_model, T=T) : sensor
_thermalized_sensor(sensor::InGaAsSensor, enabled::Bool, ::Type{T}) where {T<:AbstractFloat} =
    enabled ? InGaAsSensor(glow_rate=T(0.02), persistence_model=sensor.persistence_model, T=T) : sensor
_thermalized_sensor(sensor::HgCdTeAvalancheArraySensor, enabled::Bool, ::Type{T}) where {T<:AbstractFloat} =
    enabled ? HgCdTeAvalancheArraySensor(
        avalanche_gain=sensor.avalanche_gain,
        excess_noise_factor=sensor.excess_noise_factor,
        glow_rate=T(0.02),
        read_time=sensor.read_time,
        sampling_mode=sensor.sampling_mode,
        T=T) : sensor

function _gc_counter_delta(before, after)
    names = fieldnames(typeof(before))
    pairs = Tuple(name => (getfield(after, name) - getfield(before, name)) for name in names)
    return NamedTuple(pairs)
end

function _progress_interval(steps::Int)
    return max(1, fld(steps, 10))
end

function _timing_summary(timings::Vector{Int})
    sorted = sort(copy(timings))
    samples = length(timings)
    p50_idx = clamp(round(Int, 0.50 * samples), 1, samples)
    p95_idx = clamp(round(Int, 0.95 * samples), 1, samples)
    p99_idx = clamp(round(Int, 0.99 * samples), 1, samples)
    return (
        min_ns=sorted[1],
        max_ns=sorted[end],
        mean_ns=mean(timings),
        p50_ns=sorted[p50_idx],
        p95_ns=sorted[p95_idx],
        p99_ns=sorted[p99_idx],
    )
end

function _steady_step_alloc_bytes(ctx)
    for _ in 1:10
        revolt_like_step!(ctx)
    end
    GC.gc()
    return @allocated revolt_like_step!(ctx)
end

function run_soak(; backend_name::AbstractString="cuda", config_dir::AbstractString=_config_dir_arg,
    sensor_name::AbstractString="cmos", response_name::AbstractString="default",
    thermal_name::AbstractString="none", steps::Int=100_000, warmup::Int=1_000)
    steps > 0 || error("steps must be positive")
    warmup >= 0 || error("warmup must be >= 0")

    sensor, sensor_label = _resolve_sensor(sensor_name)
    response_model, response_label = _resolve_response_model(response_name)
    T = Float32
    thermal_model, thermal_label = _resolve_thermal_model(thermal_name, sensor_label, T)
    sensor = _thermalized_sensor(sensor, !isnothing(thermal_model), T)

    ctx = build_revolt_like_hil_context(; backend_name, config_dir, sensor, response_model, thermal_model, T)
    metadata = detector_export_metadata(ctx.det)
    backend_cfg = revolt_profile_backend(backend_name)

    for _ in 1:warmup
        revolt_like_step!(ctx)
    end

    steady_alloc_bytes = _steady_step_alloc_bytes(ctx)
    timings = Vector{Int}(undef, steps)
    interval = _progress_interval(steps)
    GC.gc()
    gc_before = Base.gc_num()
    t_run_start = time_ns()
    @inbounds for i in 1:steps
        t0 = time_ns()
        revolt_like_step!(ctx)
        timings[i] = time_ns() - t0
        if i % interval == 0 || i == steps
            println("soak_progress step=", i, "/", steps)
        end
    end
    total_wall_ns = time_ns() - t_run_start
    gc_after = Base.gc_num()

    stats = _timing_summary(timings)
    gc_delta = _gc_counter_delta(gc_before, gc_after)

    println("revolt_like_hil_soak")
    println("  backend: ", backend_cfg.label)
    println("  config_dir: ", config_dir)
    println("  sensor: ", sensor_label)
    println("  response_mode: ", response_label)
    println("  thermal_mode: ", thermal_label)
    println("  effective_frame_response: ", metadata.frame_response)
    println("  thermal_model: ", metadata.thermal_model)
    println("  detector_temperature_K: ", metadata.detector_temperature_K)
    println("  steps: ", steps)
    println("  warmup: ", warmup)
    println("  steady_step_alloc_bytes: ", steady_alloc_bytes)
    println("  total_wall_ns: ", total_wall_ns)
    println("  frame_rate_hz: ", 1.0e9 * steps / total_wall_ns)
    println("  min_ns: ", stats.min_ns)
    println("  mean_ns: ", stats.mean_ns)
    println("  p50_ns: ", stats.p50_ns)
    println("  p95_ns: ", stats.p95_ns)
    println("  p99_ns: ", stats.p99_ns)
    println("  max_ns: ", stats.max_ns)
    println("  gc_collect_delta: ", gc_delta.collect)
    println("  gc_allocd_delta: ", gc_delta.allocd)
    println("  gc_total_allocd_delta: ", gc_delta.total_allocd)
    println("  gc_malloc_delta: ", gc_delta.malloc)
    println("  gc_realloc_delta: ", gc_delta.realloc)
    println("  gc_bigalloc_delta: ", gc_delta.bigalloc)
    println("  gc_total_time_delta: ", gc_delta.total_time)
    println("  gc_pause_delta: ", gc_delta.pause)
    println("  gc_full_sweep_delta: ", gc_delta.full_sweep)
    println("  gc_max_pause_delta: ", gc_delta.max_pause)
    return nothing
end

run_soak(; backend_name=_backend_arg, config_dir=_config_dir_arg, sensor_name=_sensor_arg,
    response_name=_response_arg, thermal_name=_thermal_arg, steps=_steps_arg, warmup=_warmup_arg)

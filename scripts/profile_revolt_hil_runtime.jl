using AdaptiveOpticsSim
using Statistics: mean

include(joinpath(dirname(@__DIR__), "benchmarks", "support", "revolt_like_hil_common.jl"))

const _backend_arg = isempty(ARGS) ? "cpu" : lowercase(ARGS[1])
const _default_config_dir = joinpath(dirname(@__DIR__), "benchmarks", "assets", "revolt_like")
const _config_dir_arg = length(ARGS) >= 2 ? ARGS[2] : get(ENV, "REVOLT_CONFIG_DIR", _default_config_dir)
const _sensor_arg = length(ARGS) >= 3 ? lowercase(ARGS[3]) : "ccd"
const _response_arg = length(ARGS) >= 4 ? lowercase(ARGS[4]) : "default"
const _thermal_arg = length(ARGS) >= 5 ? lowercase(ARGS[5]) : "none"
const _samples_arg = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 6
const _warmup_arg = length(ARGS) >= 7 ? parse(Int, ARGS[7]) : 2

function _timed_stats!(f!::F; warmup::Int=2, samples::Int=6) where {F<:Function}
    for _ in 1:warmup
        f!()
    end
    GC.gc()
    timings = Vector{Int}(undef, samples)
    @inbounds for i in 1:samples
        t0 = time_ns()
        f!()
        timings[i] = time_ns() - t0
    end
    sorted = sort(timings)
    p95_idx = clamp(round(Int, 0.95 * samples), 1, samples)
    return mean(timings), sorted[p95_idx]
end

function _allocated_bytes(f!::F; warmup::Int=2, gc_before::Bool=true) where {F<:Function}
    for _ in 1:warmup
        f!()
    end
    gc_before && GC.gc()
    return @allocated f!()
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

function run_profile(; backend_name::AbstractString="cpu", config_dir::AbstractString=_config_dir_arg,
    sensor_name::AbstractString="ccd", response_name::AbstractString="default",
    thermal_name::AbstractString="none", samples::Int=6, warmup::Int=2)
    backend_cfg = revolt_profile_backend(backend_name)
    sensor, sensor_label = _resolve_sensor(sensor_name)
    response_model, response_label = _resolve_response_model(response_name)
    T = Float32
    thermal_model, thermal_label = _resolve_thermal_model(thermal_name, sensor_label, T)
    sensor = _thermalized_sensor(sensor, !isnothing(thermal_model), T)
    ctx = build_revolt_like_hil_context(; backend_name, config_dir, sensor, response_model, thermal_model, T)
    metadata = detector_export_metadata(ctx.det)
    n_active = length(ctx.active_command)
    dm_grid_command_length = length(ctx.dm.state.coefs)
    build_time_ns = time_ns()
    ctx = build_revolt_like_hil_context(; backend_name, config_dir, sensor, response_model, thermal_model, T)
    build_time_ns = time_ns() - build_time_ns

    command_map_mean_ns, command_map_p95_ns = _timed_stats!(() -> revolt_like_command_map!(ctx); warmup=warmup, samples=samples)
    dm_apply_mean_ns, dm_apply_p95_ns = _timed_stats!(() -> revolt_like_dm_apply!(ctx); warmup=warmup, samples=samples)
    sense_mean_ns, sense_p95_ns = _timed_stats!(() -> revolt_like_sense!(ctx); warmup=warmup, samples=samples)
    mosaic_mean_ns, mosaic_p95_ns = _timed_stats!(() -> revolt_like_mosaic!(ctx); warmup=warmup, samples=samples)
    total_mean_ns, total_p95_ns = _timed_stats!(() -> revolt_like_step!(ctx); warmup=warmup, samples=samples)
    command_map_alloc_bytes = _allocated_bytes(() -> revolt_like_command_map!(ctx); warmup=warmup)
    dm_apply_alloc_bytes = _allocated_bytes(() -> revolt_like_dm_apply!(ctx); warmup=warmup)
    sense_alloc_bytes = _allocated_bytes(() -> revolt_like_sense!(ctx); warmup=warmup)
    mosaic_alloc_bytes = _allocated_bytes(() -> revolt_like_mosaic!(ctx); warmup=warmup)
    total_alloc_bytes = _allocated_bytes(() -> revolt_like_step!(ctx); warmup=warmup)

    println("revolt_like_hil_runtime_profile")
    println("  backend: ", backend_cfg.label)
    println("  config_dir: ", config_dir)
    println("  sensor: ", sensor_label)
    println("  response_mode: ", response_label)
    println("  thermal_mode: ", thermal_label)
    println("  effective_frame_response: ", metadata.frame_response)
    println("  thermal_model: ", metadata.thermal_model)
    println("  detector_temperature_K: ", metadata.detector_temperature_K)
    println("  cooling_setpoint_K: ", metadata.cooling_setpoint_K)
    println("  ambient_temperature_K: ", metadata.ambient_temperature_K)
    println("  actuator_command_length: ", n_active)
    println("  extrapolated_command_length: ", length(ctx.extrapolated_command))
    println("  dm_grid_command_length: ", dm_grid_command_length)
    println("  dm_layout_shape: ", size(ctx.dm.state.coefs))
    println("  pupil_resolution: ", size(ctx.tel.state.opd, 1))
    println("  n_lenslets: ", ctx.n_lenslets)
    println("  roi_size: ", ctx.roi)
    println("  tiled_wfs_frame_shape: ", size(ctx.tiled_frame))
    println("  spot_cube_shape: ", size(ctx.wfs.state.spot_cube))
    println("  build_time_ns: ", build_time_ns)
    println("  command_map_mean_ns: ", command_map_mean_ns)
    println("  command_map_p95_ns: ", command_map_p95_ns)
    println("  dm_apply_mean_ns: ", dm_apply_mean_ns)
    println("  dm_apply_p95_ns: ", dm_apply_p95_ns)
    println("  sense_mean_ns: ", sense_mean_ns)
    println("  sense_p95_ns: ", sense_p95_ns)
    println("  mosaic_mean_ns: ", mosaic_mean_ns)
    println("  mosaic_p95_ns: ", mosaic_p95_ns)
    println("  total_mean_ns: ", total_mean_ns)
    println("  total_p95_ns: ", total_p95_ns)
    println("  frame_rate_hz: ", 1.0e9 / total_mean_ns)
    println("  command_map_alloc_bytes: ", command_map_alloc_bytes)
    println("  dm_apply_alloc_bytes: ", dm_apply_alloc_bytes)
    println("  sense_alloc_bytes: ", sense_alloc_bytes)
    println("  mosaic_alloc_bytes: ", mosaic_alloc_bytes)
    println("  total_alloc_bytes: ", total_alloc_bytes)
    return nothing
end

run_profile(; backend_name=_backend_arg, config_dir=_config_dir_arg, sensor_name=_sensor_arg, response_name=_response_arg,
    thermal_name=_thermal_arg, samples=_samples_arg, warmup=_warmup_arg)

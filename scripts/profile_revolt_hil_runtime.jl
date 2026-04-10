using AdaptiveOpticsSim
using Random
using Statistics: mean
using DelimitedFiles
using SparseArrays
using LinearAlgebra: mul!
using KernelAbstractions: @kernel, @index

const _backend_arg = isempty(ARGS) ? "cpu" : lowercase(ARGS[1])
const _default_config_dir = joinpath(dirname(@__DIR__), "benchmarks", "assets", "revolt_like")
const _config_dir_arg = length(ARGS) >= 2 ? ARGS[2] : get(ENV, "REVOLT_CONFIG_DIR", _default_config_dir)
const _sensor_arg = length(ARGS) >= 3 ? lowercase(ARGS[3]) : "ccd"
const _response_arg = length(ARGS) >= 4 ? lowercase(ARGS[4]) : "default"
const _thermal_arg = length(ARGS) >= 5 ? lowercase(ARGS[5]) : "none"
const _samples_arg = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 6
const _warmup_arg = length(ARGS) >= 7 ? parse(Int, ARGS[7]) : 2

if _backend_arg == "cuda"
    import CUDA
end

if _backend_arg == "amdgpu"
    import AMDGPU
end

@kernel function scatter_active_command_kernel!(full_command, active_command, active_indices, n_active::Int)
    i = @index(Global, Linear)
    if i <= n_active
        @inbounds full_command[active_indices[i]] = active_command[i]
    end
end

@kernel function sh_spot_mosaic_kernel!(mosaic, spot_cube, n_subap::Int, roi::Int)
    idx, u, v = @index(Global, NTuple)
    if idx <= size(spot_cube, 1) && u <= roi && v <= roi
        sub_i = mod1(idx, n_subap)
        sub_j = fld(idx - 1, n_subap) + 1
        out_i = (sub_i - 1) * roi + u
        out_j = (sub_j - 1) * roi + v
        @inbounds mosaic[out_i, out_j] = spot_cube[idx, u, v]
    end
end

function _resolve_backend(name::AbstractString)
    lowered = lowercase(name)
    if lowered == "cpu"
        return Array, nothing, "cpu"
    elseif lowered == "cuda"
        isdefined(Main, :CUDA) || error("profile_revolt_hil_runtime.jl requires CUDA.jl for backend=cuda")
        CUDA.functional() || error("profile_revolt_hil_runtime.jl requires a functional CUDA driver/device")
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.CUDABackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(AdaptiveOpticsSim.CUDABackendTag)
        backend === nothing && error("CUDA backend array type is unavailable")
        return backend, AdaptiveOpticsSim.CUDABackendTag, "cuda"
    elseif lowered == "amdgpu"
        isdefined(Main, :AMDGPU) || error("profile_revolt_hil_runtime.jl requires AMDGPU.jl for backend=amdgpu")
        AMDGPU.functional() || error("profile_revolt_hil_runtime.jl requires a functional ROCm installation and GPU")
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.AMDGPUBackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(AdaptiveOpticsSim.AMDGPUBackendTag)
        backend === nothing && error("AMDGPU backend array type is unavailable")
        return backend, AdaptiveOpticsSim.AMDGPUBackendTag, "amdgpu"
    end
    error("unsupported backend '$name'; use cpu, cuda, or amdgpu")
end

_sync_backend!(::Nothing, _) = nothing

function _sync_backend!(::Type{B}, array) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(array))
    return nothing
end

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

function _load_dm277_actuator_map(path::AbstractString)
    raw = readlines(path)
    length(raw) >= 2 || error("actuator map file is unexpectedly short")
    header = split(strip(raw[1]))
    length(header) == 4 || error("unexpected actuator map header '$((raw[1]))'")
    n_row = parse(Int, header[3])
    n_col = parse(Int, header[4])
    grid = readdlm(IOBuffer(join(raw[2:end], '\n')), ',', Int)
    size(grid) == (n_row, n_col) || error("actuator map grid size $(size(grid)) does not match header $((n_row, n_col))")
    maximum(grid) > 0 || error("actuator map has no active actuators")
    active_indices = Vector{Int}(undef, maximum(grid))
    @inbounds for j in 1:n_col, i in 1:n_row
        label = grid[i, j]
        if label > 0
            active_indices[label] = LinearIndices(grid)[i, j]
        end
    end
    return grid, active_indices
end

function _load_dm_extrapolation(path::AbstractString, T::Type{<:AbstractFloat})
    raw = readlines(path)
    length(raw) >= 2 || error("extrapolation file is unexpectedly short")
    header = split(strip(raw[1]))
    length(header) == 4 || error("unexpected extrapolation header '$((raw[1]))'")
    n_row = parse(Int, header[3])
    n_col = parse(Int, header[4])
    triplets = readdlm(IOBuffer(join(raw[2:end], '\n')), ',', Float64)
    size(triplets, 2) == 3 || error("dmExtrapolation triplets must have three columns")
    rows = Int.(triplets[:, 1]) .+ 1
    cols = Int.(triplets[:, 2]) .+ 1
    vals = T.(triplets[:, 3])
    return sparse(rows, cols, vals, n_row, n_col)
end

function _fill_active_command!(active_command::AbstractVector{T}) where {T<:AbstractFloat}
    host = Vector{T}(undef, length(active_command))
    @inbounds for i in eachindex(host)
        host[i] = T(sin(0.013 * i) + 0.25 * cos(0.031 * i))
    end
    copyto!(active_command, host)
    return active_command
end

function _scatter_active_command!(full_command::AbstractVector{T}, active_command::AbstractVector{T},
    active_indices_backend::AbstractVector{Int}) where {T<:AbstractFloat}
    fill!(full_command, zero(T))
    style = AdaptiveOpticsSim.execution_style(full_command)
    if full_command isa Array
        @inbounds for i in eachindex(active_command)
            full_command[active_indices_backend[i]] = active_command[i]
        end
    else
        AdaptiveOpticsSim.launch_kernel!(style, scatter_active_command_kernel!,
            full_command, active_command, active_indices_backend, length(active_command);
            ndrange=length(active_command))
    end
    return full_command
end

function _tile_spot_cube!(mosaic::AbstractMatrix{T}, spot_cube::AbstractArray{T,3}, n_subap::Int, roi::Int) where {T<:AbstractFloat}
    style = AdaptiveOpticsSim.execution_style(mosaic)
    if mosaic isa Array
        @inbounds for idx in 1:size(spot_cube, 1)
            sub_i = mod1(idx, n_subap)
            sub_j = fld(idx - 1, n_subap) + 1
            out_i = (sub_i - 1) * roi + 1
            out_j = (sub_j - 1) * roi + 1
            copyto!(@view(mosaic[out_i:out_i + roi - 1, out_j:out_j + roi - 1]), @view(spot_cube[idx, :, :]))
        end
    else
        AdaptiveOpticsSim.launch_kernel!(style, sh_spot_mosaic_kernel!, mosaic, spot_cube, n_subap, roi;
            ndrange=size(spot_cube))
    end
    return mosaic
end

function run_profile(; backend_name::AbstractString="cpu", config_dir::AbstractString=_config_dir_arg,
    sensor_name::AbstractString="ccd", response_name::AbstractString="default",
    thermal_name::AbstractString="none",
    samples::Int=6, warmup::Int=2)
    BackendArray, backend_tag, label = _resolve_backend(backend_name)
    sensor, sensor_label = _resolve_sensor(sensor_name)
    response_model, response_label = _resolve_response_model(response_name)
    T = Float32
    thermal_model, thermal_label = _resolve_thermal_model(thermal_name, sensor_label, T)
    sensor = _thermalized_sensor(sensor, !isnothing(thermal_model), T)
    dark_current = isnothing(thermal_model) ? zero(T) : T(0.02)
    actuator_map_path = joinpath(config_dir, "revolt_like_dmActuatorMap_277.csv")
    extrapolation_path = joinpath(config_dir, "revolt_like_dmExtrapolation.csv")
    actuator_map, active_indices_host = _load_dm277_actuator_map(actuator_map_path)
    extrapolation_host = Matrix(_load_dm_extrapolation(extrapolation_path, T))
    n_subap = 16
    roi = 22
    resolution = 352
    n_act = size(actuator_map, 1)
    n_active = length(active_indices_host)

    tel = Telescope(
        resolution=resolution,
        diameter=8.0,
        sampling_time=1e-3,
        central_obstruction=0.0,
        T=T,
        backend=BackendArray,
    )
    src = Source(band=:I, magnitude=0.0, T=T)
    dm = DeformableMirror(tel; n_act=n_act, influence_width=0.3, T=T, backend=BackendArray)
    wfs = ShackHartmann(tel; n_subap=n_subap, mode=Diffractive(), n_pix_subap=roi, diffraction_padding=2, T=T, backend=BackendArray)
    det = Detector(noise=NoiseNone(), integration_time=T(1), qe=T(1), binning=1,
        dark_current=dark_current, sensor=sensor, response_model=response_model,
        thermal_model=thermal_model, T=T, backend=BackendArray)
    active_indices_backend = BackendArray{Int}(undef, n_active)
    copyto!(active_indices_backend, active_indices_host)
    active_command = BackendArray{T}(undef, n_active)
    extrapolated_command = BackendArray{T}(undef, n_active)
    extrapolation_backend = BackendArray{T}(undef, size(extrapolation_host)...)
    copyto!(extrapolation_backend, extrapolation_host)
    _fill_active_command!(active_command)
    tiled_frame = BackendArray{T}(undef, resolution, resolution)

    build_time_ns = time_ns()
    AdaptiveOpticsSim.ensure_sh_calibration!(wfs, tel, src)
    _sync_backend!(backend_tag, wfs.state.slopes)
    build_time_ns = time_ns() - build_time_ns

    function _step!()
        reset_opd!(tel)
        mul!(extrapolated_command, extrapolation_backend, active_command)
        _scatter_active_command!(dm.state.coefs, extrapolated_command, active_indices_backend)
        apply!(dm, tel, DMReplace())
        measure!(wfs, tel, src, det; rng=Random.default_rng())
        _tile_spot_cube!(tiled_frame, wfs.state.spot_cube, n_subap, roi)
        _sync_backend!(backend_tag, dm.state.coefs)
        _sync_backend!(backend_tag, wfs.state.spot_cube)
        _sync_backend!(backend_tag, tiled_frame)
        return nothing
    end

    command_map_mean_ns, command_map_p95_ns = _timed_stats!(() -> begin
        reset_opd!(tel)
        mul!(extrapolated_command, extrapolation_backend, active_command)
        _scatter_active_command!(dm.state.coefs, extrapolated_command, active_indices_backend)
        _sync_backend!(backend_tag, dm.state.coefs)
    end; warmup=warmup, samples=samples)
    dm_apply_mean_ns, dm_apply_p95_ns = _timed_stats!(() -> begin
        mul!(extrapolated_command, extrapolation_backend, active_command)
        _scatter_active_command!(dm.state.coefs, extrapolated_command, active_indices_backend)
        apply!(dm, tel, DMReplace())
        _sync_backend!(backend_tag, tel.state.opd)
    end; warmup=warmup, samples=samples)
    sense_mean_ns, sense_p95_ns = _timed_stats!(() -> begin
        mul!(extrapolated_command, extrapolation_backend, active_command)
        _scatter_active_command!(dm.state.coefs, extrapolated_command, active_indices_backend)
        apply!(dm, tel, DMReplace())
        measure!(wfs, tel, src, det; rng=Random.default_rng())
        _sync_backend!(backend_tag, wfs.state.spot_cube)
    end; warmup=warmup, samples=samples)
    mosaic_mean_ns, mosaic_p95_ns = _timed_stats!(() -> begin
        mul!(extrapolated_command, extrapolation_backend, active_command)
        _scatter_active_command!(dm.state.coefs, extrapolated_command, active_indices_backend)
        apply!(dm, tel, DMReplace())
        measure!(wfs, tel, src, det; rng=Random.default_rng())
        _tile_spot_cube!(tiled_frame, wfs.state.spot_cube, n_subap, roi)
        _sync_backend!(backend_tag, tiled_frame)
    end; warmup=warmup, samples=samples)
    total_mean_ns, total_p95_ns = _timed_stats!(_step!; warmup=warmup, samples=samples)
    command_map_alloc_bytes = _allocated_bytes(() -> begin
        reset_opd!(tel)
        mul!(extrapolated_command, extrapolation_backend, active_command)
        _scatter_active_command!(dm.state.coefs, extrapolated_command, active_indices_backend)
        _sync_backend!(backend_tag, dm.state.coefs)
    end; warmup=warmup)
    dm_apply_alloc_bytes = _allocated_bytes(() -> begin
        mul!(extrapolated_command, extrapolation_backend, active_command)
        _scatter_active_command!(dm.state.coefs, extrapolated_command, active_indices_backend)
        apply!(dm, tel, DMReplace())
        _sync_backend!(backend_tag, tel.state.opd)
    end; warmup=warmup)
    sense_alloc_bytes = _allocated_bytes(() -> begin
        mul!(extrapolated_command, extrapolation_backend, active_command)
        _scatter_active_command!(dm.state.coefs, extrapolated_command, active_indices_backend)
        apply!(dm, tel, DMReplace())
        measure!(wfs, tel, src, det; rng=Random.default_rng())
        _sync_backend!(backend_tag, wfs.state.spot_cube)
    end; warmup=warmup)
    mosaic_alloc_bytes = _allocated_bytes(() -> begin
        mul!(extrapolated_command, extrapolation_backend, active_command)
        _scatter_active_command!(dm.state.coefs, extrapolated_command, active_indices_backend)
        apply!(dm, tel, DMReplace())
        measure!(wfs, tel, src, det; rng=Random.default_rng())
        _tile_spot_cube!(tiled_frame, wfs.state.spot_cube, n_subap, roi)
        _sync_backend!(backend_tag, tiled_frame)
    end; warmup=warmup)
    total_alloc_bytes = _allocated_bytes(_step!; warmup=warmup)
    metadata = detector_export_metadata(det)

    println("revolt_like_hil_runtime_profile")
    println("  backend: ", label)
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
    println("  extrapolated_command_length: ", length(extrapolated_command))
    println("  dm_grid_command_length: ", length(dm.state.coefs))
    println("  dm_layout_shape: ", size(actuator_map))
    println("  pupil_resolution: ", resolution)
    println("  n_subap: ", n_subap)
    println("  roi_size: ", roi)
    println("  tiled_wfs_frame_shape: ", size(tiled_frame))
    println("  spot_cube_shape: ", size(wfs.state.spot_cube))
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

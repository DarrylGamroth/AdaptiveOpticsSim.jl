using AdaptiveOpticsSim
using KernelAbstractions: @index, @kernel
using Random
using Statistics: mean

include(joinpath(@__DIR__, "revolt", "common.jl"))

const _backend_arg = isempty(ARGS) ? "cpu" : lowercase(ARGS[1])
const _model_arg = length(ARGS) >= 2 ? lowercase(ARGS[2]) : "all"
const _response_arg = length(ARGS) >= 3 ? lowercase(ARGS[3]) : "default"
const _samples_arg = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 6
const _warmup_arg = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 2

if _backend_arg == "cuda"
    import CUDA
end

if _backend_arg == "amdgpu"
    import AMDGPU
end

function _resolve_backend(name::AbstractString)
    lowered = lowercase(name)
    if lowered == "cpu"
        return Array, nothing, "cpu"
    elseif lowered == "cuda"
        isdefined(Main, :CUDA) || error("profile_revolt_hil_runtime.jl requires CUDA.jl for backend=cuda")
        CUDA.functional() || error("profile_revolt_hil_runtime.jl requires a functional CUDA driver/device")
        AdaptiveOpticsSim.disable_scalar_backend!(CUDABackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(CUDABackendTag)
        backend === nothing && error("CUDA backend array type is unavailable")
        return backend, CUDABackendTag, "cuda"
    elseif lowered == "amdgpu"
        isdefined(Main, :AMDGPU) || error("profile_revolt_hil_runtime.jl requires AMDGPU.jl for backend=amdgpu")
        AMDGPU.functional() || error("profile_revolt_hil_runtime.jl requires a functional ROCm installation and GPU")
        AdaptiveOpticsSim.disable_scalar_backend!(AMDGPUBackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(AMDGPUBackendTag)
        backend === nothing && error("AMDGPU backend array type is unavailable")
        return backend, AMDGPUBackendTag, "amdgpu"
    end
    error("unsupported backend '$name'; use cpu, cuda, or amdgpu")
end

_sync_backend!(::Nothing, _) = nothing

function _sync_backend!(::Type{B}, array) where {B<:GPUBackendTag}
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

function _resolve_models(name::AbstractString)
    lowered = lowercase(name)
    lowered == "all" && return collect(REVOLT_MODELS)
    lowered == "pwfs" && return [:pwfs]
    lowered == "pwfs_unmod" && return [:pwfs_unmod]
    lowered == "pwfs-unmod" && return [:pwfs_unmod]
    lowered == "pwfs_unmodulated" && return [:pwfs_unmod]
    lowered == "unmod" && return [:pwfs_unmod]
    lowered == "shwfs" && return [:shwfs]
    error("unsupported REVOLT model '$name'; use all, pwfs, pwfs_unmod, or shwfs")
end

function _resolve_response_mode(name::AbstractString)
    lowered = lowercase(name)
    lowered == "default" && return :default
    lowered == "null" && return :null
    error("unsupported response mode '$name'; use default or null")
end

@kernel function _fill_dm_command_kernel!(coefs, phase, n::Int)
    i = @index(Global, Linear)
    if i <= n
        T = eltype(coefs)
        @inbounds coefs[i] = T(3e-8) * (sin(T(0.017) * T(i) + phase) + T(0.25) * cos(T(0.031) * T(i) - phase))
    end
end

function _fill_dm_command!(::AdaptiveOpticsSim.ScalarCPUStyle, coefs::AbstractVector{T}, phase::T) where {T<:AbstractFloat}
    @inbounds for i in eachindex(coefs)
        coefs[i] = T(3e-8) * (sin(T(0.017) * T(i) + phase) + T(0.25) * cos(T(0.031) * T(i) - phase))
    end
    return coefs
end

function _fill_dm_command!(style::AdaptiveOpticsSim.AcceleratorStyle, coefs::AbstractVector{T}, phase::T) where {T<:AbstractFloat}
    AdaptiveOpticsSim.launch_kernel!(style, _fill_dm_command_kernel!, coefs, phase, length(coefs); ndrange=length(coefs))
    return coefs
end

function _fill_dm_command!(coefs::AbstractVector{T}, phase::T) where {T<:AbstractFloat}
    return _fill_dm_command!(AdaptiveOpticsSim.execution_style(coefs), coefs, phase)
end

function _set_dm_command!(runtime, phase::T) where {T<:AbstractFloat}
    return _fill_dm_command!(runtime.dm.state.coefs, phase)
end

function _propagate_atmosphere!(runtime)
    propagate!(runtime.atmosphere, runtime.telescope)
    return runtime.telescope
end

function _sync_profile_state!(backend_tag, runtime)
    # The real external boundary is the detector pixel output. Synchronizing the
    # final detector frame is sufficient to fence the full step on the backend
    # stream without separately synchronizing internal telescope/DM/WFS buffers.
    _sync_backend!(backend_tag, output_frame(runtime.wfs_detector))
    return nothing
end

_wfs_frame(runtime) = AdaptiveOpticsSim.wfs_output_frame(runtime.wfs, runtime.wfs_detector)

function _phase_step!(runtime, backend_tag; phase_index::Int=1)
    advance!(runtime.atmosphere, runtime.telescope; rng=runtime.rng)
    _propagate_atmosphere!(runtime)
    _set_dm_command!(runtime, eltype(runtime.dm.state.coefs)(phase_index))
    apply!(runtime.dm, runtime.telescope, DMAdditive())
    measure!(runtime.wfs, runtime.telescope, runtime.sky_source, runtime.wfs_detector; rng=runtime.rng)
    _sync_profile_state!(backend_tag, runtime)
    return nothing
end

function _report_profile(runtime, backend_label::AbstractString, response_mode::Symbol,
    build_time_ns, atmosphere_mean_ns, atmosphere_p95_ns, dm_apply_mean_ns, dm_apply_p95_ns,
    sense_mean_ns, sense_p95_ns, total_mean_ns, total_p95_ns, atmosphere_alloc_bytes,
    dm_apply_alloc_bytes, sense_alloc_bytes, total_alloc_bytes)
    wfs_metadata = detector_export_metadata(runtime.wfs_detector)
    science_metadata = detector_export_metadata(runtime.science_detector)
    gain_metadata = isnothing(runtime.gain_detector) ? nothing : detector_export_metadata(runtime.gain_detector)
    wfs_frame = _wfs_frame(runtime)
    wfs_camera = runtime.camera_configs[runtime.wfs_detector_name]
    wfs_reference_frame = detector_reference_frame(runtime.wfs_detector)
    wfs_signal_frame = detector_signal_frame(runtime.wfs_detector)
    wfs_combined_frame = detector_combined_frame(runtime.wfs_detector)
    wfs_reference_cube = detector_reference_cube(runtime.wfs_detector)
    wfs_signal_cube = detector_signal_cube(runtime.wfs_detector)
    wfs_read_cube = detector_read_cube(runtime.wfs_detector)

    println("revolt_hil_runtime_profile")
    println("  backend: ", backend_label)
    println("  model: ", runtime.label)
    println("  response_mode: ", response_mode)
    println("  wfs_family: ", runtime.wfs_family)
    println("  wfs_detector: ", runtime.wfs_detector_name)
    println("  wfs_detector_sensor: ", wfs_metadata.sensor)
    println("  wfs_detector_nominal_resolution: ", wfs_camera.resolution)
    println("  wfs_frame_response: ", wfs_metadata.frame_response)
    println("  science_detector: ", runtime.science_detector_name)
    println("  science_detector_sensor: ", science_metadata.sensor)
    println("  science_detector_nominal_resolution: ", runtime.camera_configs[runtime.science_detector_name].resolution)
    if !isnothing(gain_metadata)
        println("  gain_detector: ", runtime.gain_detector_name)
        println("  gain_detector_sensor: ", gain_metadata.sensor)
        println("  gain_detector_nominal_resolution: ", runtime.camera_configs[runtime.gain_detector_name].resolution)
    end
    println("  pupil_resolution: ", runtime.telescope.params.resolution)
    println("  n_subap: ", runtime.n_subap)
    println("  dm_grid_shape: ", (runtime.dm.params.n_act, runtime.dm.params.n_act))
    println("  dm_command_length: ", length(runtime.dm.state.coefs))
    println("  slope_length: ", length(runtime.wfs.state.slopes))
    println("  source_band: ", runtime.sky_source.params.band)
    println("  effective_sky_magnitude: ", runtime.sky_source.params.magnitude)
    println("  effective_calibration_magnitude: ", runtime.calibration_source.params.magnitude)
    println("  wfs_frame_shape: ", size(wfs_frame))
    println("  detector_output_shape: ", size(output_frame(runtime.wfs_detector)))
    println("  wfs_reference_frame_shape: ", isnothing(wfs_reference_frame) ? nothing : size(wfs_reference_frame))
    println("  wfs_signal_frame_shape: ", isnothing(wfs_signal_frame) ? nothing : size(wfs_signal_frame))
    println("  wfs_combined_frame_shape: ", isnothing(wfs_combined_frame) ? nothing : size(wfs_combined_frame))
    println("  wfs_reference_cube_shape: ", isnothing(wfs_reference_cube) ? nothing : size(wfs_reference_cube))
    println("  wfs_signal_cube_shape: ", isnothing(wfs_signal_cube) ? nothing : size(wfs_signal_cube))
    println("  wfs_read_cube_shape: ", isnothing(wfs_read_cube) ? nothing : size(wfs_read_cube))
    println("  build_time_ns: ", build_time_ns)
    println("  atmosphere_mean_ns: ", atmosphere_mean_ns)
    println("  atmosphere_p95_ns: ", atmosphere_p95_ns)
    println("  dm_apply_mean_ns: ", dm_apply_mean_ns)
    println("  dm_apply_p95_ns: ", dm_apply_p95_ns)
    println("  sense_mean_ns: ", sense_mean_ns)
    println("  sense_p95_ns: ", sense_p95_ns)
    println("  total_mean_ns: ", total_mean_ns)
    println("  total_p95_ns: ", total_p95_ns)
    println("  frame_rate_hz: ", 1.0e9 / total_mean_ns)
    println("  atmosphere_alloc_bytes: ", atmosphere_alloc_bytes)
    println("  dm_apply_alloc_bytes: ", dm_apply_alloc_bytes)
    println("  sense_alloc_bytes: ", sense_alloc_bytes)
    println("  total_alloc_bytes: ", total_alloc_bytes)
    return nothing
end

function run_profile(kind::Symbol; backend_name::AbstractString="cpu", response_name::AbstractString="default",
    samples::Int=6, warmup::Int=2)
    BackendArray, backend_tag, backend_label = _resolve_backend(backend_name)
    response_mode = _resolve_response_mode(response_name)
    runtime = revolt_setup(kind; T=Float32, backend=BackendArray, response_mode=response_mode)

    build_time_ns = time_ns()
    _phase_step!(runtime, backend_tag; phase_index=1)
    _phase_step!(runtime, backend_tag; phase_index=2)
    build_time_ns = time_ns() - build_time_ns

    phase_index = Ref(0)
    atmosphere_mean_ns, atmosphere_p95_ns = _timed_stats!(() -> begin
        phase_index[] += 1
        advance!(runtime.atmosphere, runtime.telescope; rng=runtime.rng)
        _propagate_atmosphere!(runtime)
        _sync_backend!(backend_tag, runtime.telescope.state.opd)
    end; warmup=warmup, samples=samples)
    dm_apply_mean_ns, dm_apply_p95_ns = _timed_stats!(() -> begin
        phase_index[] += 1
        advance!(runtime.atmosphere, runtime.telescope; rng=runtime.rng)
        _propagate_atmosphere!(runtime)
        _set_dm_command!(runtime, Float32(phase_index[]))
        apply!(runtime.dm, runtime.telescope, DMAdditive())
        _sync_backend!(backend_tag, runtime.telescope.state.opd)
    end; warmup=warmup, samples=samples)
    sense_mean_ns, sense_p95_ns = _timed_stats!(() -> begin
        phase_index[] += 1
        _phase_step!(runtime, backend_tag; phase_index=phase_index[])
    end; warmup=warmup, samples=samples)
    total_mean_ns, total_p95_ns = _timed_stats!(() -> begin
        phase_index[] += 1
        _phase_step!(runtime, backend_tag; phase_index=phase_index[])
    end; warmup=warmup, samples=samples)

    atmosphere_alloc_bytes = _allocated_bytes(() -> begin
        phase_index[] += 1
        advance!(runtime.atmosphere, runtime.telescope; rng=runtime.rng)
        _propagate_atmosphere!(runtime)
        _sync_backend!(backend_tag, runtime.telescope.state.opd)
    end; warmup=warmup)
    dm_apply_alloc_bytes = _allocated_bytes(() -> begin
        phase_index[] += 1
        advance!(runtime.atmosphere, runtime.telescope; rng=runtime.rng)
        _propagate_atmosphere!(runtime)
        _set_dm_command!(runtime, Float32(phase_index[]))
        apply!(runtime.dm, runtime.telescope, DMAdditive())
        _sync_backend!(backend_tag, runtime.telescope.state.opd)
    end; warmup=warmup)
    sense_alloc_bytes = _allocated_bytes(() -> begin
        phase_index[] += 1
        _phase_step!(runtime, backend_tag; phase_index=phase_index[])
    end; warmup=warmup)
    total_alloc_bytes = _allocated_bytes(() -> begin
        phase_index[] += 1
        _phase_step!(runtime, backend_tag; phase_index=phase_index[])
    end; warmup=warmup)

    _report_profile(runtime, backend_label, response_mode, build_time_ns,
        atmosphere_mean_ns, atmosphere_p95_ns, dm_apply_mean_ns, dm_apply_p95_ns,
        sense_mean_ns, sense_p95_ns, total_mean_ns, total_p95_ns, atmosphere_alloc_bytes,
        dm_apply_alloc_bytes, sense_alloc_bytes, total_alloc_bytes)
    return nothing
end

for model in _resolve_models(_model_arg)
    run_profile(model; backend_name=_backend_arg, response_name=_response_arg,
        samples=_samples_arg, warmup=_warmup_arg)
end

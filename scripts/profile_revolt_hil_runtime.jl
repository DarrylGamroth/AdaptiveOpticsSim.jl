using AdaptiveOpticsSim
using Random
using Statistics: mean
using KernelAbstractions: @kernel, @index

const _backend_arg = isempty(ARGS) ? "cpu" : lowercase(ARGS[1])

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

function _active_support_indices(n_act::Int, n_active::Int)
    total = n_act * n_act
    0 < n_active <= total || error("n_active must be in 1:$total")
    xs = range(-1.0, 1.0; length=n_act)
    ys = range(-1.0, 1.0; length=n_act)
    scores = Vector{Tuple{Float64,Int}}(undef, total)
    idx = 1
    for x in xs, y in ys
        scores[idx] = (x^2 + y^2, idx)
        idx += 1
    end
    sort!(scores; by=first)
    active = sort(map(last, scores[1:n_active]))
    return active
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

function run_profile(; backend_name::AbstractString="cpu", samples::Int=6, warmup::Int=2)
    BackendArray, backend_tag, label = _resolve_backend(backend_name)
    T = Float32
    n_subap = 16
    roi = 22
    resolution = 352
    n_act = 17
    n_active = 277

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
    det = Detector(noise=NoiseNone(), integration_time=T(1), qe=T(1), binning=1, T=T, backend=BackendArray)
    active_indices_host = _active_support_indices(n_act, n_active)
    active_indices_backend = BackendArray{Int}(undef, n_active)
    copyto!(active_indices_backend, active_indices_host)
    active_command = BackendArray{T}(undef, n_active)
    _fill_active_command!(active_command)
    tiled_frame = BackendArray{T}(undef, resolution, resolution)

    build_time_ns = time_ns()
    AdaptiveOpticsSim.ensure_sh_calibration!(wfs, tel, src)
    _sync_backend!(backend_tag, wfs.state.slopes)
    build_time_ns = time_ns() - build_time_ns

    function _step!()
        reset_opd!(tel)
        _scatter_active_command!(dm.state.coefs, active_command, active_indices_backend)
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
        _scatter_active_command!(dm.state.coefs, active_command, active_indices_backend)
        _sync_backend!(backend_tag, dm.state.coefs)
    end; warmup=warmup, samples=samples)
    dm_apply_mean_ns, dm_apply_p95_ns = _timed_stats!(() -> begin
        _scatter_active_command!(dm.state.coefs, active_command, active_indices_backend)
        apply!(dm, tel, DMReplace())
        _sync_backend!(backend_tag, tel.state.opd)
    end; warmup=warmup, samples=samples)
    sense_mean_ns, sense_p95_ns = _timed_stats!(() -> begin
        _scatter_active_command!(dm.state.coefs, active_command, active_indices_backend)
        apply!(dm, tel, DMReplace())
        measure!(wfs, tel, src, det; rng=Random.default_rng())
        _sync_backend!(backend_tag, wfs.state.spot_cube)
    end; warmup=warmup, samples=samples)
    mosaic_mean_ns, mosaic_p95_ns = _timed_stats!(() -> begin
        _scatter_active_command!(dm.state.coefs, active_command, active_indices_backend)
        apply!(dm, tel, DMReplace())
        measure!(wfs, tel, src, det; rng=Random.default_rng())
        _tile_spot_cube!(tiled_frame, wfs.state.spot_cube, n_subap, roi)
        _sync_backend!(backend_tag, tiled_frame)
    end; warmup=warmup, samples=samples)
    total_mean_ns, total_p95_ns = _timed_stats!(_step!; warmup=warmup, samples=samples)

    println("revolt_hil_runtime_profile")
    println("  backend: ", label)
    println("  actuator_command_length: ", n_active)
    println("  dm_grid_command_length: ", length(dm.state.coefs))
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
    return nothing
end

run_profile(; backend_name=_backend_arg)

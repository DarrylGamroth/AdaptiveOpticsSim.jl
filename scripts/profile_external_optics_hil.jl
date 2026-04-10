using AdaptiveOpticsSim
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

@kernel function export_phase_crop_kernel!(dest, src, col_offset::Int)
    i, j = @index(Global, NTuple)
    if i <= size(dest, 1) && j <= size(dest, 2)
        @inbounds dest[i, j] = src[i, j + col_offset]
    end
end

function _resolve_backend(name::AbstractString)
    lowered = lowercase(name)
    if lowered == "cpu"
        return Array, nothing, "cpu"
    elseif lowered == "cuda"
        isdefined(Main, :CUDA) || error("profile_external_optics_hil.jl requires CUDA.jl for backend=cuda")
        CUDA.functional() || error("profile_external_optics_hil.jl requires a functional CUDA driver/device")
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.CUDABackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(AdaptiveOpticsSim.CUDABackendTag)
        backend === nothing && error("CUDA backend array type is unavailable")
        return backend, AdaptiveOpticsSim.CUDABackendTag, "cuda"
    elseif lowered == "amdgpu"
        isdefined(Main, :AMDGPU) || error("profile_external_optics_hil.jl requires AMDGPU.jl for backend=amdgpu")
        AMDGPU.functional() || error("profile_external_optics_hil.jl requires a functional ROCm installation and GPU")
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
        host[i] = T(0.5 * sin(0.011 * i) + 0.35 * cos(0.023 * i))
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

function _export_phase_crop!(dest::AbstractMatrix{T}, src::AbstractMatrix{T}) where {T<:AbstractFloat}
    col_offset = div(size(src, 2) - size(dest, 2), 2)
    style = AdaptiveOpticsSim.execution_style(dest)
    if dest isa Array
        copyto!(dest, @view(src[:, col_offset + 1:col_offset + size(dest, 2)]))
    else
        AdaptiveOpticsSim.launch_kernel!(style, export_phase_crop_kernel!, dest, src, col_offset;
            ndrange=size(dest))
    end
    return dest
end

function run_profile(; backend_name::AbstractString="cpu", samples::Int=6, warmup::Int=2)
    BackendArray, backend_tag, label = _resolve_backend(backend_name)
    T = Float32
    phase_resolution = 640
    external_shape = (640, 512)
    n_act = 24
    n_active = 468

    tel = Telescope(
        resolution=phase_resolution,
        diameter=8.0,
        sampling_time=1e-3,
        central_obstruction=0.0,
        T=T,
        backend=BackendArray,
    )
    dm = DeformableMirror(tel; n_act=n_act, influence_width=0.25, T=T, backend=BackendArray)
    active_indices_host = _active_support_indices(n_act, n_active)
    active_indices_backend = BackendArray{Int}(undef, n_active)
    copyto!(active_indices_backend, active_indices_host)
    active_command = BackendArray{T}(undef, n_active)
    _fill_active_command!(active_command)
    export_phase = BackendArray{T}(undef, external_shape...)

    build_time_ns = 0

    function _step!()
        _scatter_active_command!(dm.state.coefs, active_command, active_indices_backend)
        apply_opd!(dm, tel)
        _export_phase_crop!(export_phase, dm.state.opd)
        _sync_backend!(backend_tag, dm.state.opd)
        _sync_backend!(backend_tag, export_phase)
        return nothing
    end

    command_map_mean_ns, command_map_p95_ns = _timed_stats!(() -> begin
        _scatter_active_command!(dm.state.coefs, active_command, active_indices_backend)
        _sync_backend!(backend_tag, dm.state.coefs)
    end; warmup=warmup, samples=samples)
    dm_phase_mean_ns, dm_phase_p95_ns = _timed_stats!(() -> begin
        _scatter_active_command!(dm.state.coefs, active_command, active_indices_backend)
        apply_opd!(dm, tel)
        _sync_backend!(backend_tag, dm.state.opd)
    end; warmup=warmup, samples=samples)
    export_mean_ns, export_p95_ns = _timed_stats!(() -> begin
        _scatter_active_command!(dm.state.coefs, active_command, active_indices_backend)
        apply_opd!(dm, tel)
        _export_phase_crop!(export_phase, dm.state.opd)
        _sync_backend!(backend_tag, export_phase)
    end; warmup=warmup, samples=samples)
    total_mean_ns, total_p95_ns = _timed_stats!(_step!; warmup=warmup, samples=samples)
    command_map_alloc_bytes = _allocated_bytes(() -> begin
        _scatter_active_command!(dm.state.coefs, active_command, active_indices_backend)
        _sync_backend!(backend_tag, dm.state.coefs)
    end; warmup=warmup)
    dm_phase_alloc_bytes = _allocated_bytes(() -> begin
        _scatter_active_command!(dm.state.coefs, active_command, active_indices_backend)
        apply_opd!(dm, tel)
        _sync_backend!(backend_tag, dm.state.opd)
    end; warmup=warmup)
    export_alloc_bytes = _allocated_bytes(() -> begin
        _scatter_active_command!(dm.state.coefs, active_command, active_indices_backend)
        apply_opd!(dm, tel)
        _export_phase_crop!(export_phase, dm.state.opd)
        _sync_backend!(backend_tag, export_phase)
    end; warmup=warmup)
    total_alloc_bytes = _allocated_bytes(_step!; warmup=warmup)

    println("external_optics_hil_profile")
    println("  backend: ", label)
    println("  actuator_command_length: ", n_active)
    println("  dm_grid_command_length: ", length(dm.state.coefs))
    println("  phase_grid_shape: ", size(dm.state.opd))
    println("  external_phase_export_shape: ", size(export_phase))
    println("  downstream_image_shape: ", external_shape)
    println("  build_time_ns: ", build_time_ns)
    println("  command_map_mean_ns: ", command_map_mean_ns)
    println("  command_map_p95_ns: ", command_map_p95_ns)
    println("  dm_phase_mean_ns: ", dm_phase_mean_ns)
    println("  dm_phase_p95_ns: ", dm_phase_p95_ns)
    println("  export_mean_ns: ", export_mean_ns)
    println("  export_p95_ns: ", export_p95_ns)
    println("  total_mean_ns: ", total_mean_ns)
    println("  total_p95_ns: ", total_p95_ns)
    println("  frame_rate_hz: ", 1.0e9 / total_mean_ns)
    println("  command_map_alloc_bytes: ", command_map_alloc_bytes)
    println("  dm_phase_alloc_bytes: ", dm_phase_alloc_bytes)
    println("  export_alloc_bytes: ", export_alloc_bytes)
    println("  total_alloc_bytes: ", total_alloc_bytes)
    return nothing
end

run_profile(; backend_name=_backend_arg)

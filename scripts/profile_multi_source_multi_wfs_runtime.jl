using AdaptiveOpticsSim
using Random
using Statistics

const _backend_arg = isempty(ARGS) ? "cpu" : lowercase(ARGS[1])
const _scale_arg = length(ARGS) >= 2 ? lowercase(ARGS[2]) : "compact"

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
        isdefined(Main, :CUDA) || (@eval import CUDA)
        CUDA.functional() || error("profile_multi_source_multi_wfs_runtime.jl requires a functional CUDA driver/device")
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.CUDABackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(AdaptiveOpticsSim.CUDABackendTag)
        backend === nothing && error("CUDA backend array type is unavailable")
        return backend, AdaptiveOpticsSim.CUDABackendTag, "cuda"
    elseif lowered == "amdgpu"
        isdefined(Main, :AMDGPU) || (@eval import AMDGPU)
        AMDGPU.functional() || error("profile_multi_source_multi_wfs_runtime.jl requires a functional ROCm installation and GPU")
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

function _timed_stats!(f!::F; warmup::Int=5, samples::Int=20) where {F<:Function}
    for _ in 1:warmup
        f!()
    end
    GC.gc()
    timings = Vector{Int}(undef, samples)
    alloc_bytes = Vector{Int}(undef, samples)
    @inbounds for i in 1:samples
        gc_enable = GC.enable(false)
        try
            t0 = time_ns()
            alloc_bytes[i] = @allocated f!()
            timings[i] = time_ns() - t0
        finally
            GC.enable(gc_enable)
        end
    end
    return mean(timings), sort(timings)[clamp(round(Int, 0.95 * samples), 1, samples)], Int(round(mean(alloc_bytes)))
end

function _resolve_scale(name::AbstractString)
    lowered = lowercase(name)
    lowered in ("compact", "medium", "representative") ||
        error("unsupported scale '$name'; use compact, medium, or representative")
    return lowered
end

function _multi_source_scale_config(name::AbstractString)
    scale = _resolve_scale(name)
    if scale == "compact"
        return (
            scale=scale,
            asterism_resolution=32,
            asterism_n_subap=4,
            source_coords=((0.0, 0.0), (1.0, 45.0)),
            runtime_resolution=16,
            runtime_n_subap=4,
            runtime_dm_n_act=4,
            branch_count=2,
            warmup=5,
            samples=20,
        )
    elseif scale == "medium"
        return (
            scale=scale,
            asterism_resolution=64,
            asterism_n_subap=8,
            source_coords=((0.0, 0.0), (1.0, 45.0), (1.5, 135.0), (2.0, 225.0)),
            runtime_resolution=32,
            runtime_n_subap=8,
            runtime_dm_n_act=8,
            branch_count=3,
            warmup=3,
            samples=12,
        )
    end
    return (
        scale=scale,
        asterism_resolution=96,
        asterism_n_subap=12,
        source_coords=((0.0, 0.0), (1.0, 45.0), (1.5, 135.0), (2.0, 225.0), (2.5, 315.0), (3.0, 90.0)),
        runtime_resolution=72,
        runtime_n_subap=12,
        runtime_dm_n_act=12,
        branch_count=4,
        warmup=2,
        samples=6,
    )
end

function _build_recon(dm, wfs, tel, src)
    imat = interaction_matrix(dm, wfs, tel, src; amplitude=eltype(dm.state.coefs)(0.05))
    return ModalReconstructor(imat; gain=eltype(dm.state.coefs)(0.5))
end

function _seed_opd!(tel::Telescope{<:Any,<:Any}, ::Type{T}) where {T<:AbstractFloat}
    host = Matrix{T}(undef, tel.params.resolution, tel.params.resolution)
    @inbounds for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        host[i, j] = T(i + j / 10)
    end
    copyto!(tel.state.opd, host)
    return tel.state.opd
end

function _source_list(coords, T::Type{<:AbstractFloat})
    return [Source(band=:I, magnitude=0.0, coordinates=(coord[1], coord[2]), T=T) for coord in coords]
end

function _build_runtime_branch(::Type{T}, BackendArray, resolution::Int, n_subap::Int, n_act::Int,
    seed::Integer; sensor::Symbol=:sh) where {T<:AbstractFloat}
    tel = Telescope(resolution=resolution, diameter=T(8.0), sampling_time=T(1e-3),
        central_obstruction=T(0.0), T=T, backend=BackendArray)
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = KolmogorovAtmosphere(tel; r0=T(0.2), L0=T(25.0), T=T, backend=BackendArray)
    dm = DeformableMirror(tel; n_act=n_act, influence_width=T(0.3), T=T, backend=BackendArray)
    wfs = sensor == :sh ?
        ShackHartmann(tel; n_subap=n_subap, mode=Diffractive(), T=T, backend=BackendArray) :
        PyramidWFS(tel; n_subap=n_subap, modulation=T(1.0), mode=Diffractive(), T=T, backend=BackendArray)
    det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1, T=T, backend=BackendArray)
    sim = AdaptiveOpticsSim.AOSimulation(tel, atm, src, dm, wfs)
    return ClosedLoopBranchConfig(
        Symbol("branch_", seed),
        sim,
        _build_recon(dm, wfs, tel, src);
        wfs_detector=det,
        rng=MersenneTwister(seed),
    )
end

function run_profile(; backend_name::AbstractString="cpu", scale_name::AbstractString="compact",
    T::Type{<:AbstractFloat}=Float32, samples::Union{Int,Nothing}=nothing,
    warmup::Union{Int,Nothing}=nothing)
    BackendArray, backend_tag, label = _resolve_backend(backend_name)
    rng = MersenneTwister(7)
    cfg = _multi_source_scale_config(scale_name)
    resolved_samples = something(samples, cfg.samples)
    resolved_warmup = something(warmup, cfg.warmup)

    tel = Telescope(resolution=cfg.asterism_resolution, diameter=T(8.0), sampling_time=T(1e-3),
        central_obstruction=T(0.0), T=T, backend=BackendArray)
    ast = Asterism(_source_list(cfg.source_coords, T))
    det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1, T=T, backend=BackendArray)

    _seed_opd!(tel, T)

    sh = ShackHartmann(tel; n_subap=cfg.asterism_n_subap, mode=Diffractive(), T=T, backend=BackendArray)
    sh_mean_ns, sh_p95_ns, sh_alloc_bytes = _timed_stats!(() -> begin
        measure!(sh, tel, ast, det; rng=rng)
        _sync_backend!(backend_tag, sh.state.slopes)
    end; warmup=resolved_warmup, samples=resolved_samples)

    pyr = PyramidWFS(tel; n_subap=cfg.asterism_n_subap, modulation=T(1.0), mode=Diffractive(), T=T, backend=BackendArray)
    pyr_mean_ns, pyr_p95_ns, pyr_alloc_bytes = _timed_stats!(() -> begin
        measure!(pyr, tel, ast, det; rng=rng)
        _sync_backend!(backend_tag, pyr.state.slopes)
    end; warmup=resolved_warmup, samples=resolved_samples)

    bio = BioEdgeWFS(tel; n_subap=cfg.asterism_n_subap, modulation=T(1.0), mode=Diffractive(), T=T, backend=BackendArray)
    bio_mean_ns, bio_p95_ns, bio_alloc_bytes = _timed_stats!(() -> begin
        measure!(bio, tel, ast, det; rng=rng)
        _sync_backend!(backend_tag, bio.state.slopes)
    end; warmup=resolved_warmup, samples=resolved_samples)

    compatible_branches = [_build_runtime_branch(T, BackendArray,
        cfg.runtime_resolution, cfg.runtime_n_subap, cfg.runtime_dm_n_act, 10 + i; sensor=:sh)
        for i in 1:cfg.branch_count]
    mixed_branches = vcat(
        [_build_runtime_branch(T, BackendArray,
            cfg.runtime_resolution, cfg.runtime_n_subap, cfg.runtime_dm_n_act, 30 + i; sensor=:sh)
            for i in 1:max(cfg.branch_count - 1, 1)],
        [_build_runtime_branch(T, BackendArray,
            cfg.runtime_resolution, cfg.runtime_n_subap, cfg.runtime_dm_n_act, 40; sensor=:pyr)],
    )
    compatible_cfg = GroupedPlatformConfig(
        Tuple(map(branch -> branch.label, compatible_branches));
        name=:compatible_grouped_runtime,
        products=GroupedRuntimeProductRequirements(wfs_frames=true, science_frames=false, wfs_stack=true, science_stack=false),
    )
    mixed_cfg = GroupedPlatformConfig(
        Tuple(map(branch -> branch.label, mixed_branches));
        name=:mixed_grouped_runtime,
        products=GroupedRuntimeProductRequirements(wfs_frames=true, science_frames=false, wfs_stack=false, science_stack=false),
    )
    composite_compatible = build_platform_scenario(compatible_cfg, compatible_branches...)
    composite_mixed = build_platform_scenario(mixed_cfg, mixed_branches...)
    prepare!(composite_compatible)
    prepare!(composite_mixed)
    step!(composite_compatible)
    step!(composite_mixed)

    comp_mean_ns, comp_p95_ns, comp_alloc_bytes = _timed_stats!(() -> begin
        step!(composite_compatible)
        _sync_backend!(backend_tag, simulation_command(composite_compatible))
    end; warmup=resolved_warmup, samples=resolved_samples)
    mixed_mean_ns, mixed_p95_ns, mixed_alloc_bytes = _timed_stats!(() -> begin
        step!(composite_mixed)
        _sync_backend!(backend_tag, simulation_command(composite_mixed))
    end; warmup=resolved_warmup, samples=resolved_samples)

    compatible_stack = simulation_grouped_wfs_stack(composite_compatible)
    mixed_stack = simulation_grouped_wfs_stack(composite_mixed)

    return (
        backend=label,
        scale=cfg.scale,
        sh_asterism_mean_ns=sh_mean_ns,
        sh_asterism_p95_ns=sh_p95_ns,
        sh_asterism_frame_rate_hz=1e9 / sh_mean_ns,
        sh_asterism_alloc_bytes=sh_alloc_bytes,
        pyramid_asterism_mean_ns=pyr_mean_ns,
        pyramid_asterism_p95_ns=pyr_p95_ns,
        pyramid_asterism_frame_rate_hz=1e9 / pyr_mean_ns,
        pyramid_asterism_alloc_bytes=pyr_alloc_bytes,
        bioedge_asterism_mean_ns=bio_mean_ns,
        bioedge_asterism_p95_ns=bio_p95_ns,
        bioedge_asterism_frame_rate_hz=1e9 / bio_mean_ns,
        bioedge_asterism_alloc_bytes=bio_alloc_bytes,
        composite_compatible_mean_ns=comp_mean_ns,
        composite_compatible_p95_ns=comp_p95_ns,
        composite_compatible_frame_rate_hz=1e9 / comp_mean_ns,
        composite_compatible_alloc_bytes=comp_alloc_bytes,
        composite_mixed_mean_ns=mixed_mean_ns,
        composite_mixed_p95_ns=mixed_p95_ns,
        composite_mixed_frame_rate_hz=1e9 / mixed_mean_ns,
        composite_mixed_alloc_bytes=mixed_alloc_bytes,
        source_count=length(cfg.source_coords),
        asterism_resolution=cfg.asterism_resolution,
        asterism_n_subap=cfg.asterism_n_subap,
        runtime_resolution=cfg.runtime_resolution,
        runtime_n_subap=cfg.runtime_n_subap,
        runtime_branch_count=cfg.branch_count,
        compatible_grouped_wfs_stack=!isnothing(compatible_stack),
        mixed_grouped_wfs_stack=!isnothing(mixed_stack),
        compatible_grouped_wfs_stack_shape=isnothing(compatible_stack) ? "none" : string(size(compatible_stack)),
        mixed_grouped_wfs_stack_shape=isnothing(mixed_stack) ? "none" : string(size(mixed_stack)),
    )
end

function print_profile(result)
    println("multi_source_multi_wfs_profile")
    for key in (
        :backend,
        :scale,
        :sh_asterism_mean_ns,
        :sh_asterism_p95_ns,
        :sh_asterism_frame_rate_hz,
        :sh_asterism_alloc_bytes,
        :pyramid_asterism_mean_ns,
        :pyramid_asterism_p95_ns,
        :pyramid_asterism_frame_rate_hz,
        :pyramid_asterism_alloc_bytes,
        :bioedge_asterism_mean_ns,
        :bioedge_asterism_p95_ns,
        :bioedge_asterism_frame_rate_hz,
        :bioedge_asterism_alloc_bytes,
        :composite_compatible_mean_ns,
        :composite_compatible_p95_ns,
        :composite_compatible_frame_rate_hz,
        :composite_compatible_alloc_bytes,
        :composite_mixed_mean_ns,
        :composite_mixed_p95_ns,
        :composite_mixed_frame_rate_hz,
        :composite_mixed_alloc_bytes,
        :source_count,
        :asterism_resolution,
        :asterism_n_subap,
        :runtime_resolution,
        :runtime_n_subap,
        :runtime_branch_count,
        :compatible_grouped_wfs_stack,
        :compatible_grouped_wfs_stack_shape,
        :mixed_grouped_wfs_stack,
        :mixed_grouped_wfs_stack_shape,
    )
        println("  ", key, ": ", getproperty(result, key))
    end
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    print_profile(run_profile(; backend_name=_backend_arg, scale_name=_scale_arg))
end

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
        isdefined(Main, :CUDA) || error("profile_multi_source_multi_wfs_runtime.jl requires CUDA.jl for backend=cuda")
        CUDA.functional() || error("profile_multi_source_multi_wfs_runtime.jl requires a functional CUDA driver/device")
        AdaptiveOpticsSim.disable_scalar_backend!(CUDABackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(CUDABackendTag)
        backend === nothing && error("CUDA backend array type is unavailable")
        return backend, CUDABackendTag, "cuda"
    elseif lowered == "amdgpu"
        isdefined(Main, :AMDGPU) || error("profile_multi_source_multi_wfs_runtime.jl requires AMDGPU.jl for backend=amdgpu")
        AMDGPU.functional() || error("profile_multi_source_multi_wfs_runtime.jl requires a functional ROCm installation and GPU")
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

function _timed_stats!(f!::F; warmup::Int=5, samples::Int=20) where {F<:Function}
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
    return mean(timings), sort(timings)[clamp(round(Int, 0.95 * samples), 1, samples)]
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
    sim = AOSimulation(tel, atm, src, dm, wfs)
    rt = ClosedLoopRuntime(sim, _build_recon(dm, wfs, tel, src); rng=MersenneTwister(seed), wfs_detector=det)
    return rt
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
    sh_mean_ns, sh_p95_ns = _timed_stats!(() -> begin
        measure!(sh, tel, ast, det; rng=rng)
        _sync_backend!(backend_tag, sh.state.slopes)
    end; warmup=resolved_warmup, samples=resolved_samples)

    pyr = PyramidWFS(tel; n_subap=cfg.asterism_n_subap, modulation=T(1.0), mode=Diffractive(), T=T, backend=BackendArray)
    pyr_mean_ns, pyr_p95_ns = _timed_stats!(() -> begin
        measure!(pyr, tel, ast, det; rng=rng)
        _sync_backend!(backend_tag, pyr.state.slopes)
    end; warmup=resolved_warmup, samples=resolved_samples)

    bio = BioEdgeWFS(tel; n_subap=cfg.asterism_n_subap, modulation=T(1.0), mode=Diffractive(), T=T, backend=BackendArray)
    bio_mean_ns, bio_p95_ns = _timed_stats!(() -> begin
        measure!(bio, tel, ast, det; rng=rng)
        _sync_backend!(backend_tag, bio.state.slopes)
    end; warmup=resolved_warmup, samples=resolved_samples)

    compatible_runtimes = [_build_runtime_branch(T, BackendArray,
        cfg.runtime_resolution, cfg.runtime_n_subap, cfg.runtime_dm_n_act, 10 + i; sensor=:sh)
        for i in 1:cfg.branch_count]
    mixed_runtimes = vcat(
        [_build_runtime_branch(T, BackendArray,
            cfg.runtime_resolution, cfg.runtime_n_subap, cfg.runtime_dm_n_act, 30 + i; sensor=:sh)
            for i in 1:max(cfg.branch_count - 1, 1)],
        [_build_runtime_branch(T, BackendArray,
            cfg.runtime_resolution, cfg.runtime_n_subap, cfg.runtime_dm_n_act, 40; sensor=:pyr)],
    )
    for runtime in compatible_runtimes
        step!(runtime)
    end
    for runtime in mixed_runtimes
        step!(runtime)
    end
    compatible_interfaces = map(SimulationInterface, compatible_runtimes)
    mixed_interfaces = map(SimulationInterface, mixed_runtimes)
    composite_compatible = CompositeSimulationInterface(compatible_interfaces...)
    composite_mixed = CompositeSimulationInterface(mixed_interfaces...)

    comp_mean_ns, comp_p95_ns = _timed_stats!(() -> begin
        step!(composite_compatible)
        _sync_backend!(backend_tag, simulation_command(composite_compatible))
    end; warmup=resolved_warmup, samples=resolved_samples)
    mixed_mean_ns, mixed_p95_ns = _timed_stats!(() -> begin
        step!(composite_mixed)
        _sync_backend!(backend_tag, simulation_command(composite_mixed))
    end; warmup=resolved_warmup, samples=resolved_samples)

    println("multi_source_multi_wfs_profile")
    println("  backend: ", label)
    println("  scale: ", cfg.scale)
    println("  sh_asterism_mean_ns: ", sh_mean_ns)
    println("  sh_asterism_p95_ns: ", sh_p95_ns)
    println("  pyramid_asterism_mean_ns: ", pyr_mean_ns)
    println("  pyramid_asterism_p95_ns: ", pyr_p95_ns)
    println("  bioedge_asterism_mean_ns: ", bio_mean_ns)
    println("  bioedge_asterism_p95_ns: ", bio_p95_ns)
    println("  composite_compatible_mean_ns: ", comp_mean_ns)
    println("  composite_compatible_p95_ns: ", comp_p95_ns)
    println("  composite_mixed_mean_ns: ", mixed_mean_ns)
    println("  composite_mixed_p95_ns: ", mixed_p95_ns)
    println("  source_count: ", length(cfg.source_coords))
    println("  asterism_resolution: ", cfg.asterism_resolution)
    println("  asterism_n_subap: ", cfg.asterism_n_subap)
    println("  runtime_resolution: ", cfg.runtime_resolution)
    println("  runtime_n_subap: ", cfg.runtime_n_subap)
    println("  runtime_branch_count: ", cfg.branch_count)
    return nothing
end

run_profile(; backend_name=_backend_arg, scale_name=_scale_arg)

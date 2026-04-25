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
        CUDA.functional() || error("profile_platform_runtime.jl requires a functional CUDA driver/device")
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.CUDABackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(AdaptiveOpticsSim.CUDABackendTag)
        backend === nothing && error("CUDA backend array type is unavailable")
        return backend, AdaptiveOpticsSim.CUDABackendTag, "cuda"
    elseif lowered == "amdgpu"
        isdefined(Main, :AMDGPU) || (@eval import AMDGPU)
        AMDGPU.functional() || error("profile_platform_runtime.jl requires a functional ROCm installation and GPU")
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
    sorted = sort(timings)
    p95 = sorted[clamp(round(Int, 0.95 * samples), 1, samples)]
    return mean(timings), p95, Int(round(mean(alloc_bytes)))
end

function _resolve_scale(name::AbstractString)
    lowered = lowercase(name)
    lowered in ("compact", "medium", "representative") ||
        error("unsupported scale '$name'; use compact, medium, or representative")
    return lowered
end

function _platform_scale_config(name::AbstractString)
    scale = _resolve_scale(name)
    if scale == "compact"
        return (
            scale=scale,
            single_resolution=24,
            single_wfs_samples=6,
            single_dm_n_act=6,
            grouped_resolution=16,
            grouped_wfs_samples=4,
            grouped_dm_n_act=4,
            grouped_branch_count=2,
            warmup=5,
            samples=20,
        )
    elseif scale == "medium"
        return (
            scale=scale,
            single_resolution=48,
            single_wfs_samples=12,
            single_dm_n_act=12,
            grouped_resolution=32,
            grouped_wfs_samples=8,
            grouped_dm_n_act=8,
            grouped_branch_count=3,
            warmup=3,
            samples=12,
        )
    end
    return (
        scale=scale,
        single_resolution=72,
        single_wfs_samples=18,
        single_dm_n_act=18,
        grouped_resolution=48,
        grouped_wfs_samples=12,
        grouped_dm_n_act=12,
        grouped_branch_count=4,
        warmup=2,
        samples=6,
    )
end

function _build_recon(dm, wfs, tel, src)
    imat = interaction_matrix(dm, wfs, tel, src; amplitude=eltype(dm.state.coefs)(0.05))
    return ModalReconstructor(imat; gain=eltype(dm.state.coefs)(0.5))
end

function _runtime_components(::Type{T}, BackendArray, resolution::Int, wfs_samples::Int, n_act::Int;
    sensor::Symbol=:sh) where {T<:AbstractFloat}
    tel = Telescope(resolution=resolution, diameter=T(8.0), sampling_time=T(1e-3),
        central_obstruction=T(0.0), T=T, backend=BackendArray)
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = KolmogorovAtmosphere(tel; r0=T(0.2), L0=T(25.0), T=T, backend=BackendArray)
    dm = DeformableMirror(tel; n_act=n_act, influence_width=T(0.3), T=T, backend=BackendArray)
    wfs = if sensor == :sh
        ShackHartmann(tel; n_lenslets=wfs_samples, mode=Diffractive(), T=T, backend=BackendArray)
    elseif sensor == :pyr
        PyramidWFS(tel; pupil_samples=wfs_samples, modulation=T(1.0), mode=Diffractive(), T=T, backend=BackendArray)
    else
        error("unsupported sensor '$sensor'")
    end
    sim = AdaptiveOpticsSim.AOSimulation(tel, atm, src, dm, wfs)
    recon = _build_recon(dm, wfs, tel, src)
    return sim, recon
end

function _build_branch(::Type{T}, BackendArray, resolution::Int, wfs_samples::Int, n_act::Int,
    seed::Integer; sensor::Symbol=:sh, science::Bool=false, label::Symbol=Symbol("branch_", seed)) where {T<:AbstractFloat}
    sim, recon = _runtime_components(T, BackendArray, resolution, wfs_samples, n_act; sensor=sensor)
    wfs_det = Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1, T=T, backend=BackendArray)
    science_det = science ?
        Detector(noise=NoiseNone(), integration_time=T(1.0), qe=T(1.0), binning=1, T=T, backend=BackendArray) :
        nothing
    return ClosedLoopBranchConfig(
        label,
        sim,
        recon;
        wfs_detector=wfs_det,
        science_detector=science_det,
        rng=MersenneTwister(seed),
    )
end

function _frame_shape(value)
    isnothing(value) && return "none"
    return string(size(value))
end

function run_profile(; backend_name::AbstractString="cpu", scale_name::AbstractString="compact",
    T::Type{<:AbstractFloat}=Float32, samples::Union{Int,Nothing}=nothing,
    warmup::Union{Int,Nothing}=nothing)
    BackendArray, backend_tag, label = _resolve_backend(backend_name)
    cfg = _platform_scale_config(scale_name)
    resolved_samples = something(samples, cfg.samples)
    resolved_warmup = something(warmup, cfg.warmup)

    single_cfg = SinglePlatformConfig(
        name=:single_platform_runtime,
        branch_label=:main,
        products=RuntimeProductRequirements(slopes=true, wfs_pixels=true, science_pixels=true),
    )
    single_branch = _build_branch(T, BackendArray,
        cfg.single_resolution, cfg.single_wfs_samples, cfg.single_dm_n_act, 11;
        sensor=:pyr, science=true, label=:main)
    build_t0 = time_ns()
    single = build_platform_scenario(single_cfg, single_branch)
    prepare!(single)
    single_build_time_ns = time_ns() - build_t0
    step!(single)
    _sync_backend!(backend_tag, simulation_command(single))
    single_mean_ns, single_p95_ns, single_alloc_bytes = _timed_stats!(() -> begin
        step!(single)
        _sync_backend!(backend_tag, simulation_command(single))
    end; warmup=resolved_warmup, samples=resolved_samples)

    compatible_branches = ntuple(i ->
        _build_branch(T, BackendArray,
            cfg.grouped_resolution, cfg.grouped_wfs_samples, cfg.grouped_dm_n_act, 100 + i;
            sensor=:sh, label=Symbol("compatible_", i)),
        cfg.grouped_branch_count,
    )
    mixed_branches = ntuple(i ->
        _build_branch(T, BackendArray,
            cfg.grouped_resolution, cfg.grouped_wfs_samples, cfg.grouped_dm_n_act, 200 + i;
            sensor=(i == cfg.grouped_branch_count ? :pyr : :sh), label=Symbol("mixed_", i)),
        cfg.grouped_branch_count,
    )
    compatible_cfg = GroupedPlatformConfig(
        compatible_branches...;
        name=:compatible_platform_runtime,
        products=GroupedRuntimeProductRequirements(wfs_frames=true, science_frames=false, wfs_stack=true, science_stack=false),
    )
    mixed_cfg = GroupedPlatformConfig(
        mixed_branches...;
        name=:mixed_platform_runtime,
        products=GroupedRuntimeProductRequirements(wfs_frames=true, science_frames=false, wfs_stack=false, science_stack=false),
    )

    build_t0 = time_ns()
    compatible = build_platform_scenario(compatible_cfg, compatible_branches...)
    prepare!(compatible)
    compatible_build_time_ns = time_ns() - build_t0
    build_t0 = time_ns()
    mixed = build_platform_scenario(mixed_cfg, mixed_branches...)
    prepare!(mixed)
    mixed_build_time_ns = time_ns() - build_t0

    step!(compatible)
    step!(mixed)
    _sync_backend!(backend_tag, simulation_command(compatible))
    _sync_backend!(backend_tag, simulation_command(mixed))

    compatible_mean_ns, compatible_p95_ns, compatible_alloc_bytes = _timed_stats!(() -> begin
        step!(compatible)
        _sync_backend!(backend_tag, simulation_command(compatible))
    end; warmup=resolved_warmup, samples=resolved_samples)
    mixed_mean_ns, mixed_p95_ns, mixed_alloc_bytes = _timed_stats!(() -> begin
        step!(mixed)
        _sync_backend!(backend_tag, simulation_command(mixed))
    end; warmup=resolved_warmup, samples=resolved_samples)

    compatible_stack = simulation_grouped_wfs_stack(compatible)
    mixed_stack = simulation_grouped_wfs_stack(mixed)
    single_science = simulation_science_frame(single)

    return (
        backend=label,
        scale=cfg.scale,
        single_wfs_family="pyramid",
        single_build_time_ns=single_build_time_ns,
        single_mean_ns=single_mean_ns,
        single_p95_ns=single_p95_ns,
        single_frame_rate_hz=1e9 / single_mean_ns,
        single_alloc_bytes=single_alloc_bytes,
        single_science_frame_shape=_frame_shape(single_science),
        compatible_build_time_ns=compatible_build_time_ns,
        compatible_mean_ns=compatible_mean_ns,
        compatible_p95_ns=compatible_p95_ns,
        compatible_frame_rate_hz=1e9 / compatible_mean_ns,
        compatible_alloc_bytes=compatible_alloc_bytes,
        compatible_branch_count=cfg.grouped_branch_count,
        compatible_grouped_wfs_stack=!isnothing(compatible_stack),
        compatible_grouped_wfs_stack_shape=_frame_shape(compatible_stack),
        mixed_build_time_ns=mixed_build_time_ns,
        mixed_mean_ns=mixed_mean_ns,
        mixed_p95_ns=mixed_p95_ns,
        mixed_frame_rate_hz=1e9 / mixed_mean_ns,
        mixed_alloc_bytes=mixed_alloc_bytes,
        mixed_branch_count=cfg.grouped_branch_count,
        mixed_grouped_wfs_stack=!isnothing(mixed_stack),
        mixed_grouped_wfs_stack_shape=_frame_shape(mixed_stack),
    )
end

function print_profile(result)
    println("platform_runtime_profile")
    for key in (
        :backend,
        :scale,
        :single_wfs_family,
        :single_build_time_ns,
        :single_mean_ns,
        :single_p95_ns,
        :single_frame_rate_hz,
        :single_alloc_bytes,
        :single_science_frame_shape,
        :compatible_build_time_ns,
        :compatible_mean_ns,
        :compatible_p95_ns,
        :compatible_frame_rate_hz,
        :compatible_alloc_bytes,
        :compatible_branch_count,
        :compatible_grouped_wfs_stack,
        :compatible_grouped_wfs_stack_shape,
        :mixed_build_time_ns,
        :mixed_mean_ns,
        :mixed_p95_ns,
        :mixed_frame_rate_hz,
        :mixed_alloc_bytes,
        :mixed_branch_count,
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

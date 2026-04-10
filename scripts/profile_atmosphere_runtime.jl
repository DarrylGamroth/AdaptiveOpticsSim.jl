using AdaptiveOpticsSim
using Random
using Statistics

const _backend_arg = isempty(ARGS) ? "cpu" : lowercase(ARGS[1])
const _scale_arg = length(ARGS) >= 2 ? lowercase(ARGS[2]) : "medium"
const _source_arg = length(ARGS) >= 3 ? lowercase(ARGS[3]) : "onaxis"

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
        isdefined(Main, :CUDA) || error("profile_atmosphere_runtime.jl requires CUDA.jl for backend=cuda")
        CUDA.functional() || error("profile_atmosphere_runtime.jl requires a functional CUDA driver/device")
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.CUDABackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(AdaptiveOpticsSim.CUDABackendTag)
        backend === nothing && error("CUDA backend array type is unavailable")
        return backend, AdaptiveOpticsSim.CUDABackendTag, "cuda"
    elseif lowered == "amdgpu"
        isdefined(Main, :AMDGPU) || error("profile_atmosphere_runtime.jl requires AMDGPU.jl for backend=amdgpu")
        AMDGPU.functional() || error("profile_atmosphere_runtime.jl requires a functional ROCm installation and GPU")
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.AMDGPUBackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(AdaptiveOpticsSim.AMDGPUBackendTag)
        backend === nothing && error("AMDGPU backend array type is unavailable")
        return backend, AdaptiveOpticsSim.AMDGPUBackendTag, "amdgpu"
    end
    error("unsupported backend '$name'; use cpu, cuda, or amdgpu")
end

function _sync_array!(::Nothing, _)
    return nothing
end

function _sync_array!(::Type{B}, A) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(A))
    return nothing
end

function _resolve_scale(name::AbstractString)
    lowered = lowercase(name)
    lowered in ("compact", "medium", "representative") ||
        error("unsupported scale '$name'; use compact, medium, or representative")
    if lowered == "compact"
        return (
            scale=lowered,
            resolution=32,
            diameter=8.0f0,
            r0=0.2f0,
            L0=25.0f0,
            fractional_cn2=Float32[0.7, 0.3],
            wind_speed=Float32[8.0, 4.0],
            wind_direction=Float32[0.0, 90.0],
            altitude=Float32[0.0, 5000.0],
            screen_resolution=65,
            stencil_size=67,
            warmup=5,
            samples=20,
        )
    elseif lowered == "medium"
        return (
            scale=lowered,
            resolution=64,
            diameter=8.0f0,
            r0=0.2f0,
            L0=25.0f0,
            fractional_cn2=Float32[0.5, 0.3, 0.2],
            wind_speed=Float32[12.0, 8.0, 4.0],
            wind_direction=Float32[0.0, 60.0, 120.0],
            altitude=Float32[0.0, 5000.0, 10_000.0],
            screen_resolution=129,
            stencil_size=131,
            warmup=4,
            samples=12,
        )
    end
    return (
        scale=lowered,
        resolution=96,
        diameter=8.0f0,
        r0=0.2f0,
        L0=25.0f0,
        fractional_cn2=Float32[0.45, 0.3, 0.15, 0.1],
        wind_speed=Float32[14.0, 10.0, 7.0, 4.0],
        wind_direction=Float32[0.0, 45.0, 90.0, 135.0],
        altitude=Float32[0.0, 4000.0, 9000.0, 14_000.0],
        screen_resolution=193,
        stencil_size=257,
        warmup=3,
        samples=8,
    )
end

function _resolve_source(name::AbstractString, T::Type{<:AbstractFloat})
    lowered = lowercase(name)
    if lowered == "onaxis"
        return Source(band=:I, magnitude=0.0, T=T), "onaxis"
    elseif lowered == "offaxis"
        return Source(band=:I, magnitude=0.0, coordinates=(45.0, 30.0), T=T), "offaxis"
    elseif lowered == "lgs"
        return LGSSource(magnitude=0.0, coordinates=(45.0, 30.0), altitude=90_000.0, T=T), "lgs"
    end
    error("unsupported source '$name'; use onaxis, offaxis, or lgs")
end

function _allocated_bytes(f!::F; warmup::Int=2, gc_before::Bool=true) where {F<:Function}
    for _ in 1:warmup
        f!()
    end
    gc_before && GC.gc()
    return @allocated f!()
end

function _step_with_sync!(atm, tel, src, backend_tag, rng)
    advance!(atm, tel; rng=rng)
    propagate!(atm, tel, src)
    _sync_array!(backend_tag, atm.state.opd)
    _sync_array!(backend_tag, tel.state.opd)
    return tel.state.opd
end

function _benchmark_model(model_name::Symbol, tel, src, backend_tag, BackendArray, cfg, T)
    rng_build = MersenneTwister(1)
    rng_step = MersenneTwister(2)
    t0 = time_ns()
    atm = if model_name === :finite
        MultiLayerAtmosphere(tel;
            r0=cfg.r0,
            L0=cfg.L0,
            fractional_cn2=cfg.fractional_cn2,
            wind_speed=cfg.wind_speed,
            wind_direction=cfg.wind_direction,
            altitude=cfg.altitude,
            T=T,
            backend=BackendArray,
        )
    elseif model_name === :infinite
        InfiniteMultiLayerAtmosphere(tel;
            r0=cfg.r0,
            L0=cfg.L0,
            fractional_cn2=cfg.fractional_cn2,
            wind_speed=cfg.wind_speed,
            wind_direction=cfg.wind_direction,
            altitude=cfg.altitude,
            screen_resolution=cfg.screen_resolution,
            stencil_size=cfg.stencil_size,
            T=T,
            backend=BackendArray,
        )
    else
        error("unsupported atmosphere model '$model_name'")
    end
    _step_with_sync!(atm, tel, src, backend_tag, rng_build)
    build_time_ns = time_ns() - t0
    timing = runtime_timing(() -> _step_with_sync!(atm, tel, src, backend_tag, rng_step);
        warmup=cfg.warmup, samples=cfg.samples, gc_before=false)
    alloc_bytes = _allocated_bytes(() -> _step_with_sync!(atm, tel, src, backend_tag, rng_step);
        warmup=cfg.warmup, gc_before=false)
    return (
        build_time_ns=build_time_ns,
        step_mean_ns=timing.mean_ns,
        step_p95_ns=timing.p95_ns,
        frame_rate_hz=1.0e9 / timing.mean_ns,
        alloc_bytes=alloc_bytes,
        sync_count_per_sample=backend_tag === nothing ? 0 : 2,
        opd_std=std(vec(Array(atm.state.opd))),
        screen_shape=model_name === :finite ? size(atm.layers[1].generator.state.opd) : size(atm.layers[1].screen.state.screen),
    )
end

function run_profile(; backend_name::AbstractString="cpu", scale_name::AbstractString="medium",
    source_name::AbstractString="onaxis")
    BackendArray, backend_tag, backend_label = _resolve_backend(backend_name)
    cfg = _resolve_scale(scale_name)
    T = Float32
    tel = Telescope(resolution=cfg.resolution, diameter=cfg.diameter, sampling_time=1.0f-3,
        central_obstruction=0.0f0, T=T, backend=BackendArray)
    src, source_label = _resolve_source(source_name, T)
    finite = _benchmark_model(:finite, tel, src, backend_tag, BackendArray, cfg, T)
    infinite = _benchmark_model(:infinite, tel, src, backend_tag, BackendArray, cfg, T)

    println("atmosphere_runtime_profile")
    println("  backend: ", backend_label)
    println("  scale: ", cfg.scale)
    println("  source: ", source_label)
    println("  pupil_resolution: ", cfg.resolution)
    println("  n_layers: ", length(cfg.fractional_cn2))
    println("  finite_screen_shape: ", finite.screen_shape)
    println("  infinite_screen_shape: ", infinite.screen_shape)
    println("  finite_build_time_ns: ", finite.build_time_ns)
    println("  infinite_build_time_ns: ", infinite.build_time_ns)
    println("  finite_runtime_step_mean_ns: ", finite.step_mean_ns)
    println("  infinite_runtime_step_mean_ns: ", infinite.step_mean_ns)
    println("  finite_runtime_step_p95_ns: ", finite.step_p95_ns)
    println("  infinite_runtime_step_p95_ns: ", infinite.step_p95_ns)
    println("  finite_frame_rate_hz: ", finite.frame_rate_hz)
    println("  infinite_frame_rate_hz: ", infinite.frame_rate_hz)
    println("  finite_runtime_alloc_bytes: ", finite.alloc_bytes)
    println("  infinite_runtime_alloc_bytes: ", infinite.alloc_bytes)
    println("  finite_sync_count_per_sample: ", finite.sync_count_per_sample)
    println("  infinite_sync_count_per_sample: ", infinite.sync_count_per_sample)
    println("  finite_opd_std: ", finite.opd_std)
    println("  infinite_opd_std: ", infinite.opd_std)
    println("  infinite_to_finite_step_ratio: ", infinite.step_mean_ns / finite.step_mean_ns)
end

run_profile(; backend_name=_backend_arg, scale_name=_scale_arg, source_name=_source_arg)

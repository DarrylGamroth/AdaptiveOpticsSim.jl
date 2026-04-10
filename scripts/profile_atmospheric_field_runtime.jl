using AdaptiveOpticsSim
using Random
using Statistics

const _backend_arg = isempty(ARGS) ? "cpu" : lowercase(ARGS[1])
const _mode_arg = length(ARGS) >= 2 ? lowercase(ARGS[2]) : "geometric"
const _atmo_arg = length(ARGS) >= 3 ? lowercase(ARGS[3]) : "finite"
const _scale_arg = length(ARGS) >= 4 ? lowercase(ARGS[4]) : "medium"

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
        isdefined(Main, :CUDA) || error("profile_atmospheric_field_runtime.jl requires CUDA.jl for backend=cuda")
        CUDA.functional() || error("profile_atmospheric_field_runtime.jl requires a functional CUDA driver/device")
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.CUDABackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(AdaptiveOpticsSim.CUDABackendTag)
        backend === nothing && error("CUDA backend array type is unavailable")
        return backend, AdaptiveOpticsSim.CUDABackendTag, "cuda"
    elseif lowered == "amdgpu"
        isdefined(Main, :AMDGPU) || error("profile_atmospheric_field_runtime.jl requires AMDGPU.jl for backend=amdgpu")
        AMDGPU.functional() || error("profile_atmospheric_field_runtime.jl requires a functional ROCm installation and GPU")
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.AMDGPUBackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(AdaptiveOpticsSim.AMDGPUBackendTag)
        backend === nothing && error("AMDGPU backend array type is unavailable")
        return backend, AdaptiveOpticsSim.AMDGPUBackendTag, "amdgpu"
    end
    error("unsupported backend '$name'; use cpu, cuda, or amdgpu")
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

function _sync_array!(::Nothing, _)
    return nothing
end

function _sync_array!(::Type{B}, A) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(A))
    return nothing
end

function _make_atmosphere(kind::Symbol, tel, cfg, T, BackendArray)
    if kind === :finite
        return MultiLayerAtmosphere(tel;
            r0=cfg.r0,
            L0=cfg.L0,
            fractional_cn2=cfg.fractional_cn2,
            wind_speed=cfg.wind_speed,
            wind_direction=cfg.wind_direction,
            altitude=cfg.altitude,
            T=T,
            backend=BackendArray,
        )
    elseif kind === :infinite
        return InfiniteMultiLayerAtmosphere(tel;
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
    end
    error("unsupported atmosphere kind '$kind'")
end

function _allocated_bytes(f!::F; warmup::Int=2, gc_before::Bool=true) where {F<:Function}
    for _ in 1:warmup
        f!()
    end
    gc_before && GC.gc()
    return @allocated f!()
end

function _profile_field_path(mode::Symbol, atmo_kind::Symbol, backend_name::AbstractString, scale_name::AbstractString)
    BackendArray, backend_tag, backend_label = _resolve_backend(backend_name)
    cfg = _resolve_scale(scale_name)
    T = Float32
    tel = Telescope(resolution=cfg.resolution, diameter=cfg.diameter, sampling_time=1.0f-3,
        central_obstruction=0.0f0, T=T, backend=BackendArray)
    src = Source(band=:I, magnitude=0.0, T=T)
    atm = _make_atmosphere(atmo_kind, tel, cfg, T, BackendArray)
    rng = MersenneTwister(3)
    advance!(atm, tel; rng=rng)

    prop = AtmosphericFieldPropagation(atm, tel, src;
        model=mode === :geometric ? GeometricAtmosphericPropagation(T=T) : LayeredFresnelAtmosphericPropagation(T=T),
        zero_padding=2,
        T=T)
    wfs = CurvatureWFS(tel; n_subap=8, T=T, backend=BackendArray)

    step! = if mode === :curvature
        () -> begin
            advance!(atm, tel; rng=rng)
            measure!(wfs, tel, src, atm)
            _sync_array!(backend_tag, wfs.state.slopes)
            return wfs.state.slopes
        end
    else
        () -> begin
            advance!(atm, tel; rng=rng)
            field = propagate_atmosphere_field!(prop, atm, tel, src)
            _sync_array!(backend_tag, field.state.field)
            return field.state.field
        end
    end

    build_t0 = time_ns()
    step!()
    build_time_ns = time_ns() - build_t0
    timing = runtime_timing(step!; warmup=cfg.warmup, samples=cfg.samples, gc_before=false)
    alloc_bytes = _allocated_bytes(step!; warmup=cfg.warmup, gc_before=false)

    println("atmospheric_field_runtime_profile")
    println("  backend: ", backend_label)
    println("  mode: ", String(mode))
    println("  atmosphere: ", String(atmo_kind))
    println("  scale: ", cfg.scale)
    println("  pupil_resolution: ", cfg.resolution)
    println("  n_layers: ", length(cfg.fractional_cn2))
    println("  build_time_ns: ", build_time_ns)
    println("  step_mean_ns: ", timing.mean_ns)
    println("  step_p95_ns: ", timing.p95_ns)
    println("  frame_rate_hz: ", 1.0e9 / timing.mean_ns)
    println("  alloc_bytes: ", alloc_bytes)
    println("  sync_count_per_sample: ", backend_tag === nothing ? 0 : 1)
end

mode = Symbol(_mode_arg)
atmo_kind = Symbol(_atmo_arg)
mode in (:geometric, :fresnel, :curvature) || error("mode must be geometric, fresnel, or curvature")
atmo_kind in (:finite, :infinite) || error("atmosphere must be finite or infinite")

_profile_field_path(mode, atmo_kind, _backend_arg, _scale_arg)

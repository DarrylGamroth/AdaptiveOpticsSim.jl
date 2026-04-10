using AdaptiveOpticsSim
using Random

const _backend_arg = isempty(ARGS) ? "cpu" : lowercase(ARGS[1])

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
        isdefined(Main, :CUDA) || error("profile_mixed_sh_asterism_runtime.jl requires CUDA.jl for backend=cuda")
        CUDA.functional() || error("profile_mixed_sh_asterism_runtime.jl requires a functional CUDA driver/device")
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.CUDABackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(AdaptiveOpticsSim.CUDABackendTag)
        backend === nothing && error("CUDA backend array type is unavailable")
        return backend, AdaptiveOpticsSim.CUDABackendTag, "cuda"
    elseif lowered == "amdgpu"
        isdefined(Main, :AMDGPU) || error("profile_mixed_sh_asterism_runtime.jl requires AMDGPU.jl for backend=amdgpu")
        AMDGPU.functional() || error("profile_mixed_sh_asterism_runtime.jl requires a functional ROCm installation and GPU")
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.AMDGPUBackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(AdaptiveOpticsSim.AMDGPUBackendTag)
        backend === nothing && error("AMDGPU backend array type is unavailable")
        return backend, AdaptiveOpticsSim.AMDGPUBackendTag, "amdgpu"
    end
    error("unsupported backend '$name'; use cpu, cuda, or amdgpu")
end

_sync_wfs!(::Nothing, _) = nothing

function _sync_wfs!(::Type{B}, wfs) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(wfs.state.slopes))
    return nothing
end

function _na_profile(T::Type{<:AbstractFloat})
    return T[
        89_500 90_000 90_500 91_000 91_500;
        0.10   0.25   0.30   0.22   0.13
    ]
end

function run_profile(; backend_name::AbstractString="cpu", samples::Int=20, warmup::Int=5)
    BackendArray, backend_tag, backend_label = _resolve_backend(backend_name)
    T = Float32

    tel = Telescope(
        resolution=112,
        diameter=8.2,
        sampling_time=1e-3,
        central_obstruction=0.30,
        T=T,
        backend=BackendArray,
    )
    ngs = Source(wavelength=T(589e-9), magnitude=0.0, coordinates=(0.0, 0.0), T=T)
    lgs = LGSSource(
        magnitude=8.0,
        wavelength=589e-9,
        altitude=90000.0,
        laser_coordinates=(5.0, 0.0),
        na_profile=_na_profile(T),
        fwhm_spot_up=1.0,
        T=T,
    )
    ast = Asterism([ngs, lgs])
    wfs = ShackHartmann(tel; n_subap=14, mode=Diffractive(), T=T, backend=BackendArray)
    det = Detector(noise=NoiseNone(), integration_time=T(1e-3), qe=T(1), binning=1, T=T, backend=BackendArray)

    rng = MersenneTwister(1)
    AdaptiveOpticsSim.randn_backend!(rng, tel.state.opd)
    tel.state.opd .*= T(5e-8)

    t0 = time_ns()
    measure!(wfs, tel, ast, det; rng=rng)
    _sync_wfs!(backend_tag, wfs)
    build_time_ns = time_ns() - t0

    timing = runtime_timing(() -> begin
        measure!(wfs, tel, ast, det; rng=rng)
        _sync_wfs!(backend_tag, wfs)
    end; warmup=warmup, samples=samples, gc_before=false)

    println("mixed_sh_asterism_runtime_profile")
    println("  backend: ", backend_label)
    println("  build_time_ns: ", build_time_ns)
    println("  measure_mean_ns: ", timing.mean_ns)
    println("  measure_p95_ns: ", timing.p95_ns)
    println("  frame_rate_hz: ", 1.0e9 / timing.mean_ns)
    println("  spot_cube_shape: ", size(wfs.state.spot_cube))
    return nothing
end

run_profile(; backend_name=_backend_arg)

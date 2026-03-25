using AdaptiveOpticsSim
using Random

include(joinpath(dirname(@__DIR__), "examples", "support", "subaru_ao188_simulation.jl"))
using .SubaruAO188Simulation

const _backend_arg = isempty(ARGS) ? "cpu" : lowercase(ARGS[1])
const _profile_arg = length(ARGS) >= 2 ? lowercase(ARGS[2]) : "none"

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
        isdefined(Main, :CUDA) || error("profile_lgs_sh_runtime.jl requires CUDA.jl for backend=cuda")
        CUDA.functional() || error("profile_lgs_sh_runtime.jl requires a functional CUDA driver/device")
        AdaptiveOpticsSim.disable_scalar_backend!(CUDABackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(CUDABackendTag)
        backend === nothing && error("CUDA backend array type is unavailable")
        return backend, CUDABackendTag, "cuda"
    elseif lowered == "amdgpu"
        isdefined(Main, :AMDGPU) || error("profile_lgs_sh_runtime.jl requires AMDGPU.jl for backend=amdgpu")
        AMDGPU.functional() || error("profile_lgs_sh_runtime.jl requires a functional ROCm installation and GPU")
        AdaptiveOpticsSim.disable_scalar_backend!(AMDGPUBackendTag)
        backend = AdaptiveOpticsSim.gpu_backend_array_type(AMDGPUBackendTag)
        backend === nothing && error("AMDGPU backend array type is unavailable")
        return backend, AMDGPUBackendTag, "amdgpu"
    end
    error("unsupported backend '$name'; use cpu, cuda, or amdgpu")
end

_sync_wfs!(::Nothing, wfs) = nothing

function _sync_wfs!(::Type{B}, wfs) where {B<:GPUBackendTag}
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(wfs.state.slopes))
    return nothing
end

function _na_profile(T::Type{<:AbstractFloat})
    return T[
        88000 89500 90500 91500 93000;
        0.10  0.25  0.30  0.25  0.10;
    ]
end

function _resolve_lgs(profile_name::AbstractString, T::Type{<:AbstractFloat})
    lowered = lowercase(profile_name)
    if lowered == "none"
        src = LGSSource(
            magnitude=8.0,
            wavelength=589e-9,
            altitude=90000.0,
            elongation_factor=1.6,
            laser_coordinates=(5.0, 0.0),
            T=T,
        )
        return src, "none"
    elseif lowered == "na" || lowered == "na_profile"
        src = LGSSource(
            magnitude=8.0,
            wavelength=589e-9,
            altitude=90000.0,
            elongation_factor=1.6,
            laser_coordinates=(5.0, 0.0),
            na_profile=_na_profile(T),
            fwhm_spot_up=1.0,
            T=T,
        )
        return src, "na_profile"
    end
    error("unsupported LGS profile '$profile_name'; use none or na")
end

function run_profile(; backend_name::AbstractString="cpu", profile_name::AbstractString="none",
    samples::Int=20, warmup::Int=5)
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
    src, lgs_label = _resolve_lgs(profile_name, T)
    wfs = ShackHartmann(tel; n_subap=14, mode=Diffractive(), T=T, backend=BackendArray)
    det = AdaptiveOpticsSim.detector_from_config(
        AO188WFSDetectorConfig(
            T=T,
            integration_time=1e-3,
            qe=0.9,
            psf_sampling=1,
            binning=1,
            gain=1.0,
            dark_current=0.0,
            noise=NoisePhotonReadout(0.3),
            sensor=CCDSensor(),
        );
        backend=BackendArray,
    )

    rng = MersenneTwister(1)
    opd_scale = T(5e-8)
    AdaptiveOpticsSim.randn_backend!(rng, tel.state.opd)
    tel.state.opd .*= opd_scale

    t0 = time_ns()
    measure!(wfs, tel, src, det; rng=rng)
    _sync_wfs!(backend_tag, wfs)
    build_time_ns = time_ns() - t0

    timing = runtime_timing(() -> begin
        measure!(wfs, tel, src, det; rng=rng)
        _sync_wfs!(backend_tag, wfs)
    end; warmup=warmup, samples=samples, gc_before=false)

    println("lgs_sh_runtime_profile")
    println("  backend: ", backend_label)
    println("  lgs_profile: ", lgs_label)
    println("  build_time_ns: ", build_time_ns)
    println("  measure_mean_ns: ", timing.mean_ns)
    println("  measure_p95_ns: ", timing.p95_ns)
    println("  frame_rate_hz: ", 1.0e9 / timing.mean_ns)
    println("  spot_cube_shape: ", size(wfs.state.spot_cube))
    println("  has_na_profile: ", AdaptiveOpticsSim.lgs_has_profile(src))
    return nothing
end

run_profile(; backend_name=_backend_arg, profile_name=_profile_arg)

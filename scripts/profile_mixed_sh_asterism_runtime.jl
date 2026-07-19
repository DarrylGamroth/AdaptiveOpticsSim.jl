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
        return CPUBackend(), nothing, "cpu"
    elseif lowered == "cuda"
        isdefined(Main, :CUDA) || error("profile_mixed_sh_asterism_runtime.jl requires CUDA.jl for backend=cuda")
        CUDA.functional() || error("profile_mixed_sh_asterism_runtime.jl requires a functional CUDA driver/device")
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.CUDABackendTag)
        return CUDABackend(), AdaptiveOpticsSim.CUDABackendTag, "cuda"
    elseif lowered == "amdgpu"
        isdefined(Main, :AMDGPU) || error("profile_mixed_sh_asterism_runtime.jl requires AMDGPU.jl for backend=amdgpu")
        AMDGPU.functional() || error("profile_mixed_sh_asterism_runtime.jl requires a functional ROCm installation and GPU")
        AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.AMDGPUBackendTag)
        return AMDGPUBackend(), AdaptiveOpticsSim.AMDGPUBackendTag, "amdgpu"
    end
    error("unsupported backend '$name'; use cpu, cuda, or amdgpu")
end

_sync_wfs!(::Nothing, _) = nothing

function _sync_wfs!(::Type{B}, wfs) where {B<:AdaptiveOpticsSim.GPUBackendTag}
    AdaptiveOpticsSim.synchronize_backend!(AdaptiveOpticsSim.execution_style(slopes(wfs)))
    return nothing
end

function _na_profile(T::Type{<:AbstractFloat})
    return T[
        89_500 90_000 90_500 91_000 91_500;
        0.10   0.25   0.30   0.22   0.13
    ]
end

function run_profile(; backend_name::AbstractString="cpu", samples::Int=20, warmup::Int=5)
    backend, backend_tag, backend_label = _resolve_backend(backend_name)
    T = Float32

    tel = Telescope(
        resolution=112,
        diameter=8.2,
        central_obstruction=0.30,
        T=T,
        backend=backend,
    )
    lgs = LGSSource(
        magnitude=8.0,
        wavelength=589e-9,
        coordinates=(0.0, 0.0),
        altitude=90000.0,
        laser_coordinates=(5.0, 0.0),
        na_profile=_na_profile(T),
        fwhm_spot_up=1.0,
        photon_irradiance=one(T),
        T=T,
    )
    second_lgs = LGSSource(
        magnitude=8.0,
        wavelength=589e-9,
        coordinates=(5.0, 90.0),
        altitude=90000.0,
        laser_coordinates=(5.0, 0.0),
        na_profile=_na_profile(T),
        fwhm_spot_up=1.0,
        photon_irradiance=T(0.75),
        T=T,
    )
    ast = Asterism([lgs, second_lgs])
    wfs = ShackHartmannWFS(tel; n_lenslets=14, mode=Diffractive(), T=T, backend=backend)
    det = Detector(noise=NoiseNone(), integration_time=T(1e-3), qe=T(1), binning=1, T=T, backend=backend)
    pupil = PupilFunction(tel; T=T, backend=backend)

    rng = runtime_rng(1)
    AdaptiveOpticsSim.randn_backend!(rng, pupil.opd)
    pupil.opd .*= T(5e-8)

    t0 = time_ns()
    measure!(wfs, pupil, ast, det; rng=rng)
    _sync_wfs!(backend_tag, wfs)
    build_time_ns = time_ns() - t0

    timing = runtime_timing(() -> begin
        measure!(wfs, pupil, ast, det; rng=rng)
        _sync_wfs!(backend_tag, wfs)
    end; warmup=warmup, samples=samples, gc_before=false)

    println("common_calibration_lgs_asterism_runtime_profile")
    println("  backend: ", backend_label)
    println("  build_time_ns: ", build_time_ns)
    println("  measure_mean_ns: ", timing.mean_ns)
    println("  measure_p95_ns: ", timing.p95_ns)
    println("  frame_rate_hz: ", 1.0e9 / timing.mean_ns)
    println("  spot_cube_shape: ", size(wfs.acquisition.spot_cube))
    return nothing
end

run_profile(; backend_name=_backend_arg)

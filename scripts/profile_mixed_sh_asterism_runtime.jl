using AdaptiveOpticsSim
using AdaptiveOpticsSim.Detectors
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.Backends
using AdaptiveOpticsSim.WavefrontSensors
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
        AdaptiveOpticsSim.Backends.disable_scalar_backend!(AdaptiveOpticsSim.Backends.CUDABackendTag)
        return CUDABackend(), AdaptiveOpticsSim.Backends.CUDABackendTag, "cuda"
    elseif lowered == "amdgpu"
        isdefined(Main, :AMDGPU) || error("profile_mixed_sh_asterism_runtime.jl requires AMDGPU.jl for backend=amdgpu")
        AMDGPU.functional() || error("profile_mixed_sh_asterism_runtime.jl requires a functional ROCm installation and GPU")
        AdaptiveOpticsSim.Backends.disable_scalar_backend!(AdaptiveOpticsSim.Backends.AMDGPUBackendTag)
        return AMDGPUBackend(), AdaptiveOpticsSim.Backends.AMDGPUBackendTag, "amdgpu"
    end
    error("unsupported backend '$name'; use cpu, cuda, or amdgpu")
end

_sync_array!(::Nothing, _) = nothing

function _sync_array!(::Type{B}, storage) where {B<:AdaptiveOpticsSim.Backends.GPUBackendTag}
    AdaptiveOpticsSim.Backends.synchronize_backend!(AdaptiveOpticsSim.Backends.execution_style(storage))
    return nothing
end

function _sodium_layer_profile(T::Type{<:AbstractFloat})
    return SodiumLayerProfile(
        T[89_500, 90_000, 90_500, 91_000, 91_500],
        T[0.10, 0.25, 0.30, 0.22, 0.13])
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
        separation_arcsec=0.0, position_angle_deg=0.0,
        altitude=90000.0,
        laser_launch_xy_m=(5.0, 0.0),
        sodium_layer_profile=_sodium_layer_profile(T),
        fwhm_spot_up=1.0,
        photon_irradiance=one(T),
        T=T,
    )
    second_lgs = LGSSource(
        magnitude=8.0,
        wavelength=589e-9,
        separation_arcsec=5.0, position_angle_deg=90.0,
        altitude=90000.0,
        laser_launch_xy_m=(5.0, 0.0),
        sodium_layer_profile=_sodium_layer_profile(T),
        fwhm_spot_up=1.0,
        photon_irradiance=T(0.75),
        T=T,
    )
    ast = Asterism([lgs, second_lgs])
    wfs = ShackHartmannWFS(tel; n_lenslets=14, T=T, backend=backend)
    det = Detector(noise=NoiseNone(), exposure_duration=T(1e-3), qe=T(1), binning=1, T=T, backend=backend)
    pupil = PupilFunction(tel; T=T, backend=backend)

    rng = runtime_rng(1)
    AdaptiveOpticsSim.Backends.randn_backend!(rng, pupil.opd)
    pupil.opd .*= T(5e-8)

    rate = shack_hartmann_rate_map(wfs, pupil, ast)
    optics = prepare_wfs_optics(shack_hartmann_optics(wfs, ast), pupil,
        rate)
    frame = similar(rate.values)
    observation = WFSObservation(frame; units=:electron_count,
        layout=:lenslet_mosaic)
    acquisition = prepare_wfs_acquisition(det, rate, observation;
        source=ast)

    step! = () -> begin
        form_wfs_optical_products!(rate, pupil, optics)
        acquire_wfs_observation!(observation, rate, acquisition, rng)
        _sync_array!(backend_tag, observation_storage(observation))
        return observation
    end

    t0 = time_ns()
    step!()
    build_time_ns = time_ns() - t0

    timing = runtime_timing(step!; warmup=warmup, samples=samples,
        gc_before=false)

    println("shack_hartmann_lgs_asterism_optical_acquisition_profile")
    println("  backend: ", backend_label)
    println("  build_time_ns: ", build_time_ns)
    println("  optical_acquisition_mean_ns: ", timing.mean_ns)
    println("  optical_acquisition_p95_ns: ", timing.p95_ns)
    println("  frame_rate_hz: ", 1.0e9 / timing.mean_ns)
    println("  photon_rate_mosaic_shape: ", size(rate.values))
    println("  observation_shape: ", size(observation_storage(observation)))
    return nothing
end

run_profile(; backend_name=_backend_arg)

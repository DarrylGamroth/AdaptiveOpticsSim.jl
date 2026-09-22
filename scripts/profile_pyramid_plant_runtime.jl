using AdaptiveOpticsSim
using AdaptiveOpticsSim.Backends
using AdaptiveOpticsSim.Detectors
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.WavefrontSensors
using Random

const PROFILE_BACKEND = lowercase(get(ENV,
    "ADAPTIVEOPTICS_PROFILE_BACKEND", "cuda"))
const PROFILE_STEPS = parse(Int, get(ENV,
    "ADAPTIVEOPTICS_PROFILE_STEPS", "1000"))

if PROFILE_BACKEND == "cuda"
    import CUDA
elseif PROFILE_BACKEND == "amdgpu"
    import AMDGPU
else
    error("ADAPTIVEOPTICS_PROFILE_BACKEND must be cuda or amdgpu")
end

function profile_backend(name::AbstractString)
    if name == "cuda"
        CUDA.functional() || error("CUDA is not functional")
        Backends.disable_scalar_backend!(Backends.CUDABackendTag)
        return Backends.CUDABackend()
    elseif name == "amdgpu"
        AMDGPU.functional() || error("AMDGPU is not functional")
        Backends.disable_scalar_backend!(Backends.AMDGPUBackendTag)
        return Backends.AMDGPUBackend()
    end
    error("ADAPTIVEOPTICS_PROFILE_BACKEND must be cuda or amdgpu")
end

function prepare_pyramid_profile(backend::AbstractArrayBackend, strategy)
    T = Float32
    resolution = 64
    telescope = Telescope(
        resolution=resolution,
        diameter=T(8),
        central_obstruction=zero(T),
        T=T,
        backend=backend,
    )
    source = Source(
        band=:custom,
        wavelength=T(750e-9),
        photon_irradiance=T(2e8),
        T=T,
    )
    pupil = PupilFunction(telescope; T=T, backend=backend)
    host_opd = Matrix{T}(undef, resolution, resolution)
    @inbounds for column in axes(host_opd, 2), row in axes(host_opd, 1)
        host_opd[row, column] = T(40e-9) * (
            T(row - 1) / T(resolution - 1) -
            T(0.4) * T(column - 1) / T(resolution - 1))
    end
    copyto!(pupil.opd, host_opd)
    sensor = PyramidWFS(
        telescope;
        pupil_samples=16,
        modulation=T(3),
        modulation_points=8,
        modulation_propagation_strategy=strategy,
        T=T,
        backend=backend,
    )
    front_end = PyramidOpticalFrontEnd(sensor, source)
    rate = pyramid_rate_map(front_end, pupil)
    optics = prepare_wfs_optics(front_end, pupil, rate)
    detector = Detector(
        noise=NoiseNone(),
        exposure_duration=T(1e-3),
        qe=one(T),
        response_model=NullFrameResponse(),
        T=T,
        backend=backend,
    )
    observation = WFSObservation(
        similar(rate.values);
        units=:electron_count,
        layout=:four_pupil_mosaic,
    )
    acquisition = prepare_wfs_acquisition(
        detector,
        rate,
        observation;
        source,
    )
    return (; pupil, rate, observation, optics, acquisition,
        rng=Xoshiro(0x50575253))
end

@inline function step_pyramid_profile!(prepared)
    form_wfs_optical_products!(
        prepared.rate,
        prepared.pupil,
        prepared.optics,
    )
    acquire_wfs_observation!(
        prepared.observation,
        prepared.rate,
        prepared.acquisition,
        prepared.rng,
    )
    Backends.synchronize_backend!(Backends.execution_style(
        observation_storage(prepared.observation)))
    return nothing
end

function run_pyramid_profile(name::AbstractString, backend, strategy)
    prepared = prepare_pyramid_profile(backend, strategy)
    for _ in 1:20
        step_pyramid_profile!(prepared)
    end
    allocation_bytes = @allocated step_pyramid_profile!(prepared)
    start_nanoseconds = time_ns()
    for _ in 1:PROFILE_STEPS
        step_pyramid_profile!(prepared)
    end
    elapsed_nanoseconds = time_ns() - start_nanoseconds
    values = Array(observation_storage(prepared.observation))
    all(isfinite, values) || error("$name produced non-finite values")
    println("pyramid_plant_profile")
    println("  backend=", PROFILE_BACKEND)
    println("  strategy=", name)
    println("  steps=", PROFILE_STEPS)
    println("  target_ready_elapsed_ns=", elapsed_nanoseconds)
    println("  target_ready_mean_ns=", elapsed_nanoseconds / PROFILE_STEPS)
    println("  warmed_host_allocation_bytes=", allocation_bytes)
    println("  output_sum=", sum(values))
    return nothing
end

backend = profile_backend(PROFILE_BACKEND)
run_pyramid_profile("pupil_tilt", backend,
    WavefrontSensors.PyramidPupilTiltStrategy())
run_pyramid_profile("shifted_mask", backend,
    WavefrontSensors.PyramidShiftedMaskStrategy())

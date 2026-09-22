"""
Profile the retained AdaptiveOpticsSim Bi-O-edge plant path:

    pupil field -> four-pupil photon-rate map -> detector observation

Set `ADAPTIVEOPTICS_PROFILE_BACKEND` to `cpu`, `cuda`, or `amdgpu` and
`ADAPTIVEOPTICS_PROFILE_STEPS` to the measured repetition count. The script
warms the prepared path and reports mean wall time and Julia heap allocation
for one complete target-ready step. CPU execution requires zero warmed heap
allocation. Direct GPU stream execution reports launch overhead; zero host
allocation is a captured CUDA/HIP Graph replay contract, measured separately.
CPU runs also print a sampling profile; GPU runs are intended to be launched
under Nsight Systems/Compute or rocprofv3 after the functional gates pass.
"""

using AdaptiveOpticsSim
using AdaptiveOpticsSim.Backends
using AdaptiveOpticsSim.Detectors
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.WavefrontSensors
using Profile
using Random

const PROFILE_BACKEND = lowercase(get(ENV,
    "ADAPTIVEOPTICS_PROFILE_BACKEND", "cpu"))
const PROFILE_STEPS = parse(Int, get(ENV,
    "ADAPTIVEOPTICS_PROFILE_STEPS", "1000"))
const PROFILE_PRINT_CPU_SAMPLES = get(ENV,
    "ADAPTIVEOPTICS_PROFILE_PRINT_CPU_SAMPLES", "1") == "1"
const PROFILE_ALLOCATION_DIAGNOSTICS = get(ENV,
    "ADAPTIVEOPTICS_PROFILE_ALLOCATION_DIAGNOSTICS", "0") == "1"
const PROFILE_WARMUPS = 20

PROFILE_STEPS > 0 || error("ADAPTIVEOPTICS_PROFILE_STEPS must be positive")

if PROFILE_BACKEND == "cuda"
    import CUDA
elseif PROFILE_BACKEND == "amdgpu"
    import AMDGPU
elseif PROFILE_BACKEND != "cpu"
    error("ADAPTIVEOPTICS_PROFILE_BACKEND must be cpu, cuda, or amdgpu")
end

function profile_backend(name::AbstractString)
    if name == "cpu"
        return CPUBackend()
    elseif name == "cuda"
        CUDA.functional() || error("CUDA is not functional")
        Backends.disable_scalar_backend!(Backends.CUDABackendTag)
        return Backends.CUDABackend()
    end
    AMDGPU.functional() || error("AMDGPU is not functional")
    Backends.disable_scalar_backend!(Backends.AMDGPUBackendTag)
    return Backends.AMDGPUBackend()
end

function prepare_bi_o_edge_profile(backend::AbstractArrayBackend)
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
    sensor = BiOEdgeWFS(
        telescope;
        pupil_samples=16,
        modulation=T(3),
        modulation_points=8,
        T=T,
        backend=backend,
    )
    front_end = BiOEdgeOpticalFrontEnd(sensor, source)
    rate = bi_o_edge_rate_map(front_end, pupil)
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
        detector, rate, observation; source)
    return (; pupil, rate, observation, optics, acquisition,
        rng=Xoshiro(0xB10E))
end

@inline function step_bi_o_edge_profile!(prepared)
    form_wfs_optical_products!(
        prepared.rate, prepared.pupil, prepared.optics)
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

@inline function replay_bi_o_edge_profile!(prepared, steps::Int)
    @inbounds for _ in 1:steps
        step_bi_o_edge_profile!(prepared)
    end
    return nothing
end

@inline _profile_region!(f::F, ::Backends.CPUBackend, ::AbstractString) where {
    F<:Function,
} = f()

@inline function _profile_region!(
    f::F, ::Backends.CUDABackend, label::AbstractString,
) where {F<:Function}
    ccall((:nvtxRangePushA, "libnvToolsExt"), Cint, (Cstring,), label)
    try
        status = ccall((:cuProfilerStart, "libcuda"), Cint, ())
        iszero(status) || error("cuProfilerStart failed with status $status")
        return f()
    finally
        status = ccall((:cuProfilerStop, "libcuda"), Cint, ())
        iszero(status) || error("cuProfilerStop failed with status $status")
        ccall((:nvtxRangePop, "libnvToolsExt"), Cint, ())
    end
end

@inline function _profile_region!(
    f::F, ::Backends.AMDGPUBackend, label::AbstractString,
) where {F<:Function}
    ccall((:roctxRangePushA, "libroctx64"), Cint, (Cstring,), label)
    try
        status = ccall(
            (:roctxProfilerResume, "librocprofiler-sdk-roctx"),
            Cint, (UInt64,), UInt64(0))
        iszero(status) || error(
            "roctxProfilerResume failed with status $status")
        return f()
    finally
        status = ccall(
            (:roctxProfilerPause, "librocprofiler-sdk-roctx"),
            Cint, (UInt64,), UInt64(0))
        iszero(status) || error(
            "roctxProfilerPause failed with status $status")
        ccall((:roctxRangePop, "libroctx64"), Cint, ())
    end
end

function run_bi_o_edge_profile(backend)
    prepared = prepare_bi_o_edge_profile(backend)
    replay_bi_o_edge_profile!(prepared, PROFILE_WARMUPS)
    allocation_bytes = @allocated step_bi_o_edge_profile!(prepared)
    if !iszero(allocation_bytes) && PROFILE_ALLOCATION_DIAGNOSTICS
        Profile.Allocs.clear()
        Profile.Allocs.start(sample_rate=1.0)
        try
            step_bi_o_edge_profile!(prepared)
        finally
            Profile.Allocs.stop()
        end
        Profile.Allocs.print(stdout, Profile.Allocs.fetch())
    end
    if backend isa Backends.CPUBackend && !iszero(allocation_bytes)
        error(
            "Bi-O-edge CPU plant step allocated " *
            "$allocation_bytes warmed Julia bytes")
    end

    start_nanoseconds = time_ns()
    _profile_region!(backend, "aos_bi_o_edge_plant") do
        replay_bi_o_edge_profile!(prepared, PROFILE_STEPS)
    end
    elapsed_nanoseconds = time_ns() - start_nanoseconds
    values = Array(observation_storage(prepared.observation))
    all(isfinite, values) || error("Bi-O-edge plant produced non-finite values")
    sum(values) > 0 || error("Bi-O-edge plant produced an empty observation")

    println("bi_o_edge_plant_profile")
    println("  backend=", PROFILE_BACKEND)
    println("  steps=", PROFILE_STEPS)
    println("  target_ready_elapsed_ns=", elapsed_nanoseconds)
    println("  target_ready_mean_ns=", elapsed_nanoseconds / PROFILE_STEPS)
    println("  warmed_host_allocation_bytes=", allocation_bytes)
    println("  output_sum=", sum(values))

    if PROFILE_BACKEND == "cpu" && PROFILE_PRINT_CPU_SAMPLES
        Profile.clear()
        Profile.@profile replay_bi_o_edge_profile!(prepared, PROFILE_STEPS)
        Profile.print(stdout; format=:flat, sortedby=:count)
    end
    return nothing
end

run_bi_o_edge_profile(profile_backend(PROFILE_BACKEND))

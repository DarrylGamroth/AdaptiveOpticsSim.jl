# CPU parity and architecture timelines for two Pyramid propagation strategies.
# AOS_PYRAMID_PROFILE=1 enables one profiler range per warmed case.
# AOS_PYRAMID_COMPLETION=cooperative|blocking overrides the default for diagnosis.
import CUDA
using LinearAlgebra, Statistics, TOML
include(joinpath(@__DIR__, "..", "scripts", "profile_pyramid_plant_runtime.jl"))

# Diagnostic process-local dispatch: preserve completion boundaries while
# comparing CUDA's cooperative wait with a direct blocking stream wait.
const _PYRAMID_COMPLETION = get(ENV, "AOS_PYRAMID_COMPLETION", "default")
_PYRAMID_COMPLETION in ("default", "cooperative", "blocking") ||
    error("unknown AOS_PYRAMID_COMPLETION policy")
if _PYRAMID_COMPLETION != "default"
    @eval function Backends.synchronize_backend!(
        ::Backends.AcceleratorStyle{CUDA.CUDABackend};
        blocking::Bool=$(_PYRAMID_COMPLETION == "blocking"))
        CUDA.synchronize(; blocking)
        return nothing
    end
end

function pyramid_efficiency_case(resolution, label, strategy)
    cpu = prepare_pyramid_profile(CPUBackend(), strategy;
        resolution, pupil_samples=div(resolution, 4))
    step_pyramid_profile!(cpu)
    reference = copy(observation_storage(cpu.observation))
    gpu = prepare_pyramid_profile(CUDABackend(), strategy;
        resolution, pupil_samples=div(resolution, 4))
    CUDA.synchronize()
    for _ in 1:32
        step_pyramid_profile!(gpu)
    end
    frame = Array(observation_storage(gpu.observation))
    @assert isapprox(frame, reference; rtol=5f-4, atol=1f-5)
    samples = zeros(UInt64, 128)
    gc_bytes = zeros(Int, length(samples))
    gc_collections = zeros(Int, length(samples))
    memory_before = CUDA.used_memory()
    for i in eachindex(samples)
        before = Base.gc_num()
        start = time_ns()
        step_pyramid_profile!(gpu)
        samples[i] = time_ns() - start
        difference = Base.GC_Diff(Base.gc_num(), before)
        gc_bytes[i] = difference.allocd
        gc_collections[i] = difference.pause
    end
    @assert CUDA.used_memory() == memory_before
    if get(ENV, "AOS_PYRAMID_PROFILE", "0") == "1"
        CUDA.@profile external=true begin
            for _ in 1:32
                step_pyramid_profile!(gpu)
            end
        end
    end
    frame = Array(observation_storage(gpu.observation))
    @assert isapprox(frame, reference; rtol=5f-4, atol=1f-5)
    return (; resolution, strategy=label, pupil_samples=div(resolution, 4),
        modulation_points=8, frame_shape=collect(size(frame)),
        maximum_cpu_error=maximum(abs, frame - reference),
        relative_cpu_error=norm(frame-reference)/norm(reference),
        samples_ns=samples, gc_bytes, gc_collections, p50_ns=median(samples),
        p99_ns=sort(samples)[ceil(Int, .99length(samples))])
end

function profile_pyramid_efficiency()
    CUDA.allowscalar(false)
    BLAS.set_num_threads(1)
    Backends.set_fft_provider_threads!(1)
    output = get(ENV, "AOS_PYRAMID_RESULTS", "/tmp/aos-pyramid-efficiency.toml")
    open(output, "w") do io
        TOML.print(io, Dict("source_revision"=>readchomp(`git rev-parse HEAD`),
            "device"=>CUDA.name(CUDA.device()), "julia"=>string(VERSION),
            "synchronization"=>_PYRAMID_COMPLETION,
            "boundary"=>"Host time_ns around complete synchronized stream optics+noiseless acquisition;128 samples,32 warmups; timings from profiler invocation are diagnostic only"))
        for n in (64, 256)
            for (label, strategy) in (("pupil_tilt", WavefrontSensors.PyramidPupilTiltStrategy()),
                ("shifted_mask", WavefrontSensors.PyramidShiftedMaskStrategy()))
                result = pyramid_efficiency_case(n, label, strategy)
                println(io, "\n[[case]]")
                TOML.print(io, Dict(string(k)=>v for (k,v) in pairs(result)))
                println("completed_pyramid_N$n strategy=$label p50_ns=$(result.p50_ns)")
                flush(stdout)
            end
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    profile_pyramid_efficiency()
end

# Focused CUDA capture selector; no optical propagation or detector RNG runs.
# JULIA_NUM_THREADS=1 julia --startup-file=no --project=<CUDA environment> test/cuda_external_exposure.jl
using Test, AdaptiveOpticsSim, CUDA
using AdaptiveOpticsSim.AlgorithmGraphs
import AdaptiveOpticsSim.Backends

CUDA.allowscalar(false)

function captured_exposure_step!(emit, execution, exposure, intensity, start, span)
    Backends._with_prepared_device_execution_context(execution.context) do
        AlgorithmGraphs._integrate_intensity!(emit,exposure,intensity,
            ModelTimestamp(start),ModelDuration(span),execution)
        Backends._synchronize_prepared_device_execution_context!(execution.context)
    end
end

@testset "CUDA prepared external exposure helpers" begin
    target = Backends.compute_device(CUDA.zeros(Float32,1))
    for T in (Float32,Float64)
        intensity = CUDA.zeros(T,17,19)
        gpu = IntensityExposure(size(intensity),ModelDuration(10);target)
        streamed = IntensityExposure(size(intensity),ModelDuration(10);target)
        cpu = IntensityExposure(size(intensity),ModelDuration(10))
        context = Backends._prepare_device_execution_context(intensity)
        execution = AlgorithmGraphs._prepare_external_exposure_execution(target,
            (;context),gpu,intensity,nothing)
        @test execution !== nothing
        @test eltype(gpu.integrated) === Float64
        @test length(execution.validation) == cld(length(intensity),256)
        @test all(graph -> isempty(graph.workspace),
            (execution.validate,execution.accumulate,execution.finish,execution.clear))
        host_frames = Matrix{Float64}[]
        device_frames = Matrix{Float64}[]
        host_emit = (image,args...) -> push!(host_frames,copy(image))
        device_emit = (image,args...) -> push!(device_frames,Array(image))
        cursor = 0
        # Non-dividing coverage crosses exposure boundaries and leaves a partial
        # exposure. Compare every output to the original ordered CPU broadcast.
        for (span,value) in ((2,T(1e-3)),(3,T(1e30)),(7,T(1e35)),(8,T(0.1)))
            host = fill(value,size(intensity))
            copyto!(intensity,host)
            integrate_intensity!(host_emit,cpu,host,ModelTimestamp(cursor),ModelDuration(span))
            integrate_intensity!((args...) -> nothing,streamed,intensity,
                ModelTimestamp(cursor),ModelDuration(span))
            captured_exposure_step!(device_emit,execution,gpu,intensity,cursor,span)
            cursor += span
            @test Array(gpu.integrated) == cpu.integrated
            @test Array(gpu.image) ≈ cpu.image rtol=2eps(Float64)
            @test Array(gpu.integrated) == Array(streamed.integrated)
            @test Array(gpu.image) == Array(streamed.image)
            @test (gpu.start_ns,gpu.covered_ns,gpu.next_ns,gpu.sequence,gpu.incomplete) ==
                (cpu.start_ns,cpu.covered_ns,cpu.next_ns,cpu.sequence,cpu.incomplete)
        end
        @test all(isapprox(device,host;rtol=2eps(Float64))
            for (device,host) in zip(device_frames,host_frames))
        before = (Array(gpu.integrated),Array(gpu.image),gpu.next_ns,gpu.sequence)
        for invalid in (T(-1),T(NaN),T(Inf))
            fill!(intensity,invalid)
            @test_throws AlgorithmGraphError captured_exposure_step!(device_emit,
                execution,gpu,intensity,cursor,1)
            @test (Array(gpu.integrated),Array(gpu.image),gpu.next_ns,gpu.sequence) == before
        end
        fill!(intensity,T(3))
        failure = ErrorException("consumer rejected image")
        callback = function (image,args...)
            @test Array(gpu.integrated) == fill(30.0,size(intensity))
            throw(failure)
        end
        @test_throws typeof(failure) captured_exposure_step!(callback,execution,gpu,
            intensity,cursor,10)
        @test gpu.sequence == 3 && gpu.covered_ns == 10 && gpu.start_ns == 20
        @test Array(gpu.integrated) == fill(30.0,size(intensity))
        Backends._with_prepared_device_execution_context(context) do
            AlgorithmGraphs._reset_intensity_exposure!(gpu,ModelTimestamp(33),execution)
            Backends._synchronize_prepared_device_execution_context!(context)
        end
        @test all(iszero,Array(gpu.integrated))
        @test Array(gpu.image) == fill(3.0,size(intensity))
        @test gpu.start_ns == 30 && gpu.next_ns == 33 && gpu.incomplete == 1
    end
end

# JULIA_NUM_THREADS=1 julia --startup-file=no --project=<CUDA environment> test/cuda_external_pupil_execution.jl
using Test, AdaptiveOpticsSim, CUDA
using AdaptiveOpticsSim.AlgorithmGraphs
import AdaptiveOpticsSim.Backends
CUDA.allowscalar(false)

include("external_exposure_execution_fixture.jl")

function delayed_exposure_writer(cycles)
    start = CUDA.clock(UInt64)
    while CUDA.clock(UInt64)-start < cycles
    end
    return nothing
end

test_external_exposure_pipeline(Backends.compute_device(CUDA.zeros(Float32,1)))

writer = CUDA.CuStream(;flags=CUDA.STREAM_NON_BLOCKING)
CUDA.@cuda stream=writer delayed_exposure_writer(UInt64(0))
CUDA.synchronize(writer)
function mutate_later_exposure!(later,unrelated_stream)
    if unrelated_stream
        CUDA.stream!(writer) do
            CUDA.@cuda stream=writer delayed_exposure_writer(UInt64(100_000_000))
            fill!(later,-1f0)
        end
    else
        fill!(later,-1f0)
    end
    return nothing
end
test_external_exposure_callback_mutation(Backends.compute_device(CUDA.zeros(Float32,1)),
    mutate_later_exposure!,() -> CUDA.synchronize(writer))

test_external_exposure_callback_failure(Backends.compute_device(CUDA.zeros(Float32,1)))

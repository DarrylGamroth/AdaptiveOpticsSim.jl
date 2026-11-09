# JULIA_NUM_THREADS=1 julia --startup-file=no --project=<AMDGPU environment> test/amdgpu_external_exposure_execution.jl
using Test, AdaptiveOpticsSim, AMDGPU
using AdaptiveOpticsSim.AlgorithmGraphs
import AdaptiveOpticsSim.Backends
AMDGPU.allowscalar(false)
include("external_exposure_execution_fixture.jl")
target = Backends.compute_device(AMDGPU.zeros(Float32,1))
test_external_exposure_pipeline(target)
test_external_exposure_callback_failure(target)

function delayed_amd_exposure_writer!(scratch,iterations)
    value = 0.9f0
    for index in 1:iterations
        value = sin(value)+0.001f0
    end
    scratch[1] = value
    return nothing
end
writer = AMDGPU.HIP.HIPStream()
scratch = AMDGPU.zeros(Float32,1)
AMDGPU.@roc stream=writer delayed_amd_exposure_writer!(scratch,1)
AMDGPU.synchronize(writer)
function mutate_later_amd_exposure!(later,unrelated_stream)
    if unrelated_stream
        AMDGPU.stream!(writer) do
            AMDGPU.@roc stream=writer delayed_amd_exposure_writer!(scratch,5_000_000)
            fill!(later,-1f0)
        end
    else
        fill!(later,-1f0)
    end
    return nothing
end
test_external_exposure_callback_mutation(target,mutate_later_amd_exposure!,
    () -> AMDGPU.synchronize(writer))

import CUDA
include("runtests_gpu_target_common.jl")

function run_graph_rng_capture_replay(
    ::Type{AdaptiveOpticsSim.Backends.CUDABackendTag},
)
    output = CUDA.zeros(Float32, 256)
    target = compute_device(output)
    rng = AdaptiveOpticsSim.Backends._prepare_graph_rng(
        target,
        UInt64(0x5eed),
    )
    style = AdaptiveOpticsSim.Backends.execution_style(output)

    # Compile both kernels before CUDA begins capture.
    AdaptiveOpticsSim.Backends.randn_backend_async!(style, rng, output)
    CUDA.synchronize()
    AdaptiveOpticsSim.Backends._reset_graph_rng!(rng, UInt64(0x5eed))
    CUDA.synchronize()

    graph = CUDA.capture() do
        AdaptiveOpticsSim.Backends.randn_backend_async!(style, rng, output)
    end
    executable = CUDA.instantiate(graph)

    CUDA.launch(executable)
    CUDA.synchronize()
    first_draw = Array(output)
    CUDA.launch(executable)
    CUDA.synchronize()
    second_draw = Array(output)
    @test second_draw != first_draw

    AdaptiveOpticsSim.Backends._reset_graph_rng!(rng, UInt64(0x5eed))
    CUDA.synchronize()
    CUDA.launch(executable)
    CUDA.synchronize()
    @test Array(output) == first_draw
    return nothing
end

run_gpu_backend_target(AdaptiveOpticsSim.Backends.CUDABackendTag)

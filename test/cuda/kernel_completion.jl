using Test
import CUDA
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Backends

function _cuda_kernel_completion_work!(output, rounds)
    index = CUDA.threadIdx().x
    if index <= length(output)
        value = 0.0f0
        for _ in 1:rounds
            value = muladd(value, 0.99999f0, 0.00001f0)
        end
        output[index] = value
    end
    return nothing
end

function _cuda_kernel_completion_launch!(output, rounds)
    CUDA.@cuda threads=4 blocks=1 _cuda_kernel_completion_work!(output, rounds)
    return nothing
end

@testset "CUDA synchronous kernel completion" begin
    CUDA.allowscalar(false)
    output = CUDA.zeros(Float32, 4)
    style = Backends.execution_style(output)
    original_stream, original_context = CUDA.stream(), CUDA.context()
    streams = (CUDA.CuStream(), CUDA.CuStream())
    for stream in streams
        CUDA.stream!(stream) do
            # Compile the kernel before exercising a pending completion.
            _cuda_kernel_completion_launch!(output, 1)
            Backends.synchronize_backend!(style)
            _cuda_kernel_completion_launch!(output, 2_000_000)
            @test Backends.synchronize_backend!(style) === nothing
            @test CUDA.isdone(stream)
            @test CUDA.stream() == stream
            @test CUDA.context() == original_context
            @test all(x -> isfinite(x) && x > 0.9f0, Array(output))

            # Explicit cooperative completion preserves the same result and
            # stream scope for low-level callers that require task cooperation.
            _cuda_kernel_completion_launch!(output, 2_000_000)
            @test Backends.synchronize_backend!(style; blocking=false) === nothing
            @test CUDA.isdone(stream)
            @test CUDA.stream() == stream
            @test all(x -> isfinite(x) && x > 0.9f0, Array(output))
        end
        @test CUDA.stream() == original_stream
    end
    @test CUDA.context() == original_context
end

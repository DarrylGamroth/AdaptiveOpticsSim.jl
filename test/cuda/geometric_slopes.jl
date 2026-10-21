using Test
using AdaptiveOpticsSim
using CUDA

CUDA.functional() || error("CUDA hardware is required for geometric slope tests")
CUDA.allowscalar(false)
include("../geometric_slopes_numerics.jl")
test_geometric_slopes_numerics(CuArray,
    AdaptiveOpticsSim.Backends.execution_style(CUDA.zeros(Float32, 1)))

@testset "CUDA geometric slope public calls" begin
    for T in (Float32, Float64), sub in (2, 16)
        opd = CuArray(T[x / 8 - y / 4 for x in 1:(2sub), y in 1:(2sub)])
        valid = CuArray(trues(2, 2))
        slopes = CUDA.zeros(T, 8)
        sampling = (T(4), T(0.25))
        AdaptiveOpticsSim.WavefrontSensors.geometric_wavefront_slopes!(
            slopes, opd, valid, sampling)
        @test Array(slopes) == T[fill(T(1 / 32), 4); fill(T(-1), 4)]
    end
end

function _complete_captured_geometric_slopes!(captured, context)
    AdaptiveOpticsSim.Backends._with_prepared_device_execution_context(context) do
        AdaptiveOpticsSim.Backends._launch_prepared_device_graph!(captured, context)
        AdaptiveOpticsSim.Backends._synchronize_prepared_device_execution_context_blocking!(context)
    end
    return nothing
end

@testset "CUDA geometric slopes fixed-storage replay" begin
    backends = AdaptiveOpticsSim.Backends
    wfs = AdaptiveOpticsSim.WavefrontSensors
    for T in (Float32, Float64)
        n, n_sub = 64, 4
        opd = T[sin(T(i) / 7) + cos(T(j) / 11) for i in 1:n, j in 1:n]
        mask = trues(n_sub, n_sub)
        reference = zeros(T, 2n_sub^2)
        wfs.geometric_slopes!(reference, opd, mask)
        output = CUDA.zeros(T, length(reference))
        device_opd, device_mask = CuArray(opd), CuArray(mask)
        style = backends.execution_style(output)
        context = backends._prepare_device_execution_context(output)
        args = (output, device_opd, device_mask, div(n, n_sub), n_sub,
            n_sub^2, one(T), one(T))
        CUDA.synchronize()
        captured = backends._with_prepared_device_execution_context(context) do
            wfs._enqueue_geometric_slopes!(style, args...)
            backends._synchronize_prepared_device_execution_context_blocking!(context)
            backends._capture_prepared_device_graph(context) do
                wfs._enqueue_geometric_slopes!(style, args...)
            end
        end
        for _ in 1:16
            _complete_captured_geometric_slopes!(captured, context)
        end
        # Function barrier keeps the exact executable/context concrete.
        bytes = @allocated _complete_captured_geometric_slopes!(captured, context)
        @test bytes == 0
        @test Array(output) ≈ reference rtol=64eps(T) atol=64eps(T)
        @test size(output) == size(reference)
    end
end

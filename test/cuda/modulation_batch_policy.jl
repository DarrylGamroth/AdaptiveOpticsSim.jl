# Run separately from modulation_batches.jl, which installs a benchmark override.
using AdaptiveOpticsSim, Test
using AdaptiveOpticsSim.Backends
import CUDA
const PolicyWFS = AdaptiveOpticsSim.WavefrontSensors
CUDA.allowscalar(false)

@testset "CUDA modulation batch preparation policy" begin
    qualified = CUDA.capability(CUDA.device()) == v"8.6" &&
        CUDA.attribute(CUDA.device(), CUDA.DEVICE_ATTRIBUTE_L2_CACHE_SIZE) == 5 * 1024^2
    for T in (ComplexF32, ComplexF64), side in (256, 272, 408, 512, 528, 560, 576)
        field = CUDA.zeros(T, side, side)
        bytes = sizeof(T) * length(field)
        for points in (1, 8, 12, 24, 30, 31, 32, 64, 96)
            baseline = PolicyWFS._pyramid_modulation_batch_size(points, bytes)
            expected = baseline
            supported = (side == 256 && points == 32) ||
                (side == 272 && points in (24, 30, 32)) ||
                (side == 408 && points == 32) ||
                (side == 512 && points in (12, 32, 96)) ||
                (side == 528 && points == 64)
            if qualified && T === ComplexF32 && supported
                cap = side in (256, 272) ? 8 : 2
                candidate = PolicyWFS._pyramid_modulation_batch_size(points, bytes, cap)
                candidate >= 2 && (expected = candidate)
            end
            selected = PolicyWFS._pyramid_modulation_batch_size(field, points)
            @test selected == expected
            @test points % selected == 0
            @test bytes * selected <= 256 * 1024^2
        end
    end
end

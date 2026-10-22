# Run separately from modulation_batches.jl, which installs a benchmark override.
using AdaptiveOpticsSim, Test
using AdaptiveOpticsSim.Backends
import CUDA
const PolicyWFS = AdaptiveOpticsSim.WavefrontSensors
CUDA.allowscalar(false)

@testset "Modulation batch scratch budget" begin
    @test PolicyWFS._pyramid_modulation_batch_size(8, 1024, 8, 2048) == 2
    @test PolicyWFS._pyramid_modulation_batch_size(31, 1024, 8, 2048) == 1
    @test PolicyWFS._pyramid_modulation_batch_size(192, 2880^2 * 8, 16) == 4
    @test PolicyWFS._pyramid_modulation_batch_size(256, 3072^2 * 8, 16) == 2
    for (side, points) in ((2880, 192), (3072, 256))
        @test PolicyWFS._pyramid_modulation_batch_size(
            points, side^2 * 8, 16, 1536 * 1024^2) == 16
    end
    for budget in (0, -1)
        @test_throws InvalidConfiguration PolicyWFS._pyramid_modulation_batch_size(
            8, 1024, 8, budget)
    end
end

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

    @testset "Fine grids and neighboring configurations" begin
        for T in (ComplexF32, ComplexF64), side in (2878, 2880, 2882, 3070, 3072, 3074)
            field = CUDA.zeros(T, side, side)
            bytes = sizeof(T) * length(field)
            style = Backends.execution_style(field)
            for points in (191, 192, 193, 255, 256, 257)
                supported = (side == 2880 && points == 192) ||
                    (side == 3072 && points == 256)
                expected_budget = qualified && T === ComplexF32 && supported ?
                    1536 * 1024^2 : 256 * 1024^2
                expected = qualified && T === ComplexF32 && supported ? 16 :
                    PolicyWFS._pyramid_modulation_batch_size(points, bytes)
                budget = PolicyWFS._pyramid_modulation_batch_budget_bytes(
                    style, field, points, bytes)
                selected = PolicyWFS._pyramid_modulation_batch_size(field, points)
                @test budget == expected_budget
                @test selected == expected
                @test points % selected == 0
                @test bytes * selected <= budget
            end
        end
        field = CUDA.zeros(ComplexF32, 3072, 3071)
        bytes = sizeof(eltype(field)) * length(field)
        @test PolicyWFS._pyramid_modulation_batch_budget_bytes(
            Backends.execution_style(field), field, 256, bytes) == 256 * 1024^2
        @test PolicyWFS._pyramid_modulation_batch_size(field, 256) ==
            PolicyWFS._pyramid_modulation_batch_size(256, bytes)
    end
end

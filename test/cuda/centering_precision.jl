using Test
using AdaptiveOpticsSim
import CUDA

@testset "CUDA pyramid centering-factor precision" begin
    CUDA.allowscalar(false)
    WFS = AdaptiveOpticsSim.WavefrontSensors
    for T in (Float32, Float64), n in (31, 512, 1024)
        host = Matrix{Complex{T}}(undef, n, n)
        device = CUDA.CuArray{Complex{T}}(undef, n, n)
        WFS.build_pyramid_phasor!(host)
        @test WFS.build_pyramid_phasor!(device) === device
        @test Array(device) == host
        for i in (1, cld(n, 2), n), j in (1, cld(n, 2), n)
            expected = Complex{T}(cispi(-BigFloat(n + 1) * (i + j - 2) / n))
            tolerance = T === Float32 ? 2eps(T) : 8n * eps(T)
            @test abs(host[i, j] - expected) <= tolerance
        end
    end
end

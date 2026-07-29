using Test
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Backends
using Metal

@testset "Metal extension owner surface" begin
    @test Sys.isapple()
    @test Sys.ARCH === :aarch64
    @test Base.get_extension(
        AdaptiveOpticsSim, :AdaptiveOpticsSimMetalExt) !== nothing

    @test Backends.gpu_backend_loaded(Backends.MetalBackendTag)
    @test Backends.gpu_backend_array_type(
        Backends.MetalBackendTag) === Metal.MtlArray
    @test Backends.gpu_backend_name(Metal.MtlArray) === :metal
    @test Backends.array_backend_selector(
        Metal.MtlArray) isa Backends.MetalBackend
    @test Backends.array_backend_type(
        Backends.MetalBackend()) === Metal.MtlArray
    @test Backends.MetalBackendTag in Backends.available_gpu_backends()

    if Metal.functional()
        uniform = Backends.backend_rand(
            Backends.MetalBackendTag, Float32, 2, 2)
        normal = Backends.backend_randn(
            Backends.MetalBackendTag, Float32, 2, 2)
        zeros_array = Backends.backend_zeros(
            Backends.MetalBackendTag, Float32, 2, 2)
        filled = Backends.backend_fill(
            Backends.MetalBackendTag, 3.0f0, 2, 2)
        @test uniform isa Metal.MtlArray
        @test normal isa Metal.MtlArray
        @test zeros_array isa Metal.MtlArray
        @test filled isa Metal.MtlArray
        @test all(isfinite, Array(uniform))
        @test all(isfinite, Array(normal))
        @test iszero(sum(Array(zeros_array)))
        @test Array(filled) == fill(3.0f0, 2, 2)
        @test Backends.compute_device_identifier(zeros_array) isa UInt
    else
        @test_skip "Metal device execution is unavailable on this host"
    end
end

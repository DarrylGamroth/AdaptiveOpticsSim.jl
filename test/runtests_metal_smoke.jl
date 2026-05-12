include("gpu_target_env.jl")

using Test
using AdaptiveOpticsSim
using LinearAlgebra
using Random

if !(Sys.isapple() && Sys.ARCH === :aarch64)
    @info "Skipping Metal hosted runner smoke: Metal target requires macOS arm64"
    @testset "Metal hosted runner smoke" begin
        @test true
    end
    exit()
end

using Metal

function run_metal_array_smoke()
    T = Float32
    a = Metal.MtlArray(T[1, 2, 3, 4])
    b = Metal.zeros(T, 4)
    b .= muladd.(T(2), a, T(1))
    @test Array(b) == T[3, 5, 7, 9]

    c = AdaptiveOpticsSim.backend_fill(AdaptiveOpticsSim.MetalBackendTag, T(3), 2, 2)
    @test c isa Metal.MtlArray
    @test Array(c) == fill(T(3), 2, 2)
    return nothing
end

include("backend_optional_common.jl")

@testset "Metal hosted runner smoke" begin
    Metal.versioninfo()
    @test Metal.functional()
    AdaptiveOpticsSim.disable_scalar_backend!(AdaptiveOpticsSim.MetalBackendTag)

    @test AdaptiveOpticsSim.gpu_backend_loaded(AdaptiveOpticsSim.MetalBackendTag)
    @test AdaptiveOpticsSim.gpu_backend_array_type(AdaptiveOpticsSim.MetalBackendTag) === Metal.MtlArray
    @test AdaptiveOpticsSim.gpu_backend_name(AdaptiveOpticsSim.MetalBackendTag) === :metal
    @test AdaptiveOpticsSim.array_backend_type(AdaptiveOpticsSim.MetalBackend()) === Metal.MtlArray
    @test AdaptiveOpticsSim.available_gpu_backends() == (AdaptiveOpticsSim.MetalBackendTag,)

    run_metal_array_smoke()
    run_optional_backend_selector_smoke(AdaptiveOpticsSim.MetalBackendTag, Metal.MtlArray)
end

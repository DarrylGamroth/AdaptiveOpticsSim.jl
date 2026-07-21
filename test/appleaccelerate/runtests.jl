using Test
using LinearAlgebra
using AppleAccelerate
using AdaptiveOpticsSim
import FFTW

function accelerate_backing_library(symbol::AbstractString)
    info = BLAS.lbt_find_backing_library(symbol, :ilp64)
    info === nothing && return nothing
    return info.libname
end

@testset "Explicit AppleAccelerate selection" begin
    @test Sys.isapple()
    @test Sys.ARCH === :aarch64
    @test Base.pkgversion(AppleAccelerate) >= v"0.7.0"
    @test AppleAccelerate.get_macos_version() >= v"15"
    @test AppleAccelerate.set_num_threads(1) == 1
    @test AppleAccelerate.get_num_threads() == 1

    blas_library = accelerate_backing_library("dgemm_")
    lapack_library = accelerate_backing_library("dgetrf_")
    @test blas_library !== nothing
    @test lapack_library !== nothing
    @test occursin("Accelerate", blas_library)
    @test occursin("Accelerate", lapack_library)

    A = [4.0 1.0; 2.0 3.0]
    B = [2.0 -1.0; 0.5 5.0]
    product = similar(A)
    mul!(product, A, B)
    @test product == [8.5 1.0; 5.5 13.0]
    rhs = [1.0, -2.0]
    @test A * (lu(A) \ rhs) ≈ rhs

    original = ComplexF64.(reshape(1:64, 8, 8))
    transformed = copy(original)
    fft_plan = AdaptiveOpticsSim.plan_fft_backend!(transformed)
    ifft_plan = AdaptiveOpticsSim.plan_ifft_backend!(transformed)
    apple_fft_extension =
        Base.get_extension(AppleAccelerate, :AppleAccelerateAbstractFFTsExt)
    @test apple_fft_extension !== nothing
    @test parentmodule(typeof(fft_plan)) === apple_fft_extension
    @test occursin("VDSPInplaceFFTPlan", string(typeof(fft_plan)))
    @test occursin("VDSPInplaceBFFTPlan", string(typeof(ifft_plan)))
    AdaptiveOpticsSim.execute_fft_plan!(transformed, fft_plan)
    AdaptiveOpticsSim.execute_fft_plan!(transformed, ifft_plan)
    @test transformed ≈ original

    # vDSP is deliberately selected only inside its documented capability
    # envelope. FFTW preserves arbitrary-size and 3D+ package behavior.
    arbitrary_size = zeros(ComplexF64, 6, 8)
    higher_dimensional = zeros(ComplexF64, 4, 4, 2)
    @test parentmodule(typeof(
        AdaptiveOpticsSim.plan_fft_backend!(arbitrary_size))) === FFTW
    @test parentmodule(typeof(
        AdaptiveOpticsSim.plan_fft_backend!(higher_dimensional))) === FFTW

    println("AppleAccelerate version: ", Base.pkgversion(AppleAccelerate))
    println("BLAS/LAPACK config: ", BLAS.get_config())
    println("AdaptiveOpticsSim FFT plan: ", typeof(fft_plan))
end

include(joinpath(dirname(@__DIR__), "runtests.jl"))

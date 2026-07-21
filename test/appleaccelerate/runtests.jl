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

function apple_fft_cycle!(buffer, fft_plan, ifft_plan)
    AdaptiveOpticsSim.execute_fft_plan!(buffer, fft_plan)
    AdaptiveOpticsSim.execute_fft_plan!(buffer, ifft_plan)
    return buffer
end

function deterministic_fft_input(::Type{T}, shape) where {T<:AbstractFloat}
    buffer = Array{Complex{T}}(undef, shape)
    @inbounds for (index, linear_index) in enumerate(eachindex(buffer))
        buffer[linear_index] = Complex{T}(
            T(index) / T(7), T(mod(index, 5) - 2) / T(11))
    end
    return buffer
end

fft_rtol(::Type{Float32}) = 2f-5
fft_rtol(::Type{Float64}) = 2e-12

function validate_apple_fft_case(
    ::Type{T}, shape::NTuple{N,Int}, apple_fft_extension,
) where {T<:AbstractFloat,N}
    original = deterministic_fft_input(T, shape)
    reference = copy(original)
    dims = ntuple(identity, N)
    reference_plan = AdaptiveOpticsSim._plan_fftw_fft!(reference, dims)
    AdaptiveOpticsSim.execute_fft_plan!(reference, reference_plan)

    transformed = copy(original)
    fft_plan = AdaptiveOpticsSim.plan_fft_backend!(transformed)
    ifft_plan = AdaptiveOpticsSim.plan_ifft_backend!(transformed)
    @test parentmodule(typeof(fft_plan)) === apple_fft_extension
    @test occursin("AdaptiveOpticsAppleFFTPlan", string(typeof(fft_plan)))
    @test occursin("AppleFFTForward", string(typeof(fft_plan)))
    @test occursin("AppleFFTInverse", string(typeof(ifft_plan)))
    @test fft_plan.setup !== nothing
    @test ifft_plan.setup !== nothing

    AdaptiveOpticsSim.execute_fft_plan!(transformed, fft_plan)
    @test isapprox(transformed, reference;
        rtol=fft_rtol(T), atol=fft_rtol(T))
    AdaptiveOpticsSim.execute_fft_plan!(transformed, ifft_plan)
    @test isapprox(transformed, original;
        rtol=fft_rtol(T), atol=fft_rtol(T))

    transformed .= original
    apple_fft_cycle!(transformed, fft_plan, ifft_plan)
    transformed .= original
    @test @allocated(apple_fft_cycle!(
        transformed, fft_plan, ifft_plan)) == 0
    @test isapprox(transformed, original;
        rtol=fft_rtol(T), atol=fft_rtol(T))
    return nothing
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

    apple_fft_extension = Base.get_extension(AdaptiveOpticsSim,
        :AdaptiveOpticsSimAppleAccelerateExt)
    @test apple_fft_extension !== nothing

    for T in (Float32, Float64), shape in ((16,), (8, 16))
        validate_apple_fft_case(T, shape, apple_fft_extension)
    end

    # vDSP is deliberately selected only inside its documented capability
    # envelope. The stable wrapper preserves FFTW behavior while prepared state
    # is resized across that boundary.
    supported_size = zeros(ComplexF64, 8, 8)
    arbitrary_original = deterministic_fft_input(Float64, (6, 8))
    arbitrary_size = copy(arbitrary_original)
    higher_dimensional = zeros(ComplexF64, 4, 4, 2)
    supported_plan = AdaptiveOpticsSim.plan_fft_backend!(supported_size)
    arbitrary_plan = AdaptiveOpticsSim.plan_fft_backend!(arbitrary_size)
    arbitrary_inverse_plan = AdaptiveOpticsSim.plan_ifft_backend!(arbitrary_size)
    @test typeof(supported_plan) === typeof(arbitrary_plan)
    @test supported_plan.setup !== nothing
    @test arbitrary_plan.setup === nothing
    @test arbitrary_inverse_plan.setup === nothing
    @test parentmodule(typeof(arbitrary_plan.fftw_plan)) === FFTW
    @test parentmodule(typeof(
        AdaptiveOpticsSim.plan_fft_backend!(higher_dimensional))) === FFTW

    arbitrary_reference = copy(arbitrary_original)
    arbitrary_reference_plan = AdaptiveOpticsSim._plan_fftw_fft!(
        arbitrary_reference, (1, 2))
    AdaptiveOpticsSim.execute_fft_plan!(
        arbitrary_reference, arbitrary_reference_plan)
    AdaptiveOpticsSim.execute_fft_plan!(arbitrary_size, arbitrary_plan)
    @test isapprox(arbitrary_size, arbitrary_reference;
        rtol=fft_rtol(Float64), atol=fft_rtol(Float64))
    AdaptiveOpticsSim.execute_fft_plan!(
        arbitrary_size, arbitrary_inverse_plan)
    @test isapprox(arbitrary_size, arbitrary_original;
        rtol=fft_rtol(Float64), atol=fft_rtol(Float64))
    arbitrary_size .= arbitrary_original
    apple_fft_cycle!(arbitrary_size, arbitrary_plan, arbitrary_inverse_plan)
    arbitrary_size .= arbitrary_original
    @test @allocated(apple_fft_cycle!(
        arbitrary_size, arbitrary_plan, arbitrary_inverse_plan)) == 0

    println("AppleAccelerate version: ", Base.pkgversion(AppleAccelerate))
    println("BLAS/LAPACK config: ", BLAS.get_config())
    println("AdaptiveOpticsSim FFT plan: ", typeof(supported_plan))
end

include(joinpath(dirname(@__DIR__), "runtests.jl"))

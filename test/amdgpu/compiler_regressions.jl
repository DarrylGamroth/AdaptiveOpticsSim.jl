using Test
using AMDGPU
using LinearAlgebra

AMDGPU.functional() || error("AMDGPU compiler regressions require a functional GPU")
AMDGPU.allowscalar(false)

# These operations reach the exception-reporting runtime through checked integer
# division, including the Cartesian indexing in generic matmul and cropped copies.
# Registry AMDGPU 2.8's precompiled runtime fails to compile them on Julia 1.13 (upstream #1075).
@testset "AMDGPU compiler regressions" begin
    @testset "checked integer division" begin
        function divide_kernel!(out, numerators, denominators)
            i = AMDGPU.workitemIdx().x
            if i <= length(out)
                @inbounds out[i] = div(numerators[i], denominators[i])
            end
            return
        end
        numerators = Int64[5, -7, 9, -9]
        denominators = Int64[2, 3, -4, -4]
        out = AMDGPU.zeros(Int64, length(numerators))
        @roc groupsize=length(out) divide_kernel!(out,
            ROCArray(numerators), ROCArray(denominators))
        @test Array(out) == div.(numerators, denominators)

        # Compiling the exception path must still preserve error reporting.
        @roc groupsize=1 divide_kernel!(out, ROCArray(Int64[1]), ROCArray(Int64[0]))
        @test_throws r"GPU Kernel Exception" AMDGPU.synchronize()
        @roc groupsize=1 divide_kernel!(out, ROCArray(Int64[6]), ROCArray(Int64[2]))
        @test Array(out)[1] == 3
    end

    for T in (Float32, Float64)
        @testset "$T transpose matrix-vector multiplication" begin
            a = reshape(T.(1:15), 5, 3)
            b = T.(1:5)
            out = AMDGPU.zeros(T, 3)
            mul!(out, transpose(ROCArray(a)), ROCArray(b))
            @test Array(out) ≈ transpose(a) * b
        end

        @testset "$T cropped matrix copy" begin
            a = reshape(Complex{T}.(1:64), 8, 8)
            source = ROCArray(a)
            out = AMDGPU.zeros(Complex{T}, 4, 4)
            copyto!(out, @view source[2:5, 3:6])
            @test Array(out) == a[2:5, 3:6]
        end
    end
end

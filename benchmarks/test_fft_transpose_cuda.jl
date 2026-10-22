# Standalone verification of the benchmark-only tiled transpose.
include("fft_transpose_cuda.jl")
using Test
import CUDA
CUDA.allowscalar(false)

@testset "Tiled FFT layout transpose" begin
    FFTTransposeExperiment.TRANSPOSE == "tiled" || error("set AOS_FFT_TRANSPOSE=tiled")
    for T in (ComplexF32, ComplexF64), shape in ((1,1,1),(17,35,3),(32,32,2),(65,33,2))
        host=reshape(T.(1:prod(shape)),shape)
        source=CUDA.CuArray(host)
        destination=similar(source,T,(shape[2],shape[1],shape[3]))
        FFTTransposeExperiment.transpose_planes!(destination,source)
        CUDA.synchronize()
        @test Array(destination) == permutedims(host,(2,1,3))
        captured=CUDA.capture() do
            FFTTransposeExperiment.transpose_planes!(destination,source)
        end
        executable=CUDA.instantiate(captured)
        fill!(source,zero(T));CUDA.synchronize()
        GC.gc(true)
        CUDA.launch(executable);CUDA.synchronize()
        @test all(iszero,Array(destination))
        copyto!(source,host);CUDA.synchronize()
        CUDA.launch(executable);CUDA.synchronize()
        @test Array(destination) == permutedims(host,(2,1,3))
    end
end

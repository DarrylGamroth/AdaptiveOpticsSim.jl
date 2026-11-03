# Optional CUDA selector. CPU execution is a numerical oracle only.
include("proper/compiled_paths_helpers.jl")
using CUDA
import AdaptiveOpticsSim.Backends
CUDA.allowscalar(false)
BLAS.set_num_threads(1)
Backends.set_fft_provider_threads!(1)
const COMPILED_TARGET=Backends.AcceleratorComputeDevice(Backends.CUDABackend(),0)
const COMPILED_ARRAY=CuArray
const COMPILED_SYNCHRONIZE=CUDA.synchronize
include("proper/compiled_paths_gpu.jl")

@testset "Compiled edge storage eligibility" begin
    field=CUDA.zeros(ComplexF32,512,512)
    target=COMPILED_TARGET
    ext=Base.get_extension(AdaptiveOpticsSim,:AdaptiveOpticsSimCUDAExt)
    coefficients=Backends._prepare_mul_coefficients(field)
    @test AdaptiveOpticsSim.Optics._prepare_compiled_paths_mul_coefficients(
        CUDA.zeros(ComplexF32,8,8),coefficients) === coefficients
    if ext._compiled_paths_device_qualified(field)
        function paths(edges)
            operations=map(edge->(;kind=CP.PathTransferFunction,auxiliary=edge),edges)
            return ((prefix=operations,closed=(),only=()),)
        end
        indices=Int32[min(i,514-i) for i in 1:512]
        quarter=reshape(Float32.(1:257^2)./Float32(257^2),257,257)
        values=complex.(quarter[indices,indices])
        edge=CuArray(values)
        cache=CP._compiled_paths_prepare_edge_multiplier(paths((edge,edge)),target,field)
        @test first(cache) === edge
        @test isequal(complex.(Array(last(cache))),values)
        @test Base.parent(last(cache)) isa CUDA.CuMatrix{Float32}
        @test CP._compiled_paths_prepare_edge_multiplier(paths((edge,copy(edge))),target,field) === nothing
        @test CP._compiled_paths_prepare_edge_multiplier(paths(()),target,field) === nothing
        asymmetric=copy(values);asymmetric[2,2]+=1f0
        @test CP._compiled_paths_prepare_edge_multiplier(paths((CuArray(asymmetric),)),target,field) === nothing
        nonreal=copy(values);nonreal[2,2]+=im
        @test CP._compiled_paths_prepare_edge_multiplier(paths((CuArray(nonreal),)),target,field) === nothing
        signed=copy(values);signed[2,2]=complex(real(signed[2,2]),-0f0)
        @test CP._compiled_paths_prepare_edge_multiplier(paths((CuArray(signed),)),target,field) === nothing
        @test hasproperty(AdaptiveOpticsSim.Optics._prepare_compiled_paths_mul_coefficients(field,coefficients),:algorithm)
        CUDA.math_mode!(CUDA.PEDANTIC_MATH)
        try
            @test AdaptiveOpticsSim.Optics._prepare_compiled_paths_mul_coefficients(field,coefficients) === coefficients
        finally
            CUDA.math_mode!(CUDA.DEFAULT_MATH)
        end
    end
end

include("proper/compiled_paths_cuda_multiplication.jl")

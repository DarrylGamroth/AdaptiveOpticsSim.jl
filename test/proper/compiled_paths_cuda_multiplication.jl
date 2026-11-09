# Isolate complex map arithmetic from the fitted-model propagation gate.
using Test, Random, LinearAlgebra, AdaptiveOpticsSim, CUDA, Proper

@testset "Captured compiled complex matrix products" begin
    CUDA.allowscalar(false)
    backend=AdaptiveOpticsSim.Backends
    field=CUDA.zeros(ComplexF32,512,512)
    original=backend._prepare_mul_coefficients(field)
    prepared=AdaptiveOpticsSim.Optics._prepare_compiled_paths_mul_coefficients(field,original)
    for (n,k) in ((512,512),(221,221),(512,72)), transposed in (false,true)
        rng=MersenneTwister(61)
        a=CuArray(randn(rng,ComplexF32,n,k).*1f-4)
        b=CuArray(randn(rng,ComplexF32,transposed ? (n,k) : (k,n)).*1f4)
        right=transposed ? transpose(b) : b
        out=CUDA.zeros(ComplexF32,n,n)
        expected=similar(out)
        backend._enqueue_mul!(expected,a,right,original)
        operation()=backend._enqueue_mul!(out,a,right,prepared)
        captured=Proper.prepare_cuda_graph(operation;warmup=1)
        Proper.launch_cuda_graph!(captured);CUDA.synchronize()
        actual,reference=Array(out),Array(expected)
        @test norm(actual-reference)/norm(reference)<=1e-5
        @test maximum(abs,actual-reference)/maximum(abs,reference)<=1e-5
        fill!(a,0)
        Proper.launch_cuda_graph!(captured);CUDA.synchronize()
        @test all(iszero,Array(out))
    end
end

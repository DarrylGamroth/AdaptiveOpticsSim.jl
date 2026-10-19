# Run explicitly in a project containing CUDA and this checkout:
# julia --project=<environment> test/cuda_prepared_multiplication_capture.jl
using Test
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Backends
using CUDA
using LinearAlgebra

CUDA.allowscalar(false)
BLAS.set_num_threads(1)

function captured_mul_node_types(graph)
    driver = CUDA.CUDACore
    count = Ref{Csize_t}(0)
    driver.cuGraphGetNodes(graph, C_NULL, count)
    nodes = Vector{driver.CUgraphNode}(undef, count[])
    driver.cuGraphGetNodes(graph, nodes, count)
    return map(nodes) do node
        kind = Ref{driver.CUgraphNodeType}()
        driver.cuGraphNodeGetType(node, kind)
        return kind[]
    end
end

function captured_mul_fixture(::Type{T}, shape) where {T}
    return reshape(T[convert(T, sin(0.41i)) +
        (T <: Complex ? convert(T, 0.3im * cos(0.27i)) : zero(T))
        for i in 1:prod(shape)], shape)
end

captured_mul_allocations(C, A, B, coeffs) =
    @allocated Backends._enqueue_mul!(C, A, B, coeffs)

captured_mul_replay_allocations(executable) = @allocated CUDA.launch(executable)

@testset "CUDA prepared multiplication capture and replay" begin
    stream = CUDA.CuStream()
    CUDA.stream!(stream) do
        for T in (Float32, Float64, ComplexF32, ComplexF64), kind in (:matrix, :transpose_right, :vector)
            host_a = captured_mul_fixture(T, (13, 11))
            host_b = captured_mul_fixture(T, kind === :matrix ? (11, 9) :
                kind === :transpose_right ? (9, 11) : (11,))
            A, raw_b = CuArray(host_a), CuArray(host_b)
            B = kind === :transpose_right ? transpose(raw_b) : raw_b
            expected_b = kind === :transpose_right ? transpose(host_b) : host_b
            C = CuArray(host_a * expected_b)
            coefficients = Backends._prepare_mul_coefficients(C)
            @test coefficients.alpha.x === coefficients.beta.x
            @test Array(coefficients.alpha.x) == T[one(T), zero(T)]
            @test Backends.compute_device(coefficients.alpha.x) == Backends.compute_device(C)
            for _ in 1:20
                Backends._enqueue_mul!(C, A, B, coefficients)
            end
            CUDA.synchronize(stream)
            captured_mul_allocations(C, A, B, coefficients)
            bytes = captured_mul_allocations(C, A, B, coefficients)
            # Installed typed cuBLAS wrappers query the library version on each
            # submission (48 host bytes); this occurs outside captured replay.
            graph = CUDA.capture() do
                Backends._enqueue_mul!(C, A, B, coefficients)
                return nothing
            end
            types = captured_mul_node_types(graph)
            @test !isempty(types)
            permitted = (CUDA.CUDACore.CU_GRAPH_NODE_TYPE_KERNEL,
                CUDA.CUDACore.CU_GRAPH_NODE_TYPE_MEM_ALLOC,
                CUDA.CUDACore.CU_GRAPH_NODE_TYPE_MEMSET,
                CUDA.CUDACore.CU_GRAPH_NODE_TYPE_MEM_FREE)
            @test all(in(permitted), types)
            @test count(==(CUDA.CUDACore.CU_GRAPH_NODE_TYPE_MEM_ALLOC), types) ==
                count(==(CUDA.CUDACore.CU_GRAPH_NODE_TYPE_MEM_FREE), types)
            executable = CUDA.instantiate(graph)
            CUDA.launch(executable)
            captured_mul_replay_allocations(executable)
            @test captured_mul_replay_allocations(executable) == 0
            CUDA.synchronize(stream)
            # Retained coefficient ownership must survive collection and reuse.
            GC.gc(true)
            for replay in 1:3
                scale = real(T)(replay)
                copyto!(A, scale .* host_a)
                CUDA.launch(executable)
                CUDA.synchronize(stream)
                expected = (scale .* host_a) * expected_b
                actual = Array(C)
                @test all(isfinite, actual)
                @test actual ≈ expected rtol=64eps(real(T)) atol=64eps(real(T))
            end
            @test Array(coefficients.alpha.x) == T[one(T), zero(T)]
            println((; T, kind, nodes=string.(types), allocation_bytes=bytes, replays=3))
            flush(stdout)
        end
    end
end

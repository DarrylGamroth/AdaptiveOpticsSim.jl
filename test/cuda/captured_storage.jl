using Test
using CUDA
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Backends

CUDA.functional() || error("CUDA hardware is required for captured storage tests")
CUDA.allowscalar(false)

function _cuda_storage_child_graph(child)
    driver = CUDA.CUDACore
    parent = CUDA.CuGraph()
    node = Ref{driver.CUgraphNode}()
    driver.cuGraphAddChildGraphNode(node, parent, C_NULL, 0, child)
    return parent
end

function _cuda_storage_handle_state(handle)
    library = CUDA.CUBLAS
    stream = Ref{CUDA.CUDACore.CUstream}()
    pointer_mode = Ref{library.cublasPointerMode_t}()
    library.cublasGetStream_v2(handle, stream)
    library.cublasGetPointerMode_v2(handle, pointer_mode)
    return (; stream=stream[], pointer_mode=pointer_mode[])
end

@testset "CUDA captured graph fixed storage admission" begin
    extension = Base.get_extension(AdaptiveOpticsSim, :AdaptiveOpticsSimCUDAExt)
    driver = CUDA.CUDACore
    empty_graph = CUDA.CuGraph()
    @test extension._require_cuda_graph_fixed_storage(empty_graph) === empty_graph
    fixed_graph = CUDA.CuGraph()
    empty_node = Ref{driver.CUgraphNode}()
    driver.cuGraphAddEmptyNode(empty_node, fixed_graph, C_NULL, 0)
    nested_fixed = _cuda_storage_child_graph(_cuda_storage_child_graph(fixed_graph))
    @test extension._require_cuda_graph_fixed_storage(nested_fixed) === nested_fixed

    stream = CUDA.CuStream()
    CUDA.stream!(stream) do
        # Capturing reserves a virtual address without executing an allocation.
        # Rejected graphs are never instantiated or launched.
        allocation = Ref{driver.CUdeviceptr}()
        allocation_graph = CUDA.capture() do
            driver.cuMemAllocAsync(allocation, 256, stream)
        end
        @test_throws InvalidConfiguration extension._require_cuda_graph_fixed_storage(
            allocation_graph)

        # A free-only graph can refer to another graph's reserved address.
        # Neither graph executes, so there is no allocation or free to clean up.
        free_graph = CUDA.CuGraph()
        free_node = Ref{driver.CUgraphNode}()
        driver.cuGraphAddMemFreeNode(free_node, free_graph, C_NULL, 0, allocation[])
        GC.@preserve allocation_graph begin
            @test_throws InvalidConfiguration extension._require_cuda_graph_fixed_storage(
                free_graph)
        end
        # CUDA currently prohibits allocation/free nodes in child graphs. The
        # two-level fixed child above exercises recursive admission on real graphs.

    end
end

@testset "CUDA captured graph workspace ownership and handle restoration" begin
    extension = Base.get_extension(AdaptiveOpticsSim, :AdaptiveOpticsSimCUDAExt)
    output = CUDA.zeros(Float32, 4)
    contexts = (Backends._prepare_device_execution_context(output),
        Backends._prepare_device_execution_context(output))
    original_stream = CUDA.stream()
    captured = map(contexts) do context
        Backends._with_prepared_device_execution_context(context) do
            handle = CUDA.CUBLAS.handle()
            before = _cuda_storage_handle_state(handle)
            result = Backends._capture_prepared_device_graph(context) do
                nothing
            end
            @test _cuda_storage_handle_state(handle) == before
            @test extension._require_cuda_graph_fixed_storage(result.graph) === result.graph
            @test result.workspace isa CuArray{UInt8,1,CUDA.DeviceMemory}
            @test Backends.compute_device(result.workspace) == Backends.compute_device(output)
            @test UInt(pointer(result.workspace)) % 256 == 0
            @test sizeof(result.workspace) ==
                (CUDA.capability(context.device) >= v"9.0" ? 32 << 20 : 4 << 20)
            @test_throws ErrorException extension._with_prepared_cublas_workspace(
                context, result.workspace) do
                error("workspace scope regression")
            end
            @test _cuda_storage_handle_state(handle) == before
            @test_throws ErrorException Backends._capture_prepared_device_graph(context) do
                error("capture callback regression")
            end
            @test _cuda_storage_handle_state(handle) == before
            @test_throws InvalidConfiguration Backends._capture_prepared_device_graph(context) do
                allocation = Ref{CUDA.CUDACore.CUdeviceptr}()
                CUDA.CUDACore.cuMemAllocAsync(allocation, 256, context.stream)
            end
            @test _cuda_storage_handle_state(handle) == before
            return result
        end
    end
    @test CUDA.stream() == original_stream
    @test contexts[1].stream != contexts[2].stream
    @test captured[1].workspace !== captured[2].workspace
    @test pointer(captured[1].workspace) != pointer(captured[2].workspace)
    addresses = map(result -> pointer(result.workspace), captured)
    GC.gc(true)
    @test map(result -> pointer(result.workspace), captured) == addresses
end

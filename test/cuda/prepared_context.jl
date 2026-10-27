_prepared_context_noop() = nothing

struct _PreparedContextCallable end
(::_PreparedContextCallable)() = :called

function _prepared_context_allocation(context)
    return @allocated Backends._with_prepared_device_execution_context(
        _prepared_context_noop, context)
end

@testset "CUDA prepared context ownership" begin
    storage = CUDA.zeros(Float32, 1)
    context = Backends._prepare_device_execution_context(storage)
    nested_context = Backends._prepare_device_execution_context(
        Backends.compute_device(storage))
    caller_context = CUDA.context()
    caller_stream = CUDA.stream()
    caller_device = CUDA.device()

    @test context.context == caller_context
    @test nested_context.context == caller_context
    result = Backends._with_prepared_device_execution_context(context) do
        @test CUDA.context() == context.context
        @test CUDA.stream() === context.stream
        Backends._with_prepared_device_execution_context(nested_context) do
            @test CUDA.stream() === nested_context.stream
        end
        @test CUDA.stream() === context.stream
        :completed
    end
    @test result === :completed
    @test Backends._with_prepared_device_execution_context(
        _PreparedContextCallable(), context) === :called
    @test CUDA.context() == caller_context
    @test CUDA.stream() === caller_stream
    @test CUDA.device() == caller_device

    @test_throws ErrorException Backends._with_prepared_device_execution_context(context) do
        error("prepared context restoration")
    end
    @test CUDA.context() == caller_context
    @test CUDA.stream() === caller_stream
    @test CUDA.device() == caller_device

    # This entry/exit check excludes synchronization and resource retirement.
    # The owner must retain its context rather than resolve it on every call.
    _prepared_context_allocation(context)
    @test _prepared_context_allocation(context) == 0
end

function _prepared_wait_regression_kernel!(output, rounds)
    index = CUDA.threadIdx().x
    if index <= length(output)
        value = 0.0f0
        for _ in 1:rounds
            value = muladd(value, 0.999999f0, 0.000001f0)
        end
        output[index] = value
    end
    return
end

function _prepared_wait_regression_allocation(context)
    return @allocated Backends._synchronize_prepared_device_execution_context_blocking!(context)
end

function _prepared_wait_regression_launch!(output, rounds)
    CUDA.@cuda threads=4 blocks=1 _prepared_wait_regression_kernel!(output, rounds)
    return nothing
end

@testset "CUDA prepared blocking wait allocation" begin
    # CUDA 6.4.1 default nonblocking policy: this contract is specific to that
    # tested configuration, not a guarantee for every CUDA.jl allocation policy.
    rounds = 20_000_000
    output = CUDA.zeros(Float32, 4)
    context = Backends._prepare_device_execution_context(output)

    Backends._with_prepared_device_execution_context(context) do
        # Warm compilation and a completed-stream wait before measuring.
        _prepared_wait_regression_launch!(output, rounds)
        Backends._synchronize_prepared_device_execution_context_blocking!(context)
        @test all(isfinite, Array(output))
        @test _prepared_wait_regression_allocation(context) == 0

        # Warm the pending-kernel path, then submit another long-running tiny
        # kernel and measure only the owner's blocking wait.
        _prepared_wait_regression_launch!(output, rounds)
        Backends._synchronize_prepared_device_execution_context_blocking!(context)
        _prepared_wait_regression_launch!(output, rounds)
        @test !CUDA.isdone(context.stream)
        @test _prepared_wait_regression_allocation(context) == 0
        copied = Array(output)
        @test all(x -> isfinite(x) && x > 0.9f0, copied)
        @test all(==(first(copied)), copied)
    end

    # Exercise CUDA's quiet-window memory-refresh boundary. The sleep is outside
    # the allocation measurement; only the warmed wait after the idle interval
    # is measured.
    sleep(10.1)
    quiet_window_bytes = Backends._with_prepared_device_execution_context(context) do
        _prepared_wait_regression_allocation(context)
    end
    @test quiet_window_bytes == 0
end

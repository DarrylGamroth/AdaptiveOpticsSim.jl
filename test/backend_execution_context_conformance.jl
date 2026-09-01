function assert_prepared_device_execution_context_conformance(
    context::C,
    target::D,
    observe::F,
    expected_active::O,
) where {
    C<:Backends._AbstractPreparedDeviceExecutionContext,
    D<:Backends.AbstractComputeDevice,
    F,
    O,
}
    @test @inferred(
        Backends._prepared_device_execution_compute_device(context)) ==
        target
    caller_context = observe()
    active_context =
        Backends._with_prepared_device_execution_context(observe, context)
    @test active_context == expected_active
    @test observe() == caller_context

    injected = ErrorException("injected prepared-context failure")
    caught = try
        Backends._with_prepared_device_execution_context(context) do
            throw(injected)
        end
        nothing
    catch error
        error
    end
    @test caught === injected
    @test observe() == caller_context
    @test @inferred(
        Backends._synchronize_prepared_device_execution_context!(context)) ===
        nothing
    @test @inferred(
        Backends._synchronize_prepared_device_execution_context_blocking!(
            context,
        )) === nothing
    return nothing
end

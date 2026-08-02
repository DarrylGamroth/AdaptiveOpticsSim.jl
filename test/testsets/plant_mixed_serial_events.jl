@inline mixed_serial_event_host_values(values::Array) = values
@inline mixed_serial_event_host_values(values::HandoffTestArray) =
    values.storage

function mixed_serial_event_result_values(fixture, id::Symbol)
    return copy(mixed_serial_event_host_values(
        path_result(path_input_publication_test_path(
            fixture.partitions, id)).values))
end

function mixed_serial_event_product_values(fixture)
    products = acquisition_products(
        fixture.prepared, :alpha_camera)
    return copy(products.observation)
end

function mixed_serial_event_rng_probe(fixture)
    acquisition_rng = only(fixture.prepared.acquisitions).rng
    path_rngs = map(fixture.prepared.execution.paths) do path
        rng_stream_state(path.rngs, Val(:provider))
    end
    return (
        rand(copy(acquisition_rng), UInt64, 4),
        map(rng -> rand(copy(rng), UInt64, 4), path_rngs),
    )
end

function mixed_serial_event_step_allocations!(fixture)
    return @allocated step_plant_events!(
        fixture.prepared, fixture.state, fixture.workspace)
end

function mixed_serial_event_run_allocations!(fixture, stop)
    return @allocated run_plant_events_until!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        stop,
    )
end

function mixed_serial_event_command_admission_allocations!(
    fixture,
    command,
    timestamp,
)
    return @allocated admit_plant_command!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        command,
        timestamp,
    )
end

function mixed_serial_event_transaction_admission_allocations!(
    fixture,
    transaction,
    timestamp,
)
    return @allocated admit_plant_command_transaction!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        transaction,
        timestamp,
    )
end

function mixed_serial_event_capture_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

@testset "Mixed serial event composition preserves target-local products" begin
    reset_handoff_test_transfer_controls!()
    mixed = mixed_serial_event_fixture()
    host = mixed_serial_event_fixture(beta_target=HostComputeDevice())
    stop = PlantTimestamp(305_000_000)

    for name in (
        :PreparedMixedResourcePlantEventLoop,
        :MixedResourcePlantEventLoopState,
        :MixedResourcePlantEventLoopWorkspace,
    )
        @test Base.ispublic(Plant, name)
        @test !Base.isexported(Plant, name)
        @test !Base.isexported(AdaptiveOpticsSim, name)
    end
    @test mixed.prepared isa Plant.PreparedMixedResourcePlantEventLoop
    @test mixed.state isa Plant.MixedResourcePlantEventLoopState
    @test mixed.workspace isa Plant.MixedResourcePlantEventLoopWorkspace
    @test plant_event_path_count(mixed.prepared) == 2
    @test plant_event_acquisition_count(mixed.prepared) == 1
    @test plant_event_generator_count(mixed.prepared) ==
        event_generator_count(mixed.prepared.scheduler)
    @test plant_event_command_endpoint_count(mixed.prepared) == 0
    @test plant_event_controllable_optic_count(mixed.prepared) == 0
    @test plant_event_sampled_aberration_count(mixed.prepared) == 0
    @test plant_event_autonomous_optic_count(mixed.prepared) == 0
    @test atmosphere_authority_target(mixed.prepared) == HostComputeDevice()
    @test command_authority_target(mixed.prepared) == HostComputeDevice()

    @test run_plant_events_until!(
        mixed.prepared, mixed.state, mixed.workspace, stop) > 0
    @test run_plant_events_until!(
        host.prepared, host.state, host.workspace, stop) > 0

    @test isconcretetype(typeof(mixed.prepared))
    @test isconcretetype(typeof(mixed.state))
    @test isconcretetype(typeof(mixed.workspace))
    @test mixed.prepared.acquisitions isa
        Memory{Plant._AbstractPreparedMixedSerialEventAcquisition}
    @test mixed.state.acquisitions isa
        Memory{Plant._AcquisitionEventLifecycleState}
    @test length(mixed.prepared.paths) ==
        length(mixed.state.path_sampled) ==
        length(mixed.workspace.execution.due_paths)
    @test length(mixed.prepared.acquisitions) ==
        length(mixed.state.acquisitions) ==
        length(mixed.state.product_sequences)
    @test mixed.state.phase == Plant._MixedSerialEventLoopIdle
    @test host.state.phase == Plant._MixedSerialEventLoopIdle
    @test scheduler_timestamp(mixed.state.scheduler) ==
        scheduler_timestamp(host.state.scheduler)
    @test atmosphere_timeline(prepared_atmosphere(
        prepared_atmosphere_authority(mixed.partitions))).sequence ==
        atmosphere_timeline(prepared_atmosphere(
            prepared_atmosphere_authority(host.partitions))).sequence
    @test mixed_serial_event_result_values(mixed, :alpha) ≈
        mixed_serial_event_result_values(host, :alpha)
    @test mixed_serial_event_result_values(mixed, :beta) ≈
        mixed_serial_event_result_values(host, :beta)
    @test mixed_serial_event_product_values(mixed) ≈
        mixed_serial_event_product_values(host)
    @test acquisition_product_sequence(
        mixed.prepared, mixed.state, :alpha_camera) == UInt64(2)
    @test acquisition_product_sequence(
        mixed.prepared, mixed.state, :alpha_camera) ==
        acquisition_product_sequence(
            host.prepared, host.state, :alpha_camera)
    @test acquisition_product_ready_timestamp(
        mixed.prepared, mixed.state, :alpha_camera) ==
        PlantTimestamp(300_000_000)
    @test compute_device(acquisition_products(
        mixed.prepared, :alpha_camera).observation) == HostComputeDevice()
    @test mixed_serial_event_rng_probe(mixed) ==
        mixed_serial_event_rng_probe(host)
    @test HANDOFF_TEST_CONTEXT_ENTRIES[] > zero(UInt64)
    @test HANDOFF_TEST_ACTIVE_DEVICE[] == zero(UInt32)
end

@testset "Mixed serial event composition closes detector lifecycle modes" begin
    for (mode, lifecycle_type) in (
        (:global, PreparedGlobalShutterAcquisition),
        (:rolling, PreparedRollingShutterAcquisition),
        (:frame_transfer, PreparedFrameTransferAcquisition),
        (:up_the_ramp, PreparedGlobalShutterAcquisition),
    )
        reset_handoff_test_transfer_controls!()
        fixture = mixed_serial_event_fixture(detector_mode=mode)
        @test run_plant_events_until!(
            fixture.prepared,
            fixture.state,
            fixture.workspace,
            PlantTimestamp(120_000_000),
        ) > 0
        acquisition = only(fixture.prepared.acquisitions)
        @test acquisition.lifecycle isa lifecycle_type
        @test acquisition_product_sequence(
            fixture.prepared, fixture.state, :alpha_camera) == UInt64(1)
        @test all(isfinite, mixed_serial_event_product_values(fixture))
        @test fixture.state.phase == Plant._MixedSerialEventLoopIdle
        if mode === :rolling
            @test rolling_opened_band_count(
                only(fixture.state.acquisitions)) == 2
            @test rolling_closed_band_count(
                only(fixture.state.acquisitions)) == 2
        elseif mode === :frame_transfer
            @test frame_transfer_product_sequence(
                only(fixture.state.acquisitions)) == UInt64(1)
        elseif mode === :up_the_ramp
            @test nondestructive_read_count(acquisition.lifecycle) == 3
        end
    end
end

@testset "Mixed serial trigger delivery starts a target-local detector" begin
    reset_handoff_test_transfer_controls!()
    fixture = mixed_serial_event_fixture(triggered=true)
    @test run_plant_events_until!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        PlantTimestamp(110_000_000),
    ) > 0
    @test acquisition_product_sequence(
        fixture.prepared, fixture.state, :alpha_camera) == UInt64(1)
    @test pending_trigger_delivery_count(fixture.state.trigger) == 0
    @test fixture.state.phase == Plant._MixedSerialEventLoopIdle
end

@testset "Mixed serial command scheduling precedes unequal-rate paths" begin
    reset_handoff_test_transfer_controls!()
    fixture = mixed_serial_event_fixture(with_command=true)
    timestamp = PlantTimestamp(50_000_000)
    command = PlantCommand(fixture.schema, 1, timestamp, 2.0)
    admission = admit_plant_command!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        command,
        zero(PlantTimestamp),
    )
    @test command_admission_status(admission) == CommandAdmittedPending
    @test run_plant_events_until!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        PlantTimestamp(160_000_000),
    ) > 0
    @test fixture.state.execution.fanout.authority.publication_sequences ==
        UInt64[1]
    @test effective_command(
        fixture.prepared,
        fixture.state,
        :effective_command_route_endpoint,
    ) == 3.0
    @test command_disposition_count(
        fixture.workspace) == 1
    disposition = command_disposition(
        fixture.workspace, 1)
    @test command_terminal_kind(disposition) == AppliedCommand
    clear_command_dispositions!(fixture.workspace)
end

@testset "Mixed serial command silence uses an independent authority target" begin
    reset_handoff_test_transfer_controls!()
    silence_policy = CommandSilencePolicy(
        ApplySafeCommand,
        AgeFromApplication;
        timeout=PlantDuration(50_000_000),
    )
    fixture = mixed_serial_event_fixture(
        beta_target=HostComputeDevice(),
        with_command=true,
        command_authority_target=HANDOFF_TEST_ACCELERATOR,
        silence_policy=silence_policy,
        safe_command=-2.0,
    )

    @test atmosphere_authority_target(fixture.prepared) ==
        HostComputeDevice()
    @test command_authority_target(fixture.prepared) ==
        HANDOFF_TEST_ACCELERATOR
    @test all(
        partition -> compute_device(partition) == HostComputeDevice(),
        prepared_partitions(fixture.partitions),
    )
    command_timestamp = PlantTimestamp(10_000_000)
    admission = admit_plant_command!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        PlantCommand(fixture.schema, 1, command_timestamp, 2.0),
        zero(PlantTimestamp),
    )
    @test command_admission_status(admission) == CommandAdmittedPending
    @test run_plant_events_until!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        PlantTimestamp(160_000_000),
    ) > 0
    @test effective_command(
        fixture.prepared,
        fixture.state,
        :effective_command_route_endpoint,
    ) == -2.0
    @test fixture.state.execution.fanout.authority.publication_sequences ==
        UInt64[2]
    @test fixture.state.execution.fanout.authority.publication_timestamps ==
        PlantTimestamp[PlantTimestamp(60_000_000)]
    @test command_disposition_count(fixture.workspace) == 1
    @test command_terminal_kind(command_disposition(
        fixture.workspace, 1)) == AppliedCommand
    @test fixture.state.phase == Plant._MixedSerialEventLoopIdle
    @test HANDOFF_TEST_ACTIVE_DEVICE[] == zero(UInt32)
end

@testset "Mixed serial event scheduling preserves atomic transactions" begin
    reset_handoff_test_transfer_controls!()
    fixture = mixed_serial_transaction_event_fixture()
    timestamp = PlantTimestamp(50_000_000)
    transaction = PlantCommandTransaction(
        PlantCommand(fixture.beta_schema, 1, timestamp, 4.0),
        PlantCommand(fixture.alpha_schema, 1, timestamp, 2.0),
    )
    admission = admit_plant_command_transaction!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        transaction,
        zero(PlantTimestamp),
    )

    @test command_admission_status(admission) == CommandAdmittedPending
    @test command_transaction_id(admission) == UInt64(1)
    @test command_transaction_member_count(admission) == UInt32(2)
    @test run_plant_events_until!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        PlantTimestamp(160_000_000),
    ) > 0
    @test effective_command(
        fixture.prepared,
        fixture.state,
        command_endpoint_id(fixture.alpha_schema),
    ) == 3.0
    @test effective_command(
        fixture.prepared,
        fixture.state,
        command_endpoint_id(fixture.beta_schema),
    ) == 5.0
    @test fixture.state.execution.fanout.authority.publication_sequences ==
        UInt64[1, 1]
    @test command_disposition_count(fixture.workspace) == 2
    @test all(
        index -> command_terminal_kind(command_disposition(
            fixture.workspace, index)) == AppliedCommand,
        1:2,
    )
    clear_command_dispositions!(fixture.workspace)
    @test HANDOFF_TEST_ACTIVE_DEVICE[] == zero(UInt32)
end

@testset "Mixed serial event failure disarms the coordinator" begin
    reset_handoff_test_transfer_controls!()
    fixture = mixed_serial_event_fixture()
    initial_alpha = mixed_serial_event_result_values(fixture, :alpha)
    initial_beta = mixed_serial_event_result_values(fixture, :beta)
    HANDOFF_TEST_FAIL_COMPLETION[] = true
    error = mixed_serial_event_capture_error() do
        step_plant_events!(
            fixture.prepared, fixture.state, fixture.workspace)
    end
    HANDOFF_TEST_FAIL_COMPLETION[] = false

    @test error isa PlantScheduleError
    error isa PlantScheduleError &&
        @test error.reason === :path_input_completion
    @test fixture.state.phase == Plant._MixedSerialEventLoopFailed
    @test fixture.state.execution.phase == Plant._MixedSerialExecutionFailed
    @test mixed_serial_event_result_values(fixture, :alpha) == initial_alpha
    @test mixed_serial_event_result_values(fixture, :beta) == initial_beta
    @test acquisition_product_sequence(
        fixture.prepared, fixture.state, :alpha_camera) == zero(UInt64)
    @test HANDOFF_TEST_ACTIVE_DEVICE[] == zero(UInt32)

    retry = mixed_serial_event_capture_error() do
        step_plant_events!(
            fixture.prepared, fixture.state, fixture.workspace)
    end
    @test retry isa PlantScheduleError
    retry isa PlantScheduleError &&
        @test retry.reason === :event_loop_failed
    reset_handoff_test_transfer_controls!()
end

@testset "Mixed serial CPU event paths are inferred and zero allocation" begin
    reset_handoff_test_transfer_controls!()
    inferred_fixture = mixed_serial_event_fixture()
    run_plant_events_until!(
        inferred_fixture.prepared,
        inferred_fixture.state,
        inferred_fixture.workspace,
        PlantTimestamp(505_000_000),
    )
    @test @inferred(Nothing, step_plant_events!(
        inferred_fixture.prepared,
        inferred_fixture.state,
        inferred_fixture.workspace,
    )) == PlantTimestamp(600_000_000)

    detector_modes = (:global, :rolling, :frame_transfer, :up_the_ramp)
    lifecycle_run_allocations = map(detector_modes) do mode
        warm_fixture = mixed_serial_event_fixture(
            beta_target=HostComputeDevice(),
            detector_mode=mode,
        )
        run_plant_events_until!(
            warm_fixture.prepared,
            warm_fixture.state,
            warm_fixture.workspace,
            PlantTimestamp(120_000_000),
        )
        measured_fixture = mixed_serial_event_fixture(
            beta_target=HostComputeDevice(),
            detector_mode=mode,
        )
        mixed_serial_event_run_allocations!(
            measured_fixture, PlantTimestamp(120_000_000))
    end
    @info(
        "Mixed serial CPU detector-lifecycle allocations",
        detector_modes,
        lifecycle_run_allocations,
    )

    warm_trigger = mixed_serial_event_fixture(
        beta_target=HostComputeDevice(), triggered=true)
    run_plant_events_until!(
        warm_trigger.prepared,
        warm_trigger.state,
        warm_trigger.workspace,
        PlantTimestamp(120_000_000),
    )
    measured_trigger = mixed_serial_event_fixture(
        beta_target=HostComputeDevice(), triggered=true)
    trigger_run_allocations = mixed_serial_event_run_allocations!(
        measured_trigger, PlantTimestamp(120_000_000))

    command_timestamp = PlantTimestamp(50_000_000)
    warm_command = mixed_serial_event_fixture(
        beta_target=HostComputeDevice(), with_command=true)
    admit_plant_command!(
        warm_command.prepared,
        warm_command.state,
        warm_command.workspace,
        PlantCommand(warm_command.schema, 1, command_timestamp, 2.0),
        zero(PlantTimestamp),
    )
    run_plant_events_until!(
        warm_command.prepared,
        warm_command.state,
        warm_command.workspace,
        PlantTimestamp(160_000_000),
    )
    measured_command = mixed_serial_event_fixture(
        beta_target=HostComputeDevice(), with_command=true)
    command_admission_allocations =
        mixed_serial_event_command_admission_allocations!(
            measured_command,
            PlantCommand(
                measured_command.schema, 1, command_timestamp, 2.0),
            zero(PlantTimestamp),
        )
    command_run_allocations = mixed_serial_event_run_allocations!(
        measured_command, PlantTimestamp(160_000_000))

    silence_policy = CommandSilencePolicy(
        ApplySafeCommand,
        AgeFromApplication;
        timeout=PlantDuration(50_000_000),
    )
    warm_silence = mixed_serial_event_fixture(
        beta_target=HostComputeDevice(),
        with_command=true,
        silence_policy=silence_policy,
        safe_command=-2.0,
    )
    admit_plant_command!(
        warm_silence.prepared,
        warm_silence.state,
        warm_silence.workspace,
        PlantCommand(
            warm_silence.schema, 1, PlantTimestamp(10_000_000), 2.0),
        zero(PlantTimestamp),
    )
    run_plant_events_until!(
        warm_silence.prepared,
        warm_silence.state,
        warm_silence.workspace,
        PlantTimestamp(160_000_000),
    )
    measured_silence = mixed_serial_event_fixture(
        beta_target=HostComputeDevice(),
        with_command=true,
        silence_policy=silence_policy,
        safe_command=-2.0,
    )
    silence_admission_allocations =
        mixed_serial_event_command_admission_allocations!(
            measured_silence,
            PlantCommand(
                measured_silence.schema,
                1,
                PlantTimestamp(10_000_000),
                2.0,
            ),
            zero(PlantTimestamp),
        )
    silence_run_allocations = mixed_serial_event_run_allocations!(
        measured_silence, PlantTimestamp(160_000_000))

    warm_transaction = mixed_serial_transaction_event_fixture(
        beta_target=HostComputeDevice())
    warm_transaction_value = PlantCommandTransaction(
        PlantCommand(
            warm_transaction.beta_schema, 1, command_timestamp, 4.0),
        PlantCommand(
            warm_transaction.alpha_schema, 1, command_timestamp, 2.0),
    )
    admit_plant_command_transaction!(
        warm_transaction.prepared,
        warm_transaction.state,
        warm_transaction.workspace,
        warm_transaction_value,
        zero(PlantTimestamp),
    )
    run_plant_events_until!(
        warm_transaction.prepared,
        warm_transaction.state,
        warm_transaction.workspace,
        PlantTimestamp(160_000_000),
    )
    measured_transaction = mixed_serial_transaction_event_fixture(
        beta_target=HostComputeDevice())
    measured_transaction_value = PlantCommandTransaction(
        PlantCommand(
            measured_transaction.beta_schema, 1, command_timestamp, 4.0),
        PlantCommand(
            measured_transaction.alpha_schema, 1, command_timestamp, 2.0),
    )
    transaction_admission_allocations =
        mixed_serial_event_transaction_admission_allocations!(
            measured_transaction,
            measured_transaction_value,
            zero(PlantTimestamp),
        )
    transaction_run_allocations = mixed_serial_event_run_allocations!(
        measured_transaction, PlantTimestamp(160_000_000))

    cpu_event_run_allocations = (
        trigger=trigger_run_allocations,
        command=command_run_allocations,
        silence=silence_run_allocations,
        transaction=transaction_run_allocations,
    )
    cpu_admission_allocations = (
        command=command_admission_allocations,
        silence=silence_admission_allocations,
        transaction=transaction_admission_allocations,
    )
    @info "Mixed serial CPU event-run allocations" cpu_event_run_allocations
    @info "Mixed serial CPU admission allocations" cpu_admission_allocations

    allocation_fixture = mixed_serial_event_fixture()
    run_plant_events_until!(
        allocation_fixture.prepared,
        allocation_fixture.state,
        allocation_fixture.workspace,
        PlantTimestamp(505_000_000),
    )
    allocations = mixed_serial_event_step_allocations!(allocation_fixture)
    @info "Mixed serial event steady-state allocations" allocations
    if !coverage_instrumented()
        # Preparation, compilation, diagnostics, and exceptional paths are
        # outside the successful steady-state CPU HIL allocation contract.
        @test all(iszero, lifecycle_run_allocations)
        @test all(iszero, values(cpu_event_run_allocations))
        @test all(iszero, values(cpu_admission_allocations))
        @test allocations == 0
    end
end

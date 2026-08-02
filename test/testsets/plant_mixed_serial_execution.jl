@inline mixed_serial_host_values(values::Array) = values
@inline mixed_serial_host_values(values::HandoffTestArray) = values.storage

function mixed_serial_fixture(;
    command::Bool=false,
)
    partitions, schema = if command
        effective_command_route_test_partitions(dimensions=())
    else
        (path_input_publication_test_partitions(), nothing)
    end
    prepared = Plant._prepare_mixed_serial_execution(partitions)
    state = Plant._prepare_mixed_serial_execution_state(prepared)
    workspace = Plant._prepare_mixed_serial_execution_workspace(
        prepared, state)
    return (; partitions, schema, prepared, state, workspace)
end

function mixed_serial_oracle_values(
    partitions::PreparedPlantPartitions,
    id::Symbol;
    surface=0.0,
)
    oracle = path_input_publication_test_oracle(partitions, id)
    atmosphere = prepared_atmosphere(
        prepared_atmosphere_authority(partitions))
    materialize_path_input!(oracle, current_epoch(atmosphere))
    execute_path!(oracle)
    values = copy(mixed_serial_host_values(path_result(oracle).values))
    @. values += surface
    return values
end

function mixed_serial_result_values(
    partitions::PreparedPlantPartitions,
    id::Symbol,
)
    path = path_input_publication_test_path(partitions, id)
    return mixed_serial_host_values(path_result(path).values)
end

function mixed_serial_capture_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

function mixed_serial_allocation_sample!(fixture, ids, timestamp)
    return @allocated Plant._execute_mixed_serial_paths!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        ids,
        timestamp,
    )
end

@testset "Prepared mixed-target serial optical execution" begin
    reset_handoff_test_transfer_controls!()
    fixture = mixed_serial_fixture()
    prepared = fixture.prepared
    state = fixture.state
    workspace = fixture.workspace

    @test prepared.paths isa Memory{Plant._AbstractPreparedMixedSerialPath}
    @test map(path -> path.id, prepared.paths) ==
        OpticalPathID[OpticalPathID(:alpha), OpticalPathID(:beta)]
    @test workspace.due_paths isa Memory{Bool}
    @test length(workspace.paths) == 2
    @test all(path -> path.controllable_optics isa Tuple &&
        path.command_dependencies isa Tuple, prepared.paths)
    @test isconcretetype(typeof(prepared))
    @test isconcretetype(typeof(state))
    @test isconcretetype(typeof(workspace))

    @test Plant._execute_mixed_serial_paths!(
        prepared,
        state,
        workspace,
        (:beta, OpticalPathID(:alpha)),
        PlantTimestamp(0),
    ) === nothing
    atmosphere = prepared_atmosphere(
        prepared_atmosphere_authority(fixture.partitions))
    @test atmosphere_timeline(atmosphere).sequence == UInt64(1)
    @test state.phase == Plant._MixedSerialExecutionIdle
    @test state.last_timestamp == PlantTimestamp(0)
    @test state.has_timestamp
    @test all(!, workspace.due_paths)
    @test HANDOFF_TEST_ACTIVE_DEVICE[] == zero(UInt32)
    @test mixed_serial_result_values(fixture.partitions, :alpha) ≈
        mixed_serial_oracle_values(fixture.partitions, :alpha)
    @test mixed_serial_result_values(fixture.partitions, :beta) ≈
        mixed_serial_oracle_values(fixture.partitions, :beta)

    @test Plant._execute_mixed_serial_paths!(
        prepared,
        state,
        workspace,
        (:alpha, :beta),
        PlantTimestamp(10),
    ) === nothing
    @test atmosphere_timeline(atmosphere).sequence == UInt64(2)
    @test state.last_timestamp == PlantTimestamp(10)
end

@testset "Mixed serial preparation selects a canonical scheduled subset" begin
    reset_handoff_test_transfer_controls!()
    partitions = path_input_publication_test_partitions()
    prepared = Plant._prepare_mixed_serial_execution(
        partitions, (:beta,))
    state = Plant._prepare_mixed_serial_execution_state(prepared)
    workspace = Plant._prepare_mixed_serial_execution_workspace(
        prepared, state)

    @test prepared.path_ids == OpticalPathID[OpticalPathID(:beta)]
    @test only(prepared.paths).id == OpticalPathID(:beta)
    workspace.due_paths[1] = true
    @test Plant._execute_selected_mixed_serial_paths!(
        prepared, state, workspace, PlantTimestamp(0)) === nothing
    @test mixed_serial_result_values(partitions, :beta) ≈
        mixed_serial_oracle_values(partitions, :beta)

    duplicate = mixed_serial_capture_error() do
        Plant._prepare_mixed_serial_execution(partitions, (:beta, :beta))
    end
    @test duplicate isa PlantPreparationError
    duplicate isa PlantPreparationError &&
        @test duplicate.reason === :duplicate_path

    unknown = mixed_serial_capture_error() do
        Plant._prepare_mixed_serial_execution(partitions, (:missing,))
    end
    @test unknown isa PlantPreparationError
    unknown isa PlantPreparationError &&
        @test unknown.reason === :unknown_path
end

@testset "Mixed serial preflight is mutation-free" begin
    reset_handoff_test_transfer_controls!()
    fixture = mixed_serial_fixture()
    atmosphere = prepared_atmosphere(
        prepared_atmosphere_authority(fixture.partitions))

    for (due, timestamp, reason) in (
        ((), PlantTimestamp(0), :empty_due_paths),
        ((:alpha, :alpha), PlantTimestamp(0), :duplicate_due_path),
        ((:unknown,), PlantTimestamp(0), :unknown_path),
    )
        error = mixed_serial_capture_error() do
            Plant._execute_mixed_serial_paths!(
                fixture.prepared,
                fixture.state,
                fixture.workspace,
                due,
                timestamp,
            )
        end
        @test error isa PlantScheduleError
        error isa PlantScheduleError && @test error.reason === reason
        @test fixture.state.phase == Plant._MixedSerialExecutionIdle
        @test atmosphere_timeline(atmosphere).sequence == UInt64(0)
        @test all(!, fixture.workspace.due_paths)
    end

    Plant._execute_mixed_serial_paths!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        (:alpha, :beta),
        PlantTimestamp(1),
    )
    error = mixed_serial_capture_error() do
        Plant._execute_mixed_serial_paths!(
            fixture.prepared,
            fixture.state,
            fixture.workspace,
            (:alpha,),
            PlantTimestamp(1),
        )
    end
    @test error isa PlantScheduleError
    error isa PlantScheduleError &&
        @test error.reason === :nonmonotonic_batch_time
    @test atmosphere_timeline(atmosphere).sequence == UInt64(1)
    @test fixture.state.phase == Plant._MixedSerialExecutionIdle
end

@testset "Mixed serial command dependency readiness" begin
    reset_handoff_test_transfer_controls!()
    fixture = mixed_serial_fixture(command=true)
    endpoint_id = command_endpoint_id(fixture.schema)
    lane = only(fixture.prepared.fanout.lanes)
    lane_state = only(fixture.state.fanout.lanes)
    disposition = only(fixture.workspace.fanout.lanes).disposition
    timestamp = PlantTimestamp(10)

    admission = admit_plant_command!(
        disposition,
        lane.authority_endpoint.endpoint,
        lane_state.endpoint,
        PlantCommand(fixture.schema, 1, timestamp, 2.0),
        zero(PlantTimestamp),
    )
    @test command_admission_status(admission) == CommandAdmittedPending
    before = atmosphere_timeline(prepared_atmosphere(
        prepared_atmosphere_authority(fixture.partitions))).sequence
    error = mixed_serial_capture_error() do
        Plant._execute_mixed_serial_paths!(
            fixture.prepared,
            fixture.state,
            fixture.workspace,
            (:alpha,),
            timestamp,
        )
    end
    @test error isa PlantScheduleError
    error isa PlantScheduleError &&
        @test error.reason === :command_not_applied
    @test atmosphere_timeline(prepared_atmosphere(
        prepared_atmosphere_authority(fixture.partitions))).sequence == before
    @test fixture.state.phase == Plant._MixedSerialExecutionIdle

    Plant._apply_next_command_fanout!(
        fixture.prepared.fanout,
        fixture.state.fanout,
        fixture.workspace.fanout,
        endpoint_id,
        timestamp,
    )
    clear_command_dispositions!(disposition)
    @test Plant._execute_mixed_serial_paths!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        (:beta, :alpha),
        timestamp,
    ) === nothing
    @test fixture.state.fanout.authority.publication_sequences == UInt64[1]
    @test mixed_serial_result_values(fixture.partitions, :alpha) ≈
        mixed_serial_oracle_values(
            fixture.partitions, :alpha; surface=3.0)
    @test mixed_serial_result_values(fixture.partitions, :beta) ≈
        mixed_serial_oracle_values(
            fixture.partitions, :beta; surface=3.0)
end

@testset "Mixed serial post-mutation failure is fail-stop" begin
    reset_handoff_test_transfer_controls!()
    fixture = mixed_serial_fixture()
    HANDOFF_TEST_FAIL_COMPLETION[] = true
    error = mixed_serial_capture_error() do
        Plant._execute_mixed_serial_paths!(
            fixture.prepared,
            fixture.state,
            fixture.workspace,
            (:alpha, :beta),
            PlantTimestamp(0),
        )
    end
    @test error isa PlantScheduleError
    error isa PlantScheduleError &&
        @test error.reason === :path_input_completion
    @test fixture.state.phase == Plant._MixedSerialExecutionFailed
    @test atmosphere_timeline(prepared_atmosphere(
        prepared_atmosphere_authority(fixture.partitions))).sequence ==
        UInt64(1)
    @test_throws PlantScheduleError Plant._execute_mixed_serial_paths!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        (:alpha,),
        PlantTimestamp(1),
    )
    reset_handoff_test_transfer_controls!()
end

@testset "Mixed serial inference and steady-state allocation" begin
    reset_handoff_test_transfer_controls!()
    fixture = mixed_serial_fixture()
    ids = (OpticalPathID(:alpha), OpticalPathID(:beta))
    @test @inferred(Plant._execute_mixed_serial_paths!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        ids,
        PlantTimestamp(0),
    )) === nothing
    allocations = mixed_serial_allocation_sample!(
        fixture, ids, PlantTimestamp(1))
    if !coverage_instrumented()
        @test allocations == 0
    end

    commanded = mixed_serial_fixture(command=true)
    endpoint_id = command_endpoint_id(commanded.schema)
    lane = only(commanded.prepared.fanout.lanes)
    lane_state = only(commanded.state.fanout.lanes)
    disposition = only(commanded.workspace.fanout.lanes).disposition
    admission = admit_plant_command!(
        disposition,
        lane.authority_endpoint.endpoint,
        lane_state.endpoint,
        PlantCommand(commanded.schema, 1, PlantTimestamp(0), 2.0),
        zero(PlantTimestamp),
    )
    @test command_admission_status(admission) == CommandAdmittedReady
    Plant._apply_next_command_fanout!(
        commanded.prepared.fanout,
        commanded.state.fanout,
        commanded.workspace.fanout,
        endpoint_id,
        PlantTimestamp(0),
    )
    clear_command_dispositions!(disposition)
    @test @inferred(Plant._execute_mixed_serial_paths!(
        commanded.prepared,
        commanded.state,
        commanded.workspace,
        ids,
        PlantTimestamp(0),
    )) === nothing
    command_allocations = mixed_serial_allocation_sample!(
        commanded, ids, PlantTimestamp(1))
    if !coverage_instrumented()
        @test command_allocations == 0
    end
end

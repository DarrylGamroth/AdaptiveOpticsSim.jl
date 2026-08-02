function command_fanout_test_fixture()
    partitions, schema = effective_command_route_test_partitions()
    prepared = Plant._prepare_command_fanout(partitions)
    state = Plant._prepare_command_fanout_state(prepared)
    workspace = Plant._prepare_command_fanout_workspace(prepared, state)
    authority = prepared_command_authority(partitions)
    endpoint_id = command_endpoint_id(schema)
    endpoint = prepared_command_authority_endpoint(authority, endpoint_id)
    endpoint_state = command_authority_endpoint_state(
        authority, state.authority, endpoint_id)
    disposition = command_authority_disposition_workspace(
        authority, workspace.authority, endpoint_id)
    return (
        partitions,
        schema,
        prepared,
        state,
        workspace,
        endpoint,
        endpoint_state,
        disposition,
    )
end

function command_fanout_test_silence_fixture(
    action::CommandSilenceAction;
    safe_command=nothing,
)
    silence_policy = action == HoldLastCommand ?
        CommandSilencePolicy() :
        CommandSilencePolicy(
            action,
            AgeFromApplication;
            timeout=PlantDuration(10),
        )
    schema = PlantCommandSchema(
        Float64,
        (2,);
        id=:command_fanout_silence_schema,
        version=1,
        endpoint=:effective_command_route_endpoint,
        units=:metre,
        sign_convention=:positive_surface_increases_opd,
        basis=CommandBasis(:actuator, :command_fanout_silence_basis),
        basis_revision=1,
        semantics=IncrementalCommand,
        bounds=UnboundedCommandValues(),
        value_policy=CommandValuePolicy(),
        sequence_policy=CommandSequencePolicy(),
        effective_time_policy=CommandEffectiveTimePolicy(),
        silence_policy,
    )
    optic = ControllableOpticDefinition(
        :effective_command_route_optic,
        EffectiveCommandRouteTestOpticModel(),
        (schema,);
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
    )
    definition = path_input_publication_test_definition(
        controllable_optics=(optic,))
    configurations = (
        CommandEndpointConfiguration(
            command_endpoint_id(schema),
            [1.0, 1.0];
            capacity=2,
            safe_command,
        ),)
    assignment = resolve_plant_partition_assignment(
        definition,
        HostComputeDevice(),
        :alpha => HostComputeDevice(),
        :beta => HANDOFF_TEST_ACCELERATOR,
    )
    partitions = with_path_input_publication_cold_scalar_indexing() do
        prepare_plant_partitions(
            definition,
            assignment;
            run_seed=0x219,
            command_authority_target=HostComputeDevice(),
            command_endpoints=configurations,
        )
    end
    prepared = Plant._prepare_command_fanout(partitions)
    state = Plant._prepare_command_fanout_state(prepared)
    workspace = Plant._prepare_command_fanout_workspace(prepared, state)
    authority = prepared_command_authority(partitions)
    endpoint_id = command_endpoint_id(schema)
    endpoint = prepared_command_authority_endpoint(authority, endpoint_id)
    endpoint_state = command_authority_endpoint_state(
        authority, state.authority, endpoint_id)
    disposition = command_authority_disposition_workspace(
        authority, workspace.authority, endpoint_id)
    return (
        partitions,
        schema,
        prepared,
        state,
        workspace,
        endpoint,
        endpoint_state,
        disposition,
    )
end

function command_fanout_test_target_values(partitions)
    endpoint_id = CommandEndpointID(:effective_command_route_endpoint)
    values = Any[]
    physical = Any[]
    for partition in prepared_partitions(partitions)
        owner = only(target_local_controllable_optic_owners(partition))
        push!(values, effective_command(
            target_local_command_endpoint_state(owner, endpoint_id)))
        push!(physical, getfield(
            target_local_controllable_optic_state(owner), :physical).active)
    end
    return values, physical
end

function command_fanout_test_admit!(endpoint, endpoint_state, disposition,
    schema, sequence, effective_timestamp, payload, admission_timestamp)
    return admit_plant_command!(
        disposition,
        endpoint,
        endpoint_state,
        PlantCommand(schema, sequence, effective_timestamp, payload),
        admission_timestamp,
    )
end

function command_fanout_test_replica_values(partitions)
    values, physical = command_fanout_test_target_values(partitions)
    return (
        map(effective_command_route_test_host_value, values),
        map(effective_command_route_test_host_value, physical),
    )
end

function command_fanout_test_no_publication(partitions)
    endpoint_id = CommandEndpointID(:effective_command_route_endpoint)
    return all(prepared_partitions(partitions)) do partition
        owner = only(target_local_controllable_optic_owners(partition))
        !has_effective_command_publication(
            target_local_command_endpoint_state(owner, endpoint_id))
    end
end

function command_fanout_test_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

function command_fanout_transaction_schema(
    id::Symbol,
    endpoint::Symbol;
    effective_time_policy=CommandEffectiveTimePolicy(),
)
    return PlantCommandSchema(
        Float64,
        (2,);
        id,
        version=1,
        endpoint,
        units=:metre,
        sign_convention=:positive_surface_increases_opd,
        basis=CommandBasis(:actuator, Symbol(id, :_basis)),
        basis_revision=1,
        semantics=IncrementalCommand,
        bounds=UnboundedCommandValues(),
        value_policy=CommandValuePolicy(),
        sequence_policy=CommandSequencePolicy(),
        effective_time_policy,
        silence_policy=CommandSilencePolicy(),
    )
end

function command_fanout_transaction_fixture(
    ;
    alpha_capacity::Integer=2,
    beta_capacity::Integer=2,
    alpha_effective_time_policy=CommandEffectiveTimePolicy(),
    beta_effective_time_policy=CommandEffectiveTimePolicy(),
)
    alpha_schema = command_fanout_transaction_schema(
        :command_fanout_transaction_alpha_schema,
        :command_fanout_transaction_alpha,
        effective_time_policy=alpha_effective_time_policy,
    )
    beta_schema = command_fanout_transaction_schema(
        :command_fanout_transaction_beta_schema,
        :command_fanout_transaction_beta,
        effective_time_policy=beta_effective_time_policy,
    )
    alpha_optic = ControllableOpticDefinition(
        :command_fanout_transaction_alpha_optic,
        EffectiveCommandRouteTestOpticModel(),
        (alpha_schema,);
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
    )
    beta_optic = ControllableOpticDefinition(
        :command_fanout_transaction_beta_optic,
        EffectiveCommandRouteTestOpticModel(),
        (beta_schema,);
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
    )
    definition = path_input_publication_test_definition(
        controllable_optics=(alpha_optic, beta_optic))
    configurations = (
        CommandEndpointConfiguration(
            command_endpoint_id(alpha_schema),
            [1.0, 1.0];
            capacity=alpha_capacity,
        ),
        CommandEndpointConfiguration(
            command_endpoint_id(beta_schema),
            [1.0, 1.0];
            capacity=beta_capacity,
        ),
    )
    assignment = resolve_plant_partition_assignment(
        definition,
        HostComputeDevice(),
        :alpha => HostComputeDevice(),
        :beta => HANDOFF_TEST_ACCELERATOR,
    )
    partitions = with_path_input_publication_cold_scalar_indexing() do
        prepare_plant_partitions(
            definition,
            assignment;
            run_seed=0x219,
            command_authority_target=HostComputeDevice(),
            command_endpoints=configurations,
        )
    end
    prepared = Plant._prepare_command_fanout(partitions)
    state = Plant._prepare_command_fanout_state(prepared)
    workspace = Plant._prepare_command_fanout_workspace(prepared, state)
    authority = prepared_command_authority(partitions)
    alpha_id = command_endpoint_id(alpha_schema)
    beta_id = command_endpoint_id(beta_schema)
    return (
        partitions=partitions,
        alpha_schema=alpha_schema,
        beta_schema=beta_schema,
        alpha_id=alpha_id,
        beta_id=beta_id,
        prepared=prepared,
        state=state,
        workspace=workspace,
        alpha_endpoint=prepared_command_authority_endpoint(authority, alpha_id),
        beta_endpoint=prepared_command_authority_endpoint(authority, beta_id),
        alpha_state=command_authority_endpoint_state(
            authority, state.authority, alpha_id),
        beta_state=command_authority_endpoint_state(
            authority, state.authority, beta_id),
        alpha_disposition=command_authority_disposition_workspace(
            authority, workspace.authority, alpha_id),
        beta_disposition=command_authority_disposition_workspace(
            authority, workspace.authority, beta_id),
    )
end

function command_fanout_transaction_target_owner(
    partition::PreparedTargetPartition,
    endpoint_id::CommandEndpointID,
)
    return only(owner for owner in target_local_controllable_optic_owners(partition)
        if any(endpoint -> command_endpoint_id(endpoint) == endpoint_id,
            owner.prepared.endpoints))
end

function command_fanout_transaction_replica_values(
    partitions,
    endpoint_id::CommandEndpointID,
)
    values = Any[]
    physical = Any[]
    for partition in prepared_partitions(partitions)
        owner = command_fanout_transaction_target_owner(partition, endpoint_id)
        push!(values, effective_command(
            target_local_command_endpoint_state(owner, endpoint_id)))
        push!(physical, getfield(
            target_local_controllable_optic_state(owner), :physical).active)
    end
    return (
        map(effective_command_route_test_host_value, values),
        map(effective_command_route_test_host_value, physical),
    )
end

function command_fanout_transaction_has_publication(
    partitions,
    endpoint_id::CommandEndpointID,
)
    return all(prepared_partitions(partitions)) do partition
        owner = command_fanout_transaction_target_owner(partition, endpoint_id)
        has_effective_command_publication(
            target_local_command_endpoint_state(owner, endpoint_id))
    end
end

function command_fanout_transaction_commands(
    fixture,
    sequence,
    timestamp::PlantTimestamp;
    reversed::Bool=false,
    alpha_payload=[1.0, 2.0],
    beta_payload=[3.0, 4.0],
)
    alpha = PlantCommand(
        fixture.alpha_schema, sequence, timestamp, alpha_payload)
    beta = PlantCommand(
        fixture.beta_schema, sequence, timestamp, beta_payload)
    return reversed ? PlantCommandTransaction(beta, alpha) :
        PlantCommandTransaction(alpha, beta)
end

function command_fanout_transaction_workspace_idle(fixture)
    return all(fixture.workspace.lanes) do lane
        !(lane.selected || lane.has_claim || lane.has_staged ||
          lane.has_finalization || lane.has_transaction_admission)
    end
end

@testset "Authority-owned command fanout" begin
    reset_handoff_test_transfer_controls!()
    (
        partitions,
        schema,
        prepared,
        state,
        workspace,
        endpoint,
        endpoint_state,
        disposition,
    ) = command_fanout_test_fixture()
    endpoint_id = command_endpoint_id(schema)
    timestamp = PlantTimestamp(10)

    admission = command_fanout_test_admit!(
        endpoint,
        endpoint_state,
        disposition,
        schema,
        1,
        timestamp,
        [2.0, 3.0],
        zero(PlantTimestamp),
    )
    @test command_admission_status(admission) == CommandAdmittedPending
    @test command_disposition_count(disposition) == 0

    @test Plant._apply_next_command_fanout!(
        prepared, state, workspace, endpoint_id, timestamp) === nothing
    values, physical = command_fanout_test_replica_values(partitions)
    @test all(==([3.0, 4.0]), values)
    @test all(==([3.0, 4.0]), physical)
    @test state.phase == Plant._CommandFanoutIdle
    @test !command_authority_failed(state.authority)
    @test state.authority.publication_sequences[1] == UInt64(1)
    @test state.authority.publication_timestamps[1] == timestamp
    @test command_disposition_count(disposition) == 1
    applied = command_disposition(disposition, 1)
    @test command_terminal_kind(applied) == AppliedCommand
    @test command_disposition_reason(applied) == CommandDispositionReason(:applied)
    @test HANDOFF_TEST_ACTIVE_DEVICE[] == zero(UInt32)

    clear_command_dispositions!(disposition)
    duplicate = command_fanout_test_admit!(
        endpoint,
        endpoint_state,
        disposition,
        schema,
        1,
        PlantTimestamp(20),
        [9.0, 9.0],
        PlantTimestamp(11),
    )
    @test command_admission_status(duplicate) == CommandTerminatedOnAdmission
    @test command_disposition_count(disposition) == 1
    rejected = command_disposition(disposition, 1)
    @test command_terminal_kind(rejected) == RejectedCommand
    @test command_disposition_reason(rejected) ==
        CommandDispositionReason(:duplicate_sequence)
    values, physical = command_fanout_test_replica_values(partitions)
    @test all(==([3.0, 4.0]), values)
    @test all(==([3.0, 4.0]), physical)
    @test state.phase == Plant._CommandFanoutIdle
    @test state.authority.publication_sequences[1] == UInt64(1)
    lane_state = only(state.lanes)
    @test all(route_state -> route_state.phase ==
        Plant._EffectiveCommandRouteIdle, lane_state.routes)
    @test HANDOFF_TEST_ACTIVE_DEVICE[] == zero(UInt32)
end

@testset "Warmed authority-owned command fanout" begin
    reset_handoff_test_transfer_controls!()
    (
        _,
        schema,
        prepared,
        state,
        workspace,
        endpoint,
        endpoint_state,
        disposition,
    ) = command_fanout_test_fixture()
    endpoint_id = command_endpoint_id(schema)

    command_fanout_test_admit!(
        endpoint,
        endpoint_state,
        disposition,
        schema,
        1,
        PlantTimestamp(1),
        [2.0, 3.0],
        zero(PlantTimestamp),
    )
    Plant._apply_next_command_fanout!(
        prepared, state, workspace, endpoint_id, PlantTimestamp(1))
    clear_command_dispositions!(disposition)

    command_fanout_test_admit!(
        endpoint,
        endpoint_state,
        disposition,
        schema,
        2,
        PlantTimestamp(2),
        [4.0, 5.0],
        PlantTimestamp(1),
    )
    allocations = @allocated Plant._apply_next_command_fanout!(
        prepared, state, workspace, endpoint_id, PlantTimestamp(2))
    clear_command_dispositions!(disposition)
    command_fanout_test_admit!(
        endpoint,
        endpoint_state,
        disposition,
        schema,
        3,
        PlantTimestamp(3),
        [6.0, 7.0],
        PlantTimestamp(2),
    )
    @test @inferred(Plant._apply_next_command_fanout!(
        prepared, state, workspace, endpoint_id, PlantTimestamp(3))) === nothing
    if !coverage_instrumented()
        @test allocations == 0
    end
end

@testset "Authority-owned command-fanout silence" begin
    reset_handoff_test_transfer_controls!()
    (
        safe_partitions,
        _,
        safe_prepared,
        safe_state,
        safe_workspace,
        _,
        _,
        safe_disposition,
    ) = command_fanout_test_silence_fixture(
        ApplySafeCommand;
        safe_command=[-2.0, -3.0],
    )
    endpoint_id = CommandEndpointID(:effective_command_route_endpoint)
    safe_transition = Plant._apply_command_silence_fanout!(
        safe_prepared,
        safe_state,
        safe_workspace,
        endpoint_id,
        PlantTimestamp(10),
    )
    @test command_silence_action(safe_transition) == ApplySafeCommand
    @test command_silence_transition_timestamp(safe_transition) ==
        PlantTimestamp(10)
    values, physical = command_fanout_test_replica_values(safe_partitions)
    @test all(==([-2.0, -3.0]), values)
    @test all(==([-2.0, -3.0]), physical)
    @test !command_fanout_test_no_publication(safe_partitions)
    @test safe_state.authority.publication_sequences[1] == UInt64(1)
    @test safe_state.authority.publication_timestamps[1] ==
        PlantTimestamp(10)
    @test safe_state.phase == Plant._CommandFanoutIdle
    @test !command_authority_failed(safe_state.authority)
    @test command_disposition_count(safe_disposition) == 0
    @test HANDOFF_TEST_ACTIVE_DEVICE[] == zero(UInt32)

    (
        failed_partitions,
        _,
        failed_prepared,
        failed_state,
        failed_workspace,
        _,
        failed_endpoint_state,
        failed_disposition,
    ) = command_fanout_test_silence_fixture(FailOnCommandSilence)
    failed_transition = Plant._apply_command_silence_fanout!(
        failed_prepared,
        failed_state,
        failed_workspace,
        endpoint_id,
        PlantTimestamp(10),
    )
    @test command_silence_action(failed_transition) == FailOnCommandSilence
    @test command_endpoint_failed(failed_endpoint_state)
    @test command_fanout_test_no_publication(failed_partitions)
    values, physical = command_fanout_test_replica_values(failed_partitions)
    @test all(==([1.0, 1.0]), values)
    @test all(==([1.0, 1.0]), physical)
    @test failed_state.authority.publication_sequences[1] == UInt64(0)
    @test failed_state.phase == Plant._CommandFanoutFailed
    @test command_authority_failed(failed_state.authority)
    @test command_disposition_count(failed_disposition) == 0
    @test HANDOFF_TEST_ACTIVE_DEVICE[] == zero(UInt32)

    (
        hold_partitions,
        _,
        hold_prepared,
        hold_state,
        hold_workspace,
        _,
        _,
        hold_disposition,
    ) = command_fanout_test_silence_fixture(HoldLastCommand)
    hold_error = command_fanout_test_error(() ->
        Plant._apply_command_silence_fanout!(
            hold_prepared,
            hold_state,
            hold_workspace,
            endpoint_id,
            PlantTimestamp(10),
        ))
    @test hold_error isa PlantCommandError
    if hold_error isa PlantCommandError
        @test hold_error.stage === :silence
        @test hold_error.reason === :silence_not_scheduled
    end
    @test command_fanout_test_no_publication(hold_partitions)
    values, physical = command_fanout_test_replica_values(hold_partitions)
    @test all(==([1.0, 1.0]), values)
    @test all(==([1.0, 1.0]), physical)
    @test hold_state.authority.publication_sequences[1] == UInt64(0)
    @test hold_state.phase == Plant._CommandFanoutIdle
    @test !command_authority_failed(hold_state.authority)
    @test command_disposition_count(hold_disposition) == 0
    @test HANDOFF_TEST_ACTIVE_DEVICE[] == zero(UInt32)
end

@testset "Warmed authority-owned command-fanout silence" begin
    reset_handoff_test_transfer_controls!()
    (
        _,
        schema,
        prepared,
        state,
        workspace,
        endpoint,
        endpoint_state,
        disposition,
    ) = command_fanout_test_silence_fixture(
        ApplySafeCommand;
        safe_command=[-2.0, -3.0],
    )
    endpoint_id = command_endpoint_id(schema)

    Plant._apply_command_silence_fanout!(
        prepared, state, workspace, endpoint_id, PlantTimestamp(10))
    command_fanout_test_admit!(
        endpoint,
        endpoint_state,
        disposition,
        schema,
        1,
        PlantTimestamp(11),
        [2.0, 3.0],
        PlantTimestamp(10),
    )
    Plant._apply_next_command_fanout!(
        prepared, state, workspace, endpoint_id, PlantTimestamp(11))
    clear_command_dispositions!(disposition)
    allocations = @allocated Plant._apply_command_silence_fanout!(
        prepared, state, workspace, endpoint_id, PlantTimestamp(21))

    command_fanout_test_admit!(
        endpoint,
        endpoint_state,
        disposition,
        schema,
        2,
        PlantTimestamp(22),
        [2.0, 3.0],
        PlantTimestamp(21),
    )
    Plant._apply_next_command_fanout!(
        prepared, state, workspace, endpoint_id, PlantTimestamp(22))
    clear_command_dispositions!(disposition)
    @test @inferred(Plant._apply_command_silence_fanout!(
        prepared, state, workspace, endpoint_id, PlantTimestamp(32))) isa
        PlantCommandSilenceTransition
    if !coverage_instrumented()
        @test allocations == 0
    end
end

@testset "Authority-owned atomic command fanout" begin
    reset_handoff_test_transfer_controls!()
    fixture = command_fanout_transaction_fixture()
    timestamp = PlantTimestamp(10)
    transaction = command_fanout_transaction_commands(
        fixture,
        1,
        timestamp;
        reversed=true,
    )
    admission = Plant._admit_command_fanout_transaction!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        transaction,
        zero(PlantTimestamp),
    )
    @test command_admission_status(admission) == CommandAdmittedPending
    @test command_transaction_id(admission) == UInt64(1)
    @test command_transaction_member_count(admission) == UInt32(2)
    @test command_scheduled_timestamp(admission) == timestamp
    @test fixture.state.phase == Plant._CommandFanoutIdle
    @test pending_command_count(fixture.alpha_state) == 1
    @test pending_command_count(fixture.beta_state) == 1

    @test Plant._apply_next_command_fanout!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        fixture.beta_id,
        timestamp,
    ) === nothing
    alpha_values, alpha_physical = command_fanout_transaction_replica_values(
        fixture.partitions, fixture.alpha_id)
    beta_values, beta_physical = command_fanout_transaction_replica_values(
        fixture.partitions, fixture.beta_id)
    @test all(==([2.0, 3.0]), alpha_values)
    @test all(==([2.0, 3.0]), alpha_physical)
    @test all(==([4.0, 5.0]), beta_values)
    @test all(==([4.0, 5.0]), beta_physical)
    @test command_fanout_transaction_has_publication(
        fixture.partitions, fixture.alpha_id)
    @test command_fanout_transaction_has_publication(
        fixture.partitions, fixture.beta_id)
    @test all(==(UInt64(1)), fixture.state.authority.publication_sequences)
    @test all(
        ==(timestamp), fixture.state.authority.publication_timestamps)
    @test fixture.state.phase == Plant._CommandFanoutIdle
    @test !command_authority_failed(fixture.state.authority)
    @test command_disposition_count(fixture.alpha_disposition) == 1
    @test command_disposition_count(fixture.beta_disposition) == 1
    @test command_terminal_kind(
        command_disposition(fixture.alpha_disposition, 1)) == AppliedCommand
    @test command_terminal_kind(
        command_disposition(fixture.beta_disposition, 1)) == AppliedCommand
    @test HANDOFF_TEST_ACTIVE_DEVICE[] == zero(UInt32)

    rejected_fixture = command_fanout_transaction_fixture()
    rejected_transaction = command_fanout_transaction_commands(
        rejected_fixture,
        1,
        timestamp;
        reversed=true,
        beta_payload=[9.0],
    )
    rejected_admission = Plant._admit_command_fanout_transaction!(
        rejected_fixture.prepared,
        rejected_fixture.state,
        rejected_fixture.workspace,
        rejected_transaction,
        zero(PlantTimestamp),
    )
    @test command_admission_status(rejected_admission) ==
        CommandTerminatedOnAdmission
    @test command_transaction_id(rejected_admission) === nothing
    @test command_scheduled_timestamp(rejected_admission) === nothing
    @test rejected_fixture.state.authority.command_transaction_sequence ==
        UInt64(0)
    @test all(endpoint_state ->
        iszero(pending_command_count(endpoint_state)) &&
            iszero(active_command_count(endpoint_state)),
        (rejected_fixture.alpha_state, rejected_fixture.beta_state))
    @test !rejected_fixture.alpha_state.has_admission
    @test !rejected_fixture.beta_state.has_admission
    @test !command_fanout_transaction_has_publication(
        rejected_fixture.partitions, rejected_fixture.alpha_id)
    @test !command_fanout_transaction_has_publication(
        rejected_fixture.partitions, rejected_fixture.beta_id)
    @test all(iszero, rejected_fixture.state.authority.publication_sequences)
    @test rejected_fixture.state.phase == Plant._CommandFanoutIdle
    @test !command_authority_failed(rejected_fixture.state.authority)
    @test command_disposition_count(rejected_fixture.alpha_disposition) == 1
    @test command_disposition_count(rejected_fixture.beta_disposition) == 1
    @test all(disposition -> command_terminal_kind(disposition) ==
        RejectedCommand, (
            command_disposition(rejected_fixture.alpha_disposition, 1),
            command_disposition(rejected_fixture.beta_disposition, 1),
        ))
    rejection_reasons = map(command_disposition_reason, (
        command_disposition(rejected_fixture.alpha_disposition, 1),
        command_disposition(rejected_fixture.beta_disposition, 1),
    ))
    @test all(reason -> startswith(String(reason.name),
        "atomic_transaction_aborted_"), rejection_reasons)

    structural_fixture = command_fanout_transaction_fixture()
    first_lane, second_lane = structural_fixture.prepared.lanes
    colliding_second_lane = Plant._PreparedCommandFanoutLane(
        second_lane.authority_endpoint,
        second_lane.endpoint_slot,
        first_lane.optic_slot,
        second_lane.routes,
    )
    colliding_prepared = Plant._PreparedCommandFanout(
        structural_fixture.prepared.binding,
        structural_fixture.prepared.authority,
        (first_lane, colliding_second_lane),
    )
    structural_error = command_fanout_test_error(() ->
        Plant._admit_command_fanout_transaction!(
            colliding_prepared,
            structural_fixture.state,
            structural_fixture.workspace,
            command_fanout_transaction_commands(
                structural_fixture, 1, timestamp),
            zero(PlantTimestamp),
        ))
    @test structural_error isa PlantCommandError
    if structural_error isa PlantCommandError
        @test structural_error.stage === :transaction
        @test structural_error.reason === :duplicate_physical_optic
    end
    @test structural_fixture.state.authority.command_transaction_sequence ==
        UInt64(0)
    @test all(endpoint_state ->
        iszero(pending_command_count(endpoint_state)) &&
            iszero(active_command_count(endpoint_state)),
        (structural_fixture.alpha_state, structural_fixture.beta_state))
    @test !command_fanout_transaction_has_publication(
        structural_fixture.partitions, structural_fixture.alpha_id)
    @test !command_fanout_transaction_has_publication(
        structural_fixture.partitions, structural_fixture.beta_id)
    @test structural_fixture.state.phase == Plant._CommandFanoutIdle
    @test HANDOFF_TEST_ACTIVE_DEVICE[] == zero(UInt32)
end

@testset "Canonical atomic command-fanout transaction precedence" begin
    reset_handoff_test_transfer_controls!()
    timestamp = PlantTimestamp(10)
    function rejected_transaction_reasons(fixture, transaction)
        admission = Plant._admit_command_fanout_transaction!(
            fixture.prepared,
            fixture.state,
            fixture.workspace,
            transaction,
            zero(PlantTimestamp),
        )
        @test command_admission_status(admission) ==
            CommandTerminatedOnAdmission
        @test command_disposition_count(fixture.alpha_disposition) == 1
        @test command_disposition_count(fixture.beta_disposition) == 1
        return map(command_disposition_reason, (
            command_disposition(fixture.alpha_disposition, 1),
            command_disposition(fixture.beta_disposition, 1),
        ))
    end

    forward_fixture = command_fanout_transaction_fixture()
    forward_reasons = rejected_transaction_reasons(
        forward_fixture,
        command_fanout_transaction_commands(
            forward_fixture,
            1,
            timestamp;
            alpha_payload=[1.0],
            beta_payload=Float32[3.0, 4.0],
        ),
    )
    reversed_fixture = command_fanout_transaction_fixture()
    reversed_reasons = rejected_transaction_reasons(
        reversed_fixture,
        command_fanout_transaction_commands(
            reversed_fixture,
            1,
            timestamp;
            reversed=true,
            alpha_payload=[1.0],
            beta_payload=Float32[3.0, 4.0],
        ),
    )
    @test forward_reasons == reversed_reasons
    @test all(reason -> startswith(String(reason.name),
        "atomic_transaction_aborted_"), forward_reasons)
    @test all(iszero, forward_fixture.state.authority.publication_sequences)
    @test all(iszero, reversed_fixture.state.authority.publication_sequences)
    @test command_fanout_transaction_workspace_idle(forward_fixture)
    @test command_fanout_transaction_workspace_idle(reversed_fixture)
    @test HANDOFF_TEST_ACTIVE_DEVICE[] == zero(UInt32)
end

@testset "Atomic command-fanout structural failure drains full calendars" begin
    reset_handoff_test_transfer_controls!()
    fail_late = CommandEffectiveTimePolicy(
        AllowFutureCommand,
        FailOnLateCommand,
        PreservePendingCommands,
    )
    fixture = command_fanout_transaction_fixture(
        alpha_effective_time_policy=fail_late)
    for (sequence, effective_timestamp) in
            ((1, PlantTimestamp(20)), (2, PlantTimestamp(30)))
        admission = command_fanout_test_admit!(
            fixture.beta_endpoint,
            fixture.beta_state,
            fixture.beta_disposition,
            fixture.beta_schema,
            sequence,
            effective_timestamp,
            [3.0, 4.0],
            zero(PlantTimestamp),
        )
        @test command_admission_status(admission) == CommandAdmittedPending
    end
    @test pending_command_count(fixture.beta_state) == 2

    late_transaction = PlantCommandTransaction(
        PlantCommand(fixture.beta_schema, 3, PlantTimestamp(0), [3.0, 4.0]),
        PlantCommand(fixture.alpha_schema, 1, PlantTimestamp(0), [1.0, 2.0]),
    )
    failure = command_fanout_test_error(() ->
        Plant._admit_command_fanout_transaction!(
            fixture.prepared,
            fixture.state,
            fixture.workspace,
            late_transaction,
            PlantTimestamp(10),
        ))
    @test failure isa PlantCommandError
    if failure isa PlantCommandError
        @test failure.stage === :transaction
        @test failure.reason === :late_command
    end
    @test command_disposition_count(fixture.alpha_disposition) == 1
    @test command_disposition_count(fixture.beta_disposition) == 3
    @test all(index -> command_terminal_kind(
            command_disposition(fixture.alpha_disposition, index)) ==
            FailedCommand, 1:command_disposition_count(fixture.alpha_disposition))
    @test all(index -> command_terminal_kind(
            command_disposition(fixture.beta_disposition, index)) ==
            FailedCommand, 1:command_disposition_count(fixture.beta_disposition))
    @test all(index -> command_disposition_reason(
            command_disposition(fixture.alpha_disposition, index)).name ==
            :atomic_transaction_aborted_late_command,
        1:command_disposition_count(fixture.alpha_disposition))
    @test all(index -> command_disposition_reason(
            command_disposition(fixture.beta_disposition, index)).name ==
            :atomic_transaction_aborted_late_command,
        1:command_disposition_count(fixture.beta_disposition))
    @test length(unique(command_presentation_id(
        command_disposition(fixture.beta_disposition, index))
        for index in 1:command_disposition_count(fixture.beta_disposition))) == 3
    @test all(endpoint_state ->
        iszero(pending_command_count(endpoint_state)) &&
            iszero(active_command_count(endpoint_state)) && endpoint_state.failed,
        (fixture.alpha_state, fixture.beta_state))
    @test command_authority_failed(fixture.state.authority)
    @test fixture.state.phase == Plant._CommandFanoutFailed
    @test command_fanout_transaction_workspace_idle(fixture)
    @test HANDOFF_TEST_ACTIVE_DEVICE[] == zero(UInt32)
end

@testset "Atomic command-fanout commit recovery" begin
    reset_handoff_test_transfer_controls!()
    fixture = command_fanout_transaction_fixture()
    timestamp = PlantTimestamp(10)
    transaction = command_fanout_transaction_commands(
        fixture, 1, timestamp; reversed=true)
    Plant._select_canonical_command_fanout_transaction_lanes!(
        fixture.prepared.lanes,
        fixture.workspace.lanes,
        transaction.commands,
    )
    Plant._with_completed_prepared_device_execution_context(
        fixture.prepared.authority.context) do
        Plant._preflight_canonical_command_fanout_transaction!(
            fixture.prepared.lanes,
            fixture.state.lanes,
            fixture.workspace.lanes,
            transaction.commands,
            zero(PlantTimestamp),
        )
    end
    fixture.state.phase = Plant._CommandFanoutCommitting
    Plant._commit_canonical_command_fanout_lane!(
        transaction.commands,
        first(fixture.prepared.lanes),
        first(fixture.state.lanes),
        first(fixture.workspace.lanes),
        zero(PlantTimestamp),
        UInt64(1),
        UInt32(2),
    )
    @test pending_command_count(fixture.alpha_state) == 1
    @test pending_command_count(fixture.beta_state) == 0
    Plant._fail_command_fanout_transaction_commit!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        transaction,
        zero(PlantTimestamp),
        UInt64(1),
    )
    @test command_disposition_count(fixture.alpha_disposition) == 1
    @test command_disposition_count(fixture.beta_disposition) == 1
    @test all(disposition ->
        command_terminal_kind(disposition) == FailedCommand &&
            command_disposition_reason(disposition).name ==
            :transaction_commit_failure, (
            command_disposition(fixture.alpha_disposition, 1),
            command_disposition(fixture.beta_disposition, 1),
        ))
    @test all(endpoint_state ->
        iszero(pending_command_count(endpoint_state)) &&
            iszero(active_command_count(endpoint_state)) && endpoint_state.failed,
        (fixture.alpha_state, fixture.beta_state))
    @test fixture.state.authority.command_transaction_sequence == UInt64(0)
    @test fixture.alpha_state.presentation_sequence == UInt64(1)
    @test fixture.beta_state.presentation_sequence == UInt64(1)
    @test command_authority_failed(fixture.state.authority)
    @test fixture.state.phase == Plant._CommandFanoutUncertain
    @test command_fanout_transaction_workspace_idle(fixture)
    @test HANDOFF_TEST_ACTIVE_DEVICE[] == zero(UInt32)
end

@testset "Atomic command-fanout transaction overflow precedence" begin
    reset_handoff_test_transfer_controls!()
    fixture = command_fanout_transaction_fixture()
    fixture.state.authority.command_transaction_sequence = typemax(UInt64)
    invalid_transaction = command_fanout_transaction_commands(
        fixture,
        1,
        PlantTimestamp(10);
        reversed=true,
        beta_payload=[9.0],
    )
    failure = command_fanout_test_error(() ->
        Plant._admit_command_fanout_transaction!(
            fixture.prepared,
            fixture.state,
            fixture.workspace,
            invalid_transaction,
            zero(PlantTimestamp),
        ))
    @test failure isa PlantCommandError
    if failure isa PlantCommandError
        @test failure.stage === :transaction
        @test failure.reason === :transaction_overflow
    end
    @test fixture.state.authority.command_transaction_sequence ==
        typemax(UInt64)
    @test command_disposition_count(fixture.alpha_disposition) == 0
    @test command_disposition_count(fixture.beta_disposition) == 0
    @test all(endpoint_state ->
        iszero(pending_command_count(endpoint_state)) &&
            iszero(active_command_count(endpoint_state)),
        (fixture.alpha_state, fixture.beta_state))
    @test fixture.state.phase == Plant._CommandFanoutIdle
    @test !command_authority_failed(fixture.state.authority)
    @test command_fanout_transaction_workspace_idle(fixture)
    @test HANDOFF_TEST_ACTIVE_DEVICE[] == zero(UInt32)
end

@testset "Warmed authority-owned atomic command fanout" begin
    reset_handoff_test_transfer_controls!()
    fixture = command_fanout_transaction_fixture()
    first_transaction = command_fanout_transaction_commands(
        fixture, 1, PlantTimestamp(1); reversed=true)
    Plant._admit_command_fanout_transaction!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        first_transaction,
        zero(PlantTimestamp),
    )
    Plant._apply_next_command_fanout!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        fixture.alpha_id,
        PlantTimestamp(1),
    )
    clear_command_dispositions!(fixture.alpha_disposition)
    clear_command_dispositions!(fixture.beta_disposition)

    second_transaction = command_fanout_transaction_commands(
        fixture, 2, PlantTimestamp(2); reversed=true)
    admission_allocations = @allocated Plant._admit_command_fanout_transaction!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        second_transaction,
        PlantTimestamp(1),
    )
    application_allocations = @allocated Plant._apply_next_command_fanout!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        fixture.beta_id,
        PlantTimestamp(2),
    )
    clear_command_dispositions!(fixture.alpha_disposition)
    clear_command_dispositions!(fixture.beta_disposition)

    third_transaction = command_fanout_transaction_commands(
        fixture, 3, PlantTimestamp(3); reversed=true)
    @test @inferred(Plant._admit_command_fanout_transaction!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        third_transaction,
        PlantTimestamp(2),
    )) isa PlantCommandTransactionAdmission
    @test @inferred(Plant._apply_next_command_fanout!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        fixture.alpha_id,
        PlantTimestamp(3),
    )) === nothing
    if !coverage_instrumented()
        @test admission_allocations == 0
        @test application_allocations == 0
    end
end

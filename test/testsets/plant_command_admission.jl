function _admission_schema(::Type{T}=Float32;
    id=:dm_schema,
    version=1,
    endpoint=:dm_command,
    dimensions=(3,),
    semantics=AbsoluteCommand,
    basis=CommandBasis(:actuator, :test_actuators),
    bounds=UniformCommandBounds(T(-1), T(1)),
    value_policy=CommandValuePolicy(),
    sequence_policy=CommandSequencePolicy(),
    effective_time_policy=CommandEffectiveTimePolicy(),
    silence_policy=CommandSilencePolicy(),
) where {T}
    return PlantCommandSchema(
        T,
        dimensions;
        id,
        version,
        endpoint,
        units=:metre,
        sign_convention=:positive_surface_increases_opd,
        basis,
        basis_revision=1,
        semantics,
        bounds,
        value_policy,
        sequence_policy,
        effective_time_policy,
        silence_policy,
    )
end

function _captured_command_admission_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

function _assert_command_admission_error(f, stage::Symbol, reason::Symbol)
    error = _captured_command_admission_error(f)
    @test error isa PlantCommandError
    if error isa PlantCommandError
        @test error.stage === stage
        @test error.reason === reason
        @test !isempty(error.msg)
    end
    return error
end

function _assert_single_disposition(workspace, kind, reason)
    @test command_disposition_count(workspace) == 1
    disposition = command_disposition(workspace, 1)
    @test command_terminal_kind(disposition) == kind
    @test command_disposition_reason(disposition) ==
        CommandDispositionReason(reason)
    return disposition
end

function _scalar_command_cycle!(workspace, endpoint, state, sequence,
    timestamp)
    command = PlantCommand(command_schema(endpoint), sequence, timestamp, 0.25)
    admission = admit_plant_command!(workspace, endpoint, state, command,
        timestamp)
    command_admission_status(admission) == CommandAdmittedReady ||
        return nothing
    claim = claim_next_application_ready_command!(endpoint, state, timestamp)
    claim === nothing && return nothing
    claimed_command_payload(endpoint, state, claim)
    mark_plant_command_applied!(workspace, endpoint, state, claim)
    clear_command_dispositions!(workspace)
    return nothing
end

function _array_command_cycle!(workspace, endpoint, state, payload, sequence,
    timestamp)
    command = PlantCommand(command_schema(endpoint), sequence, timestamp,
        payload)
    admission = admit_plant_command!(workspace, endpoint, state, command,
        timestamp)
    command_admission_status(admission) == CommandAdmittedReady ||
        return nothing
    claim = claim_next_application_ready_command!(endpoint, state, timestamp)
    claim === nothing && return nothing
    claimed_command_payload(endpoint, state, claim)
    mark_plant_command_applied!(workspace, endpoint, state, claim)
    clear_command_dispositions!(workspace)
    return nothing
end

function _multi_entry_command_cycle!(workspace, endpoint, state,
    first_sequence, base_timestamp)
    schema = command_schema(endpoint)
    timestamp_10 = base_timestamp + PlantDuration(10)
    timestamp_20 = base_timestamp + PlantDuration(20)
    timestamp_30 = base_timestamp + PlantDuration(30)
    admit_plant_command!(workspace, endpoint, state,
        PlantCommand(schema, first_sequence, timestamp_30, 0.3),
        base_timestamp)
    admit_plant_command!(workspace, endpoint, state,
        PlantCommand(schema, first_sequence + 1, timestamp_10, 0.1),
        base_timestamp)
    admit_plant_command!(workspace, endpoint, state,
        PlantCommand(schema, first_sequence + 2, timestamp_20, 0.2),
        base_timestamp)
    for timestamp in (timestamp_10, timestamp_20, timestamp_30)
        claim = claim_next_application_ready_command!(endpoint, state,
            timestamp)
        claimed_command_payload(endpoint, state, claim)
        mark_plant_command_applied!(workspace, endpoint, state, claim)
        clear_command_dispositions!(workspace)
    end
    return nothing
end

function _supersession_command_cycle!(workspace, endpoint, state,
    first_sequence, base_timestamp)
    schema = command_schema(endpoint)
    timestamp_20 = base_timestamp + PlantDuration(20)
    timestamp_25 = base_timestamp + PlantDuration(25)
    timestamp_30 = base_timestamp + PlantDuration(30)
    admit_plant_command!(workspace, endpoint, state,
        PlantCommand(schema, first_sequence + 1, timestamp_30, 0.3),
        base_timestamp)
    admit_plant_command!(workspace, endpoint, state,
        PlantCommand(schema, first_sequence, timestamp_20, 0.2),
        base_timestamp)
    admit_plant_command!(workspace, endpoint, state,
        PlantCommand(schema, first_sequence + 2, timestamp_25, 0.25),
        base_timestamp)
    clear_command_dispositions!(workspace)
    claim = claim_next_application_ready_command!(endpoint, state,
        timestamp_25)
    claimed_command_payload(endpoint, state, claim)
    mark_plant_command_applied!(workspace, endpoint, state, claim)
    clear_command_dispositions!(workspace)
    return nothing
end

function _capacity_rejection_command_cycle!(workspace, endpoint, state,
    first_sequence, base_timestamp)
    schema = command_schema(endpoint)
    timestamp_10 = base_timestamp + PlantDuration(10)
    timestamp_20 = base_timestamp + PlantDuration(20)
    admit_plant_command!(workspace, endpoint, state,
        PlantCommand(schema, first_sequence, timestamp_10, 0.1),
        base_timestamp)
    admit_plant_command!(workspace, endpoint, state,
        PlantCommand(schema, first_sequence + 1, timestamp_20, 0.2),
        base_timestamp)
    clear_command_dispositions!(workspace)
    claim = claim_next_application_ready_command!(endpoint, state,
        timestamp_10)
    claimed_command_payload(endpoint, state, claim)
    mark_plant_command_applied!(workspace, endpoint, state, claim)
    clear_command_dispositions!(workspace)
    return nothing
end

function _failure_drain_command_cycle!(workspace, endpoint, state,
    first_sequence, base_timestamp)
    schema = command_schema(endpoint)
    for (sequence_offset, time_offset) in ((0, 30), (1, 10), (2, 20))
        admit_plant_command!(workspace, endpoint, state,
            PlantCommand(schema, first_sequence + sequence_offset,
                base_timestamp + PlantDuration(time_offset), 0.0),
            base_timestamp)
    end
    fail_pending_plant_commands!(workspace, endpoint, state,
        base_timestamp + PlantDuration(5))
    clear_command_dispositions!(workspace)
    return nothing
end

function _rebuilt_application_claim(claim;
    key=command_order_key(claim),
    requested=command_requested_effective_timestamp(claim),
    admission=command_admission_timestamp(claim),
    ready=command_ready_timestamp(claim),
)
    return PlantCommandApplicationClaim(
        getfield(claim, :binding_id),
        getfield(claim, :slot),
        getfield(claim, :generation),
        command_presentation_id(claim),
        key,
        requested,
        admission,
        ready,
        AdaptiveOpticsSim._PLANT_COMMAND_CLAIM_TOKEN,
    )
end

struct AdmissionValidationFailureVector{T} <: AbstractVector{T} end

Base.size(::AdmissionValidationFailureVector) = (3,)
Base.IndexStyle(::Type{<:AdmissionValidationFailureVector}) = IndexLinear()
Base.getindex(::AdmissionValidationFailureVector, ::Int) =
    error("declared validation failure")

struct AdmissionStagingFailureVector{T} <: AbstractVector{T}
    values::Vector{T}
end

Base.size(values::AdmissionStagingFailureVector) = size(values.values)
Base.IndexStyle(::Type{<:AdmissionStagingFailureVector}) = IndexLinear()
Base.getindex(values::AdmissionStagingFailureVector, index::Int) =
    values.values[index]
Base.copyto!(::Vector{T}, ::AdmissionStagingFailureVector{T}) where {T} =
    error("declared staging failure")

@testset "Bounded plant-command admission API" begin
    for name in (
        :PlantCommandSequence,
        :CommandPresentationID,
        :CommandDispositionReason,
        :PlantCommand,
        :CommandSequenceClass,
        :CommandAdmissionStatus,
        :CommandTerminalKind,
        :PlantCommandOrderKey,
        :PlantCommandAdmission,
        :PlantCommandDisposition,
        :PreparedCommandEndpoint,
        :PlantCommandApplicationClaim,
        :prepare_command_endpoint,
        :validate_plant_command,
        :admit_plant_command!,
        :claim_next_application_ready_command!,
        :claimed_command_payload,
        :mark_plant_command_applied!,
        :fail_plant_command_application!,
        :fail_pending_plant_commands!,
        :command_requested_effective_timestamp,
        :command_scheduled_timestamp,
        :command_ready_timestamp,
        :command_terminal_timestamp,
    )
        @test Base.isexported(AdaptiveOpticsSim, name)
    end
    @test !Base.isexported(AdaptiveOpticsSim, :CommandEndpointState)
    @test !Base.isexported(AdaptiveOpticsSim, :CommandDispositionWorkspace)

    schema = _admission_schema()
    endpoint = @inferred prepare_command_endpoint(
        schema;
        capacity=3,
        sequence_window=5,
        ordinal=7,
    )
    state = @inferred CommandEndpointState(
        endpoint;
        initial_timestamp=PlantTimestamp(10),
    )
    workspace = @inferred CommandDispositionWorkspace(endpoint)

    @test endpoint isa PreparedCommandEndpoint
    @test !Base.ismutable(endpoint)
    @test Base.ismutable(state)
    @test Base.ismutable(workspace)
    @test command_schema(endpoint) === schema
    @test command_endpoint_id(endpoint) == CommandEndpointID(:dm_command)
    @test command_endpoint_capacity(endpoint) == 3
    @test command_sequence_window(endpoint) == 5
    @test command_endpoint_ordinal(endpoint) == UInt32(7)
    @test backend(endpoint) == CPUBackend()
    @test pending_command_count(state) == 0
    @test active_command_count(state) == 0
    @test command_endpoint_timestamp(state) == PlantTimestamp(10)
    @test command_disposition_count(workspace) == 0
    @test length(getfield(state, :slots)) == 3
    @test length(getfield(state, :calendar)) == 3
    @test length(getfield(state, :accepted_sequences)) == 5
    @test length(getfield(getfield(state, :payloads), :values)) == 4

    scalar_schema = _admission_schema(
        Float64;
        id=:focus_schema,
        endpoint=:focus_command,
        dimensions=(),
        basis=CommandBasis(:rigid_body, :focus),
    )
    scalar_endpoint = prepare_command_endpoint(
        scalar_schema;
        capacity=2,
        ordinal=8,
    )
    scalar_state = CommandEndpointState(scalar_endpoint)
    @test length(getfield(getfield(scalar_state, :payloads), :values)) == 2

    sequence = PlantCommandSequence(4)
    presentation = CommandPresentationID(9)
    reason = CommandDispositionReason(:late_command)
    @test sequence == PlantCommandSequence(4)
    @test presentation == CommandPresentationID(9)
    @test reason == CommandDispositionReason(:late_command)
    @test length(Set((sequence, PlantCommandSequence(4)))) == 1
    @test sprint(show, sequence) == "PlantCommandSequence(4)"
    @test isbitstype(PlantCommandSequence)
    @test isbitstype(CommandPresentationID)
    @test isbitstype(PlantCommandOrderKey)
    @test isbitstype(PlantCommandApplicationClaim)

    payload = Float32[0, 0.25, 0.5]
    command = @inferred PlantCommand(schema, sequence, PlantTimestamp(20),
        payload)
    @test command_endpoint_id(command) == CommandEndpointID(:dm_command)
    @test command_schema_id(command) == PlantCommandSchemaID(:dm_schema)
    @test command_schema_version(command) == PlantCommandSchemaVersion(1)
    @test command_sequence(command) == sequence
    @test command_requested_effective_timestamp(command) == PlantTimestamp(20)
    @test command_payload(command) === payload
    @test validate_plant_command(endpoint, command) === command

    for (operation, stage, failure_reason) in (
        (() -> PlantCommandSequence(0), :command, :invalid_sequence),
        (() -> PlantCommandSequence(true), :command, :invalid_sequence),
        (() -> CommandPresentationID(-1), :command,
            :invalid_presentation),
        (() -> CommandDispositionReason(Symbol("")), :disposition,
            :empty_reason),
        (() -> PlantCommand(schema, 1, 1, payload), :command,
            :invalid_effective_timestamp),
        (() -> prepare_command_endpoint(schema; capacity=0, ordinal=1),
            :preparation, :invalid_capacity),
        (() -> prepare_command_endpoint(schema; capacity=true, ordinal=1),
            :preparation, :invalid_capacity),
        (() -> prepare_command_endpoint(schema; capacity=1,
            sequence_window=0, ordinal=1), :preparation,
            :invalid_sequence_window),
        (() -> prepare_command_endpoint(schema; capacity=1, ordinal=false),
            :preparation, :invalid_ordinal),
        (() -> prepare_command_endpoint(scalar_schema; capacity=1, ordinal=1,
            backend=CUDABackend()), :preparation,
            :unsupported_scalar_backend),
        (() -> command_disposition(workspace, 1), :disposition,
            :invalid_index),
        (() -> command_disposition(workspace, 1.0), :disposition,
            :invalid_index),
    )
        _assert_command_admission_error(operation, stage, failure_reason)
    end

    incremental_supersession = _admission_schema(
        Float32;
        id=:incremental_schema,
        endpoint=:incremental_command,
        semantics=IncrementalCommand,
        effective_time_policy=CommandEffectiveTimePolicy(
            supersession=SupersedeOlderPendingCommands,
        ),
    )
    _assert_command_admission_error(
        () -> prepare_command_endpoint(
            incremental_supersession;
            capacity=2,
            ordinal=1,
        ),
        :preparation,
        :invalid_incremental_supersession,
    )
end

@testset "Validation outcomes and copied payload ownership" begin
    schema = _admission_schema()
    endpoint = prepare_command_endpoint(schema; capacity=2, ordinal=1)
    state = CommandEndpointState(endpoint)
    workspace = CommandDispositionWorkspace(endpoint)

    invalid_commands = (
        (PlantCommand(:other_command, :dm_schema, 1, 1,
            PlantTimestamp(1), zeros(Float32, 3)), :endpoint_mismatch),
        (PlantCommand(:dm_command, :other_schema, 1, 1,
            PlantTimestamp(2), zeros(Float32, 3)), :schema_mismatch),
        (PlantCommand(:dm_command, :dm_schema, 2, 1,
            PlantTimestamp(3), zeros(Float32, 3)),
            :schema_version_mismatch),
        (PlantCommand(schema, 1, PlantTimestamp(4), zeros(Float64, 3)),
            :element_type),
        (PlantCommand(schema, 1, PlantTimestamp(5),
            Float32[0, NaN, 0]), :nonfinite_rejected),
    )
    for (index, (command, expected_reason)) in enumerate(invalid_commands)
        admission = admit_plant_command!(workspace, endpoint, state, command,
            PlantTimestamp(index))
        @test command_admission_status(admission) ==
            CommandTerminatedOnAdmission
        @test command_sequence_class(admission) === nothing
        @test command_order_key(admission) === nothing
        disposition = _assert_single_disposition(
            workspace,
            RejectedCommand,
            expected_reason,
        )
        @test command_presentation_id(disposition) ==
            command_presentation_id(admission)
        @test command_endpoint_id(disposition) == command_endpoint_id(endpoint)
        @test command_sequence(disposition) == command_sequence(command)
        clear_command_dispositions!(workspace)
    end
    @test pending_command_count(state) == 0
    @test active_command_count(state) == 0

    validation_failure_endpoint = prepare_command_endpoint(
        schema;
        capacity=1,
        ordinal=9,
    )
    validation_failure_state = CommandEndpointState(
        validation_failure_endpoint)
    validation_failure_workspace = CommandDispositionWorkspace(
        validation_failure_endpoint)
    validation_failure_command = PlantCommand(
        schema,
        1,
        PlantTimestamp(1),
        AdmissionValidationFailureVector{Float32}(),
    )
    _assert_command_admission_error(
        () -> admit_plant_command!(validation_failure_workspace,
            validation_failure_endpoint, validation_failure_state,
            validation_failure_command, PlantTimestamp(1)),
        :admission,
        :payload_validation_failure,
    )
    _assert_single_disposition(validation_failure_workspace, FailedCommand,
        :payload_validation_failure)

    staging_failure_endpoint = prepare_command_endpoint(
        schema;
        capacity=2,
        ordinal=10,
    )
    staging_failure_state = CommandEndpointState(staging_failure_endpoint)
    staging_failure_workspace = CommandDispositionWorkspace(
        staging_failure_endpoint)
    retained_payload = Float32[0.1, 0.2, 0.3]
    retained = admit_plant_command!(staging_failure_workspace,
        staging_failure_endpoint, staging_failure_state,
        PlantCommand(schema, 1, PlantTimestamp(10), retained_payload),
        PlantTimestamp(0))
    retained_key = command_order_key(retained)
    staging_failure_command = PlantCommand(
        schema,
        2,
        PlantTimestamp(20),
        AdmissionStagingFailureVector(Float32[0.4, 0.5, 0.6]),
    )
    _assert_command_admission_error(
        () -> admit_plant_command!(staging_failure_workspace,
            staging_failure_endpoint, staging_failure_state,
            staging_failure_command, PlantTimestamp(1)),
        :admission,
        :payload_storage_failure,
    )
    _assert_single_disposition(staging_failure_workspace, FailedCommand,
        :payload_storage_failure)
    @test pending_command_count(staging_failure_state) == 1
    @test active_command_count(staging_failure_state) == 1
    @test next_command_order_key(staging_failure_endpoint,
        staging_failure_state) == retained_key
    clear_command_dispositions!(staging_failure_workspace)
    retried = admit_plant_command!(staging_failure_workspace,
        staging_failure_endpoint, staging_failure_state,
        PlantCommand(schema, 2, PlantTimestamp(20),
            Float32[0.4, 0.5, 0.6]),
        PlantTimestamp(1))
    @test command_sequence_class(retried) == InOrderCommandSequence
    @test command_presentation_id(retried) == CommandPresentationID(3)
    @test pending_command_count(staging_failure_state) == 2

    fail_schema = _admission_schema(
        Float32;
        id=:fail_schema,
        endpoint=:fail_command,
        value_policy=CommandValuePolicy(
            FailOnInvalidCommand,
            FailOnInvalidCommand,
            ValidateOnPresentation,
        ),
    )
    fail_endpoint = prepare_command_endpoint(
        fail_schema;
        capacity=1,
        ordinal=2,
    )
    fail_state = CommandEndpointState(fail_endpoint)
    fail_workspace = CommandDispositionWorkspace(fail_endpoint)
    failing_command = PlantCommand(
        fail_schema,
        1,
        PlantTimestamp(1),
        Float32[0, NaN, 0],
    )
    _assert_command_admission_error(
        () -> admit_plant_command!(fail_workspace, fail_endpoint, fail_state,
            failing_command, PlantTimestamp(1)),
        :admission,
        :nonfinite_failure,
    )
    failed = _assert_single_disposition(
        fail_workspace,
        FailedCommand,
        :nonfinite_failure,
    )
    @test command_presentation_id(failed) == CommandPresentationID(1)
    @test command_endpoint_timestamp(fail_state) == PlantTimestamp(1)

    clip_schema = _admission_schema(
        Float32;
        id=:clip_schema,
        endpoint=:clip_command,
        value_policy=CommandValuePolicy(
            RejectInvalidCommand,
            ClipInvalidCommand,
            ValidateOnPresentation,
        ),
    )
    clip_endpoint = prepare_command_endpoint(
        clip_schema;
        capacity=1,
        ordinal=3,
    )
    clip_state = CommandEndpointState(clip_endpoint)
    clip_workspace = CommandDispositionWorkspace(clip_endpoint)
    caller_payload = Float32[-2, 0.25, 2]
    command = PlantCommand(clip_schema, 1, PlantTimestamp(10), caller_payload)
    admission = @inferred admit_plant_command!(
        clip_workspace,
        clip_endpoint,
        clip_state,
        command,
        PlantTimestamp(10),
    )
    @test command_admission_status(admission) == CommandAdmittedReady
    @test command_sequence_class(admission) == InOrderCommandSequence
    @test command_disposition_count(clip_workspace) == 0
    caller_payload .= 9

    claim = claim_next_application_ready_command!(
        clip_endpoint,
        clip_state,
        PlantTimestamp(10),
    )
    @test claim isa PlantCommandApplicationClaim
    @test claimed_command_payload(clip_endpoint, clip_state, claim) ==
        Float32[-1, 0.25, 1]
    @test command_requested_effective_timestamp(claim) == PlantTimestamp(10)
    @test command_admission_timestamp(claim) == PlantTimestamp(10)
    @test command_scheduled_timestamp(claim) == PlantTimestamp(10)
    @test command_ready_timestamp(claim) == PlantTimestamp(10)
    @test iszero(command_lateness(claim))
    @test !applicable(
        PlantCommandApplicationClaim,
        getfield(claim, :binding_id),
        getfield(claim, :slot),
        getfield(claim, :generation),
        command_presentation_id(claim),
        command_order_key(claim),
        command_requested_effective_timestamp(claim),
        command_admission_timestamp(claim),
        command_ready_timestamp(claim),
    )
    forged_claims = (
        _rebuilt_application_claim(
            claim;
            key=PlantCommandOrderKey(
                PlantTimestamp(0),
                UInt32(99),
                PlantCommandSequence(999),
            ),
        ),
        _rebuilt_application_claim(
            claim;
            requested=PlantTimestamp(0),
        ),
        _rebuilt_application_claim(
            claim;
            admission=PlantTimestamp(0),
        ),
        _rebuilt_application_claim(
            claim;
            ready=PlantTimestamp(0),
        ),
    )
    _assert_command_admission_error(
        () -> mark_plant_command_applied!(
            clip_workspace,
            clip_endpoint,
            clip_state,
            forged_claims[1],
        ),
        :application,
        :stale_claim,
    )
    @test command_disposition_count(clip_workspace) == 0
    @test active_command_count(clip_state) == 1
    for index in 2:length(forged_claims)
        forged_claim = forged_claims[index]
        _assert_command_admission_error(
            () -> claimed_command_payload(
                clip_endpoint,
                clip_state,
                forged_claim,
            ),
            :application,
            :stale_claim,
        )
    end
    disposition = @inferred mark_plant_command_applied!(
        clip_workspace,
        clip_endpoint,
        clip_state,
        claim,
    )
    @test disposition === command_disposition(clip_workspace, 1)
    @test command_terminal_kind(disposition) == AppliedCommand
    @test command_disposition_reason(disposition) ==
        CommandDispositionReason(:applied)
    @test command_terminal_timestamp(disposition) == PlantTimestamp(10)
    @test superseding_command_presentation_id(disposition) === nothing
    @test active_command_count(clip_state) == 0
    _assert_command_admission_error(
        () -> claimed_command_payload(clip_endpoint, clip_state, claim),
        :application,
        :stale_claim,
    )
end

@testset "Future, late, and no-backdating policies" begin
    future_reject_schema = _admission_schema(
        Float32;
        id=:future_reject_schema,
        endpoint=:future_reject_command,
        effective_time_policy=CommandEffectiveTimePolicy(
            future=RejectFutureCommand,
        ),
    )
    future_reject_endpoint = prepare_command_endpoint(
        future_reject_schema;
        capacity=1,
        ordinal=1,
    )
    future_reject_state = CommandEndpointState(future_reject_endpoint)
    future_reject_workspace = CommandDispositionWorkspace(
        future_reject_endpoint)
    future_command = PlantCommand(
        future_reject_schema,
        1,
        PlantTimestamp(10),
        zeros(Float32, 3),
    )
    future_result = admit_plant_command!(future_reject_workspace,
        future_reject_endpoint, future_reject_state, future_command,
        PlantTimestamp(1))
    @test command_admission_status(future_result) ==
        CommandTerminatedOnAdmission
    _assert_single_disposition(future_reject_workspace, RejectedCommand,
        :future_command)

    future_schema = _admission_schema(
        Float32;
        id=:future_schema,
        endpoint=:future_command,
    )
    future_endpoint = prepare_command_endpoint(
        future_schema;
        capacity=1,
        ordinal=2,
    )
    future_state = CommandEndpointState(future_endpoint)
    future_workspace = CommandDispositionWorkspace(future_endpoint)
    pending = admit_plant_command!(future_workspace, future_endpoint,
        future_state,
        PlantCommand(future_schema, 1, PlantTimestamp(10),
            zeros(Float32, 3)),
        PlantTimestamp(1))
    @test command_admission_status(pending) == CommandAdmittedPending
    @test command_scheduled_timestamp(command_order_key(pending)) ==
        PlantTimestamp(10)
    @test claim_next_application_ready_command!(future_endpoint,
        future_state, PlantTimestamp(5)) === nothing
    @test command_endpoint_timestamp(future_state) == PlantTimestamp(5)
    future_claim = claim_next_application_ready_command!(future_endpoint,
        future_state, PlantTimestamp(10))
    fail_plant_command_application!(future_workspace, future_endpoint,
        future_state, future_claim; reason=:declared_device_failure)
    _assert_single_disposition(future_workspace, FailedCommand,
        :declared_device_failure)

    for (late_policy, expected_kind, expected_reason) in (
        (RejectLateCommand, RejectedCommand, :late_command),
        (FailOnLateCommand, FailedCommand, :late_command),
    )
        schema = _admission_schema(
            Float32;
            id=:late_policy_schema,
            endpoint=:late_policy_command,
            effective_time_policy=CommandEffectiveTimePolicy(
                late=late_policy,
            ),
        )
        endpoint = prepare_command_endpoint(schema; capacity=1, ordinal=3)
        state = CommandEndpointState(endpoint)
        workspace = CommandDispositionWorkspace(endpoint)
        command = PlantCommand(
            schema,
            1,
            PlantTimestamp(5),
            zeros(Float32, 3),
        )
        operation = () -> admit_plant_command!(workspace, endpoint, state,
            command, PlantTimestamp(10))
        if late_policy == FailOnLateCommand
            _assert_command_admission_error(operation, :admission,
                :late_command)
        else
            result = operation()
            @test command_admission_status(result) ==
                CommandTerminatedOnAdmission
        end
        disposition = _assert_single_disposition(workspace, expected_kind,
            expected_reason)
        @test command_requested_effective_timestamp(disposition) ==
            PlantTimestamp(5)
        @test command_terminal_timestamp(disposition) == PlantTimestamp(10)
        @test command_lateness(disposition) == PlantDuration(5)
    end

    late_now_schema = _admission_schema(
        Float32;
        id=:late_now_schema,
        endpoint=:late_now_command,
        effective_time_policy=CommandEffectiveTimePolicy(
            late=ApplyLateCommandNow,
        ),
    )
    late_now_endpoint = prepare_command_endpoint(
        late_now_schema;
        capacity=1,
        ordinal=4,
    )
    late_now_state = CommandEndpointState(late_now_endpoint)
    late_now_workspace = CommandDispositionWorkspace(late_now_endpoint)
    late_now = admit_plant_command!(late_now_workspace, late_now_endpoint,
        late_now_state,
        PlantCommand(late_now_schema, 1, PlantTimestamp(5),
            zeros(Float32, 3)),
        PlantTimestamp(10))
    @test command_admission_status(late_now) == CommandAdmittedReady
    @test command_scheduled_timestamp(command_order_key(late_now)) ==
        PlantTimestamp(10)
    late_claim = claim_next_application_ready_command!(late_now_endpoint,
        late_now_state, PlantTimestamp(10))
    @test command_requested_effective_timestamp(late_claim) ==
        PlantTimestamp(5)
    @test command_scheduled_timestamp(late_claim) == PlantTimestamp(10)
    @test command_ready_timestamp(late_claim) == PlantTimestamp(10)
    @test command_lateness(late_claim) == PlantDuration(5)
    mark_plant_command_applied!(late_now_workspace, late_now_endpoint,
        late_now_state, late_claim)
    late_disposition = command_disposition(late_now_workspace, 1)
    @test command_requested_effective_timestamp(late_disposition) ==
        PlantTimestamp(5)
    @test command_terminal_timestamp(late_disposition) == PlantTimestamp(10)
    @test command_lateness(late_disposition) == PlantDuration(5)
end

@testset "Sequence window, capacity, retry, and supersession" begin
    sequence_schema = _admission_schema(
        Float32;
        id=:sequence_schema,
        endpoint=:sequence_command,
        sequence_policy=CommandSequencePolicy(
            RejectSequence,
            RejectSequence,
            AcceptSequence,
            AcceptSequence,
        ),
    )
    endpoint = prepare_command_endpoint(
        sequence_schema;
        capacity=5,
        sequence_window=2,
        ordinal=1,
    )
    state = CommandEndpointState(endpoint)
    workspace = CommandDispositionWorkspace(endpoint)
    payload = zeros(Float32, 3)

    for (sequence, expected_class) in (
        (10, InOrderCommandSequence),
        (12, SkippedCommandSequence),
        (11, ReorderedCommandSequence),
    )
        result = admit_plant_command!(workspace, endpoint, state,
            PlantCommand(sequence_schema, sequence, PlantTimestamp(100),
                payload),
            PlantTimestamp(0))
        @test command_admission_status(result) == CommandAdmittedPending
        @test command_sequence_class(result) == expected_class
    end
    duplicate = admit_plant_command!(workspace, endpoint, state,
        PlantCommand(sequence_schema, 11, PlantTimestamp(100), payload),
        PlantTimestamp(0))
    @test command_sequence_class(duplicate) == DuplicateCommandSequence
    _assert_single_disposition(workspace, RejectedCommand,
        :duplicate_sequence)
    clear_command_dispositions!(workspace)
    stale = admit_plant_command!(workspace, endpoint, state,
        PlantCommand(sequence_schema, 9, PlantTimestamp(100), payload),
        PlantTimestamp(0))
    @test command_sequence_class(stale) == StaleCommandSequence
    _assert_single_disposition(workspace, RejectedCommand, :stale_sequence)

    sequence_fail_schema = _admission_schema(
        Float32;
        id=:sequence_fail_schema,
        endpoint=:sequence_fail_command,
        sequence_policy=CommandSequencePolicy(
            duplicate=FailOnSequence,
        ),
    )
    sequence_fail_endpoint = prepare_command_endpoint(
        sequence_fail_schema;
        capacity=1,
        ordinal=4,
    )
    sequence_fail_state = CommandEndpointState(sequence_fail_endpoint)
    sequence_fail_workspace = CommandDispositionWorkspace(
        sequence_fail_endpoint)
    sequence_fail_command = PlantCommand(sequence_fail_schema, 1,
        PlantTimestamp(100), payload)
    admit_plant_command!(sequence_fail_workspace, sequence_fail_endpoint,
        sequence_fail_state, sequence_fail_command, PlantTimestamp(0))
    _assert_command_admission_error(
        () -> admit_plant_command!(sequence_fail_workspace,
            sequence_fail_endpoint, sequence_fail_state,
            sequence_fail_command, PlantTimestamp(0)),
        :admission,
        :duplicate_sequence,
    )
    _assert_single_disposition(sequence_fail_workspace, FailedCommand,
        :duplicate_sequence)

    capacity_schema = _admission_schema(
        Float64;
        id=:capacity_schema,
        endpoint=:capacity_command,
        dimensions=(),
        bounds=UniformCommandBounds(-1.0, 1.0),
    )
    capacity_endpoint = prepare_command_endpoint(
        capacity_schema;
        capacity=1,
        ordinal=2,
    )
    capacity_state = CommandEndpointState(capacity_endpoint)
    capacity_workspace = CommandDispositionWorkspace(capacity_endpoint)
    first = admit_plant_command!(capacity_workspace, capacity_endpoint,
        capacity_state,
        PlantCommand(capacity_schema, 1, PlantTimestamp(100), 0.25),
        PlantTimestamp(0))
    first_key = command_order_key(first)
    rejected = admit_plant_command!(capacity_workspace, capacity_endpoint,
        capacity_state,
        PlantCommand(capacity_schema, 2, PlantTimestamp(200), 0.75),
        PlantTimestamp(0))
    @test command_sequence_class(rejected) == InOrderCommandSequence
    @test command_presentation_id(rejected) == CommandPresentationID(2)
    _assert_single_disposition(capacity_workspace, RejectedCommand,
        :calendar_capacity)
    @test pending_command_count(capacity_state) == 1
    @test active_command_count(capacity_state) == 1
    @test next_command_order_key(capacity_endpoint, capacity_state) ==
        first_key
    clear_command_dispositions!(capacity_workspace)
    first_claim = claim_next_application_ready_command!(capacity_endpoint,
        capacity_state, PlantTimestamp(100))
    @test claimed_command_payload(capacity_endpoint, capacity_state,
        first_claim) == 0.25
    mark_plant_command_applied!(capacity_workspace, capacity_endpoint,
        capacity_state, first_claim)
    clear_command_dispositions!(capacity_workspace)
    retry = admit_plant_command!(capacity_workspace, capacity_endpoint,
        capacity_state,
        PlantCommand(capacity_schema, 2, PlantTimestamp(200), 0.75),
        PlantTimestamp(100))
    @test command_sequence_class(retry) == InOrderCommandSequence
    @test command_presentation_id(retry) == CommandPresentationID(3)
    @test command_admission_status(retry) == CommandAdmittedPending

    supersession_schema = _admission_schema(
        Float32;
        id=:supersession_schema,
        endpoint=:supersession_command,
        sequence_policy=CommandSequencePolicy(reordered=AcceptSequence),
        effective_time_policy=CommandEffectiveTimePolicy(
            supersession=SupersedeOlderPendingCommands,
        ),
    )
    supersession_endpoint = prepare_command_endpoint(
        supersession_schema;
        capacity=3,
        ordinal=3,
    )
    supersession_state = CommandEndpointState(supersession_endpoint)
    supersession_workspace = CommandDispositionWorkspace(
        supersession_endpoint)
    admit_plant_command!(supersession_workspace, supersession_endpoint,
        supersession_state,
        PlantCommand(supersession_schema, 3, PlantTimestamp(100), payload),
        PlantTimestamp(0))
    admit_plant_command!(supersession_workspace, supersession_endpoint,
        supersession_state,
        PlantCommand(supersession_schema, 2, PlantTimestamp(80), payload),
        PlantTimestamp(0))
    replacing = admit_plant_command!(supersession_workspace,
        supersession_endpoint, supersession_state,
        PlantCommand(supersession_schema, 4, PlantTimestamp(90), payload),
        PlantTimestamp(0))
    @test command_presentation_id(replacing) == CommandPresentationID(3)
    @test command_disposition_count(supersession_workspace) == 2
    @test command_sequence(command_disposition(
        supersession_workspace, 1)) == PlantCommandSequence(2)
    @test command_sequence(command_disposition(
        supersession_workspace, 2)) == PlantCommandSequence(3)
    for index in 1:2
        displaced = command_disposition(supersession_workspace, index)
        @test command_terminal_kind(displaced) == SupersededCommand
        @test command_disposition_reason(displaced) ==
            CommandDispositionReason(:superseded_by_newer_sequence)
        @test superseding_command_presentation_id(displaced) ==
            command_presentation_id(replacing)
    end
    @test pending_command_count(supersession_state) == 1
    @test active_command_count(supersession_state) == 1
end

@testset "Deterministic order and terminal lifecycle" begin
    ordered_schema = _admission_schema(
        Float64;
        id=:ordered_schema,
        endpoint=:ordered_command,
        dimensions=(),
        bounds=UniformCommandBounds(-1.0, 1.0),
        sequence_policy=CommandSequencePolicy(reordered=AcceptSequence),
    )
    endpoint = prepare_command_endpoint(
        ordered_schema;
        capacity=3,
        ordinal=4,
    )
    state = CommandEndpointState(endpoint)
    workspace = CommandDispositionWorkspace(endpoint)
    admit_plant_command!(workspace, endpoint, state,
        PlantCommand(ordered_schema, 2, PlantTimestamp(50), 0.2),
        PlantTimestamp(0))
    admit_plant_command!(workspace, endpoint, state,
        PlantCommand(ordered_schema, 1, PlantTimestamp(50), 0.1),
        PlantTimestamp(0))
    admit_plant_command!(workspace, endpoint, state,
        PlantCommand(ordered_schema, 3, PlantTimestamp(40), 0.3),
        PlantTimestamp(0))
    @test command_sequence(next_command_order_key(endpoint, state)) ==
        PlantCommandSequence(3)
    claim_three = claim_next_application_ready_command!(endpoint, state,
        PlantTimestamp(40))
    @test command_sequence(claim_three) == PlantCommandSequence(3)
    mark_plant_command_applied!(workspace, endpoint, state, claim_three)
    clear_command_dispositions!(workspace)
    claim_one = claim_next_application_ready_command!(endpoint, state,
        PlantTimestamp(50))
    @test command_sequence(claim_one) == PlantCommandSequence(1)
    mark_plant_command_applied!(workspace, endpoint, state, claim_one)
    clear_command_dispositions!(workspace)
    claim_two = claim_next_application_ready_command!(endpoint, state,
        PlantTimestamp(50))
    @test command_sequence(claim_two) == PlantCommandSequence(2)
    failed = fail_plant_command_application!(workspace, endpoint, state,
        claim_two; reason=:optic_rejected_command)
    @test command_terminal_kind(failed) == FailedCommand
    @test command_disposition_reason(failed) ==
        CommandDispositionReason(:optic_rejected_command)
    @test active_command_count(state) == 0

    schema_b = _admission_schema(
        Float64;
        id=:endpoint_b_schema,
        endpoint=:endpoint_b,
        dimensions=(),
        bounds=UniformCommandBounds(-1.0, 1.0),
    )
    schema_a = _admission_schema(
        Float64;
        id=:endpoint_a_schema,
        endpoint=:endpoint_a,
        dimensions=(),
        bounds=UniformCommandBounds(-1.0, 1.0),
    )
    endpoint_b = prepare_command_endpoint(schema_b; capacity=1, ordinal=2)
    endpoint_a = prepare_command_endpoint(schema_a; capacity=1, ordinal=1)
    state_b = CommandEndpointState(endpoint_b)
    state_a = CommandEndpointState(endpoint_a)
    workspace_b = CommandDispositionWorkspace(endpoint_b)
    workspace_a = CommandDispositionWorkspace(endpoint_a)
    key_b = command_order_key(admit_plant_command!(workspace_b, endpoint_b,
        state_b, PlantCommand(schema_b, 1, PlantTimestamp(25), 0.0),
        PlantTimestamp(0)))
    key_a = command_order_key(admit_plant_command!(workspace_a, endpoint_a,
        state_a, PlantCommand(schema_a, 100, PlantTimestamp(25), 0.0),
        PlantTimestamp(0)))
    @test isless(key_a, key_b)
    @test !isless(key_b, key_a)
    @test command_endpoint_ordinal(key_a) == UInt32(1)
    @test command_endpoint_ordinal(key_b) == UInt32(2)

    drain_schema = _admission_schema(
        Float64;
        id=:drain_schema,
        endpoint=:drain_command,
        dimensions=(),
        bounds=UniformCommandBounds(-1.0, 1.0),
    )
    drain_endpoint = prepare_command_endpoint(
        drain_schema;
        capacity=3,
        ordinal=5,
    )
    drain_state = CommandEndpointState(drain_endpoint)
    drain_workspace = CommandDispositionWorkspace(drain_endpoint)
    for (sequence, timestamp) in ((1, 30), (2, 10), (3, 20))
        admit_plant_command!(drain_workspace, drain_endpoint, drain_state,
            PlantCommand(drain_schema, sequence, PlantTimestamp(timestamp),
                0.0),
            PlantTimestamp(0))
    end
    @test fail_pending_plant_commands!(drain_workspace, drain_endpoint,
        drain_state, PlantTimestamp(5); reason=:run_stopped) == 3
    @test command_sequence.(Tuple(command_disposition(
        drain_workspace, index) for index in 1:3)) == (
        PlantCommandSequence(2),
        PlantCommandSequence(3),
        PlantCommandSequence(1),
    )
    @test all(command_terminal_kind(command_disposition(
        drain_workspace, index)) == FailedCommand for index in 1:3)
    @test all(command_disposition_reason(command_disposition(
        drain_workspace, index)) == CommandDispositionReason(:run_stopped)
        for index in 1:3)
    @test pending_command_count(drain_state) == 0
    @test active_command_count(drain_state) == 0
end

@testset "Endpoint ownership and unconsumed disposition guards" begin
    schema_a = _admission_schema(
        Float64;
        id=:guard_a_schema,
        endpoint=:guard_a_command,
        dimensions=(),
        bounds=UniformCommandBounds(-1.0, 1.0),
        effective_time_policy=CommandEffectiveTimePolicy(
            future=RejectFutureCommand,
        ),
    )
    schema_b = _admission_schema(
        Float64;
        id=:guard_b_schema,
        endpoint=:guard_b_command,
        dimensions=(),
        bounds=UniformCommandBounds(-1.0, 1.0),
    )
    endpoint_a = prepare_command_endpoint(schema_a; capacity=1, ordinal=1)
    endpoint_b = prepare_command_endpoint(schema_b; capacity=1, ordinal=2)
    state_a = CommandEndpointState(endpoint_a)
    state_b = CommandEndpointState(endpoint_b)
    workspace_a = CommandDispositionWorkspace(endpoint_a)
    workspace_b = CommandDispositionWorkspace(endpoint_b)
    command_a = PlantCommand(schema_a, 1, PlantTimestamp(10), 0.0)

    _assert_command_admission_error(
        () -> admit_plant_command!(workspace_a, endpoint_a, state_b,
            command_a, PlantTimestamp(0)),
        :endpoint,
        :foreign_state,
    )
    _assert_command_admission_error(
        () -> admit_plant_command!(workspace_b, endpoint_a, state_a,
            command_a, PlantTimestamp(0)),
        :endpoint,
        :foreign_workspace,
    )
    admit_plant_command!(workspace_a, endpoint_a, state_a, command_a,
        PlantTimestamp(5))
    _assert_command_admission_error(
        () -> admit_plant_command!(workspace_a, endpoint_a, state_a,
            command_a, PlantTimestamp(5)),
        :disposition,
        :unconsumed_dispositions,
    )
    clear_command_dispositions!(workspace_a)
    _assert_command_admission_error(
        () -> claim_next_application_ready_command!(endpoint_a, state_a,
            PlantTimestamp(0)),
        :endpoint,
        :time_regression,
    )
end

@testset "Command endpoint inference and warmed allocation" begin
    scalar_schema = _admission_schema(
        Float64;
        id=:scalar_hot_schema,
        endpoint=:scalar_hot_command,
        dimensions=(),
        bounds=UniformCommandBounds(-1.0, 1.0),
    )
    scalar_endpoint = prepare_command_endpoint(
        scalar_schema;
        capacity=1,
        ordinal=1,
    )
    scalar_state = CommandEndpointState(scalar_endpoint)
    scalar_workspace = CommandDispositionWorkspace(scalar_endpoint)
    @test @inferred(_scalar_command_cycle!(scalar_workspace,
        scalar_endpoint, scalar_state, 1, PlantTimestamp(1))) === nothing

    array_schema = _admission_schema(
        Float32;
        id=:array_hot_schema,
        endpoint=:array_hot_command,
        value_policy=CommandValuePolicy(
            RejectInvalidCommand,
            ClipInvalidCommand,
            ValidateOnPresentation,
        ),
    )
    array_endpoint = prepare_command_endpoint(
        array_schema;
        capacity=1,
        ordinal=2,
    )
    array_state = CommandEndpointState(array_endpoint)
    array_workspace = CommandDispositionWorkspace(array_endpoint)
    array_payload = Float32[-2, 0.5, 2]
    @test @inferred(_array_command_cycle!(array_workspace, array_endpoint,
        array_state, array_payload, 1, PlantTimestamp(1))) === nothing

    multi_entry_schema = _admission_schema(
        Float64;
        id=:multi_entry_hot_schema,
        endpoint=:multi_entry_hot_command,
        dimensions=(),
        bounds=UniformCommandBounds(-1.0, 1.0),
    )
    multi_entry_endpoint = prepare_command_endpoint(
        multi_entry_schema;
        capacity=3,
        ordinal=3,
    )
    multi_entry_state = CommandEndpointState(multi_entry_endpoint)
    multi_entry_workspace = CommandDispositionWorkspace(multi_entry_endpoint)
    @test @inferred(_multi_entry_command_cycle!(multi_entry_workspace,
        multi_entry_endpoint, multi_entry_state, 1,
        PlantTimestamp(10))) === nothing

    supersession_schema = _admission_schema(
        Float64;
        id=:supersession_hot_schema,
        endpoint=:supersession_hot_command,
        dimensions=(),
        bounds=UniformCommandBounds(-1.0, 1.0),
        sequence_policy=CommandSequencePolicy(reordered=AcceptSequence),
        effective_time_policy=CommandEffectiveTimePolicy(
            supersession=SupersedeOlderPendingCommands,
        ),
    )
    supersession_endpoint = prepare_command_endpoint(
        supersession_schema;
        capacity=3,
        ordinal=4,
    )
    supersession_state = CommandEndpointState(supersession_endpoint)
    supersession_workspace = CommandDispositionWorkspace(
        supersession_endpoint)
    @test @inferred(_supersession_command_cycle!(supersession_workspace,
        supersession_endpoint, supersession_state, 1,
        PlantTimestamp(10))) === nothing

    capacity_schema = _admission_schema(
        Float64;
        id=:capacity_hot_schema,
        endpoint=:capacity_hot_command,
        dimensions=(),
        bounds=UniformCommandBounds(-1.0, 1.0),
    )
    capacity_endpoint = prepare_command_endpoint(
        capacity_schema;
        capacity=1,
        ordinal=5,
    )
    capacity_state = CommandEndpointState(capacity_endpoint)
    capacity_workspace = CommandDispositionWorkspace(capacity_endpoint)
    @test @inferred(_capacity_rejection_command_cycle!(capacity_workspace,
        capacity_endpoint, capacity_state, 1,
        PlantTimestamp(10))) === nothing

    drain_schema = _admission_schema(
        Float64;
        id=:drain_hot_schema,
        endpoint=:drain_hot_command,
        dimensions=(),
        bounds=UniformCommandBounds(-1.0, 1.0),
    )
    drain_endpoint = prepare_command_endpoint(
        drain_schema;
        capacity=3,
        ordinal=6,
    )
    drain_state = CommandEndpointState(drain_endpoint)
    drain_workspace = CommandDispositionWorkspace(drain_endpoint)
    @test @inferred(_failure_drain_command_cycle!(drain_workspace,
        drain_endpoint, drain_state, 1,
        PlantTimestamp(10))) === nothing

    if coverage_instrumented()
        @test_skip "command-admission allocation assertions disabled under coverage"
    else
        @test @allocated(_scalar_command_cycle!(scalar_workspace,
            scalar_endpoint, scalar_state, 2, PlantTimestamp(2))) == 0
        @test @allocated(_array_command_cycle!(array_workspace,
            array_endpoint, array_state, array_payload, 2,
            PlantTimestamp(2))) == 0
        @test @allocated(_multi_entry_command_cycle!(multi_entry_workspace,
            multi_entry_endpoint, multi_entry_state, 4,
            PlantTimestamp(50))) == 0
        @test @allocated(_supersession_command_cycle!(
            supersession_workspace, supersession_endpoint,
            supersession_state, 4, PlantTimestamp(50))) == 0
        @test @allocated(_capacity_rejection_command_cycle!(
            capacity_workspace, capacity_endpoint, capacity_state, 2,
            PlantTimestamp(50))) == 0
        @test @allocated(_failure_drain_command_cycle!(drain_workspace,
            drain_endpoint, drain_state, 4,
            PlantTimestamp(50))) == 0
    end

    scalar_slot_count = length(getfield(scalar_state, :slots))
    scalar_history_count = length(getfield(scalar_state,
        :accepted_sequences))
    for sequence in 3:1_002
        _scalar_command_cycle!(scalar_workspace, scalar_endpoint,
            scalar_state, sequence, PlantTimestamp(sequence))
    end
    @test length(getfield(scalar_state, :slots)) == scalar_slot_count
    @test length(getfield(scalar_state, :accepted_sequences)) ==
        scalar_history_count
    @test pending_command_count(scalar_state) == 0
    @test active_command_count(scalar_state) == 0
end

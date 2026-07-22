function _command_application_schema(::Type{T}=Float64;
    id=:application_schema,
    endpoint=:application_command,
    dimensions=(),
    semantics=AbsoluteCommand,
    bounds=UniformCommandBounds(T(-1), T(1)),
    value_policy=CommandValuePolicy(),
    silence_policy=CommandSilencePolicy(),
) where {T}
    basis = isempty(dimensions) ?
        CommandBasis(:rigid_body, :focus) :
        CommandBasis(:actuator, :application_test)
    return PlantCommandSchema(
        T,
        dimensions;
        id,
        version=1,
        endpoint,
        units=:metre,
        sign_convention=:positive_surface_increases_opd,
        basis,
        basis_revision=1,
        semantics,
        bounds,
        value_policy,
        sequence_policy=CommandSequencePolicy(),
        effective_time_policy=CommandEffectiveTimePolicy(),
        silence_policy,
    )
end

function _captured_command_application_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

function _assert_command_application_error(f, stage::Symbol, reason::Symbol)
    error = _captured_command_application_error(f)
    @test error isa PlantCommandError
    if error isa PlantCommandError
        @test error.stage === stage
        @test error.reason === reason
        @test !isempty(error.msg)
    end
    return error
end

function _assert_application_disposition(workspace, kind, reason)
    @test command_disposition_count(workspace) == 1
    disposition = command_disposition(workspace, 1)
    @test command_terminal_kind(disposition) == kind
    @test command_disposition_reason(disposition) ==
        CommandDispositionReason(reason)
    return disposition
end

function _apply_scalar_command_cycle!(workspace, endpoint, endpoint_state,
    application_state, sequence, timestamp, payload)
    command = PlantCommand(command_schema(endpoint), sequence, timestamp,
        payload)
    admission = admit_plant_command!(workspace, endpoint, endpoint_state,
        command, timestamp)
    command_admission_status(admission) == CommandAdmittedReady ||
        return nothing
    claim = claim_next_application_ready_command!(endpoint, endpoint_state,
        timestamp)
    claim === nothing && return nothing
    apply_claimed_plant_command!(workspace, endpoint, endpoint_state,
        application_state, claim)
    clear_command_dispositions!(workspace)
    return nothing
end

function _apply_array_command_cycle!(workspace, endpoint, endpoint_state,
    application_state, sequence, timestamp, payload)
    command = PlantCommand(command_schema(endpoint), sequence, timestamp,
        payload)
    admission = admit_plant_command!(workspace, endpoint, endpoint_state,
        command, timestamp)
    command_admission_status(admission) == CommandAdmittedReady ||
        return nothing
    claim = claim_next_application_ready_command!(endpoint, endpoint_state,
        timestamp)
    claim === nothing && return nothing
    apply_claimed_plant_command!(workspace, endpoint, endpoint_state,
        application_state, claim)
    clear_command_dispositions!(workspace)
    return nothing
end

function _safe_command_cycle!(workspace, endpoint, endpoint_state,
    application_state, sequence, application_timestamp)
    _apply_scalar_command_cycle!(workspace, endpoint, endpoint_state,
        application_state, sequence, application_timestamp, 0.25)
    deadline = next_command_silence_timestamp(endpoint, endpoint_state,
        application_state)
    apply_command_silence_transition!(workspace, endpoint, endpoint_state,
        application_state, deadline)
    return nothing
end

@testset "Effective-command application state and API" begin
    for name in (
        :PlantCommandSilenceTransition,
        :effective_command,
        :last_command_admission_timestamp,
        :last_command_application_timestamp,
        :command_endpoint_failed,
        :apply_claimed_plant_command!,
        :next_command_silence_timestamp,
        :apply_command_silence_transition!,
        :command_silence_action,
        :command_silence_age_origin,
        :command_silence_origin_timestamp,
        :command_silence_deadline,
        :command_silence_transition_timestamp,
        :command_silence_age,
    )
        @test Base.isexported(AdaptiveOpticsSim, name)
    end
    @test !Base.isexported(AdaptiveOpticsSim, :CommandApplicationState)

    schema = _command_application_schema()
    endpoint = prepare_command_endpoint(schema; capacity=2, ordinal=1)
    endpoint_state = CommandEndpointState(
        endpoint;
        initial_timestamp=PlantTimestamp(5),
    )
    application_state = @inferred CommandApplicationState(
        endpoint, endpoint_state, 0.25)
    @test Base.ismutable(application_state)
    @test effective_command(application_state) == 0.25
    @test last_command_application_timestamp(application_state) ==
        PlantTimestamp(5)
    @test last_command_admission_timestamp(endpoint_state) === nothing
    @test !command_endpoint_failed(endpoint_state)
    @test next_command_silence_timestamp(endpoint, endpoint_state,
        application_state) === nothing

    array_schema = _command_application_schema(
        Float32;
        id=:array_application_schema,
        endpoint=:array_application_command,
        dimensions=(3,),
    )
    array_endpoint = prepare_command_endpoint(
        array_schema;
        capacity=1,
        ordinal=2,
    )
    array_endpoint_state = CommandEndpointState(array_endpoint)
    initial = Float32[0.1, 0.2, 0.3]
    array_application_state = CommandApplicationState(
        array_endpoint, array_endpoint_state, initial)
    initial .= 9
    @test effective_command(array_application_state) ==
        Float32[0.1, 0.2, 0.3]
    @test backend(effective_command(array_application_state)) == CPUBackend()

    safe_array_schema = _command_application_schema(
        Float32;
        id=:safe_array_application_schema,
        endpoint=:safe_array_application_command,
        dimensions=(3,),
        silence_policy=CommandSilencePolicy(
            ApplySafeCommand,
            AgeFromApplication;
            timeout=PlantDuration(10),
        ),
    )
    safe_array_endpoint = prepare_command_endpoint(
        safe_array_schema;
        capacity=1,
        ordinal=5,
    )
    safe_array_endpoint_state = CommandEndpointState(safe_array_endpoint)
    caller_safe = fill(Float32(-0.25), 3)
    safe_array_application_state = CommandApplicationState(
        safe_array_endpoint, safe_array_endpoint_state, zeros(Float32, 3);
        safe_command=caller_safe)
    safe_array_workspace = CommandDispositionWorkspace(safe_array_endpoint)
    caller_safe .= 9
    apply_command_silence_transition!(safe_array_workspace,
        safe_array_endpoint, safe_array_endpoint_state,
        safe_array_application_state, PlantTimestamp(10))
    @test effective_command(safe_array_application_state) ==
        fill(Float32(-0.25), 3)

    safe_schema = _command_application_schema(
        Float64;
        id=:safe_configuration_schema,
        endpoint=:safe_configuration_command,
        silence_policy=CommandSilencePolicy(
            ApplySafeCommand,
            AgeFromApplication;
            timeout=PlantDuration(10),
        ),
    )
    safe_endpoint = prepare_command_endpoint(
        safe_schema;
        capacity=1,
        ordinal=3,
    )
    safe_endpoint_state = CommandEndpointState(safe_endpoint)

    other_endpoint = prepare_command_endpoint(
        _command_application_schema(
            Float64;
            id=:other_application_schema,
            endpoint=:other_application_command,
        );
        capacity=1,
        ordinal=4,
    )
    other_endpoint_state = CommandEndpointState(other_endpoint)
    other_application_state = CommandApplicationState(
        other_endpoint, other_endpoint_state, 0.0)
    sibling_endpoint_state = CommandEndpointState(
        endpoint;
        initial_timestamp=PlantTimestamp(5),
    )
    sibling_application_state = CommandApplicationState(
        endpoint, sibling_endpoint_state, 0.25)

    for (operation, stage, reason) in (
        (() -> CommandApplicationState(safe_endpoint, safe_endpoint_state,
            0.0), :preparation, :missing_safe_command),
        (() -> CommandApplicationState(endpoint,
            CommandEndpointState(endpoint), 0.0; safe_command=0.0),
            :preparation, :unexpected_safe_command),
        (() -> CommandApplicationState(endpoint,
            CommandEndpointState(endpoint), 2.0),
            :preparation, :invalid_effective_command),
        (() -> CommandApplicationState(endpoint, CommandEndpointState(endpoint),
            reshape([0.0], ())), :preparation,
            :invalid_effective_command),
        (() -> CommandApplicationState(array_endpoint,
            CommandEndpointState(array_endpoint),
            zeros(Float64, 3)), :preparation, :invalid_effective_command),
        (() -> CommandApplicationState(array_endpoint,
            CommandEndpointState(array_endpoint),
            zeros(Float32, 2)), :preparation, :invalid_effective_command),
        (() -> CommandApplicationState(endpoint, endpoint_state, 0.0),
            :preparation, :application_state_exists),
        (() -> next_command_silence_timestamp(endpoint, endpoint_state,
            other_application_state), :application,
            :foreign_application_state),
        (() -> next_command_silence_timestamp(endpoint, endpoint_state,
            sibling_application_state), :application,
            :foreign_endpoint_state),
    )
        _assert_command_application_error(operation, stage, reason)
    end

    valid_safe_application_state = CommandApplicationState(
        safe_endpoint, safe_endpoint_state, 0.0; safe_command=-0.5)
    @test effective_command(valid_safe_application_state) == 0.0

    active_schema = _command_application_schema(
        Float64;
        id=:active_application_schema,
        endpoint=:active_application_command,
    )
    active_endpoint = prepare_command_endpoint(active_schema;
        capacity=1, ordinal=6)
    active_endpoint_state = CommandEndpointState(active_endpoint)
    active_workspace = CommandDispositionWorkspace(active_endpoint)
    admit_plant_command!(active_workspace, active_endpoint,
        active_endpoint_state,
        PlantCommand(active_schema, 1, PlantTimestamp(1), 0.25),
        PlantTimestamp(1))
    _assert_command_application_error(
        () -> CommandApplicationState(
            active_endpoint, active_endpoint_state, 0.0),
        :preparation,
        :active_endpoint_state,
    )
end

@testset "Absolute, incremental, clipped, and failed application" begin
    absolute_schema = _command_application_schema()
    absolute_endpoint = prepare_command_endpoint(
        absolute_schema;
        capacity=1,
        ordinal=1,
    )
    absolute_endpoint_state = CommandEndpointState(absolute_endpoint)
    absolute_application_state = CommandApplicationState(
        absolute_endpoint, absolute_endpoint_state, 0.1)
    absolute_workspace = CommandDispositionWorkspace(absolute_endpoint)
    absolute_admission = admit_plant_command!(absolute_workspace,
        absolute_endpoint, absolute_endpoint_state,
        PlantCommand(absolute_schema, 1, PlantTimestamp(10), 0.75),
        PlantTimestamp(10))
    @test last_command_admission_timestamp(absolute_endpoint_state) ==
        PlantTimestamp(10)
    absolute_claim = claim_next_application_ready_command!(absolute_endpoint,
        absolute_endpoint_state, PlantTimestamp(10))
    absolute_disposition = @inferred apply_claimed_plant_command!(
        absolute_workspace,
        absolute_endpoint,
        absolute_endpoint_state,
        absolute_application_state,
        absolute_claim,
    )
    @test command_presentation_id(absolute_disposition) ==
        command_presentation_id(absolute_admission)
    @test command_terminal_kind(absolute_disposition) == AppliedCommand
    @test command_disposition_reason(absolute_disposition) ==
        CommandDispositionReason(:applied)
    @test effective_command(absolute_application_state) == 0.75
    @test last_command_application_timestamp(absolute_application_state) ==
        PlantTimestamp(10)
    @test active_command_count(absolute_endpoint_state) == 0

    incremental_schema = _command_application_schema(
        Float64;
        id=:incremental_application_schema,
        endpoint=:incremental_application_command,
        semantics=IncrementalCommand,
        value_policy=CommandValuePolicy(
            out_of_range=RejectInvalidCommand,
            range_stage=EnforceOnApplication,
        ),
    )
    incremental_endpoint = prepare_command_endpoint(
        incremental_schema;
        capacity=1,
        ordinal=2,
    )
    incremental_endpoint_state = CommandEndpointState(incremental_endpoint)
    incremental_application_state = CommandApplicationState(
        incremental_endpoint, incremental_endpoint_state, 0.5)
    incremental_workspace = CommandDispositionWorkspace(incremental_endpoint)
    admit_plant_command!(incremental_workspace, incremental_endpoint,
        incremental_endpoint_state,
        PlantCommand(incremental_schema, 1, PlantTimestamp(1), 0.25),
        PlantTimestamp(1))
    incremental_claim = claim_next_application_ready_command!(
        incremental_endpoint, incremental_endpoint_state, PlantTimestamp(1))
    apply_claimed_plant_command!(incremental_workspace, incremental_endpoint,
        incremental_endpoint_state, incremental_application_state,
        incremental_claim)
    @test effective_command(incremental_application_state) == 0.75
    clear_command_dispositions!(incremental_workspace)

    admit_plant_command!(incremental_workspace, incremental_endpoint,
        incremental_endpoint_state,
        PlantCommand(incremental_schema, 2, PlantTimestamp(2), 0.5),
        PlantTimestamp(2))
    rejected_claim = claim_next_application_ready_command!(
        incremental_endpoint, incremental_endpoint_state, PlantTimestamp(2))
    rejected = apply_claimed_plant_command!(incremental_workspace,
        incremental_endpoint, incremental_endpoint_state,
        incremental_application_state, rejected_claim)
    @test command_terminal_kind(rejected) == RejectedCommand
    @test command_disposition_reason(rejected) ==
        CommandDispositionReason(:out_of_range_rejected)
    @test effective_command(incremental_application_state) == 0.75

    clip_schema = _command_application_schema(
        Float64;
        id=:clip_application_schema,
        endpoint=:clip_application_command,
        semantics=IncrementalCommand,
        value_policy=CommandValuePolicy(
            out_of_range=ClipInvalidCommand,
            range_stage=EnforceOnApplication,
        ),
    )
    clip_endpoint = prepare_command_endpoint(clip_schema;
        capacity=1, ordinal=3)
    clip_endpoint_state = CommandEndpointState(clip_endpoint)
    clip_application_state = CommandApplicationState(
        clip_endpoint, clip_endpoint_state, 0.75)
    clip_workspace = CommandDispositionWorkspace(clip_endpoint)
    admit_plant_command!(clip_workspace, clip_endpoint, clip_endpoint_state,
        PlantCommand(clip_schema, 1, PlantTimestamp(1), 0.5),
        PlantTimestamp(1))
    clip_claim = claim_next_application_ready_command!(clip_endpoint,
        clip_endpoint_state, PlantTimestamp(1))
    clipped = apply_claimed_plant_command!(clip_workspace, clip_endpoint,
        clip_endpoint_state, clip_application_state, clip_claim)
    @test command_terminal_kind(clipped) == AppliedCommand
    @test command_disposition_reason(clipped) ==
        CommandDispositionReason(:applied_clipped)
    @test effective_command(clip_application_state) == 1.0

    fail_schema = _command_application_schema(
        Float64;
        id=:fail_application_schema,
        endpoint=:fail_application_command,
        semantics=IncrementalCommand,
        value_policy=CommandValuePolicy(
            out_of_range=FailOnInvalidCommand,
            range_stage=EnforceOnApplication,
        ),
    )
    fail_endpoint = prepare_command_endpoint(fail_schema;
        capacity=1, ordinal=4)
    fail_endpoint_state = CommandEndpointState(fail_endpoint)
    fail_application_state = CommandApplicationState(
        fail_endpoint, fail_endpoint_state, 0.75)
    fail_workspace = CommandDispositionWorkspace(fail_endpoint)
    admit_plant_command!(fail_workspace, fail_endpoint, fail_endpoint_state,
        PlantCommand(fail_schema, 1, PlantTimestamp(1), 0.5),
        PlantTimestamp(1))
    fail_claim = claim_next_application_ready_command!(fail_endpoint,
        fail_endpoint_state, PlantTimestamp(1))
    _assert_command_application_error(
        () -> apply_claimed_plant_command!(fail_workspace, fail_endpoint,
            fail_endpoint_state, fail_application_state, fail_claim),
        :application,
        :out_of_range_failure,
    )
    _assert_application_disposition(fail_workspace, FailedCommand,
        :out_of_range_failure)
    @test effective_command(fail_application_state) == 0.75

    finite_limit = floatmax(Float64)
    nonfinite_schema = _command_application_schema(
        Float64;
        id=:nonfinite_application_schema,
        endpoint=:nonfinite_application_command,
        semantics=IncrementalCommand,
        bounds=UniformCommandBounds(-finite_limit, finite_limit),
        value_policy=CommandValuePolicy(
            nonfinite=RejectInvalidCommand,
            out_of_range=RejectInvalidCommand,
            range_stage=EnforceOnApplication,
        ),
    )
    nonfinite_endpoint = prepare_command_endpoint(nonfinite_schema;
        capacity=1, ordinal=5)
    nonfinite_endpoint_state = CommandEndpointState(nonfinite_endpoint)
    nonfinite_application_state = CommandApplicationState(
        nonfinite_endpoint, nonfinite_endpoint_state, finite_limit)
    nonfinite_workspace = CommandDispositionWorkspace(nonfinite_endpoint)
    admit_plant_command!(nonfinite_workspace, nonfinite_endpoint,
        nonfinite_endpoint_state,
        PlantCommand(nonfinite_schema, 1, PlantTimestamp(1), finite_limit),
        PlantTimestamp(1))
    nonfinite_claim = claim_next_application_ready_command!(
        nonfinite_endpoint, nonfinite_endpoint_state, PlantTimestamp(1))
    nonfinite = apply_claimed_plant_command!(nonfinite_workspace,
        nonfinite_endpoint, nonfinite_endpoint_state,
        nonfinite_application_state, nonfinite_claim)
    @test command_terminal_kind(nonfinite) == RejectedCommand
    @test command_disposition_reason(nonfinite) ==
        CommandDispositionReason(:nonfinite_rejected)
    @test effective_command(nonfinite_application_state) == finite_limit

    missed_schema = _command_application_schema(
        Float64;
        id=:missed_application_schema,
        endpoint=:missed_application_command,
    )
    missed_endpoint = prepare_command_endpoint(missed_schema;
        capacity=1, ordinal=6)
    missed_endpoint_state = CommandEndpointState(missed_endpoint)
    missed_application_state = CommandApplicationState(
        missed_endpoint, missed_endpoint_state, 0.25)
    missed_workspace = CommandDispositionWorkspace(missed_endpoint)
    admit_plant_command!(missed_workspace, missed_endpoint,
        missed_endpoint_state,
        PlantCommand(missed_schema, 1, PlantTimestamp(5), 0.75),
        PlantTimestamp(0))
    missed_claim = claim_next_application_ready_command!(missed_endpoint,
        missed_endpoint_state, PlantTimestamp(6))
    _assert_command_application_error(
        () -> apply_claimed_plant_command!(missed_workspace,
            missed_endpoint, missed_endpoint_state,
            missed_application_state, missed_claim),
        :application,
        :missed_application_timestamp,
    )
    missed_disposition = _assert_application_disposition(
        missed_workspace, FailedCommand, :missed_application_timestamp)
    @test command_terminal_timestamp(missed_disposition) == PlantTimestamp(6)
    @test command_lateness(missed_disposition) == PlantDuration(1)
    @test effective_command(missed_application_state) == 0.25

    array_schema = _command_application_schema(
        Float32;
        id=:array_incremental_application_schema,
        endpoint=:array_incremental_application_command,
        dimensions=(3,),
        semantics=IncrementalCommand,
        value_policy=CommandValuePolicy(
            out_of_range=ClipInvalidCommand,
            range_stage=EnforceOnApplication,
        ),
    )
    array_endpoint = prepare_command_endpoint(array_schema;
        capacity=1, ordinal=7)
    array_endpoint_state = CommandEndpointState(array_endpoint)
    array_application_state = CommandApplicationState(array_endpoint,
        array_endpoint_state, Float32[0.75, 0, -0.75])
    array_workspace = CommandDispositionWorkspace(array_endpoint)
    caller_payload = Float32[0.5, 0.25, -0.5]
    admit_plant_command!(array_workspace, array_endpoint, array_endpoint_state,
        PlantCommand(array_schema, 1, PlantTimestamp(1), caller_payload),
        PlantTimestamp(1))
    caller_payload .= 9
    array_claim = claim_next_application_ready_command!(array_endpoint,
        array_endpoint_state, PlantTimestamp(1))
    array_disposition = apply_claimed_plant_command!(array_workspace,
        array_endpoint, array_endpoint_state, array_application_state,
        array_claim)
    @test command_disposition_reason(array_disposition) ==
        CommandDispositionReason(:applied_clipped)
    @test effective_command(array_application_state) ==
        Float32[1, 0.25, -1]

    absolute_array_schema = _command_application_schema(
        Float32;
        id=:absolute_array_application_schema,
        endpoint=:absolute_array_application_command,
        dimensions=(3,),
    )
    absolute_array_endpoint = prepare_command_endpoint(
        absolute_array_schema;
        capacity=1,
        ordinal=8,
    )
    absolute_array_endpoint_state = CommandEndpointState(
        absolute_array_endpoint)
    absolute_array_application_state = CommandApplicationState(
        absolute_array_endpoint, absolute_array_endpoint_state,
        zeros(Float32, 3))
    absolute_array_workspace = CommandDispositionWorkspace(
        absolute_array_endpoint)
    absolute_array_payload = Float32[0.5, 0.25, -0.5]
    admit_plant_command!(absolute_array_workspace, absolute_array_endpoint,
        absolute_array_endpoint_state,
        PlantCommand(absolute_array_schema, 1, PlantTimestamp(1),
            absolute_array_payload),
        PlantTimestamp(1))
    absolute_array_claim = claim_next_application_ready_command!(
        absolute_array_endpoint, absolute_array_endpoint_state,
        PlantTimestamp(1))
    apply_claimed_plant_command!(absolute_array_workspace,
        absolute_array_endpoint, absolute_array_endpoint_state,
        absolute_array_application_state, absolute_array_claim)
    @test effective_command(absolute_array_application_state) ==
        absolute_array_payload

    reject_array_schema = _command_application_schema(
        Float32;
        id=:reject_array_application_schema,
        endpoint=:reject_array_application_command,
        dimensions=(3,),
        semantics=IncrementalCommand,
        value_policy=CommandValuePolicy(
            out_of_range=RejectInvalidCommand,
            range_stage=EnforceOnApplication,
        ),
    )
    reject_array_endpoint = prepare_command_endpoint(reject_array_schema;
        capacity=1, ordinal=9)
    reject_array_endpoint_state = CommandEndpointState(reject_array_endpoint)
    reject_array_initial = Float32[0.75, 0, -0.75]
    reject_array_application_state = CommandApplicationState(
        reject_array_endpoint, reject_array_endpoint_state,
        reject_array_initial)
    reject_array_workspace = CommandDispositionWorkspace(reject_array_endpoint)
    admit_plant_command!(reject_array_workspace, reject_array_endpoint,
        reject_array_endpoint_state,
        PlantCommand(reject_array_schema, 1, PlantTimestamp(1),
            Float32[0.5, 0, -0.5]),
        PlantTimestamp(1))
    reject_array_claim = claim_next_application_ready_command!(
        reject_array_endpoint, reject_array_endpoint_state,
        PlantTimestamp(1))
    reject_array_disposition = apply_claimed_plant_command!(
        reject_array_workspace, reject_array_endpoint,
        reject_array_endpoint_state, reject_array_application_state,
        reject_array_claim)
    @test command_terminal_kind(reject_array_disposition) == RejectedCommand
    @test command_disposition_reason(reject_array_disposition) ==
        CommandDispositionReason(:out_of_range_rejected)
    @test effective_command(reject_array_application_state) ==
        reject_array_initial
end

@testset "Application storage-failure translation" begin
    schema = _command_application_schema(
        Float64;
        id=:application_storage_failure_schema,
        endpoint=:application_storage_failure_command,
    )
    endpoint = prepare_command_endpoint(schema; capacity=1, ordinal=1)
    endpoint_state = CommandEndpointState(endpoint)
    application_state = CommandApplicationState(
        endpoint, endpoint_state, 0.25)
    workspace = CommandDispositionWorkspace(endpoint)
    admit_plant_command!(workspace, endpoint, endpoint_state,
        PlantCommand(schema, 1, PlantTimestamp(1), 0.5),
        PlantTimestamp(1))
    claim = claim_next_application_ready_command!(endpoint, endpoint_state,
        PlantTimestamp(1))

    @test_throws InterruptException _handle_command_application_staging_error!(
        workspace, endpoint, endpoint_state, claim, InterruptException())
    @test active_command_count(endpoint_state) == 1
    _assert_command_application_error(
        () -> _handle_command_application_staging_error!(workspace,
            endpoint, endpoint_state, claim, ErrorException("copy failed")),
        :application,
        :application_storage_failure,
    )
    _assert_application_disposition(workspace, FailedCommand,
        :application_storage_failure)
    @test effective_command(application_state) == 0.25
    @test active_command_count(endpoint_state) == 0

    @test_throws InterruptException _handle_safe_command_staging_error(
        InterruptException())
    _assert_command_application_error(
        () -> _handle_safe_command_staging_error(
            ErrorException("copy failed")),
        :silence,
        :safe_storage_failure,
    )
end

@testset "Command-silence timing, equal-time order, and terminal failure" begin
    safe_schema = _command_application_schema(
        Float64;
        id=:safe_silence_schema,
        endpoint=:safe_silence_command,
        silence_policy=CommandSilencePolicy(
            ApplySafeCommand,
            AgeFromApplication;
            timeout=PlantDuration(10),
        ),
    )
    safe_endpoint = prepare_command_endpoint(safe_schema;
        capacity=2, ordinal=1)
    safe_endpoint_state = CommandEndpointState(safe_endpoint)
    safe_application_state = CommandApplicationState(safe_endpoint,
        safe_endpoint_state, 0.25; safe_command=-0.5)
    safe_workspace = CommandDispositionWorkspace(safe_endpoint)
    @test next_command_silence_timestamp(safe_endpoint, safe_endpoint_state,
        safe_application_state) == PlantTimestamp(10)

    admit_plant_command!(safe_workspace, safe_endpoint, safe_endpoint_state,
        PlantCommand(safe_schema, 1, PlantTimestamp(10), 0.75),
        PlantTimestamp(0))
    _assert_command_application_error(
        () -> apply_command_silence_transition!(safe_workspace, safe_endpoint,
            safe_endpoint_state, safe_application_state, PlantTimestamp(10)),
        :silence,
        :commands_due,
    )
    claim = claim_next_application_ready_command!(safe_endpoint,
        safe_endpoint_state, PlantTimestamp(10))
    apply_claimed_plant_command!(safe_workspace, safe_endpoint,
        safe_endpoint_state, safe_application_state, claim)
    clear_command_dispositions!(safe_workspace)
    @test effective_command(safe_application_state) == 0.75
    @test next_command_silence_timestamp(safe_endpoint, safe_endpoint_state,
        safe_application_state) == PlantTimestamp(20)

    transition = @inferred apply_command_silence_transition!(safe_workspace,
        safe_endpoint, safe_endpoint_state, safe_application_state,
        PlantTimestamp(20))
    @test command_endpoint_id(transition) == command_endpoint_id(safe_endpoint)
    @test command_silence_action(transition) == ApplySafeCommand
    @test command_silence_age_origin(transition) == AgeFromApplication
    @test command_silence_origin_timestamp(transition) == PlantTimestamp(10)
    @test command_silence_deadline(transition) == PlantTimestamp(20)
    @test command_silence_transition_timestamp(transition) ==
        PlantTimestamp(20)
    @test command_silence_age(transition) == PlantDuration(10)
    @test effective_command(safe_application_state) == -0.5
    @test command_disposition_count(safe_workspace) == 0
    @test next_command_silence_timestamp(safe_endpoint, safe_endpoint_state,
        safe_application_state) === nothing

    admit_plant_command!(safe_workspace, safe_endpoint, safe_endpoint_state,
        PlantCommand(safe_schema, 2, PlantTimestamp(30), 0.5),
        PlantTimestamp(21))
    @test next_command_silence_timestamp(safe_endpoint, safe_endpoint_state,
        safe_application_state) === nothing
    claim_two = claim_next_application_ready_command!(safe_endpoint,
        safe_endpoint_state, PlantTimestamp(30))
    apply_claimed_plant_command!(safe_workspace, safe_endpoint,
        safe_endpoint_state, safe_application_state, claim_two)
    clear_command_dispositions!(safe_workspace)
    @test next_command_silence_timestamp(safe_endpoint, safe_endpoint_state,
        safe_application_state) == PlantTimestamp(40)

    admission_schema = _command_application_schema(
        Float64;
        id=:admission_silence_schema,
        endpoint=:admission_silence_command,
        silence_policy=CommandSilencePolicy(
            ApplySafeCommand,
            AgeFromAdmission;
            timeout=PlantDuration(10),
        ),
    )
    admission_endpoint = prepare_command_endpoint(admission_schema;
        capacity=2, ordinal=2)
    admission_endpoint_state = CommandEndpointState(admission_endpoint)
    admission_application_state = CommandApplicationState(admission_endpoint,
        admission_endpoint_state, 0.0; safe_command=-0.25)
    admission_workspace = CommandDispositionWorkspace(admission_endpoint)
    @test next_command_silence_timestamp(admission_endpoint,
        admission_endpoint_state, admission_application_state) ==
        PlantTimestamp(10)
    admit_plant_command!(admission_workspace, admission_endpoint,
        admission_endpoint_state,
        PlantCommand(admission_schema, 1, PlantTimestamp(20), 0.5),
        PlantTimestamp(5))
    @test next_command_silence_timestamp(admission_endpoint,
        admission_endpoint_state, admission_application_state) ==
        PlantTimestamp(15)
    admission_transition = apply_command_silence_transition!(
        admission_workspace, admission_endpoint, admission_endpoint_state,
        admission_application_state, PlantTimestamp(15))
    @test command_silence_origin_timestamp(admission_transition) ==
        PlantTimestamp(5)
    @test effective_command(admission_application_state) == -0.25
    admission_claim = claim_next_application_ready_command!(
        admission_endpoint, admission_endpoint_state, PlantTimestamp(20))
    apply_claimed_plant_command!(admission_workspace, admission_endpoint,
        admission_endpoint_state, admission_application_state,
        admission_claim)
    clear_command_dispositions!(admission_workspace)
    @test next_command_silence_timestamp(admission_endpoint,
        admission_endpoint_state, admission_application_state) === nothing
    admit_plant_command!(admission_workspace, admission_endpoint,
        admission_endpoint_state,
        PlantCommand(admission_schema, 2, PlantTimestamp(40), 0.75),
        PlantTimestamp(25))
    @test next_command_silence_timestamp(admission_endpoint,
        admission_endpoint_state, admission_application_state) ==
        PlantTimestamp(35)

    fail_schema = _command_application_schema(
        Float64;
        id=:fail_silence_schema,
        endpoint=:fail_silence_command,
        silence_policy=CommandSilencePolicy(
            FailOnCommandSilence,
            AgeFromApplication;
            timeout=PlantDuration(10),
        ),
    )
    fail_endpoint = prepare_command_endpoint(fail_schema;
        capacity=2, ordinal=3)
    fail_endpoint_state = CommandEndpointState(fail_endpoint)
    fail_application_state = CommandApplicationState(
        fail_endpoint, fail_endpoint_state, 0.1)
    fail_workspace = CommandDispositionWorkspace(fail_endpoint)
    for (sequence, timestamp) in ((1, 30), (2, 20))
        admit_plant_command!(fail_workspace, fail_endpoint, fail_endpoint_state,
            PlantCommand(fail_schema, sequence, PlantTimestamp(timestamp), 0.0),
            PlantTimestamp(0))
    end
    _assert_command_application_error(
        () -> apply_command_silence_transition!(fail_workspace, fail_endpoint,
            fail_endpoint_state, fail_application_state, PlantTimestamp(9)),
        :silence,
        :unexpected_silence_timestamp,
    )
    fail_transition = apply_command_silence_transition!(fail_workspace,
        fail_endpoint, fail_endpoint_state, fail_application_state,
        PlantTimestamp(10))
    @test command_silence_action(fail_transition) == FailOnCommandSilence
    @test command_endpoint_failed(fail_endpoint_state)
    @test effective_command(fail_application_state) == 0.1
    @test pending_command_count(fail_endpoint_state) == 0
    @test active_command_count(fail_endpoint_state) == 0
    @test command_disposition_count(fail_workspace) == 2
    @test all(command_terminal_kind(command_disposition(fail_workspace,
        index)) == FailedCommand for index in 1:2)
    @test all(command_disposition_reason(command_disposition(fail_workspace,
        index)) == CommandDispositionReason(:command_silence)
        for index in 1:2)
    @test next_command_silence_timestamp(fail_endpoint, fail_endpoint_state,
        fail_application_state) === nothing
    _assert_command_application_error(
        () -> admit_plant_command!(fail_workspace, fail_endpoint,
            fail_endpoint_state,
            PlantCommand(fail_schema, 3, PlantTimestamp(40), 0.0),
            PlantTimestamp(10)),
        :endpoint,
        :endpoint_failed,
    )
end

@testset "Effective-command inference and warmed allocation" begin
    scalar_schema = _command_application_schema(
        Float64;
        id=:scalar_application_hot_schema,
        endpoint=:scalar_application_hot_command,
    )
    scalar_endpoint = prepare_command_endpoint(scalar_schema;
        capacity=1, ordinal=1)
    scalar_endpoint_state = CommandEndpointState(scalar_endpoint)
    scalar_application_state = CommandApplicationState(
        scalar_endpoint, scalar_endpoint_state, 0.0)
    scalar_workspace = CommandDispositionWorkspace(scalar_endpoint)
    @test @inferred(_apply_scalar_command_cycle!(scalar_workspace,
        scalar_endpoint, scalar_endpoint_state, scalar_application_state,
        1, PlantTimestamp(1), 0.25)) === nothing

    array_schema = _command_application_schema(
        Float32;
        id=:array_application_hot_schema,
        endpoint=:array_application_hot_command,
        dimensions=(3,),
        semantics=IncrementalCommand,
        value_policy=CommandValuePolicy(
            out_of_range=ClipInvalidCommand,
            range_stage=EnforceOnApplication,
        ),
    )
    array_endpoint = prepare_command_endpoint(array_schema;
        capacity=1, ordinal=2)
    array_endpoint_state = CommandEndpointState(array_endpoint)
    array_application_state = CommandApplicationState(array_endpoint,
        array_endpoint_state, zeros(Float32, 3))
    array_workspace = CommandDispositionWorkspace(array_endpoint)
    array_payload = Float32[0.1, 0.2, 0.3]
    @test @inferred(_apply_array_command_cycle!(array_workspace,
        array_endpoint, array_endpoint_state, array_application_state,
        1, PlantTimestamp(1), array_payload)) === nothing

    safe_schema = _command_application_schema(
        Float64;
        id=:safe_application_hot_schema,
        endpoint=:safe_application_hot_command,
        silence_policy=CommandSilencePolicy(
            ApplySafeCommand,
            AgeFromApplication;
            timeout=PlantDuration(10),
        ),
    )
    safe_endpoint = prepare_command_endpoint(safe_schema;
        capacity=1, ordinal=3)
    safe_endpoint_state = CommandEndpointState(safe_endpoint)
    safe_application_state = CommandApplicationState(
        safe_endpoint, safe_endpoint_state, 0.0; safe_command=-0.25)
    safe_workspace = CommandDispositionWorkspace(safe_endpoint)
    @test @inferred(_safe_command_cycle!(safe_workspace, safe_endpoint,
        safe_endpoint_state, safe_application_state, 1,
        PlantTimestamp(1))) === nothing

    if coverage_instrumented()
        @test_skip "command-application allocation assertions disabled under coverage"
    else
        @test @allocated(_apply_scalar_command_cycle!(scalar_workspace,
            scalar_endpoint, scalar_endpoint_state, scalar_application_state,
            2, PlantTimestamp(2), 0.5)) == 0
        @test @allocated(_apply_array_command_cycle!(array_workspace,
            array_endpoint, array_endpoint_state, array_application_state,
            2, PlantTimestamp(2), array_payload)) == 0
        @test @allocated(_safe_command_cycle!(safe_workspace, safe_endpoint,
            safe_endpoint_state, safe_application_state, 2,
            PlantTimestamp(21))) == 0
    end
end

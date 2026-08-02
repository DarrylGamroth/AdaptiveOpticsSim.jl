#
# Authority-owned effective-command fanout
#
# This private serial coordinator is the command-side correctness oracle for a
# prepared mixed plant. The sole authority applies command semantics once, then
# publishes the already-effective value to every target-local replica. It owns
# no scheduler, task, queue, pacing policy, or RTC transport.
#

mutable struct _PreparedCommandFanoutBinding end

struct _PreparedCommandFanoutLane{B,R<:Tuple}
    authority_endpoint::B
    endpoint_slot::UInt32
    optic_slot::UInt32
    routes::R
end

struct _PreparedCommandFanout{B,A,L<:Tuple}
    binding::B
    authority::A
    lanes::L
end

@enum _CommandFanoutPhase::UInt8 begin
    _CommandFanoutIdle = 0x01
    _CommandFanoutApplying = 0x02
    _CommandFanoutCommitting = 0x03
    _CommandFanoutFailed = 0x04
    _CommandFanoutUncertain = 0x05
end

struct _CommandFanoutLaneState{E,A,R<:Tuple}
    endpoint::E
    application::A
    routes::R
end

mutable struct _CommandFanoutState{B,A,L<:Tuple}
    binding::B
    authority::A
    lanes::L
    phase::_CommandFanoutPhase
end

mutable struct _CommandFanoutLaneWorkspace{D,F}
    disposition::D
    claim::Base.RefValue{PlantCommandApplicationClaim}
    staged::Base.RefValue{_StagedCommandApplication}
    finalization::Base.RefValue{F}
    silence::Base.RefValue{PlantCommandSilenceTransition}
    transaction_admission::Base.RefValue{_CommandTransactionAdmissionPlan}
    proposed_sequence::UInt64
    selected::Bool
    has_claim::Bool
    has_staged::Bool
    has_finalization::Bool
    has_transaction_admission::Bool
end

struct _CommandFanoutWorkspace{B,A,L<:Tuple}
    binding::B
    authority::A
    lanes::L
end

@noinline function _command_fanout_preparation_error(
    reason::Symbol,
    message::AbstractString,
)
    throw(PlantPreparationError(:command_fanout, reason, String(message)))
end

@noinline function _command_fanout_error(
    reason::Symbol,
    message::AbstractString,
)
    throw(PlantCommandError(:command_fanout, reason, String(message)))
end

function _command_fanout_endpoint_routes(
    routes::_PreparedEffectiveCommandPublicationRoutes,
    authority_endpoint::_PreparedPlantCommandEndpoint,
)
    selected = _PreparedEffectiveCommandPublicationRoute[]
    for route in routes
        route.authority_endpoint.endpoint.binding ===
            authority_endpoint.endpoint.binding || continue
        push!(selected, route)
    end
    isempty(selected) && _command_fanout_preparation_error(
        :missing_replica,
        "command endpoint $(command_endpoint_id(authority_endpoint)) has no " *
        "effective-command publication route",
    )
    return Tuple(selected)
end

function _prepare_command_fanout(partitions::PreparedPlantPartitions)
    authority = prepared_command_authority(partitions)
    routes = _prepare_effective_command_publication_routes(partitions)
    # Cold heterogeneous builder; `Tuple` below seals each concrete lane type
    # for allocation-free recursive execution.
    lanes = Any[]
    sizehint!(lanes, length(authority.endpoints))
    @inbounds for endpoint_slot in eachindex(authority.endpoints)
        authority_endpoint = authority.endpoints[endpoint_slot]
        push!(lanes, _PreparedCommandFanoutLane(
            authority_endpoint,
            UInt32(endpoint_slot),
            authority_endpoint.optic_slot,
            _command_fanout_endpoint_routes(routes, authority_endpoint),
        ))
    end
    return _PreparedCommandFanout(
        _PreparedCommandFanoutBinding(),
        authority,
        Tuple(lanes),
    )
end

function _prepare_command_fanout_lane_state(
    lane::_PreparedCommandFanoutLane,
    authority_state::CommandAuthorityState,
)
    endpoint_slot = Int(lane.endpoint_slot)
    endpoint_state = @inbounds authority_state.endpoint_states[endpoint_slot]
    application_state =
        @inbounds authority_state.application_states[endpoint_slot]
    route_states = map(
        _prepare_effective_command_publication_route_state,
        lane.routes,
    )
    return _CommandFanoutLaneState(
        endpoint_state, application_state, route_states)
end

function _prepare_command_fanout_state(prepared::_PreparedCommandFanout)
    authority_state = CommandAuthorityState(prepared.authority)
    lane_states = map(
        lane -> _prepare_command_fanout_lane_state(lane, authority_state),
        prepared.lanes,
    )
    return _CommandFanoutState(
        prepared.binding,
        authority_state,
        lane_states,
        _CommandFanoutIdle,
    )
end

function _prepare_command_fanout_lane_workspace(
    lane_state::_CommandFanoutLaneState,
    disposition::CommandDispositionWorkspace,
)
    F = _PrevalidatedCommandApplicationDisposition{
        typeof(lane_state.endpoint),
    }
    return _CommandFanoutLaneWorkspace(
        disposition,
        Ref{PlantCommandApplicationClaim}(),
        Ref{_StagedCommandApplication}(),
        Ref{F}(),
        Ref{PlantCommandSilenceTransition}(),
        Ref{_CommandTransactionAdmissionPlan}(),
        zero(UInt64),
        false,
        false,
        false,
        false,
        false,
    )
end

@inline function _select_command_fanout_endpoint!(
    ::Tuple{},
    ::Tuple{},
    endpoint::CommandEndpointID,
)
    _command_fanout_error(
        :unknown_endpoint,
        "command authority has no endpoint $endpoint",
    )
end

@inline function _select_command_fanout_endpoint!(
    lanes::Tuple,
    workspaces::Tuple,
    endpoint::CommandEndpointID,
)
    lane = first(lanes)
    if command_endpoint_id(lane.authority_endpoint) == endpoint
        first(workspaces).selected = true
        return nothing
    end
    return _select_command_fanout_endpoint!(
        Base.tail(lanes), Base.tail(workspaces), endpoint)
end

@inline _command_fanout_transaction_contains_endpoint(
    ::Tuple{}, ::CommandEndpointID) = false

@inline function _command_fanout_transaction_contains_endpoint(
    commands::Tuple,
    endpoint::CommandEndpointID,
)
    command_endpoint_id(first(commands)) == endpoint && return true
    return _command_fanout_transaction_contains_endpoint(
        Base.tail(commands), endpoint)
end

@inline function _select_canonical_command_fanout_transaction_lanes!(
    ::Tuple{},
    ::Tuple{},
    ::Tuple,
    ::Tuple,
    ::Tuple,
    count::Int,
)
    return count
end

@inline function _select_canonical_command_fanout_transaction_lanes!(
    lanes::Tuple,
    workspaces::Tuple,
    all_lanes::Tuple,
    all_workspaces::Tuple,
    commands::Tuple,
    count::Int,
)
    lane = first(lanes)
    workspace = first(workspaces)
    if _command_fanout_transaction_contains_endpoint(
        commands, command_endpoint_id(lane.authority_endpoint))
        _selected_command_fanout_has_optic(
            all_lanes, all_workspaces, lane.optic_slot) &&
            _command_admission_error(
                :transaction,
                :duplicate_physical_optic,
                "atomic transaction contains multiple endpoints for one " *
                "physical controllable optic",
            )
        workspace.selected = true
        count += 1
    end
    return _select_canonical_command_fanout_transaction_lanes!(
        Base.tail(lanes),
        Base.tail(workspaces),
        all_lanes,
        all_workspaces,
        commands,
        count,
    )
end

@inline function _select_canonical_command_fanout_transaction_lanes!(
    lanes::Tuple,
    workspaces::Tuple,
    commands::Tuple,
)
    count = _select_canonical_command_fanout_transaction_lanes!(
        lanes, workspaces, lanes, workspaces, commands, 0)
    count == length(commands) || _command_admission_error(
        :transaction,
        :unknown_endpoint,
        "atomic transaction references an endpoint outside the command " *
        "authority",
    )
    return count
end

@inline function _preflight_canonical_command_fanout_lane!(
    ::Tuple{},
    lane::_PreparedCommandFanoutLane,
    ::_CommandFanoutLaneState,
    ::_CommandFanoutLaneWorkspace,
    ::PlantTimestamp,
)
    _command_admission_error(
        :transaction,
        :transaction_invariant,
        "selected canonical command-fanout lane has no transaction member " *
        "($(command_endpoint_id(lane.authority_endpoint)))",
    )
end

@inline function _preflight_canonical_command_fanout_lane!(
    commands::Tuple,
    lane::_PreparedCommandFanoutLane,
    lane_state::_CommandFanoutLaneState,
    workspace::_CommandFanoutLaneWorkspace,
    timestamp::PlantTimestamp,
)
    command = first(commands)
    if command_endpoint_id(command) ==
            command_endpoint_id(lane.authority_endpoint)
        plan = _preflight_command_transaction_member!(
            lane.authority_endpoint.endpoint,
            lane_state.endpoint,
            workspace.disposition,
            command,
            timestamp,
        )
        workspace.transaction_admission[] = plan
        workspace.has_transaction_admission = true
        return plan.scheduled_timestamp
    end
    return _preflight_canonical_command_fanout_lane!(
        Base.tail(commands), lane, lane_state, workspace, timestamp)
end

@inline function _preflight_canonical_command_fanout_transaction!(
    ::Tuple{},
    ::Tuple{},
    ::Tuple{},
    ::Tuple,
    ::PlantTimestamp,
    scheduled::PlantTimestamp,
    has_scheduled::Bool,
)
    has_scheduled || _command_admission_error(
        :transaction,
        :transaction_invariant,
        "atomic transaction selected no canonical command-fanout lanes",
    )
    return scheduled
end

@inline function _preflight_canonical_command_fanout_transaction!(
    lanes::Tuple,
    lane_states::Tuple,
    workspaces::Tuple,
    commands::Tuple,
    timestamp::PlantTimestamp,
    scheduled::PlantTimestamp,
    has_scheduled::Bool,
)
    lane = first(lanes)
    workspace = first(workspaces)
    if workspace.selected
        member_scheduled = _preflight_canonical_command_fanout_lane!(
            commands, lane, first(lane_states), workspace, timestamp)
        has_scheduled && member_scheduled != scheduled &&
            _command_admission_error(
                :transaction,
                :scheduled_timestamp_mismatch,
                "atomic transaction members resolved to different " *
                "scheduled plant timestamps",
            )
        scheduled = member_scheduled
        has_scheduled = true
    end
    return _preflight_canonical_command_fanout_transaction!(
        Base.tail(lanes),
        Base.tail(lane_states),
        Base.tail(workspaces),
        commands,
        timestamp,
        scheduled,
        has_scheduled,
    )
end

@inline function _preflight_canonical_command_fanout_transaction!(
    lanes::Tuple,
    lane_states::Tuple,
    workspaces::Tuple,
    commands::Tuple,
    timestamp::PlantTimestamp,
)
    return _preflight_canonical_command_fanout_transaction!(
        lanes,
        lane_states,
        workspaces,
        commands,
        timestamp,
        zero(PlantTimestamp),
        false,
    )
end

@inline function _commit_canonical_command_fanout_lane!(
    ::Tuple{},
    lane::_PreparedCommandFanoutLane,
    ::_CommandFanoutLaneState,
    ::_CommandFanoutLaneWorkspace,
    ::PlantTimestamp,
    ::UInt64,
    ::UInt32,
)
    _command_admission_error(
        :transaction,
        :transaction_invariant,
        "selected canonical command-fanout lane has no commit member " *
        "($(command_endpoint_id(lane.authority_endpoint)))",
    )
end

@inline function _commit_canonical_command_fanout_lane!(
    commands::Tuple,
    lane::_PreparedCommandFanoutLane,
    lane_state::_CommandFanoutLaneState,
    workspace::_CommandFanoutLaneWorkspace,
    timestamp::PlantTimestamp,
    transaction_id::UInt64,
    member_count::UInt32,
)
    command = first(commands)
    if command_endpoint_id(command) ==
            command_endpoint_id(lane.authority_endpoint)
        workspace.has_transaction_admission || _command_admission_error(
            :transaction,
            :transaction_invariant,
            "atomic-transaction member has no prevalidated admission",
        )
        _commit_command_transaction_member!(
            lane.authority_endpoint.endpoint,
            lane_state.endpoint,
            command,
            timestamp,
            transaction_id,
            member_count,
            workspace.transaction_admission[],
        )
        workspace.has_transaction_admission = false
        return nothing
    end
    return _commit_canonical_command_fanout_lane!(
        Base.tail(commands),
        lane,
        lane_state,
        workspace,
        timestamp,
        transaction_id,
        member_count,
    )
end

@inline function _commit_canonical_command_fanout_transaction!(
    ::Tuple{},
    ::Tuple{},
    ::Tuple{},
    ::Tuple,
    ::PlantTimestamp,
    ::UInt64,
    ::UInt32,
)
    return nothing
end

@inline function _commit_canonical_command_fanout_transaction!(
    lanes::Tuple,
    lane_states::Tuple,
    workspaces::Tuple,
    commands::Tuple,
    timestamp::PlantTimestamp,
    transaction_id::UInt64,
    member_count::UInt32,
)
    lane = first(lanes)
    workspace = first(workspaces)
    workspace.selected && _commit_canonical_command_fanout_lane!(
        commands,
        lane,
        first(lane_states),
        workspace,
        timestamp,
        transaction_id,
        member_count,
    )
    return _commit_canonical_command_fanout_transaction!(
        Base.tail(lanes),
        Base.tail(lane_states),
        Base.tail(workspaces),
        commands,
        timestamp,
        transaction_id,
        member_count,
    )
end

@inline function _terminate_canonical_command_fanout_lane!(
    ::Tuple{},
    lane::_PreparedCommandFanoutLane,
    ::_CommandFanoutLaneState,
    ::_CommandFanoutLaneWorkspace,
    ::PlantTimestamp,
    ::CommandTerminalKind,
    ::CommandDispositionReason,
)
    _command_admission_error(
        :transaction,
        :transaction_invariant,
        "selected canonical command-fanout lane has no terminal member " *
        "($(command_endpoint_id(lane.authority_endpoint)))",
    )
end

@inline function _terminate_canonical_command_fanout_lane!(
    commands::Tuple,
    lane::_PreparedCommandFanoutLane,
    lane_state::_CommandFanoutLaneState,
    workspace::_CommandFanoutLaneWorkspace,
    timestamp::PlantTimestamp,
    kind::CommandTerminalKind,
    reason::CommandDispositionReason,
)
    command = first(commands)
    if command_endpoint_id(command) ==
            command_endpoint_id(lane.authority_endpoint)
        _finish_terminal_admission!(
            workspace.disposition,
            lane.authority_endpoint.endpoint,
            lane_state.endpoint,
            command,
            timestamp,
            _next_command_presentation(lane_state.endpoint),
            kind,
            reason,
            nothing,
        )
        return nothing
    end
    return _terminate_canonical_command_fanout_lane!(
        Base.tail(commands),
        lane,
        lane_state,
        workspace,
        timestamp,
        kind,
        reason,
    )
end

@inline function _terminate_canonical_command_fanout_transaction!(
    ::Tuple{},
    ::Tuple{},
    ::Tuple{},
    ::Tuple,
    ::PlantTimestamp,
    ::CommandTerminalKind,
    ::CommandDispositionReason,
)
    return nothing
end

@inline function _terminate_canonical_command_fanout_transaction!(
    lanes::Tuple,
    lane_states::Tuple,
    workspaces::Tuple,
    commands::Tuple,
    timestamp::PlantTimestamp,
    kind::CommandTerminalKind,
    reason::CommandDispositionReason,
)
    lane = first(lanes)
    workspace = first(workspaces)
    workspace.selected && _terminate_canonical_command_fanout_lane!(
        commands,
        lane,
        first(lane_states),
        workspace,
        timestamp,
        kind,
        reason,
    )
    return _terminate_canonical_command_fanout_transaction!(
        Base.tail(lanes),
        Base.tail(lane_states),
        Base.tail(workspaces),
        commands,
        timestamp,
        kind,
        reason,
    )
end

function _terminate_command_fanout_transaction_admission!(
    prepared::_PreparedCommandFanout,
    state::_CommandFanoutState,
    workspace::_CommandFanoutWorkspace,
    transaction::PlantCommandTransaction,
    timestamp::PlantTimestamp,
    kind::CommandTerminalKind,
    reason::CommandDispositionReason,
)
    commands = transaction.commands
    _terminate_canonical_command_fanout_transaction!(
        prepared.lanes,
        state.lanes,
        workspace.lanes,
        commands,
        timestamp,
        kind,
        reason,
    )
    _clear_command_fanout_selection!(workspace.lanes)
    return PlantCommandTransactionAdmission(
        zero(UInt64),
        CommandTerminatedOnAdmission,
        UInt32(length(commands)),
        nothing,
    )
end

@noinline function _handle_command_fanout_transaction_preflight_failure!(
    ::_PreparedCommandFanout,
    ::_CommandFanoutState,
    workspace::_CommandFanoutWorkspace,
    ::PlantCommandTransaction,
    ::PlantTimestamp,
    error::InterruptException,
)
    _clear_command_fanout_selection!(workspace.lanes)
    throw(error)
end

@noinline function _handle_command_fanout_transaction_preflight_failure!(
    prepared::_PreparedCommandFanout,
    state::_CommandFanoutState,
    workspace::_CommandFanoutWorkspace,
    transaction::PlantCommandTransaction,
    timestamp::PlantTimestamp,
    failure::_CommandTransactionPolicyFailure,
)
    error = failure.error
    reason = CommandDispositionReason(Symbol(
        :atomic_transaction_aborted_, error.reason))
    _terminate_command_fanout_transaction_admission!(
        prepared,
        state,
        workspace,
        transaction,
        timestamp,
        FailedCommand,
        reason,
    )
    drained = _fail_all_command_fanout_endpoints!(
        prepared.lanes,
        state.lanes,
        workspace.lanes,
        timestamp,
        reason,
    )
    state.authority.failed = true
    state.phase = drained ? _CommandFanoutFailed : _CommandFanoutUncertain
    throw(error)
end

@noinline function _handle_command_fanout_transaction_preflight_failure!(
    prepared::_PreparedCommandFanout,
    state::_CommandFanoutState,
    workspace::_CommandFanoutWorkspace,
    transaction::PlantCommandTransaction,
    timestamp::PlantTimestamp,
    error::PlantCommandError,
)
    reason = CommandDispositionReason(Symbol(
        :atomic_transaction_aborted_, error.reason))
    kind = _validation_terminal_kind(error)
    admission = _terminate_command_fanout_transaction_admission!(
        prepared,
        state,
        workspace,
        transaction,
        timestamp,
        kind,
        reason,
    )
    if kind == FailedCommand
        drained = _fail_all_command_fanout_endpoints!(
            prepared.lanes,
            state.lanes,
            workspace.lanes,
            timestamp,
            reason,
        )
        state.authority.failed = true
        state.phase = drained ?
            _CommandFanoutFailed : _CommandFanoutUncertain
        throw(error)
    end
    return admission
end

@noinline function _handle_command_fanout_transaction_preflight_failure!(
    prepared::_PreparedCommandFanout,
    state::_CommandFanoutState,
    workspace::_CommandFanoutWorkspace,
    transaction::PlantCommandTransaction,
    timestamp::PlantTimestamp,
    error,
)
    reason = CommandDispositionReason(:transaction_preflight_failure)
    _terminate_command_fanout_transaction_admission!(
        prepared,
        state,
        workspace,
        transaction,
        timestamp,
        FailedCommand,
        reason,
    )
    _fail_all_command_fanout_endpoints!(
        prepared.lanes,
        state.lanes,
        workspace.lanes,
        timestamp,
        reason,
    )
    state.authority.failed = true
    state.phase = _CommandFanoutUncertain
    _command_admission_error(
        :transaction,
        reason.name,
        "atomic transaction preflight failed unexpectedly " *
        "($(typeof(error)))",
    )
end

function _prepare_command_fanout_lane_workspaces(
    ::Tuple{},
    ::Tuple{},
    ::CommandAuthorityWorkspace,
)
    return ()
end

function _prepare_command_fanout_lane_workspaces(
    lanes::Tuple,
    lane_states::Tuple,
    authority_workspace::CommandAuthorityWorkspace,
)
    lane = first(lanes)
    disposition = @inbounds authority_workspace.dispositions[
        Int(lane.endpoint_slot)]
    workspace = _prepare_command_fanout_lane_workspace(
        first(lane_states), disposition)
    return (workspace, _prepare_command_fanout_lane_workspaces(
        Base.tail(lanes),
        Base.tail(lane_states),
        authority_workspace,
    )...)
end

function _prepare_command_fanout_authority_workspace(
    authority::PreparedCommandAuthority,
)
    dispositions = Memory{CommandDispositionWorkspace}(
        undef, length(authority.endpoints))
    @inbounds for index in eachindex(authority.endpoints)
        # A structural transaction rejection can terminalize the new member
        # before draining a full endpoint calendar.  Reserve that one extra
        # credit so fail-stop never strands an accepted presentation.
        dispositions[index] = _command_disposition_workspace(
            authority.endpoints[index].endpoint, 1)
    end
    return CommandAuthorityWorkspace(
        authority.binding,
        _partition_registry(dispositions, CommandDispositionWorkspace),
    )
end

function _prepare_command_fanout_workspace(
    prepared::_PreparedCommandFanout,
    state::_CommandFanoutState,
)
    prepared.binding === state.binding ||
        _command_fanout_preparation_error(
            :foreign_state,
            "command-fanout state belongs to another prepared fanout",
        )
    authority_workspace =
        _prepare_command_fanout_authority_workspace(prepared.authority)
    lanes = _prepare_command_fanout_lane_workspaces(
        prepared.lanes, state.lanes, authority_workspace)
    return _CommandFanoutWorkspace(
        prepared.binding, authority_workspace, lanes)
end

@inline function _require_command_fanout_binding(
    prepared::_PreparedCommandFanout,
    state::_CommandFanoutState,
    workspace::_CommandFanoutWorkspace,
)
    prepared.binding === state.binding === workspace.binding ||
        _command_fanout_error(
            :foreign_owner,
            "command-fanout plan, state, and workspace belong to different owners",
        )
    _require_command_authority_binding(prepared.authority, state.authority)
    _require_command_authority_binding(
        prepared.authority, workspace.authority)
    state.phase == _CommandFanoutIdle || _command_fanout_error(
        :fanout_busy,
        "command fanout is not idle ($(state.phase))",
    )
    command_authority_failed(state.authority) && _command_fanout_error(
        :authority_failed,
        "command authority is fail-stop",
    )
    return nothing
end

@inline _require_command_fanout_workspace_idle(::Tuple{}) = nothing

@inline function _require_command_fanout_workspace_idle(workspaces::Tuple)
    workspace = first(workspaces)
    _require_empty_command_dispositions(workspace.disposition)
    !(workspace.selected || workspace.has_claim || workspace.has_staged ||
        workspace.has_finalization || workspace.has_transaction_admission) ||
        _command_fanout_error(
        :workspace_busy,
        "command-fanout workspace retains an active transaction",
    )
    return _require_command_fanout_workspace_idle(Base.tail(workspaces))
end

@inline _clear_command_fanout_selection!(::Tuple{}) = nothing

@inline function _clear_command_fanout_selection!(workspaces::Tuple)
    workspace = first(workspaces)
    workspace.selected = false
    workspace.has_claim = false
    workspace.has_staged = false
    workspace.has_finalization = false
    workspace.has_transaction_admission = false
    workspace.proposed_sequence = zero(UInt64)
    return _clear_command_fanout_selection!(Base.tail(workspaces))
end

@inline function _command_fanout_due_metadata(
    ::Tuple{},
    ::Tuple{},
    endpoint::CommandEndpointID,
    ::PlantTimestamp,
)
    _command_fanout_error(
        :unknown_endpoint,
        "command authority has no endpoint $endpoint",
    )
end

@inline function _command_fanout_due_metadata(
    lanes::Tuple,
    lane_states::Tuple,
    endpoint::CommandEndpointID,
    timestamp::PlantTimestamp,
)
    lane = first(lanes)
    lane_state = first(lane_states)
    if command_endpoint_id(lane.authority_endpoint) == endpoint
        key = next_command_order_key(
            lane.authority_endpoint.endpoint, lane_state.endpoint)
        key === nothing && _command_fanout_error(
            :command_not_due,
            "command endpoint $endpoint has no pending command",
        )
        command_scheduled_timestamp(key) == timestamp ||
            _command_fanout_error(
                :command_not_due,
                "next command for $endpoint is scheduled at " *
                "$(command_scheduled_timestamp(key)); got $timestamp",
            )
        transaction, member_count = _next_command_transaction_metadata(
            lane.authority_endpoint.endpoint, lane_state.endpoint)
        return transaction, member_count
    end
    return _command_fanout_due_metadata(
        Base.tail(lanes),
        Base.tail(lane_states),
        endpoint,
        timestamp,
    )
end

@inline _selected_command_fanout_has_optic(
    ::Tuple{}, ::Tuple{}, ::UInt32) = false

@inline function _selected_command_fanout_has_optic(
    lanes::Tuple,
    workspaces::Tuple,
    optic_slot::UInt32,
)
    lane = first(lanes)
    workspace = first(workspaces)
    workspace.selected && lane.optic_slot == optic_slot && return true
    return _selected_command_fanout_has_optic(
        Base.tail(lanes), Base.tail(workspaces), optic_slot)
end

@inline _require_distinct_selected_command_fanout_optics(
    ::Tuple{}, ::Tuple{}) = nothing

@inline function _require_distinct_selected_command_fanout_optics(
    lanes::Tuple,
    workspaces::Tuple,
)
    lane = first(lanes)
    workspace = first(workspaces)
    if workspace.selected
        _selected_command_fanout_has_optic(
            Base.tail(lanes),
            Base.tail(workspaces),
            lane.optic_slot,
        ) && _command_fanout_error(
            :duplicate_physical_optic,
            "atomic transaction contains multiple endpoints for one " *
            "physical controllable optic",
        )
    end
    return _require_distinct_selected_command_fanout_optics(
        Base.tail(lanes), Base.tail(workspaces))
end

@inline function _select_command_fanout_members!(
    ::Tuple{},
    ::Tuple{},
    ::Tuple{},
    ::CommandEndpointID,
    ::UInt64,
    ::UInt32,
    ::PlantTimestamp,
)
    return 0
end

@inline function _select_command_fanout_members!(
    lanes::Tuple,
    lane_states::Tuple,
    workspaces::Tuple,
    endpoint::CommandEndpointID,
    transaction::UInt64,
    expected_count::UInt32,
    timestamp::PlantTimestamp,
)
    lane = first(lanes)
    lane_state = first(lane_states)
    workspace = first(workspaces)
    lane_endpoint = lane.authority_endpoint.endpoint
    member_transaction, member_count =
        _next_command_transaction_metadata(lane_endpoint, lane_state.endpoint)
    selected = if iszero(transaction)
        command_endpoint_id(lane.authority_endpoint) == endpoint
    else
        member_transaction == transaction
    end
    if selected
        if !iszero(transaction)
            member_count == expected_count || _command_fanout_error(
                :transaction_invariant,
                "atomic transaction members disagree about member count",
            )
        end
        key = next_command_order_key(lane_endpoint, lane_state.endpoint)
        key !== nothing && command_scheduled_timestamp(key) == timestamp ||
            _command_fanout_error(
                :transaction_invariant,
                "atomic transaction member is not ready at its common timestamp",
            )
        workspace.selected = true
    end
    count = _select_command_fanout_members!(
        Base.tail(lanes),
        Base.tail(lane_states),
        Base.tail(workspaces),
        endpoint,
        transaction,
        expected_count,
        timestamp,
    )
    return count + Int(selected)
end

@inline _require_command_fanout_route_idle(::Tuple{}, ::Tuple{}) = nothing

@inline function _require_command_fanout_route_idle(
    routes::Tuple,
    states::Tuple,
)
    route = first(routes)
    state = first(states)
    _require_effective_command_route_state(route, state)
    state.phase == _EffectiveCommandRouteIdle || _command_fanout_error(
        :route_busy,
        "effective-command publication route is not idle",
    )
    _require_effective_command_route_binding(route)
    _require_effective_command_route_destination_available(route)
    return _require_command_fanout_route_idle(
        Base.tail(routes), Base.tail(states))
end

@inline function _require_selected_command_fanout_ready(
    ::Tuple{}, ::Tuple{}, ::Tuple{}, ::Memory{UInt64})
    return nothing
end

@inline function _require_selected_command_fanout_ready(
    lanes::Tuple,
    lane_states::Tuple,
    workspaces::Tuple,
    sequences::Memory{UInt64},
)
    lane = first(lanes)
    lane_state = first(lane_states)
    workspace = first(workspaces)
    if workspace.selected
        _require_empty_command_dispositions(workspace.disposition)
        current_sequence = @inbounds sequences[Int(lane.endpoint_slot)]
        current_sequence != typemax(UInt64) || _command_fanout_error(
            :publication_sequence_overflow,
            "effective-command publication sequence exceeds UInt64",
        )
        workspace.proposed_sequence = current_sequence + one(UInt64)
        _require_command_fanout_route_idle(lane.routes, lane_state.routes)
    end
    return _require_selected_command_fanout_ready(
        Base.tail(lanes),
        Base.tail(lane_states),
        Base.tail(workspaces),
        sequences,
    )
end

@inline function _prevalidate_selected_command_silence!(
    ::Tuple{},
    ::Tuple{},
    ::Tuple{},
    ::Memory{UInt64},
    ::PlantTimestamp,
)
    return nothing
end

@inline function _prevalidate_selected_command_silence!(
    lanes::Tuple,
    lane_states::Tuple,
    workspaces::Tuple,
    sequences::Memory{UInt64},
    timestamp::PlantTimestamp,
)
    lane = first(lanes)
    lane_state = first(lane_states)
    workspace = first(workspaces)
    if workspace.selected
        transition = _prevalidate_command_silence_transition(
            workspace.disposition,
            lane.authority_endpoint.endpoint,
            lane_state.endpoint,
            lane_state.application,
            timestamp,
        )
        workspace.silence[] = transition
        if transition.action == ApplySafeCommand
            current_sequence = @inbounds sequences[Int(lane.endpoint_slot)]
            current_sequence != typemax(UInt64) || _command_fanout_error(
                :publication_sequence_overflow,
                "effective-command publication sequence exceeds UInt64",
            )
            workspace.proposed_sequence = current_sequence + one(UInt64)
            _require_command_fanout_route_idle(
                lane.routes, lane_state.routes)
        elseif transition.action != FailOnCommandSilence
            _command_fanout_error(
                :silence_invariant,
                "scheduled command-silence transition has no executable action",
            )
        end
    end
    return _prevalidate_selected_command_silence!(
        Base.tail(lanes),
        Base.tail(lane_states),
        Base.tail(workspaces),
        sequences,
        timestamp,
    )
end

@inline function _selected_command_silence_transition(::Tuple{})
    _command_fanout_error(
        :silence_invariant,
        "command fanout has no selected command-silence transition",
    )
end

@inline function _selected_command_silence_transition(workspaces::Tuple)
    workspace = first(workspaces)
    workspace.selected && return workspace.silence[]
    return _selected_command_silence_transition(Base.tail(workspaces))
end

@inline _stage_selected_safe_command_authority!(::Tuple{}, ::Tuple{}) = nothing

@inline function _stage_selected_safe_command_authority!(
    lane_states::Tuple,
    workspaces::Tuple,
)
    lane_state = first(lane_states)
    workspace = first(workspaces)
    if workspace.selected
        _stage_safe_command!(lane_state.application.values)
        workspace.has_staged = true
    end
    return _stage_selected_safe_command_authority!(
        Base.tail(lane_states), Base.tail(workspaces))
end

@inline function _claim_selected_command_fanout!(
    ::Tuple{}, ::Tuple{}, ::Tuple{}, ::PlantTimestamp)
    return nothing
end

@inline function _claim_selected_command_fanout!(
    lanes::Tuple,
    lane_states::Tuple,
    workspaces::Tuple,
    timestamp::PlantTimestamp,
)
    lane = first(lanes)
    lane_state = first(lane_states)
    workspace = first(workspaces)
    if workspace.selected
        claim = claim_next_application_ready_command!(
            lane.authority_endpoint.endpoint, lane_state.endpoint, timestamp)
        claim === nothing && _command_fanout_error(
            :missing_command_claim,
            "selected command endpoint did not produce its due claim",
        )
        workspace.claim[] = claim
        workspace.has_claim = true
    end
    return _claim_selected_command_fanout!(
        Base.tail(lanes),
        Base.tail(lane_states),
        Base.tail(workspaces),
        timestamp,
    )
end

@inline function _stage_selected_command_authority!(
    ::Tuple{}, ::Tuple{}, ::Tuple{})
    return nothing
end

@inline function _stage_selected_command_authority!(
    lanes::Tuple,
    lane_states::Tuple,
    workspaces::Tuple,
)
    lane = first(lanes)
    lane_state = first(lane_states)
    workspace = first(workspaces)
    if workspace.selected
        staged = _stage_claimed_plant_command!(
            lane.authority_endpoint.endpoint,
            lane_state.endpoint,
            lane_state.application,
            workspace.claim[],
        )
        workspace.staged[] = staged
        workspace.has_staged = true
    end
    return _stage_selected_command_authority!(
        Base.tail(lanes),
        Base.tail(lane_states),
        Base.tail(workspaces),
    )
end

@inline _selected_command_fanout_status(::Tuple{}) = (false, false)

@inline function _selected_command_fanout_status(workspaces::Tuple)
    workspace = first(workspaces)
    tail_rejected, tail_failed =
        _selected_command_fanout_status(Base.tail(workspaces))
    if !(workspace.selected && workspace.has_staged)
        return tail_rejected, tail_failed
    end
    decision = workspace.staged[].decision
    return (
        tail_rejected || decision != _AcceptCommandCandidate,
        tail_failed || decision == _FailCommandCandidate,
    )
end

@inline function _prevalidate_rejected_command_fanout_finalization!(
    ::Tuple{}, ::Tuple{}, ::Tuple{}, ::Bool)
    return nothing
end

@inline function _prevalidate_rejected_command_fanout_finalization!(
    lanes::Tuple,
    lane_states::Tuple,
    workspaces::Tuple,
    atomic::Bool,
)
    lane = first(lanes)
    lane_state = first(lane_states)
    workspace = first(workspaces)
    if workspace.selected
        staged = workspace.staged[]
        if staged.decision == _AcceptCommandCandidate
            kind = RejectedCommand
            reason = CommandDispositionReason(:atomic_transaction_aborted)
        else
            kind = staged.decision == _FailCommandCandidate ?
                FailedCommand : RejectedCommand
            reason = staged.reason
        end
        atomic || staged.decision != _AcceptCommandCandidate ||
            _command_fanout_error(
                :application_invariant,
                "accepted single command entered rejection finalization",
            )
        workspace.finalization[] =
            _prevalidate_command_application_disposition(
                workspace.disposition,
                lane.authority_endpoint.endpoint,
                lane_state.endpoint,
                workspace.claim[],
                kind,
                reason,
            )
        workspace.has_finalization = true
    end
    return _prevalidate_rejected_command_fanout_finalization!(
        Base.tail(lanes),
        Base.tail(lane_states),
        Base.tail(workspaces),
        atomic,
    )
end

@inline _commit_rejected_command_fanout_finalization!(::Tuple{}) = nothing

@inline function _commit_rejected_command_fanout_finalization!(
    workspaces::Tuple,
)
    workspace = first(workspaces)
    if workspace.selected
        _commit_prevalidated_command_application_disposition!(
            workspace.finalization[])
        workspace.has_claim = false
    end
    return _commit_rejected_command_fanout_finalization!(
        Base.tail(workspaces))
end

@inline function _prepare_command_fanout_routes!(
    ::Tuple{},
    ::Tuple{},
    value,
    ::PlantTimestamp,
    ::UInt64,
)
    return nothing
end

@inline function _prepare_command_fanout_routes!(
    routes::Tuple,
    route_states::Tuple,
    value,
    timestamp::PlantTimestamp,
    sequence::UInt64,
)
    route = first(routes)
    publication = _effective_command_route_publication(
        route, timestamp, sequence)
    _prepare_effective_command_publication_route!(
        route, first(route_states), publication, value)
    return _prepare_command_fanout_routes!(
        Base.tail(routes),
        Base.tail(route_states),
        value,
        timestamp,
        sequence,
    )
end

@inline function _stage_command_fanout_routes!(
    ::Tuple{},
    ::Tuple{},
    value,
    ::PlantTimestamp,
    ::UInt64,
)
    return nothing
end

@inline function _stage_command_fanout_routes!(
    routes::Tuple,
    route_states::Tuple,
    value,
    timestamp::PlantTimestamp,
    sequence::UInt64,
)
    route = first(routes)
    publication = _effective_command_route_publication(
        route, timestamp, sequence)
    _stage_effective_command_publication_route!(
        route, first(route_states), publication, value)
    return _stage_command_fanout_routes!(
        Base.tail(routes),
        Base.tail(route_states),
        value,
        timestamp,
        sequence,
    )
end

@inline _reclaim_command_fanout_routes!(::Tuple{}, ::Tuple{}) = nothing

@inline function _reclaim_command_fanout_routes!(
    routes::Tuple,
    route_states::Tuple,
)
    _reclaim_effective_command_publication_route!(
        first(routes), first(route_states))
    return _reclaim_command_fanout_routes!(
        Base.tail(routes), Base.tail(route_states))
end

@inline _commit_command_fanout_routes!(::Tuple{}, ::Tuple{}) = nothing

@inline function _commit_command_fanout_routes!(
    routes::Tuple,
    route_states::Tuple,
)
    _commit_effective_command_publication_route!(
        first(routes), first(route_states))
    return _commit_command_fanout_routes!(
        Base.tail(routes), Base.tail(route_states))
end

@inline _abandon_command_fanout_routes!(::Tuple{}, ::Tuple{}) = true

@inline function _abandon_command_fanout_routes!(
    routes::Tuple,
    route_states::Tuple,
)
    proven = try
        _abandon_effective_command_publication_route!(
            first(routes), first(route_states))
    catch
        false
    end
    return _abandon_command_fanout_routes!(
        Base.tail(routes), Base.tail(route_states)) && proven
end

@inline function _prepare_selected_command_fanout_routes!(
    ::Tuple{}, ::Tuple{}, ::Tuple{}, ::PlantTimestamp)
    return nothing
end

@inline function _prepare_selected_command_fanout_routes!(
    lanes::Tuple,
    lane_states::Tuple,
    workspaces::Tuple,
    timestamp::PlantTimestamp,
)
    lane_state = first(lane_states)
    workspace = first(workspaces)
    if workspace.selected
        _prepare_command_fanout_routes!(
            first(lanes).routes,
            lane_state.routes,
            _staged_effective_command(lane_state.application),
            timestamp,
            workspace.proposed_sequence,
        )
    end
    return _prepare_selected_command_fanout_routes!(
        Base.tail(lanes),
        Base.tail(lane_states),
        Base.tail(workspaces),
        timestamp,
    )
end

@inline function _stage_selected_command_fanout_routes!(
    ::Tuple{}, ::Tuple{}, ::Tuple{}, ::PlantTimestamp)
    return nothing
end

@inline function _stage_selected_command_fanout_routes!(
    lanes::Tuple,
    lane_states::Tuple,
    workspaces::Tuple,
    timestamp::PlantTimestamp,
)
    lane_state = first(lane_states)
    workspace = first(workspaces)
    if workspace.selected
        _stage_command_fanout_routes!(
            first(lanes).routes,
            lane_state.routes,
            _staged_effective_command(lane_state.application),
            timestamp,
            workspace.proposed_sequence,
        )
    end
    return _stage_selected_command_fanout_routes!(
        Base.tail(lanes),
        Base.tail(lane_states),
        Base.tail(workspaces),
        timestamp,
    )
end

@inline function _reclaim_selected_command_fanout_routes!(
    ::Tuple{}, ::Tuple{}, ::Tuple{})
    return nothing
end

@inline function _reclaim_selected_command_fanout_routes!(
    lanes::Tuple,
    lane_states::Tuple,
    workspaces::Tuple,
)
    first(workspaces).selected && _reclaim_command_fanout_routes!(
        first(lanes).routes, first(lane_states).routes)
    return _reclaim_selected_command_fanout_routes!(
        Base.tail(lanes),
        Base.tail(lane_states),
        Base.tail(workspaces),
    )
end

@inline function _abandon_selected_command_fanout_routes!(
    ::Tuple{}, ::Tuple{}, ::Tuple{})
    return true
end

@inline function _abandon_selected_command_fanout_routes!(
    lanes::Tuple,
    lane_states::Tuple,
    workspaces::Tuple,
)
    proven = !first(workspaces).selected ||
        _abandon_command_fanout_routes!(
            first(lanes).routes, first(lane_states).routes)
    return _abandon_selected_command_fanout_routes!(
        Base.tail(lanes),
        Base.tail(lane_states),
        Base.tail(workspaces),
    ) && proven
end

@inline function _prevalidate_selected_command_fanout_finalization!(
    ::Tuple{}, ::Tuple{}, ::Tuple{})
    return nothing
end

@inline function _prevalidate_selected_command_fanout_finalization!(
    lanes::Tuple,
    lane_states::Tuple,
    workspaces::Tuple,
)
    lane = first(lanes)
    lane_state = first(lane_states)
    workspace = first(workspaces)
    if workspace.selected
        workspace.finalization[] =
            _prevalidate_command_application_disposition(
                workspace.disposition,
                lane.authority_endpoint.endpoint,
                lane_state.endpoint,
                workspace.claim[],
                AppliedCommand,
                workspace.staged[].reason,
            )
        workspace.has_finalization = true
    end
    return _prevalidate_selected_command_fanout_finalization!(
        Base.tail(lanes),
        Base.tail(lane_states),
        Base.tail(workspaces),
    )
end

@inline function _commit_selected_command_fanout_routes!(
    ::Tuple{}, ::Tuple{}, ::Tuple{})
    return nothing
end

@inline function _commit_selected_command_fanout_routes!(
    lanes::Tuple,
    lane_states::Tuple,
    workspaces::Tuple,
)
    first(workspaces).selected && _commit_command_fanout_routes!(
        first(lanes).routes, first(lane_states).routes)
    return _commit_selected_command_fanout_routes!(
        Base.tail(lanes),
        Base.tail(lane_states),
        Base.tail(workspaces),
    )
end

@inline function _commit_selected_command_authority!(
    ::Tuple{}, ::Tuple{}, ::Tuple{}, ::Memory{UInt64})
    return nothing
end

@inline function _commit_selected_command_authority!(
    lanes::Tuple,
    lane_states::Tuple,
    workspaces::Tuple,
    sequences::Memory{UInt64},
)
    lane = first(lanes)
    lane_state = first(lane_states)
    workspace = first(workspaces)
    if workspace.selected
        _commit_staged_application!(
            lane_state.application, lane_state.endpoint)
        @inbounds sequences[Int(lane.endpoint_slot)] =
            workspace.proposed_sequence
        _commit_prevalidated_command_application_disposition!(
            workspace.finalization[])
        workspace.has_claim = false
    end
    return _commit_selected_command_authority!(
        Base.tail(lanes),
        Base.tail(lane_states),
        Base.tail(workspaces),
        sequences,
    )
end

@inline function _commit_selected_safe_command_authority!(
    ::Tuple{}, ::Tuple{}, ::Tuple{}, ::Memory{UInt64})
    return nothing
end

@inline function _commit_selected_safe_command_authority!(
    lanes::Tuple,
    lane_states::Tuple,
    workspaces::Tuple,
    sequences::Memory{UInt64},
)
    lane = first(lanes)
    lane_state = first(lane_states)
    workspace = first(workspaces)
    if workspace.selected
        _commit_prevalidated_safe_command_silence!(
            lane_state.endpoint,
            lane_state.application,
            workspace.silence[],
        )
        @inbounds sequences[Int(lane.endpoint_slot)] =
            workspace.proposed_sequence
    end
    return _commit_selected_safe_command_authority!(
        Base.tail(lanes),
        Base.tail(lane_states),
        Base.tail(workspaces),
        sequences,
    )
end

@inline function _commit_selected_failed_command_silence!(
    ::Tuple{}, ::Tuple{}, ::Tuple{}, ::PlantTimestamp)
    return nothing
end

@inline function _commit_selected_failed_command_silence!(
    lanes::Tuple,
    lane_states::Tuple,
    workspaces::Tuple,
    timestamp::PlantTimestamp,
)
    lane = first(lanes)
    lane_state = first(lane_states)
    workspace = first(workspaces)
    if workspace.selected
        fail_pending_plant_commands!(
            workspace.disposition,
            lane.authority_endpoint.endpoint,
            lane_state.endpoint,
            timestamp;
            reason=:command_silence,
        )
        _commit_prevalidated_command_silence_latch!(
            lane_state.application, workspace.silence[])
    end
    return _commit_selected_failed_command_silence!(
        Base.tail(lanes),
        Base.tail(lane_states),
        Base.tail(workspaces),
        timestamp,
    )
end

@inline function _fail_selected_command_claims!(
    ::Tuple{}, ::Tuple{}, ::Tuple{}, ::CommandDispositionReason)
    return true
end

@inline function _fail_selected_command_claims!(
    lanes::Tuple,
    lane_states::Tuple,
    workspaces::Tuple,
    reason::CommandDispositionReason,
)
    lane = first(lanes)
    lane_state = first(lane_states)
    workspace = first(workspaces)
    released = true
    if workspace.selected && workspace.has_claim
        released = try
            _finish_command_application!(
                workspace.disposition,
                lane.authority_endpoint.endpoint,
                lane_state.endpoint,
                workspace.claim[],
                FailedCommand,
                reason,
            )
            workspace.has_claim = false
            true
        catch
            false
        end
    end
    return _fail_selected_command_claims!(
        Base.tail(lanes),
        Base.tail(lane_states),
        Base.tail(workspaces),
        reason,
    ) && released
end

@inline _command_fanout_failure_reason(error::PlantCommandError) =
    CommandDispositionReason(error.reason)
@inline _command_fanout_failure_reason(::Any) =
    CommandDispositionReason(:command_fanout_failure)

@inline function _fail_all_command_fanout_endpoints!(
    ::Tuple{},
    ::Tuple{},
    ::Tuple{},
    ::PlantTimestamp,
    ::CommandDispositionReason,
)
    return true
end

@inline function _fail_all_command_fanout_endpoints!(
    lanes::Tuple,
    lane_states::Tuple,
    workspaces::Tuple,
    timestamp::PlantTimestamp,
    reason::CommandDispositionReason,
)
    lane = first(lanes)
    lane_state = first(lane_states)
    workspace = first(workspaces)
    terminal_timestamp = max(
        timestamp, command_endpoint_timestamp(lane_state.endpoint))
    drained = try
        _append_failed_pending_plant_commands!(
            workspace.disposition,
            lane.authority_endpoint.endpoint,
            lane_state.endpoint,
            terminal_timestamp;
            reason,
        )
        true
    catch
        false
    end
    lane_state.endpoint.failed = true
    return _fail_all_command_fanout_endpoints!(
        Base.tail(lanes),
        Base.tail(lane_states),
        Base.tail(workspaces),
        timestamp,
        reason,
    ) && drained
end

@inline function _command_fanout_transaction_presentation_is_pending(
    state::CommandEndpointState,
    presentation::UInt64,
    transaction_id::UInt64,
)
    @inbounds for index in 1:state.pending_count
        metadata = state.slots[Int(state.calendar[index])]
        metadata.presentation == presentation &&
            metadata.transaction == transaction_id && return true
    end
    return false
end

@inline function _append_uncommitted_command_fanout_transaction_lane!(
    ::Tuple{},
    lane::_PreparedCommandFanoutLane,
    ::_CommandFanoutLaneState,
    ::_CommandFanoutLaneWorkspace,
    ::PlantTimestamp,
    ::UInt64,
    ::CommandDispositionReason,
)
    _command_admission_error(
        :transaction,
        :transaction_invariant,
        "prevalidated canonical command-fanout lane has no transaction " *
        "member ($(command_endpoint_id(lane.authority_endpoint)))",
    )
end


@inline function _append_uncommitted_command_fanout_transaction_lane!(
    commands::Tuple,
    lane::_PreparedCommandFanoutLane,
    lane_state::_CommandFanoutLaneState,
    workspace::_CommandFanoutLaneWorkspace,
    timestamp::PlantTimestamp,
    transaction_id::UInt64,
    reason::CommandDispositionReason,
)
    command = first(commands)
    if command_endpoint_id(command) ==
            command_endpoint_id(lane.authority_endpoint)
        plan = workspace.transaction_admission[]
        if !_command_fanout_transaction_presentation_is_pending(
            lane_state.endpoint, plan.presentation.value, transaction_id)
            disposition_workspace = workspace.disposition
            destination = disposition_workspace.count + 1
            destination <= length(disposition_workspace.dispositions) ||
                _command_admission_error(
                    :disposition,
                    :disposition_capacity,
                    "transaction commit failure exceeds disposition " *
                    "workspace capacity",
                )
            @inbounds disposition_workspace.dispositions[destination] =
                _command_disposition(
                    plan.presentation,
                    lane.authority_endpoint.endpoint,
                    command.sequence,
                    FailedCommand,
                    reason,
                    command.requested_effective_timestamp,
                    timestamp,
                )
            disposition_workspace.count = destination
            lane_state.endpoint.presentation_sequence =
                plan.presentation.value
            lane_state.endpoint.current_timestamp = timestamp
        end
        workspace.has_transaction_admission = false
        return nothing
    end
    return _append_uncommitted_command_fanout_transaction_lane!(
        Base.tail(commands),
        lane,
        lane_state,
        workspace,
        timestamp,
        transaction_id,
        reason,
    )
end

@inline function _append_uncommitted_command_fanout_transaction_members!(
    ::Tuple{},
    ::Tuple{},
    ::Tuple{},
    ::Tuple,
    ::PlantTimestamp,
    ::UInt64,
    ::CommandDispositionReason,
)
    return nothing
end


@inline function _append_uncommitted_command_fanout_transaction_members!(
    lanes::Tuple,
    lane_states::Tuple,
    workspaces::Tuple,
    commands::Tuple,
    timestamp::PlantTimestamp,
    transaction_id::UInt64,
    reason::CommandDispositionReason,
)
    lane = first(lanes)
    workspace = first(workspaces)
    workspace.selected && workspace.has_transaction_admission &&
        _append_uncommitted_command_fanout_transaction_lane!(
            commands,
            lane,
            first(lane_states),
            workspace,
            timestamp,
            transaction_id,
            reason,
        )
    return _append_uncommitted_command_fanout_transaction_members!(
        Base.tail(lanes),
        Base.tail(lane_states),
        Base.tail(workspaces),
        commands,
        timestamp,
        transaction_id,
        reason,
    )
end

function _fail_command_fanout_transaction_commit!(
    prepared::_PreparedCommandFanout,
    state::_CommandFanoutState,
    workspace::_CommandFanoutWorkspace,
    transaction::PlantCommandTransaction,
    timestamp::PlantTimestamp,
    transaction_id::UInt64,
)
    reason = CommandDispositionReason(:transaction_commit_failure)
    try
        _append_uncommitted_command_fanout_transaction_members!(
            prepared.lanes,
            state.lanes,
            workspace.lanes,
            transaction.commands,
            timestamp,
            transaction_id,
            reason,
        )
    catch
        # Preserve the original commit exception; the uncertain phase below
        # records that terminal-credit recovery could not be proven.
    end
    _fail_all_command_fanout_endpoints!(
        prepared.lanes,
        state.lanes,
        workspace.lanes,
        timestamp,
        reason,
    )
    _clear_command_fanout_selection!(workspace.lanes)
    state.authority.failed = true
    # An exception after the first calendar mutation is an ownership boundary:
    # terminal credits can be recovered, but exact committed state cannot be
    # proven from the throwing call site.
    state.phase = _CommandFanoutUncertain
    return nothing
end

function _fail_command_fanout_precommit!(
    prepared::_PreparedCommandFanout,
    state::_CommandFanoutState,
    workspace::_CommandFanoutWorkspace,
    error,
    timestamp::PlantTimestamp,
)
    routes_proven = _abandon_selected_command_fanout_routes!(
        prepared.lanes, state.lanes, workspace.lanes)
    claims_released = _fail_selected_command_claims!(
        prepared.lanes,
        state.lanes,
        workspace.lanes,
        _command_fanout_failure_reason(error),
    )
    endpoints_drained = _fail_all_command_fanout_endpoints!(
        prepared.lanes,
        state.lanes,
        workspace.lanes,
        timestamp,
        _command_fanout_failure_reason(error),
    )
    state.authority.failed = true
    state.phase = routes_proven && claims_released && endpoints_drained ?
        _CommandFanoutFailed : _CommandFanoutUncertain
    return nothing
end

function _fail_command_fanout_postcommit!(
    prepared::_PreparedCommandFanout,
    state::_CommandFanoutState,
    workspace::_CommandFanoutWorkspace,
    error,
    timestamp::PlantTimestamp,
)
    reason = _command_fanout_failure_reason(error)
    _fail_selected_command_claims!(
        prepared.lanes, state.lanes, workspace.lanes, reason)
    _fail_all_command_fanout_endpoints!(
        prepared.lanes,
        state.lanes,
        workspace.lanes,
        timestamp,
        reason,
    )
    state.authority.failed = true
    state.phase = _CommandFanoutUncertain
    return nothing
end

function _admit_command_fanout_transaction!(
    prepared::_PreparedCommandFanout,
    state::_CommandFanoutState,
    workspace::_CommandFanoutWorkspace,
    transaction::PlantCommandTransaction,
    timestamp::PlantTimestamp,
)
    _require_command_fanout_binding(prepared, state, workspace)
    _require_command_fanout_workspace_idle(workspace.lanes)
    commands = transaction.commands
    length(commands) <= typemax(UInt32) || _command_admission_error(
        :transaction,
        :member_count_overflow,
        "atomic transaction member count exceeds UInt32 range",
    )
    member_count = UInt32(length(commands))
    try
        _select_canonical_command_fanout_transaction_lanes!(
            prepared.lanes, workspace.lanes, commands)
    catch
        _clear_command_fanout_selection!(workspace.lanes)
        rethrow()
    end
    transaction_id = try
        _next_command_transaction_sequence(
            state.authority.command_transaction_sequence)
    catch
        _clear_command_fanout_selection!(workspace.lanes)
        rethrow()
    end
    scheduled::PlantTimestamp = try
        _with_completed_prepared_device_execution_context(
            prepared.authority.context) do
            _preflight_canonical_command_fanout_transaction!(
                prepared.lanes,
                state.lanes,
                workspace.lanes,
                commands,
                timestamp,
            )
        end
    catch error
        return _handle_command_fanout_transaction_preflight_failure!(
            prepared,
            state,
            workspace,
            transaction,
            timestamp,
            error,
        )::PlantCommandTransactionAdmission
    end

    state.phase = _CommandFanoutCommitting
    try
        _commit_canonical_command_fanout_transaction!(
            prepared.lanes,
            state.lanes,
            workspace.lanes,
            commands,
            timestamp,
            transaction_id,
            member_count,
        )
    catch error
        _fail_command_fanout_transaction_commit!(
            prepared,
            state,
            workspace,
            transaction,
            timestamp,
            transaction_id,
        )
        rethrow()
    end
    state.authority.command_transaction_sequence = transaction_id
    _clear_command_fanout_selection!(workspace.lanes)
    state.phase = _CommandFanoutIdle
    status = scheduled <= timestamp ?
        CommandAdmittedReady : CommandAdmittedPending
    return PlantCommandTransactionAdmission(
        transaction_id, status, member_count, scheduled)
end

function _apply_next_command_fanout!(
    prepared::_PreparedCommandFanout,
    state::_CommandFanoutState,
    workspace::_CommandFanoutWorkspace,
    endpoint::CommandEndpointID,
    timestamp::PlantTimestamp,
)
    _require_command_fanout_binding(prepared, state, workspace)
    _require_command_fanout_workspace_idle(workspace.lanes)
    transaction, declared_count = _command_fanout_due_metadata(
        prepared.lanes, state.lanes, endpoint, timestamp)
    expected_count = iszero(transaction) ? UInt32(1) : declared_count
    count = try
        _select_command_fanout_members!(
            prepared.lanes,
            state.lanes,
            workspace.lanes,
            endpoint,
            transaction,
            expected_count,
            timestamp,
        )
    catch
        _clear_command_fanout_selection!(workspace.lanes)
        rethrow()
    end
    count == Int(expected_count) || begin
        _clear_command_fanout_selection!(workspace.lanes)
        _command_fanout_error(
            :transaction_invariant,
            "atomic transaction yielded $count of $(Int(expected_count)) members",
        )
    end
    try
        _require_distinct_selected_command_fanout_optics(
            prepared.lanes, workspace.lanes)
        _require_selected_command_fanout_ready(
            prepared.lanes,
            state.lanes,
            workspace.lanes,
            state.authority.publication_sequences,
        )
    catch
        _clear_command_fanout_selection!(workspace.lanes)
        rethrow()
    end

    state.phase = _CommandFanoutApplying
    try
        _claim_selected_command_fanout!(
            prepared.lanes, state.lanes, workspace.lanes, timestamp)
        _with_completed_prepared_device_execution_context(
            prepared.authority.context) do
            _stage_selected_command_authority!(
                prepared.lanes, state.lanes, workspace.lanes)
        end
    catch error
        _fail_command_fanout_precommit!(
            prepared, state, workspace, error, timestamp)
        rethrow()
    end

    rejected, failed = _selected_command_fanout_status(workspace.lanes)
    if rejected
        try
            _prevalidate_rejected_command_fanout_finalization!(
                prepared.lanes,
                state.lanes,
                workspace.lanes,
                count > 1,
            )
        catch error
            _fail_command_fanout_precommit!(
                prepared, state, workspace, error, timestamp)
            rethrow()
        end
        _commit_rejected_command_fanout_finalization!(workspace.lanes)
        _clear_command_fanout_selection!(workspace.lanes)
        if failed
            drained = _fail_all_command_fanout_endpoints!(
                prepared.lanes,
                state.lanes,
                workspace.lanes,
                timestamp,
                CommandDispositionReason(:application_policy_failure),
            )
            state.authority.failed = true
            state.phase = drained ?
                _CommandFanoutFailed : _CommandFanoutUncertain
            _command_fanout_error(
                :application_policy_failure,
                "command application policy requires fail-stop",
            )
        end
        state.phase = _CommandFanoutIdle
        return nothing
    end

    try
        _prevalidate_selected_command_fanout_finalization!(
            prepared.lanes, state.lanes, workspace.lanes)
        _prepare_selected_command_fanout_routes!(
            prepared.lanes, state.lanes, workspace.lanes, timestamp)
        _stage_selected_command_fanout_routes!(
            prepared.lanes, state.lanes, workspace.lanes, timestamp)
        _reclaim_selected_command_fanout_routes!(
            prepared.lanes, state.lanes, workspace.lanes)
    catch error
        _fail_command_fanout_precommit!(
            prepared, state, workspace, error, timestamp)
        rethrow()
    end

    state.phase = _CommandFanoutCommitting
    try
        _commit_selected_command_fanout_routes!(
            prepared.lanes, state.lanes, workspace.lanes)
        _commit_selected_command_authority!(
            prepared.lanes,
            state.lanes,
            workspace.lanes,
            state.authority.publication_sequences,
        )
    catch error
        _fail_command_fanout_postcommit!(
            prepared, state, workspace, error, timestamp)
        rethrow()
    end
    _clear_command_fanout_selection!(workspace.lanes)
    state.phase = _CommandFanoutIdle
    return nothing
end

function _apply_command_silence_fanout!(
    prepared::_PreparedCommandFanout,
    state::_CommandFanoutState,
    workspace::_CommandFanoutWorkspace,
    endpoint::CommandEndpointID,
    timestamp::PlantTimestamp,
)
    _require_command_fanout_binding(prepared, state, workspace)
    _require_command_fanout_workspace_idle(workspace.lanes)
    try
        _select_command_fanout_endpoint!(
            prepared.lanes, workspace.lanes, endpoint)
        _prevalidate_selected_command_silence!(
            prepared.lanes,
            state.lanes,
            workspace.lanes,
            state.authority.publication_sequences,
            timestamp,
        )
    catch
        _clear_command_fanout_selection!(workspace.lanes)
        rethrow()
    end
    transition = _selected_command_silence_transition(workspace.lanes)

    state.phase = _CommandFanoutApplying
    if command_silence_action(transition) == FailOnCommandSilence
        try
            _commit_selected_failed_command_silence!(
                prepared.lanes,
                state.lanes,
                workspace.lanes,
                timestamp,
            )
        catch error
            _fail_command_fanout_postcommit!(
                prepared, state, workspace, error, timestamp)
            rethrow()
        end
        drained = _fail_all_command_fanout_endpoints!(
            prepared.lanes,
            state.lanes,
            workspace.lanes,
            timestamp,
            CommandDispositionReason(:command_silence),
        )
        _clear_command_fanout_selection!(workspace.lanes)
        state.authority.failed = true
        state.phase = drained ?
            _CommandFanoutFailed : _CommandFanoutUncertain
        drained || _command_fanout_error(
            :fail_stop_drain_incomplete,
            "command-silence fail-stop could not drain every endpoint",
        )
        return transition
    end

    try
        _with_completed_prepared_device_execution_context(
            prepared.authority.context) do
            _stage_selected_safe_command_authority!(
                state.lanes, workspace.lanes)
        end
        _prepare_selected_command_fanout_routes!(
            prepared.lanes, state.lanes, workspace.lanes, timestamp)
        _stage_selected_command_fanout_routes!(
            prepared.lanes, state.lanes, workspace.lanes, timestamp)
        _reclaim_selected_command_fanout_routes!(
            prepared.lanes, state.lanes, workspace.lanes)
    catch error
        _fail_command_fanout_precommit!(
            prepared, state, workspace, error, timestamp)
        rethrow()
    end

    state.phase = _CommandFanoutCommitting
    try
        _commit_selected_command_fanout_routes!(
            prepared.lanes, state.lanes, workspace.lanes)
        _commit_selected_safe_command_authority!(
            prepared.lanes,
            state.lanes,
            workspace.lanes,
            state.authority.publication_sequences,
        )
    catch error
        _fail_command_fanout_postcommit!(
            prepared, state, workspace, error, timestamp)
        rethrow()
    end
    _clear_command_fanout_selection!(workspace.lanes)
    state.phase = _CommandFanoutIdle
    return transition
end

@inline function _apply_next_command_fanout!(
    prepared::_PreparedCommandFanout,
    state::_CommandFanoutState,
    workspace::_CommandFanoutWorkspace,
    endpoint::Symbol,
    timestamp::PlantTimestamp,
)
    return _apply_next_command_fanout!(
        prepared, state, workspace, CommandEndpointID(endpoint), timestamp)
end

@inline function _apply_command_silence_fanout!(
    prepared::_PreparedCommandFanout,
    state::_CommandFanoutState,
    workspace::_CommandFanoutWorkspace,
    endpoint::Symbol,
    timestamp::PlantTimestamp,
)
    return _apply_command_silence_fanout!(
        prepared, state, workspace, CommandEndpointID(endpoint), timestamp)
end

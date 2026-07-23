#
# Effective plant-command state and replayable command silence
#
# Admission owns presentation validation, sequence history, and the bounded
# future calendar. This layer separately owns the currently effective command,
# application-stage validation, and scheduled plant-time silence transitions.
# Physical optic mutation, multi-optic atomicity, event-loop composition, HIL
# liveness, and transports remain outside this layer.
#

abstract type _AbstractCommandApplicationValues end

mutable struct _ScalarCommandApplicationValues{T,S} <:
    _AbstractCommandApplicationValues
    effective::T
    staging::T
    safe::S
end

mutable struct _ArrayCommandApplicationValues{A<:AbstractArray,S} <:
    _AbstractCommandApplicationValues
    effective::A
    staging::A
    safe::S
end

"""
Single-writer effective-command state for one prepared endpoint.

Construct this qualified state once from a prepared endpoint, its exact
admission-state owner, and one initial effective command, before that owner has
successfully admitted a command. `safe_command` is required only by
`ApplySafeCommand`. Array values are copied into endpoint-backend storage
during construction and are borrowed read-only through `effective_command`
until the next application or safe transition.
"""
mutable struct CommandApplicationState{V<:_AbstractCommandApplicationValues}
    binding::_CommandEndpointBinding
    endpoint_state_binding::_CommandEndpointStateBinding
    values::V
    baseline_timestamp::PlantTimestamp
    last_application_timestamp::PlantTimestamp
    last_silence_origin_timestamp::PlantTimestamp
    has_silence_transition::Bool
end

@inline function _require_new_command_application_state(
    endpoint_state::CommandEndpointState)
    endpoint_state.has_application_state && _command_admission_error(
        :preparation, :application_state_exists,
        "command endpoint state already owns an application state")
    (endpoint_state.has_admission || !iszero(endpoint_state.active_count)) &&
        _command_admission_error(:preparation, :active_endpoint_state,
            "construct command application state before the endpoint's " *
            "first successful admission")
    return nothing
end

@inline function _effective_seed_uses_command_bounds(
    schema::PlantCommandSchema)
    command_semantics(schema) == AbsoluteCommand && return true
    return command_value_policy(schema).range_stage == EnforceOnApplication
end

@inline function _validate_effective_seed_bounds(
    schema::PlantCommandSchema, payload, ::UnboundedCommandValues, label)
    return payload
end

function _validate_effective_seed_bounds(schema::PlantCommandSchema,
    payload, bounds::UniformCommandBounds, label)
    _effective_seed_uses_command_bounds(schema) || return payload
    _all_command_values_in_bounds(payload, bounds) ||
        _command_admission_error(:preparation, :invalid_effective_command,
            "$label contains a value outside the effective command bounds")
    return payload
end

function _validate_effective_seed(schema::PlantCommandSchema{T,0}, payload,
    label::AbstractString) where {T}
    typeof(payload) === T || _command_admission_error(
        :preparation, :invalid_effective_command,
        "$label must have exact scalar type $T; got $(typeof(payload))")
    isfinite(payload) || _command_admission_error(
        :preparation, :invalid_effective_command,
        "$label must be finite")
    return _validate_effective_seed_bounds(schema, payload,
        command_bounds(schema), label)
end

function _validate_effective_seed(::PlantCommandSchema{T,0},
    payload::AbstractArray{T,0}, label::AbstractString) where {T}
    _command_admission_error(:preparation, :invalid_effective_command,
        "$label must have exact scalar type $T; got $(typeof(payload))")
end

function _validate_effective_seed(schema::PlantCommandSchema{T,N},
    payload::AbstractArray{T,N}, label::AbstractString) where {T,N}
    size(payload) == command_dimensions(schema) || _command_admission_error(
        :preparation, :invalid_effective_command,
        "$label has size $(size(payload)); expected " *
        "$(command_dimensions(schema))")
    _all_command_values_finite(payload) || _command_admission_error(
        :preparation, :invalid_effective_command,
        "$label must contain only finite values")
    return _validate_effective_seed_bounds(schema, payload,
        command_bounds(schema), label)
end

function _validate_effective_seed(schema::PlantCommandSchema{T,N}, payload,
    label::AbstractString) where {T,N}
    _command_admission_error(:preparation, :invalid_effective_command,
        "$label must be an AbstractArray{$T,$N}; got $(typeof(payload))")
end

@inline function _require_safe_command_configuration(
    policy::CommandSilencePolicy, ::Nothing)
    policy.action == ApplySafeCommand && _command_admission_error(
        :preparation, :missing_safe_command,
        "ApplySafeCommand requires one prepared safe command")
    return nothing
end

@inline function _require_safe_command_configuration(
    policy::CommandSilencePolicy, safe_command)
    policy.action == ApplySafeCommand || _command_admission_error(
        :preparation, :unexpected_safe_command,
        "a safe command is valid only for ApplySafeCommand")
    return safe_command
end

function _prepare_command_application_values(
    endpoint::PreparedCommandEndpoint{<:PlantCommandSchema{T,0}},
    initial_command, safe_command) where {T}
    schema = command_schema(endpoint)
    initial = _validate_effective_seed(schema, initial_command,
        "initial effective command")
    safe = safe_command === nothing ? nothing :
        _validate_effective_seed(schema, safe_command, "safe command")
    return _ScalarCommandApplicationValues(initial, initial, safe)
end

function _prepare_command_application_values(
    endpoint::PreparedCommandEndpoint{<:PlantCommandSchema{T,N}},
    initial_command, safe_command) where {T,N}
    schema = command_schema(endpoint)
    initial = _validate_effective_seed(schema, initial_command,
        "initial effective command")
    safe = safe_command === nothing ? nothing :
        _validate_effective_seed(schema, safe_command, "safe command")
    effective = allocate_array(endpoint.backend, T, schema.dimensions...)
    staging = allocate_array(endpoint.backend, T, schema.dimensions...)
    copyto!(effective, initial)
    if safe === nothing
        return _ArrayCommandApplicationValues(effective, staging, nothing)
    end
    prepared_safe = allocate_array(endpoint.backend, T, schema.dimensions...)
    copyto!(prepared_safe, safe)
    return _ArrayCommandApplicationValues(effective, staging, prepared_safe)
end

function CommandApplicationState(endpoint::PreparedCommandEndpoint,
    endpoint_state::CommandEndpointState, initial_command;
    safe_command=nothing)
    _require_command_endpoint_binding(endpoint, endpoint_state)
    _require_operational_command_endpoint(endpoint_state)
    _require_idle_command_endpoint(endpoint_state)
    _require_new_command_application_state(endpoint_state)
    policy = command_silence_policy(command_schema(endpoint))
    _require_safe_command_configuration(policy, safe_command)
    values = _prepare_command_application_values(endpoint, initial_command,
        safe_command)
    timestamp = command_endpoint_timestamp(endpoint_state)
    state = CommandApplicationState(getfield(endpoint, :binding),
        getfield(endpoint_state, :state_binding), values, timestamp,
        timestamp, timestamp, false)
    endpoint_state.has_application_state = true
    return state
end

@inline effective_command(state::CommandApplicationState) =
    _effective_command(state.values)
@inline _effective_command(values::_ScalarCommandApplicationValues) =
    values.effective
@inline _effective_command(values::_ArrayCommandApplicationValues) =
    values.effective
@inline last_command_application_timestamp(state::CommandApplicationState) =
    state.last_application_timestamp

@inline function _require_command_application_binding(
    endpoint::PreparedCommandEndpoint, endpoint_state::CommandEndpointState,
    state::CommandApplicationState)
    getfield(endpoint, :binding) === state.binding ||
        _command_admission_error(:application, :foreign_application_state,
            "effective-command state belongs to another prepared endpoint")
    getfield(endpoint_state, :state_binding) ===
        state.endpoint_state_binding || _command_admission_error(
            :application, :foreign_endpoint_state,
            "effective-command state belongs to another endpoint-state owner")
    return nothing
end

@enum _CommandCandidateDecision::UInt8 begin
    _AcceptCommandCandidate = 0x01
    _RejectCommandCandidate = 0x02
    _FailCommandCandidate = 0x03
end

struct _StagedCommandApplication
    decision::_CommandCandidateDecision
    reason::CommandDispositionReason
end

@inline function _invalid_candidate_decision(action::InvalidCommandAction)
    action == FailOnInvalidCommand && return _FailCommandCandidate
    return _RejectCommandCandidate
end

@inline function _nonfinite_application_reason(action::InvalidCommandAction)
    action == FailOnInvalidCommand &&
        return CommandDispositionReason(:nonfinite_failure)
    return CommandDispositionReason(:nonfinite_rejected)
end

@inline function _range_application_reason(action::InvalidCommandAction)
    action == FailOnInvalidCommand &&
        return CommandDispositionReason(:out_of_range_failure)
    return CommandDispositionReason(:out_of_range_rejected)
end

@inline function _validate_application_candidate_finite(
    schema::PlantCommandSchema, candidate)
    _all_command_values_finite(candidate) && return nothing
    action = command_value_policy(schema).nonfinite
    return (_invalid_candidate_decision(action),
        _nonfinite_application_reason(action))
end

@inline function _validate_scalar_application_candidate_range(
    schema::PlantCommandSchema, candidate, ::UnboundedCommandValues)
    return (_AcceptCommandCandidate, CommandDispositionReason(:applied),
        candidate)
end

@inline function _validate_scalar_application_candidate_range(
    schema::PlantCommandSchema, candidate, bounds::UniformCommandBounds)
    policy = command_value_policy(schema)
    policy.range_stage == EnforceOnApplication || return (
        _AcceptCommandCandidate, CommandDispositionReason(:applied), candidate)
    bounds.lower <= candidate <= bounds.upper && return (
        _AcceptCommandCandidate, CommandDispositionReason(:applied), candidate)
    policy.out_of_range == ClipInvalidCommand && return (
        _AcceptCommandCandidate, CommandDispositionReason(:applied_clipped),
        clamp(candidate, bounds.lower, bounds.upper))
    return (_invalid_candidate_decision(policy.out_of_range),
        _range_application_reason(policy.out_of_range), candidate)
end

@inline function _stage_application_candidate!(
    values::_ScalarCommandApplicationValues,
    schema::PlantCommandSchema, payload)
    candidate = command_semantics(schema) == AbsoluteCommand ?
        payload : values.effective + payload
    finite_result = _validate_application_candidate_finite(schema, candidate)
    if finite_result !== nothing
        values.staging = candidate
        return _StagedCommandApplication(finite_result[1], finite_result[2])
    end
    decision, reason, staged =
        _validate_scalar_application_candidate_range(schema, candidate,
            command_bounds(schema))
    values.staging = staged
    return _StagedCommandApplication(decision, reason)
end

@inline function _stage_array_application_candidate!(
    staging::AbstractArray, effective::AbstractArray, payload::AbstractArray,
    semantics::CommandValueSemantics)
    if semantics == AbsoluteCommand
        copyto!(staging, payload)
    else
        @. staging = effective + payload
    end
    return nothing
end

@inline function _validate_array_application_candidate_range!(
    schema::PlantCommandSchema, candidate::AbstractArray,
    ::UnboundedCommandValues)
    return (_AcceptCommandCandidate, CommandDispositionReason(:applied))
end

@inline function _validate_array_application_candidate_range!(
    schema::PlantCommandSchema, candidate::AbstractArray,
    bounds::UniformCommandBounds)
    policy = command_value_policy(schema)
    policy.range_stage == EnforceOnApplication || return (
        _AcceptCommandCandidate, CommandDispositionReason(:applied))
    _all_command_values_in_bounds(candidate, bounds) && return (
        _AcceptCommandCandidate, CommandDispositionReason(:applied))
    if policy.out_of_range == ClipInvalidCommand
        clamp_array!(candidate, bounds.lower, bounds.upper)
        return (_AcceptCommandCandidate,
            CommandDispositionReason(:applied_clipped))
    end
    return (_invalid_candidate_decision(policy.out_of_range),
        _range_application_reason(policy.out_of_range))
end

@inline function _stage_application_candidate!(
    values::_ArrayCommandApplicationValues,
    schema::PlantCommandSchema, payload::AbstractArray)
    _stage_array_application_candidate!(values.staging, values.effective,
        payload, command_semantics(schema))
    finite_result = _validate_application_candidate_finite(schema,
        values.staging)
    finite_result === nothing || return _StagedCommandApplication(
        finite_result[1], finite_result[2])
    decision, reason = _validate_array_application_candidate_range!(schema,
        values.staging, command_bounds(schema))
    return _StagedCommandApplication(decision, reason)
end

@inline function _commit_application_candidate!(
    values::_ScalarCommandApplicationValues)
    values.effective = values.staging
    return nothing
end

@inline function _commit_application_candidate!(
    values::_ArrayCommandApplicationValues)
    values.effective, values.staging = values.staging, values.effective
    return nothing
end

@inline _staged_effective_command(values::_ScalarCommandApplicationValues) =
    values.staging
@inline _staged_effective_command(values::_ArrayCommandApplicationValues) =
    values.staging

@inline function _commit_staged_application!(
    state::CommandApplicationState,
    endpoint_state::CommandEndpointState)
    _commit_application_candidate!(state.values)
    state.last_application_timestamp =
        command_endpoint_timestamp(endpoint_state)
    return nothing
end

function _stage_claimed_plant_command!(
    endpoint::PreparedCommandEndpoint,
    endpoint_state::CommandEndpointState,
    application_state::CommandApplicationState,
    claim::PlantCommandApplicationClaim)
    _require_command_endpoint_binding(endpoint, endpoint_state)
    _require_command_application_binding(endpoint, endpoint_state,
        application_state)
    _require_operational_command_endpoint(endpoint_state)
    slot, metadata = _require_current_command_claim(
        endpoint, endpoint_state, claim)
    metadata.scheduled_timestamp == endpoint_state.current_timestamp ||
        _command_admission_error(:application,
            :missed_application_timestamp,
            "application-ready command was claimed after its immutable " *
            "scheduled plant timestamp")
    payload = _claimed_command_payload(endpoint_state.payloads, slot)
    return _stage_application_candidate!(application_state.values,
        command_schema(endpoint), payload)
end

@inline _staged_effective_command(state::CommandApplicationState) =
    _staged_effective_command(state.values)

@noinline function _handle_command_application_staging_error!(
    workspace::CommandDispositionWorkspace,
    endpoint::PreparedCommandEndpoint,
    endpoint_state::CommandEndpointState,
    claim::PlantCommandApplicationClaim,
    error::InterruptException)
    throw(error)
end

@noinline function _fail_missed_command_application_timestamp!(
    workspace::CommandDispositionWorkspace,
    endpoint::PreparedCommandEndpoint,
    endpoint_state::CommandEndpointState,
    claim::PlantCommandApplicationClaim)
    _finish_command_application!(workspace, endpoint, endpoint_state, claim,
        FailedCommand,
        CommandDispositionReason(:missed_application_timestamp))
    _command_admission_error(:application, :missed_application_timestamp,
        "application-ready command was claimed after its immutable " *
        "scheduled plant timestamp")
end

@noinline function _handle_command_application_staging_error!(
    workspace::CommandDispositionWorkspace,
    endpoint::PreparedCommandEndpoint,
    endpoint_state::CommandEndpointState,
    claim::PlantCommandApplicationClaim,
    error)
    _finish_command_application!(workspace, endpoint, endpoint_state, claim,
        FailedCommand, CommandDispositionReason(:application_storage_failure))
    _command_admission_error(:application, :application_storage_failure,
        "failed to stage the effective plant command ($(typeof(error)))")
end

"""
    apply_claimed_plant_command!(workspace, endpoint, endpoint_state,
        application_state, claim)

Transactionally combine the outstanding claim with the separately held
effective command, enforce application-stage finite/range policy, and publish
one terminal disposition at its exact scheduled timestamp. A claim made after
that timestamp fails rather than backdating state. Absolute commands replace
the held value; incremental commands add to it. Rejection or failure leaves
the held value unchanged. Array application swaps preallocated staging and
effective buffers.
"""
function apply_claimed_plant_command!(
    workspace::CommandDispositionWorkspace,
    endpoint::PreparedCommandEndpoint,
    endpoint_state::CommandEndpointState,
    application_state::CommandApplicationState,
    claim::PlantCommandApplicationClaim)
    _require_command_endpoint_binding(endpoint, workspace)
    _require_empty_command_dispositions(workspace)
    staged = try
        _stage_claimed_plant_command!(endpoint, endpoint_state,
            application_state, claim)
    catch error
        return _handle_staged_command_application_error!(workspace, endpoint,
            endpoint_state, claim, error)
    end
    decision = staged.decision
    reason = staged.reason
    if decision == _AcceptCommandCandidate
        _commit_staged_application!(application_state, endpoint_state)
        return _finish_command_application!(workspace, endpoint,
            endpoint_state, claim, AppliedCommand, reason)
    end
    kind = decision == _FailCommandCandidate ? FailedCommand : RejectedCommand
    disposition = _finish_command_application!(workspace, endpoint,
        endpoint_state, claim, kind, reason)
    decision == _FailCommandCandidate && _command_admission_error(
        :application, reason.name,
        "application-stage command policy requires structural failure")
    return disposition
end

"""One replayable plant-time command-silence state transition."""
struct PlantCommandSilenceTransition
    endpoint::CommandEndpointID
    action::CommandSilenceAction
    age_origin::CommandAgeOrigin
    origin_timestamp::PlantTimestamp
    deadline_timestamp::PlantTimestamp
    transition_timestamp::PlantTimestamp
end

@inline command_endpoint_id(value::PlantCommandSilenceTransition) =
    value.endpoint
@inline command_silence_action(value::PlantCommandSilenceTransition) =
    value.action
@inline command_silence_age_origin(value::PlantCommandSilenceTransition) =
    value.age_origin
@inline command_silence_origin_timestamp(
    value::PlantCommandSilenceTransition) = value.origin_timestamp
@inline command_silence_deadline(value::PlantCommandSilenceTransition) =
    value.deadline_timestamp
@inline command_silence_transition_timestamp(
    value::PlantCommandSilenceTransition) = value.transition_timestamp
@inline command_silence_age(value::PlantCommandSilenceTransition) =
    value.transition_timestamp - value.origin_timestamp

@inline function _command_silence_origin_timestamp(
    age_origin::CommandAgeOrigin, endpoint_state::CommandEndpointState,
    application_state::CommandApplicationState)
    if age_origin == AgeFromAdmission
        endpoint_state.has_admission &&
            return endpoint_state.last_admission_timestamp
        return application_state.baseline_timestamp
    end
    return application_state.last_application_timestamp
end

"""
Return the next scheduled plant-time silence transition, or `nothing` for
indefinite hold, a terminally failed endpoint, or an already-latched transition
whose selected age origin has not advanced.
"""
function next_command_silence_timestamp(endpoint::PreparedCommandEndpoint,
    endpoint_state::CommandEndpointState,
    application_state::CommandApplicationState)
    _require_command_endpoint_binding(endpoint, endpoint_state)
    _require_command_application_binding(endpoint, endpoint_state,
        application_state)
    endpoint_state.failed && return nothing
    policy = command_silence_policy(command_schema(endpoint))
    policy.action == HoldLastCommand && return nothing
    origin = _command_silence_origin_timestamp(policy.age_origin,
        endpoint_state, application_state)
    application_state.has_silence_transition &&
        origin <= application_state.last_silence_origin_timestamp &&
        return nothing
    return origin + policy.timeout
end

@inline function _stage_safe_command!(
    values::_ScalarCommandApplicationValues)
    values.staging = values.safe
    return nothing
end

@inline function _stage_safe_command!(values::_ArrayCommandApplicationValues)
    copyto!(values.staging, values.safe)
    return nothing
end

@noinline _handle_safe_command_staging_error(error::InterruptException) =
    throw(error)

@noinline function _handle_safe_command_staging_error(error)
    _command_admission_error(:silence, :safe_storage_failure,
        "failed to stage the prepared safe command ($(typeof(error)))")
end

@noinline function _handle_staged_command_application_error!(
    workspace::CommandDispositionWorkspace,
    endpoint::PreparedCommandEndpoint,
    endpoint_state::CommandEndpointState,
    claim::PlantCommandApplicationClaim,
    error::PlantCommandError)
    error.reason === :missed_application_timestamp &&
        return _fail_missed_command_application_timestamp!(workspace,
            endpoint, endpoint_state, claim)
    return _handle_command_application_staging_error!(workspace, endpoint,
        endpoint_state, claim, error)
end

@noinline function _handle_staged_command_application_error!(
    workspace::CommandDispositionWorkspace,
    endpoint::PreparedCommandEndpoint,
    endpoint_state::CommandEndpointState,
    claim::PlantCommandApplicationClaim,
    error)
    return _handle_command_application_staging_error!(workspace, endpoint,
        endpoint_state, claim, error)
end

"""
    apply_command_silence_transition!(workspace, endpoint, endpoint_state,
        application_state, timestamp)

Apply the exact next safe-command or fail-on-silence event. `timestamp` must
equal `next_command_silence_timestamp`. Any command due at the same timestamp
must be claimed and resolved first. Safe transition is transactional; failure
drains every still-pending command with `:command_silence` and marks the
endpoint terminally failed.
"""
function apply_command_silence_transition!(
    workspace::CommandDispositionWorkspace,
    endpoint::PreparedCommandEndpoint,
    endpoint_state::CommandEndpointState,
    application_state::CommandApplicationState,
    timestamp::PlantTimestamp)
    _require_command_endpoint_binding(endpoint, endpoint_state)
    _require_command_endpoint_binding(endpoint, workspace)
    _require_command_application_binding(endpoint, endpoint_state,
        application_state)
    _require_operational_command_endpoint(endpoint_state)
    _require_empty_command_dispositions(workspace)
    _require_idle_command_endpoint(endpoint_state)
    _require_forward_command_timestamp(endpoint_state, timestamp)
    expected = next_command_silence_timestamp(endpoint, endpoint_state,
        application_state)
    expected === nothing && _command_admission_error(
        :silence, :silence_not_scheduled,
        "command endpoint has no pending silence transition")
    timestamp == expected || _command_admission_error(
        :silence, :unexpected_silence_timestamp,
        "command silence transition expected at $expected; got $timestamp")
    key = next_command_order_key(endpoint, endpoint_state)
    key !== nothing && command_scheduled_timestamp(key) <= timestamp &&
        _command_admission_error(:silence, :commands_due,
            "resolve application-ready commands before an equal-time " *
            "command-silence transition")

    policy = command_silence_policy(command_schema(endpoint))
    origin = _command_silence_origin_timestamp(policy.age_origin,
        endpoint_state, application_state)
    if policy.action == ApplySafeCommand
        staged = try
            _stage_safe_command!(application_state.values)
        catch error
            _handle_safe_command_staging_error(error)
        end
        _commit_application_candidate!(application_state.values)
        endpoint_state.current_timestamp = timestamp
    else
        fail_pending_plant_commands!(workspace, endpoint, endpoint_state,
            timestamp; reason=:command_silence)
        endpoint_state.failed = true
    end
    application_state.last_silence_origin_timestamp = origin
    application_state.has_silence_transition = true
    return PlantCommandSilenceTransition(command_endpoint_id(endpoint),
        policy.action, policy.age_origin, origin, expected, timestamp)
end

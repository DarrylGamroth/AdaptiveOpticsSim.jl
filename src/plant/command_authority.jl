#
# Explicitly placed sole command authority for mixed-resource plants
#
# The authority owns command admission, application semantics, silence policy,
# and effective-publication sequence. Target-local physical-optic owners remain
# role-neutral replicas. Publication routes are composed below this preparation
# boundary rather than allowing replicas to reinterpret controller commands.
#

mutable struct _PreparedCommandAuthorityBinding end
mutable struct _PreparedCommandAuthorityToken end
const _PREPARED_COMMAND_AUTHORITY_TOKEN = _PreparedCommandAuthorityToken()

"""Run-immutable preparation for the sole mixed-plant command authority."""
struct PreparedCommandAuthority{I,D,X,O,E}
    binding::_PreparedCommandAuthorityBinding
    identity::I
    target::D
    context::X
    optic_definitions::O
    endpoints::E

    function PreparedCommandAuthority(
        token::_PreparedCommandAuthorityToken,
        binding::_PreparedCommandAuthorityBinding,
        identity::I,
        target::D,
        context::X,
        optic_definitions::O,
        endpoints::E,
    ) where {I,D,X,O,E}
        token === _PREPARED_COMMAND_AUTHORITY_TOKEN || throw(
            ArgumentError("invalid internal prepared command-authority token"))
        return new{I,D,X,O,E}(
            binding,
            identity,
            target,
            context,
            optic_definitions,
            endpoints,
        )
    end
end

"""Single-writer command-admission and effective-value state."""
mutable struct CommandAuthorityState{S,A}
    binding::_PreparedCommandAuthorityBinding
    endpoint_states::S
    application_states::A
    publication_sequences::Memory{UInt64}
    command_transaction_sequence::UInt64
    failed::Bool
end

"""Caller-owned fixed disposition storage for one command authority."""
struct CommandAuthorityWorkspace{W}
    binding::_PreparedCommandAuthorityBinding
    dispositions::W
end

@inline compute_device(authority::PreparedCommandAuthority) = authority.target
@inline command_authority_target(authority::PreparedCommandAuthority) =
    authority.target
@inline command_authority_identity(authority::PreparedCommandAuthority) =
    authority.identity
@inline command_authority_failed(state::CommandAuthorityState) = state.failed

@noinline function _command_authority_preparation_error(
    reason::Symbol,
    message::AbstractString,
)
    throw(PlantPreparationError(:command_authority, reason, String(message)))
end

@inline _require_command_authority_target(target::HostComputeDevice) = target
@inline _require_command_authority_target(
    target::AcceleratorComputeDevice) = target

function _require_command_authority_target(target)
    _command_authority_preparation_error(
        :invalid_target,
        "command-authority target must be an exact HostComputeDevice or " *
        "AcceleratorComputeDevice; got $(typeof(target))",
    )
end

function _require_command_authority_target_available(
    target::AbstractComputeDevice,
)
    availability = compute_device_availability(target)
    compute_device_is_available(availability) && return target
    reason = compute_device_unavailable_reason(availability)
    _command_authority_preparation_error(
        :unavailable_target,
        "exact command-authority target $target is unavailable" *
        (isnothing(reason) ? "" : " ($reason)"),
    )
end

@inline function _require_gate9a_command_authority_target(
    assignment::ResolvedPlantPartitionAssignment,
    target::HostComputeDevice,
)
    return target
end

@inline function _require_gate9a_accelerator_pair(
    ::HostComputeDevice,
    ::AcceleratorComputeDevice,
)
    return nothing
end

function _require_gate9a_accelerator_pair(
    partition_target::AcceleratorComputeDevice,
    target::AcceleratorComputeDevice,
)
    partition_target == target && return nothing
    _command_authority_preparation_error(
        :multiple_accelerators,
        "Gate 9A permits at most one exact accelerator across path and " *
        "command-authority placement; got $partition_target and $target",
    )
end

function _require_gate9a_command_authority_target(
    assignment::ResolvedPlantPartitionAssignment,
    target::AcceleratorComputeDevice,
)
    for partition_target in partition_targets(assignment)
        _require_gate9a_accelerator_pair(partition_target, target)
    end
    return target
end

function _command_authority_context(
    target::AbstractComputeDevice,
    atmosphere_target::AbstractComputeDevice,
    atmosphere_context,
)
    target == atmosphere_target && return atmosphere_context
    context = _prepare_device_execution_context(target)
    _prepared_device_execution_compute_device(context) == target ||
        _command_authority_preparation_error(
            :execution_context_target,
            "prepared command-authority context does not match its exact target",
        )
    return context
end

function _prepare_command_authority_endpoints(
    definition::PlantDefinition,
    configurations,
    target::AbstractComputeDevice,
)
    optic_definitions = _canonical_controllable_optic_definitions(definition)
    endpoints = _prepare_plant_command_endpoints(
        definition, configurations, optic_definitions, target)
    return (
        _partition_registry(optic_definitions, ControllableOpticDefinition),
        _partition_registry(endpoints, _PreparedPlantCommandEndpoint),
    )
end

"""
    _prepare_command_authority(definition, identity, target, configurations,
        atmosphere_target, atmosphere_context)

Prepare the sole command authority on one explicit exact target. The target is
independent of atmosphere placement and need not own a path partition. Gate 9A
still permits at most one accelerator across the complete placement.

Endpoint capacity, sequence history, initial effective values, and safe-command
policy are copied from the same sealed command configurations used to prepare
every target-local replica.
"""
function _prepare_command_authority(
    definition::PlantDefinition,
    identity::CommandAuthorityIdentity,
    target,
    configurations,
    atmosphere_target::AbstractComputeDevice,
    atmosphere_context,
)
    exact_target = _require_command_authority_target(target)
    _require_command_authority_target_available(exact_target)
    context = _command_authority_context(
        exact_target, atmosphere_target, atmosphere_context)
    optic_definitions, endpoints =
        _with_completed_prepared_device_execution_context(context) do
            _prepare_command_authority_endpoints(
                definition, configurations, exact_target)
        end
    authority = PreparedCommandAuthority(
        _PREPARED_COMMAND_AUTHORITY_TOKEN,
        _PreparedCommandAuthorityBinding(),
        identity,
        exact_target,
        context,
        optic_definitions,
        endpoints,
    )
    return _require_exact_prepared_command_authority_target(
        authority, exact_target)
end

function _prepare_command_authority_state(
    authority::PreparedCommandAuthority,
    initial_timestamp::PlantTimestamp,
)
    endpoint_count = length(authority.endpoints)
    endpoint_states = Memory{CommandEndpointState}(undef, endpoint_count)
    application_states = Memory{CommandApplicationState}(
        undef, endpoint_count)
    sequences = Memory{UInt64}(undef, endpoint_count)
    fill!(sequences, zero(UInt64))
    @inbounds for index in eachindex(authority.endpoints)
        binding = authority.endpoints[index]
        endpoint = binding.endpoint
        endpoint_state = CommandEndpointState(
            endpoint; initial_timestamp)
        application_state = CommandApplicationState(
            endpoint,
            endpoint_state,
            initial_effective_command(binding);
            safe_command=safe_effective_command(binding),
        )
        endpoint_states[index] = endpoint_state
        application_states[index] = application_state
    end
    state = CommandAuthorityState(
        authority.binding,
        _partition_registry(endpoint_states, CommandEndpointState),
        _partition_registry(application_states, CommandApplicationState),
        sequences,
        zero(UInt64),
        false,
    )
    return _require_exact_command_authority_state_target(
        authority, state, authority.target)
end

function CommandAuthorityState(
    authority::PreparedCommandAuthority;
    initial_timestamp::PlantTimestamp=zero(PlantTimestamp),
)
    return _with_completed_prepared_device_execution_context(
        authority.context) do
        _prepare_command_authority_state(authority, initial_timestamp)
    end
end

function CommandAuthorityWorkspace(authority::PreparedCommandAuthority)
    dispositions = Memory{CommandDispositionWorkspace}(
        undef, length(authority.endpoints))
    @inbounds for index in eachindex(authority.endpoints)
        dispositions[index] = CommandDispositionWorkspace(
            authority.endpoints[index].endpoint)
    end
    return CommandAuthorityWorkspace(
        authority.binding,
        _partition_registry(dispositions, CommandDispositionWorkspace),
    )
end

@inline function _require_command_authority_binding(
    authority::PreparedCommandAuthority,
    state::CommandAuthorityState,
)
    authority.binding === state.binding ||
        _command_authority_preparation_error(
            :foreign_state,
            "command-authority state belongs to another prepared authority",
        )
    return nothing
end

@inline function _require_command_authority_binding(
    authority::PreparedCommandAuthority,
    workspace::CommandAuthorityWorkspace,
)
    authority.binding === workspace.binding ||
        _command_authority_preparation_error(
            :foreign_workspace,
            "command-authority workspace belongs to another prepared authority",
        )
    return nothing
end

function _command_authority_endpoint_slot(
    authority::PreparedCommandAuthority,
    id::CommandEndpointID,
)
    @inbounds for index in eachindex(authority.endpoints)
        command_endpoint_id(authority.endpoints[index]) == id && return index
    end
    _command_authority_preparation_error(
        :unknown_endpoint,
        "prepared command authority has no endpoint $id",
    )
end

@inline _command_authority_endpoint_slot(
    authority::PreparedCommandAuthority,
    name::Symbol,
) = _command_authority_endpoint_slot(authority, CommandEndpointID(name))

function _command_authority_endpoint_slot(
    ::PreparedCommandAuthority,
    id,
)
    _command_authority_preparation_error(
        :invalid_endpoint,
        "command-authority endpoint must be a CommandEndpointID or Symbol; " *
        "got $(typeof(id))",
    )
end

function prepared_command_authority_endpoint(
    authority::PreparedCommandAuthority,
    id,
)
    slot = _command_authority_endpoint_slot(authority, id)
    return @inbounds authority.endpoints[slot].endpoint
end

function command_authority_endpoint_state(
    authority::PreparedCommandAuthority,
    state::CommandAuthorityState,
    id,
)
    _require_command_authority_binding(authority, state)
    slot = _command_authority_endpoint_slot(authority, id)
    return @inbounds state.endpoint_states[slot]
end

function command_authority_application_state(
    authority::PreparedCommandAuthority,
    state::CommandAuthorityState,
    id,
)
    _require_command_authority_binding(authority, state)
    slot = _command_authority_endpoint_slot(authority, id)
    return @inbounds state.application_states[slot]
end

function command_authority_disposition_workspace(
    authority::PreparedCommandAuthority,
    workspace::CommandAuthorityWorkspace,
    id,
)
    _require_command_authority_binding(authority, workspace)
    slot = _command_authority_endpoint_slot(authority, id)
    return @inbounds workspace.dispositions[slot]
end

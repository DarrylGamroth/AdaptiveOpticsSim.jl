#
# Effective-command publications and role-neutral target-local optic owners
#
# This layer owns no command admission, authority placement, payload lease, or
# cross-domain transfer.  A later composition layer binds the sole publisher
# to CommandAuthorityIdentity and supplies already target-local payloads.
#

mutable struct _CommandAuthorityIdentityToken end
const _COMMAND_AUTHORITY_IDENTITY_TOKEN = _CommandAuthorityIdentityToken()

"""Unique run-local identity of the future sole effective-command publisher."""
mutable struct CommandAuthorityIdentity
    function CommandAuthorityIdentity(token::_CommandAuthorityIdentityToken)
        token === _COMMAND_AUTHORITY_IDENTITY_TOKEN || throw(
            ArgumentError("invalid internal command-authority token"))
        return new()
    end
end

mutable struct _EffectiveCommandPublicationSequenceToken end
const _EFFECTIVE_COMMAND_PUBLICATION_SEQUENCE_TOKEN =
    _EffectiveCommandPublicationSequenceToken()

"""Positive endpoint-local order of an effective-command publication."""
struct EffectiveCommandPublicationSequence
    value::UInt64

    function EffectiveCommandPublicationSequence(
        value::UInt64,
        token::_EffectiveCommandPublicationSequenceToken,
    )
        token === _EFFECTIVE_COMMAND_PUBLICATION_SEQUENCE_TOKEN || throw(
            ArgumentError(
                "invalid internal effective-command sequence token"))
        return new(value)
    end
end

@inline function _effective_command_publication_error(
    reason::Symbol,
    message::AbstractString,
)
    throw(PlantCommandError(
        :effective_command_publication, reason, String(message)))
end

function EffectiveCommandPublicationSequence(value::Integer)
    value > 0 || _effective_command_publication_error(
        :invalid_sequence,
        "effective-command publication sequence must be positive",
    )
    value <= typemax(UInt64) || _effective_command_publication_error(
        :invalid_sequence,
        "effective-command publication sequence exceeds UInt64",
    )
    return EffectiveCommandPublicationSequence(
        UInt64(value), _EFFECTIVE_COMMAND_PUBLICATION_SEQUENCE_TOKEN)
end

EffectiveCommandPublicationSequence(::Bool) =
    _effective_command_publication_error(
        :invalid_sequence,
        "effective-command publication sequence must be an integer, not Bool",
    )

function EffectiveCommandPublicationSequence(value)
    _effective_command_publication_error(
        :invalid_sequence,
        "effective-command publication sequence must be an Integer; got " *
        "$(typeof(value))",
    )
end

Base.:(==)(left::EffectiveCommandPublicationSequence,
    right::EffectiveCommandPublicationSequence) = left.value == right.value
Base.isequal(left::EffectiveCommandPublicationSequence,
    right::EffectiveCommandPublicationSequence) =
    isequal(left.value, right.value)
Base.hash(value::EffectiveCommandPublicationSequence, seed::UInt) =
    hash(value.value, hash(EffectiveCommandPublicationSequence, seed))
Base.isless(left::EffectiveCommandPublicationSequence,
    right::EffectiveCommandPublicationSequence) = left.value < right.value

function Base.show(io::IO, value::EffectiveCommandPublicationSequence)
    print(io, "EffectiveCommandPublicationSequence(", value.value, ")")
end

mutable struct _EffectiveCommandPublicationToken end
const _EFFECTIVE_COMMAND_PUBLICATION_TOKEN =
    _EffectiveCommandPublicationToken()

"""
Immutable metadata for one already-effective command publication.

The payload is deliberately separate.  Target-local active and staging
storage own its semantic value; Gate 9A mixed composition owns bounded payload
slots and transfer lifetime.
"""
struct EffectiveCommandPublication{
    I<:CommandAuthorityIdentity,
    D<:AbstractComputeDevice,
    S<:PlantCommandSchema,
}
    authority_identity::I
    target::D
    optic::ControllableOpticID
    schema::S
    timestamp::PlantTimestamp
    sequence::EffectiveCommandPublicationSequence

    function EffectiveCommandPublication(
        token::_EffectiveCommandPublicationToken,
        authority_identity::I,
        target::D,
        optic::ControllableOpticID,
        schema::S,
        timestamp::PlantTimestamp,
        sequence::EffectiveCommandPublicationSequence,
    ) where {
        I<:CommandAuthorityIdentity,
        D<:AbstractComputeDevice,
        S<:PlantCommandSchema,
    }
        token === _EFFECTIVE_COMMAND_PUBLICATION_TOKEN || throw(
            ArgumentError("invalid internal effective-command publication token"))
        return new{I,D,S}(
            authority_identity, target, optic, schema, timestamp, sequence)
    end
end

@inline function _effective_command_publication(
    authority_identity::CommandAuthorityIdentity,
    target::AbstractComputeDevice,
    optic::ControllableOpticID,
    schema::PlantCommandSchema,
    timestamp::PlantTimestamp,
    sequence,
)
    return EffectiveCommandPublication(
        _EFFECTIVE_COMMAND_PUBLICATION_TOKEN,
        authority_identity,
        target,
        optic,
        schema,
        timestamp,
        EffectiveCommandPublicationSequence(sequence),
    )
end

@inline command_authority_identity(
    publication::EffectiveCommandPublication) = publication.authority_identity
@inline compute_device(publication::EffectiveCommandPublication) =
    publication.target
@inline controllable_optic_id(publication::EffectiveCommandPublication) =
    publication.optic
@inline command_schema(publication::EffectiveCommandPublication) =
    publication.schema
@inline command_endpoint_id(publication::EffectiveCommandPublication) =
    command_endpoint_id(command_schema(publication))
@inline effective_command_publication_timestamp(
    publication::EffectiveCommandPublication) = publication.timestamp
@inline effective_command_publication_sequence(
    publication::EffectiveCommandPublication) = publication.sequence

mutable struct _TargetLocalCommandEndpointBinding end

"""Immutable semantic binding for one target-local effective-command endpoint."""
struct PreparedTargetLocalCommandEndpoint{
    I<:CommandAuthorityIdentity,
    D<:AbstractComputeDevice,
    S<:PlantCommandSchema,
}
    binding::_TargetLocalCommandEndpointBinding
    authority_identity::I
    target::D
    optic::ControllableOpticID
    schema::S
end

abstract type _AbstractTargetLocalEffectiveCommandValues end

mutable struct _TargetLocalScalarEffectiveCommandValues{T} <:
    _AbstractTargetLocalEffectiveCommandValues
    active::T
    staging::T
end

mutable struct _TargetLocalArrayEffectiveCommandValues{A<:AbstractArray} <:
    _AbstractTargetLocalEffectiveCommandValues
    active::A
    staging::A
end

"""Single-writer state for one role-neutral target-local command endpoint."""
mutable struct TargetLocalCommandEndpointState{
    V<:_AbstractTargetLocalEffectiveCommandValues,
}
    binding::_TargetLocalCommandEndpointBinding
    values::V
    last_timestamp::PlantTimestamp
    staged_timestamp::PlantTimestamp
    last_sequence::UInt64
    staged_sequence::UInt64
    has_publication::Bool
    has_staged_publication::Bool
end

mutable struct _TargetLocalControllableOpticBinding end

"""Run-immutable target-local physical-optic preparation and endpoint table."""
struct PreparedTargetLocalControllableOptic{
    D<:AbstractComputeDevice,
    C<:ControllableOpticDefinition,
    I,
}
    binding::_TargetLocalControllableOpticBinding
    target::D
    definition::C
    implementation::I
    endpoints::_FixedPlantRegistry{PreparedTargetLocalCommandEndpoint}
end

"""Mutable physical and effective-command state for one target-local optic."""
mutable struct TargetLocalControllableOpticState{S}
    binding::_TargetLocalControllableOpticBinding
    physical::S
    endpoints::_FixedPlantRegistry{TargetLocalCommandEndpointState}
    staged_endpoint_slot::UInt32
    has_staged_publication::Bool
end

"""Reusable physical scratch for one target-local controllable optic."""
struct TargetLocalControllableOpticWorkspace{W}
    binding::_TargetLocalControllableOpticBinding
    physical::W
end

"""Aligned prepared, state, and workspace owners for one target-local optic."""
struct TargetLocalControllableOpticOwner{P,S,W}
    prepared::P
    state::S
    workspace::W
end

"""
Cold type-restoring binding for one target-local endpoint owner.

Cardinality-neutral partition registries deliberately erase heterogeneous
endpoint types. Mixed composition prepares and retains this binding once so
the warmed local stage/commit path dispatches on concrete endpoint and state
types without searching the registry.
"""
mutable struct _TargetLocalCommandEndpointOwnerToken end
const _TARGET_LOCAL_COMMAND_ENDPOINT_OWNER_TOKEN =
    _TargetLocalCommandEndpointOwnerToken()

struct TargetLocalCommandEndpointOwner{O,E,S}
    optic::O
    endpoint::E
    state::S
    slot::UInt32

    function TargetLocalCommandEndpointOwner(
        token::_TargetLocalCommandEndpointOwnerToken,
        optic::O,
        endpoint::E,
        state::S,
        slot::UInt32,
    ) where {O,E,S}
        token === _TARGET_LOCAL_COMMAND_ENDPOINT_OWNER_TOKEN || throw(
            ArgumentError("invalid internal target-local endpoint-owner token"))
        return new{O,E,S}(optic, endpoint, state, slot)
    end
end

@inline compute_device(endpoint::PreparedTargetLocalCommandEndpoint) =
    endpoint.target
@inline controllable_optic_id(
    endpoint::PreparedTargetLocalCommandEndpoint) = endpoint.optic
@inline command_schema(endpoint::PreparedTargetLocalCommandEndpoint) =
    endpoint.schema
@inline command_endpoint_id(endpoint::PreparedTargetLocalCommandEndpoint) =
    command_endpoint_id(command_schema(endpoint))
@inline command_authority_identity(
    endpoint::PreparedTargetLocalCommandEndpoint) = endpoint.authority_identity

@inline compute_device(optic::PreparedTargetLocalControllableOptic) =
    optic.target
@inline controllable_optic_id(
    optic::PreparedTargetLocalControllableOptic) =
    controllable_optic_id(optic.definition)
@inline target_local_command_endpoints(
    optic::PreparedTargetLocalControllableOptic) = optic.endpoints
@inline prepared_target_local_controllable_optic(
    owner::TargetLocalControllableOpticOwner) = owner.prepared
@inline target_local_controllable_optic_state(
    owner::TargetLocalControllableOpticOwner) = owner.state
@inline target_local_controllable_optic_workspace(
    owner::TargetLocalControllableOpticOwner) = owner.workspace
@inline compute_device(owner::TargetLocalControllableOpticOwner) =
    compute_device(prepared_target_local_controllable_optic(owner))
@inline controllable_optic_id(owner::TargetLocalControllableOpticOwner) =
    controllable_optic_id(prepared_target_local_controllable_optic(owner))
@inline command_endpoint_id(owner::TargetLocalCommandEndpointOwner) =
    command_endpoint_id(owner.endpoint)
@inline controllable_optic_id(owner::TargetLocalCommandEndpointOwner) =
    controllable_optic_id(owner.endpoint)
@inline compute_device(owner::TargetLocalCommandEndpointOwner) =
    compute_device(owner.endpoint)

@inline _active_effective_command(
    values::_TargetLocalScalarEffectiveCommandValues) = values.active
@inline _active_effective_command(
    values::_TargetLocalArrayEffectiveCommandValues) = values.active
@inline effective_command(state::TargetLocalCommandEndpointState) =
    _active_effective_command(state.values)
@inline has_effective_command_publication(
    state::TargetLocalCommandEndpointState) = state.has_publication

@inline function last_effective_command_publication_timestamp(
    state::TargetLocalCommandEndpointState,
)
    state.has_publication || return nothing
    return state.last_timestamp
end

@inline function last_effective_command_publication_sequence(
    state::TargetLocalCommandEndpointState,
)
    state.has_publication || return nothing
    return EffectiveCommandPublicationSequence(
        state.last_sequence, _EFFECTIVE_COMMAND_PUBLICATION_SEQUENCE_TOKEN)
end

function target_local_command_endpoint_state(
    optic::PreparedTargetLocalControllableOptic,
    state::TargetLocalControllableOpticState,
    endpoint::CommandEndpointID,
)
    optic.binding === state.binding || _effective_command_publication_error(
        :foreign_optic_state,
        "target-local controllable-optic state belongs to another owner",
    )
    for index in eachindex(optic.endpoints)
        prepared_endpoint = optic.endpoints[index]
        command_endpoint_id(prepared_endpoint) == endpoint || continue
        endpoint_state = state.endpoints[index]
        prepared_endpoint.binding === endpoint_state.binding ||
            _effective_command_publication_error(
                :foreign_endpoint_state,
                "target-local command state belongs to another prepared endpoint",
            )
        return endpoint_state
    end
    _effective_command_publication_error(
        :foreign_endpoint,
        "target-local controllable optic $(controllable_optic_id(optic)) " *
        "does not own endpoint $endpoint",
    )
end

@inline function target_local_command_endpoint_state(
    owner::TargetLocalControllableOpticOwner,
    endpoint::CommandEndpointID,
)
    return target_local_command_endpoint_state(
        owner.prepared, owner.state, endpoint)
end

function _target_local_command_endpoint_slot(
    optic::PreparedTargetLocalControllableOptic,
    endpoint::CommandEndpointID,
)
    @inbounds for index in eachindex(optic.endpoints)
        command_endpoint_id(optic.endpoints[index]) == endpoint && return index
    end
    _effective_command_publication_error(
        :foreign_endpoint,
        "target-local controllable optic $(controllable_optic_id(optic)) " *
        "does not own endpoint $endpoint",
    )
end

function target_local_command_endpoint_owner(
    owner::TargetLocalControllableOpticOwner,
    endpoint_id::CommandEndpointID,
)
    slot = _target_local_command_endpoint_slot(owner.prepared, endpoint_id)
    endpoint = @inbounds owner.prepared.endpoints[slot]
    state = @inbounds owner.state.endpoints[slot]
    endpoint.binding === state.binding ||
        _effective_command_publication_error(
            :foreign_endpoint_state,
            "target-local command state belongs to another prepared endpoint",
        )
    return TargetLocalCommandEndpointOwner(
        _TARGET_LOCAL_COMMAND_ENDPOINT_OWNER_TOKEN,
        owner,
        endpoint,
        state,
        UInt32(slot),
    )
end

@inline function _require_effective_value_bounds(
    ::PlantCommandSchema,
    value,
    ::UnboundedCommandValues,
)
    return value
end

function _require_effective_value_bounds(
    schema::PlantCommandSchema,
    value,
    bounds::UniformCommandBounds,
)
    _effective_seed_uses_command_bounds(schema) || return value
    _all_command_values_in_bounds(value, bounds) ||
        _effective_command_publication_error(
            :out_of_range,
            "effective value lies outside the endpoint's command bounds",
        )
    return value
end

function _require_target_local_effective_value_layout(
    schema::PlantCommandSchema{T,0},
    value,
    ::AbstractComputeDevice,
) where {T}
    typeof(value) === T || _effective_command_publication_error(
        :numeric_type,
        "scalar effective value has type $(typeof(value)); expected $T",
    )
    return value
end

function _require_target_local_effective_value(
    schema::PlantCommandSchema{T,0},
    value,
    target::AbstractComputeDevice,
) where {T}
    _require_target_local_effective_value_layout(schema, value, target)
    isfinite(value) || _effective_command_publication_error(
        :nonfinite,
        "effective value must be finite",
    )
    return _require_effective_value_bounds(
        schema, value, command_bounds(schema))
end

function _require_target_local_effective_value_layout(
    ::PlantCommandSchema{T,0},
    value::AbstractArray{T,0},
    ::AbstractComputeDevice,
) where {T}
    _effective_command_publication_error(
        :numeric_type,
        "scalar effective value must have exact scalar type $T, not an array",
    )
end

function _require_target_local_effective_value_layout(
    schema::PlantCommandSchema{T,N},
    value::AbstractArray{T,N},
    target::AbstractComputeDevice,
) where {T,N}
    size(value) == command_dimensions(schema) ||
        _effective_command_publication_error(
            :dimensions,
            "effective value has size $(size(value)); expected " *
            "$(command_dimensions(schema))",
        )
    compute_device(value) == target || _effective_command_publication_error(
        :wrong_target,
        "effective value occupies $(compute_device(value)); expected $target",
    )
    return value
end

function _require_target_local_effective_value(
    schema::PlantCommandSchema{T,N},
    value::AbstractArray{T,N},
    target::AbstractComputeDevice,
) where {T,N}
    _require_target_local_effective_value_layout(schema, value, target)
    _all_command_values_finite(value) ||
        _effective_command_publication_error(
            :nonfinite,
            "effective value must contain only finite values",
        )
    return _require_effective_value_bounds(
        schema, value, command_bounds(schema))
end

function _require_target_local_effective_value_layout(
    schema::PlantCommandSchema{T,N},
    value,
    ::AbstractComputeDevice,
) where {T,N}
    _effective_command_publication_error(
        :numeric_type,
        "effective value must be an AbstractArray{$T,$N}; got " *
        "$(typeof(value))",
    )
end

@inline function _require_target_local_effective_value(
    schema::PlantCommandSchema,
    value,
    target::AbstractComputeDevice,
)
    return _require_target_local_effective_value_layout(schema, value, target)
end

@inline function _require_publication_schema_field(
    matches::Bool,
    reason::Symbol,
    label::AbstractString,
)
    matches || _effective_command_publication_error(
        reason, "effective-command publication $label does not match " *
        "the prepared endpoint")
    return nothing
end

@inline function _effective_command_publication_schema_matches(
    expected::PlantCommandSchema,
    supplied::PlantCommandSchema,
)
    return command_schema_id(supplied) == command_schema_id(expected) &&
        command_schema_version(supplied) ==
            command_schema_version(expected) &&
        command_endpoint_id(supplied) == command_endpoint_id(expected) &&
        command_numeric_type(supplied) === command_numeric_type(expected) &&
        command_dimensions(supplied) == command_dimensions(expected) &&
        command_units(supplied) == command_units(expected) &&
        command_sign_convention(supplied) ==
            command_sign_convention(expected) &&
        command_basis(supplied) == command_basis(expected) &&
        command_basis_revision(supplied) ==
            command_basis_revision(expected) &&
        command_semantics(supplied) == command_semantics(expected) &&
        command_bounds(supplied) == command_bounds(expected) &&
        command_value_policy(supplied) == command_value_policy(expected) &&
        command_sequence_policy(supplied) ==
            command_sequence_policy(expected) &&
        command_effective_time_policy(supplied) ==
            command_effective_time_policy(expected) &&
        command_silence_policy(supplied) ==
            command_silence_policy(expected)
end

function _require_effective_command_publication_metadata(
    endpoint::PreparedTargetLocalCommandEndpoint,
    state::TargetLocalCommandEndpointState,
    publication::EffectiveCommandPublication,
)
    endpoint.binding === state.binding || _effective_command_publication_error(
        :foreign_endpoint_state,
        "target-local command state belongs to another prepared endpoint",
    )
    command_authority_identity(publication) ===
        command_authority_identity(endpoint) ||
        _effective_command_publication_error(
            :foreign_authority,
            "effective-command publication belongs to another authority",
        )
    compute_device(publication) == compute_device(endpoint) ||
        _effective_command_publication_error(
            :wrong_target,
            "effective-command publication targets $(compute_device(publication)); " *
            "expected $(compute_device(endpoint))",
        )
    controllable_optic_id(publication) == controllable_optic_id(endpoint) ||
        _effective_command_publication_error(
            :foreign_optic,
            "effective-command publication names another controllable optic",
        )
    command_endpoint_id(publication) == command_endpoint_id(endpoint) ||
        _effective_command_publication_error(
            :foreign_endpoint,
            "effective-command publication names another endpoint",
        )

    expected = command_schema(endpoint)
    supplied = command_schema(publication)
    _require_publication_schema_field(
        command_schema_id(supplied) == command_schema_id(expected),
        :schema_id, "schema identity")
    _require_publication_schema_field(
        command_schema_version(supplied) == command_schema_version(expected),
        :schema_version, "schema version")
    _require_publication_schema_field(
        command_numeric_type(supplied) === command_numeric_type(expected),
        :numeric_type, "numeric type")
    _require_publication_schema_field(
        command_dimensions(supplied) == command_dimensions(expected),
        :dimensions, "dimensions")
    _require_publication_schema_field(
        command_units(supplied) == command_units(expected),
        :units, "units")
    _require_publication_schema_field(
        command_sign_convention(supplied) ==
            command_sign_convention(expected),
        :sign_convention, "sign convention")
    _require_publication_schema_field(
        command_basis(supplied) == command_basis(expected),
        :basis, "basis")
    _require_publication_schema_field(
        command_basis_revision(supplied) ==
            command_basis_revision(expected),
        :basis_revision, "basis revision")
    _require_publication_schema_field(
        command_semantics(supplied) == command_semantics(expected),
        :semantics, "semantics")
    _require_publication_schema_field(
        command_bounds(supplied) == command_bounds(expected),
        :bounds, "bounds")
    _require_publication_schema_field(
        command_value_policy(supplied) == command_value_policy(expected),
        :value_policy, "value policy")
    _require_publication_schema_field(
        command_sequence_policy(supplied) ==
            command_sequence_policy(expected),
        :sequence_policy, "sequence policy")
    _require_publication_schema_field(
        command_effective_time_policy(supplied) ==
            command_effective_time_policy(expected),
        :effective_time_policy, "effective-time policy")
    _require_publication_schema_field(
        command_silence_policy(supplied) ==
            command_silence_policy(expected),
        :silence_policy, "silence policy")

    state.has_staged_publication && _effective_command_publication_error(
        :stage_pending,
        "target-local endpoint already has a staged publication",
    )
    sequence = effective_command_publication_sequence(publication).value
    timestamp = effective_command_publication_timestamp(publication)
    if state.has_publication
        sequence > state.last_sequence ||
            _effective_command_publication_error(
                sequence == state.last_sequence ?
                    :duplicate_sequence : :stale_sequence,
                "effective-command publication sequence must increase",
            )
        timestamp < state.last_timestamp &&
            _effective_command_publication_error(
                :regressing_timestamp,
                "effective-command publication timestamp must not regress",
            )
    end
    return nothing
end

function _require_effective_command_publication(
    endpoint::PreparedTargetLocalCommandEndpoint,
    state::TargetLocalCommandEndpointState,
    publication::EffectiveCommandPublication,
    value,
)
    _require_effective_command_publication_metadata(
        endpoint, state, publication)
    _require_target_local_effective_value(
        command_schema(endpoint), value, endpoint.target)
    return nothing
end

# Effective values routed from the sole command authority have already passed
# finite-value and effective-bound validation while the authority staged its
# candidate. A target-local replica therefore revalidates publication identity,
# ordering, shape, and exact residency without launching a redundant device
# reduction for every replica.
function _require_routed_effective_command_publication(
    endpoint::PreparedTargetLocalCommandEndpoint,
    state::TargetLocalCommandEndpointState,
    publication::EffectiveCommandPublication,
    value,
)
    _require_effective_command_publication_metadata(
        endpoint, state, publication)
    _require_target_local_effective_value_layout(
        command_schema(endpoint), value, endpoint.target)
    return nothing
end

@inline function _copy_staged_effective_command!(
    values::_TargetLocalScalarEffectiveCommandValues,
    value,
)
    values.staging = value
    return nothing
end

@inline function _copy_staged_effective_command!(
    values::_TargetLocalArrayEffectiveCommandValues,
    value,
)
    copyto!(values.staging, value)
    return nothing
end

@inline function _commit_staged_effective_command!(
    values::_TargetLocalScalarEffectiveCommandValues,
)
    values.active = values.staging
    return nothing
end

@inline function _commit_staged_effective_command!(
    values::_TargetLocalArrayEffectiveCommandValues,
)
    previous = values.active
    values.active = values.staging
    values.staging = previous
    return nothing
end

function _stage_prevalidated_effective_command_publication!(
    endpoint_owner::TargetLocalCommandEndpointOwner,
    publication::EffectiveCommandPublication,
    value,
)
    owner = endpoint_owner.optic
    prepared = owner.prepared
    state = owner.state
    workspace = owner.workspace
    prepared.binding === state.binding === workspace.binding ||
        _effective_command_publication_error(
            :foreign_optic_state,
            "target-local controllable-optic plan, state, and workspace " *
            "belong to different owners",
        )
    endpoint = endpoint_owner.endpoint
    endpoint_state = endpoint_owner.state
    state.has_staged_publication && _effective_command_publication_error(
        :optic_stage_pending,
        "target-local controllable optic already has a staged publication",
    )
    _copy_staged_effective_command!(endpoint_state.values, value)
    stage_controllable_optic_command!(
        prepared.implementation,
        state.physical,
        workspace.physical,
        command_endpoint_id(endpoint),
        _staged_effective_command(endpoint_state.values),
        effective_command_publication_timestamp(publication),
    )
    endpoint_state.staged_timestamp =
        effective_command_publication_timestamp(publication)
    endpoint_state.staged_sequence =
        effective_command_publication_sequence(publication).value
    endpoint_state.has_staged_publication = true
    state.staged_endpoint_slot = endpoint_owner.slot
    state.has_staged_publication = true
    return nothing
end

function _stage_effective_command_publication!(
    endpoint_owner::TargetLocalCommandEndpointOwner,
    publication::EffectiveCommandPublication,
    value,
)
    _require_effective_command_publication(
        endpoint_owner.endpoint,
        endpoint_owner.state,
        publication,
        value,
    )
    return _stage_prevalidated_effective_command_publication!(
        endpoint_owner, publication, value)
end

function _stage_effective_command_publication!(
    owner::TargetLocalControllableOpticOwner,
    publication::EffectiveCommandPublication,
    value,
)
    endpoint_owner = target_local_command_endpoint_owner(
        owner, command_endpoint_id(publication))
    return _stage_effective_command_publication!(
        endpoint_owner, publication, value)
end

@inline _staged_effective_command(
    values::_TargetLocalScalarEffectiveCommandValues) = values.staging
@inline _staged_effective_command(
    values::_TargetLocalArrayEffectiveCommandValues) = values.staging

function _commit_effective_command_publication!(
    endpoint_owner::TargetLocalCommandEndpointOwner,
)
    owner = endpoint_owner.optic
    prepared = owner.prepared
    state = owner.state
    workspace = owner.workspace
    prepared.binding === state.binding === workspace.binding ||
        _effective_command_publication_error(
            :foreign_optic_state,
            "target-local controllable-optic plan, state, and workspace " *
            "belong to different owners",
        )
    endpoint = endpoint_owner.endpoint
    endpoint_state = endpoint_owner.state
    endpoint_id = command_endpoint_id(endpoint)
    endpoint_state.has_staged_publication ||
        _effective_command_publication_error(
            :no_staged_publication,
            "target-local endpoint has no staged publication to commit",
        )
    (state.has_staged_publication &&
        state.staged_endpoint_slot == endpoint_owner.slot) ||
        _effective_command_publication_error(
            :foreign_staged_endpoint,
            "target-local physical staging belongs to another endpoint",
        )
    commit_controllable_optic_command!(
        prepared.implementation,
        state.physical,
        workspace.physical,
        endpoint_id,
        endpoint_state.staged_timestamp,
    )
    _commit_staged_effective_command!(endpoint_state.values)
    endpoint_state.last_timestamp = endpoint_state.staged_timestamp
    endpoint_state.last_sequence = endpoint_state.staged_sequence
    endpoint_state.has_publication = true
    endpoint_state.has_staged_publication = false
    state.staged_endpoint_slot = UInt32(0)
    state.has_staged_publication = false
    return nothing
end


function _commit_effective_command_publication!(
    owner::TargetLocalControllableOpticOwner,
    endpoint_id::CommandEndpointID,
)
    endpoint_owner = target_local_command_endpoint_owner(owner, endpoint_id)
    return _commit_effective_command_publication!(endpoint_owner)
end

function _prepare_target_local_effective_values(
    schema::PlantCommandSchema{T,0},
    initial,
    ::AbstractComputeDevice,
) where {T}
    value = _validate_effective_seed(
        schema, initial, "initial target-local effective command")
    return _TargetLocalScalarEffectiveCommandValues(value, value)
end

function _prepare_target_local_effective_values(
    schema::PlantCommandSchema{T,N},
    initial,
    target::AbstractComputeDevice,
) where {T,N}
    validated = _validate_effective_seed(
        schema, initial, "initial target-local effective command")
    active = allocate_device_array(target, T, command_dimensions(schema)...)
    staging = allocate_device_array(target, T, command_dimensions(schema)...)
    copyto!(active, validated)
    copyto!(staging, validated)
    active === staging && throw(PlantPreparationError(
        :command_replica,
        :aliased_staging,
        "target-local active and staging command storage must be distinct",
    ))
    return _TargetLocalArrayEffectiveCommandValues(active, staging)
end

function _copy_target_local_physical_initial(
    schema::PlantCommandSchema{T,0},
    initial,
    ::AbstractComputeDevice,
) where {T}
    return _validate_effective_seed(
        schema, initial, "initial target-local physical command")
end

function _copy_target_local_physical_initial(
    schema::PlantCommandSchema{T,N},
    initial,
    target::AbstractComputeDevice,
) where {T,N}
    validated = _validate_effective_seed(
        schema, initial, "initial target-local physical command")
    copied = allocate_device_array(
        target, T, command_dimensions(schema)...)
    copyto!(copied, validated)
    return copied
end

function _canonical_target_local_command_schemas(
    definition::ControllableOpticDefinition,
)
    schemas = collect(command_schemas(definition))
    sort!(schemas; by=schema -> String(command_endpoint_id(schema).name))
    return _partition_registry(schemas, PlantCommandSchema)
end

function _partition_command_endpoint_configurations(
    definition::PlantDefinition,
    configurations,
)
    supplied = _sorted_command_endpoint_configurations(configurations)
    declared = _canonical_command_endpoint_declarations(definition)
    length(supplied) == length(declared) || throw(PlantPreparationError(
        :command_replica,
        :configuration_count,
        "plant declares $(length(declared)) command endpoints but " *
        "$(length(supplied)) target-local configurations were supplied",
    ))
    prepared = Memory{CommandEndpointConfiguration}(
        undef, length(supplied))
    @inbounds for index in eachindex(declared)
        schema = declared[index][2]
        expected = command_endpoint_id(schema)
        configuration = supplied[index]
        actual = command_endpoint_id(configuration)
        actual == expected || throw(PlantPreparationError(
            :command_replica,
            :configuration_identity,
            "target-local command configuration $actual does not match " *
            "declared endpoint $expected",
        ))
        prepared[index] = _seal_partition_command_endpoint_configuration(
            schema, configuration)
    end
    return prepared
end

@inline function _seal_partition_command_value(
    schema::PlantCommandSchema{T,0},
    value,
    label::AbstractString,
) where {T}
    return _validate_effective_seed(schema, value, label)
end

function _seal_partition_command_value(
    schema::PlantCommandSchema{T,N},
    value,
    label::AbstractString,
) where {T,N}
    validated = _validate_effective_seed(schema, value, label)
    compute_device(validated) == HostComputeDevice() || throw(
        PlantPreparationError(
            :command_replica,
            :configuration_residency,
            "$label for $(command_endpoint_id(schema)) must be host-resident " *
            "cold configuration; got $(compute_device(validated))",
        ))
    sealed = copy(validated)
    compute_device(sealed) == HostComputeDevice() || throw(
        PlantPreparationError(
            :command_replica,
            :configuration_copy,
            "copied $label for $(command_endpoint_id(schema)) is not " *
            "host-resident",
        ))
    Base.mightalias(sealed, validated) && throw(PlantPreparationError(
        :command_replica,
        :aliased_configuration,
        "copied $label for $(command_endpoint_id(schema)) still aliases " *
        "caller-owned storage",
    ))
    return sealed
end

function _seal_partition_command_endpoint_configuration(
    schema::PlantCommandSchema,
    configuration::CommandEndpointConfiguration,
)
    initial = _seal_partition_command_value(
        schema,
        initial_effective_command(configuration),
        "initial effective command",
    )
    supplied_safe = safe_effective_command(configuration)
    _require_safe_command_configuration(
        command_silence_policy(schema), supplied_safe)
    safe = supplied_safe === nothing ? nothing :
        _seal_partition_command_value(schema, supplied_safe, "safe command")
    return CommandEndpointConfiguration(
        command_endpoint_id(configuration),
        command_endpoint_capacity(configuration),
        command_sequence_window(configuration),
        initial,
        safe,
    )
end

function _prepare_target_local_command_endpoints(
    definition::ControllableOpticDefinition,
    configurations,
    authority_identity::CommandAuthorityIdentity,
    target::AbstractComputeDevice,
)
    schemas = _canonical_target_local_command_schemas(definition)
    endpoints = Memory{PreparedTargetLocalCommandEndpoint}(
        undef, length(schemas))
    states = Memory{TargetLocalCommandEndpointState}(
        undef, length(schemas))
    physical_initials = Memory{Union{Real,AbstractArray}}(
        undef, length(schemas))
    origin = zero(PlantTimestamp)
    optic_id = controllable_optic_id(definition)
    @inbounds for index in eachindex(schemas)
        schema = schemas[index]
        configuration = _command_endpoint_configuration(
            configurations, command_endpoint_id(schema))
        binding = _TargetLocalCommandEndpointBinding()
        endpoints[index] = PreparedTargetLocalCommandEndpoint(
            binding, authority_identity, target, optic_id, schema)
        values = _prepare_target_local_effective_values(
            schema, initial_effective_command(configuration), target)
        states[index] = TargetLocalCommandEndpointState(
            binding,
            values,
            origin,
            origin,
            UInt64(0),
            UInt64(0),
            false,
            false,
        )
        validate_target_local_command_endpoint_target(
            endpoints[index], states[index], target)
        physical_initials[index] = _copy_target_local_physical_initial(
            schema, initial_effective_command(configuration), target)
    end
    return (
        _partition_registry(endpoints, PreparedTargetLocalCommandEndpoint),
        _partition_registry(states, TargetLocalCommandEndpointState),
        Tuple(physical_initials),
    )
end

function _prepare_target_local_controllable_optic_owner(
    definition::ControllableOpticDefinition,
    telescope::AbstractTelescope,
    atmosphere_definition::AbstractTimedAtmosphereDefinition,
    target::AbstractComputeDevice,
    configurations,
    authority_identity::CommandAuthorityIdentity,
)
    implementation = prepare_target_local_controllable_optic(
        controllable_optic_model(definition),
        definition,
        telescope,
        atmosphere_definition,
        target,
    )
    implementation === nothing && throw(PlantPreparationError(
        :controllable_optic,
        :invalid_target_local_preparation,
        "prepare_target_local_controllable_optic returned nothing for " *
        "$(controllable_optic_id(definition))",
    ))
    ismutabletype(typeof(implementation)) && throw(PlantPreparationError(
        :controllable_optic,
        :mutable_target_local_preparation,
        "prepare_target_local_controllable_optic must return immutable " *
        "preparation data",
    ))
    validate_controllable_optic_target(implementation, target)

    endpoints, endpoint_states, physical_initials =
        _prepare_target_local_command_endpoints(
            definition, configurations, authority_identity, target)
    physical_state = prepare_controllable_optic_state(
        implementation,
        definition,
        Tuple(command_endpoint_id(endpoint) for endpoint in endpoints),
        physical_initials,
    )
    physical_workspace = prepare_controllable_optic_workspace(implementation)
    validate_controllable_optic_state_target(
        implementation, physical_state, target)
    validate_controllable_optic_workspace_target(
        implementation, physical_workspace, target)

    binding = _TargetLocalControllableOpticBinding()
    prepared = PreparedTargetLocalControllableOptic(
        binding, target, definition, implementation, endpoints)
    state = TargetLocalControllableOpticState(
        binding, physical_state, endpoint_states, UInt32(0), false)
    workspace = TargetLocalControllableOpticWorkspace(
        binding, physical_workspace)
    return TargetLocalControllableOpticOwner(prepared, state, workspace)
end

function _prepare_target_local_controllable_optic_owners(
    definitions::AbstractVector,
    telescope::AbstractTelescope,
    atmosphere_definition::AbstractTimedAtmosphereDefinition,
    target::AbstractComputeDevice,
    configurations,
    authority_identity::CommandAuthorityIdentity,
)
    owners = Memory{TargetLocalControllableOpticOwner}(
        undef, length(definitions))
    @inbounds for index in eachindex(definitions)
        owners[index] = _prepare_target_local_controllable_optic_owner(
            definitions[index],
            telescope,
            atmosphere_definition,
            target,
            configurations,
            authority_identity,
        )
    end
    return _partition_registry(owners, TargetLocalControllableOpticOwner)
end

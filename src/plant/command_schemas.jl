#
# Immutable plant-command schemas
#
# This layer defines semantic payload meaning and validation only. Mutable
# admission, effective-time calendars, application, held state, and terminal
# model dispositions are deliberately separate later layers.
#

@inline function _command_schema_error(reason::Symbol, message::AbstractString)
    throw(PlantCommandError(:schema, reason, String(message)))
end

@inline function _require_command_name(name::Symbol, reason::Symbol,
    label::AbstractString)
    isempty(String(name)) && _command_schema_error(reason,
        "$label must not be empty")
    return name
end

"""Stable declared identity of one semantic plant-command schema."""
struct PlantCommandSchemaID
    name::Symbol

    function PlantCommandSchemaID(name::Symbol)
        return new(_require_command_name(name, :empty_id,
            "plant-command schema identity"))
    end
end

PlantCommandSchemaID(value) = _command_schema_error(:invalid_id,
    "plant-command schema identity must be a Symbol; got $(typeof(value))")

struct _CommandRevisionToken end
const _COMMAND_REVISION_TOKEN = _CommandRevisionToken()

@inline function _checked_command_revision(value::Integer, reason::Symbol,
    label::AbstractString)
    value > 0 || _command_schema_error(reason,
        "$label must be positive")
    value <= typemax(UInt32) || _command_schema_error(reason,
        "$label exceeds UInt32 range")
    return UInt32(value)
end

@inline _checked_command_revision(::Bool, reason::Symbol,
    label::AbstractString) = _command_schema_error(reason,
    "$label must be a positive integer, not Bool")

"""Positive, process-stable version of a plant-command schema."""
struct PlantCommandSchemaVersion
    value::UInt32

    PlantCommandSchemaVersion(value::UInt32, ::_CommandRevisionToken) =
        new(value)
end

PlantCommandSchemaVersion(value::Integer) = PlantCommandSchemaVersion(
    _checked_command_revision(value, :invalid_version,
        "plant-command schema version"),
    _COMMAND_REVISION_TOKEN,
)

PlantCommandSchemaVersion(value) = _command_schema_error(:invalid_version,
    "plant-command schema version must be an Integer; got $(typeof(value))")

"""Positive revision identifying the basis/calibration interpretation."""
struct CommandBasisRevision
    value::UInt32

    CommandBasisRevision(value::UInt32, ::_CommandRevisionToken) = new(value)
end

CommandBasisRevision(value::Integer) = CommandBasisRevision(
    _checked_command_revision(value, :invalid_basis_revision,
        "command basis revision"),
    _COMMAND_REVISION_TOKEN,
)

CommandBasisRevision(value) = _command_schema_error(
    :invalid_basis_revision,
    "command basis revision must be an Integer; got $(typeof(value))")

"""Explicit physical unit identity for command values."""
struct CommandUnit
    name::Symbol

    function CommandUnit(name::Symbol)
        return new(_require_command_name(name, :invalid_unit,
            "command unit"))
    end
end

CommandUnit(value) = _command_schema_error(:invalid_unit,
    "command unit must be a Symbol; got $(typeof(value))")

"""Explicit physical sign-convention identity for command values."""
struct CommandSignConvention
    name::Symbol

    function CommandSignConvention(name::Symbol)
        return new(_require_command_name(name, :invalid_sign_convention,
            "command sign convention"))
    end
end

CommandSignConvention(value) = _command_schema_error(
    :invalid_sign_convention,
    "command sign convention must be a Symbol; got $(typeof(value))")

"""
    CommandBasis(kind, name)

Stable identity of the coordinate basis in which a command payload is
expressed. `kind` may describe an actuator, modal, rigid-body, or
application-defined basis; `name` identifies the particular basis.
"""
struct CommandBasis
    kind::Symbol
    name::Symbol

    function CommandBasis(kind::Symbol, name::Symbol)
        return new(
            _require_command_name(kind, :invalid_basis, "command basis kind"),
            _require_command_name(name, :invalid_basis, "command basis name"),
        )
    end
end


function CommandBasis(kind, name)
    _command_schema_error(:invalid_basis,
        "command basis kind and name must be Symbols; got $(typeof(kind)) " *
        "and $(typeof(name))")
end

const _NamedCommandSchemaValue = Union{
    PlantCommandSchemaID,
    CommandUnit,
    CommandSignConvention,
}

Base.:(==)(left::T, right::T) where {T<:_NamedCommandSchemaValue} =
    left.name == right.name
Base.isequal(left::T, right::T) where {T<:_NamedCommandSchemaValue} =
    isequal(left.name, right.name)
Base.hash(value::T, seed::UInt) where {T<:_NamedCommandSchemaValue} =
    hash(value.name, hash(T, seed))

const _CommandRevision = Union{
    PlantCommandSchemaVersion,
    CommandBasisRevision,
}

Base.:(==)(left::T, right::T) where {T<:_CommandRevision} =
    left.value == right.value
Base.isequal(left::T, right::T) where {T<:_CommandRevision} =
    isequal(left.value, right.value)
Base.hash(value::T, seed::UInt) where {T<:_CommandRevision} =
    hash(value.value, hash(T, seed))

Base.:(==)(left::CommandBasis, right::CommandBasis) =
    left.kind == right.kind && left.name == right.name
Base.isequal(left::CommandBasis, right::CommandBasis) =
    isequal(left.kind, right.kind) && isequal(left.name, right.name)
Base.hash(value::CommandBasis, seed::UInt) =
    hash(value.name, hash(value.kind, hash(CommandBasis, seed)))

function Base.show(io::IO, value::_NamedCommandSchemaValue)
    print(io, nameof(typeof(value)), "(", repr(value.name), ")")
end

function Base.show(io::IO, value::_CommandRevision)
    print(io, nameof(typeof(value)), "(", value.value, ")")
end

function Base.show(io::IO, value::CommandBasis)
    print(io, "CommandBasis(", repr(value.kind), ", ", repr(value.name), ")")
end

@inline _as_plant_command_schema_id(id::PlantCommandSchemaID) = id
@inline _as_plant_command_schema_id(name::Symbol) = PlantCommandSchemaID(name)

function _as_plant_command_schema_id(value)
    _command_schema_error(:invalid_id,
        "plant-command schema identity must be a Symbol or " *
        "PlantCommandSchemaID; got $(typeof(value))")
end

@inline _as_command_schema_version(version::PlantCommandSchemaVersion) = version
@inline _as_command_schema_version(version::Integer) =
    PlantCommandSchemaVersion(version)

function _as_command_schema_version(value)
    _command_schema_error(:invalid_version,
        "plant-command schema version must be an Integer or " *
        "PlantCommandSchemaVersion; got $(typeof(value))")
end

@inline _as_command_basis_revision(revision::CommandBasisRevision) = revision
@inline _as_command_basis_revision(revision::Integer) =
    CommandBasisRevision(revision)

function _as_command_basis_revision(value)
    _command_schema_error(:invalid_basis_revision,
        "command basis revision must be an Integer or CommandBasisRevision; " *
        "got $(typeof(value))")
end

@inline _as_command_unit(unit::CommandUnit) = unit
@inline _as_command_unit(name::Symbol) = CommandUnit(name)

function _as_command_unit(value)
    _command_schema_error(:invalid_unit,
        "command unit must be a Symbol or CommandUnit; got $(typeof(value))")
end

@inline _as_command_sign_convention(value::CommandSignConvention) = value
@inline _as_command_sign_convention(name::Symbol) = CommandSignConvention(name)

function _as_command_sign_convention(value)
    _command_schema_error(:invalid_sign_convention,
        "command sign convention must be a Symbol or " *
        "CommandSignConvention; got $(typeof(value))")
end

@inline _as_schema_command_endpoint_id(id::CommandEndpointID) = id

function _as_schema_command_endpoint_id(name::Symbol)
    _require_command_name(name, :invalid_endpoint,
        "command endpoint identity")
    return CommandEndpointID(name)
end

function _as_schema_command_endpoint_id(value)
    _command_schema_error(:invalid_endpoint,
        "command endpoint identity must be a Symbol or CommandEndpointID; " *
        "got $(typeof(value))")
end

@enum CommandValueSemantics::UInt8 begin
    AbsoluteCommand = 0x01
    IncrementalCommand = 0x02
end

@enum InvalidCommandAction::UInt8 begin
    RejectInvalidCommand = 0x01
    ClipInvalidCommand = 0x02
    FailOnInvalidCommand = 0x03
end

@enum CommandRangeStage::UInt8 begin
    ValidateOnPresentation = 0x01
    EnforceOnApplication = 0x02
end

@enum CommandSequenceAction::UInt8 begin
    AcceptSequence = 0x01
    RejectSequence = 0x02
    FailOnSequence = 0x03
end

@enum FutureCommandPolicy::UInt8 begin
    AllowFutureCommand = 0x01
    RejectFutureCommand = 0x02
end

@enum LateCommandPolicy::UInt8 begin
    RejectLateCommand = 0x01
    ApplyLateCommandNow = 0x02
    FailOnLateCommand = 0x03
end

@enum CommandSupersessionPolicy::UInt8 begin
    PreservePendingCommands = 0x01
    SupersedeOlderPendingCommands = 0x02
end

@enum CommandSilenceAction::UInt8 begin
    HoldLastCommand = 0x01
    ApplySafeCommand = 0x02
    FailOnCommandSilence = 0x03
end

@enum CommandAgeOrigin::UInt8 begin
    AgeFromAdmission = 0x01
    AgeFromApplication = 0x02
end

abstract type AbstractCommandBounds end

"""No finite lower or upper command-value bound."""
struct UnboundedCommandValues <: AbstractCommandBounds end

"""One inclusive lower/upper bound applied uniformly to all command values."""
struct UniformCommandBounds{T<:Real} <: AbstractCommandBounds
    lower::T
    upper::T

    function UniformCommandBounds(lower::T, upper::T) where {T<:Real}
        _require_command_numeric_type(T)
        isfinite(lower) && isfinite(upper) || _command_schema_error(
            :nonfinite_bounds, "uniform command bounds must be finite")
        lower <= upper || _command_schema_error(:invalid_bounds,
            "uniform command lower bound must not exceed its upper bound")
        return new{T}(lower, upper)
    end
end

function UniformCommandBounds(lower::Real, upper::Real)
    _command_schema_error(:bounds_type,
        "uniform command bounds must have one exact numeric type; got " *
        "$(typeof(lower)) and $(typeof(upper))")
end


function UniformCommandBounds(lower, upper)
    _command_schema_error(:bounds_type,
        "uniform command bounds must be real values; got $(typeof(lower)) " *
        "and $(typeof(upper))")
end

"""Finite-value and bounded-range handling for presented command payloads."""
struct CommandValuePolicy
    nonfinite::InvalidCommandAction
    out_of_range::InvalidCommandAction
    range_stage::CommandRangeStage

    function CommandValuePolicy(nonfinite::InvalidCommandAction,
        out_of_range::InvalidCommandAction,
        range_stage::CommandRangeStage)
        nonfinite == ClipInvalidCommand && _command_schema_error(
            :invalid_nonfinite_policy,
            "nonfinite command values cannot use a clipping policy")
        return new(nonfinite, out_of_range, range_stage)
    end
end

function CommandValuePolicy(nonfinite, out_of_range, range_stage)
    _command_schema_error(:invalid_value_policy,
        "command value policy requires InvalidCommandAction, " *
        "InvalidCommandAction, and CommandRangeStage values")
end

CommandValuePolicy(; nonfinite=RejectInvalidCommand,
    out_of_range=RejectInvalidCommand,
    range_stage=ValidateOnPresentation) =
    CommandValuePolicy(nonfinite, out_of_range, range_stage)

"""Duplicate, stale, reordered, and skipped sequence-number handling."""
struct CommandSequencePolicy
    duplicate::CommandSequenceAction
    stale::CommandSequenceAction
    reordered::CommandSequenceAction
    skipped::CommandSequenceAction

    function CommandSequencePolicy(duplicate::CommandSequenceAction,
        stale::CommandSequenceAction, reordered::CommandSequenceAction,
        skipped::CommandSequenceAction)
        duplicate == AcceptSequence && _command_schema_error(
            :invalid_duplicate_policy,
            "duplicate command sequences cannot be accepted as new commands")
        stale == AcceptSequence && _command_schema_error(
            :invalid_stale_policy,
            "stale command sequences cannot be accepted as new commands")
        return new(duplicate, stale, reordered, skipped)
    end
end

function CommandSequencePolicy(duplicate, stale, reordered, skipped)
    _command_schema_error(:invalid_sequence_policy,
        "command sequence policy fields must be CommandSequenceAction values")
end

CommandSequencePolicy(; duplicate=RejectSequence, stale=RejectSequence,
    reordered=RejectSequence, skipped=AcceptSequence) =
    CommandSequencePolicy(duplicate, stale, reordered, skipped)

"""Future, late, and pending-command supersession policy vocabulary."""
struct CommandEffectiveTimePolicy
    future::FutureCommandPolicy
    late::LateCommandPolicy
    supersession::CommandSupersessionPolicy

    function CommandEffectiveTimePolicy(future::FutureCommandPolicy,
        late::LateCommandPolicy,
        supersession::CommandSupersessionPolicy)
        return new(future, late, supersession)
    end
end

function CommandEffectiveTimePolicy(future, late, supersession)
    _command_schema_error(:invalid_effective_time_policy,
        "effective-time policy requires FutureCommandPolicy, " *
        "LateCommandPolicy, and CommandSupersessionPolicy values")
end

CommandEffectiveTimePolicy(; future=AllowFutureCommand,
    late=RejectLateCommand,
    supersession=PreservePendingCommands) =
    CommandEffectiveTimePolicy(future, late, supersession)

"""
    CommandSilencePolicy(action, age_origin; timeout=nothing)

Replayable plant-time behavior when a command endpoint remains silent. Holding
is indefinite and therefore has no timeout. Safe-command and failure actions
require a positive `PlantDuration`; the corresponding prepared safe value and
mutable age state belong to the later endpoint layer.
"""
struct CommandSilencePolicy{D<:Union{Nothing,PlantDuration}}
    action::CommandSilenceAction
    age_origin::CommandAgeOrigin
    timeout::D

    function CommandSilencePolicy(action::CommandSilenceAction,
        age_origin::CommandAgeOrigin, timeout::D) where {
        D<:Union{Nothing,PlantDuration}}
        _require_command_silence_timeout(action, timeout)
        return new{D}(action, age_origin, timeout)
    end
end

@inline function _require_command_silence_timeout(
    action::CommandSilenceAction, ::Nothing)
    action == HoldLastCommand || _command_schema_error(
        :missing_silence_timeout,
        "safe-command and failure silence policies require a timeout")
    return nothing
end


@inline function _require_command_silence_timeout(
    action::CommandSilenceAction, timeout::PlantDuration)
    action == HoldLastCommand && _command_schema_error(
        :invalid_silence_timeout,
        "indefinite hold must not declare a silence timeout")
    iszero(timeout) && _command_schema_error(:invalid_silence_timeout,
        "command-silence timeout must be positive")
    return nothing
end

function CommandSilencePolicy(action, age_origin, timeout)
    _command_schema_error(:invalid_silence_policy,
        "command silence policy requires CommandSilenceAction, " *
        "CommandAgeOrigin, and Nothing or PlantDuration")
end

CommandSilencePolicy(action::CommandSilenceAction,
    age_origin::CommandAgeOrigin; timeout=nothing) =
    CommandSilencePolicy(action, age_origin, timeout)

CommandSilencePolicy(; action=HoldLastCommand,
    age_origin=AgeFromApplication, timeout=nothing) =
    CommandSilencePolicy(action, age_origin, timeout)

@inline function _require_command_numeric_type(::Type{T}) where {T}
    T <: Real || _command_schema_error(:invalid_numeric_type,
        "command payload element type must be a subtype of Real; got $T")
    T === Bool && _command_schema_error(:invalid_numeric_type,
        "Bool is not a command payload numeric type")
    isconcretetype(T) || _command_schema_error(:invalid_numeric_type,
        "command payload element type must be concrete; got $T")
    isbitstype(T) || _command_schema_error(:invalid_numeric_type,
        "command payload element type must be isbits; got $T")
    return T
end

@inline _normalize_command_dimension(::Bool) = _command_schema_error(
    :invalid_dimensions, "command dimensions must not contain Bool")

function _normalize_command_dimension(dimension::Integer)
    try
        return Int(dimension)
    catch
        _command_schema_error(:invalid_dimensions,
            "command dimensions must fit Int")
    end
end

function _normalize_command_dimension(dimension)
    _command_schema_error(:invalid_dimensions,
        "command dimensions must be integer extents; got " *
        "$(typeof(dimension))")
end

function _normalize_command_dimensions(dimensions::Tuple)
    normalized = map(_normalize_command_dimension, dimensions)
    all(>(0), normalized) || _command_schema_error(:invalid_dimensions,
        "every command dimension must be positive")
    return normalized
end

@inline _require_command_basis(basis::CommandBasis) = basis

function _require_command_basis(basis)
    _command_schema_error(:invalid_basis,
        "command basis must be a CommandBasis; got $(typeof(basis))")
end

@inline _require_command_semantics(semantics::CommandValueSemantics) =
    semantics

function _require_command_semantics(semantics)
    _command_schema_error(:invalid_semantics,
        "command semantics must be a CommandValueSemantics value; got " *
        "$(typeof(semantics))")
end

@inline _require_command_value_policy(policy::CommandValuePolicy) = policy

function _require_command_value_policy(policy)
    _command_schema_error(:invalid_value_policy,
        "value_policy must be a CommandValuePolicy; got $(typeof(policy))")
end

@inline _require_command_sequence_policy(policy::CommandSequencePolicy) =
    policy

function _require_command_sequence_policy(policy)
    _command_schema_error(:invalid_sequence_policy,
        "sequence_policy must be a CommandSequencePolicy; got " *
        "$(typeof(policy))")
end

@inline _require_command_effective_time_policy(
    policy::CommandEffectiveTimePolicy) = policy

function _require_command_effective_time_policy(policy)
    _command_schema_error(:invalid_effective_time_policy,
        "effective_time_policy must be a CommandEffectiveTimePolicy; got " *
        "$(typeof(policy))")
end

@inline _require_command_silence_policy(policy::CommandSilencePolicy) = policy

function _require_command_silence_policy(policy)
    _command_schema_error(:invalid_silence_policy,
        "silence_policy must be a CommandSilencePolicy; got " *
        "$(typeof(policy))")
end

function _normalize_command_dimensions(dimensions)
    _command_schema_error(:invalid_dimensions,
        "command dimensions must be a Tuple; got $(typeof(dimensions))")
end

@inline _require_command_bounds(bounds::AbstractCommandBounds) = bounds

function _require_command_bounds(bounds)
    _command_schema_error(:invalid_bounds,
        "command bounds must be UnboundedCommandValues or " *
        "UniformCommandBounds; got $(typeof(bounds))")
end

@inline function _require_command_bounds_type(::Type{T},
    ::UnboundedCommandValues) where {T}
    return nothing
end

function _require_command_bounds_type(::Type{T},
    bounds::UniformCommandBounds{B}) where {T,B}
    T === B || _command_schema_error(:bounds_type,
        "command bounds use $B; expected exact payload element type $T")
    return nothing
end

@inline function _require_unbounded_value_policy(::UnboundedCommandValues,
    policy::CommandValuePolicy)
    policy.out_of_range == RejectInvalidCommand || _command_schema_error(
        :invalid_unbounded_policy,
        "unbounded commands must use the canonical reject range action")
    policy.range_stage == ValidateOnPresentation || _command_schema_error(
        :invalid_unbounded_policy,
        "unbounded commands must use the canonical presentation range stage")
    return nothing
end

@inline _require_unbounded_value_policy(::UniformCommandBounds,
    ::CommandValuePolicy) = nothing

struct _PlantCommandSchemaToken end
const _PLANT_COMMAND_SCHEMA_TOKEN = _PlantCommandSchemaToken()

"""
    PlantCommandSchema(T, dimensions; ...)

Immutable semantic contract for payloads presented to one core command
endpoint. `T` is the exact scalar/array element type. An empty `dimensions`
tuple denotes a scalar payload; otherwise payloads are backend-neutral
`AbstractArray{T,N}` values with that exact shape.

The schema contains no sequence state, future calendar, applied value, session,
source-clock mapping, payload lease, queue, transport, or terminal outcome.
"""
struct PlantCommandSchema{T,N,B<:AbstractCommandBounds,
    S<:CommandSilencePolicy}
    id::PlantCommandSchemaID
    version::PlantCommandSchemaVersion
    endpoint::CommandEndpointID
    dimensions::NTuple{N,Int}
    units::CommandUnit
    sign_convention::CommandSignConvention
    basis::CommandBasis
    basis_revision::CommandBasisRevision
    semantics::CommandValueSemantics
    bounds::B
    value_policy::CommandValuePolicy
    sequence_policy::CommandSequencePolicy
    effective_time_policy::CommandEffectiveTimePolicy
    silence_policy::S

    function PlantCommandSchema(::_PlantCommandSchemaToken, ::Type{T},
        dimensions::NTuple{N,Int}, id::PlantCommandSchemaID,
        version::PlantCommandSchemaVersion, endpoint::CommandEndpointID,
        units::CommandUnit, sign_convention::CommandSignConvention,
        basis::CommandBasis, basis_revision::CommandBasisRevision,
        semantics::CommandValueSemantics, bounds::B,
        value_policy::CommandValuePolicy,
        sequence_policy::CommandSequencePolicy,
        effective_time_policy::CommandEffectiveTimePolicy,
        silence_policy::S) where {T,N,B<:AbstractCommandBounds,
        S<:CommandSilencePolicy}
        _require_command_numeric_type(T)
        _require_command_bounds_type(T, bounds)
        _require_unbounded_value_policy(bounds, value_policy)
        return new{T,N,B,S}(id, version, endpoint, dimensions, units,
            sign_convention, basis, basis_revision, semantics, bounds,
            value_policy, sequence_policy, effective_time_policy,
            silence_policy)
    end
end

function PlantCommandSchema(::Type{T}, dimensions; id, version, endpoint,
    units, sign_convention, basis, basis_revision, semantics, bounds,
    value_policy, sequence_policy, effective_time_policy, silence_policy) where {T}
    normalized_dimensions = _normalize_command_dimensions(dimensions)
    resolved_bounds = _require_command_bounds(bounds)
    resolved_basis = _require_command_basis(basis)
    resolved_semantics = _require_command_semantics(semantics)
    resolved_value_policy = _require_command_value_policy(value_policy)
    resolved_sequence_policy = _require_command_sequence_policy(
        sequence_policy)
    resolved_effective_time_policy = _require_command_effective_time_policy(
        effective_time_policy)
    resolved_silence_policy = _require_command_silence_policy(silence_policy)
    return PlantCommandSchema(
        _PLANT_COMMAND_SCHEMA_TOKEN,
        T,
        normalized_dimensions,
        _as_plant_command_schema_id(id),
        _as_command_schema_version(version),
        _as_schema_command_endpoint_id(endpoint),
        _as_command_unit(units),
        _as_command_sign_convention(sign_convention),
        resolved_basis,
        _as_command_basis_revision(basis_revision),
        resolved_semantics,
        resolved_bounds,
        resolved_value_policy,
        resolved_sequence_policy,
        resolved_effective_time_policy,
        resolved_silence_policy,
    )
end

@inline command_schema_id(schema::PlantCommandSchema) = schema.id
@inline command_schema_version(schema::PlantCommandSchema) = schema.version
@inline command_endpoint_id(schema::PlantCommandSchema) = schema.endpoint
@inline command_numeric_type(::PlantCommandSchema{T}) where {T} = T
@inline command_dimensions(schema::PlantCommandSchema) = schema.dimensions
@inline command_units(schema::PlantCommandSchema) = schema.units
@inline command_sign_convention(schema::PlantCommandSchema) =
    schema.sign_convention
@inline command_basis(schema::PlantCommandSchema) = schema.basis
@inline command_basis_revision(schema::PlantCommandSchema) =
    schema.basis_revision
@inline command_semantics(schema::PlantCommandSchema) = schema.semantics
@inline command_bounds(schema::PlantCommandSchema) = schema.bounds
@inline command_value_policy(schema::PlantCommandSchema) = schema.value_policy
@inline command_sequence_policy(schema::PlantCommandSchema) =
    schema.sequence_policy
@inline command_effective_time_policy(schema::PlantCommandSchema) =
    schema.effective_time_policy
@inline command_silence_policy(schema::PlantCommandSchema) =
    schema.silence_policy

@inline function _payload_violation_reason(kind::Symbol,
    action::InvalidCommandAction)
    suffix = action == FailOnInvalidCommand ? "failure" : "rejected"
    return Symbol(kind, "_", suffix)
end

@noinline function _throw_command_payload_violation(
    schema::PlantCommandSchema, kind::Symbol,
    action::InvalidCommandAction, message::AbstractString)
    throw(PlantCommandError(:payload,
        _payload_violation_reason(kind, action),
        "command endpoint $(schema.endpoint): $message"))
end

@inline _all_command_values_finite(value::Real) = isfinite(value)
@inline _all_command_values_finite(values::AbstractArray) =
    all(isfinite, values)

@inline _all_command_values_in_bounds(value::Real,
    bounds::UniformCommandBounds) = bounds.lower <= value <= bounds.upper

@inline function _all_command_values_in_bounds(values::AbstractArray,
    bounds::UniformCommandBounds)
    lower = bounds.lower
    upper = bounds.upper
    return all(value -> lower <= value <= upper, values)
end

@inline function _validate_command_payload_finite(schema::PlantCommandSchema,
    payload)
    _all_command_values_finite(payload) && return payload
    policy = schema.value_policy
    return _throw_command_payload_violation(schema, :nonfinite,
        policy.nonfinite, "payload contains a nonfinite value")
end

@inline _validate_command_payload_range(schema::PlantCommandSchema,
    payload, ::UnboundedCommandValues) = payload

@inline function _validate_command_payload_range(schema::PlantCommandSchema,
    payload, bounds::UniformCommandBounds)
    policy = schema.value_policy
    policy.range_stage == EnforceOnApplication && return payload
    _all_command_values_in_bounds(payload, bounds) && return payload
    policy.out_of_range == ClipInvalidCommand && return payload
    return _throw_command_payload_violation(schema, :out_of_range,
        policy.out_of_range,
        "payload contains a value outside the inclusive command bounds")
end

@inline function _validate_command_payload_values(
    schema::PlantCommandSchema, payload)
    _validate_command_payload_finite(schema, payload)
    return _validate_command_payload_range(schema, payload, schema.bounds)
end

"""
    validate_plant_command_payload(schema, payload)

Validate one caller payload against the schema's exact scalar/array type,
shape, finite-value policy, and any presentation-stage range policy. The
function returns the original payload and never clips or mutates it. A clipping
policy means later admission/application may clip; application-stage range
checks are likewise deferred. This function does not admit, sequence, schedule,
or apply a command.
"""
@inline function validate_plant_command_payload(
    schema::PlantCommandSchema{T,0}, payload::Real) where {T}
    typeof(payload) === T || throw(PlantCommandError(:payload, :element_type,
        "command endpoint $(schema.endpoint): scalar type " *
        "$(typeof(payload)) does not match $T"))
    _validate_command_payload_values(schema, payload)
    return payload
end

@inline function validate_plant_command_payload(
    schema::PlantCommandSchema{T,N},
    payload::AbstractArray) where {T,N}
    eltype(payload) === T || throw(PlantCommandError(
        :payload, :element_type,
        "command endpoint $(schema.endpoint): payload element type " *
        "$(eltype(payload)) does not match $T"))
    size(payload) == schema.dimensions || throw(PlantCommandError(
        :payload, :shape,
        "command endpoint $(schema.endpoint): payload shape $(size(payload)) " *
        "does not match $(schema.dimensions)"))
    _validate_command_payload_values(schema, payload)
    return payload
end

function validate_plant_command_payload(schema::PlantCommandSchema,
    ::Real)
    throw(PlantCommandError(:payload, :shape,
        "command endpoint $(schema.endpoint): scalar payload does not match " *
        "$(schema.dimensions)"))
end

function validate_plant_command_payload(schema::PlantCommandSchema, payload)
    throw(PlantCommandError(:payload, :payload_type,
        "command endpoint $(schema.endpoint): payload must be an exact " *
        "$(command_numeric_type(schema)) " *
        "scalar or AbstractArray; got $(typeof(payload))"))
end

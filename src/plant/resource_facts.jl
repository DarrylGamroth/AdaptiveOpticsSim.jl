#
# Versioned structural resource facts
#
# These values describe explicitly owned prepared storage. Reporters enumerate
# their allocations by dispatch; aggregation never walks object fields or uses
# allocator/object-size reflection.
#

@inline function _require_resource_symbol(value::Symbol, label::AbstractString)
    isempty(String(value)) && throw(StructuralResourceError(
        :identity, :empty_identity, "$label must not be empty"))
    return value
end

"""Stable logical identity of one prepared structural-storage owner."""
struct StructuralResourceOwnerID
    category::Symbol
    component::Symbol

    function StructuralResourceOwnerID(category::Symbol, component::Symbol)
        return new(
            _require_resource_symbol(category, "resource-owner category"),
            _require_resource_symbol(component, "resource-owner component"),
        )
    end
end

Base.:(==)(left::StructuralResourceOwnerID,
    right::StructuralResourceOwnerID) =
    left.category == right.category && left.component == right.component
Base.isequal(left::StructuralResourceOwnerID,
    right::StructuralResourceOwnerID) =
    isequal(left.category, right.category) &&
    isequal(left.component, right.component)
Base.hash(owner::StructuralResourceOwnerID, seed::UInt) =
    hash(owner.component,
        hash(owner.category, hash(StructuralResourceOwnerID, seed)))
Base.isless(left::StructuralResourceOwnerID,
    right::StructuralResourceOwnerID) =
    left.category == right.category ?
    isless(left.component, right.component) :
    isless(left.category, right.category)

function Base.show(io::IO, owner::StructuralResourceOwnerID)
    print(io, "StructuralResourceOwnerID(", repr(owner.category), ", ",
        repr(owner.component), ")")
end

"""Deterministic identity and version of one resource-byte estimate method."""
struct ResourceEstimateMethod
    identity::Symbol
    version::UInt32

    function ResourceEstimateMethod(identity::Symbol,
        version::UInt32)
        iszero(version) && throw(StructuralResourceError(
            :estimate_method, :invalid_version,
            "resource-estimate method version must be positive"))
        return new(_require_resource_symbol(identity,
            "resource-estimate method identity"), version)
    end
end

function ResourceEstimateMethod(identity::Symbol, version::Integer)
    version > 0 || throw(StructuralResourceError(
        :estimate_method, :invalid_version,
        "resource-estimate method version must be positive"))
    version <= typemax(UInt32) || throw(StructuralResourceError(
        :estimate_method, :invalid_version,
        "resource-estimate method version exceeds UInt32"))
    return ResourceEstimateMethod(identity, UInt32(version))
end

function ResourceEstimateMethod(::Symbol, ::Bool)
    throw(StructuralResourceError(
        :estimate_method, :invalid_version,
        "resource-estimate method version must be an integer, not Bool"))
end

Base.:(==)(left::ResourceEstimateMethod,
    right::ResourceEstimateMethod) =
    left.identity == right.identity && left.version == right.version
Base.isequal(left::ResourceEstimateMethod,
    right::ResourceEstimateMethod) =
    isequal(left.identity, right.identity) &&
    isequal(left.version, right.version)
Base.hash(method::ResourceEstimateMethod, seed::UInt) =
    hash(method.version,
        hash(method.identity,
            hash(ResourceEstimateMethod, seed)))

function Base.show(io::IO, method::ResourceEstimateMethod)
    print(io, "ResourceEstimateMethod(", repr(method.identity),
        ", ", method.version, ")")
end

const _STRUCTURAL_RESOURCE_ESTIMATE_METHOD =
    ResourceEstimateMethod(:explicit_owned_array_storage, 1)

abstract type AbstractStructuralResourceFact end

@inline function _checked_resource_byte_count(value::Integer,
    component::Symbol)
    value >= 0 || throw(StructuralResourceError(component,
        :invalid_byte_count, "resource byte count must not be negative"))
    value <= typemax(UInt64) || throw(StructuralResourceError(component,
        :invalid_byte_count, "resource byte count exceeds UInt64"))
    return UInt64(value)
end


function _checked_resource_byte_count(::Bool, component::Symbol)
    throw(StructuralResourceError(component, :invalid_byte_count,
        "resource byte count must be an integer, not Bool"))
end

function _checked_resource_byte_count(value, component::Symbol)
    throw(StructuralResourceError(component, :invalid_byte_count,
        "resource byte count must be an Integer; got $(typeof(value))"))
end

"""Exact structural byte counts for one logical owner on one exact device."""
struct KnownStructuralResourceFact{D<:AbstractComputeDevice} <:
    AbstractStructuralResourceFact
    owner::StructuralResourceOwnerID
    device::D
    resident_bytes::UInt64
    workspace_bytes::UInt64
    method::ResourceEstimateMethod
end

function KnownStructuralResourceFact(owner::StructuralResourceOwnerID,
    device::D, resident_bytes, workspace_bytes,
    method::ResourceEstimateMethod=
        _STRUCTURAL_RESOURCE_ESTIMATE_METHOD) where {D<:AbstractComputeDevice}
    return KnownStructuralResourceFact{D}(
        owner,
        device,
        _checked_resource_byte_count(resident_bytes, :resident_bytes),
        _checked_resource_byte_count(workspace_bytes, :workspace_bytes),
        method,
    )
end

"""Explicitly unknown structural byte counts for one logical owner."""
struct UnknownStructuralResourceFact{D<:AbstractComputeDevice} <:
    AbstractStructuralResourceFact
    owner::StructuralResourceOwnerID
    device::D
    reason::Symbol
    method::ResourceEstimateMethod

    function UnknownStructuralResourceFact(owner::StructuralResourceOwnerID,
        device::D, reason::Symbol,
        method::ResourceEstimateMethod=
            _STRUCTURAL_RESOURCE_ESTIMATE_METHOD) where {
        D<:AbstractComputeDevice,
    }
        return new{D}(owner, device,
            _require_resource_symbol(reason,
                "unknown-resource reason"), method)
    end
end

"""Caller- or provider-supplied reserve for storage core cannot inspect."""
struct OpaqueResourceReserve{D<:AbstractComputeDevice}
    owner::StructuralResourceOwnerID
    device::D
    provider::Symbol
    bytes::UInt64
    method::ResourceEstimateMethod

    function OpaqueResourceReserve(owner::StructuralResourceOwnerID,
        device::D, provider::Symbol, bytes,
        method::ResourceEstimateMethod) where {
        D<:AbstractComputeDevice,
    }
        return new{D}(
            owner,
            device,
            _require_resource_symbol(provider,
                "opaque-reserve provider"),
            _checked_resource_byte_count(bytes, :opaque_reserve),
            method,
        )
    end
end

struct _StructuralResourceReportToken end
const _STRUCTURAL_RESOURCE_REPORT_TOKEN = _StructuralResourceReportToken()

struct _FixedStructuralResourceRegistry{T,N} <: AbstractVector{T}
    _storage::NTuple{N,T}
end

function _fixed_structural_resource_registry(::Type{T},
    values::AbstractVector{T}) where {T}
    count = length(values)
    return _FixedStructuralResourceRegistry{T,count}(Tuple(values))
end

Base.size(registry::_FixedStructuralResourceRegistry) =
    (length(getfield(registry, :_storage)),)
Base.axes(registry::_FixedStructuralResourceRegistry) =
    axes(getfield(registry, :_storage))
Base.length(registry::_FixedStructuralResourceRegistry) =
    length(getfield(registry, :_storage))
Base.getindex(registry::_FixedStructuralResourceRegistry, index::Int) =
    getfield(registry, :_storage)[index]
Base.IndexStyle(::Type{<:_FixedStructuralResourceRegistry}) = IndexLinear()
Base.iterate(registry::_FixedStructuralResourceRegistry, state...) =
    iterate(getfield(registry, :_storage), state...)
Base.copy(registry::_FixedStructuralResourceRegistry) =
    collect(getfield(registry, :_storage))

"""Validated exact-device resource facts and their overflow-safe totals."""
struct StructuralResourceReport{D<:AbstractComputeDevice,NF,NR}
    device::D
    facts::_FixedStructuralResourceRegistry{
        KnownStructuralResourceFact{D},NF}
    reserves::_FixedStructuralResourceRegistry{OpaqueResourceReserve{D},NR}
    resident_bytes::UInt64
    workspace_bytes::UInt64
    opaque_reserve_bytes::UInt64
    method::ResourceEstimateMethod

    function StructuralResourceReport(
        ::_StructuralResourceReportToken,
        device::D,
        facts::_FixedStructuralResourceRegistry{
            KnownStructuralResourceFact{D},NF},
        reserves::_FixedStructuralResourceRegistry{
            OpaqueResourceReserve{D},NR},
        resident_bytes::UInt64,
        workspace_bytes::UInt64,
        opaque_reserve_bytes::UInt64,
        method::ResourceEstimateMethod,
    ) where {D<:AbstractComputeDevice,NF,NR}
        return new{D,NF,NR}(device, facts, reserves, resident_bytes,
            workspace_bytes, opaque_reserve_bytes, method)
    end
end

@inline structural_resource_owner_id(fact::AbstractStructuralResourceFact) =
    fact.owner
@inline structural_resource_owner_id(reserve::OpaqueResourceReserve) =
    reserve.owner
@inline compute_device(fact::AbstractStructuralResourceFact) = fact.device
@inline compute_device(reserve::OpaqueResourceReserve) = reserve.device
@inline compute_device(report::StructuralResourceReport) = report.device
@inline structural_resource_known(::KnownStructuralResourceFact) = true
@inline structural_resource_known(::UnknownStructuralResourceFact) = false
@inline structural_resource_unknown_reason(
    fact::UnknownStructuralResourceFact) = fact.reason
@inline structural_resident_bytes(fact::KnownStructuralResourceFact) =
    fact.resident_bytes
@inline structural_resident_bytes(report::StructuralResourceReport) =
    report.resident_bytes
@inline structural_workspace_bytes(fact::KnownStructuralResourceFact) =
    fact.workspace_bytes
@inline structural_workspace_bytes(report::StructuralResourceReport) =
    report.workspace_bytes
@inline opaque_resource_reserve_bytes(reserve::OpaqueResourceReserve) =
    reserve.bytes
@inline opaque_resource_reserve_bytes(report::StructuralResourceReport) =
    report.opaque_reserve_bytes
@inline resource_estimate_method(fact::AbstractStructuralResourceFact) =
    fact.method
@inline resource_estimate_method(reserve::OpaqueResourceReserve) =
    reserve.method
@inline resource_estimate_method(report::StructuralResourceReport) =
    report.method
@inline structural_resource_facts(report::StructuralResourceReport) =
    report.facts
@inline opaque_resource_reserves(report::StructuralResourceReport) =
    report.reserves

@noinline function _structural_resource_error(component::Symbol,
    reason::Symbol, message::AbstractString)
    throw(StructuralResourceError(component, reason, String(message)))
end


@noinline function _handle_resource_add_error(::OverflowError,
    component::Symbol)
    _structural_resource_error(component, :arithmetic_overflow,
        "$component byte total exceeds UInt64")
end

@noinline _handle_resource_add_error(error, ::Symbol) = throw(error)

@inline function _checked_resource_add(left::UInt64, right::UInt64,
    component::Symbol)
    try
        return Base.checked_add(left, right)
    catch error
        _handle_resource_add_error(error, component)
    end
end

@inline function _require_fact_target(
    fact::AbstractStructuralResourceFact,
    target::AbstractComputeDevice,
)
    compute_device(fact) == target || _structural_resource_error(
        :fact, :wrong_device,
        "resource fact for $(structural_resource_owner_id(fact)) occupies " *
        "$(compute_device(fact)); expected $target")
    return fact
end

function _require_fact_target(fact, ::AbstractComputeDevice)
    _structural_resource_error(:fact, :invalid_fact,
        "structural facts must contain AbstractStructuralResourceFact " *
        "values; got $(typeof(fact))")
end

@inline function _require_known_fact(
    fact::KnownStructuralResourceFact,
    target::AbstractComputeDevice,
)
    return _require_fact_target(fact, target)
end

function _require_known_fact(
    fact::UnknownStructuralResourceFact,
    target::AbstractComputeDevice,
)
    _require_fact_target(fact, target)
    _structural_resource_error(:fact, :unknown_resource,
        "resource fact for $(structural_resource_owner_id(fact)) is " *
        "unknown: $(structural_resource_unknown_reason(fact))")
end

function _require_known_fact(fact, ::AbstractComputeDevice)
    _structural_resource_error(:fact, :invalid_fact,
        "structural facts must contain AbstractStructuralResourceFact " *
        "values; got " *
        "$(typeof(fact))")
end

@inline function _require_resource_method(
    fact::KnownStructuralResourceFact,
    method::ResourceEstimateMethod,
)
    resource_estimate_method(fact) == method || _structural_resource_error(
        :fact, :inconsistent_estimate_method,
        "resource fact for $(structural_resource_owner_id(fact)) uses " *
        "$(resource_estimate_method(fact)); expected $method")
    return fact
end

function _require_reserve_target(reserve::OpaqueResourceReserve,
    target::AbstractComputeDevice)
    compute_device(reserve) == target || _structural_resource_error(
        :opaque_reserve, :wrong_device,
        "opaque reserve for $(structural_resource_owner_id(reserve)) " *
        "occupies $(compute_device(reserve)); expected $target")
    return reserve
end

function _require_reserve_target(reserve, ::AbstractComputeDevice)
    _structural_resource_error(:opaque_reserve, :invalid_reserve,
        "opaque reserves must contain OpaqueResourceReserve values; got " *
        "$(typeof(reserve))")
end

"""
    aggregate_structural_resource_facts(facts, target; opaque_reserves=())

Validate and aggregate one complete exact-device graph. Unknown facts, duplicate
logical owners, inconsistent structural-estimate methods, wrong-device values,
and arithmetic overflow fail structurally. Opaque reserves remain a separate
total and may use provider-specific method identities.
"""
function aggregate_structural_resource_facts(
    facts,
    target::D;
    opaque_reserves=(),
) where {D<:AbstractComputeDevice}
    frozen_facts = Tuple(facts)
    isempty(frozen_facts) && _structural_resource_error(:fact,
        :empty_facts, "at least one structural resource fact is required")
    first_fact = _require_known_fact(first(frozen_facts), target)
    method = resource_estimate_method(first_fact)
    owners = Set{StructuralResourceOwnerID}()
    known_fact_values = KnownStructuralResourceFact{D}[]
    resident_bytes = UInt64(0)
    workspace_bytes = UInt64(0)
    for candidate in frozen_facts
        fact = _require_resource_method(
            _require_known_fact(candidate, target), method)
        owner = structural_resource_owner_id(fact)
        owner in owners && _structural_resource_error(:fact,
            :duplicate_owner,
            "exact device $target has duplicate resource owner $owner")
        push!(owners, owner)
        push!(known_fact_values, fact)
        resident_bytes = _checked_resource_add(resident_bytes,
            structural_resident_bytes(fact), :resident_bytes)
        workspace_bytes = _checked_resource_add(workspace_bytes,
            structural_workspace_bytes(fact), :workspace_bytes)
    end

    frozen_reserves = Tuple(opaque_reserves)
    reserve_keys = Set{Tuple{StructuralResourceOwnerID,Symbol}}()
    known_reserve_values = OpaqueResourceReserve{D}[]
    opaque_bytes = UInt64(0)
    for candidate in frozen_reserves
        reserve = _require_reserve_target(candidate, target)
        owner = structural_resource_owner_id(reserve)
        owner in owners || _structural_resource_error(:opaque_reserve,
            :unknown_owner,
            "opaque reserve owner $owner has no structural resource fact " *
            "on exact device $target")
        key = (owner, reserve.provider)
        key in reserve_keys && _structural_resource_error(:opaque_reserve,
            :duplicate_reserve,
            "exact device $target has duplicate opaque reserve $key")
        push!(reserve_keys, key)
        push!(known_reserve_values, reserve)
        opaque_bytes = _checked_resource_add(opaque_bytes,
            opaque_resource_reserve_bytes(reserve), :opaque_reserve)
    end

    sort!(known_fact_values; by=structural_resource_owner_id)
    known_facts = _fixed_structural_resource_registry(
        KnownStructuralResourceFact{D}, known_fact_values)
    sort!(known_reserve_values; by=reserve ->
        (structural_resource_owner_id(reserve), reserve.provider))
    known_reserves = _fixed_structural_resource_registry(
        OpaqueResourceReserve{D}, known_reserve_values)
    return StructuralResourceReport(
        _STRUCTURAL_RESOURCE_REPORT_TOKEN,
        target,
        known_facts,
        known_reserves,
        resident_bytes,
        workspace_bytes,
        opaque_bytes,
        method,
    )
end

"""Fail-closed extension seam for one prepared logical resource owner."""
function structural_resource_fact(owner,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice)
    return UnknownStructuralResourceFact(id, target, :unsupported_model)
end

@inline function _require_resource_owner(
    value::AbstractStructuralResourceFact,
    expected::StructuralResourceOwnerID,
)
    actual = structural_resource_owner_id(value)
    actual == expected || _structural_resource_error(:fact, :wrong_owner,
        "resource reporter returned owner $actual; expected $expected")
    return value
end

@inline function _require_resource_owner(
    value::OpaqueResourceReserve,
    expected::StructuralResourceOwnerID,
)
    actual = structural_resource_owner_id(value)
    actual == expected || _structural_resource_error(
        :opaque_reserve, :wrong_owner,
        "opaque reserve belongs to owner $actual; expected $expected")
    return value
end

"""
    require_exact_structural_resource_facts(owner, id, target;
        opaque_reserves=())

Collect and require exact facts for one prepared logical owner. Whole-graph
overloads add their facts explicitly without reflective field traversal.
"""
function require_exact_structural_resource_facts(owner,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice;
    opaque_reserves=(),
)
    fact = _require_resource_owner(
        _require_fact_target(structural_resource_fact(owner, id, target),
            target),
        id,
    )
    frozen_reserves = Tuple(opaque_reserves)
    validated_reserves = OpaqueResourceReserve{typeof(target)}[]
    for candidate in frozen_reserves
        reserve = _require_resource_owner(
            _require_reserve_target(candidate, target), id)
        push!(validated_reserves, reserve)
    end
    return aggregate_structural_resource_facts((fact,), target;
        opaque_reserves=validated_reserves)
end

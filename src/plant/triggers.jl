#
# Deterministic modeled trigger sources and finite distribution topologies
#
# This serial core is deliberately independent of wall clocks, tasks, rings,
# transports, and detector transition semantics. Exact finite fault traces are
# resolved during preparation; execution retains only per-node phase state and
# a fixed-capacity calendar of delivered edges.
#

const _TRIGGER_TOPOLOGY_COMPONENT = :trigger_topology

@noinline function _trigger_topology_error(reason::Symbol,
    message::AbstractString)
    throw(PlantScheduleError(_TRIGGER_TOPOLOGY_COMPONENT, reason,
        String(message)))
end

"""Stable declared identity of one modeled trigger source."""
struct TriggerSourceID
    name::Symbol

    TriggerSourceID(name::Symbol) = new(
        _require_component_name(name, :trigger_source))
end

"""Stable declared identity of one trigger-distribution link."""
struct TriggerLinkID
    name::Symbol

    TriggerLinkID(name::Symbol) = new(
        _require_component_name(name, :trigger_link))
end

"""Stable declared identity of one trigger consumer."""
struct TriggerConsumerID
    name::Symbol

    TriggerConsumerID(name::Symbol) = new(
        _require_component_name(name, :trigger_consumer))
end

"""Stable declared identity of one configured trigger fault."""
struct TriggerFaultID
    name::Symbol

    TriggerFaultID(name::Symbol) = new(
        _require_component_name(name, :trigger_fault))
end

const _TriggerIdentity = Union{
    TriggerSourceID,TriggerLinkID,TriggerConsumerID,TriggerFaultID}

Base.:(==)(left::T, right::T) where {T<:_TriggerIdentity} =
    left.name == right.name
Base.isequal(left::T, right::T) where {T<:_TriggerIdentity} =
    isequal(left.name, right.name)
Base.hash(id::T, seed::UInt) where {T<:_TriggerIdentity} =
    hash(id.name, hash(T, seed))

function Base.show(io::IO, id::_TriggerIdentity)
    print(io, nameof(typeof(id)), "(", repr(id.name), ")")
end

@inline _as_trigger_source_id(id::TriggerSourceID) = id
@inline _as_trigger_source_id(name::Symbol) = TriggerSourceID(name)
@inline _as_trigger_link_id(id::TriggerLinkID) = id
@inline _as_trigger_link_id(name::Symbol) = TriggerLinkID(name)
@inline _as_trigger_consumer_id(id::TriggerConsumerID) = id
@inline _as_trigger_consumer_id(name::Symbol) = TriggerConsumerID(name)
@inline _as_trigger_fault_id(id::TriggerFaultID) = id
@inline _as_trigger_fault_id(name::Symbol) = TriggerFaultID(name)

function _invalid_trigger_id(kind::Symbol, value)
    _trigger_topology_error(:invalid_id,
        "$kind identity must be a Symbol or its typed trigger identity; got " *
        "$(typeof(value))")
end

_as_trigger_source_id(value) = _invalid_trigger_id(:trigger_source, value)
_as_trigger_link_id(value) = _invalid_trigger_id(:trigger_link, value)
_as_trigger_consumer_id(value) = _invalid_trigger_id(:trigger_consumer, value)
_as_trigger_fault_id(value) = _invalid_trigger_id(:trigger_fault, value)

@enum TriggerEdgeAction::UInt8 begin
    PassTriggerEdge = 0x00
    DropTriggerEdge = 0x01
    DuplicateTriggerEdge = 0x02
end

@inline function _checked_trigger_sequence(value::Integer)
    1 <= value <= typemax(UInt64) || _trigger_topology_error(
        :invalid_sequence,
        "trigger source sequence must be in 1:typemax(UInt64)")
    return UInt64(value)
end

@inline _checked_trigger_sequence(::Bool) = _trigger_topology_error(
    :invalid_sequence,
    "trigger source sequence must be an integer count, not Bool")

@inline _checked_trigger_sequence(value) = _trigger_topology_error(
    :invalid_sequence,
    "trigger source sequence must be an integer count; got $(typeof(value))")

@inline function _next_trigger_sequence(sequence::UInt64)
    sequence != typemax(UInt64) || _trigger_topology_error(
        :sequence_overflow, "trigger source sequence exceeds UInt64 range")
    return sequence + UInt64(1)
end

"""
One exact finite-trace update, keyed by the one-based sequence of its root
trigger source. `phase_step` persists for later sequences; `jitter` and
`timestamp_label_offset` affect only this sequence. A duplicate is an explicit
second occurrence after `duplicate_delay`.
"""
struct TriggerFaultTraceEntry
    sequence::UInt64
    jitter::PlantTimeOffset
    phase_step::PlantTimeOffset
    timestamp_label_offset::PlantTimeOffset
    action::TriggerEdgeAction
    duplicate_delay::PlantDuration
    fault_id::TriggerFaultID
end

function TriggerFaultTraceEntry(sequence, fault_id;
    jitter::PlantTimeOffset=zero(PlantTimeOffset),
    phase_step::PlantTimeOffset=zero(PlantTimeOffset),
    timestamp_label_offset::PlantTimeOffset=zero(PlantTimeOffset),
    action::TriggerEdgeAction=PassTriggerEdge,
    duplicate_delay::PlantDuration=zero(PlantDuration))
    checked_sequence = _checked_trigger_sequence(sequence)
    if action != DuplicateTriggerEdge && !iszero(duplicate_delay)
        _trigger_topology_error(:invalid_fault_trace,
            "duplicate_delay must be zero unless action is " *
            "DuplicateTriggerEdge")
    end
    return TriggerFaultTraceEntry(checked_sequence, jitter, phase_step,
        timestamp_label_offset, action, duplicate_delay,
        _as_trigger_fault_id(fault_id))
end

@inline _require_trigger_fault_trace_entry(entry::TriggerFaultTraceEntry) =
    entry

function _require_trigger_fault_trace_entry(value)
    _trigger_topology_error(:invalid_fault_trace,
        "trigger fault traces must contain TriggerFaultTraceEntry values; " *
        "got $(typeof(value))")
end

"""Immutable finite exact trace used by one trigger source or link."""
struct TriggerFaultTrace
    entries::Tuple

    function TriggerFaultTrace(entries::Tuple)
        Base.@nospecialize entries
        foreach(_require_trigger_fault_trace_entry, entries)
        return new(entries)
    end
end

TriggerFaultTrace() = TriggerFaultTrace(())
TriggerFaultTrace(entries::TriggerFaultTraceEntry...) =
    TriggerFaultTrace(entries)
TriggerFaultTrace(entries::AbstractVector) = TriggerFaultTrace(Tuple(entries))

const _TriggerParentID = Union{TriggerSourceID,TriggerLinkID}

@inline _require_trigger_parent(parent::_TriggerParentID) = parent

function _require_trigger_parent(value)
    _trigger_topology_error(:invalid_parent,
        "trigger parent must be TriggerSourceID or TriggerLinkID; got " *
        "$(typeof(value))")
end

"""Cold declaration of one periodic modeled trigger source."""
struct TriggerSourceDefinition
    id::TriggerSourceID
    schedule::PeriodicSchedule
    origin::PlantTimestamp
    timing_offset::PlantTimeOffset
    timestamp_label_offset::PlantTimeOffset
    faults::TriggerFaultTrace
end

function TriggerSourceDefinition(id, schedule::PeriodicSchedule;
    origin::PlantTimestamp=zero(PlantTimestamp),
    timing_offset::PlantTimeOffset=zero(PlantTimeOffset),
    timestamp_label_offset::PlantTimeOffset=zero(PlantTimeOffset),
    faults::TriggerFaultTrace=TriggerFaultTrace())
    return TriggerSourceDefinition(_as_trigger_source_id(id), schedule,
        origin, timing_offset, timestamp_label_offset, faults)
end

"""Cold declaration of one edge in a finite trigger fan-out."""
struct TriggerLinkDefinition
    id::TriggerLinkID
    parent::_TriggerParentID
    propagation_delay::PlantDuration
    timing_skew::PlantTimeOffset
    timestamp_label_offset::PlantTimeOffset
    faults::TriggerFaultTrace

    function TriggerLinkDefinition(id::TriggerLinkID,
        parent::_TriggerParentID, propagation_delay::PlantDuration,
        timing_skew::PlantTimeOffset,
        timestamp_label_offset::PlantTimeOffset,
        faults::TriggerFaultTrace)
        return new(id, parent, propagation_delay, timing_skew,
            timestamp_label_offset, faults)
    end
end

function TriggerLinkDefinition(id, parent;
    propagation_delay::PlantDuration=zero(PlantDuration),
    timing_skew::PlantTimeOffset=zero(PlantTimeOffset),
    timestamp_label_offset::PlantTimeOffset=zero(PlantTimeOffset),
    faults::TriggerFaultTrace=TriggerFaultTrace())
    return TriggerLinkDefinition(_as_trigger_link_id(id),
        _require_trigger_parent(parent), propagation_delay, timing_skew,
        timestamp_label_offset, faults)
end

"""Cold binding from one source or distribution link to one consumer."""
struct TriggerConsumerDefinition
    id::TriggerConsumerID
    parent::_TriggerParentID

    function TriggerConsumerDefinition(id::TriggerConsumerID,
        parent::_TriggerParentID)
        return new(id, parent)
    end
end

TriggerConsumerDefinition(id, parent) = TriggerConsumerDefinition(
    _as_trigger_consumer_id(id), _require_trigger_parent(parent))

@inline trigger_source_id(definition::TriggerSourceDefinition) = definition.id
@inline trigger_link_id(definition::TriggerLinkDefinition) = definition.id
@inline trigger_consumer_id(definition::TriggerConsumerDefinition) =
    definition.id
@inline trigger_parent_id(definition::Union{
    TriggerLinkDefinition,TriggerConsumerDefinition}) = definition.parent

struct TriggerFaultSet
    domain::UInt64
    bits::UInt64
end

Base.isempty(faults::TriggerFaultSet) = iszero(faults.bits)
function Base.union(left::TriggerFaultSet, right::TriggerFaultSet)
    left.domain == right.domain || _trigger_topology_error(:foreign_fault_set,
        "trigger fault sets belong to different prepared fault domains")
    return TriggerFaultSet(left.domain, left.bits | right.bits)
end

struct _PreparedTriggerFaultSample
    sequence::UInt64
    jitter::PlantTimeOffset
    phase_step::PlantTimeOffset
    timestamp_label_offset::PlantTimeOffset
    action::TriggerEdgeAction
    duplicate_delay::PlantDuration
    faults::TriggerFaultSet
    fault_slot::UInt8
end

@inline _no_trigger_fault_sample(domain::UInt64) =
    _PreparedTriggerFaultSample(UInt64(0), zero(PlantTimeOffset),
        zero(PlantTimeOffset), zero(PlantTimeOffset), PassTriggerEdge,
        zero(PlantDuration), TriggerFaultSet(domain, UInt64(0)), UInt8(0))

struct _PreparedTriggerSource
    id::TriggerSourceID
    schedule::PeriodicSchedule
    origin::PlantTimestamp
    timing_offset::PlantTimeOffset
    timestamp_label_offset::PlantTimeOffset
    trace_start::Int
    trace_count::Int
    ordinal::UInt32
end

@enum _TriggerParentKind::UInt8 begin
    _TriggerSourceParent = 0x01
    _TriggerLinkParent = 0x02
end

struct _PreparedTriggerLink
    id::TriggerLinkID
    parent_kind::_TriggerParentKind
    parent_slot::UInt32
    root_source_slot::UInt32
    propagation_delay::PlantDuration
    timing_skew::PlantTimeOffset
    timestamp_label_offset::PlantTimeOffset
    trace_start::Int
    trace_count::Int
    ordinal::UInt32
end

struct _PreparedTriggerConsumer
    id::TriggerConsumerID
    parent_kind::_TriggerParentKind
    parent_slot::UInt32
    root_source_slot::UInt32
    ordinal::UInt32
end

mutable struct _TriggerTopologyBinding end

"""Run-immutable finite trigger graph and fixed delivery-capacity contract."""
struct PreparedTriggerTopology
    binding::_TriggerTopologyBinding
    sources::Memory{_PreparedTriggerSource}
    links::Memory{_PreparedTriggerLink}
    consumers::Memory{_PreparedTriggerConsumer}
    link_child_offsets::Memory{Int}
    link_child_slots::Memory{UInt32}
    consumer_child_offsets::Memory{Int}
    consumer_child_slots::Memory{UInt32}
    trace_entries::Memory{_PreparedTriggerFaultSample}
    fault_ids::Memory{TriggerFaultID}
    fault_domain::UInt64
    in_flight_capacity::Int
    required_in_flight_capacity::Int
    propagation_capacity::Int
    realization_capacity::Int
    fault_observation_capacity::Int
end

@inline _trigger_topology_binding_id(topology::PreparedTriggerTopology) =
    UInt64(objectid(getfield(topology, :binding)))
@inline trigger_source_count(topology::PreparedTriggerTopology) =
    length(topology.sources)
@inline trigger_link_count(topology::PreparedTriggerTopology) =
    length(topology.links)
@inline trigger_consumer_count(topology::PreparedTriggerTopology) =
    length(topology.consumers)
@inline trigger_in_flight_capacity(topology::PreparedTriggerTopology) =
    topology.in_flight_capacity
@inline required_trigger_in_flight_capacity(
    topology::PreparedTriggerTopology) = topology.required_in_flight_capacity

struct TriggerSourceHandle
    binding_id::UInt64
    slot::UInt32
end

function trigger_source_handle(topology::PreparedTriggerTopology,
    id::TriggerSourceID)
    @inbounds for slot in eachindex(topology.sources)
        topology.sources[slot].id == id && return TriggerSourceHandle(
            _trigger_topology_binding_id(topology), UInt32(slot))
    end
    _trigger_topology_error(:unknown_source,
        "prepared trigger topology has no source $id")
end

trigger_source_handle(topology::PreparedTriggerTopology, name::Symbol) =
    trigger_source_handle(topology, TriggerSourceID(name))

@inline function _require_trigger_source_handle(
    topology::PreparedTriggerTopology, handle::TriggerSourceHandle)
    handle.binding_id == _trigger_topology_binding_id(topology) ||
        _trigger_topology_error(:foreign_source,
            "trigger source handle belongs to another prepared topology")
    slot = Int(handle.slot)
    1 <= slot <= length(topology.sources) || _trigger_topology_error(
        :invalid_source, "trigger source handle contains an invalid slot")
    return slot
end

struct NominalTriggerEdge
    source_id::TriggerSourceID
    timestamp::PlantTimestamp
    sequence::UInt64
end

struct DeliveredTriggerEdge
    timestamp::PlantTimestamp
    occurrence::UInt64
    duplicate::Bool
end

struct ReportedTriggerTimestamp
    timestamp::PlantTimestamp
end

struct TriggerDelivery
    consumer_id::TriggerConsumerID
    nominal::NominalTriggerEdge
    delivered::DeliveredTriggerEdge
    reported::ReportedTriggerTimestamp
    faults::TriggerFaultSet
    consumer_ordinal::UInt32
end

@inline nominal_trigger_edge(delivery::TriggerDelivery) = delivery.nominal
@inline delivered_trigger_edge(delivery::TriggerDelivery) = delivery.delivered
@inline reported_trigger_timestamp(delivery::TriggerDelivery) =
    delivery.reported
@inline trigger_delivery_faults(delivery::TriggerDelivery) = delivery.faults
@inline trigger_delivery_consumer(delivery::TriggerDelivery) =
    delivery.consumer_id

struct TriggerSourcePreview
    handle::TriggerSourceHandle
    nominal::NominalTriggerEdge
    realized_timestamp::PlantTimestamp
    label_offset::PlantTimeOffset
    action::TriggerEdgeAction
    duplicate_delay::PlantDuration
    faults::TriggerFaultSet
    fault_slot::UInt8
end

struct TriggerSourceRealization
    preview::TriggerSourcePreview
    delivery_count::Int
end

@inline nominal_trigger_edge(preview::TriggerSourcePreview) = preview.nominal
@inline nominal_trigger_edge(realization::TriggerSourceRealization) =
    realization.preview.nominal
@inline realized_trigger_source_timestamp(preview::TriggerSourcePreview) =
    preview.realized_timestamp
@inline realized_trigger_source_timestamp(realization::TriggerSourceRealization) =
    realization.preview.realized_timestamp
@inline realized_trigger_delivery_count(realization::TriggerSourceRealization) =
    realization.delivery_count

struct TriggerFaultObservation
    binding_id::UInt64
    source_slot::UInt32
    link_slot::UInt32
    source_sequence::UInt64
    fault_slot::UInt8
    action::TriggerEdgeAction
end

struct _TriggerPropagationEdge
    parent_kind::_TriggerParentKind
    parent_slot::UInt32
    timestamp::PlantTimestamp
    label_offset::PlantTimeOffset
    faults::TriggerFaultSet
    duplicate::Bool
end

mutable struct TriggerTopologyState
    binding_id::UInt64
    source_sequences::Memory{UInt64}
    source_phase_offsets::Memory{PlantTimeOffset}
    link_phase_offsets::Memory{PlantTimeOffset}
    last_ordinary_delivery::Memory{PlantTimestamp}
    has_last_ordinary_delivery::Memory{Bool}
    pending::Memory{TriggerDelivery}
    pending_active::Memory{Bool}
    pending_count::Int
end

function TriggerTopologyState(topology::PreparedTriggerTopology)
    source_count = length(topology.sources)
    link_count = length(topology.links)
    consumer_count = length(topology.consumers)
    sequences = Memory{UInt64}(undef, source_count)
    source_offsets = Memory{PlantTimeOffset}(undef, source_count)
    link_offsets = Memory{PlantTimeOffset}(undef, link_count)
    last_delivery = Memory{PlantTimestamp}(undef, consumer_count)
    has_last_delivery = Memory{Bool}(undef, consumer_count)
    fill!(sequences, UInt64(1))
    fill!(source_offsets, zero(PlantTimeOffset))
    fill!(link_offsets, zero(PlantTimeOffset))
    fill!(last_delivery, zero(PlantTimestamp))
    fill!(has_last_delivery, false)
    pending = Memory{TriggerDelivery}(undef, topology.in_flight_capacity)
    pending_active = Memory{Bool}(undef, topology.in_flight_capacity)
    fill!(pending_active, false)
    return TriggerTopologyState(_trigger_topology_binding_id(topology),
        sequences, source_offsets, link_offsets, last_delivery,
        has_last_delivery, pending, pending_active, 0,
    )
end

mutable struct TriggerTopologyWorkspace
    binding_id::UInt64
    propagation::Memory{_TriggerPropagationEdge}
    propagation_count::Int
    staged_deliveries::Memory{TriggerDelivery}
    staged_delivery_count::Int
    consumer_occurrences::Memory{UInt64}
    next_link_phase_offsets::Memory{PlantTimeOffset}
    link_samples::Memory{_PreparedTriggerFaultSample}
    link_selected::Memory{Bool}
    next_last_ordinary_delivery::Memory{PlantTimestamp}
    next_has_last_ordinary_delivery::Memory{Bool}
    fault_observations::Memory{TriggerFaultObservation}
    fault_observation_count::Int
end

function TriggerTopologyWorkspace(topology::PreparedTriggerTopology)
    consumer_count = length(topology.consumers)
    link_count = length(topology.links)
    return TriggerTopologyWorkspace(
        _trigger_topology_binding_id(topology),
        Memory{_TriggerPropagationEdge}(undef,
            topology.propagation_capacity),
        0,
        Memory{TriggerDelivery}(undef, topology.realization_capacity),
        0,
        Memory{UInt64}(undef, consumer_count),
        Memory{PlantTimeOffset}(undef, link_count),
        Memory{_PreparedTriggerFaultSample}(undef, link_count),
        Memory{Bool}(undef, link_count),
        Memory{PlantTimestamp}(undef, consumer_count),
        Memory{Bool}(undef, consumer_count),
        Memory{TriggerFaultObservation}(undef,
            topology.fault_observation_capacity),
        0,
    )
end

@inline function _require_trigger_binding(topology::PreparedTriggerTopology,
    state::TriggerTopologyState)
    state.binding_id == _trigger_topology_binding_id(topology) ||
        _trigger_topology_error(:foreign_state,
            "trigger topology state belongs to another prepared topology")
    return nothing
end

@inline function _require_trigger_binding(topology::PreparedTriggerTopology,
    workspace::TriggerTopologyWorkspace)
    workspace.binding_id == _trigger_topology_binding_id(topology) ||
        _trigger_topology_error(:foreign_workspace,
            "trigger topology workspace belongs to another prepared topology")
    return nothing
end

@inline function _trigger_trace_sample(topology::PreparedTriggerTopology,
    start::Int, count::Int, sequence::UInt64)
    stop = start + count - 1
    @inbounds for index in start:stop
        sample = topology.trace_entries[index]
        sample.sequence == sequence && return sample
        sample.sequence > sequence && break
    end
    return _no_trigger_fault_sample(topology.fault_domain)
end

@inline function _trigger_fault_set(sample::_PreparedTriggerFaultSample)
    return sample.faults
end

function contains_trigger_fault(topology::PreparedTriggerTopology,
    faults::TriggerFaultSet, id::TriggerFaultID)
    faults.domain == topology.fault_domain || _trigger_topology_error(
        :foreign_fault_set,
        "trigger fault set belongs to another prepared fault domain")
    @inbounds for slot in eachindex(topology.fault_ids)
        if topology.fault_ids[slot] == id
            return (faults.bits & (UInt64(1) << (slot - 1))) != 0
        end
    end
    _trigger_topology_error(:unknown_fault,
        "prepared trigger topology has no fault $id")
end

contains_trigger_fault(topology::PreparedTriggerTopology,
    faults::TriggerFaultSet, name::Symbol) = contains_trigger_fault(
        topology, faults, TriggerFaultID(name))

@inline function _source_preview(topology::PreparedTriggerTopology,
    state::TriggerTopologyState, slot::Int)
    source = topology.sources[slot]
    sequence = state.source_sequences[slot]
    sample = _trigger_trace_sample(topology, source.trace_start,
        source.trace_count, sequence)
    phase_offset = state.source_phase_offsets[slot] + sample.phase_step
    realized_offset = source.timing_offset + phase_offset + sample.jitter
    nominal_timestamp = schedule_timestamp(source.schedule, sequence,
        source.origin)
    realized_timestamp = nominal_timestamp + realized_offset
    label_offset = source.timestamp_label_offset +
        sample.timestamp_label_offset
    nominal = NominalTriggerEdge(source.id, nominal_timestamp, sequence)
    return TriggerSourcePreview(TriggerSourceHandle(
        _trigger_topology_binding_id(topology), UInt32(slot)), nominal,
        realized_timestamp, label_offset, sample.action,
        sample.duplicate_delay, sample.faults, sample.fault_slot)
end

@inline function _source_preview_precedes(left::TriggerSourcePreview,
    right::TriggerSourcePreview)
    left.realized_timestamp == right.realized_timestamp ||
        return left.realized_timestamp < right.realized_timestamp
    return left.handle.slot < right.handle.slot
end

"""Return the next realized source edge without advancing trigger state."""
function next_trigger_source(topology::PreparedTriggerTopology,
    state::TriggerTopologyState)
    _require_trigger_binding(topology, state)
    preview = _source_preview(topology, state, 1)
    @inbounds for slot in 2:length(topology.sources)
        candidate = _source_preview(topology, state, slot)
        _source_preview_precedes(candidate, preview) && (preview = candidate)
    end
    return preview
end

@inline function _delivery_precedes(left::TriggerDelivery,
    right::TriggerDelivery)
    left.delivered.timestamp == right.delivered.timestamp ||
        return left.delivered.timestamp < right.delivered.timestamp
    left.consumer_ordinal == right.consumer_ordinal ||
        return left.consumer_ordinal < right.consumer_ordinal
    left.nominal.sequence == right.nominal.sequence ||
        return left.nominal.sequence < right.nominal.sequence
    return left.delivered.occurrence < right.delivered.occurrence
end

@inline function _next_pending_slot(state::TriggerTopologyState)
    best = 0
    @inbounds for slot in eachindex(state.pending_active)
        state.pending_active[slot] || continue
        if iszero(best) || _delivery_precedes(state.pending[slot],
                state.pending[best])
            best = slot
        end
    end
    return best
end

"""Return the earliest pending delivered edge without removing it."""
function next_trigger_delivery(topology::PreparedTriggerTopology,
    state::TriggerTopologyState)
    _require_trigger_binding(topology, state)
    slot = _next_pending_slot(state)
    iszero(slot) && return nothing
    return state.pending[slot]
end

"""Return the earliest pending delivery timestamp without copying its record."""
function next_trigger_delivery_timestamp(topology::PreparedTriggerTopology,
    state::TriggerTopologyState)
    _require_trigger_binding(topology, state)
    slot = _next_pending_slot(state)
    iszero(slot) && return nothing
    return state.pending[slot].delivered.timestamp
end

@inline pending_trigger_delivery_count(state::TriggerTopologyState) =
    state.pending_count

@inline function _reset_trigger_workspace!(workspace::TriggerTopologyWorkspace,
    state::TriggerTopologyState)
    workspace.propagation_count = 0
    workspace.staged_delivery_count = 0
    workspace.fault_observation_count = 0
    fill!(workspace.consumer_occurrences, UInt64(0))
    copyto!(workspace.next_link_phase_offsets, state.link_phase_offsets)
    fill!(workspace.link_selected, false)
    copyto!(workspace.next_last_ordinary_delivery,
        state.last_ordinary_delivery)
    copyto!(workspace.next_has_last_ordinary_delivery,
        state.has_last_ordinary_delivery)
    return nothing
end

@inline function _stage_fault_observation!(
    workspace::TriggerTopologyWorkspace, source_slot::Int, link_slot::Int,
    sequence::UInt64, sample::_PreparedTriggerFaultSample)
    iszero(sample.fault_slot) && return nothing
    next_count = workspace.fault_observation_count + 1
    next_count <= length(workspace.fault_observations) ||
        _trigger_topology_error(:realization_overflow,
            "trigger fault-observation workspace capacity was exceeded")
    workspace.fault_observations[next_count] = TriggerFaultObservation(
        workspace.binding_id, UInt32(source_slot), UInt32(link_slot), sequence,
        sample.fault_slot, sample.action)
    workspace.fault_observation_count = next_count
    return nothing
end

@inline function _enqueue_propagation!(workspace::TriggerTopologyWorkspace,
    edge::_TriggerPropagationEdge)
    next_count = workspace.propagation_count + 1
    next_count <= length(workspace.propagation) ||
        _trigger_topology_error(:realization_overflow,
            "trigger propagation workspace capacity was exceeded")
    workspace.propagation[next_count] = edge
    workspace.propagation_count = next_count
    return nothing
end

@inline function _next_trigger_occurrence!(
    workspace::TriggerTopologyWorkspace, consumer_slot::Int)
    occurrence = workspace.consumer_occurrences[consumer_slot]
    occurrence != typemax(UInt64) || _trigger_topology_error(
        :occurrence_overflow,
        "trigger delivery occurrence exceeds UInt64 range")
    occurrence += UInt64(1)
    workspace.consumer_occurrences[consumer_slot] = occurrence
    return occurrence
end

@inline function _stage_trigger_delivery!(
    workspace::TriggerTopologyWorkspace,
    consumer::_PreparedTriggerConsumer, consumer_slot::Int,
    nominal::NominalTriggerEdge, edge::_TriggerPropagationEdge)
    occurrence = _next_trigger_occurrence!(workspace, consumer_slot)
    reported_timestamp = edge.timestamp + edge.label_offset
    delivery = TriggerDelivery(consumer.id, nominal,
        DeliveredTriggerEdge(edge.timestamp, occurrence, edge.duplicate),
        ReportedTriggerTimestamp(reported_timestamp), edge.faults,
        consumer.ordinal)
    if !edge.duplicate
        if workspace.next_has_last_ordinary_delivery[consumer_slot] &&
                !(workspace.next_last_ordinary_delivery[consumer_slot] <
                    edge.timestamp)
            _trigger_topology_error(:ordinary_edge_overtaking,
                "ordinary trigger deliveries must remain strictly ordered " *
                "by source sequence for consumer $(consumer.id)")
        end
        workspace.next_last_ordinary_delivery[consumer_slot] = edge.timestamp
        workspace.next_has_last_ordinary_delivery[consumer_slot] = true
    end
    next_count = workspace.staged_delivery_count + 1
    next_count <= length(workspace.staged_deliveries) ||
        _trigger_topology_error(:realization_overflow,
            "trigger delivery realization capacity was exceeded")
    workspace.staged_deliveries[next_count] = delivery
    workspace.staged_delivery_count = next_count
    return nothing
end

@inline function _trigger_node_slot(topology::PreparedTriggerTopology,
    kind::_TriggerParentKind, slot::UInt32)
    kind == _TriggerSourceParent && return Int(slot)
    return length(topology.sources) + Int(slot)
end

function _prepare_link_sequence_state!(workspace::TriggerTopologyWorkspace,
    topology::PreparedTriggerTopology, state::TriggerTopologyState,
    source_slot::Int, sequence::UInt64)
    @inbounds for link_slot in eachindex(topology.links)
        link = topology.links[link_slot]
        link.root_source_slot == UInt32(source_slot) || continue
        sample = _trigger_trace_sample(topology, link.trace_start,
            link.trace_count, sequence)
        workspace.link_selected[link_slot] = true
        workspace.link_samples[link_slot] = sample
        workspace.next_link_phase_offsets[link_slot] =
            state.link_phase_offsets[link_slot] + sample.phase_step
        _stage_fault_observation!(workspace, source_slot, link_slot,
            sequence, sample)
    end
    return nothing
end

@inline function _link_output_timestamp(link::_PreparedTriggerLink,
    phase_offset::PlantTimeOffset,
    sample::_PreparedTriggerFaultSample,
    input_timestamp::PlantTimestamp)
    delayed = input_timestamp + link.propagation_delay
    return delayed + link.timing_skew + phase_offset + sample.jitter
end

function _propagate_trigger_source!(workspace::TriggerTopologyWorkspace,
    topology::PreparedTriggerTopology, preview::TriggerSourcePreview,
    source_slot::Int)
    preview.action == DropTriggerEdge && return nothing
    root_faults = preview.faults
    root = _TriggerPropagationEdge(_TriggerSourceParent,
        UInt32(source_slot), preview.realized_timestamp,
        preview.label_offset, root_faults, false)
    _enqueue_propagation!(workspace, root)
    if preview.action == DuplicateTriggerEdge
        duplicate = _TriggerPropagationEdge(_TriggerSourceParent,
            UInt32(source_slot),
            preview.realized_timestamp + preview.duplicate_delay,
            preview.label_offset, root_faults, true)
        _enqueue_propagation!(workspace, duplicate)
    end

    cursor = 1
    while cursor <= workspace.propagation_count
        edge = workspace.propagation[cursor]
        node_slot = _trigger_node_slot(topology, edge.parent_kind,
            edge.parent_slot)
        consumer_start = topology.consumer_child_offsets[node_slot]
        consumer_stop = topology.consumer_child_offsets[node_slot + 1] - 1
        @inbounds for child_index in consumer_start:consumer_stop
            consumer_slot = Int(topology.consumer_child_slots[child_index])
            consumer = topology.consumers[consumer_slot]
            _stage_trigger_delivery!(workspace, consumer, consumer_slot,
                preview.nominal, edge)
        end
        link_start = topology.link_child_offsets[node_slot]
        link_stop = topology.link_child_offsets[node_slot + 1] - 1
        @inbounds for child_index in link_start:link_stop
            link_slot = Int(topology.link_child_slots[child_index])
            link = topology.links[link_slot]
            sample = workspace.link_samples[link_slot]
            sample.action == DropTriggerEdge && continue
            output_timestamp = _link_output_timestamp(link,
                workspace.next_link_phase_offsets[link_slot], sample,
                edge.timestamp)
            label_offset = edge.label_offset +
                link.timestamp_label_offset + sample.timestamp_label_offset
            faults = union(edge.faults, sample.faults)
            output = _TriggerPropagationEdge(_TriggerLinkParent,
                UInt32(link_slot), output_timestamp, label_offset, faults,
                edge.duplicate)
            _enqueue_propagation!(workspace, output)
            if sample.action == DuplicateTriggerEdge
                duplicate = _TriggerPropagationEdge(_TriggerLinkParent,
                    UInt32(link_slot),
                    output_timestamp + sample.duplicate_delay,
                    label_offset, faults, true)
                _enqueue_propagation!(workspace, duplicate)
            end
        end
        cursor += 1
    end
    return nothing
end

@inline function _free_pending_slots(state::TriggerTopologyState)
    return length(state.pending) - state.pending_count
end

function _commit_staged_deliveries!(state::TriggerTopologyState,
    workspace::TriggerTopologyWorkspace)
    workspace.staged_delivery_count <= _free_pending_slots(state) ||
        _trigger_topology_error(:in_flight_capacity,
            "realized trigger deliveries exceed the declared in-flight " *
            "capacity; drain due deliveries before advancing sources")
    staged_slot = 1
    @inbounds for pending_slot in eachindex(state.pending_active)
        staged_slot > workspace.staged_delivery_count && break
        state.pending_active[pending_slot] && continue
        state.pending[pending_slot] = workspace.staged_deliveries[staged_slot]
        state.pending_active[pending_slot] = true
        staged_slot += 1
    end
    state.pending_count += workspace.staged_delivery_count
    return nothing
end

"""
    realize_next_trigger_source!(workspace, topology, state)

Realize the globally earliest source edge, map it through the prepared fan-out,
and append its delivered consumer edges to fixed pending storage. A delivery
strictly earlier than the next source must be popped first; equal-time sources
precede deliveries.
"""
function realize_next_trigger_source!(workspace::TriggerTopologyWorkspace,
    topology::PreparedTriggerTopology, state::TriggerTopologyState)
    _require_trigger_binding(topology, state)
    _require_trigger_binding(topology, workspace)
    preview = next_trigger_source(topology, state)
    pending_slot = _next_pending_slot(state)
    if !iszero(pending_slot) &&
            state.pending[pending_slot].delivered.timestamp <
                preview.realized_timestamp
        _trigger_topology_error(:delivery_due,
            "pop the earlier pending trigger delivery before advancing the " *
            "next trigger source")
    end
    source_slot = _require_trigger_source_handle(topology, preview.handle)
    sequence = preview.nominal.sequence
    source = topology.sources[source_slot]
    source_sample = _trigger_trace_sample(topology, source.trace_start,
        source.trace_count, sequence)
    _reset_trigger_workspace!(workspace, state)
    _stage_fault_observation!(workspace, source_slot, 0, sequence,
        source_sample)
    _prepare_link_sequence_state!(workspace, topology, state, source_slot,
        sequence)
    _propagate_trigger_source!(workspace, topology, preview, source_slot)

    workspace.staged_delivery_count <= _free_pending_slots(state) ||
        _trigger_topology_error(:in_flight_capacity,
            "realized trigger deliveries exceed the declared in-flight " *
            "capacity; drain due deliveries before advancing sources")
    next_sequence = _next_trigger_sequence(sequence)
    state.source_sequences[source_slot] = next_sequence
    state.source_phase_offsets[source_slot] =
        state.source_phase_offsets[source_slot] + source_sample.phase_step
    @inbounds for link_slot in eachindex(topology.links)
        workspace.link_selected[link_slot] || continue
        state.link_phase_offsets[link_slot] =
            workspace.next_link_phase_offsets[link_slot]
    end
    copyto!(state.last_ordinary_delivery,
        workspace.next_last_ordinary_delivery)
    copyto!(state.has_last_ordinary_delivery,
        workspace.next_has_last_ordinary_delivery)
    _commit_staged_deliveries!(state, workspace)
    return TriggerSourceRealization(preview,
        workspace.staged_delivery_count)
end

"""
Remove and return the earliest delivered trigger edge. A source edge earlier
than or equal to that delivery must be realized first.
"""
function pop_next_trigger_delivery!(topology::PreparedTriggerTopology,
    state::TriggerTopologyState)
    _require_trigger_binding(topology, state)
    slot = _next_pending_slot(state)
    iszero(slot) && return nothing
    delivery = state.pending[slot]
    source = next_trigger_source(topology, state)
    source.realized_timestamp <= delivery.delivered.timestamp &&
        _trigger_topology_error(:source_due,
            "realize the earlier or equal-time trigger source before popping " *
            "this delivery")
    state.pending_active[slot] = false
    state.pending_count -= 1
    return delivery
end

"""
    pop_next_trigger_delivery!(output, topology, state) -> Bool

Allocation-free hot-path form. Copy the earliest delivered edge into the
caller-owned `Ref{TriggerDelivery}` and return `true`, or return `false` when
none is pending. Source-before-delivery chronology is identical to the
record-returning convenience method.
"""
function pop_next_trigger_delivery!(output::Base.RefValue{TriggerDelivery},
    topology::PreparedTriggerTopology, state::TriggerTopologyState)
    _require_trigger_binding(topology, state)
    slot = _next_pending_slot(state)
    iszero(slot) && return false
    delivery = state.pending[slot]
    source = next_trigger_source(topology, state)
    source.realized_timestamp <= delivery.delivered.timestamp &&
        _trigger_topology_error(:source_due,
            "realize the earlier or equal-time trigger source before popping " *
            "this delivery")
    output[] = delivery
    state.pending_active[slot] = false
    state.pending_count -= 1
    return true
end

@inline trigger_fault_observation_count(
    workspace::TriggerTopologyWorkspace) = workspace.fault_observation_count

function trigger_fault_observation(workspace::TriggerTopologyWorkspace,
    index::Integer)
    1 <= index <= workspace.fault_observation_count ||
        _trigger_topology_error(:invalid_observation_index,
            "trigger fault-observation index is outside the latest " *
            "realization")
    return workspace.fault_observations[Int(index)]
end

trigger_fault_observation(::TriggerTopologyWorkspace, ::Bool) =
    _trigger_topology_error(:invalid_observation_index,
        "trigger fault-observation index must be an integer count, not Bool")

function trigger_fault_id(topology::PreparedTriggerTopology,
    observation::TriggerFaultObservation)
    observation.binding_id == _trigger_topology_binding_id(topology) ||
        _trigger_topology_error(:foreign_observation,
            "trigger fault observation belongs to another topology")
    slot = Int(observation.fault_slot)
    1 <= slot <= length(topology.fault_ids) || _trigger_topology_error(
        :foreign_observation,
        "trigger fault observation does not belong to this topology")
    return topology.fault_ids[slot]
end

function trigger_fault_location(topology::PreparedTriggerTopology,
    observation::TriggerFaultObservation)
    observation.binding_id == _trigger_topology_binding_id(topology) ||
        _trigger_topology_error(:foreign_observation,
            "trigger fault observation belongs to another topology")
    source_slot = Int(observation.source_slot)
    1 <= source_slot <= length(topology.sources) || _trigger_topology_error(
        :foreign_observation,
        "trigger fault observation does not belong to this topology")
    if iszero(observation.link_slot)
        return topology.sources[source_slot].id
    end
    link_slot = Int(observation.link_slot)
    1 <= link_slot <= length(topology.links) || _trigger_topology_error(
        :foreign_observation,
        "trigger fault observation does not belong to this topology")
    return topology.links[link_slot].id
end

# -- cold preparation -------------------------------------------------------

@inline _require_trigger_definition(value::TriggerSourceDefinition,
    ::Type{TriggerSourceDefinition}) = value
@inline _require_trigger_definition(value::TriggerLinkDefinition,
    ::Type{TriggerLinkDefinition}) = value
@inline _require_trigger_definition(value::TriggerConsumerDefinition,
    ::Type{TriggerConsumerDefinition}) = value

function _require_trigger_definition(value, ::Type{T}) where {T}
    _trigger_topology_error(:invalid_definition,
        "trigger topology expected $T entries; got $(typeof(value))")
end

function _collect_trigger_definition_values(source, ::Type{T}) where {T}
    result = Any[]
    sizehint!(result, length(source))
    for value in source
        push!(result, _require_trigger_definition(value, T))
    end
    sort!(result; by=value -> String(getfield(value, :id).name))
    @inbounds for index in 2:length(result)
        getfield(result[index - 1], :id) == getfield(result[index], :id) &&
            _trigger_topology_error(:duplicate_id,
                "trigger identities must be unique within each registry; " *
                "duplicate $(getfield(result[index], :id))")
    end
    return result
end

_collect_trigger_definitions(values::Tuple, type::Type) =
    _collect_trigger_definition_values(values, type)

_collect_trigger_definitions(values::AbstractVector, type::Type) =
    _collect_trigger_definition_values(values, type)

function _collect_trigger_definitions(values::NamedTuple, type::Type)
    for (name, value) in pairs(values)
        checked = _require_trigger_definition(value, type)
        name == getfield(checked, :id).name || _trigger_topology_error(
            :identity_mismatch,
            "named trigger key $(repr(name)) does not match " *
            "$(getfield(checked, :id))")
    end
    return _collect_trigger_definition_values(Base.values(values), type)
end

function _collect_trigger_definitions(values, ::Type)
    _trigger_topology_error(:invalid_registry,
        "trigger topology registries must be finite Tuple, NamedTuple, " *
        "or AbstractVector values; got $(typeof(values))")
end

@inline function _find_trigger_id(definitions, id)
    @inbounds for index in eachindex(definitions)
        getfield(definitions[index], :id) == id && return index
    end
    return 0
end

@inline function _resolve_trigger_parent(parent::TriggerSourceID, sources,
    prepared_links)
    slot = _find_trigger_id(sources, parent)
    return _TriggerSourceParent, slot, slot
end

@inline function _resolve_trigger_parent(parent::TriggerLinkID, sources,
    prepared_links)
    slot = _find_trigger_id(prepared_links, parent)
    root = iszero(slot) ? 0 :
        Int(prepared_links[slot].root_source_slot)
    return _TriggerLinkParent, slot, root
end

@inline _unresolved_trigger_parent_is_cycle(::TriggerSourceID, links) = false
@inline _unresolved_trigger_parent_is_cycle(parent::TriggerLinkID, links) =
    !iszero(_find_trigger_id(links, parent))

@inline function _prepared_parent_node_slot(source_count::Int,
    kind::_TriggerParentKind, slot::UInt32)
    kind == _TriggerSourceParent && return Int(slot)
    return source_count + Int(slot)
end

function _prepare_trigger_child_registry(source_count::Int, link_count::Int,
    links::Memory{_PreparedTriggerLink},
    consumers::Memory{_PreparedTriggerConsumer})
    node_count = source_count + link_count
    link_counts = zeros(Int, node_count)
    consumer_counts = zeros(Int, node_count)
    @inbounds for link in links
        node = _prepared_parent_node_slot(source_count, link.parent_kind,
            link.parent_slot)
        link_counts[node] += 1
    end
    @inbounds for consumer in consumers
        node = _prepared_parent_node_slot(source_count, consumer.parent_kind,
            consumer.parent_slot)
        consumer_counts[node] += 1
    end

    link_offsets = Memory{Int}(undef, node_count + 1)
    consumer_offsets = Memory{Int}(undef, node_count + 1)
    link_offsets[1] = 1
    consumer_offsets[1] = 1
    @inbounds for node in 1:node_count
        link_offsets[node + 1] = link_offsets[node] + link_counts[node]
        consumer_offsets[node + 1] =
            consumer_offsets[node] + consumer_counts[node]
    end
    link_positions = collect(link_offsets[1:node_count])
    consumer_positions = collect(consumer_offsets[1:node_count])
    link_slots = Memory{UInt32}(undef, link_count)
    consumer_slots = Memory{UInt32}(undef, length(consumers))
    @inbounds for (slot, link) in enumerate(links)
        node = _prepared_parent_node_slot(source_count, link.parent_kind,
            link.parent_slot)
        position = link_positions[node]
        link_slots[position] = UInt32(slot)
        link_positions[node] = position + 1
    end
    @inbounds for (slot, consumer) in enumerate(consumers)
        node = _prepared_parent_node_slot(source_count, consumer.parent_kind,
            consumer.parent_slot)
        position = consumer_positions[node]
        consumer_slots[position] = UInt32(slot)
        consumer_positions[node] = position + 1
    end
    return link_offsets, link_slots, consumer_offsets, consumer_slots
end

function _collect_trigger_fault_ids(sources, links)
    ids = TriggerFaultID[]
    for definition in sources
        append!(ids, (entry.fault_id for entry in definition.faults.entries))
    end
    for definition in links
        append!(ids, (entry.fault_id for entry in definition.faults.entries))
    end
    sort!(ids; by=id -> String(id.name))
    unique!(ids)
    length(ids) <= 64 || _trigger_topology_error(:fault_capacity,
        "prepared trigger topology supports at most 64 distinct fault " *
        "identities")
    return ids
end

@inline function _trigger_domain_byte(state::UInt64, byte::UInt8)
    return (state ⊻ UInt64(byte)) * UInt64(0x00000100000001b3)
end

function _trigger_fault_domain(fault_ids)
    state = UInt64(0xcbf29ce484222325)
    for byte in codeunits("AdaptiveOpticsSim.trigger_fault_domain.v1")
        state = _trigger_domain_byte(state, byte)
    end
    @inbounds for id in fault_ids
        bytes = codeunits(String(id.name))
        for byte in bytes
            state = _trigger_domain_byte(state, byte)
        end
        state = _trigger_domain_byte(state, UInt8(0))
    end
    return state
end

@inline function _trigger_fault_slot(fault_ids, id::TriggerFaultID)
    @inbounds for slot in eachindex(fault_ids)
        fault_ids[slot] == id && return slot
    end
    _trigger_topology_error(:unknown_fault,
        "internal trigger preparation lost fault identity $id")
end

function _append_trigger_trace!(prepared_entries, fault_ids,
    fault_domain::UInt64, trace::TriggerFaultTrace)
    entries = collect(trace.entries)
    sort!(entries; by=entry -> entry.sequence)
    @inbounds for index in 2:length(entries)
        entries[index - 1].sequence == entries[index].sequence &&
            _trigger_topology_error(:duplicate_fault_sequence,
                "one source or link may declare at most one fault-trace " *
                "entry per source sequence")
    end
    start = length(prepared_entries) + 1
    for entry in entries
        fault_slot = _trigger_fault_slot(fault_ids, entry.fault_id)
        bits = UInt64(1) << (fault_slot - 1)
        push!(prepared_entries, _PreparedTriggerFaultSample(
            entry.sequence, entry.jitter, entry.phase_step,
            entry.timestamp_label_offset, entry.action,
            entry.duplicate_delay, TriggerFaultSet(fault_domain, bits),
            UInt8(fault_slot)))
    end
    return start, length(entries)
end

function _checked_trigger_capacity(capacity::Integer)
    capacity >= 0 || _trigger_topology_error(:invalid_capacity,
        "trigger in-flight capacity must be nonnegative")
    capacity <= typemax(Int) || _trigger_topology_error(:invalid_capacity,
        "trigger in-flight capacity exceeds Int range")
    return Int(capacity)
end

_checked_trigger_capacity(::Bool) = _trigger_topology_error(
    :invalid_capacity,
    "trigger in-flight capacity must be an integer count, not Bool")

_checked_trigger_capacity(value) = _trigger_topology_error(
    :invalid_capacity,
    "trigger in-flight capacity must be an integer count; got " *
    "$(typeof(value))")

@inline function _trace_has_duplicate(topology::PreparedTriggerTopology,
    start::Int, count::Int)
    stop = start + count - 1
    @inbounds for index in start:stop
        topology.trace_entries[index].action == DuplicateTriggerEdge &&
            return true
    end
    return false
end

@inline function _checked_trigger_capacity_add(left::UInt128,
    right::UInt128, label::AbstractString)
    left <= typemax(UInt128) - right || _trigger_topology_error(
        :capacity_overflow, "$label exceeds UInt128 range")
    return left + right
end

@inline function _checked_trigger_capacity_double(value::UInt128,
    label::AbstractString)
    value <= typemax(UInt128) >> 1 || _trigger_topology_error(
        :capacity_overflow, "$label exceeds UInt128 range")
    return value << 1
end

function _trigger_structural_capacities(topology::PreparedTriggerTopology)
    link_multiplicity = fill(UInt128(0), length(topology.links))
    max_propagation = UInt128(0)
    max_realization = UInt128(0)
    max_observations = UInt128(0)
    per_source_outputs = fill(UInt128(0), length(topology.sources))

    @inbounds for source_slot in eachindex(topology.sources)
        source = topology.sources[source_slot]
        root_multiplicity = _trace_has_duplicate(topology,
            source.trace_start, source.trace_count) ? UInt128(2) : UInt128(1)
        propagation = root_multiplicity
        observations = UInt128(1)
        for link_slot in eachindex(topology.links)
            link = topology.links[link_slot]
            link.root_source_slot == UInt32(source_slot) || continue
            parent_multiplicity = link.parent_kind == _TriggerSourceParent ?
                root_multiplicity : link_multiplicity[Int(link.parent_slot)]
            multiplicity = _trace_has_duplicate(topology, link.trace_start,
                link.trace_count) ? _checked_trigger_capacity_double(
                    parent_multiplicity, "trigger propagation multiplicity") :
                parent_multiplicity
            link_multiplicity[link_slot] = multiplicity
            propagation = _checked_trigger_capacity_add(propagation,
                multiplicity, "trigger propagation capacity")
            observations += UInt128(1)
        end
        outputs = UInt128(0)
        for consumer in topology.consumers
            consumer.root_source_slot == UInt32(source_slot) || continue
            multiplicity = consumer.parent_kind == _TriggerSourceParent ?
                root_multiplicity :
                link_multiplicity[Int(consumer.parent_slot)]
            outputs = _checked_trigger_capacity_add(outputs, multiplicity,
                "trigger realization capacity")
        end
        iszero(outputs) && _trigger_topology_error(:unconsumed_source,
            "trigger source $(source.id) has no downstream consumer")
        per_source_outputs[source_slot] = outputs
        max_propagation = max(max_propagation, propagation)
        max_realization = max(max_realization, outputs)
        max_observations = max(max_observations, observations)
    end
    limit = UInt128(typemax(Int))
    max_propagation <= limit && max_realization <= limit &&
        max_observations <= limit || _trigger_topology_error(
            :capacity_overflow,
            "prepared trigger workspace capacity exceeds Int range")
    return Int(max_propagation), Int(max_realization),
        Int(max_observations), per_source_outputs
end

@inline function _trace_phase_at(topology::PreparedTriggerTopology,
    start::Int, count::Int, sequence::UInt64)
    phase = zero(PlantTimeOffset)
    stop = start + count - 1
    @inbounds for index in start:stop
        sample = topology.trace_entries[index]
        sample.sequence > sequence && break
        phase += sample.phase_step
    end
    return phase
end

function _path_candidate_sequences(topology::PreparedTriggerTopology,
    source::_PreparedTriggerSource, path::Vector{Int})
    candidates = UInt64[1, 2]
    function append_trace_candidates(start, count)
        stop = start + count - 1
        @inbounds for index in start:stop
            sequence = topology.trace_entries[index].sequence
            sequence > 1 && push!(candidates, sequence - UInt64(1))
            push!(candidates, sequence)
            sequence < typemax(UInt64) && push!(candidates,
                sequence + UInt64(1))
        end
    end
    append_trace_candidates(source.trace_start, source.trace_count)
    for link_slot in path
        link = topology.links[link_slot]
        append_trace_candidates(link.trace_start, link.trace_count)
    end
    sort!(candidates)
    unique!(candidates)
    return candidates
end

@inline function _source_time_at(topology::PreparedTriggerTopology,
    source::_PreparedTriggerSource, sequence::UInt64)
    sample = _trigger_trace_sample(topology, source.trace_start,
        source.trace_count, sequence)
    phase = _trace_phase_at(topology, source.trace_start,
        source.trace_count, sequence)
    nominal = schedule_timestamp(source.schedule, sequence, source.origin)
    return nominal + source.timing_offset + phase + sample.jitter, sample
end

function _consumer_link_path(topology::PreparedTriggerTopology,
    consumer::_PreparedTriggerConsumer)
    path = Int[]
    kind = consumer.parent_kind
    slot = Int(consumer.parent_slot)
    while kind == _TriggerLinkParent
        push!(path, slot)
        link = topology.links[slot]
        kind = link.parent_kind
        slot = Int(link.parent_slot)
    end
    reverse!(path)
    return path
end

function _validate_trigger_paths_and_capacity(
    topology::PreparedTriggerTopology, per_source_outputs)
    min_intervals = fill(typemax(Int64), length(topology.sources))
    max_latencies = fill(Int64(0), length(topology.sources))
    consumed_links = fill(false, length(topology.links))

    for consumer in topology.consumers
        source_slot = Int(consumer.root_source_slot)
        source = topology.sources[source_slot]
        path = _consumer_link_path(topology, consumer)
        @inbounds for link_slot in path
            consumed_links[link_slot] = true
        end
        candidates = _path_candidate_sequences(topology, source, path)
        previous_source = zero(PlantTimestamp)
        has_previous_source = false
        previous_delivery = zero(PlantTimestamp)
        has_previous_delivery = false
        for sequence in candidates
            source_time, source_sample = _source_time_at(topology, source,
                sequence)
            if has_previous_source
                previous_source < source_time || _trigger_topology_error(
                    :ordinary_edge_overtaking,
                    "ordinary source edges for $(source.id) must be strictly " *
                    "ordered")
                interval = source_time - previous_source
                min_intervals[source_slot] = min(
                    min_intervals[source_slot], plant_nanoseconds(interval))
            end
            previous_source = source_time
            has_previous_source = true
            source_sample.action == DropTriggerEdge && continue

            delivery_time = source_time
            label_offset = source.timestamp_label_offset +
                source_sample.timestamp_label_offset
            dropped = false
            extra_duplicate_delay = source_sample.action ==
                DuplicateTriggerEdge ? source_sample.duplicate_delay :
                zero(PlantDuration)
            for link_slot in path
                link = topology.links[link_slot]
                sample = _trigger_trace_sample(topology, link.trace_start,
                    link.trace_count, sequence)
                phase = _trace_phase_at(topology, link.trace_start,
                    link.trace_count, sequence)
                input_time = delivery_time
                delivery_time = _link_output_timestamp(link, phase, sample,
                    delivery_time)
                delivery_time >= input_time || _trigger_topology_error(
                    :negative_link_delay,
                    "trigger link $(link.id) realizes a negative physical " *
                    "delivery delay")
                label_offset += link.timestamp_label_offset +
                    sample.timestamp_label_offset
                if sample.action == DropTriggerEdge
                    dropped = true
                    break
                elseif sample.action == DuplicateTriggerEdge
                    extra_duplicate_delay += sample.duplicate_delay
                end
            end
            dropped && continue
            delivery_time + label_offset
            if has_previous_delivery && !(previous_delivery < delivery_time)
                _trigger_topology_error(:ordinary_edge_overtaking,
                    "ordinary delivered edges for $(consumer.id) may " *
                    "overtake under the declared trace")
            end
            previous_delivery = delivery_time
            has_previous_delivery = true
            latest_delivery = delivery_time + extra_duplicate_delay
            latency = latest_delivery - source_time
            max_latencies[source_slot] = max(max_latencies[source_slot],
                plant_nanoseconds(latency))
        end
    end

    @inbounds for link_slot in eachindex(consumed_links)
        consumed_links[link_slot] || _trigger_topology_error(
            :unconsumed_link,
            "trigger link $(topology.links[link_slot].id) has no downstream " *
            "consumer")
    end

    required = UInt128(0)
    @inbounds for source_slot in eachindex(topology.sources)
        interval = min_intervals[source_slot]
        interval == typemax(Int64) && (interval = plant_nanoseconds(
            schedule_period(topology.sources[source_slot].schedule)))
        interval > 0 || _trigger_topology_error(:ordinary_edge_overtaking,
            "ordinary source-edge interval must remain positive")
        overlapping = UInt128(max_latencies[source_slot] ÷ interval + 1)
        contribution = overlapping * per_source_outputs[source_slot]
        required = _checked_trigger_capacity_add(required, contribution,
            "required trigger in-flight capacity")
    end
    required <= UInt128(typemax(Int)) || _trigger_topology_error(
        :capacity_overflow,
        "required trigger in-flight capacity exceeds Int range")
    return Int(required)
end

"""
    prepare_trigger_topology(sources, links, consumers;
        in_flight_capacity)

Prepare a finite acyclic fan-out. Declaration order is not identity: sources,
links, consumers, traces, and fault identities are canonicalized. Preparation
validates ordinary non-overtaking delivery and a conservative fixed pending
capacity before runtime state is allocated.
"""
function prepare_trigger_topology(source_values, link_values, consumer_values;
    in_flight_capacity)
    capacity = _checked_trigger_capacity(in_flight_capacity)
    sources = _collect_trigger_definitions(source_values,
        TriggerSourceDefinition)
    links = _collect_trigger_definitions(link_values,
        TriggerLinkDefinition)
    consumers = _collect_trigger_definitions(consumer_values,
        TriggerConsumerDefinition)
    isempty(sources) && _trigger_topology_error(:empty_sources,
        "trigger topology must contain at least one source")
    isempty(consumers) && _trigger_topology_error(:empty_consumers,
        "trigger topology must contain at least one consumer")

    fault_ids = _collect_trigger_fault_ids(sources, links)
    fault_domain = _trigger_fault_domain(fault_ids)
    trace_entries = _PreparedTriggerFaultSample[]
    prepared_sources = Memory{_PreparedTriggerSource}(undef, length(sources))
    @inbounds for (slot, source) in enumerate(sources)
        start, count = _append_trigger_trace!(trace_entries, fault_ids,
            fault_domain, source.faults)
        prepared_sources[slot] = _PreparedTriggerSource(source.id,
            source.schedule, source.origin, source.timing_offset,
            source.timestamp_label_offset, start, count, UInt32(slot))
    end

    # Resolve links in stable topological order. The graph is a fan-out forest:
    # every link has exactly one parent and therefore one root source.
    prepared_link_values = _PreparedTriggerLink[]
    unresolved = copy(links)
    while !isempty(unresolved)
        resolved_any = false
        index = 1
        while index <= length(unresolved)
            link = unresolved[index]
            parent_kind, parent_slot, root_slot = _resolve_trigger_parent(
                link.parent, sources, prepared_link_values)
            if iszero(parent_slot)
                index += 1
                continue
            end
            start, count = _append_trigger_trace!(trace_entries, fault_ids,
                fault_domain, link.faults)
            push!(prepared_link_values, _PreparedTriggerLink(link.id,
                parent_kind, UInt32(parent_slot), UInt32(root_slot),
                link.propagation_delay, link.timing_skew,
                link.timestamp_label_offset, start, count,
                UInt32(length(prepared_link_values) + 1)))
            deleteat!(unresolved, index)
            resolved_any = true
        end
        resolved_any || begin
            missing = first(unresolved)
            reason = _unresolved_trigger_parent_is_cycle(
                missing.parent, links) ? :cyclic_topology : :unknown_parent
            _trigger_topology_error(reason,
                "cannot resolve parent $(missing.parent) for trigger link " *
                "$(missing.id)")
        end
    end
    prepared_links = Memory{_PreparedTriggerLink}(undef,
        length(prepared_link_values))
    copyto!(prepared_links, prepared_link_values)

    prepared_consumers = Memory{_PreparedTriggerConsumer}(undef,
        length(consumers))
    @inbounds for (slot, consumer) in enumerate(consumers)
        parent_kind, parent_slot, root_slot = _resolve_trigger_parent(
            consumer.parent, sources, prepared_link_values)
        iszero(parent_slot) && _trigger_topology_error(:unknown_parent,
            "cannot resolve parent $(consumer.parent) for trigger consumer " *
            "$(consumer.id)")
        prepared_consumers[slot] = _PreparedTriggerConsumer(consumer.id,
            parent_kind, UInt32(parent_slot), UInt32(root_slot), UInt32(slot))
    end

    link_child_offsets, link_child_slots, consumer_child_offsets,
        consumer_child_slots = _prepare_trigger_child_registry(
            length(prepared_sources), length(prepared_links), prepared_links,
            prepared_consumers)

    prepared_trace_entries = Memory{_PreparedTriggerFaultSample}(undef,
        length(trace_entries))
    copyto!(prepared_trace_entries, trace_entries)
    prepared_fault_ids = Memory{TriggerFaultID}(undef, length(fault_ids))
    copyto!(prepared_fault_ids, fault_ids)
    provisional = PreparedTriggerTopology(_TriggerTopologyBinding(),
        prepared_sources, prepared_links, prepared_consumers,
        link_child_offsets, link_child_slots, consumer_child_offsets,
        consumer_child_slots,
        prepared_trace_entries, prepared_fault_ids, fault_domain, capacity,
        0, 0, 0, 0)
    propagation_capacity, realization_capacity, observation_capacity,
        per_source_outputs = _trigger_structural_capacities(provisional)
    provisional = PreparedTriggerTopology(provisional.binding,
        prepared_sources, prepared_links, prepared_consumers,
        link_child_offsets, link_child_slots, consumer_child_offsets,
        consumer_child_slots,
        prepared_trace_entries, prepared_fault_ids, fault_domain, capacity, 0,
        propagation_capacity, realization_capacity, observation_capacity)
    required_capacity = _validate_trigger_paths_and_capacity(provisional,
        per_source_outputs)
    capacity >= required_capacity || _trigger_topology_error(
        :insufficient_in_flight_capacity,
        "declared trigger in-flight capacity $capacity is below the " *
        "prepared conservative requirement $required_capacity")
    return PreparedTriggerTopology(provisional.binding, prepared_sources,
        prepared_links, prepared_consumers, link_child_offsets,
        link_child_slots, consumer_child_offsets, consumer_child_slots,
        prepared_trace_entries,
        prepared_fault_ids, fault_domain, capacity, required_capacity,
        propagation_capacity, realization_capacity, observation_capacity)
end

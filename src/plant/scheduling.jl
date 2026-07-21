#
# Fixed-capacity deterministic plant-event scheduling
#
# The scheduler is deliberately payload- and wall-clock-free. Model-specific
# owners retain their own immutable recurrence parameters and mutable device
# state; this layer owns only the next logical event for each prepared
# generator.
#

"""Immutable declaration of one bounded plant-event generator."""
struct EventGeneratorDefinition
    first_timestamp::PlantTimestamp
    phase::PlantEventPhase
    ordinal::UInt32
    first_occurrence::UInt64
    initially_active::Bool

    function EventGeneratorDefinition(first_timestamp::PlantTimestamp,
        phase::PlantEventPhase, ordinal::UInt32,
        first_occurrence::UInt64, initially_active::Bool,
        ::_PlantTimeToken)
        return new(first_timestamp, phase, ordinal, first_occurrence,
            initially_active)
    end
end

function EventGeneratorDefinition(first_timestamp::PlantTimestamp,
    phase::PlantEventPhase, ordinal::Integer;
    occurrence::Integer=1, active::Bool=true)
    return EventGeneratorDefinition(
        first_timestamp,
        phase,
        _checked_event_ordinal(ordinal),
        _checked_event_occurrence(occurrence),
        active,
        _PLANT_TIME_TOKEN,
    )
end

function EventGeneratorDefinition(schedule::PeriodicSchedule,
    phase::PlantEventPhase, ordinal::Integer;
    origin::PlantTimestamp=zero(PlantTimestamp),
    occurrence::Integer=1, active::Bool=true)
    checked_occurrence = _checked_event_occurrence(occurrence)
    first_timestamp = schedule_timestamp(schedule, checked_occurrence,
        origin)
    return EventGeneratorDefinition(first_timestamp, phase, ordinal;
        occurrence=checked_occurrence, active=active)
end

@inline event_generator_phase(definition::EventGeneratorDefinition) =
    definition.phase
@inline event_generator_ordinal(definition::EventGeneratorDefinition) =
    definition.ordinal
@inline first_event_timestamp(definition::EventGeneratorDefinition) =
    definition.first_timestamp
@inline first_event_occurrence(definition::EventGeneratorDefinition) =
    definition.first_occurrence
@inline initially_active(definition::EventGeneratorDefinition) =
    definition.initially_active

mutable struct _EventSchedulerBinding end

"""
Fixed-capacity, run-immutable event-scheduler plan. The flat `Memory` registry
does not grow during execution; mutable generator cursors and due-scan storage
are owned separately.
"""
struct PreparedEventScheduler
    binding::_EventSchedulerBinding
    definitions::Memory{EventGeneratorDefinition}
    capacity::Int
end

@inline event_generator_count(scheduler::PreparedEventScheduler) =
    length(getfield(scheduler, :definitions))
@inline event_scheduler_capacity(scheduler::PreparedEventScheduler) =
    scheduler.capacity
@inline _event_scheduler_binding_id(scheduler::PreparedEventScheduler) =
    UInt64(objectid(getfield(scheduler, :binding)))

@inline _require_event_generator_definition(
    definition::EventGeneratorDefinition) = definition

function _require_event_generator_definition(value)
    throw(PlantScheduleError(:event_scheduler, :invalid_definition,
        "event scheduler entries must be EventGeneratorDefinition values; " *
        "got $(typeof(value))"))
end

@inline function _checked_event_scheduler_capacity(capacity::Integer)
    capacity >= 0 || throw(PlantScheduleError(:event_scheduler,
        :invalid_capacity, "event scheduler capacity must be nonnegative"))
    capacity <= min(typemax(Int), typemax(UInt32)) || throw(
        PlantScheduleError(:event_scheduler, :invalid_capacity,
            "event scheduler capacity exceeds the supported UInt32 slot range"))
    return Int(capacity)
end

@inline _checked_event_scheduler_capacity(::Bool) = throw(
    PlantScheduleError(:event_scheduler, :invalid_capacity,
        "event scheduler capacity must be an integer count, not Bool"))

function _copy_event_generator_definitions(definitions, count::Int)
    owned = Memory{EventGeneratorDefinition}(undef, count)
    index = 0
    for value in definitions
        index += 1
        owned[index] = _require_event_generator_definition(value)
    end
    return owned
end

function _require_unique_event_ordinals(
    definitions::Memory{EventGeneratorDefinition})
    @inbounds for right_index in 2:length(definitions)
        right = definitions[right_index]
        for left_index in 1:(right_index - 1)
            left = definitions[left_index]
            if left.phase == right.phase && left.ordinal == right.ordinal
                throw(PlantScheduleError(:event_scheduler,
                    :duplicate_ordinal,
                    "event generator ordinals must be unique within each " *
                    "plant event phase; duplicate ordinal $(right.ordinal) " *
                    "in $(right.phase)"))
            end
        end
    end
    return nothing
end

@inline function _event_generator_precedes(left::EventGeneratorDefinition,
    right::EventGeneratorDefinition)
    left.phase == right.phase ||
        return UInt8(left.phase) < UInt8(right.phase)
    return left.ordinal < right.ordinal
end

function _canonicalize_event_generator_definitions!(
    definitions::Memory{EventGeneratorDefinition})
    @inbounds for index in 2:length(definitions)
        definition = definitions[index]
        insertion = index - 1
        while insertion >= 1 && _event_generator_precedes(
                definition, definitions[insertion])
            definitions[insertion + 1] = definitions[insertion]
            insertion -= 1
        end
        definitions[insertion + 1] = definition
    end
    return definitions
end

function _prepare_event_scheduler(definitions, count::Int,
    capacity::Integer)
    checked_capacity = _checked_event_scheduler_capacity(capacity)
    count <= checked_capacity || throw(PlantScheduleError(
        :event_scheduler, :capacity_overflow,
        "event scheduler received $count generators for declared capacity " *
        "$checked_capacity"))
    owned = _copy_event_generator_definitions(definitions, count)
    _require_unique_event_ordinals(owned)
    _canonicalize_event_generator_definitions!(owned)
    return PreparedEventScheduler(_EventSchedulerBinding(), owned,
        checked_capacity)
end

"""
    prepare_event_scheduler(definitions; capacity=length(definitions))

Prepare a flat, fixed-capacity registry of event generators. Preparation copies
the immutable declarations, rejects capacity overflow, and requires ordinals to
be unique within each causal event phase. It does not read a wall clock,
materialize run-length event lists, or attach model callbacks.
"""
function prepare_event_scheduler(definitions::Tuple;
    capacity::Integer=length(definitions))
    return _prepare_event_scheduler(definitions, length(definitions),
        capacity)
end

function prepare_event_scheduler(definitions::AbstractVector;
    capacity::Integer=length(definitions))
    return _prepare_event_scheduler(definitions, length(definitions),
        capacity)
end

function prepare_event_scheduler(definitions; capacity::Integer=0)
    throw(PlantScheduleError(:event_scheduler, :invalid_registry,
        "event scheduler definitions must be a finite Tuple or " *
        "AbstractVector; got $(typeof(definitions))"))
end

struct EventGeneratorHandle
    binding_id::UInt64
    slot::UInt32
end

function event_generator_handle(scheduler::PreparedEventScheduler,
    phase::PlantEventPhase, ordinal::Integer)
    checked_ordinal = _checked_event_ordinal(ordinal)
    definitions = getfield(scheduler, :definitions)
    @inbounds for index in eachindex(definitions)
        definition = definitions[index]
        if definition.phase == phase && definition.ordinal == checked_ordinal
            return EventGeneratorHandle(_event_scheduler_binding_id(scheduler),
                UInt32(index))
        end
    end
    throw(PlantScheduleError(:event_scheduler, :unknown_generator,
        "event scheduler has no generator in $phase with ordinal " *
        "$checked_ordinal"))
end

@enum _EventGeneratorStatus::UInt8 begin
    _InactiveEventGenerator = 0x00
    _ScheduledEventGenerator = 0x01
    _ClaimedEventGenerator = 0x02
end

"""Compact mutable-by-replacement cursor stored in scheduler-owned memory."""
struct EventGeneratorCursor
    next_timestamp::PlantTimestamp
    next_occurrence::UInt64
    status::_EventGeneratorStatus
end

"""Single-writer mutable state for one prepared event scheduler."""
mutable struct EventSchedulerState
    binding_id::UInt64
    cursors::Memory{EventGeneratorCursor}
    current_timestamp::PlantTimestamp
    last_key::PlantEventKey
    revision::UInt64
    claim_sequence::UInt64
    outstanding_slot::UInt32
    has_last_key::Bool
end

function EventSchedulerState(scheduler::PreparedEventScheduler;
    initial_timestamp::PlantTimestamp=zero(PlantTimestamp))
    definitions = getfield(scheduler, :definitions)
    cursors = Memory{EventGeneratorCursor}(undef, length(definitions))
    @inbounds for index in eachindex(definitions)
        definition = definitions[index]
        if definition.initially_active &&
                definition.first_timestamp < initial_timestamp
            throw(PlantScheduleError(:event_scheduler, :time_regression,
                "active event generator $(definition.phase)/" *
                "$(definition.ordinal) begins before the scheduler's " *
                "initial plant timestamp"))
        end
        status = definition.initially_active ?
            _ScheduledEventGenerator : _InactiveEventGenerator
        cursors[index] = EventGeneratorCursor(definition.first_timestamp,
            definition.first_occurrence, status)
    end
    sentinel = PlantEventKey(initial_timestamp, TriggerUpdatePhase,
        UInt32(1), UInt64(1), _PLANT_TIME_TOKEN)
    return EventSchedulerState(_event_scheduler_binding_id(scheduler), cursors,
        initial_timestamp, sentinel, UInt64(0), UInt64(0), UInt32(0), false)
end

"""Caller-owned fixed-capacity scratch for one serial due-event scan."""
mutable struct EventSchedulerWorkspace
    binding_id::UInt64
    due_slots::Memory{UInt32}
    due_count::Int
    due_timestamp::PlantTimestamp
    state_revision::UInt64
    valid::Bool
end

function EventSchedulerWorkspace(scheduler::PreparedEventScheduler)
    return EventSchedulerWorkspace(_event_scheduler_binding_id(scheduler),
        Memory{UInt32}(undef, scheduler.capacity), 0,
        zero(PlantTimestamp), UInt64(0), false)
end

@inline scheduler_timestamp(state::EventSchedulerState) =
    state.current_timestamp

@inline function _require_scheduler_binding(scheduler::PreparedEventScheduler,
    state::EventSchedulerState)
    _event_scheduler_binding_id(scheduler) == state.binding_id || throw(
        PlantScheduleError(:event_scheduler, :foreign_state,
            "event scheduler state belongs to another prepared scheduler"))
    return nothing
end

@inline function _require_scheduler_binding(scheduler::PreparedEventScheduler,
    workspace::EventSchedulerWorkspace)
    _event_scheduler_binding_id(scheduler) == workspace.binding_id || throw(
        PlantScheduleError(:event_scheduler, :foreign_workspace,
            "event scheduler workspace belongs to another prepared scheduler"))
    return nothing
end

@inline function _require_scheduler_binding(scheduler::PreparedEventScheduler,
    handle::EventGeneratorHandle)
    _event_scheduler_binding_id(scheduler) == handle.binding_id || throw(
        PlantScheduleError(:event_scheduler, :foreign_generator,
            "event generator handle belongs to another prepared scheduler"))
    return nothing
end

@inline function _require_event_generator_slot(scheduler::PreparedEventScheduler,
    handle::EventGeneratorHandle)
    _require_scheduler_binding(scheduler, handle)
    slot = Int(handle.slot)
    1 <= slot <= event_generator_count(scheduler) || throw(
        PlantScheduleError(:event_scheduler, :invalid_generator,
            "event generator handle contains an invalid registry slot"))
    return slot
end

@inline function _require_idle_scheduler(state::EventSchedulerState)
    iszero(state.outstanding_slot) || throw(PlantScheduleError(
        :event_scheduler, :outstanding_claim,
        "resolve the outstanding event claim before another scheduler operation"))
    return nothing
end

@inline function _checked_scheduler_counter(value::UInt64, reason::Symbol,
    label::AbstractString)
    value != typemax(UInt64) || throw(PlantScheduleError(
        :event_scheduler, reason, "$label exceeds UInt64 range"))
    return value + UInt64(1)
end

@inline function _event_key(definition::EventGeneratorDefinition,
    cursor::EventGeneratorCursor)
    return PlantEventKey(cursor.next_timestamp, definition.phase,
        definition.ordinal, cursor.next_occurrence, _PLANT_TIME_TOKEN)
end

@inline function _require_forward_event_key(state::EventSchedulerState,
    key::PlantEventKey)
    if state.has_last_key && !isless(state.last_key, key)
        reason = key.timestamp < state.current_timestamp ?
            :time_regression : :event_order_regression
        throw(PlantScheduleError(:event_scheduler, reason,
            "event generator cannot schedule a key at or before the last " *
            "claimed plant event"))
    end
    key.timestamp < state.current_timestamp && throw(PlantScheduleError(
        :event_scheduler, :time_regression,
        "event generator cannot schedule before the current plant timestamp"))
    return nothing
end

"""
    scan_due_events!(workspace, scheduler, state)

Populate caller-owned due-slot storage for the minimum scheduled plant
timestamp and return the number of simultaneous candidates. Slots are ordered
by the complete logical event key. The operation neither claims an event nor
advances scheduler state.
"""
function scan_due_events!(workspace::EventSchedulerWorkspace,
    scheduler::PreparedEventScheduler, state::EventSchedulerState)
    _require_scheduler_binding(scheduler, state)
    _require_scheduler_binding(scheduler, workspace)
    _require_idle_scheduler(state)

    definitions = getfield(scheduler, :definitions)
    cursors = state.cursors
    found = false
    minimum_timestamp = state.current_timestamp
    count = 0
    @inbounds for index in eachindex(cursors)
        cursor = cursors[index]
        cursor.status == _ScheduledEventGenerator || continue
        timestamp = cursor.next_timestamp
        if !found || timestamp < minimum_timestamp
            found = true
            minimum_timestamp = timestamp
            count = 1
            workspace.due_slots[1] = UInt32(index)
        elseif timestamp == minimum_timestamp
            count += 1
            workspace.due_slots[count] = UInt32(index)
        end
    end

    if !found
        workspace.due_count = 0
        workspace.due_timestamp = state.current_timestamp
        workspace.state_revision = state.revision
        workspace.valid = true
        return 0
    end

    minimum_timestamp < state.current_timestamp && throw(
        PlantScheduleError(:event_scheduler, :time_regression,
            "next due event precedes the current plant timestamp"))
    # Preparation stores definitions in phase/ordinal order. Each generator
    # owns only one cursor, so the filtered due slots retain total key order
    # without a runtime sort.
    first_slot = Int(workspace.due_slots[1])
    first_key = _event_key(definitions[first_slot], cursors[first_slot])
    _require_forward_event_key(state, first_key)
    workspace.due_count = count
    workspace.due_timestamp = minimum_timestamp
    workspace.state_revision = state.revision
    workspace.valid = true
    return count
end

@inline function _require_fresh_due_scan(workspace::EventSchedulerWorkspace,
    scheduler::PreparedEventScheduler, state::EventSchedulerState)
    _require_scheduler_binding(scheduler, state)
    _require_scheduler_binding(scheduler, workspace)
    workspace.valid && workspace.state_revision == state.revision || throw(
        PlantScheduleError(:event_scheduler, :stale_due_scan,
            "due-event workspace must be rescanned after scheduler state changes"))
    return nothing
end

function due_event_count(workspace::EventSchedulerWorkspace,
    scheduler::PreparedEventScheduler, state::EventSchedulerState)
    _require_fresh_due_scan(workspace, scheduler, state)
    return workspace.due_count
end

function due_event_timestamp(workspace::EventSchedulerWorkspace,
    scheduler::PreparedEventScheduler, state::EventSchedulerState)
    _require_fresh_due_scan(workspace, scheduler, state)
    workspace.due_count > 0 || return nothing
    return workspace.due_timestamp
end

function due_event_key(workspace::EventSchedulerWorkspace,
    scheduler::PreparedEventScheduler, state::EventSchedulerState,
    index::Integer)
    _require_fresh_due_scan(workspace, scheduler, state)
    1 <= index <= workspace.due_count || throw(PlantScheduleError(
        :event_scheduler, :invalid_due_index,
        "due-event index must be within the current due scan"))
    slot = Int(workspace.due_slots[Int(index)])
    return _event_key(getfield(scheduler, :definitions)[slot],
        state.cursors[slot])
end

@inline due_event_key(::EventSchedulerWorkspace,
    ::PreparedEventScheduler, ::EventSchedulerState, ::Bool) = throw(
        PlantScheduleError(:event_scheduler, :invalid_due_index,
            "due-event index must be an integer count, not Bool"))

"""One outstanding serial claim for a prepared event generator."""
struct EventClaim
    binding_id::UInt64
    slot::UInt32
    key::PlantEventKey
    sequence::UInt64
end

@inline claimed_event_key(claim::EventClaim) = claim.key
@inline event_generator_handle(claim::EventClaim) =
    EventGeneratorHandle(claim.binding_id, claim.slot)

"""
    claim_next_event!(workspace, scheduler, state)

Claim the globally earliest logical event, or return `nothing` when every
prepared generator is inactive. Exactly one claim may be outstanding for the
single-writer state.
"""
function claim_next_event!(workspace::EventSchedulerWorkspace,
    scheduler::PreparedEventScheduler, state::EventSchedulerState)
    count = scan_due_events!(workspace, scheduler, state)
    iszero(count) && return nothing

    slot = Int(workspace.due_slots[1])
    definitions = getfield(scheduler, :definitions)
    cursor = state.cursors[slot]
    key = _event_key(definitions[slot], cursor)
    _require_forward_event_key(state, key)
    next_revision = _checked_scheduler_counter(state.revision,
        :revision_overflow, "event scheduler state revision")
    next_claim_sequence = _checked_scheduler_counter(state.claim_sequence,
        :claim_sequence_overflow, "event scheduler claim sequence")

    state.cursors[slot] = EventGeneratorCursor(cursor.next_timestamp,
        cursor.next_occurrence, _ClaimedEventGenerator)
    state.current_timestamp = key.timestamp
    state.last_key = key
    state.revision = next_revision
    state.claim_sequence = next_claim_sequence
    state.outstanding_slot = UInt32(slot)
    state.has_last_key = true
    workspace.valid = false
    return EventClaim(_event_scheduler_binding_id(scheduler), UInt32(slot), key,
        next_claim_sequence)
end

@inline function _require_current_claim(scheduler::PreparedEventScheduler,
    state::EventSchedulerState, claim::EventClaim)
    _require_scheduler_binding(scheduler, state)
    _event_scheduler_binding_id(scheduler) == claim.binding_id || throw(
        PlantScheduleError(:event_scheduler, :foreign_claim,
            "event claim belongs to another prepared scheduler"))
    state.outstanding_slot == claim.slot &&
        state.claim_sequence == claim.sequence || throw(PlantScheduleError(
            :event_scheduler, :stale_claim,
            "event claim is stale or is not the outstanding claim"))
    slot = Int(claim.slot)
    1 <= slot <= length(state.cursors) || throw(PlantScheduleError(
        :event_scheduler, :invalid_claim,
        "event claim contains an invalid registry slot"))
    cursor = state.cursors[slot]
    definition = getfield(scheduler, :definitions)[slot]
    cursor.status == _ClaimedEventGenerator &&
        _event_key(definition, cursor) == claim.key || throw(
            PlantScheduleError(:event_scheduler, :stale_claim,
                "event claim no longer matches its generator cursor"))
    return slot, cursor, definition
end

@inline function _next_event_occurrence(occurrence::UInt64)
    occurrence != typemax(UInt64) || throw(PlantScheduleError(
        :event_scheduler, :occurrence_overflow,
        "event generator occurrence exceeds UInt64 range"))
    return occurrence + UInt64(1)
end

"""
    reschedule_event!(scheduler, state, claim, next_timestamp)

Resolve the outstanding claim by publishing one next deadline for the same
generator. The next key must be strictly later in total logical order; an
equal timestamp is valid only because the per-generator occurrence advances.
"""
function reschedule_event!(scheduler::PreparedEventScheduler,
    state::EventSchedulerState, claim::EventClaim,
    next_timestamp::PlantTimestamp)
    slot, cursor, definition = _require_current_claim(scheduler, state,
        claim)
    next_occurrence = _next_event_occurrence(cursor.next_occurrence)
    next_key = PlantEventKey(next_timestamp, definition.phase,
        definition.ordinal, next_occurrence, _PLANT_TIME_TOKEN)
    _require_forward_event_key(state, next_key)
    next_revision = _checked_scheduler_counter(state.revision,
        :revision_overflow, "event scheduler state revision")

    state.cursors[slot] = EventGeneratorCursor(next_timestamp,
        next_occurrence, _ScheduledEventGenerator)
    state.revision = next_revision
    state.outstanding_slot = UInt32(0)
    return next_key
end

function reschedule_periodic_event!(scheduler::PreparedEventScheduler,
    state::EventSchedulerState, claim::EventClaim,
    schedule::PeriodicSchedule;
    origin::PlantTimestamp=zero(PlantTimestamp))
    next_occurrence = _next_event_occurrence(claim.key.occurrence)
    timestamp = schedule_timestamp(schedule, next_occurrence, origin)
    return reschedule_event!(scheduler, state, claim, timestamp)
end

"""Activate a prepared inactive generator at an explicit plant timestamp."""
function activate_event_generator!(scheduler::PreparedEventScheduler,
    state::EventSchedulerState, handle::EventGeneratorHandle,
    timestamp::PlantTimestamp)
    _require_scheduler_binding(scheduler, state)
    _require_idle_scheduler(state)
    slot = _require_event_generator_slot(scheduler, handle)
    cursor = state.cursors[slot]
    cursor.status == _InactiveEventGenerator || throw(PlantScheduleError(
        :event_scheduler, :invalid_transition,
        "only an inactive event generator can be activated"))
    definition = getfield(scheduler, :definitions)[slot]
    key = PlantEventKey(timestamp, definition.phase, definition.ordinal,
        cursor.next_occurrence, _PLANT_TIME_TOKEN)
    _require_forward_event_key(state, key)
    next_revision = _checked_scheduler_counter(state.revision,
        :revision_overflow, "event scheduler state revision")
    state.cursors[slot] = EventGeneratorCursor(timestamp,
        cursor.next_occurrence, _ScheduledEventGenerator)
    state.revision = next_revision
    return key
end

"""Cancel one scheduled, unclaimed event without consuming its occurrence."""
function deactivate_event_generator!(scheduler::PreparedEventScheduler,
    state::EventSchedulerState, handle::EventGeneratorHandle)
    _require_scheduler_binding(scheduler, state)
    _require_idle_scheduler(state)
    slot = _require_event_generator_slot(scheduler, handle)
    cursor = state.cursors[slot]
    cursor.status == _ScheduledEventGenerator || throw(PlantScheduleError(
        :event_scheduler, :invalid_transition,
        "only a scheduled event generator can be deactivated by handle"))
    next_revision = _checked_scheduler_counter(state.revision,
        :revision_overflow, "event scheduler state revision")
    state.cursors[slot] = EventGeneratorCursor(cursor.next_timestamp,
        cursor.next_occurrence, _InactiveEventGenerator)
    state.revision = next_revision
    return handle
end

"""Resolve a claimed occurrence and leave its generator inactive."""
function deactivate_event_generator!(scheduler::PreparedEventScheduler,
    state::EventSchedulerState, claim::EventClaim)
    slot, cursor, _ = _require_current_claim(scheduler, state, claim)
    next_occurrence = _next_event_occurrence(cursor.next_occurrence)
    next_revision = _checked_scheduler_counter(state.revision,
        :revision_overflow, "event scheduler state revision")
    state.cursors[slot] = EventGeneratorCursor(cursor.next_timestamp,
        next_occurrence, _InactiveEventGenerator)
    state.revision = next_revision
    state.outstanding_slot = UInt32(0)
    return EventGeneratorHandle(_event_scheduler_binding_id(scheduler),
        UInt32(slot))
end

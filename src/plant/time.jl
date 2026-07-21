const _PLANT_NANOSECONDS_PER_SECOND = Int64(1_000_000_000)

struct _PlantTimeToken end
const _PLANT_TIME_TOKEN = _PlantTimeToken()

@inline function _checked_plant_nanoseconds(value::Integer,
    component::Symbol)
    value >= 0 || throw(PlantTimeError(component, :negative,
        "$component nanoseconds must be nonnegative"))
    value <= typemax(Int64) || throw(PlantTimeError(component, :overflow,
        "$component nanoseconds exceed Int64 range"))
    return Int64(value)
end

@inline _checked_plant_nanoseconds(::Bool, component::Symbol) =
    throw(PlantTimeError(component, :invalid_type,
        "$component nanoseconds must be an integer count, not Bool"))

"""
    PlantTimestamp(nanoseconds)

One nonnegative instant on the run-local canonical plant timeline. The value is
stored as checked integer nanoseconds from the run origin and is deliberately a
different type from [`PlantDuration`](@ref).
"""
struct PlantTimestamp
    nanoseconds::Int64

    PlantTimestamp(nanoseconds::Int64, ::_PlantTimeToken) = new(nanoseconds)
end

PlantTimestamp(nanoseconds::Integer) = PlantTimestamp(
    _checked_plant_nanoseconds(nanoseconds, :plant_timestamp),
    _PLANT_TIME_TOKEN)

"""
    PlantDuration(nanoseconds)

One nonnegative elapsed interval in integer nanoseconds. A duration is not an
absolute timestamp, a sample period, or a floating-point integration time.
"""
struct PlantDuration
    nanoseconds::Int64

    PlantDuration(nanoseconds::Int64, ::_PlantTimeToken) = new(nanoseconds)
end

PlantDuration(nanoseconds::Integer) = PlantDuration(
    _checked_plant_nanoseconds(nanoseconds, :plant_duration),
    _PLANT_TIME_TOKEN)

@inline function _checked_plant_offset_nanoseconds(value::Integer)
    typemin(Int64) <= value <= typemax(Int64) || throw(PlantTimeError(
        :plant_time_offset, :overflow,
        "plant time offset nanoseconds exceed Int64 range"))
    return Int64(value)
end

@inline _checked_plant_offset_nanoseconds(::Bool) = throw(PlantTimeError(
    :plant_time_offset, :invalid_type,
    "plant time offset nanoseconds must be an integer count, not Bool"))

"""
    PlantTimeOffset(nanoseconds)

One signed displacement between a nominal and realized plant instant. Unlike
[`PlantDuration`](@ref), an offset may be negative; applying it to a timestamp
still must produce a representable nonnegative plant instant.
"""
struct PlantTimeOffset
    nanoseconds::Int64

    PlantTimeOffset(nanoseconds::Int64, ::_PlantTimeToken) = new(nanoseconds)
end

PlantTimeOffset(nanoseconds::Integer) = PlantTimeOffset(
    _checked_plant_offset_nanoseconds(nanoseconds), _PLANT_TIME_TOKEN)

@inline plant_nanoseconds(value::Union{
    PlantTimestamp,PlantDuration,PlantTimeOffset}) =
    value.nanoseconds

Base.zero(::Type{PlantTimestamp}) = PlantTimestamp(0, _PLANT_TIME_TOKEN)
Base.zero(::PlantTimestamp) = zero(PlantTimestamp)
Base.zero(::Type{PlantDuration}) = PlantDuration(0, _PLANT_TIME_TOKEN)
Base.zero(::PlantDuration) = zero(PlantDuration)
Base.zero(::Type{PlantTimeOffset}) = PlantTimeOffset(0, _PLANT_TIME_TOKEN)
Base.zero(::PlantTimeOffset) = zero(PlantTimeOffset)
Base.iszero(value::Union{PlantTimestamp,PlantDuration,PlantTimeOffset}) =
    iszero(value.nanoseconds)

Base.:(==)(left::PlantTimestamp, right::PlantTimestamp) =
    left.nanoseconds == right.nanoseconds
Base.:(==)(left::PlantDuration, right::PlantDuration) =
    left.nanoseconds == right.nanoseconds
Base.:(==)(left::PlantTimeOffset, right::PlantTimeOffset) =
    left.nanoseconds == right.nanoseconds
Base.isless(left::PlantTimestamp, right::PlantTimestamp) =
    isless(left.nanoseconds, right.nanoseconds)
Base.isless(left::PlantDuration, right::PlantDuration) =
    isless(left.nanoseconds, right.nanoseconds)
Base.:(<)(left::PlantTimestamp, right::PlantTimestamp) = isless(left, right)
Base.:(<)(left::PlantDuration, right::PlantDuration) = isless(left, right)
Base.:(<=)(left::PlantTimestamp, right::PlantTimestamp) =
    left.nanoseconds <= right.nanoseconds
Base.:(<=)(left::PlantDuration, right::PlantDuration) =
    left.nanoseconds <= right.nanoseconds
Base.hash(value::PlantTimestamp, seed::UInt) =
    hash(value.nanoseconds, hash(:PlantTimestamp, seed))
Base.hash(value::PlantDuration, seed::UInt) =
    hash(value.nanoseconds, hash(:PlantDuration, seed))
Base.hash(value::PlantTimeOffset, seed::UInt) =
    hash(value.nanoseconds, hash(:PlantTimeOffset, seed))

function Base.show(io::IO, value::PlantTimestamp)
    print(io, "PlantTimestamp(", value.nanoseconds, " ns)")
end

function Base.show(io::IO, value::PlantDuration)
    print(io, "PlantDuration(", value.nanoseconds, " ns)")
end


function Base.show(io::IO, value::PlantTimeOffset)
    print(io, "PlantTimeOffset(", value.nanoseconds, " ns)")
end

@inline function _checked_plant_add(left::Int64, right::Int64,
    operation::Symbol)
    right <= typemax(Int64) - left || throw(PlantTimeError(
        operation, :overflow, "plant-time addition exceeds Int64 range"))
    return left + right
end

@inline function _checked_plant_subtract(left::Int64, right::Int64,
    operation::Symbol)
    left >= right || throw(PlantTimeError(operation, :negative_result,
        "plant-time subtraction would produce a negative value"))
    return left - right
end


@inline function _checked_plant_offset_add(left::Int64, right::Int64,
    operation::Symbol)
    if right > 0
        left <= typemax(Int64) - right || throw(PlantTimeError(
            operation, :overflow,
            "plant-time offset addition exceeds Int64 range"))
    elseif right < 0
        left >= typemin(Int64) - right || throw(PlantTimeError(
            operation, :overflow,
            "plant-time offset addition exceeds Int64 range"))
    end
    return left + right
end

@inline function _checked_plant_scale(value::Int64, factor::Integer,
    operation::Symbol)
    factor >= 0 || throw(PlantTimeError(operation, :negative_factor,
        "plant-time scale factor must be nonnegative"))
    if iszero(value) || iszero(factor)
        return Int64(0)
    end
    factor <= typemax(Int64) ÷ value || throw(PlantTimeError(
        operation, :overflow, "plant-time multiplication exceeds Int64 range"))
    return value * Int64(factor)
end

@inline _checked_plant_scale(::Int64, ::Bool, operation::Symbol) =
    throw(PlantTimeError(operation, :invalid_type,
        "plant-time scale factor must be an integer count, not Bool"))

@inline function Base.:+(timestamp::PlantTimestamp,
    duration::PlantDuration)
    return PlantTimestamp(_checked_plant_add(timestamp.nanoseconds,
        duration.nanoseconds, :timestamp_add), _PLANT_TIME_TOKEN)
end

Base.:+(duration::PlantDuration, timestamp::PlantTimestamp) =
    timestamp + duration

@inline function Base.:-(timestamp::PlantTimestamp,
    duration::PlantDuration)
    return PlantTimestamp(_checked_plant_subtract(timestamp.nanoseconds,
        duration.nanoseconds, :timestamp_subtract), _PLANT_TIME_TOKEN)
end

@inline function Base.:-(later::PlantTimestamp, earlier::PlantTimestamp)
    return PlantDuration(_checked_plant_subtract(later.nanoseconds,
        earlier.nanoseconds, :timestamp_difference), _PLANT_TIME_TOKEN)
end

@inline function Base.:+(left::PlantDuration, right::PlantDuration)
    return PlantDuration(_checked_plant_add(left.nanoseconds,
        right.nanoseconds, :duration_add), _PLANT_TIME_TOKEN)
end

@inline function Base.:-(left::PlantDuration, right::PlantDuration)
    return PlantDuration(_checked_plant_subtract(left.nanoseconds,
        right.nanoseconds, :duration_subtract), _PLANT_TIME_TOKEN)
end

@inline function Base.:*(duration::PlantDuration, factor::Integer)
    return PlantDuration(_checked_plant_scale(duration.nanoseconds, factor,
        :duration_scale), _PLANT_TIME_TOKEN)
end

Base.:*(factor::Integer, duration::PlantDuration) = duration * factor


@inline function Base.:+(left::PlantTimeOffset, right::PlantTimeOffset)
    return PlantTimeOffset(_checked_plant_offset_add(left.nanoseconds,
        right.nanoseconds, :offset_add), _PLANT_TIME_TOKEN)
end


@inline function Base.:-(left::PlantTimeOffset, right::PlantTimeOffset)
    right.nanoseconds == typemin(Int64) && throw(PlantTimeError(
        :offset_subtract, :overflow,
        "plant-time offset subtraction exceeds Int64 range"))
    return left + PlantTimeOffset(-right.nanoseconds, _PLANT_TIME_TOKEN)
end


@inline function Base.:+(timestamp::PlantTimestamp,
    offset::PlantTimeOffset)
    nanoseconds = _checked_plant_offset_add(timestamp.nanoseconds,
        offset.nanoseconds, :timestamp_offset)
    nanoseconds >= 0 || throw(PlantTimeError(:timestamp_offset,
        :negative_result,
        "applying plant time offset would produce a negative timestamp"))
    return PlantTimestamp(nanoseconds, _PLANT_TIME_TOKEN)
end


Base.:+(offset::PlantTimeOffset, timestamp::PlantTimestamp) =
    timestamp + offset

@inline function plant_time_seconds(value::PlantTimestamp,
    ::Type{T}=Float64) where {T<:AbstractFloat}
    return T(value.nanoseconds) / T(_PLANT_NANOSECONDS_PER_SECOND)
end

@inline function plant_duration_seconds(value::PlantDuration,
    ::Type{T}=Float64) where {T<:AbstractFloat}
    return T(value.nanoseconds) / T(_PLANT_NANOSECONDS_PER_SECOND)
end

"""
    PeriodicSchedule(period; phase=PlantDuration(0))
    PeriodicSchedule(; period_ns, phase_ns=0)

Immutable nominal recurrence on plant time. The period must be positive; phase
is a nonnegative offset from a separately supplied run origin. Prepared cursor
state is introduced by the event scheduler rather than stored here.
"""
struct PeriodicSchedule
    period::PlantDuration
    phase::PlantDuration

    function PeriodicSchedule(period::PlantDuration, phase::PlantDuration,
        ::_PlantTimeToken)
        iszero(period) && throw(PlantScheduleError(:periodic_schedule,
            :zero_period, "PeriodicSchedule period must be > 0 ns"))
        return new(period, phase)
    end
end

PeriodicSchedule(period::PlantDuration;
    phase::PlantDuration=zero(PlantDuration)) =
    PeriodicSchedule(period, phase, _PLANT_TIME_TOKEN)

PeriodicSchedule(; period_ns::Integer, phase_ns::Integer=0) =
    PeriodicSchedule(PlantDuration(period_ns);
        phase=PlantDuration(phase_ns))

@inline schedule_period(schedule::PeriodicSchedule) = schedule.period
@inline schedule_phase(schedule::PeriodicSchedule) = schedule.phase

@inline function _checked_schedule_sequence(sequence::Integer)
    sequence >= 1 || throw(PlantScheduleError(:periodic_schedule,
        :invalid_sequence, "periodic schedule sequence must be >= 1"))
    return sequence
end

@inline _checked_schedule_sequence(::Bool) = throw(PlantScheduleError(
    :periodic_schedule, :invalid_sequence,
    "periodic schedule sequence must be an integer count, not Bool"))

"""
    schedule_timestamp(schedule, sequence, origin=PlantTimestamp(0))

Return the timestamp of the one-based nominal occurrence using checked
arithmetic. This is a pure cold/reference operation; mutable next-event state
belongs to a separately prepared scheduler cursor.
"""
@inline function schedule_timestamp(schedule::PeriodicSchedule,
    sequence::Integer, origin::PlantTimestamp=zero(PlantTimestamp))
    checked_sequence = _checked_schedule_sequence(sequence)
    offset = schedule.period * (checked_sequence - one(checked_sequence))
    return origin + schedule.phase + offset
end

@enum PlantEventPhase::UInt8 begin
    TriggerUpdatePhase = 0x01
    CommandApplicationPhase = 0x02
    AtmosphereEvolutionPhase = 0x03
    IntegrationBoundaryPhase = 0x04
    ExposureOpenPhase = 0x05
    OpticalSamplePhase = 0x06
    ReadoutCompletionPhase = 0x07
    AcquisitionReadyPhase = 0x08
end

@inline function _checked_event_ordinal(value::Integer)
    1 <= value <= typemax(UInt32) || throw(PlantScheduleError(
        :plant_event_key, :invalid_ordinal,
        "plant event ordinal must be in 1:typemax(UInt32)"))
    return UInt32(value)
end

@inline _checked_event_ordinal(::Bool) = throw(PlantScheduleError(
    :plant_event_key, :invalid_ordinal,
    "plant event ordinal must be an integer count, not Bool"))

@inline function _checked_event_occurrence(value::Integer)
    1 <= value <= typemax(UInt64) || throw(PlantScheduleError(
        :plant_event_key, :invalid_occurrence,
        "plant event occurrence must be in 1:typemax(UInt64)"))
    return UInt64(value)
end

@inline _checked_event_occurrence(::Bool) = throw(PlantScheduleError(
    :plant_event_key, :invalid_occurrence,
    "plant event occurrence must be an integer count, not Bool"))

"""
Deterministic logical ordering key for one plant event. This developer-facing
value is ordered by timestamp, causal phase, prepared ordinal, and occurrence;
runtime container and completion order are intentionally absent.
"""
struct PlantEventKey
    timestamp::PlantTimestamp
    phase::PlantEventPhase
    ordinal::UInt32
    occurrence::UInt64

    PlantEventKey(timestamp::PlantTimestamp, phase::PlantEventPhase,
        ordinal::UInt32, occurrence::UInt64, ::_PlantTimeToken) =
        new(timestamp, phase, ordinal, occurrence)
end

PlantEventKey(timestamp::PlantTimestamp, phase::PlantEventPhase,
    ordinal::Integer, occurrence::Integer) = PlantEventKey(
    timestamp,
    phase,
    _checked_event_ordinal(ordinal),
    _checked_event_occurrence(occurrence),
    _PLANT_TIME_TOKEN,
)

Base.:(==)(left::PlantEventKey, right::PlantEventKey) =
    left.timestamp == right.timestamp && left.phase == right.phase &&
    left.ordinal == right.ordinal && left.occurrence == right.occurrence

@inline function Base.isless(left::PlantEventKey, right::PlantEventKey)
    left.timestamp == right.timestamp ||
        return isless(left.timestamp, right.timestamp)
    left.phase == right.phase || return UInt8(left.phase) < UInt8(right.phase)
    left.ordinal == right.ordinal || return left.ordinal < right.ordinal
    return left.occurrence < right.occurrence
end

Base.:(<)(left::PlantEventKey, right::PlantEventKey) = isless(left, right)

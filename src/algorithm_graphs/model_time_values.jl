const _MODEL_NANOSECONDS_PER_SECOND = Int64(1_000_000_000)

struct _ModelTimeToken end
const _MODEL_TIME_TOKEN = _ModelTimeToken()

@inline function _checked_model_nanoseconds(value::Integer, role::Symbol)
    value >= 0 || throw(AlgorithmGraphError(
        "$role nanoseconds must be nonnegative",
    ))
    value <= typemax(Int64) || throw(AlgorithmGraphError(
        "$role nanoseconds exceed Int64 range",
    ))
    return Int64(value)
end

@inline _checked_model_nanoseconds(::Bool, role::Symbol) =
    throw(AlgorithmGraphError(
        "$role nanoseconds must be an integer count, not Bool",
    ))

"""
    ModelTimestamp(nanoseconds)

One nonnegative instant on a run-local model timeline, stored as checked
integer nanoseconds from the run origin. A timestamp is deliberately distinct
from [`ModelDuration`](@ref).
"""
struct ModelTimestamp
    nanoseconds::Int64

    ModelTimestamp(nanoseconds::Int64, ::_ModelTimeToken) = new(nanoseconds)
end

ModelTimestamp(nanoseconds::Integer) = ModelTimestamp(
    _checked_model_nanoseconds(nanoseconds, :model_timestamp),
    _MODEL_TIME_TOKEN,
)

"""
    ModelDuration(nanoseconds)

One nonnegative elapsed interval in integer nanoseconds on a model timeline.
A duration is not an absolute timestamp or a floating-point exposure duration.
"""
struct ModelDuration
    nanoseconds::Int64

    ModelDuration(nanoseconds::Int64, ::_ModelTimeToken) = new(nanoseconds)
end

ModelDuration(nanoseconds::Integer) = ModelDuration(
    _checked_model_nanoseconds(nanoseconds, :model_duration),
    _MODEL_TIME_TOKEN,
)

@inline model_nanoseconds(value::Union{ModelTimestamp,ModelDuration}) =
    value.nanoseconds

Base.zero(::Type{ModelTimestamp}) = ModelTimestamp(0, _MODEL_TIME_TOKEN)
Base.zero(::ModelTimestamp) = zero(ModelTimestamp)
Base.zero(::Type{ModelDuration}) = ModelDuration(0, _MODEL_TIME_TOKEN)
Base.zero(::ModelDuration) = zero(ModelDuration)
Base.iszero(value::Union{ModelTimestamp,ModelDuration}) =
    iszero(value.nanoseconds)

Base.:(==)(left::ModelTimestamp, right::ModelTimestamp) =
    left.nanoseconds == right.nanoseconds
Base.:(==)(left::ModelDuration, right::ModelDuration) =
    left.nanoseconds == right.nanoseconds
Base.isless(left::ModelTimestamp, right::ModelTimestamp) =
    isless(left.nanoseconds, right.nanoseconds)
Base.isless(left::ModelDuration, right::ModelDuration) =
    isless(left.nanoseconds, right.nanoseconds)
Base.:(<)(left::ModelTimestamp, right::ModelTimestamp) = isless(left, right)
Base.:(<)(left::ModelDuration, right::ModelDuration) = isless(left, right)
Base.:(<=)(left::ModelTimestamp, right::ModelTimestamp) =
    left.nanoseconds <= right.nanoseconds
Base.:(<=)(left::ModelDuration, right::ModelDuration) =
    left.nanoseconds <= right.nanoseconds
Base.hash(value::ModelTimestamp, seed::UInt) =
    hash(value.nanoseconds, hash(:ModelTimestamp, seed))
Base.hash(value::ModelDuration, seed::UInt) =
    hash(value.nanoseconds, hash(:ModelDuration, seed))

function Base.show(io::IO, value::ModelTimestamp)
    print(io, "ModelTimestamp(", value.nanoseconds, " ns)")
end

function Base.show(io::IO, value::ModelDuration)
    print(io, "ModelDuration(", value.nanoseconds, " ns)")
end

@inline function _checked_model_add(
    left::Int64,
    right::Int64,
    operation::AbstractString,
)
    right <= typemax(Int64) - left || throw(AlgorithmGraphError(
        "$operation exceeds Int64 range",
    ))
    return left + right
end

@inline function _checked_model_subtract(
    left::Int64,
    right::Int64,
    operation::AbstractString,
)
    left >= right || throw(AlgorithmGraphError(
        "$operation would produce a negative value",
    ))
    return left - right
end

@inline function _checked_model_scale(value::Int64, factor::Integer)
    factor >= 0 || throw(AlgorithmGraphError(
        "model-time scale factor must be nonnegative",
    ))
    if iszero(value) || iszero(factor)
        return Int64(0)
    end
    factor <= typemax(Int64) ÷ value || throw(AlgorithmGraphError(
        "model-time multiplication exceeds Int64 range",
    ))
    return value * Int64(factor)
end

@inline _checked_model_scale(::Int64, ::Bool) = throw(AlgorithmGraphError(
    "model-time scale factor must be an integer count, not Bool",
))

@inline function Base.:+(timestamp::ModelTimestamp, duration::ModelDuration)
    return ModelTimestamp(
        _checked_model_add(
            timestamp.nanoseconds,
            duration.nanoseconds,
            "model-time addition",
        ),
        _MODEL_TIME_TOKEN,
    )
end

Base.:+(duration::ModelDuration, timestamp::ModelTimestamp) =
    timestamp + duration

@inline function Base.:-(timestamp::ModelTimestamp, duration::ModelDuration)
    return ModelTimestamp(
        _checked_model_subtract(
            timestamp.nanoseconds,
            duration.nanoseconds,
            "model-time subtraction",
        ),
        _MODEL_TIME_TOKEN,
    )
end

@inline function Base.:-(later::ModelTimestamp, earlier::ModelTimestamp)
    return ModelDuration(
        _checked_model_subtract(
            later.nanoseconds,
            earlier.nanoseconds,
            "model timestamp difference",
        ),
        _MODEL_TIME_TOKEN,
    )
end

@inline function Base.:+(left::ModelDuration, right::ModelDuration)
    return ModelDuration(
        _checked_model_add(
            left.nanoseconds,
            right.nanoseconds,
            "model-duration addition",
        ),
        _MODEL_TIME_TOKEN,
    )
end

@inline function Base.:-(left::ModelDuration, right::ModelDuration)
    return ModelDuration(
        _checked_model_subtract(
            left.nanoseconds,
            right.nanoseconds,
            "model-duration subtraction",
        ),
        _MODEL_TIME_TOKEN,
    )
end

@inline function Base.:*(duration::ModelDuration, factor::Integer)
    return ModelDuration(
        _checked_model_scale(duration.nanoseconds, factor),
        _MODEL_TIME_TOKEN,
    )
end

Base.:*(factor::Integer, duration::ModelDuration) = duration * factor

@inline function model_time_seconds(
    value::ModelTimestamp,
    ::Type{T}=Float64,
) where {T<:AbstractFloat}
    return T(value.nanoseconds) / T(_MODEL_NANOSECONDS_PER_SECOND)
end

@inline function model_duration_seconds(
    value::ModelDuration,
    ::Type{T}=Float64,
) where {T<:AbstractFloat}
    return T(value.nanoseconds) / T(_MODEL_NANOSECONDS_PER_SECOND)
end

"""
    PeriodicSchedule(period; phase=ModelDuration(0))
    PeriodicSchedule(; period_ns, phase_ns=0)

Immutable nominal recurrence on model time. The period is positive and the
phase is a nonnegative offset from a separately supplied run origin.
"""
struct PeriodicSchedule
    period::ModelDuration
    phase::ModelDuration

    function PeriodicSchedule(
        period::ModelDuration,
        phase::ModelDuration,
        ::_ModelTimeToken,
    )
        iszero(period) && throw(AlgorithmGraphError(
            "PeriodicSchedule period must be greater than zero",
        ))
        return new(period, phase)
    end
end

PeriodicSchedule(
    period::ModelDuration;
    phase::ModelDuration=zero(ModelDuration),
) = PeriodicSchedule(period, phase, _MODEL_TIME_TOKEN)

PeriodicSchedule(; period_ns::Integer, phase_ns::Integer=0) =
    PeriodicSchedule(
        ModelDuration(period_ns);
        phase=ModelDuration(phase_ns),
    )

@inline schedule_period(schedule::PeriodicSchedule) = schedule.period
@inline schedule_phase(schedule::PeriodicSchedule) = schedule.phase

@inline function _checked_schedule_sequence(sequence::Integer)
    sequence >= 1 || throw(AlgorithmGraphError(
        "periodic schedule sequence must be at least one",
    ))
    return sequence
end

@inline _checked_schedule_sequence(::Bool) = throw(AlgorithmGraphError(
    "periodic schedule sequence must be an integer count, not Bool",
))

"""
    schedule_timestamp(schedule, sequence, origin=ModelTimestamp(0))

Return the timestamp of a one-based nominal occurrence using checked model-time
arithmetic.
"""
@inline function schedule_timestamp(
    schedule::PeriodicSchedule,
    sequence::Integer,
    origin::ModelTimestamp=zero(ModelTimestamp),
)
    checked_sequence = _checked_schedule_sequence(sequence)
    offset = schedule.period * (checked_sequence - one(checked_sequence))
    return origin + schedule.phase + offset
end

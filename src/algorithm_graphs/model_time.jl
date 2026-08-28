mutable struct _ModelTimeDriverState
    sequence::UInt64
    exhausted::Bool
end

_ModelTimeDriverState() = _ModelTimeDriverState(UInt64(0), false)

"""A deterministic fixed-step cursor on canonical plant time."""
struct FixedStepModelTimeDriver
    schedule::PeriodicSchedule
    origin::PlantTimestamp
    state::_ModelTimeDriverState
end

function FixedStepModelTimeDriver(
    schedule::PeriodicSchedule;
    origin::PlantTimestamp=zero(PlantTimestamp),
)
    return FixedStepModelTimeDriver(schedule, origin, _ModelTimeDriverState())
end

"""A deterministic cursor over one finite, sealed set of model-time boundaries."""
struct PreparedBoundaryModelTimeDriver{
    Boundaries<:FixedSizeVector{PlantTimestamp},
}
    boundaries::Boundaries
    state::_ModelTimeDriverState
end

"""One mapped run-local model timestamp with uncertainty and owned provenance."""
struct CapturedModelTimestamp{Provenance}
    timestamp::PlantTimestamp
    uncertainty::PlantDuration
    provenance::Provenance

    function CapturedModelTimestamp(
        timestamp::PlantTimestamp,
        uncertainty::PlantDuration,
        provenance::Provenance,
    ) where {Provenance}
        isbitstype(Provenance) || throw(AlgorithmGraphError(
            "captured model-time provenance must have an immutable isbits representation",
        ))
        return new{Provenance}(timestamp, uncertainty, provenance)
    end
end

@inline model_timestamp(capture::CapturedModelTimestamp) = capture.timestamp
@inline model_time_uncertainty(capture::CapturedModelTimestamp) =
    capture.uncertainty
@inline model_time_provenance(capture::CapturedModelTimestamp) =
    capture.provenance

"""A deterministic cursor over sealed captured model timestamps."""
struct PreparedCapturedModelTimeDriver{Captures<:FixedSizeVector}
    captures::Captures
    state::_ModelTimeDriverState
end

const _ModelTimeDriver = Union{
    FixedStepModelTimeDriver,
    PreparedBoundaryModelTimeDriver,
    PreparedCapturedModelTimeDriver,
}

function _seal_model_time_boundaries(boundaries)
    builder = PlantTimestamp[]
    for boundary in boundaries
        boundary isa PlantTimestamp || throw(AlgorithmGraphError(
            "model-time boundaries must be PlantTimestamp values, not " *
            "$(typeof(boundary))",
        ))
        push!(builder, boundary)
    end
    isempty(builder) && throw(AlgorithmGraphError(
        "a prepared-boundary model-time driver requires at least one boundary",
    ))
    for index in 2:length(builder)
        builder[index - 1] < builder[index] || throw(AlgorithmGraphError(
            "model-time boundaries must be strictly increasing",
        ))
    end
    return FixedSizeVectorDefault{PlantTimestamp}(builder)
end

"""
    prepare_boundary_model_time_driver(boundaries)

Copy a finite, strictly increasing sequence of canonical plant timestamps into
sealed host storage. Equal-time physical transitions must be combined by the
atomic operation that consumes the boundary; this cursor does not impose an
event-phase order.
"""
function prepare_boundary_model_time_driver(boundaries)
    prepared = _seal_model_time_boundaries(boundaries)
    return PreparedBoundaryModelTimeDriver(prepared, _ModelTimeDriverState())
end

function _seal_model_time_captures(captures)
    isempty(captures) && throw(AlgorithmGraphError(
        "a captured model-time driver requires at least one capture",
    ))
    Capture = typeof(first(captures))
    Capture <: CapturedModelTimestamp || throw(AlgorithmGraphError(
        "captured model-time entries must be CapturedModelTimestamp values, not $Capture",
    ))
    isconcretetype(Capture) || throw(AlgorithmGraphError(
        "captured model-time entries must have one concrete provenance type",
    ))
    builder = Vector{Capture}(undef, length(captures))
    prior = nothing
    index = 1
    for capture in captures
        typeof(capture) === Capture || throw(AlgorithmGraphError(
            "captured model-time entries must have one concrete provenance type",
        ))
        timestamp = model_timestamp(capture)
        prior === nothing || prior < timestamp || throw(AlgorithmGraphError(
            "captured model timestamps must be strictly increasing",
        ))
        builder[index] = capture
        prior = timestamp
        index += 1
    end
    return FixedSizeVectorDefault{Capture}(builder)
end

"""
    prepare_captured_model_time_driver(captures)

Seal a finite, strictly increasing tuple or vector of captured model timestamps.
Each entry must use the same concrete isbits provenance representation.
"""
function prepare_captured_model_time_driver(
    captures::Union{Tuple,AbstractVector},
)
    prepared = _seal_model_time_captures(captures)
    return PreparedCapturedModelTimeDriver(prepared, _ModelTimeDriverState())
end

function _checked_model_time_occurrences(occurrences::Integer)
    occurrences >= 1 || throw(AlgorithmGraphError(
        "periodic model-time occurrence count must be positive",
    ))
    occurrences <= typemax(Int) || throw(AlgorithmGraphError(
        "periodic model-time occurrence count exceeds host index range",
    ))
    return Int(occurrences)
end

_checked_model_time_occurrences(::Bool) = throw(AlgorithmGraphError(
    "periodic model-time occurrence count must be an integer, not Bool",
))

function _prepare_model_time_offsets(offsets, period::PlantDuration)
    builder = PlantDuration[]
    for offset in offsets
        offset isa PlantDuration || throw(AlgorithmGraphError(
            "periodic model-time offsets must be PlantDuration values, not " *
            "$(typeof(offset))",
        ))
        offset < period || throw(AlgorithmGraphError(
            "periodic model-time offsets must be less than the schedule period",
        ))
        push!(builder, offset)
    end
    isempty(builder) && throw(AlgorithmGraphError(
        "a periodic prepared-boundary driver requires at least one offset",
    ))
    for index in 2:length(builder)
        builder[index - 1] < builder[index] || throw(AlgorithmGraphError(
            "periodic model-time offsets must be strictly increasing",
        ))
    end
    return builder
end

"""
    prepare_boundary_model_time_driver(schedule, offsets, occurrences;
                                       origin=PlantTimestamp(0))

Expand a finite number of periodic occurrences and their strictly increasing,
half-open-cycle offsets into one sealed boundary cursor. The expansion occurs
only during preparation; repeated execution performs no calendar search.
"""
function prepare_boundary_model_time_driver(
    schedule::PeriodicSchedule,
    offsets,
    occurrences::Integer;
    origin::PlantTimestamp=zero(PlantTimestamp),
)
    count = _checked_model_time_occurrences(occurrences)
    prepared_offsets = _prepare_model_time_offsets(offsets, schedule_period(schedule))
    total = try
        Base.checked_mul(length(prepared_offsets), count)
    catch error
        error isa OverflowError || rethrow()
        throw(AlgorithmGraphError(
            "periodic model-time boundary count exceeds host index range",
        ))
    end
    boundaries = Vector{PlantTimestamp}(undef, total)
    index = 1
    for occurrence in 1:count
        base = schedule_timestamp(schedule, occurrence, origin)
        for offset in prepared_offsets
            boundaries[index] = base + offset
            index += 1
        end
    end
    return prepare_boundary_model_time_driver(boundaries)
end

@inline model_time_sequence(driver::_ModelTimeDriver) = driver.state.sequence
@inline model_time_exhausted(driver::_ModelTimeDriver) = driver.state.exhausted

@noinline function _throw_model_time_exhausted()
    throw(AlgorithmGraphError("the model-time driver is exhausted"))
end

@inline function next_model_timestamp(driver::FixedStepModelTimeDriver)
    model_time_exhausted(driver) && _throw_model_time_exhausted()
    sequence = driver.state.sequence + UInt64(1)
    return schedule_timestamp(driver.schedule, sequence, driver.origin)
end

@inline function next_model_timestamp(driver::PreparedBoundaryModelTimeDriver)
    model_time_exhausted(driver) && _throw_model_time_exhausted()
    index = Int(driver.state.sequence) + 1
    @inbounds return driver.boundaries[index]
end


"""Return the next captured record without advancing its replay cursor."""
@inline function next_model_time_capture(
    driver::PreparedCapturedModelTimeDriver,
)
    model_time_exhausted(driver) && _throw_model_time_exhausted()
    index = Int(driver.state.sequence) + 1
    @inbounds return driver.captures[index]
end

@inline next_model_timestamp(driver::PreparedCapturedModelTimeDriver) =
    model_timestamp(next_model_time_capture(driver))

@inline function _commit_model_time!(driver::FixedStepModelTimeDriver)
    state = driver.state
    state.sequence += UInt64(1)
    state.exhausted = state.sequence == typemax(UInt64)
    return driver
end

@inline function _commit_model_time!(driver::PreparedBoundaryModelTimeDriver)
    state = driver.state
    state.sequence += UInt64(1)
    state.exhausted = state.sequence == UInt64(length(driver.boundaries))
    return driver
end

@inline function _commit_model_time!(driver::PreparedCapturedModelTimeDriver)
    state = driver.state
    state.sequence += UInt64(1)
    state.exhausted = state.sequence == UInt64(length(driver.captures))
    return driver
end

"""Return and commit the next deterministic model timestamp."""
@inline function advance_model_time!(driver::_ModelTimeDriver)
    timestamp = next_model_timestamp(driver)
    _commit_model_time!(driver)
    return timestamp
end

"""Reset one model-time cursor without changing its immutable schedule."""
function reset_model_time!(driver::_ModelTimeDriver)
    driver.state.sequence = UInt64(0)
    driver.state.exhausted = false
    return driver
end

"""
    step_graph_at!(graph, driver) -> timestamp

Execute one complete-frame graph at the driver's next timestamp and commit the
time cursor only after the graph succeeds. The caller may inspect
[`next_model_timestamp`](@ref) first to prepare timestamp-dependent graph
inputs. Graph and driver sequences must remain aligned.
"""
function step_graph_at!(graph::PreparedAlgorithmGraph, driver::_ModelTimeDriver)
    graph_step_sequence(graph) == model_time_sequence(driver) || throw(
        AlgorithmGraphError(
            "algorithm graph and model-time driver sequences are not aligned",
        ),
    )
    timestamp = next_model_timestamp(driver)
    step_graph!(graph)
    _commit_model_time!(driver)
    return timestamp
end

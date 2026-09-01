mutable struct GraphHILBoundaryState
    frame_sequence::UInt64
    active_command_sequence::UInt64
    command_response_required::Bool
    failed::Bool
end

GraphHILBoundaryState() = GraphHILBoundaryState(
    UInt64(0),
    UInt64(0),
    false,
    false,
)

"""
    PreparedGraphHILBoundary

One transport-neutral, single-writer lockstep binding between a prepared
algorithm graph and host command/frame exchange buffers. Persistent sequence
and failure state is owned separately from the graph, buffers, and reset
command snapshot.
"""
struct PreparedGraphHILBoundary{
    Graph<:PreparedAlgorithmGraph,
    CommandInput<:AbstractArray,
    FrameOutput<:AbstractArray,
    CommandBuffer<:Array,
    FrameBuffer<:Array,
    InitialCommand<:Array,
}
    graph::Graph
    command_input::CommandInput
    frame_output::FrameOutput
    command_buffer::CommandBuffer
    frame_buffer::FrameBuffer
    initial_command::InitialCommand
    state::GraphHILBoundaryState
end

function _prepare_hil_host_buffer(
    values::AbstractArray{T,N},
    ::Nothing,
    role::AbstractString,
) where {T,N}
    return Array{T,N}(undef, size(values))
end

function _prepare_hil_host_buffer(
    values::AbstractArray{T,N},
    buffer,
    role::AbstractString,
) where {T,N}
    buffer isa Array{T,N} || throw(AlgorithmGraphError(
        "$role buffer must be an Array{$T,$N}, not $(typeof(buffer))",
    ))
    size(buffer) == size(values) || throw(AlgorithmGraphError(
        "$role buffer shape $(size(buffer)) does not match $(size(values))",
    ))
    return buffer
end

function _copy_hil_buffer!(graph::PreparedAlgorithmGraph, destination, source)
    _with_prepared_device_execution_context(graph.context) do
        try
            copyto!(destination, source)
        finally
            _synchronize_prepared_graph_execution_context!(
                graph.execution,
                graph.context,
            )
        end
    end
    return destination
end

"""
    prepare_graph_hil_boundary(graph; command_input, frame_output,
                               command_buffer=nothing, frame_buffer=nothing)

Prepare a lockstep external-RTC boundary around an unstepped graph. The named
command input and frame output remain the graph's exact-target arrays. The
exchange buffers are ordinary host `Array`s, allocated during preparation
unless supplied by the caller. Preparation snapshots the initial graph command
for deterministic reset. After preparation, the boundary must be the only
owner that steps the graph or mutates its exact command input; transport code
writes only [`hil_command_buffer`](@ref).
"""
function prepare_graph_hil_boundary(
    graph::PreparedAlgorithmGraph;
    command_input::Symbol,
    frame_output::Symbol,
    command_buffer=nothing,
    frame_buffer=nothing,
)
    graph_failed(graph) && throw(AlgorithmGraphError(
        "a failed graph cannot be bound to a HIL boundary",
    ))
    graph_step_pending(graph) && throw(AlgorithmGraphError(
        "a graph with a pending frame cannot be bound to a HIL boundary",
    ))
    iszero(graph_step_sequence(graph)) || throw(AlgorithmGraphError(
        "a HIL boundary requires an unstepped graph",
    ))
    graph_command = graph_input(graph, command_input)
    graph_frame = graph_output(graph, frame_output)
    eltype(graph_command) <: AbstractFloat || throw(AlgorithmGraphError(
        "a HIL command input must use an AbstractFloat element type",
    ))
    prepared_command = _prepare_hil_host_buffer(
        graph_command,
        command_buffer,
        "HIL command",
    )
    prepared_frame = _prepare_hil_host_buffer(
        graph_frame,
        frame_buffer,
        "HIL frame",
    )
    prepared_command === graph_command && throw(AlgorithmGraphError(
        "the HIL command exchange buffer must not alias the active graph command",
    ))
    prepared_frame === graph_frame && throw(AlgorithmGraphError(
        "the HIL frame exchange buffer must not alias the graph output",
    ))
    prepared_command === prepared_frame && throw(AlgorithmGraphError(
        "HIL command and frame buffers must be distinct",
    ))
    initial_command = similar(prepared_command)
    _copy_hil_buffer!(graph, initial_command, graph_command)
    _validate_hil_command(initial_command)
    copyto!(prepared_command, initial_command)
    fill!(prepared_frame, zero(eltype(prepared_frame)))
    return PreparedGraphHILBoundary(
        graph,
        graph_command,
        graph_frame,
        prepared_command,
        prepared_frame,
        initial_command,
        GraphHILBoundaryState(),
    )
end

"""Return the caller-writable host command exchange buffer."""
@inline hil_command_buffer(boundary::PreparedGraphHILBoundary) =
    boundary.command_buffer

"""Return the host buffer containing the last completed graph frame."""
@inline hil_frame_buffer(boundary::PreparedGraphHILBoundary) =
    boundary.frame_buffer

"""Return the lockstep sequence and failure state without mutating it."""
@inline function hil_boundary_status(boundary::PreparedGraphHILBoundary)
    state = boundary.state
    return (
        frame_sequence=state.frame_sequence,
        active_command_sequence=state.active_command_sequence,
        command_response_required=state.command_response_required,
        failed=state.failed || graph_failed(boundary.graph),
    )
end

@noinline function _throw_hil_boundary_failed()
    throw(AlgorithmGraphError(
        "the graph HIL boundary is stopped after a failed operation; reset it before reuse",
    ))
end

@noinline function _throw_hil_command_response_required(sequence::UInt64)
    throw(AlgorithmGraphError(
        "frame $sequence requires a same-sequence command response before the next frame",
    ))
end

function _require_hil_frame_step(boundary::PreparedGraphHILBoundary)
    state = boundary.state
    (state.failed || graph_failed(boundary.graph)) &&
        _throw_hil_boundary_failed()
    state.command_response_required &&
        _throw_hil_command_response_required(state.frame_sequence)
    graph_step_sequence(boundary.graph) == state.frame_sequence || throw(
        AlgorithmGraphError(
            "the graph and HIL boundary frame sequences are not aligned",
        ),
    )
    return nothing
end

function _stage_completed_hil_frame!(boundary::PreparedGraphHILBoundary)
    _copy_hil_buffer!(
        boundary.graph,
        boundary.frame_buffer,
        boundary.frame_output,
    )
    state = boundary.state
    state.frame_sequence = graph_step_sequence(boundary.graph)
    state.command_response_required = true
    return state.frame_sequence
end

"""
    step_hil_frame!(boundary) -> sequence

Execute one complete graph frame with the active command, copy the completed
frame to the host exchange buffer, and return its positive sequence. The next
frame is blocked until [`adopt_hil_command!`](@ref) accepts that same sequence.
"""
function step_hil_frame!(boundary::PreparedGraphHILBoundary)
    _require_hil_frame_step(boundary)
    try
        step_graph!(boundary.graph)
        return _stage_completed_hil_frame!(boundary)
    catch
        boundary.state.failed = true
        rethrow()
    end
end

"""
    step_hil_frame_at!(boundary, driver) -> (sequence, timestamp)

Execute and stage one complete HIL frame at the next model-time boundary. The
graph, HIL boundary, and model-time sequences must be aligned. The driver is
caller-owned so wall-clock pacing and external timestamp capture remain outside
the lockstep boundary.
"""
function step_hil_frame_at!(
    boundary::PreparedGraphHILBoundary,
    driver::_ModelTimeDriver,
)
    _require_hil_frame_step(boundary)
    model_time_sequence(driver) == boundary.state.frame_sequence || throw(
        AlgorithmGraphError(
            "the model-time driver and HIL boundary frame sequences are not aligned",
        ),
    )
    try
        timestamp = step_graph_at!(boundary.graph, driver)
        sequence = _stage_completed_hil_frame!(boundary)
        return (sequence=sequence, timestamp=timestamp)
    catch
        boundary.state.failed = true
        rethrow()
    end
end

function _validate_hil_command(values::Array{T}) where {T<:AbstractFloat}
    @inbounds for index in eachindex(values)
        isfinite(values[index]) || throw(AlgorithmGraphError(
            "HIL command buffer contains a non-finite value at index $index",
        ))
    end
    return nothing
end

"""
    adopt_hil_command!(boundary, sequence) -> sequence

Validate the complete caller-written host command buffer and adopt it for the
next graph frame. `sequence` must match the outstanding frame. Validation
finishes before the exact-target graph command is mutated; a failed target copy
stops the boundary because it may have partially changed that buffer.
"""
function adopt_hil_command!(
    boundary::PreparedGraphHILBoundary,
    sequence::UInt64,
)
    state = boundary.state
    (state.failed || graph_failed(boundary.graph)) &&
        _throw_hil_boundary_failed()
    state.command_response_required || throw(AlgorithmGraphError(
        "no completed HIL frame is awaiting a command response",
    ))
    sequence == state.frame_sequence || throw(AlgorithmGraphError(
        "expected command sequence $(state.frame_sequence), received $sequence",
    ))
    _validate_hil_command(boundary.command_buffer)
    try
        _copy_hil_buffer!(
            boundary.graph,
            boundary.command_input,
            boundary.command_buffer,
        )
    catch
        state.failed = true
        rethrow()
    end
    state.active_command_sequence = sequence
    state.command_response_required = false
    return sequence
end

adopt_hil_command!(boundary::PreparedGraphHILBoundary, sequence) =
    throw(AlgorithmGraphError(
        "HIL command sequence must be UInt64, not $(typeof(sequence))",
    ))

"""
    reset_hil_boundary!(boundary)

Reset the graph, restore its snapshotted sequence-zero command, clear the host
frame product, and restart lockstep sequencing. A separately owned model-time
driver must also be reset before `step_hil_frame_at!` is used again.
"""
function reset_hil_boundary!(boundary::PreparedGraphHILBoundary)
    state = boundary.state
    state.failed = true
    reset_graph!(boundary.graph)
    copyto!(boundary.command_buffer, boundary.initial_command)
    _copy_hil_buffer!(
        boundary.graph,
        boundary.command_input,
        boundary.initial_command,
    )
    fill!(boundary.frame_buffer, zero(eltype(boundary.frame_buffer)))
    state.frame_sequence = UInt64(0)
    state.active_command_sequence = UInt64(0)
    state.command_response_required = false
    state.failed = false
    return boundary
end

function reset_hil_boundary!(
    boundary::PreparedGraphHILBoundary,
    driver::_ModelTimeDriver,
)
    reset_hil_boundary!(boundary)
    reset_model_time!(driver)
    return boundary
end

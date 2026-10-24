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
    CommandInput<:Union{AbstractArray,NamedTuple},
    FrameOutput<:AbstractArray,
    CommandBuffer<:Union{Array,NamedTuple},
    FrameBuffer<:Array,
    InitialCommand<:Union{Array,NamedTuple},
}
    graph::Graph
    command_input::CommandInput
    frame_output::FrameOutput
    command_buffer::CommandBuffer
    frame_buffer::FrameBuffer
    initial_command::InitialCommand
    state::GraphHILBoundaryState
end

function _prepare_hil_command_buffers(commands::NamedTuple, ::Nothing)
    return map(command -> _prepare_hil_host_buffer(command, nothing, "HIL command"), commands)
end

function _prepare_hil_command_buffers(commands::NamedTuple{Names}, buffers::NamedTuple) where {Names}
    keys(buffers) == Names || throw(AlgorithmGraphError(
        "HIL command buffer names must match command_input names and order $Names",
    ))
    prepared = map(values(commands), values(buffers)) do command, buffer
        _prepare_hil_host_buffer(command, buffer, "HIL command")
    end
    return NamedTuple{Names}(prepared)
end

function _prepare_hil_command_buffers(commands::NamedTuple, buffers)
    throw(AlgorithmGraphError("multiple HIL commands require a NamedTuple of host buffers"))
end

function _require_distinct_hil_buffers(buffers::Tuple)
    isempty(buffers) && return nothing
    first_buffer = first(buffers)
    for buffer in Base.tail(buffers)
        Base.mightalias(first_buffer, buffer) && throw(AlgorithmGraphError(
            "HIL graph command, graph frame, and exchange buffers must not share storage",
        ))
    end
    _require_distinct_hil_buffers(Base.tail(buffers))
    return nothing
end

function _prepare_hil_boundary_buffers(
    graph::PreparedAlgorithmGraph,
    command_input::Tuple,
    frame_output::Symbol,
    command_buffer,
    frame_buffer,
)
    (1 <= length(command_input) <= 16 && all(name -> name isa Symbol, command_input)) ||
        throw(AlgorithmGraphError("HIL command_input must contain 1 to 16 Symbols"))
    allunique(command_input) || throw(AlgorithmGraphError(
        "HIL command_input names must be unique",
    ))
    graph_failed(graph) && throw(AlgorithmGraphError(
        "a failed graph cannot be bound to a HIL boundary",
    ))
    graph_step_pending(graph) && throw(AlgorithmGraphError(
        "a graph with a pending frame cannot be bound to a HIL boundary",
    ))
    iszero(graph_step_sequence(graph)) || throw(AlgorithmGraphError(
        "a HIL boundary requires an unstepped graph",
    ))
    commands = NamedTuple{command_input}(map(name -> graph_input(graph, name), command_input))
    for command in values(commands)
        eltype(command) <: AbstractFloat || throw(AlgorithmGraphError(
            "a HIL command input must use an AbstractFloat element type",
        ))
    end
    frame = graph_output(graph, frame_output)
    prepared_commands = _prepare_hil_command_buffers(commands, command_buffer)
    prepared_frame = _prepare_hil_host_buffer(frame, frame_buffer, "HIL frame")
    _require_distinct_hil_buffers((values(commands)..., frame,
        values(prepared_commands)..., prepared_frame))
    initial_commands = map(similar, prepared_commands)
    _copy_hil_commands!(graph, initial_commands, commands)
    _validate_hil_command(initial_commands)
    _restore_hil_host_commands!(prepared_commands, initial_commands)
    fill!(prepared_frame, zero(eltype(prepared_frame)))
    return (commands, frame, prepared_commands, prepared_frame, initial_commands)
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
    GC.@preserve destination source begin
        _with_prepared_device_execution_context(graph.context) do
            try
                Backends._copy_prepared_device_buffer!(graph.context, destination, source)
            finally
                _synchronize_prepared_graph_execution_context!(
                    graph.execution,
                    graph.context,
                )
            end
        end
    end
    return destination
end

function _prepare_hil_boundary_buffers(
    graph::PreparedAlgorithmGraph,
    command_input::Symbol,
    frame_output::Symbol,
    command_buffer,
    frame_buffer,
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
    (Base.mightalias(prepared_command, graph_command) ||
     Base.mightalias(prepared_command, graph_frame)) && throw(AlgorithmGraphError(
        "the HIL command exchange buffer must not alias a graph command or frame buffer",
    ))
    (Base.mightalias(prepared_frame, graph_command) ||
     Base.mightalias(prepared_frame, graph_frame)) && throw(AlgorithmGraphError(
        "the HIL frame exchange buffer must not alias a graph command or frame buffer",
    ))
    Base.mightalias(prepared_command, prepared_frame) && throw(AlgorithmGraphError(
        "HIL command and frame buffers must not share storage",
    ))
    initial_command = similar(prepared_command)
    _copy_hil_buffer!(graph, initial_command, graph_command)
    _validate_hil_command(initial_command)
    copyto!(prepared_command, initial_command)
    fill!(prepared_frame, zero(eltype(prepared_frame)))
    return (
        graph_command,
        graph_frame,
        prepared_command,
        prepared_frame,
        initial_command,
    )
end

"""
    prepare_graph_hil_boundary(graph; command_input, frame_output,
                               command_buffer=nothing, frame_buffer=nothing)

Prepare a lockstep external-RTC boundary around an unstepped graph. `command_input`
is one `Symbol` or a tuple of 1 to 16 unique `Symbol`s. For the tuple form,
`command_buffer` is `nothing` or a `NamedTuple` of host arrays with the same names
and order; each command can have a different shape and floating-point type.
The named command inputs and frame output remain the graph's exact-target arrays.
The exchange buffers are ordinary host `Array`s, allocated during preparation
unless supplied by the caller. Preparation snapshots the initial graph commands
for deterministic reset. After preparation, the boundary must be the only
owner that steps the graph or mutates its exact command input; transport code
writes only [`hil_command_buffer`](@ref).
"""
function prepare_graph_hil_boundary(
    graph::PreparedAlgorithmGraph;
    command_input::Union{Symbol,Tuple},
    frame_output::Symbol,
    command_buffer=nothing,
    frame_buffer=nothing,
)
    buffers = _prepare_hil_boundary_buffers(
        graph, command_input, frame_output, command_buffer, frame_buffer,
    )
    return PreparedGraphHILBoundary(graph, buffers..., GraphHILBoundaryState())
end

"""Return the caller-writable host command Array, or NamedTuple of Arrays for multiple ports."""
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
    graph_step_pending(boundary.graph) && throw(AlgorithmGraphError(
        "a graph with a pending frame cannot advance a HIL boundary",
    ))
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

@inline function _validate_hil_command(commands::NamedTuple)
    foreach(_validate_hil_command, values(commands))
    return nothing
end

@inline function _copy_hil_commands!(
    graph::PreparedAlgorithmGraph,
    destination::AbstractArray,
    source::AbstractArray,
)
    _copy_hil_buffer!(graph, destination, source)
    return nothing
end

@inline function _copy_hil_commands!(
    graph::PreparedAlgorithmGraph,
    destinations::NamedTuple{Names},
    sources::NamedTuple{Names},
) where {Names}
    foreach(values(destinations), values(sources)) do destination, source
        _copy_hil_buffer!(graph, destination, source)
    end
    return nothing
end

@inline function _restore_hil_host_commands!(destination::Array, source::Array)
    copyto!(destination, source)
    return nothing
end

@inline function _restore_hil_host_commands!(
    destinations::NamedTuple{Names},
    sources::NamedTuple{Names},
) where {Names}
    foreach(copyto!, values(destinations), values(sources))
    return nothing
end

"""
    adopt_hil_command!(boundary, sequence) -> sequence

Validate every complete caller-written host command buffer and adopt them for the
next graph frame. `sequence` must match the outstanding frame. Validation
finishes for all ports before any exact-target graph command is mutated; a failed
target copy stops the boundary because commands may have been partially adopted.
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
    graph_step_pending(boundary.graph) && throw(AlgorithmGraphError(
        "a graph with a pending frame cannot adopt a HIL command",
    ))
    graph_step_sequence(boundary.graph) == state.frame_sequence || throw(AlgorithmGraphError(
        "the graph and HIL boundary frame sequences are not aligned",
    ))
    _validate_hil_command(boundary.command_buffer)
    try
        _copy_hil_commands!(
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
    _restore_hil_host_commands!(boundary.command_buffer, boundary.initial_command)
    _copy_hil_commands!(
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

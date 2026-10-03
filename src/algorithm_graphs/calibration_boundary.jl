mutable struct GraphCalibrationBoundaryState
    frame_sequence::UInt64
    active_probe_sequence::UInt64
    failed::Bool
end

GraphCalibrationBoundaryState() =
    GraphCalibrationBoundaryState(UInt64(0), UInt64(0), false)

"""
    PreparedGraphCalibrationBoundary

Transport-neutral, single-writer binding for exposures under an explicitly
adopted, held probe. Probe and graph-frame sequences are independent. Persistent
sequence and failure state is separate from the exact graph bindings and host
exchange buffers.
"""
struct PreparedGraphCalibrationBoundary{
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
    state::GraphCalibrationBoundaryState
end

"""
    prepare_graph_calibration_boundary(graph; command_input, frame_output,
                                       command_buffer=nothing, frame_buffer=nothing)

Bind an unstepped graph's exact command input and frame output to distinct host
`Array` exchange buffers using the same preparation contract as
[`prepare_graph_hil_boundary`](@ref). No probe is initially adopted. After
preparation this boundary must be the sole owner that steps the graph or mutates
its exact command input. Transport writes only [`hil_command_buffer`](@ref).
"""
function prepare_graph_calibration_boundary(
    graph::PreparedAlgorithmGraph;
    command_input::Symbol,
    frame_output::Symbol,
    command_buffer=nothing,
    frame_buffer=nothing,
)
    buffers = _prepare_hil_boundary_buffers(
        graph, command_input, frame_output, command_buffer, frame_buffer,
    )
    return PreparedGraphCalibrationBoundary(
        graph, buffers..., GraphCalibrationBoundaryState(),
    )
end

@inline hil_command_buffer(boundary::PreparedGraphCalibrationBoundary) =
    boundary.command_buffer

@inline hil_frame_buffer(boundary::PreparedGraphCalibrationBoundary) =
    boundary.frame_buffer

"""Return independent exposure/probe sequences and failure state."""
@inline function hil_boundary_status(boundary::PreparedGraphCalibrationBoundary)
    state = boundary.state
    return (
        frame_sequence=state.frame_sequence,
        active_probe_sequence=state.active_probe_sequence,
        failed=state.failed || graph_failed(boundary.graph),
    )
end

function _require_calibration_boundary_available(
    boundary::PreparedGraphCalibrationBoundary,
)
    state = boundary.state
    (state.failed || graph_failed(boundary.graph)) && _throw_hil_boundary_failed()
    graph_step_pending(boundary.graph) && throw(AlgorithmGraphError(
        "a pending graph step prevents calibration probe adoption or exposure",
    ))
    graph_step_sequence(boundary.graph) == state.frame_sequence || throw(
        AlgorithmGraphError(
            "the graph and calibration boundary frame sequences are not aligned",
        ),
    )
    return nothing
end

"""
    adopt_hil_probe!(boundary, probe_sequence::UInt64) -> probe_sequence

Validate the entire host command buffer and synchronously copy it to the exact
graph command input. The positive probe sequence must strictly advance; it is
independent of exposure sequences. A failed target copy stops the boundary
because the target may have changed partially. Successful adoption completes
before any subsequent exposure and holds the command across exposures.
"""
function adopt_hil_probe!(
    boundary::PreparedGraphCalibrationBoundary,
    probe_sequence::UInt64,
)
    _require_calibration_boundary_available(boundary)
    state = boundary.state
    probe_sequence > state.active_probe_sequence || throw(AlgorithmGraphError(
        "probe sequence must be positive and advance beyond $(state.active_probe_sequence), received $probe_sequence",
    ))
    _validate_hil_command(boundary.command_buffer)
    try
        _copy_hil_buffer!(boundary.graph, boundary.command_input, boundary.command_buffer)
    catch
        state.failed = true
        rethrow()
    end
    state.active_probe_sequence = probe_sequence
    return probe_sequence
end

adopt_hil_probe!(boundary::PreparedGraphCalibrationBoundary, probe_sequence) =
    throw(AlgorithmGraphError(
        "HIL probe sequence must be UInt64, not $(typeof(probe_sequence))",
    ))

function _require_hil_exposure_step(boundary::PreparedGraphCalibrationBoundary)
    _require_calibration_boundary_available(boundary)
    iszero(boundary.state.active_probe_sequence) && throw(AlgorithmGraphError(
        "a calibration exposure requires an explicitly adopted probe",
    ))
    return nothing
end

function _stage_completed_hil_exposure!(boundary::PreparedGraphCalibrationBoundary)
    _copy_hil_buffer!(boundary.graph, boundary.frame_buffer, boundary.frame_output)
    boundary.state.frame_sequence = graph_step_sequence(boundary.graph)
    return boundary.state.frame_sequence
end

"""
    step_hil_exposure!(boundary) -> sequence

Execute one complete graph exposure under the adopted probe, synchronously stage
its frame in the host exchange buffer, and return its actual graph sequence.
The probe remains active without a per-exposure command response. The caller
must finish consuming this frame before another exposure overwrites the buffer.
"""
function step_hil_exposure!(boundary::PreparedGraphCalibrationBoundary)
    _require_hil_exposure_step(boundary)
    try
        step_graph!(boundary.graph)
        return _stage_completed_hil_exposure!(boundary)
    catch
        boundary.state.failed = true
        rethrow()
    end
end

"""
    step_hil_exposure_at!(boundary, driver) -> (sequence, timestamp)

Execute and stage exactly one exposure at the next model-time boundary, returning
the graph sequence and actual model timestamp. Graph, boundary, and driver frame
sequences must align. The caller owns transport completion, settling policy,
exposure duration, and wall-clock pacing.
"""
function step_hil_exposure_at!(
    boundary::PreparedGraphCalibrationBoundary,
    driver::_ModelTimeDriver,
)
    _require_hil_exposure_step(boundary)
    model_time_sequence(driver) == boundary.state.frame_sequence || throw(
        AlgorithmGraphError(
            "the model-time driver and calibration boundary frame sequences are not aligned",
        ),
    )
    try
        timestamp = step_graph_at!(boundary.graph, driver)
        sequence = _stage_completed_hil_exposure!(boundary)
        return (sequence=sequence, timestamp=timestamp)
    catch
        boundary.state.failed = true
        rethrow()
    end
end

"""
    reset_hil_boundary!(boundary::PreparedGraphCalibrationBoundary[, driver])

Reset the graph, restore the initial command snapshot, clear the host frame, and
clear probe and frame sequences. A new explicit probe adoption is required.
Reset the separately owned model-time driver too when supplied. Reset is an
owner operation, not confirmation of physical probe restoration by a transport.
"""
function reset_hil_boundary!(boundary::PreparedGraphCalibrationBoundary)
    state = boundary.state
    state.failed = true
    reset_graph!(boundary.graph)
    copyto!(boundary.command_buffer, boundary.initial_command)
    _copy_hil_buffer!(boundary.graph, boundary.command_input, boundary.initial_command)
    fill!(boundary.frame_buffer, zero(eltype(boundary.frame_buffer)))
    state.frame_sequence = UInt64(0)
    state.active_probe_sequence = UInt64(0)
    state.failed = false
    return boundary
end

function reset_hil_boundary!(
    boundary::PreparedGraphCalibrationBoundary,
    driver::_ModelTimeDriver,
)
    reset_hil_boundary!(boundary)
    reset_model_time!(driver)
    return boundary
end

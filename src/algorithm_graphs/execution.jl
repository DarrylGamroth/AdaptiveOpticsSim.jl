@inline function _enqueue_node!(node::PreparedGraphNode)
    enqueue_graph_node!(node.owner)
    return nothing
end

@inline _enqueue_nodes!(::Tuple{}) = nothing
@inline function _enqueue_nodes!(nodes::Tuple)
    _enqueue_node!(first(nodes))
    _enqueue_nodes!(Base.tail(nodes))
    return nothing
end

@inline _commit_delayed_links!(::Tuple{}, ::Tuple{}) = nothing
@inline function _commit_delayed_links!(links::Tuple, state_values::Tuple)
    link = first(links)
    _copy_graph_buffer!(first(state_values), link.source)
    _commit_delayed_links!(Base.tail(links), Base.tail(state_values))
    return nothing
end

@inline _reset_delayed_links!(::Tuple{}, ::Tuple{}) = nothing
@inline function _reset_delayed_links!(links::Tuple, state_values::Tuple)
    link = first(links)
    _copy_graph_buffer!(first(state_values), link.initial)
    _reset_delayed_links!(Base.tail(links), Base.tail(state_values))
    return nothing
end

@inline _reset_nodes!(::Tuple{}) = nothing
@inline function _reset_nodes!(nodes::Tuple)
    reset_graph_node!(first(nodes).owner)
    _reset_nodes!(Base.tail(nodes))
    return nothing
end

"""
    GraphStepTicket

Single-use completion ticket for one submitted graph frame. The ticket retains
the exact prepared graph whose stream and mutable execution storage remain in
flight. Call [`wait_graph_step!`](@ref) before reading graph outputs, mutating
graph inputs, resetting the graph, or submitting another frame.
"""
struct GraphStepTicket{Graph<:PreparedAlgorithmGraph}
    graph::Graph
    sequence::UInt64
end

"""Return the number of graph steps committed since preparation or reset."""
@inline graph_step_sequence(graph::PreparedAlgorithmGraph) =
    graph.state.step_sequence

"""Return whether one submitted graph frame still owns the execution storage."""
@inline graph_step_pending(graph::PreparedAlgorithmGraph) = graph.state.pending

"""Return whether a node failure has stopped this graph run."""
@inline graph_failed(graph::PreparedAlgorithmGraph) = graph.state.failed

@noinline function _throw_graph_step_pending(sequence::UInt64)
    throw(AlgorithmGraphError(
        "algorithm graph frame $sequence is still pending; wait for its " *
        "completion ticket before reuse",
    ))
end

@noinline function _throw_stale_graph_step_ticket(sequence::UInt64)
    throw(AlgorithmGraphError(
        "algorithm graph completion ticket $sequence is no longer pending",
    ))
end

function _drain_failed_graph_submission!(graph::PreparedAlgorithmGraph)
    try
        _with_prepared_device_execution_context(graph.context) do
            _drain_prepared_graph_execution_contexts!(
                graph.execution,
                graph.context,
            )
        end
    catch
        # Preserve the submission exception. The graph remains fail-stop, so
        # reset is the only legal next mutation and performs another context
        # synchronization before changing execution storage.
    end
    return nothing
end

"""
    step_graph_async!(graph) -> GraphStepTicket

Submit one complete graph frame in declaration order without synchronizing the
prepared device context at the graph boundary. Submission has bounded capacity
one because the frame owns the graph's outputs, delayed-link storage, node
state/workspaces, and RNG state until its ticket is consumed.

Some node implementations may contain a required host or backend completion
boundary and therefore block during submission. The API promises deferred
graph-boundary completion, not that every admitted node is nonblocking.
"""
function step_graph_async!(graph::PreparedAlgorithmGraph)
    state = graph.state
    state.failed && throw(AlgorithmGraphError(
        "the algorithm graph is stopped after a failed step; reset it before reuse",
    ))
    state.pending && _throw_graph_step_pending(state.pending_sequence)
    state.step_sequence == typemax(UInt64) && throw(AlgorithmGraphError(
        "algorithm graph step sequence is exhausted",
    ))
    sequence = state.step_sequence + UInt64(1)
    state.pending_sequence = sequence
    state.pending = true
    try
        _preflight_prepared_graph_execution!(graph.execution, graph.nodes)
        _with_prepared_device_execution_context(graph.context) do
            _enqueue_prepared_graph_execution!(
                graph.execution,
                graph.nodes,
                graph.delayed_links,
                state.delayed_values,
                graph.context,
            )
        end
    catch
        state.failed = true
        _drain_failed_graph_submission!(graph)
        state.pending = false
        state.pending_sequence = UInt64(0)
        rethrow()
    end
    return GraphStepTicket(graph, sequence)
end

"""
    wait_graph_step!(ticket) -> graph

Wait for the exact submitted frame, publish its sequence, and release the
graph's mutable execution storage. A ticket is single-use. An asynchronous
device failure makes the graph fail-stop until [`reset_graph!`](@ref).
"""
function wait_graph_step!(ticket::GraphStepTicket)
    graph = ticket.graph
    state = graph.state
    state.pending || _throw_stale_graph_step_ticket(ticket.sequence)
    state.pending_sequence == ticket.sequence || throw(AlgorithmGraphError(
        "algorithm graph completion ticket $(ticket.sequence) does not match " *
        "pending frame $(state.pending_sequence)",
    ))
    try
        _with_prepared_device_execution_context(graph.context) do
            _synchronize_prepared_graph_execution_context!(
                graph.execution,
                graph.context,
            )
        end
        _complete_prepared_graph_execution!(graph.execution, graph.nodes)
    catch
        state.failed = true
        state.pending = false
        state.pending_sequence = UInt64(0)
        rethrow()
    end
    state.step_sequence = ticket.sequence
    state.pending = false
    state.pending_sequence = UInt64(0)
    return graph
end

"""
    step_graph!(graph)

Invoke every prepared node once in declaration order. Delayed links and the
step sequence publish only after the execution context reports completion. A
thrown node or device error stops the graph until [`reset_graph!`](@ref). This
synchronous
compatibility surface submits one frame with [`step_graph_async!`](@ref) and
consumes its completion ticket before returning.
"""
function step_graph!(graph::PreparedAlgorithmGraph)
    return wait_graph_step!(step_graph_async!(graph))
end

"""
    reset_graph!(graph)

Reset every prepared graph node through its adapter, restore all
delayed-link initial values, and clear graph failure and sequence state. The
graph remains failed if any reset operation throws.
"""
function reset_graph!(graph::PreparedAlgorithmGraph)
    state = graph.state
    state.pending && _throw_graph_step_pending(state.pending_sequence)
    state.failed = true
    _with_prepared_device_execution_context(graph.context) do
        try
            # A prior submission or completion failure may have left backend
            # work outstanding even though no usable ticket remains.
            _synchronize_prepared_graph_execution_context!(
                graph.execution,
                graph.context,
            )
            _reset_prepared_graph_execution!(
                graph.execution,
                graph.nodes,
                graph.delayed_links,
                state.delayed_values,
                graph.context,
            )
        finally
            _drain_prepared_graph_execution_contexts!(
                graph.execution,
                graph.context,
            )
        end
    end
    state.step_sequence = UInt64(0)
    state.pending_sequence = UInt64(0)
    state.pending = false
    state.failed = false
    return graph
end

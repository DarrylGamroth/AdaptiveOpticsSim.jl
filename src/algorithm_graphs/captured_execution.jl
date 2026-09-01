"""Prepared direct-stream execution for an algorithm graph."""
struct _PreparedStreamGraphExecution end

"""One node that remains directly enqueued under captured execution."""
struct _PreparedStreamGraphNode{Node}
    node::Node
end

"""One node retained by a native CUDA Graph or HIP Graph executable."""
struct _PreparedCapturedGraphNode{Captured}
    captured::Captured
end

"""Concrete mixed sequence of direct and captured graph-node executions."""
struct _PreparedCapturedGraphExecution{Executions<:Tuple}
    executions::Executions
    captured_count::Int
end

@inline function _prepare_graph_execution(
    ::StreamGraphExecution,
    nodes::Tuple,
    context,
)
    return _PreparedStreamGraphExecution()
end

@inline function _prepare_captured_graph_node(
    node,
    context,
    ::GraphNodeCaptureUnsupported,
)
    return _PreparedStreamGraphNode(node), 0
end

function _prepare_captured_graph_node(
    node,
    context,
    ::GraphNodeCaptureSafe,
)
    # Compile and exercise the exact enqueue path before capture. Reset and
    # synchronize restore the owner's initial scientific state and ensure no
    # warm-up work can still access its storage.
    enqueue_graph_node!(node.owner)
    _synchronize_prepared_device_execution_context!(context)
    reset_graph_node!(node.owner)
    _synchronize_prepared_device_execution_context!(context)

    captured = _capture_prepared_device_graph(context) do
        enqueue_graph_node!(node.owner)
        nothing
    end
    return _PreparedCapturedGraphNode(captured), 1
end

@inline function _prepare_captured_graph_nodes(
    ::Tuple{},
    context,
)
    return (), 0
end

@inline function _prepare_captured_graph_nodes(
    nodes::Tuple,
    context,
)
    node = first(nodes)
    prepared_node, captured = _prepare_captured_graph_node(
        node,
        context,
        graph_node_capture_capability(node.owner),
    )
    tail, tail_captured = _prepare_captured_graph_nodes(
        Base.tail(nodes),
        context,
    )
    return (prepared_node, tail...), captured + tail_captured
end

function _prepare_graph_execution(
    ::CapturedGraphExecution,
    nodes::Tuple,
    context,
)
    prepared, count = _prepare_captured_graph_nodes(nodes, context)
    count > 0 || throw(AlgorithmGraphError(
        "captured graph execution requires at least one explicitly " *
        "capture-safe node on a supported accelerator target",
    ))
    return _PreparedCapturedGraphExecution(prepared, count)
end

@inline graph_execution_policy(graph::PreparedAlgorithmGraph) =
    _graph_execution_policy(graph.execution)
@inline _graph_execution_policy(::_PreparedStreamGraphExecution) =
    StreamGraphExecution()
@inline _graph_execution_policy(::_PreparedCapturedGraphExecution) =
    CapturedGraphExecution()

@inline captured_graph_node_count(graph::PreparedAlgorithmGraph) =
    _captured_graph_node_count(graph.execution)
@inline _captured_graph_node_count(::_PreparedStreamGraphExecution) = 0
@inline _captured_graph_node_count(execution::_PreparedCapturedGraphExecution) =
    execution.captured_count

@inline function _enqueue_prepared_graph_execution!(
    ::_PreparedStreamGraphExecution,
    nodes::NamedTuple,
    context,
)
    _enqueue_nodes!(values(nodes))
    return nothing
end

@inline _enqueue_prepared_captured_nodes!(::Tuple{}, context) = nothing

@inline function _enqueue_prepared_captured_nodes!(
    nodes::Tuple,
    context,
)
    _enqueue_prepared_captured_node!(first(nodes), context)
    _enqueue_prepared_captured_nodes!(Base.tail(nodes), context)
    return nothing
end

@inline function _enqueue_prepared_captured_node!(
    node::_PreparedStreamGraphNode,
    context,
)
    _enqueue_node!(node.node)
    return nothing
end

@inline function _enqueue_prepared_captured_node!(
    node::_PreparedCapturedGraphNode,
    context,
)
    _launch_prepared_device_graph!(node.captured, context)
    return nothing
end

@inline function _enqueue_prepared_graph_execution!(
    execution::_PreparedCapturedGraphExecution,
    nodes::NamedTuple,
    context,
)
    _enqueue_prepared_captured_nodes!(execution.executions, context)
    return nothing
end

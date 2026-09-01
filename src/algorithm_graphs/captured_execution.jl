"""Prepared direct-stream execution for an algorithm graph."""
struct _PreparedStreamGraphExecution end

"""One complete prepared frame retained by a native device graph."""
struct _PreparedCapturedGraphExecution{Captured}
    captured::Captured
    captured_count::Int
end

@inline function _prepare_graph_execution(
    ::StreamGraphExecution,
    nodes::Tuple,
    delayed_links::Tuple,
    delayed_values::Tuple,
    context,
)
    return _PreparedStreamGraphExecution()
end

@inline _capture_safe_node_count(::Tuple{}) = 0

@inline _require_capture_safe_node(node, ::GraphNodeCaptureSafe) = 1

@noinline function _require_capture_safe_node(
    node,
    ::GraphNodeCaptureUnsupported,
)
    throw(AlgorithmGraphError(
        "captured graph execution requires every node to be capture-safe; " *
        "node $(_node_name(node)) with owner $(typeof(node.owner)) is unsupported",
    ))
end

@inline function _capture_safe_node_count(nodes::Tuple)
    node = first(nodes)
    return _require_capture_safe_node(
        node,
        graph_node_capture_capability(node.owner),
    ) + _capture_safe_node_count(Base.tail(nodes))
end

@inline _enqueue_captured_nodes!(::Tuple{}) = nothing

@inline function _enqueue_captured_nodes!(nodes::Tuple)
    enqueue_captured_graph_node!(first(nodes).owner)
    _enqueue_captured_nodes!(Base.tail(nodes))
    return nothing
end

@inline _preflight_captured_nodes!(::Tuple{}) = nothing

@inline function _preflight_captured_nodes!(nodes::Tuple)
    _preflight_captured_graph_node!(first(nodes).owner)
    _preflight_captured_nodes!(Base.tail(nodes))
    return nothing
end

@inline _complete_captured_nodes!(::Tuple{}) = nothing

@inline function _complete_captured_nodes!(nodes::Tuple)
    _complete_captured_graph_node!(first(nodes).owner)
    _complete_captured_nodes!(Base.tail(nodes))
    return nothing
end

function _prepare_graph_execution(
    ::CapturedGraphExecution,
    nodes::Tuple,
    delayed_links::Tuple,
    delayed_values::Tuple,
    context,
)
    count = _capture_safe_node_count(nodes)
    count > 0 || throw(AlgorithmGraphError(
        "captured graph execution requires at least one capture-safe node",
    ))

    # Compile and exercise the exact complete-frame enqueue path before
    # capture. Reset and synchronize restore the graph's initial scientific
    # state and ensure no warm-up work can still access retained storage.
    _enqueue_captured_nodes!(nodes)
    _commit_delayed_links!(delayed_links, delayed_values)
    _synchronize_prepared_device_execution_context!(context)
    _reset_nodes!(nodes)
    _reset_delayed_links!(delayed_links, delayed_values)
    _synchronize_prepared_device_execution_context!(context)

    captured = _capture_prepared_device_graph(context) do
        _enqueue_captured_nodes!(nodes)
        _commit_delayed_links!(delayed_links, delayed_values)
        nothing
    end
    return _PreparedCapturedGraphExecution(captured, count)
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

@inline function _synchronize_prepared_graph_execution_context!(
    ::_PreparedStreamGraphExecution,
    context,
)
    _synchronize_prepared_device_execution_context!(context)
    return nothing
end

@inline function _synchronize_prepared_graph_execution_context!(
    ::_PreparedCapturedGraphExecution,
    context,
)
    _synchronize_prepared_device_execution_context_blocking!(context)
    return nothing
end

@inline function _preflight_prepared_graph_execution!(
    ::_PreparedStreamGraphExecution,
    nodes::NamedTuple,
)
    return nothing
end

@inline function _preflight_prepared_graph_execution!(
    ::_PreparedCapturedGraphExecution,
    nodes::NamedTuple,
)
    _preflight_captured_nodes!(values(nodes))
    return nothing
end

@inline function _complete_prepared_graph_execution!(
    ::_PreparedStreamGraphExecution,
    nodes::NamedTuple,
)
    return nothing
end

@inline function _complete_prepared_graph_execution!(
    ::_PreparedCapturedGraphExecution,
    nodes::NamedTuple,
)
    _complete_captured_nodes!(values(nodes))
    return nothing
end

@inline function _enqueue_prepared_graph_execution!(
    ::_PreparedStreamGraphExecution,
    nodes::NamedTuple,
    delayed_links::Tuple,
    delayed_values::Tuple,
    context,
)
    _enqueue_nodes!(values(nodes))
    _commit_delayed_links!(delayed_links, delayed_values)
    return nothing
end

@inline function _enqueue_prepared_graph_execution!(
    execution::_PreparedCapturedGraphExecution,
    nodes::NamedTuple,
    delayed_links::Tuple,
    delayed_values::Tuple,
    context,
)
    _launch_prepared_device_graph!(execution.captured, context)
    return nothing
end

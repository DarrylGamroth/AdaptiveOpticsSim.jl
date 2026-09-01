"""Prepared direct-stream execution for an algorithm graph."""
struct _PreparedStreamGraphExecution end

"""Prepared bounded execution groups and their retained context lanes."""
struct _PreparedGroupedStreamGraphExecution{
    Policy<:GroupedStreamGraphExecution,
    Groups<:Tuple,
    Contexts<:FixedSizeVector,
    Events<:FixedSizeVector,
}
    policy::Policy
    groups::Groups
    contexts::Contexts
    events::Events
end

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
    preparation::_GraphExecutionPreparation,
)
    return _PreparedStreamGraphExecution()
end

@inline _take_graph_execution_nodes(nodes::Tuple, ::Val{0}) = ()
@inline function _take_graph_execution_nodes(
    nodes::Tuple,
    ::Val{N},
) where {N}
    return (
        first(nodes),
        _take_graph_execution_nodes(Base.tail(nodes), Val(N - 1))...,
    )
end

@inline _drop_graph_execution_nodes(nodes::Tuple, ::Val{0}) = nodes
@inline function _drop_graph_execution_nodes(
    nodes::Tuple,
    ::Val{N},
) where {N}
    return _drop_graph_execution_nodes(Base.tail(nodes), Val(N - 1))
end

@inline _prepare_grouped_graph_nodes(::Tuple{}, ::Tuple{}) = ()
@inline function _prepare_grouped_graph_nodes(
    groups::Tuple,
    nodes::Tuple,
)
    width = length(_graph_execution_group_names(first(groups)))
    return (
        _take_graph_execution_nodes(nodes, Val(width)),
        _prepare_grouped_graph_nodes(
            Base.tail(groups),
            _drop_graph_execution_nodes(nodes, Val(width)),
        )...,
    )
end

function _prepare_grouped_graph_events(contexts::FixedSizeVector)
    first_context = contexts[1]
    first_event = _with_prepared_device_execution_context(first_context) do
        _prepare_device_execution_event(first_context)
    end
    Event = typeof(first_event)
    builder = Vector{Event}(undef, length(contexts))
    builder[1] = first_event
    for index in 2:length(contexts)
        context = contexts[index]
        candidate = _with_prepared_device_execution_context(context) do
            _prepare_device_execution_event(context)
        end
        candidate isa Event || throw(AlgorithmGraphError(
            "grouped stream execution events must have one concrete type; " *
            "got $Event and $(typeof(candidate))",
        ))
        builder[index] = candidate
    end
    return FixedSizeVectorDefault{Event}(builder)
end

function _prepare_graph_execution(
    execution::GroupedStreamGraphExecution,
    nodes::Tuple,
    delayed_links::Tuple,
    delayed_values::Tuple,
    context,
    preparation::_GraphExecutionPreparation,
)
    groups = _prepare_grouped_graph_nodes(execution.groups, nodes)
    events = _prepare_grouped_graph_events(preparation.contexts)
    return _PreparedGroupedStreamGraphExecution(
        execution,
        groups,
        preparation.contexts,
        events,
    )
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
    preparation::_GraphExecutionPreparation,
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
@inline _graph_execution_policy(
    execution::_PreparedGroupedStreamGraphExecution,
) = execution.policy
@inline _graph_execution_policy(::_PreparedCapturedGraphExecution) =
    CapturedGraphExecution()

@inline captured_graph_node_count(graph::PreparedAlgorithmGraph) =
    _captured_graph_node_count(graph.execution)
@inline _captured_graph_node_count(::_PreparedStreamGraphExecution) = 0
@inline _captured_graph_node_count(
    ::_PreparedGroupedStreamGraphExecution,
) = 0
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
    ::_PreparedGroupedStreamGraphExecution,
    context,
)
    # Successful grouped submission joins every lane into the primary context
    # before delayed-state publication. Completing the primary context is
    # therefore sufficient on the ordinary completion path.
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

@inline function _drain_prepared_graph_execution_contexts!(
    ::_PreparedStreamGraphExecution,
    context,
)
    _synchronize_prepared_device_execution_context!(context)
    return nothing
end

@inline function _drain_prepared_graph_execution_contexts!(
    ::_PreparedCapturedGraphExecution,
    context,
)
    _synchronize_prepared_device_execution_context_blocking!(context)
    return nothing
end

function _drain_prepared_graph_execution_contexts!(
    execution::_PreparedGroupedStreamGraphExecution,
    context,
)
    first_error = nothing
    for lane_context in execution.contexts
        try
            _with_prepared_device_execution_context(lane_context) do
                _synchronize_prepared_device_execution_context!(lane_context)
            end
        catch error
            isnothing(first_error) && (first_error = error)
        end
    end
    isnothing(first_error) || throw(first_error)
    return nothing
end

@inline function _preflight_prepared_graph_execution!(
    ::_PreparedStreamGraphExecution,
    nodes::NamedTuple,
)
    return nothing
end

@inline function _preflight_prepared_graph_execution!(
    ::_PreparedGroupedStreamGraphExecution,
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
    ::_PreparedGroupedStreamGraphExecution,
    nodes::NamedTuple,
)
    return nothing
end

struct _EnqueueGroupedGraphNodes end
struct _ResetGroupedGraphNodes end

@inline _invoke_grouped_graph_node!(
    ::_EnqueueGroupedGraphNodes,
    node::PreparedGraphNode,
) = _enqueue_node!(node)

@inline function _invoke_grouped_graph_node!(
    ::_ResetGroupedGraphNodes,
    node::PreparedGraphNode,
)
    reset_graph_node!(node.owner)
    return nothing
end

@inline _wait_group_predecessor_events!(context, events, ::Val{0}) = nothing
@inline function _wait_group_predecessor_events!(
    context,
    events,
    ::Val{N},
) where {N}
    _wait_group_predecessor_events!(context, events, Val(N - 1))
    _wait_prepared_device_execution_event!(events[N], context)
    return nothing
end

@inline function _enqueue_group_predecessor_waits!(
    ::Tuple{},
    contexts,
    events,
    previous_width,
    lane::Int,
)
    return nothing
end

@inline function _enqueue_group_predecessor_waits!(
    nodes::Tuple,
    contexts,
    events,
    previous_width,
    lane::Int,
)
    context = contexts[lane]
    _with_prepared_device_execution_context(context) do
        _wait_group_predecessor_events!(
            context,
            events,
            previous_width,
        )
    end
    _enqueue_group_predecessor_waits!(
        Base.tail(nodes),
        contexts,
        events,
        previous_width,
        lane + 1,
    )
    return nothing
end

@inline function _enqueue_grouped_graph_lanes!(
    operation,
    ::Tuple{},
    contexts,
    events,
    lane::Int,
)
    return nothing
end

@inline function _enqueue_grouped_graph_lanes!(
    operation,
    nodes::Tuple,
    contexts,
    events,
    lane::Int,
)
    context = contexts[lane]
    event = events[lane]
    _with_prepared_device_execution_context(context) do
        _invoke_grouped_graph_node!(operation, first(nodes))
        _record_prepared_device_execution_event!(event, context)
    end
    _enqueue_grouped_graph_lanes!(
        operation,
        Base.tail(nodes),
        contexts,
        events,
        lane + 1,
    )
    return nothing
end

@inline function _enqueue_grouped_graph_sequence!(
    operation,
    ::Tuple{},
    contexts,
    events,
    previous_width,
)
    return previous_width
end

@inline function _enqueue_grouped_graph_sequence!(
    operation,
    groups::Tuple,
    contexts,
    events,
    previous_width,
)
    nodes = first(groups)
    # Every lane must capture the preceding event generation before any lane
    # re-records its reusable event for this group. Keeping waits and records in
    # separate host-submission phases preserves actual within-group overlap.
    _enqueue_group_predecessor_waits!(
        nodes,
        contexts,
        events,
        previous_width,
        1,
    )
    _enqueue_grouped_graph_lanes!(
        operation,
        nodes,
        contexts,
        events,
        1,
    )
    return _enqueue_grouped_graph_sequence!(
        operation,
        Base.tail(groups),
        contexts,
        events,
        Val(length(nodes)),
    )
end

@inline function _join_grouped_graph_lanes!(
    execution::_PreparedGroupedStreamGraphExecution,
    final_width,
)
    primary = execution.contexts[1]
    _with_prepared_device_execution_context(primary) do
        _wait_group_predecessor_events!(
            primary,
            execution.events,
            final_width,
        )
    end
    return nothing
end

@inline function _enqueue_grouped_graph_nodes!(
    operation,
    execution::_PreparedGroupedStreamGraphExecution,
)
    final_width = _enqueue_grouped_graph_sequence!(
        operation,
        execution.groups,
        execution.contexts,
        execution.events,
        Val(0),
    )
    _join_grouped_graph_lanes!(execution, final_width)
    return nothing
end

@inline function _complete_prepared_graph_execution!(
    ::_PreparedStreamGraphExecution,
    nodes::NamedTuple,
)
    return nothing
end

@inline function _enqueue_prepared_graph_execution!(
    execution::_PreparedGroupedStreamGraphExecution,
    nodes::NamedTuple,
    delayed_links::Tuple,
    delayed_values::Tuple,
    context,
)
    _enqueue_grouped_graph_nodes!(_EnqueueGroupedGraphNodes(), execution)
    _commit_delayed_links!(delayed_links, delayed_values)
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

@inline function _reset_prepared_graph_execution!(
    ::Union{_PreparedStreamGraphExecution,_PreparedCapturedGraphExecution},
    nodes::NamedTuple,
    delayed_links::Tuple,
    delayed_values::Tuple,
    context,
)
    _reset_nodes!(values(nodes))
    _reset_delayed_links!(delayed_links, delayed_values)
    return nothing
end

@inline function _reset_prepared_graph_execution!(
    execution::_PreparedGroupedStreamGraphExecution,
    nodes::NamedTuple,
    delayed_links::Tuple,
    delayed_values::Tuple,
    context,
)
    _enqueue_grouped_graph_nodes!(_ResetGroupedGraphNodes(), execution)
    _reset_delayed_links!(delayed_links, delayed_values)
    return nothing
end

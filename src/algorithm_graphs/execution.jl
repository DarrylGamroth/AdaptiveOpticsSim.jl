@inline function _process_node!(node::PreparedAlgorithmNode)
    _process_algorithm!(node.algorithm, node.outputs, node.inputs)
    return nothing
end

@inline _process_nodes!(::Tuple{}) = nothing
@inline function _process_nodes!(nodes::Tuple)
    _process_node!(first(nodes))
    _process_nodes!(Base.tail(nodes))
    return nothing
end

@inline _commit_delayed_links!(::Tuple{}, ::Tuple{}) = nothing
@inline function _commit_delayed_links!(links::Tuple, state_values::Tuple)
    link = first(links)
    _copy_algorithm_buffer!(
        link.destination_algorithm,
        first(state_values),
        link.source,
    )
    _commit_delayed_links!(Base.tail(links), Base.tail(state_values))
    return nothing
end

@inline _reset_delayed_links!(::Tuple{}, ::Tuple{}) = nothing
@inline function _reset_delayed_links!(links::Tuple, state_values::Tuple)
    link = first(links)
    _copy_algorithm_buffer!(
        link.destination_algorithm,
        first(state_values),
        link.initial,
    )
    _reset_delayed_links!(Base.tail(links), Base.tail(state_values))
    return nothing
end

@inline _reset_nodes!(::Tuple{}) = nothing
@inline function _reset_nodes!(nodes::Tuple)
    _reset_algorithm!(first(nodes).algorithm)
    _reset_nodes!(Base.tail(nodes))
    return nothing
end

"""Return the number of graph steps committed since preparation or reset."""
@inline graph_step_sequence(graph::PreparedAlgorithmGraph) =
    graph.state.step_sequence

"""Return whether a node failure has stopped this graph run."""
@inline graph_failed(graph::PreparedAlgorithmGraph) = graph.state.failed

"""
    step_graph!(graph)

Invoke every prepared node once in declaration order. Delayed links and the
step sequence commit only after every node returns normally. A thrown node
error stops the graph until [`reset_graph!`](@ref).
"""
function step_graph!(graph::PreparedAlgorithmGraph)
    state = graph.state
    state.failed && throw(AlgorithmGraphError(
        "the algorithm graph is stopped after a failed step; reset it before reuse",
    ))
    state.step_sequence == typemax(UInt64) && throw(AlgorithmGraphError(
        "algorithm graph step sequence is exhausted",
    ))
    try
        _process_nodes!(values(graph.nodes))
        _commit_delayed_links!(graph.delayed_links, state.delayed_values)
        state.step_sequence += UInt64(1)
    catch
        state.failed = true
        rethrow()
    end
    return graph
end

"""
    reset_graph!(graph)

Reset every prepared algorithm through its declaration adapter, restore all
delayed-link initial values, and clear graph failure and sequence state. The
graph remains failed if any reset operation throws.
"""
function reset_graph!(graph::PreparedAlgorithmGraph)
    state = graph.state
    state.failed = true
    _reset_nodes!(values(graph.nodes))
    _reset_delayed_links!(graph.delayed_links, state.delayed_values)
    state.step_sequence = UInt64(0)
    state.failed = false
    return graph
end

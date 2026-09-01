"""A structured preparation or execution failure in a portable algorithm graph."""
struct AlgorithmGraphError <: AdaptiveOpticsSimError
    msg::String
end

"""
    StreamGraphExecution()

Execute every algorithm node directly on the retained host or accelerator
stream. This is the default graph execution policy.
"""
struct StreamGraphExecution end

"""
    CapturedGraphExecution()

Capture every explicitly eligible node as a backend command graph during
preparation. Nodes that do not opt in continue to execute directly on the same
retained stream, preserving declaration order. Preparation fails if the target
cannot capture any eligible node or if a claimed-eligible node cannot be
captured safely.
"""
struct CapturedGraphExecution end

"""Trait value for a graph-node owner that must execute directly."""
struct GraphNodeCaptureUnsupported end

"""
Trait value for a graph-node owner whose repeated enqueue operation is safe to
record and replay as a fixed-address accelerator command graph. Reset must
preserve every recorded storage address.
"""
struct GraphNodeCaptureSafe end

"""A typed reference to one named port on one named algorithm node."""
struct AlgorithmPortReference{Node,Port} end

@inline _node_name(::AlgorithmPortReference{Node}) where {Node} = Node
@inline _port_name(::AlgorithmPortReference{Node,Port}) where {Node,Port} = Port

function _require_graph_name(name::Symbol, role::AbstractString)
    isempty(String(name)) && throw(AlgorithmGraphError("$role name must not be empty"))
    return name
end

function AlgorithmPortReference(node::Symbol, port::Symbol)
    _require_graph_name(node, "algorithm node")
    _require_graph_name(port, "algorithm port")
    return AlgorithmPortReference{node,port}()
end

AlgorithmPortReference(reference::Pair{Symbol,Symbol}) =
    AlgorithmPortReference(first(reference), last(reference))

_port_reference(reference::AlgorithmPortReference) = reference
_port_reference(reference::Pair{Symbol,Symbol}) = AlgorithmPortReference(reference)
_port_reference(reference) = throw(AlgorithmGraphError(
    "algorithm port references use :node => :port, not $(typeof(reference))",
))

"""One cold typed graph-node declaration, configuration, and initial props."""
struct AlgorithmNodeDefinition{Name,Node,Config,Props}
    config::Config
    props::Props
end

@inline _node_name(::AlgorithmNodeDefinition{Name}) where {Name} = Name
@inline _node_type(
    ::AlgorithmNodeDefinition{Name,Node},
) where {Name,Node} = Node

"""
    algorithm_node(name, Node, config; props=NamedTuple())

Declare one named graph node backed by an AOS graph-node adapter. `config`
fixes construction and graph-rebuild values. `props` supplies scalar initial
values without changing the node's port contract. Node order in
[`algorithm_graph`](@ref) is direct-link execution order.
"""
function algorithm_node(
    name::Symbol,
    ::Type{Node},
    config::Config;
    props::Props=NamedTuple(),
) where {Node,Config,Props}
    _require_graph_name(name, "algorithm node")
    return AlgorithmNodeDefinition{name,Node,Config,Props}(config, props)
end

"""A same-step output-to-input connection."""
struct AlgorithmLink{Source<:AlgorithmPortReference,Destination<:AlgorithmPortReference}
    source::Source
    destination::Destination
end

"""
    link(source, destination)

Connect an output to an input in the same graph step. Both references use
`:node => :port`. The source node must precede the destination node.
"""
function link(source, destination)
    prepared_source = _port_reference(source)
    prepared_destination = _port_reference(destination)
    return AlgorithmLink(prepared_source, prepared_destination)
end

"""A prior-successful-step output-to-input connection and its initial value."""
struct DelayedAlgorithmLink{
    Source<:AlgorithmPortReference,
    Destination<:AlgorithmPortReference,
    Initial,
}
    source::Source
    destination::Destination
    initial::Initial
end

"""
    delayed_link(source, destination, initial)

Connect an output to an input through one successful graph-step delay. The
initial ndarray is required explicitly and is snapshotted during preparation.
"""
function delayed_link(source, destination, initial)
    prepared_source = _port_reference(source)
    prepared_destination = _port_reference(destination)
    return DelayedAlgorithmLink(
        prepared_source,
        prepared_destination,
        initial,
    )
end

"""One caller-owned graph input bound to an algorithm input port."""
struct GraphInputDefinition{
    Name,
    Destination<:AlgorithmPortReference,
    Values,
}
    destination::Destination
    values::Values
end

@inline _boundary_name(::GraphInputDefinition{Name}) where {Name} = Name

"""
    graph_input(name, destination, values)

Bind caller-owned packed ndarray storage to one graph input. The storage is
validated and retained by identity; preparation does not copy it.
"""
function graph_input(name::Symbol, destination, values)
    _require_graph_name(name, "graph input")
    prepared_destination = _port_reference(destination)
    return GraphInputDefinition{name,typeof(prepared_destination),typeof(values)}(
        prepared_destination,
        values,
    )
end

"""One caller-visible graph output aliasing an algorithm output port."""
struct GraphOutputDefinition{Name,Source<:AlgorithmPortReference}
    source::Source
end

@inline _boundary_name(::GraphOutputDefinition{Name}) where {Name} = Name

"""
    graph_output(name, source)

Expose one algorithm output as a named graph result without copying it.
"""
function graph_output(name::Symbol, source)
    _require_graph_name(name, "graph output")
    prepared_source = _port_reference(source)
    return GraphOutputDefinition{name,typeof(prepared_source)}(prepared_source)
end

"""One initial sparse-parameter replacement for a named algorithm port."""
struct SparseParameterBinding{Endpoint<:AlgorithmPortReference,Values}
    endpoint::Endpoint
    values::Values
end

"""
    sparse_parameter(endpoint, values)

Supply one startup ndarray value for a sparse-parameter port. The complete
graph is still inactive while all parameter values are validated and bound.
"""
function sparse_parameter(endpoint, values)
    prepared_endpoint = _port_reference(endpoint)
    return SparseParameterBinding(prepared_endpoint, values)
end

"""Cold static topology for one portable algorithm graph."""
struct AlgorithmGraphDefinition{Nodes,Inputs,Outputs,Links,Delays,Parameters}
    name::Symbol
    nodes::Nodes
    inputs::Inputs
    outputs::Outputs
    links::Links
    delayed_links::Delays
    parameters::Parameters
end

@inline _validate_node_definition(::AlgorithmNodeDefinition) = nothing
@inline _validate_node_definition(value) = throw(AlgorithmGraphError(
    "graph nodes must be AlgorithmNodeDefinition values, not $(typeof(value))",
))

@inline _validate_graph_input(::GraphInputDefinition) = nothing
@inline _validate_graph_input(value) = throw(AlgorithmGraphError(
    "graph inputs must be GraphInputDefinition values, not $(typeof(value))",
))

@inline _validate_graph_output(::GraphOutputDefinition) = nothing
@inline _validate_graph_output(value) = throw(AlgorithmGraphError(
    "graph outputs must be GraphOutputDefinition values, not $(typeof(value))",
))

@inline _validate_link(::AlgorithmLink) = nothing
@inline _validate_link(value) = throw(AlgorithmGraphError(
    "direct links must be AlgorithmLink values, not $(typeof(value))",
))

@inline _validate_delayed_link(::DelayedAlgorithmLink) = nothing
@inline _validate_delayed_link(value) = throw(AlgorithmGraphError(
    "delayed links must be DelayedAlgorithmLink values, not $(typeof(value))",
))

@inline _validate_parameter(::SparseParameterBinding) = nothing
@inline _validate_parameter(value) = throw(AlgorithmGraphError(
    "parameters must be SparseParameterBinding values, not $(typeof(value))",
))

@inline _validate_tuple(::Tuple{}, validator) = nothing
@inline function _validate_tuple(values::Tuple, validator)
    validator(first(values))
    _validate_tuple(Base.tail(values), validator)
    return nothing
end

@inline _node_names(::Tuple{}) = ()
@inline _node_names(nodes::Tuple) =
    (_node_name(first(nodes)), _node_names(Base.tail(nodes))...)

@inline _boundary_names(::Tuple{}) = ()
@inline _boundary_names(boundaries::Tuple) =
    (_boundary_name(first(boundaries)), _boundary_names(Base.tail(boundaries))...)

function _require_unique_names(names::Tuple, role::AbstractString)
    length(names) == length(unique(names)) ||
        throw(AlgorithmGraphError("$role names must be unique"))
    return nothing
end

"""
    algorithm_graph(nodes; name=:anonymous, inputs=(), outputs=(), links=(),
                    delayed_links=(), parameters=())

Construct a cold graph definition from concrete tuples. Direct-link execution
order is the order of `nodes`; feedback must use an explicit delayed link.
"""
function algorithm_graph(
    nodes::Tuple;
    name::Symbol=:anonymous,
    inputs::Tuple=(),
    outputs::Tuple=(),
    links::Tuple=(),
    delayed_links::Tuple=(),
    parameters::Tuple=(),
)
    _require_graph_name(name, "algorithm graph")
    isempty(nodes) && throw(AlgorithmGraphError(
        "an algorithm graph requires at least one node",
    ))
    _validate_tuple(nodes, _validate_node_definition)
    _validate_tuple(inputs, _validate_graph_input)
    _validate_tuple(outputs, _validate_graph_output)
    _validate_tuple(links, _validate_link)
    _validate_tuple(delayed_links, _validate_delayed_link)
    _validate_tuple(parameters, _validate_parameter)
    _require_unique_names(_node_names(nodes), "algorithm node")
    _require_unique_names(_boundary_names(inputs), "graph input")
    _require_unique_names(_boundary_names(outputs), "graph output")
    return AlgorithmGraphDefinition(
        name,
        nodes,
        inputs,
        outputs,
        links,
        delayed_links,
        parameters,
    )
end

"""Return the stable name of a cold or prepared algorithm graph."""
@inline graph_name(graph::AlgorithmGraphDefinition) = graph.name

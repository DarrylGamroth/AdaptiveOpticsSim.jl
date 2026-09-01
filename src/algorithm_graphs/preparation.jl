struct _GraphPortFormat{T,N}
    shape::NTuple{N,Int}
    schema::String
    layout::Symbol
end

struct _GraphPortContract{Name,Direction,Role,Format<:_GraphPortFormat}
    format::Format
end

@inline _port_name(::_GraphPortContract{Name}) where {Name} = Name
@inline _port_direction(
    ::_GraphPortContract{Name,Direction},
) where {Name,Direction} = Direction
@inline _port_role(
    ::_GraphPortContract{Name,Direction,Role},
) where {Name,Direction,Role} = Role

"""
    graph_port_contract(name, direction, role, T, shape, schema, layout)

Construct one validated graph-side port contract for a graph-node adapter.
`direction` is `:input` or `:output`; `role` is `:data` or `:parameter`; the
portable allocator currently admits `:column_major` layout.
"""
function graph_port_contract(
    name::Symbol,
    direction::Symbol,
    role::Symbol,
    ::Type{T},
    shape::NTuple{N,<:Integer},
    schema::AbstractString,
    layout::Symbol,
) where {T,N}
    _require_graph_name(name, "algorithm port")
    isconcretetype(T) && isbitstype(T) || throw(AlgorithmGraphError(
        "algorithm port $name element type $T must be concrete and isbits",
    ))
    direction in (:input, :output) || throw(AlgorithmGraphError(
        "algorithm port $name has unsupported direction $direction",
    ))
    role in (:data, :parameter) || throw(AlgorithmGraphError(
        "algorithm port $name has unsupported role $role",
    ))
    direction === :output && role === :parameter &&
        throw(AlgorithmGraphError("sparse parameter $name must be an input"))
    isempty(schema) && throw(AlgorithmGraphError(
        "algorithm port $name scientific schema must not be empty",
    ))
    converted = ntuple(Val(N)) do index
        dimension = shape[index]
        dimension isa Bool && throw(AlgorithmGraphError(
            "algorithm port $name dimensions must be integers, not Bool",
        ))
        zero(dimension) < dimension <= typemax(Int) || throw(
            AlgorithmGraphError(
                "algorithm port $name dimensions must be positive and fit " *
                "the host index range",
            ),
        )
        return Int(dimension)
    end
    format = _GraphPortFormat{T,N}(converted, String(schema), layout)
    return _GraphPortContract{name,direction,role,typeof(format)}(format)
end

@inline _validate_port_contract(::_GraphPortContract) = nothing
@inline _validate_port_contract(value) = throw(AlgorithmGraphError(
    "graph-node adapters must return graph port contracts, not $(typeof(value))",
))

@inline _port_names(::Tuple{}) = ()
@inline _port_names(ports::Tuple) =
    (_port_name(first(ports)), _port_names(Base.tail(ports))...)

function _named_tuple(names::Tuple, values::Tuple, role::AbstractString)
    length(names) == length(values) || throw(AlgorithmGraphError(
        "$role name and value counts differ",
    ))
    _require_unique_names(names, role)
    return NamedTuple{names}(values)
end

"""
    graph_node_ports(Node, config)

Return the fixed tuple of graph-side port contracts for one graph-node type and
construction configuration. Port contracts must not depend on initial or live
props.
"""
function graph_node_ports(::Type{Node}, config) where {Node}
    throw(AlgorithmGraphError(
        "no AlgorithmGraphs port adapter is loaded for node type $Node",
    ))
end

"""
    prepare_graph_node(Node, config, props, inputs, outputs, parameters, target)

Prepare one exact single-writer node owner after all graph storage, connections,
and startup sparse parameters have been admitted. The returned owner binds its
run-immutable plan, persistent state, replaceable workspace, and exact named
buffers without changing their identities.
"""
function prepare_graph_node(
    ::Type{Node},
    config,
    props,
    inputs::NamedTuple,
    outputs::NamedTuple,
    parameters::NamedTuple,
    target,
) where {Node}
    throw(AlgorithmGraphError(
        "no AlgorithmGraphs preparation adapter is loaded for node type $Node",
    ))
end

"""
    step_graph_node!(owner)

Execute one prepared graph-node owner against the exact buffers retained during
preparation. Implementations must perform bounded repeated-path work and must
not replace their input, output, parameter, state, or workspace storage.
"""
function step_graph_node!(owner)
    throw(AlgorithmGraphError(
        "the AlgorithmGraphs adapter for $(typeof(owner)) does not implement stepping",
    ))
end

"""
    enqueue_graph_node!(owner)

Submit one prepared graph-node owner without adding a graph-boundary
completion wait. The default calls [`step_graph_node!`](@ref), so adapters need
no second implementation unless they can safely preserve same-context ordering
while deferring a narrower completion boundary.
"""
@inline function enqueue_graph_node!(owner)
    step_graph_node!(owner)
    return nothing
end

"""
    enqueue_captured_graph_node!(owner)

Submit the fixed-address device operation recorded for one capture-safe graph
node. The default reuses [`enqueue_graph_node!`](@ref). Stateful adapters may
specialize this seam when ordinary stream submission also performs host-side
publication that native command-graph replay cannot repeat.
"""
@inline function enqueue_captured_graph_node!(owner)
    enqueue_graph_node!(owner)
    return nothing
end

"""
    reset_graph_node!(owner)

Reset one prepared graph-node owner at a serialized graph boundary.
"""
function reset_graph_node!(owner)
    throw(AlgorithmGraphError(
        "the AlgorithmGraphs adapter for $(typeof(owner)) does not implement reset",
    ))
end

"""
    graph_node_capture_capability(owner)

Return `GraphNodeCaptureSafe()` only when `enqueue_graph_node!` submits a
fixed-address, replay-equivalent accelerator command sequence without host
state mutation, allocation, synchronization, result queries, or pointer
replacement, and when reset preserves those addresses. The conservative
default is `GraphNodeCaptureUnsupported()`.
"""
@inline graph_node_capture_capability(owner) =
    GraphNodeCaptureUnsupported()

function _allocate_graph_buffer(
    target::AbstractComputeDevice,
    format::_GraphPortFormat{T,N},
) where {T,N}
    format.layout === :column_major || throw(AlgorithmGraphError(
        "the native graph allocator supports only column-major buffers; got " *
        "$(format.layout)",
    ))
    values = allocate_device_array(target, T, format.shape...)
    fill!(values, zero(T))
    return values
end

function _validate_graph_buffer(
    target::AbstractComputeDevice,
    format::_GraphPortFormat{T,N},
    values,
) where {T,N}
    format.layout === :column_major || throw(AlgorithmGraphError(
        "the native graph supports only column-major buffers; got " *
        "$(format.layout)",
    ))
    storage_type = array_backend_type(compute_device_backend(target))
    values isa storage_type || throw(AlgorithmGraphError(
        "buffer $(typeof(values)) is not native packed storage for graph " *
        "target $target; expected $storage_type",
    ))
    eltype(values) === T || throw(AlgorithmGraphError(
        "buffer element type $(eltype(values)) does not match declared $T",
    ))
    ndims(values) == N || throw(AlgorithmGraphError(
        "buffer rank $(ndims(values)) does not match declared rank $N",
    ))
    size(values) == format.shape || throw(AlgorithmGraphError(
        "buffer shape $(size(values)) does not match declared shape $(format.shape)",
    ))
    actual_target = compute_device(values)
    actual_target == target || throw(AlgorithmGraphError(
        "buffer target $actual_target does not match prepared graph target $target",
    ))
    return values
end

@inline function _copy_graph_buffer!(destination, source)
    copyto_backend_async!(destination, source)
    return destination
end

struct _AdmittedGraphNode{Name,Node,Config,Props,Ports<:NamedTuple}
    config::Config
    props::Props
    ports::Ports
end

@inline _node_name(::_AdmittedGraphNode{Name}) where {Name} = Name

function _admit_graph_node(
    definition::AlgorithmNodeDefinition{Name,Node},
) where {Name,Node}
    port_values = graph_node_ports(Node, definition.config)
    _validate_tuple(port_values, _validate_port_contract)
    ports = _named_tuple(_port_names(port_values), port_values, "algorithm port")
    return _AdmittedGraphNode{
        Name,
        Node,
        typeof(definition.config),
        typeof(definition.props),
        typeof(ports),
    }(
        definition.config,
        definition.props,
        ports,
    )
end

@inline _admit_graph_nodes(::Tuple{}) = ()
@inline function _admit_graph_nodes(definitions::Tuple)
    return (
        _admit_graph_node(first(definitions)),
        _admit_graph_nodes(Base.tail(definitions))...,
    )
end

function _admitted_graph_nodes(definition::AlgorithmGraphDefinition)
    values = _admit_graph_nodes(definition.nodes)
    return _named_tuple(_node_names(definition.nodes), values, "algorithm node")
end

function _require_node(nodes::NamedTuple, ::Val{Name}) where {Name}
    hasproperty(nodes, Name) || throw(AlgorithmGraphError(
        "algorithm graph has no node named $Name",
    ))
    return getproperty(nodes, Name)
end

function _require_port(node::_AdmittedGraphNode, ::Val{Name}) where {Name}
    hasproperty(node.ports, Name) || throw(AlgorithmGraphError(
        "algorithm node $(_node_name(node)) has no port named $Name",
    ))
    return getproperty(node.ports, Name)
end

function _endpoint_contract(
    nodes::NamedTuple,
    reference::AlgorithmPortReference{Node,Port},
) where {Node,Port}
    node = _require_node(nodes, Val(Node))
    return node, _require_port(node, Val(Port))
end

@inline _require_data_input(
    ::_GraphPortContract{Name,:input,:data},
    reference,
) where {Name} = nothing
@inline _require_data_input(contract::_GraphPortContract, reference) =
    throw(AlgorithmGraphError(
        "$(_node_name(reference)) => $(_port_name(reference)) is not a frame-data input",
    ))

@inline _require_data_output(
    ::_GraphPortContract{Name,:output,:data},
    reference,
) where {Name} = nothing
@inline _require_data_output(contract::_GraphPortContract, reference) =
    throw(AlgorithmGraphError(
        "$(_node_name(reference)) => $(_port_name(reference)) is not a frame-data output",
    ))

@inline _require_sparse_parameter(
    ::_GraphPortContract{Name,:input,:parameter},
    reference,
) where {Name} = nothing
@inline _require_sparse_parameter(contract::_GraphPortContract, reference) =
    throw(AlgorithmGraphError(
        "$(_node_name(reference)) => $(_port_name(reference)) is not a sparse parameter",
    ))

function _require_matching_formats(
    source::_GraphPortFormat{T,N},
    destination::_GraphPortFormat{T,N},
    context::AbstractString,
) where {T,N}
    source.shape == destination.shape || throw(AlgorithmGraphError(
        "$context shape mismatch: $(source.shape) != $(destination.shape)",
    ))
    source.schema == destination.schema || throw(AlgorithmGraphError(
        "$context scientific schema mismatch: $(source.schema) != $(destination.schema)",
    ))
    source.layout === destination.layout || throw(AlgorithmGraphError(
        "$context packed layout mismatch: $(source.layout) != $(destination.layout)",
    ))
    return nothing
end

function _require_matching_formats(source, destination, context::AbstractString)
    throw(AlgorithmGraphError(
        "$context element type or rank does not match",
    ))
end

@inline _parameter_endpoint_keys(::Tuple{}) = ()
@inline _parameter_endpoint_keys(parameters::Tuple) = (
    (_node_name(first(parameters).endpoint), _port_name(first(parameters).endpoint)),
    _parameter_endpoint_keys(Base.tail(parameters))...,
)

function _validate_parameters!(
    nodes::NamedTuple,
    parameters::Tuple,
    target::AbstractComputeDevice,
)
    keys = _parameter_endpoint_keys(parameters)
    length(keys) == length(unique(keys)) || throw(AlgorithmGraphError(
        "each sparse parameter port may be initialized at most once",
    ))
    _validate_parameters_recursive!(nodes, parameters, target)
    return nothing
end

@inline _validate_parameters_recursive!(nodes, ::Tuple{}, target) = nothing
@inline function _validate_parameters_recursive!(
    nodes,
    parameters::Tuple,
    target,
)
    binding = first(parameters)
    _, contract = _endpoint_contract(nodes, binding.endpoint)
    _require_sparse_parameter(contract, binding.endpoint)
    _validate_graph_buffer(
        target,
        contract.format,
        binding.values,
    )
    _validate_parameters_recursive!(nodes, Base.tail(parameters), target)
    return nothing
end

@inline _data_output_contracts(::Tuple{}) = ()
@inline function _data_output_contracts(ports::Tuple)
    head = first(ports)
    tail = _data_output_contracts(Base.tail(ports))
    return _prepend_data_output(head, tail)
end

@inline _prepend_data_output(
    port::_GraphPortContract{Name,:output,:data},
    tail::Tuple,
) where {Name} = (port, tail...)
@inline _prepend_data_output(port::_GraphPortContract, tail::Tuple) = tail

@inline _data_input_contracts(::Tuple{}) = ()
@inline function _data_input_contracts(ports::Tuple)
    head = first(ports)
    tail = _data_input_contracts(Base.tail(ports))
    return _prepend_data_input(head, tail)
end

@inline _prepend_data_input(
    port::_GraphPortContract{Name,:input,:data},
    tail::Tuple,
) where {Name} = (port, tail...)
@inline _prepend_data_input(port::_GraphPortContract, tail::Tuple) = tail

@inline _parameter_contracts(::Tuple{}) = ()
@inline function _parameter_contracts(ports::Tuple)
    head = first(ports)
    tail = _parameter_contracts(Base.tail(ports))
    return _prepend_parameter(head, tail)
end

@inline _prepend_parameter(
    port::_GraphPortContract{Name,:input,:parameter},
    tail::Tuple,
) where {Name} = (port, tail...)
@inline _prepend_parameter(port::_GraphPortContract, tail::Tuple) = tail

@inline _allocate_output_buffers(::Tuple{}, target) = ()
@inline function _allocate_output_buffers(contracts::Tuple, target)
    return (
        _allocate_graph_buffer(target, first(contracts).format),
        _allocate_output_buffers(
            Base.tail(contracts),
            target,
        )...,
    )
end

function _allocate_node_outputs(node::_AdmittedGraphNode, target)
    contracts = _data_output_contracts(values(node.ports))
    buffers = _allocate_output_buffers(contracts, target)
    return _named_tuple(_port_names(contracts), buffers, "algorithm output port")
end

@inline _allocate_all_node_outputs(::Tuple{}, target) = ()
@inline function _allocate_all_node_outputs(nodes::Tuple, target)
    return (
        _allocate_node_outputs(first(nodes), target),
        _allocate_all_node_outputs(Base.tail(nodes), target)...,
    )
end

function _prepared_node_outputs(nodes::NamedTuple, target)
    groups = _allocate_all_node_outputs(values(nodes), target)
    return NamedTuple{keys(nodes)}(groups)
end

function _endpoint_output(
    node_outputs::NamedTuple,
    reference::AlgorithmPortReference{Node,Port},
) where {Node,Port}
    hasproperty(node_outputs, Node) || throw(AlgorithmGraphError(
        "algorithm graph has no output owner named $Node",
    ))
    outputs = getproperty(node_outputs, Node)
    hasproperty(outputs, Port) || throw(AlgorithmGraphError(
        "algorithm node $Node has no frame-data output named $Port",
    ))
    return getproperty(outputs, Port)
end

struct _PreparedGraphInput{Name,Destination,Values}
    destination::Destination
    values::Values
end

@inline _boundary_name(::_PreparedGraphInput{Name}) where {Name} = Name

function _prepare_graph_input(
    definition::GraphInputDefinition{Name},
    nodes::NamedTuple,
    target::AbstractComputeDevice,
) where {Name}
    _, contract = _endpoint_contract(nodes, definition.destination)
    _require_data_input(contract, definition.destination)
    _validate_graph_buffer(
        target,
        contract.format,
        definition.values,
    )
    return _PreparedGraphInput{
        Name,
        typeof(definition.destination),
        typeof(definition.values),
    }(definition.destination, definition.values)
end

@inline _prepare_graph_inputs(::Tuple{}, nodes, target) = ()
@inline function _prepare_graph_inputs(definitions::Tuple, nodes, target)
    return (
        _prepare_graph_input(first(definitions), nodes, target),
        _prepare_graph_inputs(Base.tail(definitions), nodes, target)...,
    )
end

function _prepared_graph_inputs(
    definition::AlgorithmGraphDefinition,
    nodes,
    target,
)
    prepared = _prepare_graph_inputs(definition.inputs, nodes, target)
    return _named_tuple(_boundary_names(definition.inputs), prepared, "graph input")
end

function _node_ordinal(names::Tuple, sought::Symbol, ordinal::Int=1)
    isempty(names) && throw(AlgorithmGraphError(
        "algorithm graph has no node named $sought",
    ))
    first(names) === sought && return ordinal
    return _node_ordinal(Base.tail(names), sought, ordinal + 1)
end

@inline _validate_direct_links!(
    nodes,
    node_outputs,
    ::Tuple{},
    target,
) = nothing
@inline function _validate_direct_links!(
    nodes,
    node_outputs,
    links::Tuple,
    target,
)
    direct = first(links)
    _, source_contract = _endpoint_contract(nodes, direct.source)
    _, destination_contract =
        _endpoint_contract(nodes, direct.destination)
    _require_data_output(source_contract, direct.source)
    _require_data_input(destination_contract, direct.destination)
    context = "direct link $(_node_name(direct.source)) => " *
              "$(_node_name(direct.destination))"
    _require_matching_formats(
        source_contract.format,
        destination_contract.format,
        context,
    )
    names = keys(nodes)
    source_ordinal = _node_ordinal(names, _node_name(direct.source))
    destination_ordinal = _node_ordinal(names, _node_name(direct.destination))
    source_ordinal < destination_ordinal || throw(AlgorithmGraphError(
        "$context must point from an earlier node to a later node; use delayed_link for feedback",
    ))
    source_values = _endpoint_output(node_outputs, direct.source)
    _validate_graph_buffer(
        target,
        destination_contract.format,
        source_values,
    )
    _validate_direct_links!(
        nodes,
        node_outputs,
        Base.tail(links),
        target,
    )
    return nothing
end

struct _PreparedDelayedLink{Destination,Source,Initial,Format}
    destination::Destination
    source::Source
    initial::Initial
    format::Format
end

function _prepare_delayed_link(
    definition::DelayedAlgorithmLink,
    nodes,
    node_outputs,
    target,
)
    _, source_contract = _endpoint_contract(nodes, definition.source)
    _, destination_contract =
        _endpoint_contract(nodes, definition.destination)
    _require_data_output(source_contract, definition.source)
    _require_data_input(destination_contract, definition.destination)
    context = "delayed link $(_node_name(definition.source)) => " *
              "$(_node_name(definition.destination))"
    _require_matching_formats(
        source_contract.format,
        destination_contract.format,
        context,
    )
    _validate_graph_buffer(
        target,
        destination_contract.format,
        definition.initial,
    )
    source = _endpoint_output(node_outputs, definition.source)
    _validate_graph_buffer(
        target,
        destination_contract.format,
        source,
    )
    initial = _allocate_graph_buffer(
        target,
        destination_contract.format,
    )
    state = _allocate_graph_buffer(
        target,
        destination_contract.format,
    )
    _copy_graph_buffer!(
        initial,
        definition.initial,
    )
    _copy_graph_buffer!(
        state,
        definition.initial,
    )
    prepared = _PreparedDelayedLink(
        definition.destination,
        source,
        initial,
        destination_contract.format,
    )
    return prepared, state
end

@inline _prepare_delayed_links(
    ::Tuple{},
    nodes,
    node_outputs,
    target,
) = ((), ())
@inline function _prepare_delayed_links(
    definitions::Tuple,
    nodes,
    node_outputs,
    target,
)
    prepared, state = _prepare_delayed_link(
        first(definitions),
        nodes,
        node_outputs,
        target,
    )
    tail_prepared, tail_state = _prepare_delayed_links(
        Base.tail(definitions),
        nodes,
        node_outputs,
        target,
    )
    return (prepared, tail_prepared...), (state, tail_state...)
end

struct _PreparedInputBinding{Destination,Values}
    destination::Destination
    values::Values
end

@inline _input_bindings(::Tuple{}) = ()
@inline function _input_bindings(inputs::Tuple)
    input = first(inputs)
    return (
        _PreparedInputBinding(input.destination, input.values),
        _input_bindings(Base.tail(inputs))...,
    )
end

@inline _direct_bindings(::Tuple{}, node_outputs) = ()
@inline function _direct_bindings(links::Tuple, node_outputs)
    direct = first(links)
    return (
        _PreparedInputBinding(
            direct.destination,
            _endpoint_output(node_outputs, direct.source),
        ),
        _direct_bindings(Base.tail(links), node_outputs)...,
    )
end

@inline _delayed_bindings(::Tuple{}, ::Tuple{}) = ()
@inline function _delayed_bindings(links::Tuple, state::Tuple)
    delayed = first(links)
    return (
        _PreparedInputBinding(delayed.destination, first(state)),
        _delayed_bindings(Base.tail(links), Base.tail(state))...,
    )
end

@inline _binding_key(binding::_PreparedInputBinding) =
    (_node_name(binding.destination), _port_name(binding.destination))

@inline _binding_keys(::Tuple{}) = ()
@inline _binding_keys(bindings::Tuple) =
    (_binding_key(first(bindings)), _binding_keys(Base.tail(bindings))...)

function _require_unique_destinations(bindings::Tuple)
    keys = _binding_keys(bindings)
    length(keys) == length(unique(keys)) || throw(AlgorithmGraphError(
        "each frame-data input must have exactly one graph input, direct link, or delayed link",
    ))
    return nothing
end

@inline function _binding_matches(
    binding::_PreparedInputBinding,
    ::AlgorithmPortReference{Node,Port},
) where {Node,Port}
    return _node_name(binding.destination) === Node &&
           _port_name(binding.destination) === Port
end

function _require_input_binding(::Tuple{}, reference::AlgorithmPortReference)
    throw(AlgorithmGraphError(
        "frame-data input $(_node_name(reference)) => $(_port_name(reference)) is unconnected",
    ))
end

function _require_input_binding(bindings::Tuple, reference::AlgorithmPortReference)
    binding = first(bindings)
    _binding_matches(binding, reference) && return binding.values
    return _require_input_binding(Base.tail(bindings), reference)
end

@inline _bind_node_inputs(node_name::Symbol, ::Tuple{}, bindings) = ()
@inline function _bind_node_inputs(node_name::Symbol, contracts::Tuple, bindings)
    contract = first(contracts)
    reference = AlgorithmPortReference(node_name, _port_name(contract))
    return (
        _require_input_binding(bindings, reference),
        _bind_node_inputs(node_name, Base.tail(contracts), bindings)...,
    )
end

@inline function _parameter_binding_matches(
    binding::SparseParameterBinding,
    ::AlgorithmPortReference{Node,Port},
) where {Node,Port}
    return _node_name(binding.endpoint) === Node &&
           _port_name(binding.endpoint) === Port
end

function _require_parameter_binding(::Tuple{}, reference::AlgorithmPortReference)
    throw(AlgorithmGraphError(
        "sparse parameter $(_node_name(reference)) => $(_port_name(reference)) is unbound",
    ))
end

function _require_parameter_binding(
    bindings::Tuple,
    reference::AlgorithmPortReference,
)
    binding = first(bindings)
    _parameter_binding_matches(binding, reference) && return binding.values
    return _require_parameter_binding(Base.tail(bindings), reference)
end

@inline _bind_node_parameters(node_name::Symbol, ::Tuple{}, bindings) = ()
@inline function _bind_node_parameters(node_name::Symbol, contracts::Tuple, bindings)
    contract = first(contracts)
    reference = AlgorithmPortReference(node_name, _port_name(contract))
    return (
        _require_parameter_binding(bindings, reference),
        _bind_node_parameters(node_name, Base.tail(contracts), bindings)...,
    )
end

"""One prepared graph-node owner and its exact named buffer bindings."""
struct PreparedGraphNode{
    Name,
    Owner,
    Inputs<:NamedTuple,
    Outputs<:NamedTuple,
    Parameters<:NamedTuple,
}
    owner::Owner
    inputs::Inputs
    outputs::Outputs
    parameters::Parameters
end

@inline _node_name(::PreparedGraphNode{Name}) where {Name} = Name

function _prepare_exact_node(
    node::_AdmittedGraphNode{Name,Node},
    outputs::NamedTuple,
    bindings::Tuple,
    parameter_bindings::Tuple,
    target::AbstractComputeDevice,
) where {Name,Node}
    input_contracts = _data_input_contracts(values(node.ports))
    input_values = _bind_node_inputs(Name, input_contracts, bindings)
    inputs = _named_tuple(
        _port_names(input_contracts),
        input_values,
        "algorithm input port",
    )
    parameter_contracts = _parameter_contracts(values(node.ports))
    parameter_values = _bind_node_parameters(
        Name,
        parameter_contracts,
        parameter_bindings,
    )
    parameters = _named_tuple(
        _port_names(parameter_contracts),
        parameter_values,
        "sparse parameter port",
    )
    owner = prepare_graph_node(
        Node,
        node.config,
        node.props,
        inputs,
        outputs,
        parameters,
        target,
    )
    return PreparedGraphNode{
        Name,
        typeof(owner),
        typeof(inputs),
        typeof(outputs),
        typeof(parameters),
    }(
        owner,
        inputs,
        outputs,
        parameters,
    )
end

@inline _prepare_exact_nodes(
    ::Tuple{},
    ::Tuple{},
    bindings,
    parameter_bindings,
    target,
) = ()
@inline function _prepare_exact_nodes(
    nodes::Tuple,
    outputs::Tuple,
    bindings,
    parameter_bindings,
    target,
)
    return (
        _prepare_exact_node(
            first(nodes),
            first(outputs),
            bindings,
            parameter_bindings,
            target,
        ),
        _prepare_exact_nodes(
            Base.tail(nodes),
            Base.tail(outputs),
            bindings,
            parameter_bindings,
            target,
        )...,
    )
end

function _prepare_graph_outputs(
    definitions::Tuple,
    nodes::NamedTuple,
    node_outputs::NamedTuple,
)
    values = _prepare_graph_output_values(definitions, nodes, node_outputs)
    return _named_tuple(_boundary_names(definitions), values, "graph output")
end

@inline _prepare_graph_output_values(::Tuple{}, nodes, node_outputs) = ()
@inline function _prepare_graph_output_values(definitions::Tuple, nodes, node_outputs)
    definition = first(definitions)
    _, contract = _endpoint_contract(nodes, definition.source)
    _require_data_output(contract, definition.source)
    return (
        _endpoint_output(node_outputs, definition.source),
        _prepare_graph_output_values(
            Base.tail(definitions),
            nodes,
            node_outputs,
        )...,
    )
end

"""Mutable single-writer state for one exact prepared algorithm graph."""
mutable struct AlgorithmGraphState{DelayedValues}
    delayed_values::DelayedValues
    step_sequence::UInt64
    pending_sequence::UInt64
    pending::Bool
    failed::Bool
end

"""Concrete prepared owner for one portable algorithm graph run."""
struct PreparedAlgorithmGraph{
    Target<:AbstractComputeDevice,
    Context,
    Nodes<:NamedTuple,
    Execution,
    DelayedLinks<:Tuple,
    Inputs<:NamedTuple,
    Outputs<:NamedTuple,
    State<:AlgorithmGraphState,
}
    name::Symbol
    target::Target
    context::Context
    nodes::Nodes
    execution::Execution
    delayed_links::DelayedLinks
    inputs::Inputs
    outputs::Outputs
    state::State
end

"""
    prepare_algorithm_graph(definition; target=HostComputeDevice(),
                            execution=StreamGraphExecution())

Admit every node's fixed ports, validate all exact formats and connections,
allocate intermediate and delayed storage, then prepare every node against its
exact frame-data buffers and startup sparse parameters in one concrete
single-writer graph owner. Every graph buffer must be native packed storage on
the exact `target`; preparation never inserts an implicit host/device transfer.
`CapturedGraphExecution()` records the complete node and delayed-link sequence
as one native device graph. Every node owner must explicitly satisfy the
device-command-graph capture contract.
"""
function _prepare_algorithm_graph(
    definition::AlgorithmGraphDefinition,
    target::AbstractComputeDevice,
    context,
    execution,
)
    admitted_nodes = _admitted_graph_nodes(definition)
    _validate_parameters!(admitted_nodes, definition.parameters, target)

    node_outputs = _prepared_node_outputs(admitted_nodes, target)
    prepared_inputs = _prepared_graph_inputs(definition, admitted_nodes, target)
    _validate_direct_links!(
        admitted_nodes,
        node_outputs,
        definition.links,
        target,
    )
    prepared_delays, delayed_values = _prepare_delayed_links(
        definition.delayed_links,
        admitted_nodes,
        node_outputs,
        target,
    )

    bindings = (
        _input_bindings(values(prepared_inputs))...,
        _direct_bindings(definition.links, node_outputs)...,
        _delayed_bindings(prepared_delays, delayed_values)...,
    )
    _require_unique_destinations(bindings)

    node_values = _prepare_exact_nodes(
        values(admitted_nodes),
        values(node_outputs),
        bindings,
        definition.parameters,
        target,
    )
    nodes = NamedTuple{keys(admitted_nodes)}(node_values)
    outputs = _prepare_graph_outputs(
        definition.outputs,
        admitted_nodes,
        node_outputs,
    )
    state = AlgorithmGraphState(
        delayed_values,
        UInt64(0),
        UInt64(0),
        false,
        false,
    )
    prepared_execution = _prepare_graph_execution(
        execution,
        values(nodes),
        prepared_delays,
        state.delayed_values,
        context,
    )
    return PreparedAlgorithmGraph(
        definition.name,
        target,
        context,
        nodes,
        prepared_execution,
        prepared_delays,
        prepared_inputs,
        outputs,
        state,
    )
end

function prepare_algorithm_graph(
    definition::AlgorithmGraphDefinition;
    target::AbstractComputeDevice=HostComputeDevice(),
    execution::Union{StreamGraphExecution,CapturedGraphExecution}=
        StreamGraphExecution(),
)
    context = _prepare_device_execution_context(target)
    return _with_prepared_device_execution_context(context) do
        try
            return _prepare_algorithm_graph(
                definition,
                target,
                context,
                execution,
            )
        finally
            _synchronize_prepared_device_execution_context!(context)
        end
    end
end

@inline compute_device(graph::PreparedAlgorithmGraph) = graph.target
@inline graph_name(graph::PreparedAlgorithmGraph) = graph.name

@inline graph_input(
    graph::PreparedAlgorithmGraph,
    ::Val{Name},
) where {Name} = getproperty(graph.inputs, Name).values

function graph_input(graph::PreparedAlgorithmGraph, name::Symbol)
    hasproperty(graph.inputs, name) || throw(AlgorithmGraphError(
        "prepared graph has no input named $name",
    ))
    return getproperty(graph.inputs, name).values
end

@inline graph_output(
    graph::PreparedAlgorithmGraph,
    ::Val{Name},
) where {Name} = getproperty(graph.outputs, Name)

function graph_output(graph::PreparedAlgorithmGraph, name::Symbol)
    hasproperty(graph.outputs, name) || throw(AlgorithmGraphError(
        "prepared graph has no output named $name",
    ))
    return getproperty(graph.outputs, name)
end

@inline prepared_graph_node(
    graph::PreparedAlgorithmGraph,
    ::Val{Name},
) where {Name} = getproperty(graph.nodes, Name).owner

function prepared_graph_node(graph::PreparedAlgorithmGraph, name::Symbol)
    hasproperty(graph.nodes, name) || throw(AlgorithmGraphError(
        "prepared graph has no node named $name",
    ))
    return getproperty(graph.nodes, name).owner
end

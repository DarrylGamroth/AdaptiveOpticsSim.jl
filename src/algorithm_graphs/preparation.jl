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

function _graph_port_contract(
    name::Symbol,
    direction::Symbol,
    role::Symbol,
    ::Type{T},
    shape::NTuple{N,<:Integer},
    schema::AbstractString,
    layout::Symbol,
) where {T,N}
    _require_graph_name(name, "algorithm port")
    direction in (:input, :output) || throw(AlgorithmGraphError(
        "algorithm port $name has unsupported direction $direction",
    ))
    role in (:data, :parameter) || throw(AlgorithmGraphError(
        "algorithm port $name has unsupported role $role",
    ))
    direction === :output && role === :parameter &&
        throw(AlgorithmGraphError("sparse parameter $name must be an input"))
    all(>(0), shape) || throw(AlgorithmGraphError(
        "algorithm port $name dimensions must be positive",
    ))
    converted = ntuple(index -> Int(shape[index]), Val(N))
    format = _GraphPortFormat{T,N}(converted, String(schema), layout)
    return _GraphPortContract{name,direction,role,typeof(format)}(format)
end

@inline _validate_port_contract(::_GraphPortContract) = nothing
@inline _validate_port_contract(value) = throw(AlgorithmGraphError(
    "algorithm adapters must return graph port contracts, not $(typeof(value))",
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

function _prepare_algorithm_instance(::Type{Declaration}, configuration) where {Declaration}
    throw(AlgorithmGraphError(
        "no AlgorithmGraphs adapter is loaded for declaration $Declaration",
    ))
end

function _algorithm_port_contracts(::Type{Declaration}, prepared) where {Declaration}
    throw(AlgorithmGraphError(
        "the AlgorithmGraphs adapter for $Declaration did not provide port contracts",
    ))
end

function _replace_algorithm_parameter!(prepared, name::Val, values)
    throw(AlgorithmGraphError(
        "the AlgorithmGraphs adapter for $(typeof(prepared)) does not support " *
        "sparse parameter $(name)",
    ))
end

function _process_algorithm!(prepared, outputs::NamedTuple, inputs::NamedTuple)
    throw(AlgorithmGraphError(
        "the AlgorithmGraphs adapter for $(typeof(prepared)) does not implement processing",
    ))
end

function _reset_algorithm!(prepared)
    throw(AlgorithmGraphError(
        "the AlgorithmGraphs adapter for $(typeof(prepared)) does not implement reset",
    ))
end

function _allocate_algorithm_buffer(prepared, format::_GraphPortFormat{T,N}) where {T,N}
    format.layout === :column_major || throw(AlgorithmGraphError(
        "the native graph allocator supports only column-major buffers; got " *
        "$(format.layout)",
    ))
    return zeros(T, format.shape)
end

function _validate_algorithm_buffer(
    prepared,
    format::_GraphPortFormat{T,N},
    values::Array{T,N},
) where {T,N}
    size(values) == format.shape || throw(AlgorithmGraphError(
        "buffer shape $(size(values)) does not match declared shape $(format.shape)",
    ))
    format.layout === :column_major || throw(AlgorithmGraphError(
        "ordinary Julia arrays cannot satisfy $(format.layout) packed layout",
    ))
    return values
end

function _validate_algorithm_buffer(prepared, format::_GraphPortFormat{T,N}, values) where {T,N}
    throw(AlgorithmGraphError(
        "buffer $(typeof(values)) does not satisfy packed $T rank-$N storage",
    ))
end

@inline function _copy_algorithm_buffer!(prepared, destination, source)
    copyto!(destination, source)
    return destination
end

struct _PreparedNodePrototype{Name,Declaration,Algorithm,Ports<:NamedTuple}
    algorithm::Algorithm
    ports::Ports
end

@inline _node_name(::_PreparedNodePrototype{Name}) where {Name} = Name

function _prepare_node_prototype(
    definition::AlgorithmNodeDefinition{Name,Declaration},
) where {Name,Declaration}
    algorithm = _prepare_algorithm_instance(Declaration, definition.configuration)
    port_values = _algorithm_port_contracts(Declaration, algorithm)
    _validate_tuple(port_values, _validate_port_contract)
    ports = _named_tuple(_port_names(port_values), port_values, "algorithm port")
    return _PreparedNodePrototype{Name,Declaration,typeof(algorithm),typeof(ports)}(
        algorithm,
        ports,
    )
end

@inline _prepare_node_prototypes(::Tuple{}) = ()
@inline function _prepare_node_prototypes(definitions::Tuple)
    return (
        _prepare_node_prototype(first(definitions)),
        _prepare_node_prototypes(Base.tail(definitions))...,
    )
end

function _prepared_node_prototypes(definition::AlgorithmGraphDefinition)
    values = _prepare_node_prototypes(definition.nodes)
    return _named_tuple(_node_names(definition.nodes), values, "algorithm node")
end

function _require_node(prototypes::NamedTuple, ::Val{Name}) where {Name}
    hasproperty(prototypes, Name) || throw(AlgorithmGraphError(
        "algorithm graph has no node named $Name",
    ))
    return getproperty(prototypes, Name)
end

function _require_port(prototype::_PreparedNodePrototype, ::Val{Name}) where {Name}
    hasproperty(prototype.ports, Name) || throw(AlgorithmGraphError(
        "algorithm node $(_node_name(prototype)) has no port named $Name",
    ))
    return getproperty(prototype.ports, Name)
end

function _endpoint_contract(
    prototypes::NamedTuple,
    reference::AlgorithmPortReference{Node,Port},
) where {Node,Port}
    prototype = _require_node(prototypes, Val(Node))
    return prototype, _require_port(prototype, Val(Port))
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

function _validate_parameters!(prototypes::NamedTuple, parameters::Tuple)
    keys = _parameter_endpoint_keys(parameters)
    length(keys) == length(unique(keys)) || throw(AlgorithmGraphError(
        "each sparse parameter port may be initialized at most once",
    ))
    _validate_parameters_recursive!(prototypes, parameters)
    return nothing
end

@inline _validate_parameters_recursive!(prototypes, ::Tuple{}) = nothing
@inline function _validate_parameters_recursive!(prototypes, parameters::Tuple)
    binding = first(parameters)
    prototype, contract = _endpoint_contract(prototypes, binding.endpoint)
    _require_sparse_parameter(contract, binding.endpoint)
    _validate_algorithm_buffer(
        prototype.algorithm,
        contract.format,
        binding.values,
    )
    _validate_parameters_recursive!(prototypes, Base.tail(parameters))
    return nothing
end

@inline _apply_parameters!(prototypes, ::Tuple{}) = nothing
@inline function _apply_parameters!(prototypes, parameters::Tuple)
    binding = first(parameters)
    prototype = _require_node(
        prototypes,
        Val(_node_name(binding.endpoint)),
    )
    _replace_algorithm_parameter!(
        prototype.algorithm,
        Val(_port_name(binding.endpoint)),
        binding.values,
    )
    _apply_parameters!(prototypes, Base.tail(parameters))
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

@inline _allocate_output_buffers(algorithm, ::Tuple{}) = ()
@inline function _allocate_output_buffers(algorithm, contracts::Tuple)
    return (
        _allocate_algorithm_buffer(algorithm, first(contracts).format),
        _allocate_output_buffers(algorithm, Base.tail(contracts))...,
    )
end

function _allocate_node_outputs(prototype::_PreparedNodePrototype)
    contracts = _data_output_contracts(values(prototype.ports))
    buffers = _allocate_output_buffers(prototype.algorithm, contracts)
    return _named_tuple(_port_names(contracts), buffers, "algorithm output port")
end

@inline _allocate_all_node_outputs(::Tuple{}) = ()
@inline function _allocate_all_node_outputs(prototypes::Tuple)
    return (
        _allocate_node_outputs(first(prototypes)),
        _allocate_all_node_outputs(Base.tail(prototypes))...,
    )
end

function _prepared_node_outputs(prototypes::NamedTuple)
    groups = _allocate_all_node_outputs(values(prototypes))
    return NamedTuple{keys(prototypes)}(groups)
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
    prototypes::NamedTuple,
) where {Name}
    prototype, contract = _endpoint_contract(prototypes, definition.destination)
    _require_data_input(contract, definition.destination)
    _validate_algorithm_buffer(
        prototype.algorithm,
        contract.format,
        definition.values,
    )
    return _PreparedGraphInput{
        Name,
        typeof(definition.destination),
        typeof(definition.values),
    }(definition.destination, definition.values)
end

@inline _prepare_graph_inputs(::Tuple{}, prototypes) = ()
@inline function _prepare_graph_inputs(definitions::Tuple, prototypes)
    return (
        _prepare_graph_input(first(definitions), prototypes),
        _prepare_graph_inputs(Base.tail(definitions), prototypes)...,
    )
end

function _prepared_graph_inputs(definition::AlgorithmGraphDefinition, prototypes)
    prepared = _prepare_graph_inputs(definition.inputs, prototypes)
    return _named_tuple(_boundary_names(definition.inputs), prepared, "graph input")
end

function _node_ordinal(names::Tuple, sought::Symbol, ordinal::Int=1)
    isempty(names) && throw(AlgorithmGraphError(
        "algorithm graph has no node named $sought",
    ))
    first(names) === sought && return ordinal
    return _node_ordinal(Base.tail(names), sought, ordinal + 1)
end

@inline _validate_direct_links!(prototypes, node_outputs, ::Tuple{}) = nothing
@inline function _validate_direct_links!(
    prototypes,
    node_outputs,
    links::Tuple,
)
    direct = first(links)
    _, source_contract = _endpoint_contract(prototypes, direct.source)
    destination_prototype, destination_contract =
        _endpoint_contract(prototypes, direct.destination)
    _require_data_output(source_contract, direct.source)
    _require_data_input(destination_contract, direct.destination)
    context = "direct link $(_node_name(direct.source)) => " *
              "$(_node_name(direct.destination))"
    _require_matching_formats(
        source_contract.format,
        destination_contract.format,
        context,
    )
    names = keys(prototypes)
    source_ordinal = _node_ordinal(names, _node_name(direct.source))
    destination_ordinal = _node_ordinal(names, _node_name(direct.destination))
    source_ordinal < destination_ordinal || throw(AlgorithmGraphError(
        "$context must point from an earlier node to a later node; use delayed_link for feedback",
    ))
    source_values = _endpoint_output(node_outputs, direct.source)
    _validate_algorithm_buffer(
        destination_prototype.algorithm,
        destination_contract.format,
        source_values,
    )
    _validate_direct_links!(
        prototypes,
        node_outputs,
        Base.tail(links),
    )
    return nothing
end

struct _PreparedDelayedLink{Destination,Algorithm,Source,Initial,Format}
    destination::Destination
    destination_algorithm::Algorithm
    source::Source
    initial::Initial
    format::Format
end

function _prepare_delayed_link(
    definition::DelayedAlgorithmLink,
    prototypes,
    node_outputs,
)
    _, source_contract = _endpoint_contract(prototypes, definition.source)
    destination_prototype, destination_contract =
        _endpoint_contract(prototypes, definition.destination)
    _require_data_output(source_contract, definition.source)
    _require_data_input(destination_contract, definition.destination)
    context = "delayed link $(_node_name(definition.source)) => " *
              "$(_node_name(definition.destination))"
    _require_matching_formats(
        source_contract.format,
        destination_contract.format,
        context,
    )
    _validate_algorithm_buffer(
        destination_prototype.algorithm,
        destination_contract.format,
        definition.initial,
    )
    source = _endpoint_output(node_outputs, definition.source)
    _validate_algorithm_buffer(
        destination_prototype.algorithm,
        destination_contract.format,
        source,
    )
    initial = _allocate_algorithm_buffer(
        destination_prototype.algorithm,
        destination_contract.format,
    )
    state = _allocate_algorithm_buffer(
        destination_prototype.algorithm,
        destination_contract.format,
    )
    _copy_algorithm_buffer!(
        destination_prototype.algorithm,
        initial,
        definition.initial,
    )
    _copy_algorithm_buffer!(
        destination_prototype.algorithm,
        state,
        definition.initial,
    )
    prepared = _PreparedDelayedLink(
        definition.destination,
        destination_prototype.algorithm,
        source,
        initial,
        destination_contract.format,
    )
    return prepared, state
end

@inline _prepare_delayed_links(::Tuple{}, prototypes, node_outputs) = ((), ())
@inline function _prepare_delayed_links(definitions::Tuple, prototypes, node_outputs)
    prepared, state = _prepare_delayed_link(
        first(definitions),
        prototypes,
        node_outputs,
    )
    tail_prepared, tail_state = _prepare_delayed_links(
        Base.tail(definitions),
        prototypes,
        node_outputs,
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

"""One prepared algorithm and its exact named graph buffer bindings."""
struct PreparedAlgorithmNode{Name,Algorithm,Inputs<:NamedTuple,Outputs<:NamedTuple}
    algorithm::Algorithm
    inputs::Inputs
    outputs::Outputs
end

@inline _node_name(::PreparedAlgorithmNode{Name}) where {Name} = Name

function _bind_prepared_node(
    prototype::_PreparedNodePrototype{Name},
    outputs::NamedTuple,
    bindings::Tuple,
) where {Name}
    contracts = _data_input_contracts(values(prototype.ports))
    input_values = _bind_node_inputs(Name, contracts, bindings)
    inputs = _named_tuple(_port_names(contracts), input_values, "algorithm input port")
    return PreparedAlgorithmNode{Name,typeof(prototype.algorithm),typeof(inputs),typeof(outputs)}(
        prototype.algorithm,
        inputs,
        outputs,
    )
end

@inline _bind_prepared_nodes(::Tuple{}, ::Tuple{}, bindings) = ()
@inline function _bind_prepared_nodes(prototypes::Tuple, outputs::Tuple, bindings)
    return (
        _bind_prepared_node(first(prototypes), first(outputs), bindings),
        _bind_prepared_nodes(
            Base.tail(prototypes),
            Base.tail(outputs),
            bindings,
        )...,
    )
end

function _prepare_graph_outputs(
    definitions::Tuple,
    prototypes::NamedTuple,
    node_outputs::NamedTuple,
)
    values = _prepare_graph_output_values(definitions, prototypes, node_outputs)
    return _named_tuple(_boundary_names(definitions), values, "graph output")
end

@inline _prepare_graph_output_values(::Tuple{}, prototypes, node_outputs) = ()
@inline function _prepare_graph_output_values(definitions::Tuple, prototypes, node_outputs)
    definition = first(definitions)
    _, contract = _endpoint_contract(prototypes, definition.source)
    _require_data_output(contract, definition.source)
    return (
        _endpoint_output(node_outputs, definition.source),
        _prepare_graph_output_values(
            Base.tail(definitions),
            prototypes,
            node_outputs,
        )...,
    )
end

"""Mutable single-writer state for one exact prepared algorithm graph."""
mutable struct AlgorithmGraphState{DelayedValues}
    delayed_values::DelayedValues
    step_sequence::UInt64
    failed::Bool
end

"""Concrete prepared owner for one portable algorithm graph run."""
struct PreparedAlgorithmGraph{
    Nodes<:NamedTuple,
    DelayedLinks<:Tuple,
    Inputs<:NamedTuple,
    Outputs<:NamedTuple,
    State<:AlgorithmGraphState,
}
    nodes::Nodes
    delayed_links::DelayedLinks
    inputs::Inputs
    outputs::Outputs
    state::State
end

"""
    prepare_algorithm_graph(definition)

Prepare every algorithm, validate all exact formats and connections, install
startup sparse parameters, allocate intermediate and delayed storage, and bind
one concrete single-writer graph owner.
"""
function prepare_algorithm_graph(definition::AlgorithmGraphDefinition)
    prototypes = _prepared_node_prototypes(definition)
    _validate_parameters!(prototypes, definition.parameters)
    _apply_parameters!(prototypes, definition.parameters)

    node_outputs = _prepared_node_outputs(prototypes)
    prepared_inputs = _prepared_graph_inputs(definition, prototypes)
    _validate_direct_links!(
        prototypes,
        node_outputs,
        definition.links,
    )
    prepared_delays, delayed_values = _prepare_delayed_links(
        definition.delayed_links,
        prototypes,
        node_outputs,
    )

    bindings = (
        _input_bindings(values(prepared_inputs))...,
        _direct_bindings(definition.links, node_outputs)...,
        _delayed_bindings(prepared_delays, delayed_values)...,
    )
    _require_unique_destinations(bindings)

    node_values = _bind_prepared_nodes(
        values(prototypes),
        values(node_outputs),
        bindings,
    )
    nodes = NamedTuple{keys(prototypes)}(node_values)
    outputs = _prepare_graph_outputs(
        definition.outputs,
        prototypes,
        node_outputs,
    )
    state = AlgorithmGraphState(delayed_values, UInt64(0), false)
    return PreparedAlgorithmGraph(
        nodes,
        prepared_delays,
        prepared_inputs,
        outputs,
        state,
    )
end

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

@inline prepared_algorithm(
    graph::PreparedAlgorithmGraph,
    ::Val{Name},
) where {Name} = getproperty(graph.nodes, Name).algorithm

function prepared_algorithm(graph::PreparedAlgorithmGraph, name::Symbol)
    hasproperty(graph.nodes, name) || throw(AlgorithmGraphError(
        "prepared graph has no algorithm node named $name",
    ))
    return getproperty(graph.nodes, name).algorithm
end

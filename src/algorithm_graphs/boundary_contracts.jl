"""
    graph_input_contract(definition, name::Symbol)

Describe a declared graph input without allocating graph arrays or preparing
numerical execution. Return an immutable named tuple with the boundary `name`,
port `direction` and `role`, `eltype`, `shape`, scientific `schema` and `layout`.
The contract comes from the destination node adapter, not from the supplied
input buffer. Unknown or ambiguous boundaries and invalid endpoints throw
`AlgorithmGraphError`.
"""
function graph_input_contract(definition::AlgorithmGraphDefinition, name::Symbol)
    boundary = _require_boundary(definition.inputs, name)
    _, port = _endpoint_contract(_admitted_graph_nodes(definition), boundary.destination)
    _require_data_input(port, boundary.destination)
    return _boundary_contract(port, name)
end

"""
    graph_output_contract(definition, name::Symbol)

Describe a declared graph output using the source node adapter's contract.
Return the same immutable fields as [`graph_input_contract`](@ref), without
preparing or stepping numerical owners. Unknown or ambiguous boundaries and
invalid endpoints throw `AlgorithmGraphError`.
"""
function graph_output_contract(definition::AlgorithmGraphDefinition, name::Symbol)
    boundary = _require_boundary(definition.outputs, name)
    _, port = _endpoint_contract(_admitted_graph_nodes(definition), boundary.source)
    _require_data_output(port, boundary.source)
    return _boundary_contract(port, name)
end

function _require_boundary(boundaries::Tuple, name::Symbol)
    matches = filter(boundary -> _boundary_name(boundary) === name, boundaries)
    length(matches) == 1 || throw(AlgorithmGraphError(
        "algorithm graph requires exactly one boundary named $name",
    ))
    return only(matches)
end

function _boundary_contract(
        port::_GraphPortContract{Name, Direction, Role, _GraphPortFormat{T, N}}, name::Symbol
) where {Name, Direction, Role, T, N}
    return (;
        name, direction = Direction, role = Role, eltype = T,
        shape = port.format.shape, schema = port.format.schema, layout = port.format.layout
    )
end

const ALGORITHM_GRAPH_TOML_SCHEMA_VERSION = 1

@inline function _file_error(context::AbstractString, message::AbstractString)
    return AlgorithmGraphError("$context: $message")
end

function _require_table_keys(
    table::AbstractDict,
    allowed::Tuple,
    required::Tuple,
    context::AbstractString,
)
    for key in keys(table)
        key in allowed || throw(_file_error(context, "unknown field '$key'"))
    end
    for key in required
        haskey(table, key) || throw(_file_error(context, "missing field '$key'"))
    end
    return nothing
end

function _require_named_fields(
    values::NamedTuple,
    allowed::Tuple,
    required::Tuple,
    context::AbstractString,
)
    for key in keys(values)
        key in allowed || throw(_file_error(context, "unknown field '$key'"))
    end
    for key in required
        hasproperty(values, key) || throw(_file_error(
            context,
            "missing field '$key'",
        ))
    end
    return nothing
end

function _file_string(value, context::AbstractString)
    value isa AbstractString || throw(_file_error(
        context,
        "expected a string, got $(typeof(value))",
    ))
    return String(value)
end

function _file_integer(value, context::AbstractString)
    value isa Integer && !(value isa Bool) || throw(_file_error(
        context,
        "expected an integer, got $(typeof(value))",
    ))
    return Int(value)
end

function _file_real(value, context::AbstractString)
    value isa Real && !(value isa Bool) || throw(_file_error(
        context,
        "expected a real number, got $(typeof(value))",
    ))
    isfinite(value) || throw(_file_error(context, "value must be finite"))
    return value
end

function _file_identifier(value, context::AbstractString)
    text = _file_string(value, context)
    Base.isidentifier(text) || throw(_file_error(
        context,
        "'$text' is not a valid identifier",
    ))
    return Symbol(text)
end

function _file_endpoint(value, context::AbstractString)
    text = _file_string(value, context)
    parts = split(text, '.'; keepempty=true)
    length(parts) == 2 || throw(_file_error(
        context,
        "endpoint '$text' must use node.port syntax",
    ))
    node = _file_identifier(parts[1], "$context node")
    port = _file_identifier(parts[2], "$context port")
    return AlgorithmPortReference(node, port)
end

_normalize_toml_value(value::Union{Bool,Integer,AbstractFloat,AbstractString}) =
    value

function _normalize_toml_value(values::AbstractVector)
    builder = Vector{Any}(undef, length(values))
    for index in eachindex(values)
        builder[index] = _normalize_toml_value(values[index])
    end
    return Tuple(builder)
end


function _normalize_toml_value(table::AbstractDict)
    names = sort!(collect(String, keys(table)))
    symbols = Tuple(Symbol(name) for name in names)
    values = Tuple(_normalize_toml_value(table[name]) for name in names)
    return NamedTuple{symbols}(values)
end

function _normalize_toml_value(value)
    throw(AlgorithmGraphError(
        "unsupported TOML value type $(typeof(value)) in graph config or props",
    ))
end

function _table_or_empty(table::AbstractDict, key::String, context::AbstractString)
    value = get(table, key, Dict{String,Any}())
    value isa AbstractDict || throw(_file_error(
        context,
        "field '$key' must be a table",
    ))
    return value
end

function _table_array_or_empty(
    table::AbstractDict,
    key::String,
    context::AbstractString,
)
    value = get(table, key, Any[])
    value isa AbstractVector || throw(_file_error(
        context,
        "field '$key' must be an array of tables",
    ))
    for (index, entry) in pairs(value)
        entry isa AbstractDict || throw(_file_error(
            context,
            "field '$key' entry $index must be a table",
        ))
    end
    return value
end

function _require_binding(
    bindings::NamedTuple,
    value,
    context::AbstractString,
)
    name = _file_identifier(value, context)
    hasproperty(bindings, name) || throw(_file_error(
        context,
        "no external binding named '$name' was supplied",
    ))
    return getproperty(bindings, name)
end

function _discrete_integrator_f32_node(
    name::Symbol,
    config::NamedTuple,
    props::NamedTuple,
)
    config_fields = (
        :extent,
        :input_schema,
        :output_schema,
        :sample_period_s,
    )
    prop_fields = (:gain, :tau_s)
    _require_named_fields(
        config,
        config_fields,
        config_fields,
        "node '$name' config",
    )
    _require_named_fields(
        props,
        prop_fields,
        prop_fields,
        "node '$name' props",
    )
    return discrete_integrator_node(
        name;
        extent=_file_integer(config.extent, "node '$name' config.extent"),
        sample_period_s=_file_real(
            config.sample_period_s,
            "node '$name' config.sample_period_s",
        ),
        input_schema=_file_string(
            config.input_schema,
            "node '$name' config.input_schema",
        ),
        output_schema=_file_string(
            config.output_schema,
            "node '$name' config.output_schema",
        ),
        gain=_file_real(props.gain, "node '$name' props.gain"),
        tau_s=_file_real(props.tau_s, "node '$name' props.tau_s"),
        T=Float32,
    )
end

function _modal_opd_expansion_f32_node(
    name::Symbol,
    config::NamedTuple,
    props::NamedTuple,
)
    config_fields = (
        :basis_schema,
        :coefficients_schema,
        :mode_count,
        :opd_schema,
        :pupil_columns,
        :pupil_rows,
        :pupil_support_schema,
    )
    _require_named_fields(
        config,
        config_fields,
        config_fields,
        "node '$name' config",
    )
    _require_named_fields(props, (), (), "node '$name' props")
    return modal_opd_expansion_node(
        name;
        pupil_rows=_file_integer(
            config.pupil_rows,
            "node '$name' config.pupil_rows",
        ),
        pupil_columns=_file_integer(
            config.pupil_columns,
            "node '$name' config.pupil_columns",
        ),
        mode_count=_file_integer(
            config.mode_count,
            "node '$name' config.mode_count",
        ),
        coefficients_schema=_file_string(
            config.coefficients_schema,
            "node '$name' config.coefficients_schema",
        ),
        opd_schema=_file_string(
            config.opd_schema,
            "node '$name' config.opd_schema",
        ),
        basis_schema=_file_string(
            config.basis_schema,
            "node '$name' config.basis_schema",
        ),
        pupil_support_schema=_file_string(
            config.pupil_support_schema,
            "node '$name' config.pupil_support_schema",
        ),
        T=Float32,
    )
end

"""
    builtin_graph_node_types()

Return the immutable node-type map accepted by the maintained TOML loader.
Optional packages may merge additional cold factories into this named tuple
and pass the result explicitly to [`load_algorithm_graph`](@ref). No mutable
global registry or dynamic Julia evaluation is used.
"""
function builtin_graph_node_types()
    return (
        discrete_integrator_f32=_discrete_integrator_f32_node,
        modal_opd_expansion_f32=_modal_opd_expansion_f32_node,
    )
end

function _parse_file_node(entry::AbstractDict, node_types::NamedTuple, index::Int)
    context = "nodes[$index]"
    _require_table_keys(
        entry,
        ("name", "type", "config", "props"),
        ("name", "type", "config"),
        context,
    )
    name = _file_identifier(entry["name"], "$context.name")
    node_type = _file_identifier(entry["type"], "$context.type")
    hasproperty(node_types, node_type) || throw(_file_error(
        "$context.type",
        "unknown graph node type '$node_type'",
    ))
    config = _normalize_toml_value(
        _table_or_empty(entry, "config", context),
    )
    props = _normalize_toml_value(
        _table_or_empty(entry, "props", context),
    )
    factory = getproperty(node_types, node_type)
    definition = factory(name, config, props)
    definition isa AlgorithmNodeDefinition || throw(_file_error(
        context,
        "node type '$node_type' did not return an AlgorithmNodeDefinition",
    ))
    return definition
end

function _parse_file_input(entry::AbstractDict, bindings::NamedTuple, index::Int)
    context = "inputs[$index]"
    fields = ("name", "destination", "binding")
    _require_table_keys(entry, fields, fields, context)
    name = _file_identifier(entry["name"], "$context.name")
    destination = _file_endpoint(entry["destination"], "$context.destination")
    values = _require_binding(bindings, entry["binding"], "$context.binding")
    return graph_input(name, destination, values)
end

function _parse_file_output(entry::AbstractDict, index::Int)
    context = "outputs[$index]"
    fields = ("name", "source")
    _require_table_keys(entry, fields, fields, context)
    name = _file_identifier(entry["name"], "$context.name")
    source = _file_endpoint(entry["source"], "$context.source")
    return graph_output(name, source)
end

function _parse_file_link(entry::AbstractDict, index::Int)
    context = "links[$index]"
    fields = ("source", "destination")
    _require_table_keys(entry, fields, fields, context)
    source = _file_endpoint(entry["source"], "$context.source")
    destination = _file_endpoint(entry["destination"], "$context.destination")
    return link(source, destination)
end

function _parse_file_delayed_link(
    entry::AbstractDict,
    bindings::NamedTuple,
    index::Int,
)
    context = "delayed_links[$index]"
    fields = ("source", "destination", "initial")
    _require_table_keys(entry, fields, fields, context)
    source = _file_endpoint(entry["source"], "$context.source")
    destination = _file_endpoint(entry["destination"], "$context.destination")
    initial = _require_binding(bindings, entry["initial"], "$context.initial")
    return delayed_link(source, destination, initial)
end

function _parse_file_parameter(
    entry::AbstractDict,
    bindings::NamedTuple,
    index::Int,
)
    context = "parameters[$index]"
    fields = ("destination", "binding")
    _require_table_keys(entry, fields, fields, context)
    destination = _file_endpoint(entry["destination"], "$context.destination")
    values = _require_binding(bindings, entry["binding"], "$context.binding")
    return sparse_parameter(destination, values)
end

function _parse_entries(entries::AbstractVector, parse_entry)
    builder = Vector{Any}(undef, length(entries))
    for (index, entry) in pairs(entries)
        builder[index] = parse_entry(entry, Int(index))
    end
    return Tuple(builder)
end

"""
    load_algorithm_graph(path; node_types=builtin_graph_node_types(),
                         bindings=NamedTuple())

Load a versioned, declarative TOML graph and compile it into the same concrete
[`AlgorithmGraphDefinition`](@ref) used by direct Julia construction.
`node_types` is an explicit immutable named tuple of cold factories. `bindings`
supplies caller-owned frame inputs, delayed-link initial values, and large
sparse parameters by name; the TOML file never embeds or implicitly loads
calibration arrays.
"""
function load_algorithm_graph(
    path::AbstractString;
    node_types::NamedTuple=builtin_graph_node_types(),
    bindings::NamedTuple=NamedTuple(),
)
    root = try
        TOML.parsefile(path)
    catch error
        throw(AlgorithmGraphError(
            "could not parse algorithm graph '$path': $(sprint(showerror, error))",
        ))
    end
    _require_table_keys(
        root,
        (
            "schema_version",
            "name",
            "nodes",
            "inputs",
            "outputs",
            "links",
            "delayed_links",
            "parameters",
        ),
        ("schema_version", "name", "nodes"),
        "algorithm graph '$path'",
    )
    version = _file_integer(root["schema_version"], "schema_version")
    version == ALGORITHM_GRAPH_TOML_SCHEMA_VERSION || throw(AlgorithmGraphError(
        "unsupported algorithm graph schema_version $version; expected " *
        "$ALGORITHM_GRAPH_TOML_SCHEMA_VERSION",
    ))
    name = _file_identifier(root["name"], "name")

    node_entries = _table_array_or_empty(root, "nodes", "algorithm graph '$path'")
    input_entries = _table_array_or_empty(root, "inputs", "algorithm graph '$path'")
    output_entries = _table_array_or_empty(root, "outputs", "algorithm graph '$path'")
    link_entries = _table_array_or_empty(root, "links", "algorithm graph '$path'")
    delayed_entries = _table_array_or_empty(
        root,
        "delayed_links",
        "algorithm graph '$path'",
    )
    parameter_entries = _table_array_or_empty(
        root,
        "parameters",
        "algorithm graph '$path'",
    )

    nodes = _parse_entries(
        node_entries,
        (entry, index) -> _parse_file_node(entry, node_types, index),
    )
    inputs = _parse_entries(
        input_entries,
        (entry, index) -> _parse_file_input(entry, bindings, index),
    )
    outputs = _parse_entries(output_entries, _parse_file_output)
    links = _parse_entries(link_entries, _parse_file_link)
    delayed_links = _parse_entries(
        delayed_entries,
        (entry, index) -> _parse_file_delayed_link(entry, bindings, index),
    )
    parameters = _parse_entries(
        parameter_entries,
        (entry, index) -> _parse_file_parameter(entry, bindings, index),
    )
    return algorithm_graph(
        nodes;
        name=name,
        inputs=inputs,
        outputs=outputs,
        links=links,
        delayed_links=delayed_links,
        parameters=parameters,
    )
end

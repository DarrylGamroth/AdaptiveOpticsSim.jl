"""
    compose_proper_camera_graph(graph::AbstractDict; bindings=NamedTuple(),
        output_shape, output_sampling_m, array_center_m=(0, 0), output_schema,
        name, pixel_node_name="camera_pixels")

Return an independent declarative graph that integrates one scalar Float32 Proper
train's intensity onto rectangular camera pixels. `graph` must have exactly one
`proper_optical_train_f32` node, one `intensity` output from that node's `output`
port, and no links, delayed links or parameters. Caller-owned `bindings` supply
the optical factory's construction calibrations; frame inputs are preserved and
need not be supplied during composition. Load `Proper` before calling.

Geometry is in logical `(row, column)` order. `output_shape` gives the detector
ROI dimensions and `output_sampling_m` gives its physical pixel pitches in metres.
`array_center_m` explicitly places the geometric center of that ROI relative to
the optical axis. It defaults to zero and is never inferred from an ROI size.
The factory plan supplies native optical sampling and the FFT-centered optical
axis at `resolution ÷ 2 + 1`. The detector's fractional one-based pixel coordinate
at physical zero is `(output_shape + 1)/2 - array_center_m/output_sampling_m`.

The appended `pixel_integration_f32` node preserves cell-integrated flux, using
the source's intensity schema and caller's `output_schema`. Composition performs
no file I/O, accelerator preparation, radiometric normalization, exposure or
detector response. The source graph and supplied arrays are not modified.
"""
function compose_proper_camera_graph(graph::AbstractDict;
        bindings::NamedTuple=NamedTuple(), output_shape, output_sampling_m,
        array_center_m=(0, 0), output_schema::AbstractString,
        name::AbstractString, pixel_node_name::AbstractString="camera_pixels")
    context = "Proper camera graph"
    fields = ("schema_version", "name", "nodes", "inputs", "outputs", "links",
        "delayed_links", "parameters")
    _require_table_keys(graph, fields, ("schema_version", "name", "nodes", "outputs"), context)
    _file_integer(graph["schema_version"], "schema_version") == ALGORITHM_GRAPH_SCHEMA_VERSION ||
        throw(AlgorithmGraphError("unsupported camera graph schema_version"))
    _file_identifier(graph["name"], "source graph name")
    _file_identifier(name, "camera graph name")
    pixel_name = _file_identifier(pixel_node_name, "camera pixel node name")
    nodes = _table_array_or_empty(graph, "nodes", context)
    length(nodes) == 1 || throw(AlgorithmGraphError("camera source requires exactly one optical node"))
    node = only(nodes)
    get(node, "type", nothing) == "proper_optical_train_f32" ||
        throw(AlgorithmGraphError("camera source requires a scalar Float32 Proper optical train"))
    source_name = _file_identifier(get(node, "name", ""), "optical node name")
    source_name == pixel_name && throw(AlgorithmGraphError("camera pixel node name collides with optical node"))
    for key in ("links", "delayed_links", "parameters")
        isempty(_table_array_or_empty(graph, key, context)) ||
            throw(AlgorithmGraphError("camera source must not contain $key"))
    end
    inputs = _table_array_or_empty(graph, "inputs", context)
    outputs = _table_array_or_empty(graph, "outputs", context)
    length(outputs) == 1 || throw(AlgorithmGraphError("camera source requires one intensity output"))
    output = only(outputs)
    _require_table_keys(output, ("name", "source"), ("name", "source"), "camera source output")
    output["name"] == "intensity" || throw(AlgorithmGraphError("camera source output must be named intensity"))
    reference = _file_endpoint(output["source"], "camera source output")
    _node_name(reference) == source_name && _port_name(reference) == :output ||
        throw(AlgorithmGraphError("camera source output must select the optical node's output port"))
    config = _table_or_empty(node, "config", context)
    any(key -> haskey(config, key), ("channels", "output_sampling_m")) &&
        throw(AlgorithmGraphError("camera source requires native scalar Proper sampling"))
    get(config, "output_quantity", "intensity") == "intensity" ||
        throw(AlgorithmGraphError("camera source requires intensity output"))
    definition = _parse_file_node(node,
        (; proper_optical_train_f32=proper_optical_train_node_factory(; T=Float32)), bindings, 1)
    ports = graph_node_ports(_node_type(definition), definition.config)
    input_names = Tuple(map(inputs) do entry
        _require_table_keys(entry, ("name", "destination", "binding"),
            ("name", "destination", "binding"), "camera source input")
        input_name = _file_identifier(entry["name"], "camera source input name")
        _file_identifier(entry["binding"], "camera source input binding")
        reference = _file_endpoint(entry["destination"], "camera source input destination")
        _node_name(reference) == source_name && any(port ->
            _port_direction(port) == :input && _port_name(port) == _port_name(reference), ports) ||
            throw(AlgorithmGraphError("camera source input must select an optical input port"))
        return input_name
    end)
    _require_unique_names(input_names, "camera source input")
    destinations = Tuple(_port_name(_file_endpoint(entry["destination"], "camera source input")) for entry in inputs)
    _require_unique_names(destinations, "camera source input destination")
    Set(destinations) == Set(_port_name(port) for port in ports if _port_direction(port) == :input) ||
        throw(AlgorithmGraphError("camera source must bind every optical input port"))
    source_port = only(filter(port -> _port_direction(port) == :output, ports))
    source_shape = source_port.format.shape
    pitch = Float64(optical_sampling_m(definition.config))
    shape = _pixel_integration_shape(output_shape isa AbstractVector ? Tuple(output_shape) : output_shape,
        "camera output_shape")
    sampling = _pixel_integration_pair(output_sampling_m isa AbstractVector ? Tuple(output_sampling_m) : output_sampling_m,
        "camera output_sampling_m"; positive=true)
    center = _pixel_integration_pair(array_center_m isa AbstractVector ? Tuple(array_center_m) : array_center_m,
        "camera array_center_m")
    source_center = map(n -> Float64(fld(n, 2) + 1), source_shape)
    destination_center = ntuple(axis -> (shape[axis] + 1)/2 - center[axis]/sampling[axis], 2)
    pixel_config = (; input_shape=source_shape,
        output_shape=shape, input_sampling_m=(pitch, pitch), output_sampling_m=sampling,
        input_center_pixel=source_center, output_center_pixel=destination_center,
        input_schema=source_port.format.schema, output_schema=String(output_schema))
    pixel_integration_node(pixel_name; pixel_config..., T=Float32)
    pixel_entry = Dict("name" => String(pixel_node_name), "type" => "pixel_integration_f32",
        "config" => Dict(String(key) => (value isa Tuple ? collect(value) : value)
            for (key, value) in pairs(pixel_config)))
    composed = Dict(key => deepcopy(value) for (key, value) in graph)
    composed["name"] = String(name)
    composed["nodes"] = [deepcopy(node), pixel_entry]
    composed["links"] = [Dict("source" => String(output["source"]),
        "destination" => "$(pixel_node_name).input")]
    composed["outputs"] = [Dict("name" => "intensity", "source" => "$(pixel_node_name).output")]
    return composed
end

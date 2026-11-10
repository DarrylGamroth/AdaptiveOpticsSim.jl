function camera_source(definition, node, port)
    destination = AOG.AlgorithmPortReference(node, port)
    sources = filter(link -> link.destination == destination, definition.links)
    length(sources) == 1 ||
        throw(ArgumentError("PWFS camera requires one direct link into $node:$port"))
    return only(sources).source
end

function camera_node(definition, name)
    matches = filter(node -> AOG._node_name(node) === name, definition.nodes)
    length(matches) == 1 || throw(ArgumentError("PWFS camera has no node named $name"))
    return only(matches)
end

function pwfs_camera_nodes(definition)
    outputs = filter(output -> AOG._boundary_name(output) === :frame, definition.outputs)
    length(outputs) == 1 || throw(ArgumentError("PWFS camera requires a frame output"))
    frame = only(outputs).source
    layout = camera_node(definition, AOG._node_name(frame))
    detector_frame = camera_source(definition, AOG._node_name(layout), :detector_frame)
    detector = camera_node(definition, AOG._node_name(detector_frame))
    photon_rate = camera_source(definition, AOG._node_name(detector), :photon_rate)
    pyramid = camera_node(definition, AOG._node_name(photon_rate))
    map(AOG._port_name, (photon_rate, detector_frame, frame)) == (:photon_rate, :frame, :frame) ||
        throw(ArgumentError("PWFS camera requires Pyramid photon rate, EMCCD frame, and raw frame layout"))
    pyramid_rate_config(pyramid)
    emccd_config(detector)
    detector_frame_layout_config(layout)
    return pyramid, detector, layout
end

function O.pwfs_camera_filter_graph(definition::AOG.AlgorithmGraphDefinition)
    nodes = pwfs_camera_nodes(definition)
    names = map(AOG._node_name, nodes)
    pyramid, detector, layout = names
    algorithms = (PyramidRateF32, EMCCDF32, DetectorFrameLayoutF32)
    declarations = map(nodes, names, algorithms) do original, name, algorithm
        (; name, label=FGA.algorithm_label(algorithm), config=(; definition=original))
    end
    document = (; nodes=declarations, inputs=(; opd="$pyramid:opd"),
        links=("$pyramid:photon_rate" => "$detector:photon_rate",
            "$detector:frame" => "$layout:detector_frame"),
        outputs=(; frame="$layout:frame"))
    return (; document, algorithms)
end

function O.with_filtergraphs_pwfs_camera(definition::AOG.AlgorithmGraphDefinition)
    names = map(AOG._node_name, pwfs_camera_nodes(definition))
    return replace_camera_nodes(definition, names)
end

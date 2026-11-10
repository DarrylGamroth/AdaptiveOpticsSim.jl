function camera_source(definition, node, port)
    destination = AOG.AlgorithmPortReference(node, port)
    sources = filter(link -> link.destination == destination, definition.links)
    length(sources) == 1 ||
        throw(ArgumentError("camera requires one direct link into $node:$port"))
    return only(sources).source
end

function camera_node(definition, name)
    matches = filter(node -> AOG._node_name(node) === name, definition.nodes)
    length(matches) == 1 || throw(ArgumentError("camera has no node named $name"))
    return only(matches)
end

function camera_chain_nodes(definition, frame_output::Symbol)
    outputs = filter(output -> AOG._boundary_name(output) === frame_output, definition.outputs)
    length(outputs) == 1 || throw(ArgumentError("camera requires a raw frame output named $frame_output"))
    frame = only(outputs).source
    layout = camera_node(definition, AOG._node_name(frame))
    detector_frame = camera_source(definition, AOG._node_name(layout), :detector_frame)
    detector = camera_node(definition, AOG._node_name(detector_frame))
    photon_rate = camera_source(definition, AOG._node_name(detector), :photon_rate)
    optics = camera_node(definition, AOG._node_name(photon_rate))
    map(AOG._port_name, (photon_rate, detector_frame, frame)) == (:photon_rate, :frame, :frame) ||
        throw(ArgumentError("camera requires optical photon rate, a detector frame, and raw frame layout"))
    names = map(AOG._node_name, (optics, detector, layout))
    length(unique(names)) == 3 ||
        throw(ArgumentError("camera requires distinct optics, detector, and layout nodes"))
    return optics, detector, layout
end

function pwfs_camera_nodes(definition)
    pyramid, detector, layout = camera_chain_nodes(definition, :frame)
    pyramid_rate_config(pyramid)
    emccd_config(detector)
    detector_frame_layout_config(layout)
    return pyramid, detector, layout
end

function O.pwfs_camera_filter_graph(definition::AOG.AlgorithmGraphDefinition)
    return pwfs_camera_description(pwfs_camera_nodes(definition))
end

pwfs_camera_description(nodes) =
    camera_description(nodes, (PyramidRateF32, EMCCDF32, DetectorFrameLayoutF32))

function camera_description(nodes, algorithms)
    names = map(AOG._node_name, nodes)
    optics, detector, layout = names
    declarations = map(nodes, names, algorithms) do original, name, algorithm
        (; name, label=FGA.algorithm_label(algorithm), config=(; definition=original))
    end
    document = (; nodes=declarations, inputs=(; opd="$optics:opd"),
        links=("$optics:photon_rate" => "$detector:photon_rate",
            "$detector:frame" => "$layout:detector_frame"),
        outputs=(; frame="$layout:frame"))
    return (; document, algorithms)
end

function O.with_filtergraphs_pwfs_camera(definition::AOG.AlgorithmGraphDefinition)
    names = map(AOG._node_name, pwfs_camera_nodes(definition))
    return replace_camera_nodes(definition, names)
end

function shack_hartmann_camera_nodes(definition, frame_output::Symbol)
    optics, detector, layout = camera_chain_nodes(definition, frame_output)
    shack_hartmann_rate_config(optics)
    cmos_config(detector)
    detector_frame_layout_config(layout)
    return optics, detector, layout
end

function O.shack_hartmann_camera_filter_graph(
        definition::AOG.AlgorithmGraphDefinition; frame_output::Symbol=:frame)
    nodes = shack_hartmann_camera_nodes(definition, frame_output)
    return camera_description(nodes, (ShackHartmannRateF32, CMOSF32, DetectorFrameLayoutF32))
end

function O.with_filtergraphs_shack_hartmann_camera(
        definition::AOG.AlgorithmGraphDefinition; frame_output::Symbol=:frame)
    names = map(AOG._node_name, shack_hartmann_camera_nodes(definition, frame_output))
    return replace_camera_nodes(definition, names)
end

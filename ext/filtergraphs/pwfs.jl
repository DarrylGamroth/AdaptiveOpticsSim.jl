function O.pwfs_filter_graph(definition::AOG.AlgorithmNodeDefinition)
    pyramid_rate_config(definition)
    name = AOG._node_name(definition)
    node = (; name, label=FGA.algorithm_label(PyramidRateF32), config=(; definition))
    document = (; nodes=(node,), inputs=(; opd="$name:opd"),
        outputs=(; photon_rate="$name:photon_rate"))
    return (; document, algorithms=(PyramidRateF32,))
end

# AOS retains its surrounding schedule and binds the Algorithms to its ports.
struct FilterGraphsCameraNode{N} end

function O.with_filtergraphs_pwfs(definition::AOG.AlgorithmGraphDefinition; node::Symbol=:pyramid)
    pyramid_rate_config(camera_node(definition, node))
    return replace_camera_nodes(definition, (node,))
end

function replace_camera_nodes(definition, names)
    nodes = map(definition.nodes) do original
        name = AOG._node_name(original)
        name in names || return original
        return AOG.algorithm_node(name, FilterGraphsCameraNode{AOG._node_type(original)}, original)
    end
    return AOG.algorithm_graph(nodes; definition.name, definition.inputs, definition.outputs,
        definition.links, definition.delayed_links, definition.parameters)
end

function AOG.graph_node_ports(::Type{FilterGraphsCameraNode{N}}, definition::AOG.AlgorithmNodeDefinition) where {N}
    return AOG.graph_node_ports(N, definition.config)
end

struct BoundCameraAlgorithm{P,I,O}
    algorithm::P
    input::I
    output::O
end

camera_algorithm(::Type{AOG.PyramidRateNode{Float32}}) = (PyramidRateF32, PyramidRateConfiguration)
camera_algorithm(::Type{AOG.EMCCDDetectorAcquisitionNode{Float32}}) = (EMCCDF32, EMCCDConfiguration)
camera_algorithm(::Type{AOG.DetectorFrameLayoutNode}) = (DetectorFrameLayoutF32, DetectorFrameLayoutConfiguration)

function AOG.prepare_graph_node(::Type{FilterGraphsCameraNode{N}}, definition::AOG.AlgorithmNodeDefinition,
        ::NamedTuple{()}, inputs::NamedTuple, outputs::NamedTuple,
        ::NamedTuple{()}, ::B.HostComputeDevice) where {N}
    declaration, configuration = camera_algorithm(N)
    algorithm = FGA.prepare_algorithm(declaration, configuration(definition))
    return BoundCameraAlgorithm(algorithm, only(inputs), only(outputs))
end

captured_camera_algorithm(::Type{AOG.PyramidRateNode{Float32}}) = prepare_captured_pyramid
captured_camera_algorithm(::Type{AOG.EMCCDDetectorAcquisitionNode{Float32}}) = prepare_captured_emccd
captured_camera_algorithm(::Type{AOG.DetectorFrameLayoutNode}) = prepare_captured_detector_frame_layout

function AOG.prepare_graph_node(::Type{FilterGraphsCameraNode{N}}, definition::AOG.AlgorithmNodeDefinition,
        ::NamedTuple{()}, inputs::NamedTuple, outputs::NamedTuple,
        ::NamedTuple{()}, ::B.AcceleratorComputeDevice) where {N}
    return captured_camera_algorithm(N)(definition, only(inputs), only(outputs))
end

function AOG.step_graph_node!(node::BoundCameraAlgorithm)
    FGA.process!(node.output, node.algorithm, node.input)
    return nothing
end

AOG.reset_graph_node!(node::BoundCameraAlgorithm) = FGA.reset!(node.algorithm)

const CapturedCameraAlgorithm = Union{CapturedPyramidRate,CapturedEMCCD,CapturedDetectorFrameLayout}
AOG.graph_node_capture_capability(::CapturedCameraAlgorithm) = AOG.GraphNodeCaptureSafe()
AOG.enqueue_captured_graph_node!(node::CapturedCameraAlgorithm) = FGA.enqueue!(node)
AOG.enqueue_graph_node!(node::CapturedPyramidRate) = FGA.enqueue!(node)
AOG.step_graph_node!(node::CapturedCameraAlgorithm) = AOG.step_graph_node!(node.owner)
AOG.reset_graph_node!(node::CapturedCameraAlgorithm) = FGA.reset!(node)

AOG._preflight_captured_graph_node!(node::CapturedEMCCD) = FGA.preflight(node)
AOG._complete_captured_graph_node!(node::CapturedEMCCD) = FGA.complete!(node, Int32(0))

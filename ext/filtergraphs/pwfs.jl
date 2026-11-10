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

function AOG.prepare_graph_node(::Type{FilterGraphsCameraNode{N}}, definition::AOG.AlgorithmNodeDefinition,
        ::NamedTuple{()}, inputs::NamedTuple, outputs::NamedTuple,
        ::NamedTuple{()}, ::B.HostComputeDevice) where {N}
    input, output = only(inputs), only(outputs)
    algorithm = prepare_cpu_camera_algorithm(N, definition, input, output)
    return BoundCameraAlgorithm(algorithm, input, output)
end

# AOS admits the routes and allocates their arrays before preparing the nodes.
# The ordinary FilterGraphs CPU API can receive different arrays on each call.
function prepare_cpu_camera_algorithm(::Type{AOG.PyramidRateNode{Float32}}, definition, input, output)
    return prepare_bound_camera_algorithm(PyramidRateF32, definition,
        prepare_host_pyramid(definition, input, output))
end

function prepare_cpu_camera_algorithm(::Type{AOG.EMCCDDetectorAcquisitionNode{Float32}}, definition, input, output)
    return prepare_bound_camera_algorithm(EMCCDF32, definition,
        prepare_host_emccd(definition, input, output))
end

function prepare_cpu_camera_algorithm(::Type{AOG.DetectorFrameLayoutNode}, definition, input, output)
    return FGA.prepare_algorithm(DetectorFrameLayoutF32, DetectorFrameLayoutConfiguration(definition))
end

function prepare_bound_camera_algorithm(declaration, definition, prepared)
    ports = AOG.graph_node_ports(AOG._node_type(definition), definition.config)
    formats = map(port -> camera_algorithm_format(port.format), ports)
    return FGA.PreparedAlgorithm(declaration, prepared.plan, prepared.workspace, formats)
end

camera_algorithm_format(format::AOG._GraphPortFormat{T}) where {T} =
    FGA.AlgorithmFormat(T, format.shape, format.schema)

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

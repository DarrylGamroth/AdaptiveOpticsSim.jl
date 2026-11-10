function O.pwfs_filter_graph(definition::AOG.AlgorithmNodeDefinition)
    pyramid_rate_config(definition)
    name = AOG._node_name(definition)
    node = (; name, label=FGA.algorithm_label(PyramidRateF32), config=(; definition))
    document = (; nodes=(node,), inputs=(; opd="$name:opd"),
        outputs=(; photon_rate="$name:photon_rate"))
    return (; document, algorithms=(PyramidRateF32,))
end

# The surrounding AOS graph still owns source, detector, and model time. Bind
# the FilterGraphAlgorithms node at the same point in the existing schedule.
struct FilterGraphsPyramidRate end

function O.with_filtergraphs_pwfs(definition::AOG.AlgorithmGraphDefinition; node::Symbol=:pyramid)
    any(candidate -> AOG._node_name(candidate) === node, definition.nodes) ||
        throw(ArgumentError("algorithm graph has no PWFS node named $node"))
    nodes = map(definition.nodes) do original
        AOG._node_name(original) === node || return original
        pyramid_rate_config(original)
        return AOG.algorithm_node(node, FilterGraphsPyramidRate, original)
    end
    return AOG.algorithm_graph(nodes; definition.name, definition.inputs, definition.outputs,
        definition.links, definition.delayed_links, definition.parameters)
end

function AOG.graph_node_ports(::Type{FilterGraphsPyramidRate}, definition::AOG.AlgorithmNodeDefinition)
    return AOG.graph_node_ports(AOG.PyramidRateNode{Float32}, pyramid_rate_config(definition))
end

struct BoundPyramidRate{P,I,O}
    algorithm::P
    opd::I
    photon_rate::O
end

function AOG.prepare_graph_node(::Type{FilterGraphsPyramidRate}, definition::AOG.AlgorithmNodeDefinition,
        ::NamedTuple{()}, inputs::NamedTuple{(:opd,)}, outputs::NamedTuple{(:photon_rate,)},
        ::NamedTuple{()}, ::B.HostComputeDevice)
    algorithm = FGA.prepare_algorithm(PyramidRateF32, PyramidRateConfiguration(definition))
    return BoundPyramidRate(algorithm, inputs.opd, outputs.photon_rate)
end

function AOG.prepare_graph_node(::Type{FilterGraphsPyramidRate}, definition::AOG.AlgorithmNodeDefinition,
        ::NamedTuple{()}, inputs::NamedTuple{(:opd,)}, outputs::NamedTuple{(:photon_rate,)},
        ::NamedTuple{()}, ::B.AcceleratorComputeDevice)
    return prepare_captured_pyramid(definition, inputs.opd, outputs.photon_rate)
end

function AOG.step_graph_node!(node::BoundPyramidRate)
    FGA.process!(node.photon_rate, node.algorithm, node.opd)
    return nothing
end

AOG.reset_graph_node!(node::BoundPyramidRate) = FGA.reset!(node.algorithm)

AOG.graph_node_capture_capability(::CapturedPyramidRate) = AOG.GraphNodeCaptureSafe()
AOG.enqueue_graph_node!(node::CapturedPyramidRate) = FGA.enqueue!(node)
AOG.step_graph_node!(node::CapturedPyramidRate) = AOG.step_graph_node!(node.owner)
AOG.reset_graph_node!(node::CapturedPyramidRate) = FGA.reset!(node)

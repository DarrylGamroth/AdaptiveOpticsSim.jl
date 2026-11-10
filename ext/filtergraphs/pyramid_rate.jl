function pyramid_rate_config(definition)
    config = definition.config
    AOG._node_type(definition) === AOG.PyramidRateNode{Float32} &&
        config isa AOG.PyramidRateNodeConfig{Float32} &&
        config.input_quantity === :opd && isempty(definition.props) ||
        throw(ArgumentError("pyramid rate requires a Float32 OPD-input definition with no properties"))
    return config
end

function pyramid_rate_format(definition, name)
    ports = AOG.graph_node_ports(AOG.PyramidRateNode{Float32}, pyramid_rate_config(definition))
    return only(port for port in ports if AOG._port_name(port) === name).format
end

function prepare_pyramid_rate(definition, opd, photon_rate)
    config = pyramid_rate_config(definition)
    return AOG.prepare_graph_node(AOG.PyramidRateNode{Float32}, config,
        (;), (; opd), (; photon_rate), (;), B.compute_device(photon_rate))
end

struct PyramidRateWorkspace{O,I,R}
    owner::O
    opd::I
    photon_rate::R
end

FGA.commit!(::PyramidRateWorkspace) = nothing
FGA.discard!(::PyramidRateWorkspace) = nothing
FGA.reset!(workspace::PyramidRateWorkspace) = AOG.reset_graph_node!(workspace.owner)

function prepare_host_pyramid(definition)
    opd = zeros(Float32, pyramid_rate_format(definition, :opd).shape)
    photon_rate = zeros(Float32, pyramid_rate_format(definition, :photon_rate).shape)
    owner = prepare_pyramid_rate(definition, opd, photon_rate)
    return (; plan=definition, workspace=PyramidRateWorkspace(owner, opd, photon_rate))
end

function process_pyramid_rate!(photon_rate, workspace::PyramidRateWorkspace, opd)
    copyto!(workspace.opd, opd)
    AOG.step_graph_node!(workspace.owner)
    copyto!(photon_rate, workspace.photon_rate)
    return nothing
end

@algorithm PyramidRateF32 begin
    label = "aos-pyramid-rate-f32"

    @configuration PyramidRateConfiguration begin
        definition::AOG.AlgorithmNodeDefinition => graph_rebuild
    end

    ports(c) = (
        input(:opd, Float32, pyramid_rate_format(c.definition, :opd).shape,
            pyramid_rate_format(c.definition, :opd).schema),
        output(:photon_rate, Float32, pyramid_rate_format(c.definition, :photon_rate).shape,
            pyramid_rate_format(c.definition, :photon_rate).schema; metadata_from=:opd),
    )

    prepare(c) = prepare_host_pyramid(c.definition)
    process!(plan, workspace) = process_pyramid_rate!(photon_rate, workspace, opd)
end

struct CapturedPyramidRate{O}
    owner::O
end

FGA.capture_support(::Type{PyramidRateF32}, device) = optical_capture_support(device)

function prepare_captured_pyramid(definition, opd, photon_rate)
    owner = prepare_pyramid_rate(definition, opd, photon_rate)
    AOG.graph_node_capture_capability(owner) isa AOG.GraphNodeCaptureSafe ||
        throw(ArgumentError("pyramid rate does not support native capture on the selected device"))
    return CapturedPyramidRate(owner)
end

function FGA.prepare_capture(::FGA.CompleteTransition, device,
        prepared::FGA.PreparedAlgorithm{PyramidRateF32},
        outputs::Tuple, inputs::Tuple, status, index::Integer)
    return prepare_captured_pyramid(prepared.plan, only(inputs), only(outputs))
end

FGA.preflight(::CapturedPyramidRate) = nothing
FGA.preserves_inputs(::CapturedPyramidRate) = true
FGA.enqueue!(node::CapturedPyramidRate) = AOG.enqueue_captured_graph_node!(node.owner)
FGA.complete!(::CapturedPyramidRate, ::Int32) = nothing
FGA.reset!(node::CapturedPyramidRate) = AOG.reset_graph_node!(node.owner)

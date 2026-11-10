function shack_hartmann_rate_config(definition)
    config = definition.config
    AOG._node_type(definition) === AOG.ShackHartmannRateNode{Float32} &&
        config isa AOG.ShackHartmannRateNodeConfig{Float32} &&
        isempty(definition.props) ||
        throw(ArgumentError("Shack–Hartmann rate requires a Float32 definition with no properties"))
    return config
end

function shack_hartmann_rate_format(definition, name)
    ports = AOG.graph_node_ports(AOG.ShackHartmannRateNode{Float32},
        shack_hartmann_rate_config(definition))
    return only(port for port in ports if AOG._port_name(port) === name).format
end

function prepare_shack_hartmann_rate(definition, opd, photon_rate)
    config = shack_hartmann_rate_config(definition)
    return AOG.prepare_graph_node(AOG.ShackHartmannRateNode{Float32}, config,
        (;), (; opd), (; photon_rate), (;), B.compute_device(photon_rate))
end

struct ShackHartmannRateWorkspace{O,I,R}
    owner::O
    opd::I
    photon_rate::R
end

FGA.commit!(::ShackHartmannRateWorkspace) = nothing
FGA.discard!(::ShackHartmannRateWorkspace) = nothing
FGA.reset!(workspace::ShackHartmannRateWorkspace) =
    AOG.reset_graph_node!(workspace.owner)

function prepare_host_shack_hartmann(definition)
    opd = zeros(Float32, shack_hartmann_rate_format(definition, :opd).shape)
    photon_rate = zeros(Float32,
        shack_hartmann_rate_format(definition, :photon_rate).shape)
    return prepare_host_shack_hartmann(definition, opd, photon_rate)
end

function prepare_host_shack_hartmann(definition, opd, photon_rate)
    owner = prepare_shack_hartmann_rate(definition, opd, photon_rate)
    return (; plan=definition, workspace=ShackHartmannRateWorkspace(owner, opd, photon_rate))
end

function process_shack_hartmann_rate!(photon_rate,
        workspace::ShackHartmannRateWorkspace, opd)
    opd === workspace.opd || copyto!(workspace.opd, opd)
    AOG.step_graph_node!(workspace.owner)
    photon_rate === workspace.photon_rate || copyto!(photon_rate, workspace.photon_rate)
    return nothing
end

@algorithm ShackHartmannRateF32 begin
    label = "aos-shack-hartmann-rate-f32"

    @configuration ShackHartmannRateConfiguration begin
        definition::AOG.AlgorithmNodeDefinition => graph_rebuild
    end

    ports(c) = (
        input(:opd, Float32, shack_hartmann_rate_format(c.definition, :opd).shape,
            shack_hartmann_rate_format(c.definition, :opd).schema),
        output(:photon_rate, Float32,
            shack_hartmann_rate_format(c.definition, :photon_rate).shape,
            shack_hartmann_rate_format(c.definition, :photon_rate).schema;
            metadata_from=:opd),
    )

    prepare(c) = prepare_host_shack_hartmann(c.definition)
    process!(plan, workspace) =
        process_shack_hartmann_rate!(photon_rate, workspace, opd)
end

struct CapturedShackHartmannRate{O}
    owner::O
end

FGA.capture_support(::Type{ShackHartmannRateF32}, device) =
    optical_capture_support(device)

function prepare_captured_shack_hartmann(definition, opd, photon_rate)
    owner = prepare_shack_hartmann_rate(definition, opd, photon_rate)
    AOG.graph_node_capture_capability(owner) isa AOG.GraphNodeCaptureSafe ||
        throw(ArgumentError(
            "Shack–Hartmann rate does not support native capture on the selected device"))
    return CapturedShackHartmannRate(owner)
end

function FGA.prepare_capture(::FGA.CompleteTransition, device,
        prepared::FGA.PreparedAlgorithm{ShackHartmannRateF32},
        outputs::Tuple, inputs::Tuple, status, index::Integer)
    return prepare_captured_shack_hartmann(prepared.plan, only(inputs), only(outputs))
end

FGA.preflight(::CapturedShackHartmannRate) = nothing
FGA.preserves_inputs(::CapturedShackHartmannRate) = true
FGA.enqueue!(node::CapturedShackHartmannRate) =
    AOG.enqueue_captured_graph_node!(node.owner)
FGA.complete!(::CapturedShackHartmannRate, ::Int32) = nothing
FGA.reset!(node::CapturedShackHartmannRate) =
    AOG.reset_graph_node!(node.owner)

function O.shack_hartmann_rate_filter_graph(definition::AOG.AlgorithmNodeDefinition)
    shack_hartmann_rate_config(definition)
    name = AOG._node_name(definition)
    node = (; name, label=FGA.algorithm_label(ShackHartmannRateF32),
        config=(; definition))
    document = (; nodes=(node,), inputs=(; opd="$name:opd"),
        outputs=(; photon_rate="$name:photon_rate"))
    return (; document, algorithms=(ShackHartmannRateF32,))
end

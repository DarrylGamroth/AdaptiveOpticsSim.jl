function cmos_config(definition::AOG.AlgorithmNodeDefinition)
    config = definition.config
    AOG._node_type(definition) === AOG.CMOSDetectorAcquisitionNode{Float32} &&
        config isa AOG.CMOSDetectorAcquisitionNodeConfig{Float32} &&
        isempty(definition.props) || throw(
        ArgumentError("CMOS acquisition requires a Float32 definition with no properties"))
    return config
end

function cmos_ports(definition)
    return AOG.graph_node_ports(AOG.CMOSDetectorAcquisitionNode{Float32}, cmos_config(definition))
end

function prepare_cmos(definition, photon_rate, frame)
    return AOG.prepare_graph_node(AOG.CMOSDetectorAcquisitionNode{Float32},
        cmos_config(definition), (;), (; photon_rate), (; frame), (;), B.compute_device(frame))
end

struct CMOSWorkspace{O}
    owner::O
end

FGA.commit!(::CMOSWorkspace) = nothing
FGA.discard!(::CMOSWorkspace) = nothing
FGA.reset!(workspace::CMOSWorkspace) = AOG.reset_graph_node!(workspace.owner)

function prepare_host_cmos(definition)
    input, output = cmos_ports(definition)
    return prepare_host_cmos(definition, zeros(Float32, input.format.shape),
        zeros(Float32, output.format.shape))
end

function prepare_host_cmos(definition, photon_rate, frame)
    return (; plan=definition, workspace=CMOSWorkspace(prepare_cmos(definition, photon_rate, frame)))
end

function process_cmos!(frame, workspace::CMOSWorkspace, photon_rate)
    owner = workspace.owner
    bound_rate = owner.prepared.input.values
    photon_rate === bound_rate || copyto!(bound_rate, photon_rate)
    AOG.step_graph_node!(owner)
    frame === owner.output || copyto!(frame, owner.output)
    return nothing
end

@algorithm CMOSF32 begin
    label = "aos-cmos-f32"

    @configuration CMOSConfiguration begin
        definition::AOG.AlgorithmNodeDefinition => graph_rebuild
    end

    ports(c) = (
        input(:photon_rate, Float32, cmos_ports(c.definition)[1].format.shape,
            cmos_ports(c.definition)[1].format.schema),
        output(:frame, Float32, cmos_ports(c.definition)[2].format.shape,
            cmos_ports(c.definition)[2].format.schema; metadata_from=:photon_rate),
    )

    prepare(c) = prepare_host_cmos(c.definition)
    process!(plan, workspace) = process_cmos!(frame, workspace, photon_rate)
end

function O.cmos_filter_graph(definition::AOG.AlgorithmNodeDefinition)
    cmos_config(definition)
    name = AOG._node_name(definition)
    node = (; name, label=FGA.algorithm_label(CMOSF32), config=(; definition))
    document = (; nodes=(node,), inputs=(; photon_rate="$name:photon_rate"),
        outputs=(; frame="$name:frame"))
    return (; document, algorithms=(CMOSF32,))
end

function prepare_cpu_camera_algorithm(::Type{AOG.CMOSDetectorAcquisitionNode{Float32}},
        definition, input, output)
    return prepare_bound_camera_algorithm(CMOSF32, definition,
        prepare_host_cmos(definition, input, output))
end

struct CapturedCMOS{O,S,R}
    owner::O
    status::S
    candidate_rng::R
end

FGA.capture_support(::Type{CMOSF32}, device) = optical_capture_support(device)

function prepare_captured_cmos(definition, photon_rate, frame, status=nothing)
    owner = prepare_cmos(definition, photon_rate, frame)
    AOG.graph_node_capture_capability(owner) isa AOG.GraphNodeCaptureSafe ||
        throw(ArgumentError("CMOS acquisition does not support native capture on the selected device"))
    candidate_rng = isnothing(status) ? nothing :
        B._prepare_counter_rng(B.compute_device(frame), owner.rng_seed)
    return CapturedCMOS(owner, status, candidate_rng)
end

function FGA.prepare_capture(::FGA.CompleteTransition, device,
        prepared::FGA.PreparedAlgorithm{CMOSF32}, outputs::Tuple, inputs::Tuple,
        status, index::Integer)
    return prepare_captured_cmos(prepared.plan, only(inputs), only(outputs), status)
end

FGA.preflight(node::CapturedCMOS) = AOG._preflight_captured_graph_node!(node.owner)
FGA.preserves_inputs(::CapturedCMOS) = true
FGA.reset!(node::CapturedCMOS) = AOG.reset_graph_node!(node.owner)
FGA.complete!(node::CapturedCMOS, ::Int32) = AOG._complete_captured_graph_node!(node.owner)

FGA.enqueue!(node::CapturedCMOS{O,Nothing}) where {O} =
    AOG.enqueue_captured_graph_node!(node.owner)

function FGA.enqueue!(node::CapturedCMOS)
    owner = node.owner
    candidate = node.candidate_rng.state.draw_sequence
    committed = owner.rng.state.draw_sequence
    B.copyto_backend_async!(candidate, committed)
    AOG._enqueue_detector_acquisition!(owner.prepared, node.candidate_rng, owner.output)
    B.launch_kernel_async!(B.execution_style(committed), commit_detector_rng!,
        committed, candidate, node.status; ndrange=1)
    return nothing
end

captured_camera_algorithm(::Type{AOG.CMOSDetectorAcquisitionNode{Float32}}) = prepare_captured_cmos
AOG._preflight_captured_graph_node!(node::CapturedCMOS) = FGA.preflight(node)
AOG._complete_captured_graph_node!(node::CapturedCMOS) = FGA.complete!(node, Int32(0))

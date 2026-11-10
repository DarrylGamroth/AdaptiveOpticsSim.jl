function emccd_config(definition)
    config = definition.config
    AOG._node_type(definition) === AOG.EMCCDDetectorAcquisitionNode{Float32} &&
        config isa AOG.EMCCDDetectorAcquisitionNodeConfig{Float32} &&
        isempty(definition.props) ||
        throw(ArgumentError("EMCCD acquisition requires a Float32 definition with no properties"))
    return config
end

function emccd_format(definition, name)
    ports = AOG.graph_node_ports(AOG.EMCCDDetectorAcquisitionNode{Float32}, emccd_config(definition))
    return only(port for port in ports if AOG._port_name(port) === name).format
end

function prepare_emccd(definition, photon_rate, frame)
    return AOG.prepare_graph_node(AOG.EMCCDDetectorAcquisitionNode{Float32},
        emccd_config(definition), (;), (; photon_rate), (; frame), (;), B.compute_device(frame))
end

struct EMCCDWorkspace{O,I,R}
    owner::O
    photon_rate::I
    frame::R
end

FGA.commit!(::EMCCDWorkspace) = nothing
FGA.discard!(::EMCCDWorkspace) = nothing
FGA.reset!(workspace::EMCCDWorkspace) = AOG.reset_graph_node!(workspace.owner)

function prepare_host_emccd(definition)
    photon_rate = zeros(Float32, emccd_format(definition, :photon_rate).shape)
    frame = zeros(Float32, emccd_format(definition, :frame).shape)
    return prepare_host_emccd(definition, photon_rate, frame)
end

function prepare_host_emccd(definition, photon_rate, frame)
    owner = prepare_emccd(definition, photon_rate, frame)
    return (; plan=definition, workspace=EMCCDWorkspace(owner, photon_rate, frame))
end

function process_emccd!(frame, workspace::EMCCDWorkspace, photon_rate)
    photon_rate === workspace.photon_rate || copyto!(workspace.photon_rate, photon_rate)
    AOG.step_graph_node!(workspace.owner)
    frame === workspace.frame || copyto!(frame, workspace.frame)
    return nothing
end

@algorithm EMCCDF32 begin
    label = "aos-emccd-f32"

    @configuration EMCCDConfiguration begin
        definition::AOG.AlgorithmNodeDefinition => graph_rebuild
    end

    ports(c) = (
        input(:photon_rate, Float32, emccd_format(c.definition, :photon_rate).shape,
            emccd_format(c.definition, :photon_rate).schema),
        output(:frame, Float32, emccd_format(c.definition, :frame).shape,
            emccd_format(c.definition, :frame).schema; metadata_from=:photon_rate),
    )

    prepare(c) = prepare_host_emccd(c.definition)
    process!(plan, workspace) = process_emccd!(frame, workspace, photon_rate)
end

struct CapturedEMCCD{O,S,R}
    owner::O
    status::S
    candidate_rng::R
end

FGA.capture_support(::Type{EMCCDF32}, device) = optical_capture_support(device)

function prepare_captured_emccd(definition, photon_rate, frame, status=nothing)
    owner = prepare_emccd(definition, photon_rate, frame)
    AOG.graph_node_capture_capability(owner) isa AOG.GraphNodeCaptureSafe ||
        throw(ArgumentError("EMCCD acquisition does not support native capture on the selected device"))
    candidate_rng = isnothing(status) ? nothing :
        B._prepare_counter_rng(B.compute_device(frame), owner.rng_seed)
    return CapturedEMCCD(owner, status, candidate_rng)
end

function FGA.prepare_capture(::FGA.CompleteTransition, device,
        prepared::FGA.PreparedAlgorithm{EMCCDF32}, outputs::Tuple, inputs::Tuple,
        status, index::Integer)
    return prepare_captured_emccd(prepared.plan, only(inputs), only(outputs), status)
end

FGA.preflight(node::CapturedEMCCD) = AOG._preflight_captured_graph_node!(node.owner)
FGA.preserves_inputs(::CapturedEMCCD) = true
FGA.reset!(node::CapturedEMCCD) = AOG.reset_graph_node!(node.owner)

function FGA.complete!(node::CapturedEMCCD, status::Int32)
    # A preceding failure suppresses this node's completion. A later failure
    # leaves a successful detector acquisition committed.
    AOG._complete_captured_graph_node!(node.owner)
    return nothing
end

FGA.enqueue!(node::CapturedEMCCD{O,Nothing}) where {O} =
    AOG.enqueue_captured_graph_node!(node.owner)

@kernel function commit_detector_rng!(committed, candidate, status)
    if iszero(@inbounds status[1])
        @inbounds committed[1] = candidate[1]
    end
end

function FGA.enqueue!(node::CapturedEMCCD)
    owner = node.owner
    candidate = node.candidate_rng.state.draw_sequence
    committed = owner.rng.state.draw_sequence
    # Draw into a private counter so a rejected preceding node cannot advance
    # the detector stream. The shared status is ordered before this commit.
    B.copyto_backend_async!(candidate, committed)
    AOG._enqueue_detector_acquisition!(owner.prepared, node.candidate_rng, owner.output)
    B.launch_kernel_async!(B.execution_style(committed), commit_detector_rng!,
        committed, candidate, node.status; ndrange=1)
    return nothing
end

struct PWFSLoanCamera end

function loan_camera_definition(definition, camera_nodes, frame_slots)
    pyramid, detector, layout = map(AOG._node_name, camera_nodes)
    camera_names = (pyramid, detector, layout)
    within_camera(port) = AOG._node_name(port) in camera_names
    any(link -> within_camera(link.source) || within_camera(link.destination), definition.delayed_links) &&
        throw(ArgumentError("PWFS camera loans require direct camera links"))
    any(parameter -> within_camera(parameter.endpoint), definition.parameters) &&
        throw(ArgumentError("PWFS camera loans require fixed camera configuration"))

    camera_input = AOG.AlgorithmPortReference(pyramid, :opd)
    loan_input = AOG.AlgorithmPortReference(layout, :opd)
    function redirect_input(port)
        within_camera(port) || return port
        port == camera_input || throw(ArgumentError("unsupported PWFS camera input"))
        return loan_input
    end
    function redirect_output(port)
        within_camera(port) || return port
        port == AOG.AlgorithmPortReference(pyramid, :photon_rate) &&
            return AOG.AlgorithmPortReference(layout, :photon_rate)
        port == AOG.AlgorithmPortReference(detector, :frame) &&
            return AOG.AlgorithmPortReference(layout, :detector_frame)
        port == AOG.AlgorithmPortReference(layout, :frame) && return port
        throw(ArgumentError("unsupported PWFS camera output"))
    end

    retained = filter(node -> !(AOG._node_name(node) in (pyramid, detector)), definition.nodes)
    nodes = map(retained) do node
        AOG._node_name(node) === layout || return node
        AOG.algorithm_node(layout, PWFSLoanCamera, (; camera_nodes, frame_slots))
    end
    external_links = filter(link -> !(within_camera(link.source) && within_camera(link.destination)), definition.links)
    length(definition.links) - length(external_links) == 2 ||
        throw(ArgumentError("unexpected internal PWFS camera links"))
    raw_frame = AOG.AlgorithmPortReference(layout, :frame)
    any(link -> link.source == raw_frame, external_links) &&
        throw(ArgumentError("the raw PWFS frame must terminate at the HIL boundary"))
    links = map(link -> AOG.link(redirect_output(link.source), redirect_input(link.destination)), external_links)
    inputs = map(definition.inputs) do input
        AOG.graph_input(AOG._boundary_name(input), redirect_input(input.destination), input.values)
    end
    outputs = map(definition.outputs) do output
        name = AOG._boundary_name(output)
        output.source == raw_frame && name !== :frame &&
            throw(ArgumentError("PWFS camera loans require one raw frame export named frame"))
        AOG.graph_output(name, redirect_output(output.source))
    end
    return AOG.algorithm_graph(nodes; definition.name, inputs, outputs, links,
        definition.delayed_links, definition.parameters)
end

function AOG.graph_node_ports(::Type{PWFSLoanCamera}, config)
    pyramid, detector, layout = config.camera_nodes
    opd, rate = AOG.graph_node_ports(AOG.PyramidRateNode{Float32}, pyramid.config)
    _, electrons = AOG.graph_node_ports(AOG.EMCCDDetectorAcquisitionNode{Float32}, detector.config)
    _, frame = AOG.graph_node_ports(AOG.DetectorFrameLayoutNode, layout.config)
    detector_frame = AOG.graph_port_contract(:detector_frame, :output, :data, Float32,
        electrons.format.shape, electrons.format.schema, electrons.format.layout)
    return (opd, rate, detector_frame, frame)
end

mutable struct PWFSLoanState
    slot::Int
    published::Int
    metadata::Union{Nothing,FGA.SampleMetadata}
    destination::Union{Nothing,Matrix{UInt16}}
end

struct PreparedPWFSLoanCamera{C}
    camera::C
    exposure_ns::Int
    state::PWFSLoanState
end

function AOG.prepare_graph_node(::Type{PWFSLoanCamera}, config, ::NamedTuple{()},
        inputs::NamedTuple{(:opd,)}, outputs::NamedTuple, ::NamedTuple{()}, ::B.HostComputeDevice)
    recipe = pwfs_camera_description(config.camera_nodes)
    graph = FG.prepare_graph(recipe.document; algorithms=recipe.algorithms)
    # The compiled camera order is Pyramid, EMCCD, layout. AOS diagnostics
    # use the same two buffers that connect those Algorithms.
    buffers = ((outputs.photon_rate,), (outputs.detector_frame,), last(graph.node_outputs))
    slots = loan_camera_slots(config.frame_slots, outputs.frame)
    camera = FG._prepare_cpu_loans(graph, ((; inputs.opd),), slots, buffers)
    exposure_ns = round(Int, Float64(config.camera_nodes[2].config.exposure_duration_s)*1e9)
    return PreparedPWFSLoanCamera(camera, exposure_ns, PWFSLoanState(0, 0, nothing, nothing))
end

loan_camera_slots(slots::Tuple, frame) = map(slot -> (; frame=slot), slots)
loan_camera_slots(::Nothing, frame) = ((; frame),)

# Keep the nested camera executor separate from the enclosing AOS schedule.
@noinline function AOG.step_graph_node!(owner::PreparedPWFSLoanCamera)
    state = owner.state
    metadata = (; opd=state.metadata)
    destination = state.destination
    result = if destination === nothing
        state.slot > 0 || throw(AOG.AlgorithmGraphError("PWFS camera requires an admitted HIL frame slot"))
        FG.process_cpu_loan!(owner.camera, 1, state.slot, metadata)
    else
        FG._process_cpu_output_loan!(owner.camera, PermutedDimsArray(destination, (2, 1)), metadata)
    end
    result === (; frame=true) || throw(AOG.AlgorithmGraphError("PWFS camera rejected the frame"))
    return nothing
end

function AOG.reset_graph_node!(owner::PreparedPWFSLoanCamera)
    FG.reset!(owner.camera)
    owner.state.slot = 0
    owner.state.published = 0
    owner.state.metadata = nothing
    owner.state.destination = nothing
    return nothing
end

struct PreparedPWFSLoanHIL{H,C,S}
    boundary::H
    camera::C
    frame_slots::S
end

struct PreparedPWFSBorrowedHIL{H,C} <: AOG.PreparedBorrowedFrameHILBoundary
    boundary::H
    camera::C
end

const PWFSHIL = Union{PreparedPWFSLoanHIL,PreparedPWFSBorrowedHIL}

function O.prepare_filtergraphs_pwfs_hil(definition::AOG.AlgorithmGraphDefinition;
        command_input, frame_slots::Union{Tuple,Nothing}=nothing)
    frame_slots === nothing && !isdefined(FG, :_process_cpu_output_loan!) &&
        throw(ArgumentError("borrowed camera frames require transient output loan support in FilterGraphs"))
    camera_nodes = pwfs_camera_nodes(definition)
    graph = AOG.prepare_algorithm_graph(loan_camera_definition(definition, camera_nodes, frame_slots))
    camera = AOG.prepared_graph_node(graph, AOG._node_name(last(camera_nodes)))
    boundary = AOG.prepare_graph_hil_boundary(graph; command_input, frame_output=:frame)
    frame_slots === nothing && return PreparedPWFSBorrowedHIL(boundary, camera)
    # Check external storage that lives outside the camera subgraph as well.
    for input in values(graph.inputs), frame in frame_slots
        input isa Array && FG._cpu_loan_overlap(input, frame) &&
            throw(ArgumentError("a PWFS output slot overlaps a simulator input"))
    end
    return PreparedPWFSLoanHIL(boundary, camera, frame_slots)
end

AOG.hil_command_buffer(owner::PWFSHIL) = AOG.hil_command_buffer(owner.boundary)
AOG.hil_boundary_status(owner::PWFSHIL) = AOG.hil_boundary_status(owner.boundary)
AOG.adopt_hil_command!(owner::PWFSHIL, sequence) =
    AOG.adopt_hil_command!(owner.boundary, sequence)

function AOG.hil_frame_buffer(owner::PreparedPWFSLoanHIL)
    slot = owner.camera.state.published
    return slot == 0 ? AOG.hil_frame_buffer(owner.boundary) : owner.frame_slots[slot]
end

AOG.hil_frame_buffer(owner::PreparedPWFSBorrowedHIL) = AOG.hil_frame_buffer(owner.boundary)

AOG.hil_frame_exposure_nanoseconds(owner::PreparedPWFSBorrowedHIL) = UInt64(owner.camera.exposure_ns)

function AOG.graph_output(owner::PreparedPWFSLoanHIL, name::Symbol)
    name === :frame && return AOG.hil_frame_buffer(owner)
    return AOG.graph_output(owner.boundary.graph, name)
end

function AOG.graph_output(owner::PreparedPWFSBorrowedHIL, name::Symbol)
    name === :frame && throw(ArgumentError("the raw frame belongs to the caller's output loan"))
    return AOG.graph_output(owner.boundary.graph, name)
end

function AOG.step_hil_frame_at!(owner::PreparedPWFSLoanHIL, driver::AOG._ModelTimeDriver, slot::Int)
    1 <= slot <= length(owner.frame_slots) || throw(BoundsError(owner.frame_slots, slot))
    return step_pwfs_hil!(owner, driver, slot, nothing)
end

function AOG.step_hil_frame_at!(owner::PreparedPWFSBorrowedHIL, driver::AOG._ModelTimeDriver,
        destination::Matrix{UInt16})
    size(destination) == reverse(size(AOG.hil_frame_buffer(owner))) ||
        throw(DimensionMismatch("the row-major frame destination has the wrong storage shape"))
    return step_pwfs_hil!(owner, driver, 0, destination)
end

function step_pwfs_hil!(owner, driver, slot, destination)
    boundary = owner.boundary
    AOG._require_hil_frame_step(boundary)
    model_time = AOG.model_time_sequence(driver)
    model_time == boundary.state.frame_sequence || throw(
        AOG.AlgorithmGraphError("the model-time driver and HIL boundary frame sequences are not aligned"),
    )
    state = owner.camera.state
    try
        next = AOG.next_model_timestamp(driver)
        state.metadata = FGA.SampleMetadata(model_time+UInt64(1); timestamp_nanoseconds=AOG.model_nanoseconds(next),
            exposure_duration_nanoseconds=owner.camera.exposure_ns, terminal=true)
        state.slot = slot
        state.destination = destination
        timestamp = AOG.step_graph_at!(boundary.graph, driver)
        sequence = AOG._publish_completed_hil_frame!(boundary)
        state.published = slot
        return (; sequence, timestamp)
    catch
        boundary.state.failed = true
        rethrow()
    finally
        state.slot = 0
        state.metadata = nothing
        state.destination = nothing
    end
end

function AOG.reset_hil_boundary!(owner::PWFSHIL, driver::AOG._ModelTimeDriver)
    AOG.reset_hil_boundary!(owner.boundary, driver)
    return owner
end

function gpi2_ports(definition, direction)
    ports = AOG.graph_node_ports(AOG._node_type(definition), definition.config)
    return filter(port -> AOG._port_direction(port) === direction &&
        AOG._port_role(port) === :data, ports)
end

gpi2_port_values(ports, arrays) = NamedTuple{map(AOG._port_name, ports)}(arrays)

function gpi2_parameters(parameters, ::B.HostComputeDevice)
    return parameters
end

function gpi2_parameters(parameters, target::B.AcceleratorComputeDevice)
    return map(parameters) do values
        storage = B.allocate_device_array(target, eltype(values), size(values)...)
        copyto!(storage, values)
    end
end

function prepare_gpi2_node(definition, inputs::Tuple, outputs::Tuple, parameters=(;))
    target = B.compute_device(first(outputs))
    return AOG.prepare_graph_node(AOG._node_type(definition), definition.config,
        definition.props, gpi2_port_values(gpi2_ports(definition, :input), inputs),
        gpi2_port_values(gpi2_ports(definition, :output), outputs),
        gpi2_parameters(parameters, target), target)
end

struct CapturedGPI2Node{O}
    owner::O
end

FGA.preflight(node::CapturedGPI2Node) = AOG._preflight_captured_graph_node!(node.owner)
FGA.preserves_inputs(::CapturedGPI2Node) = true
FGA.enqueue!(node::CapturedGPI2Node) = AOG.enqueue_captured_graph_node!(node.owner)
FGA.complete!(node::CapturedGPI2Node, ::Int32) = AOG._complete_captured_graph_node!(node.owner)
FGA.reset!(node::CapturedGPI2Node) = AOG.reset_graph_node!(node.owner)

gpi2_algorithm(::Type{AOG.PyramidRateNode{Float32}}) = PyramidRateF32
gpi2_algorithm(::Type{AOG.EMCCDDetectorAcquisitionNode{Float32}}) = EMCCDF32
gpi2_algorithm(::Type{AOG.DetectorFrameLayoutNode}) = DetectorFrameLayoutF32

function gpi2_algorithm(type)
    algorithm = gpi2_science_relay_algorithm(type)
    algorithm === nothing && throw(ArgumentError("unsupported GPI2 plant node $type"))
    return algorithm
end

gpi2_node_configuration(node, ::NamedTuple{()}) = (; definition=node)
gpi2_node_configuration(node, parameters) = (; definition=node, parameters)

function gpi2_node_configuration(node::AOG.AlgorithmNodeDefinition{
        Name,AOG.PupilOPDCompositionNode{Float32}}, parameters::NamedTuple{()}) where {Name}
    return (; definition=node, parameters)
end

function gpi2_parameter_bindings(node, bindings)
    name = AOG._node_name(node)
    contracts = AOG._parameter_contracts(values(node.ports))
    parameters = AOG._bind_node_parameters(name, contracts, bindings)
    return NamedTuple{map(AOG._port_name, contracts)}(parameters)
end

gpi2_endpoint(port) = string(AOG._node_name(port), ':', AOG._port_name(port))

function O.gpi2_plant_filter_graph(definition::AOG.AlgorithmGraphDefinition)
    isempty(definition.delayed_links) ||
        throw(ArgumentError("the GPI2 plant recipe requires direct links"))
    admitted = AOG._admitted_graph_nodes(definition)
    bindings = map(definition.parameters) do binding
        AOG.SparseParameterBinding(binding.endpoint, Array(binding.values))
    end
    AOG._validate_parameters!(admitted, bindings, B.HostComputeDevice())
    algorithms = map(node -> gpi2_algorithm(AOG._node_type(node)), definition.nodes)
    nodes = map(definition.nodes, algorithms, values(admitted)) do node, algorithm, contract
        (; name=AOG._node_name(node), label=FGA.algorithm_label(algorithm),
            config=gpi2_node_configuration(node, gpi2_parameter_bindings(contract, bindings)))
    end
    links = map(link -> gpi2_endpoint(link.source) => gpi2_endpoint(link.destination),
        definition.links)
    inputs = NamedTuple{map(AOG._boundary_name, definition.inputs)}(
        map(input -> gpi2_endpoint(input.destination), definition.inputs))
    outputs = NamedTuple{map(AOG._boundary_name, definition.outputs)}(
        map(output -> gpi2_endpoint(output.source), definition.outputs))
    return (; document=(; nodes, links, inputs, outputs), algorithms=Tuple(unique(algorithms)))
end

struct PreparedGPI2HIL{G,I,O,C,F}
    graph::G
    inputs::I
    outputs::O
    active_commands::C
    command_buffer::C
    initial_commands::C
    frame_buffer::F
    exposure_ns::Int
    state::AOG.GraphHILBoundaryState
end

function gpi2_outputs(graph::FG.PreparedGraph)
    outputs = map(graph.output_routes) do route
        FG._route_value(route, (;), graph.node_outputs)
    end
    return NamedTuple{keys(graph.output_formats)}(outputs)
end

gpi2_outputs(graph::FG.PreparedCapturedGraph) = graph.device_outputs

function O.prepare_filtergraphs_gpi2_hil(definition::AOG.AlgorithmGraphDefinition;
        command_input::Tuple, target=B.HostComputeDevice(), exposure_ns=nothing)
    _, detector, _ = pwfs_camera_nodes(definition)
    duration = round(Int, Float64(detector.config.exposure_duration_s)*1e9)
    exposure_ns === nothing || exposure_ns == duration ||
        throw(ArgumentError("the HIL exposure does not match the prepared EMCCD exposure"))
    !isempty(command_input) && allunique(command_input) ||
        throw(ArgumentError("GPI2 command names must be nonempty and distinct"))
    inputs = NamedTuple{map(AOG._boundary_name, definition.inputs)}(
        map(input -> Array(input.values), definition.inputs))
    recipe = O.gpi2_plant_filter_graph(definition)
    checked = FG.prepare_graph(recipe.document; algorithms=recipe.algorithms)
    inputs = FG._validate_external_arrays(:input, checked.input_formats, inputs)
    active = NamedTuple{command_input}(map(name -> getproperty(inputs, name), command_input))
    all(values(active)) do command
        eltype(command) === Float32
    end || throw(ArgumentError("GPI2 HIL commands require Float32 arrays"))
    AOG._validate_hil_command(active)
    graph = prepare_execution(checked, target)
    outputs = gpi2_outputs(graph)
    frame = zeros(UInt16, size(outputs.frame))
    return PreparedGPI2HIL(graph, inputs, outputs, active, map(copy, active),
        map(copy, active), frame, duration, AOG.GraphHILBoundaryState())
end

AOG.hil_frame_buffer(owner::PreparedGPI2HIL) = owner.frame_buffer
AOG.hil_command_buffer(owner::PreparedGPI2HIL) = owner.command_buffer
AOG.graph_input(owner::PreparedGPI2HIL, name::Symbol) = getproperty(owner.inputs, name)

function AOG.graph_output(owner::PreparedGPI2HIL, name::Symbol)
    return name === :frame ? owner.frame_buffer : getproperty(owner.outputs, name)
end

function AOG.hil_boundary_status(owner::PreparedGPI2HIL)
    state = owner.state
    return (; state.frame_sequence, state.active_command_sequence,
        state.command_response_required, state.failed)
end

function run_gpi2_plant!(owner::PreparedGPI2HIL{<:FG.PreparedGraph}, metadata)
    graph = owner.graph
    admitted = FG._validate_external_metadata(graph.input_formats, metadata)
    result = FG._process_candidate_admitted!(nothing, graph, owner.inputs, admitted, graph.result_types)
    result isa NamedTuple && result.terminal && all(values(result.availability)) ||
        throw(AOG.AlgorithmGraphError("the GPI2 plant did not complete a frame"))
    copyto!(owner.frame_buffer, owner.outputs.frame)
    return nothing
end

function run_gpi2_plant!(owner::PreparedGPI2HIL{<:FG.PreparedCapturedGraph}, metadata)
    result = FG.process_device!(owner.graph, owner.inputs, metadata) do outputs
        copyto!(owner.frame_buffer, outputs.frame)
    end
    result isa NamedTuple && all(values(result)) ||
        throw(AOG.AlgorithmGraphError("the captured GPI2 plant did not complete a frame"))
    return nothing
end

function AOG.step_hil_frame_at!(owner::PreparedGPI2HIL, driver::AOG._ModelTimeDriver)
    state = owner.state
    state.failed && AOG._throw_hil_boundary_failed()
    state.command_response_required && AOG._throw_hil_command_response_required(state.frame_sequence)
    AOG.model_time_sequence(driver) == state.frame_sequence ||
        throw(AOG.AlgorithmGraphError("the model-time driver and HIL boundary frame sequences are not aligned"))
    try
        timestamp = AOG.next_model_timestamp(driver)
        sequence = state.frame_sequence + UInt64(1)
        sample = FGA.SampleMetadata(sequence; timestamp_nanoseconds=AOG.model_nanoseconds(timestamp),
            exposure_duration_nanoseconds=owner.exposure_ns, terminal=true)
        metadata = map(_ -> sample, owner.inputs)
        run_gpi2_plant!(owner, metadata)
        AOG._commit_model_time!(driver)
        state.frame_sequence = sequence
        state.command_response_required = true
        return (; sequence, timestamp)
    catch
        state.failed = true
        rethrow()
    end
end

function AOG.adopt_hil_command!(owner::PreparedGPI2HIL, sequence::UInt64)
    state = owner.state
    state.failed && AOG._throw_hil_boundary_failed()
    AOG._require_hil_command_sequence(state, sequence)
    AOG._validate_hil_command_shape(owner.command_buffer, owner.active_commands)
    AOG._validate_hil_command(owner.command_buffer)
    try
        AOG._restore_hil_host_commands!(owner.active_commands, owner.command_buffer)
    catch
        state.failed = true
        rethrow()
    end
    return AOG._commit_hil_command!(state, sequence)
end

AOG.adopt_hil_command!(::PreparedGPI2HIL, sequence) =
    throw(AOG.AlgorithmGraphError("HIL command sequence must be UInt64, not $(typeof(sequence))"))

function AOG.reset_hil_boundary!(owner::PreparedGPI2HIL)
    state = owner.state
    state.failed = true
    AOG._validate_hil_command_shape(owner.command_buffer, owner.initial_commands)
    FG.reset!(owner.graph)
    AOG._restore_hil_host_commands!(owner.command_buffer, owner.initial_commands)
    AOG._restore_hil_host_commands!(owner.active_commands, owner.initial_commands)
    fill!(owner.frame_buffer, 0)
    state.frame_sequence = UInt64(0)
    state.active_command_sequence = UInt64(0)
    state.command_response_required = false
    state.failed = false
    return owner
end

function AOG.reset_hil_boundary!(owner::PreparedGPI2HIL, driver::AOG._ModelTimeDriver)
    AOG.reset_hil_boundary!(owner)
    AOG.reset_model_time!(driver)
    return owner
end

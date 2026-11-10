# Keep the complex fields and Proper FFT plans inside the train. The real-valued
# graph boundary exposes the pupil, the applied DM surface and raw intensity.
const OPTICAL_TRAIN_INPUTS = (:pupil_opd, :pupil_amplitude, :dm_opd)

function optical_train_port(definition, name)
    return only(port for port in AOG.graph_node_ports(AOG._node_type(definition), definition.config)
        if AOG._port_name(port) === name).format
end

function prepare_optical_train(definition, inputs, output, target)
    return AOG.prepare_graph_node(AOG._node_type(definition), definition.config,
        definition.props, inputs, (; output), NamedTuple(), target)
end

struct TrainWorkspace{O,I,A}
    owner::O
    inputs::I
    output::A
end

FGA.commit!(::TrainWorkspace) = nothing
FGA.discard!(::TrainWorkspace) = nothing
FGA.reset!(workspace::TrainWorkspace) = AOG.reset_graph_node!(workspace.owner)

function prepare_host_train(definition)
    proper = proper_extension()
    AOG._node_type(definition) === proper.OpticalTrainNode{Float32} ||
        throw(ArgumentError("optical train requires a scalar Float32 Proper definition"))
    ports = AOG.graph_node_ports(AOG._node_type(definition), definition.config)
    Tuple(AOG._port_name(port) for port in ports if AOG._port_direction(port) === :input) ==
        OPTICAL_TRAIN_INPUTS || throw(ArgumentError("optical train requires pupil and DM OPD inputs"))
    definition.config.output_quantity === Val(:intensity) ||
        throw(ArgumentError("optical train requires real Float32 intensity output"))

    target = B.HostComputeDevice()
    inputs = NamedTuple{OPTICAL_TRAIN_INPUTS}(map(OPTICAL_TRAIN_INPUTS) do name
        zeros(Float32, optical_train_port(definition, name).shape)
    end)
    output = zeros(Float32, optical_train_port(definition, :output).shape)
    owner = prepare_optical_train(definition, inputs, output, target)
    return (; plan=definition, workspace=TrainWorkspace(owner, inputs, output))
end

function process_optical_train!(output, workspace::TrainWorkspace, inputs)
    foreach(copyto!, values(workspace.inputs), inputs)
    AOG.step_graph_node!(workspace.owner)
    copyto!(output, workspace.output)
    return nothing
end

@algorithm OpticalTrainF32 begin
    label = "aos-optical-train-f32"

    @configuration OpticalTrainConfiguration begin
        definition::AOG.AlgorithmNodeDefinition => graph_rebuild
    end

    ports(c) = (
        input(:pupil_opd, Float32, optical_train_port(c.definition, :pupil_opd).shape,
            optical_train_port(c.definition, :pupil_opd).schema),
        input(:pupil_amplitude, Float32, optical_train_port(c.definition, :pupil_amplitude).shape,
            optical_train_port(c.definition, :pupil_amplitude).schema),
        input(:dm_opd, Float32, optical_train_port(c.definition, :dm_opd).shape,
            optical_train_port(c.definition, :dm_opd).schema),
        output(:output, Float32, optical_train_port(c.definition, :output).shape,
            optical_train_port(c.definition, :output).schema; metadata_from=:pupil_opd),
    )

    prepare(c) = prepare_host_train(c.definition)
    process!(plan, workspace) = process_optical_train!(
        output, workspace, (pupil_opd, pupil_amplitude, dm_opd))
end

struct CapturedOpticalTrain{O}
    owner::O
end

FGA.capture_support(::Type{OpticalTrainF32}, device) = optical_capture_support(device)

function FGA.prepare_capture(::FGA.CompleteTransition, device,
        prepared::FGA.PreparedAlgorithm{OpticalTrainF32},
        outputs::Tuple, inputs::Tuple, status, index::Integer)
    output = only(outputs)
    target = B.compute_device(output)
    owner = prepare_optical_train(prepared.plan, NamedTuple{OPTICAL_TRAIN_INPUTS}(inputs),
        output, target)
    AOG.graph_node_capture_capability(owner) isa AOG.GraphNodeCaptureSafe ||
        throw(ArgumentError("optical train does not support native capture on the selected device"))
    return CapturedOpticalTrain(owner)
end

FGA.preflight(::CapturedOpticalTrain) = nothing
FGA.preserves_inputs(::CapturedOpticalTrain) = true
FGA.enqueue!(node::CapturedOpticalTrain) = AOG.enqueue_captured_graph_node!(node.owner)
FGA.complete!(::CapturedOpticalTrain, ::Int32) = nothing
FGA.reset!(node::CapturedOpticalTrain) = AOG.reset_graph_node!(node.owner)

struct AtmosphereWorkspace{O,A}
    owner::O
    output::A
end

FGA.commit!(::AtmosphereWorkspace) = nothing
FGA.discard!(::AtmosphereWorkspace) = nothing
FGA.reset!(workspace::AtmosphereWorkspace) = AOG.reset_graph_node!(workspace.owner)

function gpi2_atmosphere_format(definition)
    AOG._node_type(definition) === AOG.MultiLayerAtmosphereOPDNode{Float32} ||
        throw(ArgumentError("GPI2 atmosphere requires a Float32 multilayer atmosphere"))
    return only(gpi2_ports(definition, :output)).format
end

function prepare_host_gpi2_atmosphere(definition)
    output = zeros(Float32, gpi2_atmosphere_format(definition).shape)
    return prepare_host_gpi2_atmosphere(definition, output)
end

function prepare_host_gpi2_atmosphere(definition, output)
    owner = prepare_gpi2_node(definition, (), (output,))
    return (; plan=definition, workspace=AtmosphereWorkspace(owner, output))
end

function process_gpi2_atmosphere!(output, workspace::AtmosphereWorkspace)
    AOG.step_graph_node!(workspace.owner)
    output === workspace.output || copyto!(output, workspace.output)
    return nothing
end

@algorithm MultiLayerAtmosphereOPDF32 begin
    label = "aos-multilayer-atmosphere-opd-f32"

    @configuration MultiLayerAtmosphereOPDConfiguration begin
        definition::AOG.AlgorithmNodeDefinition => graph_rebuild
    end

    ports(c) = (output(:atmosphere_opd, Float32, gpi2_atmosphere_format(c.definition).shape,
        gpi2_atmosphere_format(c.definition).schema),)

    prepare(c) = prepare_host_gpi2_atmosphere(c.definition)
    process!(plan, workspace) = process_gpi2_atmosphere!(atmosphere_opd, workspace)
end

gpi2_algorithm(::Type{AOG.MultiLayerAtmosphereOPDNode{Float32}}) = MultiLayerAtmosphereOPDF32

function gpi2_science_relay_algorithm(type)
    proper = proper_extension()
    proper === nothing && return nothing
    if type === proper.OpticalTrainNode{Float32} ||
            type === proper.PolychromaticOpticalTrainNode{Float32}
        return ScienceRelayF32
    end
    return nothing
end

function prepare_host_gpi2_science_relay(definition)
    gpi2_science_relay_algorithm(AOG._node_type(definition)) === ScienceRelayF32 ||
        throw(ArgumentError("GPI2 science relay requires a Float32 Proper optical train"))
    inputs = gpi2_ports(definition, :input)
    outputs = gpi2_ports(definition, :output)
    map(AOG._port_name, inputs) == (:pupil_opd, :pupil_amplitude) &&
        map(AOG._port_name, outputs) == (:output,) &&
        only(outputs).format isa AOG._GraphPortFormat{Float32} ||
        throw(ArgumentError("GPI2 science relay requires pupil OPD and amplitude and a real intensity output"))
    pupil = map(port -> zeros(Float32, port.format.shape), inputs)
    output = zeros(Float32, only(outputs).format.shape)
    return prepare_host_gpi2_science_relay(definition, pupil, output)
end

function prepare_host_gpi2_science_relay(definition, inputs::Tuple, output)
    owner = prepare_gpi2_node(definition, inputs, (output,))
    return (; plan=definition, workspace=TrainWorkspace(owner, inputs, output))
end

function process_gpi2_science_relay!(output, workspace::TrainWorkspace, opd, amplitude)
    bound_opd, bound_amplitude = workspace.inputs
    opd === bound_opd || copyto!(bound_opd, opd)
    amplitude === bound_amplitude || copyto!(bound_amplitude, amplitude)
    AOG.step_graph_node!(workspace.owner)
    output === workspace.output || copyto!(output, workspace.output)
    return nothing
end

@algorithm ScienceRelayF32 begin
    label = "aos-science-relay-f32"

    @configuration ScienceRelayConfiguration begin
        definition::AOG.AlgorithmNodeDefinition => graph_rebuild
    end

    ports(c) = (
        input(:pupil_opd, Float32, optical_train_port(c.definition, :pupil_opd).shape,
            optical_train_port(c.definition, :pupil_opd).schema),
        input(:pupil_amplitude, Float32, optical_train_port(c.definition, :pupil_amplitude).shape,
            optical_train_port(c.definition, :pupil_amplitude).schema),
        output(:output, Float32, optical_train_port(c.definition, :output).shape,
            optical_train_port(c.definition, :output).schema; metadata_from=:pupil_opd),
    )

    prepare(c) = prepare_host_gpi2_science_relay(c.definition)
    process!(plan, workspace) = process_gpi2_science_relay!(
        output, workspace, pupil_opd, pupil_amplitude)
end

const EnvironmentAlgorithm = Union{MultiLayerAtmosphereOPDF32,ScienceRelayF32}

FGA.capture_support(::Type{<:EnvironmentAlgorithm}, device) = optical_capture_support(device)

function FGA.prepare_capture(::FGA.CompleteTransition, device,
        prepared::FGA.PreparedAlgorithm{D}, outputs::Tuple, inputs::Tuple,
        status, index::Integer) where {D<:EnvironmentAlgorithm}
    owner = prepare_gpi2_node(prepared.plan, inputs, outputs)
    AOG.graph_node_capture_capability(owner) isa AOG.GraphNodeCaptureSafe ||
        throw(ArgumentError("GPI2 environment node does not support native capture on the selected device"))
    return CapturedGPI2Node(owner)
end

gpi2_opd_format(definition, direction, index=1) =
    gpi2_ports(definition, direction)[index].format

struct GPI2OPDWorkspace{O}
    owner::O
end

FGA.commit!(::GPI2OPDWorkspace) = nothing
FGA.discard!(::GPI2OPDWorkspace) = nothing
FGA.reset!(workspace::GPI2OPDWorkspace) = AOG.reset_graph_node!(workspace.owner)

function prepare_host_gpi2_opd(configuration, parameters=configuration.parameters)
    definition = configuration.definition
    inputs = map(port -> zeros(Float32, port.format.shape), gpi2_ports(definition, :input))
    outputs = map(port -> zeros(Float32, port.format.shape), gpi2_ports(definition, :output))
    owner = prepare_gpi2_node(definition, inputs, outputs, parameters)
    return (; plan=(; definition, parameters), workspace=GPI2OPDWorkspace(owner))
end

function process_grid_gaussian_mirror!(surface, workspace, command)
    owner = workspace.owner
    bound = AOG._GridGaussianDeformableMirrorSurfaceOwner(
        owner.deformable_mirror, command, owner.actuator_grid_indices, surface,
    )
    return AOG.step_graph_node!(bound)
end

@algorithm GridGaussianMirrorF32 begin
    label = "aos-grid-gaussian-deformable-mirror-surface-f32"

    @configuration GridGaussianMirrorConfiguration begin
        definition::AOG.AlgorithmNodeDefinition => graph_rebuild
        parameters::NamedTuple => graph_rebuild
    end

    ports(c) = (
        input(:pdm_command, Float32, gpi2_opd_format(c.definition, :input).shape,
            gpi2_opd_format(c.definition, :input).schema),
        output(:surface_opd, Float32, gpi2_opd_format(c.definition, :output).shape,
            gpi2_opd_format(c.definition, :output).schema; metadata_from=:pdm_command),
    )

    prepare(c) = prepare_host_gpi2_opd(c)
    process!(plan, workspace) = process_grid_gaussian_mirror!(surface_opd, workspace, pdm_command)
end

function process_second_order_response!(output, workspace, input)
    response = workspace.owner.prepared
    bound = O.PreparedSecondOrderResponse(response.plan, response.state, input, output)
    return AOG.step_graph_node!(AOG._SecondOrderResponseOwner(bound))
end

@algorithm SecondOrderResponseF32 begin
    label = "aos-second-order-response-f32"

    @configuration SecondOrderResponseConfiguration begin
        definition::AOG.AlgorithmNodeDefinition => graph_rebuild
    end

    ports(c) = (
        input(:input, Float32, gpi2_opd_format(c.definition, :input).shape,
            gpi2_opd_format(c.definition, :input).schema),
        output(:output, Float32, gpi2_opd_format(c.definition, :output).shape,
            gpi2_opd_format(c.definition, :output).schema; metadata_from=:input),
    )

    prepare(c) = prepare_host_gpi2_opd(c, (;))
    process!(plan, workspace) = process_second_order_response!(output, workspace, input)
end

@algorithm ModalOPDExpansionF32 begin
    label = "aos-modal-opd-expansion-f32"

    @configuration ModalOPDExpansionConfiguration begin
        definition::AOG.AlgorithmNodeDefinition => graph_rebuild
        parameters::NamedTuple => graph_rebuild
    end

    ports(c) = (
        input(:coefficients, Float32, gpi2_opd_format(c.definition, :input).shape,
            gpi2_opd_format(c.definition, :input).schema),
        output(:opd, Float32, gpi2_opd_format(c.definition, :output).shape,
            gpi2_opd_format(c.definition, :output).schema; metadata_from=:coefficients),
    )

    prepare(c) = prepare_host_gpi2_opd(c)
    process!(plan, workspace) = AOG.step_graph_node!(
        AOG._ModalOPDExpansionOwner(workspace.owner.plan, coefficients, opd))
end

@algorithm PupilOPDCompositionF32 begin
    label = "aos-pupil-opd-composition-f32"

    @configuration PupilOPDCompositionConfiguration begin
        definition::AOG.AlgorithmNodeDefinition => graph_rebuild
        parameters::NamedTuple => graph_rebuild
    end

    ports(c) = (
        input(:uncompensated_opd, Float32, gpi2_opd_format(c.definition, :input, 1).shape,
            gpi2_opd_format(c.definition, :input, 1).schema),
        input(:surface_opd, Float32, gpi2_opd_format(c.definition, :input, 2).shape,
            gpi2_opd_format(c.definition, :input, 2).schema),
        output(:pupil_opd, Float32, gpi2_opd_format(c.definition, :output).shape,
            gpi2_opd_format(c.definition, :output).schema; metadata_from=:surface_opd),
    )

    prepare(c) = c
    process!(plan) = AOG.step_graph_node!(
        AOG._PupilOPDCompositionOwner(uncompensated_opd, surface_opd, pupil_opd))
end

gpi2_algorithm(::Type{AOG.GridGaussianDeformableMirrorSurfaceNode{Float32}}) = GridGaussianMirrorF32
gpi2_algorithm(::Type{AOG.SecondOrderResponseNode{Float32}}) = SecondOrderResponseF32
gpi2_algorithm(::Type{AOG.ModalOPDExpansionNode{Float32}}) = ModalOPDExpansionF32
gpi2_algorithm(::Type{AOG.PupilOPDCompositionNode{Float32}}) = PupilOPDCompositionF32

const GPI2OPDAlgorithm = Union{
    GridGaussianMirrorF32,SecondOrderResponseF32,ModalOPDExpansionF32,PupilOPDCompositionF32,
}
FGA.capture_support(::Type{<:GPI2OPDAlgorithm}, device) = optical_capture_support(device)

function FGA.prepare_capture(::FGA.CompleteTransition, device,
        prepared::FGA.PreparedAlgorithm{D}, outputs::Tuple, inputs::Tuple,
        status, index::Integer) where {D<:GPI2OPDAlgorithm}
    configuration = prepared.plan
    owner = prepare_gpi2_node(configuration.definition, inputs, outputs, configuration.parameters)
    AOG.graph_node_capture_capability(owner) isa AOG.GraphNodeCaptureSafe ||
        throw(ArgumentError("OPD stage does not support native capture on the selected device"))
    return CapturedGPI2Node(owner)
end

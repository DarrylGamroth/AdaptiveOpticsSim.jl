function O.llowfs_filter_graph(definition::AOG.AlgorithmGraphDefinition)
    length(definition.nodes) == 2 && isempty(definition.delayed_links) &&
        isempty(definition.parameters) && length(definition.links) == 1 &&
        length(definition.outputs) == 1 ||
        throw(ArgumentError("LLOWFS requires one optical train feeding one pixel integration node"))
    pixel_index = findfirst(node -> AOG._node_type(node) === AOG.PixelIntegrationNode{Float32},
        definition.nodes)
    pixel_index === nothing && throw(ArgumentError("LLOWFS requires Float32 pixel integration"))
    train, pixels = definition.nodes[3 - pixel_index], definition.nodes[pixel_index]
    train_name, pixel_name = AOG._node_name(train), AOG._node_name(pixels)
    link = only(definition.links)
    output = only(definition.outputs)
    link.source == AOG.AlgorithmPortReference(train_name, :output) &&
        link.destination == AOG.AlgorithmPortReference(pixel_name, :input) &&
        output.source == AOG.AlgorithmPortReference(pixel_name, :output) &&
        AOG._boundary_name(output) === :intensity ||
        throw(ArgumentError("LLOWFS must export the integrated intensity"))
    input_names = (:pupil_opd, :pupil_amplitude, :dm_opd)
    length(definition.inputs) == length(input_names) && all(input_names) do name
        any(definition.inputs) do input
            AOG._boundary_name(input) === name &&
                input.destination == AOG.AlgorithmPortReference(train_name, name)
        end
    end || throw(ArgumentError("LLOWFS must bind pupil OPD, pupil amplitude, and DM OPD"))
    nodes = (
        (; name=train_name, label=FGA.algorithm_label(OpticalTrainF32), config=(; definition=train)),
        (; name=pixel_name, label=FGA.algorithm_label(PixelIntegrationF32), config=(; definition=pixels)),
    )
    links = ("$train_name:output" => "$pixel_name:input",)
    inputs = NamedTuple{input_names}(map(name -> "$train_name:$name", input_names))
    document = (; nodes, links, inputs, outputs=(; intensity="$pixel_name:output"))
    return (; document, algorithms=(OpticalTrainF32, PixelIntegrationF32))
end

prepare_direct_llowfs(::Nothing, observation, target) = nothing

function prepare_direct_llowfs(definition, observation, target)
    description = O.llowfs_filter_graph(definition)
    pixels = last(description.document.nodes)
    owner = AOG.prepared_graph_node(observation, pixels.name)
    config = pixels.config.definition.config
    owner isa AOG._PixelIntegrationOwner && all(fieldnames(typeof(config))) do name
        getfield(owner.plan.config, name) == getfield(config, name)
    end || throw(ArgumentError("direct LLOWFS pixel integration differs from the AOS observation"))
    graph = FG.prepare_graph(description.document; algorithms=description.algorithms)
    for (name, format) in pairs(FG.input_formats(graph))
        check_llowfs_format(format, AOG.graph_input_contract(definition, name),
            AOG.graph_input(observation, name))
    end
    check_llowfs_format(only(FG.output_formats(graph)), AOG.graph_output_contract(definition, :intensity),
        AOG.graph_output(observation, :intensity))
    return (; graph=prepare_execution(graph, target),
        destinations=(; intensity=AOG.graph_output(observation, :intensity)))
end

function check_llowfs_format(format, contract, buffer)
    (format.element_type, format.shape, format.schema) ==
        (contract.eltype, contract.shape, contract.schema) &&
        (format.element_type, format.shape) == (eltype(buffer), size(buffer)) ||
        throw(ArgumentError("direct LLOWFS format differs from the AOS observation"))
    return nothing
end

run_direct_llowfs!(::Nothing, prepared) = nothing

function run_direct_llowfs!(selection, prepared)
    observation = prepared.graphs.llowfs
    AOG._require_external_pupil_graph_ready(observation)
    # DM updates use the AOS observation stream. Complete them before
    # FilterGraphs copies that input on its capture stream.
    B._synchronize_prepared_device_execution_context!(observation.context)
    inputs = (pupil_opd=prepared.pupil_opd, pupil_amplitude=prepared.pupil_amplitude,
        dm_opd=AOG.graph_input(observation, :dm_opd))
    result = run_observations!(selection.graph, inputs, selection.destinations, observation.context)
    result isa NamedTuple && result.intensity || throw(AOG.AlgorithmGraphError("direct LLOWFS rejected an optical sample"))
    return nothing
end

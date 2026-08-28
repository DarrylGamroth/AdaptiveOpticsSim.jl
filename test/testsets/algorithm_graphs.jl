const AOG = AdaptiveOpticsSim.AlgorithmGraphs

struct GraphTestGainDeclaration end

struct GraphTestGainConfiguration
    extent::Int
end

mutable struct GraphTestGainPrepared{T}
    gain::Vector{T}
end

function AOG._prepare_algorithm_instance(
    ::Type{GraphTestGainDeclaration},
    configuration::GraphTestGainConfiguration,
)
    configuration.extent > 0 || throw(ArgumentError("extent must be positive"))
    return GraphTestGainPrepared(ones(Float32, configuration.extent))
end

function AOG._algorithm_port_contracts(
    ::Type{GraphTestGainDeclaration},
    prepared::GraphTestGainPrepared{T},
) where {T}
    shape = (length(prepared.gain),)
    schema = "test.graph.signal.f32/1"
    return (
        AOG._graph_port_contract(
            :input,
            :input,
            :data,
            T,
            shape,
            schema,
            :column_major,
        ),
        AOG._graph_port_contract(
            :output,
            :output,
            :data,
            T,
            shape,
            schema,
            :column_major,
        ),
        AOG._graph_port_contract(
            :gain,
            :input,
            :parameter,
            T,
            shape,
            "test.graph.gain.f32/1",
            :column_major,
        ),
    )
end

function AOG._process_algorithm!(
    prepared::GraphTestGainPrepared,
    outputs::NamedTuple,
    inputs::NamedTuple,
)
    @inbounds for index in eachindex(outputs.output, inputs.input, prepared.gain)
        outputs.output[index] = inputs.input[index] * prepared.gain[index]
    end
    return nothing
end

AOG._reset_algorithm!(::GraphTestGainPrepared) = nothing

function AOG._replace_algorithm_parameter!(
    prepared::GraphTestGainPrepared,
    ::Val{:gain},
    values,
)
    copyto!(prepared.gain, values)
    return nothing
end

struct GraphTestAddDeclaration end

struct GraphTestAddConfiguration
    extent::Int
    schema::String
    increment::Float32
end

struct GraphTestAddPrepared
    extent::Int
    schema::String
    increment::Float32
end

function AOG._prepare_algorithm_instance(
    ::Type{GraphTestAddDeclaration},
    configuration::GraphTestAddConfiguration,
)
    return GraphTestAddPrepared(
        configuration.extent,
        configuration.schema,
        configuration.increment,
    )
end

function AOG._algorithm_port_contracts(
    ::Type{GraphTestAddDeclaration},
    prepared::GraphTestAddPrepared,
)
    format(name, direction) = AOG._graph_port_contract(
        name,
        direction,
        :data,
        Float32,
        (prepared.extent,),
        prepared.schema,
        :column_major,
    )
    return (format(:input, :input), format(:output, :output))
end

function AOG._process_algorithm!(
    prepared::GraphTestAddPrepared,
    outputs::NamedTuple,
    inputs::NamedTuple,
)
    @inbounds for index in eachindex(outputs.output, inputs.input)
        outputs.output[index] = inputs.input[index] + prepared.increment
    end
    return nothing
end

AOG._reset_algorithm!(::GraphTestAddPrepared) = nothing

struct GraphTestSourceDeclaration end

struct GraphTestSourceConfiguration
    extent::Int
    value::Float32
end

struct GraphTestSourcePrepared
    extent::Int
    value::Float32
end

AOG._prepare_algorithm_instance(
    ::Type{GraphTestSourceDeclaration},
    configuration::GraphTestSourceConfiguration,
) = GraphTestSourcePrepared(configuration.extent, configuration.value)

function AOG._algorithm_port_contracts(
    ::Type{GraphTestSourceDeclaration},
    prepared::GraphTestSourcePrepared,
)
    return (
        AOG._graph_port_contract(
            :output,
            :output,
            :data,
            Float32,
            (prepared.extent,),
            "test.graph.signal.f32/1",
            :column_major,
        ),
    )
end

function AOG._process_algorithm!(
    prepared::GraphTestSourcePrepared,
    outputs::NamedTuple,
    ::NamedTuple{()},
)
    fill!(outputs.output, prepared.value)
    return nothing
end

AOG._reset_algorithm!(::GraphTestSourcePrepared) = nothing

struct GraphTestSinkDeclaration end

struct GraphTestSinkConfiguration
    extent::Int
end

mutable struct GraphTestSinkPrepared
    extent::Int
    sum::Float32
end

AOG._prepare_algorithm_instance(
    ::Type{GraphTestSinkDeclaration},
    configuration::GraphTestSinkConfiguration,
) = GraphTestSinkPrepared(configuration.extent, 0.0f0)

function AOG._algorithm_port_contracts(
    ::Type{GraphTestSinkDeclaration},
    prepared::GraphTestSinkPrepared,
)
    return (
        AOG._graph_port_contract(
            :input,
            :input,
            :data,
            Float32,
            (prepared.extent,),
            "test.graph.signal.f32/1",
            :column_major,
        ),
    )
end


function AOG._process_algorithm!(
    prepared::GraphTestSinkPrepared,
    ::NamedTuple{()},
    inputs::NamedTuple,
)
    total = 0.0f0
    @inbounds for value in inputs.input
        total += value
    end
    prepared.sum = total
    return nothing
end


function AOG._reset_algorithm!(prepared::GraphTestSinkPrepared)
    prepared.sum = 0.0f0
    return nothing
end

struct GraphTestFailureDeclaration end

struct GraphTestFailurePrepared
    extent::Int
end

AOG._prepare_algorithm_instance(
    ::Type{GraphTestFailureDeclaration},
    extent::Int,
) = GraphTestFailurePrepared(extent)

function AOG._algorithm_port_contracts(
    ::Type{GraphTestFailureDeclaration},
    prepared::GraphTestFailurePrepared,
)
    format(name, direction) = AOG._graph_port_contract(
        name,
        direction,
        :data,
        Float32,
        (prepared.extent,),
        "test.graph.signal.f32/1",
        :column_major,
    )
    return (format(:input, :input), format(:output, :output))
end

function AOG._process_algorithm!(
    ::GraphTestFailurePrepared,
    outputs::NamedTuple,
    inputs::NamedTuple,
)
    any(<(0.0f0), inputs.input) && throw(DomainError(
        inputs.input,
        "negative graph fixture input",
    ))
    copyto!(outputs.output, inputs.input)
    return nothing
end

AOG._reset_algorithm!(::GraphTestFailurePrepared) = nothing

@testset "portable algorithm graph direct links and sparse parameters" begin
    input = Float32[1, 2, 3]
    definition = algorithm_graph(
        (
            algorithm_node(
                :gain,
                GraphTestGainDeclaration,
                GraphTestGainConfiguration(3),
            ),
            algorithm_node(
                :sink,
                GraphTestSinkDeclaration,
                GraphTestSinkConfiguration(3),
            ),
        );
        inputs=(graph_input(:residual, :gain => :input, input),),
        outputs=(graph_output(:command, :gain => :output),),
        links=(link(:gain => :output, :sink => :input),),
        parameters=(sparse_parameter(:gain => :gain, Float32[2, 3, 4]),),
    )
    graph = prepare_algorithm_graph(definition)

    @test graph_input(graph, Val(:residual)) === input
    @test graph_input(graph, :residual) === input
    @test prepared_algorithm(graph, Val(:gain)) isa GraphTestGainPrepared
    @test iszero(graph_step_sequence(graph))
    @test all(iszero, graph_output(graph, Val(:command)))
    @test step_graph!(graph) === graph
    @test graph_output(graph, Val(:command)) == Float32[2, 6, 12]
    @test graph_output(graph, :command) === graph_output(graph, Val(:command))
    @test prepared_algorithm(graph, Val(:sink)).sum == 20.0f0
    @test graph_step_sequence(graph) == UInt64(1)
    @test !graph_failed(graph)

    step_graph!(graph)
    @test @allocated(step_graph!(graph)) == 0
    @test @inferred(step_graph!(graph)) === graph
    @test all(isconcretetype, fieldtypes(typeof(graph)))
end

@testset "portable algorithm graph delayed feedback and reset" begin
    configuration = GraphTestAddConfiguration(
        2,
        "test.graph.signal.f32/1",
        1.0f0,
    )
    definition = algorithm_graph(
        (
            algorithm_node(:first, GraphTestAddDeclaration, configuration),
            algorithm_node(:second, GraphTestAddDeclaration, configuration),
        );
        outputs=(graph_output(:result, :second => :output),),
        links=(link(:first => :output, :second => :input),),
        delayed_links=(delayed_link(
            :second => :output,
            :first => :input,
            zeros(Float32, 2),
        ),),
    )
    graph = prepare_algorithm_graph(definition)

    step_graph!(graph)
    @test graph_output(graph, Val(:result)) == Float32[2, 2]
    step_graph!(graph)
    @test graph_output(graph, Val(:result)) == Float32[4, 4]
    @test graph_step_sequence(graph) == UInt64(2)
    @test @allocated(step_graph!(graph)) == 0

    @test reset_graph!(graph) === graph
    @test graph_step_sequence(graph) == UInt64(0)
    @test !graph_failed(graph)
    step_graph!(graph)
    @test graph_output(graph, Val(:result)) == Float32[2, 2]
end

@testset "portable algorithm graph admits pure sources and sinks" begin
    source = prepare_algorithm_graph(algorithm_graph(
        (algorithm_node(
            :source,
            GraphTestSourceDeclaration,
            GraphTestSourceConfiguration(2, 7.0f0),
        ),);
        outputs=(graph_output(:sample, :source => :output),),
    ))
    step_graph!(source)
    @test graph_output(source, Val(:sample)) == Float32[7, 7]

    input = Float32[2, 5]
    sink = prepare_algorithm_graph(algorithm_graph(
        (algorithm_node(
            :sink,
            GraphTestSinkDeclaration,
            GraphTestSinkConfiguration(2),
        ),);
        inputs=(graph_input(:sample, :sink => :input, input),),
    ))
    step_graph!(sink)
    @test isempty(sink.outputs)
    @test prepared_algorithm(sink, Val(:sink)).sum == 7.0f0
end

@testset "portable algorithm graph failure is fail-stop" begin
    input = Float32[1, -1]
    graph = prepare_algorithm_graph(algorithm_graph(
        (algorithm_node(:checked, GraphTestFailureDeclaration, 2),);
        inputs=(graph_input(:input, :checked => :input, input),),
        outputs=(graph_output(:output, :checked => :output),),
    ))

    @test_throws DomainError step_graph!(graph)
    @test graph_failed(graph)
    @test graph_step_sequence(graph) == UInt64(0)
    @test_throws AlgorithmGraphError step_graph!(graph)
    fill!(input, 3.0f0)
    reset_graph!(graph)
    step_graph!(graph)
    @test graph_output(graph, Val(:output)) == input
end

@testset "portable algorithm graph rejects invalid topology" begin
    configuration = GraphTestAddConfiguration(
        2,
        "test.graph.signal.f32/1",
        1.0f0,
    )
    node = algorithm_node(:add, GraphTestAddDeclaration, configuration)

    @test_throws AlgorithmGraphError algorithm_graph(())
    @test_throws AlgorithmGraphError algorithm_graph((node, node))
    @test_throws AlgorithmGraphError prepare_algorithm_graph(algorithm_graph((node,)))
    @test_throws AlgorithmGraphError prepare_algorithm_graph(algorithm_graph(
        (node,);
        inputs=(graph_input(:input, :add => :missing, zeros(Float32, 2)),),
    ))
    @test_throws AlgorithmGraphError prepare_algorithm_graph(algorithm_graph(
        (node,);
        inputs=(graph_input(:input, :add => :input, zeros(Float32, 3)),),
    ))

    backward = algorithm_graph(
        (
            algorithm_node(:first, GraphTestAddDeclaration, configuration),
            algorithm_node(:second, GraphTestAddDeclaration, configuration),
        );
        links=(link(:second => :output, :first => :input),),
        delayed_links=(delayed_link(
            :first => :output,
            :second => :input,
            zeros(Float32, 2),
        ),),
    )
    @test_throws AlgorithmGraphError prepare_algorithm_graph(backward)

    mismatched = GraphTestAddConfiguration(
        2,
        "test.graph.other.f32/1",
        1.0f0,
    )
    schema_mismatch = algorithm_graph(
        (
            algorithm_node(:first, GraphTestAddDeclaration, configuration),
            algorithm_node(:second, GraphTestAddDeclaration, mismatched),
        );
        inputs=(graph_input(:input, :first => :input, zeros(Float32, 2)),),
        links=(link(:first => :output, :second => :input),),
    )
    @test_throws AlgorithmGraphError prepare_algorithm_graph(schema_mismatch)

    duplicate_destination = algorithm_graph(
        (node,);
        inputs=(graph_input(:input, :add => :input, zeros(Float32, 2)),),
        delayed_links=(delayed_link(
            :add => :output,
            :add => :input,
            zeros(Float32, 2),
        ),),
    )
    @test_throws AlgorithmGraphError prepare_algorithm_graph(duplicate_destination)
end

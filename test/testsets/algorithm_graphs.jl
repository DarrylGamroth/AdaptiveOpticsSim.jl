using YAML

const AOG = AdaptiveOpticsSim.AlgorithmGraphs

include(joinpath(
    dirname(dirname(@__DIR__)),
    "examples",
    "support",
    "hil_reference_systems.jl",
))

function warmed_graph_step_allocation_bytes(graph)
    step_graph!(graph)
    return @allocated step_graph!(graph)
end

function warmed_async_graph_step_allocation_bytes(graph)
    wait_graph_step!(step_graph_async!(graph))
    return @allocated wait_graph_step!(step_graph_async!(graph))
end

struct GraphTestGainDeclaration end

struct GraphTestGainConfiguration
    extent::Int
end

struct GraphTestGainOwner{G,I,O}
    gain::G
    input::I
    output::O
end

function AOG.graph_node_ports(
    ::Type{GraphTestGainDeclaration},
    configuration::GraphTestGainConfiguration,
)
    configuration.extent > 0 || throw(ArgumentError("extent must be positive"))
    shape = (configuration.extent,)
    schema = "test.graph.signal.f32/1"
    return (
        AOG.graph_port_contract(
            :input,
            :input,
            :data,
            Float32,
            shape,
            schema,
            :column_major,
        ),
        AOG.graph_port_contract(
            :output,
            :output,
            :data,
            Float32,
            shape,
            schema,
            :column_major,
        ),
        AOG.graph_port_contract(
            :gain,
            :input,
            :parameter,
            Float32,
            shape,
            "test.graph.gain.f32/1",
            :column_major,
        ),
    )
end

function AOG.prepare_graph_node(
    ::Type{GraphTestGainDeclaration},
    ::GraphTestGainConfiguration,
    props,
    inputs::NamedTuple,
    outputs::NamedTuple,
    parameters::NamedTuple,
    target,
)
    return GraphTestGainOwner(parameters.gain, inputs.input, outputs.output)
end

function AOG.step_graph_node!(owner::GraphTestGainOwner)
    @inbounds for index in eachindex(owner.output, owner.input, owner.gain)
        owner.output[index] = owner.input[index] * owner.gain[index]
    end
    return nothing
end

AOG.reset_graph_node!(::GraphTestGainOwner) = nothing

struct GraphTestAddDeclaration end

struct GraphTestAddConfiguration
    extent::Int
    schema::String
    increment::Float32
end

struct GraphTestAddOwner{I,O}
    increment::Float32
    input::I
    output::O
end

function AOG.graph_node_ports(
    ::Type{GraphTestAddDeclaration},
    configuration::GraphTestAddConfiguration,
)
    format(name, direction) = AOG.graph_port_contract(
        name,
        direction,
        :data,
        Float32,
        (configuration.extent,),
        configuration.schema,
        :column_major,
    )
    return (format(:input, :input), format(:output, :output))
end

function AOG.prepare_graph_node(
    ::Type{GraphTestAddDeclaration},
    configuration::GraphTestAddConfiguration,
    props,
    inputs::NamedTuple,
    outputs::NamedTuple,
    parameters::NamedTuple,
    target,
)
    return GraphTestAddOwner(
        configuration.increment,
        inputs.input,
        outputs.output,
    )
end

function AOG.step_graph_node!(owner::GraphTestAddOwner)
    @inbounds for index in eachindex(owner.output, owner.input)
        owner.output[index] = owner.input[index] + owner.increment
    end
    return nothing
end

AOG.reset_graph_node!(::GraphTestAddOwner) = nothing

struct GraphTestSumDeclaration end

struct GraphTestSumConfiguration
    extent::Int
    schema::String
end

struct GraphTestSumOwner{L,R,O}
    left::L
    right::R
    output::O
end

function AOG.graph_node_ports(
    ::Type{GraphTestSumDeclaration},
    configuration::GraphTestSumConfiguration,
)
    format(name, direction) = AOG.graph_port_contract(
        name,
        direction,
        :data,
        Float32,
        (configuration.extent,),
        configuration.schema,
        :column_major,
    )
    return (
        format(:left, :input),
        format(:right, :input),
        format(:output, :output),
    )
end

function AOG.prepare_graph_node(
    ::Type{GraphTestSumDeclaration},
    ::GraphTestSumConfiguration,
    props,
    inputs::NamedTuple,
    outputs::NamedTuple,
    parameters::NamedTuple,
    target,
)
    return GraphTestSumOwner(inputs.left, inputs.right, outputs.output)
end

function AOG.step_graph_node!(owner::GraphTestSumOwner)
    @inbounds for index in eachindex(
        owner.output,
        owner.left,
        owner.right,
    )
        owner.output[index] = owner.left[index] + owner.right[index]
    end
    return nothing
end

AOG.reset_graph_node!(::GraphTestSumOwner) = nothing

struct GraphTestSourceDeclaration end

struct GraphTestSourceConfiguration
    extent::Int
    value::Float32
end

struct GraphTestSourceOwner{O}
    value::Float32
    output::O
end

function AOG.graph_node_ports(
    ::Type{GraphTestSourceDeclaration},
    configuration::GraphTestSourceConfiguration,
)
    return (
        AOG.graph_port_contract(
            :output,
            :output,
            :data,
            Float32,
            (configuration.extent,),
            "test.graph.signal.f32/1",
            :column_major,
        ),
    )
end

function AOG.prepare_graph_node(
    ::Type{GraphTestSourceDeclaration},
    configuration::GraphTestSourceConfiguration,
    props,
    ::NamedTuple{()},
    outputs::NamedTuple,
    parameters::NamedTuple,
    target,
)
    return GraphTestSourceOwner(configuration.value, outputs.output)
end

function AOG.step_graph_node!(owner::GraphTestSourceOwner)
    fill!(owner.output, owner.value)
    return nothing
end

AOG.reset_graph_node!(::GraphTestSourceOwner) = nothing

@testset "graph port contracts reject execution-unsafe declarations" begin
    @test_throws AlgorithmGraphError AOG.graph_port_contract(
        :empty_schema,
        :input,
        :data,
        Float32,
        (1,),
        "",
        :column_major,
    )
    @test_throws AlgorithmGraphError AOG.graph_port_contract(
        :abstract_element,
        :input,
        :data,
        AbstractFloat,
        (1,),
        "test.graph.signal/1",
        :column_major,
    )
    @test_throws AlgorithmGraphError AOG.graph_port_contract(
        :boolean_dimension,
        :input,
        :data,
        Float32,
        (true,),
        "test.graph.signal.f32/1",
        :column_major,
    )
    @test_throws AlgorithmGraphError AOG.graph_port_contract(
        :oversized_dimension,
        :input,
        :data,
        Float32,
        (big(typemax(Int)) + 1,),
        "test.graph.signal.f32/1",
        :column_major,
    )
end

struct GraphTestExactBindingDeclaration end

struct GraphTestExactBindingOwner{I,O}
    input::I
    output::O
end

function AOG.graph_node_ports(
    ::Type{GraphTestExactBindingDeclaration},
    extent::Int,
)
    format(name, direction) = AOG.graph_port_contract(
        name,
        direction,
        :data,
        Float32,
        (extent,),
        "test.graph.signal.f32/1",
        :column_major,
    )
    return (format(:input, :input), format(:output, :output))
end

function AOG.prepare_graph_node(
    ::Type{GraphTestExactBindingDeclaration},
    extent::Int,
    props,
    inputs::NamedTuple,
    outputs::NamedTuple,
    parameters::NamedTuple,
    target,
)
    return GraphTestExactBindingOwner(inputs.input, outputs.output)
end

function AOG.step_graph_node!(owner::GraphTestExactBindingOwner)
    copyto!(owner.output, owner.input)
    return nothing
end

AOG.reset_graph_node!(::GraphTestExactBindingOwner) = nothing

struct GraphTestSinkDeclaration end

struct GraphTestSinkConfiguration
    extent::Int
end

mutable struct GraphTestSinkOwner{I}
    sum::Float32
    input::I
end

function AOG.graph_node_ports(
    ::Type{GraphTestSinkDeclaration},
    configuration::GraphTestSinkConfiguration,
)
    return (
        AOG.graph_port_contract(
            :input,
            :input,
            :data,
            Float32,
            (configuration.extent,),
            "test.graph.signal.f32/1",
            :column_major,
        ),
    )
end

function AOG.prepare_graph_node(
    ::Type{GraphTestSinkDeclaration},
    configuration::GraphTestSinkConfiguration,
    props,
    inputs::NamedTuple,
    ::NamedTuple{()},
    parameters::NamedTuple,
    target,
)
    return GraphTestSinkOwner(0.0f0, inputs.input)
end

function AOG.step_graph_node!(owner::GraphTestSinkOwner)
    total = 0.0f0
    @inbounds for value in owner.input
        total += value
    end
    owner.sum = total
    return nothing
end

function AOG.reset_graph_node!(owner::GraphTestSinkOwner)
    owner.sum = 0.0f0
    return nothing
end

struct GraphTestFailureDeclaration end

struct GraphTestFailureOwner{I,O}
    input::I
    output::O
end

function AOG.graph_node_ports(
    ::Type{GraphTestFailureDeclaration},
    extent::Int,
)
    format(name, direction) = AOG.graph_port_contract(
        name,
        direction,
        :data,
        Float32,
        (extent,),
        "test.graph.signal.f32/1",
        :column_major,
    )
    return (format(:input, :input), format(:output, :output))
end

function AOG.prepare_graph_node(
    ::Type{GraphTestFailureDeclaration},
    extent::Int,
    props,
    inputs::NamedTuple,
    outputs::NamedTuple,
    parameters::NamedTuple,
    target,
)
    return GraphTestFailureOwner(inputs.input, outputs.output)
end

function AOG.step_graph_node!(owner::GraphTestFailureOwner)
    any(<(0.0f0), owner.input) && throw(DomainError(
        owner.input,
        "negative graph fixture input",
    ))
    copyto!(owner.output, owner.input)
    return nothing
end

AOG.reset_graph_node!(::GraphTestFailureOwner) = nothing

@testset "algorithm graph dependency-derived execution order" begin
    schema = "test.graph.signal.f32/1"
    source = algorithm_node(:source, GraphTestSourceDeclaration,
        GraphTestSourceConfiguration(2, 2.0f0))
    add = algorithm_node(:add, GraphTestAddDeclaration,
        GraphTestAddConfiguration(2, schema, 1.0f0))
    response = second_order_response_node(:response; natural_frequency_hz=40,
        damping_ratio=1, sample_period_s=0.001, count=2, schema)
    sink = algorithm_node(:sink, GraphTestSinkDeclaration, GraphTestSinkConfiguration(2))
    other = algorithm_node(:other, GraphTestSourceDeclaration,
        GraphTestSourceConfiguration(2, 7.0f0))
    isolated = algorithm_node(:isolated, GraphTestSourceDeclaration,
        GraphTestSourceConfiguration(2, 9.0f0))
    definition = algorithm_graph((sink, response, other, add, source, isolated);
        outputs=(graph_output(:position, :response => :output),
            graph_output(:isolated, :isolated => :output)),
        links=(link(:response => :output, :sink => :input),
            link(:add => :output, :response => :input),
            link(:source => :output, :add => :input)))
    graph = prepare_algorithm_graph(definition)
    @test AOG._node_names(definition.nodes) == (:sink, :response, :other, :add, :source, :isolated)
    # Select the earliest declared ready node again after every dependency release.
    @test keys(graph.nodes) == (:other, :source, :add, :response, :sink, :isolated)
    source_owner = prepared_graph_node(graph, :source)
    add_owner = prepared_graph_node(graph, :add)
    response_owner = prepared_graph_node(graph, :response)
    sink_owner = prepared_graph_node(graph, :sink)
    @test add_owner.input === source_owner.output
    @test response_owner.prepared.input === add_owner.output
    @test sink_owner.input === response_owner.prepared.output
    position = graph_output(graph, :position)
    @test response_owner.prepared.output === position
    step_graph!(graph)
    omega_dt = 2π * 40 * 0.001
    expected = 3 * (1 - exp(-omega_dt) * (1 + omega_dt))
    @test position ≈ fill(expected, 2) rtol=5e-6
    @test sink_owner.sum ≈ 2expected rtol=5e-6
    @test graph_output(graph, :isolated) == fill(9.0f0, 2)
    first_position = copy(position)
    @test warmed_graph_step_allocation_bytes(graph) == 0
    @test @inferred(step_graph!(graph)) === graph
    reset_graph!(graph)
    @test prepared_graph_node(graph, :response) === response_owner
    @test graph_output(graph, :position) === position
    step_graph!(graph)
    @test position == first_position

    disconnected = prepare_algorithm_graph(algorithm_graph((other, source, isolated)))
    @test keys(disconnected.nodes) == (:other, :source, :isolated)

    # The explicit grouped policy retains its declaration-order admission rule.
    @test_throws AlgorithmGraphError prepare_algorithm_graph(definition;
        execution=GroupedStreamGraphExecution((:other, :source), (:add,),
            (:response,), (:sink, :isolated)))
    @test_throws AlgorithmGraphError prepare_algorithm_graph(definition;
        execution=GroupedStreamGraphExecution((:sink,), (:response, :other),
            (:add,), (:source, :isolated)))
end

@testset "algorithm graph dependency fan-in and repeated node pairs" begin
    schema = "test.graph.signal.f32/1"
    source = algorithm_node(:source, GraphTestSourceDeclaration,
        GraphTestSourceConfiguration(2, 2.0f0))
    left = algorithm_node(:left, GraphTestAddDeclaration,
        GraphTestAddConfiguration(2, schema, 1.0f0))
    right = algorithm_node(:right, GraphTestAddDeclaration,
        GraphTestAddConfiguration(2, schema, 2.0f0))
    sum_node = algorithm_node(:sum, GraphTestSumDeclaration,
        GraphTestSumConfiguration(2, schema))
    diamond = prepare_algorithm_graph(algorithm_graph((sum_node, right, left, source);
        outputs=(graph_output(:result, :sum => :output),),
        links=(link(:source => :output, :left => :input),
            link(:source => :output, :right => :input),
            link(:left => :output, :sum => :left),
            link(:right => :output, :sum => :right))))
    @test keys(diamond.nodes) == (:source, :right, :left, :sum)
    sum_owner = prepared_graph_node(diamond, :sum)
    @test sum_owner.left === prepared_graph_node(diamond, :left).output
    @test sum_owner.right === prepared_graph_node(diamond, :right).output
    step_graph!(diamond)
    @test graph_output(diamond, :result) == fill(7.0f0, 2)

    shared = prepare_algorithm_graph(algorithm_graph((sum_node, source);
        outputs=(graph_output(:result, :sum => :output),),
        links=(link(:source => :output, :sum => :left),
            link(:source => :output, :sum => :right))))
    @test keys(shared.nodes) == (:source, :sum)
    shared_sum = prepared_graph_node(shared, :sum)
    @test shared_sum.left === shared_sum.right === prepared_graph_node(shared, :source).output
    step_graph!(shared)
    @test graph_output(shared, :result) == fill(4.0f0, 2)
end

@testset "algorithm graph direct cycles and delayed feedback" begin
    configuration = GraphTestAddConfiguration(2, "test.graph.signal.f32/1", 1.0f0)
    first = algorithm_node(:first, GraphTestAddDeclaration, configuration)
    second = algorithm_node(:second, GraphTestAddDeclaration, configuration)
    for definition in (
        algorithm_graph((first,); links=(link(:first => :output, :first => :input),)),
        algorithm_graph((first, second); links=(
            link(:first => :output, :second => :input),
            link(:second => :output, :first => :input))),
    )
        error = try
            prepare_algorithm_graph(definition)
            nothing
        catch exception
            exception
        end
        @test error isa AlgorithmGraphError
        @test occursin("cycle", sprint(showerror, error))
        @test occursin("first", sprint(showerror, error))
        @test occursin("delayed_link", sprint(showerror, error))
    end

    graph = prepare_algorithm_graph(algorithm_graph((second, first);
        outputs=(graph_output(:result, :second => :output),),
        links=(link(:first => :output, :second => :input),),
        delayed_links=(delayed_link(:second => :output, :first => :input,
            zeros(Float32, 2)),)))
    @test keys(graph.nodes) == (:first, :second)
    initial = copy(graph.state.delayed_values[1])
    ticket = step_graph_async!(graph)
    @test graph_step_sequence(graph) == UInt64(0)
    @test graph_step_pending(graph)
    wait_graph_step!(ticket)
    @test graph_output(graph, :result) == Float32[2, 2]
    @test graph_step_sequence(graph) == UInt64(1)
    step_graph!(graph)
    @test graph_output(graph, :result) == Float32[4, 4]
    reset_graph!(graph)
    @test graph.state.delayed_values[1] == initial
    step_graph!(graph)
    @test graph_output(graph, :result) == Float32[2, 2]

    input = Float32[-3, -3]
    failing = prepare_algorithm_graph(algorithm_graph((
        algorithm_node(:checked, GraphTestFailureDeclaration, 2),
        algorithm_node(:feedback, GraphTestAddDeclaration, configuration),
        algorithm_node(:add, GraphTestAddDeclaration, configuration),
        algorithm_node(:source, GraphTestExactBindingDeclaration, 2));
        inputs=(graph_input(:input, :source => :input, input),),
        links=(link(:source => :output, :add => :input),
            link(:add => :output, :checked => :input)),
        delayed_links=(delayed_link(:checked => :output, :feedback => :input,
            fill(5.0f0, 2)),)))
    @test keys(failing.nodes) == (:feedback, :source, :add, :checked)
    @test_throws DomainError step_graph!(failing)
    @test graph_failed(failing)
    @test graph_step_sequence(failing) == UInt64(0)
    @test failing.state.delayed_values[1] == fill(5.0f0, 2)
    @test_throws AlgorithmGraphError step_graph!(failing)
    fill!(input, 3.0f0)
    reset_graph!(failing)
    step_graph!(failing)
    @test !graph_failed(failing)
    @test graph_step_sequence(failing) == UInt64(1)
    @test failing.state.delayed_values[1] == fill(4.0f0, 2)
    @test prepared_graph_node(failing, :feedback).output == fill(6.0f0, 2)
end

@testset "graph files preserve identities with dependency ordering" begin
    document = Dict("schema_version" => 1, "name" => "scrambled",
        "nodes" => [Dict("name" => "sink", "type" => "test_sink", "config" => Dict()),
            Dict("name" => "response", "type" => "second_order_response_f32",
                "config" => Dict("natural_frequency_hz" => 40, "damping_ratio" => 1,
                    "sample_period_s" => 0.001, "count" => 2,
                    "schema" => "test.graph.signal.f32/1")),
            Dict("name" => "source", "type" => "test_source", "config" => Dict())],
        "links" => [Dict("source" => "source.output", "destination" => "response.input"),
            Dict("source" => "response.output", "destination" => "sink.input")],
        "outputs" => [Dict("name" => "position", "source" => "response.output")])
    node_types = merge(builtin_graph_node_types(), (
        test_source=(name, config, props) -> algorithm_node(name, GraphTestSourceDeclaration,
            GraphTestSourceConfiguration(2, 2.0f0)),
        test_sink=(name, config, props) -> algorithm_node(name, GraphTestSinkDeclaration,
            GraphTestSinkConfiguration(2)),))
    definition = load_algorithm_graph(document; node_types)
    graph = prepare_algorithm_graph(definition)
    @test AOG._node_names(definition.nodes) == (:sink, :response, :source)
    @test keys(graph.nodes) == (:source, :response, :sink)
    step_graph!(graph)
    @test prepared_graph_node(graph, :sink).sum ≈ sum(graph_output(graph, :position))
    @test all(>(0), graph_output(graph, :position))
end

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
    target = AdaptiveOpticsSim.Backends.HostComputeDevice()
    graph = prepare_algorithm_graph(definition; target=target)

    @test AdaptiveOpticsSim.Backends.compute_device(graph) == target
    @test AdaptiveOpticsSim.Backends.compute_device(
        graph_output(graph, Val(:command)),
    ) ==
        target
    @test graph_input(graph, Val(:residual)) === input
    @test graph_input(graph, :residual) === input
    @test prepared_graph_node(graph, Val(:gain)) isa GraphTestGainOwner
    @test iszero(graph_step_sequence(graph))
    @test all(iszero, graph_output(graph, Val(:command)))
    @test step_graph!(graph) === graph
    @test graph_output(graph, Val(:command)) == Float32[2, 6, 12]
    @test graph_output(graph, :command) === graph_output(graph, Val(:command))
    @test prepared_graph_node(graph, Val(:sink)).sum == 20.0f0
    @test graph_step_sequence(graph) == UInt64(1)
    @test !graph_failed(graph)

    @test warmed_graph_step_allocation_bytes(graph) == 0
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

@testset "portable algorithm graph asynchronous completion ownership" begin
    input = Float32[1, 2, 3]
    graph = prepare_algorithm_graph(algorithm_graph(
        (algorithm_node(
            :gain,
            GraphTestGainDeclaration,
            GraphTestGainConfiguration(3),
        ),);
        inputs=(graph_input(:input, :gain => :input, input),),
        outputs=(graph_output(:output, :gain => :output),),
        parameters=(sparse_parameter(:gain => :gain, Float32[2, 3, 4]),),
    ))

    ticket = @inferred step_graph_async!(graph)
    @test ticket isa GraphStepTicket{typeof(graph)}
    @test graph_step_pending(graph)
    @test iszero(graph_step_sequence(graph))
    @test_throws AlgorithmGraphError step_graph_async!(graph)
    @test_throws AlgorithmGraphError reset_graph!(graph)
    @test_throws AlgorithmGraphError wait_graph_step!(
        GraphStepTicket(graph, UInt64(2)),
    )

    @test @inferred(wait_graph_step!(ticket)) === graph
    @test !graph_step_pending(graph)
    @test graph_step_sequence(graph) == UInt64(1)
    @test graph_output(graph, Val(:output)) == Float32[2, 6, 12]
    @test_throws AlgorithmGraphError wait_graph_step!(ticket)

    @test warmed_async_graph_step_allocation_bytes(graph) == 0
    @test !graph_step_pending(graph)
end

@testset "bounded grouped stream graph execution" begin
    schema = "test.graph.signal.f32/1"
    source_configuration_a = GraphTestSourceConfiguration(3, 2.0f0)
    source_configuration_b = GraphTestSourceConfiguration(3, 5.0f0)
    add_configuration_a = GraphTestAddConfiguration(3, schema, 1.0f0)
    add_configuration_b = GraphTestAddConfiguration(3, schema, 3.0f0)
    sum_configuration = GraphTestSumConfiguration(3, schema)
    definition = algorithm_graph(
        (
            algorithm_node(
                :source_a,
                GraphTestSourceDeclaration,
                source_configuration_a,
            ),
            algorithm_node(
                :source_b,
                GraphTestSourceDeclaration,
                source_configuration_b,
            ),
            algorithm_node(
                :add_a,
                GraphTestAddDeclaration,
                add_configuration_a,
            ),
            algorithm_node(
                :add_b,
                GraphTestAddDeclaration,
                add_configuration_b,
            ),
            algorithm_node(
                :sum,
                GraphTestSumDeclaration,
                sum_configuration,
            ),
        );
        name=:grouped_stream_diamond,
        outputs=(graph_output(:result, :sum => :output),),
        links=(
            link(:source_a => :output, :add_a => :input),
            link(:source_b => :output, :add_b => :input),
            link(:add_a => :output, :sum => :left),
            link(:add_b => :output, :sum => :right),
        ),
    )
    execution = GroupedStreamGraphExecution(
        (:source_a, :source_b),
        (:add_a, :add_b),
        (:sum,),
    )
    graph = prepare_algorithm_graph(definition; execution)

    @test graph_execution_policy(graph) === execution
    @test captured_graph_node_count(graph) == 0
    @test length(graph.execution.contexts) == 2
    @test graph.execution.contexts isa FixedSizeVector
    @test graph.execution.events isa FixedSizeVector
    @test @inferred(step_graph!(graph)) === graph
    @test graph_output(graph, Val(:result)) == fill(11.0f0, 3)
    @test graph_step_sequence(graph) == UInt64(1)
    @test warmed_async_graph_step_allocation_bytes(graph) == 0

    @test reset_graph!(graph) === graph
    @test iszero(graph_step_sequence(graph))
    @test step_graph!(graph) === graph
    @test graph_output(graph, Val(:result)) == fill(11.0f0, 3)

    @test_throws AlgorithmGraphError GroupedStreamGraphExecution()
    @test_throws AlgorithmGraphError GroupedStreamGraphExecution(())
    @test_throws AlgorithmGraphError GroupedStreamGraphExecution(
        (:source_a, :source_a),
    )
    @test_throws AlgorithmGraphError prepare_algorithm_graph(
        definition;
        execution=GroupedStreamGraphExecution(
            (:source_b, :source_a),
            (:add_a, :add_b),
            (:sum,),
        ),
    )
    @test_throws AlgorithmGraphError prepare_algorithm_graph(
        definition;
        execution=GroupedStreamGraphExecution(
            (:source_a, :source_b, :add_a, :add_b),
            (:sum,),
        ),
    )
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
    @test prepared_graph_node(sink, Val(:sink)).sum == 7.0f0
end

@testset "portable algorithm graph binds exact node execution owners" begin
    input = Float32[3, 5]
    graph = prepare_algorithm_graph(algorithm_graph(
        (algorithm_node(:bound, GraphTestExactBindingDeclaration, 2),);
        inputs=(graph_input(:sample, :bound => :input, input),),
        outputs=(graph_output(:copy, :bound => :output),),
    ))

    owner = prepared_graph_node(graph, Val(:bound))
    @test owner isa GraphTestExactBindingOwner
    @test owner.input === input
    @test owner.output === graph_output(graph, Val(:copy))
    @test step_graph!(graph) === graph
    @test graph_output(graph, Val(:copy)) == input
    @test @allocated(step_graph!(graph)) == 0
    @test @inferred(step_graph!(graph)) === graph
end

@testset "graph HIL lockstep boundary" begin
    initial_command = Float32[1, 2]
    graph = prepare_algorithm_graph(algorithm_graph(
        (algorithm_node(:bound, GraphTestExactBindingDeclaration, 2),);
        inputs=(graph_input(:pdm_command, :bound => :input, initial_command),),
        outputs=(graph_output(:shwfs_frame, :bound => :output),),
    ))
    command_buffer = zeros(Float32, 2)
    frame_buffer = zeros(Float32, 2)
    boundary = prepare_graph_hil_boundary(
        graph;
        command_input=:pdm_command,
        frame_output=:shwfs_frame,
        command_buffer,
        frame_buffer,
    )

    @test boundary isa PreparedGraphHILBoundary
    @test hil_command_buffer(boundary) === command_buffer
    @test hil_frame_buffer(boundary) === frame_buffer
    @test command_buffer == initial_command
    @test all(iszero, frame_buffer)
    @test @inferred(hil_boundary_status(boundary)) == (
        frame_sequence=UInt64(0),
        active_command_sequence=UInt64(0),
        command_response_required=false,
        failed=false,
    )
    @test_throws AlgorithmGraphError adopt_hil_command!(boundary, UInt64(1))

    @test @inferred(step_hil_frame!(boundary)) == UInt64(1)
    @test frame_buffer == initial_command
    @test hil_boundary_status(boundary).command_response_required
    @test_throws AlgorithmGraphError step_hil_frame!(boundary)
    @test_throws AlgorithmGraphError adopt_hil_command!(boundary, 1)

    command_buffer .= Float32[3, 4]
    @test_throws AlgorithmGraphError adopt_hil_command!(boundary, UInt64(2))
    @test graph_input(graph, Val(:pdm_command)) == initial_command
    command_buffer[1] = NaN32
    @test_throws AlgorithmGraphError adopt_hil_command!(boundary, UInt64(1))
    @test graph_input(graph, Val(:pdm_command)) == initial_command
    @test !hil_boundary_status(boundary).failed

    command_buffer .= Float32[3, 4]
    @test @inferred(adopt_hil_command!(boundary, UInt64(1))) == UInt64(1)
    @test graph_input(graph, Val(:pdm_command)) == Float32[3, 4]
    @test step_hil_frame!(boundary) == UInt64(2)
    @test frame_buffer == Float32[3, 4]

    command_buffer .= Float32[5, 6]
    @test @allocated(adopt_hil_command!(boundary, UInt64(2))) == 0
    @test @allocated(step_hil_frame!(boundary)) == 0
    @test frame_buffer == Float32[5, 6]
    command_buffer .= Float32[7, 8]
    @test @allocated(adopt_hil_command!(boundary, UInt64(3))) == 0
    @test hil_boundary_status(boundary) == (
        frame_sequence=UInt64(3),
        active_command_sequence=UInt64(3),
        command_response_required=false,
        failed=false,
    )

    @test reset_hil_boundary!(boundary) === boundary
    @test graph_input(graph, Val(:pdm_command)) == initial_command
    @test command_buffer == initial_command
    @test all(iszero, frame_buffer)
    @test hil_boundary_status(boundary).frame_sequence == UInt64(0)

    timed_graph = prepare_algorithm_graph(algorithm_graph(
        (algorithm_node(:bound, GraphTestExactBindingDeclaration, 2),);
        inputs=(graph_input(
            :pdm_command,
            :bound => :input,
            Float32[9, 10],
        ),),
        outputs=(graph_output(:shwfs_frame, :bound => :output),),
    ))
    timed_boundary = prepare_graph_hil_boundary(
        timed_graph;
        command_input=:pdm_command,
        frame_output=:shwfs_frame,
    )
    driver = FixedStepModelTimeDriver(AlgorithmGraphs.PeriodicSchedule(ModelDuration(50)))
    @test step_hil_frame_at!(timed_boundary, driver) == (
        sequence=UInt64(1),
        timestamp=ModelTimestamp(0),
    )
    hil_command_buffer(timed_boundary) .= Float32[11, 12]
    adopt_hil_command!(timed_boundary, UInt64(1))
    @test @allocated(step_hil_frame_at!(timed_boundary, driver)) == 0
    @test hil_frame_buffer(timed_boundary) == Float32[11, 12]
    @test reset_hil_boundary!(timed_boundary, driver) === timed_boundary
    @test model_time_sequence(driver) == UInt64(0)

    unbound_graph = prepare_algorithm_graph(algorithm_graph(
        (algorithm_node(:bound, GraphTestExactBindingDeclaration, 2),);
        inputs=(graph_input(
            :pdm_command,
            :bound => :input,
            Float32[1, 2],
        ),),
        outputs=(graph_output(:shwfs_frame, :bound => :output),),
    ))
    @test_throws AlgorithmGraphError prepare_graph_hil_boundary(
        unbound_graph;
        command_input=:pdm_command,
        frame_output=:shwfs_frame,
        command_buffer=graph_input(unbound_graph, Val(:pdm_command)),
    )
    @test_throws AlgorithmGraphError prepare_graph_hil_boundary(
        unbound_graph;
        command_input=:pdm_command,
        frame_output=:shwfs_frame,
        frame_buffer=zeros(Float32, 3),
    )

    nonfinite_graph = prepare_algorithm_graph(algorithm_graph(
        (algorithm_node(:bound, GraphTestExactBindingDeclaration, 2),);
        inputs=(graph_input(
            :pdm_command,
            :bound => :input,
            Float32[NaN, 0],
        ),),
        outputs=(graph_output(:shwfs_frame, :bound => :output),),
    ))
    @test_throws AlgorithmGraphError prepare_graph_hil_boundary(
        nonfinite_graph;
        command_input=:pdm_command,
        frame_output=:shwfs_frame,
    )

    failing_graph = prepare_algorithm_graph(algorithm_graph(
        (algorithm_node(:checked, GraphTestFailureDeclaration, 2),);
        inputs=(graph_input(
            :pdm_command,
            :checked => :input,
            Float32[-1, 1],
        ),),
        outputs=(graph_output(:shwfs_frame, :checked => :output),),
    ))
    failing_boundary = prepare_graph_hil_boundary(
        failing_graph;
        command_input=:pdm_command,
        frame_output=:shwfs_frame,
    )
    @test_throws DomainError step_hil_frame!(failing_boundary)
    @test hil_boundary_status(failing_boundary).failed
    @test_throws AlgorithmGraphError step_hil_frame!(failing_boundary)
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

@testset "native modal OPD graph node" begin
    coefficients = Float32[0.1, 0.2]
    basis = zeros(Float32, 2, 2, 2)
    fill!(@view(basis[:, :, 1]), 1.0f0)
    fill!(@view(basis[:, :, 2]), 2.0f0)
    pupil_support = Bool[true false; true true]
    coefficient_schema = "test.graph.modal-coefficients.f32/1"

    definition = algorithm_graph(
        (
            modal_opd_expansion_node(
                :modal_opd;
                pupil_rows=2,
                pupil_columns=2,
                mode_count=2,
                coefficients_schema=coefficient_schema,
                opd_schema="test.graph.opd.f32/1",
                basis_schema="test.graph.modal-basis.f32/1",
                pupil_support_schema="test.graph.pupil-support.bool/1",
            ),
        );
        name=:native_modal_opd,
        inputs=(graph_input(
            :coefficients,
            :modal_opd => :coefficients,
            coefficients,
        ),),
        outputs=(graph_output(:opd, :modal_opd => :opd),),
        parameters=(
            sparse_parameter(:modal_opd => :basis, basis),
            sparse_parameter(:modal_opd => :pupil_support, pupil_support),
        ),
    )
    graph = prepare_algorithm_graph(definition)

    @test graph_name(definition) === :native_modal_opd
    @test graph_name(graph) === :native_modal_opd
    @test step_graph!(graph) === graph
    @test graph_output(graph, Val(:opd)) ≈ Float32[0.5 0.0; 0.5 0.5]
    coefficients .= Float32[0.3, -0.1]
    step_graph!(graph)
    @test graph_output(graph, Val(:opd)) ≈ Float32[0.1 0.0; 0.1 0.1]
    @test @allocated(step_graph!(graph)) == 0
    @test @inferred(step_graph!(graph)) === graph
end

@testset "measured deformable-mirror surface graph node" begin
    pdm_command = Float32[2, -3]
    actuator_coordinates = Float32[-0.5 0.5; 0 0]
    influence_functions = Float32[
        1 0
        0 1
        0 1
        1 0
    ]
    node = deformable_mirror_surface_node(
        :dm;
        resolution=2,
        telescope_diameter_m=1.22,
        actuator_count=2,
        actuator_model=ClippedActuators(-1.0f0, 1.0f0),
        pdm_command_schema="test.graph.pdm-command.f32/1",
        surface_opd_schema="test.graph.dm-surface-opd.f32/1",
        actuator_coordinates_schema=
            "test.graph.dm-actuator-coordinates.f32/1",
        influence_functions_schema=
            "test.graph.dm-influence-functions.f32/1",
    )
    definition = algorithm_graph(
        (node,);
        name=:measured_deformable_mirror,
        inputs=(
            graph_input(:pdm_command, :dm => :pdm_command, pdm_command),
        ),
        outputs=(
            graph_output(:surface_opd, :dm => :surface_opd),
        ),
        parameters=(
            sparse_parameter(
                :dm => :actuator_coordinates,
                actuator_coordinates,
            ),
            sparse_parameter(
                :dm => :influence_functions,
                influence_functions,
            ),
        ),
    )
    graph = prepare_algorithm_graph(definition)
    owner = prepared_graph_node(graph, Val(:dm))
    output = graph_output(graph, Val(:surface_opd))

    @test owner.pdm_command === pdm_command
    @test owner.output === output
    @test influence_model(owner.deformable_mirror).modes !==
        influence_functions
    @test influence_model(owner.deformable_mirror).modes ==
        influence_functions
    @test @inferred(step_graph!(graph)) === graph
    @test output == Float32[1 -1; -1 1]
    pdm_command .= Float32[0.25, 0.5]
    step_graph!(graph)
    @test output == Float32[0.25 0.5; 0.5 0.25]
    @test warmed_graph_step_allocation_bytes(graph) == 0

    reset_graph!(graph)
    @test all(iszero, output)
    step_graph!(graph)
    @test output == Float32[0.25 0.5; 0.5 0.25]

    @test_throws AlgorithmGraphError deformable_mirror_surface_node(
        :invalid_dm;
        resolution=2,
        telescope_diameter_m=1.22,
        actuator_count=0,
        pdm_command_schema="test.graph.pdm-command.f32/1",
        surface_opd_schema="test.graph.dm-surface-opd.f32/1",
        actuator_coordinates_schema=
            "test.graph.dm-actuator-coordinates.f32/1",
        influence_functions_schema=
            "test.graph.dm-influence-functions.f32/1",
    )

    nonfinite_influence_functions = copy(influence_functions)
    nonfinite_influence_functions[1] = NaN32
    invalid_definition = algorithm_graph(
        (node,);
        inputs=(
            graph_input(:pdm_command, :dm => :pdm_command, pdm_command),
        ),
        parameters=(
            sparse_parameter(
                :dm => :actuator_coordinates,
                actuator_coordinates,
            ),
            sparse_parameter(
                :dm => :influence_functions,
                nonfinite_influence_functions,
            ),
        ),
    )
    @test_throws AlgorithmGraphError prepare_algorithm_graph(
        invalid_definition,
    )
end

@testset "analytic Gaussian deformable-mirror surface graph node" begin
    pdm_command = Float32[2.0f-8, -1.0f-8]
    actuator_coordinates = Float32[-0.25 0.25; -0.25 0.25]
    node = gaussian_deformable_mirror_surface_node(
        :dm;
        resolution=4,
        telescope_diameter_m=1.22,
        actuator_count=2,
        influence_width=0.2,
        pdm_command_schema="test.graph.pdm-actuator-opd-m.f32/1",
        surface_opd_schema="test.graph.dm-surface-opd-m.f32/1",
        actuator_coordinates_schema=
            "test.graph.dm-normalized-pupil-coordinates.f32/1",
    )
    definition = algorithm_graph(
        (node,);
        name=:analytic_gaussian_deformable_mirror,
        inputs=(
            graph_input(:pdm_command, :dm => :pdm_command, pdm_command),
        ),
        outputs=(
            graph_output(:surface_opd, :dm => :surface_opd),
        ),
        parameters=(
            sparse_parameter(
                :dm => :actuator_coordinates,
                actuator_coordinates,
            ),
        ),
    )
    graph = prepare_algorithm_graph(definition)
    owner = prepared_graph_node(graph, Val(:dm))
    output = graph_output(graph, Val(:surface_opd))

    @test influence_model(owner.deformable_mirror) ==
        GaussianInfluenceWidth(0.2f0)
    @test @inferred(step_graph!(graph)) === graph
    diagonal_response = exp(-0.5f0 / (2 * 0.2f0^2))
    cross_response = exp(-0.25f0 / (2 * 0.2f0^2))
    @test output[2, 2] ≈ 2.0f-8 - 1.0f-8 * diagonal_response
    @test output[3, 3] ≈ 2.0f-8 * diagonal_response - 1.0f-8
    @test output[2, 3] ≈ 1.0f-8 * cross_response
    @test iszero(output[1, 1])
    @test warmed_graph_step_allocation_bytes(graph) == 0

    reset_graph!(graph)
    @test all(iszero, output)
    step_graph!(graph)
    @test output[2, 2] ≈ 2.0f-8 - 1.0f-8 * diagonal_response

    @test_throws AlgorithmGraphError gaussian_deformable_mirror_surface_node(
        :invalid_dm;
        resolution=4,
        telescope_diameter_m=1.22,
        actuator_count=2,
        influence_width=0.0,
        pdm_command_schema="test.graph.pdm-actuator-opd-m.f32/1",
        surface_opd_schema="test.graph.dm-surface-opd-m.f32/1",
        actuator_coordinates_schema=
            "test.graph.dm-normalized-pupil-coordinates.f32/1",
    )
end

@testset "separable grid Gaussian deformable-mirror graph node" begin
    actuator_axis_count = 3
    actuator_pitch = 0.4f0
    influence_width = 0.2f0
    actuator_grid_indices = Int32[1, 3, 5, 7, 9]
    pdm_command = Float32[2, -1, 0.5, 1.5, -0.25] .* 1.0f-8
    grid_coordinates = Float32[-0.4, 0, 0.4]
    actuator_coordinates = Matrix{Float32}(
        undef,
        2,
        length(actuator_grid_indices),
    )
    for command_index in eachindex(actuator_grid_indices)
        grid_index = CartesianIndices((3, 3))[
            actuator_grid_indices[command_index]
        ]
        actuator_coordinates[1, command_index] =
            grid_coordinates[grid_index[1]]
        actuator_coordinates[2, command_index] =
            grid_coordinates[grid_index[2]]
    end
    general_node = gaussian_deformable_mirror_surface_node(
        :general_dm;
        resolution=16,
        telescope_diameter_m=1.22,
        actuator_count=length(pdm_command),
        influence_width,
        pdm_command_schema="test.graph.pdm-actuator-opd-m.f32/1",
        surface_opd_schema="test.graph.dm-surface-opd-m.f32/1",
        actuator_coordinates_schema=
            "test.graph.dm-normalized-pupil-coordinates.f32/1",
    )
    grid_node = grid_gaussian_deformable_mirror_surface_node(
        :grid_dm;
        resolution=16,
        telescope_diameter_m=1.22,
        actuator_count=length(pdm_command),
        actuator_axis_count,
        actuator_pitch,
        influence_width,
        pdm_command_schema="test.graph.pdm-actuator-opd-m.f32/1",
        surface_opd_schema="test.graph.dm-surface-opd-m.f32/1",
        actuator_grid_indices_schema=
            "test.graph.dm-actuator-grid-indices.i32/1",
    )
    definition = algorithm_graph(
        (general_node, grid_node);
        name=:separable_grid_gaussian_deformable_mirror,
        inputs=(
            graph_input(
                :general_command,
                :general_dm => :pdm_command,
                copy(pdm_command),
            ),
            graph_input(
                :grid_command,
                :grid_dm => :pdm_command,
                copy(pdm_command),
            ),
        ),
        outputs=(
            graph_output(:general_surface, :general_dm => :surface_opd),
            graph_output(:grid_surface, :grid_dm => :surface_opd),
        ),
        parameters=(
            sparse_parameter(
                :general_dm => :actuator_coordinates,
                actuator_coordinates,
            ),
            sparse_parameter(
                :grid_dm => :actuator_grid_indices,
                actuator_grid_indices,
            ),
        ),
    )
    graph = prepare_algorithm_graph(definition)
    @test graph_execution_policy(graph) isa StreamGraphExecution
    @test captured_graph_node_count(graph) == 0
    step_graph!(graph)
    general_surface = graph_output(graph, Val(:general_surface))
    grid_surface = graph_output(graph, Val(:grid_surface))
    owner = prepared_graph_node(graph, Val(:grid_dm))
    @test owner.deformable_mirror.state.separable_x !== nothing
    @test grid_surface ≈ general_surface rtol = 2.0f-5 atol = 1.0f-12
    @test warmed_graph_step_allocation_bytes(graph) == 0

    reset_graph!(graph)
    @test all(iszero, grid_surface)

    @test_throws AlgorithmGraphError prepare_algorithm_graph(
        definition;
        execution=CapturedGraphExecution(),
    )

    duplicate_mapping = algorithm_graph(
        (grid_node,);
        inputs=(graph_input(
            :grid_command,
            :grid_dm => :pdm_command,
            copy(pdm_command),
        ),),
        parameters=(sparse_parameter(
            :grid_dm => :actuator_grid_indices,
            Int32[1, 1, 5, 7, 9],
        ),),
    )
    @test_throws AlgorithmGraphError prepare_algorithm_graph(
        duplicate_mapping,
    )
    @test_throws AlgorithmGraphError grid_gaussian_deformable_mirror_surface_node(
        :invalid_grid_dm;
        resolution=16,
        telescope_diameter_m=1.22,
        actuator_count=5,
        actuator_axis_count=1,
        actuator_pitch,
        influence_width,
        pdm_command_schema="test.graph.pdm-actuator-opd-m.f32/1",
        surface_opd_schema="test.graph.dm-surface-opd-m.f32/1",
        actuator_grid_indices_schema=
            "test.graph.dm-actuator-grid-indices.i32/1",
    )
end

@testset "finite multilayer-atmosphere OPD graph node" begin
    node = multilayer_atmosphere_opd_node(
        :atmosphere;
        resolution=16,
        telescope_diameter_m=1.22,
        r0=0.15,
        reference_wavelength_m=500e-9,
        L0=30.0,
        fractional_cn2=(0.7, 0.3),
        wind_speed=(5.0, 9.0),
        wind_direction_deg=(0.0, 90.0),
        altitude=(0.0, 5000.0),
        layer_ids=(:ground, :high),
        atmosphere_step=0.001,
        rng_seed=17,
        atmosphere_opd_schema="test.graph.atmosphere-opd-m.f32/1",
    )
    definition = algorithm_graph(
        (node,);
        name=:multilayer_atmosphere_opd,
        outputs=(
            graph_output(:atmosphere_opd, :atmosphere => :atmosphere_opd),
        ),
    )
    graph = prepare_algorithm_graph(definition)
    output = graph_output(graph, Val(:atmosphere_opd))

    @test @inferred(step_graph!(graph)) === graph
    @test size(output) == (16, 16)
    @test all(isfinite, output)
    @test !all(iszero, output)
    first_frame = copy(output)

    @test warmed_graph_step_allocation_bytes(graph) == 0
    @test output != first_frame

    reset_graph!(graph)
    @test all(iszero, output)
    step_graph!(graph)
    @test output == first_frame

    @test_throws AlgorithmGraphError multilayer_atmosphere_opd_node(
        :invalid_atmosphere;
        resolution=16,
        telescope_diameter_m=1.22,
        r0=0.15,
        reference_wavelength_m=500e-9,
        fractional_cn2=(1.0,),
        wind_speed=(5.0,),
        wind_direction_deg=(0.0,),
        altitude=(0.0,),
        layer_ids=(:ground,),
        atmosphere_step=0.0,
        rng_seed=17,
        atmosphere_opd_schema="test.graph.atmosphere-opd-m.f32/1",
    )
end

@testset "atmosphere YAML follows the declared source direction" begin
    config = Dict{String,Any}(
        "resolution" => 16, "telescope_diameter_m" => 1.22,
        "r0" => 0.15, "reference_wavelength_m" => 500e-9,
        "fractional_cn2" => [0.5, 0.5], "wind_speed" => [5.0, 9.0],
        "wind_direction_deg" => [0.0, 90.0], "altitude" => [0.0, 5000.0],
        "layer_ids" => ["ground", "high"], "atmosphere_step" => 0.001,
        "rng_seed" => 17, "atmosphere_opd_schema" => "test.atmosphere/1",
    )
    document = Dict("schema_version" => 1, "name" => "directional_atmosphere",
        "nodes" => [Dict("name" => "atmosphere",
            "type" => "multilayer_atmosphere_opd_f32", "config" => config)],
        "outputs" => [Dict("name" => "opd", "source" => "atmosphere.atmosphere_opd")])
    for angle in (0.0, 90.0)
        config["source_separation_arcsec"] = 10.0
        config["source_position_angle_deg"] = angle
        graph = prepare_algorithm_graph(load_algorithm_graph(document))
        owner = prepared_graph_node(graph, Val(:atmosphere))
        source = Source(; separation_arcsec=10, position_angle_deg=angle, T=Float32)
        renderer = prepare_atmosphere_renderer(owner.atmosphere, owner.telescope, source)
        onaxis = prepare_atmosphere_renderer(owner.atmosphere, owner.telescope)
        reference = PupilFunction(owner.telescope)
        unshifted = PupilFunction(owner.telescope)
        @test renderer.shift_x[1] == 0
        @test renderer.shift_y[1] == 0
        displacement = Float32(10 * pi / (180 * 3600) * 5000 / (1.22 / 16))
        @test renderer.shift_x[2] ≈ displacement * cosd(Float32(angle)) atol=1e-6
        @test renderer.shift_y[2] ≈ displacement * sind(Float32(angle)) atol=1e-6
        first_frame = zeros(Float32, 16, 16)
        for frame in 1:2
            step_graph!(graph)
            epoch = current_epoch(owner.atmosphere)
            render_atmosphere!(reference, renderer, owner.atmosphere, epoch)
            render_atmosphere!(unshifted, onaxis, owner.atmosphere, epoch)
            @test graph_output(graph, Val(:opd)) == reference.opd
            @test graph_output(graph, Val(:opd)) != unshifted.opd
            frame == 1 && copyto!(first_frame, reference.opd)
        end
        @test warmed_graph_step_allocation_bytes(graph) == 0
        reset_graph!(graph)
        step_graph!(graph)
        @test graph_output(graph, Val(:opd)) == first_frame
    end
    delete!(config, "source_separation_arcsec")
    delete!(config, "source_position_angle_deg")
    graph = prepare_algorithm_graph(load_algorithm_graph(document))
    owner = prepared_graph_node(graph, Val(:atmosphere))
    onaxis = prepare_atmosphere_renderer(owner.atmosphere, owner.telescope)
    reference = PupilFunction(owner.telescope)
    step_graph!(graph)
    render_atmosphere!(reference, onaxis, owner.atmosphere, current_epoch(owner.atmosphere))
    @test graph_output(graph, Val(:opd)) == reference.opd
    for field in ("source_separation_arcsec", "source_position_angle_deg")
        for invalid in (Inf, NaN, true, "invalid")
            config[field] = invalid
            @test_throws AlgorithmGraphError load_algorithm_graph(document)
        end
        delete!(config, field)
    end
end

@testset "additive pupil-OPD composition graph node" begin
    uncompensated_opd = Float32[1 2; 3 4]
    surface_opd = Float32[-0.25 0.5; 0.75 -1]
    node = pupil_opd_composition_node(
        :composition;
        resolution=2,
        uncompensated_opd_schema="test.graph.uncompensated-opd.f32/1",
        surface_opd_schema="test.graph.surface-opd.f32/1",
        pupil_opd_schema="test.graph.pupil-opd.f32/1",
    )
    definition = algorithm_graph(
        (node,);
        name=:pupil_opd_composition,
        inputs=(
            graph_input(
                :uncompensated_opd,
                :composition => :uncompensated_opd,
                uncompensated_opd,
            ),
            graph_input(
                :surface_opd,
                :composition => :surface_opd,
                surface_opd,
            ),
        ),
        outputs=(
            graph_output(:pupil_opd, :composition => :pupil_opd),
        ),
    )
    graph = prepare_algorithm_graph(definition)
    output = graph_output(graph, Val(:pupil_opd))

    @test @inferred(step_graph!(graph)) === graph
    @test output == Float32[0.75 2.5; 3.75 3]
    @test warmed_graph_step_allocation_bytes(graph) == 0
    reset_graph!(graph)
    @test all(iszero, output)
end

@testset "diffractive Shack-Hartmann photon-rate graph node" begin
    pupil_opd = zeros(Float32, 16, 16)
    node = shack_hartmann_rate_node(
        :shwfs;
        resolution=16,
        telescope_diameter_m=8.0,
        n_lenslets=4,
        n_pix_subap=4,
        pixel_scale_arcsec=0.1,
        source_wavelength_m=0.75e-6,
        source_photon_irradiance_m2_s=1.0,
        opd_schema="test.graph.pupil-opd.f32/1",
        photon_rate_schema="test.graph.shwfs-photon-rate.f32/1",
    )
    definition = algorithm_graph(
        (node,);
        name=:shwfs_photon_rate,
        inputs=(graph_input(:pupil_opd, :shwfs => :opd, pupil_opd),),
        outputs=(
            graph_output(:shwfs_photon_rate, :shwfs => :photon_rate),
        ),
    )
    graph = prepare_algorithm_graph(definition)
    owner = prepared_graph_node(graph, Val(:shwfs))
    output = graph_output(graph, Val(:shwfs_photon_rate))

    @test owner.prepared.input.opd === pupil_opd
    @test owner.prepared.output.values === output
    @test size(output) == (16, 16)
    step_graph!(graph)
    @test all(isfinite, output)
    @test sum(output) > 0
    @test warmed_graph_step_allocation_bytes(graph) == 0
    @test @inferred(step_graph!(graph)) === graph

    @test_throws AlgorithmGraphError shack_hartmann_rate_node(
        :invalid_shwfs;
        resolution=16,
        telescope_diameter_m=8.0,
        n_lenslets=4,
        n_pix_subap=3,
        pixel_scale_arcsec=0.1,
        source_wavelength_m=0.75e-6,
        source_photon_irradiance_m2_s=1.0,
        opd_schema="test.graph.pupil-opd.f32/1",
        photon_rate_schema="test.graph.shwfs-photon-rate.f32/1",
    )
end

@testset "diffractive Pyramid photon-rate graph node" begin
    pupil_opd = zeros(Float32, 8, 8)
    node = pyramid_rate_node(
        :pwfs;
        resolution=8,
        telescope_diameter_m=1.22,
        pupil_samples=4,
        modulation=2,
        modulation_points=4,
        n_pix_separation=2,
        n_pix_edge=1,
        source_wavelength_m=0.55e-6,
        source_photon_irradiance_m2_s=1.0,
        opd_schema="test.graph.pupil-opd.f32/1",
        photon_rate_schema="test.graph.pwfs-photon-rate.f32/1",
    )
    definition = algorithm_graph(
        (node,);
        name=:pwfs_photon_rate,
        inputs=(graph_input(:pupil_opd, :pwfs => :opd, pupil_opd),),
        outputs=(graph_output(:pwfs_photon_rate, :pwfs => :photon_rate),),
    )
    graph = prepare_algorithm_graph(definition)
    owner = prepared_graph_node(graph, Val(:pwfs))
    output = graph_output(graph, Val(:pwfs_photon_rate))

    @test owner.prepared.input.opd === pupil_opd
    @test owner.prepared.output.values === output
    @test size(output) == (12, 12)
    step_graph!(graph)
    @test all(isfinite, output)
    @test sum(output) > 0
    @test warmed_graph_step_allocation_bytes(graph) == 0
    @test @inferred(step_graph!(graph)) === graph

    @test_throws AlgorithmGraphError pyramid_rate_node(
        :invalid_pwfs;
        resolution=8,
        telescope_diameter_m=1.22,
        pupil_samples=4,
        n_pix_edge=1,
        source_wavelength_m=0.55e-6,
        source_photon_irradiance_m2_s=1.0,
        opd_schema="test.graph.pupil-opd.f32/1",
        photon_rate_schema="test.graph.pwfs-photon-rate.f32/1",
    )
end

@testset "fractional Pyramid pupil shifts in graph files" begin
    opd = zeros(Float32, 8, 8)
    options = (resolution=8, telescope_diameter_m=1.22, pupil_samples=4,
        modulation=0, modulation_points=4,
        modulation_propagation_strategy=PyramidPupilTiltStrategy(),
        source_wavelength_m=0.55e-6, source_photon_irradiance_m2_s=1.0,
        opd_schema="test.opd/1", photon_rate_schema="test.rate/1")
    prepare_node(node) = prepare_algorithm_graph(algorithm_graph((node,);
        inputs=(graph_input(:opd, :pwfs => :opd, opd),),
        outputs=(graph_output(:rate, :pwfs => :photon_rate),)))
    shifted = prepare_node(pyramid_rate_node(:pwfs; options...,
        pupil_shift_x_pixels=0.25,
        pupil_shift_y_pixels=(0.1, -0.3, 0.15, -0.2)))
    ideal = prepare_node(pyramid_rate_node(:pwfs; options...))
    step_graph!(shifted)
    step_graph!(ideal)
    @test all(isfinite, graph_output(shifted, Val(:rate)))
    @test !isapprox(graph_output(shifted, Val(:rate)),
        graph_output(ideal, Val(:rate)); rtol=1f-4)
    @test warmed_graph_step_allocation_bytes(shifted) == 0

    text = """
    schema_version: 1
    name: fractional_pupils
    nodes:
    - name: pwfs
      type: pyramid_rate_f32
      config:
        resolution: 8
        telescope_diameter_m: 1.22
        pupil_samples: 4
        modulation: 0.0
        modulation_points: 4
        modulation_propagation_strategy: pupil_tilt
        source_wavelength_m: 5.5e-07
        source_photon_irradiance_m2_s: 1.0
        opd_schema: test.opd/1
        photon_rate_schema: test.rate/1
        pupil_shift_x_pixels: 0.25
        pupil_shift_y_pixels:
        - 0.1
        - -0.3
        - 0.15
        - -0.2
    inputs:
    - name: opd
      destination: pwfs.opd
      binding: opd
    outputs:
    - name: rate
      source: pwfs.photon_rate
    """
    mktemp() do path, io
        write(io, text)
        close(io)
        loaded = prepare_algorithm_graph(load_algorithm_graph(path; bindings=(; opd)))
        step_graph!(loaded)
        @test graph_output(loaded, Val(:rate)) == graph_output(shifted, Val(:rate))
        write(path, replace(text,
            "modulation_propagation_strategy: pupil_tilt" =>
                "modulation_propagation_strategy: shifted_mask"))
        @test_throws AlgorithmGraphError load_algorithm_graph(path; bindings=(; opd))
        for invalid in ("[0.1, 0.2]", "[0.1, 0.2, inf, 0.4]", "true", "1e40")
            write(path, replace(text,
                "pupil_shift_x_pixels: 0.25" => "pupil_shift_x_pixels: $invalid"))
            @test_throws AlgorithmGraphError load_algorithm_graph(path; bindings=(; opd))
        end
    end
end

@testset "single-read CCD detector-acquisition graph node" begin
    photon_rate = fill(1_000.0f0, 4, 4)
    node = ccd_detector_acquisition_node(
        :detector;
        rows=4,
        columns=4,
        pixel_scale_arcsec=0.1,
        wavelength_m=0.75e-6,
        exposure_duration_s=0.01,
        quantum_efficiency=0.5,
        rng_seed=0x1234,
        photon_rate_schema="test.graph.shwfs-photon-rate.f32/1",
        frame_schema="test.graph.shwfs-frame.f32/1",
        photon_noise=true,
        readout_noise=true,
        readout_noise_e=2.7,
    )
    definition = algorithm_graph(
        (node,);
        name=:ccd_detector_acquisition,
        inputs=(
            graph_input(
                :photon_rate,
                :detector => :photon_rate,
                photon_rate,
            ),
        ),
        outputs=(graph_output(:frame, :detector => :frame),),
    )
    graph = prepare_algorithm_graph(definition)
    owner = prepared_graph_node(graph, Val(:detector))
    frame = graph_output(graph, Val(:frame))

    @test owner.prepared.input.values === photon_rate
    @test owner.output === frame
    @test size(frame) == (4, 4)
    step_graph!(graph)
    first_frame = copy(frame)
    @test all(isfinite, frame)
    @test warmed_graph_step_allocation_bytes(graph) == 0
    @test @inferred(step_graph!(graph)) === graph
    @test frame != first_frame
    reset_graph!(graph)
    @test all(iszero, frame)
    step_graph!(graph)
    @test frame == first_frame

    deterministic = ccd_detector_acquisition_node(
        :detector;
        rows=4,
        columns=4,
        binning=2,
        pixel_scale_arcsec=0.1,
        wavelength_m=0.75e-6,
        exposure_duration_s=0.01,
        quantum_efficiency=0.5,
        rng_seed=1,
        photon_rate_schema="test.graph.shwfs-photon-rate.f32/1",
        frame_schema="test.graph.shwfs-frame.f32/1",
        photon_noise=false,
    )
    deterministic_graph = prepare_algorithm_graph(algorithm_graph(
        (deterministic,);
        inputs=(
            graph_input(
                :photon_rate,
                :detector => :photon_rate,
                photon_rate,
            ),
        ),
        outputs=(graph_output(:frame, :detector => :frame),),
    ))
    step_graph!(deterministic_graph)
    @test graph_output(deterministic_graph, Val(:frame)) ==
        fill(20.0f0, 2, 2)

    @test_throws AlgorithmGraphError ccd_detector_acquisition_node(
        :invalid_detector;
        rows=4,
        columns=4,
        pixel_scale_arcsec=0.1,
        wavelength_m=0.75e-6,
        exposure_duration_s=0.01,
        quantum_efficiency=0.5,
        rng_seed=1,
        photon_rate_schema="test.graph.shwfs-photon-rate.f32/1",
        frame_schema="test.graph.shwfs-frame.f32/1",
        photon_noise=false,
        readout_noise=false,
        readout_noise_e=2.7,
    )
    invalid_rate = copy(photon_rate)
    invalid_rate[1, 1] = -1
    invalid_definition = algorithm_graph(
        (node,);
        inputs=(
            graph_input(
                :photon_rate,
                :detector => :photon_rate,
                invalid_rate,
            ),
        ),
    )
    @test_throws InvalidConfiguration prepare_algorithm_graph(
        invalid_definition,
    )
end

@testset "global-shutter CMOS detector-acquisition graph node" begin
    photon_rate = fill(1_000.0f0, 4, 4)
    node = cmos_detector_acquisition_node(
        :detector;
        rows=4,
        columns=4,
        pixel_scale_arcsec=0.1,
        wavelength_m=0.8e-6,
        exposure_duration_s=0.01,
        quantum_efficiency=0.5,
        gain=2,
        dark_current_e_per_pixel_s=0,
        bits=12,
        full_well_e=1_000,
        photon_noise=false,
        readout_noise=false,
        readout_noise_e=0,
        column_readout_noise_e=0,
        row_readout_noise_e=0,
        rng_seed=0x1834,
        photon_rate_schema="test.graph.shwfs-photon-rate.f32/1",
        frame_schema="test.graph.shwfs-frame.f32/1",
    )
    definition = algorithm_graph(
        (node,);
        name=:cmos_detector_acquisition,
        inputs=(
            graph_input(
                :photon_rate,
                :detector => :photon_rate,
                photon_rate,
            ),
        ),
        outputs=(graph_output(:frame, :detector => :frame),),
    )
    graph = prepare_algorithm_graph(definition)
    owner = prepared_graph_node(graph, Val(:detector))
    frame = graph_output(graph, Val(:frame))

    @test owner.prepared.input.values === photon_rate
    @test owner.output === frame
    @test owner.prepared.detector.params.sensor isa
        AdaptiveOpticsSim.Detectors.CMOSSensor
    @test owner.prepared.detector.params.timing_model isa
        AdaptiveOpticsSim.Detectors.GlobalShutter
    step_graph!(graph)
    @test frame == fill(41.0f0, 4, 4)
    @test warmed_graph_step_allocation_bytes(graph) == 0
    @test @inferred(step_graph!(graph)) === graph
    reset_graph!(graph)
    @test all(iszero, frame)
    step_graph!(graph)
    @test frame == fill(41.0f0, 4, 4)

    @test_throws AlgorithmGraphError cmos_detector_acquisition_node(
        :invalid_detector;
        rows=4,
        columns=4,
        pixel_scale_arcsec=0.1,
        wavelength_m=0.8e-6,
        exposure_duration_s=0.01,
        quantum_efficiency=0.5,
        bits=12,
        rng_seed=1,
        photon_rate_schema="test.graph.shwfs-photon-rate.f32/1",
        frame_schema="test.graph.shwfs-frame.f32/1",
    )
end

@testset "complete-frame EMCCD detector-acquisition graph node" begin
    photon_rate = fill(1_000.0f0, 4, 4)
    node = emccd_detector_acquisition_node(
        :detector;
        rows=4,
        columns=4,
        normalized_pupil_sampling=0.05,
        wavelength_m=0.55e-6,
        exposure_duration_s=0.01,
        quantum_efficiency=0.95,
        gain=1,
        dark_current_e_per_pixel_s=20,
        bits=14,
        full_well_e=16_000,
        photon_noise=true,
        readout_noise=true,
        readout_noise_e=1,
        rng_seed=0x2234,
        photon_rate_schema="test.graph.pwfs-photon-rate.f32/1",
        frame_schema="test.graph.pwfs-frame.f32/1",
    )
    definition = algorithm_graph(
        (node,);
        name=:emccd_detector_acquisition,
        inputs=(
            graph_input(
                :photon_rate,
                :detector => :photon_rate,
                photon_rate,
            ),
        ),
        outputs=(graph_output(:frame, :detector => :frame),),
    )
    graph = prepare_algorithm_graph(definition)
    owner = prepared_graph_node(graph, Val(:detector))
    frame = graph_output(graph, Val(:frame))

    @test owner.prepared.input.values === photon_rate
    @test owner.output === frame
    @test size(frame) == (4, 4)
    step_graph!(graph)
    first_frame = copy(frame)
    @test all(isfinite, frame)
    @test sum(frame) > 0
    @test warmed_graph_step_allocation_bytes(graph) == 0
    @test @inferred(step_graph!(graph)) === graph
    @test frame != first_frame
    reset_graph!(graph)
    @test all(iszero, frame)
    step_graph!(graph)
    @test frame == first_frame

    @test_throws AlgorithmGraphError emccd_detector_acquisition_node(
        :invalid_detector;
        rows=4,
        columns=4,
        normalized_pupil_sampling=0.05,
        wavelength_m=0.55e-6,
        exposure_duration_s=0.01,
        quantum_efficiency=0.95,
        gain=10,
        em_gain_min=1,
        em_gain_max=5,
        rng_seed=1,
        photon_rate_schema="test.graph.pwfs-photon-rate.f32/1",
        frame_schema="test.graph.pwfs-frame.f32/1",
    )
end


@testset "declarative detector selectors preserve native acquisitions" begin
    photon_rate = fill(1_000.0f0, 2, 3)
    common = (
        rows=2,
        columns=3,
        wavelength_m=0.75e-6,
        exposure_duration_s=0.01,
        quantum_efficiency=0.9,
        photon_noise=true,
        readout_noise=true,
        readout_noise_e=1.0,
        rng_seed=42,
        photon_rate_schema="test.detector.photon-rate/1",
        frame_schema="test.detector.frame/1",
    )
    families = (
        ("ccd", ccd_detector_acquisition_node, (; pixel_scale_arcsec=0.1)),
        ("cmos", cmos_detector_acquisition_node, (; pixel_scale_arcsec=0.1)),
        ("emccd", emccd_detector_acquisition_node,
            (; normalized_pupil_sampling=0.05, gain=10, em_gain_max=10)),
    )
    for (detector, constructor, response) in families
        config = merge(common, response)
        entry = Dict{String,Any}(
            "name" => "camera",
            "type" => "detector",
            "detector" => detector,
            "config" => Dict{String,Any}(string(k) => v for (k, v) in pairs(config)),
        )
        document = Dict{String,Any}(
            "schema_version" => 1,
            "name" => "detector_selection",
            "nodes" => [entry],
            "inputs" => [Dict("name" => "photons",
                "destination" => "camera.photon_rate", "binding" => "photons")],
            "outputs" => [Dict("name" => "frame", "source" => "camera.frame")],
        )
        native = algorithm_graph(
            (constructor(:camera; config...),);
            inputs=(graph_input(:photons, :camera => :photon_rate, photon_rate),),
            outputs=(graph_output(:frame, :camera => :frame),),
        )
        for explicit_shutter in (false, true)
            explicit_shutter && (entry["config"]["shutter"] = "global")
            loaded = load_algorithm_graph(document; bindings=(; photons=photon_rate))
            actual = prepare_algorithm_graph(loaded)
            expected = prepare_algorithm_graph(native)
            for _ in 1:2
                step_graph!(actual)
                step_graph!(expected)
                @test graph_output(actual, Val(:frame)) == graph_output(expected, Val(:frame))
            end
            @test size(graph_output(actual, Val(:frame))) == (2, 3)
            reset_graph!(actual)
            reset_graph!(expected)
            step_graph!(actual)
            step_graph!(expected)
            @test graph_output(actual, Val(:frame)) == graph_output(expected, Val(:frame))
        end
        for shutter in ("rolling", "frame_transfer", "invalid", 1)
            entry["config"]["shutter"] = shutter
            @test_throws AlgorithmGraphError load_algorithm_graph(document; bindings=(; photons=photon_rate))
        end
        entry["config"]["shutter"] = "global"
        for unsupported in ("hgcdte", "invalid", 1)
            entry["detector"] = unsupported
            @test_throws AlgorithmGraphError load_algorithm_graph(document; bindings=(; photons=photon_rate))
        end
        delete!(entry, "detector")
        @test_throws AlgorithmGraphError load_algorithm_graph(document; bindings=(; photons=photon_rate))
        entry["detector"] = detector
        entry["type"] = "$(detector)_detector_acquisition_f32"
        @test_throws AlgorithmGraphError load_algorithm_graph(document; bindings=(; photons=photon_rate))
        entry["type"] = "modal_opd_expansion_f32"
        @test_throws AlgorithmGraphError load_algorithm_graph(document; bindings=(; photons=photon_rate))
    end
end

@testset "YAML InGaAs acquisition and calibrated pixel MTF" begin
    photons = reshape(Float32.(20:34), 3, 5)
    calibration = pixel_mtf_transfer((7, 9); box_width_pixels=1.03,
        diffusion_variance_pixels2=0.09, ipc_per_neighbor=0.03, T=Float32)
    path = joinpath(dirname(@__DIR__), "graph_files", "ingaas_pixel_mtf.yaml")
    document = YAML.load_file(path)
    bindings = (; photons, calibrated_transfer=copy(calibration))
    loaded = load_algorithm_graph(path; bindings)
    native = algorithm_graph((
        pixel_mtf_node(:response; input_shape=(3, 5), padding=2, T=Float32,
            input_schema="test.ingaas.photon-rate/1",
            output_schema="test.ingaas.photon-rate/1", transfer=calibration),
        ingaas_detector_acquisition_node(:camera; rows=3, columns=5,
            pixel_sampling_m=(12e-6, 14e-6), spectral_channel=:test_ingaas_band,
            exposure_duration_s=0.5, quantum_efficiency=0.75, gain=2,
            dark_current_e_per_pixel_s=2, glow_e_per_pixel_s=3,
            bits=12, full_well_e=1000, photon_noise=true, readout_noise=true,
            readout_noise_e=0.7, persistence_model=ExponentialPersistence(0.25, 0.5),
            nonlinearity_model=SaturatingFrameNonlinearity(0.01), rng_seed=43969,
            photon_rate_schema="test.ingaas.photon-rate/1",
            frame_schema="test.ingaas.frame/1")); name=:ingaas_pixel_mtf,
        inputs=(graph_input(:photons, :response => :input, photons),),
        links=(link(:response => :output, :camera => :photon_rate),),
        outputs=(graph_output(:filtered, :response => :output),
            graph_output(:frame, :camera => :frame)))
    @test typeof(loaded) === typeof(native)
    # Both constructors snapshot the binding before graph preparation.
    fill!(bindings.calibrated_transfer, Float32(NaN))
    actual, expected = prepare_algorithm_graph(loaded), prepare_algorithm_graph(native)
    owner = prepared_graph_node(actual, Val(:camera))
    @test owner.prepared.input.metadata.sampling == (12f-6, 14f-6)
    @test owner.prepared.input.metadata.spectral == IntegratedSpectralChannel(:test_ingaas_band)
    @test size(graph_output(actual, Val(:frame))) == (3, 5)
    first_frame = nothing
    for factor in (1, 2, 0)
        photons .= reshape(Float32.(20:34), 3, 5) .* factor
        step_graph!(actual)
        step_graph!(expected)
        @test graph_output(actual, Val(:filtered)) == graph_output(expected, Val(:filtered))
        @test graph_output(actual, Val(:frame)) == graph_output(expected, Val(:frame))
        @test owner.latent_buffer == prepared_graph_node(expected, Val(:camera)).latent_buffer
        factor == 1 && (first_frame = copy(graph_output(actual, Val(:frame))))
    end
    @test any(>(0), owner.latent_buffer)
    reset_graph!(actual)
    reset_graph!(expected)
    @test all(iszero, owner.latent_buffer)
    @test all(iszero, graph_output(actual, Val(:frame)))
    @test all(iszero, graph_output(actual, Val(:filtered)))
    photons .= reshape(Float32.(20:34), 3, 5)
    step_graph!(actual)
    step_graph!(expected)
    @test graph_output(actual, Val(:frame)) == first_frame
    @test graph_output(actual, Val(:frame)) == graph_output(expected, Val(:frame))

    valid_bindings = (; photons, calibrated_transfer=calibration)
    # Null models may be omitted or declared explicitly, with identical frames.
    null_document = deepcopy(document)
    config = null_document["nodes"][2]["config"]
    config["persistence_model"] = Dict("mode" => "none")
    config["nonlinearity_model"] = Dict("mode" => "none")
    explicit_null = prepare_algorithm_graph(load_algorithm_graph(null_document; bindings=valid_bindings))
    delete!(config, "persistence_model")
    delete!(config, "nonlinearity_model")
    delete!(config, "shutter")
    implicit_null = prepare_algorithm_graph(load_algorithm_graph(null_document; bindings=valid_bindings))
    step_graph!(explicit_null)
    step_graph!(implicit_null)
    @test graph_output(explicit_null, Val(:frame)) == graph_output(implicit_null, Val(:frame))

    for (field, invalid) in (
            ("pixel_sampling_m", [12e-6]), ("pixel_sampling_m", [12e-6, 0]),
            ("spectral_channel", 1), ("shutter", "rolling"),
            ("photon_noise", 1), ("bits", true), ("full_well_e", -1),
            ("quantum_efficiency", 1.1), ("glow_e_per_pixel_s", -1),
            ("persistence_model", "exponential"),
            ("persistence_model", Dict("mode" => "unknown")),
            ("persistence_model", Dict("mode" => "exponential", "coupling" => 0.2)),
            ("persistence_model", Dict("mode" => "none", "decay" => 0.2)),
            ("persistence_model", Dict("mode" => "exponential", "coupling" => -0.1, "decay" => 0.5)),
            ("persistence_model", Dict("mode" => "exponential", "coupling" => 0.1, "decay" => 1.1)),
            ("persistence_model", Dict("mode" => "exponential", "coupling" => true, "decay" => 0.5)),
            ("nonlinearity_model", "saturating"),
            ("nonlinearity_model", Dict("mode" => "unknown")),
            ("nonlinearity_model", Dict("mode" => "saturating")),
            ("nonlinearity_model", Dict("mode" => "none", "coefficient" => 0)),
            ("nonlinearity_model", Dict("mode" => "saturating", "coefficient" => -0.1)),
            ("nonlinearity_model", Dict("mode" => "saturating", "coefficient" => Inf)),
            ("wavelength_m", 1.5e-6))
        invalid_document = deepcopy(document)
        invalid_document["nodes"][2]["config"][field] = invalid
        @test_throws AlgorithmGraphError load_algorithm_graph(invalid_document; bindings=valid_bindings)
    end
    @test_throws AlgorithmGraphError load_algorithm_graph(document; bindings=(; photons))
    for transfer in (zeros(Float32, 9, 7), ones(Float32, 7),
            fill(2f0, 7, 9), fill(Float32(NaN), 7, 9), ones(ComplexF32, 7, 9))
        @test_throws AlgorithmGraphError load_algorithm_graph(document;
            bindings=(; photons, calibrated_transfer=transfer))
    end
    for (field, invalid) in (("input_shape", [5, 3]), ("input_shape", [3]),
            ("padding", -1), ("transfer", "calibrated_transfer"),
            ("file", "calibration.fits"))
        invalid_document = deepcopy(document)
        invalid_document["nodes"][1]["config"][field] = invalid
        @test_throws AlgorithmGraphError load_algorithm_graph(invalid_document; bindings=valid_bindings)
    end
    for invalid in ("missing", [[1, 1], [1, 1]], "bad.name")
        invalid_document = deepcopy(document)
        invalid_document["nodes"][1]["bindings"]["transfer"] = invalid
        @test_throws AlgorithmGraphError load_algorithm_graph(invalid_document; bindings=valid_bindings)
    end
    invalid_document = deepcopy(document)
    delete!(invalid_document["nodes"][1], "bindings")
    @test_throws AlgorithmGraphError load_algorithm_graph(invalid_document; bindings=valid_bindings)
end

@testset "cold node factories declare configuration binding names" begin
    factories = builtin_graph_node_types()
    config = (input_shape=(3, 5), padding=2,
        input_schema="test.photon-rate/1", output_schema="test.photon-rate/1",
        transfer="calibrated_transfer")
    @test AOG.graph_node_bindings(factories.pixel_mtf, config) == (:calibrated_transfer,)
    renamed = merge(config, (; transfer="measured_response"))
    @test AOG.graph_node_bindings(factories.pixel_mtf, renamed) == (:measured_response,)
    @test AOG.graph_node_bindings(factories.detector, (;)) == ()
    custom(name, config, props) = nothing
    @test AOG.graph_node_bindings(custom, (;)) == ()
    unbound = Base.structdiff(config, NamedTuple{(:transfer,)})
    @test AOG.graph_node_configuration(factories.pixel_mtf, unbound,
        (; transfer="calibrated_transfer")) == config
    @test AOG.graph_node_configuration(factories.detector, (;)) == (;)
    @test_throws AlgorithmGraphError AOG.graph_node_configuration(factories.pixel_mtf, config)
    @test_throws AlgorithmGraphError AOG.graph_node_configuration(factories.pixel_mtf, config,
        (; transfer="calibrated_transfer"))
    @test_throws AlgorithmGraphError AOG.graph_node_configuration(factories.pixel_mtf, unbound,
        (; transfer="calibrated_transfer", extra="unused"))
    @test_throws AlgorithmGraphError AOG.graph_node_configuration(factories.detector, (;), (; extra="unused"))
    for invalid in (1, "invalid.name", (1, 2))
        @test_throws AlgorithmGraphError AOG.graph_node_bindings(
            factories.pixel_mtf, merge(config, (; transfer=invalid)))
    end
    @test_throws AlgorithmGraphError AOG.graph_node_bindings(
        factories.pixel_mtf, Base.structdiff(config, NamedTuple{(:transfer,)}))
    @test_throws AlgorithmGraphError AOG.graph_node_bindings(
        factories.pixel_mtf, merge(config, (; unknown=true)))
end

@testset "YAML graph files compile to native graph definitions" begin
    normalized_yaml =
        AdaptiveOpticsSim.AlgorithmGraphs._normalize_graph_value(
            Any[1, "two", Any[true, 3.0]],
        )
    @test normalized_yaml == (1, "two", (true, 3.0))
    @test isconcretetype(typeof(normalized_yaml))

    coefficients = Float32[0.1, 0.2]
    basis = zeros(Float32, 2, 2, 2)
    fill!(@view(basis[:, :, 1]), 1.0f0)
    fill!(@view(basis[:, :, 2]), 2.0f0)
    pupil_support = Bool[true false; true true]
    path = joinpath(
        dirname(@__DIR__),
        "graph_files",
        "native_modal_opd.yaml",
    )

    definition = load_algorithm_graph(
        path;
        bindings=(;
            coefficients,
            basis,
            pupil_support,
        ),
    )
    @test graph_name(definition) === :native_modal_opd
    @test isconcretetype(typeof(definition))
    @test all(isconcretetype, fieldtypes(typeof(definition)))
    graph = prepare_algorithm_graph(definition)
    step_graph!(graph)
    @test graph_output(graph, Val(:opd)) ≈ Float32[0.5 0.0; 0.5 0.5]
    step_graph!(graph)
    @test @allocated(step_graph!(graph)) == 0
    @test @inferred(step_graph!(graph)) === graph

    @test_throws AlgorithmGraphError load_algorithm_graph(
        path;
        bindings=(; coefficients, basis),
    )
    valid_yaml = read(path, String)
    merged_yaml = replace(valid_yaml, "  config:\n" =>
        "  config:\n    <<: &defaults\n      mode_count: 1\n")
    mktemp() do merged_path, io
        write(io, merged_yaml)
        close(io)
        merged_definition = load_algorithm_graph(
            merged_path;
            bindings=(; coefficients, basis, pupil_support),
        )
        @test typeof(merged_definition) === typeof(definition)
        merged_graph = prepare_algorithm_graph(merged_definition)
        step_graph!(merged_graph)
        @test graph_output(merged_graph, Val(:opd)) == graph_output(graph, Val(:opd))
    end
    for invalid in (
        replace(valid_yaml, "schema_version: 1" =>
            "schema_version: 1\nschema_version: 1"),
        replace(valid_yaml, "- name: modal_opd" =>
            "- name: modal_opd\n  name: modal_opd"),
        replace(merged_yaml, "    mode_count: 2" =>
            "    mode_count: 2\n    mode_count: 2"),
        valid_yaml * "---\n" * valid_yaml,
    )
        mktemp() do invalid_path, io
            write(io, invalid)
            close(io)
            @test_throws AlgorithmGraphError load_algorithm_graph(
                invalid_path;
                bindings=(; coefficients, basis, pupil_support),
            )
        end
    end

    @test keys(builtin_graph_node_types()) == (
        :proper_optical_train_f64,
        :proper_optical_train_f32,
        :detector,
        :pixel_mtf,
        :deformable_mirror_surface_f32,
        :gaussian_deformable_mirror_surface_f32,
        :grid_gaussian_deformable_mirror_surface_f32,
        :modal_opd_expansion_f32,
        :second_order_response_f32,
        :multilayer_atmosphere_opd_f32,
        :spatial_ao_residual_f32,
        :pupil_opd_composition_f32,
        :pyramid_rate_f32,
        :ideal_pupil_relay_f32,
        :ideal_pupil_relay_f64,
        :shack_hartmann_rate_f32,
    )

    for invalid in (
        "schema_version: 2\nname: invalid_version\nnodes: []\n",
        "schema_version: 1\nname: unknown_field\nunexpected: true\nnodes: []\n",
        "schema_version: 1\nname: unknown_type\nnodes:\n" *
            "  - name: node\n    type: not_registered\n    config:\n      extent: 2\n",
        "schema_version: [unterminated\n",
        "- not_a_root_mapping\n",
        "null\n",
        "# empty document\n",
    )
        mktemp() do invalid_path, io
            write(io, invalid)
            close(io)
            @test_throws AlgorithmGraphError load_algorithm_graph(invalid_path)
        end
    end

    mktemp() do custom_path, io
        write(io, """
        schema_version: 1
        name: custom_type
        nodes:
          - name: source
            type: test_source
            config:
              extent: 2
              value: 7.0
        outputs:
          - name: sample
            source: source.output
        """)
        close(io)
        test_source = function (name, config, props)
            isempty(props) || error("test-source props must be empty")
            return algorithm_node(
                name,
                GraphTestSourceDeclaration,
                GraphTestSourceConfiguration(
                    Int(config.extent),
                    Float32(config.value),
                ),
            )
        end
        node_types = merge(
            builtin_graph_node_types(),
            (; test_source),
        )
        root = Dict{String,Any}(
            "schema_version" => 1,
            "name" => "custom_type",
            "nodes" => [Dict{String,Any}(
                "name" => "source",
                "type" => "test_source",
                "config" => Dict{String,Any}("extent" => 2, "value" => 7.0),
            )],
            "outputs" => [Dict{String,Any}(
                "name" => "sample",
                "source" => "source.output",
            )],
        )
        file_definition = load_algorithm_graph(custom_path; node_types)
        dictionary_definition = load_algorithm_graph(root; node_types, bindings=(; unused=42))
        @test typeof(dictionary_definition) === typeof(file_definition)
        @test graph_name(dictionary_definition) === graph_name(file_definition)

        custom = prepare_algorithm_graph(dictionary_definition)
        file_graph = prepare_algorithm_graph(file_definition)
        step_graph!(custom)
        step_graph!(file_graph)
        @test graph_output(custom, Val(:sample)) == Float32[7, 7]
        @test graph_output(custom, Val(:sample)) ==
            graph_output(file_graph, Val(:sample))

        bound_source(name, config, props, bindings) = test_source(name,
            merge(config, (; value=bindings.initial_value)), props)
        bound_definition = load_algorithm_graph(root;
            node_types=merge(node_types, (; test_source=bound_source)),
            bindings=(; initial_value=9.0))
        bound_graph = prepare_algorithm_graph(bound_definition)
        step_graph!(bound_graph)
        @test graph_output(bound_graph, Val(:sample)) == Float32[9, 9]

        invalid_root = deepcopy(root)
        invalid_root["schema_version"] = 2
        @test_throws AlgorithmGraphError load_algorithm_graph(
            invalid_root;
            node_types,
        )
        invalid_root = deepcopy(root)
        invalid_root["unexpected"] = true
        @test_throws AlgorithmGraphError load_algorithm_graph(
            invalid_root;
            node_types,
        )
        invalid_root = deepcopy(root)
        invalid_root["nodes"] = ["not a node table"]
        @test_throws AlgorithmGraphError load_algorithm_graph(
            invalid_root;
            node_types,
        )
    end
end

@testset "YAML graph binds a measured deformable mirror" begin
    pdm_command = Float32[2, 3]
    actuator_coordinates = Float32[-0.5 0.5; 0 0]
    influence_functions = Float32[
        1 0
        0 1
        0 1
        1 0
    ]
    path = joinpath(
        dirname(@__DIR__),
        "graph_files",
        "measured_deformable_mirror.yaml",
    )
    definition = load_algorithm_graph(
        path;
        bindings=(;
            pdm_command,
            actuator_coordinates,
            influence_functions,
        ),
    )
    graph = prepare_algorithm_graph(definition)
    step_graph!(graph)
    @test graph_output(graph, Val(:surface_opd)) == Float32[2 3; 3 2]
    @test warmed_graph_step_allocation_bytes(graph) == 0
    @test @inferred(step_graph!(graph)) === graph
end

@testset "fully declared HIL reference systems" begin
    @test HILReferenceSystems.supported_wavefront_sensors() ==
        (:shack_hartmann, :pyramid)
    @test HILReferenceSystems.actuator_count() == 25
    @test basename(HILReferenceSystems.graph_path(:shack_hartmann)) ==
        "shack_hartmann_hil_reference.yaml"
    @test basename(HILReferenceSystems.graph_path(:pyramid)) ==
        "pyramid_hil_reference.yaml"
    @test_throws ArgumentError HILReferenceSystems.graph_path(:unknown)

    coordinates = HILReferenceSystems.actuator_coordinates(Float32)
    @test size(coordinates) == (2, 25)
    @test coordinates[:, 1] == Float32[-0.8, -0.8]
    @test coordinates[:, 13] == Float32[0, 0]
    @test coordinates[:, 25] == Float32[0.8, 0.8]
    pitch = coordinates[1, 2] - coordinates[1, 1]
    width = HILReferenceSystems.influence_width(Float32)
    @test exp(-(pitch * pitch) / (2 * width * width)) ≈ 0.3f0

    valid_subapertures =
        HILReferenceSystems.shack_hartmann_valid_subapertures()
    lenslet_order = HILReferenceSystems.shack_hartmann_lenslet_order()
    @test size(valid_subapertures) == (8, 8)
    @test count(valid_subapertures) == 52
    @test lenslet_order == UInt32.(findall(vec(valid_subapertures)))
    @test size(HILReferenceSystems.shack_hartmann_reference_signal()) ==
        (64, 2)

    cases = (
        (
            wavefront_sensor=:shack_hartmann,
            graph_name=:shack_hartmann_hil_reference,
            frame_shape=(64, 64),
        ),
        (
            wavefront_sensor=:pyramid,
            graph_name=:pyramid_hil_reference,
            frame_shape=(36, 36),
        ),
    )
    for case in cases
        prepared = HILReferenceSystems.prepare_hil_reference_system(
            case.wavefront_sensor,
        )
        graph = prepared.graph
        boundary = prepared.boundary
        @test graph_name(graph) === case.graph_name
        @test graph_input(graph, Val(:uncompensated_opd)) ===
            prepared.uncompensated_opd
        @test all(iszero, prepared.uncompensated_opd)

        sequence = step_hil_frame!(boundary)
        flat_frame = copy(hil_frame_buffer(boundary))
        @test sequence == UInt64(1)
        @test size(flat_frame) == case.frame_shape
        @test all(isfinite, flat_frame)
        @test sum(flat_frame) > 0
        @test all(iszero, graph_output(graph, Val(:dm_surface_opd)))

        fill!(hil_command_buffer(boundary), 0.0f0)
        hil_command_buffer(boundary)[13] = 2.0f-8
        adopt_hil_command!(boundary, sequence)
        sequence = step_hil_frame!(boundary)
        surface_opd = graph_output(graph, Val(:dm_surface_opd))
        @test maximum(surface_opd) ≈ 2.0f-8 rtol = 5.0f-3
        @test minimum(surface_opd) >= 0.0f0
        @test graph_output(graph, Val(:pupil_opd)) ≈ surface_opd
        @test hil_frame_buffer(boundary) != flat_frame

        fill!(hil_command_buffer(boundary), 0.0f0)
        adopt_hil_command!(boundary, sequence)
        @test @allocated(step_hil_frame!(boundary)) == 0
        @test all(iszero, graph_output(graph, Val(:dm_surface_opd)))
    end

    for case in cases
        prepared =
            HILReferenceSystems.prepare_atmospheric_hil_reference_system(
                case.wavefront_sensor;
                atmosphere_step=1.0e-3,
                rng_seed=7,
            )
        graph = prepared.graph
        boundary = prepared.boundary
        @test graph_name(graph) === Symbol(
            case.wavefront_sensor,
            "_atmosphere_hil_reference",
        )
        @test size(hil_frame_buffer(boundary)) == case.frame_shape

        sequence = step_hil_frame!(boundary)
        first_atmosphere = copy(graph_output(
            graph,
            Val(:atmosphere_opd),
        ))
        first_frame = copy(hil_frame_buffer(boundary))
        @test all(isfinite, first_atmosphere)
        @test !all(iszero, first_atmosphere)
        @test all(isfinite, first_frame)
        @test sum(first_frame) > 0

        fill!(hil_command_buffer(boundary), 0.0f0)
        adopt_hil_command!(boundary, sequence)
        @test @allocated(step_hil_frame!(boundary)) == 0
        @test graph_output(graph, Val(:atmosphere_opd)) != first_atmosphere

        reset_hil_boundary!(boundary)
        sequence = step_hil_frame!(boundary)
        @test graph_output(graph, Val(:atmosphere_opd)) == first_atmosphere
        @test hil_frame_buffer(boundary) == first_frame
        fill!(hil_command_buffer(boundary), 0.0f0)
        adopt_hil_command!(boundary, sequence)
    end

    atmospheric =
        HILReferenceSystems.prepare_atmospheric_hil_reference_system(
            :shack_hartmann;
            atmosphere_step=1.0e-3,
            rng_seed=11,
        )
    diagnostics =
        HILReferenceSystems.prepare_hil_science_diagnostics()
    @test size(HILReferenceSystems.open_loop_psf(diagnostics)) == (128, 128)
    @test size(HILReferenceSystems.closed_loop_psf(diagnostics)) == (128, 128)
    @test HILReferenceSystems.open_loop_on_axis_strehl(diagnostics) ≈ 1.0f0
    @test HILReferenceSystems.closed_loop_on_axis_strehl(diagnostics) ≈ 1.0f0
    @test maximum(HILReferenceSystems.open_loop_psf(diagnostics)) < 1.0f0
    @test maximum(HILReferenceSystems.closed_loop_psf(diagnostics)) < 1.0f0

    sequence = step_hil_frame!(atmospheric.boundary)
    atmosphere_opd = graph_output(atmospheric.graph, Val(:atmosphere_opd))
    residual_opd = graph_output(atmospheric.graph, Val(:pupil_opd))
    HILReferenceSystems.update_hil_science_diagnostics!(
        diagnostics,
        atmosphere_opd,
        residual_opd,
    )
    open_loop_psf = HILReferenceSystems.open_loop_psf(diagnostics)
    closed_loop_psf = HILReferenceSystems.closed_loop_psf(diagnostics)
    @test all(isfinite, open_loop_psf)
    @test all(isfinite, closed_loop_psf)
    @test open_loop_psf ≈ closed_loop_psf
    @test maximum(open_loop_psf) <= 1.001f0
    @test maximum(closed_loop_psf) <= 1.001f0
    @test 0.0f0 <
        HILReferenceSystems.open_loop_on_axis_strehl(diagnostics) < 1.0f0
    @test HILReferenceSystems.closed_loop_on_axis_strehl(diagnostics) ≈
        HILReferenceSystems.open_loop_on_axis_strehl(diagnostics)
    @test @allocated(
        HILReferenceSystems.update_hil_science_diagnostics!(
            diagnostics,
            atmosphere_opd,
            residual_opd,
        ),
    ) == 0
    fill!(hil_command_buffer(atmospheric.boundary), 0.0f0)
    adopt_hil_command!(atmospheric.boundary, sequence)
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
    packed_parent = zeros(Float32, 3)
    wrapped_input = @view packed_parent[1:2]
    @test_throws AlgorithmGraphError prepare_algorithm_graph(algorithm_graph(
        (node,);
        inputs=(graph_input(:input, :add => :input, wrapped_input),),
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
    reordered = prepare_algorithm_graph(backward)
    @test keys(reordered.nodes) == (:second, :first)
    step_graph!(reordered)
    @test prepared_graph_node(reordered, :first).output == Float32[2, 2]

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

@testset "fixed-step model-time driver" begin
    schedule = AlgorithmGraphs.PeriodicSchedule(
        ModelDuration(100);
        phase=ModelDuration(25),
    )
    driver = FixedStepModelTimeDriver(
        schedule;
        origin=ModelTimestamp(1_000),
    )

    @test model_nanoseconds(ModelTimestamp(1_000)) == 1_000
    @test model_nanoseconds(ModelDuration(100)) == 100
    @test model_time_seconds(ModelTimestamp(1_000)) == 1.0e-6
    @test model_duration_seconds(ModelDuration(100)) == 1.0e-7
    @test schedule_period(schedule) == ModelDuration(100)
    @test schedule_phase(schedule) == ModelDuration(25)
    @test schedule_timestamp(schedule, 2, ModelTimestamp(1_000)) ==
          ModelTimestamp(1_125)
    @test_throws AlgorithmGraphError ModelTimestamp(-1)
    @test_throws AlgorithmGraphError ModelDuration(true)
    @test_throws AlgorithmGraphError AlgorithmGraphs.PeriodicSchedule(
        ModelDuration(0),
    )

    @test model_time_sequence(driver) == UInt64(0)
    @test !model_time_exhausted(driver)
    @test next_model_timestamp(driver) == ModelTimestamp(1_025)
    @test advance_model_time!(driver) == ModelTimestamp(1_025)
    @test model_time_sequence(driver) == UInt64(1)
    @test next_model_timestamp(driver) == ModelTimestamp(1_125)
    @test @inferred(advance_model_time!(driver)) == ModelTimestamp(1_125)
    @test @allocated(advance_model_time!(driver)) == 0

    @test reset_model_time!(driver) === driver
    @test model_time_sequence(driver) == UInt64(0)
    @test next_model_timestamp(driver) == ModelTimestamp(1_025)
end

@testset "prepared-boundary model-time driver" begin
    driver = prepare_boundary_model_time_driver((
        ModelTimestamp(10),
        ModelTimestamp(20),
        ModelTimestamp(35),
    ))

    @test next_model_timestamp(driver) == ModelTimestamp(10)
    @test advance_model_time!(driver) == ModelTimestamp(10)
    @test @inferred(advance_model_time!(driver)) == ModelTimestamp(20)
    @test @allocated(advance_model_time!(driver)) == 0
    @test model_time_sequence(driver) == UInt64(3)
    @test model_time_exhausted(driver)
    @test_throws AlgorithmGraphError next_model_timestamp(driver)
    @test_throws AlgorithmGraphError advance_model_time!(driver)
    reset_model_time!(driver)
    @test next_model_timestamp(driver) == ModelTimestamp(10)

    @test_throws AlgorithmGraphError prepare_boundary_model_time_driver(())
    @test_throws AlgorithmGraphError prepare_boundary_model_time_driver((
        ModelTimestamp(10),
        ModelTimestamp(10),
    ))
    @test_throws AlgorithmGraphError prepare_boundary_model_time_driver((
        ModelTimestamp(20),
        ModelTimestamp(10),
    ))
    @test_throws AlgorithmGraphError prepare_boundary_model_time_driver((1, 2))
end

@testset "periodic prepared-boundary model time" begin
    schedule = AlgorithmGraphs.PeriodicSchedule(
        ModelDuration(1_000);
        phase=ModelDuration(100),
    )
    driver = prepare_boundary_model_time_driver(
        schedule,
        (ModelDuration(0), ModelDuration(250), ModelDuration(750)),
        2;
        origin=ModelTimestamp(10),
    )
    observed = ntuple(_ -> advance_model_time!(driver), 6)
    @test observed == (
        ModelTimestamp(110),
        ModelTimestamp(360),
        ModelTimestamp(860),
        ModelTimestamp(1_110),
        ModelTimestamp(1_360),
        ModelTimestamp(1_860),
    )
    @test model_time_exhausted(driver)

    @test_throws AlgorithmGraphError prepare_boundary_model_time_driver(
        schedule,
        (),
        1,
    )
    @test_throws AlgorithmGraphError prepare_boundary_model_time_driver(
        schedule,
        (ModelDuration(0), ModelDuration(1_000)),
        1,
    )
    @test_throws AlgorithmGraphError prepare_boundary_model_time_driver(
        schedule,
        (ModelDuration(250), ModelDuration(100)),
        1,
    )
    @test_throws AlgorithmGraphError prepare_boundary_model_time_driver(
        schedule,
        (ModelDuration(0),),
        0,
    )
    @test_throws AlgorithmGraphError prepare_boundary_model_time_driver(
        schedule,
        (ModelDuration(0),),
        true,
    )
end

@testset "captured model-time replay" begin
    captures = (
        CapturedModelTimestamp(
            ModelTimestamp(15),
            ModelDuration(2),
            (source=UInt32(7), sequence=UInt64(41)),
        ),
        CapturedModelTimestamp(
            ModelTimestamp(28),
            ModelDuration(3),
            (source=UInt32(7), sequence=UInt64(42)),
        ),
    )
    driver = prepare_captured_model_time_driver(captures)

    @test next_model_time_capture(driver) === captures[1]
    @test model_timestamp(next_model_time_capture(driver)) ==
        ModelTimestamp(15)
    @test model_time_uncertainty(next_model_time_capture(driver)) ==
        ModelDuration(2)
    @test model_time_provenance(next_model_time_capture(driver)) ==
        (source=UInt32(7), sequence=UInt64(41))
    @test advance_model_time!(driver) == ModelTimestamp(15)
    @test @inferred(next_model_time_capture(driver)) === captures[2]
    @test @allocated(advance_model_time!(driver)) == 0
    @test model_time_exhausted(driver)

    reset_model_time!(driver)
    @test next_model_timestamp(driver) == ModelTimestamp(15)
    @test_throws AlgorithmGraphError prepare_captured_model_time_driver(())
    @test_throws AlgorithmGraphError prepare_captured_model_time_driver((
        captures[1],
        CapturedModelTimestamp(
            ModelTimestamp(15),
            ModelDuration(3),
            (source=UInt32(7), sequence=UInt64(42)),
        ),
    ))
    @test_throws AlgorithmGraphError prepare_captured_model_time_driver((
        captures[1],
        CapturedModelTimestamp(
            ModelTimestamp(28),
            ModelDuration(3),
            (stream=UInt32(7), sequence=UInt64(42)),
        ),
    ))
    @test_throws AlgorithmGraphError CapturedModelTimestamp(
        ModelTimestamp(1),
        ModelDuration(0),
        "borrowed provenance",
    )
end

@testset "model-time graph stepping commits together" begin
    graph = prepare_algorithm_graph(algorithm_graph(
        (algorithm_node(
            :source,
            GraphTestSourceDeclaration,
            GraphTestSourceConfiguration(2, 4.0f0),
        ),);
        outputs=(graph_output(:sample, :source => :output),),
    ))
    driver = FixedStepModelTimeDriver(AlgorithmGraphs.PeriodicSchedule(ModelDuration(50)))

    @test step_graph_at!(graph, driver) == ModelTimestamp(0)
    @test graph_output(graph, Val(:sample)) == Float32[4, 4]
    @test graph_step_sequence(graph) == model_time_sequence(driver) == UInt64(1)
    @test @allocated(step_graph_at!(graph, driver)) == 0

    step_graph!(graph)
    @test_throws AlgorithmGraphError step_graph_at!(graph, driver)

    reset_graph!(graph)
    reset_model_time!(driver)
    @test step_graph_at!(graph, driver) == ModelTimestamp(0)

    reset_graph!(graph)
    captured_driver = prepare_captured_model_time_driver((
        CapturedModelTimestamp(
            ModelTimestamp(75),
            ModelDuration(4),
            (source=UInt32(1), sequence=UInt64(1)),
        ),
    ))
    @test step_graph_at!(graph, captured_driver) == ModelTimestamp(75)
    @test model_time_exhausted(captured_driver)

    failing_input = Float32[-1, 1]
    failing_graph = prepare_algorithm_graph(algorithm_graph(
        (algorithm_node(:checked, GraphTestFailureDeclaration, 2),);
        inputs=(graph_input(:input, :checked => :input, failing_input),),
        outputs=(graph_output(:output, :checked => :output),),
    ))
    failing_driver = FixedStepModelTimeDriver(AlgorithmGraphs.PeriodicSchedule(ModelDuration(50)))
    @test_throws DomainError step_graph_at!(failing_graph, failing_driver)
    @test model_time_sequence(failing_driver) == UInt64(0)
end


# Fault injection wraps the exact command backing array; a copy changes one
# target element and then fails. The graph keeps its original exact binding.
struct CalibrationTestFailingCommand{A<:AbstractVector{Float32}} <: AbstractVector{Float32}
    values::A
    fail::Base.RefValue{Bool}
end
Base.size(values::CalibrationTestFailingCommand) = size(values.values)
Base.IndexStyle(::Type{<:CalibrationTestFailingCommand}) = IndexLinear()
Base.getindex(values::CalibrationTestFailingCommand, index::Int) = values.values[index]
function Base.copyto!(destination::CalibrationTestFailingCommand, source::Array{Float32})
    if destination.fail[]
        destination.values[1] = source[1]
        throw(ErrorException("injected partial target copy"))
    end
    copyto!(destination.values, source)
    return destination
end

function calibration_test_graph(
    command=Float32[1, 2]; declaration=GraphTestExactBindingDeclaration,
)
    return prepare_algorithm_graph(algorithm_graph(
        (algorithm_node(:plant, declaration, 2),);
        inputs=(graph_input(:command, :plant => :input, command),),
        outputs=(graph_output(:frame, :plant => :output),),
    ))
end

function calibration_test_boundary(graph=calibration_test_graph(); kwargs...)
    return prepare_graph_calibration_boundary(
        graph; command_input=:command, frame_output=:frame, kwargs...,
    )
end

@testset "graph calibration held-probe boundary" begin
    graph = calibration_test_graph()
    command_buffer = zeros(Float32, 2)
    frame_buffer = zeros(Float32, 2)
    boundary = calibration_test_boundary(graph; command_buffer, frame_buffer)
    @test boundary isa PreparedGraphCalibrationBoundary
    @test hil_command_buffer(boundary) === command_buffer
    @test hil_frame_buffer(boundary) === frame_buffer
    @test boundary.command_input === graph_input(graph, Val(:command))
    @test boundary.frame_output === graph_output(graph, Val(:frame))
    @test @inferred(hil_boundary_status(boundary)) == (
        frame_sequence=UInt64(0), active_probe_sequence=UInt64(0), failed=false,
    )
    @test command_buffer == Float32[1, 2]
    @test all(iszero, frame_buffer)
    @test_throws AlgorithmGraphError step_hil_exposure!(boundary)
    @test_throws AlgorithmGraphError adopt_hil_probe!(boundary, UInt64(0))
    @test_throws AlgorithmGraphError adopt_hil_probe!(boundary, 1)
    @test_throws AlgorithmGraphError adopt_hil_probe!(boundary, true)
    @test graph_step_sequence(graph) == UInt64(0)

    command_buffer .= Float32[3, 4]
    @test @inferred(adopt_hil_probe!(boundary, UInt64(7))) == UInt64(7)
    @test graph_input(graph, Val(:command)) == Float32[3, 4]
    @test graph_step_sequence(graph) == UInt64(0) # Adoption is not an exposure.
    @test @inferred(step_hil_exposure!(boundary)) == UInt64(1)
    @test frame_buffer == Float32[3, 4]
    command_buffer .= Float32[20, 30] # Mailbox mutation cannot alter the held probe.
    @test @inferred(step_hil_exposure!(boundary)) == UInt64(2)
    @test frame_buffer == Float32[3, 4]
    @test hil_boundary_status(boundary).active_probe_sequence == UInt64(7)
    saved_status = hil_boundary_status(boundary)
    for sequence in (UInt64(0), UInt64(6), UInt64(7))
        @test_throws AlgorithmGraphError adopt_hil_probe!(boundary, sequence)
        @test hil_boundary_status(boundary) == saved_status
        @test graph_input(graph, Val(:command)) == Float32[3, 4]
    end
    for invalid in (NaN32, Inf32, -Inf32)
        command_buffer .= Float32[99, invalid] # The finite prefix is never copied.
        @test_throws AlgorithmGraphError adopt_hil_probe!(boundary, UInt64(8))
        @test graph_input(graph, Val(:command)) == Float32[3, 4]
        @test hil_boundary_status(boundary) == saved_status
    end
    command_buffer .= Float32[5, 6]
    @test @allocated(adopt_hil_probe!(boundary, UInt64(8))) == 0
    @test @allocated(step_hil_exposure!(boundary)) == 0
    @test frame_buffer == Float32[5, 6]
    @test @allocated(step_hil_exposure!(boundary)) == 0
    @test hil_boundary_status(boundary).frame_sequence == UInt64(4)
    @test @inferred(reset_hil_boundary!(boundary)) === boundary
    @test graph_input(graph, Val(:command)) == Float32[1, 2]
    @test command_buffer == Float32[1, 2]
    @test all(iszero, frame_buffer)
    @test hil_boundary_status(boundary) == (
        frame_sequence=UInt64(0), active_probe_sequence=UInt64(0), failed=false,
    )
    @test_throws AlgorithmGraphError step_hil_exposure!(boundary)
    @test adopt_hil_probe!(boundary, UInt64(1)) == UInt64(1)
    @test step_hil_exposure!(boundary) == UInt64(1)
    @test adopt_hil_probe!(boundary, typemax(UInt64)) == typemax(UInt64)
    @test_throws AlgorithmGraphError adopt_hil_probe!(boundary, typemax(UInt64))
    @test step_hil_exposure!(boundary) == UInt64(2)
end

@testset "graph calibration model-time exposures" begin
    boundary = calibration_test_boundary()
    driver = FixedStepModelTimeDriver(
        AlgorithmGraphs.PeriodicSchedule(ModelDuration(50)); origin=ModelTimestamp(100),
    )
    @test_throws AlgorithmGraphError step_hil_exposure_at!(boundary, driver)
    @test model_time_sequence(driver) == UInt64(0)
    adopt_hil_probe!(boundary, UInt64(9))
    @test @inferred(step_hil_exposure_at!(boundary, driver)) == (
        sequence=UInt64(1), timestamp=ModelTimestamp(100),
    )
    @test @inferred(step_hil_exposure_at!(boundary, driver)) == (
        sequence=UInt64(2), timestamp=ModelTimestamp(150),
    )
    @test @allocated(step_hil_exposure_at!(boundary, driver)) == 0
    @test model_time_sequence(driver) == UInt64(3)
    @test next_model_timestamp(driver) == ModelTimestamp(250)
    @test hil_boundary_status(boundary).active_probe_sequence == UInt64(9)
    @test @inferred(reset_hil_boundary!(boundary, driver)) === boundary
    @test model_time_sequence(driver) == UInt64(0)
    @test hil_boundary_status(boundary).active_probe_sequence == UInt64(0)
    @test_throws AlgorithmGraphError step_hil_exposure_at!(boundary, driver)
    adopt_hil_probe!(boundary, UInt64(2))
    advance_model_time!(driver)
    @test_throws AlgorithmGraphError step_hil_exposure_at!(boundary, driver)
    @test graph_step_sequence(boundary.graph) == UInt64(0)
    @test !hil_boundary_status(boundary).failed
    reset_hil_boundary!(boundary, driver)
    adopt_hil_probe!(boundary, UInt64(3))
    @test step_hil_exposure_at!(boundary, driver).timestamp == ModelTimestamp(100)

    finite_boundary = calibration_test_boundary()
    finite_driver = prepare_boundary_model_time_driver((ModelTimestamp(5), ModelTimestamp(20)))
    adopt_hil_probe!(finite_boundary, UInt64(100))
    @test step_hil_exposure_at!(finite_boundary, finite_driver).timestamp == ModelTimestamp(5)
    @test step_hil_exposure_at!(finite_boundary, finite_driver).timestamp == ModelTimestamp(20)
    @test model_time_exhausted(finite_driver)
    @test_throws AlgorithmGraphError step_hil_exposure_at!(finite_boundary, finite_driver)
    @test hil_boundary_status(finite_boundary).failed
    @test graph_step_sequence(finite_boundary.graph) == UInt64(2)
    reset_hil_boundary!(finite_boundary, finite_driver)
    @test !model_time_exhausted(finite_driver)
    @test_throws AlgorithmGraphError step_hil_exposure_at!(finite_boundary, finite_driver)
end

@testset "HIL boundary rejects shared Array storage before initialization" begin
    # Reshaping via a matrix makes distinct Array objects with shared storage.
    alias_vector(values) = reshape(reshape(values, 1, length(values)), length(values))
    for prepare in (prepare_graph_hil_boundary, prepare_graph_calibration_boundary)
        for role in (:command_command, :frame_frame, :command_frame, :frame_command, :exchange)
            graph = calibration_test_graph()
            graph_command = graph_input(graph, Val(:command))
            graph_frame = graph_output(graph, Val(:frame))
            fill!(graph_frame, 42.0f0)
            command_buffer = Float32[8, 9]
            frame_buffer = Float32[10, 11]
            if role === :command_command
                command_buffer = alias_vector(graph_command)
                @test command_buffer !== graph_command
            elseif role === :frame_frame
                frame_buffer = alias_vector(graph_frame)
                @test frame_buffer !== graph_frame
            elseif role === :command_frame
                command_buffer = alias_vector(graph_frame)
            elseif role === :frame_command
                frame_buffer = alias_vector(graph_command)
            else
                frame_buffer = alias_vector(command_buffer)
                @test command_buffer !== frame_buffer
            end
            saved_command = copy(command_buffer)
            saved_frame = copy(frame_buffer)
            @test_throws AlgorithmGraphError prepare(
                graph; command_input=:command, frame_output=:frame,
                command_buffer, frame_buffer,
            )
            @test graph_command == Float32[1, 2]
            @test graph_frame == Float32[42, 42]
            @test command_buffer == saved_command
            @test frame_buffer == saved_frame
        end
    end
end

@testset "graph calibration ownership and fail-stop" begin
    graph = calibration_test_graph()
    @test_throws AlgorithmGraphError calibration_test_boundary(
        graph; command_buffer=graph_input(graph, Val(:command)),
    )
    @test_throws AlgorithmGraphError calibration_test_boundary(
        graph; frame_buffer=graph_output(graph, Val(:frame)),
    )
    shared = zeros(Float32, 2)
    @test_throws AlgorithmGraphError calibration_test_boundary(
        graph; command_buffer=shared, frame_buffer=shared,
    )
    @test_throws AlgorithmGraphError calibration_test_boundary(graph; frame_buffer=zeros(Float32, 3))
    @test_throws AlgorithmGraphError calibration_test_boundary(graph; command_buffer=zeros(Float64, 2))
    @test_throws AlgorithmGraphError calibration_test_boundary(calibration_test_graph(Float32[NaN, 0]))
    ticket = step_graph_async!(graph)
    @test_throws AlgorithmGraphError calibration_test_boundary(graph)
    wait_graph_step!(ticket)
    @test_throws AlgorithmGraphError calibration_test_boundary(graph)

    pending = calibration_test_boundary()
    adopt_hil_probe!(pending, UInt64(1))
    ticket = step_graph_async!(pending.graph) # Deliberate caller ownership violation.
    @test_throws AlgorithmGraphError adopt_hil_probe!(pending, UInt64(2))
    @test_throws AlgorithmGraphError step_hil_exposure!(pending)
    @test graph_step_pending(pending.graph)
    wait_graph_step!(ticket)
    @test_throws AlgorithmGraphError step_hil_exposure!(pending)
    @test_throws AlgorithmGraphError adopt_hil_probe!(pending, UInt64(2))
    reset_hil_boundary!(pending)
    @test_throws AlgorithmGraphError step_hil_exposure!(pending)

    failing = calibration_test_boundary(calibration_test_graph(; declaration=GraphTestFailureDeclaration))
    hil_command_buffer(failing) .= Float32[1, -1]
    adopt_hil_probe!(failing, UInt64(4))
    @test_throws DomainError step_hil_exposure!(failing)
    @test hil_boundary_status(failing).failed
    @test hil_boundary_status(failing).frame_sequence == UInt64(0)
    @test all(iszero, hil_frame_buffer(failing))
    @test_throws AlgorithmGraphError adopt_hil_probe!(failing, UInt64(5))
    @test_throws AlgorithmGraphError step_hil_exposure!(failing)
    @test_throws AlgorithmGraphError calibration_test_boundary(failing.graph)
    reset_hil_boundary!(failing)
    @test_throws AlgorithmGraphError step_hil_exposure!(failing)
    adopt_hil_probe!(failing, UInt64(1))
    @test step_hil_exposure!(failing) == UInt64(1)

    prepared = calibration_test_boundary()
    fault = CalibrationTestFailingCommand(prepared.command_input, Ref(true))
    copy_failure = PreparedGraphCalibrationBoundary(
        prepared.graph, fault, prepared.frame_output, prepared.command_buffer,
        prepared.frame_buffer, prepared.initial_command, prepared.state,
    )
    hil_command_buffer(copy_failure) .= Float32[8, 9]
    @test_throws ErrorException adopt_hil_probe!(copy_failure, UInt64(2))
    @test fault.values == Float32[8, 2]
    @test hil_boundary_status(copy_failure).active_probe_sequence == UInt64(0)
    @test hil_boundary_status(copy_failure).failed
    @test_throws AlgorithmGraphError adopt_hil_probe!(copy_failure, UInt64(3))
    @test_throws AlgorithmGraphError step_hil_exposure!(copy_failure)
    @test_throws ErrorException reset_hil_boundary!(copy_failure)
    @test hil_boundary_status(copy_failure).failed
    fault.fail[] = false
    reset_hil_boundary!(copy_failure)
    @test fault.values == Float32[1, 2]
    @test !hil_boundary_status(copy_failure).failed
    @test_throws AlgorithmGraphError step_hil_exposure!(copy_failure)
end

function multi_command_hil_test_graph(;
    commands=(dm=Float32[1, 2], tiptilt=Float32[3], offload=Float32[4, 5, 6]),
)
    nodes = map(keys(commands), values(commands)) do name, command
        algorithm_node(name, GraphTestExactBindingDeclaration, length(command))
    end
    inputs = map(keys(commands), values(commands)) do name, command
        graph_input(name, name => :input, command)
    end
    return prepare_algorithm_graph(algorithm_graph(
        nodes; inputs, outputs=(graph_output(:frame, :dm => :output),),
    ))
end

function multi_command_hil_test_boundary(graph=multi_command_hil_test_graph(); kwargs...)
    return prepare_graph_hil_boundary(
        graph; command_input=(:dm, :tiptilt, :offload), frame_output=:frame, kwargs...,
    )
end

@testset "graph calibration holds complete named probes" begin
    graph = multi_command_hil_test_graph()
    prepare(graph; kwargs...) = prepare_graph_calibration_boundary(
        graph; command_input=(:dm, :tiptilt, :offload), frame_output=:frame, kwargs...,
    )
    buffers = (dm=zeros(Float32, 2), tiptilt=zeros(Float32, 1), offload=zeros(Float32, 3))
    boundary = prepare(graph; command_buffer=buffers)
    initial = map(copy, boundary.command_input)
    @test keys(hil_command_buffer(boundary)) == (:dm, :tiptilt, :offload)
    @test hil_command_buffer(boundary) === buffers
    @test buffers == initial
    @test all(isconcretetype, fieldtypes(typeof(boundary)))
    @test_throws AlgorithmGraphError step_hil_exposure!(boundary)
    foreach(buffer -> fill!(buffer, 9.0f0), buffers)
    for invalid in (NaN32, Inf32, -Inf32)
        buffers.offload[end] = invalid
        @test_throws AlgorithmGraphError adopt_hil_probe!(boundary, UInt64(7))
        @test boundary.command_input == initial
        @test !hil_boundary_status(boundary).failed
        @test hil_boundary_status(boundary).active_probe_sequence == 0
    end
    buffers.offload[end] = 9.0f0
    @test @inferred(adopt_hil_probe!(boundary, UInt64(7))) == UInt64(7)
    @test boundary.command_input == buffers
    @test graph_step_sequence(graph) == 0
    @test @inferred(step_hil_exposure!(boundary)) == UInt64(1)
    @test hil_frame_buffer(boundary) == Float32[9, 9]
    foreach(buffer -> fill!(buffer, 3.0f0), buffers)
    @test @allocated(step_hil_exposure!(boundary)) == 0
    @test hil_frame_buffer(boundary) == Float32[9, 9]
    @test hil_boundary_status(boundary).active_probe_sequence == 7
    @test_throws AlgorithmGraphError adopt_hil_probe!(boundary, UInt64(7))
    @test @allocated(adopt_hil_probe!(boundary, UInt64(8))) == 0
    @test @inferred(reset_hil_boundary!(boundary)) === boundary
    @test boundary.command_input == buffers == initial
    @test all(iszero, hil_frame_buffer(boundary))
    @test_throws AlgorithmGraphError step_hil_exposure!(boundary)

    # A malformed last buffer cannot copy either earlier command or reset the graph.
    resize!(buffers.offload, 2)
    foreach(buffer -> fill!(buffer, 8.0f0), buffers)
    @test_throws AlgorithmGraphError adopt_hil_probe!(boundary, UInt64(1))
    @test boundary.command_input == initial
    @test_throws AlgorithmGraphError reset_hil_boundary!(boundary)
    @test hil_boundary_status(boundary).failed
    @test boundary.command_input == initial
    resize!(buffers.offload, 3)
    reset_hil_boundary!(boundary)

    # Reuse ordinary HIL's exact names/order, shapes and distinct-storage admission.
    for invalid in (
        (dm=buffers.dm, tiptilt=buffers.tiptilt),
        (tiptilt=buffers.tiptilt, dm=buffers.dm, offload=buffers.offload),
        merge(buffers, (foreign=zeros(Float32, 1),)),
        merge(buffers, (offload=zeros(Float32, 2),)),
        merge(buffers, (offload=zeros(Float64, 3),)),
        merge(buffers, (offload=view(buffers.offload, :),)),
        merge(buffers, (dm=buffers.offload,)),
    )
        @test_throws AlgorithmGraphError prepare(multi_command_hil_test_graph(); command_buffer=invalid)
    end
    for names in ((), (:dm, :dm), (:dm, :missing), (:dm, 1), ntuple(_ -> :dm, 17))
        @test_throws AlgorithmGraphError prepare_graph_calibration_boundary(
            multi_command_hil_test_graph(); command_input=names, frame_output=:frame,
        )
    end
    one = prepare_graph_calibration_boundary(
        calibration_test_graph(); command_input=(:command,), frame_output=:frame,
    )
    @test keys(hil_command_buffer(one)) == (:command,)
    @test adopt_hil_probe!(one, UInt64(1)) == UInt64(1)

    # A target-copy failure after an earlier command changed fences every exposure.
    fault = CalibrationTestFailingCommand(boundary.command_input.offload, Ref(true))
    failing = PreparedGraphCalibrationBoundary(
        boundary.graph, merge(boundary.command_input, (offload=fault,)),
        boundary.frame_output, boundary.command_buffer, boundary.frame_buffer,
        boundary.initial_command, boundary.state,
    )
    foreach(buffer -> fill!(buffer, 8.0f0), hil_command_buffer(failing))
    @test_throws ErrorException adopt_hil_probe!(failing, UInt64(2))
    @test failing.command_input.dm == Float32[8, 8]
    @test fault.values == Float32[8, 5, 6]
    @test hil_boundary_status(failing).failed
    @test hil_boundary_status(failing).active_probe_sequence == 0
    @test_throws AlgorithmGraphError step_hil_exposure!(failing)
    @test_throws AlgorithmGraphError adopt_hil_probe!(failing, UInt64(3))
    @test_throws ErrorException reset_hil_boundary!(failing)
    fault.fail[] = false
    reset_hil_boundary!(failing)
    @test failing.command_input == initial
    @test hil_command_buffer(failing) == initial
    @test_throws AlgorithmGraphError step_hil_exposure!(failing)
end

@testset "graph HIL complete multi-command response" begin
    graph = multi_command_hil_test_graph()
    buffers = (dm=zeros(Float32, 2), tiptilt=zeros(Float32, 1), offload=zeros(Float32, 3))
    boundary = multi_command_hil_test_boundary(graph; command_buffer=buffers)
    initial = map(copy, boundary.command_input)
    @test hil_command_buffer(boundary) === buffers
    @test hil_frame_buffer(boundary) isa Array
    @test buffers == initial
    @test all(isconcretetype, fieldtypes(typeof(boundary)))
    @test_throws AlgorithmGraphError adopt_hil_command!(boundary, UInt64(1))
    @test @inferred(step_hil_frame!(boundary)) == UInt64(1)
    @test hil_frame_buffer(boundary) == initial.dm
    @test_throws AlgorithmGraphError step_hil_frame!(boundary)
    foreach(buffer -> fill!(buffer, 9.0f0), buffers)
    @test_throws AlgorithmGraphError adopt_hil_command!(boundary, UInt64(2))
    @test_throws AlgorithmGraphError adopt_hil_command!(boundary, 1)
    @test boundary.command_input == initial
    for invalid in (NaN32, Inf32, -Inf32)
        buffers.offload[end] = invalid
        @test_throws AlgorithmGraphError adopt_hil_command!(boundary, UInt64(1))
        @test boundary.command_input == initial # Neither finite earlier port was copied.
        @test !hil_boundary_status(boundary).failed
        @test hil_boundary_status(boundary).command_response_required
    end
    buffers.offload[end] = 9.0f0
    @test @inferred(adopt_hil_command!(boundary, UInt64(1))) == UInt64(1)
    @test boundary.command_input == buffers
    @test_throws AlgorithmGraphError adopt_hil_command!(boundary, UInt64(1))
    @test step_hil_frame!(boundary) == UInt64(2)
    @test hil_frame_buffer(boundary) == Float32[9, 9]
    foreach(buffer -> fill!(buffer, 7.0f0), buffers)
    @test @allocated(adopt_hil_command!(boundary, UInt64(2))) == 0
    @test @allocated(step_hil_frame!(boundary)) == 0
    @test @allocated(adopt_hil_command!(boundary, UInt64(3))) == 0
    @test @inferred(reset_hil_boundary!(boundary)) === boundary
    @test boundary.command_input == initial
    @test buffers == initial
    @test all(iszero, hil_frame_buffer(boundary))
    @test hil_boundary_status(boundary) == (
        frame_sequence=UInt64(0), active_command_sequence=UInt64(0),
        command_response_required=false, failed=false,
    )

    one = prepare_graph_hil_boundary(
        calibration_test_graph(); command_input=(:command,), frame_output=:frame,
    )
    @test hil_command_buffer(one) isa NamedTuple{(:command,)}
    @test step_hil_frame!(one) == UInt64(1)
    @test adopt_hil_command!(one, UInt64(1)) == UInt64(1)
end

@testset "graph HIL multi-command preparation contracts" begin
    graph = multi_command_hil_test_graph()
    for names in ((), (:dm, :dm), (:dm, :missing), (:dm, 1), ntuple(_ -> :dm, 17))
        @test_throws AlgorithmGraphError prepare_graph_hil_boundary(
            graph; command_input=names, frame_output=:frame,
        )
    end
    valid = (dm=zeros(Float32, 2), tiptilt=zeros(Float32, 1), offload=zeros(Float32, 3))
    for invalid in (
        (dm=valid.dm, tiptilt=valid.tiptilt),
        (tiptilt=valid.tiptilt, dm=valid.dm, offload=valid.offload),
        merge(valid, (extra=zeros(Float32, 1),)),
        merge(valid, (offload=zeros(Float32, 2),)),
        merge(valid, (offload=zeros(Float64, 3),)),
        merge(valid, (offload=view(zeros(Float32, 3), :),)),
        zeros(Float32, 6),
    )
        @test_throws AlgorithmGraphError multi_command_hil_test_boundary(
            graph; command_buffer=invalid,
        )
    end
    nonfinite = multi_command_hil_test_graph(;
        commands=(dm=Float32[1, 2], tiptilt=Float32[3], offload=Float32[4, NaN, 6]),
    )
    foreach(buffer -> fill!(buffer, 42.0f0), valid)
    @test_throws AlgorithmGraphError multi_command_hil_test_boundary(
        nonfinite; command_buffer=valid,
    )
    @test all(buffer -> all(==(42.0f0), buffer), valid)
    ticket = step_graph_async!(graph)
    @test_throws AlgorithmGraphError multi_command_hil_test_boundary(graph)
    wait_graph_step!(ticket)
    @test_throws AlgorithmGraphError multi_command_hil_test_boundary(graph)
end

@testset "graph HIL multi-command storage ownership" begin
    alias_vector(values) = reshape(reshape(values, 1, length(values)), length(values))
    for role in (:graph_commands, :graph_frame, :exchange_commands,
                 :command_other_graph, :command_graph_frame,
                 :frame_other_graph, :frame_exchange)
        commands = (dm=Float32[1, 2], tiptilt=Float32[3, 4], offload=Float32[5, 6])
        role === :graph_commands && (commands = merge(commands, (offload=alias_vector(commands.dm),)))
        graph = multi_command_hil_test_graph(; commands)
        graph_frame = graph_output(graph, Val(:frame))
        role === :graph_frame && (graph_frame = commands.offload)
        # Substitute only the published frame for fault injection; exact node owners remain intact.
        role === :graph_frame && (graph = AOG.PreparedAlgorithmGraph(
            graph.name, graph.target, graph.context, graph.nodes, graph.execution,
            graph.delayed_links, graph.inputs, (frame=graph_frame,), graph.state,
        ))
        buffers = (dm=Float32[8, 9], tiptilt=Float32[10, 11], offload=Float32[12, 13])
        frame_buffer = Float32[14, 15]
        if role === :exchange_commands
            buffers = merge(buffers, (offload=alias_vector(buffers.dm),))
        elseif role === :command_other_graph
            buffers = merge(buffers, (dm=alias_vector(commands.offload),))
        elseif role === :command_graph_frame
            buffers = merge(buffers, (offload=alias_vector(graph_frame),))
        elseif role === :frame_other_graph
            frame_buffer = alias_vector(commands.offload)
        elseif role === :frame_exchange
            frame_buffer = alias_vector(buffers.offload)
        end
        saved_commands, saved_buffers, saved_frame = map(copy, commands), map(copy, buffers), copy(frame_buffer)
        @test_throws AlgorithmGraphError multi_command_hil_test_boundary(
            graph; command_buffer=buffers, frame_buffer,
        )
        @test commands == saved_commands
        @test buffers == saved_buffers
        @test frame_buffer == saved_frame
    end
end

@testset "graph HIL multi-command pending frames and partial copy failure" begin
    pending = multi_command_hil_test_boundary()
    step_hil_frame!(pending)
    ticket = step_graph_async!(pending.graph) # Deliberate violation of sole execution ownership.
    saved = map(copy, pending.command_input)
    foreach(buffer -> fill!(buffer, 8.0f0), hil_command_buffer(pending))
    @test_throws AlgorithmGraphError adopt_hil_command!(pending, UInt64(1))
    @test pending.command_input == saved
    wait_graph_step!(ticket)
    @test_throws AlgorithmGraphError adopt_hil_command!(pending, UInt64(1))
    reset_hil_boundary!(pending)
    @test step_hil_frame!(pending) == UInt64(1)

    prepared = multi_command_hil_test_boundary()
    fault = CalibrationTestFailingCommand(prepared.command_input.offload, Ref(true))
    failing = PreparedGraphHILBoundary(
        prepared.graph, merge(prepared.command_input, (offload=fault,)),
        prepared.frame_output, prepared.command_buffer, prepared.frame_buffer,
        prepared.initial_command, prepared.state,
    )
    step_hil_frame!(failing)
    foreach(buffer -> fill!(buffer, 8.0f0), hil_command_buffer(failing))
    @test_throws ErrorException adopt_hil_command!(failing, UInt64(1))
    @test failing.command_input.dm == Float32[8, 8]
    @test failing.command_input.offload.values == Float32[8, 5, 6]
    @test hil_boundary_status(failing).failed
    @test hil_boundary_status(failing).active_command_sequence == UInt64(0)
    @test_throws AlgorithmGraphError step_hil_frame!(failing)
    @test_throws AlgorithmGraphError adopt_hil_command!(failing, UInt64(1))
    @test_throws ErrorException reset_hil_boundary!(failing)
    @test hil_boundary_status(failing).failed
    fault.fail[] = false
    reset_hil_boundary!(failing)
    @test failing.command_input.dm == Float32[1, 2]
    @test failing.command_input.offload.values == Float32[4, 5, 6]
    @test hil_command_buffer(failing) == failing.initial_command
    @test !hil_boundary_status(failing).failed
end

@testset "HIL single-command and calibration buffers retain complete shapes" begin
    for prepare in (prepare_graph_hil_boundary, prepare_graph_calibration_boundary), delta in (-1, 1)
        graph = calibration_test_graph()
        boundary = prepare(graph; command_input=:command, frame_output=:frame)
        adopt! = prepare === prepare_graph_hil_boundary ? adopt_hil_command! : adopt_hil_probe!
        step! = prepare === prepare_graph_hil_boundary ? step_hil_frame! : step_hil_exposure!
        prepare === prepare_graph_hil_boundary && step!(boundary)
        command = hil_command_buffer(boundary)
        resize!(command, 2 + delta)
        fill!(command, 8.0f0)
        saved_status = hil_boundary_status(boundary)
        @test_throws AlgorithmGraphError adopt!(boundary, UInt64(1))
        @test graph_input(graph, Val(:command)) == Float32[1, 2]
        @test hil_boundary_status(boundary) == saved_status
        resize!(command, 2)
        fill!(command, 8.0f0)
        @test adopt!(boundary, UInt64(1)) == UInt64(1)
        step!(boundary)
        saved_sequence = graph_step_sequence(graph)
        saved_frame = copy(hil_frame_buffer(boundary))
        resize!(command, 2 + delta)
        fill!(command, 9.0f0)
        saved_command = copy(command)
        @test_throws AlgorithmGraphError reset_hil_boundary!(boundary)
        @test hil_boundary_status(boundary).failed
        @test graph_step_sequence(graph) == saved_sequence
        @test graph_input(graph, Val(:command)) == Float32[8, 8]
        @test command == saved_command
        @test hil_frame_buffer(boundary) == saved_frame
        @test_throws AlgorithmGraphError adopt!(boundary, UInt64(2))
        @test_throws AlgorithmGraphError step!(boundary)
        resize!(command, 2)
        @test reset_hil_boundary!(boundary) === boundary
        @test command == Float32[1, 2]
        @test graph_input(graph, Val(:command)) == Float32[1, 2]
        @test !hil_boundary_status(boundary).failed
        @test graph_step_sequence(graph) == UInt64(0)
    end
end

@testset "HIL multi-command shape checks complete before mutation" begin
    for name in (:dm, :offload), delta in (-1, 1)
        boundary = multi_command_hil_test_boundary()
        initial = map(copy, boundary.command_input)
        step_hil_frame!(boundary)
        buffers = hil_command_buffer(boundary)
        command = buffers[name]
        extent = length(command)
        resize!(command, extent + delta)
        foreach(buffer -> fill!(buffer, 8.0f0), buffers)
        saved_status = hil_boundary_status(boundary)
        @test_throws AlgorithmGraphError adopt_hil_command!(boundary, UInt64(1))
        @test boundary.command_input == initial
        @test hil_boundary_status(boundary) == saved_status
        resize!(command, extent)
        fill!(command, 8.0f0)
        @test adopt_hil_command!(boundary, UInt64(1)) == UInt64(1)
        step_hil_frame!(boundary)
        saved_sequence = graph_step_sequence(boundary.graph)
        saved_commands = map(copy, boundary.command_input)
        saved_frame = copy(hil_frame_buffer(boundary))
        resize!(command, extent + delta)
        foreach(buffer -> fill!(buffer, 9.0f0), buffers)
        saved_buffers = map(copy, buffers)
        @test_throws AlgorithmGraphError reset_hil_boundary!(boundary)
        @test hil_boundary_status(boundary).failed
        @test graph_step_sequence(boundary.graph) == saved_sequence
        @test boundary.command_input == saved_commands
        @test buffers == saved_buffers
        @test hil_frame_buffer(boundary) == saved_frame
        @test_throws AlgorithmGraphError adopt_hil_command!(boundary, UInt64(2))
        @test_throws AlgorithmGraphError step_hil_frame!(boundary)
        resize!(command, extent)
        @test reset_hil_boundary!(boundary) === boundary
        @test boundary.command_input == initial
        @test buffers == initial
        @test !hil_boundary_status(boundary).failed
        @test graph_step_sequence(boundary.graph) == UInt64(0)
    end
end

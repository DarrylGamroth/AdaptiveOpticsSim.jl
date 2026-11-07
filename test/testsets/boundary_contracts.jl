using Test, AdaptiveOpticsSim
using AdaptiveOpticsSim.AlgorithmGraphs
using AdaptiveOpticsSim.AlgorithmGraphs: graph_input_contract, graph_output_contract

@testset "cold graph boundary contracts describe adapter ports" begin
    node = second_order_response_node(:response;
        natural_frequency_hz = 20, damping_ratio = 0.7, sample_period_s = 0.001,
        count = 2, schema = "test.response/1")
    command = zeros(Float32, 2)
    definition = algorithm_graph((node,);
        inputs = (graph_input(:command, :response => :input, command),),
        outputs = (graph_output(:frame, :response => :output),))
    input = graph_input_contract(definition, :command)
    output = graph_output_contract(definition, :frame)
    @test input == (; name = :command, direction = :input, role = :data,
        eltype = Float32, shape = (2,), schema = "test.response/1", layout = :column_major)
    @test output == merge(input, (; name = :frame, direction = :output))
    @test all(iszero, command)
    mismatched = algorithm_graph((node,);
        inputs = (graph_input(:command, :response => :input, zeros(Float64, 3)),))
    @test graph_input_contract(mismatched, :command) == input
    @test_throws AlgorithmGraphError graph_input_contract(definition, :missing)
    @test_throws AlgorithmGraphError graph_output_contract(definition, :missing)
    reversed = algorithm_graph((node,);
        inputs = (graph_input(:command, :response => :output, command),),
        outputs = (graph_output(:frame, :response => :input),))
    @test_throws AlgorithmGraphError graph_input_contract(reversed, :command)
    @test_throws AlgorithmGraphError graph_output_contract(reversed, :frame)
    aliased = algorithm_graph((node,);
        outputs = (graph_output(:first, :response => :output),
            graph_output(:second, :response => :output)))
    @test graph_output_contract(aliased, :first).schema == "test.response/1"
    @test graph_output_contract(aliased, :second).shape == (2,)
end

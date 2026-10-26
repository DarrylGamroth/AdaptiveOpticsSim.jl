# Run explicitly with CUDA and this checkout in the active environment:
# JULIA_NUM_THREADS=1 julia --startup-file=no --project=<environment> test/cuda_second_order_response.jl
using Test
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Backends
using AdaptiveOpticsSim.AlgorithmGraphs
import AdaptiveOpticsSim.AlgorithmGraphs: prepared_graph_node
using CUDA

CUDA.allowscalar(false)

function cuda_response_definition(input)
    node = second_order_response_node(:response; natural_frequency_hz=40,
        damping_ratio=0.7, sample_period_s=0.001, count=length(input),
        schema="test.response-coordinate/1", T=eltype(input))
    return algorithm_graph((node,); name=:cuda_response,
        inputs=(graph_input(:command, :response => :input, input),),
        outputs=(graph_output(:position, :response => :output),))
end

@testset "CUDA second-order response capture and reset" begin
    for T in (Float32, Float64)
        host_input = zeros(T, 129)
        cpu = prepare_algorithm_graph(cuda_response_definition(host_input))
        stream_input, captured_input = CuArray(host_input), CuArray(host_input)
        target = compute_device(stream_input)
        stream = prepare_algorithm_graph(cuda_response_definition(stream_input); target)
        captured = prepare_algorithm_graph(cuda_response_definition(captured_input);
            target, execution=CapturedGraphExecution())
        @test captured_graph_node_count(captured) == 1
        owner = prepared_graph_node(captured, Val(:response)).prepared
        @test AdaptiveOpticsSim.AlgorithmGraphs.graph_node_capture_capability(
            prepared_graph_node(captured, Val(:response))) isa
            AdaptiveOpticsSim.AlgorithmGraphs.GraphNodeCaptureSafe
        addresses = pointer.((owner.state.position, owner.state.velocity, owner.output))
        first_output = zeros(T, 129)
        tolerance = T === Float32 ? 5e-6 : 5e-13
        for repetition in 1:2
            reset_graph!(cpu)
            reset_graph!(stream)
            reset_graph!(captured)
            @test pointer.((owner.state.position, owner.state.velocity, owner.output)) == addresses
            @test all(iszero, Array(owner.state.position))
            @test all(iszero, Array(owner.state.velocity))
            @test all(iszero, Array(owner.output))
            for n in 1:100
                host_input .= T(sin(0.13n)) .* T.(1:129) ./ T(129)
                copyto!(stream_input, host_input)
                copyto!(captured_input, host_input)
                step_graph!(cpu)
                step_graph!(stream)
                step_graph!(captured)
                reference = graph_output(cpu, Val(:position))
                @test Array(graph_output(stream, Val(:position))) ≈ reference rtol=tolerance atol=tolerance
                @test Array(graph_output(captured, Val(:position))) ≈ reference rtol=tolerance atol=tolerance
                if n == 1
                    repetition == 1 ? copyto!(first_output, Array(owner.output)) :
                        (@test Array(owner.output) == first_output)
                end
            end
        end
    end
end

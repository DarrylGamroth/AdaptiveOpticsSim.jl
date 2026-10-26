function response_step_reference(t, frequency, damping)
    omega = 2π * frequency
    if damping < 1
        q = sqrt(1 - damping^2)
        return 1 - exp(-damping*omega*t) *
            (cos(omega*q*t) + damping/q*sin(omega*q*t))
    elseif damping == 1
        return 1 - exp(-omega*t) * (1 + omega*t)
    end
    r1 = -omega * (damping - sqrt(damping^2 - 1))
    r2 = -omega * (damping + sqrt(damping^2 - 1))
    return 1 + (r2*exp(r1*t) - r1*exp(r2*t)) / (r1-r2)
end

function response_test_definition(input; damping_ratio=0.7)
    node = second_order_response_node(:response; natural_frequency_hz=40,
        damping_ratio, sample_period_s=0.001, count=length(input),
        schema="test.response-coordinate/1", T=eltype(input))
    return algorithm_graph((node,); name=:response_test,
        inputs=(graph_input(:command, :response => :input, input),),
        outputs=(graph_output(:position, :response => :output),))
end

function response_allocation_bytes(owner)
    advance_response!(owner)
    return @allocated advance_response!(owner)
end

@testset "exact held-input second-order response" begin
    optics = AdaptiveOpticsSim.Optics
    for T in (Float32, Float64), damping in (0.2, 1.0, 2.0)
        input = T[1, -2, 0]
        output = similar(input)
        plan = optics.SecondOrderResponsePlan(; natural_frequency_hz=40,
            damping_ratio=damping, sample_period_s=0.001, T)
        owner = prepare_second_order_response(input, plan; output)
        @test owner.input === input && owner.output === output
        @test owner.state.position !== output && owner.state.velocity !== output
        @test owner.state.position !== owner.state.velocity
        @test owner.state.position isa FixedSizeVector{T}
        @test owner.state.velocity isa FixedSizeVector{T}
        @test all(iszero, output)
        tolerance = T === Float32 ? 5e-6 : 5e-13
        first_output = similar(output)
        for n in 1:1000
            advance_response!(owner)
            if n in (1, 2, 10, 50, 100, 500, 1000)
                expected = response_step_reference(n*0.001, 40, damping)
                @test output ≈ input .* expected rtol=tolerance atol=tolerance
            end
            n == 1 && copyto!(first_output, output)
        end
        @test output ≈ input rtol=tolerance atol=tolerance
        @test reset_response!(owner) === owner
        @test owner.input === input && owner.output === output
        @test all(iszero, owner.state.position) && all(iszero, owner.state.velocity)
        @test all(iszero, output)
        advance_response!(owner)
        @test output == first_output
        if !coverage_instrumented()
            @test response_allocation_bytes(owner) == 0
        end
    end

    # Piecewise-constant sinusoidal forcing: independently derive the complex
    # sampled frequency response from differences of the continuous step curve.
    for signal_frequency in (5.0, 40.0, 150.0)
        input = [0.0]
        plan = optics.SecondOrderResponsePlan(; natural_frequency_hz=40,
            damping_ratio=0.7, sample_period_s=0.001, T=Float64)
        owner = prepare_second_order_response(input, plan)
        phase_step = 2π * signal_frequency * 0.001
        transfer = sum((response_step_reference(k*0.001, 40, 0.7) -
            response_step_reference((k-1)*0.001, 40, 0.7)) * cis(-phase_step*(k-1))
            for k in 1:1000)
        for n in 1:600
            input[1] = sin(phase_step*(n-1))
            advance_response!(owner)
            if n > 400
                @test owner.output[1] ≈ imag(transfer*cis(phase_step*(n-1))) atol=1e-12
            end
        end
    end

    # One-based strided host views follow the same independent-coordinate law.
    storage = [1.0, 99.0, -2.0, 99.0]
    plan = optics.SecondOrderResponsePlan(; natural_frequency_hz=40,
        damping_ratio=1, sample_period_s=0.001, T=Float64)
    owner = prepare_second_order_response(view(storage, 1:2:4), plan)
    @test advance_response!(owner) ≈ [1, -2] .* response_step_reference(0.001, 40, 1)
    @test storage == [1, 99, -2, 99]
end

@testset "second-order response validation and graph adapter" begin
    optics = AdaptiveOpticsSim.Optics
    defaults = (; natural_frequency_hz=40, damping_ratio=0.7, sample_period_s=0.001)
    for key in keys(defaults), value in (0, -1, Inf, NaN, true)
        @test_throws AdaptiveOpticsSim.InvalidConfiguration optics.SecondOrderResponsePlan(;
            merge(defaults, NamedTuple{(key,)}((value,)))...)
    end
    @test_throws AdaptiveOpticsSim.InvalidConfiguration optics.SecondOrderResponsePlan(; defaults..., T=Float16)
    @test_throws AdaptiveOpticsSim.InvalidConfiguration optics.SecondOrderResponsePlan(;
        merge(defaults, (; sample_period_s=1e-30))...)
    @test_throws AdaptiveOpticsSim.InvalidConfiguration optics.SecondOrderResponsePlan(;
        merge(defaults, (; natural_frequency_hz=1e300))...)
    # Very long intervals can settle entirely within one hold without instability.
    settled = optics.SecondOrderResponsePlan(; merge(defaults, (; sample_period_s=100.0))...)
    @test advance_response!(prepare_second_order_response(Float32[1], settled)) ≈ Float32[1]
    plan = optics.SecondOrderResponsePlan(; defaults...)
    @test_throws AdaptiveOpticsSim.DimensionMismatchError prepare_second_order_response(Float32[], plan)
    @test_throws AdaptiveOpticsSim.DimensionMismatchError prepare_second_order_response(Float32[1], plan; output=zeros(Float32, 2))
    input = Float32[1, 2]
    @test_throws AdaptiveOpticsSim.InvalidConfiguration prepare_second_order_response(input, plan; output=input)
    overlap = Float32[1, 2, 3]
    @test_throws AdaptiveOpticsSim.InvalidConfiguration prepare_second_order_response(
        view(overlap, 1:2), plan; output=view(overlap, 2:3))
    for count in (0, -1, true)
        @test_throws AlgorithmGraphError second_order_response_node(:response;
            defaults..., count, schema="test.response/1")
    end
    @test_throws AlgorithmGraphError second_order_response_node(:response;
        defaults..., count=2, schema="")
    graph = prepare_algorithm_graph(response_test_definition(input))
    owner = prepared_graph_node(graph, Val(:response))
    @test AdaptiveOpticsSim.AlgorithmGraphs.graph_node_capture_capability(owner) isa
        AdaptiveOpticsSim.AlgorithmGraphs.GraphNodeCaptureUnsupported
    step_graph!(graph)
    first = copy(graph_output(graph, Val(:position)))
    @test first ≈ input .* response_step_reference(0.001, 40, 0.7)
    step_graph!(graph)
    reset_graph!(graph)
    @test all(iszero, graph_output(graph, Val(:position)))
    step_graph!(graph)
    @test graph_output(graph, Val(:position)) == first

    document = Dict("schema_version" => 1, "name" => "response_file",
        "nodes" => [Dict("name" => "response", "type" => "second_order_response_f32",
            "config" => Dict("natural_frequency_hz" => 40, "damping_ratio" => 0.7,
                "sample_period_s" => 0.001, "count" => 2, "schema" => "test.response-coordinate/1"))],
        "inputs" => [Dict("name" => "command", "destination" => "response.input", "binding" => "command")],
        "outputs" => [Dict("name" => "position", "source" => "response.output")])
    loaded = prepare_algorithm_graph(load_algorithm_graph(document; bindings=(command=input,)))
    step_graph!(loaded)
    @test graph_output(loaded, Val(:position)) == first
    document["nodes"][1]["config"]["unknown"] = 1
    @test_throws AlgorithmGraphError load_algorithm_graph(document; bindings=(command=input,))
end

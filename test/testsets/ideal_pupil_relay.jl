using Test, AdaptiveOpticsSim
using AdaptiveOpticsSim.AlgorithmGraphs

function prepare_test_pupil_relay(input, output_size; permutation=(2,1), overrides=(;))
    T = real(eltype(input))
    options = (input_size=size(input), output_size, axes=permutation,
        input_sampling_m=(T(0.002),T(0.003)),
        output_sampling_m=ntuple(k -> T(2)*(T(0.002),T(0.003))[permutation[k]],2),
        wavelength_m=T(750e-9), input_schema="test.dimensionless-input-field/1",
        output_schema="test.dimensionless-pupil-field/1", T)
    node = ideal_pupil_relay_node(:relay; merge(options,overrides)...)
    return prepare_algorithm_graph(algorithm_graph((node,);
        inputs=(graph_input(:field,:relay=>:input,input),),
        outputs=(graph_output(:field,:relay=>:output),)))
end

@testset "Matched-grid ideal pupil relay" begin
    for T in (Float32,Float64), permutation in ((1,2),(2,1)),
            (source_size,output_size) in (((8,10),(5,6)),((9,11),(6,7)),((8,8),(8,8)))
        input = Complex{T}[complex(T(10i+j),T(i-2j)) for i in 1:source_size[1],j in 1:source_size[2]]
        original = copy(input)
        graph = prepare_test_pupil_relay(input,output_size;permutation)
        output = graph_output(graph,:field)
        owner = AlgorithmGraphs.prepared_graph_node(graph,:relay)
        @test owner.input === input
        @test owner.output === output
        @test @inferred(step_graph!(graph)) === graph
        expected = similar(output)
        source_center = ntuple(k -> fld(source_size[k],2)+1,2)
        output_center = ntuple(k -> fld(output_size[k],2)+1,2)
        for j in axes(expected,2), i in axes(expected,1)
            a = i-output_center[1]+source_center[permutation[1]]
            b = j-output_center[2]+source_center[permutation[2]]
            expected[i,j] = permutation == (1,2) ? original[a,b] : original[b,a]
        end
        @test output == expected
        @test output[output_center...] == original[source_center...]
        @test input == original
        @test sum(abs2,output) <= sum(abs2,input)
        if output_size != source_size
            @test sum(abs2,output) < sum(abs2,input)
        else
            @test sum(abs2,output) ≈ sum(abs2,input)
        end
        input .*= T(0.5)
        step_graph!(graph)
        @test output == T(0.5).*expected
        reset_graph!(graph)
        @test graph_output(graph,:field) === output
        @test all(iszero,output)
        step_graph!(graph)
        @test output == T(0.5).*expected
        @test @allocated(step_graph!(graph)) == 0
    end

    input = zeros(ComplexF32,8,10)
    for bad in ((output_size=(11,8),), (axes=(1,1),), (axes=(true,2),),
            (input_sampling_m=(0.0,0.003),), (output_sampling_m=(0.006,0.008),),
            (wavelength_m=0.0,), (input_schema="",), (output_schema="",),
            (input_size=(true,10),))
        @test_throws AlgorithmGraphError prepare_test_pupil_relay(input,(5,6);overrides=bad)
    end
end

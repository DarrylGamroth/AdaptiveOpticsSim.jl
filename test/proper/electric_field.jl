using AdaptiveOpticsSim, Proper, Test
using AdaptiveOpticsSim.AlgorithmGraphs

function electric_field_train_graph(config; T=Float64,
        opd=zeros(T,config.pupil_resolution,config.pupil_resolution),
        amplitude=ones(T,size(opd)), runtime_inputs=(;))
    node = proper_optical_train_node_factory(;T)(:optics, config, (;))
    inputs = (graph_input(:opd, :optics=>:pupil_opd, opd),
        graph_input(:amplitude, :optics=>:pupil_amplitude, amplitude),
        map(keys(runtime_inputs)) do name
            graph_input(name, :optics=>name, getproperty(runtime_inputs,name))
        end...)
    return prepare_algorithm_graph(algorithm_graph((node,); inputs,
        outputs=(graph_output(:field, :optics=>:output),)))
end

@testset "Proper planar complex field output" begin
    for T in (Float32,Float64)
        config = (resolution=16, pupil_resolution=8, diameter_m=0.008,
            wavelength_um=0.75, steps=(), output_quantity=:electric_field)
        amplitude = reshape(T.(1:64),8,8)/T(64)
        opd = reshape(T.(1:64),8,8)*T(1e-9)
        graph = electric_field_train_graph(config;T,opd,amplitude)
        output = graph_output(graph,:field)
        owner = AlgorithmGraphs.prepared_graph_node(graph,:optics)
        @test eltype(output) === Complex{T}
        @test optical_sampling_m(owner) ≈ T(0.001)
        node = proper_optical_train_node_factory(;T)(:cold,config,(;))
        @test optical_sampling_m(node.config) == optical_sampling_m(owner)
        @test owner.output === output
        @test @inferred(step_graph!(graph)) === graph
        expected = zeros(Complex{T},16,16)
        expected[5:12,5:12] .= amplitude .* cis.(T(2pi)/T(0.75e-6).*opd)
        @test output ≈ expected rtol=20eps(T)
        intensity = electric_field_train_graph(merge(config,(output_quantity="intensity",));T,opd,amplitude)
        step_graph!(intensity)
        @test abs2.(output) ≈ graph_output(intensity,:field) rtol=20eps(T)
        saved = copy(output)
        reset_graph!(graph)
        @test graph_output(graph,:field) === output
        @test all(iszero,output)
        step_graph!(graph)
        @test output == saved
        amplitude .*= T(0.5)
        step_graph!(graph)
        @test output ≈ T(0.5).*saved rtol=20eps(T)

        screen = (op="add_opd", name="separate pupil screen",input="correction",
            size=(8,8),axes=(1,2),sampling_m=(0.001,0.001),center_pixel=(5,5),
            schema="test.runtime-opd/1",outside_opd_m=0)
        runtime = electric_field_train_graph(merge(config,(steps=(screen,),));T,opd,amplitude,
            runtime_inputs=(correction=opd,))
        step_graph!(runtime)
        expected[5:12,5:12] .= amplitude .* cis.(T(4pi)/T(0.75e-6).*opd)
        @test graph_output(runtime,:field) ≈ expected rtol=30eps(T)

        propagated_config = merge(config,(steps=((op="lens",focal_length_m=0.2),
            (op="propagate",distance_m=0.2)),))
        propagated = electric_field_train_graph(propagated_config;T,opd,amplitude)
        propagated_intensity = electric_field_train_graph(
            merge(propagated_config,(output_quantity=:intensity,));T,opd,amplitude)
        pitch_before = optical_sampling_m(AlgorithmGraphs.prepared_graph_node(propagated,:optics))
        @test pitch_before > zero(T)
        step_graph!(propagated)
        step_graph!(propagated_intensity)
        @test abs2.(graph_output(propagated,:field)) ≈
            graph_output(propagated_intensity,:field) rtol=30eps(T)
        @test optical_sampling_m(AlgorithmGraphs.prepared_graph_node(propagated,:optics)) == pitch_before
        @test pitch_before != T(0.001)

        @test_throws AlgorithmGraphError electric_field_train_graph(
            merge(config,(steps=((op="lens",focal_length_m=0.2),),));T)
        spherical = proper_optical_train_node_factory(;T)(:spherical,
            merge(config,(steps=((op="lens",focal_length_m=0.2),),)),(;))
        @test_throws AlgorithmGraphError optical_sampling_m(spherical.config)
        @test_throws AlgorithmGraphError electric_field_train_graph(
            merge(config,(channels=((wavelength_um=0.75,weight=1.0,steps=()),),));T)
        @test_throws AlgorithmGraphError electric_field_train_graph(
            merge(config,(output_sampling_m=0.001,));T)
        @test_throws AlgorithmGraphError electric_field_train_graph(
            merge(config,(output_quantity=:phase,));T)
    end
end

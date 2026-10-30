using Test, AdaptiveOpticsSim
using AdaptiveOpticsSim.AlgorithmGraphs
using YAML

function pixel_integration_test_node(input_shape, output_shape; overrides=(;), T=Float64)
    geometry = (input_shape, output_shape,
        input_sampling_m=(1.0,1.0), output_sampling_m=(1.0,1.0),
        input_center_pixel=ntuple(k -> (input_shape[k]+1)/2,2),
        output_center_pixel=ntuple(k -> (output_shape[k]+1)/2,2),
        input_schema="test.relative-power.optical-grid/1",
        output_schema="test.relative-power.detector-grid/1", T)
    return pixel_integration_node(:pixels; merge(geometry,overrides)...)
end

function pixel_integration_test_graph(input, output_shape; overrides=(;))
    node = pixel_integration_test_node(size(input), output_shape;
        T=eltype(input), overrides)
    graph = prepare_algorithm_graph(algorithm_graph((node,);
        inputs=(graph_input(:flux,:pixels=>:input,input),),
        outputs=(graph_output(:flux,:pixels=>:output),)))
    return graph
end

pixel_integration_step_allocations(graph) = @allocated step_graph!(graph)

@testset "Rectangular pixel integration" begin
    for T in (Float32,Float64)
        # Exact two-by-three binning preserves each source cell's integrated
        # flux and does not swap array axes on a nonsquare grid.
        input = T[10y+x for y in 1:4,x in 1:6]
        source = copy(input)
        graph = pixel_integration_test_graph(input,(2,2);
            overrides=(output_sampling_m=(2.0,3.0),))
        expected = T[sum(input[2y-1:2y,3x-2:3x]) for y in 1:2,x in 1:2]
        output = graph_output(graph,:flux)
        @test @inferred(step_graph!(graph)) === graph
        @test output == expected
        @test sum(output) == sum(input)
        @test input == source
        owner = AlgorithmGraphs.prepared_graph_node(graph,:pixels)
        @test owner.input === input
        @test owner.output === output
        @test_throws AlgorithmGraphError optical_sampling_m(owner)
        @test AlgorithmGraphs.graph_node_capture_capability(owner) isa
            AlgorithmGraphs.GraphNodeCaptureUnsupported
        input .*= T(0.5)
        step_graph!(graph)
        @test output == T(0.5).*expected
        reset_graph!(graph)
        @test graph_output(graph,:flux) === output
        @test all(iszero,output)
        step_graph!(graph)
        @test output == T(0.5).*expected
        pixel_integration_step_allocations(graph)
        @test pixel_integration_step_allocations(graph) == 0

        # Uniform source density: each source cell has flux 6 and area 6 m².
        # A complete output pixel of area 0.75 m² must contain flux 0.75.
        graph = pixel_integration_test_graph(fill(T(6),4,6),(3,5); overrides=(
            input_sampling_m=(2.0,3.0), output_sampling_m=(0.5,1.5)))
        step_graph!(graph)
        @test graph_output(graph,:flux) ≈ fill(T(0.75),3,5)

        # One source pixel split into four equal area contributions. Changing
        # only its x registration preserves row identity.
        graph = pixel_integration_test_graph(fill(T(8),1,1),(2,2);
            overrides=(output_sampling_m=(0.5,0.5),))
        step_graph!(graph)
        @test graph_output(graph,:flux) == fill(T(2),2,2)
        @test optical_sampling_m(AlgorithmGraphs.prepared_graph_node(graph,:pixels)) == 0.5
        impulse = zeros(T,3,5)
        impulse[1,4] = T(8)
        graph = pixel_integration_test_graph(impulse,(3,5);
            overrides=(output_center_pixel=(2.0,3.25),))
        step_graph!(graph)
        expected = zeros(T,3,5)
        expected[1,4], expected[1,5] = T(6), T(2)
        @test graph_output(graph,:flux) == expected
        @test sum(graph_output(graph,:flux)) == sum(impulse)

        # A fractional two-axis displacement retains exactly 3/4 × 1/2 of
        # this source pixel. The missing flux remains outside the output.
        graph = pixel_integration_test_graph(fill(T(8),1,1),(1,1);
            overrides=(output_center_pixel=(0.75,0.5),))
        step_graph!(graph)
        @test graph_output(graph,:flux) == fill(T(3),1,1)
        @test sum(graph_output(graph,:flux)) < T(8)
        graph = pixel_integration_test_graph(fill(T(8),1,1),(1,1);
            overrides=(output_center_pixel=(-5.0,0.5),))
        step_graph!(graph)
        @test iszero(only(graph_output(graph,:flux)))

        # A larger containing detector grid preserves total flux, including
        # nonintegral pixel boundaries and noncentral input registration.
        graph = pixel_integration_test_graph(source,(8,9); overrides=(
            input_sampling_m=(0.7,0.9), output_sampling_m=(0.6,0.8),
            input_center_pixel=(2.1,3.8)))
        step_graph!(graph)
        @test sum(graph_output(graph,:flux)) ≈ sum(source) rtol=8eps(T)
    end
end

@testset "Pixel integration explicit geometry and boundary validation" begin
    declare(; kwargs...) = pixel_integration_test_node((3,5),(2,4); overrides=(;kwargs...))
    for bad in ((input_shape=(true,5),), (output_shape=(0,4),),
            (input_shape=(typemax(Int),2),), (output_shape=(3.0,4),),
            (input_sampling_m=(0.0,1.0),), (output_sampling_m=(-1.0,1.0),),
            (input_sampling_m=(NaN,1.0),), (output_sampling_m=(Inf,1.0),),
            (input_sampling_m=(true,1.0),), (output_sampling_m=(1.0,),),
            (input_center_pixel=(Inf,1.0),), (output_center_pixel=(NaN,1.0),),
            (output_center_pixel=(1e300,1.0),),
            (output_sampling_m=(1e-300,1.0),),
            (input_schema="",), (output_schema="",))
        @test_throws AlgorithmGraphError declare(; bad...)
    end
    node = declare(output_sampling_m=(2.0,2.0))
    @test optical_sampling_m(node.config) == 2.0
    ports = AlgorithmGraphs.graph_node_ports(AlgorithmGraphs.PixelIntegrationNode{Float64},node.config)
    @test ports[1].format.shape == (3,5)
    @test ports[2].format.shape == (2,4)
    @test ports[1].format.schema == "test.relative-power.optical-grid/1"
    @test ports[2].format.schema == "test.relative-power.detector-grid/1"
    @test_throws AlgorithmGraphError AlgorithmGraphs.prepare_graph_node(
        AlgorithmGraphs.PixelIntegrationNode{Float64},node.config,(;),
        (input=zeros(3,5),),(output=zeros(2,4),),(;),
        AdaptiveOpticsSim.Backends.AcceleratorComputeDevice(AdaptiveOpticsSim.Backends.CUDABackend(),0))
    @test_throws AlgorithmGraphError prepare_algorithm_graph(algorithm_graph((node,);
        inputs=(graph_input(:flux,:pixels=>:input,zeros(Float32,3,5)),)))
    downstream = pixel_integration_test_node((2,4),(1,2);
        overrides=(input_schema="test.incompatible-quantity/1",))
    downstream = algorithm_node(:downstream,AlgorithmGraphs.PixelIntegrationNode{Float64},
        downstream.config)
    @test_throws AlgorithmGraphError prepare_algorithm_graph(algorithm_graph((node,downstream);
        inputs=(graph_input(:flux,:pixels=>:input,ones(3,5)),),
        links=(link(:pixels=>:output,:downstream=>:input),)))

    identity_node = pixel_integration_test_node((2,4),(2,4))
    aliased = ones(2,4)
    @test_throws AlgorithmGraphError AlgorithmGraphs.prepare_graph_node(
        AlgorithmGraphs.PixelIntegrationNode{Float64},identity_node.config,(;),
        (input=aliased,),(output=aliased,),(;),AdaptiveOpticsSim.Backends.HostComputeDevice())

    # Invalid flux anywhere in the source is rejected before output mutation,
    # even when that source cell lies outside the output crop.
    for bad in (-1.0,NaN,Inf)
        input = ones(3,5)
        graph = pixel_integration_test_graph(input,(1,1))
        step_graph!(graph)
        owner = AlgorithmGraphs.prepared_graph_node(graph,:pixels)
        before = copy(owner.output)
        input[1,1] = bad
        @test_throws AlgorithmGraphError AlgorithmGraphs.step_graph_node!(owner)
        @test isequal(owner.output,before)
    end
end

@testset "Pixel integration YAML factories" begin
    for (T,type) in ((Float32,"pixel_integration_f32"),(Float64,"pixel_integration_f64"))
        document = """
        schema_version: 1
        name: pixel_integration
        nodes:
          - name: pixels
            type: $type
            config:
              input_shape: [4, 6]
              output_shape: [2, 3]
              input_sampling_m: [0.00001, 0.00001]
              output_sampling_m: [0.00002, 0.00002]
              input_center_pixel: [2.5, 3.5]
              output_center_pixel: [1.5, 2.0]
              input_schema: test.relative-power.optical-grid/1
              output_schema: test.relative-power.detector-grid/1
        inputs:
          - name: flux
            destination: pixels.input
            binding: flux
        outputs:
          - name: flux
            source: pixels.output
        """
        root = YAML.load(document)
        graph = prepare_algorithm_graph(load_algorithm_graph(root; bindings=(flux=ones(T,4,6),)))
        step_graph!(graph)
        @test graph_output(graph,:flux) == fill(T(4),2,3)
        config = AlgorithmGraphs.prepared_graph_node(graph,:pixels).plan.config
        @test optical_sampling_m(config) == 2e-5
        factory = getproperty(builtin_graph_node_types(),Symbol(type))
        settings = (input_shape=(4,6),output_shape=(2,3),input_sampling_m=(1.,1.),
            output_sampling_m=(2.,2.),input_center_pixel=(2.5,3.5),
            output_center_pixel=(1.5,2.),input_schema="a",output_schema="b")
        @test_throws AlgorithmGraphError factory(:pixels,merge(settings,(unexpected=true,)),(;))
        @test_throws AlgorithmGraphError factory(:pixels,
            Base.structdiff(settings,NamedTuple{(:input_center_pixel,)}),(;))
        @test_throws AlgorithmGraphError factory(:pixels,settings,(gain=1.,))
    end
end

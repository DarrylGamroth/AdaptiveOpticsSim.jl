using AdaptiveOpticsSim, Proper, Test, YAML
using AdaptiveOpticsSim.AlgorithmGraphs

function train_graph(config; opd=zeros(config.pupil_resolution,config.pupil_resolution), amplitude=ones(size(opd)), T=Float64, bindings=(;))
    node = proper_optical_train_node_factory(; T)(:optics,config,(;),bindings)
    graph = algorithm_graph((node,);name=:train,
        inputs=(graph_input(:opd,:optics=>:pupil_opd,opd),graph_input(:amplitude,:optics=>:pupil_amplitude,amplitude)),
        outputs=(graph_output(:intensity,:optics=>:output),))
    prepare_algorithm_graph(graph;target=AdaptiveOpticsSim.Backends.HostComputeDevice())
end

@testset "YAML scalar optical train" begin
    config = (resolution=16,pupil_resolution=8,diameter_m=0.008,wavelength_um=1.55,steps=())
    amplitude = reshape(collect(1.0:64.0),8,8)./64
    opd = zeros(8,8)
    graph = train_graph(config;opd,amplitude)
    step_graph!(graph)
    expected = zeros(16,16); expected[5:12,5:12] .= abs2.(amplitude)
    @test graph_output(graph,:intensity) ≈ expected
    @test sum(graph_output(graph,:intensity)) ≈ sum(abs2,amplitude)
    reset_graph!(graph); step_graph!(graph)
    @test graph_output(graph,:intensity) ≈ expected
    amplitude .*= 0.5
    step_graph!(graph)
    @test graph_output(graph,:intensity) ≈ expected.*0.25

    odd_pupil_config = merge(config, (resolution=6, pupil_resolution=3))
    odd_pupil = reshape(collect(1.0:9.0), 3, 3)
    odd_pupil_graph = train_graph(odd_pupil_config; opd=zeros(3,3), amplitude=odd_pupil)
    step_graph!(odd_pupil_graph)
    odd_pupil_expected = zeros(6,6)
    odd_pupil_expected[3:5,3:5] .= abs2.(odd_pupil)
    @test graph_output(odd_pupil_graph,:intensity) ≈ odd_pupil_expected

    node_types = builtin_graph_node_types()
    yaml_graph = prepare_algorithm_graph(load_algorithm_graph(
        joinpath(@__DIR__,"graph_files","proper_optical_train.yaml");node_types,
        bindings=(pupil_opd=opd,pupil_amplitude=amplitude)))
    step_graph!(yaml_graph)
    @test graph_output(yaml_graph,:intensity) ≈ expected.*0.25
    @test haskey(builtin_graph_node_types(), :proper_optical_train_f32)

    steps = ((op="circular_aperture",radius_m=0.003,name="entrance clear aperture"),
        (op="circular_obscuration",radius_m=0.0005),
        (op="multiply",amplitude=0.7),
        (op="propagate",distance_m=-0.01),
        (op="lens",focal_length_m=0.2),(op="propagate",distance_m=0.2))
    propagation = merge(config,(;steps))
    graph = train_graph(propagation;opd,amplitude)
    step_graph!(graph)
    context = Proper.RunContext(Matrix{ComplexF64})
    wf = prop_begin(0.008,1.55e-6,16;beam_diam_fraction=0.5,context)
    padded = zeros(16,16); padded[5:12,5:12] .= amplitude
    prop_multiply(wf,padded)
    prop_circular_aperture(wf,0.003)
    prop_circular_obscuration(wf,0.0005)
    prop_multiply(wf,0.7)
    prop_propagate(wf,-0.01,context,"propagate")
    prop_lens(wf,0.2,context,"lens")
    prop_propagate(wf,0.2,context,"propagate")
    direct,_ = prop_end(wf)
    owner = AlgorithmGraphs.prepared_graph_node(graph,:optics)
    @test optical_sampling_m(owner) ≈ prop_get_sampling(wf)
    @test owner.plan.steps[1].name == "entrance clear aperture"
    @test @inferred(step_graph!(graph)) === graph
    @test graph_output(graph,:intensity) ≈ direct
    before = copy(graph_output(graph,:intensity))
    opd .= [1e-7*(j-4) for i in 1:8,j in 1:8]
    step_graph!(graph)
    @test !isapprox(graph_output(graph,:intensity),before;rtol=1e-3)

    circles = ((radius_m=0.0011,x_m=0.002,y_m=-0.001),)
    transmission = train_graph(merge(config,(steps=((op="circular_mask",circles=circles,mode="transmit"),),)))
    reflection = train_graph(merge(config,(steps=((op="circular_mask",circles=circles,mode="reflect"),),)))
    step_graph!(transmission); step_graph!(reflection)
    transmitted = graph_output(transmission,:intensity)
    reflected = graph_output(reflection,:intensity)
    @test transmitted+reflected ≈ [i in 5:12 && j in 5:12 ? 1.0 : 0.0 for i in 1:16,j in 1:16]
    @test transmitted[8,11] ≈ 1 # y=-1 mm, x=+2 mm; row then column
    @test transmitted[11,8] ≈ 0

    f32 = train_graph(config; T=Float32, opd=zeros(Float32,8,8), amplitude=ones(Float32,8,8))
    step_graph!(f32)
    @test eltype(graph_output(f32,:intensity)) === Float32
    @test sum(graph_output(f32,:intensity)) ≈ 64f0
    schema = "test.optical-intensity/1;sha256=synthetic"
    identified = train_graph(merge(config,(output_schema=schema,)))
    @test AlgorithmGraphs.prepared_graph_node(identified,:optics).plan.schemas[3] == schema

    for circle in ((radius_m=1e-12,x_m=0.0005,y_m=0.0005),
            (radius_m=0.001,x_m=1.0,y_m=0.0))
        unsupported = merge(config,(steps=((op="circular_mask",mode="transmit",circles=(circle,)),),))
        @test_throws AlgorithmGraphError train_graph(unsupported)
    end

    factory = proper_optical_train_node_factory()
    for bad in ((steps=((op="eval",),),), (diameter_m=true,),
        (pupil_resolution=17,), (resolution=true,), (resolution=17,), (wavelength_um=Inf,),
        (unexpected=1,), (steps=((op="multiply",amplitude=1.1),),),
        (steps=((op="propagate",distance_m=NaN),),),
        (steps=((op="lens",focal_length_m=0.0),),),
        (steps=((op="circular_mask",circles=circles,mode="other"),),),
        (steps=((op="circular_aperture",radius_m=true),),))
        @test_throws AlgorithmGraphError factory(:bad,merge(config,bad),(;))
    end
end

include("proper/sampled_masks.jl")

include("proper/polychromatic_train.jl")

# Independent tensor-product tent-basis reference, including the outside field.
function sampled_mask_reference(amplitude, phase, n, pitch, sampling, center, outside)
    source = amplitude .* cis.(phase)
    return [begin
        y = (row-fld(n,2)-1)*pitch/sampling[1]+center[1]
        x = (col-fld(n,2)-1)*pitch/sampling[2]+center[2]
        outside + sum(max(0, 1-abs(y-i))*max(0, 1-abs(x-j))*(source[i,j]-outside)
            for i in axes(source,1), j in axes(source,2))
    end for row in 1:n, col in 1:n]
end

struct OffsetMaskForTest <: AbstractMatrix{Float64} end
Base.size(::OffsetMaskForTest) = (2,2)
Base.axes(::OffsetMaskForTest) = (0:1,0:1)

@testset "Sampled complex optical masks" begin
    factory = proper_optical_train_node_factory()
    base = (resolution=8,pupil_resolution=8,diameter_m=0.008,wavelength_um=1.55)
    step = (op="sampled_mask",name="sampled apodizer",amplitude="transmission",
        sampling_m=(0.001,0.001),center_pixel=(5.0,5.0),
        outside_amplitude=0.0,outside_phase_rad=0.0)
    config = merge(base,(steps=(step,),transmission="calibration"))
    @test AlgorithmGraphs.graph_node_bindings(factory, config) == (:calibration,)
    @test AlgorithmGraphs.graph_node_configuration(factory, merge(base,(steps=(step,),)),
        (transmission="calibration",)) == config
    @test_throws AlgorithmGraphError factory(:optics,config,(;))
    @test_throws AlgorithmGraphError AlgorithmGraphs.graph_node_bindings(factory,merge(base,(steps=(step,),)))
    @test_throws AlgorithmGraphError AlgorithmGraphs.graph_node_configuration(factory,
        merge(base,(steps=(step,),)),(transmission="calibration",unused="calibration"))
    legacy = factory(:legacy,merge(base,(steps=(),)),(;))
    @test legacy !== nothing

    for T in (Float32,Float64)
        graph = train_graph(config; T, opd=zeros(T,8,8), amplitude=ones(T,8,8),
            bindings=(calibration=fill(0.5,8,8),))
        step_graph!(graph)
        @test graph_output(graph,:intensity) ≈ fill(T(0.25),8,8)
        @test eltype(AlgorithmGraphs.prepared_graph_node(graph,:optics).workspace.masks[1]) === Complex{T}
        @test @inferred(step_graph!(graph)) === graph
    end

    amplitude = [0.1 0.3 0.7 0.2; 0.6 0.8 0.4 1.0; 0.9 0.2 0.5 0.0]
    phase = [0.0 0.3 -0.4 0.8; -0.6 1.1 0.7 0.0; 1.8 -0.9 0.2 -1.2]
    asymmetric = merge(step,(phase="retardance",sampling_m=(0.0015,0.0007),
        center_pixel=(1.25,2.7),outside_amplitude=0.2,outside_phase_rad=0.4))
    config = merge(base,(steps=(asymmetric,),transmission="calibration",retardance="phase_calibration"))
    expected = sampled_mask_reference(amplitude,phase,8,0.001,asymmetric.sampling_m,
        asymmetric.center_pixel,0.2cis(0.4))
    graph = train_graph(config;bindings=(calibration=amplitude,phase_calibration=phase))
    owner = AlgorithmGraphs.prepared_graph_node(graph,:optics)
    @test owner.workspace.masks[1] ≈ expected
    step_graph!(graph)
    @test graph_output(graph,:intensity) ≈ abs2.(expected)
    @test expected[5,5] != expected[4,6] # An asymmetric fixture detects axis exchange.

    # At a propagated plane, phase changes intensity even for fixed amplitude.
    propagation = ((op="lens",focal_length_m=0.1),(op="propagate",distance_m=0.1))
    propagated = train_graph(merge(config,(steps=(asymmetric,propagation...),));
        bindings=(calibration=amplitude,phase_calibration=phase))
    step_graph!(propagated)
    context = Proper.RunContext(Matrix{ComplexF64})
    wf = prop_begin(0.008,1.55e-6,8;beam_diam_fraction=1.0,context)
    prop_multiply(wf,expected)
    prop_lens(wf,0.1,context,"lens")
    prop_propagate(wf,0.1,context,"propagate")
    direct,_ = prop_end(wf)
    @test graph_output(propagated,:intensity) ≈ direct
    no_phase = train_graph(merge(config,(steps=(asymmetric,propagation...),));
        bindings=(calibration=amplitude,phase_calibration=zeros(size(phase))))
    step_graph!(no_phase)
    @test !isapprox(graph_output(propagated,:intensity),graph_output(no_phase,:intensity);rtol=1e-3)

    # Resampling at a focal plane uses its changed physical pitch.
    focal = train_graph(merge(config,(steps=(propagation...,asymmetric),));
        bindings=(calibration=amplitude,phase_calibration=phase))
    focal_mask = AlgorithmGraphs.prepared_graph_node(focal,:optics).workspace.masks[3]
    focal_expected = sampled_mask_reference(amplitude,phase,8,prop_get_sampling(wf),
        asymmetric.sampling_m,asymmetric.center_pixel,0.2cis(0.4))
    @test focal_mask ≈ focal_expected
    @test !isapprox(focal_mask,expected)

    # The middle phasor remains near -1 across the ±π angle boundary.
    wrap = merge(step,(phase="retardance",center_pixel=(1.0,1.5)))
    for T in (Float32,Float64)
        wrapped = train_graph(merge(base,(steps=(wrap,),transmission="calibration",retardance="phase_calibration"));
            T,opd=zeros(T,8,8),amplitude=ones(T,8,8),
            bindings=(calibration=ones(1,2),phase_calibration=[π-0.1 -π+0.1]))
        wrapped_mask = AlgorithmGraphs.prepared_graph_node(wrapped,:optics).workspace.masks[1]
        @test wrapped_mask[5,5] ≈ -T(cos(0.1))+0im
        step_graph!(wrapped)
        @test graph_output(wrapped,:intensity)[5,5] ≈ T(cos(0.1)^2)
    end

    # Two roles can refer to the same external array without duplicate loading.
    shared = merge(base,(steps=(merge(step,(phase="retardance",)),),transmission="same",retardance="same"))
    @test AlgorithmGraphs.graph_node_bindings(factory,shared) == (:same,)
    @test AlgorithmGraphs.graph_node_configuration(factory,
        merge(base,(steps=shared.steps,)),(transmission="same",retardance="same")) == shared

    # Constructor snapshots views before preparation and all subsequent stepping.
    input_amplitude, input_phase = copy(amplitude), copy(phase)
    node = factory(:optics,config,(;),(calibration=view(input_amplitude,:,:),phase_calibration=view(input_phase,:,:)))
    fill!(input_amplitude,0); fill!(input_phase,2)
    snapshot_graph = prepare_algorithm_graph(algorithm_graph((node,);name=:snapshot,
        inputs=(graph_input(:opd,:optics=>:pupil_opd,zeros(8,8)),
            graph_input(:amplitude,:optics=>:pupil_amplitude,ones(8,8))),
        outputs=(graph_output(:intensity,:optics=>:output),)))
    step_graph!(snapshot_graph)
    @test graph_output(snapshot_graph,:intensity) ≈ abs2.(expected)

    # Ordinary YAML bindings resolve role aliases; graph source is cold input.
    yaml = """
    schema_version: 1
    name: sampled_train
    nodes:
    - name: optics
      type: proper_optical_train_f64
      bindings: {transmission: calibrated_amplitude}
      config:
        resolution: 8
        pupil_resolution: 8
        diameter_m: 0.008
        wavelength_um: 1.55
        steps:
        - op: sampled_mask
          name: sampled apodizer
          amplitude: transmission
          sampling_m: [0.001, 0.001]
          center_pixel: [5.0, 5.0]
          outside_amplitude: 0.0
          outside_phase_rad: 0.0
    inputs:
    - {name: opd, destination: optics.pupil_opd, binding: opd}
    - {name: amplitude, destination: optics.pupil_amplitude, binding: amplitude}
    outputs:
    - {name: intensity, source: optics.output}
    """
    yaml_graph = mktempdir() do directory
        path = joinpath(directory,"sampled.yaml")
        write(path,yaml)
        graph = prepare_algorithm_graph(load_algorithm_graph(path;
            bindings=(opd=zeros(8,8),amplitude=ones(8,8),calibrated_amplitude=fill(0.5,8,8))))
        rm(path)
        graph
    end
    step_graph!(yaml_graph); step_graph!(yaml_graph)
    @test graph_output(yaml_graph,:intensity) ≈ fill(0.25,8,8)
    step_graph!(snapshot_graph)
    @test graph_output(snapshot_graph,:intensity) ≈ abs2.(expected)

    for bad in ((sampling_m=(0.0,0.001),),(sampling_m=(-1.0,1.0),),
            (sampling_m=(Inf,1.0),),(sampling_m=(1.0,),),(sampling_m=(true,1.0),),
            (center_pixel=(NaN,1.0),),(center_pixel=(1.0,Inf),),
            (center_pixel=(1.0,),),(outside_amplitude=-0.1,),(outside_amplitude=1.1,),
            (outside_amplitude=NaN,),(outside_phase_rad=Inf,),(unexpected=1,),
            (amplitude=1,),(phase=1,),(op="fits_mask",))
        @test_throws AlgorithmGraphError factory(:bad,merge(base,
            (steps=(merge(step,bad),),transmission="calibration")),(;),(calibration=amplitude,))
    end
    for omitted in (:name,:amplitude,:sampling_m,:center_pixel,:outside_amplitude,:outside_phase_rad)
        incomplete = NamedTuple{Tuple(k for k in keys(step) if k != omitted)}(
            Tuple(getproperty(step,k) for k in keys(step) if k != omitted))
        @test_throws AlgorithmGraphError factory(:bad,merge(base,
            (steps=(incomplete,),transmission="calibration")),(;),(calibration=amplitude,))
    end
    amplitude_config = merge(base,(steps=(step,),transmission="calibration"))
    for invalid in (fill(-0.1,2,2),fill(1.1,2,2),fill(NaN,2,2),fill(Inf,2,2),
            ones(ComplexF64,2,2),ones(2),zeros(0,2),OffsetMaskForTest())
        @test_throws AlgorithmGraphError factory(:bad,amplitude_config,(;),(calibration=invalid,))
    end
    for invalid in (zeros(2,2),fill(NaN,size(amplitude)),fill(Inf,size(amplitude)),
            ones(ComplexF64,size(amplitude)),ones(length(amplitude)))
        @test_throws AlgorithmGraphError factory(:bad,config,(;),(calibration=amplitude,phase_calibration=invalid))
    end
    @test_throws AlgorithmGraphError proper_optical_train_node_factory(T=Float32)(:bad,
        config,(;),(calibration=amplitude,phase_calibration=fill(1e100,size(amplitude))))
end

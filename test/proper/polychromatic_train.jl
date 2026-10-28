# Independent dense cell-intersection reference; production uses a sparse sweep.
function cell_overlap_reference(image, source_pitch, destination_pitch)
    n = size(image,1)
    center = fld(n,2)+1
    bounds(i,p) = ((i-center-0.5)*p, (i-center+0.5)*p)
    fractions = [begin
        slo,shi = bounds(s,source_pitch)
        dlo,dhi = bounds(d,destination_pitch)
        max(0.0,min(shi,dhi)-max(slo,dlo))/source_pitch
    end for d in 1:n, s in 1:n]
    fractions*image*transpose(fractions)
end

struct UnsupportedProperDevice <: AdaptiveOpticsSim.Backends.AbstractComputeDevice end

@testset "Polychromatic Proper optical train" begin
    extension = Base.get_extension(AdaptiveOpticsSim,:AdaptiveOpticsSimProperExt)
    config = (resolution=16,pupil_resolution=8,diameter_m=0.008,wavelength_um=1.5,
        steps=((op="lens",focal_length_m=0.2),(op="propagate",distance_m=0.2)))
    channel(λ,w;steps=config.steps,transmission=1.0) =
        (wavelength_um=λ,weight=w,transmission=transmission,steps=steps)
    for T in (Float32,Float64)
        opd = T[2e-8*(row-col) for row in 1:8,col in 1:8]
        amplitude = T[0.2+0.8*(row+2col)/24 for row in 1:8,col in 1:8]
        scalar = train_graph(config;T,opd,amplitude)
        single = train_graph(merge(config,(channels=(channel(1.5,7.0),),));T,opd,amplitude)
        step_graph!(scalar); step_graph!(single)
        @test graph_output(single,:intensity) == graph_output(scalar,:intensity)
        @test optical_sampling_m(AlgorithmGraphs.prepared_graph_node(single,:optics)) ≈
            optical_sampling_m(AlgorithmGraphs.prepared_graph_node(scalar,:optics))
        @test @inferred(step_graph!(single)) === single

        channels = (channel(1.0,2.0;steps=((op="multiply",amplitude=0.6),config.steps...)),
            channel(2.0,3.0;transmission=0.4))
        mixed = train_graph(merge(config,(;channels));T,opd,amplitude)
        owner = AlgorithmGraphs.prepared_graph_node(mixed,:optics)
        reference_pitch = optical_sampling_m(owner) # Prepared output pitch, before first step.
        step_graph!(mixed)
        expected = zeros(16,16)
        raw = Matrix{T}[]
        for c in channels
            scalar_channel = train_graph(merge(config,(wavelength_um=c.wavelength_um,steps=c.steps));T,opd,amplitude)
            step_graph!(scalar_channel)
            image = graph_output(scalar_channel,:intensity)
            push!(raw,copy(image))
            pitch = optical_sampling_m(AlgorithmGraphs.prepared_graph_node(scalar_channel,:optics))
            expected .+= (c.weight/5)*c.transmission.*cell_overlap_reference(image,pitch,reference_pitch)
        end
        @test !isapprox(raw[1],raw[2];rtol=1e-3)
        @test graph_output(mixed,:intensity) ≈ expected rtol=100eps(T)
        @test !isapprox(graph_output(mixed,:intensity),0.4raw[1]+0.24raw[2];rtol=1e-3)
        before = copy(graph_output(mixed,:intensity))
        reset_graph!(mixed)
        @test all(iszero,graph_output(mixed,:intensity))
        step_graph!(mixed)
        @test graph_output(mixed,:intensity) == before
        @test optical_sampling_m(owner) == reference_pitch
        @test @inferred(step_graph!(mixed)) === mixed
        opd .*= 3
        step_graph!(mixed)
        @test !isapprox(graph_output(mixed,:intensity),before;rtol=1e-3)
        @test_throws AlgorithmGraphError AlgorithmGraphs.prepare_graph_node(
            extension.PolychromaticOpticalTrainNode{T},owner.plan,(;),
            (pupil_opd=opd,pupil_amplitude=amplitude),(output=zeros(T,16,16),),(;),UnsupportedProperDevice())
    end

    # Exact power tests independent of propagation: expansion conserves total;
    # cropping removes edge-cell power without redistributing it.
    source = reshape(collect(1.0:64.0),8,8)
    for (sp,dp) in ((1.0,1.0),(1.0,2.0),(2.0,1.0),(1.25,1.0))
        overlaps = extension._cell_overlaps(8,sp,dp)
        output = zeros(8,8)
        extension._add_channel_intensity!(output,source,zeros(8,8),overlaps,0.3)
        @test output ≈ 0.3cell_overlap_reference(source,sp,dp)
        if dp >= sp
            @test sum(output) ≈ 0.3sum(source)
        else
            @test sum(output) < 0.3sum(source)
        end
    end
    edge = zeros(8,8); edge[1,1] = 1
    output = zeros(8,8)
    extension._add_channel_intensity!(output,edge,zeros(8,8),extension._cell_overlaps(8,2.0,1.0),1.0)
    @test iszero(sum(output))

    # A longer-wavelength PSF occupies more physical detector area after mapping.
    wide = merge(config,(resolution=64,pupil_resolution=16,wavelength_um=2.0))
    short = train_graph(merge(wide,(channels=(channel(1.0,1.0),),)))
    long = train_graph(merge(wide,(channels=(channel(2.0,1.0),),)))
    step_graph!(short); step_graph!(long)
    radius2 = [(row-33)^2+(col-33)^2 for row in 1:64,col in 1:64]
    extent(g) = sum(radius2.*graph_output(g,:intensity))/sum(graph_output(g,:intensity))
    @test extent(long) > 1.5extent(short)
    @test sum(graph_output(short,:intensity)) ≈ sum(graph_output(long,:intensity))
    filtered = train_graph(merge(wide,(channels=(channel(1.0,1.0;transmission=0.0),
        channel(2.0,1.0;transmission=0.3)),)))
    step_graph!(filtered)
    @test graph_output(filtered,:intensity) ≈ 0.15graph_output(long,:intensity)

    # Mask bindings are the union over the reference and all channel trains.
    mask(role) = (op="sampled_mask",name=role,amplitude=role,sampling_m=(0.001,0.001),
        center_pixel=(5.0,5.0),outside_amplitude=0.0,outside_phase_rad=0.0)
    masked_config = merge(config,(steps=(mask("reference_mask"),),
        channels=(channel(1.0,1.0;steps=(mask("blue_mask"),)),
            channel(2.0,3.0;steps=(mask("red_mask"),))),
        reference_mask="reference",blue_mask="blue",red_mask="red"))
    @test AlgorithmGraphs.graph_node_bindings(proper_optical_train_node_factory(),masked_config) ==
        (:reference,:blue,:red)
    bindings = (reference=ones(8,8),blue=fill(0.2,8,8),red=fill(0.6,8,8))
    masked = train_graph(masked_config;bindings)
    fill!(bindings.blue,0); fill!(bindings.red,0) # Cold snapshots retain the masks.
    step_graph!(masked)
    @test sum(graph_output(masked,:intensity)) ≈ 64*(0.25*0.2^2+0.75*0.6^2)
    reset_graph!(masked); step_graph!(masked)
    @test sum(graph_output(masked,:intensity)) ≈ 64*(0.25*0.2^2+0.75*0.6^2)

    yaml = """
    schema_version: 1
    name: polychromatic_train
    nodes:
    - name: optics
      type: proper_optical_train_f64
      bindings: {blue_mask: blue}
      config:
        resolution: 8
        pupil_resolution: 8
        diameter_m: 0.008
        wavelength_um: 1.5
        steps: []
        channels:
        - wavelength_um: 1.0
          weight: 2.0
          steps:
          - op: sampled_mask
            name: blue apodizer
            amplitude: blue_mask
            sampling_m: [0.001, 0.001]
            center_pixel: [5.0, 5.0]
            outside_amplitude: 0.0
            outside_phase_rad: 0.0
        - wavelength_um: 2.0
          weight: 6.0
          transmission: 0.5
          steps: [{op: multiply, amplitude: 0.25}]
    inputs:
    - {name: opd, destination: optics.pupil_opd, binding: opd}
    - {name: amplitude, destination: optics.pupil_amplitude, binding: amplitude}
    outputs:
    - {name: intensity, source: optics.output}
    """
    yaml_graph = mktempdir() do directory
        path = joinpath(directory,"channels.yaml")
        write(path,yaml)
        prepare_algorithm_graph(load_algorithm_graph(path;
            bindings=(opd=zeros(8,8),amplitude=ones(8,8),blue=fill(0.5,8,8))))
    end
    step_graph!(yaml_graph)
    @test graph_output(yaml_graph,:intensity) ≈ fill(0.25*0.5^2+0.75*0.5*0.25^2,8,8)
    @test_throws AlgorithmGraphError train_graph(masked_config;
        bindings=(blue=ones(8,8),red=ones(8,8)))

    factory = proper_optical_train_node_factory()
    good = channel(1.5,1.0)
    for channels in ((),ntuple(_ -> good,65),[good],(1,),
            ((wavelength_um=1.5,weight=1.0),),
            (merge(good,(unknown=1,)),),
            (merge(good,(wavelength_um=0.0,)),),
            (merge(good,(wavelength_um=Inf,)),),
            (merge(good,(wavelength_um=true,)),),
            (merge(good,(weight=-1.0,)),),
            (merge(good,(weight=NaN,)),),
            (merge(good,(weight=Inf,)),),
            (merge(good,(weight=true,)),),
            (merge(good,(weight=0.0,)),),
            (merge(good,(transmission=-0.1,)),),
            (merge(good,(transmission=1.1,)),),
            (merge(good,(transmission=Inf,)),),
            (merge(good,(transmission=NaN,)),),
            (merge(good,(transmission=true,)),),
            (merge(good,(steps=((op="unsupported",),),)),))
        @test_throws AlgorithmGraphError factory(:bad,merge(config,(;channels)),(;))
    end
    @test_throws AlgorithmGraphError proper_optical_train_node_factory(T=Float32)(:bad,
        merge(config,(channels=(merge(good,(wavelength_um=1e-100,)),),)),(;))
    @test_throws AlgorithmGraphError proper_optical_train_node_factory(T=Float32)(:bad,
        merge(config,(channels=(merge(good,(transmission=1+eps(Float64),)),),)),(;))
    # Proportional source weights may include dark channels; transmission defaults to one.
    default_response = (wavelength_um=1.5,weight=1.0,steps=config.steps)
    accepted = train_graph(merge(config,(channels=(default_response,merge(good,(weight=0.0,))),)))
    step_graph!(accepted)
    @test sum(graph_output(accepted,:intensity)) > 0
end

@testset "Fixed Proper detector sampling" begin
    extension = Base.get_extension(AdaptiveOpticsSim,:AdaptiveOpticsSimProperExt)
    entrance = (resolution=16,pupil_resolution=8,diameter_m=0.008,wavelength_um=1.5,steps=())
    focus = merge(entrance,(resolution=32,
        steps=((op="lens",focal_length_m=0.2),(op="propagate",distance_m=0.2))))
    for T in (Float32,Float64)
        opd = zeros(T,8,8)
        amplitude = T[0.1+0.9*(row+col)/16 for row in 1:8,col in 1:8]
        scalar = train_graph(entrance;T,opd,amplitude)
        fixed = train_graph(merge(entrance,(output_sampling_m=0.002,));T,opd,amplitude)
        scalar_owner = AlgorithmGraphs.prepared_graph_node(scalar,:optics)
        owner = AlgorithmGraphs.prepared_graph_node(fixed,:optics)
        @test scalar_owner isa extension.OpticalTrainOwner
        @test owner isa extension.PolychromaticOpticalTrainOwner
        @test optical_sampling_m(owner) === T(0.002)
        step_graph!(scalar); step_graph!(fixed)
        expected = cell_overlap_reference(graph_output(scalar,:intensity),
            optical_sampling_m(scalar_owner),T(0.002))
        @test graph_output(fixed,:intensity) ≈ expected rtol=100eps(T)
        @test sum(graph_output(fixed,:intensity)) ≈ sum(graph_output(scalar,:intensity))
        @test @inferred(step_graph!(fixed)) === fixed
        reset_graph!(fixed); step_graph!(fixed)
        @test graph_output(fixed,:intensity) ≈ expected rtol=100eps(T)

        # Both native focal-plane pitches differ from the fixed detector pitch.
        # Changing the selected central wavelength must never move that grid.
        pitch = T(4e-5)
        images = Matrix{T}[]
        native_pitches = T[]
        for wavelength in (1.0,2.0)
            config = merge(focus,(wavelength_um=wavelength,))
            native = train_graph(config;T,opd,amplitude)
            detector = train_graph(merge(config,(output_sampling_m=pitch,));T,opd,amplitude)
            detector_owner = AlgorithmGraphs.prepared_graph_node(detector,:optics)
            @test optical_sampling_m(detector_owner) === pitch
            step_graph!(native); step_graph!(detector)
            native_pitch = optical_sampling_m(AlgorithmGraphs.prepared_graph_node(native,:optics))
            push!(native_pitches,native_pitch)
            @test native_pitch != pitch
            mapped = cell_overlap_reference(graph_output(native,:intensity),native_pitch,pitch)
            @test graph_output(detector,:intensity) ≈ mapped rtol=100eps(T)
            @test sum(graph_output(detector,:intensity)) ≈ sum(graph_output(native,:intensity))
            push!(images,copy(graph_output(detector,:intensity)))
        end
        @test native_pitches[1] != native_pitches[2]
        @test !isapprox(images[1],images[2];rtol=1e-3)
        channels = ((wavelength_um=1.0,weight=1.0,steps=focus.steps),
            (wavelength_um=2.0,weight=3.0,steps=focus.steps))
        outputs = map((1.0,2.0)) do reference_wavelength
            graph = train_graph(merge(focus,(wavelength_um=reference_wavelength,
                output_sampling_m=pitch,channels=channels));T,opd,amplitude)
            @test optical_sampling_m(AlgorithmGraphs.prepared_graph_node(graph,:optics)) === pitch
            step_graph!(graph)
            copy(graph_output(graph,:intensity))
        end
        @test outputs[1] == outputs[2]
        @test outputs[1] ≈ T(0.25).*images[1]+T(0.75).*images[2] rtol=100eps(T)
    end
    for T in (Float32,Float64), invalid in (0.0,-0.001,NaN,Inf,true,"0.001",(0.001,))
        @test_throws AlgorithmGraphError proper_optical_train_node_factory(;T)(:bad,
            merge(entrance,(output_sampling_m=invalid,)),(;))
    end
    for invalid in (1e-100,1e100)
        @test_throws AlgorithmGraphError proper_optical_train_node_factory(T=Float32)(:bad,
            merge(entrance,(output_sampling_m=invalid,)),(;))
    end
    @test_throws AlgorithmGraphError proper_optical_train_node_factory()(:bad,
        merge(entrance,(output_sampling_mm=1.0,)),(;))
end

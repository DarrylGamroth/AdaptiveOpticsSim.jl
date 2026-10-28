# Include after external_pupil.jl, which supplies small optical/response graphs.
@testset "whole SCC exposures select static fringed and unfringed graphs" begin
    fringed=external_intensity_test_graph(Ref(8.0))
    unfringed=external_intensity_test_graph(Ref(2.0))
    llowfs=external_intensity_test_graph(Ref(3.0))
    response=external_response_graph()
    specification=(graph=response,input=:photon_rate,output=:frame,
        normalized_to_photon_rate=1e9,pixel_sampling_m=(1.0,1.0),
        exposure_duration_s=4e-9,frame_schema="test.frame/1")
    owner=prepare_external_pupil((;scc=fringed,llowfs);
        exposure_durations=(scc=ModelDuration(4),llowfs=ModelDuration(2)),
        observations=(scc=(output=:intensity,sampling_node=:optics,response=specification,
                chopping=(;unfringed,initial_phase=:fringed)),
            llowfs=(output=:intensity,sampling_node=:optics)))
    products=NamedTuple[]
    emit(name,image,metadata)=push!(products,(;name,image=copy(image),metadata))
    domain=ntuple(_->UInt8(1),16)
    for (sequence,(time,span)) in enumerate(((0,1),(1,2),(3,2),(5,2),(7,2)))
        @test step_external_pupil!(emit,owner,zeros(4,4),ones(4,4),
            ModelTimestamp(time),ModelDuration(span);domain,generation=UInt64(1),sequence=UInt64(sequence))
    end
    scc=filter(p->p.name==:scc,products)
    ll=filter(p->p.name==:llowfs,products)
    @test [get(p.metadata,:chopped_image_kind,nothing) for p in scc]==[:fringed,:unfringed]
    @test [p.image[1] for p in scc]≈[16f0,4f0]
    @test [p.metadata.pts for p in scc]==[0,4]
    @test all(p->p.image==fill(3.0,4,4) && !haskey(p.metadata,:chopped_image_kind),ll)
    @test length(ll)==4 && llowfs.state.step_sequence==5
    @test response.state.step_sequence==2
    @test unfringed.state.step_sequence==3
end

function chopped_test_owner(; initial_phase=:fringed, fringed=external_intensity_test_graph(Ref(8.0)),
        unfringed=external_intensity_test_graph(Ref(2.0)))
    prepare_external_pupil((scc=fringed,);exposure_durations=(scc=ModelDuration(4),),
        observations=(scc=(output=:intensity,sampling_node=:optics,
            chopping=(;unfringed,initial_phase)),))
end

@testset "chopping follows model time through gaps, generations and reset" begin
    owner=chopped_test_owner()
    products=NamedTuple[]
    emit(name,image,metadata)=push!(products,(;image=copy(image),metadata))
    domain=ntuple(_->UInt8(1),16)
    feed(time,span,sequence; generation=1,lost=0)=step_external_pupil!(emit,owner,
        zeros(4,4),ones(4,4),ModelTimestamp(time),ModelDuration(span);
        domain,generation=UInt64(generation),sequence=UInt64(sequence),lost_samples=UInt64(lost))
    @test feed(0,4,1)
    @test feed(8,4,3;lost=1) # Missing unfringed exposure must not invert the phase.
    @test [p.metadata.chopped_image_kind for p in products]==[:fringed,:fringed]
    @test [p.metadata.sequence for p in products]==[1,2]
    @test feed(14,2,5;lost=2) # Incomplete unfringed interval cannot be published.
    @test length(products)==2
    @test feed(16,4,6;lost=2)
    @test last(products).metadata.pts==16 && last(products).metadata.chopped_image_kind===:fringed
    @test feed(0,4,1;generation=2)
    @test last(products).metadata.sequence==1 && last(products).metadata.chopped_image_kind===:fringed
    @test owner.graphs.scc.state.step_sequence==1
    @test owner.observations.scc.chopping.unfringed.state.step_sequence==0
    reset_external_pupil!(owner)
    @test owner.observations.scc.chopping.unfringed.state.step_sequence==0
    @test feed(0,4,1;generation=3)
    @test last(products).metadata.sequence==1 && last(products).image==fill(8.0,4,4)
    reversed=chopped_test_owner(;initial_phase=:unfringed)
    phases=Symbol[]
    for (sequence,time) in enumerate((0,4))
        @test step_external_pupil!((name,image,metadata)->push!(phases,metadata.chopped_image_kind),
            reversed,zeros(4,4),ones(4,4),ModelTimestamp(time),ModelDuration(4);
            domain,generation=UInt64(1),sequence=UInt64(sequence))
    end
    @test phases==[:unfringed,:fringed]
end

@testset "chopping validates both variants before exposure mutation" begin
    value=Ref(2.0)
    owner=chopped_test_owner(;unfringed=external_intensity_test_graph(value))
    domain=ntuple(_->UInt8(1),16)
    products=Any[]
    feed(time,span,sequence)=step_external_pupil!((name,image,metadata)->push!(products,metadata),
        owner,zeros(4,4),ones(4,4),ModelTimestamp(time),ModelDuration(span);
        domain,generation=UInt64(1),sequence=UInt64(sequence))
    @test feed(0,3,1)
    prior=exposure_state(owner.exposures.scc)
    for invalid in (-1.0,NaN,Inf)
        value[]=invalid
        @test_throws AlgorithmGraphError feed(3,2,2)
        @test isequal(exposure_state(owner.exposures.scc),prior)
        @test isempty(products) && owner.state.sequence==1 && owner.state.last_end_ns==3
    end
    value[]=2
    @test feed(3,2,2)
    @test only(products).chopped_image_kind===:fringed

    graph=external_intensity_test_graph(Ref(8.0))
    for chopping in ((;unfringed=graph,initial_phase=:fringed),
            (;unfringed=nothing,initial_phase=:fringed),
            (;unfringed=external_intensity_test_graph(Ref(2.0)),initial_phase=:other),
            (;unfringed=external_intensity_test_graph(Ref(2.0))),
            (;unfringed=external_intensity_test_graph(Ref(2.0)),initial_phase=:fringed,other=1))
        @test_throws AlgorithmGraphError prepare_external_pupil((scc=graph,);
            exposure_durations=(scc=ModelDuration(4),),
            observations=(scc=(output=:intensity,sampling_node=:optics,chopping),))
    end
    shared=graph_input(graph,:pupil_opd)
    @test_throws AlgorithmGraphError chopped_test_owner(;fringed=graph,
        unfringed=external_test_graph(opd=shared))
    @test_throws AlgorithmGraphError chopped_test_owner(;fringed=graph,
        unfringed=external_test_graph()) # Same shape but incompatible physical sampling.
    config=(resolution=8,pupil_resolution=4,diameter_m=0.004,wavelength_um=1.55,steps=())
    node=proper_optical_train_node_factory()(:optics,config,(;))
    wrong_shape=prepare_algorithm_graph(algorithm_graph((node,);
        inputs=(graph_input(:pupil_opd,:optics=>:pupil_opd,zeros(4,4)),
            graph_input(:pupil_amplitude,:optics=>:pupil_amplitude,ones(4,4))),
        outputs=(graph_output(:intensity,:optics=>:output),)))
    @test_throws AlgorithmGraphError chopped_test_owner(;unfringed=wrong_shape)
end

function chopped_spectral_graph(reference_open)
    main=(radius_m=0.002,x_m=0.0,y_m=0.0)
    reference=(radius_m=0.0006,x_m=0.003,y_m=0.0)
    circles=reference_open ? (main,reference) : (main,)
    steps=((op="circular_mask",mode="transmit",circles),
        (op="lens",focal_length_m=0.1),(op="propagate",distance_m=0.1))
    config=(resolution=16,pupil_resolution=8,diameter_m=0.008,wavelength_um=1.55,
        steps,channels=((wavelength_um=1.5,weight=1.0,steps),
            (wavelength_um=1.6,weight=2.0,steps)),output_sampling_m=1e-5)
    node=proper_optical_train_node_factory()(:optics,config,(;))
    prepare_algorithm_graph(algorithm_graph((node,);
        inputs=(graph_input(:pupil_opd,:optics=>:pupil_opd,zeros(8,8)),
            graph_input(:pupil_amplitude,:optics=>:pupil_amplitude,ones(8,8))),
        outputs=(graph_output(:intensity,:optics=>:output),)))
end

@testset "static spectral variants gate the reference field before intensity" begin
    fringed=chopped_spectral_graph(true)
    unfringed=chopped_spectral_graph(false)
    expected=map((true,false)) do open
        graph=chopped_spectral_graph(open)
        step_graph!(graph)
        copy(graph_output(graph,:intensity))
    end
    @test !isapprox(expected[1],expected[2];rtol=1e-3)
    owner=chopped_test_owner(;fringed,unfringed)
    images=Matrix{Float64}[];kinds=Symbol[]
    domain=ntuple(_->UInt8(1),16)
    for (sequence,time) in enumerate((0,4))
        @test step_external_pupil!(owner,zeros(8,8),ones(8,8),ModelTimestamp(time),ModelDuration(4);
                domain,generation=UInt64(1),sequence=UInt64(sequence)) do name,image,metadata
            push!(images,copy(image));push!(kinds,metadata.chopped_image_kind)
        end
    end
    @test kinds==[:fringed,:unfringed]
    @test images[1]≈expected[1] rtol=1e-12
    @test images[2]≈expected[2] rtol=1e-12
end

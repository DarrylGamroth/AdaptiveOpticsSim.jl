struct ExposureExecutionTestNode end
struct ExposureExecutionTestOwner{V,O}
    value::V
    output::O
end
function AlgorithmGraphs.graph_node_ports(::Type{ExposureExecutionTestNode}, value)
    port(name,direction) = AlgorithmGraphs.graph_port_contract(name,direction,:data,
        Float32,(4,4),"test.external-intensity/1",:column_major)
    return (port(:pupil_opd,:input),port(:pupil_amplitude,:input),port(:output,:output))
end
AlgorithmGraphs.prepare_graph_node(::Type{ExposureExecutionTestNode},value,props,
    inputs::NamedTuple,outputs::NamedTuple,parameters::NamedTuple,target) =
    ExposureExecutionTestOwner(value,outputs.output)
AlgorithmGraphs.step_graph_node!(owner::ExposureExecutionTestOwner) = (fill!(owner.output,owner.value[]);nothing)
AlgorithmGraphs.reset_graph_node!(owner::ExposureExecutionTestOwner) = (fill!(owner.output,0);nothing)
AlgorithmGraphs.optical_sampling_m(::ExposureExecutionTestOwner) = 1.0

function exposure_execution_optical(target,value)
    node = algorithm_node(:optics,ExposureExecutionTestNode,value)
    opd = Backends.allocate_device_array(target,Float32,4,4); fill!(opd,0)
    amplitude = Backends.allocate_device_array(target,Float32,4,4); fill!(amplitude,1)
    return prepare_algorithm_graph(algorithm_graph((node,);
        inputs=(graph_input(:pupil_opd,:optics=>:pupil_opd,opd),
            graph_input(:pupil_amplitude,:optics=>:pupil_amplitude,amplitude)),
        outputs=(graph_output(:intensity,:optics=>:output),));target)
end

function exposure_execution_owner(target;captured=true,scale=1e9)
    values = (first=Ref(2f0),fringed=Ref(3f0),unfringed=Ref(7f0))
    first = exposure_execution_optical(target,values.first)
    second = exposure_execution_optical(target,values.fringed)
    alternate = exposure_execution_optical(target,values.unfringed)
    detector = ingaas_detector_acquisition_node(:detector;rows=4,columns=4,
        pixel_sampling_m=(1.0,1.0),spectral_channel=:reference,
        exposure_duration_s=4e-9,quantum_efficiency=0.5,photon_noise=false,
        rng_seed=17,photon_rate_schema="test.rate/1",frame_schema="test.frame/1",T=Float32)
    rate = Backends.allocate_device_array(target,Float32,4,4); fill!(rate,0)
    response = prepare_algorithm_graph(algorithm_graph((detector,);
        inputs=(graph_input(:photon_rate,:detector=>:photon_rate,rate),),
        outputs=(graph_output(:frame,:detector=>:frame),));target)
    spec = (graph=response,input=:photon_rate,output=:frame,
        normalized_to_photon_rate=scale,pixel_sampling_m=(1.0,1.0),
        exposure_duration_s=4e-9,frame_schema="test.frame/1")
    owner = prepare_external_pupil((;first,second);
        exposure_durations=(first=ModelDuration(4),second=ModelDuration(4)),
        observations=(first=(output=:intensity,sampling_node=:optics,response=spec),
            second=(output=:intensity,sampling_node=:optics,
                chopping=(unfringed=alternate,initial_phase=:fringed))))
    captured && return owner,values
    # The original device path is the numerical and callback-order oracle.
    execution = map(variants -> map(_ -> nothing,variants),owner.exposure_execution)
    return PreparedExternalPupil(owner.graphs,owner.exposures,owner.pupil_opd,
        owner.pupil_amplitude,owner.observations,owner.input_names,owner.propagation,
        execution,owner.state),values
end

function execution_snapshot(owner)
    exposures = map(owner.exposures) do e
        (e.start_ns,e.covered_ns,e.next_ns,Array(e.integrated),Array(e.image),e.sequence,e.incomplete)
    end
    s = owner.state
    return (exposures,(s.last_end_ns,s.generation,s.sequence,s.domain,s.lost_samples,s.discontinuities,s.rejected))
end
function execution_feed!(emit,owner,start,span,sequence;generation=UInt64(1))
    step_external_pupil!(emit,owner,zeros(Float32,4,4),ones(Float32,4,4),
        ModelTimestamp(start),ModelDuration(span);domain=ntuple(_ -> UInt8(1),16),
        generation,sequence=UInt64(sequence))
end
function test_external_exposure_pipeline(target)
    @testset "external exposure response, chopping, admission and reset" begin
        captured,values = exposure_execution_owner(target)
        original,reference_values = exposure_execution_owner(target;captured=false)
        products = NamedTuple[];reference = NamedTuple[]
        collect_product(destination) = (name,image,metadata) -> push!(destination,(;name,image=Array(image),metadata))
        for (sequence,(start,span,value)) in enumerate(((0,1,2f0),(1,2,4f0),(3,2,6f0),(5,2,8f0),(7,1,10f0)))
            values.first[] = reference_values.first[] = value
            @test execution_feed!(collect_product(products),captured,start,span,sequence)
            @test execution_feed!(collect_product(reference),original,start,span,sequence)
            @test execution_snapshot(captured) == execution_snapshot(original)
            @test products == reference
        end
        @test [p.image[1] for p in products if p.name == :second] == [3.0,7.0]
        @test [p.metadata.chopped_image_kind for p in products if p.name == :second] == [:fringed,:unfringed]
        @test Array(graph_input(captured.observations.first.response.graph,:photon_rate)) ==
            Array(graph_input(original.observations.first.response.graph,:photon_rate))
        before = execution_snapshot(captured)
        values.first[] = 1f35
        @test_throws AlgorithmGraphError execution_feed!((args...) -> nothing,captured,8,2,6)
        @test execution_snapshot(captured) == before
        @test length(products) == 4
        reset_external_pupil!(captured);reset_external_pupil!(original)
        @test execution_snapshot(captured) == execution_snapshot(original)
        @test all(e -> all(iszero,Array(e.integrated)) && all(iszero,Array(e.image)),Base.values(captured.exposures))
    end
end

function test_external_exposure_callback_failure(target)
    @testset "external callback throw preserves original partial commit and reset" begin
        snapshots = map((false,true)) do captured
            owner,_ = exposure_execution_owner(target;captured)
            @test execution_feed!((args...) -> nothing,owner,0,2,1)
            before_later = execution_snapshot(owner)[1].second
            @test_throws ErrorException execution_feed!((args...) -> error("consumer rejected frame"),owner,2,2,2)
            after = execution_snapshot(owner)
            @test after[1].second == before_later
            @test owner.exposures.first.covered_ns == 4 && owner.exposures.first.start_ns == 0
            @test owner.exposures.first.sequence == 1 && owner.state.sequence == 1
            reset_external_pupil!(owner)
            @test all(e -> e.sequence == 0 && all(iszero,Array(e.image)),values(owner.exposures))
            return after
        end
        @test snapshots[1] == snapshots[2]
    end
end

function test_external_exposure_callback_mutation(target,mutate!,finish!)
    @testset "callback writes to a later observation retain managed stream visibility" begin
        for unrelated_stream in (false,true)
            snapshots = map((false,true)) do captured
                owner,_ = exposure_execution_owner(target;captured)
                @test execution_feed!((args...) -> nothing,owner,0,2,1)
                before_later = execution_snapshot(owner)[1].second
                callback = function (name,image,metadata)
                    name == :first || return nothing
                    later = graph_output(owner.graphs.second,:intensity)
                    mutate!(later,unrelated_stream)
                    return nothing
                end
                @test_throws AlgorithmGraphError execution_feed!(callback,owner,2,2,2)
                finish!()
                after = execution_snapshot(owner)
                @test after[1].second == before_later
                @test owner.state.sequence == 1 && owner.state.last_end_ns == 2
                @test owner.exposures.first.sequence == 1 && owner.exposures.first.next_ns == 4
                return after
            end
            @test snapshots[1] == snapshots[2]
        end
    end
end

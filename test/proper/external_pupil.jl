using AdaptiveOpticsSim.AlgorithmGraphs
using Test

struct ExternalIntensityTestNode end
struct ExternalIntensityTestOwner{V,O}
    value::V
    output::O
end
function AlgorithmGraphs.graph_node_ports(::Type{ExternalIntensityTestNode}, value)
    port(name,direction) = AlgorithmGraphs.graph_port_contract(name,direction,:data,
        Float64,(4,4),"test.external-intensity/1",:column_major)
    return (port(:pupil_opd,:input),port(:pupil_amplitude,:input),port(:output,:output))
end
AlgorithmGraphs.prepare_graph_node(::Type{ExternalIntensityTestNode},value,props,
        inputs::NamedTuple,outputs::NamedTuple,parameters::NamedTuple,target) =
    ExternalIntensityTestOwner(value,outputs.output)
AlgorithmGraphs.step_graph_node!(owner::ExternalIntensityTestOwner) = (fill!(owner.output,owner.value[]); nothing)
AlgorithmGraphs.reset_graph_node!(::ExternalIntensityTestOwner) = nothing
AlgorithmGraphs.optical_sampling_m(::ExternalIntensityTestOwner) = 1.0

function external_intensity_test_graph(value)
    node = algorithm_node(:optics,ExternalIntensityTestNode,value)
    return prepare_algorithm_graph(algorithm_graph((node,);
        inputs=(graph_input(:pupil_opd,:optics=>:pupil_opd,zeros(4,4)),
            graph_input(:pupil_amplitude,:optics=>:pupil_amplitude,ones(4,4))),
        outputs=(graph_output(:intensity,:optics=>:output),)))
end

exposure_state(exposure) = (exposure.start_ns,exposure.covered_ns,exposure.next_ns,
    copy(exposure.integrated),copy(exposure.image),exposure.sequence,exposure.incomplete)

function external_test_graph(; opd=zeros(4,4), amplitude=zeros(4,4))
    config = (resolution=4,pupil_resolution=4,diameter_m=0.004,wavelength_um=1.55,steps=())
    node = proper_optical_train_node_factory()(:optics,config,(;))
    return prepare_algorithm_graph(algorithm_graph((node,); name=:external_test,
        inputs=(graph_input(:pupil_opd,:optics=>:pupil_opd,opd),
            graph_input(:pupil_amplitude,:optics=>:pupil_amplitude,amplitude)),
        outputs=(graph_output(:intensity,:optics=>:output),)))
end

function external_test_owner(graphs=(slow=external_test_graph(),fast=external_test_graph()))
    return prepare_external_pupil(graphs;
        exposure_durations=(slow=ModelDuration(10),fast=ModelDuration(4)),
        observations=(slow=(output=:intensity,sampling_node=:optics),
            fast=(output=:intensity,sampling_node=:optics)))
end

@testset "Transport-neutral external pupil exposures" begin
    owner = external_test_owner()
    samples = NamedTuple[]
    emit = (name,image,metadata) -> push!(samples,(;name,image=copy(image),metadata))
    opd = zeros(4,4)
    amplitude = ones(4,4)
    domain = ntuple(_->UInt8(1),16)
    function feed(time,sequence,value; generation=UInt64(1),lost_samples=UInt64(0))
        fill!(amplitude,sqrt(value))
        step_external_pupil!(emit,owner,opd,amplitude,ModelTimestamp(time),ModelDuration(2);
            domain,generation,sequence=UInt64(sequence),lost_samples)
    end
    for (sequence,value) in enumerate((1.0,3.0,5.0,7.0,9.0))
        @test feed(2(sequence-1),sequence,value)
    end
    slow = filter(s->s.name==:slow,samples)
    fast = filter(s->s.name==:fast,samples)
    @test length(slow)==1 && length(fast)==2
    @test slow[1].image ≈ fill(5.0,4,4)
    @test fast[1].image ≈ fill(2.0,4,4)
    @test fast[2].image ≈ fill(6.0,4,4)
    @test slow[1].metadata.pts==0 && slow[1].metadata.duration==10
    @test fast[2].metadata.pts==4 && fast[2].metadata.duration==4
    @test slow[1].metadata.sampling_m ≈ 0.001
    @test !feed(8,5,9.0) # duplicate
    @test !feed(9,6,9.0) # overlapping interval
    emitted = length(samples)
    @test feed(12,7,11.0) # missing interval and sequence discard partial coverage
    @test length(samples)==emitted
    @test feed(14,8,13.0)
    @test last(samples).name==:fast
    @test last(samples).metadata.pts==12
    @test last(samples).image ≈ fill(12.0,4,4)
    @test feed(0,1,2.0;generation=UInt64(2))
    @test feed(2,2,4.0;generation=UInt64(2))
    @test last(samples).metadata.generation==2 && last(samples).metadata.sequence==1
    @test last(samples).image ≈ fill(3.0,4,4)
    @test !feed(100,100,1.0;generation=UInt64(1))
    @test !step_external_pupil!(emit,owner,opd,amplitude,ModelTimestamp(4),ModelDuration(2);
        domain=ntuple(_->UInt8(2),16),generation=UInt64(2),sequence=UInt64(3))
    fill!(amplitude,NaN)
    @test_throws AlgorithmGraphError step_external_pupil!(emit,owner,opd,amplitude,
        ModelTimestamp(4),ModelDuration(2);domain,generation=UInt64(2),sequence=UInt64(3))
    @test owner.state.sequence==2 && owner.state.last_end_ns==4
    @test owner.exposures.slow.incomplete>0
    @test feed(4,3,6.0;generation=UInt64(2),lost_samples=UInt64(1))
    @test owner.exposures.slow.covered_ns==2 # counter-only loss discards earlier partial coverage
    @test feed(6,4,8.0;generation=UInt64(2),lost_samples=UInt64(1))
    @test last(samples).image ≈ fill(7.0,4,4)

    # Cross-quantity aliases would destroy the coherent OPD/amplitude pair.
    shared = zeros(4,4)
    @test_throws AlgorithmGraphError external_test_owner((
        slow=external_test_graph(opd=shared,amplitude=shared),fast=external_test_graph()))
    @test_throws AlgorithmGraphError external_test_owner((
        slow=external_test_graph(amplitude=shared),fast=external_test_graph(opd=shared)))
    first_graph = external_test_graph()
    @test_throws AlgorithmGraphError external_test_owner((slow=first_graph,
        fast=external_test_graph(opd=graph_output(first_graph,:intensity))))

    # A later invalid graph must not advance any observation's integration.
    value = Ref(1.0)
    invalid_owner = external_test_owner((slow=external_test_graph(),
        fast=external_intensity_test_graph(value)))
    accept = (name,image,metadata) -> nothing
    @test step_external_pupil!(accept,invalid_owner,opd,ones(4,4),ModelTimestamp(0),ModelDuration(2);
        domain,generation=UInt64(1),sequence=UInt64(1))
    prior = map(exposure_state,invalid_owner.exposures)
    for invalid in (-1.0,NaN,Inf)
        value[] = invalid
        @test_throws AlgorithmGraphError step_external_pupil!(accept,invalid_owner,opd,ones(4,4),
            ModelTimestamp(2),ModelDuration(2);domain,generation=UInt64(1),sequence=UInt64(2))
        @test isequal(map(exposure_state,invalid_owner.exposures),prior)
        @test invalid_owner.state.sequence==1 && invalid_owner.state.last_end_ns==2
    end
end

@testset "Exactly covered intensity means" begin
    images = Matrix{Float64}[]
    emit = (image,args...) -> push!(images,copy(image))
    exposure = IntensityExposure((1,1),ModelDuration(10))
    integrate_intensity!(emit,exposure,fill(1.0,1,1),ModelTimestamp(0),ModelDuration(2))
    integrate_intensity!(emit,exposure,fill(3.0,1,1),ModelTimestamp(2),ModelDuration(8))
    @test only(images) ≈ fill(2.6,1,1)
    reset_intensity_exposure!(exposure,ModelTimestamp(13))
    integrate_intensity!(emit,exposure,fill(9.0,1,1),ModelTimestamp(13),ModelDuration(7))
    @test length(images)==1
    @test exposure.incomplete==1
    integrate_intensity!(emit,exposure,fill(5.0,1,1),ModelTimestamp(20),ModelDuration(10))
    @test length(images)==2 && last(images)==fill(5.0,1,1)

    # Covered duration alone cannot prove that the interval had no overlap/gap.
    broken = IntensityExposure((1,1),ModelDuration(10))
    integrate_intensity!(emit,broken,ones(1,1),ModelTimestamp(0),ModelDuration(4))
    @test_throws AlgorithmGraphError integrate_intensity!(emit,broken,ones(1,1),
        ModelTimestamp(2),ModelDuration(4))
    @test_throws AlgorithmGraphError integrate_intensity!(emit,broken,ones(1,1),
        ModelTimestamp(8),ModelDuration(2))

    validated = IntensityExposure((1,1),ModelDuration(2))
    integrate_intensity!(emit,validated,ones(1,1),ModelTimestamp(0),ModelDuration(1))
    prior = exposure_state(validated)
    for invalid in (-1.0,NaN,Inf,1.0+0.0im,big"1e4000")
        @test_throws AlgorithmGraphError integrate_intensity!(emit,validated,fill(invalid,1,1),
            ModelTimestamp(1),ModelDuration(1))
        @test isequal(exposure_state(validated),prior)
    end
    # Emitting into image would overwrite an aliased source before later overlaps.
    validated.image[1] = 3
    prior = exposure_state(validated)
    for input in (validated.image,validated.integrated,view(validated.image,:,:),
            view(validated.integrated,:,:))
        @test_throws AlgorithmGraphError integrate_intensity!(emit,validated,input,
            ModelTimestamp(1),ModelDuration(3))
        @test isequal(exposure_state(validated),prior)
    end
    integrate_intensity!(emit,validated,fill(3.0,1,1),ModelTimestamp(1),ModelDuration(3))
    @test images[end-1:end] == [fill(2.0,1,1),fill(3.0,1,1)]
end

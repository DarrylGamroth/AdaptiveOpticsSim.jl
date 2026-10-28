# Explicit optional selector: julia --project=<AMDGPU/Proper environment> test/amdgpu_external_pupil.jl
using Test, AdaptiveOpticsSim, Proper, AMDGPU
using AdaptiveOpticsSim.AlgorithmGraphs
import AdaptiveOpticsSim.Backends

AMDGPU.allowscalar(false)
const TARGET = Backends.AcceleratorComputeDevice(Backends.AMDGPUBackend(),1)
@test Backends.compute_device_is_available(Backends.compute_device_availability(TARGET))

function optical_graph(target; spectrum=false)
    steps=((op="lens",focal_length_m=0.14),(op="propagate",distance_m=0.14))
    config=(resolution=32,pupil_resolution=16,diameter_m=0.008,wavelength_um=1.55,steps)
    if spectrum
        config=merge(config,(channels=((wavelength_um=1.5,weight=1.0,steps),
            (wavelength_um=1.6,weight=2.0,transmission=0.4,steps)),))
    end
    node=proper_optical_train_node_factory()(:optics,config,(;))
    opd=Backends.allocate_device_array(target,Float64,16,16);fill!(opd,0)
    amplitude=similar(opd);fill!(amplitude,0)
    return prepare_algorithm_graph(algorithm_graph((node,);
        inputs=(graph_input(:pupil_opd,:optics=>:pupil_opd,opd),
            graph_input(:pupil_amplitude,:optics=>:pupil_amplitude,amplitude)),
        outputs=(graph_output(:intensity,:optics=>:output),));target)
end

function optical_owner(target; spectrum=false, response=false)
    graphs=(slow=optical_graph(target;spectrum),fast=optical_graph(target;spectrum))
    durations=(slow=ModelDuration(10),fast=ModelDuration(4))
    observations=map(keys(graphs)) do name
        selection=(output=:intensity,sampling_node=:optics)
        response || return selection
        graph=getproperty(graphs,name)
        step_graph!(graph)
        pitch=optical_sampling_m(AlgorithmGraphs.prepared_graph_node(graph,:optics))
        reset_graph!(graph)
        duration=AlgorithmGraphs.model_nanoseconds(getproperty(durations,name))*1e-9
        detector=ingaas_detector_acquisition_node(:detector;rows=32,columns=32,
            pixel_sampling_m=(pitch,pitch),spectral_channel=:reference,
            exposure_duration_s=duration,quantum_efficiency=0.5,
            photon_noise=false,rng_seed=17,photon_rate_schema="test.rate/1",
            frame_schema="test.frame/1",T=Float32)
        input=Backends.allocate_device_array(target,Float32,32,32);fill!(input,0)
        detector_graph=prepare_algorithm_graph(algorithm_graph((detector,);
            inputs=(graph_input(:photon_rate,:detector=>:photon_rate,input),),
            outputs=(graph_output(:frame,:detector=>:frame),));target)
        merge(selection,(response=(graph=detector_graph,input=:photon_rate,output=:frame,
            normalized_to_photon_rate=1e9,pixel_sampling_m=(pitch,pitch),
            exposure_duration_s=duration,frame_schema="test.frame/1"),))
    end
    return prepare_external_pupil(graphs;exposure_durations=durations,
        observations=NamedTuple{keys(graphs)}(observations))
end

@testset "AMDGPU external-pupil propagation and exposures" begin
    for spectrum in (false,true), response in (false,true)
        cpu=optical_owner(Backends.HostComputeDevice();spectrum,response)
        gpu=optical_owner(TARGET;spectrum,response)
        @test all(e->Backends.compute_device(e.integrated)==TARGET,values(gpu.exposures))
        @test all(e->Backends.compute_device(e.image)==TARGET,values(gpu.exposures))
        host_products=NamedTuple[];device_products=NamedTuple[]
        opd=[3e-8sin(0.3i+0.7j) for i in 1:16,j in 1:16]
        amplitude=[0.2+0.6sin(0.1i+0.3j)^2 for i in 1:16,j in 1:16]
        domain=ntuple(_->UInt8(1),16)
        function feed(owner,products,time,sequence,generation)
            return step_external_pupil!(owner,opd,amplitude,ModelTimestamp(time),ModelDuration(2);
                    domain,generation=UInt64(generation),sequence=UInt64(sequence)) do name,image,metadata
                @test Backends.compute_device(image)==(owner===gpu ? TARGET : Backends.HostComputeDevice())
                push!(products,(;name,image=Array(image),metadata))
            end
        end
        for (time,sequence,generation) in ((0,1,1),(2,2,1),(4,3,1),(6,4,1),(8,5,1),
                (12,7,1),(14,8,1),(0,1,2),(2,2,2))
            amplitude .*= 0.95
            @test feed(cpu,host_products,time,sequence,generation)
            @test feed(gpu,device_products,time,sequence,generation)
        end
        @test length(host_products)==length(device_products)==5
        for (host,device) in zip(host_products,device_products)
            @test host.name==device.name
            @test host.metadata.pts==device.metadata.pts
            @test host.metadata.duration==device.metadata.duration
            @test host.metadata.generation==device.metadata.generation
            @test host.metadata.sequence==device.metadata.sequence
            @test host.metadata.sampling_m≈device.metadata.sampling_m
            @test host.image≈device.image rtol=(response ? 2e-6 : 2e-10) atol=1e-11
            @test all(isfinite,device.image) && sum(device.image)>0
        end
        @test !feed(gpu,device_products,2,2,2) # duplicate
        @test !feed(gpu,device_products,100,100,1) # old generation
        @test gpu.state.discontinuities==cpu.state.discontinuities
        @test gpu.exposures.slow.incomplete==cpu.exposures.slow.incomplete
        @test_throws AlgorithmGraphError step_external_pupil!((args...)->nothing,gpu,
            AMDGPU.ROCArray(opd),AMDGPU.ROCArray(amplitude),ModelTimestamp(4),ModelDuration(2);
            domain,generation=UInt64(2),sequence=UInt64(3))
    end
    @test_throws AlgorithmGraphError prepare_external_pupil(
        (slow=optical_graph(TARGET),fast=optical_graph(Backends.HostComputeDevice()));
        exposure_durations=(slow=ModelDuration(10),fast=ModelDuration(4)),
        observations=(slow=(output=:intensity,sampling_node=:optics),fast=(output=:intensity,sampling_node=:optics)))
end

@testset "AMDGPU intensity validation precedes exposure mutation" begin
    exposure=IntensityExposure((4,4),ModelDuration(4);target=TARGET)
    products=Matrix{Float64}[]
    emit=(image,args...)->push!(products,Array(image))
    a=AMDGPU.ROCArray(fill(1.0,4,4))
    integrate_intensity!(emit,exposure,a,ModelTimestamp(0),ModelDuration(2))
    fill!(a,3)
    integrate_intensity!(emit,exposure,a,ModelTimestamp(2),ModelDuration(2))
    @test only(products)==fill(2.0,4,4)
    for value in (-1.0,NaN,Inf)
        fill!(a,value)
        @test_throws AlgorithmGraphError integrate_intensity!(emit,exposure,a,ModelTimestamp(4),ModelDuration(2))
        @test exposure.next_ns==4 && exposure.covered_ns==0 && exposure.sequence==1
    end
    @test_throws AlgorithmGraphError integrate_intensity!(emit,exposure,zeros(4,4),ModelTimestamp(4),ModelDuration(2))
end
@test !any(p.name=="CUDA" for p in keys(Base.loaded_modules))

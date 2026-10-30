# Explicit optional selector: julia --project=<source AMDGPU environment> test/amdgpu_pixel_integration.jl
using Test, AdaptiveOpticsSim, AMDGPU
using AdaptiveOpticsSim.AlgorithmGraphs
import AdaptiveOpticsSim.Backends
AMDGPU.allowscalar(false)
const TARGET = Backends.AcceleratorComputeDevice(Backends.AMDGPUBackend(),1)

function integration_graph(input, target; T=eltype(input), output_shape=(3,5),
        output_sampling=(1.0,1.0), center=(2.0,3.25))
    node=pixel_integration_node(:pixels;input_shape=size(input),output_shape,
        input_sampling_m=(1.0,1.0),output_sampling_m=output_sampling,
        input_center_pixel=((size(input,1)+1)/2,(size(input,2)+1)/2),
        output_center_pixel=center,input_schema="test.pixel-integration.input/1",
        output_schema="test.pixel-integration.output/1",T)
    execution=target===Backends.HostComputeDevice() ? StreamGraphExecution() : CapturedGraphExecution()
    return prepare_algorithm_graph(algorithm_graph((node,);
        inputs=(graph_input(:flux,:pixels=>:input,input),),
        outputs=(graph_output(:flux,:pixels=>:output),));target,execution)
end

@testset "captured AMDGPU pixel integration retains geometry and current inputs" begin
    for T in (Float32,Float64), (shape,out,sampling,center) in (
            ((3,5),(3,5),(1.,1.),(2.,3.25)),
            ((4,6),(2,2),(2.,3.),(1.5,1.5)),
            ((1,1),(2,2),(.5,.5),(1.5,1.5)),
            ((1,1),(1,1),(1.,1.),(-5.,.5)),
            ((3,5),(8,9),(.6,.8),(4.1,5.2)))
        host=reshape(T.(1:prod(shape)),shape)
        device=AMDGPU.ROCArray(host)
        cpu=integration_graph(host,Backends.HostComputeDevice();output_shape=out,output_sampling=sampling,center)
        gpu=integration_graph(device,TARGET;output_shape=out,output_sampling=sampling,center)
        @test graph_execution_policy(gpu) isa CapturedGraphExecution
        @test captured_graph_node_count(gpu)==1
        for scale in (T(1),T(.5))
            host .*= scale
            copyto!(device,host);AMDGPU.synchronize()
            step_graph!(cpu);step_graph!(gpu)
            @test Array(graph_output(gpu,:flux)) ≈ graph_output(cpu,:flux) rtol=8eps(T) atol=eps(T)
        end
        reset_graph!(gpu)
        @test all(iszero,Array(graph_output(gpu,:flux)))
        step_graph!(gpu)
        @test Array(graph_output(gpu,:flux)) ≈ graph_output(cpu,:flux) rtol=8eps(T) atol=eps(T)
    end
end

@testset "captured AMDGPU pixel integration fails before publishing invalid input" begin
    host=ones(Float32,3,5)
    device=AMDGPU.ROCArray(host)
    # The invalid corner is outside this cropped output. CPU validation still
    # rejects it, so the captured input check must cover the entire source.
    cpu=integration_graph(host,Backends.HostComputeDevice();output_shape=(1,1),center=(1.,1.))
    gpu=integration_graph(device,TARGET;output_shape=(1,1),center=(1.,1.))
    for value in (-1f0,Float32(NaN),Float32(Inf))
        reset_graph!(cpu);reset_graph!(gpu)
        step_graph!(cpu);step_graph!(gpu)
        host[1,1]=value;copyto!(device,host);AMDGPU.synchronize()
        @test_throws AlgorithmGraphError step_graph!(cpu)
        ticket=step_graph_async!(gpu)
        @test_throws AlgorithmGraphError wait_graph_step!(ticket)
        @test gpu.state.step_sequence==1
        @test gpu.state.failed
        @test_throws AlgorithmGraphError step_graph!(gpu)
        host[1,1]=1f0;copyto!(device,host);AMDGPU.synchronize()
        reset_graph!(gpu)
        step_graph!(gpu)
        @test gpu.state.step_sequence==1
        @test only(Array(graph_output(gpu,:flux)))==1f0
    end
end
@test !any(id->id.name=="CUDA",keys(Base.loaded_modules))

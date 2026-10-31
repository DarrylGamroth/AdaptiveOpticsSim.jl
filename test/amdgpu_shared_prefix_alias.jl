# Explicit functional selector; no throughput or deadline claim.
# HSA_OVERRIDE_GFX_VERSION=<host setting> julia --project=<environment> test/amdgpu_shared_prefix_alias.jl
using Test, LinearAlgebra, AdaptiveOpticsSim, Proper, AMDGPU
using AdaptiveOpticsSim.AlgorithmGraphs
using AdaptiveOpticsSim.Backends

AMDGPU.allowscalar(false)
BLAS.set_num_threads(1)
Backends.set_fft_provider_threads!(1)

function alias_optical_graph(opd, amplitude, scale; execution=StreamGraphExecution())
    config = (;resolution=4,pupil_resolution=4,diameter_m=0.004,wavelength_um=1.55,
        steps=((op="lens",focal_length_m=0.1),(op="propagate",distance_m=0.02),
            (op="multiply",amplitude=scale),(op="propagate",distance_m=0.03)))
    node = proper_optical_train_node_factory(;T=Float32)(:optics,config,(;))
    return prepare_algorithm_graph(algorithm_graph((node,);
        inputs=(graph_input(:pupil_opd,:optics=>:pupil_opd,opd),
            graph_input(:pupil_amplitude,:optics=>:pupil_amplitude,amplitude)),
        outputs=(graph_output(:intensity,:optics=>:output),));
        target=compute_device(opd),execution)
end

function alias_external_pupil(a,b)
    return prepare_external_pupil((;a,b);
        exposure_durations=(a=ModelDuration(2),b=ModelDuration(2)),
        observations=(a=(output=:intensity,sampling_node=:optics),
            b=(output=:intensity,sampling_node=:optics)))
end

@testset "Native GPU pupil views retain independent staging" begin
    backing = AMDGPU.zeros(Float32,4,12)
    first_pupil, second_pupil = view(backing,:,1:4), view(backing,:,3:6)
    disjoint = view(backing,:,9:12)
    println("Julia=",VERSION," AMDGPU=",Base.pkgversion(AMDGPU))
    println("AMDGPU source=",pathof(AMDGPU))
    println("view types=",typeof(first_pupil)," / ",typeof(second_pupil))
    println("mightalias method=",which(Base.mightalias,(typeof(first_pupil),typeof(second_pupil))))
    println("dataids=",Base.dataids(first_pupil)," / ",Base.dataids(second_pupil))
    println("pointer delta bytes=",UInt(pointer(second_pupil))-UInt(pointer(first_pupil)),
        " first extent bytes=",sizeof(first_pupil))
    @test first_pupil isa ROCArray
    @test second_pupil isa ROCArray
    @test first_pupil !== second_pupil
    @test Base.mightalias(first_pupil,second_pupil)
    @test !Base.mightalias(first_pupil,disjoint)
    copyto!(first_pupil,ones(Float32,4,4))
    @test Array(second_pupil)[:,1:2] == ones(Float32,4,2)
    @test Array(second_pupil)[:,3:4] == zeros(Float32,4,2)
    fill!(backing,0)

    a = alias_optical_graph(first_pupil,AMDGPU.ones(Float32,4,4),1.0;
        execution=CapturedGraphExecution())
    b = alias_optical_graph(second_pupil,AMDGPU.ones(Float32,4,4),0.5;
        execution=CapturedGraphExecution())
    owner = alias_external_pupil(a,b)
    @test owner.propagation === nothing
    @test all(graph -> graph_execution_policy(graph) isa CapturedGraphExecution,(a,b))
    @test all(graph -> captured_graph_node_count(graph) == 1,(a,b))
    pupil = Float32[1e-7sin(0.3i+0.8j) for i in 1:4,j in 1:4]
    reference = alias_optical_graph(copy(pupil),ones(Float32,4,4),1.0)
    step_graph!(reference)
    @test step_external_pupil!((name,image,metadata)->nothing,owner,pupil,ones(Float32,4,4),
        ModelTimestamp(0),ModelDuration(2);domain=ntuple(_->UInt8(9),16),
        generation=UInt64(1),sequence=UInt64(1))
    @test Array(graph_output(a,:intensity)) ≈ graph_output(reference,:intensity) rtol=8e-4 atol=8e-6
    @test Array(first_pupil)[:,3:4] == pupil[:,1:2]
end

@testset "Exact shared native pupil bindings remain capture eligible" begin
    opd, amplitude = AMDGPU.zeros(Float32,4,4), AMDGPU.ones(Float32,4,4)
    a = alias_optical_graph(opd,amplitude,1.0;execution=CapturedGraphExecution())
    b = alias_optical_graph(opd,amplitude,0.5;execution=CapturedGraphExecution())
    owner = alias_external_pupil(a,b)
    @test owner.propagation !== nothing
    joint = AlgorithmGraphs.optical_execution_graph(owner,:a)
    @test joint === AlgorithmGraphs.optical_execution_graph(owner,:b)
    @test graph_execution_policy(joint) isa CapturedGraphExecution
    @test captured_graph_node_count(joint) == length(joint.nodes) == 1
    @test step_external_pupil!((name,image,metadata)->nothing,owner,zeros(Float32,4,4),ones(Float32,4,4),
        ModelTimestamp(0),ModelDuration(2);domain=ntuple(_->UInt8(10),16),
        generation=UInt64(1),sequence=UInt64(1))
    @test Array(graph_output(b,:intensity)) ≈ 0.25f0 .* Array(graph_output(a,:intensity)) rtol=8e-4 atol=8e-6
end

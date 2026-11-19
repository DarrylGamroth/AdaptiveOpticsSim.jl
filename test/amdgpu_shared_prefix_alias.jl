# Explicit functional selector; no throughput or deadline claim.
# HSA_OVERRIDE_GFX_VERSION=<host setting> julia --project=<environment> test/amdgpu_shared_prefix_alias.jl
using Test, LinearAlgebra, AdaptiveOpticsSim, Proper, AMDGPU
using AdaptiveOpticsSim.AlgorithmGraphs
using AdaptiveOpticsSim.Backends

AMDGPU.allowscalar(false)
BLAS.set_num_threads(1)
Backends.set_fft_provider_threads!(1)

function alias_optical_graph(opd, amplitude, scale; execution=StreamGraphExecution(), calibration=nothing)
    prefix = ((op="lens",focal_length_m=0.1),(op="propagate",distance_m=0.02))
    branch = (op="multiply",amplitude=scale)
    finish = (op="propagate",distance_m=0.03)
    config = (;resolution=4,pupil_resolution=4,diameter_m=0.004,wavelength_um=1.55,
        steps=(prefix...,branch,finish))
    if calibration !== nothing
        sampled = (op="sampled_mask",name="prefix_mask",amplitude="mask_amplitude",phase="mask_phase",
            sampling_m=(0.0009,0.0007),center_pixel=(2.8,3.2),outside_amplitude=0.31,
            outside_phase_rad=0.0,resampling="band_limited")
        config = merge(config,(;steps=(prefix...,sampled,branch,
            merge(sampled,(;name="suffix_mask")),finish),mask_amplitude="mask_amplitude",
            mask_phase="mask_phase"))
    end
    node = proper_optical_train_node_factory(;T=Float32)(:optics,config,(;),something(calibration,(;)))
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

@testset "Captured shared paths reuse sampled masks across replay and reset" begin
    calibration = (mask_amplitude=Float32[0.21+0.08i+0.02j for i in 1:4,j in 1:4],
        mask_phase=Float32[0.30sin(i+2j) for i in 1:4,j in 1:4])
    opd, amplitude = AMDGPU.zeros(Float32,4,4), AMDGPU.ones(Float32,4,4)
    paths = map((1.0,0.5)) do scale
        alias_optical_graph(opd,amplitude,scale;execution=CapturedGraphExecution(),calibration)
    end
    originals = map(path -> AlgorithmGraphs.prepared_graph_node(path,:optics),paths)
    masks = Tuple(train.workspace.masks[index] for train in originals for index in (3,5))
    coefficients = map(Array,masks)
    @test coefficients[1][2,3] != coefficients[1][3,2]
    plant = alias_external_pupil(paths...)
    joint = AlgorithmGraphs.optical_execution_graph(plant,:a)
    @test joint !== paths[1]
    @test joint === AlgorithmGraphs.optical_execution_graph(plant,:b)
    @test graph_execution_policy(joint) isa CapturedGraphExecution
    @test captured_graph_node_count(joint) == length(joint.nodes) == 1
    shared = AlgorithmGraphs.prepared_graph_node(joint,:propagation)
    @test shared.prefix.workspace.masks[3] === originals[1].workspace.masks[3]
    @test shared.geometry === originals[1].surfaces[4]
    for (branch, original) in zip(shared.branches,originals)
        @test branch.masks[2] === original.workspace.masks[5]
    end
    trains = (originals...,shared.prefix,map(branch -> branch.owner,shared.branches)...)
    for index in eachindex(trains), other in 1:index-1
        a, b = trains[index], trains[other]
        @test a.workspace.context !== b.workspace.context
        @test !Base.mightalias(a.workspace.wavefront.field,b.workspace.wavefront.field)
        @test !Base.mightalias(a.workspace.shifted_mask,b.workspace.shifted_mask)
        @test !Base.mightalias(a.workspace.padded_opd,b.workspace.padded_opd)
        @test !Base.mightalias(a.workspace.padded_amplitude,b.workspace.padded_amplitude)
        @test !Base.mightalias(a.output,b.output)
    end
    domain = ntuple(_->UInt8(11),16)
    frames = Tuple{Matrix{Float32},Matrix{Float32}}[]
    for (generation, sequences) in ((UInt64(1),1:2),(UInt64(2),1:1))
        generation == 2 && reset_external_pupil!(plant)
        for sequence in sequences
            pupil = Float32[sequence*1e-7sin(0.3i+0.8j) for i in 1:4,j in 1:4]
            transmission = Float32[0.7+0.1cos(sequence+i-j) for i in 1:4,j in 1:4]
            @test step_external_pupil!((name,image,metadata)->nothing,plant,pupil,transmission,
                ModelTimestamp(2(sequence-1)),ModelDuration(2);domain,generation,sequence=UInt64(sequence))
            @test graph_step_sequence(joint) == sequence
            frames_current = map(paths) do path
                captured = Array(graph_output(path,:intensity))
                step_graph!(path)
                @test captured ≈ Array(graph_output(path,:intensity)) rtol=8e-4 atol=8e-6
                captured
            end
            push!(frames,frames_current)
        end
    end
    @test frames[1] != frames[2]
    @test frames[1] == frames[3]
    @test map(Array,masks) == coefficients
end

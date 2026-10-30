# Explicit optional selector: julia --project=<source AMDGPU/Proper environment> test/amdgpu_proper_runtime_opd.jl
using Test, AdaptiveOpticsSim, Proper, AMDGPU, LinearAlgebra
using AdaptiveOpticsSim.AlgorithmGraphs
using AdaptiveOpticsSim.Backends
AMDGPU.allowscalar(false)
BLAS.set_num_threads(1)
Backends.set_fft_provider_threads!(1)

function sampled_opd_definition(config, opd, amplitude, screen)
    node=proper_optical_train_node_factory(T=eltype(opd))(:optics,config,(;))
    return algorithm_graph((node,);inputs=(graph_input(:opd,:optics=>:pupil_opd,opd),
        graph_input(:amplitude,:optics=>:pupil_amplitude,amplitude),
        graph_input(:screen,:optics=>:correction,screen)),
        outputs=(graph_output(:intensity,:optics=>:output),))
end

@testset "AMDGPU captured sampled OPD and entrance registration" begin
    for T in (Float32,Float64)
        screen_step=(op="add_opd",name="sampled correction",input="correction",size=(9,7),axes=(2,1),
            sampling_m=(0.0007,0.0004),center_pixel=(4.3,3.8),
            schema="test.sampled-opd-metres/1",outside_opd_m=-2e-8)
        config=(resolution=32,pupil_resolution=16,diameter_m=.008,wavelength_um=1.55,
            pupil_grid=(size=(16,16),axes=(2,1),sampling_m=(.0006,.0004),center_pixel=(8.5,8.5)),
            steps=(screen_step,(op="lens",focal_length_m=.12),(op="propagate",distance_m=.12)))
        opd=T[2e-8sin(.3i-.7j) for i in 1:16,j in 1:16]
        amplitude=T[.2+.6sin(.5i+.2j)^2 for i in 1:16,j in 1:16]
        screen=T[2e-7sin(.3i-.4j) for i in 1:9,j in 1:7]
        device=map(ROCArray,(opd,amplitude,screen))
        cpu=prepare_algorithm_graph(sampled_opd_definition(config,opd,amplitude,screen))
        gpu=prepare_algorithm_graph(sampled_opd_definition(config,device...);
            target=compute_device(first(device)),execution=CapturedGraphExecution())
        @test captured_graph_node_count(gpu)==1
        owner=AlgorithmGraphs.prepared_graph_node(gpu,:optics)
        maps,entrance=owner.interpolation,owner.entrance
        @test compute_device(entrance.row_indices)==compute_device(first(device))
        first_image=nothing
        rtol,atol=T===Float32 ? (8e-4,8e-6) : (8e-11,8e-12)
        for factor in (T(1),T(-1))
            screen .*= factor
            copyto!(device[3],screen);AMDGPU.synchronize()
            step_graph!(cpu);step_graph!(gpu)
            image=Array(graph_output(gpu,:intensity))
            @test image ≈ graph_output(cpu,:intensity) rtol=rtol atol=atol
            @test all(isfinite,image) && sum(image)>0
            @test owner.interpolation===maps && owner.entrance===entrance
            if first_image===nothing
                first_image=copy(image)
            else
                @test !isapprox(image,first_image;rtol=1e-3,atol=atol)
            end
        end
        reset_graph!(gpu)
        @test all(iszero,Array(graph_output(gpu,:intensity)))
        step_graph!(gpu)
        @test Array(graph_output(gpu,:intensity)) ≈ graph_output(cpu,:intensity) rtol=rtol atol=atol
    end
end
@test !any(id->id.name=="CUDA",keys(Base.loaded_modules))

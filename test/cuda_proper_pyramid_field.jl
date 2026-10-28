# Run explicitly with this checkout, Proper and CUDA in the selected environment.
# JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --startup-file=no --project=<environment> test/cuda_proper_pyramid_field.jl
using Test, AdaptiveOpticsSim, Proper, CUDA, LinearAlgebra
using AdaptiveOpticsSim.AlgorithmGraphs
using AdaptiveOpticsSim.Backends

CUDA.allowscalar(false)
BLAS.set_num_threads(1)
Backends.set_fft_provider_threads!(1)

function cuda_pyramid_field_definition(config, opd, amplitude, screens)
    T = eltype(opd)
    optics = proper_optical_train_node_factory(;T)(:optics,config,(;))
    pitch = optical_sampling_m(optics.config)
    p = config.pupil_resolution
    relay = ideal_pupil_relay_node(:relay; input_size=(config.resolution,config.resolution),
        output_size=(p,p), axes=(2,1), input_sampling_m=(pitch,pitch),
        output_sampling_m=(T(2)/p,T(2)/p), wavelength_m=T(750e-9),
        input_schema="test.proper-field/1", output_schema="test.receiving-pupil-field/1",T)
    pyramid = pyramid_rate_node(:pyramid; resolution=p,telescope_diameter_m=2,
        pupil_samples=p,source_wavelength_m=750e-9,source_photon_irradiance_m2_s=9,
        input_quantity=:electric_field, field_schema="test.receiving-pupil-field/1",
        photon_rate_schema="test.pyramid-rate/1",modulation=2,modulation_points=4,T)
    inputs = (graph_input(:opd,:optics=>:pupil_opd,opd),
        graph_input(:amplitude,:optics=>:pupil_amplitude,amplitude),
        map(keys(screens)) do role
            graph_input(role,:optics=>role,getproperty(screens,role))
        end...)
    return algorithm_graph((optics,relay,pyramid);name=:proper_pyramid_field,inputs,
        links=(link(:optics=>:output,:relay=>:input),
            link(:relay=>:output,:pyramid=>:electric_field)),
        outputs=(graph_output(:proper_field,:optics=>:output),
            graph_output(:pupil_field,:relay=>:output),
            graph_output(:rate,:pyramid=>:photon_rate)))
end

cuda_pyramid_field_host_allocations(graph) = @allocated step_graph!(graph)
cuda_pyramid_field_device_allocations(graph) = CUDA.@allocated step_graph!(graph)

function compare_cuda_pyramid_field(::Type{T}) where {T}
    p, n = 16, 32
    first_screen = (op="add_opd",name="upstream DM",input="upstream",size=(13,17),
        axes=(2,1),sampling_m=(0.0007,0.00055),center_pixel=(8.7,6.4),
        schema="test.dm-opd-metres/1",outside_opd_m=0)
    second_screen = (op="add_opd",name="downstream DM",input="downstream",size=(9,7),
        axes=(1,2),sampling_m=(0.0011,0.0013),center_pixel=(4.8,3.2),
        schema="test.dm-opd-metres/1",outside_opd_m=0)
    config = (resolution=n,pupil_resolution=p,diameter_m=0.008,wavelength_um=0.75,
        output_quantity=:electric_field,output_schema="test.proper-field/1",
        steps=(first_screen,(op="propagate",distance_m=0.8),second_screen,
            (op="propagate",distance_m=0.6)))
    initial_opd = T[3e-8*sin(0.4i-0.7j) for i in 1:p,j in 1:p]
    initial_amplitude = T[0.25+0.65sin(0.3i+0.6j)^2 for i in 1:p,j in 1:p]
    initial_screens = (upstream=T[1e-7*sin(0.3i+0.5j) for i in 1:13,j in 1:17],
        downstream=T[8e-8*cos(0.4i-0.6j) for i in 1:9,j in 1:7])
    host_opd, host_amplitude = copy(initial_opd), copy(initial_amplitude)
    host_screens = map(copy,initial_screens)
    opd, amplitude = CuArray(initial_opd), CuArray(initial_amplitude)
    screens = map(CuArray,initial_screens)
    host = prepare_algorithm_graph(cuda_pyramid_field_definition(config,host_opd,
        host_amplitude,host_screens);target=Backends.HostComputeDevice())
    captured = prepare_algorithm_graph(cuda_pyramid_field_definition(config,opd,
        amplitude,screens);target=compute_device(opd),execution=CapturedGraphExecution())
    @test captured_graph_node_count(captured) == 3
    proper_owner = AlgorithmGraphs.prepared_graph_node(captured,Val(:optics))
    relay_owner = AlgorithmGraphs.prepared_graph_node(captured,Val(:relay))
    pyramid_owner = AlgorithmGraphs.prepared_graph_node(captured,Val(:pyramid))
    for owner in (proper_owner,relay_owner,pyramid_owner)
        @test AlgorithmGraphs.graph_node_capture_capability(owner) isa AlgorithmGraphs.GraphNodeCaptureSafe
    end
    proper_field = graph_output(captured,Val(:proper_field))
    pupil_field = graph_output(captured,Val(:pupil_field))
    rate = graph_output(captured,Val(:rate))
    @test relay_owner.input === proper_field
    @test relay_owner.output === pupil_field
    @test pyramid_owner.field_input === pupil_field
    @test pyramid_owner.prepared.input.values !== pupil_field
    workspace = proper_owner.train.workspace
    retained = (workspace.wavefront.field,workspace.context.workspace.fft.scratch,
        proper_field,pupil_field,rate,pyramid_owner.prepared.input.values,
        pyramid_owner.prepared.workspace.modulation_batch.field_stack)
    addresses = map(pointer,retained)
    rtol, atol = T === Float32 ? (1e-3,1e-5) : (2e-10,2e-11)
    function compare_outputs!()
        CUDA.synchronize() # Publish caller writes before retained streams consume them.
        step_graph!(host)
        step_graph!(captured)
        for name in (:proper_field,:pupil_field,:rate)
            @test Array(graph_output(captured,name)) ≈ graph_output(host,name) rtol=rtol atol=atol
        end
        @test map(pointer,retained) == addresses
        @test workspace.wavefront.field === retained[1]
        @test workspace.context.workspace.fft.scratch === retained[2]
        return nothing
    end
    compare_outputs!()
    first_rate = Array(rate)
    full = Array(proper_field)
    pupil = Array(pupil_field)
    # Proper rows are y; the ideal relay explicitly maps them to output (x,y).
    @test pupil ≈ permutedims(full[9:24,9:24],(2,1)) rtol=rtol atol=atol
    @test sum(abs2,full) ≈ sum(abs2,initial_amplitude) rtol=rtol
    @test sum(abs2,pupil) < sum(abs2,full)
    scale_squared = T(9)*(T(2)/p)^2
    @test sum(first_rate) ≈ sum(abs2,pupil)*scale_squared rtol=rtol
    @test all(isfinite,first_rate) && sum(first_rate)>0

    # A change in one separated DM must reach the same captured sensor frame.
    host_screens.downstream .+= T[1.8e-7*(i-0.3j)/9 for i in 1:9,j in 1:7]
    copyto!(screens.downstream,host_screens.downstream)
    compare_outputs!()
    changed_rate = Array(rate)
    changed_field = Array(proper_field)
    @test !isapprox(changed_rate,first_rate;rtol=1e-3,atol=atol)

    # Piston rotates the complex field while preserving the measured rate.
    angle = T(0.73)
    host_opd .+= T(750e-9)*angle/T(2pi)
    copyto!(opd,host_opd)
    compare_outputs!()
    @test Array(proper_field) ≈ cis(angle).*changed_field rtol=rtol atol=atol
    @test Array(rate) ≈ changed_rate rtol=rtol atol=atol

    # The field crop and WFS must retain actual flux rather than renormalize it.
    host_amplitude .*= T(0.5)
    copyto!(amplitude,host_amplitude)
    compare_outputs!()
    @test Array(rate) ≈ T(0.25).*changed_rate rtol=rtol atol=atol
    @test sum(Array(rate)) ≈ sum(abs2,Array(pupil_field))*scale_squared rtol=rtol

    for _ in 1:32
        step_graph!(captured)
    end
    cuda_pyramid_field_host_allocations(captured)
    cuda_pyramid_field_device_allocations(captured)
    @test cuda_pyramid_field_host_allocations(captured) == 0
    @test cuda_pyramid_field_device_allocations(captured) == 0

    rate_before_reset = Array(rate)
    host_opd .= initial_opd
    host_amplitude .= initial_amplitude
    copyto!(opd,initial_opd)
    copyto!(amplitude,initial_amplitude)
    for role in keys(screens)
        getproperty(host_screens,role) .= getproperty(initial_screens,role)
        copyto!(getproperty(screens,role),getproperty(initial_screens,role))
    end
    CUDA.synchronize()
    reset_graph!(host)
    reset_graph!(captured)
    for name in (:proper_field,:pupil_field)
        @test all(iszero,Array(graph_output(captured,name)))
    end
    # The stateless Pyramid adapter retains its last product until the next step.
    @test Array(rate) == rate_before_reset
    @test map(pointer,retained) == addresses
    compare_outputs!()
    @test Array(rate) ≈ first_rate rtol=rtol atol=atol
    return nothing
end

@testset "CUDA Proper field to analytic Pyramid capture" begin
    for T in (Float32,Float64)
        compare_cuda_pyramid_field(T)
    end
end

@testset "CUDA rectangular ideal pupil relay capture" begin
    for T in (Float32,Float64)
        source = Complex{T}[complex(10i+j,i-2j) for i in 1:8,j in 1:10]
        input = CuArray(source)
        node = ideal_pupil_relay_node(:relay;input_size=(8,10),output_size=(5,6),
            axes=(2,1),input_sampling_m=(0.002,0.003),
            output_sampling_m=(0.006,0.004),wavelength_m=750e-9,
            input_schema="test.relative-field/1",output_schema="test.pupil-field/1",T)
        graph = prepare_algorithm_graph(algorithm_graph((node,);
            inputs=(graph_input(:input,:relay=>:input,input),),
            outputs=(graph_output(:output,:relay=>:output),));
            target=compute_device(input),execution=CapturedGraphExecution())
        @test captured_graph_node_count(graph) == 1
        step_graph!(graph)
        @test Array(graph_output(graph,:output)) == permutedims(source[2:7,4:8],(2,1))
    end
end

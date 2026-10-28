# Run explicitly with this checkout, Proper and CUDA in the selected environment.
# JULIA_NUM_THREADS=1 julia --startup-file=no --project=<environment> test/cuda_proper_runtime_opd.jl
using Test, AdaptiveOpticsSim, Proper, CUDA, LinearAlgebra
using AdaptiveOpticsSim.AlgorithmGraphs
using AdaptiveOpticsSim.Backends

CUDA.allowscalar(false)
BLAS.set_num_threads(1)
Backends.set_fft_provider_threads!(1)
const CUDA_OPD_EXT = Base.get_extension(AdaptiveOpticsSim, :AdaptiveOpticsSimProperExt)

function cuda_opd_definition(config, opd, amplitude, screens; bindings=(;))
    node = proper_optical_train_node_factory(T=eltype(opd))(:optics, config, (;), bindings)
    inputs = (graph_input(:opd, :optics=>:pupil_opd, opd),
        graph_input(:amplitude, :optics=>:pupil_amplitude, amplitude),
        map(keys(screens)) do role
            graph_input(role, :optics=>role, getproperty(screens,role))
        end...)
    return algorithm_graph((node,);name=:cuda_runtime_opd,inputs,
        outputs=(graph_output(:intensity,:optics=>:output),))
end

cuda_opd_host_allocations(graph) = @allocated step_graph!(graph)
cuda_opd_device_allocations(graph) = CUDA.@allocated step_graph!(graph)

function compare_cuda_opd_capture(::Type{T}, config, initial_screens; bindings=(;)) where {T}
    p = config.pupil_resolution
    initial_opd = T[2e-8sin(0.3i-0.7j) for i in 1:p,j in 1:p]
    initial_amplitude = T[0.2+0.6sin(0.5i+0.2j)^2 for i in 1:p,j in 1:p]
    host_opd, host_amplitude = copy(initial_opd), copy(initial_amplitude)
    host_screens = map(copy, initial_screens)
    opd, amplitude = CuArray(initial_opd), CuArray(initial_amplitude)
    screens = map(CuArray, initial_screens)
    target = compute_device(opd)
    definition = cuda_opd_definition(config,opd,amplitude,screens;bindings)
    streamed = prepare_algorithm_graph(definition;target,execution=StreamGraphExecution())
    captured = prepare_algorithm_graph(definition;target,execution=CapturedGraphExecution())
    host = prepare_algorithm_graph(cuda_opd_definition(config,host_opd,host_amplitude,
        host_screens;bindings);target=Backends.HostComputeDevice())
    owner = AlgorithmGraphs.prepared_graph_node(captured,:optics)
    @test captured_graph_node_count(captured) == 1
    @test AlgorithmGraphs.graph_node_capture_capability(owner) isa AlgorithmGraphs.GraphNodeCaptureSafe
    workspace = owner.train.workspace
    field, scratch, maps, entrance = workspace.wavefront.field,
        workspace.context.workspace.fft.scratch, owner.interpolation, owner.entrance
    padded_opd, padded_amplitude = workspace.padded_opd, workspace.padded_amplitude
    fft_plans = (workspace.context.workspace.fft.forward_plan,
        workspace.context.workspace.fft.backward_plan)
    for plan in (maps...,entrance)
        plan === nothing && continue
        @test compute_device(plan.row_indices) == target
        @test compute_device(plan.column_weights) == target
    end
    rtol, atol = T === Float32 ? (8e-4,8e-6) : (8e-11,8e-12)
    first_image = zeros(T,config.resolution,config.resolution)
    for frame in 1:4
        host_opd .= initial_opd .+ T(frame-1) .* T[5e-8*(i-0.3j)/p for i in 1:p,j in 1:p]
        host_amplitude .= initial_amplitude .* T(1-0.13*(frame-1))
        for role in keys(host_screens)
            values, initial = getproperty(host_screens,role), getproperty(initial_screens,role)
            values .= initial .* T(1+0.23*(frame-1))
            copyto!(getproperty(screens,role),values)
        end
        copyto!(opd,host_opd); copyto!(amplitude,host_amplitude)
        CUDA.synchronize()
        step_graph!(host); step_graph!(streamed); step_graph!(captured)
        actual = Array(graph_output(captured,:intensity))
        @test actual ≈ Array(graph_output(streamed,:intensity)) rtol=rtol atol=atol
        @test actual ≈ graph_output(host,:intensity) rtol=rtol atol=atol
        @test all(isfinite,actual) && sum(actual)>0
        @test workspace.wavefront.field === field
        @test workspace.context.workspace.fft.scratch === scratch
        @test workspace.context.workspace.fft.forward_plan === fft_plans[1]
        @test workspace.context.workspace.fft.backward_plan === fft_plans[2]
        @test workspace.padded_opd === padded_opd && workspace.padded_amplitude === padded_amplitude
        @test owner.interpolation === maps && owner.entrance === entrance
        if frame == 1
            first_image .= actual
        else
            @test !isapprox(actual,first_image;rtol=1e-3,atol=atol)
        end
    end
    # A screen-only mutation must reach the retained executable, independently
    # of changes to entrance amplitude or OPD.
    if !isempty(keys(host_screens))
        before = copy(graph_output(host,:intensity))
        role = first(keys(host_screens))
        values = getproperty(host_screens,role)
        values .*= -1
        copyto!(getproperty(screens,role),values)
        CUDA.synchronize()
        step_graph!(host); step_graph!(captured)
        @test Array(graph_output(captured,:intensity)) ≈ graph_output(host,:intensity) rtol=rtol atol=atol
        @test !isapprox(graph_output(host,:intensity),before;rtol=1e-3,atol=atol)
    end
    for _ in 1:32
        step_graph!(captured)
    end
    cuda_opd_host_allocations(captured)
    @test cuda_opd_host_allocations(captured) == 0
    @test cuda_opd_device_allocations(captured) == 0
    @test cuda_opd_device_allocations(streamed) == 0
    copyto!(opd,initial_opd); copyto!(amplitude,initial_amplitude)
    for role in keys(screens)
        copyto!(getproperty(screens,role),getproperty(initial_screens,role))
    end
    CUDA.synchronize()
    reset_graph!(streamed); reset_graph!(captured)
    @test all(iszero,Array(graph_output(captured,:intensity)))
    @test workspace.wavefront.field === field
    step_graph!(streamed); step_graph!(captured)
    @test Array(graph_output(captured,:intensity)) ≈ first_image rtol=rtol atol=atol
    @test Array(graph_output(captured,:intensity)) ≈ Array(graph_output(streamed,:intensity)) rtol=rtol atol=atol
    return nothing
end

# Independent tent-basis interpolation of real samples, including outside OPD.
function cuda_opd_grid_reference(input, grid, n, pitch, outside)
    source = permutedims(input,grid.axes)
    return [begin
        y = (row-fld(n,2)-1)*pitch/grid.sampling_m[1]+grid.center_pixel[1]
        x = (column-fld(n,2)-1)*pitch/grid.sampling_m[2]+grid.center_pixel[2]
        outside+sum(max(0,1-abs(y-i))*max(0,1-abs(x-j))*(source[i,j]-outside)
            for i in axes(source,1),j in axes(source,2))
    end for row in 1:n,column in 1:n]
end

@testset "CUDA real sampled OPD interpolation" begin
    for T in (Float32,Float64), permutation in ((1,2),(2,1))
        input = T[2e-7*(i^2-2j) for i in 1:3,j in 1:5]
        grid = (size=size(input),axes=permutation,sampling_m=(T(0.0013),T(0.0007)),
            center_pixel=(T(2.2),T(3.4)))
        device_input = CuArray(input)
        plan = CUDA_OPD_EXT._bilinear_interpolation_plan(grid,16,T(0.001),T(3e-8))
        device_plan = CUDA_OPD_EXT._prepare_grid_interpolation(plan,compute_device(device_input))
        output = CUDA.zeros(T,16,16)
        CUDA_OPD_EXT._interpolate_grid!(output,device_input,device_plan)
        @test Array(output) ≈ cuda_opd_grid_reference(input,grid,16,T(0.001),T(3e-8)) rtol=20eps(T)
    end
end

@testset "CUDA runtime OPD and entrance-grid capture" begin
    for T in (Float32,Float64)
        base = (resolution=32,pupil_resolution=16,diameter_m=0.008,wavelength_um=1.55)
        pupil_grid = (size=(16,16),axes=(2,1),sampling_m=(0.0006,0.0004),center_pixel=(8.5,8.5))
        first_step = (op="add_opd",name="upstream corrector",input="upstream",size=(13,17),
            axes=(2,1),sampling_m=(0.0007,0.00035),center_pixel=(8.4,6.7),
            schema="test.opd-metres/1",outside_opd_m=0.0)
        second_step = (op="add_opd",name="downstream corrector",input="downstream",size=(9,7),
            axes=(1,2),sampling_m=(0.000015,0.00001),center_pixel=(4.4,3.8),
            schema="test.opd-metres/1",outside_opd_m=-3e-8)
        mask = (op="sampled_mask",name="sampled focal mask",amplitude="transmission",
            sampling_m=(0.000012,0.00001),center_pixel=(3.2,4.4),
            outside_amplitude=0.8,outside_phase_rad=0.0)
        propagation = ((op="lens",focal_length_m=0.12),(op="propagate",distance_m=0.12),
            mask,(op="propagate",distance_m=0.03))
        calibration = T[0.5+0.4sin(0.2i+0.4j)^2 for i in 1:5,j in 1:7]
        initial = (upstream=T[2e-7sin(0.3i-0.4j) for i in 1:13,j in 1:17],
            downstream=T[1e-7cos(0.5i+0.2j) for i in 1:9,j in 1:7])
        steps = (first_step,propagation[1:3]...,second_step,propagation[4])
        compare_cuda_opd_capture(T,merge(base,(;steps,transmission="calibration")),initial;
            bindings=(;calibration))
        compare_cuda_opd_capture(T,merge(base,(;steps,pupil_grid,transmission="calibration")),initial;
            bindings=(;calibration))
        compare_cuda_opd_capture(T,merge(base,(;steps=propagation,pupil_grid,transmission="calibration")),(;);
            bindings=(;calibration))
    end
end

using AdaptiveOpticsSim, Proper, Test
using AdaptiveOpticsSim.AlgorithmGraphs

const RUNTIME_OPD_BASE = (resolution=16, pupil_resolution=16,
    diameter_m=0.016, wavelength_um=1.55)

runtime_opd_step(input; size=(16,16), axes=(1,2), sampling_m=(0.001,0.001),
    center_pixel=(9.0,9.0), outside_opd_m=0.0) =
    (; op="add_opd", name=string(input), input=string(input), size, axes,
        sampling_m, center_pixel, schema="test.opd-metres/1", outside_opd_m)

function runtime_opd_graph(steps, screens; T=Float64, base=RUNTIME_OPD_BASE,
        opd=zeros(T,base.pupil_resolution,base.pupil_resolution),
        amplitude=ones(T,size(opd)))
    node = proper_optical_train_node_factory(; T)(:optics,merge(base,(;steps)),(;))
    inputs = (graph_input(:opd,:optics=>:pupil_opd,opd),
        graph_input(:amplitude,:optics=>:pupil_amplitude,amplitude),
        map(keys(screens)) do role
            graph_input(role,:optics=>role,getproperty(screens,role))
        end...)
    return prepare_algorithm_graph(algorithm_graph((node,);name=:runtime_opd,
        inputs,outputs=(graph_output(:intensity,:optics=>:output),));
        target=AdaptiveOpticsSim.Backends.HostComputeDevice())
end

# Independent tensor-product tent reference, including each outside neighbor.
function runtime_opd_reference(input, step, n, pitch)
    source = permutedims(input, step.axes)
    outside = step.outside_opd_m
    return [begin
        y = (row-fld(n,2)-1)*pitch/step.sampling_m[1]+step.center_pixel[1]
        x = (column-fld(n,2)-1)*pitch/step.sampling_m[2]+step.center_pixel[2]
        outside + sum(max(0,1-abs(y-i))*max(0,1-abs(x-j))*(source[i,j]-outside)
            for i in axes(source,1), j in axes(source,2))
    end for row in 1:n, column in 1:n]
end

function runtime_opd_direct(steps, screens; base=RUNTIME_OPD_BASE,
        amplitude=ones(base.resolution,base.resolution),opd=zeros(size(amplitude)))
    context = Proper.RunContext(Matrix{ComplexF64})
    wf = prop_begin(base.diameter_m,base.wavelength_um*1e-6,base.resolution;
        beam_diam_fraction=base.pupil_resolution/base.resolution,context)
    prop_multiply(wf,amplitude)
    prop_add_phase(wf,opd)
    for step in steps
        if step.op == "add_opd"
            opd = runtime_opd_reference(getproperty(screens,Symbol(step.input)),step,
                base.resolution,prop_get_sampling(wf))
            prop_add_phase(wf,opd)
        elseif step.op == "propagate"
            prop_propagate(wf,step.distance_m,context)
        elseif step.op == "lens"
            prop_lens(wf,step.focal_length_m,context)
        elseif step.op == "circular_aperture"
            prop_circular_aperture(wf,step.radius_m)
        elseif step.op == "multiply"
            prop_multiply(wf,step.amplitude)
        end
    end
    return first(prop_end(wf))
end

function runtime_opd_allocations(graph)
    step_graph!(graph)
    step_graph!(graph)
    return @allocated step_graph!(graph)
end

@testset "Long runtime OPD train" begin
    segment = ((op="lens",focal_length_m=0.8),(op="propagate",distance_m=0.02),
        (op="circular_aperture",radius_m=0.006),(op="multiply",amplitude=0.999))
    optics = Tuple(step for _ in 1:9 for step in segment)
    steps = (runtime_opd_step(:first),runtime_opd_step(:second),optics...,
        runtime_opd_step(:third),runtime_opd_step(:fourth))
    @test length(steps) == 40
    screens = (first=[1e-7sin(0.3i-0.2j) for i in 1:16,j in 1:16],
        second=[2e-8cos(0.2i+0.5j) for i in 1:16,j in 1:16],
        third=fill(1e-8,16,16),fourth=zeros(16,16))
    reference = runtime_opd_direct(steps,screens)
    for T in (Float32,Float64)
        graph = runtime_opd_graph(steps,map(x->T.(x),screens);T)
        step_graph!(graph)
        @test graph_output(graph,:intensity) ≈ reference
        @test runtime_opd_allocations(graph) == 0
    end
end

@testset "Sampled entrance pupil grid" begin
    grid = (size=(8,8),axes=(2,1),sampling_m=(0.0011,0.0009),center_pixel=(4.5,4.5))
    base = (resolution=16,pupil_resolution=8,diameter_m=0.008,wavelength_um=1.55)
    sampled = merge(base,(pupil_grid=grid,))
    amplitude = [0.2+0.5i/8+0.3j/8 for i in 1:8,j in 1:8]
    opd = [1.5e-7sin(0.8i-0.3j) for i in 1:8,j in 1:8]
    geometry = merge(grid,(outside_opd_m=0.0,))
    mapped_amplitude = runtime_opd_reference(amplitude,geometry,16,0.001)
    mapped_opd = runtime_opd_reference(opd,geometry,16,0.001)
    propagation = ((op="lens",focal_length_m=0.1),(op="propagate",distance_m=0.1))
    reference = runtime_opd_direct(propagation,(;);base,amplitude=mapped_amplitude,opd=mapped_opd)
    for T in (Float32,Float64)
        pupil = runtime_opd_graph((),(;);T,base=sampled,amplitude=T.(amplitude),opd=T.(opd))
        step_graph!(pupil)
        @test graph_output(pupil,:intensity) ≈ abs2.(mapped_amplitude)
        @test runtime_opd_allocations(pupil) == 0

        propagated = runtime_opd_graph(propagation,(;);T,base=sampled,
            amplitude=T.(amplitude),opd=T.(opd))
        step_graph!(propagated)
        @test graph_output(propagated,:intensity) ≈ reference
        @test runtime_opd_allocations(propagated) == 0
        @test @inferred(step_graph!(propagated)) === propagated

        # Entrance and screen registration agree, so their OPD maps cancel.
        screen = runtime_opd_step(:correction;size=(8,8),axes=grid.axes,
            sampling_m=grid.sampling_m,center_pixel=grid.center_pixel)
        corrected = runtime_opd_graph((screen,propagation...),(correction=T.(-opd),);
            T,base=sampled,amplitude=T.(amplitude),opd=T.(opd))
        step_graph!(corrected)
        @test graph_output(corrected,:intensity) ≈ runtime_opd_direct(propagation,(;);
            base,amplitude=mapped_amplitude)
        @test runtime_opd_allocations(corrected) == 0
    end
    # The legacy padding convention reproduces the registration mismatch when
    # a physical source grid is omitted; the opt-in path changes that result.
    legacy = runtime_opd_graph(propagation,(;);base,amplitude,opd)
    step_graph!(legacy)
    @test !isapprox(graph_output(legacy,:intensity),reference;rtol=1e-3)

    factory = proper_optical_train_node_factory()
    for invalid in (0,merge(grid,(size=(7,8),)),merge(grid,(axes=(2,2),)),
            merge(grid,(sampling_m=(0.0,0.001),)),
            merge(grid,(center_pixel=(Inf,4.5),)),
            merge(grid,(outside_opd_m=0.0,)))
        @test_throws AlgorithmGraphError factory(:bad,merge(base,(pupil_grid=invalid,steps=())),(;))
    end
    @test_throws AlgorithmGraphError factory(:bad,
        merge(sampled,(steps=(),output_sampling_m=0.001)),(;))
    @test_throws AlgorithmGraphError factory(:bad,
        merge(sampled,(steps=(),channels=((wavelength_um=1.55,weight=1.0,steps=()),))), (;))
    node = factory(:optics,merge(sampled,(steps=(),)),(;))
    device = AdaptiveOpticsSim.Backends.AcceleratorComputeDevice(
        AdaptiveOpticsSim.Backends.CUDABackend(),0)
    @test_throws AlgorithmGraphError AlgorithmGraphs.prepare_graph_node(
        typeof(node).parameters[2],node.config,(;),
        (pupil_opd=zeros(8,8),pupil_amplitude=ones(8,8)),
        (output=zeros(16,16),),(;),device)
end

@testset "Runtime sampled OPD steps" begin
    a = [2e-7sin(0.3i+0.2j) for i in 1:16,j in 1:16]
    b = [1e-7cos(0.5i-0.3j) for i in 1:16,j in 1:16]
    amplitude = [exp(-((i-9)^2+(j-9)^2)/30) for i in 1:16,j in 1:16]
    propagation = ((op="lens",focal_length_m=0.1),(op="propagate",distance_m=0.1))
    first_step, second_step = runtime_opd_step(:first), runtime_opd_step(:second)
    for T in (Float32,Float64)
        empty = runtime_opd_graph(propagation,(;);T,amplitude=T.(amplitude))
        zero_screen = runtime_opd_graph((first_step,propagation...),
            (first=zeros(T,16,16),);T,amplitude=T.(amplitude))
        step_graph!(empty); step_graph!(zero_screen)
        @test graph_output(zero_screen,:intensity) ≈ graph_output(empty,:intensity)

        addition = runtime_opd_graph((first_step,second_step,propagation...),
            (first=T.(a),second=T.(b));T,amplitude=T.(amplitude))
        entrance = runtime_opd_graph(propagation,(;);T,opd=T.(a.+b),amplitude=T.(amplitude))
        step_graph!(addition); step_graph!(entrance)
        @test graph_output(addition,:intensity) ≈ graph_output(entrance,:intensity)
        cancellation = runtime_opd_graph((first_step,second_step,propagation...),
            (first=T.(a),second=T.(-a));T,amplitude=T.(amplitude))
        step_graph!(cancellation)
        @test graph_output(cancellation,:intensity) ≈ graph_output(empty,:intensity)

        piston = runtime_opd_graph((first_step,propagation...),
            (first=fill(T(2e-7),16,16),);T,amplitude=T.(amplitude))
        step_graph!(piston)
        @test graph_output(piston,:intensity) ≈ graph_output(empty,:intensity)
        @test @inferred(step_graph!(addition)) === addition
        @test runtime_opd_allocations(addition) == 0
    end

    # Distinct planes retain diffraction between OPD applications. The second
    # plane follows a focal propagation, so its mapping uses the changed pitch.
    rectangular = [1e-7*(i^2-2j) for i in 1:5,j in 1:3]
    transposed = runtime_opd_step(:second;size=(5,3),axes=(2,1),
        sampling_m=(0.000012,0.000008),center_pixel=(1.7,2.8),outside_opd_m=3e-8)
    steps = (first_step,propagation...,transposed,(op="propagate",distance_m=0.03))
    graph = runtime_opd_graph(steps,(first=a,second=rectangular);amplitude)
    step_graph!(graph)
    @test graph_output(graph,:intensity) ≈ runtime_opd_direct(steps,
        (first=a,second=rectangular);amplitude)
    colocated = runtime_opd_graph((first_step,transposed,propagation...,
        (op="propagate",distance_m=0.03)),(first=a,second=rectangular);amplitude)
    step_graph!(colocated)
    @test !isapprox(graph_output(graph,:intensity),graph_output(colocated,:intensity);rtol=1e-3)
    @test runtime_opd_allocations(graph) == 0

    # Equivalent explicit axis declarations agree, while swapping an asymmetric
    # rectangle without the matching permutation yields a different wavefront.
    proper_axes = merge(transposed,(size=(3,5),axes=(1,2)))
    same = runtime_opd_graph((proper_axes,propagation...),
        (second=permutedims(rectangular),);amplitude)
    swapped = runtime_opd_graph((transposed,propagation...),(second=rectangular,);amplitude)
    wrong_axes = runtime_opd_graph((merge(transposed,(axes=(1,2),)),propagation...),
        (second=rectangular,);amplitude)
    step_graph!(same); step_graph!(swapped); step_graph!(wrong_axes)
    @test graph_output(same,:intensity) ≈ graph_output(swapped,:intensity)
    @test !isapprox(graph_output(same,:intensity),graph_output(wrong_axes,:intensity);rtol=1e-3)

    # Input values remain live; reset clears products/scratch, preserving inputs.
    live = copy(a)
    changing = runtime_opd_graph((first_step,propagation...),(first=live,);amplitude)
    step_graph!(changing)
    before = copy(graph_output(changing,:intensity))
    fill!(live,0)
    reset_graph!(changing)
    @test iszero(graph_output(changing,:intensity))
    step_graph!(changing)
    @test !isapprox(graph_output(changing,:intensity),before;rtol=1e-3)
    @test graph_output(changing,:intensity) ≈ runtime_opd_direct(propagation,(;);amplitude)

    # Native even AOS pupils center between samples. Interpolating OPD across
    # the half-pixel origin preserves unit amplitude, including phase wraps.
    half_center = runtime_opd_step(:first;center_pixel=(8.5,8.5))
    half = runtime_opd_graph((half_center,propagation...),(first=a,);amplitude)
    step_graph!(half)
    @test graph_output(half,:intensity) ≈ runtime_opd_direct((half_center,propagation...),
        (first=a,);amplitude)
    wrapping = runtime_opd_step(:first;size=(1,2),center_pixel=(1.0,1.5))
    phase_only = runtime_opd_graph((wrapping,),(first=[0.7e-6 -0.7e-6],))
    step_graph!(phase_only)
    @test graph_output(phase_only,:intensity) ≈ ones(16,16)

    factory = proper_optical_train_node_factory()
    for invalid in (merge(first_step,(input="pupil_opd",)),
            merge(first_step,(size=(0,16),)),merge(first_step,(size=(16,),)),
            merge(first_step,(axes=(1,1),)),merge(first_step,(axes=(1,2,3),)),
            merge(first_step,(sampling_m=(0.0,0.001),)),
            merge(first_step,(center_pixel=(NaN,9.0),)),
            merge(first_step,(outside_opd_m=Inf,)),merge(first_step,(schema="",)))
        @test_throws AlgorithmGraphError factory(:bad,merge(RUNTIME_OPD_BASE,(steps=(invalid,),)),(;))
    end
    @test_throws AlgorithmGraphError factory(:bad,
        merge(RUNTIME_OPD_BASE,(steps=(first_step,first_step),)),(;))
    @test_throws AlgorithmGraphError factory(:bad,
        merge(RUNTIME_OPD_BASE,(steps=ntuple(i->runtime_opd_step(Symbol("s",i)),65),)),(;))
    @test_throws AlgorithmGraphError factory(:bad,
        merge(RUNTIME_OPD_BASE,(steps=(first_step,
            ntuple(_->(op="multiply",amplitude=1.0),256)...),)),(;))
    @test_throws AlgorithmGraphError factory(:bad,
        merge(RUNTIME_OPD_BASE,(steps=ntuple(_->(op="multiply",amplitude=1.0),257),
            pupil_grid=(size=(16,16),axes=(1,2),sampling_m=(0.001,0.001),center_pixel=(9.0,9.0)))),(;))
    @test_throws AlgorithmGraphError factory(:bad,
        merge(RUNTIME_OPD_BASE,(steps=(first_step,),output_sampling_m=0.001)),(;))
    @test_throws AlgorithmGraphError factory(:bad,
        merge(RUNTIME_OPD_BASE,(steps=(),channels=((wavelength_um=1.55,weight=1.0,
            steps=(first_step,)),))), (;))
    @test_throws AlgorithmGraphError runtime_opd_graph((first_step,),(first=zeros(3,4),))
    @test_throws AlgorithmGraphError runtime_opd_graph((first_step,),(first=zeros(Float32,16,16),))
    node = factory(:optics,merge(RUNTIME_OPD_BASE,(steps=(first_step,),)),(;))
    device = AdaptiveOpticsSim.Backends.AcceleratorComputeDevice(
        AdaptiveOpticsSim.Backends.CUDABackend(),0)
    @test_throws AlgorithmGraphError AlgorithmGraphs.prepare_graph_node(
        typeof(node).parameters[2],node.config,(;),
        (pupil_opd=zeros(16,16),pupil_amplitude=ones(16,16),first=zeros(16,16)),
        (output=zeros(16,16),),(;),device)
end

# Run in an environment with this checkout, Proper, FilterGraphAlgorithms and FilterGraphs.
# Select cpu (default), cuda, or amdgpu as the first argument.
using AdaptiveOpticsSim, Proper, FilterGraphs, FilterGraphAlgorithms, Test
using AdaptiveOpticsSim.AlgorithmGraphs
const B = AdaptiveOpticsSim.Backends
const O = AdaptiveOpticsSim.Optics

function optical_train(target, policy, wavelength, radius)
    steps = ((; op="circular_aperture", name="aperture", radius_m=radius),)
    config = (; resolution=8, pupil_resolution=8, diameter_m=0.004,
        wavelength_um=wavelength, steps)
    node = proper_optical_train_node_factory(T=Float32)(:optics, config, (;))
    buffer(value) = fill!(B.allocate_device_array(target, Float32, 8, 8), value)
    definition = algorithm_graph((node,);
        inputs=(graph_input(:pupil_opd, :optics=>:pupil_opd, buffer(0)),
            graph_input(:pupil_amplitude, :optics=>:pupil_amplitude, buffer(1))),
        outputs=(graph_output(:intensity, :optics=>:output),))
    return prepare_algorithm_graph(definition; target, execution=policy)
end

function exposure_owner(target, policy)
    llowfs = optical_train(target, policy, 1.4, 0.0016)
    scc = optical_train(target, policy, 1.55, 0.0019)
    unfringed = optical_train(target, policy, 1.55, 0.0012)
    return prepare_external_pupil((; llowfs, scc);
        exposure_durations=(llowfs=ModelDuration(2), scc=ModelDuration(4)),
        observations=(llowfs=(output=:intensity, sampling_node=:optics),
            scc=(output=:intensity, sampling_node=:optics,
                chopping=(; unfringed, initial_phase=:fringed))))
end

function collect_exposures(owner)
    events = NamedTuple[]
    emit(name, image, metadata) = push!(events, (; name, image=Array(image), metadata))
    domain = ntuple(_ -> UInt8(1), 16)
    opd = zeros(Float32, 8, 8)
    amplitude = ones(Float32, 8, 8)
    feed(time, duration, seq; generation=1) = step_external_pupil!(
        emit, owner, opd, amplitude, ModelTimestamp(time), ModelDuration(duration);
        domain, generation=UInt64(generation), sequence=UInt64(seq))
    for (seq, (time, duration)) in enumerate(((0,1), (1,2), (3,2), (5,2), (7,2)))
        amplitude .*= 0.99f0
        @test feed(time, duration, seq)
    end
    @test !feed(7, 2, 5)
    @test feed(0, 2, 1; generation=2)
    @test feed(2, 2, 2; generation=2)
    reset_external_pupil!(owner)
    @test feed(0, 2, 1; generation=3)
    return events
end

selection = get(ARGS, 1, "cpu")
if selection == "cuda"
    using CUDA
    @test CUDA.functional()
    CUDA.allowscalar(false)
    target = B.AcceleratorComputeDevice(B.CUDABackend(), 0)
elseif selection == "amdgpu"
    using AMDGPU
    @test AMDGPU.functional()
    AMDGPU.allowscalar(false)
    target = B.AcceleratorComputeDevice(B.AMDGPUBackend(), AMDGPU.device_id(AMDGPU.device()))
else
    @assert selection == "cpu"
    target = B.HostComputeDevice()
end
policy = selection == "cpu" ? StreamGraphExecution() : CapturedGraphExecution()

@testset "FilterGraphs independent exposures and SCC chopping" begin
    oracle = exposure_owner(target, policy)
    prepared = exposure_owner(target, policy)
    pilot = O.prepare_filtergraphs_external_pupil(prepared)
    if selection == "cpu"
        @test AlgorithmGraphs.optical_execution_graph(pilot, :scc) isa FilterGraphs.PreparedGraph
    else
        @test all(g -> g isa FilterGraphs.PreparedCapturedGraph, (
            AlgorithmGraphs.optical_execution_graph(pilot, :scc),
            AlgorithmGraphs.optical_execution_graph(pilot, :scc; unfringed=true),
            pilot.propagation.crossing.graph,
        ))
    end
    @test_throws ArgumentError O.prepare_filtergraphs_external_pupil(pilot)
    actual = collect_exposures(pilot)
    expected = collect_exposures(oracle)
    @test length(actual) == length(expected) == 10
    for (a, e) in zip(actual, expected)
        @test (a.name, a.metadata) == (e.name, e.metadata)
        @test a.image ≈ e.image rtol=2e-6
    end
    llowfs = filter(x -> x.name === :llowfs, actual)
    scc = filter(x -> x.name === :scc, actual)
    @test [x.metadata.pts for x in llowfs] == [0,2,4,6,0,2,0]
    @test [x.metadata.pts for x in scc] == [0,4,0]
    @test [x.metadata.chopped_image_kind for x in scc] == [:fringed,:unfringed,:fringed]
    @test all(x -> !haskey(x.metadata, :chopped_image_kind), llowfs)
end

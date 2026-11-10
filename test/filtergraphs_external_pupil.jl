# Run in an environment with this checkout, Proper, FilterGraphAlgorithms and FilterGraphs.
# Select cpu (default), cuda, or amdgpu as the first argument.
using AdaptiveOpticsSim, Proper, FilterGraphs, FilterGraphAlgorithms, LinearAlgebra, Test
using AdaptiveOpticsSim.AlgorithmGraphs
const B = AdaptiveOpticsSim.Backends
const O = AdaptiveOpticsSim.Optics

BLAS.set_num_threads(1)
B.set_fft_provider_threads!(1)

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

function llowfs_definition(target; output_sampling_m=(0.0008, 0.0007))
    shape = (8, 8)
    pitch = (0.0005, 0.0005)
    center = (5.0, 5.0)
    schema = "test.optical-cell-flux/1"
    steps = (
        (; op="add_opd", name="mirror", input="dm_opd", size=shape,
            axes=(1,2), sampling_m=pitch, center_pixel=center,
            schema="test.dm-opd/1", outside_opd_m=0.0),
        (; op="lens", focal_length_m=0.1),
        (; op="propagate", distance_m=0.07),
    )
    config = (; resolution=8, pupil_resolution=8, diameter_m=0.004, wavelength_um=1.4,
        output_schema=schema, steps)
    train = proper_optical_train_node_factory(T=Float32)(:optics, config, (;))
    pixels = pixel_integration_node(:pixels; T=Float32, input_shape=shape, output_shape=(5,6),
        input_sampling_m=pitch, output_sampling_m, input_center_pixel=center,
        output_center_pixel=(3.0, 3.5), input_schema=schema, output_schema="test.pixel-flux/1")
    buffer(value) = fill!(B.allocate_device_array(target, Float32, shape...), value)
    return algorithm_graph((train, pixels);
        inputs=(graph_input(:pupil_opd, :optics=>:pupil_opd, buffer(0)),
            graph_input(:pupil_amplitude, :optics=>:pupil_amplitude, buffer(1)),
            graph_input(:dm_opd, :optics=>:dm_opd, buffer(0))),
        links=(link(:optics=>:output, :pixels=>:input),),
        outputs=(graph_output(:intensity, :pixels=>:output),))
end

function exposure_owner(target, policy; definition=nothing)
    llowfs = definition === nothing ? optical_train(target, policy, 1.4, 0.0016) :
        prepare_algorithm_graph(definition; target, execution=policy)
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
        if hasproperty(owner.graphs.llowfs.inputs, :dm_opd)
            dm = Float32[1e-8 * seq * sin(i + j) for i in 1:8, j in 1:8]
            AlgorithmGraphs.copyto_graph_input!(owner.graphs.llowfs, :dm_opd, dm)
        end
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

@testset "Direct optical train and pixel integration preserve AOS exposure timing" begin
    oracle = exposure_owner(target, policy; definition=llowfs_definition(target))
    definition = llowfs_definition(target)
    original = exposure_owner(target, policy; definition)
    mismatched = llowfs_definition(target; output_sampling_m=(0.0009, 0.0007))
    @test_throws ArgumentError O.prepare_filtergraphs_external_pupil(original; llowfs_definition=mismatched)
    direct = O.prepare_filtergraphs_external_pupil(original; llowfs_definition=definition)
    execution = AlgorithmGraphs.optical_execution_graph(direct, :llowfs)
    @test execution isa (selection == "cpu" ? FilterGraphs.PreparedGraph : FilterGraphs.PreparedCapturedGraph)
    @test length(execution.nodes) == 2
    expected, actual = collect_exposures(oracle), collect_exposures(direct)
    @test length(expected) == length(actual) == 10
    for (a, b) in zip(expected, actual)
        @test (a.name, a.metadata) == (b.name, b.metadata)
        @test a.image ≈ b.image rtol=2e-6
    end
end

# Run with AdaptiveOpticsSim, Proper, FilterGraphAlgorithms and FilterGraphs.
# The first argument selects cpu (default), cuda, or amdgpu.
using AdaptiveOpticsSim, Proper, FilterGraphAlgorithms, FilterGraphs, Test

const A = AdaptiveOpticsSim.AlgorithmGraphs
const B = AdaptiveOpticsSim.Backends
const Ext = Base.get_extension(AdaptiveOpticsSim, :AdaptiveOpticsSimFilterGraphsExt)

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

@testset "FilterGraphs global-shutter CMOS acquisition" begin
    shape = (8, 10)
    sampling = (1e-3, 1e-3)
    center = (4.5, 5.5)
    rate_schema = "test.photon-rate/1"
    frame_schema = "test.cmos-frame/1"
    binning = selection == "cpu" ? 2 : 1
    detector = A.cmos_detector_acquisition_node(:camera; rows=shape[1], columns=shape[2],
        binning, pixel_scale_arcsec=0.697, wavelength_m=0.8e-6, exposure_duration_s=0.001896,
        quantum_efficiency=0.27, gain=15.848932, dark_current_e_per_pixel_s=0.96,
        bits=12, full_well_e=94_000, photon_noise=true, readout_noise=true,
        readout_noise_e=2.33, column_readout_noise_e=0.25, row_readout_noise_e=0.5,
        rng_seed=0x1811, photon_rate_schema=rate_schema, frame_schema)
    function detector_definition()
        rate = fill!(B.allocate_device_array(target, Float32, shape...), 0)
        return A.algorithm_graph((detector,);
            inputs=(A.graph_input(:photon_rate, :camera=>:photon_rate, rate),),
            outputs=(A.graph_output(:frame, :camera=>:frame),))
    end
    definition = detector_definition()
    substituted = Ext.replace_camera_nodes(detector_definition(), (:camera,))
    @test A.graph_input_contract(substituted, :photon_rate) == A.graph_input_contract(definition, :photon_rate)
    @test A.graph_output_contract(substituted, :frame) == A.graph_output_contract(definition, :frame)
    execution = selection == "cpu" ? A.StreamGraphExecution() : A.CapturedGraphExecution()
    baseline = A.prepare_algorithm_graph(definition; target, execution)
    bound = A.prepare_algorithm_graph(substituted; target, execution)

    recipe = AdaptiveOpticsSim.Optics.cmos_filter_graph(detector)
    pixels = A.pixel_integration_node(:pixels; input_shape=shape, output_shape=shape,
        input_sampling_m=sampling, output_sampling_m=sampling, input_center_pixel=center,
        output_center_pixel=center, input_schema=rate_schema, output_schema=rate_schema)
    pixel = (; name=:pixels, label=FilterGraphAlgorithms.algorithm_label(Ext.PixelIntegrationF32),
        config=(; definition=pixels))
    document = (; nodes=(pixel, recipe.document.nodes...), inputs=(; photon_rate="pixels:input"),
        links=("pixels:output"=>"camera:photon_rate",), recipe.document.outputs)
    checked = FilterGraphs.prepare_graph(document;
        algorithms=(Ext.PixelIntegrationF32, recipe.algorithms...))
    graph = if selection == "cuda"
        captured = FilterGraphs.prepare_captured_graph(checked, CUDA.device())
        @test captured.executable.executable isa CUDA.CuGraphExec
        captured
    elseif selection == "amdgpu"
        FilterGraphs.prepare_captured_graph(checked, AMDGPU.device())
    else
        checked
    end
    if selection != "cpu"
        @test graph isa FilterGraphs.PreparedCapturedGraph
        @test A.captured_graph_node_count(bound) == 1
    end

    rates = Float32[2e5 + 100 * (i + 3j) for i in 1:shape[1], j in 1:shape[2]]
    bad_rates = copy(rates)
    bad_rates[3, 5] = -1
    frame = zeros(Float32, div.(shape, binning))
    metadata(sequence) = (; photon_rate=FilterGraphAlgorithms.SampleMetadata(sequence;
        timestamp_nanoseconds=(sequence-1)*2_000_000,
        exposure_duration_nanoseconds=1_896_000, terminal=true))
    for aos in (baseline, bound)
        A.copyto_graph_input!(aos, :photon_rate, rates)
        A.step_graph!(aos)
    end
    @test FilterGraphs.process!((; frame), graph, (; photon_rate=rates), metadata(1)) === (; frame=true)
    first_frame = copy(frame)
    @test first_frame == Array(A.graph_output(baseline, :frame)) == Array(A.graph_output(bound, :frame))
    @test all(value -> isinteger(value) && 0 <= value <= 4095, first_frame)

    @test FilterGraphs.process!((; frame), graph, (; photon_rate=bad_rates), metadata(2)) ===
        Ext.InvalidOpticalInput
    @test frame == first_frame
    @test FilterGraphs.output_metadata(graph).frame == metadata(1).photon_rate
    A.step_graph!(baseline)
    A.step_graph!(bound)
    @test FilterGraphs.process!((; frame), graph, (; photon_rate=rates), metadata(3)) === (; frame=true)
    second_frame = copy(frame)
    @test second_frame == Array(A.graph_output(baseline, :frame)) == Array(A.graph_output(bound, :frame))
    @test second_frame != first_frame
    @test FilterGraphs.output_metadata(graph).frame == metadata(3).photon_rate

    A.reset_graph!(bound)
    A.step_graph!(bound)
    @test Array(A.graph_output(bound, :frame)) == first_frame
    FilterGraphs.reset!(graph)
    @test FilterGraphs.process!((; frame), graph, (; photon_rate=rates), metadata(1)) === (; frame=true)
    @test frame == first_frame
    @test FilterGraphs.process!((; frame), graph, (; photon_rate=rates), metadata(2)) === (; frame=true)
    @test frame == second_frame
end

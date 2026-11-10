# Run with AdaptiveOpticsSim, Proper, FilterGraphAlgorithms and FilterGraphs.
# The first argument selects cpu (default), cuda, or amdgpu.
using AdaptiveOpticsSim, Proper, FilterGraphAlgorithms, FilterGraphs, Test

const A = AdaptiveOpticsSim.AlgorithmGraphs
const Ext = Base.get_extension(AdaptiveOpticsSim, :AdaptiveOpticsSimFilterGraphsExt)

selection = get(ARGS, 1, "cpu")
if selection == "cuda"
    using CUDA
    @test CUDA.functional()
    CUDA.allowscalar(false)
elseif selection == "amdgpu"
    using AMDGPU
    @test AMDGPU.functional()
    AMDGPU.allowscalar(false)
else
    @assert selection == "cpu"
end

@testset "Upstream rejection preserves the FilterGraphs EMCCD noise sequence" begin
    shape = (8, 8)
    sampling = (1e-3, 1e-3)
    center = (4.5, 4.5)
    rate_schema = "test.photon-rate/1"
    pixels = A.pixel_integration_node(:pixels; input_shape=shape, output_shape=shape,
        input_sampling_m=sampling, output_sampling_m=sampling, input_center_pixel=center,
        output_center_pixel=center, input_schema=rate_schema, output_schema=rate_schema)
    detector = A.emccd_detector_acquisition_node(:camera; rows=8, columns=8,
        normalized_pupil_sampling=0.05, wavelength_m=0.55e-6, exposure_duration_s=500e-6,
        quantum_efficiency=0.82, rng_seed=12, photon_rate_schema=rate_schema,
        frame_schema="test.emccd/1", gain=70, readout_noise=true, readout_noise_e=5,
        excess_noise_factor=sqrt(2), clock_induced_charge_e_per_pixel_frame=0.001)
    document = (;
        nodes=(
            (; name=:pixels, label=FilterGraphAlgorithms.algorithm_label(Ext.PixelIntegrationF32), config=(; definition=pixels)),
            (; name=:camera, label=FilterGraphAlgorithms.algorithm_label(Ext.EMCCDF32), config=(; definition=detector)),
        ),
        links=("pixels:output"=>"camera:photon_rate",),
        inputs=(; photon_rate="pixels:input"),
        outputs=(; frame="camera:frame"),
    )
    checked = FilterGraphs.prepare_graph(document; algorithms=(Ext.PixelIntegrationF32, Ext.EMCCDF32))
    graph = if selection == "cuda"
        captured = FilterGraphs.prepare_captured_graph(checked, CUDA.device())
        @test captured.executable.executable isa CUDA.CuGraphExec
        captured
    elseif selection == "amdgpu"
        captured = FilterGraphs.prepare_captured_graph(checked, AMDGPU.device())
        @test captured isa FilterGraphs.PreparedCapturedGraph
        captured
    else
        checked
    end

    rates = Float32[1e4 + 100 * (i + 3j) for i in 1:shape[1], j in 1:shape[2]]
    rejected = copy(rates)
    rejected[3, 5] = -1
    frame = zeros(Float32, shape)
    metadata(seq) = (; photon_rate=FilterGraphAlgorithms.SampleMetadata(seq; timestamp_nanoseconds=(seq-1)*500_000,
        exposure_duration_nanoseconds=500_000, terminal=true))
    @test FilterGraphs.process!((; frame), graph, (; photon_rate=rates), metadata(1)) === (; frame=true)
    first_frame = copy(frame)
    @test FilterGraphs.process!((; frame), graph, (; photon_rate=rejected), metadata(2)) ===
        Ext.InvalidOpticalInput
    @test frame == first_frame
    @test FilterGraphs.output_metadata(graph).frame == metadata(1).photon_rate
    @test FilterGraphs.process!((; frame), graph, (; photon_rate=rates), metadata(3)) === (; frame=true)
    second_frame = copy(frame)
    @test second_frame != first_frame

    FilterGraphs.reset!(graph)
    @test FilterGraphs.process!((; frame), graph, (; photon_rate=rates), metadata(1)) === (; frame=true)
    @test frame == first_frame
    @test FilterGraphs.process!((; frame), graph, (; photon_rate=rates), metadata(2)) === (; frame=true)
    @test frame == second_frame
end

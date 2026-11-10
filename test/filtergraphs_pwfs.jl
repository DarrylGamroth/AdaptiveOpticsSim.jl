# Run with AdaptiveOpticsSim, Proper, FilterGraphAlgorithms and FilterGraphs.
# The first argument selects cpu (default), cuda, or amdgpu.
using AdaptiveOpticsSim, Proper, FilterGraphs, FilterGraphAlgorithms, LinearAlgebra, Test
using AdaptiveOpticsSim.AlgorithmGraphs
const A = AdaptiveOpticsSim.AlgorithmGraphs
const B = AdaptiveOpticsSim.Backends
const O = AdaptiveOpticsSim.Optics

BLAS.set_num_threads(1)
B.set_fft_provider_threads!(1)

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

@testset "Direct and AOS-bound FilterGraphAlgorithms PWFS" begin
    node = pyramid_rate_node(:pyramid; T=Float32, resolution=8, telescope_diameter_m=1.22,
        pupil_samples=4, modulation=2, modulation_points=4, n_pix_separation=2, n_pix_edge=1,
        source_wavelength_m=0.55e-6, source_photon_irradiance_m2_s=1.0,
        opd_schema="test.pupil-opd/1", photon_rate_schema="test.photon-rate/1")
    opd = fill!(B.allocate_device_array(target, Float32, 8, 8), 0)
    definition = algorithm_graph((node,); inputs=(graph_input(:opd, :pyramid=>:opd, opd),),
        outputs=(graph_output(:photon_rate, :pyramid=>:photon_rate),))
    baseline = prepare_algorithm_graph(definition; target, execution=policy)
    substituted = O.with_filtergraphs_pwfs(definition)
    @test A.graph_input_contract(substituted, :opd) == A.graph_input_contract(definition, :opd)
    @test A.graph_output_contract(substituted, :photon_rate) == A.graph_output_contract(definition, :photon_rate)
    @test_throws ArgumentError O.with_filtergraphs_pwfs(definition; node=:missing)
    bridge = prepare_algorithm_graph(substituted; target, execution=policy)
    recipe = O.pwfs_filter_graph(node)
    checked = FilterGraphs.prepare_graph(recipe.document; algorithms=recipe.algorithms)
    direct = if selection == "cuda"
        result = FilterGraphs.prepare_captured_graph(checked, CUDA.device())
        @test result.executable.executable isa CUDA.CuGraphExec
        result
    elseif selection == "amdgpu"
        FilterGraphs.prepare_captured_graph(checked, AMDGPU.device())
    else
        checked
    end
    if selection != "cpu"
        @test direct isa FilterGraphs.PreparedCapturedGraph
        @test captured_graph_node_count(bridge) == 1
    end
    photon_rate = zeros(Float32, size(graph_output(baseline, :photon_rate)))
    cases = (
        zeros(Float32, 8, 8),
        Float32[(i + 2j) * 1e-8 for i in 1:8, j in 1:8],
        Float32[2e-8 * sin(i + j) for i in 1:8, j in 1:8],
    )
    rates = map(cases) do host_opd
        A.copyto_graph_input!(baseline, :opd, host_opd)
        step_graph!(baseline)
        A.copyto_graph_input!(bridge, :opd, host_opd)
        step_graph!(bridge)
        @test FilterGraphs.process!((; photon_rate), direct, (; opd=host_opd)) === (; photon_rate=true)
        expected = Array(graph_output(baseline, :photon_rate))
        @test all(isfinite, photon_rate) && any(>(0), photon_rate)
        @test photon_rate == expected
        @test Array(graph_output(bridge, :photon_rate)) == expected
        copy(photon_rate)
    end
    @test rates[1] != rates[2] != rates[3]
    reset_graph!(bridge)
    FilterGraphs.reset!(direct)
    @test FilterGraphs.process!((; photon_rate), direct, (; opd=first(cases))) === (; photon_rate=true)
    @test photon_rate == first(rates)
end

@testset "FilterGraphs PWFS camera exposure and detector noise" begin
    rate_schema = "test.photon-rate/1"
    detector_schema = "test.detector-electrons/1"
    pyramid = pyramid_rate_node(:optics; T=Float32, resolution=8, telescope_diameter_m=1.22,
        pupil_samples=4, modulation=2, modulation_points=4, n_pix_separation=2, n_pix_edge=1,
        source_wavelength_m=0.55e-6, source_photon_irradiance_m2_s=1e6,
        opd_schema="test.pupil-opd/1", photon_rate_schema=rate_schema)
    detector = A.emccd_detector_acquisition_node(:camera; rows=12, columns=12,
        normalized_pupil_sampling=0.25, wavelength_m=0.55e-6, exposure_duration_s=0.0005,
        quantum_efficiency=0.9, rng_seed=2, gain=5, readout_noise=true, readout_noise_e=2,
        excess_noise_factor=sqrt(2), clock_induced_charge_e_per_pixel_frame=0.01,
        photon_rate_schema=rate_schema, frame_schema=detector_schema)
    layout = A.builtin_graph_node_types().detector_frame_layout_f32(:layout,
        (; input_size=(12,12), output_size=(14,12), origin=(0,0), scale=1, bias=0, output_element_type="uint16",
            input_schema=detector_schema, output_schema="test.raw-frame/1",
            range=(0,16383), rounding="nearest_even"), (;))
    opd = fill!(B.allocate_device_array(target, Float32, 8, 8), 0)
    definition = algorithm_graph((pyramid, detector, layout);
        inputs=(graph_input(:opd, :optics=>:opd, opd),),
        links=(link(:optics=>:photon_rate, :camera=>:photon_rate),
            link(:camera=>:frame, :layout=>:detector_frame)),
        outputs=(graph_output(:frame, :layout=>:frame),))
    baseline = prepare_algorithm_graph(definition; target, execution=policy)
    bridged = prepare_algorithm_graph(O.with_filtergraphs_pwfs_camera(definition); target, execution=policy)
    if selection == "cpu"
        optics = A.prepared_graph_node(bridged, :optics)
        camera = A.prepared_graph_node(bridged, :camera)
        @test optics.algorithm.workspace.opd === graph_input(bridged, :opd)
        @test optics.algorithm.workspace.photon_rate === optics.output === camera.input
        @test camera.algorithm.workspace.photon_rate === camera.input
        @test camera.algorithm.workspace.frame === camera.output
    end
    recipe = O.pwfs_camera_filter_graph(definition)
    checked = FilterGraphs.prepare_graph(recipe.document; algorithms=recipe.algorithms)
    direct = if selection == "cuda"
        graph = FilterGraphs.prepare_captured_graph(checked, CUDA.device())
        @test graph.executable.executable isa CUDA.CuGraphExec
        graph
    elseif selection == "amdgpu"
        FilterGraphs.prepare_captured_graph(checked, AMDGPU.device())
    else
        checked
    end
    @test length(direct.nodes) == 3
    if selection != "cpu"
        @test direct isa FilterGraphs.PreparedCapturedGraph
        @test captured_graph_node_count(bridged) == 3
    end
    frame = zeros(UInt16, 14, 12)
    host_opd = Float32[(i + 2j)*1e-8 for i in 1:8, j in 1:8]
    original_frames = map(1:3) do sequence
        A.copyto_graph_input!(baseline, :opd, host_opd)
        step_graph!(baseline)
        A.copyto_graph_input!(bridged, :opd, host_opd)
        step_graph!(bridged)
        sample = FilterGraphAlgorithms.SampleMetadata(sequence;
            timestamp_nanoseconds=(sequence-1)*500_000, exposure_duration_nanoseconds=500_000,
            terminal=true)
        @test FilterGraphs.process!((; frame), direct, (; opd=host_opd), (; opd=sample)) === (; frame=true)
        @test FilterGraphs.output_metadata(direct).frame == sample
        @test all(iszero, @view frame[13:14, :])
        @test Array(graph_output(bridged, :frame)) == Array(graph_output(baseline, :frame)) == frame
        copy(frame)
    end
    @test original_frames[1] != original_frames[2]
    FilterGraphs.reset!(direct)
    @test FilterGraphs.process!((; frame), direct, (; opd=host_opd)) === (; frame=true)
    @test frame == first(original_frames)

    if selection == "cpu" && isdefined(FilterGraphs, :prepare_cpu_loans)
        inputs = ntuple(_ -> (; opd=copy(host_opd)), 2)
        outputs = ntuple(_ -> (; frame=similar(frame)), 3)
        loans = FilterGraphs.prepare_cpu_loans(recipe.document; algorithms=recipe.algorithms,
            input_slots=inputs, output_slots=outputs)
        first_slot, second_slot = loans.nodes_by_input
        first_pyramid = first_slot[1].prepared.workspace
        second_pyramid = second_slot[1].prepared.workspace
        @test first_pyramid.opd === inputs[1].opd
        @test second_pyramid.opd === inputs[2].opd
        @test first_pyramid.owner.prepared.workspace === second_pyramid.owner.prepared.workspace
        @test first_slot[2].prepared === second_slot[2].prepared
        detector = first_slot[2].prepared.workspace
        @test first_pyramid.photon_rate === second_pyramid.photon_rate === detector.photon_rate
        @test detector.frame === only(loans.graph.node_outputs[2])
        @test FilterGraphs.process_cpu_loan!(loans, 1, 3) === (; frame=true)
        @test FilterGraphs.process_cpu_loan!(loans, 2, 1) === (; frame=true)
        @test outputs[3].frame == original_frames[1]
        @test outputs[1].frame == original_frames[2]
        FilterGraphs.reset!(loans)
        @test FilterGraphs.process_cpu_loan!(loans, 2, 2) === (; frame=true)
        @test outputs[2].frame == original_frames[1]
    end
end

import KernelAbstractions

function ingaas_test_graph_definition(input; kwargs...)
    config = merge((rows=size(input, 1), columns=size(input, 2),
            pixel_sampling_m=(12e-6, 14e-6), spectral_channel=:test_ingaas_band,
            exposure_duration_s=0.5, quantum_efficiency=0.75, gain=2,
            rng_seed=0xabc1, photon_rate_schema="test.ingaas.photon-rate/1",
            frame_schema="test.ingaas.frame/1", photon_noise=false, T=eltype(input)),
        (; kwargs...))
    node = ingaas_detector_acquisition_node(:detector; config...)
    return algorithm_graph((node,); name=:ingaas_test,
        inputs=(graph_input(:rate, :detector => :photon_rate, input),),
        outputs=(graph_output(:frame, :detector => :frame),))
end

function ingaas_graph_step_allocations(graph)
    step_graph!(graph)
    return @allocated step_graph!(graph)
end

@testset "InGaAs complete-frame graph node" begin
    detectors = AdaptiveOpticsSim.Detectors
    graphs = AdaptiveOpticsSim.AlgorithmGraphs
    backend_module = AdaptiveOpticsSim.Backends
    style = backend_module.AcceleratorStyle(KernelAbstractions.CPU())
    for T in (Float32, Float64)
        input = reshape(T.(1:15), 3, 5)
        graph = prepare_algorithm_graph(ingaas_test_graph_definition(input))
        owner = prepared_graph_node(graph, Val(:detector))
        output = graph_output(graph, Val(:frame))
        @test owner.prepared.input.values === input
        @test owner.output === output
        @test owner.prepared.detector.params.sensor isa InGaAsSensor
        @test owner.prepared.input.metadata.coordinate_domain isa MetricCoordinates
        @test owner.prepared.input.metadata.sampling == (T(12e-6), T(14e-6))
        @test owner.prepared.input.metadata.spectral ==
            IntegratedSpectralChannel(:test_ingaas_band)
        @test size(output) == size(input)
        @test @inferred(step_graph!(graph)) === graph
        @test output == input .* T(0.75)
        input .*= T(2)
        step_graph!(graph)
        @test output == input .* T(0.75)
        @test detectors._supports_enqueued_detector_acquisition(style,
            owner.prepared, owner.output)
        @test graphs.graph_node_capture_capability(owner) isa
            graphs.GraphNodeCaptureUnsupported # Host graph has no device counter RNG.
        if !coverage_instrumented()
            @test ingaas_graph_step_allocations(graph) == 0
        end
        reset_graph!(graph)
        @test all(iszero, output)
        @test all(iszero, owner.latent_buffer)
        @test graph_step_sequence(graph) == 0

        # Saturation occurs in charge before conversion gain; ADC quantization
        # then follows the existing Detector convention, with no axis mapping.
        saturated = prepare_algorithm_graph(ingaas_test_graph_definition(
            fill(T(20), 3, 5); full_well_e=5))
        step_graph!(saturated)
        @test graph_output(saturated, Val(:frame)) == fill(T(10), 3, 5)
        adc = prepare_algorithm_graph(ingaas_test_graph_definition(
            fill(T(4), 3, 5); bits=8, full_well_e=5))
        step_graph!(adc)
        @test graph_output(adc, Val(:frame)) == fill(T(153), 3, 5)
        adc_owner = prepared_graph_node(adc, Val(:detector))
        @test detectors._supports_enqueued_detector_acquisition(style,
            adc_owner.prepared, adc_owner.output)

        # A direct general InGaAs acquisition is the independent adapter oracle:
        # identical seeds and input sequences must produce identical frames.
        noisy_input = fill(T(20), 3, 5)
        noisy = prepare_algorithm_graph(ingaas_test_graph_definition(noisy_input;
            photon_noise=true, readout_noise=true, readout_noise_e=0.7,
            dark_current_e_per_pixel_s=2, glow_e_per_pixel_s=3))
        noisy_owner = prepared_graph_node(noisy, Val(:detector))
        reference_detector = Detector(exposure_duration=T(0.5), qe=T(0.75), gain=T(2),
            dark_current=T(2), noise=NoisePhotonReadout(T(0.7)),
            sensor=InGaAsSensor(glow_rate=3, T=T), response_model=NullFrameResponse(), T=T)
        reference = prepare_detector_acquisition(reference_detector, noisy_owner.prepared.input)
        reference_rng = backend_module.runtime_rng(0xabc1)
        step_graph!(noisy)
        first = copy(graph_output(noisy, Val(:frame)))
        @test first == capture!(reference, reference_rng)
        noisy_input .*= T(2)
        step_graph!(noisy)
        @test graph_output(noisy, Val(:frame)) == capture!(reference, reference_rng)
        @test graph_output(noisy, Val(:frame)) != first
        noisy_input ./= T(2)
        reset_graph!(noisy)
        step_graph!(noisy)
        @test graph_output(noisy, Val(:frame)) == first
        @test detectors._supports_enqueued_detector_acquisition(style,
            noisy_owner.prepared, noisy_owner.output)
        # Exercise the new asynchronous sensor dispatches with KA's CPU backend
        # and the genuine counter-state representation, without a GPU runtime.
        enqueued_detector = noisy_owner.prepared.detector
        frame = enqueued_detector.products.frame
        fill!(frame, zero(T))
        counter_rng() = backend_module._PreparedCounterRNG(
            backend_module._CounterRandomPlan(UInt64(42)),
            backend_module._CounterRandomState(UInt64[0]))
        enqueued_rng, expected_rng = counter_rng(), counter_rng()
        detectors._enqueue_detector_sensor_statistics!(style,
            enqueued_detector.params.sensor, enqueued_detector, enqueued_rng)
        expected_glow = fill(T(1.5), size(frame))
        backend_module.poisson_noise_async!(style, expected_rng, expected_glow)
        KernelAbstractions.synchronize(style.backend)
        @test frame == expected_glow
        @test enqueued_rng.state.draw_sequence == UInt64[1]
        @test detectors._enqueue_detector_pre_readout_gain!(style,
            enqueued_detector.params.sensor, enqueued_detector, enqueued_rng) === frame
        @test detectors._enqueue_detector_sensor_readout_noise!(style,
            enqueued_detector.params.sensor, enqueued_detector, enqueued_rng) === frame
        detectors._enqueue_detector_post_readout_gain!(style,
            enqueued_detector.params.sensor, enqueued_detector)
        KernelAbstractions.synchronize(style.backend)
        @test frame == T(2) .* expected_glow
        quiet_detector = owner.prepared.detector
        quiet_rng = counter_rng()
        fill!(quiet_detector.products.frame, zero(T))
        detectors._enqueue_detector_sensor_statistics!(style,
            quiet_detector.params.sensor, quiet_detector, quiet_rng)
        KernelAbstractions.synchronize(style.backend)
        @test all(iszero, quiet_detector.products.frame)
        @test quiet_rng.state.draw_sequence == UInt64[0]
        if !coverage_instrumented()
            @test ingaas_graph_step_allocations(noisy) == 0
        end

        # Non-null models stay effective in stream execution and fail the
        # existing captured acquisition admission, rather than disappearing.
        persistence_input = fill(T(8), 3, 5)
        persistence = prepare_algorithm_graph(ingaas_test_graph_definition(persistence_input;
            exposure_duration_s=1, quantum_efficiency=1, gain=4,
            persistence_model=ExponentialPersistence(0.25, 0.5)))
        persistent_owner = prepared_graph_node(persistence, Val(:detector))
        step_graph!(persistence)
        @test graph_output(persistence, Val(:frame)) == fill(T(32), 3, 5)
        @test persistent_owner.latent_buffer == fill(T(2), 3, 5)
        fill!(persistence_input, zero(T))
        step_graph!(persistence)
        @test graph_output(persistence, Val(:frame)) == fill(T(8), 3, 5)
        @test persistent_owner.latent_buffer == fill(T(1.5), 3, 5)
        @test !detectors._supports_enqueued_detector_acquisition(style,
            persistent_owner.prepared, persistent_owner.output)
        reset_graph!(persistence)
        @test all(iszero, persistent_owner.latent_buffer)

        nonlinear = prepare_algorithm_graph(ingaas_test_graph_definition(input;
            nonlinearity_model=SaturatingFrameNonlinearity(0.1)))
        nonlinear_owner = prepared_graph_node(nonlinear, Val(:detector))
        step_graph!(nonlinear)
        charge = input .* T(0.375)
        @test graph_output(nonlinear, Val(:frame)) ≈ T(2) .* charge ./ (1 .+ T(0.1) .* charge)
        @test !detectors._supports_enqueued_detector_acquisition(style,
            nonlinear_owner.prepared, nonlinear_owner.output)
        defects = PixelResponseNonuniformity(fill(T(1.1), 3, 5); T=T)
        defect_detector = Detector(noise=NoiseNone(), sensor=InGaAsSensor(T=T),
            response_model=NullFrameResponse(), defect_model=defects, T=T)
        defect_prepared = prepare_detector_acquisition(defect_detector, owner.prepared.input)
        @test !detectors._supports_enqueued_detector_acquisition(style,
            defect_prepared, similar(input))
    end

    input = ones(Float32, 3, 5)
    for bad in ((quantum_efficiency=1.1,), (exposure_duration_s=0,),
            (gain=0,), (dark_current_e_per_pixel_s=-1,), (glow_e_per_pixel_s=-1,),
            (glow_e_per_pixel_s=Inf,), (pixel_sampling_m=(0, 1e-6),),
            (spectral_channel=Symbol(""),), (bits=8,), (bits=0, full_well_e=10),
            (readout_noise_e=2,), (readout_noise=true, readout_noise_e=0),
            (rng_seed=-1,))
        @test_throws AlgorithmGraphError ingaas_test_graph_definition(input; bad...)
    end
    @test_throws MethodError ingaas_test_graph_definition(input;
        defect_model=NullDetectorDefectModel())
    @test_throws MethodError ingaas_test_graph_definition(input; binning=2)
    @test Base.isexported(graphs, :ingaas_detector_acquisition_node)
    @test Base.ispublic(graphs, :InGaAsDetectorAcquisitionNode)
    @test !Base.ispublic(AdaptiveOpticsSim, :InGaAsDetectorAcquisitionNode)
end

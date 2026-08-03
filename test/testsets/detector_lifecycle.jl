@testset "Detector lifecycle, readout, and output behavior" begin
    incremental_detector = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse())
    capture!(incremental_detector, ones(2, 2); rng=MersenneTwister(206),
        integration_duration=0.6)
    @test_throws InvalidConfiguration capture!(incremental_detector,
        ones(2, 2); rng=MersenneTwister(206), integration_duration=0.5)
    @test incremental_detector.state.integrated_time == 0.6
    @test !readout_ready(incremental_detector)
    incremental_pending_accum = copy(incremental_detector.state.accum_buffer)
    @test_throws DimensionMismatchError capture!(incremental_detector,
        ones(4, 4); rng=MersenneTwister(206), integration_duration=0.1)
    @test incremental_detector.state.integrated_time == 0.6
    @test !readout_ready(incremental_detector)
    @test incremental_detector.state.accum_buffer == incremental_pending_accum
    @test_throws InvalidConfiguration Detector(integration_time=0.0)
    @test_throws InvalidConfiguration Detector(integration_time=Inf)

    transition_values = fill(2.0, 2, 2)
    transition_map = detector_test_intensity_map(transition_values)
    transition_detector = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse())
    transition_plan = prepare_detector_acquisition(transition_detector,
        transition_map)
    transition_source = Source(band=:custom, wavelength=0.6e-6,
        photon_irradiance=1.0)
    temporal_calls = Ref(0)
    transition_temporal = FunctionFrameSource(t -> begin
        temporal_calls[] += 1
        transition_values
    end)

    capture!(transition_detector, transition_map, transition_plan;
        rng=MersenneTwister(208), integration_duration=0.25)
    pending_time = transition_detector.state.integrated_time
    pending_frame = copy(transition_detector.state.frame)
    pending_accum = copy(transition_detector.state.accum_buffer)
    @test pending_time == 0.25
    @test !readout_ready(transition_detector)

    @test_throws InvalidConfiguration capture!(transition_detector,
        transition_values; rng=MersenneTwister(209))
    @test_throws InvalidConfiguration capture!(transition_detector,
        transition_values, transition_source; rng=MersenneTwister(210))
    @test_throws InvalidConfiguration capture!(transition_detector,
        transition_temporal; rng=MersenneTwister(211))
    @test_throws InvalidConfiguration capture!(transition_detector,
        transition_map, transition_plan; rng=MersenneTwister(212))
    @test_throws InvalidConfiguration capture_with_quantum_efficiency!(
        transition_detector, transition_values, 0.25,
        MersenneTwister(213))
    pending_stack = fill(3.0, 1, 2, 2)
    @test_throws InvalidConfiguration AdaptiveOpticsSim.Detectors.capture_stack!(
        transition_detector, pending_stack, similar(pending_stack);
        rng=MersenneTwister(213))
    pending_generalized_output = fill(UInt8(7), 1, 2, 2)
    pending_generalized_input = fill(3.0, 1, 2, 2)
    @test_throws InvalidConfiguration AdaptiveOpticsSim.Detectors.capture_stack!(
        transition_detector, pending_generalized_output,
        pending_generalized_input; rng=MersenneTwister(213))
    @test all(==(UInt8(7)), pending_generalized_output)
    @test temporal_calls[] == 0
    @test transition_detector.state.integrated_time == pending_time
    @test !readout_ready(transition_detector)
    @test transition_detector.state.frame == pending_frame
    @test transition_detector.state.accum_buffer == pending_accum

    pending_prepare_values = fill(4.0, 4, 4)
    pending_prepare_map = detector_test_intensity_map(pending_prepare_values)
    pending_frame_shape = size(transition_detector.state.frame)
    pending_accum_shape = size(transition_detector.state.accum_buffer)
    @test_throws InvalidConfiguration prepare_detector_acquisition(
        transition_detector, pending_prepare_map)
    @test transition_detector.state.integrated_time == pending_time
    @test !readout_ready(transition_detector)
    @test size(transition_detector.state.frame) == pending_frame_shape
    @test size(transition_detector.state.accum_buffer) == pending_accum_shape
    @test transition_detector.state.frame == pending_frame
    @test transition_detector.state.accum_buffer == pending_accum

    capture!(transition_detector, transition_map, transition_plan;
        rng=MersenneTwister(214), integration_duration=0.75)
    @test readout_ready(transition_detector)
    @test iszero(transition_detector.state.integrated_time)
    @test all(iszero, transition_detector.state.accum_buffer)

    capture!(transition_detector, transition_values;
        rng=MersenneTwister(215))
    @test readout_ready(transition_detector)
    @test iszero(transition_detector.state.integrated_time)
    capture!(transition_detector, transition_values, transition_source;
        rng=MersenneTwister(216))
    @test readout_ready(transition_detector)
    @test iszero(transition_detector.state.integrated_time)
    capture!(transition_detector, transition_temporal;
        rng=MersenneTwister(217))
    @test temporal_calls[] > 0
    @test readout_ready(transition_detector)
    @test iszero(transition_detector.state.integrated_time)
    capture!(transition_detector, transition_map, transition_plan;
        rng=MersenneTwister(218))
    @test readout_ready(transition_detector)
    @test iszero(transition_detector.state.integrated_time)
    @test all(iszero, transition_detector.state.accum_buffer)
    explicit_qe_frame = capture_with_quantum_efficiency!(
        transition_detector, transition_values, 0.25,
        MersenneTwister(219))
    @test explicit_qe_frame == transition_values .* 0.25
    @test readout_ready(transition_detector)
    @test iszero(transition_detector.state.integrated_time)

    glow_full_duration = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, response_model=NullFrameResponse(),
        sensor=InGaAsSensor(glow_rate=20.0))
    glow_short_duration = Detector(integration_time=0.25,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
        sensor=InGaAsSensor(glow_rate=20.0))
    prepare_detector_buffers!(glow_full_duration, (2, 2))
    prepare_detector_buffers!(glow_short_duration, (2, 2))
    fill!(glow_full_duration.state.frame, 0.0)
    fill!(glow_short_duration.state.frame, 0.0)
    apply_sensor_statistics!(glow_full_duration.params.sensor,
        glow_full_duration, MersenneTwister(207), 0.25)
    apply_sensor_statistics!(glow_short_duration.params.sensor,
        glow_short_duration, MersenneTwister(207), 0.25)
    @test glow_full_duration.state.frame == glow_short_duration.state.frame

    qe_curve = SampledQuantumEfficiency([0.50e-6, 0.60e-6], [0.2, 0.8])
    @test qe_at(qe_curve, 0.55e-6) ≈ 0.5
    @test qe_at(qe_curve, 0.70e-6) == 0.0
    det_qe_curve = Detector(integration_time=1.0, noise=NoiseNone(), qe=qe_curve, binning=1,
        response_model=NullFrameResponse())
    @test det_qe_curve.params.qe ≈ 0.8
    @test capture!(det_qe_curve, ones(2, 2); rng=MersenneTwister(30)) ≈ fill(0.8, 2, 2)
    src_qe = Source(wavelength=0.55e-6)
    @test effective_qe(det_qe_curve, src_qe) ≈ 0.5
    @test capture!(det_qe_curve, ones(2, 2), src_qe; rng=MersenneTwister(31)) ≈ fill(0.5, 2, 2)
    generalized_qe_detector = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=qe_curve, psf_sampling=2,
        response_model=NullFrameResponse())
    generalized_qe_input = ones(2, 4, 4)
    generalized_qe_output = zeros(2, 2, 2)
    AdaptiveOpticsSim.Detectors.capture_stack!(generalized_qe_detector,
        generalized_qe_output, generalized_qe_input, src_qe;
        rng=MersenneTwister(31))
    @test generalized_qe_output ≈ fill(2.0, 2, 2, 2)
    @test readout_ready(generalized_qe_detector)
    @test iszero(generalized_qe_detector.state.integrated_time)
    spectral_qe = with_spectrum(Source(wavelength=0.55e-6),
        SpectralBundle([0.50e-6, 0.60e-6], [0.25, 0.75]))
    @test effective_qe(det_qe_curve, spectral_qe) ≈ 0.65
    @test capture!(det_qe_curve, ones(2, 2), spectral_qe; rng=MersenneTwister(32)) ≈ fill(0.65, 2, 2)
    @test_throws InvalidConfiguration ScalarQuantumEfficiency(1.5)
    @test_throws InvalidConfiguration SampledQuantumEfficiency([0.60e-6, 0.50e-6], [0.8, 0.2])
    @test_throws InvalidConfiguration SampledQuantumEfficiency([0.50e-6, 0.60e-6], [0.2, 1.2])

    det_tuple = Detector(integration_time=1.0, noise=(NoisePhoton(), NoiseReadout(0.5)),
        qe=1.0, binning=1)
    @test det_tuple.noise isa NoisePhotonReadout
    @test AdaptiveOpticsSim.Detectors.detector_execution_strategy(AdaptiveOpticsSim.Backends.execution_style(det_tuple.state.frame), det_tuple) isa AdaptiveOpticsSim.Detectors.DetectorDirectStrategy

    det_sat = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1, full_well=5.0)
    frame_sat = capture!(det_sat, fill(10.0, 4, 4); rng=MersenneTwister(2))
    @test maximum(frame_sat) == 5.0

    det_adc = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1, bits=8, full_well=10.0)
    frame_adc = capture!(det_adc, fill(10.0, 4, 4); rng=MersenneTwister(2))
    @test frame_adc isa Matrix{UInt8}
    @test output_frame(det_adc) === frame_adc
    @test maximum(frame_adc) == 0xff
    @test minimum(frame_adc) >= 0x00
    @test eltype(det_adc.state.frame) == Float64
    metadata_adc = detector_export_metadata(det_adc)
    @test metadata_adc.noise == :none
    @test metadata_adc.sensor == :ccd
    @test metadata_adc.output_type == UInt8
    @test metadata_adc.frame_size == (4, 4)
    @test metadata_adc.output_size == (4, 4)
    @test_throws InvalidConfiguration Detector(noise=NoiseNone(), bits=8)
    @test_throws InvalidConfiguration Detector(noise=NoiseNone(), bits=0,
        full_well=10.0)

    det_adc_float = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        bits=8, full_well=10.0, output_type=Float32)
    frame_adc_float = capture!(det_adc_float, fill(10.0, 4, 4); rng=MersenneTwister(2))
    @test frame_adc_float isa Matrix{Float32}
    @test maximum(frame_adc_float) == Float32(255.0)

    det_adc_window_corr = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        bits=8, full_well=100.0, readout_window=FrameWindow(2:5, 3:7), output_type=UInt16,
        sensor=HgCdTeSensor(sampling_mode=SingleRead()),
        response_model=NullFrameResponse(),
        correction_model=CompositeFrameReadoutCorrection((
            ReferenceRowCommonModeCorrection(1),
            ReferenceColumnCommonModeCorrection(1),
        )))
    adc_window_in = reshape(collect(1.0:96.0), 2, 6, 8)
    adc_window_out = Array{UInt16}(undef, 2, 4, 5)
    generalized_adc_window = AdaptiveOpticsSim.Detectors.capture_stack!(det_adc_window_corr, adc_window_out, copy(adc_window_in);
        rng=MersenneTwister(10))
    @test size(generalized_adc_window) == (2, 4, 5)
    @test generalized_adc_window[1, :, :] == capture!(Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
            bits=8, full_well=100.0, readout_window=FrameWindow(2:5, 3:7), output_type=UInt16,
            sensor=HgCdTeSensor(sampling_mode=SingleRead()),
            response_model=NullFrameResponse(),
            correction_model=CompositeFrameReadoutCorrection((
                ReferenceRowCommonModeCorrection(1),
                ReferenceColumnCommonModeCorrection(1),
            ))),
        @view(adc_window_in[1, :, :]); rng=MersenneTwister(10))
    @test generalized_adc_window[2, :, :] == capture!(Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
            bits=8, full_well=100.0, readout_window=FrameWindow(2:5, 3:7), output_type=UInt16,
            sensor=HgCdTeSensor(sampling_mode=SingleRead()),
            response_model=NullFrameResponse(),
            correction_model=CompositeFrameReadoutCorrection((
                ReferenceRowCommonModeCorrection(1),
                ReferenceColumnCommonModeCorrection(1),
            ))),
        @view(adc_window_in[2, :, :]); rng=MersenneTwister(10))

    det_window = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        readout_window=FrameWindow(2:3, 2:4))
    psf_window = reshape(collect(1.0:16.0), 4, 4)
    frame_window = copy(capture!(det_window, psf_window; rng=MersenneTwister(2)))
    @test size(frame_window) == (2, 3)
    @test frame_window == psf_window[2:3, 2:4]
    meta_window = detector_export_metadata(det_window)
    @test meta_window.window_rows == (2, 3)
    @test meta_window.window_cols == (2, 4)
    det_window_oob = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        readout_window=FrameWindow(2:5, 1:2))
    @test_throws DimensionMismatchError capture!(det_window_oob, psf_window; rng=MersenneTwister(2))
    @test_throws InvalidConfiguration FrameWindow(0:1, 1:2)

    det_dark = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1, dark_current=100.0)
    frame_dark = capture!(det_dark, zeros(4, 4); rng=MersenneTwister(2))
    @test sum(frame_dark) > 0

    det_buffered = Detector(integration_time=2.0, noise=NoiseNone(), qe=1.0, binning=1)
    frame_partial = copy(capture!(det_buffered, fill(1.0, 4, 4); rng=MersenneTwister(2), integration_duration=1.0))
    @test !readout_ready(det_buffered)
    @test sum(frame_partial) == 16.0
    frame_buffered = copy(capture!(det_buffered, fill(1.0, 4, 4); rng=MersenneTwister(2), integration_duration=1.0))
    @test readout_ready(det_buffered)
    @test sum(frame_buffered) == 32.0
    reset_integration!(det_buffered)
    @test readout_ready(det_buffered)
    @test det_buffered.state.integrated_time == 0.0
    metadata_buffered = detector_export_metadata(det_buffered)
    @test metadata_buffered.output_size == size(output_frame(det_buffered))
    @test metadata_buffered.psf_sampling == 1
    @test metadata_buffered.binning == 1

    det_background_flux = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        background_flux=2.0)
    frame_background_flux = capture!(det_background_flux, zeros(4, 4); rng=MersenneTwister(2))
    @test sum(frame_background_flux) > 0

    det_background_map = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        background_map=fill(1.0, 4, 4))
    frame_background_map = capture!(det_background_map, zeros(4, 4); rng=MersenneTwister(2))
    @test frame_background_map == fill(-1.0, 4, 4)

    cube = Array{Float64}(undef, 2, 4, 4)
    cube[1, :, :] .= fill(1.0, 4, 4)
    cube[2, :, :] .= fill(2.0, 4, 4)
    scratch = similar(cube)
    det_stack = Detector(integration_time=1.0, noise=NoiseNone(), qe=0.5, binning=1)
    AdaptiveOpticsSim.Detectors.capture_stack!(det_stack, cube, scratch; rng=MersenneTwister(10))
    @test cube[1, :, :] ≈ fill(0.5, 4, 4)
    @test cube[2, :, :] ≈ fill(1.0, 4, 4)

    allocation_stack_detector = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse())
    allocation_stack_cube = ones(2, 4, 4)
    allocation_stack_scratch = similar(allocation_stack_cube)
    @test_detector_allocation fixed_stack_capture_allocations(
        allocation_stack_detector,
        allocation_stack_cube, allocation_stack_scratch, Xoshiro(10)) == 0

    psf = reshape(Float64.(1:256), 16, 16)
    det_fused = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, psf_sampling=2, binning=2)
    frame_fused = copy(AdaptiveOpticsSim.Detectors.fill_frame!(det_fused, psf, 1.0))
    manual_mid = zeros(Float64, 8, 8)
    manual_out = zeros(Float64, 4, 4)
    AdaptiveOpticsSim.bin2d!(manual_mid, psf, 2)
    AdaptiveOpticsSim.bin2d!(manual_out, manual_mid, 2)
    @test frame_fused == manual_out
end

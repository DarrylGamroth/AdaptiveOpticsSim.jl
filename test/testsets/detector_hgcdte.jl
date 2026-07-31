@testset "HgCdTe detector families" begin
    zero_psf = zeros(4, 4)
    det_hgcdte_single = Detector(integration_time=1.0, noise=NoiseReadout(4.0), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeSensor())
    frame_hgcdte_single = copy(capture!(det_hgcdte_single, zero_psf; rng=MersenneTwister(16)))
    single_products = readout_products(det_hgcdte_single)
    @test single_products isa HgCdTeReadoutProducts
    @test single_products isa MultiReadFrameReadoutProducts
    @test detector_reference_frame(det_hgcdte_single) === nothing
    @test detector_signal_frame(det_hgcdte_single) !== nothing
    @test detector_combined_frame(det_hgcdte_single) == frame_hgcdte_single
    @test detector_reference_cube(det_hgcdte_single) === nothing
    @test detector_signal_cube(det_hgcdte_single) !== nothing
    @test detector_read_cube(det_hgcdte_single) === nothing
    @test detector_read_times(det_hgcdte_single) === nothing

    conventional_detector = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, binning=2)
    @test fieldtype(typeof(conventional_detector.state),
        :readout_products) === NoFrameReadoutProducts
    @test fieldtype(typeof(det_hgcdte_single.state), :readout_products) !== FrameReadoutProducts

    struct DummyReadoutProducts{A,V} <: FrameReadoutProducts
        signal_frame::A
        read_times::V
    end
    AdaptiveOpticsSim.Detectors.detector_signal_frame(products::DummyReadoutProducts) = products.signal_frame
    AdaptiveOpticsSim.Detectors.detector_read_times(products::DummyReadoutProducts) = products.read_times

    dummy_products = DummyReadoutProducts(fill(3.0, 2, 2), [0.25, 0.5])
    @test detector_reference_frame(dummy_products) === nothing
    @test detector_signal_frame(dummy_products) == fill(3.0, 2, 2)
    @test detector_combined_frame(dummy_products) === nothing
    @test detector_reference_cube(dummy_products) === nothing
    @test detector_signal_cube(dummy_products) === nothing
    @test detector_read_cube(dummy_products) === nothing
    @test detector_read_times(dummy_products) == [0.25, 0.5]

    det_hgcdte_ndr = Detector(integration_time=1.0, noise=NoiseReadout(4.0), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeSensor(
            sampling_mode=AveragedNonDestructiveReads(4)))
    frame_hgcdte_ndr = copy(capture!(det_hgcdte_ndr, zero_psf; rng=MersenneTwister(16)))
    @test std(vec(frame_hgcdte_ndr)) < std(vec(frame_hgcdte_single))
    @test supports_nondestructive_reads(det_hgcdte_ndr.params.sensor)
    @test supports_readout_correction(det_hgcdte_ndr.params.sensor)
    @test supports_read_cube(det_hgcdte_ndr.params.sensor)
    hgcdte_meta = detector_export_metadata(det_hgcdte_ndr)
    @test hgcdte_meta.sampling_mode == :averaged_non_destructive_reads
    @test hgcdte_meta.sampling_reads == 4
    @test hgcdte_meta.sampling_reference_reads == 0
    @test hgcdte_meta.sampling_signal_reads == 4
    @test hgcdte_meta.readout_sigma == 2.0
    @test hgcdte_meta.provides_signal_frame
    @test !hgcdte_meta.provides_reference_frame
    @test hgcdte_meta.provides_combined_frame
    @test !hgcdte_meta.provides_reference_cube
    @test hgcdte_meta.provides_signal_cube
    @test hgcdte_meta.provides_read_cube
    @test hgcdte_meta.signal_cube_reads == 4
    @test hgcdte_meta.read_cube_reads == 4
    @test detector_combined_frame(det_hgcdte_ndr) == frame_hgcdte_ndr
    @test size(detector_signal_cube(det_hgcdte_ndr)) == (4, 4, 4)
    @test length(detector_read_times(det_hgcdte_ndr)) == 4
    det_hgcdte_cds = Detector(integration_time=1.0, noise=NoiseReadout(4.0), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeSensor(
            sampling_mode=CorrelatedDoubleSampling()))
    frame_hgcdte_cds = copy(capture!(det_hgcdte_cds, zero_psf; rng=MersenneTwister(16)))
    @test std(vec(frame_hgcdte_cds)) > std(vec(frame_hgcdte_single))
    cds_meta = detector_export_metadata(det_hgcdte_cds)
    @test cds_meta.sampling_mode == :correlated_double_sampling
    @test cds_meta.sampling_reference_reads == 1
    @test cds_meta.sampling_signal_reads == 1
    @test cds_meta.readout_sigma ≈ 4.0 * sqrt(2.0)
    @test cds_meta.provides_reference_cube
    @test cds_meta.provides_signal_cube
    @test cds_meta.reference_cube_reads == 1
    @test cds_meta.signal_cube_reads == 1
    @test size(detector_reference_cube(det_hgcdte_cds)) == (4, 4, 1)
    @test size(detector_signal_cube(det_hgcdte_cds)) == (4, 4, 1)
    @test detector_combined_frame(det_hgcdte_cds) ≈
        detector_signal_frame(det_hgcdte_cds) .- detector_reference_frame(det_hgcdte_cds)
    det_hgcdte_fowler = Detector(integration_time=1.0, noise=NoiseReadout(4.0), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeSensor(
            sampling_mode=FowlerSampling(8)))
    frame_hgcdte_fowler = copy(capture!(det_hgcdte_fowler, zero_psf; rng=MersenneTwister(16)))
    fowler_meta = detector_export_metadata(det_hgcdte_fowler)
    @test fowler_meta.sampling_mode == :fowler_sampling
    @test fowler_meta.sampling_reads == 16
    @test fowler_meta.sampling_reference_reads == 8
    @test fowler_meta.sampling_signal_reads == 8
    @test fowler_meta.readout_sigma == 2.0
    @test fowler_meta.reference_cube_reads == 8
    @test fowler_meta.signal_cube_reads == 8
    @test detector_combined_frame(det_hgcdte_fowler) ≈
        detector_signal_frame(det_hgcdte_fowler) .- detector_reference_frame(det_hgcdte_fowler)

    ramp_input = fill(5.0, 4, 4)
    ramp_detector = Detector(integration_time=2.0, noise=NoiseNone(),
        qe=1.0, gain=1.0, response_model=NullFrameResponse(),
        sensor=HgCdTeSensor(read_time=0.1,
            sampling_mode=UpTheRampSampling(5)))
    ramp_rng = MersenneTwister(161)
    ramp_frame = copy(capture!(ramp_detector, ramp_input, ramp_rng))
    @test ramp_frame == fill(10.0, 4, 4)
    @test detector_ramp_slope(ramp_detector) == fill(5.0, 4, 4)
    @test detector_ramp_intercept(ramp_detector) == zeros(4, 4)
    @test detector_signal_frame(ramp_detector) == ramp_frame
    @test detector_combined_frame(ramp_detector) == ramp_frame
    @test detector_reference_frame(ramp_detector) === nothing
    @test detector_reference_cube(ramp_detector) === nothing
    @test size(detector_signal_cube(ramp_detector)) == (4, 4, 5)
    @test detector_signal_cube(ramp_detector) ===
        detector_read_cube(ramp_detector)
    @test detector_ramp_cube(ramp_detector) ===
        detector_read_cube(ramp_detector)
    @test detector_ramp_times(ramp_detector) ===
        detector_read_times(ramp_detector)
    @test AdaptiveOpticsSim.Detectors.detector_ramp_acquisition(
        ramp_detector) == :synthesized_final_charge
    @test_throws InvalidConfiguration begin
        AdaptiveOpticsSim.Detectors.ensure_up_the_ramp_products!(
            ramp_detector, 5; acquisition=:invalid)
    end
    @test detector_read_times(ramp_detector) == [0.0, 0.5, 1.0, 1.5, 2.0]
    @test vec(Array(detector_read_cube(ramp_detector)[1, 1, :])) ==
        [0.0, 2.5, 5.0, 7.5, 10.0]
    ramp_meta = detector_export_metadata(ramp_detector)
    @test ramp_meta.sampling_mode == :up_the_ramp
    @test ramp_meta.sampling_reads == 5
    @test ramp_meta.sampling_reference_reads == 0
    @test ramp_meta.sampling_signal_reads == 5
    @test ramp_meta.sampling_read_time == 0.1
    @test ramp_meta.sampling_wallclock_time == 2.1
    @test ramp_meta.provides_signal_frame
    @test ramp_meta.provides_combined_frame
    @test ramp_meta.provides_signal_cube
    @test ramp_meta.provides_read_cube
    @test ramp_meta.signal_cube_reads == 5
    @test ramp_meta.read_cube_reads == 5
    @test supports_up_the_ramp(ramp_detector.params.sensor)
    capture!(ramp_detector, ramp_input, ramp_rng)
    if coverage_instrumented()
        @test_skip "up-the-ramp allocation assertions are disabled under coverage instrumentation"
    else
        @test @allocated(capture!(ramp_detector, ramp_input, ramp_rng)) == 0
    end
    ramp_stack = fill(5.0, 2, 4, 4)
    @test_throws InvalidConfiguration capture_stack!(ramp_detector,
        ramp_stack, similar(ramp_stack), MersenneTwister(161))

    ramp_window_detector = Detector(integration_time=2.0,
        noise=NoiseNone(), qe=1.0, gain=1.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeSensor(read_time=0.1,
            sampling_mode=UpTheRampSampling(5)),
        readout_window=FrameWindow(2:3, 2:4))
    ramp_window_rng = MersenneTwister(162)
    ramp_window_frame = capture!(ramp_window_detector, ramp_input,
        ramp_window_rng)
    @test size(ramp_window_frame) == (2, 3)
    @test ramp_window_frame == fill(10.0, 2, 3)
    @test size(detector_ramp_slope(ramp_window_detector)) == (2, 3)
    @test size(detector_read_cube(ramp_window_detector)) == (2, 3, 5)
    capture!(ramp_window_detector, ramp_input, ramp_window_rng)
    if coverage_instrumented()
        @test_skip "up-the-ramp allocation assertions are disabled under coverage instrumentation"
    else
        @test @allocated(capture!(ramp_window_detector, ramp_input,
            ramp_window_rng)) == 0
    end

    ramp_noise_detector = Detector(integration_time=1.0,
        noise=NoiseReadout(4.0), qe=1.0, gain=1.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeSensor(read_time=0.01,
            sampling_mode=UpTheRampSampling(16)))
    ramp_noise_frame = copy(capture!(ramp_noise_detector, zeros(64, 64);
        rng=MersenneTwister(163)))
    expected_ramp_sigma = 4.0 * sqrt(12 * 15 / (16 * 17))
    @test isapprox(std(ramp_noise_frame), expected_ramp_sigma; rtol=0.12)
    @test detector_export_metadata(ramp_noise_detector).readout_sigma ≈
        expected_ramp_sigma

    @test_throws InvalidConfiguration validate_frame_sampling_mode(
        UpTheRampSampling(1))
    @test_throws InvalidConfiguration CCDSensor(
        sampling_mode=UpTheRampSampling(4))
    @test_throws InvalidConfiguration HgCdTeSensor(
        sampling_mode=SkipperSampling(4))
    invalid_ramp_schedule = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0,
        sensor=HgCdTeSensor(read_time=0.3,
            sampling_mode=UpTheRampSampling(5)))
    @test_throws InvalidConfiguration capture!(invalid_ramp_schedule,
        ones(4, 4); rng=MersenneTwister(164))

    multiread_noise_fixture = zeros(64, 64)
    single_std = std(vec(copy(capture!(det_hgcdte_single, multiread_noise_fixture;
        rng=MersenneTwister(160)))))
    ndr_std = std(vec(copy(capture!(det_hgcdte_ndr, multiread_noise_fixture;
        rng=MersenneTwister(160)))))
    cds_std = std(vec(copy(capture!(det_hgcdte_cds, multiread_noise_fixture;
        rng=MersenneTwister(160)))))
    fowler_std = std(vec(copy(capture!(det_hgcdte_fowler, multiread_noise_fixture;
        rng=MersenneTwister(160)))))
    @test isapprox(single_std, 4.0; rtol=0.15)
    @test isapprox(ndr_std, 2.0; rtol=0.15)
    @test isapprox(cds_std, 4.0 * sqrt(2.0); rtol=0.15)
    @test isapprox(fowler_std, 2.0; rtol=0.15)
    det_hgcdte_timed_single = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        dark_current=1000.0, gain=1.0,
        sensor=HgCdTeSensor(read_time=1.0))
    frame_hgcdte_timed_single = copy(capture!(det_hgcdte_timed_single, zero_psf; rng=MersenneTwister(17)))
    det_hgcdte_timed_cds = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        dark_current=1000.0, gain=1.0,
        sensor=HgCdTeSensor(read_time=1.0,
            sampling_mode=CorrelatedDoubleSampling()))
    frame_hgcdte_timed_cds = copy(capture!(det_hgcdte_timed_cds, zero_psf; rng=MersenneTwister(17)))
    @test sum(frame_hgcdte_timed_cds) > sum(frame_hgcdte_timed_single)
    det_hgcdte_timed_glow = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        dark_current=3.0, gain=1.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeSensor(glow_rate=2.0, read_time=1.0,
            sampling_mode=CorrelatedDoubleSampling()))
    frame_hgcdte_timed_glow = copy(capture!(det_hgcdte_timed_glow, zero_psf; rng=MersenneTwister(17)))
    det_hgcdte_timed_noglow = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        dark_current=3.0, gain=1.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeSensor(read_time=1.0,
            sampling_mode=CorrelatedDoubleSampling()))
    frame_hgcdte_timed_noglow = copy(capture!(det_hgcdte_timed_noglow, zero_psf; rng=MersenneTwister(17)))
    @test sum(frame_hgcdte_timed_glow) > sum(frame_hgcdte_timed_noglow)
    timed_meta = detector_export_metadata(det_hgcdte_timed_cds)
    @test timed_meta.sampling_read_time == 1.0
    @test timed_meta.sampling_wallclock_time == 3.0
    det_hgcdte_windowed = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeSensor(read_time=1.0,
            sampling_mode=CorrelatedDoubleSampling()),
        readout_window=FrameWindow(2:3, 2:3))
    capture!(det_hgcdte_windowed, fill(10.0, 4, 4); rng=MersenneTwister(18))
    windowed_meta = detector_export_metadata(det_hgcdte_windowed)
    @test windowed_meta.sampling_read_time == 1.0
    @test windowed_meta.sampling_wallclock_time == 3.0
    @test detector_combined_frame(det_hgcdte_windowed) !== nothing
    @test size(detector_reference_cube(det_hgcdte_windowed)) == (2, 2, 1)
    @test size(detector_signal_cube(det_hgcdte_windowed)) == (2, 2, 1)
    @test size(detector_signal_frame(det_hgcdte_windowed)) == (2, 2)
    @test size(detector_read_cube(det_hgcdte_windowed)) == (2, 2, 2)
    @test detector_read_times(det_hgcdte_windowed) == [1.0, 2.0]
    @test detector_combined_frame(det_hgcdte_windowed) ≈
        detector_signal_frame(det_hgcdte_windowed) .- detector_reference_frame(det_hgcdte_windowed)
    det_hgcdte_windowed_corrected = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeSensor(read_time=1.0,
            sampling_mode=CorrelatedDoubleSampling()),
        readout_window=FrameWindow(2:3, 2:3),
        correction_model=CompositeFrameReadoutCorrection((
            ReferenceRowCommonModeCorrection(1),
            ReferenceColumnCommonModeCorrection(1),
        )))
    row_pattern = repeat(reshape([1.0, 2.0, 3.0, 4.0], :, 1), 1, 4)
    col_pattern = repeat(reshape([1.0, 2.0, 3.0, 4.0], 1, :), 4, 1)
    windowed_corrected_input = row_pattern .+ col_pattern
    windowed_corrected = capture!(det_hgcdte_windowed_corrected, windowed_corrected_input; rng=MersenneTwister(18))
    @test maximum(abs, windowed_corrected) < 1e-6
    @test detector_combined_frame(det_hgcdte_windowed_corrected) ≈
        detector_signal_frame(det_hgcdte_windowed_corrected) .- detector_reference_frame(det_hgcdte_windowed_corrected)
    @test size(detector_read_cube(det_hgcdte_windowed_corrected)) == (2, 2, 2)
    det_hgcdte_corrected = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeSensor(),
        response_model=NullFrameResponse(),
        correction_model=ReferencePixelCommonModeCorrection(1, 1))
    corrected_frame = capture!(det_hgcdte_corrected, fill(5.0, 4, 4); rng=MersenneTwister(19))
    corrected_meta = detector_export_metadata(det_hgcdte_corrected)
    @test corrected_meta.readout_correction == :reference_pixel_common_mode
    @test corrected_meta.correction_edge_rows == 1
    @test corrected_meta.correction_edge_cols == 1
    @test abs(mean(corrected_frame)) < 1e-6

    calibration_input = reshape(collect(1.0:64.0), 8, 8)
    calibration_sensor = HgCdTeSensor(
        read_time=0.1, sampling_mode=CorrelatedDoubleSampling())
    calibration_correction = ReferencePixelCommonModeCorrection(1, 1)
    calibration_detector = Detector(noise=NoiseNone(), qe=1.0, gain=3.0,
        sensor=calibration_sensor, response_model=NullFrameResponse(),
        correction_model=calibration_correction)
    capture_detector = Detector(noise=NoiseNone(), qe=1.0, gain=3.0,
        sensor=calibration_sensor, response_model=NullFrameResponse(),
        correction_model=calibration_correction)
    deterministic_reference = copy(
        AdaptiveOpticsSim.WavefrontSensors.detector_calibration_frame!(
            calibration_detector, calibration_input, 1.0))
    noiseless_cds = copy(capture!(capture_detector, calibration_input;
        rng=MersenneTwister(1901)))
    @test deterministic_reference == noiseless_cds
    @test iszero(calibration_detector.state.integrated_time)
    @test readout_ready(calibration_detector)
    @test readout_products(calibration_detector) isa NoFrameReadoutProducts

    ramp_calibration_sensor = HgCdTeSensor(
        read_time=0.1,
        sampling_mode=UpTheRampSampling(4))
    ramp_calibration_detector = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, gain=3.0,
        sensor=ramp_calibration_sensor, response_model=NullFrameResponse(),
        correction_model=calibration_correction)
    ramp_capture_detector = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, gain=3.0,
        sensor=ramp_calibration_sensor, response_model=NullFrameResponse(),
        correction_model=calibration_correction)
    ramp_reference = copy(
        AdaptiveOpticsSim.WavefrontSensors.detector_calibration_frame!(
            ramp_calibration_detector, calibration_input, 1.0))
    noiseless_ramp = copy(capture!(ramp_capture_detector,
        calibration_input; rng=MersenneTwister(1902)))
    @test ramp_reference ≈ noiseless_ramp atol=1e-12 rtol=1e-12
    @test iszero(ramp_calibration_detector.state.integrated_time)
    @test readout_ready(ramp_calibration_detector)
    @test readout_products(ramp_calibration_detector) isa NoFrameReadoutProducts

    invalid_ramp_calibration = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, gain=3.0,
        sensor=HgCdTeSensor(
            read_time=0.4, sampling_mode=UpTheRampSampling(4)),
        response_model=NullFrameResponse(),
        correction_model=calibration_correction)
    @test_throws InvalidConfiguration begin
        AdaptiveOpticsSim.WavefrontSensors.detector_calibration_signature(
            invalid_ramp_calibration, UInt(0))
    end
    @test_throws InvalidConfiguration begin
        AdaptiveOpticsSim.WavefrontSensors.detector_calibration_frame!(
            invalid_ramp_calibration, calibration_input, 1.0)
    end

    det_hgcdte_row_corrected = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeSensor(),
        response_model=NullFrameResponse(),
        correction_model=ReferenceRowCommonModeCorrection(1))
    row_corrected = capture!(det_hgcdte_row_corrected, row_pattern; rng=MersenneTwister(20))
    @test maximum(abs, row_corrected) < 1e-6
    det_hgcdte_col_corrected = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeSensor(),
        response_model=NullFrameResponse(),
        correction_model=ReferenceColumnCommonModeCorrection(1))
    col_corrected = capture!(det_hgcdte_col_corrected, col_pattern; rng=MersenneTwister(21))
    @test maximum(abs, col_corrected) < 1e-6
    det_hgcdte_output_corrected = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeSensor(),
        response_model=NullFrameResponse(),
        correction_model=ReferenceOutputCommonModeCorrection(2; edge_rows=1, edge_cols=1))
    output_pattern = hcat(fill(5.0, 4, 2), fill(10.0, 4, 2))
    output_corrected = capture!(det_hgcdte_output_corrected, output_pattern; rng=MersenneTwister(22))
    output_meta = detector_export_metadata(det_hgcdte_output_corrected)
    @test output_meta.readout_correction == :reference_output_common_mode
    @test output_meta.correction_group_cols == 2
    @test maximum(abs, output_corrected) < 1e-6
    det_hgcdte_composite = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeSensor(),
        response_model=NullFrameResponse(),
        correction_model=CompositeFrameReadoutCorrection((
            ReferenceRowCommonModeCorrection(1),
            ReferenceColumnCommonModeCorrection(1))))
    composite_pattern = row_pattern .+ col_pattern
    composite_corrected = capture!(det_hgcdte_composite, composite_pattern; rng=MersenneTwister(23))
    composite_meta = detector_export_metadata(det_hgcdte_composite)
    @test composite_meta.readout_correction == :composite
    @test composite_meta.correction_stage_count == 2
    @test maximum(abs, composite_corrected) < 1e-6

    shared_readout = AdaptiveOpticsSim.Detectors.HgCdTeReadout(
        glow_rate=0.25, read_time=0.01,
        sampling_mode=FowlerSampling(2))
    @test HgCdTeSensor(shared_readout).readout === shared_readout
    @test HgCdTeAvalancheArraySensor(shared_readout;
        avalanche_gain=3.0).readout === shared_readout
    @test !supports_avalanche_gain(HgCdTeSensor(shared_readout))
    @test supports_avalanche_gain(HgCdTeAvalancheArraySensor(
        shared_readout; avalanche_gain=3.0))
    @test !hasfield(typeof(shared_readout), :persistence_model)
    shared_readout_persistent = HgCdTeSensor(shared_readout;
        persistence_model=ExponentialPersistence(0.25, 0.0))
    shared_readout_avalanche = HgCdTeAvalancheArraySensor(shared_readout;
        avalanche_gain=3.0)
    @test shared_readout_persistent.persistence_model isa
        ExponentialPersistence
    @test shared_readout_avalanche.persistence_model isa
        AdaptiveOpticsSim.Detectors.NullPersistence

    conventional_saturation = Detector(
        integration_time=1.0, noise=NoiseNone(), qe=1.0,
        full_well=100.0, sensor=HgCdTeSensor())
    @test capture!(conventional_saturation, fill(150.0, 2, 2);
        rng=Xoshiro(2301)) == fill(100.0, 2, 2)

    conventional_nonlinearity = Detector(
        integration_time=1.0, noise=NoiseNone(), qe=1.0,
        response_model=NullFrameResponse(),
        nonlinearity_model=SaturatingFrameNonlinearity(0.1),
        sensor=HgCdTeSensor())
    @test capture!(conventional_nonlinearity, fill(10.0, 2, 2);
        rng=Xoshiro(2302)) == fill(5.0, 2, 2)

    conventional_persistence = Detector(
        integration_time=1.0, noise=NoiseNone(), qe=1.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeSensor(
            persistence_model=ExponentialPersistence(0.5, 0.0)))
    capture!(conventional_persistence, fill(4.0, 2, 2);
        rng=Xoshiro(2303))
    @test capture!(conventional_persistence, zeros(2, 2);
        rng=Xoshiro(2304)) == fill(2.0, 2, 2)

    @test_throws InvalidConfiguration HgCdTeAvalancheArraySensor(avalanche_gain=0.5)
    @test_throws InvalidConfiguration HgCdTeAvalancheArraySensor(excess_noise_factor=0.5)
    @test_throws InvalidConfiguration HgCdTeAvalancheArraySensor(
        avalanche_gain=NaN)
    @test_throws InvalidConfiguration HgCdTeAvalancheArraySensor(
        excess_noise_factor=Inf)
    @test_throws InvalidConfiguration HgCdTeSensor(glow_rate=-1.0)
    @test_throws InvalidConfiguration HgCdTeSensor(read_time=-1.0)
    @test_throws InvalidConfiguration HgCdTeSensor(glow_rate=NaN)
    @test_throws InvalidConfiguration HgCdTeSensor(read_time=Inf)
    @test_throws InvalidConfiguration HgCdTeSensor(
        sampling_mode=AveragedNonDestructiveReads(0))
    @test_throws InvalidConfiguration HgCdTeSensor(
        sampling_mode=FowlerSampling(0))
    @test_throws InvalidConfiguration ReferencePixelCommonModeCorrection(0, 0)
    @test_throws InvalidConfiguration ReferenceRowCommonModeCorrection(0)
    @test_throws InvalidConfiguration ReferenceColumnCommonModeCorrection(0)
    @test_throws InvalidConfiguration ReferenceOutputCommonModeCorrection(0)
    @test_throws InvalidConfiguration CompositeFrameReadoutCorrection(())

end

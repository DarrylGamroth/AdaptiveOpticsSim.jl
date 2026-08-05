@testset "Sampled detector parameter ownership" begin
    canonical_detector = Detector(exposure_duration=1.0)
    @test :exposure_duration in fieldnames(typeof(canonical_detector.params))
    @test :integration_time ∉ fieldnames(typeof(canonical_detector.params))
    @test_throws MethodError Detector(integration_time=1.0)
    @test_throws MethodError RollingShutter(line_time=1.0)
    @test_throws MethodError FrameTransferAcquisition(transfer_time=1.0)
    @test_throws MethodError HgCdTeSensor(read_time=1.0)

    raw_kernel = [0.0 0.1 0.0; 0.1 0.6 0.1; 0.0 0.1 0.0]
    sampled_response = SampledFrameResponse(raw_kernel)
    sampled_snapshot = copy(sampled_response.kernel)
    raw_kernel[2, 2] = 0.2
    @test sampled_response.kernel == sampled_snapshot

    typed_kernel = copy(sampled_snapshot)
    typed_response = SampledFrameResponse{
        Float64,typeof(typed_kernel)}(typed_kernel)
    typed_kernel[2, 2] = 0.2
    @test typed_response.kernel == sampled_snapshot

    gaussian_kernel = [0.25, 0.5, 0.25]
    gaussian_response = GaussianPixelResponse{
        Float64,typeof(gaussian_kernel)}(0.75, gaussian_kernel)
    gaussian_kernel[2] = 0.25
    @test gaussian_response.kernel == [0.25, 0.5, 0.25]

    rectangular_kernel_x = [0.25, 0.5, 0.25]
    rectangular_kernel_y = [0.125, 0.75, 0.125]
    rectangular_response = RectangularPixelAperture{
        Float64,typeof(rectangular_kernel_x),typeof(rectangular_kernel_y)}(
            2.0, 2.0, 0.75, 0.5, rectangular_kernel_x,
            rectangular_kernel_y)
    rectangular_kernel_x[2] = 0.25
    rectangular_kernel_y[2] = 0.5
    @test rectangular_response.kernel_x == [0.25, 0.5, 0.25]
    @test rectangular_response.kernel_y == [0.125, 0.75, 0.125]

    raw_wavelengths = [0.5e-6, 0.6e-6, 0.7e-6]
    raw_qe = [0.2, 0.8, 0.4]
    sampled_qe = SampledQuantumEfficiency(raw_wavelengths, raw_qe)
    wavelength_snapshot = copy(sampled_qe.wavelengths)
    qe_snapshot = copy(sampled_qe.values)
    raw_wavelengths[1] = 0.4e-6
    raw_qe[1] = 0.9
    @test sampled_qe.wavelengths == wavelength_snapshot
    @test sampled_qe.values == qe_snapshot

    typed_wavelengths = copy(wavelength_snapshot)
    typed_values = copy(qe_snapshot)
    typed_qe = AdaptiveOpticsSim.Detectors.SampledQuantumEfficiency{
        Float64,typeof(typed_wavelengths)}(
            typed_wavelengths, typed_values, 0.0)
    typed_wavelengths[1] = 0.4e-6
    typed_values[1] = 0.9
    @test typed_qe.wavelengths == wavelength_snapshot
    @test typed_qe.values == qe_snapshot

    raw_cmos_sigma = [0.1 0.2; 0.3 0.4]
    cmos_noise = CMOSReadNoiseMap(raw_cmos_sigma)
    raw_cmos_sigma[1, 1] = 9.0
    @test cmos_noise.sigma == [0.1 0.2; 0.3 0.4]

    typed_cmos_sigma = [0.1 0.2; 0.3 0.4]
    typed_cmos_noise = CMOSReadNoiseMap{
        Float64,typeof(typed_cmos_sigma)}(typed_cmos_sigma)
    typed_cmos_sigma[1, 1] = 9.0
    @test typed_cmos_noise.sigma == [0.1 0.2; 0.3 0.4]

    raw_output_gains = [1.0, 1.1]
    raw_output_offsets = [0.0, 0.2]
    output_pattern = StaticCMOSOutputPattern(2, raw_output_gains,
        raw_output_offsets)
    typed_output_pattern = StaticCMOSOutputPattern{
        Float64,typeof(raw_output_gains),typeof(raw_output_offsets)}(
            2, raw_output_gains, raw_output_offsets)
    raw_output_gains[1] = 9.0
    raw_output_offsets[1] = 9.0
    @test output_pattern.gains == [1.0, 1.1]
    @test output_pattern.offsets == [0.0, 0.2]
    @test typed_output_pattern.gains == [1.0, 1.1]
    @test typed_output_pattern.offsets == [0.0, 0.2]

    raw_background = [0.1 0.2; 0.3 0.4]
    typed_background = BackgroundFrame{
        Float64,typeof(raw_background)}(raw_background)
    raw_background[1, 1] = 9.0
    @test typed_background.map == [0.1 0.2; 0.3 0.4]

    inferred_background_values = [0.4 0.3; 0.2 0.1]
    inferred_background = BackgroundFrame(inferred_background_values)
    inferred_background_values[1, 1] = 9.0
    @test inferred_background.map == [0.4 0.3; 0.2 0.1]

    sensor_sigma = zeros(2, 2)
    sensor_gains = [1.0, 2.0]
    sensor_offsets = [0.0, 10.0]
    external_cmos_sensor = CMOSSensor(
        column_readout_sigma=0.25,
        row_readout_sigma=0.5,
        readout_noise_model=CMOSReadNoiseMap(sensor_sigma),
        output_model=StaticCMOSOutputPattern(1, sensor_gains,
            sensor_offsets),
        timing_model=RollingShutter(0.01),
    )
    owned_cmos_detector = Detector(noise=NoiseNone(), qe=1.0,
        response_model=NullFrameResponse(), sensor=external_cmos_sensor,
        T=Float32)
    owned_cmos_sensor = owned_cmos_detector.params.sensor
    @test owned_cmos_sensor !== external_cmos_sensor
    @test owned_cmos_sensor.readout_noise_model.sigma !==
        external_cmos_sensor.readout_noise_model.sigma
    @test owned_cmos_sensor.output_model.gains !==
        external_cmos_sensor.output_model.gains
    @test owned_cmos_sensor.output_model.offsets !==
        external_cmos_sensor.output_model.offsets
    @test eltype(owned_cmos_sensor.readout_noise_model.sigma) === Float32
    @test eltype(owned_cmos_sensor.output_model.gains) === Float32
    @test owned_cmos_sensor.column_readout_sigma == 0.25f0
    @test owned_cmos_sensor.row_readout_sigma == 0.5f0
    @test owned_cmos_sensor.timing_model.line_duration == 0.01f0
    owned_cmos_before = copy(capture!(owned_cmos_detector,
        ones(Float32, 2, 2); rng=MersenneTwister(2300)))
    fill!(external_cmos_sensor.readout_noise_model.sigma, 9.0)
    fill!(external_cmos_sensor.output_model.gains, 9.0)
    fill!(external_cmos_sensor.output_model.offsets, 9.0)
    owned_cmos_after = copy(capture!(owned_cmos_detector,
        ones(Float32, 2, 2); rng=MersenneTwister(2300)))
    @test owned_cmos_after == owned_cmos_before
    @test all(iszero, owned_cmos_sensor.readout_noise_model.sigma)
    @test owned_cmos_sensor.output_model.gains == Float32[1, 2]
    @test owned_cmos_sensor.output_model.offsets == Float32[0, 10]

    # Ordinary detector construction accepts an extension response after the
    # extension supplies conversion and validation, but WFS reference caching
    # must not key instance-dependent behavior by type alone.
    unkeyed_response = Detector(noise=NoiseNone(),
        response_model=UnkeyedCalibrationFrameResponse(0.25))
    @test_throws InvalidConfiguration begin
        AdaptiveOpticsSim.WavefrontSensors.detector_calibration_signature(
            unkeyed_response, UInt(0))
    end

    detector = Detector(exposure_duration=1.0, noise=NoiseNone(),
        qe=sampled_qe, response_model=sampled_response)
    @test detector.params.response_model.kernel !== sampled_response.kernel
    @test detector.params.quantum_efficiency_model.wavelengths !==
        sampled_qe.wavelengths
    @test detector.params.quantum_efficiency_model.values !==
        sampled_qe.values
    src = Source(band=:custom, wavelength=0.6e-6,
        photon_irradiance=1.0)
    impulse = zeros(5, 5)
    impulse[3, 3] = 1.0
    before = copy(capture!(detector, impulse, src;
        rng=MersenneTwister(2301)))
    signature_before = AdaptiveOpticsSim.WavefrontSensors.detector_calibration_signature(
        detector, UInt(7))

    sampled_response.kernel[2, 2] = 0.2
    sampled_qe.values[2] = 0.1
    after = copy(capture!(detector, impulse, src;
        rng=MersenneTwister(2301)))
    @test after == before
    @test AdaptiveOpticsSim.WavefrontSensors.detector_calibration_signature(
        detector, UInt(7)) == signature_before

    raw_gain_map = [1.0 0.8; 0.9 1.1]
    raw_dark_map = [0.1 0.2; 0.3 0.4]
    raw_bad_mask = Bool[false true; false false]
    @test_throws InvalidConfiguration PixelResponseNonuniformity(
        [1.0 NaN; 0.9 1.1])
    @test_throws InvalidConfiguration PixelResponseNonuniformity(
        [1.0 -0.1; 0.9 1.1])
    @test_throws InvalidConfiguration DarkSignalNonuniformity(
        [0.1 Inf; 0.3 0.4])
    @test_throws InvalidConfiguration DarkSignalNonuniformity(
        [0.1 -0.2; 0.3 0.4])
    prnu = PixelResponseNonuniformity(raw_gain_map)
    dsnu = DarkSignalNonuniformity(raw_dark_map)
    bad_pixels = BadPixelMask(raw_bad_mask; throughput=0.0)
    raw_gain_map[1, 1] = 0.1
    raw_dark_map[1, 1] = 0.9
    raw_bad_mask[1, 1] = true
    @test prnu.gain_map == [1.0 0.8; 0.9 1.1]
    @test dsnu.dark_map == [0.1 0.2; 0.3 0.4]
    @test bad_pixels.mask == Bool[false true; false false]

    typed_gain_map = [1.0 0.8; 0.9 1.1]
    typed_dark_map = [0.1 0.2; 0.3 0.4]
    typed_bad_mask = Bool[false true; false false]
    typed_prnu = PixelResponseNonuniformity{
        Float64,typeof(typed_gain_map)}(typed_gain_map)
    typed_dsnu = DarkSignalNonuniformity{
        Float64,typeof(typed_dark_map)}(typed_dark_map)
    typed_bad_pixels = BadPixelMask{
        Float64,typeof(typed_bad_mask)}(typed_bad_mask, 0.0)
    typed_gain_map[1, 1] = 0.2
    typed_dark_map[1, 1] = 0.9
    typed_bad_mask[1, 1] = true
    @test typed_prnu.gain_map[1, 1] == 1.0
    @test typed_dsnu.dark_map[1, 1] == 0.1
    @test !typed_bad_pixels.mask[1, 1]

    composite_defects = CompositeDetectorDefectModel(prnu, bad_pixels)
    @test composite_defects.stages[1].gain_map !== prnu.gain_map
    @test composite_defects.stages[2].mask !== bad_pixels.mask

    defect_detector = Detector(noise=NoiseNone(), sensor=CMOSSensor(),
        defect_model=composite_defects)
    defect_signature = AdaptiveOpticsSim.WavefrontSensors.detector_calibration_signature(
        defect_detector, UInt(8))
    @test defect_detector.params.defect_model.stages[1].gain_map !==
        prnu.gain_map
    @test defect_detector.params.defect_model.stages[2].mask !==
        bad_pixels.mask
    prnu.gain_map[1, 1] = 0.2
    bad_pixels.mask[2, 2] = true
    @test composite_defects.stages[1].gain_map[1, 1] == 1.0
    @test !composite_defects.stages[2].mask[2, 2]
    @test AdaptiveOpticsSim.WavefrontSensors.detector_calibration_signature(
        defect_detector, UInt(8)) == defect_signature

    replacement_detector = Detector(noise=NoiseNone(), sensor=CMOSSensor(),
        defect_model=CompositeDetectorDefectModel(
            PixelResponseNonuniformity([0.7 0.8; 0.9 1.1]),
            BadPixelMask(Bool[false true; true false]; throughput=0.0)))
    @test AdaptiveOpticsSim.WavefrontSensors.detector_calibration_signature(
        replacement_detector, UInt(8)) != defect_signature
end

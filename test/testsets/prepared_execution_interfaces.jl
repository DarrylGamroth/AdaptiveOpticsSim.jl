function pe01_detector_rate_map(values::AbstractMatrix{T}) where {
    T<:AbstractFloat,
}
    metadata = OpticalPlaneMetadata(DetectorPlane(), values;
        coordinate_domain=AngularCoordinates(),
        sampling=(one(T), one(T)),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition(),
        spectral=MonochromaticChannel(T(0.55e-6)))
    return IntensityMap(metadata, values)
end

@testset "PE-01 calibration command-plan interface" begin
    actuator_commands = Matrix{Float64}(I, 3, 3)
    assert_calibration_command_plan_conformance(
        Calibration.ActuatorCalibrationCommands(3),
        actuator_commands,
        0.25,
    )

    matrix_commands = [1.0 3.0; 2.0 4.0; 5.0 6.0]
    assert_calibration_command_plan_conformance(
        Calibration.MatrixCalibrationCommands(matrix_commands),
        matrix_commands,
        0.5,
    )
end

@testset "PE-01 four-pupil LGS plan interface" begin
    target = HostComputeDevice()
    intensity = ones(4, 4)
    scratch = similar(intensity)
    fft_buffer = zeros(ComplexF64, size(intensity))
    ifft_buffer = similar(fft_buffer)
    fft_plan = Backends.plan_fft_backend!(fft_buffer)
    ifft_plan = Backends.plan_ifft_backend!(ifft_buffer)

    no_lgs = WavefrontSensors.NoPreparedFourPupilLGS()
    baseline = copy(intensity)
    assert_prepared_four_pupil_lgs_conformance!(no_lgs, intensity,
        scratch, fft_buffer, fft_plan, ifft_buffer, ifft_plan, target)
    @test intensity == baseline

    elongation = WavefrontSensors.PreparedFourPupilElongation(
        [0.25, 0.5, 0.25], 1)
    fill!(intensity, 1.0)
    assert_prepared_four_pupil_lgs_conformance!(elongation, intensity,
        scratch, fft_buffer, fft_plan, ifft_buffer, ifft_plan, target)
    @test intensity == ones(4, 4)

    sodium = WavefrontSensors.PreparedFourPupilSodiumLayerProfile(
        ones(ComplexF64, size(intensity)))
    intensity .= reshape(collect(1.0:16.0), 4, 4)
    baseline = copy(intensity)
    assert_prepared_four_pupil_lgs_conformance!(sodium, intensity,
        scratch, fft_buffer, fft_plan, ifft_buffer, ifft_plan, target)
    @test intensity ≈ baseline atol=1e-12 rtol=1e-12

    unavailable = AcceleratorComputeDevice(CUDABackend(), 0)
    @test WavefrontSensors._require_exact_prepared_four_pupil_lgs_target(
        no_lgs, unavailable) === nothing
    @test_throws WFSPreparationError begin
        WavefrontSensors._require_exact_prepared_four_pupil_lgs_target(
            elongation, unavailable)
    end
    @test_throws WFSPreparationError begin
        WavefrontSensors._require_exact_prepared_four_pupil_lgs_target(
            sodium, unavailable)
    end
end

@testset "PE-01 prepared device-context interface" begin
    target = HostComputeDevice()
    context = Backends._prepare_device_execution_context(target)
    observe_host() = HostComputeDevice()
    @test @inferred(Backends._with_prepared_device_execution_context(
        observe_host, context)) == target
    assert_prepared_device_execution_context_conformance(
        context, target, observe_host, target)

    @test_throws ComputeDeviceError begin
        Backends._prepare_device_execution_context(
            AcceleratorComputeDevice(CUDABackend(), 0))
    end
end

@testset "PE-01 acquisition-lifecycle interfaces" begin
    target = HostComputeDevice()
    definition = GlobalShutterAcquisitionDefinition(
        PlantDuration(1_000_000_000))
    map = pe01_detector_rate_map(ones(2, 2))
    detector = Detector(exposure_duration=1.0, noise=NoiseNone())
    foreign_detector = Detector(exposure_duration=1.0, noise=NoiseNone())
    prepared = prepare_global_shutter_acquisition(
        detector, map, definition)
    foreign_prepared = prepare_global_shutter_acquisition(
        foreign_detector, map, definition)
    assert_prepared_acquisition_lifecycle_conformance(
        prepared,
        GlobalShutterAcquisitionState(prepared),
        GlobalShutterAcquisitionState(foreign_prepared),
        target,
        DetectorAcquisitionError;
        detector_backed=true,
    )

    rolling_definition = RollingShutterAcquisitionDefinition(
        PlantDuration(1_000_000_000))
    rolling_sensor = CMOSSensor(
        timing_model=RollingShutter(0.01; row_group_size=1))
    foreign_rolling_sensor = CMOSSensor(
        timing_model=RollingShutter(0.01; row_group_size=1))
    rolling_detector = Detector(exposure_duration=1.0, noise=NoiseNone(),
        sensor=rolling_sensor)
    foreign_rolling_detector = Detector(exposure_duration=1.0,
        noise=NoiseNone(), sensor=foreign_rolling_sensor)
    rolling = prepare_rolling_shutter_acquisition(
        rolling_detector, map, rolling_definition)
    foreign_rolling = prepare_rolling_shutter_acquisition(
        foreign_rolling_detector, map, rolling_definition)
    assert_prepared_acquisition_lifecycle_conformance(
        rolling,
        RollingShutterAcquisitionState(rolling),
        RollingShutterAcquisitionState(foreign_rolling),
        target,
        DetectorAcquisitionError;
        detector_backed=true,
    )

    transfer_definition = FrameTransferAcquisitionDefinition(
        PlantDuration(1_000_000_000))
    transfer_sensor = EMCCDSensor(excess_noise_factor=1.0,
        acquisition_mode=FrameTransferAcquisition(transfer_duration=0.1))
    foreign_transfer_sensor = EMCCDSensor(excess_noise_factor=1.0,
        acquisition_mode=FrameTransferAcquisition(transfer_duration=0.1))
    transfer_detector = Detector(exposure_duration=1.0, noise=NoiseNone(),
        sensor=transfer_sensor)
    foreign_transfer_detector = Detector(exposure_duration=1.0,
        noise=NoiseNone(), sensor=foreign_transfer_sensor)
    transfer = prepare_frame_transfer_acquisition(
        transfer_detector, map, transfer_definition)
    foreign_transfer = prepare_frame_transfer_acquisition(
        foreign_transfer_detector, map, transfer_definition)
    assert_prepared_acquisition_lifecycle_conformance(
        transfer,
        FrameTransferAcquisitionState(transfer),
        FrameTransferAcquisitionState(foreign_transfer),
        target,
        DetectorAcquisitionError;
        detector_backed=true,
    )

    direct_definition = DirectMeasurementAcquisitionDefinition(
        PlantDuration(1_000_000_000))
    measurement = WFSMeasurement(zeros(2);
        units=:metre, kind=:modal_residual)
    foreign_measurement = WFSMeasurement(zeros(2);
        units=:metre, kind=:modal_residual)
    direct = prepare_direct_measurement_acquisition(
        measurement, direct_definition)
    foreign_direct = prepare_direct_measurement_acquisition(
        foreign_measurement, direct_definition)
    assert_prepared_acquisition_lifecycle_conformance(
        direct,
        DirectMeasurementAcquisitionState(direct),
        DirectMeasurementAcquisitionState(foreign_direct),
        target,
        PlantScheduleError;
        detector_backed=false,
    )
end

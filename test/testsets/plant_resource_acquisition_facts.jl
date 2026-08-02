@inline acquisition_resource_array_bytes(array) =
    UInt64(length(array) * sizeof(eltype(array)))

struct AcquisitionResourceUnknownQuantumEfficiency <:
    AbstractQuantumEfficiencyModel
    values::Vector{Float64}
end

struct AcquisitionResourceFakeBackend <: AbstractArrayBackend end
struct AcquisitionResourceUnsupportedThermalState <:
    AbstractDetectorThermalState end
struct AcquisitionResourceUnsupportedReadout <: FrameReadoutProducts end
struct AcquisitionResourceUnsupportedLifecycle <:
    AbstractPreparedAcquisitionLifecycle end
struct AcquisitionResourceUnsupportedLifecycleState <:
    AbstractAcquisitionLifecycleState end

@inline function acquisition_resource_sum_bytes(arrays)
    total = UInt64(0)
    for array in arrays
        total += acquisition_resource_array_bytes(array)
    end
    return total
end

function acquisition_resource_intensity_map(
    values::AbstractMatrix{T}) where {T<:AbstractFloat}
    metadata = OpticalPlaneMetadata(DetectorPlane(), values;
        coordinate_domain=AngularCoordinates(),
        sampling=(one(T), one(T)),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition(),
        spectral=MonochromaticChannel(T(0.55e-6)))
    return IntensityMap(metadata, values)
end

@testset "Acquisition product structural resource facts" begin
    target = HostComputeDevice()
    id = StructuralResourceOwnerID(:acquisition, :science)
    observation = zeros(Float32, 3, 4)
    products = AcquisitionProducts(observation;
        metadata=(kind=:detector_frame, units=:electron))
    fact = structural_resource_fact(products, id, target)

    @test structural_resource_known(fact)
    @test structural_resident_bytes(fact) ==
        acquisition_resource_array_bytes(observation)
    @test iszero(structural_workspace_bytes(fact))

    measurement = WFSMeasurement(zeros(Float64, 5);
        units=:metre, kind=:modal_residual)
    direct_products = AcquisitionProducts(nothing, measurement;
        metadata=(kind=:modal_residual, units=:metre))
    direct_fact = structural_resource_fact(direct_products, id, target)
    @test structural_resident_bytes(direct_fact) ==
        acquisition_resource_array_bytes(measurement.storage)

    wfs_observation = WFSObservation(zeros(Float32, 2, 3);
        units=:detected_electrons, layout=:lenslet_mosaic)
    observation_products = AcquisitionProducts(wfs_observation, nothing;
        metadata=(kind=:wfs_observation,))
    observation_fact = structural_resource_fact(
        observation_products, id, target)
    @test structural_resource_known(observation_fact)
    @test structural_resident_bytes(observation_fact) ==
        acquisition_resource_array_bytes(wfs_observation.storage)

    unsupported = AcquisitionProducts(Ref(1.0);
        metadata=(kind=:unsupported,))
    unsupported_fact = structural_resource_fact(unsupported, id, target)
    @test !structural_resource_known(unsupported_fact)
    @test structural_resource_unknown_reason(unsupported_fact) ==
        :unsupported_acquisition_products
end


@testset "Detector parameter-model structural dispatch" begin
    target = HostComputeDevice()
    fixed_temperature = FixedTemperature(temperature_K=80.0)
    dynamic_temperature = FirstOrderThermalModel(
        ambient_temperature_K=290.0,
        setpoint_temperature_K=80.0,
        time_constant_s=10.0,
    )
    gaussian_response = GaussianPixelResponse(response_width_px=0.8)
    aperture_response = RectangularPixelAperture(
        pitch_x_px=1.0,
        pitch_y_px=1.0,
        fill_factor_x=0.9,
        fill_factor_y=0.8,
    )
    coupling = InterpixelCapacitance(
        [0.0 0.1 0.0; 0.1 0.6 0.1; 0.0 0.1 0.0])
    dark_nonuniformity = DarkSignalNonuniformity(ones(2, 2))
    correction = CompositeFrameReadoutCorrection(
        ReferenceRowCommonModeCorrection(1),
        ReferenceColumnCommonModeCorrection(1),
    )
    sampled_qe = AdaptiveOpticsSim.Detectors.SampledQuantumEfficiency(
        [0.5e-6, 0.6e-6], [0.8, 0.9])
    background = BackgroundFrame(ones(2, 2))
    read_noise = CMOSReadNoiseMap(ones(2, 2))
    output_pattern = AdaptiveOpticsSim.Detectors.StaticCMOSOutputPattern(
        2, [1.0], [0.0])

    cases = (
        (fixed_temperature, ()),
        (dynamic_temperature, ()),
        (gaussian_response, (gaussian_response.kernel,)),
        (aperture_response,
            (aperture_response.kernel_x, aperture_response.kernel_y)),
        (coupling, (coupling.response.kernel,)),
        (dark_nonuniformity, (dark_nonuniformity.dark_map,)),
        (correction, ()),
        (sampled_qe, (sampled_qe.wavelengths, sampled_qe.values)),
        (background, (background.map,)),
        (read_noise, (read_noise.sigma,)),
        (output_pattern, (output_pattern.gains, output_pattern.offsets)),
        (CCDSensor(sampling_mode=SkipperSampling(3)), ()),
        (InGaAsSensor(), ()),
        (HgCdTeAvalancheArraySensor(), ()),
    )

    for (model, arrays) in cases
        bytes = Plant._detector_parameter_bytes(model, target)
        expected = acquisition_resource_sum_bytes(arrays)
        @test bytes.known
        @test bytes.storage.present == !iszero(expected)
        @test bytes.storage.bytes == expected
    end
end


@testset "Detector readout structural dispatch" begin
    target = HostComputeDevice()
    wrong_target = AcceleratorComputeDevice(
        AcquisitionResourceFakeBackend(), UInt32(1))
    window = FrameWindow(1:2, 1:2)

    memory_bytes = Plant._acquisition_optional_array_bytes(
        Memory{Float64}(undef, 2), wrong_target, :resident_bytes)
    @test !memory_bytes.present
    @test iszero(memory_bytes.bytes)

    mean_frame = zeros(Float64, 2, 2)
    skipper = SkipperReadoutProducts(mean_frame, 3)
    skipper_storage = Plant._detector_readout_storage(
        skipper, nothing, target)
    @test skipper_storage.known
    @test skipper_storage.resident.bytes ==
        acquisition_resource_array_bytes(mean_frame)
    @test iszero(skipper_storage.workspace.bytes)

    reference_frame = zeros(Float64, 2, 2)
    signal_frame = zeros(Float64, 2, 2)
    sampled_cube = zeros(Float64, 2, 2, 3)
    sampled = SampledFrameReadoutProducts(
        reference_frame, signal_frame, sampled_cube)
    sampled_storage = Plant._detector_readout_storage(
        sampled, nothing, target)
    @test sampled_storage.known
    @test sampled_storage.resident.bytes ==
        acquisition_resource_sum_bytes(
            (reference_frame, signal_frame, sampled_cube))

    combined_frame = zeros(Float64, 2, 2)
    reference_cube = zeros(Float64, 2, 2, 2)
    signal_cube = zeros(Float64, 2, 2, 2)
    read_cube = zeros(Float64, 2, 2, 4)
    read_times = zeros(Float64, 4)
    workspace_reference_cube = zeros(Float64, 2, 2, 2)
    workspace_signal_cube = zeros(Float64, 2, 2, 2)
    multi_read = MultiReadFrameReadoutProducts(
        reference_frame,
        signal_frame,
        combined_frame,
        reference_cube,
        signal_cube,
        read_cube,
        read_times,
        workspace_reference_cube,
        workspace_signal_cube,
    )
    multi_read_resident = acquisition_resource_sum_bytes((
        reference_frame,
        signal_frame,
        combined_frame,
        reference_cube,
        signal_cube,
        read_cube,
        read_times,
    ))
    unwindowed = Plant._detector_readout_storage(
        multi_read, nothing, target)
    @test unwindowed.known
    @test unwindowed.resident.bytes == multi_read_resident
    @test iszero(unwindowed.workspace.bytes)
    windowed = Plant._detector_readout_storage(
        multi_read, window, target)
    @test windowed.resident.bytes == multi_read_resident
    @test windowed.workspace.bytes == acquisition_resource_sum_bytes((
        workspace_reference_cube, workspace_signal_cube))

    slope_frame = zeros(Float64, 2, 2)
    intercept_frame = zeros(Float64, 2, 2)
    integrated_frame = zeros(Float64, 2, 2)
    ramp_cube = zeros(Float64, 2, 2, 3)
    ramp_times = zeros(Float64, 3)
    workspace_slope = zeros(Float64, 2, 2)
    workspace_intercept = zeros(Float64, 2, 2)
    workspace_integrated = zeros(Float64, 2, 2)
    workspace_cube = zeros(Float64, 2, 2, 3)
    ramp = UpTheRampReadoutProducts(
        slope_frame,
        intercept_frame,
        integrated_frame,
        ramp_cube,
        ramp_times,
        workspace_slope,
        workspace_intercept,
        workspace_integrated,
        workspace_cube,
        SynthesizedFinalChargeRamp,
    )
    ramp_storage = Plant._detector_readout_storage(ramp, window, target)
    @test ramp_storage.known
    @test ramp_storage.workspace.bytes == acquisition_resource_sum_bytes((
        workspace_slope,
        workspace_intercept,
        workspace_integrated,
        workspace_cube,
    ))

    @test Plant._detector_thermal_state_storage(NoThermalState())
    @test Plant._detector_thermal_state_storage(DetectorThermalState(80.0))
    @test !Plant._detector_thermal_state_storage(
        AcquisitionResourceUnsupportedThermalState())
    @test !Plant._detector_readout_storage(
        AcquisitionResourceUnsupportedReadout(), nothing, target).known

    detector = Detector(
        integration_time=0.1,
        noise=NoiseNone(),
        qe=1.0,
        response_model=NullFrameResponse(),
    )
    standalone_multi_read = Plant._standalone_detector_state_resource_fact(
        detector.state, multi_read,
        StructuralResourceOwnerID(:acquisition, :multi_read), target)
    @test !structural_resource_known(standalone_multi_read)
    @test structural_resource_unknown_reason(standalone_multi_read) ==
        :detector_readout_window_required

    unsupported_lifecycle = structural_resource_fact(
        AcquisitionResourceUnsupportedLifecycle(),
        AcquisitionResourceUnsupportedLifecycleState(),
        StructuralResourceOwnerID(:acquisition, :unsupported),
        target,
    )
    @test !structural_resource_known(unsupported_lifecycle)
    @test structural_resource_unknown_reason(unsupported_lifecycle) ==
        :unsupported_acquisition_lifecycle
end

@testset "Conventional detector structural resource facts" begin
    target = HostComputeDevice()
    id = StructuralResourceOwnerID(:acquisition, :detector)
    response = SampledFrameResponse(
        Float64[0 0.1 0; 0.1 0.6 0.1; 0 0.1 0])
    defect = CompositeDetectorDefectModel(
        PixelResponseNonuniformity(fill(0.9, 4, 4)),
        BadPixelMask(falses(4, 4)))
    detector = Detector(
        integration_time=0.25,
        noise=NoiseNone(),
        qe=1.0,
        response_model=response,
        defect_model=defect,
        sensor=CMOSSensor(timing_model=GlobalShutter()),
    )
    map = acquisition_resource_intensity_map(fill(2.0, 4, 4))
    prepared = prepare_global_shutter_acquisition(
        detector, map,
        GlobalShutterAcquisitionDefinition(PlantDuration(250_000_000)))
    state = GlobalShutterAcquisitionState(prepared)
    fact = structural_resource_fact(prepared, state, id, target)

    detector_resident = acquisition_resource_sum_bytes((
        detector.state.frame,
        detector.state.accum_buffer,
        detector.state.latent_buffer,
        detector.params.response_model.kernel,
        detector.params.defect_model.stages[1].gain_map,
        detector.params.defect_model.stages[2].mask,
    ))
    detector_workspace = acquisition_resource_sum_bytes((
        detector.state.presampling_buffer,
        detector.state.presampling_scratch,
        detector.state.response_buffer,
        detector.state.bin_buffer,
        detector.state.temporal_buffer,
        detector.state.noise_buffer,
        detector.state.noise_buffer_host,
        detector.state.batched_buffer_host,
    ))
    @test structural_resource_known(fact)
    @test structural_resident_bytes(fact) == detector_resident
    @test structural_workspace_bytes(fact) == detector_workspace

    params = detector.params
    unsupported_params = DetectorParams(
        params.integration_time,
        params.qe,
        params.psf_sampling,
        params.binning,
        params.gain,
        params.dark_current,
        params.bits,
        params.full_well,
        params.sensor,
        AcquisitionResourceUnknownQuantumEfficiency([0.8]),
        params.response_model,
        params.charge_coupling_model,
        params.defect_model,
        params.timing_model,
        params.correction_model,
        params.nonlinearity_model,
        params.thermal_model,
        params.readout_window,
        params.output_type,
    )
    unsupported_detector = Detector{
        typeof(detector.noise),
        typeof(unsupported_params),
        typeof(detector.state),
        typeof(detector.background_flux),
        typeof(detector.background_map),
        CPUBackend,
    }(
        detector.noise,
        unsupported_params,
        detector.state,
        detector.background_flux,
        detector.background_map,
    )
    unsupported_fact = structural_resource_fact(
        unsupported_detector, id, target)
    @test !structural_resource_known(unsupported_fact)
    @test structural_resource_unknown_reason(unsupported_fact) ==
        :unsupported_detector_parameter_model
end

@testset "Acquisition lifecycle structural resource facts" begin
    target = HostComputeDevice()
    id = StructuralResourceOwnerID(:acquisition, :lifecycle)

    measurement = WFSMeasurement(zeros(Float64, 6);
        units=:metre, kind=:modal_residual)
    prepared_direct = prepare_direct_measurement_acquisition(
        measurement,
        DirectMeasurementAcquisitionDefinition(PlantDuration(1_000_000)))
    direct_state = DirectMeasurementAcquisitionState(prepared_direct)
    direct_fact = structural_resource_fact(
        prepared_direct, direct_state, id, target)
    @test structural_resident_bytes(direct_fact) ==
        acquisition_resource_array_bytes(prepared_direct.instantaneous_sample)
    @test structural_workspace_bytes(direct_fact) ==
        acquisition_resource_array_bytes(prepared_direct.integrated_sample)

    detector = Detector(
        integration_time=0.1,
        noise=NoiseNone(),
        qe=1.0,
        response_model=NullFrameResponse(),
        sensor=EMCCDSensor(
            acquisition_mode=FrameTransferAcquisition(0.01)),
    )
    map = acquisition_resource_intensity_map(fill(1.0, 3, 3))
    prepared_transfer = prepare_frame_transfer_acquisition(
        detector, map,
        FrameTransferAcquisitionDefinition(PlantDuration(100_000_000);
            readout_duration=PlantDuration(20_000_000)))
    transfer_state = FrameTransferAcquisitionState(prepared_transfer)
    transfer_fact = structural_resource_fact(
        prepared_transfer, transfer_state, id, target)
    detector_fact = structural_resource_fact(detector, id, target)
    @test structural_resident_bytes(transfer_fact) ==
        structural_resident_bytes(detector_fact) +
        acquisition_resource_array_bytes(prepared_transfer.storage_frame)
    @test structural_workspace_bytes(transfer_fact) ==
        structural_workspace_bytes(detector_fact)

    rolling_detector = Detector(
        integration_time=0.1,
        noise=NoiseNone(),
        qe=1.0,
        response_model=NullFrameResponse(),
        sensor=CMOSSensor(timing_model=RollingShutter(0.01)),
    )
    rolling_map = acquisition_resource_intensity_map(fill(1.0, 3, 3))
    prepared_rolling = prepare_rolling_shutter_acquisition(
        rolling_detector, rolling_map,
        RollingShutterAcquisitionDefinition(PlantDuration(100_000_000)))
    rolling_state = RollingShutterAcquisitionState(prepared_rolling)
    rolling_fact = structural_resource_fact(
        prepared_rolling, rolling_state, id, target)
    rolling_detector_fact = structural_resource_fact(
        rolling_detector, id, target)
    @test structural_resident_bytes(rolling_fact) ==
        structural_resident_bytes(rolling_detector_fact)
    @test structural_workspace_bytes(rolling_fact) ==
        structural_workspace_bytes(rolling_detector_fact)

    ramp_detector = Detector(
        integration_time=1.0,
        noise=NoiseNone(),
        qe=1.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeSensor(
            read_time=0.1,
            sampling_mode=UpTheRampSampling(3)),
    )
    ramp_map = acquisition_resource_intensity_map(fill(1.0, 2, 2))
    prepared_ramp = prepare_global_shutter_acquisition(
        ramp_detector, ramp_map,
        GlobalShutterAcquisitionDefinition(PlantDuration(1_000_000_000);
            readout_duration=PlantDuration(100_000_000)))
    ramp_state = GlobalShutterAcquisitionState(prepared_ramp)
    ramp_fact = structural_resource_fact(
        prepared_ramp, ramp_state, id, target)
    ramp_products = prepared_ramp.readout_products
    ramp_resident = acquisition_resource_sum_bytes((
        ramp_detector.state.frame,
        ramp_detector.state.accum_buffer,
        ramp_detector.state.latent_buffer,
        ramp_products.slope_frame,
        ramp_products.intercept_frame,
        ramp_products.integrated_frame,
        ramp_products.read_cube,
        ramp_products.read_times,
        prepared_ramp.read_offsets,
    ))
    ramp_workspace = acquisition_resource_sum_bytes((
        ramp_detector.state.presampling_buffer,
        ramp_detector.state.presampling_scratch,
        ramp_detector.state.response_buffer,
        ramp_detector.state.bin_buffer,
        ramp_detector.state.temporal_buffer,
        ramp_detector.state.noise_buffer,
        ramp_detector.state.noise_buffer_host,
        ramp_detector.state.batched_buffer_host,
    ))
    @test structural_resident_bytes(ramp_fact) == ramp_resident
    @test structural_workspace_bytes(ramp_fact) == ramp_workspace
end

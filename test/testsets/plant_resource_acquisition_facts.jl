@inline acquisition_resource_array_bytes(array) =
    UInt64(length(array) * sizeof(eltype(array)))

struct AcquisitionResourceUnknownQuantumEfficiency <:
    AbstractQuantumEfficiencyModel
    values::Vector{Float64}
end

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

    unsupported = AcquisitionProducts(Ref(1.0);
        metadata=(kind=:unsupported,))
    unsupported_fact = structural_resource_fact(unsupported, id, target)
    @test !structural_resource_known(unsupported_fact)
    @test structural_resource_unknown_reason(unsupported_fact) ==
        :unsupported_acquisition_products
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

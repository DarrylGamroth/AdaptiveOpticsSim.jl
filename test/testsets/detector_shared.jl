@testset "Shared frame-detector behavior" begin
    cadence_free_tel = Telescope(resolution=8, diameter=1.0,
        central_obstruction=0.0)
    @test !hasfield(typeof(cadence_free_tel.params), :sampling_time)

    psf = fill(1.0, 8, 8)
    det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=2)
    frame = capture!(det, psf; rng=MersenneTwister(2))
    @test size(frame) == (4, 4)
    @test sum(frame) == sum(psf)

    det_sampling = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, psf_sampling=2, binning=2)
    frame_sampling = capture!(det_sampling, psf; rng=MersenneTwister(2))
    @test size(frame_sampling) == (2, 2)
    @test sum(frame_sampling) == sum(psf)

    rate_values = reshape(Float64.(1:16), 4, 4)
    shared_rate = detector_test_intensity_map(rate_values)
    shared_rate_before = copy(shared_rate.values)
    short_exposure = Detector(integration_time=0.25, noise=NoiseNone(),
        qe=0.5, response_model=NullFrameResponse())
    long_exposure = Detector(integration_time=1.5, noise=NoiseNone(),
        qe=0.5, response_model=NullFrameResponse())
    short_plan = prepare_detector_acquisition(short_exposure, shared_rate)
    long_plan = prepare_detector_acquisition(long_exposure, shared_rate)
    host_target = AdaptiveOpticsSim.Backends.HostComputeDevice()
    @test AdaptiveOpticsSim.Detectors._require_exact_detector_target(
        short_exposure, host_target) === short_exposure
    @test AdaptiveOpticsSim.Detectors._require_exact_detector_acquisition_target(
        short_exposure, short_plan, host_target) === short_plan
    @test_throws InvalidConfiguration begin
        AdaptiveOpticsSim.Detectors._require_exact_detector_acquisition_target(
            short_exposure, long_plan, host_target)
    end
    wrong_target = AdaptiveOpticsSim.Backends.AcceleratorComputeDevice(
        CUDABackend(), 0)
    wrong_target_error = try
        AdaptiveOpticsSim.Detectors._require_exact_detector_acquisition_target(
            short_exposure, short_plan, wrong_target)
        nothing
    catch error
        error
    end
    @test wrong_target_error isa AdaptiveOpticsSim.Backends.ComputeDeviceError
    @test wrong_target_error.operation ===
        :validate_detector_acquisition_target

    for counting_detector in (
        SPADArrayDetector((4, 4); noise=NoiseNone()),
        MKIDArrayDetector(noise=NoiseNone()),
    )
        @test AdaptiveOpticsSim.Detectors._require_exact_counting_detector_target(
            counting_detector, host_target) === counting_detector
        counting_target_error = try
            AdaptiveOpticsSim.Detectors._require_exact_counting_detector_target(
                counting_detector, wrong_target)
            nothing
        catch error
            error
        end
        @test counting_target_error isa
            AdaptiveOpticsSim.Backends.ComputeDeviceError
        @test counting_target_error.operation ===
            :validate_detector_acquisition_target
    end
    @test wrong_target_error.reason === :wrong_device
    @test wrong_target_error.device == wrong_target
    @test_throws InvalidConfiguration begin
        AdaptiveOpticsSim.Detectors._require_exact_detector_readout_products_target(
            UnsupportedExactTargetReadoutProducts(), host_target)
    end
    sampled_products = AdaptiveOpticsSim.Detectors.SampledFrameReadoutProducts(
        nothing, zeros(2, 2), zeros(2, 2, 2))
    @test AdaptiveOpticsSim.Detectors._require_exact_detector_readout_products_target(
        sampled_products, host_target) === sampled_products
    rich_kernel = [0.0 0.1 0.0; 0.1 0.6 0.1; 0.0 0.1 0.0]
    rich_detector = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=SampledQuantumEfficiency([0.5e-6, 0.6e-6], [0.4, 0.8]),
        sensor=CMOSSensor(
            readout_noise_model=CMOSReadNoiseMap(ones(4, 4)),
            output_model=StaticCMOSOutputPattern(
                2, [1.0, 1.0], [0.0, 0.0])),
        response_model=SampledFrameResponse(rich_kernel),
        charge_coupling_model=InterpixelCapacitance(rich_kernel),
        defect_model=CompositeDetectorDefectModel(
            PixelResponseNonuniformity(ones(4, 4)),
            DarkSignalNonuniformity(zeros(4, 4)),
            BadPixelMask(falses(4, 4))),
        background_flux=fill(0.1, 4, 4),
        background_map=fill(0.05, 4, 4))
    rich_plan = prepare_detector_acquisition(rich_detector, shared_rate)
    @test AdaptiveOpticsSim.Detectors._require_exact_detector_acquisition_target(
        rich_detector, rich_plan, host_target) === rich_plan
    @test !applicable(DetectorAcquisitionPlan, short_exposure.params,
        shared_rate.metadata, 1.0, 1.0)
    @test !applicable(typeof(short_plan),
        short_exposure.params, shared_rate.metadata, 1.0, 1.0)
    short_frame = copy(capture!(short_exposure, shared_rate, short_plan;
        rng=MersenneTwister(200)))
    long_frame = copy(capture!(long_exposure, shared_rate, long_plan;
        rng=MersenneTwister(200)))
    @test short_frame == rate_values .* 0.125
    @test long_frame == rate_values .* 0.75
    @test long_frame == short_frame .* 6
    @test shared_rate.values == shared_rate_before

    for invalid_values in (
            [-1.0 1.0; 1.0 1.0],
            [NaN 1.0; 1.0 1.0],
            [Inf 1.0; 1.0 1.0],
            [-Inf 1.0; 1.0 1.0],
            Matrix{Float64}(undef, 0, 0),
        )
        @test_throws InvalidConfiguration prepare_detector_acquisition(
            short_exposure, detector_test_intensity_map(invalid_values))
    end

    trusted_values = ones(2, 2)
    trusted_rate = detector_test_intensity_map(trusted_values)
    trusted_detector = Detector(integration_time=1.0, noise=NoiseNone(),
        response_model=NullFrameResponse())
    trusted_plan = prepare_detector_acquisition(trusted_detector,
        trusted_rate)
    trusted_values[1, 1] = -1.0
    @test capture!(trusted_detector, trusted_rate, trusted_plan;
        rng=MersenneTwister(201))[1, 1] == -1.0

    spectral_rate = detector_test_intensity_map(ones(2, 2);
        spectral=MonochromaticChannel(0.55e-6))
    spectral_detector = Detector(integration_time=2.0, noise=NoiseNone(),
        qe=SampledQuantumEfficiency([0.50e-6, 0.60e-6], [0.2, 0.8]),
        response_model=NullFrameResponse())
    spectral_plan = prepare_detector_acquisition(spectral_detector,
        spectral_rate)
    @test capture!(spectral_detector, spectral_rate, spectral_plan;
        rng=MersenneTwister(200)) ≈ fill(1.0, 2, 2)
    @test_throws InvalidConfiguration prepare_detector_acquisition(
        spectral_detector, detector_test_intensity_map(ones(2, 2);
            spectral=UnspecifiedSpectralCoordinate()))
    @test_throws InvalidConfiguration prepare_detector_acquisition(
        Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0,
            response_model=NullFrameResponse()),
        detector_test_intensity_map(ones(2, 2);
            spectral=UnspecifiedSpectralCoordinate()))
    @test_throws InvalidConfiguration prepare_detector_acquisition(
        spectral_detector, detector_test_intensity_map(ones(2, 2);
            spectral=IntegratedSpectralChannel(:science_passband)))

    for (qe, spectral) in (
        (1.0, UnspecifiedSpectralCoordinate()),
        (SampledQuantumEfficiency([0.50e-6, 0.60e-6], [0.2, 0.8]),
            IntegratedSpectralChannel(:science_passband)),
    )
        rejecting_detector = Detector(integration_time=1.0,
            noise=NoiseNone(), qe=qe, response_model=NullFrameResponse())
        state_before = detector_state_snapshot(rejecting_detector)
        @test_throws InvalidConfiguration prepare_detector_acquisition(
            rejecting_detector, detector_test_intensity_map(ones(2, 2);
                spectral=spectral))
        @test detector_state_matches_snapshot(rejecting_detector,
            state_before)
    end

    invalid_window_detector = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, readout_window=FrameWindow(1:3, 1:3),
        response_model=NullFrameResponse())
    invalid_window_state_before = detector_state_snapshot(
        invalid_window_detector)
    @test_throws DimensionMismatchError prepare_detector_acquisition(
        invalid_window_detector, detector_test_intensity_map(ones(2, 2)))
    @test detector_state_matches_snapshot(invalid_window_detector,
        invalid_window_state_before)

    response_kernel = [0.0 0.1 0.0; 0.1 0.6 0.1; 0.0 0.1 0.0]
    density_values = zeros(9, 9)
    density_values[3, 5] = 8.0
    density_rate = detector_test_intensity_map(density_values;
        sampling=(0.5, 0.25), spatial_measure=SpatialDensityMeasure())
    density_detector = Detector(integration_time=2.0, noise=NoiseNone(),
        qe=0.5, binning=3,
        response_model=SampledFrameResponse(response_kernel))
    density_plan = prepare_detector_acquisition(density_detector, density_rate)
    density_frame = copy(capture!(density_detector, density_rate,
        density_plan; rng=MersenneTwister(201)))
    manual_response = zeros(9, 9)
    manual_response[3, 5] = 4.8
    manual_response[2, 5] = 0.8
    manual_response[4, 5] = 0.8
    manual_response[3, 4] = 0.8
    manual_response[3, 6] = 0.8
    manual_binned = zeros(3, 3)
    bin2d!(manual_binned, manual_response, 3)
    @test density_frame ≈ manual_binned .* 0.125
    @test sum(density_frame) ≈ sum(density_values) * 0.125
    @test_detector_allocation prepared_detector_capture_allocations(
        density_detector,
        density_rate, density_plan, Xoshiro(202)) == 0

    cell_rate = detector_test_intensity_map(copy(density_values);
        sampling=(0.5, 0.25), spatial_measure=CellIntegratedMeasure())
    cell_detector = Detector(integration_time=2.0, noise=NoiseNone(),
        qe=0.5, binning=3,
        response_model=SampledFrameResponse(response_kernel))
    cell_plan = prepare_detector_acquisition(cell_detector, cell_rate)
    cell_frame = capture!(cell_detector, cell_rate, cell_plan;
        rng=MersenneTwister(202))
    @test cell_frame ≈ manual_binned
    @test sum(cell_frame) ≈ sum(density_values)

    normalized_map = detector_test_intensity_map(ones(2, 2);
        normalization=DimensionlessNormalization())
    normalized_detector = Detector(integration_time=0.5,
        noise=NoiseNone(), qe=0.25, response_model=NullFrameResponse())
    @test_throws InvalidConfiguration prepare_detector_acquisition(
        normalized_detector, normalized_map)
    normalized_plan = prepare_detector_acquisition(normalized_detector,
        normalized_map; normalized_to_photon_rate=40.0)
    @test capture!(normalized_detector, normalized_map, normalized_plan;
        rng=MersenneTwister(203)) == fill(5.0, 2, 2)
    @test_throws InvalidConfiguration prepare_detector_acquisition(
        normalized_detector, normalized_map; normalized_to_photon_rate=0.0)
    @test_throws InvalidConfiguration prepare_detector_acquisition(
        normalized_detector, shared_rate; normalized_to_photon_rate=40.0)

    # An external optical executor, including Proper.jl, enters through the
    # same declared product boundary. The package does not infer semantics from
    # an otherwise bare array.
    external_values = fill(3.0, 2, 2)
    external_metadata = OpticalPlaneMetadata(DetectorPlane(), external_values;
        coordinate_domain=AngularCoordinates(), sampling=(0.2, 0.2),
        spectral=IntegratedSpectralChannel(:proper_science_passband),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    external_product = IntensityMap(external_metadata, external_values)
    external_detector = Detector(integration_time=0.4, noise=NoiseNone(),
        qe=0.5, response_model=NullFrameResponse())
    external_plan = prepare_detector_acquisition(external_detector,
        external_product)
    @test capture!(external_detector, external_product, external_plan;
        rng=MersenneTwister(2031)) ≈ fill(0.6, 2, 2)
    undeclared_external_metadata = OpticalPlaneMetadata(DetectorPlane(),
        external_values; coordinate_domain=AngularCoordinates(),
        sampling=(0.2, 0.2))
    undeclared_external_product = IntensityMap(
        undeclared_external_metadata, external_values)
    @test_throws InvalidConfiguration prepare_detector_acquisition(
        external_detector, undeclared_external_product)

    pupil_rate = detector_test_intensity_map(ones(2, 2); kind=PupilPlane())
    point_rate = detector_test_intensity_map(ones(2, 2);
        spatial_measure=PointSampledMeasure())
    coherent_rate = detector_test_intensity_map(ones(2, 2);
        coherence=CoherentFieldCombination())
    @test_throws InvalidConfiguration prepare_detector_acquisition(
        normalized_detector, pupil_rate)
    @test_throws InvalidConfiguration prepare_detector_acquisition(
        normalized_detector, point_rate)
    @test_throws InvalidConfiguration prepare_detector_acquisition(
        normalized_detector, coherent_rate)
    float32_detector = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, response_model=NullFrameResponse(), T=Float32)
    @test_throws InvalidConfiguration prepare_detector_acquisition(
        float32_detector, shared_rate)

    oversampled_response_detector = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, psf_sampling=2,
        response_model=SampledFrameResponse(response_kernel))
    @test_throws InvalidConfiguration prepare_detector_acquisition(
        oversampled_response_detector, shared_rate)
    @test_throws InvalidConfiguration capture!(oversampled_response_detector,
        rate_values; rng=MersenneTwister(204))

    allocation_detector = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse())
    allocation_plan = prepare_detector_acquisition(allocation_detector,
        shared_rate)
    @test_detector_allocation prepared_detector_capture_allocations(
        allocation_detector,
        shared_rate, allocation_plan, Xoshiro(205)) == 0
    @test_detector_allocation prepared_detector_readiness_allocations(
        allocation_detector,
        shared_rate, allocation_plan) == 0

    busy_prepared_detector = Detector(integration_time=2.0,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse())
    busy_prepared_plan = prepare_detector_acquisition(
        busy_prepared_detector, shared_rate)
    capture!(busy_prepared_detector, shared_rate, busy_prepared_plan;
        rng=MersenneTwister(2399), integration_duration=0.25)
    busy_integrated_time = busy_prepared_detector.state.integrated_time
    busy_accumulation = copy(busy_prepared_detector.state.accum_buffer)
    @test_throws InvalidConfiguration begin
        AdaptiveOpticsSim.Detectors._require_prepared_whole_acquisition(
            busy_prepared_detector, shared_rate, busy_prepared_plan)
    end
    @test busy_prepared_detector.state.integrated_time == busy_integrated_time
    @test busy_prepared_detector.state.accum_buffer == busy_accumulation

    prepared_readout_rate = detector_test_intensity_map(ones(4, 4))
    prepared_readout_builders = (
        :skipper => (() -> Detector(integration_time=2.0,
            noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
            sensor=CCDSensor(sampling_mode=SkipperSampling(4)))),
        :hgcdte_single => (() -> Detector(integration_time=2.0,
            noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
            sensor=HgCdTeSensor())),
        :hgcdte_ndr => (() -> Detector(integration_time=2.0,
            noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
            sensor=HgCdTeSensor(
                sampling_mode=AveragedNonDestructiveReads(4)))),
        :hgcdte_cds => (() -> Detector(integration_time=2.0,
            noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
            sensor=HgCdTeSensor(
                sampling_mode=CorrelatedDoubleSampling()))),
        :hgcdte_fowler => (() -> Detector(integration_time=2.0,
            noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
            sensor=HgCdTeSensor(
                sampling_mode=FowlerSampling(2)))),
        :hgcdte_ramp => (() -> Detector(integration_time=2.0,
            noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
            sensor=HgCdTeSensor(read_time=0.1,
                sampling_mode=UpTheRampSampling(5)))),
        :hgcdte_windowed_cds => (() -> Detector(integration_time=2.0,
            noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
            sensor=HgCdTeSensor(
                sampling_mode=CorrelatedDoubleSampling()),
            readout_window=FrameWindow(2:3, 2:3))),
    )
    for (label, builder) in prepared_readout_builders
        prepared_detector = builder()
        prepared_plan = prepare_detector_acquisition(prepared_detector,
            prepared_readout_rate)
        @test AdaptiveOpticsSim.Detectors._require_exact_detector_acquisition_target(
            prepared_detector, prepared_plan, host_target) === prepared_plan
        @test readout_products(prepared_detector) isa FrameReadoutProducts
        @test !(readout_products(prepared_detector) isa
            NoFrameReadoutProducts)
        @test prepared_detector_exposed_storage_is_zero(prepared_detector)
        prepared_metadata = detector_export_metadata(prepared_detector)
        @test prepared_metadata.provides_signal_frame
        @test prepared_metadata.provides_combined_frame
        products = readout_products(prepared_detector)
        if products isa MultiReadFrameReadoutProducts
            @test isnothing(products.workspace_reference_cube) ||
                all(iszero, products.workspace_reference_cube)
            @test all(iszero, products.workspace_signal_cube)
        elseif products isa UpTheRampReadoutProducts
            @test all(iszero, products.workspace_slope)
            @test all(iszero, products.workspace_intercept)
            @test all(iszero, products.workspace_integrated)
            @test all(iszero, products.workspace_cube)
        end
        if coverage_instrumented()
            @test_skip "first prepared detector capture allocation assertion is disabled under coverage instrumentation: $label"
        else
            @test prepared_first_detector_capture_allocations(builder,
                prepared_readout_rate) == 0
        end
    end

    invalid_prepared_defect = Detector(noise=NoiseNone(),
        response_model=NullFrameResponse(), sensor=CMOSSensor(),
        defect_model=PixelResponseNonuniformity(ones(2, 2)))
    invalid_prepared_defect_state = detector_state_snapshot(
        invalid_prepared_defect)
    @test_throws DimensionMismatchError prepare_detector_acquisition(
        invalid_prepared_defect, prepared_readout_rate)
    @test detector_state_matches_snapshot(invalid_prepared_defect,
        invalid_prepared_defect_state)

    invalid_prepared_flux = Detector(noise=NoiseNone(),
        response_model=NullFrameResponse(), background_flux=ones(2, 2))
    invalid_prepared_flux_state = detector_state_snapshot(
        invalid_prepared_flux)
    @test_throws DimensionMismatchError prepare_detector_acquisition(
        invalid_prepared_flux, prepared_readout_rate)
    @test detector_state_matches_snapshot(invalid_prepared_flux,
        invalid_prepared_flux_state)

    invalid_prepared_background = Detector(noise=NoiseNone(),
        response_model=NullFrameResponse(), background_map=ones(2, 2))
    invalid_prepared_background_state = detector_state_snapshot(
        invalid_prepared_background)
    @test_throws DimensionMismatchError prepare_detector_acquisition(
        invalid_prepared_background, prepared_readout_rate)
    @test detector_state_matches_snapshot(invalid_prepared_background,
        invalid_prepared_background_state)

    invalid_prepared_cmos_noise = Detector(noise=NoiseNone(),
        response_model=NullFrameResponse(), sensor=CMOSSensor(
            readout_noise_model=CMOSReadNoiseMap(ones(2, 2))))
    invalid_prepared_cmos_noise_state = detector_state_snapshot(
        invalid_prepared_cmos_noise)
    @test_throws DimensionMismatchError prepare_detector_acquisition(
        invalid_prepared_cmos_noise, prepared_readout_rate)
    @test detector_state_matches_snapshot(invalid_prepared_cmos_noise,
        invalid_prepared_cmos_noise_state)

    invalid_prepared_cmos_output = Detector(noise=NoiseNone(),
        response_model=NullFrameResponse(), sensor=CMOSSensor(
            output_model=StaticCMOSOutputPattern(1, [1.0], [0.0])))
    invalid_prepared_cmos_output_state = detector_state_snapshot(
        invalid_prepared_cmos_output)
    @test_throws DimensionMismatchError prepare_detector_acquisition(
        invalid_prepared_cmos_output, prepared_readout_rate)
    @test detector_state_matches_snapshot(invalid_prepared_cmos_output,
        invalid_prepared_cmos_output_state)

    invalid_prepared_ramp = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
        sensor=HgCdTeSensor(read_time=0.3,
            sampling_mode=UpTheRampSampling(5)))
    invalid_prepared_ramp_state = detector_state_snapshot(
        invalid_prepared_ramp)
    @test_throws InvalidConfiguration prepare_detector_acquisition(
        invalid_prepared_ramp, prepared_readout_rate)
    @test detector_state_matches_snapshot(invalid_prepared_ramp,
        invalid_prepared_ramp_state)

    full_cds = Detector(integration_time=2.0, noise=NoiseReadout(2.0),
        qe=1.0, response_model=NullFrameResponse(),
        sensor=HgCdTeSensor(read_time=0.1,
            sampling_mode=CorrelatedDoubleSampling()))
    windowed_cds = Detector(integration_time=2.0,
        noise=NoiseReadout(2.0), qe=1.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeSensor(read_time=0.1,
            sampling_mode=CorrelatedDoubleSampling()),
        readout_window=FrameWindow(2:3, 2:3))
    full_cds_plan = prepare_detector_acquisition(full_cds,
        prepared_readout_rate)
    windowed_cds_plan = prepare_detector_acquisition(windowed_cds,
        prepared_readout_rate)
    full_cds_frame = copy(capture!(full_cds, prepared_readout_rate,
        full_cds_plan, MersenneTwister(2402)))
    windowed_cds_frame = copy(capture!(windowed_cds,
        prepared_readout_rate, windowed_cds_plan, MersenneTwister(2402)))
    @test windowed_cds_frame == full_cds_frame[2:3, 2:3]
    @test detector_reference_frame(windowed_cds) ==
        detector_reference_frame(full_cds)[2:3, 2:3]
    @test detector_signal_frame(windowed_cds) ==
        detector_signal_frame(full_cds)[2:3, 2:3]
    @test detector_combined_frame(windowed_cds) ==
        detector_combined_frame(full_cds)[2:3, 2:3]
    @test detector_reference_cube(windowed_cds) ==
        detector_reference_cube(full_cds)[2:3, 2:3, :]
    @test detector_signal_cube(windowed_cds) ==
        detector_signal_cube(full_cds)[2:3, 2:3, :]
    @test detector_read_cube(windowed_cds) ==
        detector_read_cube(full_cds)[2:3, 2:3, :]

    binned_ndr = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, binning=2, response_model=NullFrameResponse(),
        sensor=HgCdTeSensor(
            sampling_mode=AveragedNonDestructiveReads(4)))
    binned_ndr_plan = prepare_detector_acquisition(binned_ndr,
        prepared_readout_rate)
    @test capture!(binned_ndr, prepared_readout_rate, binned_ndr_plan,
        MersenneTwister(2403)) == fill(4.0, 2, 2)

    replacement_storage = copy(shared_rate.values)
    replacement_storage_map = IntensityMap(shared_rate.metadata,
        replacement_storage)
    @test replacement_storage_map.metadata === shared_rate.metadata
    @test replacement_storage_map.values !== shared_rate.values
    @test_throws InvalidConfiguration capture!(allocation_detector,
        replacement_storage_map, allocation_plan; rng=MersenneTwister(205))
    @test_detector_allocation prepared_detector_capture_allocations(
        allocation_detector,
        shared_rate, allocation_plan, Xoshiro(205)) == 0
    identical_allocation_detector = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse())
    @test identical_allocation_detector.params === allocation_detector.params
    @test identical_allocation_detector.state !== allocation_detector.state
    @test_throws InvalidConfiguration capture!(identical_allocation_detector,
        shared_rate, allocation_plan; rng=MersenneTwister(205))
    incremental_allocation_detector = Detector(integration_time=4.0,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse())
    incremental_allocation_plan = prepare_detector_acquisition(
        incremental_allocation_detector, shared_rate)
    @test_detector_allocation prepared_incremental_capture_allocations(
        incremental_allocation_detector, shared_rate,
        incremental_allocation_plan, Xoshiro(206), 0.5) == 0

    mismatched_rate = detector_test_intensity_map(copy(rate_values);
        sampling=(2.0, 1.0))
    @test_throws InvalidConfiguration capture!(allocation_detector,
        mismatched_rate, allocation_plan; rng=MersenneTwister(207))
    @test_throws InvalidConfiguration capture!(long_exposure, shared_rate,
        allocation_plan; rng=MersenneTwister(207))

end

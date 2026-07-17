function detector_test_intensity_map(values::AbstractMatrix{T};
    kind::AbstractOpticalPlaneKind=FocalPlane(),
    sampling::NTuple{2,T}=(one(T), one(T)),
    normalization::AbstractOpticalNormalization=PhotonRateNormalization(),
    spatial_measure::AbstractSpatialMeasure=CellIntegratedMeasure(),
    coherence::AbstractCombinationPolicy=IncoherentIntensityAddition(),
    spectral::AbstractSpectralCoordinate=MonochromaticChannel(0.55e-6)) where {T<:AbstractFloat}
    metadata = OpticalPlaneMetadata(kind, values;
        coordinate_domain=AngularCoordinates(), sampling=sampling,
        normalization=normalization, spatial_measure=spatial_measure,
        coherence=coherence, spectral=spectral)
    return IntensityMap(metadata, values)
end

function detector_state_snapshot(det::Detector)
    names = fieldnames(typeof(det.state))
    values = map(name -> getfield(det.state, name), names)
    return NamedTuple{names}(values)
end

function detector_state_matches_snapshot(det::Detector, snapshot::NamedTuple)
    return all(name -> getfield(det.state, name) === getfield(snapshot, name),
        keys(snapshot))
end

struct UnkeyedCalibrationFrameResponse{T<:AbstractFloat} <:
    AdaptiveOpticsSim.AbstractFrameResponse
    alpha::T
end

function AdaptiveOpticsSim.convert_frame_response_model(
    model::UnkeyedCalibrationFrameResponse, ::Type{T}, backend) where
    {T<:AbstractFloat}
    return UnkeyedCalibrationFrameResponse{T}(T(model.alpha))
end

AdaptiveOpticsSim.validate_frame_response_model(
    model::UnkeyedCalibrationFrameResponse) = model

@testset "Sampled detector parameter ownership" begin
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
    typed_qe = AdaptiveOpticsSim.SampledQuantumEfficiency{
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
    @test owned_cmos_sensor.timing_model.line_time == 0.01f0
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
        AdaptiveOpticsSim.detector_calibration_signature(
            unkeyed_response, UInt(0))
    end

    detector = Detector(integration_time=1.0, noise=NoiseNone(),
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
    signature_before = AdaptiveOpticsSim.detector_calibration_signature(
        detector, UInt(7))

    sampled_response.kernel[2, 2] = 0.2
    sampled_qe.values[2] = 0.1
    after = copy(capture!(detector, impulse, src;
        rng=MersenneTwister(2301)))
    @test after == before
    @test AdaptiveOpticsSim.detector_calibration_signature(
        detector, UInt(7)) == signature_before

    raw_gain_map = [1.0 0.8; 0.9 1.1]
    raw_dark_map = [0.1 0.2; 0.3 0.4]
    raw_bad_mask = Bool[false true; false false]
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
    defect_signature = AdaptiveOpticsSim.detector_calibration_signature(
        defect_detector, UInt(8))
    @test defect_detector.params.defect_model.stages[1].gain_map !==
        prnu.gain_map
    @test defect_detector.params.defect_model.stages[2].mask !==
        bad_pixels.mask
    prnu.gain_map[1, 1] = 0.2
    bad_pixels.mask[2, 2] = true
    @test composite_defects.stages[1].gain_map[1, 1] == 1.0
    @test !composite_defects.stages[2].mask[2, 2]
    @test AdaptiveOpticsSim.detector_calibration_signature(
        defect_detector, UInt(8)) == defect_signature

    replacement_detector = Detector(noise=NoiseNone(), sensor=CMOSSensor(),
        defect_model=CompositeDetectorDefectModel(
            PixelResponseNonuniformity([0.7 0.8; 0.9 1.1]),
            BadPixelMask(Bool[false true; true false]; throughput=0.0)))
    @test AdaptiveOpticsSim.detector_calibration_signature(
        replacement_detector, UInt(8)) != defect_signature
end

function prepared_detector_capture_allocations(det, map, plan, rng)
    capture!(det, map, plan, rng)
    return @allocated capture!(det, map, plan, rng)
end

function prepared_detector_readiness_allocations(det, map, plan)
    AdaptiveOpticsSim._require_prepared_whole_acquisition(det, map, plan)
    return @allocated AdaptiveOpticsSim._require_prepared_whole_acquisition(
        det, map, plan)
end

function prepared_first_detector_capture_allocations(builder, map)
    warm_detector = builder()
    warm_plan = prepare_detector_acquisition(warm_detector, map)
    capture!(warm_detector, map, warm_plan, Xoshiro(2400))

    detector = builder()
    plan = prepare_detector_acquisition(detector, map)
    rng = Xoshiro(2401)
    return @allocated capture!(detector, map, plan, rng)
end

function prepared_detector_exposed_storage_is_zero(det::Detector)
    products = (
        output_frame(det),
        detector_reference_frame(det),
        detector_signal_frame(det),
        detector_combined_frame(det),
        detector_reference_cube(det),
        detector_signal_cube(det),
        detector_read_cube(det),
        detector_read_times(det),
        detector_ramp_slope(det),
        detector_ramp_intercept(det),
    )
    return all(product -> isnothing(product) || all(iszero, product),
        products)
end

function prepared_incremental_capture_allocations(det, map, plan, rng,
    sample_time)
    capture!(det, map, plan; rng=rng, sample_time=sample_time)
    return @allocated capture!(det, map, plan; rng=rng,
        sample_time=sample_time)
end

function prepared_first_incremental_capture_allocations(det, map, plan, rng,
    sample_time)
    capture!(det, map, plan; rng=rng, sample_time=sample_time)
    reset_integration!(det)
    return @allocated capture!(det, map, plan; rng=rng,
        sample_time=sample_time)
end

function fixed_stack_capture_allocations(det, cube, scratch, rng)
    AdaptiveOpticsSim.capture_stack!(det, cube, scratch, rng)
    fill!(cube, one(eltype(cube)))
    return @allocated AdaptiveOpticsSim.capture_stack!(det, cube, scratch,
        rng)
end

@testset "Detector" begin
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
    @test prepared_detector_capture_allocations(density_detector,
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
    @test prepared_detector_capture_allocations(allocation_detector,
        shared_rate, allocation_plan, Xoshiro(205)) == 0
    @test prepared_detector_readiness_allocations(allocation_detector,
        shared_rate, allocation_plan) == 0

    busy_prepared_detector = Detector(integration_time=2.0,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse())
    busy_prepared_plan = prepare_detector_acquisition(
        busy_prepared_detector, shared_rate)
    capture!(busy_prepared_detector, shared_rate, busy_prepared_plan;
        rng=MersenneTwister(2399), sample_time=0.25)
    busy_integrated_time = busy_prepared_detector.state.integrated_time
    busy_accumulation = copy(busy_prepared_detector.state.accum_buffer)
    @test_throws InvalidConfiguration begin
        AdaptiveOpticsSim._require_prepared_whole_acquisition(
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
            sensor=HgCdTeAvalancheArraySensor())),
        :hgcdte_ndr => (() -> Detector(integration_time=2.0,
            noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
            sensor=HgCdTeAvalancheArraySensor(
                sampling_mode=AveragedNonDestructiveReads(4)))),
        :hgcdte_cds => (() -> Detector(integration_time=2.0,
            noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
            sensor=HgCdTeAvalancheArraySensor(
                sampling_mode=CorrelatedDoubleSampling()))),
        :hgcdte_fowler => (() -> Detector(integration_time=2.0,
            noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
            sensor=HgCdTeAvalancheArraySensor(
                sampling_mode=FowlerSampling(2)))),
        :hgcdte_ramp => (() -> Detector(integration_time=2.0,
            noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
            sensor=HgCdTeAvalancheArraySensor(read_time=0.1,
                sampling_mode=UpTheRampSampling(5)))),
        :hgcdte_windowed_cds => (() -> Detector(integration_time=2.0,
            noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
            sensor=HgCdTeAvalancheArraySensor(
                sampling_mode=CorrelatedDoubleSampling()),
            readout_window=FrameWindow(2:3, 2:3))),
    )
    for (label, builder) in prepared_readout_builders
        prepared_detector = builder()
        prepare_detector_acquisition(prepared_detector,
            prepared_readout_rate)
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
        sensor=HgCdTeAvalancheArraySensor(read_time=0.3,
            sampling_mode=UpTheRampSampling(5)))
    invalid_prepared_ramp_state = detector_state_snapshot(
        invalid_prepared_ramp)
    @test_throws InvalidConfiguration prepare_detector_acquisition(
        invalid_prepared_ramp, prepared_readout_rate)
    @test detector_state_matches_snapshot(invalid_prepared_ramp,
        invalid_prepared_ramp_state)

    full_cds = Detector(integration_time=2.0, noise=NoiseReadout(2.0),
        qe=1.0, response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(read_time=0.1,
            sampling_mode=CorrelatedDoubleSampling()))
    windowed_cds = Detector(integration_time=2.0,
        noise=NoiseReadout(2.0), qe=1.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(read_time=0.1,
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
        sensor=HgCdTeAvalancheArraySensor(
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
    @test prepared_detector_capture_allocations(allocation_detector,
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
    @test prepared_incremental_capture_allocations(
        incremental_allocation_detector, shared_rate,
        incremental_allocation_plan, Xoshiro(206), 0.5) == 0

    mismatched_rate = detector_test_intensity_map(copy(rate_values);
        sampling=(2.0, 1.0))
    @test_throws InvalidConfiguration capture!(allocation_detector,
        mismatched_rate, allocation_plan; rng=MersenneTwister(207))
    @test_throws InvalidConfiguration capture!(long_exposure, shared_rate,
        allocation_plan; rng=MersenneTwister(207))

    incremental_detector = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse())
    capture!(incremental_detector, ones(2, 2); rng=MersenneTwister(206),
        sample_time=0.6)
    @test_throws InvalidConfiguration capture!(incremental_detector,
        ones(2, 2); rng=MersenneTwister(206), sample_time=0.5)
    @test incremental_detector.state.integrated_time == 0.6
    @test !readout_ready(incremental_detector)
    incremental_pending_accum = copy(incremental_detector.state.accum_buffer)
    @test_throws DimensionMismatchError capture!(incremental_detector,
        ones(4, 4); rng=MersenneTwister(206), sample_time=0.1)
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
        rng=MersenneTwister(208), sample_time=0.25)
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
    @test_throws InvalidConfiguration AdaptiveOpticsSim.capture_stack!(
        transition_detector, pending_stack, similar(pending_stack);
        rng=MersenneTwister(213))
    pending_generalized_output = fill(UInt8(7), 1, 2, 2)
    pending_generalized_input = fill(3.0, 1, 2, 2)
    @test_throws InvalidConfiguration AdaptiveOpticsSim.capture_stack!(
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
        rng=MersenneTwister(214), sample_time=0.75)
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
    AdaptiveOpticsSim.capture_stack!(generalized_qe_detector,
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
    @test AdaptiveOpticsSim.detector_execution_plan(AdaptiveOpticsSim.execution_style(det_tuple.state.frame), det_tuple) isa AdaptiveOpticsSim.DetectorDirectPlan

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
        sensor=HgCdTeAvalancheArraySensor(sampling_mode=SingleRead()),
        response_model=NullFrameResponse(),
        correction_model=CompositeFrameReadoutCorrection((
            ReferenceRowCommonModeCorrection(1),
            ReferenceColumnCommonModeCorrection(1),
        )))
    adc_window_in = reshape(collect(1.0:96.0), 2, 6, 8)
    adc_window_out = Array{UInt16}(undef, 2, 4, 5)
    generalized_adc_window = AdaptiveOpticsSim.capture_stack!(det_adc_window_corr, adc_window_out, copy(adc_window_in);
        rng=MersenneTwister(10))
    @test size(generalized_adc_window) == (2, 4, 5)
    @test generalized_adc_window[1, :, :] == capture!(Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
            bits=8, full_well=100.0, readout_window=FrameWindow(2:5, 3:7), output_type=UInt16,
            sensor=HgCdTeAvalancheArraySensor(sampling_mode=SingleRead()),
            response_model=NullFrameResponse(),
            correction_model=CompositeFrameReadoutCorrection((
                ReferenceRowCommonModeCorrection(1),
                ReferenceColumnCommonModeCorrection(1),
            ))),
        @view(adc_window_in[1, :, :]); rng=MersenneTwister(10))
    @test generalized_adc_window[2, :, :] == capture!(Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
            bits=8, full_well=100.0, readout_window=FrameWindow(2:5, 3:7), output_type=UInt16,
            sensor=HgCdTeAvalancheArraySensor(sampling_mode=SingleRead()),
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

    zero_psf = zeros(4, 4)
    rng_ccd = MersenneTwister(7)
    rng_emccd = MersenneTwister(7)
    det_ccd = Detector(integration_time=1.0, noise=NoiseReadout(1.0), qe=1.0, binning=1,
        gain=10.0, sensor=CCDSensor())
    det_emccd = Detector(integration_time=1.0, noise=NoiseReadout(1.0), qe=1.0, binning=1,
        gain=10.0, sensor=EMCCDSensor())
    frame_ccd = copy(capture!(det_ccd, zero_psf; rng=rng_ccd))
    frame_emccd = copy(capture!(det_emccd, zero_psf; rng=rng_emccd))
    @test frame_ccd ≈ 10 .* frame_emccd
    @test_throws InvalidConfiguration EMCCDSensor(excess_noise_factor=0.5)

    uniform_signal = fill(50.0, 8, 8)
    det_emccd_base = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=EMCCDSensor())
    det_emccd_excess = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=EMCCDSensor(excess_noise_factor=sqrt(2.0)))
    frame_emccd_base = copy(capture!(det_emccd_base, uniform_signal; rng=MersenneTwister(8)))
    frame_emccd_excess = copy(capture!(det_emccd_excess, uniform_signal; rng=MersenneTwister(8)))
    @test frame_emccd_base == uniform_signal
    @test std(vec(frame_emccd_excess)) > 0
    det_emccd_stochastic = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=5.0, sensor=EMCCDSensor(excess_noise_factor=1.4,
            multiplication_model=StochasticMultiplicationRegister(0.6)))
    frame_emccd_stochastic = copy(capture!(det_emccd_stochastic, uniform_signal; rng=MersenneTwister(124)))
    @test std(vec(frame_emccd_stochastic)) > 0
    @test isapprox(mean(frame_emccd_stochastic), 250.0; rtol=0.1)
    @test isapprox(var(frame_emccd_stochastic), 1200.0; rtol=0.35)

    det_ccd_cic = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=CCDSensor(clock_induced_charge_per_frame=5.0))
    frame_ccd_cic = capture!(det_ccd_cic, zero_psf; rng=MersenneTwister(11))
    @test sum(frame_ccd_cic) > 0
    @test supports_clock_induced_charge(det_ccd_cic.params.sensor)
    det_ccd_cic_long = Detector(integration_time=10.0, noise=NoiseNone(),
        qe=1.0, sensor=CCDSensor(clock_induced_charge_per_frame=5.0))
    @test capture!(det_ccd_cic_long, zero_psf; rng=MersenneTwister(11)) ==
        frame_ccd_cic
    @test_throws InvalidConfiguration CCDSensor(
        clock_induced_charge_per_frame=-1.0)

    skipper_input = zeros(64, 64)
    skipper_single = Detector(noise=NoiseReadout(4.0), qe=1.0,
        sensor=CCDSensor(sampling_mode=SkipperSampling(1), read_time=2e-6),
        response_model=NullFrameResponse())
    skipper_many = Detector(noise=NoiseReadout(4.0), qe=1.0,
        sensor=CCDSensor(sampling_mode=SkipperSampling(16), read_time=2e-6),
        response_model=NullFrameResponse())
    skipper_single_frame = copy(capture!(skipper_single, skipper_input;
        rng=MersenneTwister(91)))
    skipper_many_rng = MersenneTwister(91)
    skipper_many_frame = copy(capture!(skipper_many, skipper_input,
        skipper_many_rng))
    @test std(skipper_many_frame) < 0.35 * std(skipper_single_frame)
    @test detector_signal_frame(skipper_many) == skipper_many_frame
    @test detector_combined_frame(skipper_many) == skipper_many_frame
    @test detector_read_cube(skipper_many) === nothing
    skipper_meta = detector_export_metadata(skipper_many)
    @test skipper_meta.sampling_mode == :skipper
    @test skipper_meta.sampling_reads == 16
    @test skipper_meta.sampling_signal_reads == 16
    @test skipper_meta.sampling_read_time == 2e-6
    @test skipper_meta.sampling_wallclock_time == 1.0 + 32e-6
    @test skipper_meta.provides_signal_frame
    @test !skipper_meta.provides_read_cube
    @test supports_nondestructive_reads(skipper_many.params.sensor)
    @test_throws InvalidConfiguration CCDSensor(
        sampling_mode=SkipperSampling(0))
    capture!(skipper_many, skipper_input, skipper_many_rng)
    @test @allocated(capture!(skipper_many, skipper_input,
        skipper_many_rng)) == 0
    skipper_transition = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
        sensor=CCDSensor(sampling_mode=SkipperSampling(2)))
    capture!(skipper_transition, ones(2, 2); rng=MersenneTwister(92))
    @test all(==(2.0), skipper_transition.state.accum_buffer)
    capture!(skipper_transition, zeros(2, 2); rng=MersenneTwister(93),
        sample_time=0.5)
    skipper_incremental_frame = copy(capture!(skipper_transition,
        zeros(2, 2); rng=MersenneTwister(94), sample_time=0.5))
    @test all(iszero, skipper_incremental_frame)
    @test readout_ready(skipper_transition)
    det_emccd_cic = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=EMCCDSensor(clock_induced_charge_per_frame=3.0))
    frame_emccd_cic = copy(capture!(det_emccd_cic, zero_psf;
        rng=MersenneTwister(125)))
    @test sum(frame_emccd_cic) > 0
    det_emccd_cic_long = Detector(integration_time=10.0, noise=NoiseNone(),
        qe=1.0, sensor=EMCCDSensor(clock_induced_charge_per_frame=3.0))
    @test capture!(det_emccd_cic_long, zero_psf;
        rng=MersenneTwister(125)) == frame_emccd_cic
    emccd_cic_whole = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, gain=5.0, response_model=NullFrameResponse(),
        sensor=EMCCDSensor(clock_induced_charge_per_frame=3.0))
    emccd_cic_split = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, gain=5.0, response_model=NullFrameResponse(),
        sensor=EMCCDSensor(clock_induced_charge_per_frame=3.0))
    emccd_cic_whole_frame = copy(capture!(emccd_cic_whole, zero_psf;
        rng=MersenneTwister(130)))
    capture!(emccd_cic_split, zero_psf; rng=MersenneTwister(130),
        sample_time=0.5)
    emccd_cic_split_frame = copy(capture!(emccd_cic_split, zero_psf;
        rng=MersenneTwister(130), sample_time=0.5))
    @test emccd_cic_split_frame == emccd_cic_whole_frame
    det_emccd_sat = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=5.0, sensor=EMCCDSensor(register_full_well=100.0))
    @test maximum(capture!(det_emccd_sat, fill(50.0, 4, 4); rng=MersenneTwister(126))) == 100.0
    emccd_saturated_stack = fill(50.0, 2, 4, 4)
    AdaptiveOpticsSim.capture_stack!(det_emccd_sat, emccd_saturated_stack,
        similar(emccd_saturated_stack); rng=MersenneTwister(126))
    @test all(==(100.0), emccd_saturated_stack)
    @test_throws InvalidConfiguration EMCCDSensor(
        clock_induced_charge_per_frame=-1.0)
    @test_throws InvalidConfiguration EMCCDSensor(register_full_well=0.0)
    @test_throws InvalidConfiguration EMCCDSensor(multiplication_model=StochasticMultiplicationRegister(-1.0))
    @test_throws InvalidConfiguration EMCCDSensor(em_gain_range=(10.0, 1.0))
    @test_throws InvalidConfiguration EMCCDSensor(readout_rate_hz=-1.0)
    @test_throws InvalidConfiguration FrameTransferAcquisition(
        transfer_time=-1.0)
    @test_throws InvalidConfiguration PhotonCountingEMMode(threshold=-1.0)
    @test_throws InvalidConfiguration PhotonCountingEMMode(threshold=1.0, detection_efficiency=1.5)

    emccd_timing_input = fill(2.0, 4, 4)
    emccd_sequential = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, gain=1.0, response_model=NullFrameResponse(),
        sensor=EMCCDSensor(readout_rate_hz=1000.0,
            acquisition_mode=SequentialAcquisition()))
    emccd_frame_transfer = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, gain=1.0,
        response_model=NullFrameResponse(),
        sensor=EMCCDSensor(readout_rate_hz=1000.0,
            acquisition_mode=FrameTransferAcquisition(
                transfer_time=0.002)))
    sequential_frame = copy(capture!(emccd_sequential, emccd_timing_input;
        rng=MersenneTwister(126)))
    frame_transfer_frame = copy(capture!(emccd_frame_transfer,
        emccd_timing_input; rng=MersenneTwister(126)))
    @test frame_transfer_frame == sequential_frame
    sequential_meta = detector_export_metadata(emccd_sequential)
    frame_transfer_meta = detector_export_metadata(emccd_frame_transfer)
    @test sequential_meta.acquisition_mode == :sequential
    @test sequential_meta.frame_transfer_time === nothing
    @test sequential_meta.sampling_read_time == 0.016
    @test sequential_meta.sampling_wallclock_time == 1.016
    @test sequential_meta.steady_state_frame_period == 1.016
    @test frame_transfer_meta.acquisition_mode == :frame_transfer
    @test frame_transfer_meta.frame_transfer_time == 0.002
    @test frame_transfer_meta.sampling_read_time == 0.016
    @test frame_transfer_meta.sampling_wallclock_time == 1.018
    @test frame_transfer_meta.steady_state_frame_period == 1.002

    emccd_readout_limited = Detector(integration_time=0.001,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
        sensor=EMCCDSensor(readout_rate_hz=1000.0,
            acquisition_mode=FrameTransferAcquisition(
                transfer_time=0.002)))
    capture!(emccd_readout_limited, emccd_timing_input;
        rng=MersenneTwister(126))
    readout_limited_meta = detector_export_metadata(emccd_readout_limited)
    @test readout_limited_meta.steady_state_frame_period ≈ 0.018

    emccd_unknown_timing = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0,
        sensor=EMCCDSensor(acquisition_mode=FrameTransferAcquisition()))
    capture!(emccd_unknown_timing, emccd_timing_input;
        rng=MersenneTwister(126))
    unknown_timing_meta = detector_export_metadata(emccd_unknown_timing)
    @test unknown_timing_meta.sampling_read_time === nothing
    @test unknown_timing_meta.sampling_wallclock_time === nothing
    @test unknown_timing_meta.steady_state_frame_period === nothing

    det_emccd_conventional = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=10.0, sensor=EMCCDSensor(output_path=ConventionalOutput()))
    @test capture!(det_emccd_conventional, uniform_signal; rng=MersenneTwister(127)) == uniform_signal

    det_emccd_pc = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=10.0,
        sensor=EMCCDSensor(operating_mode=PhotonCountingEMMode(threshold=5.0)))
    pc_frame = capture!(det_emccd_pc, [0.0 0.4; 0.6 1.0]; rng=MersenneTwister(128))
    @test pc_frame == [0.0 0.0; 1.0 1.0]
    @test supports_photon_number_resolving(det_emccd_pc.params.sensor)

    det_emccd_efficiency = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, gain=10.0,
        sensor=EMCCDSensor(operating_mode=PhotonCountingEMMode(
            threshold=5.0, detection_efficiency=0.8)))
    efficiency_frame = capture!(det_emccd_efficiency, ones(100, 100);
        rng=MersenneTwister(128))
    @test all(x -> x == 0.0 || x == 1.0, efficiency_frame)
    @test isapprox(mean(efficiency_frame), 0.8; atol=0.025)

    det_emccd_pc_batched = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=10.0, response_model=NullFrameResponse(),
        sensor=EMCCDSensor(operating_mode=PhotonCountingEMMode(threshold=5.0)))
    pc_cube = reshape(Float64[0.0, 0.4, 0.6, 1.0], 1, 2, 2)
    pc_scratch = similar(pc_cube)
    capture_stack!(det_emccd_pc_batched, pc_cube, pc_scratch, MersenneTwister(129))
    @test pc_cube == reshape(Float64[0.0, 0.0, 1.0, 1.0], 1, 2, 2)

    @test emccd_snr(1.0; readout_noise=20.0, gain=100.0, excess_noise_factor=1.0) >
        emccd_snr(1.0; readout_noise=20.0, gain=100.0, excess_noise_factor=sqrt(2.0))
    @test emccd_snr(1.0; readout_noise=20.0, gain=100.0, operating_mode=PhotonCountingEMMode(threshold=0.5)) >
        emccd_snr(1.0; readout_noise=20.0, gain=100.0, excess_noise_factor=sqrt(2.0))
    @test emccd_snr(1.0; readout_noise=20.0, gain=100.0, output_path=EMOutput()) >
        emccd_snr(1.0; readout_noise=20.0, gain=100.0, output_path=ConventionalOutput())

    det_cmos = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=CMOSSensor(column_readout_sigma=1.0))
    frame_cmos = copy(capture!(det_cmos, zeros(8, 8); rng=MersenneTwister(12)))
    @test !all(iszero, frame_cmos)
    @test all(j -> isapprox(std(frame_cmos[:, j]), 0.0; atol=1e-8), axes(frame_cmos, 2))
    @test std(vec(frame_cmos[1, :])) > 0
    @test supports_column_readout_noise(det_cmos.params.sensor)
    @test detector_export_metadata(det_cmos).frame_response == :none
    @test_throws InvalidConfiguration CMOSSensor(column_readout_sigma=-1.0)

    det_cmos_rows = Detector(noise=NoiseNone(), qe=1.0,
        sensor=CMOSSensor(row_readout_sigma=1.0))
    frame_cmos_rows = copy(capture!(det_cmos_rows, zeros(8, 8);
        rng=MersenneTwister(1212)))
    @test all(i -> isapprox(std(frame_cmos_rows[i, :]), 0.0; atol=1e-8),
        axes(frame_cmos_rows, 1))
    @test std(frame_cmos_rows[:, 1]) > 0

    sigma_map = zeros(4, 4)
    sigma_map[2, 3] = 2.0
    det_cmos_map = Detector(noise=NoiseNone(), qe=1.0,
        sensor=CMOSSensor(readout_noise_model=CMOSReadNoiseMap(sigma_map)))
    frame_cmos_map = copy(capture!(det_cmos_map, zeros(4, 4);
        rng=MersenneTwister(1213)))
    @test count(x -> !iszero(x), frame_cmos_map) == 1
    @test frame_cmos_map[2, 3] != 0
    @test_throws InvalidConfiguration CMOSSensor(row_readout_sigma=-1.0)
    @test_throws InvalidConfiguration CMOSReadNoiseMap(fill(-1.0, 2, 2))
    prnu_map = [1.0 0.5 1.0 0.5; 1.0 0.5 1.0 0.5; 1.0 0.5 1.0 0.5; 1.0 0.5 1.0 0.5]
    dsnu_map = fill(0.25, 4, 4)
    bad_mask = falses(4, 4)
    bad_mask[2, 3] = true
    det_cmos_structured = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=CMOSSensor(output_model=StaticCMOSOutputPattern(2, [1.0, 2.0], [0.0, 10.0]),
            timing_model=RollingShutter(1e-3)),
        response_model=NullFrameResponse(),
        defect_model=CompositeDetectorDefectModel(
            PixelResponseNonuniformity(prnu_map),
            DarkSignalNonuniformity(dsnu_map),
            BadPixelMask(bad_mask; throughput=0.0)))
    structured_frame = capture!(det_cmos_structured, fill(2.0, 4, 4); rng=MersenneTwister(120))
    @test structured_frame[1, 1] ≈ 2.25
    @test structured_frame[1, 2] ≈ 1.25
    @test structured_frame[1, 3] ≈ 14.5
    @test structured_frame[2, 3] ≈ 10.5
    structured_meta = detector_export_metadata(det_cmos_structured)
    @test structured_meta.detector_defects == :composite
    @test structured_meta.has_prnu
    @test structured_meta.has_dsnu
    @test structured_meta.has_bad_pixels
    @test structured_meta.timing_model == :rolling_shutter
    @test structured_meta.timing_line_time == 1e-3
    @test structured_meta.sampling_wallclock_time == 1.004
    @test supports_detector_defect_maps(det_cmos_structured.params.sensor)
    @test supports_shutter_timing(det_cmos_structured.params.sensor)
    @test_throws InvalidConfiguration CMOSSensor(timing_model=RollingShutter(-1.0))
    @test_throws InvalidConfiguration CMOSSensor(timing_model=RollingShutter(1e-3; row_group_size=0))
    @test_throws InvalidConfiguration CMOSSensor(output_model=StaticCMOSOutputPattern(2, [1.0], [0.0, 1.0]))

    rolling_det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=CMOSSensor(timing_model=RollingShutter(0.25)),
        response_model=NullFrameResponse())
    rolling_source = InPlaceFrameSource((out, t) -> fill!(out, t), (4, 4))
    rolling_rng = MersenneTwister(127)
    rolling_frame = capture!(rolling_det, rolling_source, rolling_rng)
    @test rolling_frame == repeat(reshape([0.0, 0.25, 0.5, 0.75], :, 1), 1, 4)
    @test detector_export_metadata(rolling_det).sampling_wallclock_time == 2.0
    @test @allocated(capture!(rolling_det, rolling_source, rolling_rng)) == 0

    global_exposure_calls = Tuple{Float64,Float64}[]
    global_exposure_source = FunctionExposureFrameSource(
        (start_time, exposure_time) -> begin
            push!(global_exposure_calls, (start_time, exposure_time))
            fill(exposure_time, 2, 2)
        end)
    global_exposure_det = Detector(integration_time=2.0,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse())
    global_exposure_frame = capture!(global_exposure_det,
        global_exposure_source; rng=MersenneTwister(1271))
    @test global_exposure_calls == [(0.0, 2.0)]
    @test global_exposure_frame == fill(4.0, 2, 2)

    rolling_exposure_calls = Tuple{Float64,Float64}[]
    rolling_exposure_source = InPlaceExposureFrameSource(
        (out, start_time, exposure_time) -> begin
            push!(rolling_exposure_calls, (start_time, exposure_time))
            fill!(out, start_time + exposure_time)
        end, (4, 4))
    rolling_exposure_det = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0,
        sensor=CMOSSensor(timing_model=RollingShutter(0.25;
            row_group_size=2)), response_model=NullFrameResponse())
    rolling_exposure_frame = capture!(rolling_exposure_det,
        rolling_exposure_source; rng=MersenneTwister(1272))
    @test rolling_exposure_calls == [(0.0, 1.0), (0.25, 1.0)]
    @test rolling_exposure_frame ==
        repeat(reshape([1.0, 1.0, 1.25, 1.25], :, 1), 1, 4)

    global_reset_det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=CMOSSensor(timing_model=RollingShutter(0.25; exposure_mode=GlobalResetExposure())),
        response_model=NullFrameResponse())
    constant_source = FunctionFrameSource(t -> ones(4, 4))
    global_reset_frame = capture!(global_reset_det, constant_source; rng=MersenneTwister(129))
    @test global_reset_frame == repeat(reshape([1.0, 1.25, 1.5, 1.75], :, 1), 1, 4)
    @test global_reset_det.params.timing_model.exposure_mode == GlobalResetExposure()

    interval_source = FunctionExposureFrameSource((start_time, exposure_time) ->
        fill(start_time <= 1.4 < start_time + exposure_time ? 10.0 : 0.0, 4, 4))
    rolling_interval_frame = capture!(rolling_det, interval_source; rng=MersenneTwister(132))
    @test rolling_interval_frame[1:2, :] == zeros(2, 4)
    @test rolling_interval_frame[3:4, :] == fill(10.0, 2, 4)
    global_reset_interval_frame = capture!(global_reset_det, interval_source; rng=MersenneTwister(133))
    @test global_reset_interval_frame[1:2, :] == zeros(2, 4)
    @test global_reset_interval_frame[3, :] == fill(15.0, 4)
    @test global_reset_interval_frame[4, :] == fill(17.5, 4)

    pulse_source = FunctionFrameSource(t -> fill(t >= 0.5 ? 10.0 : 0.0, 4, 4))
    pulse_frame = capture!(rolling_det, pulse_source; rng=MersenneTwister(128))
    @test pulse_frame[1:2, :] == zeros(2, 4)
    @test pulse_frame[3:4, :] == fill(10.0, 2, 4)


    det_ingaas = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=InGaAsSensor(glow_rate=3.0))
    frame_ingaas = capture!(det_ingaas, zero_psf; rng=MersenneTwister(13))
    @test sum(frame_ingaas) > 0
    @test supports_sensor_glow(det_ingaas.params.sensor)
    @test detector_export_metadata(det_ingaas).frame_response == :gaussian
    @test_throws InvalidConfiguration InGaAsSensor(glow_rate=-1.0)
    det_ingaas_persist = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=NullFrameResponse(),
        sensor=InGaAsSensor(persistence_model=ExponentialPersistence(0.5, 0.0)))
    capture!(det_ingaas_persist, fill(4.0, 4, 4); rng=MersenneTwister(121))
    persisted = capture!(det_ingaas_persist, zeros(4, 4); rng=MersenneTwister(122))
    @test sum(persisted) ≈ 32.0

    persistence_sensor = InGaAsSensor(
        persistence_model=ExponentialPersistence(0.5, 0.0))
    persistence_whole = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, binning=1, response_model=NullFrameResponse(),
        sensor=persistence_sensor)
    persistence_split = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, binning=1, response_model=NullFrameResponse(),
        sensor=persistence_sensor)
    persistence_seed = fill(4.0, 4, 4)
    persistence_dark = zeros(4, 4)
    capture!(persistence_whole, persistence_seed; rng=MersenneTwister(123))
    capture!(persistence_split, persistence_seed; rng=MersenneTwister(123))
    persistence_whole_frame = copy(capture!(persistence_whole,
        persistence_dark; rng=MersenneTwister(124)))
    persistence_partial_frame = copy(capture!(persistence_split,
        persistence_dark; rng=MersenneTwister(124), sample_time=0.5))
    persistence_split_frame = copy(capture!(persistence_split,
        persistence_dark; rng=MersenneTwister(124), sample_time=0.5))
    @test persistence_partial_frame == persistence_whole_frame
    @test persistence_split_frame == persistence_whole_frame
    @test persistence_split.state.latent_buffer ==
        persistence_whole.state.latent_buffer

    persistence_prepared_whole = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, binning=1,
        response_model=NullFrameResponse(), sensor=persistence_sensor)
    persistence_prepared_split = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, binning=1,
        response_model=NullFrameResponse(), sensor=persistence_sensor)
    persistence_rate_map = detector_test_intensity_map(persistence_dark)
    persistence_whole_plan = prepare_detector_acquisition(
        persistence_prepared_whole, persistence_rate_map)
    persistence_split_plan = prepare_detector_acquisition(
        persistence_prepared_split, persistence_rate_map)
    capture!(persistence_prepared_whole, persistence_seed;
        rng=MersenneTwister(125))
    capture!(persistence_prepared_split, persistence_seed;
        rng=MersenneTwister(125))
    persistence_prepared_whole_frame = copy(capture!(
        persistence_prepared_whole, persistence_rate_map,
        persistence_whole_plan; rng=MersenneTwister(126)))
    capture!(persistence_prepared_split, persistence_rate_map,
        persistence_split_plan; rng=MersenneTwister(126), sample_time=0.5)
    persistence_prepared_split_frame = copy(capture!(
        persistence_prepared_split, persistence_rate_map,
        persistence_split_plan; rng=MersenneTwister(126), sample_time=0.5))
    @test persistence_prepared_split_frame ==
        persistence_prepared_whole_frame
    @test persistence_prepared_split.state.latent_buffer ==
        persistence_prepared_whole.state.latent_buffer

    persistence_allocation_detector = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, binning=1,
        response_model=NullFrameResponse(), sensor=persistence_sensor)
    persistence_allocation_plan = prepare_detector_acquisition(
        persistence_allocation_detector, persistence_rate_map)
    capture!(persistence_allocation_detector, persistence_seed;
        rng=Xoshiro(127))
    @test prepared_first_incremental_capture_allocations(
        persistence_allocation_detector, persistence_rate_map,
        persistence_allocation_plan, Xoshiro(128), 0.5) == 0
    persist_meta = detector_export_metadata(det_ingaas_persist)
    @test persist_meta.persistence_model == :exponential
    det_ingaas_nonlinear = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=NullFrameResponse(),
        nonlinearity_model=SaturatingFrameNonlinearity(0.1),
        sensor=InGaAsSensor())
    nonlinear_frame = capture!(det_ingaas_nonlinear, fill(10.0, 2, 2); rng=MersenneTwister(123))
    @test nonlinear_frame == fill(5.0, 2, 2)
    nonlinear_meta = detector_export_metadata(det_ingaas_nonlinear)
    @test nonlinear_meta.nonlinearity_model == :saturating
    @test supports_detector_persistence(det_ingaas_persist.params.sensor)
    @test supports_detector_nonlinearity(det_ingaas_nonlinear.params.sensor)
    @test_throws InvalidConfiguration InGaAsSensor(persistence_model=ExponentialPersistence(1.1, 0.0))
    @test_throws InvalidConfiguration Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        nonlinearity_model=SaturatingFrameNonlinearity(-0.1), sensor=InGaAsSensor())

    arrhenius = ArrheniusRateLaw(300.0, 6000.0)
    linear = LinearTemperatureLaw(300.0, 0.01)
    exp_law = ExponentialTemperatureLaw(300.0, 0.01)
    @test evaluate_temperature_law(arrhenius, 10.0, 80.0) < 10.0
    @test evaluate_temperature_law(linear, 2.0, 250.0) ≈ 1.0
    @test evaluate_temperature_law(exp_law, 2.0, 250.0) < 2.0

    thermal_det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        dark_current=10.0,
        response_model=NullFrameResponse(),
        thermal_model=FixedTemperature(temperature_K=80.0, dark_current_law=arrhenius),
        sensor=CCDSensor())
    thermal_meta = detector_export_metadata(thermal_det)
    @test supports_detector_thermal_model(thermal_det)
    @test !supports_dynamic_thermal_state(thermal_det.params.thermal_model)
    @test supports_temperature_dependent_dark_current(thermal_det)
    @test detector_temperature(thermal_det) == 80.0
    @test thermal_meta.thermal_model == :fixed_temperature
    @test thermal_meta.detector_temperature_K == 80.0
    @test thermal_meta.cooling_setpoint_K == 80.0
    @test thermal_meta.dark_current_law == :arrhenius
    @test effective_dark_current(thermal_det) < thermal_det.params.dark_current
    @test thermal_state(thermal_det) isa NoThermalState
    @test advance_thermal!(thermal_det, 1.0) === thermal_det

    dynamic_model = FirstOrderThermalModel(
        ambient_temperature_K=295.0,
        setpoint_temperature_K=120.0,
        initial_temperature_K=300.0,
        time_constant_s=2.0,
        min_temperature_K=80.0,
        max_temperature_K=320.0,
        dark_current_law=arrhenius,
    )
    dynamic_det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        dark_current=10.0,
        response_model=NullFrameResponse(),
        thermal_model=dynamic_model,
        sensor=CCDSensor())
    dynamic_meta = detector_export_metadata(dynamic_det)
    @test supports_detector_thermal_model(dynamic_det)
    @test supports_dynamic_thermal_state(dynamic_det.params.thermal_model)
    @test thermal_state(dynamic_det) isa DetectorThermalState
    @test detector_temperature(dynamic_det) == 300.0
    @test dynamic_meta.thermal_model == :first_order
    @test dynamic_meta.detector_temperature_K == 300.0
    @test dynamic_meta.ambient_temperature_K == 295.0
    @test dynamic_meta.cooling_setpoint_K == 120.0
    @test dynamic_meta.thermal_time_constant_s == 2.0
    dark_current_initial = effective_dark_current(dynamic_det)
    @test advance_thermal!(dynamic_det, 2.0) === dynamic_det
    @test detector_temperature(dynamic_det) ≈ 120.0 + 180.0 * exp(-1.0)
    @test effective_dark_current(dynamic_det) < dark_current_initial
    reset_integration!(dynamic_det)
    capture!(dynamic_det, fill(1.0f0, 4, 4); rng=MersenneTwister(24))
    @test detector_temperature(dynamic_det) < 120.0 + 180.0 * exp(-1.0)
    @test_throws InvalidConfiguration advance_thermal!(dynamic_det, -1.0)
    @test_throws InvalidConfiguration FirstOrderThermalModel(
        ambient_temperature_K=295.0,
        setpoint_temperature_K=120.0,
        initial_temperature_K=60.0,
        time_constant_s=2.0,
        min_temperature_K=80.0,
        max_temperature_K=320.0,
    )

    dynamic_stack_det = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, response_model=NullFrameResponse(),
        thermal_model=dynamic_model, sensor=CCDSensor())
    dynamic_stack_cube = zeros(2, 2, 2)
    AdaptiveOpticsSim.capture_stack!(dynamic_stack_det, dynamic_stack_cube,
        similar(dynamic_stack_cube); rng=MersenneTwister(24))
    expected_stack_temperature = 120.0 + 180.0 * exp(-0.5)
    @test detector_temperature(dynamic_stack_det) ≈ expected_stack_temperature
    @test readout_ready(dynamic_stack_det)

    dynamic_generalized_det = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, psf_sampling=2,
        response_model=NullFrameResponse(), thermal_model=dynamic_model,
        sensor=CCDSensor())
    dynamic_generalized_input = zeros(2, 4, 4)
    dynamic_generalized_output = zeros(2, 2, 2)
    AdaptiveOpticsSim.capture_stack!(dynamic_generalized_det,
        dynamic_generalized_output, dynamic_generalized_input;
        rng=MersenneTwister(24))
    @test detector_temperature(dynamic_generalized_det) ≈
        expected_stack_temperature
    @test readout_ready(dynamic_generalized_det)

    incremental_rate_law = LinearTemperatureLaw(120.0, 0.005)
    incremental_rate_input = zeros(128, 128)
    static_rate_model = FixedTemperature(temperature_K=250.0,
        dark_current_law=incremental_rate_law)
    static_rate_whole = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, dark_current=100.0, response_model=NullFrameResponse(),
        thermal_model=static_rate_model, sensor=CCDSensor())
    static_rate_split = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, dark_current=100.0, response_model=NullFrameResponse(),
        thermal_model=static_rate_model, sensor=CCDSensor())
    static_rate_whole_frame = copy(capture!(static_rate_whole,
        incremental_rate_input; rng=MersenneTwister(25)))
    capture!(static_rate_split, incremental_rate_input;
        rng=MersenneTwister(26), sample_time=0.5)
    static_rate_split_frame = copy(capture!(static_rate_split,
        incremental_rate_input; rng=MersenneTwister(27), sample_time=0.5))
    static_rate_expected = evaluate_temperature_law(incremental_rate_law,
        100.0, 250.0)
    @test mean(static_rate_whole_frame) ≈ static_rate_expected rtol=0.01
    @test mean(static_rate_split_frame) ≈ static_rate_expected rtol=0.01
    @test mean(static_rate_split_frame) ≈
        mean(static_rate_whole_frame) rtol=0.01

    incremental_dynamic_model = FirstOrderThermalModel(
        ambient_temperature_K=300.0,
        setpoint_temperature_K=120.0,
        initial_temperature_K=300.0,
        time_constant_s=1.0,
        min_temperature_K=100.0,
        max_temperature_K=320.0,
        dark_current_law=incremental_rate_law,
    )
    dynamic_rate_whole = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, dark_current=100.0, response_model=NullFrameResponse(),
        thermal_model=incremental_dynamic_model, sensor=CCDSensor())
    dynamic_rate_split = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, dark_current=100.0, response_model=NullFrameResponse(),
        thermal_model=incremental_dynamic_model, sensor=CCDSensor())
    dynamic_rate_whole_frame = copy(capture!(dynamic_rate_whole,
        incremental_rate_input; rng=MersenneTwister(28)))
    capture!(dynamic_rate_split, incremental_rate_input;
        rng=MersenneTwister(29), sample_time=0.5)
    dynamic_rate_split_frame = copy(capture!(dynamic_rate_split,
        incremental_rate_input; rng=MersenneTwister(30), sample_time=0.5))
    half_exposure_temperature = 120.0 + 180.0 * exp(-0.5)
    dynamic_whole_expected = evaluate_temperature_law(
        incremental_rate_law, 100.0, 300.0)
    dynamic_split_expected = 0.5 * dynamic_whole_expected +
        0.5 * evaluate_temperature_law(incremental_rate_law, 100.0,
            half_exposure_temperature)
    @test mean(dynamic_rate_whole_frame) ≈ dynamic_whole_expected rtol=0.01
    @test mean(dynamic_rate_split_frame) ≈ dynamic_split_expected rtol=0.01
    @test detector_temperature(dynamic_rate_split) ≈
        120.0 + 180.0 * exp(-1.0)

    incremental_glow_model = FirstOrderThermalModel(
        ambient_temperature_K=300.0,
        setpoint_temperature_K=120.0,
        initial_temperature_K=300.0,
        time_constant_s=1.0,
        min_temperature_K=100.0,
        max_temperature_K=320.0,
        glow_rate_law=incremental_rate_law,
    )
    dynamic_glow_split = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, response_model=NullFrameResponse(),
        thermal_model=incremental_glow_model,
        sensor=InGaAsSensor(glow_rate=100.0))
    capture!(dynamic_glow_split, incremental_rate_input;
        rng=MersenneTwister(31), sample_time=0.5)
    dynamic_glow_frame = copy(capture!(dynamic_glow_split,
        incremental_rate_input; rng=MersenneTwister(32), sample_time=0.5))
    @test mean(dynamic_glow_frame) ≈ dynamic_split_expected rtol=0.01

    hgcdte_glow_sensor = HgCdTeAvalancheArraySensor(glow_rate=60.0,
        read_time=0.1, sampling_mode=CorrelatedDoubleSampling())
    hgcdte_glow_whole = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, dark_current=40.0, response_model=NullFrameResponse(),
        sensor=hgcdte_glow_sensor)
    hgcdte_glow_split = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, dark_current=40.0, response_model=NullFrameResponse(),
        sensor=hgcdte_glow_sensor)
    hgcdte_glow_whole_frame = copy(capture!(hgcdte_glow_whole,
        incremental_rate_input; rng=MersenneTwister(33)))
    capture!(hgcdte_glow_split, incremental_rate_input;
        rng=MersenneTwister(34), sample_time=0.5)
    hgcdte_glow_split_frame = copy(capture!(hgcdte_glow_split,
        incremental_rate_input; rng=MersenneTwister(35), sample_time=0.5))
    hgcdte_glow_expected = (40.0 + 60.0) * (1.0 + 2 * 0.1)
    @test mean(hgcdte_glow_whole_frame) ≈ hgcdte_glow_expected rtol=0.01
    @test mean(hgcdte_glow_split_frame) ≈ hgcdte_glow_expected rtol=0.01

    thermal_ingaas = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=NullFrameResponse(),
        thermal_model=FixedTemperature(temperature_K=250.0, glow_rate_law=linear),
        sensor=InGaAsSensor(glow_rate=2.0))
    @test supports_temperature_dependent_glow(thermal_ingaas)
    @test effective_glow_rate(thermal_ingaas) ≈ 1.0

    thermal_emccd = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=NullFrameResponse(),
        thermal_model=FixedTemperature(temperature_K=250.0, cic_rate_law=linear),
        sensor=EMCCDSensor(clock_induced_charge_per_frame=2.0))
    @test effective_cic_rate(thermal_emccd) ≈ 1.0

    det_saphira = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(avalanche_gain=5.0))
    frame_saphira = copy(capture!(det_saphira, uniform_signal; rng=MersenneTwister(14)))
    @test frame_saphira == 5.0 .* uniform_signal
    @test supports_avalanche_gain(det_saphira.params.sensor)
    @test supports_sensor_glow(det_saphira.params.sensor)
    @test detector_export_metadata(det_saphira).frame_response == :none
    @test detector_export_metadata(det_saphira).charge_coupling == :none
    saphira_impulse = zeros(5, 5)
    saphira_impulse[3, 3] = 100.0
    det_saphira_ipc = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, sensor=HgCdTeAvalancheArraySensor(),
        charge_coupling_model=InterpixelCapacitance(
            [0.0 0.01 0.0; 0.01 0.96 0.01; 0.0 0.01 0.0]))
    frame_saphira_ipc = capture!(det_saphira_ipc, saphira_impulse;
        rng=MersenneTwister(14))
    @test frame_saphira_ipc[3, 3] == 96.0
    @test frame_saphira_ipc[2, 3] == 1.0
    @test detector_export_metadata(det_saphira_ipc).charge_coupling ==
        :interpixel_capacitance
    det_saphira_excess = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(avalanche_gain=1.0, excess_noise_factor=sqrt(2.0)))
    frame_saphira_excess = copy(capture!(det_saphira_excess, uniform_signal; rng=MersenneTwister(14)))
    @test std(vec(frame_saphira_excess)) > 0
    det_saphira_sat = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, full_well=100.0, sensor=HgCdTeAvalancheArraySensor(avalanche_gain=5.0))
    frame_saphira_sat = copy(capture!(det_saphira_sat, uniform_signal; rng=MersenneTwister(15)))
    @test maximum(frame_saphira_sat) == 100.0
    saphira_saturated_stack = fill(50.0, 2, 8, 8)
    AdaptiveOpticsSim.capture_stack!(det_saphira_sat,
        saphira_saturated_stack, similar(saphira_saturated_stack);
        rng=MersenneTwister(15))
    @test all(==(100.0), saphira_saturated_stack)
    det_saphira_single = Detector(integration_time=1.0, noise=NoiseReadout(4.0), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor())
    frame_saphira_single = copy(capture!(det_saphira_single, zero_psf; rng=MersenneTwister(16)))
    single_products = readout_products(det_saphira_single)
    @test single_products isa HgCdTeReadoutProducts
    @test single_products isa MultiReadFrameReadoutProducts
    @test detector_reference_frame(det_saphira_single) === nothing
    @test detector_signal_frame(det_saphira_single) !== nothing
    @test detector_combined_frame(det_saphira_single) == frame_saphira_single
    @test detector_reference_cube(det_saphira_single) === nothing
    @test detector_signal_cube(det_saphira_single) !== nothing
    @test detector_read_cube(det_saphira_single) === nothing
    @test detector_read_times(det_saphira_single) === nothing

    @test fieldtype(typeof(det.state), :readout_products) === NoFrameReadoutProducts
    @test fieldtype(typeof(det_saphira_single.state), :readout_products) !== FrameReadoutProducts

    struct DummyReadoutProducts{A,V} <: FrameReadoutProducts
        signal_frame::A
        read_times::V
    end
    AdaptiveOpticsSim.detector_signal_frame(products::DummyReadoutProducts) = products.signal_frame
    AdaptiveOpticsSim.detector_read_times(products::DummyReadoutProducts) = products.read_times

    dummy_products = DummyReadoutProducts(fill(3.0, 2, 2), [0.25, 0.5])
    @test detector_reference_frame(dummy_products) === nothing
    @test detector_signal_frame(dummy_products) == fill(3.0, 2, 2)
    @test detector_combined_frame(dummy_products) === nothing
    @test detector_reference_cube(dummy_products) === nothing
    @test detector_signal_cube(dummy_products) === nothing
    @test detector_read_cube(dummy_products) === nothing
    @test detector_read_times(dummy_products) == [0.25, 0.5]

    det_saphira_ndr = Detector(integration_time=1.0, noise=NoiseReadout(4.0), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(sampling_mode=AveragedNonDestructiveReads(4)))
    frame_saphira_ndr = copy(capture!(det_saphira_ndr, zero_psf; rng=MersenneTwister(16)))
    @test std(vec(frame_saphira_ndr)) < std(vec(frame_saphira_single))
    @test supports_nondestructive_reads(det_saphira_ndr.params.sensor)
    @test supports_readout_correction(det_saphira_ndr.params.sensor)
    @test supports_read_cube(det_saphira_ndr.params.sensor)
    saphira_meta = detector_export_metadata(det_saphira_ndr)
    @test saphira_meta.sampling_mode == :averaged_non_destructive_reads
    @test saphira_meta.sampling_reads == 4
    @test saphira_meta.sampling_reference_reads == 0
    @test saphira_meta.sampling_signal_reads == 4
    @test saphira_meta.readout_sigma == 2.0
    @test saphira_meta.provides_signal_frame
    @test !saphira_meta.provides_reference_frame
    @test saphira_meta.provides_combined_frame
    @test !saphira_meta.provides_reference_cube
    @test saphira_meta.provides_signal_cube
    @test saphira_meta.provides_read_cube
    @test saphira_meta.signal_cube_reads == 4
    @test saphira_meta.read_cube_reads == 4
    @test detector_combined_frame(det_saphira_ndr) == frame_saphira_ndr
    @test size(detector_signal_cube(det_saphira_ndr)) == (4, 4, 4)
    @test length(detector_read_times(det_saphira_ndr)) == 4
    det_saphira_cds = Detector(integration_time=1.0, noise=NoiseReadout(4.0), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(sampling_mode=CorrelatedDoubleSampling()))
    frame_saphira_cds = copy(capture!(det_saphira_cds, zero_psf; rng=MersenneTwister(16)))
    @test std(vec(frame_saphira_cds)) > std(vec(frame_saphira_single))
    cds_meta = detector_export_metadata(det_saphira_cds)
    @test cds_meta.sampling_mode == :correlated_double_sampling
    @test cds_meta.sampling_reference_reads == 1
    @test cds_meta.sampling_signal_reads == 1
    @test cds_meta.readout_sigma ≈ 4.0 * sqrt(2.0)
    @test cds_meta.provides_reference_cube
    @test cds_meta.provides_signal_cube
    @test cds_meta.reference_cube_reads == 1
    @test cds_meta.signal_cube_reads == 1
    @test size(detector_reference_cube(det_saphira_cds)) == (4, 4, 1)
    @test size(detector_signal_cube(det_saphira_cds)) == (4, 4, 1)
    @test detector_combined_frame(det_saphira_cds) ≈
        detector_signal_frame(det_saphira_cds) .- detector_reference_frame(det_saphira_cds)
    det_saphira_fowler = Detector(integration_time=1.0, noise=NoiseReadout(4.0), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(sampling_mode=FowlerSampling(8)))
    frame_saphira_fowler = copy(capture!(det_saphira_fowler, zero_psf; rng=MersenneTwister(16)))
    fowler_meta = detector_export_metadata(det_saphira_fowler)
    @test fowler_meta.sampling_mode == :fowler_sampling
    @test fowler_meta.sampling_reads == 16
    @test fowler_meta.sampling_reference_reads == 8
    @test fowler_meta.sampling_signal_reads == 8
    @test fowler_meta.readout_sigma == 2.0
    @test fowler_meta.reference_cube_reads == 8
    @test fowler_meta.signal_cube_reads == 8
    @test detector_combined_frame(det_saphira_fowler) ≈
        detector_signal_frame(det_saphira_fowler) .- detector_reference_frame(det_saphira_fowler)

    ramp_input = fill(5.0, 4, 4)
    ramp_detector = Detector(integration_time=2.0, noise=NoiseNone(),
        qe=1.0, gain=1.0, response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(read_time=0.1,
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
        sensor=HgCdTeAvalancheArraySensor(read_time=0.1,
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
        sensor=HgCdTeAvalancheArraySensor(read_time=0.01,
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
    @test_throws InvalidConfiguration HgCdTeAvalancheArraySensor(
        sampling_mode=SkipperSampling(4))
    invalid_ramp_schedule = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0,
        sensor=HgCdTeAvalancheArraySensor(read_time=0.3,
            sampling_mode=UpTheRampSampling(5)))
    @test_throws InvalidConfiguration capture!(invalid_ramp_schedule,
        ones(4, 4); rng=MersenneTwister(164))

    multiread_noise_fixture = zeros(64, 64)
    single_std = std(vec(copy(capture!(det_saphira_single, multiread_noise_fixture;
        rng=MersenneTwister(160)))))
    ndr_std = std(vec(copy(capture!(det_saphira_ndr, multiread_noise_fixture;
        rng=MersenneTwister(160)))))
    cds_std = std(vec(copy(capture!(det_saphira_cds, multiread_noise_fixture;
        rng=MersenneTwister(160)))))
    fowler_std = std(vec(copy(capture!(det_saphira_fowler, multiread_noise_fixture;
        rng=MersenneTwister(160)))))
    @test isapprox(single_std, 4.0; rtol=0.15)
    @test isapprox(ndr_std, 2.0; rtol=0.15)
    @test isapprox(cds_std, 4.0 * sqrt(2.0); rtol=0.15)
    @test isapprox(fowler_std, 2.0; rtol=0.15)
    det_saphira_timed_single = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        dark_current=1000.0, gain=1.0, sensor=HgCdTeAvalancheArraySensor(read_time=1.0))
    frame_saphira_timed_single = copy(capture!(det_saphira_timed_single, zero_psf; rng=MersenneTwister(17)))
    det_saphira_timed_cds = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        dark_current=1000.0, gain=1.0,
        sensor=HgCdTeAvalancheArraySensor(read_time=1.0, sampling_mode=CorrelatedDoubleSampling()))
    frame_saphira_timed_cds = copy(capture!(det_saphira_timed_cds, zero_psf; rng=MersenneTwister(17)))
    @test sum(frame_saphira_timed_cds) > sum(frame_saphira_timed_single)
    det_saphira_timed_glow = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        dark_current=3.0, gain=1.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(glow_rate=2.0, read_time=1.0, sampling_mode=CorrelatedDoubleSampling()))
    frame_saphira_timed_glow = copy(capture!(det_saphira_timed_glow, zero_psf; rng=MersenneTwister(17)))
    det_saphira_timed_noglow = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        dark_current=3.0, gain=1.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(read_time=1.0, sampling_mode=CorrelatedDoubleSampling()))
    frame_saphira_timed_noglow = copy(capture!(det_saphira_timed_noglow, zero_psf; rng=MersenneTwister(17)))
    @test sum(frame_saphira_timed_glow) > sum(frame_saphira_timed_noglow)
    timed_meta = detector_export_metadata(det_saphira_timed_cds)
    @test timed_meta.sampling_read_time == 1.0
    @test timed_meta.sampling_wallclock_time == 3.0
    det_saphira_windowed = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(read_time=1.0, sampling_mode=CorrelatedDoubleSampling()),
        readout_window=FrameWindow(2:3, 2:3))
    capture!(det_saphira_windowed, fill(10.0, 4, 4); rng=MersenneTwister(18))
    windowed_meta = detector_export_metadata(det_saphira_windowed)
    @test windowed_meta.sampling_read_time == 0.5
    @test windowed_meta.sampling_wallclock_time == 2.0
    @test detector_combined_frame(det_saphira_windowed) !== nothing
    @test size(detector_reference_cube(det_saphira_windowed)) == (2, 2, 1)
    @test size(detector_signal_cube(det_saphira_windowed)) == (2, 2, 1)
    @test size(detector_signal_frame(det_saphira_windowed)) == (2, 2)
    @test size(detector_read_cube(det_saphira_windowed)) == (2, 2, 2)
    @test detector_read_times(det_saphira_windowed) == [0.5, 1.0]
    @test detector_combined_frame(det_saphira_windowed) ≈
        detector_signal_frame(det_saphira_windowed) .- detector_reference_frame(det_saphira_windowed)
    det_saphira_windowed_corrected = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(read_time=1.0, sampling_mode=CorrelatedDoubleSampling()),
        readout_window=FrameWindow(2:3, 2:3),
        correction_model=CompositeFrameReadoutCorrection((
            ReferenceRowCommonModeCorrection(1),
            ReferenceColumnCommonModeCorrection(1),
        )))
    row_pattern = repeat(reshape([1.0, 2.0, 3.0, 4.0], :, 1), 1, 4)
    col_pattern = repeat(reshape([1.0, 2.0, 3.0, 4.0], 1, :), 4, 1)
    windowed_corrected_input = row_pattern .+ col_pattern
    windowed_corrected = capture!(det_saphira_windowed_corrected, windowed_corrected_input; rng=MersenneTwister(18))
    @test maximum(abs, windowed_corrected) < 1e-6
    @test detector_combined_frame(det_saphira_windowed_corrected) ≈
        detector_signal_frame(det_saphira_windowed_corrected) .- detector_reference_frame(det_saphira_windowed_corrected)
    @test size(detector_read_cube(det_saphira_windowed_corrected)) == (2, 2, 2)
    det_saphira_corrected = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(),
        response_model=NullFrameResponse(),
        correction_model=ReferencePixelCommonModeCorrection(1, 1))
    corrected_frame = capture!(det_saphira_corrected, fill(5.0, 4, 4); rng=MersenneTwister(19))
    corrected_meta = detector_export_metadata(det_saphira_corrected)
    @test corrected_meta.readout_correction == :reference_pixel_common_mode
    @test corrected_meta.correction_edge_rows == 1
    @test corrected_meta.correction_edge_cols == 1
    @test abs(mean(corrected_frame)) < 1e-6

    calibration_input = reshape(collect(1.0:64.0), 8, 8)
    calibration_sensor = HgCdTeAvalancheArraySensor(avalanche_gain=2.0,
        read_time=0.1, sampling_mode=CorrelatedDoubleSampling())
    calibration_correction = ReferencePixelCommonModeCorrection(1, 1)
    calibration_detector = Detector(noise=NoiseNone(), qe=1.0, gain=3.0,
        sensor=calibration_sensor, response_model=NullFrameResponse(),
        correction_model=calibration_correction)
    capture_detector = Detector(noise=NoiseNone(), qe=1.0, gain=3.0,
        sensor=calibration_sensor, response_model=NullFrameResponse(),
        correction_model=calibration_correction)
    deterministic_reference = copy(
        AdaptiveOpticsSim.detector_calibration_frame!(
            calibration_detector, calibration_input, 1.0))
    noiseless_cds = copy(capture!(capture_detector, calibration_input;
        rng=MersenneTwister(1901)))
    @test deterministic_reference == noiseless_cds
    @test iszero(calibration_detector.state.integrated_time)
    @test readout_ready(calibration_detector)
    @test readout_products(calibration_detector) isa NoFrameReadoutProducts

    ramp_calibration_sensor = HgCdTeAvalancheArraySensor(
        avalanche_gain=2.0, read_time=0.1,
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
        AdaptiveOpticsSim.detector_calibration_frame!(
            ramp_calibration_detector, calibration_input, 1.0))
    noiseless_ramp = copy(capture!(ramp_capture_detector,
        calibration_input; rng=MersenneTwister(1902)))
    @test ramp_reference ≈ noiseless_ramp atol=1e-12 rtol=1e-12
    @test iszero(ramp_calibration_detector.state.integrated_time)
    @test readout_ready(ramp_calibration_detector)
    @test readout_products(ramp_calibration_detector) isa NoFrameReadoutProducts

    invalid_ramp_calibration = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, gain=3.0,
        sensor=HgCdTeAvalancheArraySensor(avalanche_gain=2.0,
            read_time=0.4, sampling_mode=UpTheRampSampling(4)),
        response_model=NullFrameResponse(),
        correction_model=calibration_correction)
    @test_throws InvalidConfiguration begin
        AdaptiveOpticsSim.detector_calibration_signature(
            invalid_ramp_calibration, UInt(0))
    end
    @test_throws InvalidConfiguration begin
        AdaptiveOpticsSim.detector_calibration_frame!(
            invalid_ramp_calibration, calibration_input, 1.0)
    end

    det_saphira_row_corrected = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(),
        response_model=NullFrameResponse(),
        correction_model=ReferenceRowCommonModeCorrection(1))
    row_corrected = capture!(det_saphira_row_corrected, row_pattern; rng=MersenneTwister(20))
    @test maximum(abs, row_corrected) < 1e-6
    det_saphira_col_corrected = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(),
        response_model=NullFrameResponse(),
        correction_model=ReferenceColumnCommonModeCorrection(1))
    col_corrected = capture!(det_saphira_col_corrected, col_pattern; rng=MersenneTwister(21))
    @test maximum(abs, col_corrected) < 1e-6
    det_saphira_output_corrected = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(),
        response_model=NullFrameResponse(),
        correction_model=ReferenceOutputCommonModeCorrection(2; edge_rows=1, edge_cols=1))
    output_pattern = hcat(fill(5.0, 4, 2), fill(10.0, 4, 2))
    output_corrected = capture!(det_saphira_output_corrected, output_pattern; rng=MersenneTwister(22))
    output_meta = detector_export_metadata(det_saphira_output_corrected)
    @test output_meta.readout_correction == :reference_output_common_mode
    @test output_meta.correction_group_cols == 2
    @test maximum(abs, output_corrected) < 1e-6
    det_saphira_composite = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(),
        response_model=NullFrameResponse(),
        correction_model=CompositeFrameReadoutCorrection((
            ReferenceRowCommonModeCorrection(1),
            ReferenceColumnCommonModeCorrection(1))))
    composite_pattern = row_pattern .+ col_pattern
    composite_corrected = capture!(det_saphira_composite, composite_pattern; rng=MersenneTwister(23))
    composite_meta = detector_export_metadata(det_saphira_composite)
    @test composite_meta.readout_correction == :composite
    @test composite_meta.correction_stage_count == 2
    @test maximum(abs, composite_corrected) < 1e-6
    @test_throws InvalidConfiguration HgCdTeAvalancheArraySensor(avalanche_gain=0.5)
    @test_throws InvalidConfiguration HgCdTeAvalancheArraySensor(excess_noise_factor=0.5)
    @test_throws InvalidConfiguration HgCdTeAvalancheArraySensor(glow_rate=-1.0)
    @test_throws InvalidConfiguration HgCdTeAvalancheArraySensor(read_time=-1.0)
    @test_throws InvalidConfiguration HgCdTeAvalancheArraySensor(sampling_mode=AveragedNonDestructiveReads(0))
    @test_throws InvalidConfiguration HgCdTeAvalancheArraySensor(sampling_mode=FowlerSampling(0))
    @test_throws InvalidConfiguration ReferencePixelCommonModeCorrection(0, 0)
    @test_throws InvalidConfiguration ReferenceRowCommonModeCorrection(0)
    @test_throws InvalidConfiguration ReferenceColumnCommonModeCorrection(0)
    @test_throws InvalidConfiguration ReferenceOutputCommonModeCorrection(0)
    @test_throws InvalidConfiguration CompositeFrameReadoutCorrection(())

    @test_throws InvalidConfiguration Detector(integration_time=1.0, noise=NoisePhoton(), qe=1.0, binning=1,
        sensor=APDSensor())

    det_default_ccd = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1, sensor=CCDSensor())
    @test detector_export_metadata(det_default_ccd).frame_response == :none

    apd = APDDetector(integration_time=1.0, qe=0.5, gain=2.0, dark_count_rate=0.0, noise=NoiseNone())
    channels = fill(4.0, 2, 8)
    apd_out = capture!(apd, channels; rng=MersenneTwister(9))
    @test apd_out == fill(4.0, 2, 8)
    @test channel_output(apd) === apd_out
    apd_meta = detector_export_metadata(apd)
    @test apd_meta isa CountingDetectorExportMetadata
    @test apd_meta.sensor == :apd
    @test apd_meta.readout.output_size == (2, 8)
    @test apd_meta.readout.n_channels == 16
    @test apd_meta.dead_time_model == :none
    @test apd_meta.dead_time === nothing
    @test supports_counting_noise(apd)
    @test !supports_dead_time(apd)
    @test !supports_channel_gain_map(apd)

    apd_gain_map = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), channel_gain_map=fill(0.5, 2, 8))
    @test capture!(apd_gain_map, fill(2.0, 2, 8); rng=MersenneTwister(9)) == fill(1.0, 2, 8)
    @test supports_channel_gain_map(apd_gain_map)

    linear_apd = LinearAPDDetector(integration_time=0.5, qe=0.5,
        avalanche_gain=4.0, conversion_gain=2.0, noise=NoiseNone())
    linear_apd_out = capture!(linear_apd, 10.0; rng=MersenneTwister(90))
    @test linear_apd_out == [20.0]
    @test channel_output(linear_apd) === linear_apd_out
    @test supports_avalanche_gain(linear_apd)
    linear_apd_meta = detector_export_metadata(linear_apd)
    @test linear_apd_meta.topology == :single_element
    @test linear_apd_meta.n_channels == 1
    linear_apd_rng = MersenneTwister(90)
    capture!(linear_apd, 10.0; rng=linear_apd_rng)
    @test @allocated(capture!(linear_apd, 10.0;
        rng=linear_apd_rng)) == 0

    linear_apd_bank = LinearAPDDetector(topology=APDChannelBank(4),
        integration_time=1.0, qe=0.5, avalanche_gain=2.0,
        dark_current=1.0, noise=NoiseNone())
    @test capture!(linear_apd_bank, fill(3.0, 4);
        rng=MersenneTwister(91)) == fill(5.0, 4)
    @test detector_export_metadata(linear_apd_bank).topology == :channel_bank
    @test_throws DimensionMismatchError capture!(linear_apd_bank, 3.0)
    @test_throws DimensionMismatchError capture!(linear_apd_bank, fill(3.0, 3))

    linear_apd_noisy = LinearAPDDetector(topology=APDChannelBank(4096),
        noise=NoisePhotonReadout(2.0), avalanche_gain=3.0,
        excess_noise_factor=1.3)
    noisy_apd_out = capture!(linear_apd_noisy, fill(20.0, 4096);
        rng=MersenneTwister(92))
    @test std(noisy_apd_out) > 0
    @test isapprox(mean(noisy_apd_out), 60.0; rtol=0.05)

    @test_throws InvalidConfiguration APDChannelBank(1)
    @test_throws InvalidConfiguration LinearAPDDetector(integration_time=0.0)
    @test_throws InvalidConfiguration LinearAPDDetector(qe=1.1)
    @test_throws InvalidConfiguration LinearAPDDetector(avalanche_gain=0.5)
    @test_throws InvalidConfiguration LinearAPDDetector(excess_noise_factor=0.5)
    @test_throws InvalidConfiguration LinearAPDDetector(dark_current=-1.0)
    @test_throws InvalidConfiguration LinearAPDDetector(conversion_gain=0.0)

    spad_sensor = SPADArraySensor(pde=0.5, dark_count_rate=0.0, fill_factor=0.8)
    spad = SPADArrayDetector(integration_time=1.0, noise=NoiseNone(), sensor=spad_sensor)
    spad_out = capture!(spad, fill(10.0, 2, 8); rng=MersenneTwister(9))
    @test spad_out == fill(4.0, 2, 8)
    @test output_frame(spad) === spad_out
    @test channel_output(spad) === spad_out
    spad_meta = detector_export_metadata(spad)
    @test spad_meta isa CountingDetectorExportMetadata
    @test spad_meta.sensor == :spad_array
    @test spad_meta.qe == 0.5
    @test spad_meta.fill_factor == 0.8
    @test spad_meta.gain == 1.0
    @test spad_meta.readout.output_size == (2, 8)
    @test spad_meta.readout.n_channels == 16
    @test !supports_channel_gain_map(spad)
    @test supports_counting_noise(spad)
    @test !supports_dead_time(spad)

    spad_dead = SPADArrayDetector(
        integration_time=1.0,
        noise=NoiseNone(),
        sensor=SPADArraySensor(pde=0.5, dark_count_rate=0.0, fill_factor=0.8,
            dead_time_model=NonParalyzableDeadTime(0.5)),
    )
    @test capture!(spad_dead, fill(10.0, 2, 8); rng=MersenneTwister(9)) ≈ fill(4.0 / 3.0, 2, 8)
    @test supports_dead_time(spad_dead)

    spad_gate = SPADArrayDetector(
        integration_time=1.0,
        noise=NoiseNone(),
        gate_model=DutyCycleGate(0.5),
        sensor=SPADArraySensor(pde=0.5, dark_count_rate=0.0, fill_factor=0.8),
    )
    @test capture!(spad_gate, fill(10.0, 2, 8); rng=MersenneTwister(9)) == fill(2.0, 2, 8)
    @test supports_counting_gating(spad_gate)

    spad_afterpulse = SPADArrayDetector(
        integration_time=1.0,
        noise=NoiseNone(),
        sensor=SPADArraySensor(pde=1.0, dark_count_rate=0.0, fill_factor=1.0,
            correlation_model=AfterpulsingModel(0.25)),
    )
    @test capture!(spad_afterpulse, fill(4.0, 2, 8); rng=MersenneTwister(9)) == fill(5.0, 2, 8)
    @test supports_afterpulsing(spad_afterpulse)

    spad_crosstalk = SPADArrayDetector(
        integration_time=1.0,
        noise=NoiseNone(),
        sensor=SPADArraySensor(pde=1.0, dark_count_rate=0.0, fill_factor=1.0,
            correlation_model=ChannelCrosstalkModel(0.4)),
    )
    spad_crosstalk_in = zeros(3, 3)
    spad_crosstalk_in[2, 2] = 10.0
    spad_crosstalk_out = capture!(spad_crosstalk, spad_crosstalk_in; rng=MersenneTwister(9))
    @test spad_crosstalk_out[2, 2] ≈ 6.0
    @test spad_crosstalk_out[1, 2] ≈ 1.0
    @test spad_crosstalk_out[2, 1] ≈ 1.0
    @test spad_crosstalk_out[2, 3] ≈ 1.0
    @test spad_crosstalk_out[3, 2] ≈ 1.0
    @test sum(spad_crosstalk_out) ≈ sum(spad_crosstalk_in)
    @test supports_channel_crosstalk(spad_crosstalk)

    spad_single = SPADArrayDetector(
        integration_time=1.0,
        noise=NoiseNone(),
        sensor=SPADArraySensor(pde=1.0, dark_count_rate=0.0, fill_factor=1.0,
            correlation_model=ChannelCrosstalkModel(0.4)),
    )
    @test capture!(spad_single, fill(10.0, 1, 1), MersenneTwister(9)) == fill(10.0, 1, 1)

    spad_dynamic = SPADArrayDetector(
        integration_time=1.0,
        noise=NoiseNone(),
        sensor=SPADArraySensor(pde=1.0, dark_count_rate=10.0, fill_factor=1.0),
        thermal_model=FixedTemperature(temperature_K=80.0, dark_count_law=arrhenius),
    )
    spad_dynamic_meta = detector_export_metadata(spad_dynamic)
    @test supports_detector_thermal_model(spad_dynamic)
    @test supports_temperature_dependent_dark_counts(spad_dynamic)
    @test detector_temperature(spad_dynamic) == 80.0
    @test spad_dynamic_meta.thermal_model == :fixed_temperature
    @test spad_dynamic_meta.dark_count_law == :arrhenius
    @test effective_dark_count_rate(spad_dynamic) < spad_dynamic.params.sensor.dark_count_rate

    @test_throws InvalidConfiguration SPADArraySensor(pde=1.5)
    @test_throws InvalidConfiguration SPADArraySensor(dark_count_rate=-1.0)
    @test_throws InvalidConfiguration SPADArraySensor(fill_factor=0.0)
    @test_throws InvalidConfiguration SPADArrayDetector(noise=NoiseReadout(1.0))

    mkid_sensor = MKIDArraySensor(qe=0.7, dark_count_rate=0.0, fill_factor=0.9,
        energy_resolution=12.0, timing_jitter_s=2e-6, wavelength_range_m=(0.8e-6, 1.4e-6))
    mkid = MKIDArrayDetector(integration_time=2.0, noise=NoiseNone(), sensor=mkid_sensor,
        output_type=UInt16)
    mkid_out = capture!(mkid, fill(10.0, 2, 8); rng=MersenneTwister(9))
    @test mkid_out == fill(UInt16(13), 2, 8)
    @test output_frame(mkid) === mkid_out
    @test supports_photon_number_resolving(mkid.params.sensor)
    @test supports_energy_resolving(mkid.params.sensor)
    @test !supports_dead_time(mkid)
    mkid_meta = detector_export_metadata(mkid)
    @test mkid_meta isa CountingDetectorExportMetadata
    @test mkid_meta.sensor == :mkid_array
    @test mkid_meta.qe == 0.7
    @test mkid_meta.fill_factor == 0.9
    @test mkid_meta.energy_resolution == 12.0
    @test mkid_meta.timing_jitter_s == 2e-6
    @test mkid_meta.wavelength_min_m == 0.8e-6
    @test mkid_meta.wavelength_max_m == 1.4e-6
    @test mkid_meta.readout.output_size == (2, 8)
    @test mkid_meta.readout.n_channels == 16

    mkid_band = MKIDArrayDetector(
        integration_time=1.0,
        noise=NoiseNone(),
        sensor=MKIDArraySensor(qe=1.0, dark_count_rate=0.0, fill_factor=1.0,
            wavelength_range_m=(0.8e-6, 1.4e-6)),
    )
    inside_band = Source(band=:custom, wavelength=1.0e-6)
    outside_band = Source(band=:custom, wavelength=0.55e-6)
    @test capture!(mkid_band, fill(2.0, 2, 2), inside_band, MersenneTwister(10)) == fill(2.0, 2, 2)
    @test capture!(mkid_band, fill(2.0, 2, 2),
        Source(band=:custom, wavelength=0.8e-6), MersenneTwister(10)) == fill(2.0, 2, 2)
    @test capture!(mkid_band, fill(2.0, 2, 2),
        Source(band=:custom, wavelength=1.4e-6), MersenneTwister(10)) == fill(2.0, 2, 2)
    @test capture!(mkid_band, fill(2.0, 2, 2), outside_band, MersenneTwister(10)) == zeros(2, 2)
    spectral_band = with_spectrum(
        inside_band,
        SpectralBundle([0.6e-6, 1.0e-6, 1.6e-6], [0.2, 0.3, 0.5]),
    )
    @test capture!(mkid_band, fill(10.0, 2, 2), spectral_band, MersenneTwister(10)) ≈ fill(3.0, 2, 2)

    mkid_dead = MKIDArrayDetector(
        integration_time=1.0,
        noise=NoiseNone(),
        sensor=MKIDArraySensor(qe=1.0, dark_count_rate=0.0, fill_factor=1.0,
            dead_time_model=ParalyzableDeadTime(0.5)),
    )
    @test capture!(mkid_dead, fill(4.0, 2, 8); rng=MersenneTwister(9)) ≈ fill(4.0 * exp(-2.0), 2, 8)
    @test supports_dead_time(mkid_dead)
    @test supports_paralyzable_dead_time(mkid_dead)

    mkid_gate = MKIDArrayDetector(
        integration_time=1.0,
        noise=NoiseNone(),
        gate_model=DutyCycleGate(0.25),
        sensor=MKIDArraySensor(qe=1.0, dark_count_rate=0.0, fill_factor=1.0),
    )
    @test capture!(mkid_gate, fill(8.0, 2, 8); rng=MersenneTwister(9)) == fill(2.0, 2, 8)
    @test supports_counting_gating(mkid_gate)

    @test_throws InvalidConfiguration MKIDArraySensor(qe=1.5)
    @test_throws InvalidConfiguration MKIDArraySensor(dark_count_rate=-1.0)
    @test_throws InvalidConfiguration MKIDArraySensor(fill_factor=0.0)
    @test_throws InvalidConfiguration MKIDArraySensor(energy_resolution=0.0)
    @test_throws InvalidConfiguration MKIDArraySensor(timing_jitter_s=-1.0)
    @test_throws InvalidConfiguration MKIDArraySensor(wavelength_range_m=(1.4e-6, 0.8e-6))
    @test_throws InvalidConfiguration MKIDArraySensor(wavelength_range_m=(NaN, 1.4e-6))
    @test_throws InvalidConfiguration MKIDArraySensor(energy_resolution=Inf)
    @test_throws InvalidConfiguration MKIDArraySensor(timing_jitter_s=NaN)
    @test_throws InvalidConfiguration MKIDArrayDetector(noise=NoiseReadout(1.0))

    detector_artifact_path = normpath(joinpath(@__DIR__, "..", "..", "benchmarks", "results",
        "detectors", "2026-07-12-detector-mkid-validation.toml"))
    @test isfile(detector_artifact_path)
    detector_artifact = TOML.parsefile(detector_artifact_path)
    @test all(values(detector_artifact["interpretation"]))
    @test issubset(Set(["apd", "spad_array", "mkid_array"]),
        Set(detector_artifact["scope"]["families"]))
    detector_manifest = TOML.parsefile(joinpath(dirname(detector_artifact_path),
        "manifest.toml"))
    detector_entries = Dict(entry["id"] => entry for entry in
        detector_manifest["artifacts"])
    @test detector_entries["DET-VAL-2026-07-12"]["status"] == "active"
    @test detector_entries["DET-VAL-2026-04-23"]["status"] == "superseded"
    @test detector_entries["DET-VAL-2026-04-23"]["superseded_by"] ==
        "DET-VAL-2026-07-12"

    function impulse_transfer_magnitude(frame, spatial_frequency_x,
        spatial_frequency_y)
        center_i = fld(size(frame, 1), 2) + 1
        center_j = fld(size(frame, 2), 2) + 1
        response = zero(Complex{eltype(frame)})
        @inbounds for j in axes(frame, 2), i in axes(frame, 1)
            phase = -2pi * (spatial_frequency_y * (i - center_i) +
                            spatial_frequency_x * (j - center_j))
            response += frame[i, j] * cis(phase)
        end
        return abs(response) / sum(frame)
    end

    det_mtf = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=GaussianPixelResponse(response_width_px=0.75))
    impulse = zeros(9, 9)
    impulse[5, 5] = 1.0
    frame_mtf = capture!(det_mtf, impulse; rng=MersenneTwister(3))
    @test sum(frame_mtf) ≈ 1.0 atol=1e-6
    @test frame_mtf[5, 5] < 1.0
    @test frame_mtf[5, 4] > 0
    @test supports_detector_mtf(det_mtf)
    @test detector_mtf(det_mtf, 0.0, 0.0) ≈ 1.0
    @test detector_mtf(det_mtf, 0.5, 0.0) < 1.0
    @test detector_mtf(det_mtf, 0.5, 0.0) ≈
        impulse_transfer_magnitude(frame_mtf, 0.5, 0.0)
    mtf_meta = detector_export_metadata(det_mtf)
    @test mtf_meta.frame_response == :gaussian
    @test mtf_meta.response_width_px == 0.75
    @test mtf_meta.response_application_domain == :image
    @test mtf_meta.response_is_separable
    @test !mtf_meta.response_is_shift_invariant
    @test mtf_meta.response_support_rows == mtf_meta.response_support_cols
    @test mtf_meta.pitch_x_px === nothing
    @test mtf_meta.aperture_shape === nothing
    @test_throws InvalidConfiguration GaussianPixelResponse(response_width_px=0.0)

    sampled_kernel = [0.0 0.1 0.0; 0.1 0.6 0.1; 0.0 0.1 0.0]
    sampled_det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=SampledFrameResponse(sampled_kernel))
    sampled_frame = capture!(sampled_det, impulse; rng=MersenneTwister(6))
    @test sum(sampled_frame) ≈ 1.0 atol=1e-6
    @test sampled_frame[5, 5] ≈ 0.6 atol=1e-6
    @test sampled_frame[5, 4] ≈ 0.1 atol=1e-6
    sampled_meta = detector_export_metadata(sampled_det)
    @test sampled_meta.frame_response == :sampled
    @test sampled_meta.response_application_domain == :image
    @test !sampled_meta.response_is_separable
    @test sampled_meta.response_support_rows == 3
    @test sampled_meta.response_support_cols == 3
    @test sampled_meta.aperture_shape == :sampled
    @test supports_detector_mtf(sampled_det)
    @test detector_mtf(sampled_det, 0.0, 0.0) ≈ 1.0
    @test_throws InvalidConfiguration SampledFrameResponse(zeros(3, 3))
    @test_throws InvalidConfiguration SampledFrameResponse(ones(2, 3))
    negative_kernel = copy(sampled_kernel)
    negative_kernel[1, 1] = -0.1
    @test_throws InvalidConfiguration SampledFrameResponse(negative_kernel)
    nan_kernel = copy(sampled_kernel)
    nan_kernel[1, 1] = NaN
    @test_throws InvalidConfiguration SampledFrameResponse(nan_kernel)
    infinite_kernel = copy(sampled_kernel)
    infinite_kernel[1, 1] = Inf
    @test_throws InvalidConfiguration SampledFrameResponse(infinite_kernel)
    @test_throws InvalidConfiguration SampledFrameResponse(
        fill(0.2, 3, 3); normalize=false)

    asymmetric_kernel = [0.0 0.0 0.0; 0.1 0.2 0.7; 0.0 0.0 0.0]
    asymmetric_det = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, binning=1,
        response_model=SampledFrameResponse(asymmetric_kernel))
    center_impulse = zeros(9, 9)
    center_impulse[5, 5] = 1.0
    center_frame = copy(capture!(asymmetric_det, center_impulse;
        rng=MersenneTwister(61)))
    left_impulse = zeros(9, 9)
    left_impulse[5, 1] = 1.0
    left_frame = copy(capture!(asymmetric_det, left_impulse;
        rng=MersenneTwister(62)))
    right_impulse = zeros(9, 9)
    right_impulse[5, end] = 1.0
    right_frame = copy(capture!(asymmetric_det, right_impulse;
        rng=MersenneTwister(63)))
    @test sum(center_frame) ≈ 1.0
    @test sum(left_frame) ≈ 0.3
    @test sum(right_frame) ≈ 0.9
    @test sum(left_frame) <= sum(left_impulse)
    @test sum(right_frame) <= sum(right_impulse)
    @test minimum(left_frame) >= 0
    @test minimum(right_frame) >= 0
    @test detector_mtf(asymmetric_det, 0.37, 0.0) ≈
        impulse_transfer_magnitude(center_frame, 0.37, 0.0)
    @test !detector_export_metadata(asymmetric_det).response_is_shift_invariant

    cube_mtf = Array{Float64}(undef, 2, size(impulse, 1), size(impulse, 2))
    cube_mtf[1, :, :] .= impulse
    cube_mtf[2, :, :] .= impulse
    scratch_mtf = similar(cube_mtf)
    stack_mtf = AdaptiveOpticsSim.capture_stack!(det_mtf, cube_mtf, scratch_mtf; rng=MersenneTwister(10))
    @test size(stack_mtf) == size(cube_mtf)
    @test all(isfinite, stack_mtf)
    cube_sampled = Array{Float64}(undef, 2, size(impulse, 1), size(impulse, 2))
    cube_sampled[1, :, :] .= impulse
    cube_sampled[2, :, :] .= impulse
    scratch_sampled = similar(cube_sampled)
    stack_sampled = AdaptiveOpticsSim.capture_stack!(sampled_det, cube_sampled, scratch_sampled; rng=MersenneTwister(10))
    @test size(stack_sampled) == size(cube_sampled)
    @test all(isfinite, stack_sampled)
    @test stack_sampled[1, :, :] ≈ sampled_frame atol=1e-6
    @test stack_sampled[2, :, :] ≈ sampled_frame atol=1e-6
    det_stack_adc = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        bits=8, full_well=10.0, output_type=UInt16)
    cube_stack_adc = fill(10.0, 2, 4, 4)
    stack_adc = AdaptiveOpticsSim.capture_stack!(det_stack_adc, cube_stack_adc, similar(cube_stack_adc);
        rng=MersenneTwister(10))
    @test stack_adc === cube_stack_adc
    @test all(stack_adc .== 255.0)

    rect_det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=RectangularPixelAperture(pitch_x_px=2.0, pitch_y_px=2.0,
            fill_factor_x=0.6, fill_factor_y=0.8))
    rect_frame = capture!(rect_det, impulse; rng=MersenneTwister(4))
    @test sum(rect_frame) ≈ 1.0 atol=1e-6
    @test rect_frame[5, 5] < 1.0
    rect_meta = detector_export_metadata(rect_det)
    @test rect_meta.frame_response == :rectangular_aperture
    @test rect_meta.pitch_x_px == 2.0
    @test rect_meta.pitch_y_px == 2.0
    @test rect_meta.fill_factor_x == 0.6
    @test rect_meta.fill_factor_y == 0.8
    @test rect_meta.aperture_shape == :rectangular
    @test rect_meta.response_application_domain == :image
    @test supports_detector_mtf(rect_det)
    @test detector_mtf(rect_det, 0.5, 0.0) ≈
        impulse_transfer_magnitude(rect_frame, 0.5, 0.0)
    unit_rectangular = RectangularPixelAperture()
    half_fill_rectangular = RectangularPixelAperture(fill_factor_x=0.5,
        fill_factor_y=0.5)
    @test detector_mtf(unit_rectangular, 0.5, 0.0) ≈ 1.0
    @test unit_rectangular.kernel_x == half_fill_rectangular.kernel_x
    @test unit_rectangular.kernel_y == half_fill_rectangular.kernel_y
    @test detector_mtf(unit_rectangular, 0.5, 0.0) ==
        detector_mtf(half_fill_rectangular, 0.5, 0.0)
    @test !AdaptiveOpticsSim.supports_subpixel_geometry(unit_rectangular)
    @test !AdaptiveOpticsSim.supports_subpixel_geometry(
        half_fill_rectangular)

    @test_throws InvalidConfiguration RectangularPixelAperture(fill_factor_x=0.0)
    @test_throws InvalidConfiguration RectangularPixelAperture(fill_factor_y=1.5)
    @test_throws InvalidConfiguration RectangularPixelAperture(pitch_x_px=Inf)
    @test_throws InvalidConfiguration RectangularPixelAperture(fill_factor_y=NaN)
    @test_throws InvalidConfiguration GaussianPixelResponse(response_width_px=Inf)
    @test_throws InvalidConfiguration GaussianPixelResponse(truncate_at=Inf)

    ipc_kernel = [0.0 0.01 0.0; 0.01 0.96 0.01; 0.0 0.01 0.0]
    ipc_det = Detector(integration_time=1.0, noise=NoisePhoton(), qe=1.0,
        response_model=NullFrameResponse(),
        charge_coupling_model=InterpixelCapacitance(ipc_kernel))
    ipc_input = zeros(9, 9)
    ipc_input[5, 5] = 100.0
    ipc_frame = capture!(ipc_det, ipc_input; rng=MersenneTwister(45))
    @test ipc_frame[5, 4] > 0
    @test !isinteger(ipc_frame[5, 4])
    @test sum(ipc_frame) ≈ round(sum(ipc_frame)) atol=1e-10
    ipc_meta = detector_export_metadata(ipc_det)
    @test ipc_meta.charge_coupling == :interpixel_capacitance
    @test ipc_meta.charge_coupling_support_rows == 3
    @test ipc_meta.charge_coupling_support_cols == 3

    det_window_stack = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        readout_window=FrameWindow(2:8, 2:8))
    cube_window = Array{Float64}(undef, 2, size(impulse, 1), size(impulse, 2))
    cube_window[1, :, :] .= impulse
    cube_window[2, :, :] .= impulse
    scratch_window = similar(cube_window)
    @test_throws InvalidConfiguration AdaptiveOpticsSim.capture_stack!(det_window_stack, cube_window, scratch_window; rng=MersenneTwister(10))
    input_window_stack = copy(cube_window)
    output_window_stack = Array{Float64}(undef, 2, 7, 7)
    generalized_window = AdaptiveOpticsSim.capture_stack!(det_window_stack, output_window_stack, input_window_stack; rng=MersenneTwister(10))
    @test size(generalized_window) == (2, 7, 7)
    @test generalized_window[1, :, :] ≈ capture!(det_window_stack, impulse; rng=MersenneTwister(10))

    corrected_stack_models = (
        ReferencePixelCommonModeCorrection(1, 1),
        ReferenceRowCommonModeCorrection(1),
        ReferenceColumnCommonModeCorrection(1),
        ReferenceOutputCommonModeCorrection(2; edge_rows=1, edge_cols=1),
        CompositeFrameReadoutCorrection((
            ReferenceRowCommonModeCorrection(1),
            ReferenceColumnCommonModeCorrection(1),
        )),
    )
    for correction_model in corrected_stack_models
        corrected_stack_det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
            sensor=HgCdTeAvalancheArraySensor(sampling_mode=SingleRead()),
            correction_model=correction_model)
        corrected_stack_in = Array{Float64}(undef, 2, 5, 5)
        corrected_stack_in[1, :, :] .= reshape(collect(1.0:25.0), 5, 5)
        corrected_stack_in[2, :, :] .= reshape(collect(26.0:50.0), 5, 5)
        corrected_stack_ref = copy(corrected_stack_in)
        corrected_stack = AdaptiveOpticsSim.capture_stack!(corrected_stack_det, corrected_stack_in,
            similar(corrected_stack_in); rng=MersenneTwister(10))
        @test size(corrected_stack) == size(corrected_stack_in)
        corrected_frame_1 = capture!(Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
                sensor=HgCdTeAvalancheArraySensor(sampling_mode=SingleRead()),
                correction_model=correction_model),
            @view(corrected_stack_ref[1, :, :]); rng=MersenneTwister(10))
        corrected_frame_2 = capture!(Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
                sensor=HgCdTeAvalancheArraySensor(sampling_mode=SingleRead()),
                correction_model=correction_model),
            @view(corrected_stack_ref[2, :, :]); rng=MersenneTwister(10))
        @test corrected_stack[1, :, :] ≈ corrected_frame_1 atol=1e-12 rtol=1e-12
        @test corrected_stack[2, :, :] ≈ corrected_frame_2 atol=1e-12 rtol=1e-12
    end

    det_cmos_batched = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=CMOSSensor(column_readout_sigma=1.0))
    @test_throws InvalidConfiguration AdaptiveOpticsSim.capture_stack!(det_cmos_batched, cube_mtf, scratch_mtf; rng=MersenneTwister(10))

    det_generalized = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, psf_sampling=2, binning=2,
        bits=8, full_well=10.0, output_type=UInt8)
    input_generalized = zeros(Float64, 2, 8, 8)
    input_generalized[1, 4, 4] = 10.0
    input_generalized[2, 5, 5] = 10.0
    output_generalized = Array{UInt8}(undef, 2, 2, 2)
    generalized_stack = AdaptiveOpticsSim.capture_stack!(det_generalized, output_generalized, input_generalized; rng=MersenneTwister(10))
    @test size(generalized_stack) == (2, 2, 2)
    @test readout_ready(det_generalized)
    @test iszero(det_generalized.state.integrated_time)
    @test generalized_stack[1, :, :] == capture!(det_generalized, @view(input_generalized[1, :, :]); rng=MersenneTwister(10))
    @test generalized_stack[2, :, :] == capture!(det_generalized, @view(input_generalized[2, :, :]); rng=MersenneTwister(10))

    apd_dead_time = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), dead_time_model=NonParalyzableDeadTime(0.5))
    dead_time_out = capture!(apd_dead_time, fill(4.0, 2, 8); rng=MersenneTwister(9))
    @test dead_time_out ≈ fill(4.0 / 3.0, 2, 8)
    dead_time_meta = detector_export_metadata(apd_dead_time)
    @test dead_time_meta.dead_time_model == :nonparalyzable
    @test dead_time_meta.dead_time == 0.5
    @test supports_dead_time(apd_dead_time)
    apd_paralyzable = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), dead_time_model=ParalyzableDeadTime(0.5))
    @test capture!(apd_paralyzable, fill(4.0, 2, 8); rng=MersenneTwister(9)) ≈ fill(4.0 * exp(-2.0), 2, 8)
    @test supports_paralyzable_dead_time(apd_paralyzable)
    apd_gated = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), gate_model=DutyCycleGate(0.5))
    @test capture!(apd_gated, fill(4.0, 2, 8); rng=MersenneTwister(9)) == fill(2.0, 2, 8)
    gated_meta = detector_export_metadata(apd_gated)
    @test gated_meta.gate_model == :duty_cycle
    @test gated_meta.duty_cycle == 0.5
    @test supports_counting_gating(apd_gated)
    apd_afterpulse = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), correlation_model=AfterpulsingModel(0.25))
    @test capture!(apd_afterpulse, fill(4.0, 2, 8); rng=MersenneTwister(9)) == fill(5.0, 2, 8)
    @test supports_afterpulsing(apd_afterpulse)
    apd_crosstalk = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), correlation_model=ChannelCrosstalkModel(0.4))
    crosstalk_in = zeros(3, 3)
    crosstalk_in[2, 2] = 10.0
    crosstalk_out = capture!(apd_crosstalk, crosstalk_in; rng=MersenneTwister(9))
    @test crosstalk_out[2, 2] ≈ 6.0
    @test crosstalk_out[1, 2] ≈ 1.0
    @test crosstalk_out[2, 1] ≈ 1.0
    @test crosstalk_out[2, 3] ≈ 1.0
    @test crosstalk_out[3, 2] ≈ 1.0
    @test supports_channel_crosstalk(apd_crosstalk)
    apd_composite_corr = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), correlation_model=CompositeCountingCorrelation(AfterpulsingModel(0.1), ChannelCrosstalkModel(0.2)))
    composite_meta_apd = detector_export_metadata(apd_composite_corr)
    @test composite_meta_apd.correlation_model == :composite
    @test composite_meta_apd.afterpulse_probability == 0.1
    @test composite_meta_apd.crosstalk == 0.2

    @test_throws InvalidConfiguration APDDetector(noise=NoiseReadout(1.0))
    @test_throws InvalidConfiguration APDDetector(dead_time_model=NonParalyzableDeadTime(-1.0))
    @test_throws InvalidConfiguration APDDetector(dead_time_model=ParalyzableDeadTime(-1.0))
    @test_throws InvalidConfiguration APDDetector(gate_model=DutyCycleGate(0.0))
    @test_throws InvalidConfiguration APDDetector(correlation_model=AfterpulsingModel(1.5))
    @test_throws InvalidConfiguration APDDetector(correlation_model=ChannelCrosstalkModel(-0.1))
    apd_thermal = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=10.0,
        noise=NoiseNone(), thermal_model=FixedTemperature(temperature_K=80.0, dark_count_law=arrhenius))
    apd_thermal_meta = detector_export_metadata(apd_thermal)
    @test supports_detector_thermal_model(apd_thermal)
    @test supports_temperature_dependent_dark_counts(apd_thermal)
    @test detector_temperature(apd_thermal) == 80.0
    @test apd_thermal_meta.thermal_model == :fixed_temperature
    @test apd_thermal_meta.dark_count_law == :arrhenius
    @test effective_dark_count_rate(apd_thermal) < apd_thermal.params.dark_count_rate
    apd_dynamic = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=10.0,
        noise=NoiseNone(), thermal_model=FirstOrderThermalModel(
            ambient_temperature_K=295.0,
            setpoint_temperature_K=120.0,
            initial_temperature_K=300.0,
            time_constant_s=2.0,
            min_temperature_K=80.0,
            max_temperature_K=320.0,
            dark_count_law=arrhenius))
    apd_dynamic_initial = effective_dark_count_rate(apd_dynamic)
    @test supports_detector_thermal_model(apd_dynamic)
    @test supports_dynamic_thermal_state(apd_dynamic.params.thermal_model)
    @test thermal_state(apd_dynamic) isa DetectorThermalState
    @test advance_thermal!(apd_dynamic, 2.0) === apd_dynamic
    @test detector_temperature(apd_dynamic) ≈ 120.0 + 180.0 * exp(-1.0)
    @test effective_dark_count_rate(apd_dynamic) < apd_dynamic_initial
    capture!(apd_dynamic, fill(1.0, 2, 2); rng=MersenneTwister(25))
    @test detector_temperature(apd_dynamic) < 120.0 + 180.0 * exp(-1.0)

    det_buffered = Detector(integration_time=2.0, noise=NoiseNone(), qe=1.0, binning=1)
    frame_partial = copy(capture!(det_buffered, fill(1.0, 4, 4); rng=MersenneTwister(2), sample_time=1.0))
    @test !readout_ready(det_buffered)
    @test sum(frame_partial) == 16.0
    frame_buffered = copy(capture!(det_buffered, fill(1.0, 4, 4); rng=MersenneTwister(2), sample_time=1.0))
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
    AdaptiveOpticsSim.capture_stack!(det_stack, cube, scratch; rng=MersenneTwister(10))
    @test cube[1, :, :] ≈ fill(0.5, 4, 4)
    @test cube[2, :, :] ≈ fill(1.0, 4, 4)

    allocation_stack_detector = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse())
    allocation_stack_cube = ones(2, 4, 4)
    allocation_stack_scratch = similar(allocation_stack_cube)
    @test fixed_stack_capture_allocations(allocation_stack_detector,
        allocation_stack_cube, allocation_stack_scratch, Xoshiro(10)) == 0

    psf = reshape(Float64.(1:256), 16, 16)
    det_fused = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, psf_sampling=2, binning=2)
    frame_fused = copy(AdaptiveOpticsSim.fill_frame!(det_fused, psf, 1.0))
    manual_mid = zeros(Float64, 8, 8)
    manual_out = zeros(Float64, 4, 4)
    AdaptiveOpticsSim.bin2d!(manual_mid, psf, 2)
    AdaptiveOpticsSim.bin2d!(manual_out, manual_mid, 2)
    @test frame_fused == manual_out
end

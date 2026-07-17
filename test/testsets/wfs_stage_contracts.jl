include(joinpath(@__DIR__, "..", "wfs_stage_contract_fixtures.jl"))

function contract_rate_map(values::AbstractMatrix{T};
    sampling::NTuple{2,T}=(one(T), one(T)),
    origin::NTuple{2,T}=AdaptiveOpticsSim.centered_grid_origin(size(values),
        sampling),
    coordinate_domain::AbstractPlaneCoordinateDomain=AngularCoordinates(),
    spectral::AbstractSpectralCoordinate=MonochromaticChannel(T(0.75e-6)),
    normalization::AbstractOpticalNormalization=PhotonRateNormalization(),
    spatial_measure::AbstractSpatialMeasure=CellIntegratedMeasure(),
    coherence::AbstractCombinationPolicy=IncoherentIntensityAddition(),
    kind::AbstractOpticalPlaneKind=DetectorPlane(),
    orientation::PlaneAxisOrientation=PlaneAxisOrientation(),
) where {T<:AbstractFloat}
    metadata = OpticalPlaneMetadata(kind, values;
        coordinate_domain, sampling, origin, orientation, spectral,
        normalization, spatial_measure, coherence)
    return IntensityMap(metadata, values)
end

function contract_stage_allocation_bytes(optical_output, optical_input,
    optical_plan, observation, acquisition_plan, rng, measurement,
    estimator_plan)
    run_contract_stages!(optical_output, optical_input, optical_plan,
        observation, acquisition_plan, rng, measurement, estimator_plan)
    return @allocated run_contract_stages!(optical_output, optical_input,
        optical_plan, observation, acquisition_plan, rng, measurement,
        estimator_plan)
end

function contract_direct_allocation_bytes(measurement, input, plan)
    estimate_wfs_measurement!(measurement, input, plan)
    return @allocated estimate_wfs_measurement!(measurement, input, plan)
end

function contract_sampled_response(input::AbstractMatrix{T},
    kernel::AbstractMatrix{T}) where {T<:AbstractFloat}
    output = zeros(T, size(input))
    radius_i = fld(size(kernel, 1), 2)
    radius_j = fld(size(kernel, 2), 2)
    @inbounds for ki in axes(kernel, 1), kj in axes(kernel, 2)
        offset_i = ki - radius_i - 1
        offset_j = kj - radius_j - 1
        weight = kernel[ki, kj]
        for j in axes(input, 2), i in axes(input, 1)
            ii = i + offset_i
            jj = j + offset_j
            if checkbounds(Bool, input, ii, jj)
                output[i, j] += weight * input[ii, jj]
            end
        end
    end
    return output
end

@testset "Prepared WFS stage products and protocols" begin
    @test Base.isexported(AdaptiveOpticsSim, :WFSObservation)
    @test Base.isexported(AdaptiveOpticsSim, :WFSMeasurement)
    @test Base.isexported(AdaptiveOpticsSim, :prepare_wfs_optical_formation)
    @test Base.isexported(AdaptiveOpticsSim, :acquire_wfs_observation!)
    @test Base.isexported(AdaptiveOpticsSim, :DirectMeasurementPath)

    scalar_observation = WFSObservation(Ref(1.0); units=:electron_count,
        layout=:scalar_channel)
    vector_observation = WFSObservation(zeros(Float32, 3);
        units=:adu, layout=:channel_vector)
    matrix_observation = WFSObservation(zeros(Float32, 2, 3);
        units=:electron_count)
    stack_observation = WFSObservation(zeros(Float32, 2, 3, 4);
        units=:electron_count, layout=:read_stack)
    measurement = WFSMeasurement(zeros(Float32, 5);
        units=:radian, kind=:centroid_slopes)

    @test observation_storage(scalar_observation) isa Base.RefValue{Float64}
    @test observation_metadata(scalar_observation).dimensions == ()
    @test observation_metadata(vector_observation).dimensions == (3,)
    @test observation_metadata(matrix_observation).dimensions == (2, 3)
    @test observation_metadata(stack_observation).dimensions == (2, 3, 4)
    @test observation_units(stack_observation) === :electron_count
    @test measurement_storage(measurement) isa Vector{Float32}
    @test measurement_metadata(measurement).kind === :centroid_slopes
    @test measurement_units(measurement) === :radian
    @test isconcretetype(typeof(stack_observation))
    @test isconcretetype(typeof(measurement))

    missing_units = try
        WFSObservation(zeros(2); units=nothing)
        nothing
    catch err
        err
    end
    @test missing_units isa WFSPreparationError
    @test missing_units.stage === :acquisition
    @test missing_units.reason === :units

    metadata_storage = zeros(2)
    declared_observation_metadata = WFSObservationMetadata(metadata_storage;
        layout=:packed)
    observation_with_metadata = WFSObservation(metadata_storage, :adu,
        declared_observation_metadata)
    @test observation_metadata(observation_with_metadata) ===
        declared_observation_metadata

    measurement_metadata_storage = zeros(2)
    declared_measurement_metadata = WFSMeasurementMetadata(
        measurement_metadata_storage; kind=:phase)
    measurement_with_metadata = WFSMeasurement(measurement_metadata_storage,
        :radian, declared_measurement_metadata)
    @test measurement_metadata(measurement_with_metadata) ===
        declared_measurement_metadata

    wrong_metadata = WFSObservationMetadata(zeros(2); layout=:dense)
    mismatched = try
        WFSObservation(zeros(3), :adu, wrong_metadata)
        nothing
    catch err
        err
    end
    @test mismatched isa WFSPreparationError
    @test mismatched.reason === :shape

    wrong_measurement_metadata = WFSMeasurementMetadata(zeros(2);
        kind=:phase)
    mismatched_measurement = try
        WFSMeasurement(zeros(3), :radian, wrong_measurement_metadata)
        nothing
    catch err
        err
    end
    @test mismatched_measurement isa WFSPreparationError
    @test mismatched_measurement.stage === :estimation
    @test mismatched_measurement.reason === :shape

    first_device_storage = ContractDeviceArray(zeros(2),
        ContractPlaneDevice(1))
    second_device_storage = ContractDeviceArray(zeros(2),
        ContractPlaneDevice(2))
    first_device_metadata = WFSObservationMetadata(first_device_storage;
        layout=:device_channels)
    @test WFSObservation(first_device_storage, :adu,
        first_device_metadata).metadata.device == ContractPlaneDevice(1)
    @test typeof(backend(first_device_storage)) ===
        typeof(backend(second_device_storage))
    device_mismatch = try
        WFSObservation(second_device_storage, :adu, first_device_metadata)
        nothing
    catch err
        err
    end
    @test device_mismatch isa WFSPreparationError
    @test device_mismatch.stage === :acquisition
    @test device_mismatch.reason === :device

    for invalid_storage in (Any[1.0], Real[1.0])
        invalid_numeric_type = try
            WFSObservation(invalid_storage; units=:adu)
            nothing
        catch err
            err
        end
        @test invalid_numeric_type isa WFSPreparationError
        @test invalid_numeric_type.reason === :numeric_type
    end
    invalid_measurement_numeric_type = try
        WFSMeasurement(Real[1.0]; units=:radian, kind=:phase)
        nothing
    catch err
        err
    end
    @test invalid_measurement_numeric_type isa WFSPreparationError
    @test invalid_measurement_numeric_type.reason === :numeric_type
end

@testset "Explicit WFS optical formation" begin
    T = Float64
    tel = Telescope(resolution=4, diameter=T(2), central_obstruction=zero(T),
        T=T)
    src = Source(band=:custom, wavelength=T(0.75e-6),
        photon_irradiance=T(8), T=T)
    pupil = PupilFunction(tel; T=T)
    pupil.opd .= reshape(T.(1:16), 4, 4) .* T(1e-9)
    pupil_opd = copy(pupil.opd)
    pupil_amplitude = copy(pupil.amplitude)

    rate = contract_rate_map(zeros(T, 4, 4))
    model = ContractRateModel(T(3), T(1e6), rate.metadata)
    plan = prepare_wfs_optical_formation(model, pupil, rate)
    @test @inferred(form_wfs_optical_products!(rate, pupil, plan)) === rate
    expected = @. T(3) * (abs2(pupil_amplitude) + T(1e6) * abs(pupil_opd))
    @test rate.values == expected
    @test pupil.opd == pupil_opd
    @test pupil.amplitude == pupil_amplitude

    original_rate = copy(rate.values)
    tel.state.opd .= T(99)
    fill!(rate.values, zero(T))
    form_wfs_optical_products!(rate, pupil, plan)
    @test rate.values == original_rate

    field = ElectricField(pupil, src; zero_padding=1, T=T)
    field_plan = prepare_pupil_field(tel, pupil, src, field)
    fill_electric_field!(field, pupil, field_plan)
    field_values = copy(field.values)
    field_rate = contract_rate_map(zeros(T, 4, 4))
    field_model = ContractRateModel(T(2), zero(T), field_rate.metadata)
    prepared_field = prepare_wfs_optical_formation(field_model, field,
        field_rate)
    @test @inferred(form_wfs_optical_products!(field_rate, field,
        prepared_field)) === field_rate
    @test field_rate.values ≈ T(2) .* abs2.(field_values) atol=0 rtol=0
    @test field.values == field_values

    bad_maps = (
        contract_rate_map(zeros(T, 4, 4);
            kind=FocalPlane()),
        contract_rate_map(zeros(T, 4, 4);
            normalization=DimensionlessNormalization()),
        contract_rate_map(zeros(T, 4, 4);
            spatial_measure=PointSampledMeasure()),
        contract_rate_map(zeros(T, 4, 4);
            coherence=NonCombinableProduct()),
        contract_rate_map(zeros(T, 4, 4);
            spectral=UnspecifiedSpectralCoordinate()),
    )
    for bad_map in bad_maps
        err = try
            prepare_wfs_optical_formation(
                ContractRateModel(one(T), zero(T), bad_map.metadata),
                pupil, bad_map)
            nothing
        catch caught
            caught
        end
        @test err isa WFSPreparationError
        @test err.stage === :optical_formation
    end

    incompatible_products = (
        contract_rate_map(zeros(T, 4, 4);
            sampling=(T(2), T(2))),
        contract_rate_map(zeros(T, 4, 4);
            coordinate_domain=MetricCoordinates()),
        contract_rate_map(zeros(T, 4, 4);
            origin=(T(-1), T(-1))),
        contract_rate_map(zeros(T, 4, 4);
            orientation=PlaneAxisOrientation((:y, :x), (1, -1))),
        contract_rate_map(zeros(T, 4, 4);
            spectral=MonochromaticChannel(T(0.9e-6))),
    )
    for incompatible in incompatible_products
        err = try
            prepare_wfs_optical_formation(model, pupil, incompatible)
            nothing
        catch caught
            caught
        end
        @test err isa WFSPreparationError
        @test err.reason === :plane_metadata
    end

    replacement = contract_rate_map(zeros(T, 4, 4);
        sampling=rate.metadata.sampling,
        origin=rate.metadata.origin)
    replacement_before = copy(replacement.values)
    @test_throws WFSPreparationError form_wfs_optical_products!(replacement,
        pupil, plan)
    @test replacement.values == replacement_before
end

@testset "Single-plane WFS composition and exact-once exposure" begin
    T = Float64
    tel = Telescope(resolution=6, diameter=T(3), central_obstruction=zero(T),
        T=T)
    pupil = PupilFunction(tel; T=T)
    pupil.opd .= T(2e-9)
    rate = contract_rate_map(zeros(T, 6, 6))
    optical_model = ContractRateModel(T(5), T(1e6), rate.metadata)
    optical_plan = prepare_wfs_optical_formation(optical_model, pupil, rate)

    detector = Detector(noise=NoiseNone(), integration_time=T(0.4),
        qe=T(0.5), response_model=NullFrameResponse(), T=T)
    binding = contract_detector_binding(detector, rate)
    observation = binding.observation
    acquisition_plan = prepare_wfs_acquisition(
        ContractDetectorAcquisitionModel((binding,)), rate, observation)
    measurement = WFSMeasurement(Ref(zero(T)); units=:electron_count,
        kind=:total_signal)
    sum_estimator = ContractSumEstimator(:electron_count, :total_signal)
    estimator_plan = prepare_wfs_estimation(sum_estimator,
        observation, measurement)
    rng = Xoshiro(0x51a7)

    @test wfs_measurement_path(estimator_plan) isa AcquiredObservationPath
    @test @inferred(form_wfs_optical_products!(rate, pupil,
        optical_plan)) === rate
    @test @inferred(acquire_wfs_observation!(observation, rate,
        acquisition_plan, rng)) === observation
    @test @inferred(estimate_wfs_measurement!(measurement, observation,
        estimator_plan)) === measurement
    expected = sum(rate.values) * T(0.4) * T(0.5)
    @test measurement.storage[] ≈ expected atol=eps(T) * expected
    @test observation.storage === output_frame(detector)
    @test contract_stage_allocation_bytes(rate, pupil, optical_plan,
        observation, acquisition_plan, rng, measurement,
        estimator_plan) == 0

    incompatible_measurement = WFSMeasurement(Ref(zero(T));
        units=:electron_count, kind=:phase)
    incompatible_estimator = try
        prepare_wfs_estimation(sum_estimator, observation,
            incompatible_measurement)
        nothing
    catch err
        err
    end
    @test incompatible_estimator isa WFSPreparationError
    @test incompatible_estimator.stage === :estimation
    @test incompatible_estimator.reason === :estimator

    incompatible_numeric_measurement = WFSMeasurement(zeros(Int, 1);
        units=:electron_count, kind=:total_signal)
    numeric_estimator_error = try
        prepare_wfs_estimation(sum_estimator, observation,
            incompatible_numeric_measurement)
        nothing
    catch err
        err
    end
    @test numeric_estimator_error isa WFSPreparationError
    @test numeric_estimator_error.reason === :estimator
    @test incompatible_numeric_measurement.storage == zeros(Int, 1)

    response_rate_values = zeros(T, 6, 6)
    response_rate_values[3, 2] = T(4)
    response_rate = contract_rate_map(response_rate_values)
    kernel = T[0 0 0; 0.1 0.2 0.7; 0 0 0]
    response_detector = Detector(noise=NoiseNone(), integration_time=T(1.5),
        qe=T(0.4), binning=2,
        response_model=SampledFrameResponse(kernel; normalize=false, T=T),
        T=T)
    response_binding = contract_detector_binding(response_detector,
        response_rate)
    response_plan = prepare_wfs_acquisition(
        ContractDetectorAcquisitionModel((response_binding,)), response_rate,
        response_binding.observation)
    acquire_wfs_observation!(response_binding.observation, response_rate,
        response_plan, Xoshiro(0x45))
    presampled = contract_sampled_response(response_rate_values, kernel)
    expected_response = zeros(T, 3, 3)
    AdaptiveOpticsSim.bin2d!(expected_response, presampled, 2)
    expected_response .*= T(1.5) * T(0.4)
    @test response_binding.observation.storage ≈ expected_response atol=eps(T)
    binned_first = zeros(T, 3, 3)
    AdaptiveOpticsSim.bin2d!(binned_first, response_rate_values, 2)
    reordered = contract_sampled_response(binned_first, kernel)
    reordered .*= T(1.5) * T(0.4)
    @test !isapprox(response_binding.observation.storage, reordered;
        atol=T(1e-12), rtol=T(1e-12))
    response_twice = contract_sampled_response(presampled, kernel)
    response_twice_binned = zeros(T, 3, 3)
    AdaptiveOpticsSim.bin2d!(response_twice_binned, response_twice, 2)
    response_twice_binned .*= T(1.5) * T(0.4)
    @test !isapprox(response_binding.observation.storage,
        response_twice_binned; atol=T(1e-12), rtol=T(1e-12))
end

@testset "Prepared physical Shack-Hartmann stages" begin
    T = Float64
    tel = Telescope(resolution=16, diameter=T(8),
        central_obstruction=zero(T), T=T)
    src = Source(band=:custom, wavelength=T(0.75e-6),
        photon_irradiance=T(10), T=T)
    pupil = PupilFunction(tel; T=T)
    pupil.opd .= reshape(T.(1:256), 16, 16) .* T(1e-10)
    pupil_before = copy(pupil.opd)

    mla = MicrolensArray(; n_lenslets=4, diffraction_padding=2,
        n_pix_subap=4, T=T)
    @test mla.params isa MicrolensArrayParams{T}
    @test mla.params.n_lenslets == 4
    @test_throws InvalidConfiguration MicrolensArray(; n_lenslets=4,
        n_pix_subap=3, T=T)
    configured_extraction = CenterOfGravityExtraction(T(0.125); T=T)
    configured_sensor = ShackHartmannWFS(tel; n_lenslets=4,
        mode=Diffractive(), n_pix_subap=4,
        slope_extraction=configured_extraction, T=T)
    @test slope_extraction_model(configured_sensor).threshold == T(0.125)

    geometric = ShackHartmannWFS(tel; n_lenslets=4, mode=Geometric(),
        T=T)
    @test microlens_array(geometric).params.n_lenslets == 4
    @test geometric.optical_workspace === nothing
    @test geometric.acquisition === nothing
    measure!(geometric, tel)
    @test @allocated(measure!(geometric, tel)) == 0
    direct_measurement = WFSMeasurement(similar(slopes(geometric));
        units=:radian, kind=:geometric_slopes)
    direct_plan = prepare_wfs_estimation(geometric, pupil,
        direct_measurement)
    @test wfs_measurement_path(direct_plan) isa DirectMeasurementPath
    estimate_wfs_measurement!(direct_measurement, pupil, direct_plan)
    @test all(isfinite, direct_measurement.storage)
    @test any(!iszero, direct_measurement.storage)

    copyto!(tel.state.opd, pupil.opd)
    legacy = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(),
        n_pix_subap=4, T=T)
    staged = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(),
        n_pix_subap=4, T=T)
    prepare_sampling!(legacy, tel, src)
    sampled_spots_peak!(legacy, tel, src)
    expected_rate = shack_hartmann_detector_image(
        AdaptiveOpticsSim.sh_sampled_spot_cube(legacy), 4)

    rate = shack_hartmann_rate_map(staged, pupil, src)
    front_end = ShackHartmannOpticalFrontEnd(staged, src)
    optical_plan = prepare_wfs_optical_formation(front_end, pupil, rate)
    form_wfs_optical_products!(rate, pupil, optical_plan)
    @test rate.values == expected_rate
    @test pupil.opd == pupil_before
    @test rate.metadata.normalization isa PhotonRateNormalization
    @test rate.metadata.spatial_measure isa CellIntegratedMeasure
    @test rate.metadata.coherence isa IncoherentIntensityAddition
    @test @allocated(form_wfs_optical_products!(rate, pupil,
        optical_plan)) == 0

    field = ElectricField(pupil, src; zero_padding=1, T=T)
    field_formation = prepare_pupil_field(tel, pupil, src, field;
        center_even_grid=false)
    fill_electric_field!(field, pupil, field_formation)
    field_sensor = ShackHartmannWFS(tel; n_lenslets=4,
        mode=Diffractive(), n_pix_subap=4, T=T)
    field_rate = shack_hartmann_rate_map(field_sensor, field)
    field_plan = prepare_wfs_optical_formation(
        ShackHartmannOpticalFrontEnd(field_sensor), field, field_rate)
    form_wfs_optical_products!(field_rate, field, field_plan)
    @test field_rate.values ≈ rate.values atol=T(2e-12) rtol=T(2e-12)
    @test @allocated(form_wfs_optical_products!(field_rate, field,
        field_plan)) == 0

    rate_before = copy(rate.values)
    tel.state.opd .= T(99)
    form_wfs_optical_products!(rate, pupil, optical_plan)
    @test rate.values == rate_before
    @test pupil.opd == pupil_before

    short_detector = Detector(noise=NoiseNone(), integration_time=T(0.25),
        qe=one(T), response_model=NullFrameResponse(), T=T)
    long_detector = Detector(noise=NoiseNone(), integration_time=T(0.75),
        qe=one(T), response_model=NullFrameResponse(), T=T)
    short_observation = WFSObservation(similar(rate.values);
        units=:electron_count, layout=:lenslet_mosaic)
    long_observation = WFSObservation(similar(rate.values);
        units=:electron_count, layout=:lenslet_mosaic)
    short_plan = prepare_wfs_acquisition(short_detector, rate,
        short_observation)
    long_plan = prepare_wfs_acquisition(long_detector, rate,
        long_observation)
    short_rng = Xoshiro(0x51)
    long_rng = Xoshiro(0x52)
    acquire_wfs_observation!(short_observation, rate, short_plan, short_rng)
    acquire_wfs_observation!(long_observation, rate, long_plan, long_rng)
    @test short_observation.storage ≈ rate.values .* T(0.25) atol=0 rtol=0
    @test long_observation.storage ≈ rate.values .* T(0.75) atol=0 rtol=0
    @test long_observation.storage ≈ T(3) .* short_observation.storage
    @test rate.values == rate_before
    @test @allocated(acquire_wfs_observation!(short_observation, rate,
        short_plan, short_rng)) == 0

    fill!(staged.calibration.reference_signal_2d, zero(T))
    staged.calibration.slopes_units = one(T)
    measurement = WFSMeasurement(similar(slopes(staged));
        units=:pixel, kind=:centroid_slopes)
    estimator_plan = prepare_wfs_estimation(staged, short_observation,
        measurement)
    @test wfs_measurement_path(estimator_plan) isa AcquiredObservationPath
    estimate_wfs_measurement!(measurement, short_observation,
        estimator_plan)
    @test all(isfinite, measurement.storage)
    @test @allocated(estimate_wfs_measurement!(measurement,
        short_observation, estimator_plan)) == 0

    response_values = zeros(T, size(rate.values))
    response_values[7, 6] = T(4)
    response_rate = contract_rate_map(response_values;
        sampling=rate.metadata.sampling,
        spectral=rate.metadata.spectral)
    kernel = T[0 0 0; 0.1 0.2 0.7; 0 0 0]
    response_detector = Detector(noise=NoiseNone(), integration_time=T(1.5),
        qe=T(0.4), binning=2,
        response_model=SampledFrameResponse(kernel; normalize=false, T=T),
        T=T)
    response_observation = WFSObservation(zeros(T, 8, 8);
        units=:electron_count, layout=:lenslet_mosaic)
    response_plan = prepare_wfs_acquisition(response_detector,
        response_rate, response_observation)
    acquire_wfs_observation!(response_observation, response_rate,
        response_plan, Xoshiro(0x53))
    presampled = contract_sampled_response(response_values, kernel)
    expected_response = zeros(T, 8, 8)
    AdaptiveOpticsSim.bin2d!(expected_response, presampled, 2)
    expected_response .*= T(1.5) * T(0.4)
    @test response_observation.storage ≈ expected_response atol=eps(T)
    binned_first = zeros(T, 8, 8)
    AdaptiveOpticsSim.bin2d!(binned_first, response_values, 2)
    reordered = contract_sampled_response(binned_first, kernel)
    reordered .*= T(1.5) * T(0.4)
    @test !isapprox(response_observation.storage, reordered;
        atol=T(1e-12), rtol=T(1e-12))

    spectral = with_spectrum(src, SpectralBundle(
        T[0.9 * wavelength(src), 1.1 * wavelength(src)], T[0.4, 0.6];
        T=T))
    spectral_sensor = ShackHartmannWFS(tel; n_lenslets=4,
        mode=Diffractive(), n_pix_subap=4, T=T)
    spectral_rates = shack_hartmann_rate_map(spectral_sensor, pupil,
        spectral)
    @test spectral_rates isa OpticalProductBundle
    @test length(spectral_rates) == 2
    @test spectral_rates[1].metadata.sampling !=
        spectral_rates[2].metadata.sampling
    spectral_plan = prepare_wfs_optical_formation(
        ShackHartmannOpticalFrontEnd(spectral_sensor, spectral), pupil,
        spectral_rates)
    form_wfs_optical_products!(spectral_rates, pupil, spectral_plan)
    @test sum(spectral_rates[1].values) /
        sum(spectral_rates[2].values) ≈ T(2 / 3) rtol=T(2e-6)
    @test_throws WFSPreparationError prepare_wfs_optical_formation(
        ShackHartmannOpticalFrontEnd(spectral_sensor, spectral), pupil,
        spectral_rates[1])
    incompatible_spectral_sensor = ShackHartmannWFS(tel; n_lenslets=4,
        mode=Diffractive(), n_pix_subap=4, pixel_scale=T(0.04), T=T)
    incompatible_spectral_rates = shack_hartmann_rate_map(
        incompatible_spectral_sensor, pupil, spectral)
    @test_throws WFSPreparationError prepare_wfs_optical_formation(
        ShackHartmannOpticalFrontEnd(incompatible_spectral_sensor, spectral),
        pupil, incompatible_spectral_rates)

    lgs = LGSSource(wavelength=wavelength(src),
        photon_irradiance=T(6), elongation_factor=T(1.8), T=T)
    copyto!(tel.state.opd, pupil.opd)
    legacy_lgs = ShackHartmannWFS(tel; n_lenslets=4,
        mode=Diffractive(), n_pix_subap=4, T=T)
    prepare_sampling!(legacy_lgs, tel, lgs)
    AdaptiveOpticsSim.sampled_spots_peak!(legacy_lgs, tel, lgs)
    expected_lgs = shack_hartmann_detector_image(
        legacy_lgs.acquisition.spot_cube, 4)
    staged_lgs = ShackHartmannWFS(tel; n_lenslets=4,
        mode=Diffractive(), n_pix_subap=4, T=T)
    lgs_rate = shack_hartmann_rate_map(staged_lgs, pupil, lgs)
    lgs_plan = prepare_wfs_optical_formation(
        ShackHartmannOpticalFrontEnd(staged_lgs, lgs), pupil, lgs_rate)
    form_wfs_optical_products!(lgs_rate, pupil, lgs_plan)
    @test lgs_rate.values ≈ expected_lgs rtol=T(2e-12) atol=T(2e-12)

    sodium_lgs = LGSSource(wavelength=wavelength(src),
        photon_irradiance=T(6),
        na_profile=T[80_000 90_000 100_000; 0.2 0.6 0.2],
        laser_coordinates=(T(1), T(-0.5)), fwhm_spot_up=T(0.8), T=T)
    legacy_sodium = ShackHartmannWFS(tel; n_lenslets=4,
        mode=Diffractive(), n_pix_subap=4, T=T)
    prepare_sampling!(legacy_sodium, tel, sodium_lgs)
    AdaptiveOpticsSim.sampled_spots_peak!(legacy_sodium, tel, sodium_lgs)
    expected_sodium = shack_hartmann_detector_image(
        legacy_sodium.acquisition.spot_cube, 4)
    staged_sodium = ShackHartmannWFS(tel; n_lenslets=4,
        mode=Diffractive(), n_pix_subap=4, T=T)
    sodium_rate = shack_hartmann_rate_map(staged_sodium, pupil, sodium_lgs)
    sodium_plan = prepare_wfs_optical_formation(
        ShackHartmannOpticalFrontEnd(staged_sodium, sodium_lgs), pupil,
        sodium_rate)
    form_wfs_optical_products!(sodium_rate, pupil, sodium_plan)
    @test sodium_rate.values ≈ expected_sodium rtol=T(2e-12) atol=T(2e-12)

    asterism = Asterism([
        Source(band=:custom, wavelength=wavelength(src),
            photon_irradiance=T(3), T=T),
        Source(band=:custom, wavelength=wavelength(src),
            photon_irradiance=T(7), T=T),
    ])
    legacy_asterism = ShackHartmannWFS(tel; n_lenslets=4,
        mode=Diffractive(), n_pix_subap=4, T=T)
    prepare_sampling!(legacy_asterism, tel, first(asterism.sources))
    AdaptiveOpticsSim.sampled_spots_peak_asterism_stacked!(
        AdaptiveOpticsSim.ScalarCPUStyle(), legacy_asterism, tel, asterism)
    expected_asterism = shack_hartmann_detector_image(
        legacy_asterism.acquisition.spot_cube, 4)
    staged_asterism = ShackHartmannWFS(tel; n_lenslets=4,
        mode=Diffractive(), n_pix_subap=4, T=T)
    asterism_rate = shack_hartmann_rate_map(staged_asterism, pupil,
        asterism)
    asterism_plan = prepare_wfs_optical_formation(
        ShackHartmannOpticalFrontEnd(staged_asterism, asterism), pupil,
        asterism_rate)
    form_wfs_optical_products!(asterism_rate, pupil, asterism_plan)
    @test asterism_rate.values ≈ expected_asterism rtol=T(2e-12) atol=T(2e-12)
end

@testset "Rate fan-out, bundles, and packed observations" begin
    T = Float64
    tel = Telescope(resolution=4, diameter=T(2), central_obstruction=zero(T),
        T=T)
    pupil = PupilFunction(tel; T=T)
    pupil.opd .= reshape(T.(1:16), 4, 4) .* T(1e-9)

    shared_rate = contract_rate_map(zeros(T, 4, 4))
    shared_model = ContractRateModel(T(4), T(1e6), shared_rate.metadata)
    shared_plan = prepare_wfs_optical_formation(shared_model, pupil,
        shared_rate)
    form_wfs_optical_products!(shared_rate, pupil, shared_plan)
    shared_before = copy(shared_rate.values)
    short_detector = Detector(noise=NoiseNone(), integration_time=T(0.25),
        qe=T(0.8), T=T)
    long_detector = Detector(noise=NoiseNone(), integration_time=T(0.75),
        qe=T(0.8), T=T)
    short_binding = contract_detector_binding(short_detector, shared_rate)
    long_binding = contract_detector_binding(long_detector, shared_rate)
    observations = (short_binding.observation, long_binding.observation)
    fanout_plan = prepare_wfs_acquisition(
        ContractDetectorAcquisitionModel((short_binding, long_binding)),
        shared_rate, observations)
    fanout_measurement = WFSMeasurement(zeros(T, 2);
        units=:electron_count, kind=:channel_totals)
    fanout_estimator = prepare_wfs_estimation(
        ContractSumEstimator(:electron_count, :channel_totals),
        observations, fanout_measurement)
    rngs = (Xoshiro(1), Xoshiro(2))
    short_before = copy(output_frame(short_detector))
    long_before = copy(output_frame(long_detector))
    short_rng_error = try
        acquire_wfs_observation!(observations, shared_rate, fanout_plan,
            (Xoshiro(1),))
        nothing
    catch err
        err
    end
    @test short_rng_error isa WFSPreparationError
    @test short_rng_error.reason === :detector_mapping
    @test output_frame(short_detector) == short_before
    @test output_frame(long_detector) == long_before
    wrong_rng_error = try
        acquire_wfs_observation!(observations, shared_rate, fanout_plan,
            (Xoshiro(1), 1))
        nothing
    catch err
        err
    end
    @test wrong_rng_error isa WFSPreparationError
    @test wrong_rng_error.reason === :detector_mapping
    @test output_frame(short_detector) == short_before
    @test output_frame(long_detector) == long_before
    @test @inferred(acquire_wfs_observation!(observations, shared_rate,
        fanout_plan, rngs)) === observations
    @test shared_rate.values == shared_before
    estimate_wfs_measurement!(fanout_measurement, observations,
        fanout_estimator)
    @test fanout_measurement.storage[2] ≈
        T(3) * fanout_measurement.storage[1] atol=T(1e-12)

    first_rate = contract_rate_map(zeros(T, 4, 4);
        sampling=(T(0.5), T(0.5)),
        spectral=MonochromaticChannel(T(0.6e-6)))
    second_rate = contract_rate_map(zeros(T, 4, 4);
        sampling=(T(0.25), T(0.25)),
        spectral=MonochromaticChannel(T(0.9e-6)))
    bundle = OpticalProductBundle(first_rate, second_rate)
    bundle_model = ContractBundleRateModel((
        ContractRateModel(T(2), T(1e6), first_rate.metadata),
        ContractRateModel(T(7), T(1e6), second_rate.metadata),
    ))
    bundle_plan = prepare_wfs_optical_formation(bundle_model, pupil, bundle)
    @test @inferred(form_wfs_optical_products!(bundle, pupil,
        bundle_plan)) === bundle
    @test first_rate.metadata.spectral == MonochromaticChannel(T(0.6e-6))
    @test second_rate.metadata.spectral == MonochromaticChannel(T(0.9e-6))
    @test first_rate.values != second_rate.values

    plane_count_error = try
        prepare_wfs_optical_formation(bundle_model, pupil,
            OpticalProductBundle(first_rate))
        nothing
    catch err
        err
    end
    @test plane_count_error isa WFSPreparationError
    @test plane_count_error.stage === :optical_formation
    @test plane_count_error.reason === :plane_count

    first_detector = Detector(noise=NoiseNone(), integration_time=T(0.4),
        qe=T(0.5), T=T)
    second_detector = Detector(noise=NoiseNone(), integration_time=T(0.6),
        qe=T(0.25), T=T)
    first_binding = contract_detector_binding(first_detector, first_rate)
    second_binding = contract_detector_binding(second_detector, second_rate)
    bundle_observations = (first_binding.observation,
        second_binding.observation)
    bundle_acquisition = prepare_wfs_acquisition(
        ContractDetectorAcquisitionModel((first_binding, second_binding)),
        bundle, bundle_observations)
    bundle_measurement = WFSMeasurement(zeros(T, 2);
        units=:electron_count, kind=:branch_totals)
    bundle_estimator = prepare_wfs_estimation(
        ContractSumEstimator(:electron_count, :branch_totals),
        bundle_observations, bundle_measurement)
    @test @inferred(acquire_wfs_observation!(bundle_observations, bundle,
        bundle_acquisition, (Xoshiro(3), Xoshiro(4)))) ===
        bundle_observations
    @test @inferred(estimate_wfs_measurement!(bundle_measurement,
        bundle_observations, bundle_estimator)) === bundle_measurement

    first_snapshot = copy(output_frame(first_detector))
    second_snapshot = copy(output_frame(second_detector))
    missing_branch_error = try
        prepare_wfs_acquisition(
            ContractDetectorAcquisitionModel((first_binding, second_binding)),
            OpticalProductBundle(first_rate), bundle_observations)
        nothing
    catch err
        err
    end
    @test missing_branch_error isa WFSPreparationError
    @test missing_branch_error.reason === :detector_mapping
    @test output_frame(first_detector) == first_snapshot
    @test output_frame(second_detector) == second_snapshot

    duplicate_error = try
        prepare_wfs_acquisition(
            ContractDetectorAcquisitionModel((first_binding, first_binding)),
            first_rate, (first_binding.observation,
                first_binding.observation))
        nothing
    catch err
        err
    end
    @test duplicate_error isa WFSPreparationError
    @test duplicate_error.reason === :detector_mapping

    packed_first = contract_rate_map(fill(T(2), 2, 2))
    packed_second = contract_rate_map(fill(T(5), 2, 2);
        spectral=MonochromaticChannel(T(0.8e-6)))
    packed_bundle = OpticalProductBundle(packed_first, packed_second)
    regions = ((1:2, 1:2), (3:4, 1:2))
    packed = WFSObservation(zeros(T, 4, 2); units=:electron_count,
        layout=regions)
    packed_model = ContractPackedAcquisition(regions, T(0.3))
    packed_plan = prepare_wfs_acquisition(packed_model, packed_bundle, packed)
    integer_packed = WFSObservation(zeros(Int, 4, 2);
        units=:electron_count, layout=regions)
    integer_packed_before = copy(integer_packed.storage)
    integer_packed_error = try
        prepare_wfs_acquisition(packed_model, packed_bundle, integer_packed)
        nothing
    catch err
        err
    end
    @test integer_packed_error isa WFSPreparationError
    @test integer_packed_error.reason === :detector_mapping
    @test integer_packed.storage == integer_packed_before
    density_first = contract_rate_map(fill(T(2), 2, 2);
        sampling=(T(0.25), T(0.5)),
        spatial_measure=SpatialDensityMeasure())
    density_bundle = OpticalProductBundle(density_first, packed_second)
    density_before = copy(packed.storage)
    density_error = try
        prepare_wfs_acquisition(packed_model, density_bundle, packed)
        nothing
    catch err
        err
    end
    @test density_error isa WFSPreparationError
    @test density_error.reason === :radiometry
    @test packed.storage == density_before
    packed_rng = Xoshiro(9)
    @test @inferred(acquire_wfs_observation!(packed, packed_bundle,
        packed_plan, packed_rng)) === packed
    @test packed.storage[1:2, :] == fill(T(0.6), 2, 2)
    @test packed.storage[3:4, :] == fill(T(1.5), 2, 2)
    @test @allocated(acquire_wfs_observation!(packed, packed_bundle,
        packed_plan, packed_rng)) == 0

    for (bad_regions, duration, reason) in (
        (((1:2, 1:2), (2:3, 1:2)), T(0.3), :detector_mapping),
        (((1:2, 1:2), (4:4, 1:2)), T(0.3), :detector_mapping),
        (((1:2, 1:2), (3:5, 1:2)), T(0.3), :detector_mapping),
        (regions, zero(T), :duration),
        (regions, -one(T), :duration),
        (regions, T(Inf), :duration),
        (regions, T(NaN), :duration),
    )
        before = copy(packed.storage)
        bad = try
            bad_observation = WFSObservation(packed.storage;
                units=:electron_count, layout=bad_regions)
            prepare_wfs_acquisition(
                ContractPackedAcquisition(bad_regions, duration),
                packed_bundle, bad_observation)
            nothing
        catch err
            err
        end
        @test bad isa WFSPreparationError
        @test bad.reason === reason
        @test packed.storage == before
    end
end

@testset "Intentional direct WFS measurements" begin
    T = Float64
    tel = Telescope(resolution=4, diameter=T(2), central_obstruction=zero(T),
        T=T)
    src = Source(band=:custom, wavelength=T(0.75e-6),
        photon_irradiance=one(T), T=T)
    pupil = PupilFunction(tel; T=T)
    pupil.opd .= reshape(T.(1:16), 4, 4) .* T(1e-9)
    direct_measurement = WFSMeasurement(zeros(T, 2);
        units=:metre, kind=:direct_geometric_signal)
    direct_plan = prepare_wfs_estimation(ContractDirectEstimator(T(2),
        :metre, :direct_geometric_signal),
        pupil, direct_measurement)
    @test wfs_measurement_path(direct_plan) isa DirectMeasurementPath
    @test measurement_units(direct_measurement) === :metre
    @test !applicable(form_wfs_optical_products!, direct_measurement, pupil,
        direct_plan)
    @test !applicable(acquire_wfs_observation!, direct_measurement, pupil,
        direct_plan, Xoshiro(1))
    @test @inferred(estimate_wfs_measurement!(direct_measurement, pupil,
        direct_plan)) === direct_measurement
    @test direct_measurement.storage[1] == T(2) * sum(pupil.opd)
    @test direct_measurement.storage[2] == T(2) * sum(pupil.amplitude)
    @test contract_direct_allocation_bytes(direct_measurement, pupil,
        direct_plan) == 0

    integer_measurement = WFSMeasurement(zeros(Int, 2);
        units=:metre, kind=:direct_geometric_signal)
    integer_measurement_error = try
        prepare_wfs_estimation(ContractDirectEstimator(T(2), :metre,
            :direct_geometric_signal), pupil, integer_measurement)
        nothing
    catch err
        err
    end
    @test integer_measurement_error isa WFSPreparationError
    @test integer_measurement_error.reason === :estimator
    @test all(iszero, integer_measurement.storage)

    field = ElectricField(pupil, src; zero_padding=1, T=T)
    field_formation = prepare_pupil_field(tel, pupil, src, field)
    fill_electric_field!(field, pupil, field_formation)
    field_measurement = WFSMeasurement(zeros(T, 2);
        units=:field_sum, kind=:direct_field_signal)
    field_direct = prepare_wfs_estimation(ContractDirectEstimator(one(T),
        :field_sum, :direct_field_signal),
        field, field_measurement)
    @test wfs_measurement_path(field_direct) isa DirectMeasurementPath
    @test @inferred(estimate_wfs_measurement!(field_measurement, field,
        field_direct)) === field_measurement

    integer_field_measurement = WFSMeasurement(zeros(Int, size(field.values));
        units=:field_intensity, kind=:direct_field_measurement)
    integer_field_before = copy(integer_field_measurement.storage)
    integer_field_error = try
        prepare_wfs_estimation(ContractDirectCopyEstimator(
            :field_intensity, :direct_field_measurement), field,
            integer_field_measurement)
        nothing
    catch err
        err
    end
    @test integer_field_error isa WFSPreparationError
    @test integer_field_error.reason === :estimator
    @test integer_field_measurement.storage == integer_field_before

    replacement = WFSMeasurement(zeros(T, 2);
        units=:metre, kind=:direct_geometric_signal)
    replacement_before = copy(replacement.storage)
    @test_throws WFSPreparationError estimate_wfs_measurement!(replacement,
        pupil, direct_plan)
    @test replacement.storage == replacement_before
end

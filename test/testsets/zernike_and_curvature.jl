@inline wfs_optical_rate_storage(wfs::ZernikeWFS) =
    wfs.acquisition.products.frame
@inline wfs_optical_rate_storage(wfs) =
    wfs.acquisition.products.frame

@testset "Zernike optical front end" begin
    tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.0)
    pupil = PupilFunction(tel)
    source = Source(band=:custom, wavelength=0.75e-6, photon_irradiance=10.0)
    wfs = ZernikeWFS(tel; pupil_samples=8, diffraction_padding=2)
    front_end = ZernikeOpticalFrontEnd(wfs, source)
    rate = zernike_rate_map(front_end, pupil)
    plan = prepare_wfs_optics(front_end, pupil, rate)
    form_wfs_optical_products!(rate, pupil, plan)
    @test all(isfinite, rate.values)
    @test all(>=(0.0), rate.values)
    if coverage_instrumented()
        @test_skip "Zernike optical allocation assertion is disabled under coverage instrumentation"
    else
        @test @allocated(form_wfs_optical_products!(rate, pupil, plan)) == 0
    end

    attenuated_tel = Telescope(resolution=32, diameter=8.0,
        central_obstruction=0.0, pupil_reflectivity=0.25)
    attenuated_pupil = PupilFunction(attenuated_tel)
    attenuated_wfs = ZernikeWFS(attenuated_tel; pupil_samples=8,
        diffraction_padding=2)
    attenuated_front_end = ZernikeOpticalFrontEnd(attenuated_wfs, source)
    attenuated_rate = zernike_rate_map(attenuated_front_end, attenuated_pupil)
    attenuated_plan = prepare_wfs_optics(attenuated_front_end, attenuated_pupil,
        attenuated_rate)
    form_wfs_optical_products!(attenuated_rate, attenuated_pupil, attenuated_plan)
    @test sum(attenuated_rate.values) ≈ 0.25 * sum(rate.values) rtol=1e-12
end

@testset "Zernike and Curvature source-composition support boundaries" begin
    telescope = Telescope(resolution=16, diameter=8.0,
        central_obstruction=0.0)
    pupil = PupilFunction(telescope)
    source = Source(band=:custom, wavelength=0.75e-6,
        photon_irradiance=1.0)
    spectral = with_spectrum(source,
        SpectralBundle([0.70e-6, 0.80e-6], [0.5, 0.5]))
    extended = with_extended_source(source,
        PointCloudSourceModel([(0.0, 0.0)], [1.0]))

    for expanded in (spectral, extended)
        zernike = ZernikeWFS(telescope; pupil_samples=2)
        front_end = ZernikeOpticalFrontEnd(zernike, expanded)
        rate = zernike_rate_map(front_end, pupil)
        @test_throws UnsupportedAlgorithm prepare_wfs_optics(
            front_end, pupil, rate)
        @test_throws UnsupportedAlgorithm measure!(
            CurvatureWFS(telescope; pupil_samples=2), pupil, expanded)
    end
end

@testset "Curvature pupil-reflectivity throughput" begin
    transmission = 0.25
    source = Source(band=:custom, wavelength=0.75e-6,
        photon_irradiance=1.0)
    full_telescope = Telescope(resolution=16, diameter=8.0,
        central_obstruction=0.0)
    attenuated_telescope = Telescope(resolution=16, diameter=8.0,
        central_obstruction=0.0, pupil_reflectivity=transmission)
    full_pupil = PupilFunction(full_telescope)
    attenuated_pupil = PupilFunction(attenuated_telescope)

    for style in (ScalarCPUStyle(), KA_CPU_STYLE)
        full_sensor = CurvatureWFS(full_telescope; pupil_samples=2,
            diffraction_padding=2)
        attenuated_sensor = CurvatureWFS(attenuated_telescope;
            pupil_samples=2, diffraction_padding=2)
        curvature_intensity!(style, full_sensor, full_pupil, source)
        curvature_intensity!(style, attenuated_sensor, attenuated_pupil,
            source)
        full_rate = sum(wfs_optical_rate_storage(full_sensor))
        @test full_rate > 0
        @test sum(wfs_optical_rate_storage(attenuated_sensor)) ≈
            transmission * full_rate rtol=1e-12
    end
end

@testset "Curvature diffractive photon-rate conservation" begin
    tel = Telescope(resolution=16, diameter=8.0,
        central_obstruction=0.0, pupil_reflectivity=0.25)
    pupil = PupilFunction(tel)
    src = Source(band=:custom, wavelength=0.75e-6,
        photon_irradiance=1.0)
    expected_two_branch_rate = 2 * sum(pupil_photon_rate_map(tel, src))

    for style in (ScalarCPUStyle(), KA_CPU_STYLE), padding in (1, 2, 3)
        wfs = CurvatureWFS(tel; pupil_samples=2,
            diffraction_padding=padding)
        curvature_intensity!(style, wfs, pupil, src)
        @test sum(wfs.front_end.propagation.workspace.intensity_stack) ≈
            expected_two_branch_rate atol=1e-10 rtol=1e-12
        if padding == 1
            @test sum(wfs_optical_rate_storage(wfs)) ≈
                expected_two_branch_rate atol=1e-10 rtol=1e-12
        else
            @test 0 < sum(wfs_optical_rate_storage(wfs)) <=
                expected_two_branch_rate
        end
    end
end

@testset "Curvature KernelAbstractions CPU stage parity" begin
    T = Float64
    tel = Telescope(resolution=8, diameter=T(8),
        central_obstruction=zero(T), T=T)
    pupil = PupilFunction(tel; T=T)
    pupil.opd .= reshape(T.(1:64), 8, 8) .* T(1e-10)
    source = Source(band=:custom, wavelength=T(0.75e-6),
        photon_irradiance=T(10), T=T)

    scalar_sensor = CurvatureWFS(tel; pupil_samples=2, T=T)
    accelerated_sensor = CurvatureWFS(tel; pupil_samples=2, T=T)
    @test backend(accelerated_sensor) isa CPUBackend
    scalar_front_end = CurvatureOpticalFrontEnd(scalar_sensor, source)
    accelerated_front_end = CurvatureOpticalFrontEnd(
        accelerated_sensor, source)

    scalar_fields = copy(
        WavefrontSensors._form_curvature_branch_fields!(
            ScalarCPUStyle(), scalar_front_end, pupil))
    accelerated_fields = copy(
        WavefrontSensors._form_curvature_branch_fields!(
            KA_CPU_STYLE, accelerated_front_end, pupil))
    @test accelerated_fields ≈ scalar_fields rtol=T(2e-12) atol=T(2e-12)

    scalar_rates = curvature_rate_maps(scalar_front_end, pupil)
    accelerated_rates = curvature_rate_maps(accelerated_front_end, pupil)
    scalar_workspace = scalar_sensor.front_end.propagation.workspace
    accelerated_workspace = accelerated_sensor.front_end.propagation.workspace
    intensity = reshape(T.(1:length(scalar_workspace.intensity_stack)),
        size(scalar_workspace.intensity_stack))
    copyto!(scalar_workspace.intensity_stack, intensity)
    copyto!(accelerated_workspace.intensity_stack, intensity)
    WavefrontSensors._sample_curvature_rate_planes!(
        ScalarCPUStyle(), scalar_rates, scalar_front_end)
    WavefrontSensors._sample_curvature_rate_planes!(
        KA_CPU_STYLE, accelerated_rates, accelerated_front_end)
    @test accelerated_rates[1].values ≈ scalar_rates[1].values
    @test accelerated_rates[2].values ≈ scalar_rates[2].values

    padding = scalar_sensor.front_end.propagation.plan.defocus_pair.diffraction_padding
    field = ElectricField(pupil, source; zero_padding=padding, T=T)
    field_plan = prepare_pupil_field(pupil, source, field;
        center_even_grid=false)
    fill_electric_field!(field, pupil, field_plan)
    scalar_field_stack = copy(
        WavefrontSensors._form_curvature_field_input!(
            ScalarCPUStyle(), scalar_workspace, field))
    accelerated_field_stack = copy(
        WavefrontSensors._form_curvature_field_input!(
            KA_CPU_STYLE, accelerated_workspace, field))
    @test accelerated_field_stack ≈ scalar_field_stack

    set_curvature_calibration!(accelerated_sensor, zeros(T, 2, 2);
        wavelength_m=wavelength(source), signature=UInt(0x4b414350))
    observation = WFSObservation(zeros(T, 2, 4);
        units=:photon_count, layout=:curvature_branch_channels)
    measurement = WFSMeasurement(zeros(T, 4);
        units=:dimensionless, kind=:curvature_signal)
    estimator = prepare_wfs_estimation(
        accelerated_sensor, observation, measurement)
    @test validate_wfs_target(estimator, HostComputeDevice()) === estimator
end

@testset "Curvature WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.0)
    pupil = PupilFunction(tel)
    src = Source(band=:I, magnitude=0.0)
    wfs = CurvatureWFS(tel; pupil_samples=8, defocus_rms_nm=500.0)

    @test size(wfs_optical_rate_storage(wfs)) == (16, 8)
    @test length(slopes(wfs)) == 64
    @test_throws InvalidConfiguration measure!(wfs, pupil)
    @test_throws InvalidConfiguration measure!(wfs, pupil,
        Asterism([src, Source(band=:I, magnitude=0.0)]))

    flat_slopes = copy(measure!(wfs, pupil, src))
    @test wfs.estimator.state.calibrated
    @test all(isfinite, flat_slopes)
    @test all(>=(0.0), wfs_optical_rate_storage(wfs))
    @test flat_slopes ≈ zero.(flat_slopes) atol=1e-10

    det = Detector(noise=NoiseNone(), binning=1)
    det_slopes = copy(measure!(wfs, pupil, src, det))
    @test det_slopes ≈ flat_slopes atol=1e-10
    @test size(output_frame(det)) == size(wfs_optical_rate_storage(wfs))
    @test wfs_detector_image(wfs, det) === output_frame(det)

    zb = ZernikeBasis(tel, 5)
    compute_zernike!(zb, tel)
    focus = @view zb.modes[:, :, 5]
    @. pupil.opd = 5e-8 * focus
    slopes_plus = copy(measure!(wfs, pupil, src))
    @. pupil.opd = -5e-8 * focus
    slopes_minus = copy(measure!(wfs, pupil, src))
    fill!(pupil.opd, 0.0)

    @test norm(slopes_plus) > 1e-6
    @test norm(slopes_minus) > 1e-6
    @test dot(slopes_plus, slopes_minus) < 0

    counting = CurvatureWFS(tel; pupil_samples=8,
        defocus_rms_nm=500.0, readout_model=CurvatureChannelReadout())
    counting_flat = copy(measure!(counting, pupil, src))
    @test size(wfs_optical_rate_storage(counting)) == (2, 64)
    @test counting_flat ≈ zero.(counting_flat) atol=1e-10
    @test_throws InvalidConfiguration measure!(counting, pupil, src, det)
    apd = LinearAPDDetector(topology=LinearAPDChannelBank(128),
        exposure_duration=1.0, qe=1.0, avalanche_gain=1.0,
        dark_current=0.0, noise=NoiseNone())
    counting_apd = copy(measure!(counting, pupil, src, apd))
    @test counting_apd ≈ counting_flat atol=1e-10
    @test detector_export_metadata(apd).n_channels ==
        length(wfs_optical_rate_storage(counting))
    spad = SPADArrayDetector(size(wfs_optical_rate_storage(counting));
        exposure_duration=1.0,
        noise=NoiseNone(),
        sensor=SPADArraySensor(active_area_detection_efficiency=1.0, dark_count_rate=0.0, fill_factor=1.0),
    )
    counting_spad = copy(measure!(counting, pupil, src, spad))
    @test counting_spad ≈ counting_flat atol=1e-10
    @test detector_export_metadata(spad).readout.output_size ==
        size(wfs_optical_rate_storage(counting))
    mkid = MKIDArrayDetector(
        exposure_duration=1.0,
        noise=NoiseNone(),
        sensor=MKIDArraySensor(qe=1.0, dark_count_rate=0.0, fill_factor=1.0,
            characteristics=MKIDArrayCharacteristics(
                wavelength_passband_m=(
                    0.9 * wavelength(src), 1.1 * wavelength(src)))),
    )
    counting_mkid = copy(measure!(counting, pupil, src, mkid))
    @test counting_mkid ≈ counting_flat atol=1e-10
    outside_mkid_band = Source(band=:custom, magnitude=0.0,
        wavelength=2 * wavelength(src), photon_irradiance=1.0)
    counting_mkid_outside = copy(measure!(counting, pupil,
        outside_mkid_band, mkid))
    @test all(iszero, output_frame(mkid))
    @test all(iszero, counting_mkid_outside)
    wrong_bank = LinearAPDDetector(topology=LinearAPDChannelBank(64),
        noise=NoiseNone())
    @test_throws InvalidConfiguration measure!(counting, pupil, src,
        wrong_bank)
    @test_throws InvalidConfiguration CurvatureWFS(tel; pupil_samples=8, readout_model=CurvatureChannelReadout(),
        readout_pixels_per_sample=2)

    response = CurvatureBranchResponse(T=Float64, plus_throughput=1.2, minus_throughput=0.8,
        plus_background=5.0, minus_background=1.0)
    imbalanced = CurvatureWFS(tel; pupil_samples=8, defocus_rms_nm=500.0, branch_response=response)
    imbalanced_flat = copy(measure!(imbalanced, pupil, src))
    @test imbalanced_flat ≈ zero.(imbalanced_flat) atol=1e-10
    plus_mean = mean(@view wfs_optical_rate_storage(imbalanced)[
        1:imbalanced.estimator.params.pupil_samples, :])
    minus_mean = mean(@view wfs_optical_rate_storage(imbalanced)[
        imbalanced.estimator.params.pupil_samples+1:end, :])
    @test plus_mean > minus_mean
    @test_throws InvalidConfiguration CurvatureBranchResponse(plus_throughput=-1.0)

    oversampled = CurvatureWFS(tel; pupil_samples=8, readout_crop_resolution=16, readout_pixels_per_sample=2)
    oversampled_flat = copy(measure!(oversampled, pupil, src))
    @test size(wfs_optical_rate_storage(oversampled)) == (32, 16)
    @test size(oversampled.front_end.propagation.workspace.frame_plus) == (16, 16)
    @test size(oversampled.estimator.workspace.reduced_plus) == (8, 8)
    @test oversampled_flat ≈ zero.(oversampled_flat) atol=1e-10
    @test_throws InvalidConfiguration CurvatureWFS(tel; pupil_samples=8, readout_crop_resolution=18, readout_pixels_per_sample=2)

    atm = MultiLayerAtmosphere(tel;
        r0=0.2,
        reference_wavelength_m=TEST_ATMOSPHERE_REFERENCE_WAVELENGTH_M,
        L0=25.0,
        fractional_cn2=[0.7, 0.3],
        wind_speed=[8.0, 4.0],
        wind_direction_deg=[0.0, 90.0],
        altitude=[0.0, 5000.0],
    )
    advance_by!(atm, TEST_ATMOSPHERE_STEP; rng=MersenneTwister(3))
    atm_slopes = copy(measure!(wfs, pupil, src, atm))
    @test all(isfinite, atm_slopes)
    @test norm(atm_slopes) > 0

    det_atm = Detector(noise=NoiseNone(), binning=1)
    det_atm_slopes = copy(measure!(wfs, pupil, src, atm, det_atm))
    @test all(isfinite, det_atm_slopes)

    ast = Asterism([src, Source(band=:I, magnitude=0.0, coordinates=(1.0, 90.0))])
    ast_slopes = copy(measure!(wfs, pupil, ast, atm))
    @test length(ast_slopes) == length(slopes(wfs))
    @test all(isfinite, ast_slopes)
    @test norm(ast_slopes) > 0

    common_qe_wavelength = 550e-9
    common_qe_ast = Asterism([
        Source(band=:custom, wavelength=common_qe_wavelength,
            photon_irradiance=1.0, coordinates=(0.0, 0.0)),
        Source(band=:custom, wavelength=common_qe_wavelength,
            photon_irradiance=2.0, coordinates=(0.5, 90.0)),
    ])
    wavelength_dependent_qe = SampledQuantumEfficiency(
        [500e-9, common_qe_wavelength, 600e-9], [0.1, 0.35, 0.9])
    sampled_qe_wfs = CurvatureWFS(tel; pupil_samples=8,
        defocus_rms_nm=500.0)
    sampled_qe_det = Detector(noise=NoiseNone(),
        qe=wavelength_dependent_qe, exposure_duration=1.0, binning=1)
    sampled_qe_slopes = copy(measure!(sampled_qe_wfs, pupil,
        common_qe_ast, atm, sampled_qe_det; rng=MersenneTwister(23)))
    sampled_qe_frame = copy(output_frame(sampled_qe_det))

    scalar_qe_wfs = CurvatureWFS(tel; pupil_samples=8,
        defocus_rms_nm=500.0)
    scalar_qe_det = Detector(noise=NoiseNone(), qe=0.35,
        exposure_duration=1.0, binning=1)
    scalar_qe_slopes = copy(measure!(scalar_qe_wfs, pupil,
        common_qe_ast, atm, scalar_qe_det; rng=MersenneTwister(23)))
    @test sum(sampled_qe_frame) > 0
    @test sampled_qe_frame ≈ output_frame(scalar_qe_det)
    @test sampled_qe_slopes ≈ scalar_qe_slopes

    mixed_qe_ast = Asterism([
        Source(band=:custom, wavelength=common_qe_wavelength,
            photon_irradiance=1.0),
        Source(band=:custom, wavelength=600e-9,
            photon_irradiance=1.0),
    ])
    @test_throws InvalidConfiguration measure!(
        CurvatureWFS(tel; pupil_samples=8, defocus_rms_nm=500.0),
        pupil, mixed_qe_ast, atm)
    @test_throws InvalidConfiguration measure!(
        CurvatureWFS(tel; pupil_samples=8, defocus_rms_nm=500.0),
        pupil, mixed_qe_ast, atm, sampled_qe_det)

    mixed_ngs_lgs = Asterism(AdaptiveOpticsSim.Optics.AbstractSource[
        Source(band=:custom, wavelength=589e-9,
            photon_irradiance=1.0),
        LGSSource(wavelength=589e-9, elongation_factor=1.3,
            photon_irradiance=1.0),
    ])
    @test_throws InvalidConfiguration measure!(
        CurvatureWFS(tel; pupil_samples=8, defocus_rms_nm=500.0),
        pupil, mixed_ngs_lgs, atm)
    @test_throws InvalidConfiguration measure!(
        CurvatureWFS(tel; pupil_samples=8, defocus_rms_nm=500.0),
        pupil, mixed_ngs_lgs, atm, sampled_qe_det)
end

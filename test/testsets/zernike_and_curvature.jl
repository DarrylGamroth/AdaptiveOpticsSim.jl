function zernike_normalization_allocations(normalization, wfs, tel, src,
    frame, normalization_scale)
    return @allocated zernike_normalization(normalization, wfs, tel, src,
        frame, normalization_scale)
end

function zernike_signal_allocations(wfs, tel, src, frame,
    normalization_scale)
    return @allocated zernike_signal!(wfs, tel, frame, src,
        normalization_scale)
end

@testset "Zernike WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    wfs = ZernikeWFS(tel; pupil_samples=8, diffraction_padding=2)

    @test size(wfs.state.camera_frame) == (8, 8)
    @test length(slopes(wfs)) == count(wfs.state.valid_mask)
    @test_throws InvalidConfiguration measure!(wfs, tel)
    @test_throws InvalidConfiguration measure!(wfs, tel, Asterism([src, Source(band=:I, magnitude=0.0)]))

    flat_slopes = copy(measure!(wfs, tel, src))
    @test wfs.state.calibrated
    @test all(isfinite, flat_slopes)
    @test all(>=(0.0), wfs.state.camera_frame)
    @test flat_slopes ≈ zero.(flat_slopes) atol=1e-10

    det = Detector(noise=NoiseNone(), binning=1)
    det_slopes = copy(measure!(wfs, tel, src, det))
    @test det_slopes ≈ flat_slopes atol=1e-10
    @test size(output_frame(det)) == size(wfs.state.camera_frame)
    @test wfs_detector_image(wfs, det) === output_frame(det)

    zb = ZernikeBasis(tel, 5)
    compute_zernike!(zb, tel)
    focus = @view zb.modes[:, :, 5]
    @. tel.state.opd = 5e-8 * focus
    slopes_plus = copy(measure!(wfs, tel, src))
    @. tel.state.opd = -5e-8 * focus
    slopes_minus = copy(measure!(wfs, tel, src))
    fill!(tel.state.opd, 0.0)

    @test norm(slopes_plus) > 1e-6
    @test norm(slopes_minus) > 1e-6
    @test dot(slopes_plus, slopes_minus) < 0
end

@testset "Zernike normalization hot-path contract" begin
    tel = Telescope(resolution=32, diameter=8.0,
        central_obstruction=0.0)
    src = Source(band=:custom, wavelength=0.75e-6,
        photon_irradiance=10.0)
    frame = reshape(collect(range(0.25, 2.0; length=64)), 8, 8)
    normalization_scale = 0.375

    for normalization in (MeanValidFluxNormalization(),
            IncidenceFluxNormalization())
        wfs = ZernikeWFS(tel; pupil_samples=8,
            normalization=normalization)
        fill!(wfs.state.reference_signal_2d, 0.0)
        valid = Array(wfs.state.valid_mask)
        expected = if normalization isa MeanValidFluxNormalization
            max(sum(frame[valid]) / count(valid), eps(eltype(frame)))
        else
            photon_rate = pupil_photon_rate_map(tel, src)
            nominal = similar(wfs.state.nominal_frame)
            sampled = similar(wfs.state.normalization_frame)
            sample_zernike_frame!(sampled, nominal, wfs, photon_rate, tel)
            sum(sampled[valid]) / count(valid) * normalization_scale
        end

        actual = zernike_normalization(normalization, wfs, tel, src,
            frame, normalization_scale)
        @test actual ≈ expected rtol=2e-15
        @test zernike_normalization_allocations(normalization, wfs, tel,
            src, frame, normalization_scale) == 0

        zernike_signal!(wfs, tel, frame, src, normalization_scale)
        @test zernike_signal_allocations(wfs, tel, src, frame,
            normalization_scale) == 0
        @test all(isfinite, slopes(wfs))
    end
end

@testset "Zernike incidence-normalization contract" begin
    tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.0)
    @inbounds for j in axes(tel.state.opd, 2), i in axes(tel.state.opd, 1)
        tel.state.opd[i, j] = 2e-8 * sinpi(2 * i / 32) * cospi(2 * j / 32)
    end
    src = Source(band=:custom, wavelength=0.75e-6,
        photon_irradiance=10.0)
    reference_wfs = ZernikeWFS(tel; pupil_samples=8,
        normalization=IncidenceFluxNormalization())
    asymmetric_response = [0.0 0.0 0.0;
                           0.0 0.2 0.8;
                           0.0 0.0 0.0]
    reference = copy(measure!(reference_wfs, tel, src,
        Detector(noise=NoiseNone(), integration_time=1.0, qe=1.0,
            response_model=SampledFrameResponse(asymmetric_response))))
    scaled_wfs = ZernikeWFS(tel; pupil_samples=8,
        normalization=IncidenceFluxNormalization())
    scaled = copy(measure!(scaled_wfs, tel, src,
        Detector(noise=NoiseNone(), integration_time=0.5, qe=0.25,
            response_model=SampledFrameResponse(asymmetric_response))))
    @test norm(reference) > 1e-6
    @test scaled ≈ reference atol=1e-12 rtol=1e-12

    zero_src = Source(band=:custom, wavelength=wavelength(src),
        photon_irradiance=0.0)
    for normalization in (MeanValidFluxNormalization(),
            IncidenceFluxNormalization())
        zero_wfs = ZernikeWFS(tel; pupil_samples=8,
            normalization=normalization)
        @test all(iszero, measure!(zero_wfs, tel, zero_src))
        @test all(isfinite, slopes(zero_wfs))
    end
end

@testset "Zernike detector-response references" begin
    tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.0)
    src = Source(band=:custom, wavelength=0.75e-6,
        photon_irradiance=10.0)
    asymmetric_response = [0.0 0.0 0.0;
                           0.0 0.2 0.8;
                           0.0 0.0 0.0]
    responses = (
        GaussianPixelResponse(response_width_px=0.75),
        SampledFrameResponse(asymmetric_response),
    )

    for response in responses
        wfs = ZernikeWFS(tel; pupil_samples=8)
        det = Detector(noise=NoiseNone(), integration_time=0.4, qe=0.3,
            response_model=response)
        flat = copy(measure!(wfs, tel, src, det))
        @test all(iszero, flat)
        @test wfs.state.calibration_signature ==
            detector_calibration_signature(det,
                telescope_aperture_calibration_signature(tel,
                    calibration_signature(src)))
        @test measure!(wfs, tel, src, det) == flat
    end

    probe_wfs = ZernikeWFS(tel; pupil_samples=8)
    probe_detector = Detector(noise=NoiseNone(), sensor=CMOSSensor())
    measure!(probe_wfs, tel, src, probe_detector;
        rng=MersenneTwister(705))
    detector_size = size(output_frame(probe_detector))
    gain_map = ones(detector_size)
    gain_map[1:2:end] .= 0.6
    bad_mask = falses(detector_size)
    bad_mask[2, 2] = true
    prnu = PixelResponseNonuniformity(gain_map)
    bad_pixels = BadPixelMask(bad_mask; throughput=0.0)
    detector = Detector(noise=NoiseNone(), sensor=CMOSSensor(),
        defect_model=CompositeDetectorDefectModel(prnu, bad_pixels))
    wfs = ZernikeWFS(tel; pupil_samples=8)

    flat = copy(measure!(wfs, tel, src, detector;
        rng=MersenneTwister(706)))
    @test all(iszero, flat)
    @test all(iszero, measure!(wfs, tel, src, detector;
        rng=MersenneTwister(707)))
    signature = wfs.state.calibration_signature

    prnu.gain_map[1, 1] = 0.2
    bad_pixels.mask[end, end] = true
    @test detector_calibration_signature(detector,
        telescope_aperture_calibration_signature(tel,
            calibration_signature(src))) == signature

    replacement_gain = copy(gain_map)
    replacement_gain[1, 1] = 0.35
    replacement_mask = copy(bad_mask)
    replacement_mask[end, end] = true
    replacement = Detector(noise=NoiseNone(), sensor=CMOSSensor(),
        defect_model=CompositeDetectorDefectModel(
            PixelResponseNonuniformity(replacement_gain),
            BadPixelMask(replacement_mask; throughput=0.0)))
    replacement_signature = detector_calibration_signature(replacement,
        telescope_aperture_calibration_signature(tel,
            calibration_signature(src)))
    @test replacement_signature != signature
    @test all(iszero, measure!(wfs, tel, src, replacement;
        rng=MersenneTwister(708)))
    @test wfs.state.calibration_signature == replacement_signature

    fill!(tel.state.opd, 1e-8)
    opd_before_failure = copy(tel.state.opd)
    @test_throws DimensionMismatchError measure!(
        ZernikeWFS(tel; pupil_samples=8), tel, src,
        Detector(noise=NoiseNone(), psf_sampling=3))
    @test tel.state.opd == opd_before_failure
end

@testset "Zernike and Curvature source-composition support boundaries" begin
    tel = Telescope(resolution=16, diameter=8.0,
        central_obstruction=0.0)
    src = Source(band=:custom, wavelength=0.75e-6,
        photon_irradiance=1.0)
    spectral = with_spectrum(src,
        SpectralBundle([0.70e-6, 0.80e-6], [0.5, 0.5]))
    extended = with_extended_source(src,
        PointCloudSourceModel([(0.0, 0.0)], [1.0]))

    for expanded in (spectral, extended)
        @test_throws UnsupportedAlgorithm measure!(
            ZernikeWFS(tel; pupil_samples=2), tel, expanded)
        @test_throws UnsupportedAlgorithm measure!(
            CurvatureWFS(tel; pupil_samples=2), tel, expanded)
    end
end

@testset "Zernike and curvature pupil-reflectivity throughput" begin
    transmission = 0.25
    src = Source(band=:custom, wavelength=0.75e-6,
        photon_irradiance=1.0)
    full_tel = Telescope(resolution=16, diameter=8.0,
        central_obstruction=0.0)
    attenuated_tel = Telescope(resolution=16, diameter=8.0,
        central_obstruction=0.0, pupil_reflectivity=transmission)

    full_zernike = ZernikeWFS(full_tel; pupil_samples=2,
        diffraction_padding=2)
    attenuated_zernike = ZernikeWFS(attenuated_tel; pupil_samples=2,
        diffraction_padding=2)
    zernike_pupil_intensity!(full_zernike, full_tel, src)
    zernike_pupil_intensity!(attenuated_zernike, attenuated_tel, src)
    full_zernike_rate = sum(full_zernike.state.pupil_intensity)
    @test full_zernike_rate > 0
    @test sum(attenuated_zernike.state.pupil_intensity) ≈
        transmission * full_zernike_rate rtol=1e-12

    for style in (ScalarCPUStyle(), KA_CPU_STYLE)
        full_curvature = CurvatureWFS(full_tel; pupil_samples=2,
            diffraction_padding=2)
        attenuated_curvature = CurvatureWFS(attenuated_tel; pupil_samples=2,
            diffraction_padding=2)
        curvature_intensity!(style, full_curvature, full_tel, src)
        curvature_intensity!(style, attenuated_curvature, attenuated_tel, src)
        full_curvature_rate = sum(full_curvature.state.camera_frame)
        @test full_curvature_rate > 0
        @test sum(attenuated_curvature.state.camera_frame) ≈
            transmission * full_curvature_rate rtol=1e-12
    end
end

@testset "Curvature diffractive photon-rate conservation" begin
    tel = Telescope(resolution=16, diameter=8.0,
        central_obstruction=0.0, pupil_reflectivity=0.25)
    src = Source(band=:custom, wavelength=0.75e-6,
        photon_irradiance=1.0)
    expected_two_branch_rate = 2 * sum(pupil_photon_rate_map(tel, src))

    for style in (ScalarCPUStyle(), KA_CPU_STYLE), padding in (1, 2, 3)
        wfs = CurvatureWFS(tel; pupil_samples=2,
            diffraction_padding=padding)
        curvature_intensity!(style, wfs, tel, src)
        @test sum(wfs.state.intensity_stack) ≈
            expected_two_branch_rate atol=1e-10 rtol=1e-12
        if padding == 1
            @test sum(wfs.state.camera_frame) ≈
                expected_two_branch_rate atol=1e-10 rtol=1e-12
        else
            @test 0 < sum(wfs.state.camera_frame) <= expected_two_branch_rate
        end
    end
end

@testset "Curvature WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    wfs = CurvatureWFS(tel; pupil_samples=8, defocus_rms_nm=500.0)

    @test size(wfs.state.camera_frame) == (16, 8)
    @test length(slopes(wfs)) == 64
    @test_throws InvalidConfiguration measure!(wfs, tel)
    @test_throws InvalidConfiguration measure!(wfs, tel, Asterism([src, Source(band=:I, magnitude=0.0)]))

    flat_slopes = copy(measure!(wfs, tel, src))
    @test wfs.state.calibrated
    @test all(isfinite, flat_slopes)
    @test all(>=(0.0), wfs.state.camera_frame)
    @test flat_slopes ≈ zero.(flat_slopes) atol=1e-10

    det = Detector(noise=NoiseNone(), binning=1)
    det_slopes = copy(measure!(wfs, tel, src, det))
    @test det_slopes ≈ flat_slopes atol=1e-10
    @test size(output_frame(det)) == size(wfs.state.camera_frame)
    @test wfs_detector_image(wfs, det) === output_frame(det)

    zb = ZernikeBasis(tel, 5)
    compute_zernike!(zb, tel)
    focus = @view zb.modes[:, :, 5]
    @. tel.state.opd = 5e-8 * focus
    slopes_plus = copy(measure!(wfs, tel, src))
    @. tel.state.opd = -5e-8 * focus
    slopes_minus = copy(measure!(wfs, tel, src))
    fill!(tel.state.opd, 0.0)

    @test norm(slopes_plus) > 1e-6
    @test norm(slopes_minus) > 1e-6
    @test dot(slopes_plus, slopes_minus) < 0

    counting = CurvatureWFS(tel; pupil_samples=8, defocus_rms_nm=500.0, readout_model=CurvatureCountingReadout())
    counting_flat = copy(measure!(counting, tel, src))
    @test size(counting.state.camera_frame) == (2, 64)
    @test counting_flat ≈ zero.(counting_flat) atol=1e-10
    @test_throws InvalidConfiguration measure!(counting, tel, src, det)
    apd = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0, noise=NoiseNone())
    counting_apd = copy(measure!(counting, tel, src, apd))
    @test counting_apd ≈ counting_flat atol=1e-10
    @test detector_export_metadata(apd).readout.output_size == size(counting.state.camera_frame)
    spad = SPADArrayDetector(
        integration_time=1.0,
        noise=NoiseNone(),
        sensor=SPADArraySensor(pde=1.0, dark_count_rate=0.0, fill_factor=1.0),
    )
    counting_spad = copy(measure!(counting, tel, src, spad))
    @test counting_spad ≈ counting_flat atol=1e-10
    @test detector_export_metadata(spad).readout.output_size == size(counting.state.camera_frame)
    mkid = MKIDArrayDetector(
        integration_time=1.0,
        noise=NoiseNone(),
        sensor=MKIDArraySensor(qe=1.0, dark_count_rate=0.0, fill_factor=1.0,
            wavelength_range_m=(0.9 * wavelength(src), 1.1 * wavelength(src))),
    )
    counting_mkid = copy(measure!(counting, tel, src, mkid))
    @test counting_mkid ≈ counting_flat atol=1e-10
    outside_mkid_band = Source(band=:custom, magnitude=0.0,
        wavelength=2 * wavelength(src), photon_irradiance=1.0)
    counting_mkid_outside = copy(measure!(counting, tel, outside_mkid_band, mkid))
    @test all(iszero, output_frame(mkid))
    @test all(iszero, counting_mkid_outside)
    apd_dead = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), dead_time_model=NonParalyzableDeadTime(0.25))
    counting_dead = copy(measure!(counting, tel, src, apd_dead))
    @test counting_dead ≈ counting_flat atol=1e-10
    @test_throws InvalidConfiguration CurvatureWFS(tel; pupil_samples=8, readout_model=CurvatureCountingReadout(),
        readout_pixels_per_sample=2)

    response = CurvatureBranchResponse(T=Float64, plus_throughput=1.2, minus_throughput=0.8,
        plus_background=5.0, minus_background=1.0)
    imbalanced = CurvatureWFS(tel; pupil_samples=8, defocus_rms_nm=500.0, branch_response=response)
    imbalanced_flat = copy(measure!(imbalanced, tel, src))
    @test imbalanced_flat ≈ zero.(imbalanced_flat) atol=1e-10
    plus_mean = mean(@view imbalanced.state.camera_frame[1:imbalanced.params.pupil_samples, :])
    minus_mean = mean(@view imbalanced.state.camera_frame[imbalanced.params.pupil_samples+1:end, :])
    @test plus_mean > minus_mean
    @test_throws InvalidConfiguration CurvatureBranchResponse(plus_throughput=-1.0)

    oversampled = CurvatureWFS(tel; pupil_samples=8, readout_crop_resolution=16, readout_pixels_per_sample=2)
    oversampled_flat = copy(measure!(oversampled, tel, src))
    @test size(oversampled.state.camera_frame) == (32, 16)
    @test size(oversampled.state.frame_plus) == (16, 16)
    @test size(oversampled.state.reduced_plus) == (8, 8)
    @test oversampled_flat ≈ zero.(oversampled_flat) atol=1e-10
    @test_throws InvalidConfiguration CurvatureWFS(tel; pupil_samples=8, readout_crop_resolution=18, readout_pixels_per_sample=2)

    atm = MultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[0.7, 0.3],
        wind_speed=[8.0, 4.0],
        wind_direction=[0.0, 90.0],
        altitude=[0.0, 5000.0],
    )
    advance_by!(atm, TEST_ATMOSPHERE_STEP; rng=MersenneTwister(3))
    atm_slopes = copy(measure!(wfs, tel, src, atm))
    @test all(isfinite, atm_slopes)
    @test norm(atm_slopes) > 0

    det_atm = Detector(noise=NoiseNone(), binning=1)
    det_atm_slopes = copy(measure!(wfs, tel, src, atm, det_atm))
    @test all(isfinite, det_atm_slopes)

    ast = Asterism([src, Source(band=:I, magnitude=0.0, coordinates=(1.0, 90.0))])
    ast_slopes = copy(measure!(wfs, tel, ast, atm))
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
        qe=wavelength_dependent_qe, integration_time=1.0, binning=1)
    sampled_qe_slopes = copy(measure!(sampled_qe_wfs, tel,
        common_qe_ast, atm, sampled_qe_det; rng=MersenneTwister(23)))
    sampled_qe_frame = copy(output_frame(sampled_qe_det))

    scalar_qe_wfs = CurvatureWFS(tel; pupil_samples=8,
        defocus_rms_nm=500.0)
    scalar_qe_det = Detector(noise=NoiseNone(), qe=0.35,
        integration_time=1.0, binning=1)
    scalar_qe_slopes = copy(measure!(scalar_qe_wfs, tel,
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
        tel, mixed_qe_ast, atm)
    @test_throws InvalidConfiguration measure!(
        CurvatureWFS(tel; pupil_samples=8, defocus_rms_nm=500.0),
        tel, mixed_qe_ast, atm, sampled_qe_det)

    mixed_ngs_lgs = Asterism(AdaptiveOpticsSim.AbstractSource[
        Source(band=:custom, wavelength=589e-9,
            photon_irradiance=1.0),
        LGSSource(wavelength=589e-9, elongation_factor=1.3,
            photon_irradiance=1.0),
    ])
    @test_throws InvalidConfiguration measure!(
        CurvatureWFS(tel; pupil_samples=8, defocus_rms_nm=500.0),
        tel, mixed_ngs_lgs, atm)
    @test_throws InvalidConfiguration measure!(
        CurvatureWFS(tel; pupil_samples=8, defocus_rms_nm=500.0),
        tel, mixed_ngs_lgs, atm, sampled_qe_det)
end

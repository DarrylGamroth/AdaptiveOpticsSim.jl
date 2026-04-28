@testset "Shack-Hartmann valid subaperture policies" begin
    tel = Telescope(resolution=352, diameter=1.22, sampling_time=1 / 500)
    sh_geom = ShackHartmannWFS(tel; n_lenslets=16, mode=Diffractive(), T=Float32)
    sh_flux = ShackHartmannWFS(tel;
        n_lenslets=16,
        mode=Diffractive(),
        valid_subaperture_policy=FluxThresholdValidSubapertures(light_ratio=0.5, T=Float32),
        T=Float32)

    geom_mask = copy(sh_geom.state.valid_mask_host)
    flux_mask = copy(sh_flux.state.valid_mask_host)

    @test sum(geom_mask) == 216
    @test sum(flux_mask) == 208
    @test flux_mask != geom_mask

    missing = sort!([[idx.I[1] - 1, idx.I[2] - 1] for idx in findall(geom_mask .& .!flux_mask)])
    @test missing == [[0, 4], [0, 11], [4, 0], [4, 15], [11, 0], [11, 15], [15, 4], [15, 11]]
    @test isempty(findall(flux_mask .& .!geom_mask))
end

@testset "Shack-Hartmann resized detector mosaic" begin
    tel = Telescope(resolution=64, diameter=8.0, sampling_time=1e-3, central_obstruction=0.1)
    src = Source(band=:I, magnitude=7.0)
    sh = ShackHartmannWFS(tel; n_lenslets=16, mode=Diffractive(), pixel_scale=0.06, n_pix_subap=8)
    prepare_runtime_wfs!(sh, tel, src)
    measure!(sh, tel, src)
    image = wfs_detector_image(sh)
    @test size(image) == (128, 128)
    @test size(AdaptiveOpticsSim.sh_exported_spot_cube(sh)) == (16 * 16, 8, 8)
end

@testset "Shack-Hartmann signal extraction branches" begin
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)

    zero_intensity = zeros(Float64, 3, 3)
    @test AdaptiveOpticsSim.centroid_from_intensity!(AdaptiveOpticsSim.ScalarCPUStyle(), zero_intensity, 0.1) ==
          (0.0, 0.0, 0.0)

    centroid_input = [1.0 0.0 0.0; 0.0 4.0 0.0; 0.0 0.0 9.0]
    total, sx, sy = AdaptiveOpticsSim.centroid_from_intensity_cutoff!(
        AdaptiveOpticsSim.ScalarCPUStyle(), copy(centroid_input), 3.0)
    @test total == 13.0
    @test sx ≈ 22 / 13
    @test sy ≈ 22 / 13
    @test AdaptiveOpticsSim.centroid_from_intensity!(copy(centroid_input), 0.25) ==
          AdaptiveOpticsSim.centroid_from_intensity!(AdaptiveOpticsSim.ScalarCPUStyle(), copy(centroid_input), 0.25)
    @test AdaptiveOpticsSim.centroid_from_intensity_cutoff!(copy(centroid_input), 3.0) ==
          AdaptiveOpticsSim.centroid_from_intensity_cutoff!(AdaptiveOpticsSim.ScalarCPUStyle(), copy(centroid_input), 3.0)
    @test AdaptiveOpticsSim.centroid_from_intensity!(KA_CPU_STYLE, copy(centroid_input), 0.25) ==
          AdaptiveOpticsSim.centroid_from_intensity!(AdaptiveOpticsSim.ScalarCPUStyle(), copy(centroid_input), 0.25)
    @test AdaptiveOpticsSim.centroid_from_intensity_cutoff!(KA_CPU_STYLE, copy(centroid_input), 3.0) ==
          AdaptiveOpticsSim.centroid_from_intensity_cutoff!(AdaptiveOpticsSim.ScalarCPUStyle(), copy(centroid_input), 3.0)
    @test AdaptiveOpticsSim.centroid_from_intensity_cutoff!(
        AdaptiveOpticsSim.ScalarCPUStyle(), fill(1.0, 2, 2), 2.0) == (0.0, 0.0, 0.0)

    sh = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), n_pix_subap=4, threshold_cog=0.2)
    @test AdaptiveOpticsSim.centroid_from_spot!(sh, copy(centroid_input), 0.25) ==
          AdaptiveOpticsSim.centroid_from_intensity!(AdaptiveOpticsSim.ScalarCPUStyle(), copy(centroid_input), 0.25)
    @test AdaptiveOpticsSim.centroid_from_spot!(sh, copy(centroid_input)) ==
          AdaptiveOpticsSim.centroid_from_intensity!(
              AdaptiveOpticsSim.ScalarCPUStyle(), copy(centroid_input), AdaptiveOpticsSim.centroid_threshold(sh))
    @test AdaptiveOpticsSim.centroid_from_spot_cutoff!(sh, copy(centroid_input), 3.0) ==
          AdaptiveOpticsSim.centroid_from_intensity_cutoff!(AdaptiveOpticsSim.ScalarCPUStyle(), copy(centroid_input), 3.0)
    @test AdaptiveOpticsSim.centroid_from_spot!(KA_CPU_STYLE, sh, copy(centroid_input), 0.25) ==
          AdaptiveOpticsSim.centroid_from_intensity!(AdaptiveOpticsSim.ScalarCPUStyle(), copy(centroid_input), 0.25)
    @test AdaptiveOpticsSim.centroid_from_spot_cutoff!(KA_CPU_STYLE, sh, copy(centroid_input), 3.0) ==
          AdaptiveOpticsSim.centroid_from_intensity_cutoff!(AdaptiveOpticsSim.ScalarCPUStyle(), copy(centroid_input), 3.0)
    spot_view = @view centroid_input[:, :]
    @test AdaptiveOpticsSim.sync_sh_staged_spot!(KA_CPU_STYLE, centroid_input) === centroid_input
    @test AdaptiveOpticsSim.sync_sh_staged_view!(KA_CPU_STYLE, spot_view) === spot_view

    fill!(sh.state.valid_mask, true)
    fill!(sh.state.valid_mask_host, true)
    sh.state.valid_mask[1, 3] = false
    sh.state.valid_mask_host[1, 3] = false
    fill!(sh.state.spot_cube, 0.0)
    sh.state.spot_cube[1, 2, 3] = 10.0
    sh.state.spot_cube[2, 1, 1] = 0.1
    sh.state.spot_cube[4, 4, 1] = 8.0

    scalar_slopes = AdaptiveOpticsSim.sh_signal_from_spots!(AdaptiveOpticsSim.ScalarCPUStyle(), sh, 0.5)
    offset = sh.params.n_lenslets * sh.params.n_lenslets
    @test scalar_slopes[1] == 2.0
    @test scalar_slopes[offset + 1] == 1.0
    @test scalar_slopes[2] == 0.0
    @test scalar_slopes[offset + 2] == 0.0
    @test scalar_slopes[3] == 0.0
    @test scalar_slopes[offset + 3] == 0.0
    @test AdaptiveOpticsSim.sh_signal_from_spots!(sh, 10.0, 0.05)[1] == 2.0
    @test AdaptiveOpticsSim.sh_signal_from_spots!(sh, 10.0, slope_extraction_model(sh))[1] == 2.0

    sh_accel = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    fill!(sh_accel.state.valid_mask, true)
    fill!(sh_accel.state.spot_cube, 0.0)
    sh_accel.state.spot_cube[1, 2, 3] = 10.0
    accel_slopes = copy(AdaptiveOpticsSim.sh_signal_from_spots!(KA_CPU_STYLE, sh_accel, 0.5))
    @test accel_slopes[1] == 2.0
    @test accel_slopes[offset + 1] == 1.0

    sh.state.reference_signal_2d .= reshape(collect(range(0.1, 3.2; length=2 * offset)), 2 * sh.params.n_lenslets, :)
    sh.state.slopes_units = 2.0
    fill!(sh.state.spot_cube, 0.0)
    sh.state.spot_cube[1, 2, 3] = 10.0
    calibrated = copy(AdaptiveOpticsSim.sh_signal_from_spots_calibrated!(AdaptiveOpticsSim.ScalarCPUStyle(), sh, 0.5))
    reference = vec(sh.state.reference_signal_2d)
    @test calibrated[1] ≈ (2.0 - reference[1]) / 2
    @test calibrated[offset + 1] ≈ (1.0 - reference[offset + 1]) / 2
    @test calibrated[3] ≈ -reference[3] / 2
    @test calibrated[offset + 3] ≈ -reference[offset + 3] / 2
    @test AdaptiveOpticsSim.sh_signal_from_spots_calibrated!(sh, 10.0, 0.05)[1] ≈
          (2.0 - reference[1]) / 2
    @test AdaptiveOpticsSim.sh_signal_from_spots_calibrated!(sh, 10.0, slope_extraction_model(sh))[1] ≈
          (2.0 - reference[1]) / 2

    sh_accel.state.reference_signal_2d .= 0.25
    sh_accel.state.slopes_units = 2.0
    fill!(sh_accel.state.spot_cube, 0.0)
    sh_accel.state.spot_cube[1, 2, 3] = 10.0
    accel_calibrated = copy(AdaptiveOpticsSim.sh_signal_from_spots_calibrated!(KA_CPU_STYLE, sh_accel, 0.5))
    @test accel_calibrated[1] ≈ (2.0 - 0.25) / 2
    @test accel_calibrated[offset + 1] ≈ (1.0 - 0.25) / 2

    slopes_for_invalid = ones(Float64, 2 * offset)
    AdaptiveOpticsSim.zero_invalid_sh_slopes!(AdaptiveOpticsSim.ScalarCPUStyle(), slopes_for_invalid, sh.state.valid_mask)
    @test slopes_for_invalid[3] == 0.0
    @test slopes_for_invalid[offset + 3] == 0.0
    @test slopes_for_invalid[1] == 1.0
    slopes_for_invalid_accel = ones(Float64, 2 * offset)
    AdaptiveOpticsSim.zero_invalid_sh_slopes!(KA_CPU_STYLE, slopes_for_invalid_accel, sh.state.valid_mask)
    @test slopes_for_invalid_accel[3] == 0.0
    @test slopes_for_invalid_accel[offset + 3] == 0.0
    @test slopes_for_invalid_accel[1] == 1.0

    cube_for_invalid = ones(Float64, offset, 2, 2)
    AdaptiveOpticsSim.zero_invalid_sh_spot_cube!(AdaptiveOpticsSim.ScalarCPUStyle(), cube_for_invalid, sh.state.valid_mask)
    @test all(iszero, cube_for_invalid[3, :, :])
    @test all(==(1.0), cube_for_invalid[1, :, :])
    cube_for_invalid_accel = ones(Float64, offset, 2, 2)
    AdaptiveOpticsSim.zero_invalid_sh_spot_cube!(KA_CPU_STYLE, cube_for_invalid_accel, sh.state.valid_mask)
    @test all(iszero, cube_for_invalid_accel[3, :, :])
    @test all(==(1.0), cube_for_invalid_accel[1, :, :])

    mean_signal = collect(1.0:(2 * offset))
    @test AdaptiveOpticsSim.mean_valid_signal(mean_signal, sh.state.valid_mask) ≈
          AdaptiveOpticsSim.packed_valid_pair_mean(AdaptiveOpticsSim.ScalarCPUStyle(), mean_signal, sh.state.valid_mask)
    @test AdaptiveOpticsSim.mean_valid_signal(KA_CPU_STYLE, mean_signal, sh.state.valid_mask) ≈
          AdaptiveOpticsSim.mean_valid_signal(mean_signal, sh.state.valid_mask)

    scalar_ramp = zeros(Float64, 8, 8)
    ka_ramp = similar(scalar_ramp)
    AdaptiveOpticsSim.fill_calibration_ramp!(AdaptiveOpticsSim.ScalarCPUStyle(), scalar_ramp, 1e-3, 8)
    AdaptiveOpticsSim.fill_calibration_ramp!(KA_CPU_STYLE, ka_ramp, 1e-3, 8)
    @test ka_ramp ≈ scalar_ramp

    sh_scalar = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    sh_ka = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    for wfs in (sh_scalar, sh_ka)
        fill!(wfs.state.valid_mask, true)
        fill!(wfs.state.spot_cube, 0.0)
        wfs.state.spot_cube[1, 2, 3] = 10.0
        wfs.state.spot_cube[4, 4, 1] = 8.0
    end
    scalar_device_reference = copy(AdaptiveOpticsSim.sh_signal_from_spots!(AdaptiveOpticsSim.ScalarCPUStyle(), sh_scalar, 0.5))
    ka_device_stats = copy(AdaptiveOpticsSim.sh_signal_from_spots_device_stats!(KA_CPU_STYLE, sh_ka, 0.5))
    @test ka_device_stats ≈ scalar_device_reference

    sh_scalar.state.reference_signal_2d .= 0.25
    sh_ka.state.reference_signal_2d .= 0.25
    sh_scalar.state.slopes_units = 2.0
    sh_ka.state.slopes_units = 2.0
    fill!(sh_scalar.state.spot_cube, 0.0)
    fill!(sh_ka.state.spot_cube, 0.0)
    sh_scalar.state.spot_cube[1, 2, 3] = 10.0
    sh_ka.state.spot_cube[1, 2, 3] = 10.0
    scalar_calibrated = copy(AdaptiveOpticsSim.sh_signal_from_spots_calibrated!(
        AdaptiveOpticsSim.ScalarCPUStyle(), sh_scalar, 0.5))
    ka_calibrated = copy(AdaptiveOpticsSim.sh_signal_from_spots_calibrated_device_stats!(KA_CPU_STYLE, sh_ka, 0.5))
    @test ka_calibrated ≈ scalar_calibrated

    det = Detector(noise=NoiseNone(), binning=1)
    sh_det = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    point_peak = AdaptiveOpticsSim.sampled_spots_peak!(sh_det, tel, src, det, MersenneTwister(21))
    @test isfinite(point_peak)
    @test point_peak > 0
    sh_det_accel = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    AdaptiveOpticsSim.prepare_sampling!(sh_det_accel, tel, src)
    point_peak_accel = AdaptiveOpticsSim.sampled_spots_peak!(KA_CPU_STYLE, sh_det_accel, tel, src, det, MersenneTwister(21))
    @test point_peak_accel ≈ point_peak

    poly = with_spectrum(src, SpectralBundle([SpectralSample(0.95 * wavelength(src), 0.5),
        SpectralSample(1.05 * wavelength(src), 0.5)]))
    sh_poly = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    poly_peak = AdaptiveOpticsSim.sampled_spots_peak!(sh_poly, tel, poly, det, MersenneTwister(22))
    @test isfinite(poly_peak)
    @test poly_peak > 0
    sh_poly_accel = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    AdaptiveOpticsSim.prepare_sampling!(sh_poly_accel, tel, src)
    poly_peak_accel = AdaptiveOpticsSim.sampled_spots_peak!(KA_CPU_STYLE, sh_poly_accel, tel, poly, det, MersenneTwister(22))
    @test poly_peak_accel ≈ poly_peak

    ext_single = with_extended_source(src, PointCloudSourceModel([(0.0, 0.0)], [1.0]))
    sh_ext_single = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    ext_single_peak = AdaptiveOpticsSim.sampled_spots_peak!(sh_ext_single, tel, ext_single, det, MersenneTwister(23))
    @test ext_single_peak ≈ point_peak atol=1e-8 rtol=1e-8
    ext_double = with_extended_source(src, PointCloudSourceModel([(0.0, 0.0), (0.1, 0.0)], [0.5, 0.5]))
    sh_ext_double_accel = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    AdaptiveOpticsSim.prepare_sampling!(sh_ext_double_accel, tel, src)
    ext_double_peak_accel = AdaptiveOpticsSim.sampled_spots_peak!(
        KA_CPU_STYLE, sh_ext_double_accel, tel, ext_double, det, MersenneTwister(23))
    @test isfinite(ext_double_peak_accel)
    @test ext_double_peak_accel > 0

    ast_batched = Asterism([src, Source(band=:I, magnitude=0.0, coordinates=(0.1, 0.0))])
    sh_ast_batched = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    AdaptiveOpticsSim.prepare_sampling!(sh_ast_batched, tel, src)
    AdaptiveOpticsSim.ensure_sh_calibration!(sh_ast_batched, tel, src)
    ast_batched_slopes = AdaptiveOpticsSim.measure_sh_asterism_batched!(KA_CPU_STYLE, sh_ast_batched, tel, ast_batched)
    @test length(ast_batched_slopes) == 2 * 4 * 4
    @test all(isfinite, ast_batched_slopes)
    sh_ast_batched_det = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    AdaptiveOpticsSim.prepare_sampling!(sh_ast_batched_det, tel, src)
    AdaptiveOpticsSim.ensure_sh_calibration!(sh_ast_batched_det, tel, src)
    ast_batched_det_slopes = AdaptiveOpticsSim.measure_sh_asterism_batched!(
        KA_CPU_STYLE, sh_ast_batched_det, tel, ast_batched, det, MersenneTwister(26))
    @test length(ast_batched_det_slopes) == 2 * 4 * 4
    @test all(isfinite, ast_batched_det_slopes)

    lgs = LGSSource(elongation_factor=1.3)
    sh_lgs = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    lgs_peak = AdaptiveOpticsSim.sampled_spots_peak!(sh_lgs, tel, lgs, det, MersenneTwister(24))
    @test isfinite(lgs_peak)
    @test lgs_peak >= 0
    sh_lgs_accel = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    AdaptiveOpticsSim.prepare_sampling!(sh_lgs_accel, tel, lgs)
    lgs_peak_accel = AdaptiveOpticsSim.sampled_spots_peak!(KA_CPU_STYLE, sh_lgs_accel, tel, lgs, det, MersenneTwister(24))
    @test isfinite(lgs_peak_accel)
    @test lgs_peak_accel >= 0

    na_profile = [80000.0 90000.0 100000.0; 0.2 0.6 0.2]
    lgs_profile = LGSSource(elongation_factor=1.2, na_profile=na_profile, fwhm_spot_up=1.0)
    sh_lgs_profile_accel = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    AdaptiveOpticsSim.prepare_sampling!(sh_lgs_profile_accel, tel, lgs_profile)
    @test AdaptiveOpticsSim.sampled_spots_peak!(KA_CPU_STYLE, sh_lgs_profile_accel, tel, lgs_profile) >= 0
    sh_lgs_profile_det_accel = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    AdaptiveOpticsSim.prepare_sampling!(sh_lgs_profile_det_accel, tel, lgs_profile)
    @test AdaptiveOpticsSim.sampled_spots_peak!(
        KA_CPU_STYLE, sh_lgs_profile_det_accel, tel, lgs_profile, det, MersenneTwister(25)) >= 0
end

@testset "Asterism PSF" begin
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src1 = Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0))
    src2 = Source(band=:I, magnitude=0.0, coordinates=(1.0, 90.0))
    @test coordinates_xy_arcsec(src1) == (0.0, 0.0)
    @test coordinates_xy_arcsec(src2)[1] ≈ 0.0 atol=1e-12
    @test coordinates_xy_arcsec(src2)[2] ≈ 1.0
    ast = Asterism([src1, src2])
    psf = compute_psf!(tel, ast; zero_padding=2)
    @test size(tel.state.psf_stack, 3) == 2
    @test size(psf) == (32, 32)
    @test sum(psf) >= sum(@view tel.state.psf_stack[:, :, 1])
end

@testset "Polychromatic WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    λ0 = wavelength(src)
    bundle_single = SpectralBundle([SpectralSample(λ0, 1.0)])
    bundle_broad = SpectralBundle([SpectralSample(0.9 * λ0, 0.4), SpectralSample(1.1 * λ0, 0.6)])
    poly_single = with_spectrum(src, bundle_single)
    poly_broad = with_spectrum(src, bundle_broad)

    @test sum(sample.weight for sample in bundle_broad) ≈ 1.0
    @test weighted_wavelength(bundle_broad) ≈ (0.9 * λ0 * 0.4 + 1.1 * λ0 * 0.6)
    @test has_spectral_bundle(poly_broad)
    @test is_polychromatic(poly_broad)
    @test !is_polychromatic(poly_single)
    @test spectral_reference_source(poly_broad) === src

    zb = ZernikeBasis(tel, 5)
    compute_zernike!(zb, tel)
    focus = @view zb.modes[:, :, 5]
    @. tel.state.opd = 5e-8 * focus

    sh_mono = ShackHartmannWFS(tel; n_lenslets=8, mode=Diffractive())
    sh_single = ShackHartmannWFS(tel; n_lenslets=8, mode=Diffractive())
    sh_broad = ShackHartmannWFS(tel; n_lenslets=8, mode=Diffractive())

    mono_slopes = copy(measure!(sh_mono, tel, src))
    single_slopes = copy(measure!(sh_single, tel, poly_single))
    broad_slopes_1 = copy(measure!(sh_broad, tel, poly_broad))
    broad_slopes_2 = copy(measure!(sh_broad, tel, poly_broad))

    @test single_slopes ≈ mono_slopes atol=1e-10 rtol=1e-10
    @test broad_slopes_1 ≈ broad_slopes_2 atol=1e-10 rtol=1e-10
    @test norm(broad_slopes_1 - mono_slopes) > 1e-8
    @test supports_stacked_sources(sh_broad, poly_broad)

    pyr_mono = PyramidWFS(tel; pupil_samples=8, mode=Diffractive(), modulation=1.0)
    pyr_single = PyramidWFS(tel; pupil_samples=8, mode=Diffractive(), modulation=1.0)
    pyr_broad = PyramidWFS(tel; pupil_samples=8, mode=Diffractive(), modulation=1.0)

    mono_pyr = copy(measure!(pyr_mono, tel, src))
    single_pyr = copy(measure!(pyr_single, tel, poly_single))
    broad_pyr_1 = copy(measure!(pyr_broad, tel, poly_broad))
    broad_pyr_2 = copy(measure!(pyr_broad, tel, poly_broad))

    @test single_pyr ≈ mono_pyr atol=1e-10 rtol=1e-10
    @test broad_pyr_1 ≈ broad_pyr_2 atol=1e-10 rtol=1e-10
    @test norm(broad_pyr_1 - mono_pyr) > 1e-8
    @test supports_stacked_sources(pyr_broad, poly_broad)

    det = Detector(noise=NoiseNone(), binning=1)
    spectral_frame = measure!(sh_broad, tel, poly_broad, det)
    @test size(sh_broad.state.detector_noise_cube) == size(sh_broad.state.spot_cube)
    @test spectral_frame ≈ broad_slopes_1 atol=1e-10 rtol=1e-10

    pyr_det = Detector(noise=NoiseNone(), binning=1)
    pyr_det_slopes = measure!(pyr_broad, tel, poly_broad, pyr_det)
    @test size(output_frame(pyr_det)) == size(pyr_broad.state.camera_frame)
    @test pyr_det_slopes ≈ broad_pyr_1 atol=1e-10 rtol=1e-10

    fill!(tel.state.opd, 0.0)
end

@testset "Extended-source WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    point_model = PointCloudSourceModel([(0.0, 0.0)], [1.0])
    gaussian_model = GaussianDiskSourceModel(sigma_arcsec=0.35, n_side=5)
    image_model = SampledImageSourceModel([0.0 1.0 0.0; 1.0 2.0 1.0; 0.0 1.0 0.0], pixel_scale_arcsec=0.2)
    ext_point = with_extended_source(src, point_model)
    ext_gauss = with_extended_source(src, gaussian_model)
    ext_image = with_extended_source(src, image_model)

    point_ast = extended_source_asterism(ext_point)
    image_ast = extended_source_asterism(ext_image)
    @test has_extended_source_model(ext_gauss)
    @test !has_extended_source_model(src)
    @test length(point_ast) == 1
    @test length(image_ast) == 5
    @test AdaptiveOpticsSim.photon_flux(point_ast.sources[1]) ≈ AdaptiveOpticsSim.photon_flux(src)
    @test sum(AdaptiveOpticsSim.photon_flux(sample) for sample in image_ast.sources) ≈ AdaptiveOpticsSim.photon_flux(src)

    zb = ZernikeBasis(tel, 5)
    compute_zernike!(zb, tel)
    focus = @view zb.modes[:, :, 5]
    @. tel.state.opd = 5e-8 * focus

    sh_point = ShackHartmannWFS(tel; n_lenslets=8, mode=Diffractive())
    sh_ext_point = ShackHartmannWFS(tel; n_lenslets=8, mode=Diffractive())
    sh_ext = ShackHartmannWFS(tel; n_lenslets=8, mode=Diffractive())
    point_slopes = copy(measure!(sh_point, tel, src))
    ext_point_slopes = copy(measure!(sh_ext_point, tel, ext_point))
    point_peak = AdaptiveOpticsSim.sampled_spots_peak!(sh_point, tel, src)
    ext_peak = AdaptiveOpticsSim.sampled_spots_peak!(sh_ext, tel, ext_gauss)
    ext_slopes_1 = copy(measure!(sh_ext, tel, ext_gauss))
    ext_slopes_2 = copy(measure!(sh_ext, tel, ext_gauss))

    @test ext_point_slopes ≈ point_slopes atol=1e-10 rtol=1e-10
    @test ext_slopes_1 ≈ ext_slopes_2 atol=1e-10 rtol=1e-10
    @test ext_peak < point_peak
    @test norm(sh_ext.state.spot_cube - sh_point.state.spot_cube) > 1e-8
    @test supports_stacked_sources(sh_ext, ext_gauss)

    pyr_point = PyramidWFS(tel; pupil_samples=8, mode=Diffractive(), modulation=1.0)
    pyr_ext_point = PyramidWFS(tel; pupil_samples=8, mode=Diffractive(), modulation=1.0)
    pyr_ext = PyramidWFS(tel; pupil_samples=8, mode=Diffractive(), modulation=1.0)
    pyr_point_slopes = copy(measure!(pyr_point, tel, src))
    pyr_ext_point_slopes = copy(measure!(pyr_ext_point, tel, ext_point))
    pyr_ext_slopes_1 = copy(measure!(pyr_ext, tel, ext_gauss))
    pyr_ext_slopes_2 = copy(measure!(pyr_ext, tel, ext_gauss))

    @test pyr_ext_point_slopes ≈ pyr_point_slopes atol=1e-10 rtol=1e-10
    @test pyr_ext_slopes_1 ≈ pyr_ext_slopes_2 atol=1e-10 rtol=1e-10
    @test norm(pyr_ext.state.intensity - pyr_point.state.intensity) > 1e-10
    @test supports_stacked_sources(pyr_ext, ext_gauss)

    det = Detector(noise=NoiseNone(), binning=1)
    sh_det_slopes = measure!(sh_ext, tel, ext_gauss, det)
    @test size(sh_ext.state.detector_noise_cube) == size(sh_ext.state.spot_cube)
    @test sh_det_slopes ≈ ext_slopes_1 atol=1e-10 rtol=1e-10

    pyr_det = Detector(noise=NoiseNone(), binning=1)
    pyr_det_slopes = measure!(pyr_ext, tel, ext_gauss, pyr_det)
    @test size(output_frame(pyr_det)) == size(pyr_ext.state.camera_frame)
    @test pyr_det_slopes ≈ pyr_ext_slopes_1 atol=1e-10 rtol=1e-10

    fill!(tel.state.opd, 0.0)
end

@testset "Pupil masks and misregistration" begin
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    base_sum = sum(tel.state.pupil)
    apply_spiders!(tel; thickness=0.5, angles=[0.0, 90.0])
    @test sum(tel.state.pupil) < base_sum

    custom = trues(16, 16)
    custom[:, 9:end] .= false
    set_pupil!(tel, custom)
    @test sum(tel.state.pupil) == sum(custom)
    @test tel.state.pupil_reflectivity == Float64.(custom)

    reflectivity = fill(0.5, 16, 16)
    set_pupil_reflectivity!(tel, reflectivity)
    @test tel.state.pupil_reflectivity[:, 1:8] == fill(0.5, 16, 8)
    @test tel.state.pupil_reflectivity[:, 9:end] == fill(0.0, 16, 8)

    tel2 = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm1 = DeformableMirror(tel2; n_act=2, influence_width=0.3)
    mis = Misregistration(shift_x=0.1, shift_y=0.0, rotation_deg=5.0, T=Float64)
    @test rotation_deg(mis) ≈ 5.0
    @test rotation_rad(mis) ≈ deg2rad(5.0)
    dm2 = DeformableMirror(tel2; n_act=2, influence_width=0.3, misregistration=mis)
    @test dm1.state.modes != dm2.state.modes
end

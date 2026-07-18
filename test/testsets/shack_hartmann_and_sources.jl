function sh_spectral_grid_guard_allocation_bytes(wfs, src)
    AdaptiveOpticsSim.require_sh_common_spectral_grid(wfs, src)
    return @allocated AdaptiveOpticsSim.require_sh_common_spectral_grid(
        wfs, src)
end

@testset "Shack-Hartmann valid subaperture policies" begin
    tel = Telescope(resolution=352, diameter=1.22)
    sh_geom = ShackHartmannWFS(tel; n_lenslets=16, mode=Diffractive(), T=Float32)
    sh_flux = ShackHartmannWFS(tel;
        n_lenslets=16,
        mode=Diffractive(),
        valid_subaperture_policy=FluxThresholdValidSubapertures(light_ratio=0.5, T=Float32),
        T=Float32)

    geom_mask = copy(sh_geom.front_end.layout.valid_mask_host)
    flux_mask = copy(sh_flux.front_end.layout.valid_mask_host)

    @test sum(geom_mask) == 216
    @test sum(flux_mask) == 208
    @test flux_mask != geom_mask

    missing = sort!([[idx.I[1] - 1, idx.I[2] - 1] for idx in findall(geom_mask .& .!flux_mask)])
    @test missing == [[0, 4], [0, 11], [4, 0], [4, 15], [11, 0], [11, 15], [15, 4], [15, 11]]
    @test isempty(findall(flux_mask .& .!geom_mask))
end

@testset "Shack-Hartmann resized detector mosaic" begin
    tel = Telescope(resolution=64, diameter=8.0, central_obstruction=0.1)
    src = Source(band=:I, magnitude=7.0)
    sh = ShackHartmannWFS(tel; n_lenslets=16, mode=Diffractive(), pixel_scale_arcsec=0.06, n_pix_subap=8)
    prepare_runtime_wfs!(sh, tel, src)
    measure!(sh, tel, src)
    image = wfs_detector_image(sh)
    @test size(image) == (128, 128)
    @test size(AdaptiveOpticsSim.sh_exported_spot_cube(sh)) == (16 * 16, 8, 8)
end

@testset "Shack-Hartmann signal extraction branches" begin
    tel = Telescope(resolution=16, diameter=8.0, central_obstruction=0.0)
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

    fill!(sh.front_end.layout.valid_mask, true)
    fill!(sh.front_end.layout.valid_mask_host, true)
    sh.front_end.layout.valid_mask[1, 3] = false
    sh.front_end.layout.valid_mask_host[1, 3] = false
    invalid_index = LinearIndices(sh.front_end.layout.valid_mask)[1, 3]
    fill!(sh.acquisition.spot_cube, 0.0)
    sh.acquisition.spot_cube[1, 2, 3] = 10.0
    sh.acquisition.spot_cube[2, 1, 1] = 0.1
    sh.acquisition.spot_cube[4, 4, 1] = 8.0

    scalar_slopes = AdaptiveOpticsSim.sh_signal_from_spots!(AdaptiveOpticsSim.ScalarCPUStyle(), sh, 0.5)
    n_lenslets = microlens_array(sh.front_end).params.n_lenslets
    offset = n_lenslets * n_lenslets
    @test scalar_slopes[1] == 1.0
    @test scalar_slopes[offset + 1] == 2.0
    @test scalar_slopes[2] == 0.0
    @test scalar_slopes[offset + 2] == 0.0
    @test scalar_slopes[invalid_index] == 0.0
    @test scalar_slopes[offset + invalid_index] == 0.0
    @test AdaptiveOpticsSim.sh_signal_from_spots!(sh, 10.0, 0.05)[1] == 1.0
    @test AdaptiveOpticsSim.sh_signal_from_spots!(sh, 10.0, slope_extraction_model(sh))[1] == 1.0

    sh_accel = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    fill!(sh_accel.front_end.layout.valid_mask, true)
    fill!(sh_accel.acquisition.spot_cube, 0.0)
    sh_accel.acquisition.spot_cube[1, 2, 3] = 10.0
    accel_slopes = copy(AdaptiveOpticsSim.sh_signal_from_spots!(KA_CPU_STYLE, sh_accel, 0.5))
    @test accel_slopes[1] == 1.0
    @test accel_slopes[offset + 1] == 2.0

    sh.calibration.reference_signal_2d .= reshape(
        collect(range(0.1, 3.2; length=2 * offset)), offset, 2)
    sh.calibration.centroid_response = 2.0
    fill!(sh.acquisition.spot_cube, 0.0)
    sh.acquisition.spot_cube[1, 2, 3] = 10.0
    calibrated = copy(AdaptiveOpticsSim.sh_signal_from_spots_calibrated!(AdaptiveOpticsSim.ScalarCPUStyle(), sh, 0.5))
    reference = vec(sh.calibration.reference_signal_2d)
    @test calibrated[1] ≈ (1.0 - reference[1]) / 2
    @test calibrated[offset + 1] ≈ (2.0 - reference[offset + 1]) / 2
    @test calibrated[invalid_index] ≈ -reference[invalid_index] / 2
    @test calibrated[offset + invalid_index] ≈
        -reference[offset + invalid_index] / 2
    @test AdaptiveOpticsSim.sh_signal_from_spots_calibrated!(sh, 10.0, 0.05)[1] ≈
          (1.0 - reference[1]) / 2
    @test AdaptiveOpticsSim.sh_signal_from_spots_calibrated!(sh, 10.0, slope_extraction_model(sh))[1] ≈
          (1.0 - reference[1]) / 2

    sh_accel.calibration.reference_signal_2d .= 0.25
    sh_accel.calibration.centroid_response = 2.0
    fill!(sh_accel.acquisition.spot_cube, 0.0)
    sh_accel.acquisition.spot_cube[1, 2, 3] = 10.0
    accel_calibrated = copy(AdaptiveOpticsSim.sh_signal_from_spots_calibrated!(KA_CPU_STYLE, sh_accel, 0.5))
    @test accel_calibrated[1] ≈ (1.0 - 0.25) / 2
    @test accel_calibrated[offset + 1] ≈ (2.0 - 0.25) / 2

    slopes_for_invalid = ones(Float64, 2 * offset)
    AdaptiveOpticsSim.zero_invalid_sh_slopes!(AdaptiveOpticsSim.ScalarCPUStyle(), slopes_for_invalid, sh.front_end.layout.valid_mask)
    @test slopes_for_invalid[invalid_index] == 0.0
    @test slopes_for_invalid[offset + invalid_index] == 0.0
    @test slopes_for_invalid[1] == 1.0
    slopes_for_invalid_accel = ones(Float64, 2 * offset)
    AdaptiveOpticsSim.zero_invalid_sh_slopes!(KA_CPU_STYLE, slopes_for_invalid_accel, sh.front_end.layout.valid_mask)
    @test slopes_for_invalid_accel[invalid_index] == 0.0
    @test slopes_for_invalid_accel[offset + invalid_index] == 0.0
    @test slopes_for_invalid_accel[1] == 1.0

    cube_for_invalid = ones(Float64, offset, 2, 2)
    AdaptiveOpticsSim.zero_invalid_sh_spot_cube!(AdaptiveOpticsSim.ScalarCPUStyle(), cube_for_invalid, sh.front_end.layout.valid_mask)
    @test all(iszero, cube_for_invalid[invalid_index, :, :])
    @test all(==(1.0), cube_for_invalid[1, :, :])
    cube_for_invalid_accel = ones(Float64, offset, 2, 2)
    AdaptiveOpticsSim.zero_invalid_sh_spot_cube!(KA_CPU_STYLE, cube_for_invalid_accel, sh.front_end.layout.valid_mask)
    @test all(iszero, cube_for_invalid_accel[invalid_index, :, :])
    @test all(==(1.0), cube_for_invalid_accel[1, :, :])

    mean_signal = collect(1.0:(2 * offset))
    @test AdaptiveOpticsSim.mean_valid_signal(mean_signal, sh.front_end.layout.valid_mask) ≈
          AdaptiveOpticsSim.packed_valid_pair_mean(AdaptiveOpticsSim.ScalarCPUStyle(), mean_signal, sh.front_end.layout.valid_mask)
    @test AdaptiveOpticsSim.mean_valid_signal(KA_CPU_STYLE, mean_signal, sh.front_end.layout.valid_mask) ≈
          AdaptiveOpticsSim.mean_valid_signal(mean_signal, sh.front_end.layout.valid_mask)

    scalar_ramp = zeros(Float64, 8, 8)
    ka_ramp = similar(scalar_ramp)
    AdaptiveOpticsSim.fill_calibration_ramp!(AdaptiveOpticsSim.ScalarCPUStyle(), scalar_ramp, 1e-3, 8)
    AdaptiveOpticsSim.fill_calibration_ramp!(KA_CPU_STYLE, ka_ramp, 1e-3, 8)
    @test ka_ramp ≈ scalar_ramp

    sh_scalar = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    sh_ka = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    for wfs in (sh_scalar, sh_ka)
        fill!(wfs.front_end.layout.valid_mask, true)
        fill!(wfs.acquisition.spot_cube, 0.0)
        wfs.acquisition.spot_cube[1, 2, 3] = 10.0
        wfs.acquisition.spot_cube[4, 4, 1] = 8.0
    end
    scalar_device_reference = copy(AdaptiveOpticsSim.sh_signal_from_spots!(AdaptiveOpticsSim.ScalarCPUStyle(), sh_scalar, 0.5))
    ka_device_stats = copy(AdaptiveOpticsSim.sh_signal_from_spots_device_stats!(KA_CPU_STYLE, sh_ka, 0.5))
    @test ka_device_stats ≈ scalar_device_reference

    sh_scalar.calibration.reference_signal_2d .= 0.25
    sh_ka.calibration.reference_signal_2d .= 0.25
    sh_scalar.calibration.centroid_response = 2.0
    sh_ka.calibration.centroid_response = 2.0
    fill!(sh_scalar.acquisition.spot_cube, 0.0)
    fill!(sh_ka.acquisition.spot_cube, 0.0)
    sh_scalar.acquisition.spot_cube[1, 2, 3] = 10.0
    sh_ka.acquisition.spot_cube[1, 2, 3] = 10.0
    scalar_calibrated = copy(AdaptiveOpticsSim.sh_signal_from_spots_calibrated!(
        AdaptiveOpticsSim.ScalarCPUStyle(), sh_scalar, 0.5))
    ka_calibrated = copy(AdaptiveOpticsSim.sh_signal_from_spots_calibrated_device_stats!(KA_CPU_STYLE, sh_ka, 0.5))
    @test ka_calibrated ≈ scalar_calibrated

    det = Detector(noise=NoiseNone(), binning=1)
    sh_det = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    point_peak = AdaptiveOpticsSim.sampled_spots_peak!(sh_det, tel, src, det, MersenneTwister(21))
    point_spots = copy(sh_det.acquisition.spot_cube)
    @test isfinite(point_peak)
    @test point_peak > 0
    sh_det_accel = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    AdaptiveOpticsSim.prepare_sampling!(sh_det_accel, tel, src)
    point_peak_accel = AdaptiveOpticsSim.sampled_spots_peak!(KA_CPU_STYLE, sh_det_accel, tel, src, det, MersenneTwister(21))
    @test point_peak_accel ≈ point_peak

    distinct_spectral = with_spectrum(src, SpectralBundle([
        SpectralSample(0.95 * wavelength(src), 0.5),
        SpectralSample(1.05 * wavelength(src), 0.5),
    ]))
    sh_distinct = ShackHartmannWFS(tel; n_lenslets=4,
        mode=Diffractive(), n_pix_subap=4)
    @test_throws InvalidConfiguration AdaptiveOpticsSim.sampled_spots_peak!(
        sh_distinct, tel, distinct_spectral, det, MersenneTwister(22))
    AdaptiveOpticsSim.prepare_sampling!(sh_distinct, tel, src)
    @test_throws InvalidConfiguration AdaptiveOpticsSim.sampled_spots_peak!(
        AdaptiveOpticsSim.ScalarCPUStyle(), sh_distinct, tel,
        distinct_spectral, det, MersenneTwister(22))
    @test_throws InvalidConfiguration AdaptiveOpticsSim.sampled_spots_peak!(
        KA_CPU_STYLE, sh_distinct, tel, distinct_spectral, det,
        MersenneTwister(22))

    common_spectral = with_spectrum(src, SpectralBundle([
        SpectralSample(wavelength(src), 0.5),
        SpectralSample(wavelength(src), 0.5),
    ]))
    sh_spectral = ShackHartmannWFS(tel; n_lenslets=4,
        mode=Diffractive(), n_pix_subap=4)
    spectral_peak = AdaptiveOpticsSim.sampled_spots_peak!(sh_spectral, tel,
        common_spectral, det, MersenneTwister(22))
    spectral_spots = copy(sh_spectral.acquisition.spot_cube)
    @test isfinite(spectral_peak)
    @test spectral_peak > 0
    @test spectral_peak ≈ point_peak
    @test spectral_spots ≈ point_spots
    sh_spectral_accel = ShackHartmannWFS(tel; n_lenslets=4,
        mode=Diffractive(), n_pix_subap=4)
    AdaptiveOpticsSim.prepare_sampling!(sh_spectral_accel, tel,
        common_spectral)
    spectral_peak_accel = AdaptiveOpticsSim.sampled_spots_peak!(
        KA_CPU_STYLE, sh_spectral_accel, tel, common_spectral, det,
        MersenneTwister(22))
    @test spectral_peak_accel ≈ spectral_peak
    @test sh_spectral_accel.acquisition.spot_cube ≈ spectral_spots

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

    lgs = LGSSource(elongation_factor=1.3, photon_irradiance=1.0)
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
    lgs_profile = LGSSource(elongation_factor=1.2, na_profile=na_profile,
        fwhm_spot_up=1.0, photon_irradiance=1.0)
    sh_lgs_profile_accel = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    AdaptiveOpticsSim.prepare_sampling!(sh_lgs_profile_accel, tel, lgs_profile)
    @test AdaptiveOpticsSim.sampled_spots_peak!(KA_CPU_STYLE, sh_lgs_profile_accel, tel, lgs_profile) >= 0
    sh_lgs_profile_det_accel = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    AdaptiveOpticsSim.prepare_sampling!(sh_lgs_profile_det_accel, tel, lgs_profile)
    @test AdaptiveOpticsSim.sampled_spots_peak!(
        KA_CPU_STYLE, sh_lgs_profile_det_accel, tel, lgs_profile, det, MersenneTwister(25)) >= 0
end

@testset "Shack-Hartmann LGS calibration cache identity" begin
    tel = Telescope(resolution=16, diameter=8.0,
        central_obstruction=0.0)
    profile = [80000.0 90000.0 100000.0; 0.2 0.6 0.2]
    shifted_profile = [80000.0 90000.0 100000.0; 0.1 0.6 0.3]
    base = LGSSource(na_profile=profile, laser_coordinates=(0.0, 0.0),
        fwhm_spot_up=0.8, elongation_factor=1.2,
        photon_irradiance=1.0)
    changed_elongation = LGSSource(na_profile=profile,
        laser_coordinates=(0.0, 0.0), fwhm_spot_up=0.8,
        elongation_factor=1.3, photon_irradiance=1.0)
    changed_launch = LGSSource(na_profile=profile,
        laser_coordinates=(1.5, -0.5), fwhm_spot_up=0.8,
        elongation_factor=1.2, photon_irradiance=1.0)
    changed_uplink_width = LGSSource(na_profile=profile,
        laser_coordinates=(0.0, 0.0), fwhm_spot_up=1.1,
        elongation_factor=1.2, photon_irradiance=1.0)
    changed_geometry = LGSSource(na_profile=profile,
        laser_coordinates=(1.5, -0.5), fwhm_spot_up=1.1,
        elongation_factor=1.4, photon_irradiance=1.0)
    changed_profile = LGSSource(na_profile=shifted_profile,
        laser_coordinates=(1.5, -0.5), fwhm_spot_up=1.1,
        elongation_factor=1.4, photon_irradiance=1.0)
    reshaped_profile = LGSSource(
        na_profile=[80000.0 90000.0 90000.0 100000.0;
            0.2 0.3 0.3 0.2],
        laser_coordinates=(0.0, 0.0), fwhm_spot_up=0.8,
        elongation_factor=1.2, photon_irradiance=1.0)

    base_signature = AdaptiveOpticsSim.calibration_signature(base)
    for changed_source in (changed_elongation, changed_launch,
        changed_uplink_width, changed_profile, reshaped_profile)
        @test AdaptiveOpticsSim.calibration_signature(changed_source) !=
            base_signature
    end
    altitude_base = LGSSource(altitude=90000.0, elongation_factor=1.2,
        photon_irradiance=1.0)
    altitude_changed = LGSSource(altitude=95000.0, elongation_factor=1.2,
        photon_irradiance=1.0)
    @test AdaptiveOpticsSim.calibration_signature(altitude_base) !=
        AdaptiveOpticsSim.calibration_signature(altitude_changed)

    signatures = map((base, changed_geometry, changed_profile)) do src
        AdaptiveOpticsSim.telescope_aperture_calibration_signature(tel,
            AdaptiveOpticsSim.calibration_signature(src))
    end

    wfs = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(),
        n_pix_subap=4)
    AdaptiveOpticsSim.prepare_sampling!(wfs, tel, base)
    for (src, signature) in zip(
        (base, changed_geometry, changed_profile), signatures)
        AdaptiveOpticsSim.ensure_sh_calibration!(wfs, tel, src)
        @test wfs.calibration.signature == signature
        @test wfs.calibration.calibrated
    end

    bioedge = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive())
    for (src, signature) in zip(
        (base, changed_geometry, changed_profile), signatures)
        AdaptiveOpticsSim.ensure_bioedge_calibration!(bioedge, tel, src)
        @test bioedge.estimator.state.calibration_signature == signature
        @test bioedge.estimator.state.calibrated
    end
end

@testset "Shack-Hartmann incoherent detector acquisition order" begin
    tel = Telescope(resolution=16, diameter=8.0, central_obstruction=0.0)
    src = Source(band=:I, magnitude=10.0)
    extended = with_extended_source(src,
        PointCloudSourceModel([(0.0, 0.0), (0.2, 0.0)], [0.5, 0.5]))
    detector_source = first(extended_source_asterism(extended).sources)

    for style in (ScalarCPUStyle(), KA_CPU_STYLE)
        expected_wfs = ShackHartmannWFS(tel;
            n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
        prepare_sampling!(expected_wfs, tel, src)
        sampled_spots_peak!(style, expected_wfs, tel, extended)
        expected_frame = copy(expected_wfs.acquisition.spot_cube)
        expected_detector = Detector(noise=NoiseReadout(0.25), qe=0.8,
            integration_time=0.01, binning=1)
        capture_stack!(expected_detector, expected_frame, similar(expected_frame),
            detector_source, MersenneTwister(777))
        zero_invalid_sh_spot_cube!(style, expected_frame,
            expected_wfs.front_end.layout.valid_mask)

        actual_wfs = ShackHartmannWFS(tel;
            n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
        prepare_sampling!(actual_wfs, tel, src)
        actual_detector = Detector(noise=NoiseReadout(0.25), qe=0.8,
            integration_time=0.01, binning=1)
        sampled_spots_peak!(style, actual_wfs, tel, extended,
            actual_detector, MersenneTwister(777))

        @test actual_wfs.acquisition.spot_cube ≈ expected_frame atol=1e-10 rtol=1e-12
    end

    saturation_wfs = ShackHartmannWFS(tel;
        n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    saturation_detector = Detector(noise=NoiseNone(), qe=1.0,
        integration_time=1.0, binning=1, full_well=100.0)
    sampled_spots_peak!(saturation_wfs, tel, extended,
        saturation_detector, MersenneTwister(778))
    @test maximum(saturation_wfs.acquisition.spot_cube) == 100.0

    qe_source = Source(band=:custom, wavelength=0.55e-6,
        photon_irradiance=1.0)
    sampled_qe = SampledQuantumEfficiency(
        [0.50e-6, 0.60e-6], [0.2, 0.8])
    for style in (ScalarCPUStyle(), KA_CPU_STYLE)
        sampled_qe_wfs = ShackHartmannWFS(tel;
            n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
        prepare_sampling!(sampled_qe_wfs, tel, qe_source)
        sampled_qe_detector = Detector(noise=NoiseNone(), qe=sampled_qe,
            integration_time=1.0, binning=1)
        sampled_spots_peak!(style, sampled_qe_wfs, tel, qe_source,
            sampled_qe_detector, MersenneTwister(779))

        scalar_qe_wfs = ShackHartmannWFS(tel;
            n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
        prepare_sampling!(scalar_qe_wfs, tel, qe_source)
        scalar_qe_detector = Detector(noise=NoiseNone(), qe=0.5,
            integration_time=1.0, binning=1)
        sampled_spots_peak!(style, scalar_qe_wfs, tel, qe_source,
            scalar_qe_detector, MersenneTwister(779))

        @test sampled_qe_wfs.acquisition.spot_cube ≈
            scalar_qe_wfs.acquisition.spot_cube atol=1e-10 rtol=1e-12
    end

    mixed_wavelength = Asterism([
        Source(band=:custom, wavelength=700e-9, photon_irradiance=1.0),
        Source(band=:custom, wavelength=800e-9, photon_irradiance=1.0),
    ])
    prepared_mixed_wfs = ShackHartmannWFS(tel;
        n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
    @test_throws InvalidConfiguration begin
        prepare_runtime_wfs!(prepared_mixed_wfs, tel, mixed_wavelength)
    end
    @test !prepared_mixed_wfs.calibration.calibrated
    prepared_mixed_pyramid = PyramidWFS(tel;
        pupil_samples=4, mode=Diffractive())
    @test_throws InvalidConfiguration begin
        prepare_runtime_wfs!(prepared_mixed_pyramid, tel, mixed_wavelength)
    end
    @test !prepared_mixed_pyramid.estimator.state.calibrated
    for style in (ScalarCPUStyle(), KA_CPU_STYLE)
        mixed_wfs = ShackHartmannWFS(tel;
            n_lenslets=4, mode=Diffractive(), n_pix_subap=4)
        prepare_sampling!(mixed_wfs, tel, first(mixed_wavelength.sources))
        @test_throws InvalidConfiguration begin
            AdaptiveOpticsSim.sampled_spots_peak_asterism_stacked!(
                style, mixed_wfs, tel, mixed_wavelength)
        end
    end
end

@testset "Shack-Hartmann diffractive photon-rate conservation" begin
    tel = Telescope(resolution=16, diameter=8.0,
        central_obstruction=0.0, pupil_reflectivity=0.25)
    src = Source(band=:custom, wavelength=0.75e-6,
        photon_irradiance=1.0)
    expected_rate = sum(pupil_photon_rate_map(tel, src))

    for style in (ScalarCPUStyle(), KA_CPU_STYLE), padding in (1, 2, 3)
        wfs = ShackHartmannWFS(tel; n_lenslets=1, mode=Diffractive(),
            n_pix_subap=16, diffraction_padding=padding)
        prepare_sampling!(wfs, tel, src)
        compute_intensity_stack!(style, wfs, tel, src)
        @test sum(wfs.front_end.propagation.intensity_stack) ≈ expected_rate atol=1e-10 rtol=1e-12
    end

    ast = Asterism([
        src,
        Source(band=:custom, wavelength=wavelength(src),
            photon_irradiance=2.0),
    ])
    wfs = ShackHartmannWFS(tel; n_lenslets=1, mode=Diffractive(),
        n_pix_subap=16, diffraction_padding=2)
    prepare_sampling!(wfs, tel, src)
    compute_intensity_asterism_stack!(KA_CPU_STYLE, wfs, tel, ast)
    @test sum(wfs.front_end.propagation.intensity_stack) ≈ 3 * expected_rate atol=1e-10 rtol=1e-12
end

@testset "Asterism direct imaging" begin
    tel = Telescope(resolution=16, diameter=8.0, central_obstruction=0.0)
    src1 = Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0))
    src2 = Source(band=:I, magnitude=0.0, coordinates=(1.0, 90.0))
    @test coordinates_xy_arcsec(src1) == (0.0, 0.0)
    @test coordinates_xy_arcsec(src2)[1] ≈ 0.0 atol=1e-12
    @test coordinates_xy_arcsec(src2)[2] ≈ 1.0
    ast = Asterism([src1, src2])
    prepared = prepare_direct_imaging(tel, PupilFunction(tel), ast;
        zero_padding=2)
    combined = form_direct_image!(prepared)
    components = map(direct_imaging_output,
        direct_imaging_components(prepared))
    @test length(components) == 2
    @test size(intensity_values(combined)) == (32, 32)
    @test sum(intensity_values(combined)) >=
        sum(intensity_values(components[1]))

    spectral = with_spectrum(src1,
        SpectralBundle([0.9 * wavelength(src1), 1.1 * wavelength(src1)],
            [0.5, 0.5]))
    extended = with_extended_source(src1,
        PointCloudSourceModel([(0.0, 0.0), (0.1, 0.0)], [0.5, 0.5]))
    for nested in (spectral, extended, ast)
        @test_throws UnsupportedAlgorithm Asterism([nested])
        heterogeneous = AdaptiveOpticsSim.AbstractSource[src1, nested]
        @test_throws UnsupportedAlgorithm Asterism(heterogeneous)
    end
end

@testset "Polychromatic WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    λ0 = wavelength(src)
    λshift = 1.05 * λ0
    bundle_single = SpectralBundle([SpectralSample(λ0, 1.0)])
    bundle_shifted_single = SpectralBundle([
        SpectralSample(λshift, 1.0),
    ])
    bundle_common = SpectralBundle([
        SpectralSample(λshift, 0.4),
        SpectralSample(λshift, 0.6),
    ])
    bundle_broad = SpectralBundle([SpectralSample(0.9 * λ0, 0.4), SpectralSample(1.1 * λ0, 0.6)])
    poly_single = with_spectrum(src, bundle_single)
    poly_shifted_single = with_spectrum(src, bundle_shifted_single)
    poly_common = with_spectrum(src, bundle_common)
    poly_broad = with_spectrum(src, bundle_broad)
    shifted_mono = source_with_wavelength_and_radiometric_value(src,
        λshift, photon_irradiance(src))

    @test sum(sample.weight for sample in bundle_broad) ≈ 1.0
    @test weighted_wavelength(bundle_broad) ≈ (0.9 * λ0 * 0.4 + 1.1 * λ0 * 0.6)
    @test has_spectral_bundle(poly_broad)
    @test is_polychromatic(poly_broad)
    @test is_polychromatic(poly_common)
    @test !is_polychromatic(poly_single)
    @test !is_polychromatic(poly_shifted_single)
    @test spectral_reference_source(poly_broad) === src

    zb = ZernikeBasis(tel, 5)
    compute_zernike!(zb, tel)
    focus = @view zb.modes[:, :, 5]
    @. tel.state.opd = 5e-8 * focus

    sh_mono = ShackHartmannWFS(tel; n_lenslets=8, mode=Diffractive())
    sh_single = ShackHartmannWFS(tel; n_lenslets=8, mode=Diffractive())
    sh_shifted_mono = ShackHartmannWFS(tel; n_lenslets=8,
        mode=Diffractive())
    sh_shifted_single = ShackHartmannWFS(tel; n_lenslets=8,
        mode=Diffractive())
    sh_common = ShackHartmannWFS(tel; n_lenslets=8, mode=Diffractive())
    sh_broad = ShackHartmannWFS(tel; n_lenslets=8, mode=Diffractive())

    mono_slopes = copy(measure!(sh_mono, tel, src))
    single_slopes = copy(measure!(sh_single, tel, poly_single))
    shifted_mono_slopes = copy(measure!(sh_shifted_mono, tel,
        shifted_mono))
    shifted_single_slopes = copy(measure!(sh_shifted_single, tel,
        poly_shifted_single))
    common_slopes_1 = copy(measure!(sh_common, tel, poly_common))
    common_slopes_2 = copy(measure!(sh_common, tel, poly_common))

    @test single_slopes ≈ mono_slopes atol=1e-10 rtol=1e-10
    @test shifted_single_slopes ≈ shifted_mono_slopes atol=1e-10 rtol=1e-10
    @test common_slopes_1 ≈ common_slopes_2 atol=1e-10 rtol=1e-10
    @test common_slopes_1 ≈ shifted_mono_slopes atol=1e-10 rtol=1e-10
    @test sh_common.acquisition.spot_cube ≈
        sh_shifted_mono.acquisition.spot_cube atol=1e-8 rtol=1e-12
    @test sum(sh_common.acquisition.spot_cube) ≈
        sum(sh_shifted_mono.acquisition.spot_cube) atol=1e-8 rtol=1e-12
    @test sh_shifted_single.calibration.wavelength == λshift
    @test sh_common.calibration.wavelength == λshift
    @test AdaptiveOpticsSim.sh_has_common_spectral_grid(sh_common,
        poly_common)
    @test !AdaptiveOpticsSim.sh_has_common_spectral_grid(sh_broad,
        poly_broad)
    @test sh_spectral_grid_guard_allocation_bytes(sh_common,
        poly_common) == 0
    @test supports_prepared_runtime(sh_common, poly_common)
    @test supports_stacked_sources(sh_common, poly_common)
    @test supports_grouped_execution(sh_common, poly_common)
    @test !supports_prepared_runtime(sh_broad, poly_broad)
    @test !supports_stacked_sources(sh_broad, poly_broad)
    @test !supports_grouped_execution(sh_broad, poly_broad)
    @test_throws InvalidConfiguration measure!(sh_broad, tel, poly_broad)

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
    spectral_frame = copy(measure!(sh_common, tel, poly_common, det))
    @test size(sh_common.acquisition.detector_noise_cube) ==
        size(sh_common.acquisition.spot_cube)
    @test spectral_frame ≈ common_slopes_1 atol=1e-10 rtol=1e-10
    @test sh_common.acquisition.spot_cube ≈
        sh_shifted_mono.acquisition.spot_cube atol=1e-8 rtol=1e-12

    sampled_qe = SampledQuantumEfficiency(
        [0.9 * λshift, 1.1 * λshift], [0.25, 0.75])
    effective_qe = qe_at(sampled_qe, λshift)
    exposure = 2.5
    for style in (ScalarCPUStyle(), KA_CPU_STYLE)
        spectral_wfs = ShackHartmannWFS(tel; n_lenslets=8,
            mode=Diffractive())
        spectral_detector = Detector(noise=NoiseNone(), qe=sampled_qe,
            integration_time=exposure, binning=1,
            response_model=NullFrameResponse())
        prepare_sampling!(spectral_wfs, tel, poly_common)
        sampled_spots_peak!(style, spectral_wfs, tel, poly_common,
            spectral_detector, MersenneTwister(793))

        monochromatic_wfs = ShackHartmannWFS(tel; n_lenslets=8,
            mode=Diffractive())
        monochromatic_detector = Detector(noise=NoiseNone(), qe=sampled_qe,
            integration_time=exposure, binning=1,
            response_model=NullFrameResponse())
        prepare_sampling!(monochromatic_wfs, tel, shifted_mono)
        sampled_spots_peak!(style, monochromatic_wfs, tel, shifted_mono,
            monochromatic_detector, MersenneTwister(793))

        optical_wfs = ShackHartmannWFS(tel; n_lenslets=8,
            mode=Diffractive())
        prepare_sampling!(optical_wfs, tel, shifted_mono)
        sampled_spots_peak!(style, optical_wfs, tel, shifted_mono)

        @test spectral_wfs.acquisition.spot_cube ≈
            monochromatic_wfs.acquisition.spot_cube atol=1e-8 rtol=1e-11
        @test spectral_wfs.acquisition.spot_cube ≈ (
            optical_wfs.acquisition.spot_cube .* (effective_qe * exposure)
        ) atol=1e-8 rtol=1e-11
    end

    mixed_wavelengths = [
        4.000015249516764e-7,
        4.0000149709840115e-7,
    ]
    @test Float32(mixed_wavelengths[1]) == Float32(mixed_wavelengths[2])
    mixed_source = with_spectrum(src, SpectralBundle(
        mixed_wavelengths, [0.5, 0.5]; T=Float64))
    mixed_scalar_wfs = ShackHartmannWFS(tel; n_lenslets=8,
        mode=Diffractive(), T=Float32)
    mixed_ka_wfs = ShackHartmannWFS(tel; n_lenslets=8,
        mode=Diffractive(), T=Float32)
    @test AdaptiveOpticsSim.sh_has_common_spectral_grid(mixed_scalar_wfs,
        mixed_source)
    prepare_sampling!(mixed_scalar_wfs, tel, mixed_source)
    prepare_sampling!(mixed_ka_wfs, tel, mixed_source)
    sampled_spots_peak!(ScalarCPUStyle(), mixed_scalar_wfs, tel,
        mixed_source)
    sampled_spots_peak!(KA_CPU_STYLE, mixed_ka_wfs, tel, mixed_source)
    @test mixed_ka_wfs.front_end.propagation.opd_to_cycles_host[1] ==
        mixed_ka_wfs.front_end.propagation.opd_to_cycles_host[2]
    @test mixed_ka_wfs.acquisition.spot_cube ≈
        mixed_scalar_wfs.acquisition.spot_cube atol=2e-5 rtol=2e-5

    guard_wfs = ShackHartmannWFS(tel; n_lenslets=8,
        mode=Diffractive())
    guard_detector = Detector(noise=NoisePhotonReadout(0.1), binning=1)
    measure!(guard_wfs, tel, poly_common, guard_detector;
        rng=MersenneTwister(791))
    guard_state_before = (
        slopes=copy(guard_wfs.estimator.slopes),
        intensity=copy(guard_wfs.front_end.propagation.intensity),
        spot_cube=copy(guard_wfs.acquisition.spot_cube),
        exported_spot_cube=copy(guard_wfs.acquisition.exported_spot_cube),
        reference_signal=copy(guard_wfs.calibration.reference_signal_2d),
        effective_padding=guard_wfs.front_end.propagation.effective_padding,
        binning_pixel_scale=guard_wfs.front_end.propagation.binning_pixel_scale,
        sampled_n_pix_subap=guard_wfs.front_end.propagation.sampled_n_pix_subap,
        phasor_ratio=guard_wfs.front_end.propagation.phasor_ratio,
        calibrated=guard_wfs.calibration.calibrated,
        calibration_wavelength=guard_wfs.calibration.wavelength,
        calibration_signature=guard_wfs.calibration.signature,
    )
    detector_frame_before = copy(output_frame(guard_detector))
    detector_integrated_time_before = guard_detector.state.integrated_time
    detector_readout_ready_before = guard_detector.state.readout_ready
    rejection_rng = MersenneTwister(792)
    rejection_rng_reference = copy(rejection_rng)
    @test_throws InvalidConfiguration measure!(guard_wfs, tel, poly_broad,
        guard_detector; rng=rejection_rng)
    @test isequal(guard_wfs.estimator.slopes, guard_state_before.slopes)
    @test isequal(guard_wfs.front_end.propagation.intensity, guard_state_before.intensity)
    @test isequal(guard_wfs.acquisition.spot_cube,
        guard_state_before.spot_cube)
    @test isequal(guard_wfs.acquisition.exported_spot_cube,
        guard_state_before.exported_spot_cube)
    @test isequal(guard_wfs.calibration.reference_signal_2d,
        guard_state_before.reference_signal)
    @test guard_wfs.front_end.propagation.effective_padding ==
        guard_state_before.effective_padding
    @test guard_wfs.front_end.propagation.binning_pixel_scale ==
        guard_state_before.binning_pixel_scale
    @test guard_wfs.front_end.propagation.sampled_n_pix_subap ==
        guard_state_before.sampled_n_pix_subap
    @test guard_wfs.front_end.propagation.phasor_ratio == guard_state_before.phasor_ratio
    @test guard_wfs.calibration.calibrated == guard_state_before.calibrated
    @test guard_wfs.calibration.wavelength ==
        guard_state_before.calibration_wavelength
    @test guard_wfs.calibration.signature ==
        guard_state_before.calibration_signature
    @test isequal(output_frame(guard_detector), detector_frame_before)
    @test guard_detector.state.integrated_time ==
        detector_integrated_time_before
    @test guard_detector.state.readout_ready == detector_readout_ready_before
    @test rand(rejection_rng, UInt64) ==
        rand(rejection_rng_reference, UInt64)

    pyr_det = Detector(noise=NoiseNone(), binning=1)
    pyr_det_slopes = measure!(pyr_broad, tel, poly_broad, pyr_det)
    @test size(output_frame(pyr_det)) == size(pyr_broad.acquisition.state.camera_frame)
    @test pyr_det_slopes ≈ broad_pyr_1 atol=1e-10 rtol=1e-10

    selective_bundle = SpectralBundle(
        [0.8 * λ0, 1.2 * λ0], [0.5, 0.5])
    selective_source = with_spectrum(src, selective_bundle)
    selective_qe = SampledQuantumEfficiency(
        [0.8 * λ0, 1.2 * λ0], [0.0, 1.0])
    transmitted_sample = selective_bundle[2]
    transmitted_source = source_with_wavelength_and_radiometric_value(
        selective_source, transmitted_sample.wavelength,
        photon_irradiance(selective_source) * transmitted_sample.weight)

    for style in (ScalarCPUStyle(), KA_CPU_STYLE)
        selective_sh = ShackHartmannWFS(tel;
            n_lenslets=8, mode=Diffractive())
        prepare_sampling!(selective_sh, tel, src)
        fill!(selective_sh.acquisition.spot_cube, 17.0)
        spot_cube_before = copy(selective_sh.acquisition.spot_cube)
        selective_detector = Detector(noise=NoiseNone(), qe=selective_qe,
            integration_time=1.0, binning=1,
            response_model=NullFrameResponse())
        selective_rng = MersenneTwister(780)
        selective_rng_reference = copy(selective_rng)
        @test_throws InvalidConfiguration sampled_spots_peak!(style,
            selective_sh, tel, selective_source, selective_detector,
            selective_rng)
        @test isequal(selective_sh.acquisition.spot_cube, spot_cube_before)
        @test rand(selective_rng, UInt64) ==
            rand(selective_rng_reference, UInt64)
    end

    selective_pyramid = PyramidWFS(tel;
        pupil_samples=8, mode=Diffractive(), modulation=1.0)
    selective_pyramid_detector = Detector(noise=NoiseNone(), qe=selective_qe,
        integration_time=1.0, binning=1,
        response_model=NullFrameResponse())
    measure!(selective_pyramid, tel, selective_source,
        selective_pyramid_detector; rng=MersenneTwister(781))
    selective_pyramid_frame = copy(output_frame(selective_pyramid_detector))

    transmitted_pyramid = PyramidWFS(tel;
        pupil_samples=8, mode=Diffractive(), modulation=1.0)
    transmitted_pyramid_detector = Detector(noise=NoiseNone(), qe=1.0,
        integration_time=1.0, binning=1,
        response_model=NullFrameResponse())
    measure!(transmitted_pyramid, tel, transmitted_source,
        transmitted_pyramid_detector; rng=MersenneTwister(781))
    @test selective_pyramid_frame ≈
        output_frame(transmitted_pyramid_detector) atol=1e-8 rtol=1e-11

    spectral_spad_pyramid = PyramidWFS(tel;
        pupil_samples=8, mode=Diffractive(), modulation=1.0)
    spectral_spad = SPADArrayDetector(integration_time=1.0,
        noise=NoiseNone(),
        sensor=SPADArraySensor(pde=1.0, dark_count_rate=0.0,
            fill_factor=1.0))
    spectral_spad_slopes = measure!(spectral_spad_pyramid, tel,
        selective_source, spectral_spad; rng=MersenneTwister(782))
    @test all(isfinite, spectral_spad_slopes)
    @test size(output_frame(spectral_spad)) ==
        size(spectral_spad_pyramid.acquisition.state.camera_frame)

    fill!(tel.state.opd, 0.0)
end

@testset "Extended-source model constructor validation" begin
    point = PointCloudSourceModel([(0.0, 0.0), (0.2, -0.1)], [0.25, 0.75])
    point_offsets = copy(point.offsets_xy_arcsec)
    point_weights = copy(point.weights)
    point_parameterized = PointCloudSourceModel{
        Float64,typeof(point_offsets),typeof(point_weights),
    }(point_offsets, point_weights)
    @test point_parameterized.offsets_xy_arcsec == point.offsets_xy_arcsec
    @test_throws InvalidConfiguration PointCloudSourceModel{
        Float64,Vector{NTuple{2,Float64}},Vector{Float64},
    }([(NaN, 0.0)], [1.0])
    @test_throws InvalidConfiguration PointCloudSourceModel{
        Float64,Vector{NTuple{2,Float64}},Vector{Float64},
    }([(0.0, 0.0)], [2.0])

    gaussian = GaussianDiskSourceModel(sigma_arcsec=0.35, n_side=5)
    gaussian_offsets = copy(gaussian.offsets_xy_arcsec)
    gaussian_weights = copy(gaussian.weights)
    gaussian_parameterized = GaussianDiskSourceModel{
        Float64,typeof(gaussian_offsets),typeof(gaussian_weights),
    }(gaussian.sigma_arcsec, gaussian.radius_sigma, gaussian.n_side,
        gaussian_offsets, gaussian_weights)
    @test gaussian_parameterized.weights == gaussian.weights
    invalid_gaussian_offsets = copy(gaussian_offsets)
    invalid_gaussian_offsets[1] =
        (invalid_gaussian_offsets[1][1] + 0.1,
         invalid_gaussian_offsets[1][2])
    @test_throws InvalidConfiguration GaussianDiskSourceModel{
        Float64,typeof(invalid_gaussian_offsets),typeof(gaussian_weights),
    }(gaussian.sigma_arcsec, gaussian.radius_sigma, gaussian.n_side,
        invalid_gaussian_offsets, gaussian_weights)

    sampled = SampledImageSourceModel(
        [0.0 1.0 0.0; 1.0 2.0 1.0; 0.0 1.0 0.0];
        pixel_scale_arcsec=0.2)
    sampled_image = copy(sampled.image)
    sampled_offsets = copy(sampled.offsets_xy_arcsec)
    sampled_weights = copy(sampled.weights)
    sampled_parameterized = SampledImageSourceModel{
        Float64,typeof(sampled_image),typeof(sampled_offsets),
        typeof(sampled_weights),
    }(sampled_image, sampled.pixel_scale_arcsec, sampled_offsets,
        sampled_weights)
    @test sampled_parameterized.offsets_xy_arcsec ==
        sampled.offsets_xy_arcsec
    invalid_sampled_image = copy(sampled_image)
    invalid_sampled_image[2, 2] += 1.0
    @test_throws InvalidConfiguration SampledImageSourceModel{
        Float64,typeof(invalid_sampled_image),typeof(sampled_offsets),
        typeof(sampled_weights),
    }(invalid_sampled_image, sampled.pixel_scale_arcsec, sampled_offsets,
        sampled_weights)

    src = Source(band=:custom, wavelength=1.0e-6,
        photon_irradiance=1.0)
    extended = with_extended_source(src, point)
    bundle = SpectralBundle([0.9e-6, 1.1e-6], [0.5, 0.5])
    @test_throws UnsupportedAlgorithm with_spectrum(extended, bundle)
end

@testset "Extended-source WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    point_model = PointCloudSourceModel([(0.0, 0.0)], [1.0])
    gaussian_model = GaussianDiskSourceModel(sigma_arcsec=0.35, n_side=5)
    image_model = SampledImageSourceModel([0.0 1.0 0.0; 1.0 2.0 1.0; 0.0 1.0 0.0], pixel_scale_arcsec=0.2)
    ext_point = with_extended_source(src, point_model)
    ext_gauss = with_extended_source(src, gaussian_model)
    ext_image = with_extended_source(src, image_model)
    normalized_src = Source(band=:custom, wavelength=wavelength(src),
        normalized_power=3.0)
    normalized_extended = with_extended_source(normalized_src,
        PointCloudSourceModel([(0.0, 0.0), (0.1, 0.0)], [0.25, 0.75]))
    normalized_ast = extended_source_asterism(normalized_extended)
    @test all(source_radiometry(sample) isa NormalizedTestSource
        for sample in normalized_ast.sources)
    @test sum(source_radiometric_value(sample)
        for sample in normalized_ast.sources) ≈ 3.0
    @test_throws InvalidConfiguration photon_irradiance(normalized_extended)

    point_ast = extended_source_asterism(ext_point)
    image_ast = extended_source_asterism(ext_image)
    @test has_extended_source_model(ext_gauss)
    @test !has_extended_source_model(src)
    @test length(point_ast) == 1
    @test length(image_ast) == 5
    @test extended_source_asterism(ext_image) === image_ast
    @test @allocated(extended_source_asterism(ext_image)) == 0
    @test AdaptiveOpticsSim.photon_irradiance(point_ast.sources[1]) ≈ AdaptiveOpticsSim.photon_irradiance(src)
    @test sum(AdaptiveOpticsSim.photon_irradiance(sample) for sample in image_ast.sources) ≈ AdaptiveOpticsSim.photon_irradiance(src)
    frozen_coords = coordinates_xy_arcsec(image_ast.sources[1])
    old_offset = image_model.offsets_xy_arcsec[1]
    image_model.offsets_xy_arcsec[1] = (old_offset[1] + 0.1, old_offset[2])
    frozen_ast = extended_source_asterism(ext_image)
    @test coordinates_xy_arcsec(frozen_ast.sources[1]) == frozen_coords

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
    point_spots = copy(sh_point.acquisition.spot_cube)
    ext_peak = AdaptiveOpticsSim.sampled_spots_peak!(sh_ext, tel, ext_gauss)
    extended_spots = copy(sh_ext.acquisition.spot_cube)
    ext_slopes_1 = copy(measure!(sh_ext, tel, ext_gauss))
    ext_slopes_2 = copy(measure!(sh_ext, tel, ext_gauss))

    @test ext_point_slopes ≈ point_slopes atol=1e-10 rtol=1e-10
    @test ext_slopes_1 ≈ ext_slopes_2 atol=1e-10 rtol=1e-10
    @test ext_peak <= point_peak * (1 + 1e-12)
    @test sum(extended_spots) ≈ sum(point_spots) rtol=1e-12
    sh_relative_morphology = norm(extended_spots - point_spots) /
                             norm(point_spots)
    # The quadrature currently preserves rates, but direct WFS propagation does
    # not yet apply each component's angular offset to the focal-plane shape.
    @test_broken sh_relative_morphology > 1e-6
    @test supports_stacked_sources(sh_ext, ext_gauss)

    pyr_point = PyramidWFS(tel; pupil_samples=8, mode=Diffractive(), modulation=1.0)
    pyr_ext_point = PyramidWFS(tel; pupil_samples=8, mode=Diffractive(), modulation=1.0)
    pyr_ext = PyramidWFS(tel; pupil_samples=8, mode=Diffractive(), modulation=1.0)
    pyr_point_slopes = copy(measure!(pyr_point, tel, src))
    pyramid_point_intensity = copy(pyr_point.front_end.propagation.intensity)
    pyr_ext_point_slopes = copy(measure!(pyr_ext_point, tel, ext_point))
    pyr_ext_slopes_1 = copy(measure!(pyr_ext, tel, ext_gauss))
    pyramid_extended_intensity = copy(pyr_ext.front_end.propagation.intensity)
    pyr_ext_slopes_2 = copy(measure!(pyr_ext, tel, ext_gauss))

    @test pyr_ext_point_slopes ≈ pyr_point_slopes atol=1e-10 rtol=1e-10
    @test pyr_ext_slopes_1 ≈ pyr_ext_slopes_2 atol=1e-10 rtol=1e-10
    @test sum(pyramid_extended_intensity) ≈ sum(pyramid_point_intensity) rtol=1e-12
    pyramid_relative_morphology = norm(pyramid_extended_intensity - pyramid_point_intensity) /
                                  norm(pyramid_point_intensity)
    @test_broken pyramid_relative_morphology > 1e-6
    @test supports_stacked_sources(pyr_ext, ext_gauss)

    det = Detector(noise=NoiseNone(), binning=1)
    sh_det_slopes = measure!(sh_ext, tel, ext_gauss, det)
    @test size(sh_ext.acquisition.detector_noise_cube) == size(sh_ext.acquisition.spot_cube)
    @test sh_det_slopes ≈ ext_slopes_1 atol=1e-10 rtol=1e-10

    pyr_det = Detector(noise=NoiseNone(), binning=1)
    pyr_det_slopes = measure!(pyr_ext, tel, ext_gauss, pyr_det)
    @test size(output_frame(pyr_det)) == size(pyr_ext.acquisition.state.camera_frame)
    @test pyr_det_slopes ≈ pyr_ext_slopes_1 atol=1e-10 rtol=1e-10

    fill!(tel.state.opd, 0.0)
end

@testset "Pupil masks and misregistration" begin
    tel = Telescope(resolution=16, diameter=8.0, central_obstruction=0.0)
    revision = AdaptiveOpticsSim.aperture_revision(tel)
    base_sum = sum(pupil_mask(tel))
    apply_spiders!(tel; thickness=0.5, angles=[0.0, 90.0])
    @test sum(pupil_mask(tel)) < base_sum
    @test AdaptiveOpticsSim.aperture_revision(tel) == revision + 1
    revision = AdaptiveOpticsSim.aperture_revision(tel)

    custom = trues(16, 16)
    custom[:, 9:end] .= false
    set_pupil!(tel, custom)
    @test sum(pupil_mask(tel)) == sum(custom)
    @test pupil_reflectivity(tel) == Float64.(custom)
    @test AdaptiveOpticsSim.aperture_revision(tel) == revision + 1
    custom[1, 1] = false
    @test pupil_mask(tel)[1, 1]
    revision = AdaptiveOpticsSim.aperture_revision(tel)

    reflectivity = fill(0.5, 16, 16)
    set_pupil_reflectivity!(tel, reflectivity)
    @test pupil_reflectivity(tel)[:, 1:8] == fill(0.5, 16, 8)
    @test pupil_reflectivity(tel)[:, 9:end] == fill(0.0, 16, 8)
    @test AdaptiveOpticsSim.aperture_revision(tel) == revision + 1
    reflectivity[1, 1] = 0.25
    @test pupil_reflectivity(tel)[1, 1] == 0.5

    tel2 = Telescope(resolution=16, diameter=8.0, central_obstruction=0.0)
    dm1 = DeformableMirror(tel2; n_act=2, influence_width=0.3)
    mis = Misregistration(shift_x=0.1, shift_y=0.0, rotation_deg=5.0, T=Float64)
    @test rotation_deg(mis) ≈ 5.0
    @test rotation_rad(mis) ≈ deg2rad(5.0)
    dm2 = DeformableMirror(tel2; n_act=2, influence_width=0.3, misregistration=mis)
    @test dm1.state.modes != dm2.state.modes
end

@testset "Telescope aperture calibration revisions" begin
    tel = Telescope(resolution=16, diameter=8.0,
        central_obstruction=0.0)
    src = Source(band=:custom, wavelength=0.75e-6,
        photon_irradiance=1.0)
    sh = ShackHartmannWFS(tel; n_lenslets=4, mode=Diffractive(),
        n_pix_subap=4)
    pyramid = PyramidWFS(tel; pupil_samples=4, mode=Diffractive(),
        modulation=1.0)
    bioedge = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive(),
        modulation=1.0)
    zernike = ZernikeWFS(tel; pupil_samples=4)
    curvature = CurvatureWFS(tel; pupil_samples=4)
    sensors = (sh, pyramid, bioedge, zernike, curvature)

    prepare_sampling!(sh, tel, src)
    ensure_sh_calibration!(sh, tel, src)
    ensure_pyramid_calibration!(pyramid, tel, src)
    ensure_bioedge_calibration!(bioedge, tel, src)
    ensure_zernike_calibration!(zernike, tel, src)
    ensure_curvature_calibration!(curvature, tel, src)
    initial_signature = telescope_aperture_calibration_signature(tel,
        calibration_signature(src))
    @test sh.calibration.signature == initial_signature
    @test all(wfs_calibration_signature(wfs) == initial_signature
        for wfs in sensors)

    for wfs in sensors
        fill!(reference_signal(wfs), NaN)
    end
    initial_revision = AdaptiveOpticsSim.aperture_revision(tel)
    set_pupil_reflectivity!(tel, 0.7)
    @test AdaptiveOpticsSim.aperture_revision(tel) == initial_revision + 1
    ensure_sh_calibration!(sh, tel, src)
    ensure_pyramid_calibration!(pyramid, tel, src)
    ensure_bioedge_calibration!(bioedge, tel, src)
    ensure_zernike_calibration!(zernike, tel, src)
    ensure_curvature_calibration!(curvature, tel, src)
    reflectivity_signature = telescope_aperture_calibration_signature(tel,
        calibration_signature(src))
    @test reflectivity_signature != initial_signature
    @test sh.calibration.signature == reflectivity_signature
    @test all(wfs_calibration_signature(wfs) == reflectivity_signature
        for wfs in sensors)
    @test all(all(isfinite, reference_signal(wfs)) for wfs in sensors)

    for wfs in sensors
        fill!(reference_signal(wfs), NaN)
    end
    half_pupil = trues(16, 16)
    half_pupil[:, 9:end] .= false
    set_pupil!(tel, half_pupil)
    ensure_sh_calibration!(sh, tel, src)
    ensure_pyramid_calibration!(pyramid, tel, src)
    ensure_bioedge_calibration!(bioedge, tel, src)
    ensure_zernike_calibration!(zernike, tel, src)
    ensure_curvature_calibration!(curvature, tel, src)
    pupil_signature = telescope_aperture_calibration_signature(tel,
        calibration_signature(src))
    @test pupil_signature != reflectivity_signature
    @test sh.calibration.signature == pupil_signature
    @test all(wfs_calibration_signature(wfs) == pupil_signature
        for wfs in sensors)
    @test all(all(isfinite, reference_signal(wfs)) for wfs in sensors)
end

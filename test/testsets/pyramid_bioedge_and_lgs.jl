@testset "Pyramid, BioEdge, and LGS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        tel.state.opd[i, j] = i
    end

    pyr = PyramidWFS(tel; pupil_samples=4, modulation=1.0)
    pyr_slopes = measure!(pyr, tel)
    @test length(pyr_slopes) == 2 * 4 * 4

    bio = BioEdgeWFS(tel; pupil_samples=4)
    bio_slopes = measure!(bio, tel)
    @test length(bio_slopes) == 2 * 4 * 4

    ngs = Source(band=:I, magnitude=0.0)

    pyr_direct = PyramidWFS(tel; pupil_samples=2, mode=Diffractive(), modulation=0.0)
    pyr_direct.state.nominal_detector_resolution = 4
    AdaptiveOpticsSim.resize_pyramid_signal_buffers!(pyr_direct, 4)
    pyr_direct.state.valid_i4q .= Bool[1 0; 1 1]
    AdaptiveOpticsSim.update_pyramid_valid_signal!(pyr_direct)
    @test AdaptiveOpticsSim.update_pyramid_valid_signal_indices!(pyr_direct) == 3
    AdaptiveOpticsSim.resize_pyramid_slope_buffers!(pyr_direct)
    fill!(pyr_direct.state.reference_signal_2d, 0.0)
    pyr_frame = [4.0 4.0 1.0 1.0;
                 4.0 4.0 1.0 1.0;
                 3.0 3.0 2.0 2.0;
                 3.0 3.0 2.0 2.0]
    pyr_direct_slopes = AdaptiveOpticsSim.pyramid_signal!(pyr_direct, tel, pyr_frame)
    @test length(pyr_direct_slopes) == 6
    @test pyr_direct_slopes[1:3] ≈ fill(0.4, 3)
    @test pyr_direct_slopes[4:6] ≈ zeros(3)

    pyr_invalid = PyramidWFS(tel; pupil_samples=2, mode=Diffractive(), modulation=0.0)
    pyr_invalid.state.nominal_detector_resolution = 4
    AdaptiveOpticsSim.resize_pyramid_signal_buffers!(pyr_invalid, 4)
    fill!(pyr_invalid.state.valid_i4q, false)
    AdaptiveOpticsSim.update_pyramid_valid_signal!(pyr_invalid)
    @test AdaptiveOpticsSim.update_pyramid_valid_signal_indices!(pyr_invalid) == 0
    @test_throws InvalidConfiguration AdaptiveOpticsSim.resize_pyramid_slope_buffers!(pyr_invalid)

    pyr_incidence = PyramidWFS(tel; pupil_samples=2, mode=Diffractive(), modulation=0.0,
        normalization=IncidenceFluxNormalization())
    expected_pyr_norm = AdaptiveOpticsSim.photon_flux(ngs) * tel.params.sampling_time *
                        (tel.params.diameter / pyr_incidence.params.pupil_samples)^2
    @test AdaptiveOpticsSim.pyramid_normalization(pyr_incidence.params.normalization,
        pyr_incidence, tel, ngs, 3, 10.0) ≈ expected_pyr_norm
    @test AdaptiveOpticsSim.pyramid_normalization(pyr_incidence.params.normalization,
        pyr_incidence, tel, nothing, 3, 10.0) == 1.0

    bio_direct = BioEdgeWFS(tel; pupil_samples=2, mode=Diffractive())
    bio_direct.state.nominal_detector_resolution = 4
    AdaptiveOpticsSim.resize_bioedge_signal_buffers!(bio_direct, 4)
    bio_direct.state.valid_i4q .= Bool[1 0; 1 1]
    AdaptiveOpticsSim.update_bioedge_valid_signal!(bio_direct)
    @test AdaptiveOpticsSim.update_bioedge_valid_signal_indices!(bio_direct) == 3
    AdaptiveOpticsSim.resize_bioedge_slope_buffers!(bio_direct)
    fill!(bio_direct.state.reference_signal_2d, 0.0)
    fill!(bio_direct.state.optical_gain, 2.0)
    bio_frame = copy(pyr_frame)
    bio_direct_slopes = AdaptiveOpticsSim.bioedge_signal!(bio_direct, tel, bio_frame)
    @test length(bio_direct_slopes) == 6
    @test bio_direct_slopes[1:3] ≈ fill(0.8, 3)
    @test bio_direct_slopes[4:6] ≈ zeros(3)

    bio_invalid = BioEdgeWFS(tel; pupil_samples=2, mode=Diffractive())
    bio_invalid.state.nominal_detector_resolution = 4
    AdaptiveOpticsSim.resize_bioedge_signal_buffers!(bio_invalid, 4)
    fill!(bio_invalid.state.valid_i4q, false)
    AdaptiveOpticsSim.update_bioedge_valid_signal!(bio_invalid)
    @test AdaptiveOpticsSim.update_bioedge_valid_signal_indices!(bio_invalid) == 0
    @test_throws InvalidConfiguration AdaptiveOpticsSim.resize_bioedge_slope_buffers!(bio_invalid)

    bio_incidence = BioEdgeWFS(tel; pupil_samples=2, mode=Diffractive(),
        normalization=IncidenceFluxNormalization())
    expected_bio_norm = AdaptiveOpticsSim.photon_flux(ngs) * tel.params.sampling_time *
                        (tel.params.diameter / bio_incidence.params.pupil_samples)^2
    @test AdaptiveOpticsSim.bioedge_normalization(bio_incidence.params.normalization,
        bio_incidence, tel, ngs, 3, 10.0) ≈ expected_bio_norm
    @test AdaptiveOpticsSim.bioedge_normalization(bio_incidence.params.normalization,
        bio_incidence, tel, nothing, 3, 10.0) == 1.0

    sh = ShackHartmann(tel; n_lenslets=4)
    lgs = LGSSource(elongation_factor=2.0)
    slopes_ngs = measure!(sh, tel, ngs)
    slopes_lgs = measure!(sh, tel, lgs)
    n = sh.params.n_lenslets * sh.params.n_lenslets
    @test slopes_lgs[n+1:end] ≈ slopes_ngs[n+1:end] .* 2.0
end

@testset "Diffractive WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        tel.state.opd[i, j] = i
    end
    ngs = Source(band=:I, magnitude=0.0)
    lgs = LGSSource(elongation_factor=1.5)

    sh = ShackHartmann(tel; n_lenslets=4, mode=Diffractive())
    @test_throws InvalidConfiguration measure!(sh, tel)
    sh_slopes = measure!(sh, tel, ngs)
    @test length(sh_slopes) == 2 * 4 * 4
    @test all(isfinite, sh_slopes)
    sh_lgs = measure!(sh, tel, lgs)
    @test all(isfinite, sh_lgs)

    na_profile = [80000.0 90000.0 100000.0; 0.2 0.6 0.2]
    lgs_profile = LGSSource(elongation_factor=1.2, na_profile=na_profile, fwhm_spot_up=1.0)
    sh_profile = ShackHartmann(tel; n_lenslets=4, mode=Diffractive())
    sh_profile_slopes = measure!(sh_profile, tel, lgs_profile)
    @test all(isfinite, sh_profile_slopes)

    sh_sampled = ShackHartmann(tel; n_lenslets=4, mode=Diffractive(), pixel_scale=0.06, n_pix_subap=8)
    sh_sampled_slopes = measure!(sh_sampled, tel, ngs)
    @test length(sh_sampled_slopes) == 2 * 4 * 4

    pyr_sampled = PyramidWFS(tel; pupil_samples=4, mode=Diffractive(), n_pix_separation=4, binning=2)
    pyr_sampled_slopes = measure!(pyr_sampled, tel, ngs)
    @test length(pyr_sampled_slopes) == 2 * count(pyr_sampled.state.valid_i4q)
    pyr_intensity = reshape(Float64.(1:size(tel.state.opd, 1)^2), size(tel.state.opd))
    pyr_frame = copy(AdaptiveOpticsSim.sample_pyramid_intensity!(pyr_sampled, tel, pyr_intensity))
    pyr_camera = zeros(Float64, 4, 4)
    pyr_manual = zeros(Float64, 2, 2)
    AdaptiveOpticsSim.bin2d!(pyr_camera, pyr_intensity, 8)
    AdaptiveOpticsSim.bin2d!(pyr_manual, pyr_camera, 2)
    @test pyr_frame == pyr_manual

    bio_sampled = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive(), binning=2)
    bio_sampled_slopes = measure!(bio_sampled, tel, ngs)
    @test length(bio_sampled_slopes) == 2 * count(bio_sampled.state.valid_i4q)
    bio_intensity = reshape(Float64.(1:size(tel.state.opd, 1)^2), size(tel.state.opd))
    bio_frame = copy(AdaptiveOpticsSim.sample_bioedge_intensity!(bio_sampled, tel, bio_intensity))
    bio_camera = zeros(Float64, 4, 4)
    bio_manual = similar(bio_frame)
    AdaptiveOpticsSim.bin2d!(bio_camera, bio_intensity, 8)
    AdaptiveOpticsSim.bin2d!(bio_manual, bio_camera, div(size(bio_camera, 1), size(bio_frame, 1)))
    @test bio_frame == bio_manual

    pyr_profile = PyramidWFS(tel; pupil_samples=4, mode=Diffractive())
    pyr_profile_slopes = measure!(pyr_profile, tel, lgs_profile)
    @test all(isfinite, pyr_profile_slopes)

    bio_profile = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive())
    bio_profile_slopes = measure!(bio_profile, tel, lgs_profile)
    @test all(isfinite, bio_profile_slopes)

    pyr = PyramidWFS(tel; pupil_samples=4, mode=Diffractive())
    @test_throws InvalidConfiguration measure!(pyr, tel)
    pyr_slopes = measure!(pyr, tel, ngs)
    @test length(pyr_slopes) == 2 * 4 * 4

    bio = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive())
    @test_throws InvalidConfiguration measure!(bio, tel)
    bio_slopes = measure!(bio, tel, ngs)
    @test length(bio_slopes) == 2 * 4 * 4

    det = Detector(noise=NoiseNone(), binning=1)
    sh_det = ShackHartmann(tel; n_lenslets=4, mode=Diffractive())
    sh_det_slopes = measure!(sh_det, tel, ngs, det)
    @test length(sh_det_slopes) == 2 * 4 * 4
    sh_det_image = wfs_detector_image(sh_det, det; gap=1)
    @test ndims(sh_det_image) == 2
    sh_adu_det = Detector(noise=NoiseNone(), binning=1, full_well=30_000.0, bits=12, output_type=UInt16)
    sh_adu = ShackHartmann(tel; n_lenslets=4, mode=Diffractive())
    measure!(sh_adu, tel, ngs, sh_adu_det; rng=MersenneTwister(15))
    sh_adu_image = wfs_detector_image(sh_adu, sh_adu_det; gap=1)
    @test sh_adu_image isa Matrix{UInt16}
    @test maximum(sh_adu_image) <= 0x0fff
    pyr_det = PyramidWFS(tel; pupil_samples=4, mode=Diffractive())
    pyr_det_slopes = measure!(pyr_det, tel, ngs, det)
    @test length(pyr_det_slopes) == 2 * 4 * 4
    @test wfs_detector_image(pyr_det, det) === output_frame(det)
    bio_det = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive())
    bio_det_slopes = measure!(bio_det, tel, ngs, det)
    @test length(bio_det_slopes) == 2 * 4 * 4
    @test wfs_detector_image(bio_det, det) === output_frame(det)

    ast = Asterism([ngs, Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0))])
    sh_ast = ShackHartmann(tel; n_lenslets=4, mode=Diffractive())
    sh_ast_slopes = copy(measure!(sh_ast, tel, ast))
    @test length(sh_ast_slopes) == 2 * 4 * 4
    sh_ast_serial = ShackHartmann(tel; n_lenslets=4, mode=Diffractive())
    AdaptiveOpticsSim.prepare_sampling!(sh_ast_serial, tel, ast.sources[1])
    AdaptiveOpticsSim.ensure_sh_calibration!(sh_ast_serial, tel, ast.sources[1])
    fill!(sh_ast_serial.state.detector_noise_cube, zero(eltype(sh_ast_serial.state.detector_noise_cube)))
    for src in ast.sources
        AdaptiveOpticsSim.sampled_spots_peak!(sh_ast_serial, tel, src)
        sh_ast_serial.state.detector_noise_cube .+= sh_ast_serial.state.spot_cube
    end
    copyto!(sh_ast_serial.state.spot_cube, sh_ast_serial.state.detector_noise_cube)
    sh_ast_serial_peak = maximum(sh_ast_serial.state.spot_cube)
    AdaptiveOpticsSim.sh_signal_from_spots!(sh_ast_serial, sh_ast_serial_peak, slope_extraction_model(sh_ast_serial))
    AdaptiveOpticsSim.subtract_reference_and_scale!(sh_ast_serial)
    sh_ast_serial_slopes = copy(sh_ast_serial.state.slopes)
    @test norm(sh_ast_slopes - sh_ast_serial_slopes) / norm(sh_ast_slopes) < 0.07
    mixed_ngs = Source(wavelength=wavelength(lgs_profile), magnitude=0.0, coordinates=(0.0, 0.0))
    mixed_ast = Asterism([mixed_ngs, lgs_profile])
    sh_mixed_det = ShackHartmann(tel; n_lenslets=4, mode=Diffractive())
    sh_mixed_det_slopes = copy(measure!(sh_mixed_det, tel, mixed_ast, det; rng=MersenneTwister(14)))
    sh_mixed_det_frame = copy(sh_mixed_det.state.spot_cube)
    sh_mixed_det_manual = ShackHartmann(tel; n_lenslets=4, mode=Diffractive())
    AdaptiveOpticsSim.prepare_sampling!(sh_mixed_det_manual, tel, mixed_ast.sources[1])
    AdaptiveOpticsSim.ensure_sh_calibration!(sh_mixed_det_manual, tel, mixed_ast.sources[1])
    fill!(sh_mixed_det_manual.state.detector_noise_cube, zero(eltype(sh_mixed_det_manual.state.detector_noise_cube)))
    for src in mixed_ast.sources
        AdaptiveOpticsSim.sampled_spots_peak!(sh_mixed_det_manual, tel, src, det, MersenneTwister(14))
        sh_mixed_det_manual.state.detector_noise_cube .+= sh_mixed_det_manual.state.spot_cube
    end
    copyto!(sh_mixed_det_manual.state.spot_cube, sh_mixed_det_manual.state.detector_noise_cube)
    @test sh_mixed_det_frame ≈ sh_mixed_det_manual.state.spot_cube
    @test length(sh_mixed_det_slopes) == 2 * 4 * 4
    pyr_ast = PyramidWFS(tel; pupil_samples=4, mode=Diffractive())
    pyr_ast_slopes = copy(measure!(pyr_ast, tel, ast))
    @test length(pyr_ast_slopes) == 2 * 4 * 4
    pyr_ast_serial = PyramidWFS(tel; pupil_samples=4, mode=Diffractive())
    AdaptiveOpticsSim.ensure_pyramid_calibration!(pyr_ast_serial, tel, ast.sources[1])
    pyr_ast_stack = @view AdaptiveOpticsSim.ensure_pyramid_asterism_stack!(pyr_ast_serial, length(ast.sources))[:, :, 1:length(ast.sources)]
    fill!(pyr_ast_serial.state.intensity, zero(eltype(pyr_ast_serial.state.intensity)))
    for (src_idx, src) in pairs(ast.sources)
        AdaptiveOpticsSim.pyramid_intensity!(@view(pyr_ast_stack[:, :, src_idx]), pyr_ast_serial, tel, src)
        pyr_ast_serial.state.intensity .+= @view(pyr_ast_stack[:, :, src_idx])
    end
    pyr_ast_intensity = AdaptiveOpticsSim.sample_pyramid_intensity!(pyr_ast_serial, tel, pyr_ast_serial.state.intensity)
    AdaptiveOpticsSim.pyramid_signal!(pyr_ast_serial, tel, pyr_ast_intensity)
    @. pyr_ast_serial.state.slopes *= pyr_ast_serial.state.optical_gain
    @test pyr_ast_slopes ≈ pyr_ast_serial.state.slopes
    bio_ast = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive())
    bio_ast_slopes = copy(measure!(bio_ast, tel, ast))
    @test length(bio_ast_slopes) == 2 * 4 * 4
    bio_ast_serial = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive())
    AdaptiveOpticsSim.ensure_bioedge_calibration!(bio_ast_serial, tel, ast.sources[1])
    fill!(bio_ast_serial.state.binned_intensity, zero(eltype(bio_ast_serial.state.binned_intensity)))
    for src in ast.sources
        AdaptiveOpticsSim.bioedge_intensity!(bio_ast_serial.state.intensity, bio_ast_serial, tel, src)
        bio_ast_serial.state.binned_intensity .+= bio_ast_serial.state.intensity
    end
    bio_ast_intensity = AdaptiveOpticsSim.sample_bioedge_intensity!(bio_ast_serial, tel, bio_ast_serial.state.binned_intensity)
    AdaptiveOpticsSim.bioedge_signal!(bio_ast_serial, tel, bio_ast_intensity)
    @test bio_ast_slopes ≈ bio_ast_serial.state.slopes

    pyr_ast_det = PyramidWFS(tel; pupil_samples=4, mode=Diffractive())
    pyr_ast_det_slopes = copy(measure!(pyr_ast_det, tel, ast, det))
    pyr_ast_det_serial = PyramidWFS(tel; pupil_samples=4, mode=Diffractive())
    AdaptiveOpticsSim.ensure_pyramid_calibration!(pyr_ast_det_serial, tel, ast.sources[1])
    pyr_ast_det_stack = @view AdaptiveOpticsSim.ensure_pyramid_asterism_stack!(pyr_ast_det_serial, length(ast.sources))[:, :, 1:length(ast.sources)]
    fill!(pyr_ast_det_serial.state.intensity, zero(eltype(pyr_ast_det_serial.state.intensity)))
    for (src_idx, src) in pairs(ast.sources)
        AdaptiveOpticsSim.pyramid_intensity!(@view(pyr_ast_det_stack[:, :, src_idx]), pyr_ast_det_serial, tel, src)
        pyr_ast_det_serial.state.intensity .+= @view(pyr_ast_det_stack[:, :, src_idx])
    end
    pyr_ast_det_intensity = AdaptiveOpticsSim.sample_pyramid_intensity!(pyr_ast_det_serial, tel, pyr_ast_det_serial.state.intensity)
    pyr_ast_det_frame = capture!(det, pyr_ast_det_intensity; rng=MersenneTwister(12))
    AdaptiveOpticsSim.resize_pyramid_signal_buffers!(pyr_ast_det_serial, size(pyr_ast_det_frame, 1))
    AdaptiveOpticsSim.pyramid_signal!(pyr_ast_det_serial, tel, pyr_ast_det_frame)
    @. pyr_ast_det_serial.state.slopes *= pyr_ast_det_serial.state.optical_gain
    @test pyr_ast_det_slopes ≈ pyr_ast_det_serial.state.slopes

    bio_ast_det = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive())
    bio_ast_det_slopes = copy(measure!(bio_ast_det, tel, ast, det; rng=MersenneTwister(13)))
    bio_ast_det_serial = BioEdgeWFS(tel; pupil_samples=4, mode=Diffractive())
    AdaptiveOpticsSim.ensure_bioedge_calibration!(bio_ast_det_serial, tel, ast.sources[1])
    fill!(bio_ast_det_serial.state.binned_intensity, zero(eltype(bio_ast_det_serial.state.binned_intensity)))
    for src in ast.sources
        AdaptiveOpticsSim.bioedge_intensity!(bio_ast_det_serial.state.intensity, bio_ast_det_serial, tel, src)
        bio_ast_det_serial.state.binned_intensity .+= bio_ast_det_serial.state.intensity
    end
    bio_ast_det_intensity = AdaptiveOpticsSim.sample_bioedge_intensity!(bio_ast_det_serial, tel, bio_ast_det_serial.state.binned_intensity)
    bio_ast_det_frame = capture!(det, bio_ast_det_intensity; rng=MersenneTwister(13))
    AdaptiveOpticsSim.resize_bioedge_signal_buffers!(bio_ast_det_serial, size(bio_ast_det_frame, 1))
    AdaptiveOpticsSim.bioedge_signal!(bio_ast_det_serial, tel, bio_ast_det_frame)
    @test bio_ast_det_slopes ≈ bio_ast_det_serial.state.slopes
end

@testset "Shack-Hartmann subapertures" begin
    tel = Telescope(resolution=24, diameter=8.0, sampling_time=1e-3, central_obstruction=0.1)
    src = Source(band=:I, magnitude=0.0)
    sh = ShackHartmann(tel; n_lenslets=6, mode=Diffractive(), pixel_scale=0.06, n_pix_subap=8, threshold_cog=0.02)

    layout = subaperture_layout(sh)
    calibration = subaperture_calibration(sh)
    @test layout isa SubapertureLayout
    @test calibration isa SubapertureCalibration
    @test layout.n_subap == 6
    @test layout.subap_pixels == 4
    @test layout.pitch_m ≈ tel.params.diameter / 6
    @test !calibration.calibrated
    @test slope_extraction_model(sh) isa CenterOfGravityExtraction
    @test slope_extraction_model(sh).threshold ≈ 0.02
    @test n_valid_subapertures(layout) == count(layout.valid_mask_host)

    prepare_runtime_wfs!(sh, tel, src)
    @test calibration.calibrated
    @test calibration.slopes_units == sh.state.slopes_units
    @test calibration.wavelength == sh.state.calibration_wavelength
    @test calibration.signature == sh.state.calibration_signature
    @test calibration.reference_signal_2d === sh.state.reference_signal_2d
    @test calibration.reference_signal_host === sh.state.reference_signal_host
    @test length(valid_subaperture_indices(layout)) == n_valid_subapertures(layout)

    slopes = measure!(sh, tel, src)
    @test all(isfinite, slopes)
    meta = AdaptiveOpticsSim.wfs_output_metadata(sh)
    @test meta.n_valid_subap == n_valid_subapertures(layout)
    @test meta.subap_pixels == layout.subap_pixels
    @test meta.calibrated

    dm = DeformableMirror(tel; n_act=5)
    imat = interaction_matrix(dm, sh, tel, src; amplitude=1e-8)
    @test size(imat.matrix, 1) == length(sh.state.slopes)
    @test size(imat.matrix, 2) == length(dm.state.coefs)
end

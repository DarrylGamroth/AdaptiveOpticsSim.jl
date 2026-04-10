@testset "Shack-Hartmann valid subaperture policies" begin
    tel = Telescope(resolution=352, diameter=1.22, sampling_time=1 / 500)
    sh_geom = ShackHartmann(tel; n_subap=16, mode=Diffractive(), T=Float32)
    sh_flux = ShackHartmann(tel;
        n_subap=16,
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

    sh_mono = ShackHartmann(tel; n_subap=8, mode=Diffractive())
    sh_single = ShackHartmann(tel; n_subap=8, mode=Diffractive())
    sh_broad = ShackHartmann(tel; n_subap=8, mode=Diffractive())

    mono_slopes = copy(measure!(sh_mono, tel, src))
    single_slopes = copy(measure!(sh_single, tel, poly_single))
    broad_slopes_1 = copy(measure!(sh_broad, tel, poly_broad))
    broad_slopes_2 = copy(measure!(sh_broad, tel, poly_broad))

    @test single_slopes ≈ mono_slopes atol=1e-10 rtol=1e-10
    @test broad_slopes_1 ≈ broad_slopes_2 atol=1e-10 rtol=1e-10
    @test norm(broad_slopes_1 - mono_slopes) > 1e-8
    @test supports_stacked_sources(sh_broad, poly_broad)

    pyr_mono = PyramidWFS(tel; n_subap=8, mode=Diffractive(), modulation=1.0)
    pyr_single = PyramidWFS(tel; n_subap=8, mode=Diffractive(), modulation=1.0)
    pyr_broad = PyramidWFS(tel; n_subap=8, mode=Diffractive(), modulation=1.0)

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

    sh_point = ShackHartmann(tel; n_subap=8, mode=Diffractive())
    sh_ext_point = ShackHartmann(tel; n_subap=8, mode=Diffractive())
    sh_ext = ShackHartmann(tel; n_subap=8, mode=Diffractive())
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

    pyr_point = PyramidWFS(tel; n_subap=8, mode=Diffractive(), modulation=1.0)
    pyr_ext_point = PyramidWFS(tel; n_subap=8, mode=Diffractive(), modulation=1.0)
    pyr_ext = PyramidWFS(tel; n_subap=8, mode=Diffractive(), modulation=1.0)
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

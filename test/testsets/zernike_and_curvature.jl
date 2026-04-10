@testset "Zernike WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    wfs = ZernikeWFS(tel; n_subap=8, diffraction_padding=2)

    @test size(wfs.state.camera_frame) == (8, 8)
    @test length(wfs.state.slopes) == count(wfs.state.valid_mask)
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

@testset "Curvature WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    wfs = CurvatureWFS(tel; n_subap=8, defocus_rms_nm=500.0)

    @test size(wfs.state.camera_frame) == (16, 8)
    @test length(wfs.state.slopes) == 64
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

    counting = CurvatureWFS(tel; n_subap=8, defocus_rms_nm=500.0, readout_model=CurvatureCountingReadout())
    counting_flat = copy(measure!(counting, tel, src))
    @test size(counting.state.camera_frame) == (2, 64)
    @test counting_flat ≈ zero.(counting_flat) atol=1e-10
    @test_throws InvalidConfiguration measure!(counting, tel, src, det)
    apd = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0, noise=NoiseNone())
    counting_apd = copy(measure!(counting, tel, src, apd))
    @test counting_apd ≈ counting_flat atol=1e-10
    @test detector_export_metadata(apd).readout.output_size == size(counting.state.camera_frame)
    apd_dead = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), dead_time_model=NonParalyzableDeadTime(0.25))
    counting_dead = copy(measure!(counting, tel, src, apd_dead))
    @test counting_dead ≈ counting_flat atol=1e-10
    @test_throws InvalidConfiguration CurvatureWFS(tel; n_subap=8, readout_model=CurvatureCountingReadout(),
        readout_pixels_per_subap=2)

    response = CurvatureBranchResponse(T=Float64, plus_throughput=1.2, minus_throughput=0.8,
        plus_background=5.0, minus_background=1.0)
    imbalanced = CurvatureWFS(tel; n_subap=8, defocus_rms_nm=500.0, branch_response=response)
    imbalanced_flat = copy(measure!(imbalanced, tel, src))
    @test imbalanced_flat ≈ zero.(imbalanced_flat) atol=1e-10
    plus_mean = mean(@view imbalanced.state.camera_frame[1:imbalanced.params.n_subap, :])
    minus_mean = mean(@view imbalanced.state.camera_frame[imbalanced.params.n_subap+1:end, :])
    @test plus_mean > minus_mean
    @test_throws InvalidConfiguration CurvatureBranchResponse(plus_throughput=-1.0)

    oversampled = CurvatureWFS(tel; n_subap=8, readout_crop_resolution=16, readout_pixels_per_subap=2)
    oversampled_flat = copy(measure!(oversampled, tel, src))
    @test size(oversampled.state.camera_frame) == (32, 16)
    @test size(oversampled.state.frame_plus) == (16, 16)
    @test size(oversampled.state.reduced_plus) == (8, 8)
    @test oversampled_flat ≈ zero.(oversampled_flat) atol=1e-10
    @test_throws InvalidConfiguration CurvatureWFS(tel; n_subap=8, readout_crop_resolution=18, readout_pixels_per_subap=2)

    atm = MultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[0.7, 0.3],
        wind_speed=[8.0, 4.0],
        wind_direction=[0.0, 90.0],
        altitude=[0.0, 5000.0],
    )
    advance!(atm, tel; rng=MersenneTwister(3))
    atm_slopes = copy(measure!(wfs, tel, src, atm))
    @test all(isfinite, atm_slopes)
    @test norm(atm_slopes) > 0

    det_atm = Detector(noise=NoiseNone(), binning=1)
    det_atm_slopes = copy(measure!(wfs, tel, src, atm, det_atm))
    @test all(isfinite, det_atm_slopes)

    ast = Asterism([src, Source(band=:I, magnitude=0.0, coordinates=(1.0, 90.0))])
    ast_slopes = copy(measure!(wfs, tel, ast, atm))
    @test length(ast_slopes) == length(wfs.state.slopes)
    @test all(isfinite, ast_slopes)
    @test norm(ast_slopes) > 0
end

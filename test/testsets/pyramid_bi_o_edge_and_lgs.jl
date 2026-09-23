@testset "LGS convolution normalization" begin
    n = 8
    expected = reshape(collect(range(0.25, 2.0; length=n * n)), n, n)
    identity_kernel_fft = ones(ComplexF64, n, n)
    fft_buffer = zeros(ComplexF64, n, n)
    fft_plan = AdaptiveOpticsSim.Backends.plan_fft_backend!(fft_buffer)
    ifft_plan = AdaptiveOpticsSim.Backends.plan_ifft_backend!(fft_buffer)
    actual = copy(expected)
    WavefrontSensors.apply_lgs_convolution!(actual, identity_kernel_fft,
        fft_buffer, fft_plan, ifft_plan)
    @test actual ≈ expected rtol=1e-12 atol=1e-12
    @test sum(actual) ≈ sum(expected) rtol=1e-12

    expected_stack = cat(expected, reverse(expected; dims=2); dims=3)
    stack_fft = zeros(ComplexF64, n, n, 2)
    stack_plan = AdaptiveOpticsSim.Backends.plan_fft_backend!(stack_fft, (1, 2))
    stack_ifft = AdaptiveOpticsSim.Backends.plan_ifft_backend!(stack_fft, (1, 2))
    actual_stack = copy(expected_stack)
    WavefrontSensors.apply_lgs_convolution_stack!(actual_stack,
        ones(ComplexF64, n, n, 2), stack_fft, stack_plan, stack_ifft)
    @test actual_stack ≈ expected_stack rtol=1e-12 atol=1e-12
end

@testset "Pyramid physical four-pupil propagation" begin
    T = Float64
    tel = Telescope(resolution=16, diameter=T(8), central_obstruction=zero(T), T=T)
    pupil = PupilFunction(tel; T=T)
    pupil.opd .= reshape(T.(1:256), 16, 16) .* T(1e-10)
    src = Source(band=:custom, wavelength=T(0.75e-6), photon_irradiance=T(10), T=T)

    pyramid = PyramidWFS(tel; pupil_samples=4, modulation=T(1), T=T)
    @test @inferred(pyramid_focal_mask(pyramid)) ===
        pyramid_propagation_workspace(pyramid).pyramid_mask
    @test @allocated(pyramid_focal_mask(pyramid)) == 0
    @test !hasfield(typeof(pyramid), :estimator)
    @test !applicable(measure!, pyramid, pupil, src)
    @test !applicable(slopes, pyramid)
    @test !supports_valid_subaperture_mask(pyramid)
    @test valid_subaperture_mask(pyramid) === nothing
    @test !supports_reference_signal(pyramid)
    @test reference_signal(pyramid) === nothing

    front_end = PyramidOpticalFrontEnd(pyramid, src)
    rate = pyramid_rate_map(front_end, pupil)
    plan = prepare_wfs_optics(front_end, pupil, rate)
    @test @inferred(form_wfs_optical_products!(rate, pupil, plan)) === rate
    @test all(isfinite, rate.values)
    @test sum(rate.values) > zero(T)
    @test @allocated(form_wfs_optical_products!(rate, pupil, plan)) == 0

    intensity = pyramid_propagation_workspace(pyramid).intensity
    frame = WavefrontSensors.sample_pyramid_intensity!(pyramid, pupil, intensity)
    @test frame === pyramid_acquisition_products(pyramid).frame
    @test size(frame) == (8, 8)
    detector = Detector(noise=NoiseNone(), exposure_duration=T(0.25), qe=T(0.5), T=T)
    captured = capture!(detector, frame, src; rng=MersenneTwister(19))
    @test captured === output_frame(detector)
    @test size(captured) == size(frame)
    @test captured ≈ frame .* T(0.125)

    shifted = PyramidWFS(tel; pupil_samples=4, modulation=T(1), T=T,
        modulation_propagation_strategy=PyramidShiftedMaskStrategy())
    shifted_front_end = PyramidOpticalFrontEnd(shifted, src)
    shifted_rate = pyramid_rate_map(shifted_front_end, pupil)
    shifted_plan = prepare_wfs_optics(shifted_front_end, pupil, shifted_rate)
    form_wfs_optical_products!(shifted_rate, pupil, shifted_plan)
    @test all(isfinite, shifted_rate.values)
    @test sum(shifted_rate.values) > zero(T)

    unit_mask = PyramidWFS(tel; pupil_samples=4, modulation=zero(T), T=T, mask_scale=one(T))
    scaled_mask = PyramidWFS(tel; pupil_samples=4, modulation=zero(T), T=T, mask_scale=T(1.5))
    @test pyramid_propagation_workspace(unit_mask).pyramid_mask !=
          pyramid_propagation_workspace(scaled_mask).pyramid_mask
    face_shifted_mask = PyramidWFS(tel; pupil_samples=4,
        modulation=zero(T), T=T,
        pupil_shift_x_pixels=(1, -1, -1, 1),
        pupil_shift_y_pixels=(1, 1, -1, -1))
    @test pyramid_propagation_workspace(unit_mask).pyramid_mask !=
          pyramid_propagation_workspace(face_shifted_mask).pyramid_mask
    uniformly_shifted_mask = PyramidWFS(tel; pupil_samples=4,
        modulation=zero(T), T=T, pupil_shift_x_pixels=1,
        pupil_shift_y_pixels=-1)
    @test uniformly_shifted_mask.front_end.phase_mask.pupil_shift_x_pixels ==
          (1, 1, 1, 1)
    @test uniformly_shifted_mask.front_end.phase_mask.pupil_shift_y_pixels ==
          (-1, -1, -1, -1)
    old_unit_mask = PyramidWFS(tel; pupil_samples=4, modulation=zero(T), T=T, old_mask=true, mask_scale=one(T))
    old_scaled_mask = PyramidWFS(tel; pupil_samples=4, modulation=zero(T), T=T, old_mask=true, mask_scale=T(1.5))
    @test pyramid_propagation_workspace(old_unit_mask).pyramid_mask !=
          pyramid_propagation_workspace(old_scaled_mask).pyramid_mask
    for invalid_scale in (0.0, -1.0, NaN, Inf)
        @test_throws InvalidConfiguration PyramidWFS(tel; pupil_samples=4, mask_scale=invalid_scale)
    end
    for invalid_shift in ((1, 2, 3), (1, 2, 3, Inf), NaN, "one")
        @test_throws InvalidConfiguration PyramidWFS(tel; pupil_samples=4,
            pupil_shift_x_pixels=invalid_shift)
    end
    @test_throws InvalidConfiguration PyramidWFS(tel; pupil_samples=4, n_pix_separation=-1)
    @test_throws InvalidConfiguration PyramidWFS(tel; pupil_samples=4, n_pix_separation=2, n_pix_edge=-1)

    lgs = LGSSource(wavelength=T(589e-9), photon_irradiance=T(4),
        sodium_layer_profile=SodiumLayerProfile(T[80_000, 90_000, 100_000], T[0.2, 0.6, 0.2]),
        laser_launch_xy_m=(T(1), T(-0.5)), fwhm_spot_up=T(0.8), T=T)
    @test WavefrontSensors.ensure_lgs_kernel!(pyramid, pupil, lgs) === pyramid
    lgs_tag = pyramid_propagation_workspace(pyramid).lgs_kernel_tag
    lgs_front_end = PyramidOpticalFrontEnd(pyramid, lgs)
    lgs_rate = pyramid_rate_map(lgs_front_end, pupil)
    lgs_plan = prepare_wfs_optics(lgs_front_end, pupil, lgs_rate)
    form_wfs_optical_products!(lgs_rate, pupil, lgs_plan)
    @test all(isfinite, lgs_rate.values)
    lgs.params.sodium_layer_profile.relative_weights .= T[0.8, 0.1, 0.1]
    WavefrontSensors.ensure_lgs_kernel!(pyramid, pupil, lgs)
    @test pyramid_propagation_workspace(pyramid).lgs_kernel_tag != lgs_tag

    spectral = with_spectrum(src, SpectralBundle(T[0.70e-6, 0.80e-6], T[0.25, 0.75]; T=T))
    spectral_front_end = PyramidOpticalFrontEnd(pyramid, spectral)
    spectral_rates = pyramid_rate_map(spectral_front_end, pupil)
    spectral_plan = prepare_wfs_optics(spectral_front_end, pupil, spectral_rates)
    form_wfs_optical_products!(spectral_rates, pupil, spectral_plan)
    @test all(product -> all(isfinite, product.values), spectral_rates)

    path_source = Asterism([src, Source(band=:custom, wavelength=wavelength(src), photon_irradiance=T(3), T=T)])
    second_pupil = PupilFunction(tel; T=T)
    copyto!(second_pupil.opd, pupil.opd)
    path_front_end = PyramidOpticalFrontEnd(pyramid, path_source)
    path_rates = pyramid_rate_map(path_front_end, (pupil, second_pupil))
    path_plan = prepare_wfs_optics(path_front_end, (pupil, second_pupil), path_rates)
    form_wfs_optical_products!(path_rates, (pupil, second_pupil), path_plan)
    @test sum(path_rates[1].values) / sum(path_rates[2].values) ≈ T(10 / 3)
end

@testset "Bi-O-edge physical four-pupil propagation" begin
    T = Float64
    tel = Telescope(resolution=16, diameter=T(8),
        central_obstruction=zero(T), T=T)
    pupil = PupilFunction(tel; T=T)
    pupil.opd .= reshape(T.(1:256), 16, 16) .* T(1e-10)
    src = Source(band=:custom, wavelength=T(0.75e-6),
        photon_irradiance=T(10), T=T)
    bio = BiOEdgeWFS(tel; pupil_samples=2, modulation=zero(T), T=T)

    @test !hasfield(typeof(bio), :estimator)
    @test !applicable(measure!, bio, pupil, src)
    @test !applicable(slopes, bio)
    @test !supports_valid_subaperture_mask(bio)
    @test valid_subaperture_mask(bio) === nothing
    @test !supports_reference_signal(bio)
    @test reference_signal(bio) === nothing

    front_end = BiOEdgeOpticalFrontEnd(bio, src)
    rate = bi_o_edge_rate_map(front_end, pupil)
    plan = prepare_wfs_optics(front_end, pupil, rate)
    @test @inferred(form_wfs_optical_products!(rate, pupil, plan)) === rate
    @test all(isfinite, rate.values)
    @test sum(rate.values) > zero(T)
    @test @allocated(form_wfs_optical_products!(rate, pupil, plan)) == 0

    detector = Detector(noise=NoiseNone(), exposure_duration=T(0.25),
        qe=T(0.5), T=T)
    captured = capture!(detector, rate.values, src; rng=MersenneTwister(19))
    @test captured === output_frame(detector)
    @test size(captured) == size(rate.values)
    @test captured ≈ rate.values .* T(0.125)

    spectral = with_spectrum(src,
        SpectralBundle(T[0.70e-6, 0.80e-6], T[0.25, 0.75]; T=T))
    spectral_front_end = BiOEdgeOpticalFrontEnd(bio, spectral)
    spectral_rates = bi_o_edge_rate_map(spectral_front_end, pupil)
    spectral_plan = prepare_wfs_optics(
        spectral_front_end, pupil, spectral_rates)
    form_wfs_optical_products!(spectral_rates, pupil, spectral_plan)
    @test all(product -> all(isfinite, product.values), spectral_rates)

    path_source = Asterism([src, Source(band=:custom,
        wavelength=wavelength(src), photon_irradiance=T(3), T=T)])
    second_pupil = PupilFunction(tel; T=T)
    copyto!(second_pupil.opd, pupil.opd)
    path_front_end = BiOEdgeOpticalFrontEnd(bio, path_source)
    path_rates = bi_o_edge_rate_map(
        path_front_end, (pupil, second_pupil))
    path_plan = prepare_wfs_optics(
        path_front_end, (pupil, second_pupil), path_rates)
    form_wfs_optical_products!(
        path_rates, (pupil, second_pupil), path_plan)
    @test sum(path_rates[1].values) / sum(path_rates[2].values) ≈ T(10 / 3)
end

@testset "Four-pupil pupil-reflectivity throughput" begin
    transmission = 0.25
    src = Source(band=:custom, wavelength=0.75e-6, photon_irradiance=1.0)
    full_tel = Telescope(resolution=16, diameter=8.0, central_obstruction=0.0)
    attenuated_tel = Telescope(resolution=16, diameter=8.0, central_obstruction=0.0, pupil_reflectivity=transmission)
    full_pupil = PupilFunction(full_tel)
    attenuated_pupil = PupilFunction(attenuated_tel)

    full_pyramid = PyramidWFS(full_tel; pupil_samples=2, modulation=0.0)
    attenuated_pyramid = PyramidWFS(attenuated_tel; pupil_samples=2, modulation=0.0)
    full_intensity = pyramid_propagation_workspace(full_pyramid).intensity
    attenuated_intensity = pyramid_propagation_workspace(attenuated_pyramid).intensity
    pyramid_intensity!(full_intensity, full_pyramid, full_pupil, src)
    pyramid_intensity!(attenuated_intensity, attenuated_pyramid, attenuated_pupil, src)
    @test sum(full_intensity) > 0
    @test sum(attenuated_intensity) ≈ transmission * sum(full_intensity)

    full_bio = BiOEdgeWFS(full_tel; pupil_samples=2, modulation=0.0)
    attenuated_bio = BiOEdgeWFS(attenuated_tel; pupil_samples=2, modulation=0.0)
    bi_o_edge_intensity!(full_bio.front_end.propagation.workspace.intensity, full_bio, full_pupil, src)
    bi_o_edge_intensity!(attenuated_bio.front_end.propagation.workspace.intensity, attenuated_bio, attenuated_pupil, src)
    @test sum(attenuated_bio.front_end.propagation.workspace.intensity) ≈
          transmission * sum(full_bio.front_end.propagation.workspace.intensity)
end

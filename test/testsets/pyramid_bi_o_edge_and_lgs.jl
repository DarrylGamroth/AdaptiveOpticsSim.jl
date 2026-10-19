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

# Independent block integration oracle for the detector's native pixel scale.
function pyramid_shift_reference_sampling(intensity, factor)
    side = div(size(intensity, 1), factor)
    return [sum(@view intensity[(i-1)*factor+1:i*factor,
        (j-1)*factor+1:j*factor]) for i in 1:side, j in 1:side]
end

@testset "Pyramid fractional pupil shifts" begin
    for T in (Float32, Float64)
        @testset "$T" begin
            tol = T === Float32 ? 3e-5 : 2e-12
            tel = Telescope(resolution=16, diameter=T(8),
                central_obstruction=zero(T), T=T)
            pupil = PupilFunction(tel; T=T)
            pupil.opd .= reshape(T.(1:256), 16, 16) .* T(1e-10)
            src = Source(band=:custom, wavelength=T(0.75e-6),
                photon_irradiance=T(10), T=T)
            baseline = PyramidWFS(tel; pupil_samples=4, modulation=0, T=T)
            zero_shift = PyramidWFS(tel; pupil_samples=4, modulation=0, T=T,
                pupil_shift_x_pixels=(0, 0, 0, 0), pupil_shift_y_pixels=0.0)
            @test baseline.front_end.phase_mask.pupil_shift_x_pixels ===
                ntuple(_ -> zero(T), 4)
            @test baseline.front_end.phase_mask.pupil_shift_y_pixels ===
                ntuple(_ -> zero(T), 4)
            @test pyramid_focal_mask(baseline) == pyramid_focal_mask(zero_shift)
            for old_mask in (false, true)
                integer = PyramidWFS(tel; pupil_samples=4, modulation=0, T=T,
                    old_mask=old_mask, pupil_shift_x_pixels=(1, -1, -1, 1),
                    pupil_shift_y_pixels=(1, 1, -1, -1))
                real_integer = PyramidWFS(tel; pupil_samples=4, modulation=0, T=T,
                    old_mask=old_mask, pupil_shift_x_pixels=T.((1, -1, -1, 1)),
                    pupil_shift_y_pixels=T.((1, 1, -1, -1)))
                @test pyramid_focal_mask(integer) == pyramid_focal_mask(real_integer)
            end
            for invalid in (NaN, Inf, -Inf, (0, 0, NaN, 0),
                    (0, 0, 0, Inf), "quarter")
                @test_throws InvalidConfiguration PyramidWFS(tel;
                    pupil_samples=4, T=T, pupil_shift_x_pixels=invalid)
                @test_throws InvalidConfiguration PyramidWFS(tel;
                    pupil_samples=4, T=T, pupil_shift_y_pixels=invalid)
            end
            if T === Float32
                @test_throws InvalidConfiguration PyramidWFS(tel;
                    pupil_samples=4, T=T, pupil_shift_x_pixels=floatmax(Float64))
                @test_throws InvalidConfiguration PyramidWFS(tel;
                    pupil_samples=4, T=T,
                    pupil_shift_y_pixels=(0, 0, floatmax(Float64), 0))
            end

            for (dx, dy) in ((T(0.25), T(-0.375)), (T(-0.25), T(0.375)))
                shifted = PyramidWFS(tel; pupil_samples=4, modulation=0, T=T,
                    pupil_shift_x_pixels=dx, pupil_shift_y_pixels=dy)
                @test shifted.front_end.phase_mask.pupil_shift_x_pixels ===
                    ntuple(_ -> dx, 4)
                @test shifted.front_end.phase_mask.pupil_shift_y_pixels ===
                    ntuple(_ -> dy, 4)
                # Uniform face shifts add a common linear phase to every plane.
                # The Fourier shift theorem gives displacement (4dx, -4dy)
                # in propagated array samples, or (dx, -dy) native pixels.
                mask = pyramid_focal_mask(baseline)
                side = size(mask, 1)
                coords = [T(2pi) * (T(i-1) / side - T(0.5) + inv(T(2side)))
                    for i in 1:side]
                ramp = [cis(-T(4) * (x * dx - y * dy))
                    for x in coords, y in coords]
                @test pyramid_focal_mask(shifted) ≈ mask .* ramp rtol=tol atol=tol
                field = zeros(Complex{T}, side, side)
                offset = div(side - 16, 2)
                field[offset+1:offset+16, offset+1:offset+16] .=
                    sqrt(T(10) * T(0.5)^2) .* pupil.amplitude .*
                    cispi.(T(2) .* pupil.opd ./ wavelength(src))
                # Independent FFT propagation of the analytic phase ramp.
                phasor = [cis(-T(pi) * T(side+1) / T(side) * (i+j-2))
                    for i in 1:side, j in 1:side]
                expected = abs2.(AdaptiveOpticsSim.FFTW.ifft(
                    AdaptiveOpticsSim.FFTW.fft(field .* phasor) .* mask .* ramp))
                intensity = pyramid_propagation_workspace(shifted).intensity
                pyramid_intensity!(intensity, shifted, pupil, src)
                @test intensity ≈ expected rtol=tol atol=tol
                frame = WavefrontSensors.sample_pyramid_intensity!(shifted, pupil,
                    intensity)
                @test frame ≈ pyramid_shift_reference_sampling(expected, 4) rtol=tol atol=tol
                base_intensity = pyramid_propagation_workspace(baseline).intensity
                pyramid_intensity!(base_intensity, baseline, pupil, src)
                base_frame = WavefrontSensors.sample_pyramid_intensity!(baseline,
                    pupil, base_intensity)
                @test norm(frame - base_frame) > T(0.01) * norm(base_frame)
                @test sum(intensity) ≈ sum(base_intensity) rtol=tol
                # Final binning integrates the same displaced diffraction field.
                binned = PyramidWFS(tel; pupil_samples=4, modulation=0, T=T,
                    binning=2, pupil_shift_x_pixels=dx, pupil_shift_y_pixels=dy)
                front_end = PyramidOpticalFrontEnd(binned, src)
                rate = pyramid_rate_map(front_end, pupil)
                plan = prepare_wfs_optics(front_end, pupil, rate)
                form_wfs_optical_products!(rate, pupil, plan)
                @test rate.values ≈ pyramid_shift_reference_sampling(expected, 8) rtol=tol atol=tol
                @test rate.values ≈ pyramid_shift_reference_sampling(frame, 2) rtol=tol atol=tol
                @test WavefrontSensors.sample_pyramid_intensity!(binned, pupil,
                    pyramid_propagation_workspace(binned).intensity) ≈ rate.values
                detector = Detector(noise=NoiseNone(), exposure_duration=T(0.25),
                    qe=T(0.5), T=T)
                @test capture!(detector, rate.values, src) ≈ rate.values .* T(0.125)
                @test @allocated(form_wfs_optical_products!(rate, pupil, plan)) == 0
            end

            paired_x = T.((0.25, -0.375, -0.375, 0.25))
            paired_y = T.((-0.125, -0.125, 0.375, 0.375))
            for (sx, sy, paired) in ((paired_x, paired_y, true),
                    (T.((0.25, -0.375, 0.125, -0.25)),
                        T.((-0.125, 0.25, 0.375, -0.25)), false))
                shifted = PyramidWFS(tel; pupil_samples=4, modulation=T(1), T=T,
                    modulation_propagation_strategy=PyramidShiftedMaskStrategy(),
                    pupil_shift_x_pixels=sx, pupil_shift_y_pixels=sy)
                front_end = PyramidOpticalFrontEnd(shifted, src)
                rate = pyramid_rate_map(front_end, pupil)
                plan = prepare_wfs_optics(front_end, pupil, rate)
                form_wfs_optical_products!(rate, pupil, plan)
                workspace = pyramid_propagation_workspace(shifted)
                batch = workspace.modulation_batch
                @test (batch isa PyramidSeparableShiftedMaskModulationWorkspace) == paired
                side = size(workspace.intensity, 1)
                points = length(batch.operating_weights)
                # Build the four physical planes independently at each shifted
                # focal coordinate, then take their upper envelope.
                masks = Array{Complex{T}}(undef, side, side, points)
                for point in 1:points, j in 1:side, i in 1:side
                    x = T(2pi) * (T(i-1) / side - T(0.5) + inv(T(2side))) +
                        batch.axis_1_shifts_rad[point]
                    y = T(2pi) * (T(j-1) / side - T(0.5) + inv(T(2side))) +
                        batch.axis_2_shifts_rad[point]
                    planes = ((2+sx[1])*x + (2-sy[1])*y,
                        (-2+sx[2])*x + (2-sy[2])*y,
                        (-2+sx[3])*x + (-2-sy[3])*y,
                        (2+sx[4])*x + (-2-sy[4])*y)
                    masks[i, j, point] = cis(-T(4) * maximum(planes))
                end
                if paired
                    for point in 1:points
                        @test batch.axis_1_factors[:, point] *
                            transpose(batch.axis_2_factors[:, point]) ≈
                            masks[:, :, point] rtol=tol atol=tol
                    end
                else
                    @test batch.shifted_masks ≈ masks rtol=tol atol=tol
                end
                full_batch = PyramidShiftedMaskModulationWorkspace(
                    similar(batch.field_stack), masks, batch.operating_weights,
                    batch.axis_1_shifts_rad, batch.axis_2_shifts_rad,
                    batch.bfft_plan, batch.batch_size)
                expected = zero(workspace.intensity)
                WavefrontSensors._propagate_pyramid_shifted_masks!(ScalarCPUStyle(),
                    expected, front_end, full_batch)
                @test workspace.intensity ≈ expected rtol=tol atol=tol
                @test rate.values ≈ pyramid_shift_reference_sampling(expected, 4) rtol=tol atol=tol
                @test @allocated(form_wfs_optical_products!(rate, pupil, plan)) == 0
            end
        end
    end
end

@testset "Fractional Pyramid mask in AOC gain sensing" begin
    tel = Telescope(resolution=16, diameter=8.0, central_obstruction=0.0)
    pupil = PupilFunction(tel)
    src = Source(band=:custom, wavelength=0.75e-6, photon_irradiance=10.0)
    sensor = PyramidWFS(tel; pupil_samples=4, modulation=1.0,
        pupil_shift_x_pixels=(0.25, -0.375, 0.125, -0.25),
        pupil_shift_y_pixels=(-0.125, 0.25, 0.375, -0.25))
    mask = pyramid_focal_mask(sensor)
    frame = pyramid_modulation_frame(sensor, pupil, src)
    side = size(mask, 1)
    coords = collect(range(-1.0, 1.0; length=side))
    basis = cat([x for x in coords, y in coords],
        [x^2 + y^2 for x in coords, y in coords]; dims=3)
    optical_gains = AdaptiveOpticsCalibration.OpticalGains
    specification = optical_gains.GainSensingSpecification(mask, basis, frame)
    plan = AdaptiveOpticsCalibration.prepare(optical_gains.GainSensing(), specification)
    @test plan.mask == mask
    @test plan.mask !== mask
    product = AdaptiveOpticsCalibration.allocate_result(plan)
    workspace = AdaptiveOpticsCalibration.allocate_workspace(plan)
    AdaptiveOpticsCalibration.process!(product, workspace, plan,
        optical_gains.GainSensingInputs(frame))
    @test !any(optical_gains.weak_mode_mask(product))
    @test all(isfinite, optical_gains.current_sensitivities(product))
    @test optical_gains.optical_gains(product) ≈ ones(2) rtol=1e-12 atol=1e-12
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

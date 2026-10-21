@testset "Shack-Hartmann physical subaperture geometry" begin
    telescope = Telescope(resolution=352, diameter=1.22)
    geometric = ShackHartmannWFS(telescope; n_lenslets=16, T=Float32)
    flux = ShackHartmannWFS(telescope;
        n_lenslets=16,
        valid_subaperture_policy=
            RelativeIlluminationValidSubapertures(peak_fraction=0.5f0, T=Float32),
        T=Float32)

    geometric_mask = copy(geometric.front_end.layout.valid_mask_host)
    flux_mask = copy(flux.front_end.layout.valid_mask_host)
    @test count(geometric_mask) == 216
    @test count(flux_mask) == 208
    @test flux_mask != geometric_mask
    @test isempty(findall(flux_mask .& .!geometric_mask))

    sensor = ShackHartmannWFS(
        Telescope(resolution=4, diameter=2.0);
        n_lenslets=2,
        n_pix_subap=2,
        T=Float32,
    )
    layout = sensor.front_end.layout
    revision = subaperture_layout_revision(layout)
    @test_throws DimensionMismatchError set_valid_subapertures!(sensor,
        fill(true, 1, 1))
    parent = fill(false, 4, 4)
    parent[2:3, 2:3] .= Bool[true false; false true]
    @test @inferred(set_valid_subapertures!(sensor, @view(parent[2:3, 2:3]))) === sensor
    @test subaperture_layout_revision(layout) == revision + UInt(1)
    @test n_valid_subapertures(layout) == 2
    @test valid_subaperture_indices(layout) == CartesianIndex{2}[
        CartesianIndex(1, 1), CartesianIndex(2, 2)]
end

@testset "Shack-Hartmann physical optics and acquisition" begin
    T = Float64
    telescope = Telescope(resolution=16, diameter=T(8),
        central_obstruction=zero(T), T=T)
    pupil = PupilFunction(telescope; T=T)
    source = Source(band=:custom, wavelength=T(0.75e-6),
        photon_irradiance=T(10), T=T)
    sensor = ShackHartmannWFS(telescope; n_lenslets=4,
        n_pix_subap=4, T=T)
    optics = shack_hartmann_optics(sensor, source)
    rate = shack_hartmann_rate_map(optics, pupil)
    prepared = @inferred prepare_wfs_optics(optics, pupil, rate)
    @test prepared isa WavefrontSensors.PreparedShackHartmannOptics
    @test @inferred(form_wfs_optical_products!(rate, pupil, prepared)) === rate
    @test all(isfinite, rate.values)
    @test sum(rate.values) > zero(T)
    @test rate.metadata.normalization isa PhotonRateNormalization
    @test rate.metadata.spatial_measure isa CellIntegratedMeasure
    @test rate.metadata.coherence isa IncoherentIntensityAddition

    detector = Detector(noise=NoiseNone(), exposure_duration=T(0.25),
        qe=one(T), response_model=NullFrameResponse(), T=T)
    observation = WFSObservation(similar(rate.values);
        units=:electron_count, layout=:lenslet_mosaic)
    acquisition = @inferred prepare_wfs_acquisition(detector, rate,
        observation)
    @test @inferred(acquire_wfs_observation!(observation, rate,
        acquisition, Xoshiro(0x309))) === observation
    @test observation.storage ≈ rate.values .* T(0.25) atol=0 rtol=0
    @test_throws WFSPreparationError prepare_wfs_estimation(sensor,
        observation,
        WFSMeasurement(zeros(T, 2); units=:pixel, kind=:centroid_slopes))

    if !coverage_instrumented()
        @test @allocated(form_wfs_optical_products!(rate, pupil,
            prepared)) == 0
        rng = Xoshiro(0x309)
        acquire_wfs_observation!(observation, rate, acquisition, rng)
        @test @allocated(acquire_wfs_observation!(observation, rate,
            acquisition, rng)) == 0
    end
end

@testset "Shack-Hartmann physical source composition" begin
    T = Float64
    telescope = Telescope(resolution=16, diameter=T(8),
        central_obstruction=zero(T), T=T)
    pupil = PupilFunction(telescope; T=T)
    source = Source(band=:custom, wavelength=T(0.75e-6),
        photon_irradiance=T(6), T=T)

    asterism = Asterism([
        source,
        Source(band=:custom, wavelength=wavelength(source),
            photon_irradiance=T(4), separation_arcsec=T(0.1), position_angle_deg=T(-0.05), T=T),
    ])
    asterism_sensor = ShackHartmannWFS(telescope; n_lenslets=4,
        n_pix_subap=4, T=T)
    asterism_rate = shack_hartmann_rate_map(asterism_sensor, pupil,
        asterism)
    asterism_plan = prepare_wfs_optics(
        shack_hartmann_optics(asterism_sensor, asterism), pupil,
        asterism_rate)
    form_wfs_optical_products!(asterism_rate, pupil, asterism_plan)
    @test all(isfinite, asterism_rate.values)
    @test supports_stacked_sources(asterism_sensor, asterism)
    @test supports_grouped_execution(asterism_sensor, asterism)

    lgs = LGSSource(wavelength=wavelength(source), photon_irradiance=T(6),
        elongation_factor=T(1.8), T=T)
    lgs_sensor = ShackHartmannWFS(telescope; n_lenslets=4,
        n_pix_subap=4, T=T)
    lgs_rate = shack_hartmann_rate_map(lgs_sensor, pupil, lgs)
    lgs_plan = prepare_wfs_optics(shack_hartmann_optics(lgs_sensor, lgs),
        pupil, lgs_rate)
    form_wfs_optical_products!(lgs_rate, pupil, lgs_plan)
    @test all(isfinite, lgs_rate.values)
    @test sum(lgs_rate.values) > zero(T)

    spectral = with_spectrum(source, SpectralBundle(
        T[0.9 * wavelength(source), 1.1 * wavelength(source)],
        T[0.4, 0.6]; T=T))
    spectral_sensor = ShackHartmannWFS(telescope; n_lenslets=4,
        n_pix_subap=4, T=T)
    spectral_rates = shack_hartmann_rate_map(spectral_sensor, pupil,
        spectral)
    @test spectral_rates isa OpticalProductBundle
    @test length(spectral_rates) == 2
    spectral_plan = prepare_wfs_optics(
        shack_hartmann_optics(spectral_sensor, spectral), pupil,
        spectral_rates)
    form_wfs_optical_products!(spectral_rates, pupil, spectral_plan)
    @test all(all(isfinite, product.values) for product in spectral_rates)

    extended = with_extended_source(source,
        PointCloudSourceModel([(0.0, 0.0), (0.2, 0.0)], [0.5, 0.5]))
    extended_sensor = ShackHartmannWFS(telescope; n_lenslets=4,
        n_pix_subap=4, T=T)
    extended_rate = shack_hartmann_rate_map(extended_sensor, pupil,
        extended)
    extended_plan = prepare_wfs_optics(
        shack_hartmann_optics(extended_sensor, extended), pupil,
        extended_rate)
    form_wfs_optical_products!(extended_rate, pupil, extended_plan)
    @test all(isfinite, extended_rate.values)
    @test sum(extended_rate.values) > zero(T)
end

@testset "Shack-Hartmann geometric wavefront truth" begin
    T = Float64
    telescope = Telescope(resolution=8, diameter=T(8),
        central_obstruction=zero(T), T=T)
    pupil = PupilFunction(telescope; T=T)
    layout_sensor = ShackHartmannWFS(telescope; n_lenslets=2,
        n_pix_subap=4, T=T)
    pupil.opd .= reshape(T.(1:64), 8, 8) .* T(1e-9)
    truth = zeros(T, 8)
    layout = layout_sensor.front_end.layout
    @test @inferred(geometric_wavefront_slopes!(truth, pupil.opd,
        layout.valid_mask, pupil.metadata.sampling)) === truth
    @test all(isfinite, truth)
    @test any(!iszero, truth)
    if !coverage_instrumented()
        @test @allocated(geometric_wavefront_slopes!(truth, pupil.opd,
            layout.valid_mask, pupil.metadata.sampling)) == 0
    end
    @test !isdefined(WavefrontSensors, :measure!)
    @test !isdefined(WavefrontSensors, :slopes)
end

@testset "Shack-Hartmann centered detector-cell sampling" begin
    T = Float64
    λ = T(0.8e-6)
    diameter = T(1.22 / 16)
    irradiance = T(10)
    source = Source(band=:custom, wavelength=λ,
        photon_irradiance=irradiance, T=T)
    arcsec_per_rad = T(180 * 3600 / π)
    # Integrate a piecewise-constant native FFT intensity by its cumulative
    # integral. This oracle uses the analytic finite-square-aperture DFT,
    # rather than the propagation FFT or the detector sampler.
    function integrated_square_dft(sub, pad, ratio, n_out, tilt_cycles, shift)
        native = [abs2(sum(cispi(2 * k * (tilt_cycles -
            (i - (pad + 1) / 2 + shift * ratio / 2) / pad))
            for k in 0:(sub - 1))) for i in 1:pad]
        cumulative = [zero(T); cumsum(native)]
        function integral(edge)
            clipped = clamp(edge, zero(T), T(pad))
            cell = min(floor(Int, clipped), pad - 1)
            return cumulative[cell + 1] +
                (clipped - cell) * native[cell + 1]
        end
        edges = [(pad - n_out * ratio) / 2 + j * ratio
            for j in 0:n_out]
        return [integral(edges[j + 1]) - integral(edges[j])
            for j in 1:n_out]
    end
    for sub in (15, 16), shannon in (true, false), shift in (false, true),
            requested_scale in (nothing, T(0.937), T(0.697))
        telescope = Telescope(resolution=sub, diameter=diameter,
            central_obstruction=zero(T), T=T)
        pupil = PupilFunction(telescope; T=T)
        fill!(pupil.amplitude, one(T))
        sensor = ShackHartmannWFS(telescope; n_lenslets=1,
            n_pix_subap=22, diffraction_padding=3,
            shannon_sampling=shannon, half_pixel_shift=shift,
            pixel_scale_arcsec=requested_scale, T=T)
        set_valid_subapertures!(sensor, fill(true, 1, 1))
        optics = shack_hartmann_optics(sensor, source)
        rate = shack_hartmann_rate_map(optics, pupil)
        prepared = prepare_wfs_optics(optics, pupil, rate)
        workspace = Optics.microlens_propagation_workspace(optics.propagation)
        pad = size(workspace.fft_stack, 1)
        native_scale = λ / diameter / workspace.effective_padding * arcsec_per_rad
        detector_scale = rate.metadata.sampling[1] * arcsec_per_rad
        if requested_scale !== nothing
            @test detector_scale ≈ requested_scale rtol=4eps(T)
        end
        ratio = detector_scale / native_scale
        for (tilt_x, tilt_y) in ((zero(T), zero(T)), (T(0.007), T(-0.011)))
            # A pupil OPD ramp with cycles per pupil sample known analytically.
            pupil.opd .= [λ * (tilt_x * (i - 1) + tilt_y * (j - 1))
                for i in 1:sub, j in 1:sub]
            form_wfs_optical_products!(rate, pupil, prepared)
            x = integrated_square_dft(sub, pad, ratio, 22, tilt_x, shift)
            y = integrated_square_dft(sub, pad, ratio, 22, tilt_y, shift)
            expected = (irradiance * (diameter / sub)^2 / pad^2) .* (x * y')
            @test rate.values ≈ expected rtol=2e-12 atol=2e-14
            coordinates = T.(1:22) .- T(11.5)
            actual_centroid = (sum(rate.values; dims=2)[:]' * coordinates) / sum(rate.values)
            expected_centroid = dot(x, coordinates) / sum(x)
            @test actual_centroid ≈ expected_centroid atol=2e-12
            if iszero(tilt_x) && iszero(tilt_y) && !shift
                @test rate.values ≈ reverse(rate.values; dims=1) rtol=2e-12 atol=2e-14
                @test abs(actual_centroid) < 2e-12
            end
        end
        form_wfs_optical_products!(rate, pupil, prepared)
        if !coverage_instrumented()
            @test @allocated(form_wfs_optical_products!(rate, pupil, prepared)) == 0
        end
        # A window covering the entire native grid preserves Parseval's
        # incident rate, also for fractional detector/native pixel ratios.
        wide_sensor = ShackHartmannWFS(telescope; n_lenslets=1,
            n_pix_subap=128, diffraction_padding=3,
            shannon_sampling=shannon, half_pixel_shift=shift,
            pixel_scale_arcsec=requested_scale, T=T)
        set_valid_subapertures!(wide_sensor, fill(true, 1, 1))
        wide_optics = shack_hartmann_optics(wide_sensor, source)
        wide_rate = shack_hartmann_rate_map(wide_optics, pupil)
        wide_prepared = prepare_wfs_optics(wide_optics, pupil, wide_rate)
        form_wfs_optical_products!(wide_rate, pupil, wide_prepared)
        @test sum(wide_rate.values) ≈ irradiance * diameter^2 rtol=2e-13
        @test sum(rate.values) <= sum(wide_rate.values) * (1 + 2e-13)
    end
end

@testset "Shack-Hartmann coarse detector sampling" begin
    T = Float64
    telescope = Telescope(resolution=15, diameter=T(1.22 / 16), T=T)
    pupil = PupilFunction(telescope; T=T)
    fill!(pupil.amplitude, one(T))
    source = Source(band=:custom, wavelength=T(0.8e-6),
        photon_irradiance=T(10), T=T)
    # A finite footprint can be much larger than the integer index range.
    # Its four central detector cells still partition the native FFT plane.
    sensor = ShackHartmannWFS(telescope; n_lenslets=1,
        n_pix_subap=22, diffraction_padding=3,
        pixel_scale_arcsec=T(1e300), T=T)
    set_valid_subapertures!(sensor, fill(true, 1, 1))
    optics = shack_hartmann_optics(sensor, source)
    rate = shack_hartmann_rate_map(optics, pupil)
    prepared = prepare_wfs_optics(optics, pupil, rate)
    form_wfs_optical_products!(rate, pupil, prepared)
    incident_rate = T(10) * telescope.params.diameter^2
    @test sum(rate.values) ≈ incident_rate rtol=2e-13
    @test rate.values[11:12, 11:12] ≈ fill(incident_rate / 4, 2, 2) rtol=2e-13
    @test count(!iszero, rate.values) == 4
    invalid_sensor = ShackHartmannWFS(telescope; n_lenslets=1,
        n_pix_subap=22, diffraction_padding=3,
        pixel_scale_arcsec=floatmax(T), T=T)
    @test_throws InvalidConfiguration shack_hartmann_rate_map(
        shack_hartmann_optics(invalid_sensor, source), pupil)
    shifted_sensor = ShackHartmannWFS(telescope; n_lenslets=1,
        n_pix_subap=22, diffraction_padding=3, half_pixel_shift=true,
        pixel_scale_arcsec=floatmax(T) / 4, T=T)
    shifted_optics = shack_hartmann_optics(shifted_sensor, source)
    shifted_rate = shack_hartmann_rate_map(shifted_optics, pupil)
    @test_throws InvalidConfiguration prepare_wfs_optics(
        shifted_optics, pupil, shifted_rate)
end

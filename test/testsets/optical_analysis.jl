@testset "Spatial filter" begin
    tel = Telescope(resolution=8, diameter=8.0, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    sf = SpatialFilter(tel; shape=CircularFilter(), diameter=4, zero_padding=2)
    wavefront = PupilFunction(tel)
    field = ElectricField(wavefront, src; zero_padding=2,
        normalization=DimensionlessNormalization(),
        spatial_measure=PointSampledMeasure(),
        coherence=CoherentFieldCombination())
    formation = prepare_pupil_field(wavefront, src, field;
        center_even_grid=false, amplitude_scale=1)
    fill_electric_field!(field, wavefront, formation)
    filtered = PupilFunction(tel)
    prepared = prepare_spatial_filter(tel, sf, field, filtered)
    filter!(prepared)
    @test size(filtered.opd) == (8, 8)
    @test size(filtered.amplitude) == (8, 8)
end

@testset "Phase statistics" begin
    tel = Telescope(resolution=8, diameter=8.0, central_obstruction=0.0)
    atm = KolmogorovAtmosphere(tel;
        r0=0.2,
        reference_wavelength_m=TEST_ATMOSPHERE_REFERENCE_WAVELENGTH_M,
        L0=25.0,
    )
    rho = [0.0, 1e-6, 0.1, 1.0]
    cov = phase_covariance(rho, atm)
    @test length(cov) == length(rho)
    @test all(isfinite, cov)
    @test cov[1] >= cov[2]
    @test cov[2] <= cov[1]
    @test cov[3] <= cov[2]
    @test cov[4] <= cov[3]
    @test abs(cov[1] - cov[2]) / cov[1] < 1e-3
    var = phase_variance(atm)
    @test var > 0
    psd = phase_spectrum([0.1], atm)
    @test length(psd) == 1
    @test psd[1] > 0
    delta = tel.params.diameter / tel.params.resolution
    screen = ft_phase_screen(atm, 8, 0.1; rng=MersenneTwister(1))
    @test size(screen) == (8, 8)
    ensure_psd!(atm, delta)
    runtime_screen_rng = MersenneTwister(7)
    helper_screen_rng = MersenneTwister(7)
    advance_by!(atm, TEST_ATMOSPHERE_STEP; rng=runtime_screen_rng)
    helper_screen, helper_psd = ft_phase_screen(atm, tel.params.resolution, delta; rng=helper_screen_rng, return_psd=true)
    @test helper_screen ≈ atm.state.phase_rad
    @test helper_psd ≈ atm.state.psd

    for z in (1e-6, 1e-3, 0.1, 1.0, 4.0, 10.0, 40.0, 140.0)
        ref = SpecialFunctions.besselk(5 / 6, z)
        approx = AdaptiveOpticsSim._kv56_scalar(z)
        scaled_ref = z^(5 / 6) * ref
        scaled_approx = AdaptiveOpticsSim._scaled_kv56_scalar(z)
        scaled_cpu = AdaptiveOpticsSim._scaled_kv56_cpu(z)
        @test isapprox(real(approx), ref; rtol=2e-4, atol=1e-10)
        @test isapprox(scaled_approx, scaled_ref; rtol=1e-7, atol=1e-12)
        @test isapprox(scaled_cpu, scaled_ref; rtol=8eps(Float64), atol=1e-12)
    end
end

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
    plan = prepare_spatial_filter(tel, sf, field, filtered)
    workspace = SpatialFilterWorkspace(sf)
    filter!(filtered, field, sf, plan, workspace)
    @test size(filtered.opd) == (8, 8)
    @test size(filtered.amplitude) == (8, 8)
end

@testset "Gain sensing camera" begin
    mask = ones(8, 8)
    basis = rand(8, 8, 3)
    gsc = GainSensingCamera(mask, basis)
    frame = abs.(randn(8, 8))
    calibrate!(gsc, frame)
    og = compute_optical_gains!(gsc, frame)
    @test length(og) == 3
    @test length(AdaptiveOpticsSim.weak_mode_mask(gsc)) == 3
    @test all(isfinite, og)
    @test AdaptiveOpticsSim.detector_metadata(gsc) === nothing

    weak_gsc = GainSensingCamera(mask, zeros(8, 8, 2); sensitivity_floor=1e-6)
    calibrate!(weak_gsc, frame)
    weak_og = compute_optical_gains!(weak_gsc, frame)
    @test all(AdaptiveOpticsSim.weak_mode_mask(weak_gsc))
    @test weak_og == ones(2)

    det = Detector(noise=NoiseReadout(1e-3), integration_time=2.0, qe=0.8, psf_sampling=2, binning=4)
    gsc_with_det = GainSensingCamera(mask, basis; detector=det)
    metadata = AdaptiveOpticsSim.detector_metadata(gsc_with_det)
    @test metadata isa AdaptiveOpticsSim.GSCDetectorMetadata
    @test metadata.integration_time == 2.0
    @test metadata.qe == 0.8
    @test metadata.psf_sampling == 2
    @test metadata.binning == 4
    @test metadata.noise == :readout
    @test metadata.readout_sigma == 1e-3
    @test occursin("psf_sampling=2", sprint(show, MIME"text/plain"(), gsc_with_det))

    AdaptiveOpticsSim.detach_detector!(gsc_with_det)
    @test AdaptiveOpticsSim.detector_metadata(gsc_with_det) === nothing
    AdaptiveOpticsSim.attach_detector!(gsc_with_det, det)
    @test AdaptiveOpticsSim.detector_metadata(gsc_with_det) isa AdaptiveOpticsSim.GSCDetectorMetadata
end

@testset "Phase statistics" begin
    tel = Telescope(resolution=8, diameter=8.0, central_obstruction=0.0)
    atm = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
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
    @test helper_screen ≈ atm.state.opd
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

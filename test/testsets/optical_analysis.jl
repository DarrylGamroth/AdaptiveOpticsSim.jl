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

@testset "Gain sensing camera" begin
    mask = ones(8, 8)
    basis = rand(8, 8, 3)
    gsc = GainSensingCamera(mask, basis)
    frame = abs.(randn(8, 8))
    calibrate!(gsc, frame)
    og = compute_optical_gains!(gsc, frame)
    @test length(og) == 3
    @test length(Calibration.weak_mode_mask(gsc)) == 3
    @test all(isfinite, og)
    @test Calibration.detector_metadata(gsc) === nothing

    mixed_precision_frame = Float32.(frame)
    mixed_precision_gsc = GainSensingCamera(mask, basis)
    calibrate!(mixed_precision_gsc, mixed_precision_frame)
    @test all(isfinite,
        compute_optical_gains!(mixed_precision_gsc, mixed_precision_frame))

    weak_gsc = GainSensingCamera(mask, zeros(8, 8, 2); sensitivity_floor=1e-6)
    calibrate!(weak_gsc, frame)
    weak_og = compute_optical_gains!(weak_gsc, frame)
    @test all(Calibration.weak_mode_mask(weak_gsc))
    @test weak_og == ones(2)

    det = Detector(noise=NoiseReadout(1e-3), exposure_duration=2.0, qe=0.8, psf_sampling=2, binning=4)
    gsc_with_det = GainSensingCamera(mask, basis; detector=det)
    metadata = Calibration.detector_metadata(gsc_with_det)
    @test metadata isa Calibration.GSCDetectorMetadata
    @test metadata.exposure_duration == 2.0
    @test metadata.qe == 0.8
    @test metadata.psf_sampling == 2
    @test metadata.binning == 4
    @test metadata.noise == :readout
    @test metadata.readout_sigma == 1e-3
    @test occursin("psf_sampling=2", sprint(show, MIME"text/plain"(), gsc_with_det))

    Calibration.detach_detector!(gsc_with_det)
    @test Calibration.detector_metadata(gsc_with_det) === nothing
    Calibration.attach_detector!(gsc_with_det, det)
    @test Calibration.detector_metadata(gsc_with_det) isa
        Calibration.GSCDetectorMetadata
end

@testset "Frozen S6 gain-sensing camera CPU source characterization" begin
    fixture = TOML.parsefile(joinpath(@__DIR__, "..", "fixtures",
        "aos_s6_gsc_cpu.toml"))
    @test fixture["schema"] == "test.adaptive-optics-sim/gain-sensing-camera-cpu/1"
    @test fixture["source_repository"] == "AdaptiveOpticsSim.jl"
    @test fixture["source_revision"] ==
        "4dec9c2469cd6f6c3d0a6dfe275ca777901c51c5"
    @test fixture["source_paths"] == ["src/calibration/gain_sensing_camera.jl"]
    @test fixture["source_test_path"] == "test/testsets/optical_analysis.jl"
    @test fixture["array_backend"] == "CPU"
    @test fixture["float_type"] == "Float64"
    @test fixture["storage_order"] == "column-major"
    @test fixture["fft_centering"] == "fftshift(fft(fftshift(x))) / N"

    fixture_array(section) = reshape(Float64.(section["values"]),
        Tuple(Int.(section["shape"])))
    fixture_complex_array(section) = complex.(
        fixture_array(section["real"]), fixture_array(section["imag"]),
    )
    fixture_complex_vector(section) = complex.(Float64.(section["real"]),
        Float64.(section["imag"]))

    mask = fixture_complex_array(fixture["mask"])
    basis = fixture_array(fixture["basis"])
    reference_frame = fixture_array(fixture["frames"]["reference"])
    current_frame = fixture_array(fixture["frames"]["current"])
    reference_ir = fixture_array(fixture["impulse_responses"]["reference"])
    current_ir = fixture_array(fixture["impulse_responses"]["current"])
    reference_sensitivity = fixture_complex_vector(
        fixture["sensitivities"]["reference"])
    current_sensitivity = fixture_complex_vector(
        fixture["sensitivities"]["current"])
    expected_weak = Bool.(fixture["weak_mode_mask"])
    expected_gains = Float64.(fixture["optical_gains"])

    @test size(mask) == (8, 8)
    @test size(basis) == (8, 8, 3)
    @test iseven(size(mask, 1))
    @test fixture["basis"]["alignment"] ==
        "already aligned with the mask and detector-frame grid; no padding"
    @test fixture["mask"]["unit"] ==
        "dimensionless complex focal-plane transmission"
    @test fixture["frames"]["unit"] ==
        "detector intensity sample; normalized spatial distribution is used"
    @test fixture["basis"]["unit"] == "caller-defined aligned modal phase basis"
    @test fixture["sensitivities"]["unit"] == "basis-unit²"
    @test fixture["optical_gains_unit"] == "dimensionless signed ratio"

    gsc = GainSensingCamera(mask, basis;
        sensitivity_floor=fixture["sensitivity_floor"])
    @test gsc.mask == mask
    @test gsc.basis == basis
    @test_throws InvalidConfiguration compute_optical_gains!(gsc, current_frame)
    @test_throws InvalidConfiguration calibrate!(gsc, zeros(size(reference_frame)))

    calibrate!(gsc, reference_frame; n_jobs=1)
    @test sum(reference_frame) ≈ fixture["frames"]["reference"]["total_flux"]
    @test sum(gsc.frame_buffer) ≈ fixture["frames"]["reference"]["normalized_total"]
    @test gsc.ir_calib ≈ reference_ir rtol=1e-12 atol=1e-14
    @test gsc.sensi_calib ≈ reference_sensitivity rtol=1e-12 atol=1e-14
    @test Calibration.weak_mode_mask(gsc) == expected_weak
    @test any(.!Calibration.weak_mode_mask(gsc))
    @test any(Calibration.weak_mode_mask(gsc))

    gains = compute_optical_gains!(gsc, current_frame)
    @test sum(current_frame) ≈ fixture["frames"]["current"]["total_flux"]
    @test sum(gsc.frame_buffer) ≈ fixture["frames"]["current"]["normalized_total"]
    @test gsc.ir_buffer ≈ current_ir rtol=1e-12 atol=1e-14
    @test gsc.sensi_buffer ≈ current_sensitivity rtol=1e-12 atol=1e-14
    @test gains ≈ expected_gains rtol=1e-12 atol=1e-14
    @test gains[end] == 1.0
    @test_throws InvalidConfiguration compute_optical_gains!(gsc,
        zeros(size(current_frame)))

    # The compact fixture characterizes the numerical estimator. Retain the
    # full physical modulation-frame and signed-gain references separately so
    # its transfer to a calibration package cannot silently alter the plant.
    branch = load_reference_bundle(default_reference_root())
    physical_frame = only(filter(case ->
        case.id == fixture["physical_branch_reference"]["current_image_case"],
        branch.cases))
    signed_gains = only(filter(case ->
        case.id == fixture["physical_branch_reference"]["signed_gain_case"],
        branch.cases))
    @test physical_frame.shape == (48, 48)
    @test physical_frame.config["storage_order"] == "C"
    @test size(load_reference_array(physical_frame)) == (48, 48)
    @test all(<(0), load_reference_array(signed_gains))
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

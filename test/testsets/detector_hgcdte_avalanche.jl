function hgcdte_avalanche_samples(detector, input, rng, frame_count::Int)
    samples = Vector{Float64}(undef, frame_count * length(input))
    offset = 0
    for _ in 1:frame_count
        frame = capture!(detector, input, rng)
        copyto!(samples, offset + 1, vec(frame), 1, length(frame))
        offset += length(frame)
    end
    return samples
end

function test_hgcdte_avalanche_moments(samples, expected_mean,
    expected_variance, expected_fourth_central; sigma_limit=6.0)
    sample_count = length(samples)
    mean_limit = sigma_limit * sqrt(expected_variance / sample_count)
    variance_se = sqrt((expected_fourth_central -
        ((sample_count - 3) / (sample_count - 1)) *
        expected_variance^2) / sample_count)
    @test abs(mean(samples) - expected_mean) <= mean_limit
    @test abs(var(samples) - expected_variance) <=
        sigma_limit * variance_se
end

function run_scheduled_hgcdte_avalanche_ramp!(
    prepared, state, values, rng)
    start = PlantTimestamp(0)
    middle = PlantTimestamp(500_000_000)
    close = PlantTimestamp(1_000_000_000)
    begin_exposure!(prepared, state, start)
    take_nondestructive_read!(prepared, state, start, rng)
    accumulate_exposure_interval!(
        prepared, state, start, middle, rng)
    take_nondestructive_read!(prepared, state, middle, rng)
    fill!(values, zero(eltype(values)))
    accumulate_exposure_interval!(
        prepared, state, middle, close, rng)
    take_nondestructive_read!(prepared, state, close, rng)
    cube = detector_ramp_cube(prepared.detector)
    retained = true
    @inbounds for column in axes(cube, 2), row in axes(cube, 1)
        retained &= cube[row, column, 3] == cube[row, column, 2]
    end
    close_exposure!(prepared, state, close)
    complete_readout!(prepared, state, close, rng)
    mark_acquisition_ready!(prepared, state, close)
    return retained
end

function scheduled_hgcdte_avalanche_fixture(seed)
    values = fill(128.0, 8, 8)
    map = detector_test_intensity_map(values;
        kind=DetectorPlane())
    detector = Detector(
        integration_time=1.0,
        qe=1.0,
        noise=NoiseNone(),
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(
            avalanche_gain=2.0,
            excess_noise_factor=1.5,
            multiplication_model=
                ClippedGaussianAvalancheMultiplicationApproximation(),
            read_time=0.0,
            sampling_mode=UpTheRampSampling(3)),
    )
    definition = GlobalShutterAcquisitionDefinition(
        PlantDuration(1_000_000_000))
    prepared = prepare_global_shutter_acquisition(
        detector, map, definition)
    state = GlobalShutterAcquisitionState(prepared)
    return (; detector, prepared, state, values, rng=Xoshiro(seed))
end

@testset "HgCdTe linear-avalanche multiplication qualification" begin
    gamma_model = ConditionalGammaAvalancheMultiplication()
    approximate_model =
        ClippedGaussianAvalancheMultiplicationApproximation()
    @test !Base.isexported(
        AdaptiveOpticsSim, :ConditionalGammaAvalancheMultiplication)
    @test isdefined(
        AdaptiveOpticsSim.Detectors,
        :ConditionalGammaAvalancheMultiplication)
    @test !supports_photon_counting(HgCdTeAvalancheArraySensor())
    @test detector_sensor_symbol(HgCdTeAvalancheArraySensor()) ==
        :hgcdte_linear_avalanche_array

    exact_detector = Detector(
        integration_time=1.0,
        qe=1.0,
        noise=NoiseNone(),
        gain=3.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(
            avalanche_gain=5.0,
            excess_noise_factor=1.0,
            multiplication_model=gamma_model),
    )
    exact_input = fill(7.0, 4, 4)
    @test capture!(exact_detector, exact_input, Xoshiro(6101)) ==
        fill(105.0, 4, 4)

    frame_count = 16
    frame_size = 32
    for (case_index, (charge, avalanche_gain, excess_noise_factor)) in
        enumerate(((1.0, 10.0, 1.5), (8.0, 5.0, 1.25),
            (20.0, 3.0, 1.1)))
        factor_minus_one = excess_noise_factor - 1
        shape = charge / factor_minus_one
        expected_mean = charge * avalanche_gain
        expected_variance =
            charge * avalanche_gain^2 * factor_minus_one
        expected_fourth =
            3expected_variance^2 * (1 + 2 / shape)
        detector = Detector(
            integration_time=1.0,
            qe=1.0,
            noise=NoiseNone(),
            response_model=NullFrameResponse(),
            sensor=HgCdTeAvalancheArraySensor(
                avalanche_gain=avalanche_gain,
                excess_noise_factor=excess_noise_factor,
                multiplication_model=gamma_model),
        )
        samples = hgcdte_avalanche_samples(
            detector, fill(charge, frame_size, frame_size),
            Xoshiro(6110 + case_index), frame_count)
        @test all(>=(0), samples)
        test_hgcdte_avalanche_moments(
            samples, expected_mean, expected_variance,
            expected_fourth)
    end

    approximate_charge = 64.0
    approximate_gain = 5.0
    approximate_factor = 1.5
    approximate_variance = approximate_charge *
        approximate_gain^2 * (approximate_factor - 1)
    approximate_detector = Detector(
        integration_time=1.0,
        qe=1.0,
        noise=NoiseNone(),
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(
            avalanche_gain=approximate_gain,
            excess_noise_factor=approximate_factor,
            multiplication_model=approximate_model),
    )
    approximate_samples = hgcdte_avalanche_samples(
        approximate_detector,
        fill(approximate_charge, frame_size, frame_size),
        Xoshiro(6120), frame_count)
    @test approximate_charge / (approximate_factor - 1) >= 25
    @test all(>=(0), approximate_samples)
    test_hgcdte_avalanche_moments(
        approximate_samples,
        approximate_charge * approximate_gain,
        approximate_variance,
        3approximate_variance^2)

    @test_throws InvalidConfiguration begin
        validate_hgcdte_avalanche_backend(gamma_model, CUDABackend())
    end
    @test isnothing(validate_hgcdte_avalanche_backend(
        approximate_model, CUDABackend()))
end

@testset "HgCdTe linear-avalanche ordering and response" begin
    input = fill(10.0, 4, 4)
    signal_with_read_noise = Detector(
        integration_time=1.0,
        qe=1.0,
        noise=NoiseReadout(2.0),
        gain=3.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(
            avalanche_gain=4.0,
            excess_noise_factor=1.0),
    )
    read_noise_only = Detector(
        integration_time=1.0,
        qe=1.0,
        noise=NoiseReadout(2.0),
        gain=3.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(
            avalanche_gain=1.0,
            excess_noise_factor=1.0),
    )
    signal_output = copy(capture!(
        signal_with_read_noise, input, Xoshiro(6130)))
    noise_output = copy(capture!(
        read_noise_only, zeros(4, 4), Xoshiro(6130)))
    @test signal_output ≈ fill(120.0, 4, 4) + noise_output rtol=0.0 atol=256eps(Float64)

    glow_base = Detector(
        integration_time=1.0,
        qe=1.0,
        noise=NoiseNone(),
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(
            avalanche_gain=1.0, glow_rate=8.0),
    )
    glow_gained = Detector(
        integration_time=1.0,
        qe=1.0,
        noise=NoiseNone(),
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(
            avalanche_gain=5.0, glow_rate=8.0),
    )
    base_glow = copy(capture!(glow_base, zeros(4, 4), Xoshiro(6131)))
    gained_glow = copy(capture!(
        glow_gained, zeros(4, 4), Xoshiro(6131)))
    @test gained_glow == 5base_glow

    saturated = Detector(
        integration_time=1.0,
        qe=1.0,
        noise=NoiseNone(),
        full_well=100.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(avalanche_gain=5.0),
    )
    @test capture!(saturated, fill(50.0, 4, 4), Xoshiro(6132)) ==
        fill(100.0, 4, 4)
    saturated_cube = fill(50.0, 2, 4, 4)
    capture_stack!(saturated, saturated_cube, similar(saturated_cube),
        Xoshiro(6132))
    @test all(==(100.0), saturated_cube)

    response = GaussianPixelResponse(response_width_px=0.7)
    response_detector = Detector(
        qe=1.0,
        noise=NoiseNone(),
        response_model=response,
        sensor=HgCdTeAvalancheArraySensor(avalanche_gain=1.0),
    )
    mtf = detector_mtf(response_detector, 0.17, -0.09)
    capture!(response_detector, ones(4, 4), Xoshiro(6133))
    @test detector_mtf(response_detector, 0.17, -0.09) == mtf

    ipc_detector = Detector(
        qe=1.0,
        noise=NoiseNone(),
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(avalanche_gain=2.0),
        charge_coupling_model=InterpixelCapacitance(
            [0.0 0.01 0.0; 0.01 0.96 0.01; 0.0 0.01 0.0]),
    )
    impulse = zeros(5, 5)
    impulse[3, 3] = 50.0
    ipc_output = capture!(ipc_detector, impulse, Xoshiro(6134))
    @test ipc_output[3, 3] == 96.0
    @test ipc_output[2, 3] == 1.0
    @test sum(ipc_output) == 100.0

    windowed = Detector(
        integration_time=1.0,
        qe=1.0,
        noise=NoiseNone(),
        response_model=NullFrameResponse(),
        readout_window=FrameWindow(2:3, 2:3),
        sensor=HgCdTeAvalancheArraySensor(
            avalanche_gain=2.0,
            read_time=0.25,
            sampling_mode=CorrelatedDoubleSampling()),
    )
    @test capture!(windowed, fill(4.0, 4, 4), Xoshiro(6135)) ==
        fill(8.0, 2, 2)
    metadata = detector_export_metadata(windowed)
    @test metadata.sampling_wallclock_time == 1.5
    @test metadata.frame_response == :none
    @test metadata.charge_coupling == :none
end

@testset "HgCdTe linear-avalanche prepared execution" begin
    values = fill(64.0, 16, 16)
    map = detector_test_intensity_map(values)
    detector = Detector(
        integration_time=0.25,
        qe=1.0,
        noise=NoiseNone(),
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(
            avalanche_gain=5.0,
            excess_noise_factor=1.5,
            multiplication_model=
                ClippedGaussianAvalancheMultiplicationApproximation()),
    )
    plan = prepare_detector_acquisition(detector, map)
    rng = Xoshiro(6140)
    @test @inferred(capture!(detector, map, plan, rng)) isa
        Matrix{Float64}
    @test_detector_allocation prepared_detector_capture_allocations(
        detector, map, plan, rng) == 0

    gamma_detector = Detector(
        integration_time=0.25,
        qe=1.0,
        noise=NoiseNone(),
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(
            avalanche_gain=5.0,
            excess_noise_factor=1.5,
            multiplication_model=
                ConditionalGammaAvalancheMultiplication()),
    )
    gamma_plan = prepare_detector_acquisition(gamma_detector, map)
    gamma_rng = Xoshiro(6141)
    @test @inferred(capture!(
        gamma_detector, map, gamma_plan, gamma_rng)) isa Matrix{Float64}
    @test_detector_allocation prepared_detector_capture_allocations(
        gamma_detector, map, gamma_plan, gamma_rng) == 0

    warm_scheduled = scheduled_hgcdte_avalanche_fixture(6142)
    @test run_scheduled_hgcdte_avalanche_ramp!(
        warm_scheduled.prepared, warm_scheduled.state,
        warm_scheduled.values, warm_scheduled.rng)
    scheduled = scheduled_hgcdte_avalanche_fixture(6143)
    if coverage_instrumented()
        @test_skip "scheduled avalanche allocation assertion is disabled under coverage instrumentation"
    else
        @test @allocated(run_scheduled_hgcdte_avalanche_ramp!(
            scheduled.prepared, scheduled.state,
            scheduled.values, scheduled.rng)) == 0
    end
    @test detector_ramp_acquisition(scheduled.detector) ==
        :scheduled_evolving_charge
end

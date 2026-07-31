function ccd_qualification_samples(detector, input, rng, frame_count::Int)
    samples = Vector{Float64}(undef, frame_count * length(input))
    offset = 0
    for _ in 1:frame_count
        frame = capture!(detector, input, rng)
        copyto!(samples, offset + 1, vec(frame), 1, length(frame))
        offset += length(frame)
    end
    return samples
end

function test_ccd_moments(samples, expected_mean, expected_variance,
    expected_fourth_central; sigma_limit=6.0)
    sample_count = length(samples)
    mean_limit = sigma_limit * sqrt(expected_variance / sample_count)
    variance_se = sqrt((expected_fourth_central -
        ((sample_count - 3) / (sample_count - 1)) * expected_variance^2) /
        sample_count)
    @test abs(mean(samples) - expected_mean) <= mean_limit
    @test abs(var(samples) - expected_variance) <=
        sigma_limit * variance_se
end

function test_ccd_poisson_moments(samples, expected_mean)
    variance = expected_mean
    fourth_central = expected_mean + 3 * expected_mean^2
    return test_ccd_moments(samples, expected_mean, variance,
        fourth_central)
end

function test_ccd_gaussian_moments(samples, expected_sigma)
    variance = expected_sigma^2
    return test_ccd_moments(samples, 0.0, variance, 3 * variance^2)
end

@testset "Conventional CCD contract" begin
    @test_throws InvalidConfiguration CCDSensor(
        clock_induced_charge_per_frame=-1.0)
    @test_throws InvalidConfiguration CCDSensor(
        clock_induced_charge_per_frame=Inf)
    @test_throws InvalidConfiguration CCDSensor(
        clock_induced_charge_per_frame=NaN)
    @test_throws InvalidConfiguration CCDSensor(sample_duration=-1.0)
    @test_throws InvalidConfiguration CCDSensor(sample_duration=Inf)
    @test_throws InvalidConfiguration CCDSensor(sample_duration=1e-3)
    @test_throws InvalidConfiguration CCDSensor(
        sampling_mode=UpTheRampSampling(2))

    sensor = CCDSensor(clock_induced_charge_per_frame=0.25)
    detector = Detector(integration_time=2.0, noise=NoiseNone(), qe=1.0,
        sensor=sensor)
    @test isimmutable(sensor)
    @test detector.params.sensor isa CCDSensor{Float64,SingleRead}
    @test detector.params.sensor.clock_induced_charge_per_frame == 0.25
    @test supports_clock_induced_charge(detector.params.sensor)
    @test !supports_nondestructive_reads(detector.params.sensor)
    @test !supports_multi_read_readout_products(detector.params.sensor)
    metadata = detector_export_metadata(detector)
    @test metadata.sensor == :ccd
    @test metadata.sampling_mode == :single_read
    @test isnothing(metadata.sampling_read_time)
    @test metadata.sampling_wallclock_time == 2.0
    @test metadata.cic_per_frame_law == :none

    default_detector = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, sensor=CCDSensor())
    @test detector_export_metadata(default_detector).frame_response == :none
    @test !supports_detector_mtf(default_detector)
    response = SampledFrameResponse([
        0.00 0.05 0.00
        0.10 0.70 0.05
        0.00 0.10 0.00
    ])
    response_detector = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, sensor=CCDSensor(), response_model=response)
    @test supports_detector_mtf(response_detector)
    @test detector_mtf(response_detector, 0.21, 0.17) ==
        detector_mtf(response, 0.21, 0.17)

    deterministic_input = reshape(collect(1.0:16.0), 4, 4)
    prnu = [0.8 1.0; 1.2 0.9]
    deterministic_detector = Detector(
        integration_time=2.0,
        qe=0.5,
        binning=2,
        noise=NoiseNone(),
        gain=2.0,
        full_well=20.0,
        bits=8,
        output_type=UInt16,
        sensor=CCDSensor(),
        defect_model=PixelResponseNonuniformity(prnu),
    )
    expected = [
        sum(@view deterministic_input[1:2, 1:2]) *
            prnu[1, 1]
        sum(@view deterministic_input[3:4, 1:2]) *
            prnu[2, 1]
        sum(@view deterministic_input[1:2, 3:4]) *
            prnu[1, 2]
        sum(@view deterministic_input[3:4, 3:4]) *
            prnu[2, 2]
    ]
    expected = reshape(expected, 2, 2)
    clamp!(expected, 0.0, 20.0)
    expected .*= 2.0 * 255 / 20
    expected_output = round.(UInt16, clamp.(expected, 0.0, 255.0))
    @test capture!(deterministic_detector, deterministic_input,
        Xoshiro(3001)) == expected_output

    zero_input = zeros(16, 16)
    short_cic = Detector(integration_time=0.1, noise=NoiseNone(), qe=1.0,
        sensor=CCDSensor(clock_induced_charge_per_frame=3.0))
    long_cic = Detector(integration_time=10.0, noise=NoiseNone(), qe=1.0,
        sensor=CCDSensor(clock_induced_charge_per_frame=3.0))
    short_frame = copy(capture!(short_cic, zero_input, Xoshiro(3002)))
    long_frame = copy(capture!(long_cic, zero_input, Xoshiro(3002)))
    @test short_frame == long_frame

    first_forward = Detector(integration_time=0.5, noise=NoiseNone(), qe=1.0,
        sensor=CCDSensor(clock_induced_charge_per_frame=2.0))
    second_forward = Detector(integration_time=4.0, noise=NoiseNone(), qe=1.0,
        sensor=CCDSensor(clock_induced_charge_per_frame=5.0))
    forward = (
        copy(capture!(first_forward, zero_input, Xoshiro(3003))),
        copy(capture!(second_forward, zero_input, Xoshiro(3004))),
    )
    first_reverse = Detector(integration_time=0.5, noise=NoiseNone(), qe=1.0,
        sensor=CCDSensor(clock_induced_charge_per_frame=2.0))
    second_reverse = Detector(integration_time=4.0, noise=NoiseNone(), qe=1.0,
        sensor=CCDSensor(clock_induced_charge_per_frame=5.0))
    reverse_second = copy(capture!(second_reverse, zero_input, Xoshiro(3004)))
    reverse_first = copy(capture!(first_reverse, zero_input, Xoshiro(3003)))
    @test forward == (reverse_first, reverse_second)

    prepared_values = fill(12.0, 16, 16)
    prepared_map = detector_test_intensity_map(prepared_values)
    prepared_detector = Detector(integration_time=0.25,
        noise=NoisePhotonReadout(1.5), qe=0.8, dark_current=4.0,
        sensor=CCDSensor(clock_induced_charge_per_frame=0.75))
    prepared_plan = prepare_detector_acquisition(prepared_detector,
        prepared_map)
    prepared_rng = Xoshiro(3005)
    @test @inferred(capture!(prepared_detector, prepared_map, prepared_plan,
        prepared_rng)) isa Matrix{Float64}
    @test_detector_allocation prepared_detector_capture_allocations(
        prepared_detector, prepared_map, prepared_plan, prepared_rng) == 0

    valid_output = copy(output_frame(prepared_detector))
    mismatched_map = detector_test_intensity_map(copy(prepared_values);
        sampling=(2.0, 1.0))
    @test_throws InvalidConfiguration capture!(prepared_detector,
        mismatched_map, prepared_plan, Xoshiro(3006))
    @test output_frame(prepared_detector) == valid_output
end

@testset "Conventional CCD statistical qualification" begin
    # Each contract uses 16 independent frames of 32x32 independent pixels.
    # Six-standard-error limits are fixed before sampling.
    frame_count = 16
    zero_input = zeros(32, 32)
    signal_input = fill(50.0, 32, 32)

    shot_mean = 20.0
    shot_detector = Detector(integration_time=0.5, qe=0.8,
        noise=NoisePhoton(), sensor=CCDSensor())
    shot_samples = ccd_qualification_samples(shot_detector, signal_input,
        Xoshiro(3101), frame_count)
    test_ccd_poisson_moments(shot_samples, shot_mean)

    dark_mean = 24.0
    dark_detector = Detector(integration_time=2.0, qe=1.0,
        dark_current=12.0, noise=NoiseNone(), sensor=CCDSensor())
    dark_samples = ccd_qualification_samples(dark_detector, zero_input,
        Xoshiro(3102), frame_count)
    test_ccd_poisson_moments(dark_samples, dark_mean)

    cic_mean = 3.5
    cic_detector = Detector(integration_time=7.0, qe=1.0,
        noise=NoiseNone(),
        sensor=CCDSensor(clock_induced_charge_per_frame=cic_mean))
    cic_samples = ccd_qualification_samples(cic_detector, zero_input,
        Xoshiro(3103), frame_count)
    test_ccd_poisson_moments(cic_samples, cic_mean)

    read_sigma = 3.0
    read_detector = Detector(integration_time=1.0, qe=1.0,
        noise=NoiseReadout(read_sigma), sensor=CCDSensor())
    read_samples = ccd_qualification_samples(read_detector, zero_input,
        Xoshiro(3104), frame_count)
    test_ccd_gaussian_moments(read_samples, read_sigma)

    combined_poisson_mean = 20.0 + 6.0 + 2.0
    combined_sigma = 2.5
    combined_variance = combined_poisson_mean + combined_sigma^2
    combined_fourth = combined_poisson_mean + 3 * combined_variance^2
    combined_detector = Detector(integration_time=0.5, qe=0.8,
        dark_current=12.0, noise=NoisePhotonReadout(combined_sigma),
        sensor=CCDSensor(clock_induced_charge_per_frame=2.0))
    combined_samples = ccd_qualification_samples(combined_detector,
        signal_input, Xoshiro(3105), frame_count)
    test_ccd_moments(combined_samples, combined_poisson_mean,
        combined_variance, combined_fourth)
end

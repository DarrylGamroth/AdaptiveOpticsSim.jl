function test_linear_apd_moments(samples, expected_mean, expected_variance;
    sigma_limit=6.0)
    sample_count = length(samples)
    mean_limit = sigma_limit * sqrt(expected_variance / sample_count)
    variance_limit = sigma_limit * expected_variance *
        sqrt(2 / (sample_count - 1))
    @test abs(mean(samples) - expected_mean) <= mean_limit
    @test abs(var(samples) - expected_variance) <= variance_limit
end

@testset "Linear-mode APD ownership and topology" begin
    @test !isdefined(Detectors, :APDDetector)
    @test !isdefined(Detectors, :APDSensor)

    single = LinearAPDDetector(integration_time=0.5, qe=0.5,
        avalanche_gain=4.0, conversion_gain=2.0, noise=NoiseNone())
    output = capture!(single, 10.0; rng=Xoshiro(9001))
    @test output == [20.0]
    @test channel_output(single) === output
    @test supports_avalanche_gain(single)
    metadata = detector_export_metadata(single)
    @test metadata.topology == :single_element
    @test metadata.n_channels == 1
    @test metadata.avalanche_gain == 4.0
    @test metadata.excess_noise_factor == 1.0
    @test !applicable(capture!, single, fill(10.0, 1, 1))

    bank = LinearAPDDetector(topology=LinearAPDChannelBank(4),
        integration_time=1.0, qe=0.5, avalanche_gain=2.0,
        dark_current=1.0, noise=NoiseNone())
    @test capture!(bank, fill(3.0, 4); rng=Xoshiro(9002)) ==
        fill(5.0, 4)
    @test detector_export_metadata(bank).topology == :channel_bank
    @test_throws DimensionMismatchError capture!(bank, 3.0)
    @test_throws DimensionMismatchError capture!(bank, fill(3.0, 3))
    @test !applicable(capture!, bank, fill(3.0, 2, 2))
end

@testset "Linear-mode APD deterministic signal order" begin
    detector = LinearAPDDetector(
        topology=LinearAPDChannelBank(4),
        integration_time=2.0,
        qe=0.5,
        avalanche_gain=4.0,
        excess_noise_factor=1.0,
        dark_current=1.0,
        conversion_gain=2.0,
        noise=NoiseNone(),
    )
    # (3 photons/s * 0.5 * 2 s + 1 electron/s * 2 s) * 4 * 2
    @test capture!(detector, fill(3.0, 4), Xoshiro(9010)) ==
        fill(40.0, 4)
end


@testset "Linear-mode APD statistical moments" begin
    sample_count = 16_384

    multiplied_shot = LinearAPDDetector(
        topology=LinearAPDChannelBank(sample_count),
        integration_time=1.0,
        qe=1.0,
        avalanche_gain=3.0,
        excess_noise_factor=1.4,
        conversion_gain=2.0,
        noise=NoisePhoton(),
    )
    shot_samples = copy(capture!(multiplied_shot,
        fill(80.0, sample_count), Xoshiro(9020)))
    test_linear_apd_moments(shot_samples, 480.0, 4032.0)

    multiplication_only = LinearAPDDetector(
        topology=LinearAPDChannelBank(sample_count),
        avalanche_gain=3.0,
        excess_noise_factor=1.4,
        noise=NoiseNone(),
    )
    multiplication_samples = copy(capture!(multiplication_only,
        fill(80.0, sample_count), Xoshiro(9021)))
    test_linear_apd_moments(multiplication_samples, 240.0, 288.0)

    read_only = LinearAPDDetector(
        topology=LinearAPDChannelBank(sample_count),
        avalanche_gain=3.0,
        conversion_gain=2.0,
        noise=NoiseReadout(2.0),
    )
    read_samples = copy(capture!(read_only,
        fill(10.0, sample_count), Xoshiro(9022)))
    test_linear_apd_moments(read_samples, 60.0, 16.0)
end


@testset "Linear-mode APD failures, replay, and hot path" begin
    @test_throws InvalidConfiguration LinearAPDChannelBank(1)
    for bad in (0.0, -1.0, Inf, NaN)
        @test_throws InvalidConfiguration LinearAPDDetector(
            integration_time=bad)
    end
    for bad in (-0.1, 1.1, Inf, NaN)
        @test_throws InvalidConfiguration LinearAPDDetector(qe=bad)
    end
    for bad in (0.5, Inf, NaN)
        @test_throws InvalidConfiguration LinearAPDDetector(
            avalanche_gain=bad)
        @test_throws InvalidConfiguration LinearAPDDetector(
            excess_noise_factor=bad)
    end
    for bad in (-1.0, Inf, NaN)
        @test_throws InvalidConfiguration LinearAPDDetector(
            dark_current=bad)
    end
    for bad in (0.0, -1.0, Inf, NaN)
        @test_throws InvalidConfiguration LinearAPDDetector(
            conversion_gain=bad)
    end

    function replay_detector()
        return LinearAPDDetector(topology=LinearAPDChannelBank(64),
            integration_time=0.5, qe=0.75, avalanche_gain=5.0,
            excess_noise_factor=1.2, dark_current=0.5,
            conversion_gain=2.0, noise=NoisePhotonReadout(0.25))
    end
    input = fill(40.0, 64)
    @test capture!(replay_detector(), input, Xoshiro(9030)) ==
        capture!(replay_detector(), input, Xoshiro(9030))

    allocation_detector = LinearAPDDetector(
        topology=LinearAPDChannelBank(64), noise=NoiseNone())
    allocation_input = fill(2.0, 64)
    rng = Xoshiro(9031)
    capture!(allocation_detector, allocation_input, rng)
    @test @inferred(capture!(allocation_detector,
        allocation_input, rng)) isa Vector{Float64}
    @test_detector_allocation @allocated(capture!(allocation_detector,
        allocation_input, rng)) == 0
end

function test_mkid_poisson_moments(samples, expected_mean;
    sigma_limit=6.0)
    sample_count = length(samples)
    mean_limit = sigma_limit * sqrt(expected_mean / sample_count)
    variance_limit = sigma_limit * sqrt(
        (expected_mean + 2 * expected_mean^2) / (sample_count - 1))
    @test abs(mean(samples) - expected_mean) <= mean_limit
    @test abs(var(samples) - expected_mean) <= variance_limit
end

@testset "MKID characteristics and observable boundary" begin
    default_characteristics = MKIDArrayCharacteristics()
    @test default_characteristics.energy_resolving_power === nothing
    @test default_characteristics.photon_arrival_time_resolution_s === nothing
    @test default_characteristics.wavelength_passband_m === nothing
    @test isconcretetype(typeof(default_characteristics))
    @test !ismutabletype(typeof(default_characteristics))

    characteristics = MKIDArrayCharacteristics(
        energy_resolving_power=12.0,
        photon_arrival_time_resolution_s=2e-6,
        wavelength_passband_m=(0.8e-6, 1.4e-6))
    sensor = MKIDArraySensor(qe=0.7, dark_count_rate=0.0,
        fill_factor=0.9, characteristics=characteristics)
    detector = MKIDArrayDetector(integration_time=2.0, noise=NoiseNone(),
        sensor=sensor, output_type=UInt16)
    output = capture!(detector, fill(10.0, 2, 8), Xoshiro(9201))

    @test output == fill(UInt16(13), 2, 8)
    @test output_frame(detector) === output
    @test channel_output(detector) === output
    @test isconcretetype(typeof(sensor))
    @test !ismutabletype(typeof(sensor))
    @test !ismutabletype(typeof(detector.params))
    @test supports_photon_counting(sensor)
    @test !supports_energy_resolving(sensor)
    @test !supports_photon_number_resolving(sensor)
    @test !supports_first_order_afterpulse_mean_response(detector)
    @test !supports_nearest_neighbor_count_redistribution(detector)
    @test counting_mean_response_model(detector) isa NullCountingMeanResponse
    @test counting_layout(detector) === :pixel_counts

    metadata = detector_export_metadata(detector)
    @test metadata isa MKIDArrayExportMetadata
    @test metadata.observable === :accumulated_count_image
    @test !metadata.provides_energy_estimates
    @test !metadata.provides_photon_arrival_timestamps
    @test metadata.characteristics.energy_resolving_power == 12.0
    @test metadata.characteristics.photon_arrival_time_resolution_s == 2e-6
    @test metadata.characteristics.wavelength_passband_m ==
        (0.8e-6, 1.4e-6)
    @test metadata.counting.sensor === :mkid_array
    @test metadata.counting.detection_efficiency == 0.7
    @test metadata.counting.fill_factor == 0.9
    @test metadata.counting.mean_response_model === :none
    @test metadata.counting.mean_afterpulses_per_detection === nothing
    @test metadata.counting.redistribution_fraction === nothing
    @test metadata.counting.readout.layout === :pixel_counts
    @test metadata.counting.readout.output_size == (2, 8)
    @test metadata.counting.readout.n_channels == 16

    metadata32 = detector_export_metadata(detector; T=Float32)
    @test metadata32 isa MKIDArrayExportMetadata{Float32}
    @test metadata32.characteristics.energy_resolving_power === 12.0f0
    @test metadata32.characteristics.photon_arrival_time_resolution_s === 2.0f-6
    @test metadata32.characteristics.wavelength_passband_m ==
        (0.8f-6, 1.4f-6)

    @test !applicable(detector_mtf, detector, 0.1, 0.1)
end

@testset "MKID deterministic radiometry and passband" begin
    detector = MKIDArrayDetector(integration_time=2.0, noise=NoiseNone(),
        gate_model=DutyCycleGate(0.25), sensor=MKIDArraySensor(
            qe=0.5, dark_count_rate=4.0, fill_factor=0.8))
    @test counting_exposure_time(detector) == 0.5
    @test capture!(detector, fill(10.0, 2, 8), Xoshiro(9210)) ==
        fill(4.0, 2, 8)
    @test supports_counting_gating(detector)
    @test !supports_dead_time(detector)

    passband = (0.8e-6, 1.4e-6)
    band_detector = MKIDArrayDetector(noise=NoiseNone(),
        sensor=MKIDArraySensor(qe=1.0, dark_count_rate=0.0,
            fill_factor=1.0, characteristics=MKIDArrayCharacteristics(
                wavelength_passband_m=passband)))
    inside = Source(band=:custom, wavelength=1.0e-6)
    lower_boundary = Source(band=:custom, wavelength=first(passband))
    upper_boundary = Source(band=:custom, wavelength=last(passband))
    below = Source(band=:custom, wavelength=0.7e-6)
    above = Source(band=:custom, wavelength=1.5e-6)
    input = fill(2.0, 2, 2)

    @test capture!(band_detector, input, inside, Xoshiro(9211)) == input
    @test capture!(band_detector, input, lower_boundary, Xoshiro(9212)) ==
        input
    @test capture!(band_detector, input, upper_boundary, Xoshiro(9213)) ==
        input
    @test capture!(band_detector, input, below, Xoshiro(9214)) ==
        zeros(2, 2)
    @test capture!(band_detector, input, above, Xoshiro(9215)) ==
        zeros(2, 2)
    # Matrix-only capture is an explicitly prefiltered contract.
    @test capture!(band_detector, input, Xoshiro(9216)) == input

    spectral_source = with_spectrum(inside,
        SpectralBundle([0.6e-6, 0.8e-6, 1.0e-6, 1.4e-6, 1.6e-6],
            [0.1, 0.2, 0.3, 0.15, 0.25]))
    @test counting_source_throughput(
        band_detector, spectral_source) == 0.65
    @test capture!(band_detector, fill(10.0, 2, 2), spectral_source,
        Xoshiro(9217)) ≈ fill(6.5, 2, 2)

    unfiltered_detector = MKIDArrayDetector(noise=NoiseNone(),
        sensor=MKIDArraySensor(qe=1.0, dark_count_rate=0.0,
            fill_factor=1.0))
    @test counting_source_throughput(
        unfiltered_detector, spectral_source) == 1.0
    @test capture!(unfiltered_detector, input, below, Xoshiro(9218)) == input
end

@testset "MKID dead-time mean laws" begin
    dead_time = 0.1
    dimensionless_rates = (0.0, 1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0)
    for (model, expected_response) in (
        (NoDeadTime(), x -> x / dead_time),
        (NonParalyzableDeadTime(dead_time),
            x -> (x / dead_time) / (1 + x)),
        (ParalyzableDeadTime(dead_time),
            x -> (x / dead_time) * exp(-x)),
    )
        detector = MKIDArrayDetector(integration_time=1.0,
            noise=NoiseNone(), sensor=MKIDArraySensor(qe=1.0,
                dark_count_rate=0.0, fill_factor=1.0,
                dead_time_model=model))
        for x in dimensionless_rates
            observed = only(capture!(detector,
                fill(x / dead_time, 1, 1), Xoshiro(9220)))
            @test observed ≈ expected_response(x) rtol=1e-13 atol=1e-13
        end
    end

    nonparalyzable = MKIDArrayDetector(noise=NoiseNone(),
        sensor=MKIDArraySensor(qe=1.0,
            dead_time_model=NonParalyzableDeadTime(dead_time)))
    @test only(capture!(nonparalyzable, fill(1e8, 1, 1),
        Xoshiro(9221))) ≈ inv(dead_time) rtol=1e-6

    paralyzable = MKIDArrayDetector(noise=NoiseNone(),
        sensor=MKIDArraySensor(qe=1.0,
            dead_time_model=ParalyzableDeadTime(dead_time)))
    @test only(capture!(paralyzable, fill(inv(dead_time), 1, 1),
        Xoshiro(9222))) ≈ exp(-1) / dead_time
    @test only(capture!(paralyzable, fill(100 / dead_time, 1, 1),
        Xoshiro(9223))) < 1e-38
end

@testset "MKID Poisson surrogate moments" begin
    dimensions = (128, 128)
    cases = (
        (MKIDArrayDetector(noise=NoisePhoton(),
            sensor=MKIDArraySensor(qe=1.0)),
            fill(20.0, dimensions), 20.0),
        (MKIDArrayDetector(integration_time=2.0, noise=NoisePhoton(),
            sensor=MKIDArraySensor(qe=0.0, dark_count_rate=6.0)),
            zeros(dimensions), 12.0),
        (MKIDArrayDetector(noise=NoisePhoton(),
            sensor=MKIDArraySensor(qe=0.5, dark_count_rate=3.0)),
            fill(10.0, dimensions), 8.0),
        (MKIDArrayDetector(noise=NoisePhoton(),
            sensor=MKIDArraySensor(qe=1.0,
                dead_time_model=NonParalyzableDeadTime(0.05))),
            fill(100.0, dimensions), 100 / 6),
    )
    for (index, (detector, input, expected_mean)) in enumerate(cases)
        samples = vec(copy(capture!(detector, input,
            Xoshiro(9230 + index))))
        test_mkid_poisson_moments(samples, expected_mean)
    end
end

@testset "MKID failures, replay, and hot path" begin
    detector = MKIDArrayDetector(noise=NoiseNone(),
        sensor=MKIDArraySensor(qe=1.0), output_type=UInt16)
    input = fill(2.6, 2, 3)
    @test capture!(detector, input, Xoshiro(9240)) ==
        fill(UInt16(3), 2, 3)
    @test capture!(detector, fill(1e9, 2, 3), Xoshiro(9241)) ==
        fill(typemax(UInt16), 2, 3)

    arrays = (detector.state.counts, detector.state.noise_buffer,
        detector.state.host_buffer, detector.state.output_buffer,
        detector.state.output_buffer_host)
    snapshots = map(copy, arrays)
    for invalid_input in (
        fill(-1.0, 3, 2), fill(Inf, 3, 2), fill(NaN, 3, 2))
        @test_throws InvalidConfiguration capture!(detector,
            invalid_input, Xoshiro(9242))
        current_arrays = (detector.state.counts,
            detector.state.noise_buffer, detector.state.host_buffer,
            detector.state.output_buffer, detector.state.output_buffer_host)
        @test all(current_arrays[index] === arrays[index]
            for index in eachindex(arrays))
        @test all(isequal(current_arrays[index], snapshots[index])
            for index in eachindex(arrays))
    end
    @test size(capture!(detector, fill(1.0, 3, 2), Xoshiro(9243))) ==
        (3, 2)

    for bad in (0.0, -1.0, Inf, NaN)
        @test_throws InvalidConfiguration MKIDArrayDetector(
            integration_time=bad)
    end
    for bad in (-0.1, 1.1, Inf, NaN)
        @test_throws InvalidConfiguration MKIDArraySensor(qe=bad)
    end
    for bad in (-1.0, Inf, NaN)
        @test_throws InvalidConfiguration MKIDArraySensor(
            dark_count_rate=bad)
    end
    for bad in (0.0, -0.1, 1.1, Inf, NaN)
        @test_throws InvalidConfiguration MKIDArraySensor(fill_factor=bad)
    end
    for bad in (-1.0, Inf, NaN)
        @test_throws InvalidConfiguration MKIDArraySensor(
            dead_time_model=NonParalyzableDeadTime(bad))
        @test_throws InvalidConfiguration MKIDArraySensor(
            dead_time_model=ParalyzableDeadTime(bad))
    end
    for bad in (0.0, -1.0, Inf, NaN)
        @test_throws InvalidConfiguration MKIDArrayCharacteristics(
            energy_resolving_power=bad)
        @test_throws InvalidConfiguration MKIDArrayCharacteristics(
            photon_arrival_time_resolution_s=bad)
    end
    for passband in ((1.4e-6, 0.8e-6), (NaN, 1.4e-6),
        (0.8e-6, Inf), (0.0, 1.4e-6))
        @test_throws InvalidConfiguration MKIDArrayCharacteristics(
            wavelength_passband_m=passband)
    end
    for bad in (0.0, -0.1, 1.1, Inf, NaN)
        @test_throws InvalidConfiguration MKIDArrayDetector(
            gate_model=DutyCycleGate(bad))
    end
    @test_throws InvalidConfiguration MKIDArrayDetector(
        noise=NoiseReadout(1.0))

    function replay_detector()
        return MKIDArrayDetector(integration_time=0.5,
            noise=NoisePhoton(), gate_model=DutyCycleGate(0.75),
            sensor=MKIDArraySensor(qe=0.7, dark_count_rate=0.25,
                dead_time_model=NonParalyzableDeadTime(1e-3)))
    end
    replay_input = fill(40.0, 16, 16)
    @test capture!(replay_detector(), replay_input, Xoshiro(9250)) ==
        capture!(replay_detector(), replay_input, Xoshiro(9250))

    allocation_detector = MKIDArrayDetector(noise=NoiseNone(),
        sensor=MKIDArraySensor(qe=0.5))
    allocation_input = fill(2.0, 8, 8)
    rng = Xoshiro(9251)
    capture!(allocation_detector, allocation_input, rng)
    @test @inferred(capture!(allocation_detector,
        allocation_input, rng)) isa Matrix{Float64}
    @test_detector_allocation @allocated(capture!(allocation_detector,
        allocation_input, rng)) == 0
end

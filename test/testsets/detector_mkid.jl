@testset "MKID detector" begin
    mkid_sensor = MKIDArraySensor(qe=0.7, dark_count_rate=0.0, fill_factor=0.9,
        energy_resolution=12.0, timing_jitter_s=2e-6, wavelength_range_m=(0.8e-6, 1.4e-6))
    mkid = MKIDArrayDetector(integration_time=2.0, noise=NoiseNone(), sensor=mkid_sensor,
        output_type=UInt16)
    mkid_out = capture!(mkid, fill(10.0, 2, 8); rng=MersenneTwister(9))
    @test mkid_out == fill(UInt16(13), 2, 8)
    @test output_frame(mkid) === mkid_out
    @test supports_photon_number_resolving(mkid.params.sensor)
    @test supports_energy_resolving(mkid.params.sensor)
    @test !supports_dead_time(mkid)
    mkid_meta = detector_export_metadata(mkid)
    @test mkid_meta isa CountingDetectorExportMetadata
    @test mkid_meta.sensor == :mkid_array
    @test mkid_meta.detection_efficiency == 0.7
    @test mkid_meta.fill_factor == 0.9
    @test mkid_meta.energy_resolution == 12.0
    @test mkid_meta.timing_jitter_s == 2e-6
    @test mkid_meta.wavelength_min_m == 0.8e-6
    @test mkid_meta.wavelength_max_m == 1.4e-6
    @test mkid_meta.readout.output_size == (2, 8)
    @test mkid_meta.readout.n_channels == 16

    mkid_band = MKIDArrayDetector(
        integration_time=1.0,
        noise=NoiseNone(),
        sensor=MKIDArraySensor(qe=1.0, dark_count_rate=0.0, fill_factor=1.0,
            wavelength_range_m=(0.8e-6, 1.4e-6)),
    )
    inside_band = Source(band=:custom, wavelength=1.0e-6)
    outside_band = Source(band=:custom, wavelength=0.55e-6)
    @test capture!(mkid_band, fill(2.0, 2, 2), inside_band, MersenneTwister(10)) == fill(2.0, 2, 2)
    @test capture!(mkid_band, fill(2.0, 2, 2),
        Source(band=:custom, wavelength=0.8e-6), MersenneTwister(10)) == fill(2.0, 2, 2)
    @test capture!(mkid_band, fill(2.0, 2, 2),
        Source(band=:custom, wavelength=1.4e-6), MersenneTwister(10)) == fill(2.0, 2, 2)
    @test capture!(mkid_band, fill(2.0, 2, 2), outside_band, MersenneTwister(10)) == zeros(2, 2)
    spectral_band = with_spectrum(
        inside_band,
        SpectralBundle([0.6e-6, 1.0e-6, 1.6e-6], [0.2, 0.3, 0.5]),
    )
    @test capture!(mkid_band, fill(10.0, 2, 2), spectral_band, MersenneTwister(10)) ≈ fill(3.0, 2, 2)

    mkid_dead = MKIDArrayDetector(
        integration_time=1.0,
        noise=NoiseNone(),
        sensor=MKIDArraySensor(qe=1.0, dark_count_rate=0.0, fill_factor=1.0,
            dead_time_model=ParalyzableDeadTime(0.5)),
    )
    @test capture!(mkid_dead, fill(4.0, 2, 8); rng=MersenneTwister(9)) ≈ fill(4.0 * exp(-2.0), 2, 8)
    @test supports_dead_time(mkid_dead)
    @test supports_paralyzable_dead_time(mkid_dead)

    stochastic_dimensions = (128, 128)
    stochastic_dead = MKIDArrayDetector(
        integration_time=1.0,
        noise=NoisePhoton(),
        sensor=MKIDArraySensor(qe=1.0, dark_count_rate=0.0,
            fill_factor=1.0,
            dead_time_model=NonParalyzableDeadTime(0.05)),
    )
    stochastic_samples = vec(copy(capture!(stochastic_dead,
        fill(100.0, stochastic_dimensions), Xoshiro(11))))
    expected_stochastic_mean = 100 / 6
    stochastic_sample_count = length(stochastic_samples)
    @test abs(mean(stochastic_samples) - expected_stochastic_mean) <=
        6 * sqrt(expected_stochastic_mean / stochastic_sample_count)
    @test abs(var(stochastic_samples) - expected_stochastic_mean) <=
        6 * sqrt((expected_stochastic_mean +
            2 * expected_stochastic_mean^2) /
            (stochastic_sample_count - 1))

    mkid_gate = MKIDArrayDetector(
        integration_time=1.0,
        noise=NoiseNone(),
        gate_model=DutyCycleGate(0.25),
        sensor=MKIDArraySensor(qe=1.0, dark_count_rate=0.0, fill_factor=1.0),
    )
    @test capture!(mkid_gate, fill(8.0, 2, 8); rng=MersenneTwister(9)) == fill(2.0, 2, 8)
    @test supports_counting_gating(mkid_gate)

    @test_throws InvalidConfiguration MKIDArraySensor(qe=1.5)
    @test_throws InvalidConfiguration MKIDArraySensor(dark_count_rate=-1.0)
    @test_throws InvalidConfiguration MKIDArraySensor(fill_factor=0.0)
    @test_throws InvalidConfiguration MKIDArraySensor(energy_resolution=0.0)
    @test_throws InvalidConfiguration MKIDArraySensor(timing_jitter_s=-1.0)
    @test_throws InvalidConfiguration MKIDArraySensor(wavelength_range_m=(1.4e-6, 0.8e-6))
    @test_throws InvalidConfiguration MKIDArraySensor(wavelength_range_m=(NaN, 1.4e-6))
    @test_throws InvalidConfiguration MKIDArraySensor(energy_resolution=Inf)
    @test_throws InvalidConfiguration MKIDArraySensor(timing_jitter_s=NaN)
    @test_throws InvalidConfiguration MKIDArrayDetector(noise=NoiseReadout(1.0))

end

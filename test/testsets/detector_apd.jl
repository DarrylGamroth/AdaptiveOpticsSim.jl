@testset "APD detectors" begin
    arrhenius = ArrheniusRateLaw(300.0, 6000.0)
    @test_throws InvalidConfiguration Detector(integration_time=1.0, noise=NoisePhoton(), qe=1.0, binning=1,
        sensor=APDSensor())

    apd = APDDetector(integration_time=1.0, qe=0.5, gain=2.0, dark_count_rate=0.0, noise=NoiseNone())
    channels = fill(4.0, 2, 8)
    apd_out = capture!(apd, channels; rng=MersenneTwister(9))
    @test apd_out == fill(4.0, 2, 8)
    @test channel_output(apd) === apd_out
    apd_meta = detector_export_metadata(apd)
    @test apd_meta isa CountingDetectorExportMetadata
    @test apd_meta.sensor == :apd
    @test apd_meta.readout.output_size == (2, 8)
    @test apd_meta.readout.n_channels == 16
    @test apd_meta.dead_time_model == :none
    @test apd_meta.dead_time === nothing
    @test supports_counting_noise(apd)
    @test !supports_dead_time(apd)
    @test !supports_channel_gain_map(apd)

    apd_gain_map = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), channel_gain_map=fill(0.5, 2, 8))
    @test capture!(apd_gain_map, fill(2.0, 2, 8); rng=MersenneTwister(9)) == fill(1.0, 2, 8)
    @test supports_channel_gain_map(apd_gain_map)

    linear_apd = LinearAPDDetector(integration_time=0.5, qe=0.5,
        avalanche_gain=4.0, conversion_gain=2.0, noise=NoiseNone())
    linear_apd_out = capture!(linear_apd, 10.0; rng=MersenneTwister(90))
    @test linear_apd_out == [20.0]
    @test channel_output(linear_apd) === linear_apd_out
    @test supports_avalanche_gain(linear_apd)
    linear_apd_meta = detector_export_metadata(linear_apd)
    @test linear_apd_meta.topology == :single_element
    @test linear_apd_meta.n_channels == 1
    linear_apd_rng = MersenneTwister(90)
    capture!(linear_apd, 10.0; rng=linear_apd_rng)
    @test_detector_allocation @allocated(capture!(linear_apd, 10.0;
        rng=linear_apd_rng)) == 0

    linear_apd_bank = LinearAPDDetector(topology=APDChannelBank(4),
        integration_time=1.0, qe=0.5, avalanche_gain=2.0,
        dark_current=1.0, noise=NoiseNone())
    @test capture!(linear_apd_bank, fill(3.0, 4);
        rng=MersenneTwister(91)) == fill(5.0, 4)
    @test detector_export_metadata(linear_apd_bank).topology == :channel_bank
    @test_throws DimensionMismatchError capture!(linear_apd_bank, 3.0)
    @test_throws DimensionMismatchError capture!(linear_apd_bank, fill(3.0, 3))

    linear_apd_noisy = LinearAPDDetector(topology=APDChannelBank(4096),
        noise=NoisePhotonReadout(2.0), avalanche_gain=3.0,
        excess_noise_factor=1.3)
    noisy_apd_out = capture!(linear_apd_noisy, fill(20.0, 4096);
        rng=MersenneTwister(92))
    @test std(noisy_apd_out) > 0
    @test isapprox(mean(noisy_apd_out), 60.0; rtol=0.05)

    @test_throws InvalidConfiguration APDChannelBank(1)
    @test_throws InvalidConfiguration LinearAPDDetector(integration_time=0.0)
    @test_throws InvalidConfiguration LinearAPDDetector(qe=1.1)
    @test_throws InvalidConfiguration LinearAPDDetector(avalanche_gain=0.5)
    @test_throws InvalidConfiguration LinearAPDDetector(excess_noise_factor=0.5)
    @test_throws InvalidConfiguration LinearAPDDetector(dark_current=-1.0)
    @test_throws InvalidConfiguration LinearAPDDetector(conversion_gain=0.0)

    apd_dead_time = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), dead_time_model=NonParalyzableDeadTime(0.5))
    dead_time_out = capture!(apd_dead_time, fill(4.0, 2, 8); rng=MersenneTwister(9))
    @test dead_time_out ≈ fill(4.0 / 3.0, 2, 8)
    dead_time_meta = detector_export_metadata(apd_dead_time)
    @test dead_time_meta.dead_time_model == :nonparalyzable
    @test dead_time_meta.dead_time == 0.5
    @test supports_dead_time(apd_dead_time)
    apd_paralyzable = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), dead_time_model=ParalyzableDeadTime(0.5))
    @test capture!(apd_paralyzable, fill(4.0, 2, 8); rng=MersenneTwister(9)) ≈ fill(4.0 * exp(-2.0), 2, 8)
    @test supports_paralyzable_dead_time(apd_paralyzable)
    apd_gated = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), gate_model=DutyCycleGate(0.5))
    @test capture!(apd_gated, fill(4.0, 2, 8); rng=MersenneTwister(9)) == fill(2.0, 2, 8)
    gated_meta = detector_export_metadata(apd_gated)
    @test gated_meta.gate_model == :duty_cycle
    @test gated_meta.duty_cycle == 0.5
    @test supports_counting_gating(apd_gated)
    apd_afterpulse = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), correlation_model=AfterpulsingModel(0.25))
    @test capture!(apd_afterpulse, fill(4.0, 2, 8); rng=MersenneTwister(9)) == fill(5.0, 2, 8)
    @test supports_afterpulsing(apd_afterpulse)
    apd_crosstalk = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), correlation_model=ChannelCrosstalkModel(0.4))
    crosstalk_in = zeros(3, 3)
    crosstalk_in[2, 2] = 10.0
    crosstalk_out = capture!(apd_crosstalk, crosstalk_in; rng=MersenneTwister(9))
    @test crosstalk_out[2, 2] ≈ 6.0
    @test crosstalk_out[1, 2] ≈ 1.0
    @test crosstalk_out[2, 1] ≈ 1.0
    @test crosstalk_out[2, 3] ≈ 1.0
    @test crosstalk_out[3, 2] ≈ 1.0
    @test supports_channel_crosstalk(apd_crosstalk)
    apd_composite_corr = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), correlation_model=CompositeCountingCorrelation(AfterpulsingModel(0.1), ChannelCrosstalkModel(0.2)))
    composite_meta_apd = detector_export_metadata(apd_composite_corr)
    @test composite_meta_apd.correlation_model == :composite
    @test composite_meta_apd.afterpulse_probability == 0.1
    @test composite_meta_apd.crosstalk == 0.2

    @test_throws InvalidConfiguration APDDetector(noise=NoiseReadout(1.0))
    @test_throws InvalidConfiguration APDDetector(dead_time_model=NonParalyzableDeadTime(-1.0))
    @test_throws InvalidConfiguration APDDetector(dead_time_model=ParalyzableDeadTime(-1.0))
    @test_throws InvalidConfiguration APDDetector(gate_model=DutyCycleGate(0.0))
    @test_throws InvalidConfiguration APDDetector(correlation_model=AfterpulsingModel(1.5))
    @test_throws InvalidConfiguration APDDetector(correlation_model=ChannelCrosstalkModel(-0.1))
    apd_thermal = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=10.0,
        noise=NoiseNone(), thermal_model=FixedTemperature(temperature_K=80.0, dark_count_law=arrhenius))
    apd_thermal_meta = detector_export_metadata(apd_thermal)
    @test supports_detector_thermal_model(apd_thermal)
    @test supports_temperature_dependent_dark_counts(apd_thermal)
    @test detector_temperature(apd_thermal) == 80.0
    @test apd_thermal_meta.thermal_model == :fixed_temperature
    @test apd_thermal_meta.dark_count_law == :arrhenius
    @test effective_dark_count_rate(apd_thermal) < apd_thermal.params.dark_count_rate
    apd_dynamic = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=10.0,
        noise=NoiseNone(), thermal_model=FirstOrderThermalModel(
            ambient_temperature_K=295.0,
            setpoint_temperature_K=120.0,
            initial_temperature_K=300.0,
            time_constant_s=2.0,
            min_temperature_K=80.0,
            max_temperature_K=320.0,
            dark_count_law=arrhenius))
    apd_dynamic_initial = effective_dark_count_rate(apd_dynamic)
    @test supports_detector_thermal_model(apd_dynamic)
    @test supports_dynamic_thermal_state(apd_dynamic.params.thermal_model)
    @test thermal_state(apd_dynamic) isa DetectorThermalState
    @test advance_thermal!(apd_dynamic, 2.0) === apd_dynamic
    @test detector_temperature(apd_dynamic) ≈ 120.0 + 180.0 * exp(-1.0)
    @test effective_dark_count_rate(apd_dynamic) < apd_dynamic_initial
    capture!(apd_dynamic, fill(1.0, 2, 2); rng=MersenneTwister(25))
    @test detector_temperature(apd_dynamic) < 120.0 + 180.0 * exp(-1.0)

end

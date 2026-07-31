@testset "SPAD detector" begin
    arrhenius = ArrheniusRateLaw(300.0, 6000.0)
    spad_sensor = SPADArraySensor(pde=0.5, dark_count_rate=0.0, fill_factor=0.8)
    spad = SPADArrayDetector(integration_time=1.0, noise=NoiseNone(), sensor=spad_sensor)
    spad_out = capture!(spad, fill(10.0, 2, 8); rng=MersenneTwister(9))
    @test spad_out == fill(4.0, 2, 8)
    @test output_frame(spad) === spad_out
    @test channel_output(spad) === spad_out
    spad_meta = detector_export_metadata(spad)
    @test spad_meta isa CountingDetectorExportMetadata
    @test spad_meta.sensor == :spad_array
    @test spad_meta.qe == 0.5
    @test spad_meta.fill_factor == 0.8
    @test spad_meta.gain == 1.0
    @test spad_meta.readout.output_size == (2, 8)
    @test spad_meta.readout.n_channels == 16
    @test !supports_channel_gain_map(spad)
    @test supports_counting_noise(spad)
    @test !supports_dead_time(spad)

    spad_dead = SPADArrayDetector(
        integration_time=1.0,
        noise=NoiseNone(),
        sensor=SPADArraySensor(pde=0.5, dark_count_rate=0.0, fill_factor=0.8,
            dead_time_model=NonParalyzableDeadTime(0.5)),
    )
    @test capture!(spad_dead, fill(10.0, 2, 8); rng=MersenneTwister(9)) ≈ fill(4.0 / 3.0, 2, 8)
    @test supports_dead_time(spad_dead)

    spad_gate = SPADArrayDetector(
        integration_time=1.0,
        noise=NoiseNone(),
        gate_model=DutyCycleGate(0.5),
        sensor=SPADArraySensor(pde=0.5, dark_count_rate=0.0, fill_factor=0.8),
    )
    @test capture!(spad_gate, fill(10.0, 2, 8); rng=MersenneTwister(9)) == fill(2.0, 2, 8)
    @test supports_counting_gating(spad_gate)

    spad_afterpulse = SPADArrayDetector(
        integration_time=1.0,
        noise=NoiseNone(),
        sensor=SPADArraySensor(pde=1.0, dark_count_rate=0.0, fill_factor=1.0,
            correlation_model=AfterpulsingModel(0.25)),
    )
    @test capture!(spad_afterpulse, fill(4.0, 2, 8); rng=MersenneTwister(9)) == fill(5.0, 2, 8)
    @test supports_afterpulsing(spad_afterpulse)

    spad_crosstalk = SPADArrayDetector(
        integration_time=1.0,
        noise=NoiseNone(),
        sensor=SPADArraySensor(pde=1.0, dark_count_rate=0.0, fill_factor=1.0,
            correlation_model=ChannelCrosstalkModel(0.4)),
    )
    spad_crosstalk_in = zeros(3, 3)
    spad_crosstalk_in[2, 2] = 10.0
    spad_crosstalk_out = capture!(spad_crosstalk, spad_crosstalk_in; rng=MersenneTwister(9))
    @test spad_crosstalk_out[2, 2] ≈ 6.0
    @test spad_crosstalk_out[1, 2] ≈ 1.0
    @test spad_crosstalk_out[2, 1] ≈ 1.0
    @test spad_crosstalk_out[2, 3] ≈ 1.0
    @test spad_crosstalk_out[3, 2] ≈ 1.0
    @test sum(spad_crosstalk_out) ≈ sum(spad_crosstalk_in)
    @test supports_channel_crosstalk(spad_crosstalk)

    spad_single = SPADArrayDetector(
        integration_time=1.0,
        noise=NoiseNone(),
        sensor=SPADArraySensor(pde=1.0, dark_count_rate=0.0, fill_factor=1.0,
            correlation_model=ChannelCrosstalkModel(0.4)),
    )
    @test capture!(spad_single, fill(10.0, 1, 1), MersenneTwister(9)) == fill(10.0, 1, 1)

    spad_dynamic = SPADArrayDetector(
        integration_time=1.0,
        noise=NoiseNone(),
        sensor=SPADArraySensor(pde=1.0, dark_count_rate=10.0, fill_factor=1.0),
        thermal_model=FixedTemperature(temperature_K=80.0, dark_count_law=arrhenius),
    )
    spad_dynamic_meta = detector_export_metadata(spad_dynamic)
    @test supports_detector_thermal_model(spad_dynamic)
    @test supports_temperature_dependent_dark_counts(spad_dynamic)
    @test detector_temperature(spad_dynamic) == 80.0
    @test spad_dynamic_meta.thermal_model == :fixed_temperature
    @test spad_dynamic_meta.dark_count_law == :arrhenius
    @test effective_dark_count_rate(spad_dynamic) < spad_dynamic.params.sensor.dark_count_rate

    @test_throws InvalidConfiguration SPADArraySensor(pde=1.5)
    @test_throws InvalidConfiguration SPADArraySensor(dark_count_rate=-1.0)
    @test_throws InvalidConfiguration SPADArraySensor(fill_factor=0.0)
    @test_throws InvalidConfiguration SPADArrayDetector(noise=NoiseReadout(1.0))

end

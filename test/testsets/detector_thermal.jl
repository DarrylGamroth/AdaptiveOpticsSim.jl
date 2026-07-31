@testset "Detector thermal models" begin
    arrhenius = ArrheniusRateLaw(300.0, 6000.0)
    linear = LinearTemperatureLaw(300.0, 0.01)
    exp_law = ExponentialTemperatureLaw(300.0, 0.01)
    @test evaluate_temperature_law(arrhenius, 10.0, 80.0) < 10.0
    @test evaluate_temperature_law(linear, 2.0, 250.0) ≈ 1.0
    @test evaluate_temperature_law(exp_law, 2.0, 250.0) < 2.0

    thermal_det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        dark_current=10.0,
        response_model=NullFrameResponse(),
        thermal_model=FixedTemperature(temperature_K=80.0, dark_current_law=arrhenius),
        sensor=CCDSensor())
    thermal_meta = detector_export_metadata(thermal_det)
    @test supports_detector_thermal_model(thermal_det)
    @test !supports_dynamic_thermal_state(thermal_det.params.thermal_model)
    @test supports_temperature_dependent_dark_current(thermal_det)
    @test detector_temperature(thermal_det) == 80.0
    @test thermal_meta.thermal_model == :fixed_temperature
    @test thermal_meta.detector_temperature_K == 80.0
    @test thermal_meta.cooling_setpoint_K == 80.0
    @test thermal_meta.dark_current_law == :arrhenius
    @test effective_dark_current(thermal_det) < thermal_det.params.dark_current
    @test thermal_state(thermal_det) isa NoThermalState
    @test advance_thermal!(thermal_det, 1.0) === thermal_det

    dynamic_model = FirstOrderThermalModel(
        ambient_temperature_K=295.0,
        setpoint_temperature_K=120.0,
        initial_temperature_K=300.0,
        time_constant_s=2.0,
        min_temperature_K=80.0,
        max_temperature_K=320.0,
        dark_current_law=arrhenius,
    )
    dynamic_det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        dark_current=10.0,
        response_model=NullFrameResponse(),
        thermal_model=dynamic_model,
        sensor=CCDSensor())
    dynamic_meta = detector_export_metadata(dynamic_det)
    @test supports_detector_thermal_model(dynamic_det)
    @test supports_dynamic_thermal_state(dynamic_det.params.thermal_model)
    @test thermal_state(dynamic_det) isa DetectorThermalState
    @test detector_temperature(dynamic_det) == 300.0
    @test dynamic_meta.thermal_model == :first_order
    @test dynamic_meta.detector_temperature_K == 300.0
    @test dynamic_meta.ambient_temperature_K == 295.0
    @test dynamic_meta.cooling_setpoint_K == 120.0
    @test dynamic_meta.thermal_time_constant_s == 2.0
    dark_current_initial = effective_dark_current(dynamic_det)
    @test advance_thermal!(dynamic_det, 2.0) === dynamic_det
    @test detector_temperature(dynamic_det) ≈ 120.0 + 180.0 * exp(-1.0)
    @test effective_dark_current(dynamic_det) < dark_current_initial
    reset_integration!(dynamic_det)
    capture!(dynamic_det, fill(1.0f0, 4, 4); rng=MersenneTwister(24))
    @test detector_temperature(dynamic_det) < 120.0 + 180.0 * exp(-1.0)
    @test_throws InvalidConfiguration advance_thermal!(dynamic_det, -1.0)
    @test_throws InvalidConfiguration FirstOrderThermalModel(
        ambient_temperature_K=295.0,
        setpoint_temperature_K=120.0,
        initial_temperature_K=60.0,
        time_constant_s=2.0,
        min_temperature_K=80.0,
        max_temperature_K=320.0,
    )

    dynamic_stack_det = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, response_model=NullFrameResponse(),
        thermal_model=dynamic_model, sensor=CCDSensor())
    dynamic_stack_cube = zeros(2, 2, 2)
    AdaptiveOpticsSim.Detectors.capture_stack!(dynamic_stack_det, dynamic_stack_cube,
        similar(dynamic_stack_cube); rng=MersenneTwister(24))
    expected_stack_temperature = 120.0 + 180.0 * exp(-0.5)
    @test detector_temperature(dynamic_stack_det) ≈ expected_stack_temperature
    @test readout_ready(dynamic_stack_det)

    dynamic_generalized_det = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, psf_sampling=2,
        response_model=NullFrameResponse(), thermal_model=dynamic_model,
        sensor=CCDSensor())
    dynamic_generalized_input = zeros(2, 4, 4)
    dynamic_generalized_output = zeros(2, 2, 2)
    AdaptiveOpticsSim.Detectors.capture_stack!(dynamic_generalized_det,
        dynamic_generalized_output, dynamic_generalized_input;
        rng=MersenneTwister(24))
    @test detector_temperature(dynamic_generalized_det) ≈
        expected_stack_temperature
    @test readout_ready(dynamic_generalized_det)

    incremental_rate_law = LinearTemperatureLaw(120.0, 0.005)
    incremental_rate_input = zeros(128, 128)
    static_rate_model = FixedTemperature(temperature_K=250.0,
        dark_current_law=incremental_rate_law)
    static_rate_whole = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, dark_current=100.0, response_model=NullFrameResponse(),
        thermal_model=static_rate_model, sensor=CCDSensor())
    static_rate_split = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, dark_current=100.0, response_model=NullFrameResponse(),
        thermal_model=static_rate_model, sensor=CCDSensor())
    static_rate_whole_frame = copy(capture!(static_rate_whole,
        incremental_rate_input; rng=MersenneTwister(25)))
    capture!(static_rate_split, incremental_rate_input;
        rng=MersenneTwister(26), integration_duration=0.5)
    static_rate_split_frame = copy(capture!(static_rate_split,
        incremental_rate_input; rng=MersenneTwister(27), integration_duration=0.5))
    static_rate_expected = evaluate_temperature_law(incremental_rate_law,
        100.0, 250.0)
    @test mean(static_rate_whole_frame) ≈ static_rate_expected rtol=0.01
    @test mean(static_rate_split_frame) ≈ static_rate_expected rtol=0.01
    @test mean(static_rate_split_frame) ≈
        mean(static_rate_whole_frame) rtol=0.01

    incremental_dynamic_model = FirstOrderThermalModel(
        ambient_temperature_K=300.0,
        setpoint_temperature_K=120.0,
        initial_temperature_K=300.0,
        time_constant_s=1.0,
        min_temperature_K=100.0,
        max_temperature_K=320.0,
        dark_current_law=incremental_rate_law,
    )
    dynamic_rate_whole = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, dark_current=100.0, response_model=NullFrameResponse(),
        thermal_model=incremental_dynamic_model, sensor=CCDSensor())
    dynamic_rate_split = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, dark_current=100.0, response_model=NullFrameResponse(),
        thermal_model=incremental_dynamic_model, sensor=CCDSensor())
    dynamic_rate_whole_frame = copy(capture!(dynamic_rate_whole,
        incremental_rate_input; rng=MersenneTwister(28)))
    capture!(dynamic_rate_split, incremental_rate_input;
        rng=MersenneTwister(29), integration_duration=0.5)
    dynamic_rate_split_frame = copy(capture!(dynamic_rate_split,
        incremental_rate_input; rng=MersenneTwister(30), integration_duration=0.5))
    half_exposure_temperature = 120.0 + 180.0 * exp(-0.5)
    dynamic_whole_expected = evaluate_temperature_law(
        incremental_rate_law, 100.0, 300.0)
    dynamic_split_expected = 0.5 * dynamic_whole_expected +
        0.5 * evaluate_temperature_law(incremental_rate_law, 100.0,
            half_exposure_temperature)
    @test mean(dynamic_rate_whole_frame) ≈ dynamic_whole_expected rtol=0.01
    @test mean(dynamic_rate_split_frame) ≈ dynamic_split_expected rtol=0.01
    @test detector_temperature(dynamic_rate_split) ≈
        120.0 + 180.0 * exp(-1.0)

    incremental_glow_model = FirstOrderThermalModel(
        ambient_temperature_K=300.0,
        setpoint_temperature_K=120.0,
        initial_temperature_K=300.0,
        time_constant_s=1.0,
        min_temperature_K=100.0,
        max_temperature_K=320.0,
        glow_rate_law=incremental_rate_law,
    )
    dynamic_glow_split = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, response_model=NullFrameResponse(),
        thermal_model=incremental_glow_model,
        sensor=InGaAsSensor(glow_rate=100.0))
    capture!(dynamic_glow_split, incremental_rate_input;
        rng=MersenneTwister(31), integration_duration=0.5)
    dynamic_glow_frame = copy(capture!(dynamic_glow_split,
        incremental_rate_input; rng=MersenneTwister(32), integration_duration=0.5))
    @test mean(dynamic_glow_frame) ≈ dynamic_split_expected rtol=0.01

    hgcdte_glow_sensor = HgCdTeSensor(glow_rate=60.0,
        read_time=0.1, sampling_mode=CorrelatedDoubleSampling())
    hgcdte_glow_whole = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, dark_current=40.0, response_model=NullFrameResponse(),
        sensor=hgcdte_glow_sensor)
    hgcdte_glow_split = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, dark_current=40.0, response_model=NullFrameResponse(),
        sensor=hgcdte_glow_sensor)
    hgcdte_glow_whole_frame = copy(capture!(hgcdte_glow_whole,
        incremental_rate_input; rng=MersenneTwister(33)))
    capture!(hgcdte_glow_split, incremental_rate_input;
        rng=MersenneTwister(34), integration_duration=0.5)
    hgcdte_glow_split_frame = copy(capture!(hgcdte_glow_split,
        incremental_rate_input; rng=MersenneTwister(35), integration_duration=0.5))
    hgcdte_glow_expected = (40.0 + 60.0) * (1.0 + 2 * 0.1)
    @test mean(hgcdte_glow_whole_frame) ≈ hgcdte_glow_expected rtol=0.01
    @test mean(hgcdte_glow_split_frame) ≈ hgcdte_glow_expected rtol=0.01

    thermal_ingaas = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=NullFrameResponse(),
        thermal_model=FixedTemperature(temperature_K=250.0, glow_rate_law=linear),
        sensor=InGaAsSensor(glow_rate=2.0))
    @test supports_temperature_dependent_glow(thermal_ingaas)
    @test effective_glow_rate(thermal_ingaas) ≈ 1.0

    thermal_emccd = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=NullFrameResponse(),
        thermal_model=FixedTemperature(temperature_K=250.0,
            cic_per_frame_law=linear),
        sensor=EMCCDSensor(clock_induced_charge_per_frame=2.0))
    @test effective_cic_per_frame(thermal_emccd) ≈ 1.0

end

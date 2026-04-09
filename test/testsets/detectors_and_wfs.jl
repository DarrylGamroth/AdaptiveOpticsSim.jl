@testset "Detector" begin
    psf = fill(1.0, 8, 8)
    det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=2)
    frame = capture!(det, psf; rng=MersenneTwister(2))
    @test size(frame) == (4, 4)
    @test sum(frame) == sum(psf)

    det_sampling = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, psf_sampling=2, binning=2)
    frame_sampling = capture!(det_sampling, psf; rng=MersenneTwister(2))
    @test size(frame_sampling) == (2, 2)
    @test sum(frame_sampling) == sum(psf)

    det_tuple = Detector(integration_time=1.0, noise=(NoisePhoton(), NoiseReadout(0.5)),
        qe=1.0, binning=1)
    @test det_tuple.noise isa NoisePhotonReadout
    @test AdaptiveOpticsSim.detector_execution_plan(AdaptiveOpticsSim.execution_style(det_tuple.state.frame), det_tuple) isa AdaptiveOpticsSim.DetectorDirectPlan

    det_sat = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1, full_well=5.0)
    frame_sat = capture!(det_sat, fill(10.0, 4, 4); rng=MersenneTwister(2))
    @test maximum(frame_sat) == 5.0

    det_adc = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1, bits=8, full_well=10.0)
    frame_adc = capture!(det_adc, fill(10.0, 4, 4); rng=MersenneTwister(2))
    @test frame_adc isa Matrix{UInt8}
    @test output_frame(det_adc) === frame_adc
    @test maximum(frame_adc) == 0xff
    @test minimum(frame_adc) >= 0x00
    @test eltype(det_adc.state.frame) == Float64
    metadata_adc = detector_export_metadata(det_adc)
    @test metadata_adc.noise == :none
    @test metadata_adc.sensor == :ccd
    @test metadata_adc.output_precision == UInt8
    @test metadata_adc.frame_size == (4, 4)
    @test metadata_adc.output_size == (4, 4)

    det_adc_float = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        bits=8, full_well=10.0, output_precision=Float32)
    frame_adc_float = capture!(det_adc_float, fill(10.0, 4, 4); rng=MersenneTwister(2))
    @test frame_adc_float isa Matrix{Float32}
    @test maximum(frame_adc_float) == Float32(255.0)

    det_window = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        readout_window=FrameWindow(2:3, 2:4))
    psf_window = reshape(collect(1.0:16.0), 4, 4)
    frame_window = copy(capture!(det_window, psf_window; rng=MersenneTwister(2)))
    @test size(frame_window) == (2, 3)
    @test frame_window == psf_window[2:3, 2:4]
    meta_window = detector_export_metadata(det_window)
    @test meta_window.window_rows == (2, 3)
    @test meta_window.window_cols == (2, 4)
    det_window_oob = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        readout_window=FrameWindow(2:5, 1:2))
    @test_throws DimensionMismatchError capture!(det_window_oob, psf_window; rng=MersenneTwister(2))
    @test_throws InvalidConfiguration FrameWindow(0:1, 1:2)

    det_dark = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1, dark_current=100.0)
    frame_dark = capture!(det_dark, zeros(4, 4); rng=MersenneTwister(2))
    @test sum(frame_dark) > 0

    zero_psf = zeros(4, 4)
    rng_ccd = MersenneTwister(7)
    rng_emccd = MersenneTwister(7)
    det_ccd = Detector(integration_time=1.0, noise=NoiseReadout(1.0), qe=1.0, binning=1,
        gain=10.0, sensor=CCDSensor())
    det_emccd = Detector(integration_time=1.0, noise=NoiseReadout(1.0), qe=1.0, binning=1,
        gain=10.0, sensor=EMCCDSensor())
    frame_ccd = copy(capture!(det_ccd, zero_psf; rng=rng_ccd))
    frame_emccd = copy(capture!(det_emccd, zero_psf; rng=rng_emccd))
    @test frame_ccd ≈ 10 .* frame_emccd
    @test_throws InvalidConfiguration EMCCDSensor(excess_noise_factor=0.5)

    uniform_signal = fill(50.0, 8, 8)
    det_emccd_base = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=EMCCDSensor())
    det_emccd_excess = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=EMCCDSensor(excess_noise_factor=sqrt(2.0)))
    frame_emccd_base = copy(capture!(det_emccd_base, uniform_signal; rng=MersenneTwister(8)))
    frame_emccd_excess = copy(capture!(det_emccd_excess, uniform_signal; rng=MersenneTwister(8)))
    @test frame_emccd_base == uniform_signal
    @test std(vec(frame_emccd_excess)) > 0
    det_emccd_stochastic = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=5.0, sensor=EMCCDSensor(excess_noise_factor=1.4,
            multiplication_model=StochasticMultiplicationRegister(0.6)))
    frame_emccd_stochastic = copy(capture!(det_emccd_stochastic, uniform_signal; rng=MersenneTwister(124)))
    @test std(vec(frame_emccd_stochastic)) > 0

    det_ccd_cic = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=CCDSensor(clock_induced_charge_rate=5.0))
    frame_ccd_cic = capture!(det_ccd_cic, zero_psf; rng=MersenneTwister(11))
    @test sum(frame_ccd_cic) > 0
    @test supports_clock_induced_charge(det_ccd_cic.params.sensor)
    @test_throws InvalidConfiguration CCDSensor(clock_induced_charge_rate=-1.0)
    det_emccd_cic = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=EMCCDSensor(cic_rate=3.0))
    @test sum(capture!(det_emccd_cic, zero_psf; rng=MersenneTwister(125))) > 0
    det_emccd_sat = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=5.0, sensor=EMCCDSensor(register_full_well=100.0))
    @test maximum(capture!(det_emccd_sat, fill(50.0, 4, 4); rng=MersenneTwister(126))) == 100.0
    @test_throws InvalidConfiguration EMCCDSensor(cic_rate=-1.0)
    @test_throws InvalidConfiguration EMCCDSensor(register_full_well=0.0)
    @test_throws InvalidConfiguration EMCCDSensor(multiplication_model=StochasticMultiplicationRegister(-1.0))

    det_cmos = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=CMOSSensor(column_readout_sigma=1.0))
    frame_cmos = copy(capture!(det_cmos, zeros(8, 8); rng=MersenneTwister(12)))
    @test !all(iszero, frame_cmos)
    @test all(j -> isapprox(std(frame_cmos[:, j]), 0.0; atol=1e-8), axes(frame_cmos, 2))
    @test std(vec(frame_cmos[1, :])) > 0
    @test supports_column_readout_noise(det_cmos.params.sensor)
    @test detector_export_metadata(det_cmos).frame_response == :gaussian
    @test_throws InvalidConfiguration CMOSSensor(column_readout_sigma=-1.0)
    prnu_map = [1.0 0.5 1.0 0.5; 1.0 0.5 1.0 0.5; 1.0 0.5 1.0 0.5; 1.0 0.5 1.0 0.5]
    dsnu_map = fill(0.25, 4, 4)
    bad_mask = falses(4, 4)
    bad_mask[2, 3] = true
    det_cmos_structured = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=CMOSSensor(output_model=StaticCMOSOutputPattern(2, [1.0, 2.0], [0.0, 10.0]),
            timing_model=RollingShutter(1e-3)),
        response_model=NullFrameResponse(),
        defect_model=CompositeDetectorDefectModel(
            PixelResponseNonuniformity(prnu_map),
            DarkSignalNonuniformity(dsnu_map),
            BadPixelMask(bad_mask; throughput=0.0)))
    structured_frame = capture!(det_cmos_structured, fill(2.0, 4, 4); rng=MersenneTwister(120))
    @test structured_frame[1, 1] ≈ 2.25
    @test structured_frame[1, 2] ≈ 1.25
    @test structured_frame[1, 3] ≈ 14.5
    @test structured_frame[2, 3] ≈ 10.5
    structured_meta = detector_export_metadata(det_cmos_structured)
    @test structured_meta.detector_defects == :composite
    @test structured_meta.has_prnu
    @test structured_meta.has_dsnu
    @test structured_meta.has_bad_pixels
    @test structured_meta.timing_model == :rolling_shutter
    @test structured_meta.timing_line_time == 1e-3
    @test structured_meta.sampling_wallclock_time == 1.004
    @test supports_detector_defect_maps(det_cmos_structured.params.sensor)
    @test supports_shutter_timing(det_cmos_structured.params.sensor)
    @test_throws InvalidConfiguration CMOSSensor(timing_model=RollingShutter(-1.0))
    @test_throws InvalidConfiguration CMOSSensor(output_model=StaticCMOSOutputPattern(2, [1.0], [0.0, 1.0]))

    det_ingaas = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=InGaAsSensor(glow_rate=3.0))
    frame_ingaas = capture!(det_ingaas, zero_psf; rng=MersenneTwister(13))
    @test sum(frame_ingaas) > 0
    @test supports_sensor_glow(det_ingaas.params.sensor)
    @test detector_export_metadata(det_ingaas).frame_response == :gaussian
    @test_throws InvalidConfiguration InGaAsSensor(glow_rate=-1.0)
    det_ingaas_persist = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=NullFrameResponse(),
        sensor=InGaAsSensor(persistence_model=ExponentialPersistence(0.5, 0.0)))
    capture!(det_ingaas_persist, fill(4.0, 4, 4); rng=MersenneTwister(121))
    persisted = capture!(det_ingaas_persist, zeros(4, 4); rng=MersenneTwister(122))
    @test sum(persisted) ≈ 32.0
    persist_meta = detector_export_metadata(det_ingaas_persist)
    @test persist_meta.persistence_model == :exponential
    det_ingaas_nonlinear = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=NullFrameResponse(),
        nonlinearity_model=SaturatingFrameNonlinearity(0.1),
        sensor=InGaAsSensor())
    nonlinear_frame = capture!(det_ingaas_nonlinear, fill(10.0, 2, 2); rng=MersenneTwister(123))
    @test nonlinear_frame == fill(5.0, 2, 2)
    nonlinear_meta = detector_export_metadata(det_ingaas_nonlinear)
    @test nonlinear_meta.nonlinearity_model == :saturating
    @test supports_detector_persistence(det_ingaas_persist.params.sensor)
    @test supports_detector_nonlinearity(det_ingaas_nonlinear.params.sensor)
    @test_throws InvalidConfiguration InGaAsSensor(persistence_model=ExponentialPersistence(1.1, 0.0))
    @test_throws InvalidConfiguration Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        nonlinearity_model=SaturatingFrameNonlinearity(-0.1), sensor=InGaAsSensor())

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

    thermal_ingaas = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=NullFrameResponse(),
        thermal_model=FixedTemperature(temperature_K=250.0, glow_rate_law=linear),
        sensor=InGaAsSensor(glow_rate=2.0))
    @test supports_temperature_dependent_glow(thermal_ingaas)
    @test effective_glow_rate(thermal_ingaas) ≈ 1.0

    thermal_emccd = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=NullFrameResponse(),
        thermal_model=FixedTemperature(temperature_K=250.0, cic_rate_law=linear),
        sensor=EMCCDSensor(cic_rate=2.0))
    @test effective_cic_rate(thermal_emccd) ≈ 1.0

    det_saphira = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(avalanche_gain=5.0))
    frame_saphira = copy(capture!(det_saphira, uniform_signal; rng=MersenneTwister(14)))
    @test frame_saphira == 5.0 .* uniform_signal
    @test supports_avalanche_gain(det_saphira.params.sensor)
    @test supports_sensor_glow(det_saphira.params.sensor)
    @test detector_export_metadata(det_saphira).frame_response == :sampled
    det_saphira_excess = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(avalanche_gain=1.0, excess_noise_factor=sqrt(2.0)))
    frame_saphira_excess = copy(capture!(det_saphira_excess, uniform_signal; rng=MersenneTwister(14)))
    @test std(vec(frame_saphira_excess)) > 0
    det_saphira_sat = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, full_well=100.0, sensor=HgCdTeAvalancheArraySensor(avalanche_gain=5.0))
    frame_saphira_sat = copy(capture!(det_saphira_sat, uniform_signal; rng=MersenneTwister(15)))
    @test maximum(frame_saphira_sat) == 100.0
    det_saphira_single = Detector(integration_time=1.0, noise=NoiseReadout(4.0), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor())
    frame_saphira_single = copy(capture!(det_saphira_single, zero_psf; rng=MersenneTwister(16)))
    single_products = readout_products(det_saphira_single)
    @test single_products isa HgCdTeReadoutProducts
    @test detector_reference_frame(det_saphira_single) === nothing
    @test detector_signal_frame(det_saphira_single) !== nothing
    @test detector_combined_frame(det_saphira_single) == frame_saphira_single
    @test detector_reference_cube(det_saphira_single) === nothing
    @test detector_signal_cube(det_saphira_single) !== nothing
    @test detector_read_cube(det_saphira_single) === nothing
    @test detector_read_times(det_saphira_single) === nothing
    det_saphira_ndr = Detector(integration_time=1.0, noise=NoiseReadout(4.0), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(sampling_mode=AveragedNonDestructiveReads(4)))
    frame_saphira_ndr = copy(capture!(det_saphira_ndr, zero_psf; rng=MersenneTwister(16)))
    @test std(vec(frame_saphira_ndr)) < std(vec(frame_saphira_single))
    @test supports_nondestructive_reads(det_saphira_ndr.params.sensor)
    @test supports_readout_correction(det_saphira_ndr.params.sensor)
    @test supports_read_cube(det_saphira_ndr.params.sensor)
    saphira_meta = detector_export_metadata(det_saphira_ndr)
    @test saphira_meta.sampling_mode == :averaged_non_destructive_reads
    @test saphira_meta.sampling_reads == 4
    @test saphira_meta.sampling_reference_reads == 0
    @test saphira_meta.sampling_signal_reads == 4
    @test saphira_meta.readout_sigma == 2.0
    @test saphira_meta.provides_signal_frame
    @test !saphira_meta.provides_reference_frame
    @test saphira_meta.provides_combined_frame
    @test !saphira_meta.provides_reference_cube
    @test saphira_meta.provides_signal_cube
    @test saphira_meta.provides_read_cube
    @test saphira_meta.signal_cube_reads == 4
    @test saphira_meta.read_cube_reads == 4
    @test detector_combined_frame(det_saphira_ndr) == frame_saphira_ndr
    @test size(detector_signal_cube(det_saphira_ndr)) == (4, 4, 4)
    @test length(detector_read_times(det_saphira_ndr)) == 4
    det_saphira_cds = Detector(integration_time=1.0, noise=NoiseReadout(4.0), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(sampling_mode=CorrelatedDoubleSampling()))
    frame_saphira_cds = copy(capture!(det_saphira_cds, zero_psf; rng=MersenneTwister(16)))
    @test std(vec(frame_saphira_cds)) > std(vec(frame_saphira_single))
    cds_meta = detector_export_metadata(det_saphira_cds)
    @test cds_meta.sampling_mode == :correlated_double_sampling
    @test cds_meta.sampling_reference_reads == 1
    @test cds_meta.sampling_signal_reads == 1
    @test cds_meta.readout_sigma ≈ 4.0 * sqrt(2.0)
    @test cds_meta.provides_reference_cube
    @test cds_meta.provides_signal_cube
    @test cds_meta.reference_cube_reads == 1
    @test cds_meta.signal_cube_reads == 1
    @test size(detector_reference_cube(det_saphira_cds)) == (4, 4, 1)
    @test size(detector_signal_cube(det_saphira_cds)) == (4, 4, 1)
    @test detector_combined_frame(det_saphira_cds) ≈
        detector_signal_frame(det_saphira_cds) .- detector_reference_frame(det_saphira_cds)
    det_saphira_fowler = Detector(integration_time=1.0, noise=NoiseReadout(4.0), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(sampling_mode=FowlerSampling(8)))
    frame_saphira_fowler = copy(capture!(det_saphira_fowler, zero_psf; rng=MersenneTwister(16)))
    fowler_meta = detector_export_metadata(det_saphira_fowler)
    @test fowler_meta.sampling_mode == :fowler_sampling
    @test fowler_meta.sampling_reads == 16
    @test fowler_meta.sampling_reference_reads == 8
    @test fowler_meta.sampling_signal_reads == 8
    @test fowler_meta.readout_sigma == 2.0
    @test fowler_meta.reference_cube_reads == 8
    @test fowler_meta.signal_cube_reads == 8
    @test detector_combined_frame(det_saphira_fowler) ≈
        detector_signal_frame(det_saphira_fowler) .- detector_reference_frame(det_saphira_fowler)
    det_saphira_timed_single = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        dark_current=1000.0, gain=1.0, sensor=HgCdTeAvalancheArraySensor(read_time=1.0))
    frame_saphira_timed_single = copy(capture!(det_saphira_timed_single, zero_psf; rng=MersenneTwister(17)))
    det_saphira_timed_cds = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        dark_current=1000.0, gain=1.0,
        sensor=HgCdTeAvalancheArraySensor(read_time=1.0, sampling_mode=CorrelatedDoubleSampling()))
    frame_saphira_timed_cds = copy(capture!(det_saphira_timed_cds, zero_psf; rng=MersenneTwister(17)))
    @test sum(frame_saphira_timed_cds) > sum(frame_saphira_timed_single)
    timed_meta = detector_export_metadata(det_saphira_timed_cds)
    @test timed_meta.sampling_read_time == 1.0
    @test timed_meta.sampling_wallclock_time == 3.0
    det_saphira_windowed = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(read_time=1.0, sampling_mode=CorrelatedDoubleSampling()),
        readout_window=FrameWindow(2:3, 2:3))
    capture!(det_saphira_windowed, fill(10.0, 4, 4); rng=MersenneTwister(18))
    windowed_meta = detector_export_metadata(det_saphira_windowed)
    @test windowed_meta.sampling_read_time == 0.5
    @test windowed_meta.sampling_wallclock_time == 2.0
    @test detector_combined_frame(det_saphira_windowed) !== nothing
    @test size(detector_reference_cube(det_saphira_windowed)) == (2, 2, 1)
    @test size(detector_signal_cube(det_saphira_windowed)) == (2, 2, 1)
    @test size(detector_signal_frame(det_saphira_windowed)) == (2, 2)
    @test size(detector_read_cube(det_saphira_windowed)) == (2, 2, 2)
    @test detector_read_times(det_saphira_windowed) == [0.5, 1.0]
    det_saphira_corrected = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(),
        response_model=NullFrameResponse(),
        correction_model=ReferencePixelCommonModeCorrection(1, 1))
    corrected_frame = capture!(det_saphira_corrected, fill(5.0, 4, 4); rng=MersenneTwister(19))
    corrected_meta = detector_export_metadata(det_saphira_corrected)
    @test corrected_meta.readout_correction == :reference_pixel_common_mode
    @test corrected_meta.correction_edge_rows == 1
    @test corrected_meta.correction_edge_cols == 1
    @test abs(mean(corrected_frame)) < 1e-6
    det_saphira_row_corrected = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(),
        response_model=NullFrameResponse(),
        correction_model=ReferenceRowCommonModeCorrection(1))
    row_pattern = repeat(reshape([1.0, 2.0, 3.0, 4.0], :, 1), 1, 4)
    row_corrected = capture!(det_saphira_row_corrected, row_pattern; rng=MersenneTwister(20))
    @test maximum(abs, row_corrected) < 1e-6
    det_saphira_col_corrected = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(),
        response_model=NullFrameResponse(),
        correction_model=ReferenceColumnCommonModeCorrection(1))
    col_pattern = repeat(reshape([1.0, 2.0, 3.0, 4.0], 1, :), 4, 1)
    col_corrected = capture!(det_saphira_col_corrected, col_pattern; rng=MersenneTwister(21))
    @test maximum(abs, col_corrected) < 1e-6
    det_saphira_output_corrected = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(),
        response_model=NullFrameResponse(),
        correction_model=ReferenceOutputCommonModeCorrection(2; edge_rows=1, edge_cols=1))
    output_pattern = hcat(fill(5.0, 4, 2), fill(10.0, 4, 2))
    output_corrected = capture!(det_saphira_output_corrected, output_pattern; rng=MersenneTwister(22))
    output_meta = detector_export_metadata(det_saphira_output_corrected)
    @test output_meta.readout_correction == :reference_output_common_mode
    @test output_meta.correction_group_cols == 2
    @test maximum(abs, output_corrected) < 1e-6
    det_saphira_composite = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(),
        response_model=NullFrameResponse(),
        correction_model=CompositeFrameReadoutCorrection((
            ReferenceRowCommonModeCorrection(1),
            ReferenceColumnCommonModeCorrection(1))))
    composite_pattern = row_pattern .+ col_pattern
    composite_corrected = capture!(det_saphira_composite, composite_pattern; rng=MersenneTwister(23))
    composite_meta = detector_export_metadata(det_saphira_composite)
    @test composite_meta.readout_correction == :composite
    @test composite_meta.correction_stage_count == 2
    @test maximum(abs, composite_corrected) < 1e-6
    @test_throws InvalidConfiguration HgCdTeAvalancheArraySensor(avalanche_gain=0.5)
    @test_throws InvalidConfiguration HgCdTeAvalancheArraySensor(excess_noise_factor=0.5)
    @test_throws InvalidConfiguration HgCdTeAvalancheArraySensor(glow_rate=-1.0)
    @test_throws InvalidConfiguration HgCdTeAvalancheArraySensor(read_time=-1.0)
    @test_throws InvalidConfiguration HgCdTeAvalancheArraySensor(sampling_mode=AveragedNonDestructiveReads(0))
    @test_throws InvalidConfiguration HgCdTeAvalancheArraySensor(sampling_mode=FowlerSampling(0))
    @test_throws InvalidConfiguration ReferencePixelCommonModeCorrection(0, 0)
    @test_throws InvalidConfiguration ReferenceRowCommonModeCorrection(0)
    @test_throws InvalidConfiguration ReferenceColumnCommonModeCorrection(0)
    @test_throws InvalidConfiguration ReferenceOutputCommonModeCorrection(0)
    @test_throws InvalidConfiguration CompositeFrameReadoutCorrection(())

    @test_throws InvalidConfiguration Detector(integration_time=1.0, noise=NoisePhoton(), qe=1.0, binning=1,
        sensor=APDSensor())

    det_default_ccd = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1, sensor=CCDSensor())
    @test detector_export_metadata(det_default_ccd).frame_response == :none

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

    det_mtf = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=GaussianPixelResponse(response_width_px=0.75))
    impulse = zeros(9, 9)
    impulse[5, 5] = 1.0
    frame_mtf = capture!(det_mtf, impulse; rng=MersenneTwister(3))
    @test sum(frame_mtf) ≈ 1.0 atol=1e-6
    @test frame_mtf[5, 5] < 1.0
    @test frame_mtf[5, 4] > 0
    @test supports_detector_mtf(det_mtf)
    mtf_meta = detector_export_metadata(det_mtf)
    @test mtf_meta.frame_response == :gaussian
    @test mtf_meta.response_width_px == 0.75
    @test mtf_meta.response_application_domain == :image
    @test mtf_meta.response_is_separable
    @test mtf_meta.response_is_shift_invariant
    @test mtf_meta.response_support_rows == mtf_meta.response_support_cols
    @test mtf_meta.pitch_x_px === nothing
    @test mtf_meta.aperture_shape === nothing
    @test_throws InvalidConfiguration GaussianPixelResponse(response_width_px=0.0)

    sampled_kernel = [0.0 0.1 0.0; 0.1 0.6 0.1; 0.0 0.1 0.0]
    sampled_det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=SampledFrameResponse(sampled_kernel))
    sampled_frame = capture!(sampled_det, impulse; rng=MersenneTwister(6))
    @test sum(sampled_frame) ≈ 1.0 atol=1e-6
    @test sampled_frame[5, 5] ≈ 0.6 atol=1e-6
    @test sampled_frame[5, 4] ≈ 0.1 atol=1e-6
    sampled_meta = detector_export_metadata(sampled_det)
    @test sampled_meta.frame_response == :sampled
    @test sampled_meta.response_application_domain == :image
    @test !sampled_meta.response_is_separable
    @test sampled_meta.response_support_rows == 3
    @test sampled_meta.response_support_cols == 3
    @test sampled_meta.aperture_shape == :sampled
    @test supports_detector_mtf(sampled_det)
    @test_throws InvalidConfiguration SampledFrameResponse(zeros(3, 3))
    @test_throws InvalidConfiguration SampledFrameResponse(ones(2, 3))

    cube_mtf = Array{Float64}(undef, 2, size(impulse, 1), size(impulse, 2))
    cube_mtf[1, :, :] .= impulse
    cube_mtf[2, :, :] .= impulse
    scratch_mtf = similar(cube_mtf)
    stack_mtf = AdaptiveOpticsSim.capture_stack!(det_mtf, cube_mtf, scratch_mtf; rng=MersenneTwister(10))
    @test size(stack_mtf) == size(cube_mtf)
    @test all(isfinite, stack_mtf)
    cube_sampled = Array{Float64}(undef, 2, size(impulse, 1), size(impulse, 2))
    cube_sampled[1, :, :] .= impulse
    cube_sampled[2, :, :] .= impulse
    scratch_sampled = similar(cube_sampled)
    stack_sampled = AdaptiveOpticsSim.capture_stack!(sampled_det, cube_sampled, scratch_sampled; rng=MersenneTwister(10))
    @test size(stack_sampled) == size(cube_sampled)
    @test all(isfinite, stack_sampled)
    @test stack_sampled[1, :, :] ≈ sampled_frame atol=1e-6
    @test stack_sampled[2, :, :] ≈ sampled_frame atol=1e-6

    rect_det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=RectangularPixelAperture(pitch_x_px=2.0, pitch_y_px=2.0,
            fill_factor_x=0.6, fill_factor_y=0.8))
    rect_frame = capture!(rect_det, impulse; rng=MersenneTwister(4))
    @test sum(rect_frame) ≈ 1.0 atol=1e-6
    @test rect_frame[5, 5] < 1.0
    rect_meta = detector_export_metadata(rect_det)
    @test rect_meta.frame_response == :rectangular_aperture
    @test rect_meta.pitch_x_px == 2.0
    @test rect_meta.pitch_y_px == 2.0
    @test rect_meta.fill_factor_x == 0.6
    @test rect_meta.fill_factor_y == 0.8
    @test rect_meta.aperture_shape == :rectangular
    @test rect_meta.response_application_domain == :image
    @test supports_detector_mtf(rect_det)

    mtf_det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        response_model=SeparablePixelMTF(pitch_x_px=1.0, pitch_y_px=1.0,
            fill_factor_x=0.7, fill_factor_y=0.7))
    mtf_frame = capture!(mtf_det, impulse; rng=MersenneTwister(5))
    @test sum(mtf_frame) ≈ 1.0 atol=1e-6
    mtf_meta2 = detector_export_metadata(mtf_det)
    @test mtf_meta2.frame_response == :separable_mtf
    @test mtf_meta2.response_is_separable
    @test mtf_meta2.response_application_domain == :image
    @test mtf_meta2.aperture_shape == :rectangular
    @test_throws InvalidConfiguration RectangularPixelAperture(fill_factor_x=0.0)
    @test_throws InvalidConfiguration SeparablePixelMTF(fill_factor_y=1.5)

    det_window_stack = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        readout_window=FrameWindow(2:8, 2:8))
    cube_window = Array{Float64}(undef, 2, size(impulse, 1), size(impulse, 2))
    cube_window[1, :, :] .= impulse
    cube_window[2, :, :] .= impulse
    scratch_window = similar(cube_window)
    @test_throws InvalidConfiguration AdaptiveOpticsSim.capture_stack!(det_window_stack, cube_window, scratch_window; rng=MersenneTwister(10))
    input_window_stack = copy(cube_window)
    output_window_stack = Array{Float64}(undef, 2, 7, 7)
    generalized_window = AdaptiveOpticsSim.capture_stack!(det_window_stack, output_window_stack, input_window_stack; rng=MersenneTwister(10))
    @test size(generalized_window) == (2, 7, 7)
    @test generalized_window[1, :, :] ≈ capture!(det_window_stack, impulse; rng=MersenneTwister(10))

    corrected_stack_det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=HgCdTeAvalancheArraySensor(sampling_mode=SingleRead()),
        correction_model=ReferencePixelCommonModeCorrection(1, 1))
    corrected_stack_in = Array{Float64}(undef, 2, 5, 5)
    corrected_stack_in[1, :, :] .= reshape(collect(1.0:25.0), 5, 5)
    corrected_stack_in[2, :, :] .= reshape(collect(26.0:50.0), 5, 5)
    corrected_stack_ref = copy(corrected_stack_in)
    corrected_stack = AdaptiveOpticsSim.capture_stack!(corrected_stack_det, corrected_stack_in,
        similar(corrected_stack_in); rng=MersenneTwister(10))
    @test size(corrected_stack) == size(corrected_stack_in)
    corrected_frame_1 = capture!(Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
            sensor=HgCdTeAvalancheArraySensor(sampling_mode=SingleRead()),
            correction_model=ReferencePixelCommonModeCorrection(1, 1)),
        @view(corrected_stack_ref[1, :, :]); rng=MersenneTwister(10))
    corrected_frame_2 = capture!(Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
            sensor=HgCdTeAvalancheArraySensor(sampling_mode=SingleRead()),
            correction_model=ReferencePixelCommonModeCorrection(1, 1)),
        @view(corrected_stack_ref[2, :, :]); rng=MersenneTwister(10))
    @test corrected_stack[1, :, :] ≈ corrected_frame_1
    @test corrected_stack[2, :, :] ≈ corrected_frame_2

    det_cmos_batched = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=CMOSSensor(column_readout_sigma=1.0))
    @test_throws InvalidConfiguration AdaptiveOpticsSim.capture_stack!(det_cmos_batched, cube_mtf, scratch_mtf; rng=MersenneTwister(10))

    det_generalized = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, psf_sampling=2, binning=2,
        bits=8, full_well=10.0, output_precision=UInt8)
    input_generalized = zeros(Float64, 2, 8, 8)
    input_generalized[1, 4, 4] = 10.0
    input_generalized[2, 5, 5] = 10.0
    output_generalized = Array{UInt8}(undef, 2, 2, 2)
    generalized_stack = AdaptiveOpticsSim.capture_stack!(det_generalized, output_generalized, input_generalized; rng=MersenneTwister(10))
    @test size(generalized_stack) == (2, 2, 2)
    @test generalized_stack[1, :, :] == capture!(det_generalized, @view(input_generalized[1, :, :]); rng=MersenneTwister(10))
    @test generalized_stack[2, :, :] == capture!(det_generalized, @view(input_generalized[2, :, :]); rng=MersenneTwister(10))

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

    det_buffered = Detector(integration_time=2.0, noise=NoiseNone(), qe=1.0, binning=1)
    frame_partial = copy(capture!(det_buffered, fill(1.0, 4, 4); rng=MersenneTwister(2), sample_time=1.0))
    @test !readout_ready(det_buffered)
    @test sum(frame_partial) == 16.0
    frame_buffered = copy(capture!(det_buffered, fill(1.0, 4, 4); rng=MersenneTwister(2), sample_time=1.0))
    @test readout_ready(det_buffered)
    @test sum(frame_buffered) == 32.0
    reset_integration!(det_buffered)
    @test readout_ready(det_buffered)
    @test det_buffered.state.integrated_time == 0.0
    metadata_buffered = detector_export_metadata(det_buffered)
    @test metadata_buffered.output_size == size(output_frame(det_buffered))
    @test metadata_buffered.psf_sampling == 1
    @test metadata_buffered.binning == 1

    det_background_flux = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        background_flux=2.0)
    frame_background_flux = capture!(det_background_flux, zeros(4, 4); rng=MersenneTwister(2))
    @test sum(frame_background_flux) > 0

    det_background_map = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        background_map=fill(1.0, 4, 4))
    frame_background_map = capture!(det_background_map, zeros(4, 4); rng=MersenneTwister(2))
    @test frame_background_map == fill(-1.0, 4, 4)

    cube = Array{Float64}(undef, 2, 4, 4)
    cube[1, :, :] .= fill(1.0, 4, 4)
    cube[2, :, :] .= fill(2.0, 4, 4)
    scratch = similar(cube)
    det_stack = Detector(integration_time=1.0, noise=NoiseNone(), qe=0.5, binning=1)
    AdaptiveOpticsSim.capture_stack!(det_stack, cube, scratch; rng=MersenneTwister(10))
    @test cube[1, :, :] ≈ fill(0.5, 4, 4)
    @test cube[2, :, :] ≈ fill(1.0, 4, 4)

    psf = reshape(Float64.(1:256), 16, 16)
    det_fused = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, psf_sampling=2, binning=2)
    frame_fused = copy(AdaptiveOpticsSim.fill_frame!(det_fused, psf, 1.0))
    manual_mid = zeros(Float64, 8, 8)
    manual_out = zeros(Float64, 4, 4)
    AdaptiveOpticsSim.bin2d!(manual_mid, psf, 2)
    AdaptiveOpticsSim.bin2d!(manual_out, manual_mid, 2)
    @test frame_fused == manual_out
end

@testset "Shack-Hartmann valid subaperture policies" begin
    tel = Telescope(resolution=352, diameter=1.22, sampling_time=1 / 500)
    sh_geom = ShackHartmann(tel; n_subap=16, mode=Diffractive(), T=Float32)
    sh_flux = ShackHartmann(tel;
        n_subap=16,
        mode=Diffractive(),
        valid_subaperture_policy=FluxThresholdValidSubapertures(light_ratio=0.5, T=Float32),
        T=Float32)

    geom_mask = copy(sh_geom.state.valid_mask_host)
    flux_mask = copy(sh_flux.state.valid_mask_host)

    @test sum(geom_mask) == 216
    @test sum(flux_mask) == 208
    @test flux_mask != geom_mask

    missing = sort!([[idx.I[1] - 1, idx.I[2] - 1] for idx in findall(geom_mask .& .!flux_mask)])
    @test missing == [[0, 4], [0, 11], [4, 0], [4, 15], [11, 0], [11, 15], [15, 4], [15, 11]]
    @test isempty(findall(flux_mask .& .!geom_mask))
end

@testset "Asterism PSF" begin
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src1 = Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0))
    src2 = Source(band=:I, magnitude=0.0, coordinates=(1.0, 90.0))
    @test coordinates_xy_arcsec(src1) == (0.0, 0.0)
    @test coordinates_xy_arcsec(src2)[1] ≈ 0.0 atol=1e-12
    @test coordinates_xy_arcsec(src2)[2] ≈ 1.0
    ast = Asterism([src1, src2])
    psf = compute_psf!(tel, ast; zero_padding=2)
    @test size(tel.state.psf_stack, 3) == 2
    @test size(psf) == (32, 32)
    @test sum(psf) >= sum(@view tel.state.psf_stack[:, :, 1])
end

@testset "Polychromatic WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    λ0 = wavelength(src)
    bundle_single = SpectralBundle([SpectralSample(λ0, 1.0)])
    bundle_broad = SpectralBundle([SpectralSample(0.9 * λ0, 0.4), SpectralSample(1.1 * λ0, 0.6)])
    poly_single = with_spectrum(src, bundle_single)
    poly_broad = with_spectrum(src, bundle_broad)

    @test sum(sample.weight for sample in bundle_broad) ≈ 1.0
    @test weighted_wavelength(bundle_broad) ≈ (0.9 * λ0 * 0.4 + 1.1 * λ0 * 0.6)
    @test has_spectral_bundle(poly_broad)
    @test is_polychromatic(poly_broad)
    @test !is_polychromatic(poly_single)
    @test spectral_reference_source(poly_broad) === src

    zb = ZernikeBasis(tel, 5)
    compute_zernike!(zb, tel)
    focus = @view zb.modes[:, :, 5]
    @. tel.state.opd = 5e-8 * focus

    sh_mono = ShackHartmann(tel; n_subap=8, mode=Diffractive())
    sh_single = ShackHartmann(tel; n_subap=8, mode=Diffractive())
    sh_broad = ShackHartmann(tel; n_subap=8, mode=Diffractive())

    mono_slopes = copy(measure!(sh_mono, tel, src))
    single_slopes = copy(measure!(sh_single, tel, poly_single))
    broad_slopes_1 = copy(measure!(sh_broad, tel, poly_broad))
    broad_slopes_2 = copy(measure!(sh_broad, tel, poly_broad))

    @test single_slopes ≈ mono_slopes atol=1e-10 rtol=1e-10
    @test broad_slopes_1 ≈ broad_slopes_2 atol=1e-10 rtol=1e-10
    @test norm(broad_slopes_1 - mono_slopes) > 1e-8
    @test supports_stacked_sources(sh_broad, poly_broad)

    pyr_mono = PyramidWFS(tel; n_subap=8, mode=Diffractive(), modulation=1.0)
    pyr_single = PyramidWFS(tel; n_subap=8, mode=Diffractive(), modulation=1.0)
    pyr_broad = PyramidWFS(tel; n_subap=8, mode=Diffractive(), modulation=1.0)

    mono_pyr = copy(measure!(pyr_mono, tel, src))
    single_pyr = copy(measure!(pyr_single, tel, poly_single))
    broad_pyr_1 = copy(measure!(pyr_broad, tel, poly_broad))
    broad_pyr_2 = copy(measure!(pyr_broad, tel, poly_broad))

    @test single_pyr ≈ mono_pyr atol=1e-10 rtol=1e-10
    @test broad_pyr_1 ≈ broad_pyr_2 atol=1e-10 rtol=1e-10
    @test norm(broad_pyr_1 - mono_pyr) > 1e-8
    @test supports_stacked_sources(pyr_broad, poly_broad)

    det = Detector(noise=NoiseNone(), binning=1)
    spectral_frame = measure!(sh_broad, tel, poly_broad, det)
    @test size(sh_broad.state.detector_noise_cube) == size(sh_broad.state.spot_cube)
    @test spectral_frame ≈ broad_slopes_1 atol=1e-10 rtol=1e-10

    pyr_det = Detector(noise=NoiseNone(), binning=1)
    pyr_det_slopes = measure!(pyr_broad, tel, poly_broad, pyr_det)
    @test size(output_frame(pyr_det)) == size(pyr_broad.state.camera_frame)
    @test pyr_det_slopes ≈ broad_pyr_1 atol=1e-10 rtol=1e-10

    fill!(tel.state.opd, 0.0)
end

@testset "Extended-source WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    point_model = PointCloudSourceModel([(0.0, 0.0)], [1.0])
    gaussian_model = GaussianDiskSourceModel(sigma_arcsec=0.35, n_side=5)
    image_model = SampledImageSourceModel([0.0 1.0 0.0; 1.0 2.0 1.0; 0.0 1.0 0.0], pixel_scale_arcsec=0.2)
    ext_point = with_extended_source(src, point_model)
    ext_gauss = with_extended_source(src, gaussian_model)
    ext_image = with_extended_source(src, image_model)

    point_ast = extended_source_asterism(ext_point)
    image_ast = extended_source_asterism(ext_image)
    @test has_extended_source_model(ext_gauss)
    @test !has_extended_source_model(src)
    @test length(point_ast) == 1
    @test length(image_ast) == 5
    @test AdaptiveOpticsSim.photon_flux(point_ast.sources[1]) ≈ AdaptiveOpticsSim.photon_flux(src)
    @test sum(AdaptiveOpticsSim.photon_flux(sample) for sample in image_ast.sources) ≈ AdaptiveOpticsSim.photon_flux(src)

    zb = ZernikeBasis(tel, 5)
    compute_zernike!(zb, tel)
    focus = @view zb.modes[:, :, 5]
    @. tel.state.opd = 5e-8 * focus

    sh_point = ShackHartmann(tel; n_subap=8, mode=Diffractive())
    sh_ext_point = ShackHartmann(tel; n_subap=8, mode=Diffractive())
    sh_ext = ShackHartmann(tel; n_subap=8, mode=Diffractive())
    point_slopes = copy(measure!(sh_point, tel, src))
    ext_point_slopes = copy(measure!(sh_ext_point, tel, ext_point))
    point_peak = AdaptiveOpticsSim.sampled_spots_peak!(sh_point, tel, src)
    ext_peak = AdaptiveOpticsSim.sampled_spots_peak!(sh_ext, tel, ext_gauss)
    ext_slopes_1 = copy(measure!(sh_ext, tel, ext_gauss))
    ext_slopes_2 = copy(measure!(sh_ext, tel, ext_gauss))

    @test ext_point_slopes ≈ point_slopes atol=1e-10 rtol=1e-10
    @test ext_slopes_1 ≈ ext_slopes_2 atol=1e-10 rtol=1e-10
    @test ext_peak < point_peak
    @test norm(sh_ext.state.spot_cube - sh_point.state.spot_cube) > 1e-8
    @test supports_stacked_sources(sh_ext, ext_gauss)

    pyr_point = PyramidWFS(tel; n_subap=8, mode=Diffractive(), modulation=1.0)
    pyr_ext_point = PyramidWFS(tel; n_subap=8, mode=Diffractive(), modulation=1.0)
    pyr_ext = PyramidWFS(tel; n_subap=8, mode=Diffractive(), modulation=1.0)
    pyr_point_slopes = copy(measure!(pyr_point, tel, src))
    pyr_ext_point_slopes = copy(measure!(pyr_ext_point, tel, ext_point))
    pyr_ext_slopes_1 = copy(measure!(pyr_ext, tel, ext_gauss))
    pyr_ext_slopes_2 = copy(measure!(pyr_ext, tel, ext_gauss))

    @test pyr_ext_point_slopes ≈ pyr_point_slopes atol=1e-10 rtol=1e-10
    @test pyr_ext_slopes_1 ≈ pyr_ext_slopes_2 atol=1e-10 rtol=1e-10
    @test norm(pyr_ext.state.intensity - pyr_point.state.intensity) > 1e-10
    @test supports_stacked_sources(pyr_ext, ext_gauss)

    det = Detector(noise=NoiseNone(), binning=1)
    sh_det_slopes = measure!(sh_ext, tel, ext_gauss, det)
    @test size(sh_ext.state.detector_noise_cube) == size(sh_ext.state.spot_cube)
    @test sh_det_slopes ≈ ext_slopes_1 atol=1e-10 rtol=1e-10

    pyr_det = Detector(noise=NoiseNone(), binning=1)
    pyr_det_slopes = measure!(pyr_ext, tel, ext_gauss, pyr_det)
    @test size(output_frame(pyr_det)) == size(pyr_ext.state.camera_frame)
    @test pyr_det_slopes ≈ pyr_ext_slopes_1 atol=1e-10 rtol=1e-10

    fill!(tel.state.opd, 0.0)
end

@testset "Pupil masks and misregistration" begin
    tel = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    base_sum = sum(tel.state.pupil)
    apply_spiders!(tel; thickness=0.5, angles=[0.0, 90.0])
    @test sum(tel.state.pupil) < base_sum

    custom = trues(16, 16)
    custom[:, 9:end] .= false
    set_pupil!(tel, custom)
    @test sum(tel.state.pupil) == sum(custom)
    @test tel.state.pupil_reflectivity == Float64.(custom)

    reflectivity = fill(0.5, 16, 16)
    set_pupil_reflectivity!(tel, reflectivity)
    @test tel.state.pupil_reflectivity[:, 1:8] == fill(0.5, 16, 8)
    @test tel.state.pupil_reflectivity[:, 9:end] == fill(0.0, 16, 8)

    tel2 = Telescope(resolution=16, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    dm1 = DeformableMirror(tel2; n_act=2, influence_width=0.3)
    mis = Misregistration(shift_x=0.1, shift_y=0.0, rotation_deg=5.0, T=Float64)
    @test rotation_deg(mis) ≈ 5.0
    @test rotation_rad(mis) ≈ deg2rad(5.0)
    dm2 = DeformableMirror(tel2; n_act=2, influence_width=0.3, misregistration=mis)
    @test dm1.state.modes != dm2.state.modes
end

@testset "Pyramid, BioEdge, and LGS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        tel.state.opd[i, j] = i
    end

    pyr = PyramidWFS(tel; n_subap=4, modulation=1.0)
    pyr_slopes = measure!(pyr, tel)
    @test length(pyr_slopes) == 2 * 4 * 4

    bio = BioEdgeWFS(tel; n_subap=4)
    bio_slopes = measure!(bio, tel)
    @test length(bio_slopes) == 2 * 4 * 4

    sh = ShackHartmann(tel; n_subap=4)
    ngs = Source(band=:I, magnitude=0.0)
    lgs = LGSSource(elongation_factor=2.0)
    slopes_ngs = measure!(sh, tel, ngs)
    slopes_lgs = measure!(sh, tel, lgs)
    n = sh.params.n_subap * sh.params.n_subap
    @test slopes_lgs[n+1:end] ≈ slopes_ngs[n+1:end] .* 2.0
end

@testset "Diffractive WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    for i in 1:tel.params.resolution, j in 1:tel.params.resolution
        tel.state.opd[i, j] = i
    end
    ngs = Source(band=:I, magnitude=0.0)
    lgs = LGSSource(elongation_factor=1.5)

    sh = ShackHartmann(tel; n_subap=4, mode=Diffractive())
    @test_throws InvalidConfiguration measure!(sh, tel)
    sh_slopes = measure!(sh, tel, ngs)
    @test length(sh_slopes) == 2 * 4 * 4
    @test all(isfinite, sh_slopes)
    sh_lgs = measure!(sh, tel, lgs)
    @test all(isfinite, sh_lgs)

    na_profile = [80000.0 90000.0 100000.0; 0.2 0.6 0.2]
    lgs_profile = LGSSource(elongation_factor=1.2, na_profile=na_profile, fwhm_spot_up=1.0)
    sh_profile = ShackHartmann(tel; n_subap=4, mode=Diffractive())
    sh_profile_slopes = measure!(sh_profile, tel, lgs_profile)
    @test all(isfinite, sh_profile_slopes)

    sh_sampled = ShackHartmann(tel; n_subap=4, mode=Diffractive(), pixel_scale=0.06, n_pix_subap=8)
    sh_sampled_slopes = measure!(sh_sampled, tel, ngs)
    @test length(sh_sampled_slopes) == 2 * 4 * 4

    pyr_sampled = PyramidWFS(tel; n_subap=4, mode=Diffractive(), n_pix_separation=4, binning=2)
    pyr_sampled_slopes = measure!(pyr_sampled, tel, ngs)
    @test length(pyr_sampled_slopes) == 2 * count(pyr_sampled.state.valid_i4q)
    pyr_intensity = reshape(Float64.(1:size(tel.state.opd, 1)^2), size(tel.state.opd))
    pyr_frame = copy(AdaptiveOpticsSim.sample_pyramid_intensity!(pyr_sampled, tel, pyr_intensity))
    pyr_camera = zeros(Float64, 4, 4)
    pyr_manual = zeros(Float64, 2, 2)
    AdaptiveOpticsSim.bin2d!(pyr_camera, pyr_intensity, 8)
    AdaptiveOpticsSim.bin2d!(pyr_manual, pyr_camera, 2)
    @test pyr_frame == pyr_manual

    bio_sampled = BioEdgeWFS(tel; n_subap=4, mode=Diffractive(), binning=2)
    bio_sampled_slopes = measure!(bio_sampled, tel, ngs)
    @test length(bio_sampled_slopes) == 2 * count(bio_sampled.state.valid_i4q)
    bio_intensity = reshape(Float64.(1:size(tel.state.opd, 1)^2), size(tel.state.opd))
    bio_frame = copy(AdaptiveOpticsSim.sample_bioedge_intensity!(bio_sampled, tel, bio_intensity))
    bio_camera = zeros(Float64, 4, 4)
    bio_manual = similar(bio_frame)
    AdaptiveOpticsSim.bin2d!(bio_camera, bio_intensity, 8)
    AdaptiveOpticsSim.bin2d!(bio_manual, bio_camera, div(size(bio_camera, 1), size(bio_frame, 1)))
    @test bio_frame == bio_manual

    pyr_profile = PyramidWFS(tel; n_subap=4, mode=Diffractive())
    pyr_profile_slopes = measure!(pyr_profile, tel, lgs_profile)
    @test all(isfinite, pyr_profile_slopes)

    bio_profile = BioEdgeWFS(tel; n_subap=4, mode=Diffractive())
    bio_profile_slopes = measure!(bio_profile, tel, lgs_profile)
    @test all(isfinite, bio_profile_slopes)

    pyr = PyramidWFS(tel; n_subap=4, mode=Diffractive())
    @test_throws InvalidConfiguration measure!(pyr, tel)
    pyr_slopes = measure!(pyr, tel, ngs)
    @test length(pyr_slopes) == 2 * 4 * 4

    bio = BioEdgeWFS(tel; n_subap=4, mode=Diffractive())
    @test_throws InvalidConfiguration measure!(bio, tel)
    bio_slopes = measure!(bio, tel, ngs)
    @test length(bio_slopes) == 2 * 4 * 4

    det = Detector(noise=NoiseNone(), binning=1)
    sh_det = ShackHartmann(tel; n_subap=4, mode=Diffractive())
    sh_det_slopes = measure!(sh_det, tel, ngs, det)
    @test length(sh_det_slopes) == 2 * 4 * 4
    pyr_det = PyramidWFS(tel; n_subap=4, mode=Diffractive())
    pyr_det_slopes = measure!(pyr_det, tel, ngs, det)
    @test length(pyr_det_slopes) == 2 * 4 * 4

    ast = Asterism([ngs, Source(band=:I, magnitude=0.0, coordinates=(0.0, 0.0))])
    sh_ast = ShackHartmann(tel; n_subap=4, mode=Diffractive())
    sh_ast_slopes = copy(measure!(sh_ast, tel, ast))
    @test length(sh_ast_slopes) == 2 * 4 * 4
    sh_ast_serial = ShackHartmann(tel; n_subap=4, mode=Diffractive())
    AdaptiveOpticsSim.prepare_sampling!(sh_ast_serial, tel, ast.sources[1])
    AdaptiveOpticsSim.ensure_sh_calibration!(sh_ast_serial, tel, ast.sources[1])
    fill!(sh_ast_serial.state.detector_noise_cube, zero(eltype(sh_ast_serial.state.detector_noise_cube)))
    for src in ast.sources
        AdaptiveOpticsSim.sampled_spots_peak!(sh_ast_serial, tel, src)
        sh_ast_serial.state.detector_noise_cube .+= sh_ast_serial.state.spot_cube
    end
    copyto!(sh_ast_serial.state.spot_cube, sh_ast_serial.state.detector_noise_cube)
    sh_ast_serial_peak = maximum(sh_ast_serial.state.spot_cube)
    AdaptiveOpticsSim.sh_signal_from_spots!(sh_ast_serial, sh_ast_serial_peak, slope_extraction_model(sh_ast_serial))
    AdaptiveOpticsSim.subtract_reference_and_scale!(sh_ast_serial)
    sh_ast_serial_slopes = copy(sh_ast_serial.state.slopes)
    @test sh_ast_slopes ≈ sh_ast_serial_slopes
    mixed_ngs = Source(wavelength=wavelength(lgs_profile), magnitude=0.0, coordinates=(0.0, 0.0))
    mixed_ast = Asterism([mixed_ngs, lgs_profile])
    sh_mixed_det = ShackHartmann(tel; n_subap=4, mode=Diffractive())
    sh_mixed_det_slopes = copy(measure!(sh_mixed_det, tel, mixed_ast, det; rng=MersenneTwister(14)))
    sh_mixed_det_frame = copy(sh_mixed_det.state.spot_cube)
    sh_mixed_det_manual = ShackHartmann(tel; n_subap=4, mode=Diffractive())
    AdaptiveOpticsSim.prepare_sampling!(sh_mixed_det_manual, tel, mixed_ast.sources[1])
    AdaptiveOpticsSim.ensure_sh_calibration!(sh_mixed_det_manual, tel, mixed_ast.sources[1])
    fill!(sh_mixed_det_manual.state.detector_noise_cube, zero(eltype(sh_mixed_det_manual.state.detector_noise_cube)))
    for src in mixed_ast.sources
        AdaptiveOpticsSim.sampled_spots_peak!(sh_mixed_det_manual, tel, src, det, MersenneTwister(14))
        sh_mixed_det_manual.state.detector_noise_cube .+= sh_mixed_det_manual.state.spot_cube
    end
    copyto!(sh_mixed_det_manual.state.spot_cube, sh_mixed_det_manual.state.detector_noise_cube)
    @test sh_mixed_det_frame ≈ sh_mixed_det_manual.state.spot_cube
    @test length(sh_mixed_det_slopes) == 2 * 4 * 4
    pyr_ast = PyramidWFS(tel; n_subap=4, mode=Diffractive())
    pyr_ast_slopes = copy(measure!(pyr_ast, tel, ast))
    @test length(pyr_ast_slopes) == 2 * 4 * 4
    pyr_ast_serial = PyramidWFS(tel; n_subap=4, mode=Diffractive())
    AdaptiveOpticsSim.ensure_pyramid_calibration!(pyr_ast_serial, tel, ast.sources[1])
    pyr_ast_stack = @view AdaptiveOpticsSim.ensure_pyramid_asterism_stack!(pyr_ast_serial, length(ast.sources))[:, :, 1:length(ast.sources)]
    fill!(pyr_ast_serial.state.intensity, zero(eltype(pyr_ast_serial.state.intensity)))
    for (src_idx, src) in pairs(ast.sources)
        AdaptiveOpticsSim.pyramid_intensity!(@view(pyr_ast_stack[:, :, src_idx]), pyr_ast_serial, tel, src)
        pyr_ast_serial.state.intensity .+= @view(pyr_ast_stack[:, :, src_idx])
    end
    pyr_ast_intensity = AdaptiveOpticsSim.sample_pyramid_intensity!(pyr_ast_serial, tel, pyr_ast_serial.state.intensity)
    AdaptiveOpticsSim.pyramid_signal!(pyr_ast_serial, tel, pyr_ast_intensity)
    @. pyr_ast_serial.state.slopes *= pyr_ast_serial.state.optical_gain
    @test pyr_ast_slopes ≈ pyr_ast_serial.state.slopes
    bio_ast = BioEdgeWFS(tel; n_subap=4, mode=Diffractive())
    bio_ast_slopes = copy(measure!(bio_ast, tel, ast))
    @test length(bio_ast_slopes) == 2 * 4 * 4
    bio_ast_serial = BioEdgeWFS(tel; n_subap=4, mode=Diffractive())
    AdaptiveOpticsSim.ensure_bioedge_calibration!(bio_ast_serial, tel, ast.sources[1])
    fill!(bio_ast_serial.state.binned_intensity, zero(eltype(bio_ast_serial.state.binned_intensity)))
    for src in ast.sources
        AdaptiveOpticsSim.bioedge_intensity!(bio_ast_serial.state.intensity, bio_ast_serial, tel, src)
        bio_ast_serial.state.binned_intensity .+= bio_ast_serial.state.intensity
    end
    bio_ast_intensity = AdaptiveOpticsSim.sample_bioedge_intensity!(bio_ast_serial, tel, bio_ast_serial.state.binned_intensity)
    AdaptiveOpticsSim.bioedge_signal!(bio_ast_serial, tel, bio_ast_intensity)
    @test bio_ast_slopes ≈ bio_ast_serial.state.slopes

    pyr_ast_det = PyramidWFS(tel; n_subap=4, mode=Diffractive())
    pyr_ast_det_slopes = copy(measure!(pyr_ast_det, tel, ast, det))
    pyr_ast_det_serial = PyramidWFS(tel; n_subap=4, mode=Diffractive())
    AdaptiveOpticsSim.ensure_pyramid_calibration!(pyr_ast_det_serial, tel, ast.sources[1])
    pyr_ast_det_stack = @view AdaptiveOpticsSim.ensure_pyramid_asterism_stack!(pyr_ast_det_serial, length(ast.sources))[:, :, 1:length(ast.sources)]
    fill!(pyr_ast_det_serial.state.intensity, zero(eltype(pyr_ast_det_serial.state.intensity)))
    for (src_idx, src) in pairs(ast.sources)
        AdaptiveOpticsSim.pyramid_intensity!(@view(pyr_ast_det_stack[:, :, src_idx]), pyr_ast_det_serial, tel, src)
        pyr_ast_det_serial.state.intensity .+= @view(pyr_ast_det_stack[:, :, src_idx])
    end
    pyr_ast_det_intensity = AdaptiveOpticsSim.sample_pyramid_intensity!(pyr_ast_det_serial, tel, pyr_ast_det_serial.state.intensity)
    pyr_ast_det_frame = capture!(det, pyr_ast_det_intensity; rng=MersenneTwister(12))
    AdaptiveOpticsSim.resize_pyramid_signal_buffers!(pyr_ast_det_serial, size(pyr_ast_det_frame, 1))
    AdaptiveOpticsSim.pyramid_signal!(pyr_ast_det_serial, tel, pyr_ast_det_frame)
    @. pyr_ast_det_serial.state.slopes *= pyr_ast_det_serial.state.optical_gain
    @test pyr_ast_det_slopes ≈ pyr_ast_det_serial.state.slopes

    bio_ast_det = BioEdgeWFS(tel; n_subap=4, mode=Diffractive())
    bio_ast_det_slopes = copy(measure!(bio_ast_det, tel, ast, det; rng=MersenneTwister(13)))
    bio_ast_det_serial = BioEdgeWFS(tel; n_subap=4, mode=Diffractive())
    AdaptiveOpticsSim.ensure_bioedge_calibration!(bio_ast_det_serial, tel, ast.sources[1])
    fill!(bio_ast_det_serial.state.binned_intensity, zero(eltype(bio_ast_det_serial.state.binned_intensity)))
    for src in ast.sources
        AdaptiveOpticsSim.bioedge_intensity!(bio_ast_det_serial.state.intensity, bio_ast_det_serial, tel, src)
        bio_ast_det_serial.state.binned_intensity .+= bio_ast_det_serial.state.intensity
    end
    bio_ast_det_intensity = AdaptiveOpticsSim.sample_bioedge_intensity!(bio_ast_det_serial, tel, bio_ast_det_serial.state.binned_intensity)
    bio_ast_det_frame = capture!(det, bio_ast_det_intensity; rng=MersenneTwister(13))
    AdaptiveOpticsSim.resize_bioedge_signal_buffers!(bio_ast_det_serial, size(bio_ast_det_frame, 1))
    AdaptiveOpticsSim.bioedge_signal!(bio_ast_det_serial, tel, bio_ast_det_frame)
    @test bio_ast_det_slopes ≈ bio_ast_det_serial.state.slopes
end

@testset "Shack-Hartmann subapertures" begin
    tel = Telescope(resolution=24, diameter=8.0, sampling_time=1e-3, central_obstruction=0.1)
    src = Source(band=:I, magnitude=0.0)
    sh = ShackHartmann(tel; n_subap=6, mode=Diffractive(), pixel_scale=0.06, n_pix_subap=8, threshold_cog=0.02)

    layout = subaperture_layout(sh)
    calibration = subaperture_calibration(sh)
    @test layout isa SubapertureLayout
    @test calibration isa SubapertureCalibration
    @test layout.n_subap == 6
    @test layout.subap_pixels == 4
    @test layout.pitch_m ≈ tel.params.diameter / 6
    @test !calibration.calibrated
    @test slope_extraction_model(sh) isa CenterOfGravityExtraction
    @test slope_extraction_model(sh).threshold ≈ 0.02
    @test n_valid_subapertures(layout) == count(layout.valid_mask_host)

    prepare_runtime_wfs!(sh, tel, src)
    @test calibration.calibrated
    @test calibration.slopes_units == sh.state.slopes_units
    @test calibration.wavelength == sh.state.calibration_wavelength
    @test calibration.signature == sh.state.calibration_signature
    @test calibration.reference_signal_2d === sh.state.reference_signal_2d
    @test calibration.reference_signal_host === sh.state.reference_signal_host
    @test length(valid_subaperture_indices(layout)) == n_valid_subapertures(layout)

    slopes = measure!(sh, tel, src)
    @test all(isfinite, slopes)
    meta = AdaptiveOpticsSim.wfs_output_metadata(sh)
    @test meta.n_valid_subap == n_valid_subapertures(layout)
    @test meta.subap_pixels == layout.subap_pixels
    @test meta.calibrated

    dm = DeformableMirror(tel; n_act=5)
    imat = interaction_matrix(dm, sh, tel, src; amplitude=1e-8)
    @test size(imat.matrix, 1) == length(sh.state.slopes)
    @test size(imat.matrix, 2) == length(dm.state.coefs)
end

@testset "Zernike WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    wfs = ZernikeWFS(tel; n_subap=8, diffraction_padding=2)

    @test size(wfs.state.camera_frame) == (8, 8)
    @test length(wfs.state.slopes) == count(wfs.state.valid_mask)
    @test_throws InvalidConfiguration measure!(wfs, tel)
    @test_throws InvalidConfiguration measure!(wfs, tel, Asterism([src, Source(band=:I, magnitude=0.0)]))

    flat_slopes = copy(measure!(wfs, tel, src))
    @test wfs.state.calibrated
    @test all(isfinite, flat_slopes)
    @test all(>=(0.0), wfs.state.camera_frame)
    @test flat_slopes ≈ zero.(flat_slopes) atol=1e-10

    det = Detector(noise=NoiseNone(), binning=1)
    det_slopes = copy(measure!(wfs, tel, src, det))
    @test det_slopes ≈ flat_slopes atol=1e-10
    @test size(output_frame(det)) == size(wfs.state.camera_frame)

    zb = ZernikeBasis(tel, 5)
    compute_zernike!(zb, tel)
    focus = @view zb.modes[:, :, 5]
    @. tel.state.opd = 5e-8 * focus
    slopes_plus = copy(measure!(wfs, tel, src))
    @. tel.state.opd = -5e-8 * focus
    slopes_minus = copy(measure!(wfs, tel, src))
    fill!(tel.state.opd, 0.0)

    @test norm(slopes_plus) > 1e-6
    @test norm(slopes_minus) > 1e-6
    @test dot(slopes_plus, slopes_minus) < 0
end

@testset "Curvature WFS" begin
    tel = Telescope(resolution=32, diameter=8.0, sampling_time=1e-3, central_obstruction=0.0)
    src = Source(band=:I, magnitude=0.0)
    wfs = CurvatureWFS(tel; n_subap=8, defocus_rms_nm=500.0)

    @test size(wfs.state.camera_frame) == (16, 8)
    @test length(wfs.state.slopes) == 64
    @test_throws InvalidConfiguration measure!(wfs, tel)
    @test_throws InvalidConfiguration measure!(wfs, tel, Asterism([src, Source(band=:I, magnitude=0.0)]))

    flat_slopes = copy(measure!(wfs, tel, src))
    @test wfs.state.calibrated
    @test all(isfinite, flat_slopes)
    @test all(>=(0.0), wfs.state.camera_frame)
    @test flat_slopes ≈ zero.(flat_slopes) atol=1e-10

    det = Detector(noise=NoiseNone(), binning=1)
    det_slopes = copy(measure!(wfs, tel, src, det))
    @test det_slopes ≈ flat_slopes atol=1e-10
    @test size(output_frame(det)) == size(wfs.state.camera_frame)

    zb = ZernikeBasis(tel, 5)
    compute_zernike!(zb, tel)
    focus = @view zb.modes[:, :, 5]
    @. tel.state.opd = 5e-8 * focus
    slopes_plus = copy(measure!(wfs, tel, src))
    @. tel.state.opd = -5e-8 * focus
    slopes_minus = copy(measure!(wfs, tel, src))
    fill!(tel.state.opd, 0.0)

    @test norm(slopes_plus) > 1e-6
    @test norm(slopes_minus) > 1e-6
    @test dot(slopes_plus, slopes_minus) < 0

    counting = CurvatureWFS(tel; n_subap=8, defocus_rms_nm=500.0, readout_model=CurvatureCountingReadout())
    counting_flat = copy(measure!(counting, tel, src))
    @test size(counting.state.camera_frame) == (2, 64)
    @test counting_flat ≈ zero.(counting_flat) atol=1e-10
    @test_throws InvalidConfiguration measure!(counting, tel, src, det)
    apd = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0, noise=NoiseNone())
    counting_apd = copy(measure!(counting, tel, src, apd))
    @test counting_apd ≈ counting_flat atol=1e-10
    @test detector_export_metadata(apd).readout.output_size == size(counting.state.camera_frame)
    apd_dead = APDDetector(integration_time=1.0, qe=1.0, gain=1.0, dark_count_rate=0.0,
        noise=NoiseNone(), dead_time_model=NonParalyzableDeadTime(0.25))
    counting_dead = copy(measure!(counting, tel, src, apd_dead))
    @test counting_dead ≈ counting_flat atol=1e-10
    @test_throws InvalidConfiguration CurvatureWFS(tel; n_subap=8, readout_model=CurvatureCountingReadout(),
        readout_pixels_per_subap=2)

    response = CurvatureBranchResponse(T=Float64, plus_throughput=1.2, minus_throughput=0.8,
        plus_background=5.0, minus_background=1.0)
    imbalanced = CurvatureWFS(tel; n_subap=8, defocus_rms_nm=500.0, branch_response=response)
    imbalanced_flat = copy(measure!(imbalanced, tel, src))
    @test imbalanced_flat ≈ zero.(imbalanced_flat) atol=1e-10
    plus_mean = mean(@view imbalanced.state.camera_frame[1:imbalanced.params.n_subap, :])
    minus_mean = mean(@view imbalanced.state.camera_frame[imbalanced.params.n_subap+1:end, :])
    @test plus_mean > minus_mean
    @test_throws InvalidConfiguration CurvatureBranchResponse(plus_throughput=-1.0)

    oversampled = CurvatureWFS(tel; n_subap=8, readout_crop_resolution=16, readout_pixels_per_subap=2)
    oversampled_flat = copy(measure!(oversampled, tel, src))
    @test size(oversampled.state.camera_frame) == (32, 16)
    @test size(oversampled.state.frame_plus) == (16, 16)
    @test size(oversampled.state.reduced_plus) == (8, 8)
    @test oversampled_flat ≈ zero.(oversampled_flat) atol=1e-10
    @test_throws InvalidConfiguration CurvatureWFS(tel; n_subap=8, readout_crop_resolution=18, readout_pixels_per_subap=2)

    atm = MultiLayerAtmosphere(tel;
        r0=0.2,
        L0=25.0,
        fractional_cn2=[0.7, 0.3],
        wind_speed=[8.0, 4.0],
        wind_direction=[0.0, 90.0],
        altitude=[0.0, 5000.0],
    )
    advance!(atm, tel; rng=MersenneTwister(3))
    atm_slopes = copy(measure!(wfs, tel, src, atm))
    @test all(isfinite, atm_slopes)
    @test norm(atm_slopes) > 0

    det_atm = Detector(noise=NoiseNone(), binning=1)
    det_atm_slopes = copy(measure!(wfs, tel, src, atm, det_atm))
    @test all(isfinite, det_atm_slopes)

    ast = Asterism([src, Source(band=:I, magnitude=0.0, coordinates=(1.0, 90.0))])
    ast_slopes = copy(measure!(wfs, tel, ast, atm))
    @test length(ast_slopes) == length(wfs.state.slopes)
    @test all(isfinite, ast_slopes)
    @test norm(ast_slopes) > 0
end

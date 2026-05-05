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
    @test metadata_adc.output_type == UInt8
    @test metadata_adc.frame_size == (4, 4)
    @test metadata_adc.output_size == (4, 4)

    det_adc_float = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        bits=8, full_well=10.0, output_type=Float32)
    frame_adc_float = capture!(det_adc_float, fill(10.0, 4, 4); rng=MersenneTwister(2))
    @test frame_adc_float isa Matrix{Float32}
    @test maximum(frame_adc_float) == Float32(255.0)

    det_adc_window_corr = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        bits=8, full_well=100.0, readout_window=FrameWindow(2:5, 3:7), output_type=UInt16,
        sensor=HgCdTeAvalancheArraySensor(sampling_mode=SingleRead()),
        response_model=NullFrameResponse(),
        correction_model=CompositeFrameReadoutCorrection((
            ReferenceRowCommonModeCorrection(1),
            ReferenceColumnCommonModeCorrection(1),
        )))
    adc_window_in = reshape(collect(1.0:96.0), 2, 6, 8)
    adc_window_out = Array{UInt16}(undef, 2, 4, 5)
    generalized_adc_window = AdaptiveOpticsSim.capture_stack!(det_adc_window_corr, adc_window_out, copy(adc_window_in);
        rng=MersenneTwister(10))
    @test size(generalized_adc_window) == (2, 4, 5)
    @test generalized_adc_window[1, :, :] == capture!(Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
            bits=8, full_well=100.0, readout_window=FrameWindow(2:5, 3:7), output_type=UInt16,
            sensor=HgCdTeAvalancheArraySensor(sampling_mode=SingleRead()),
            response_model=NullFrameResponse(),
            correction_model=CompositeFrameReadoutCorrection((
                ReferenceRowCommonModeCorrection(1),
                ReferenceColumnCommonModeCorrection(1),
            ))),
        @view(adc_window_in[1, :, :]); rng=MersenneTwister(10))
    @test generalized_adc_window[2, :, :] == capture!(Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
            bits=8, full_well=100.0, readout_window=FrameWindow(2:5, 3:7), output_type=UInt16,
            sensor=HgCdTeAvalancheArraySensor(sampling_mode=SingleRead()),
            response_model=NullFrameResponse(),
            correction_model=CompositeFrameReadoutCorrection((
                ReferenceRowCommonModeCorrection(1),
                ReferenceColumnCommonModeCorrection(1),
            ))),
        @view(adc_window_in[2, :, :]); rng=MersenneTwister(10))

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
    @test_throws InvalidConfiguration EMCCDSensor(em_gain_range=(10.0, 1.0))
    @test_throws InvalidConfiguration EMCCDSensor(readout_rate_hz=-1.0)
    @test_throws InvalidConfiguration PhotonCountingEMMode(threshold=-1.0)
    @test_throws InvalidConfiguration PhotonCountingEMMode(threshold=1.0, detection_efficiency=1.5)

    det_emccd_conventional = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=10.0, sensor=EMCCDSensor(output_path=ConventionalOutput()))
    @test capture!(det_emccd_conventional, uniform_signal; rng=MersenneTwister(127)) == uniform_signal

    det_emccd_pc = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=10.0,
        sensor=EMCCDSensor(operating_mode=PhotonCountingEMMode(threshold=5.0, detection_efficiency=0.8)))
    pc_frame = capture!(det_emccd_pc, [0.0 0.4; 0.6 1.0]; rng=MersenneTwister(128))
    @test pc_frame == [0.0 0.0; 0.8 0.8]
    @test supports_photon_number_resolving(det_emccd_pc.params.sensor)

    det_emccd_pc_batched = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=10.0, response_model=NullFrameResponse(),
        sensor=EMCCDSensor(operating_mode=PhotonCountingEMMode(threshold=5.0)))
    pc_cube = reshape(Float64[0.0, 0.4, 0.6, 1.0], 1, 2, 2)
    pc_scratch = similar(pc_cube)
    capture_stack!(det_emccd_pc_batched, pc_cube, pc_scratch, MersenneTwister(129))
    @test pc_cube == reshape(Float64[0.0, 0.0, 1.0, 1.0], 1, 2, 2)

    @test emccd_snr(1.0; readout_noise=20.0, gain=100.0, excess_noise_factor=1.0) >
        emccd_snr(1.0; readout_noise=20.0, gain=100.0, excess_noise_factor=sqrt(2.0))
    @test emccd_snr(1.0; readout_noise=20.0, gain=100.0, operating_mode=PhotonCountingEMMode(threshold=0.5)) >
        emccd_snr(1.0; readout_noise=20.0, gain=100.0, excess_noise_factor=sqrt(2.0))
    @test emccd_snr(1.0; readout_noise=20.0, gain=100.0, output_path=EMOutput()) >
        emccd_snr(1.0; readout_noise=20.0, gain=100.0, output_path=ConventionalOutput())

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
    @test_throws InvalidConfiguration CMOSSensor(timing_model=RollingShutter(1e-3; row_group_size=0))
    @test_throws InvalidConfiguration CMOSSensor(output_model=StaticCMOSOutputPattern(2, [1.0], [0.0, 1.0]))

    rolling_det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=CMOSSensor(timing_model=RollingShutter(0.25)),
        response_model=NullFrameResponse())
    rolling_source = InPlaceFrameSource((out, t) -> fill!(out, t), (4, 4))
    rolling_frame = capture!(rolling_det, rolling_source; rng=MersenneTwister(127))
    @test rolling_frame == repeat(reshape([0.0, 0.25, 0.5, 0.75], :, 1), 1, 4)
    @test detector_export_metadata(rolling_det).sampling_wallclock_time == 2.0

    global_reset_det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=CMOSSensor(timing_model=RollingShutter(0.25; exposure_mode=GlobalResetExposure())),
        response_model=NullFrameResponse())
    constant_source = FunctionFrameSource(t -> ones(4, 4))
    global_reset_frame = capture!(global_reset_det, constant_source; rng=MersenneTwister(129))
    @test global_reset_frame == repeat(reshape([1.0, 1.25, 1.5, 1.75], :, 1), 1, 4)
    @test global_reset_det.params.timing_model.exposure_mode == GlobalResetExposure()

    interval_source = FunctionExposureFrameSource((start_time, exposure_time) ->
        fill(start_time <= 1.4 < start_time + exposure_time ? 10.0 : 0.0, 4, 4))
    rolling_interval_frame = capture!(rolling_det, interval_source; rng=MersenneTwister(132))
    @test rolling_interval_frame[1:2, :] == zeros(2, 4)
    @test rolling_interval_frame[3:4, :] == fill(10.0, 2, 4)
    global_reset_interval_frame = capture!(global_reset_det, interval_source; rng=MersenneTwister(133))
    @test global_reset_interval_frame[1:2, :] == zeros(2, 4)
    @test global_reset_interval_frame[3, :] == fill(15.0, 4)
    @test global_reset_interval_frame[4, :] == fill(17.5, 4)

    pulse_source = FunctionFrameSource(t -> fill(t >= 0.5 ? 10.0 : 0.0, 4, 4))
    pulse_frame = capture!(rolling_det, pulse_source; rng=MersenneTwister(128))
    @test pulse_frame[1:2, :] == zeros(2, 4)
    @test pulse_frame[3:4, :] == fill(10.0, 2, 4)

    qcmos_sensor = QCMOSSensor(T=Float32)
    @test detector_sensor_symbol(qcmos_sensor) == :qcmos
    @test qcmos_camera_symbol(qcmos_camera_model(qcmos_sensor)) == :orca_quest_2
    @test qcmos_scan_mode_symbol(qcmos_scan_mode(qcmos_sensor)) == :ultra_quiet
    @test qcmos_readout_noise_sigma(qcmos_sensor) ≈ 0.30f0
    @test qcmos_readout_noise_median(qcmos_sensor) ≈ 0.25f0
    @test qcmos_dsnu_sigma(qcmos_sensor) ≈ 0.06f0
    @test qcmos_prnu_sigma(qcmos_sensor) ≈ 0.001f0
    @test qcmos_conversion_factor_e_per_count(qcmos_sensor) ≈ 0.107f0
    @test qcmos_pixel_size_um(qcmos_sensor) ≈ 4.6f0
    @test qcmos_peak_qe(qcmos_sensor) ≈ 0.85f0
    @test qcmos_full_well(qcmos_sensor) ≈ 7000.0f0
    @test qcmos_dark_current(qcmos_sensor) ≈ 0.016f0
    @test qcmos_frame_size(qcmos_sensor) == (2304, 4096)
    @test qcmos_frame_rate_hz(qcmos_sensor) ≈ 25.4f0
    @test qcmos_min_exposure_time(qcmos_sensor) ≈ 33.9f-6
    @test qcmos_sensor.timing_model isa RollingShutter{Float32}
    @test qcmos_sensor.timing_model.line_time ≈ 33.9f-6
    @test qcmos_sensor.timing_model.row_group_size == 2
    @test qcmos_sensor.timing_model.exposure_mode == RollingExposure()
    @test supports_detector_defect_maps(qcmos_sensor)
    @test supports_shutter_timing(qcmos_sensor)
    @test supports_photon_number_resolving(qcmos_sensor)
    @test !supports_raw_digital_output(qcmos_sensor)
    @test supports_raw_digital_output(QCMOSSensor(scan_mode=QCMOSRawScan()))

    det_qcmos = ORCAQuest2Detector(integration_time=1e-3, response_model=NullFrameResponse(), T=Float32)
    @test det_qcmos.noise isa NoisePhotonReadout
    @test det_qcmos.noise.sigma ≈ 0.30f0
    @test det_qcmos.params.sensor isa QCMOSSensor
    @test det_qcmos.params.qe ≈ 0.85f0
    @test det_qcmos.params.full_well ≈ 7000.0f0
    @test det_qcmos.params.bits == 16
    @test det_qcmos.params.output_type == UInt16
    frame_qcmos = capture!(det_qcmos, fill(100.0f0, 4, 4); rng=MersenneTwister(130))
    @test frame_qcmos isa Matrix{UInt16}
    meta_qcmos = detector_export_metadata(det_qcmos)
    @test meta_qcmos.sensor == :qcmos
    @test meta_qcmos.noise == :photon_readout
    @test meta_qcmos.readout_sigma ≈ 0.30f0
    @test meta_qcmos.frame_response == :none
    @test meta_qcmos.timing_model == :rolling_shutter
    @test meta_qcmos.timing_line_time ≈ 33.9f-6
    @test meta_qcmos.sampling_wallclock_time ≈ inv(25.4f0)

    qcmos_standard = ORCAQuest2Sensor(scan_mode=QCMOSStandardScan(), T=Float32)
    @test qcmos_readout_noise_sigma(qcmos_standard) ≈ 0.43f0
    @test qcmos_readout_noise_median(qcmos_standard) ≈ 0.39f0
    @test qcmos_frame_rate_hz(qcmos_standard) ≈ 120.0f0
    @test qcmos_min_exposure_time(qcmos_standard) ≈ 7.2f-6
    @test qcmos_standard.timing_model.line_time ≈ 7.2f-6
    qcmos_legacy = ORCAQuestSensor(scan_mode=QCMOSUltraQuietScan(), T=Float32)
    @test qcmos_readout_noise_sigma(qcmos_legacy) ≈ 0.27f0
    @test qcmos_frame_rate_hz(qcmos_legacy) ≈ 5.0f0
    @test qcmos_min_exposure_time(qcmos_legacy) ≈ 199.9f-3
    @test qcmos_legacy.timing_model.line_time ≈ 172.8f-6
    qcmos_iq = ORCAQuestIQSensor(T=Float32)
    @test qcmos_dark_current(qcmos_iq) ≈ 0.032f0
    det_qcmos_custom = ORCAQuestIQDetector(scan_mode=QCMOSStandardScan(), integration_time=2e-3,
        response_model=NullFrameResponse(), noise=NoiseNone(), T=Float32)
    @test det_qcmos_custom.noise isa NoiseNone
    @test det_qcmos_custom.params.integration_time ≈ 2.0f-3
    @test det_qcmos_custom.params.sensor.scan_mode isa QCMOSStandardScan
    det_qcmos_temporal = ORCAQuest2Detector(integration_time=1.0, response_model=NullFrameResponse(),
        noise=NoiseNone(), qe=1.0, dark_current=0.0, bits=nothing, output_type=nothing, T=Float32)
    qcmos_temporal_source = InPlaceFrameSource((out, t) -> fill!(out, t * 1.0f6), (4, 4))
    qcmos_temporal = capture!(det_qcmos_temporal, qcmos_temporal_source; rng=MersenneTwister(131))
    @test qcmos_temporal[1:2, :] == zeros(Float32, 2, 4)
    @test qcmos_temporal[3:4, :] == fill(33.9f0, 2, 4)
    @test_throws InvalidConfiguration QCMOSSensor(readout_noise_sigma=-1.0)
    @test_throws InvalidConfiguration QCMOSSensor(peak_qe=1.5)
    @test_throws InvalidConfiguration QCMOSSensor(full_well=0.0)
    @test_throws InvalidConfiguration QCMOSSensor(frame_size=(0, 4096))

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
    @test single_products isa MultiReadFrameReadoutProducts
    @test detector_reference_frame(det_saphira_single) === nothing
    @test detector_signal_frame(det_saphira_single) !== nothing
    @test detector_combined_frame(det_saphira_single) == frame_saphira_single
    @test detector_reference_cube(det_saphira_single) === nothing
    @test detector_signal_cube(det_saphira_single) !== nothing
    @test detector_read_cube(det_saphira_single) === nothing
    @test detector_read_times(det_saphira_single) === nothing

    @test fieldtype(typeof(det.state), :readout_products) === NoFrameReadoutProducts
    @test fieldtype(typeof(det_saphira_single.state), :readout_products) !== FrameReadoutProducts

    struct DummyReadoutProducts{A,V} <: FrameReadoutProducts
        signal_frame::A
        read_times::V
    end
    AdaptiveOpticsSim.detector_signal_frame(products::DummyReadoutProducts) = products.signal_frame
    AdaptiveOpticsSim.detector_read_times(products::DummyReadoutProducts) = products.read_times

    dummy_products = DummyReadoutProducts(fill(3.0, 2, 2), [0.25, 0.5])
    @test detector_reference_frame(dummy_products) === nothing
    @test detector_signal_frame(dummy_products) == fill(3.0, 2, 2)
    @test detector_combined_frame(dummy_products) === nothing
    @test detector_reference_cube(dummy_products) === nothing
    @test detector_signal_cube(dummy_products) === nothing
    @test detector_read_cube(dummy_products) === nothing
    @test detector_read_times(dummy_products) == [0.25, 0.5]

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
    det_saphira_timed_glow = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        dark_current=3.0, gain=1.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(glow_rate=2.0, read_time=1.0, sampling_mode=CorrelatedDoubleSampling()))
    frame_saphira_timed_glow = copy(capture!(det_saphira_timed_glow, zero_psf; rng=MersenneTwister(17)))
    det_saphira_timed_noglow = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        dark_current=3.0, gain=1.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(read_time=1.0, sampling_mode=CorrelatedDoubleSampling()))
    frame_saphira_timed_noglow = copy(capture!(det_saphira_timed_noglow, zero_psf; rng=MersenneTwister(17)))
    @test sum(frame_saphira_timed_glow) > sum(frame_saphira_timed_noglow)
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
    @test detector_combined_frame(det_saphira_windowed) ≈
        detector_signal_frame(det_saphira_windowed) .- detector_reference_frame(det_saphira_windowed)
    det_saphira_windowed_corrected = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(read_time=1.0, sampling_mode=CorrelatedDoubleSampling()),
        readout_window=FrameWindow(2:3, 2:3),
        correction_model=CompositeFrameReadoutCorrection((
            ReferenceRowCommonModeCorrection(1),
            ReferenceColumnCommonModeCorrection(1),
        )))
    row_pattern = repeat(reshape([1.0, 2.0, 3.0, 4.0], :, 1), 1, 4)
    col_pattern = repeat(reshape([1.0, 2.0, 3.0, 4.0], 1, :), 4, 1)
    windowed_corrected_input = row_pattern .+ col_pattern
    windowed_corrected = capture!(det_saphira_windowed_corrected, windowed_corrected_input; rng=MersenneTwister(18))
    @test maximum(abs, windowed_corrected) < 1e-6
    @test detector_combined_frame(det_saphira_windowed_corrected) ≈
        detector_signal_frame(det_saphira_windowed_corrected) .- detector_reference_frame(det_saphira_windowed_corrected)
    @test size(detector_read_cube(det_saphira_windowed_corrected)) == (2, 2, 2)
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
    row_corrected = capture!(det_saphira_row_corrected, row_pattern; rng=MersenneTwister(20))
    @test maximum(abs, row_corrected) < 1e-6
    det_saphira_col_corrected = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=1.0, sensor=HgCdTeAvalancheArraySensor(),
        response_model=NullFrameResponse(),
        correction_model=ReferenceColumnCommonModeCorrection(1))
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
    @test supports_channel_crosstalk(spad_crosstalk)

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
    det_stack_adc = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        bits=8, full_well=10.0, output_type=UInt16)
    cube_stack_adc = fill(10.0, 2, 4, 4)
    stack_adc = AdaptiveOpticsSim.capture_stack!(det_stack_adc, cube_stack_adc, similar(cube_stack_adc);
        rng=MersenneTwister(10))
    @test stack_adc === cube_stack_adc
    @test all(stack_adc .== 255.0)

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

    corrected_stack_models = (
        ReferencePixelCommonModeCorrection(1, 1),
        ReferenceRowCommonModeCorrection(1),
        ReferenceColumnCommonModeCorrection(1),
        ReferenceOutputCommonModeCorrection(2; edge_rows=1, edge_cols=1),
        CompositeFrameReadoutCorrection((
            ReferenceRowCommonModeCorrection(1),
            ReferenceColumnCommonModeCorrection(1),
        )),
    )
    for correction_model in corrected_stack_models
        corrected_stack_det = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
            sensor=HgCdTeAvalancheArraySensor(sampling_mode=SingleRead()),
            correction_model=correction_model)
        corrected_stack_in = Array{Float64}(undef, 2, 5, 5)
        corrected_stack_in[1, :, :] .= reshape(collect(1.0:25.0), 5, 5)
        corrected_stack_in[2, :, :] .= reshape(collect(26.0:50.0), 5, 5)
        corrected_stack_ref = copy(corrected_stack_in)
        corrected_stack = AdaptiveOpticsSim.capture_stack!(corrected_stack_det, corrected_stack_in,
            similar(corrected_stack_in); rng=MersenneTwister(10))
        @test size(corrected_stack) == size(corrected_stack_in)
        corrected_frame_1 = capture!(Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
                sensor=HgCdTeAvalancheArraySensor(sampling_mode=SingleRead()),
                correction_model=correction_model),
            @view(corrected_stack_ref[1, :, :]); rng=MersenneTwister(10))
        corrected_frame_2 = capture!(Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
                sensor=HgCdTeAvalancheArraySensor(sampling_mode=SingleRead()),
                correction_model=correction_model),
            @view(corrected_stack_ref[2, :, :]); rng=MersenneTwister(10))
        @test corrected_stack[1, :, :] ≈ corrected_frame_1 atol=1e-12 rtol=1e-12
        @test corrected_stack[2, :, :] ≈ corrected_frame_2 atol=1e-12 rtol=1e-12
    end

    det_cmos_batched = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=CMOSSensor(column_readout_sigma=1.0))
    @test_throws InvalidConfiguration AdaptiveOpticsSim.capture_stack!(det_cmos_batched, cube_mtf, scratch_mtf; rng=MersenneTwister(10))

    det_generalized = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, psf_sampling=2, binning=2,
        bits=8, full_well=10.0, output_type=UInt8)
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

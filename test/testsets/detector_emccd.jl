@testset "EMCCD detector" begin
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
    @test isapprox(mean(frame_emccd_stochastic), 250.0; rtol=0.1)
    @test isapprox(var(frame_emccd_stochastic), 1200.0; rtol=0.35)

    det_emccd_cic = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=EMCCDSensor(clock_induced_charge_per_frame=3.0))
    frame_emccd_cic = copy(capture!(det_emccd_cic, zero_psf;
        rng=MersenneTwister(125)))
    @test sum(frame_emccd_cic) > 0
    det_emccd_cic_long = Detector(integration_time=10.0, noise=NoiseNone(),
        qe=1.0, sensor=EMCCDSensor(clock_induced_charge_per_frame=3.0))
    @test capture!(det_emccd_cic_long, zero_psf;
        rng=MersenneTwister(125)) == frame_emccd_cic
    emccd_cic_whole = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, gain=5.0, response_model=NullFrameResponse(),
        sensor=EMCCDSensor(clock_induced_charge_per_frame=3.0))
    emccd_cic_split = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, gain=5.0, response_model=NullFrameResponse(),
        sensor=EMCCDSensor(clock_induced_charge_per_frame=3.0))
    emccd_cic_whole_frame = copy(capture!(emccd_cic_whole, zero_psf;
        rng=MersenneTwister(130)))
    capture!(emccd_cic_split, zero_psf; rng=MersenneTwister(130),
        integration_duration=0.5)
    emccd_cic_split_frame = copy(capture!(emccd_cic_split, zero_psf;
        rng=MersenneTwister(130), integration_duration=0.5))
    @test emccd_cic_split_frame == emccd_cic_whole_frame
    det_emccd_sat = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=5.0, sensor=EMCCDSensor(register_full_well=100.0))
    @test maximum(capture!(det_emccd_sat, fill(50.0, 4, 4); rng=MersenneTwister(126))) == 100.0
    emccd_saturated_stack = fill(50.0, 2, 4, 4)
    AdaptiveOpticsSim.Detectors.capture_stack!(det_emccd_sat, emccd_saturated_stack,
        similar(emccd_saturated_stack); rng=MersenneTwister(126))
    @test all(==(100.0), emccd_saturated_stack)
    @test_throws InvalidConfiguration EMCCDSensor(
        clock_induced_charge_per_frame=-1.0)
    @test_throws InvalidConfiguration EMCCDSensor(register_full_well=0.0)
    @test_throws InvalidConfiguration EMCCDSensor(multiplication_model=StochasticMultiplicationRegister(-1.0))
    @test_throws InvalidConfiguration EMCCDSensor(em_gain_range=(10.0, 1.0))
    @test_throws InvalidConfiguration EMCCDSensor(readout_rate_hz=-1.0)
    @test_throws InvalidConfiguration FrameTransferAcquisition(
        transfer_time=-1.0)
    @test_throws InvalidConfiguration PhotonCountingEMMode(threshold=-1.0)
    @test_throws InvalidConfiguration PhotonCountingEMMode(threshold=1.0, detection_efficiency=1.5)

    emccd_timing_input = fill(2.0, 4, 4)
    emccd_sequential = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, gain=1.0, response_model=NullFrameResponse(),
        sensor=EMCCDSensor(readout_rate_hz=1000.0,
            acquisition_mode=SequentialAcquisition()))
    emccd_frame_transfer = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, gain=1.0,
        response_model=NullFrameResponse(),
        sensor=EMCCDSensor(readout_rate_hz=1000.0,
            acquisition_mode=FrameTransferAcquisition(
                transfer_time=0.002)))
    sequential_frame = copy(capture!(emccd_sequential, emccd_timing_input;
        rng=MersenneTwister(126)))
    frame_transfer_frame = copy(capture!(emccd_frame_transfer,
        emccd_timing_input; rng=MersenneTwister(126)))
    @test frame_transfer_frame == sequential_frame
    sequential_meta = detector_export_metadata(emccd_sequential)
    frame_transfer_meta = detector_export_metadata(emccd_frame_transfer)
    @test sequential_meta.acquisition_mode == :sequential
    @test sequential_meta.frame_transfer_time === nothing
    @test sequential_meta.sampling_read_time == 0.016
    @test sequential_meta.sampling_wallclock_time == 1.016
    @test sequential_meta.steady_state_frame_period == 1.016
    @test frame_transfer_meta.acquisition_mode == :frame_transfer
    @test frame_transfer_meta.frame_transfer_time == 0.002
    @test frame_transfer_meta.sampling_read_time == 0.016
    @test frame_transfer_meta.sampling_wallclock_time == 1.018
    @test frame_transfer_meta.steady_state_frame_period == 1.002

    emccd_readout_limited = Detector(integration_time=0.001,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
        sensor=EMCCDSensor(readout_rate_hz=1000.0,
            acquisition_mode=FrameTransferAcquisition(
                transfer_time=0.002)))
    capture!(emccd_readout_limited, emccd_timing_input;
        rng=MersenneTwister(126))
    readout_limited_meta = detector_export_metadata(emccd_readout_limited)
    @test readout_limited_meta.steady_state_frame_period ≈ 0.018

    emccd_unknown_timing = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0,
        sensor=EMCCDSensor(acquisition_mode=FrameTransferAcquisition()))
    capture!(emccd_unknown_timing, emccd_timing_input;
        rng=MersenneTwister(126))
    unknown_timing_meta = detector_export_metadata(emccd_unknown_timing)
    @test unknown_timing_meta.sampling_read_time === nothing
    @test unknown_timing_meta.sampling_wallclock_time === nothing
    @test unknown_timing_meta.steady_state_frame_period === nothing

    det_emccd_conventional = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=10.0, sensor=EMCCDSensor(output_path=ConventionalOutput()))
    @test capture!(det_emccd_conventional, uniform_signal; rng=MersenneTwister(127)) == uniform_signal

    det_emccd_pc = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        gain=10.0,
        sensor=EMCCDSensor(operating_mode=PhotonCountingEMMode(threshold=5.0)))
    pc_frame = capture!(det_emccd_pc, [0.0 0.4; 0.6 1.0]; rng=MersenneTwister(128))
    @test pc_frame == [0.0 0.0; 1.0 1.0]
    @test supports_photon_number_resolving(det_emccd_pc.params.sensor)

    det_emccd_efficiency = Detector(integration_time=1.0, noise=NoiseNone(),
        qe=1.0, gain=10.0,
        sensor=EMCCDSensor(operating_mode=PhotonCountingEMMode(
            threshold=5.0, detection_efficiency=0.8)))
    efficiency_frame = capture!(det_emccd_efficiency, ones(100, 100);
        rng=MersenneTwister(128))
    @test all(x -> x == 0.0 || x == 1.0, efficiency_frame)
    @test isapprox(mean(efficiency_frame), 0.8; atol=0.025)

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

end

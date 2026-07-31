function emccd_qualification_samples(detector, input, rng,
    frame_count::Int)
    samples = Vector{Float64}(undef, frame_count * length(input))
    offset = 0
    for _ in 1:frame_count
        frame = capture!(detector, input, rng)
        copyto!(samples, offset + 1, vec(frame), 1, length(frame))
        offset += length(frame)
    end
    return samples
end

function test_emccd_moments(samples, expected_mean, expected_variance,
    expected_fourth_central; sigma_limit=6.0)
    sample_count = length(samples)
    mean_limit = sigma_limit * sqrt(expected_variance / sample_count)
    variance_se = sqrt((expected_fourth_central -
        ((sample_count - 3) / (sample_count - 1)) *
        expected_variance^2) / sample_count)
    @test abs(mean(samples) - expected_mean) <= mean_limit
    @test abs(var(samples) - expected_variance) <=
        sigma_limit * variance_se
end

function test_emccd_cdf(samples, threshold, expected_probability;
    sigma_limit=6.0)
    sample_count = length(samples)
    observed_probability =
        count(value -> value <= threshold, samples) / sample_count
    probability_se = sqrt(expected_probability *
        (1 - expected_probability) / sample_count)
    @test abs(observed_probability - expected_probability) <=
        sigma_limit * probability_se
end

struct TestUnsupportedEMGainModel <: AbstractEMGainModel end
struct TestUnsupportedEMCCDOperatingMode <:
    AbstractEMCCDOperatingMode end
struct TestUnsupportedEMCCDOutputPath <: AbstractEMCCDOutputPath end
struct TestUnsupportedEMCCDAcquisitionMode <:
    AbstractEMCCDAcquisitionMode end

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
            multiplication_model=ConditionalGammaMultiplication(
                minimum_conditional_noise_factor=0.6)))
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
    @test_throws InvalidConfiguration ConditionalGammaMultiplication(
        minimum_conditional_noise_factor=-1.0)
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
    @test supports_photon_counting(det_emccd_pc.params.sensor)
    @test !supports_photon_number_resolving(det_emccd_pc.params.sensor)

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

@testset "EMCCD parameter, model, and saturation contracts" begin
    @test Base.ispublic(AdaptiveOpticsSim.Detectors,
        :ClippedGaussianMultiplicationApproximation)
    @test Base.ispublic(AdaptiveOpticsSim.Detectors,
        :ConditionalGammaMultiplication)
    @test_throws InvalidConfiguration EMCCDSensor(excess_noise_factor=Inf)
    @test_throws InvalidConfiguration EMCCDSensor(
        clock_induced_charge_per_frame=NaN)
    @test_throws InvalidConfiguration EMCCDSensor(register_full_well=Inf)
    @test_throws InvalidConfiguration EMCCDSensor(
        em_gain_range=(Inf, Inf))
    @test_throws InvalidConfiguration EMCCDSensor(
        em_gain_range=(1.0, NaN))
    @test_throws InvalidConfiguration EMCCDSensor(readout_rate_hz=Inf)
    @test_throws InvalidConfiguration FrameTransferAcquisition(
        transfer_time=Inf)
    @test_throws InvalidConfiguration PhotonCountingEMMode(threshold=Inf)
    @test_throws InvalidConfiguration PhotonCountingEMMode(
        threshold=1.0, detection_efficiency=NaN)
    @test_throws InvalidConfiguration ClippedGaussianMultiplicationApproximation(
        minimum_conditional_noise_factor=Inf)
    @test_throws InvalidConfiguration ConditionalGammaMultiplication(
        minimum_conditional_noise_factor=NaN)
    approximate_model = ClippedGaussianMultiplicationApproximation()
    conditional_model = ConditionalGammaMultiplication()
    @test is_approximate_em_gain_model(approximate_model)
    @test !is_approximate_em_gain_model(conditional_model)
    @test em_gain_model_symbol(approximate_model) ==
        :clipped_gaussian_approximation
    @test em_gain_model_symbol(conditional_model) == :conditional_gamma
    @test supports_photon_counting(EMCCDSensor(
        operating_mode=PhotonCountingEMMode(threshold=1.0)))
    @test !supports_photon_counting(CCDSensor())

    unsupported_gain = TestUnsupportedEMGainModel()
    unsupported_mode = TestUnsupportedEMCCDOperatingMode()
    unsupported_output = TestUnsupportedEMCCDOutputPath()
    unsupported_acquisition = TestUnsupportedEMCCDAcquisitionMode()
    @test convert_em_gain_model(unsupported_gain, Float32) ===
        unsupported_gain
    @test convert_emccd_operating_mode(unsupported_mode, Float32) ===
        unsupported_mode
    @test convert_emccd_acquisition_mode(unsupported_acquisition,
        Float32) === unsupported_acquisition
    @test_throws InvalidConfiguration validate_em_gain_model(
        unsupported_gain)
    @test_throws InvalidConfiguration validate_emccd_operating_mode(
        unsupported_mode)
    @test_throws InvalidConfiguration validate_emccd_output_path(
        unsupported_output)
    @test_throws InvalidConfiguration validate_emccd_acquisition_mode(
        unsupported_acquisition)

    sensor32 = EMCCDSensor(
        multiplication_model=ConditionalGammaMultiplication(
            minimum_conditional_noise_factor=0.75),
        register_full_well=100.0,
        em_gain_range=(2.0, 8.0),
        readout_rate_hz=10_000.0,
    )
    detector32 = Detector(sensor=sensor32, gain=4.0, T=Float32)
    @test detector32.params.sensor isa
        EMCCDSensor{Float32,<:ConditionalGammaMultiplication}
    @test detector32.params.sensor.multiplication_model isa
        ConditionalGammaMultiplication{Float32}
    @test_throws InvalidConfiguration Detector(sensor=sensor32, gain=1.0)
    @test_throws InvalidConfiguration Detector(sensor=sensor32, gain=9.0)
    @test_throws InvalidConfiguration Detector(sensor=sensor32, gain=Inf)
    @test_throws InvalidConfiguration begin
        AdaptiveOpticsSim.Detectors.validate_em_gain_backend(
            EMOutput(), detector32.params.sensor.multiplication_model,
            CUDABackend())
    end
    @test isnothing(
        AdaptiveOpticsSim.Detectors.validate_em_gain_backend(
            ConventionalOutput(),
            detector32.params.sensor.multiplication_model,
            CUDABackend()))
    @test_throws InvalidConfiguration begin
        _apply_conditional_gamma_multiplication!(
            AcceleratorStyle(KernelAbstractions.CPU()),
            detector32.params.sensor.multiplication_model,
            detector32.params.sensor,
            zeros(Float32, 2, 2),
            zeros(Float32, 2, 2),
            detector32.params.gain,
            Xoshiro(3300))
    end

    input_limited = Detector(integration_time=1.0, qe=1.0,
        noise=NoiseNone(), gain=5.0, full_well=10.0,
        sensor=EMCCDSensor(excess_noise_factor=1.0,
            register_full_well=100.0))
    @test capture!(input_limited, fill(50.0, 4, 4),
        Xoshiro(3301)) == fill(50.0, 4, 4)

    register_limited = Detector(integration_time=1.0, qe=1.0,
        noise=NoiseNone(), gain=5.0, full_well=30.0,
        sensor=EMCCDSensor(excess_noise_factor=1.0,
            register_full_well=100.0))
    @test capture!(register_limited, fill(50.0, 4, 4),
        Xoshiro(3302)) == fill(100.0, 4, 4)

    conventional = Detector(integration_time=1.0, qe=1.0,
        noise=NoiseNone(), gain=20.0, full_well=10.0,
        sensor=EMCCDSensor(output_path=ConventionalOutput(),
            em_gain_range=(1.0, 2.0)))
    @test capture!(conventional, fill(50.0, 4, 4),
        Xoshiro(3303)) == fill(10.0, 4, 4)

    batched_gamma = Detector(integration_time=1.0, qe=1.0,
        noise=NoiseNone(), gain=5.0,
        response_model=NullFrameResponse(),
        sensor=EMCCDSensor(excess_noise_factor=sqrt(2.0),
            multiplication_model=ConditionalGammaMultiplication()))
    gamma_cube = fill(8.0, 2, 4, 4)
    capture_stack!(batched_gamma, gamma_cube, similar(gamma_cube),
        Xoshiro(3304))
    @test all(isfinite, gamma_cube)
    @test all(>(0.0), gamma_cube)
end

@testset "EMCCD multiplication qualification" begin
    frame_count = 16
    frame_size = 32

    exponential_gain = 12.0
    exponential_detector = Detector(integration_time=1.0, qe=1.0,
        noise=NoiseNone(), gain=exponential_gain,
        sensor=EMCCDSensor(excess_noise_factor=sqrt(2.0),
            multiplication_model=ConditionalGammaMultiplication()))
    exponential_samples = emccd_qualification_samples(
        exponential_detector, ones(frame_size, frame_size),
        Xoshiro(3401), frame_count)
    exponential_variance = exponential_gain^2
    test_emccd_moments(exponential_samples, exponential_gain,
        exponential_variance, 9 * exponential_gain^4)
    test_emccd_cdf(exponential_samples, exponential_gain * log(2), 0.5)
    test_emccd_cdf(exponential_samples, exponential_gain * log(10), 0.9)

    erlang_charge = 8.0
    erlang_gain = 5.0
    erlang_detector = Detector(integration_time=1.0, qe=1.0,
        noise=NoiseNone(), gain=erlang_gain,
        sensor=EMCCDSensor(excess_noise_factor=sqrt(2.0),
            multiplication_model=ConditionalGammaMultiplication()))
    erlang_samples = emccd_qualification_samples(erlang_detector,
        fill(erlang_charge, frame_size, frame_size),
        Xoshiro(3402), frame_count)
    erlang_variance = erlang_charge * erlang_gain^2
    erlang_fourth = 3 * erlang_charge * (erlang_charge + 2) *
        erlang_gain^4
    test_emccd_moments(erlang_samples, erlang_charge * erlang_gain,
        erlang_variance, erlang_fourth)
    erlang_mean_cdf = 1 - exp(-erlang_charge) *
        sum(erlang_charge^order / factorial(order)
            for order in 0:(Int(erlang_charge) - 1))
    test_emccd_cdf(erlang_samples, erlang_charge * erlang_gain,
        erlang_mean_cdf)

    fractional_charge = 10.0
    fractional_gain = 5.0
    fractional_factor2 = 1.4^2 - 1
    fractional_shape = fractional_charge / fractional_factor2
    fractional_variance =
        fractional_charge * fractional_gain^2 * fractional_factor2
    fractional_fourth = (3 + 6 / fractional_shape) *
        fractional_variance^2
    fractional_detector = Detector(integration_time=1.0, qe=1.0,
        noise=NoiseNone(), gain=fractional_gain,
        sensor=EMCCDSensor(excess_noise_factor=1.4,
            multiplication_model=ConditionalGammaMultiplication()))
    fractional_samples = emccd_qualification_samples(
        fractional_detector,
        fill(fractional_charge, frame_size, frame_size),
        Xoshiro(3403), frame_count)
    test_emccd_moments(fractional_samples,
        fractional_charge * fractional_gain,
        fractional_variance, fractional_fourth)

    approximate_charge = 64.0
    approximate_gain = 5.0
    approximate_variance = approximate_charge * approximate_gain^2
    approximate_detector = Detector(integration_time=1.0, qe=1.0,
        noise=NoiseNone(), gain=approximate_gain,
        sensor=EMCCDSensor(excess_noise_factor=sqrt(2.0),
            multiplication_model=
                ClippedGaussianMultiplicationApproximation()))
    approximate_samples = emccd_qualification_samples(
        approximate_detector,
        fill(approximate_charge, frame_size, frame_size),
        Xoshiro(3404), frame_count)
    @test all(>=(0), approximate_samples)
    test_emccd_moments(approximate_samples,
        approximate_charge * approximate_gain,
        approximate_variance, 3 * approximate_variance^2)

    cic_mean = 8.0
    cic_gain = 5.0
    cic_variance = 2 * cic_mean * cic_gain^2
    cic_fourth = 24 * cic_mean * cic_gain^4 +
        3 * cic_variance^2
    cic_detector = Detector(integration_time=0.25, qe=1.0,
        noise=NoiseNone(), gain=cic_gain,
        sensor=EMCCDSensor(excess_noise_factor=sqrt(2.0),
            clock_induced_charge_per_frame=cic_mean,
            multiplication_model=ConditionalGammaMultiplication()))
    cic_samples = emccd_qualification_samples(cic_detector,
        zeros(frame_size, frame_size), Xoshiro(3405), frame_count)
    test_emccd_moments(cic_samples, cic_mean * cic_gain,
        cic_variance, cic_fourth)

    prepared_values = fill(8.0, 16, 16)
    prepared_map = detector_test_intensity_map(prepared_values)
    prepared_detector = Detector(integration_time=0.25, qe=1.0,
        noise=NoiseNone(), gain=5.0,
        sensor=EMCCDSensor(excess_noise_factor=sqrt(2.0),
            multiplication_model=ConditionalGammaMultiplication()))
    prepared_plan = prepare_detector_acquisition(
        prepared_detector, prepared_map)
    prepared_rng = Xoshiro(3406)
    @test @inferred(capture!(prepared_detector, prepared_map,
        prepared_plan, prepared_rng)) isa Matrix{Float64}
    @test_detector_allocation prepared_detector_capture_allocations(
        prepared_detector, prepared_map, prepared_plan,
        prepared_rng) == 0
end

@testset "EMCCD photon-counting and SNR qualification" begin
    frame_count = 16
    frame_size = 32
    efficiency = 0.75
    for (case_index, incident_mean) in enumerate((0.05, 0.5, 2.0))
        detector = Detector(integration_time=1.0, qe=1.0,
            noise=NoisePhoton(), gain=10.0,
            sensor=EMCCDSensor(excess_noise_factor=1.0,
                multiplication_model=
                    ClippedGaussianMultiplicationApproximation(),
                operating_mode=PhotonCountingEMMode(
                    threshold=5.0,
                    detection_efficiency=efficiency)))
        samples = emccd_qualification_samples(detector,
            fill(incident_mean, frame_size, frame_size),
            Xoshiro(3500 + case_index), frame_count)
        expected_probability =
            efficiency * (1 - exp(-incident_mean))
        expected_variance =
            expected_probability * (1 - expected_probability)
        expected_fourth = expected_variance *
            (1 - 3 * expected_variance)
        @test all(value -> value == 0.0 || value == 1.0, samples)
        test_emccd_moments(samples, expected_probability,
            expected_variance, expected_fourth)
        @test mean(samples) < efficiency * incident_mean
    end

    signal = 10.0
    dark = 2.0
    cic = 1.0
    read_noise = 30.0
    gain = 100.0
    excess_noise = sqrt(2.0)
    @test emccd_snr(signal; dark_electrons=dark,
        cic_electrons=cic, readout_noise=read_noise, gain,
        excess_noise_factor=excess_noise) ≈
        signal / sqrt(excess_noise^2 * (signal + dark + cic) +
            (read_noise / gain)^2)
    @test emccd_snr(signal; dark_electrons=dark,
        cic_electrons=cic, readout_noise=read_noise, gain,
        output_path=ConventionalOutput()) ≈
        signal / sqrt(signal + dark + cic + read_noise^2)
    counting_mode = PhotonCountingEMMode(
        threshold=5.0, detection_efficiency=efficiency)
    @test emccd_snr(signal; dark_electrons=dark,
        cic_electrons=cic, readout_noise=read_noise, gain,
        operating_mode=counting_mode) ≈
        efficiency * signal /
        sqrt(efficiency * (signal + dark + cic))
end

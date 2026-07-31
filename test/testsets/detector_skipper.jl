@testset "Skipper CCD detector" begin
    skipper_input = zeros(64, 64)
    skipper_single = Detector(noise=NoiseReadout(4.0), qe=1.0,
        sensor=CCDSensor(sampling_mode=SkipperSampling(1), read_time=2e-6),
        response_model=NullFrameResponse())
    skipper_many = Detector(noise=NoiseReadout(4.0), qe=1.0,
        sensor=CCDSensor(sampling_mode=SkipperSampling(16), read_time=2e-6),
        response_model=NullFrameResponse())
    skipper_single_frame = copy(capture!(skipper_single, skipper_input;
        rng=MersenneTwister(91)))
    skipper_many_rng = MersenneTwister(91)
    skipper_many_frame = copy(capture!(skipper_many, skipper_input,
        skipper_many_rng))
    @test std(skipper_many_frame) < 0.35 * std(skipper_single_frame)
    @test detector_signal_frame(skipper_many) == skipper_many_frame
    @test detector_combined_frame(skipper_many) == skipper_many_frame
    @test detector_read_cube(skipper_many) === nothing
    skipper_meta = detector_export_metadata(skipper_many)
    @test skipper_meta.sampling_mode == :skipper
    @test skipper_meta.sampling_reads == 16
    @test skipper_meta.sampling_signal_reads == 16
    @test skipper_meta.sampling_read_time == 2e-6
    @test skipper_meta.sampling_wallclock_time == 1.0 + 32e-6
    @test skipper_meta.provides_signal_frame
    @test !skipper_meta.provides_read_cube
    @test supports_nondestructive_reads(skipper_many.params.sensor)
    @test_throws InvalidConfiguration CCDSensor(
        sampling_mode=SkipperSampling(0))
    capture!(skipper_many, skipper_input, skipper_many_rng)
    @test_detector_allocation @allocated(capture!(skipper_many, skipper_input,
        skipper_many_rng)) == 0
    skipper_transition = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
        sensor=CCDSensor(sampling_mode=SkipperSampling(2)))
    capture!(skipper_transition, ones(2, 2); rng=MersenneTwister(92))
    @test all(==(2.0), skipper_transition.state.accum_buffer)
    capture!(skipper_transition, zeros(2, 2); rng=MersenneTwister(93),
        integration_duration=0.5)
    skipper_incremental_frame = copy(capture!(skipper_transition,
        zeros(2, 2); rng=MersenneTwister(94), integration_duration=0.5))
    @test all(iszero, skipper_incremental_frame)
    @test readout_ready(skipper_transition)
end

@testset "CMOS-family detector" begin
    det_cmos = Detector(integration_time=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=CMOSSensor(column_readout_sigma=1.0))
    frame_cmos = copy(capture!(det_cmos, zeros(8, 8); rng=MersenneTwister(12)))
    @test !all(iszero, frame_cmos)
    @test all(j -> isapprox(std(frame_cmos[:, j]), 0.0; atol=1e-8), axes(frame_cmos, 2))
    @test std(vec(frame_cmos[1, :])) > 0
    @test supports_column_readout_noise(det_cmos.params.sensor)
    @test detector_export_metadata(det_cmos).frame_response == :none
    @test_throws InvalidConfiguration CMOSSensor(column_readout_sigma=-1.0)

    det_cmos_rows = Detector(noise=NoiseNone(), qe=1.0,
        sensor=CMOSSensor(row_readout_sigma=1.0))
    frame_cmos_rows = copy(capture!(det_cmos_rows, zeros(8, 8);
        rng=MersenneTwister(1212)))
    @test all(i -> isapprox(std(frame_cmos_rows[i, :]), 0.0; atol=1e-8),
        axes(frame_cmos_rows, 1))
    @test std(frame_cmos_rows[:, 1]) > 0

    sigma_map = zeros(4, 4)
    sigma_map[2, 3] = 2.0
    det_cmos_map = Detector(noise=NoiseNone(), qe=1.0,
        sensor=CMOSSensor(readout_noise_model=CMOSReadNoiseMap(sigma_map)))
    frame_cmos_map = copy(capture!(det_cmos_map, zeros(4, 4);
        rng=MersenneTwister(1213)))
    @test count(x -> !iszero(x), frame_cmos_map) == 1
    @test frame_cmos_map[2, 3] != 0
    @test_throws InvalidConfiguration CMOSSensor(row_readout_sigma=-1.0)
    @test_throws InvalidConfiguration CMOSReadNoiseMap(fill(-1.0, 2, 2))
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
    rolling_rng = MersenneTwister(127)
    rolling_frame = capture!(rolling_det, rolling_source, rolling_rng)
    @test rolling_frame == repeat(reshape([0.0, 0.25, 0.5, 0.75], :, 1), 1, 4)
    @test detector_export_metadata(rolling_det).sampling_wallclock_time == 2.0
    @test_detector_allocation @allocated(
        capture!(rolling_det, rolling_source, rolling_rng)) == 0

    global_exposure_calls = Tuple{Float64,Float64}[]
    global_exposure_source = FunctionExposureFrameSource(
        (start_time, exposure_time) -> begin
            push!(global_exposure_calls, (start_time, exposure_time))
            fill(exposure_time, 2, 2)
        end)
    global_exposure_det = Detector(integration_time=2.0,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse())
    global_exposure_frame = capture!(global_exposure_det,
        global_exposure_source; rng=MersenneTwister(1271))
    @test global_exposure_calls == [(0.0, 2.0)]
    @test global_exposure_frame == fill(4.0, 2, 2)

    rolling_exposure_calls = Tuple{Float64,Float64}[]
    rolling_exposure_source = InPlaceExposureFrameSource(
        (out, start_time, exposure_time) -> begin
            push!(rolling_exposure_calls, (start_time, exposure_time))
            fill!(out, start_time + exposure_time)
        end, (4, 4))
    rolling_exposure_det = Detector(integration_time=1.0,
        noise=NoiseNone(), qe=1.0,
        sensor=CMOSSensor(timing_model=RollingShutter(0.25;
            row_group_size=2)), response_model=NullFrameResponse())
    rolling_exposure_frame = capture!(rolling_exposure_det,
        rolling_exposure_source; rng=MersenneTwister(1272))
    @test rolling_exposure_calls == [(0.0, 1.0), (0.25, 1.0)]
    @test rolling_exposure_frame ==
        repeat(reshape([1.0, 1.0, 1.25, 1.25], :, 1), 1, 4)

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


end

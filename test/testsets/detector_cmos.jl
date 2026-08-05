struct UnsupportedCMOSReadNoise <: AbstractCMOSReadNoiseModel end
struct UnsupportedCMOSOutput <: AbstractCMOSOutputModel end

function cmos_covariance_samples(detector, rng, frame_count::Int,
    frame_size::Int)
    iseven(frame_size) || error("CMOS covariance frame size must be even")
    samples_per_frame = frame_size ÷ 2
    sample_count = frame_count * samples_per_frame
    x11 = Vector{Float64}(undef, sample_count)
    x12 = similar(x11)
    x21 = similar(x11)
    x22 = similar(x11)
    input = zeros(frame_size, frame_size)
    offset = 0
    for _ in 1:frame_count
        frame = capture!(detector, input, rng)
        @inbounds for block in 1:samples_per_frame
            row = 2block - 1
            col = 2block - 1
            index = offset + block
            x11[index] = frame[row, col]
            x12[index] = frame[row, col + 1]
            x21[index] = frame[row + 1, col]
            x22[index] = frame[row + 1, col + 1]
        end
        offset += samples_per_frame
    end
    return (; x11, x12, x21, x22)
end

function cmos_sample_covariance(x, y)
    return sum((x .- mean(x)) .* (y .- mean(y))) / (length(x) - 1)
end

function cmos_moment_bounds(sample_count, variance;
    covariance=variance, sigma_limit=6.0)
    mean_limit = sigma_limit * sqrt(variance / sample_count)
    variance_limit = sigma_limit * variance *
        sqrt(2 / (sample_count - 1))
    covariance_limit = sigma_limit *
        sqrt((variance^2 + covariance^2) / sample_count)
    return (; mean_limit, variance_limit, covariance_limit)
end

@testset "CMOS-family detector" begin
    det_cmos = Detector(exposure_duration=1.0, noise=NoiseNone(), qe=1.0, binning=1,
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
    @test_throws InvalidConfiguration CMOSReadNoiseMap(fill(NaN, 2, 2))
    prnu_map = [1.0 0.5 1.0 0.5; 1.0 0.5 1.0 0.5; 1.0 0.5 1.0 0.5; 1.0 0.5 1.0 0.5]
    dsnu_map = fill(0.25, 4, 4)
    bad_mask = falses(4, 4)
    bad_mask[2, 3] = true
    det_cmos_structured = Detector(exposure_duration=1.0, noise=NoiseNone(), qe=1.0, binning=1,
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
    @test structured_meta.timing_line_duration == 1e-3
    @test structured_meta.sampling_acquisition_duration == 1.004
    @test supports_detector_defect_maps(det_cmos_structured.params.sensor)
    @test supports_shutter_timing(det_cmos_structured.params.sensor)
    @test_throws InvalidConfiguration CMOSSensor(timing_model=RollingShutter(-1.0))
    @test_throws InvalidConfiguration CMOSSensor(timing_model=RollingShutter(Inf))
    @test_throws InvalidConfiguration CMOSSensor(timing_model=RollingShutter(1e-3; row_group_size=0))
    @test_throws InvalidConfiguration CMOSSensor(output_model=StaticCMOSOutputPattern(2, [1.0], [0.0, 1.0]))
    @test_throws InvalidConfiguration CMOSSensor(column_readout_sigma=Inf)
    @test_throws InvalidConfiguration CMOSSensor(row_readout_sigma=NaN)
    @test_throws InvalidConfiguration StaticCMOSOutputPattern(
        1, [Inf], [0.0])
    @test_throws InvalidConfiguration StaticCMOSOutputPattern(
        1, [1.0], [NaN])
    @test_throws UnsupportedAlgorithm CMOSSensor(
        readout_noise_model=UnsupportedCMOSReadNoise())
    @test_throws UnsupportedAlgorithm CMOSSensor(
        output_model=UnsupportedCMOSOutput())

    rolling_det = Detector(exposure_duration=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=CMOSSensor(timing_model=RollingShutter(0.25)),
        response_model=NullFrameResponse())
    rolling_source = InPlaceFrameSource((out, t) -> fill!(out, t), (4, 4))
    rolling_rng = MersenneTwister(127)
    rolling_frame = capture!(rolling_det, rolling_source, rolling_rng)
    @test rolling_frame == repeat(reshape([0.0, 0.25, 0.5, 0.75], :, 1), 1, 4)
    @test detector_export_metadata(rolling_det).sampling_acquisition_duration == 2.0
    @test_detector_allocation @allocated(
        capture!(rolling_det, rolling_source, rolling_rng)) == 0

    global_exposure_calls = Tuple{Float64,Float64}[]
    global_exposure_source = FunctionExposureFrameSource(
        (start_offset_s, exposure_duration) -> begin
            push!(global_exposure_calls, (start_offset_s, exposure_duration))
            fill(exposure_duration, 2, 2)
        end)
    global_exposure_det = Detector(exposure_duration=2.0,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse())
    global_exposure_frame = capture!(global_exposure_det,
        global_exposure_source; rng=MersenneTwister(1271))
    @test global_exposure_calls == [(0.0, 2.0)]
    @test global_exposure_frame == fill(4.0, 2, 2)

    rolling_exposure_calls = Tuple{Float64,Float64}[]
    rolling_exposure_source = InPlaceExposureFrameSource(
        (out, start_offset_s, exposure_duration) -> begin
            push!(rolling_exposure_calls, (start_offset_s, exposure_duration))
            fill!(out, start_offset_s + exposure_duration)
        end, (4, 4))
    rolling_exposure_det = Detector(exposure_duration=1.0,
        noise=NoiseNone(), qe=1.0,
        sensor=CMOSSensor(timing_model=RollingShutter(0.25;
            row_group_size=2)), response_model=NullFrameResponse())
    rolling_exposure_frame = capture!(rolling_exposure_det,
        rolling_exposure_source; rng=MersenneTwister(1272))
    @test rolling_exposure_calls == [(0.0, 1.0), (0.25, 1.0)]
    @test rolling_exposure_frame ==
        repeat(reshape([1.0, 1.0, 1.25, 1.25], :, 1), 1, 4)

    global_reset_det = Detector(exposure_duration=1.0, noise=NoiseNone(), qe=1.0, binning=1,
        sensor=CMOSSensor(timing_model=RollingShutter(0.25; exposure_mode=GlobalResetExposure())),
        response_model=NullFrameResponse())
    constant_source = FunctionFrameSource(t -> ones(4, 4))
    global_reset_frame = capture!(global_reset_det, constant_source; rng=MersenneTwister(129))
    @test global_reset_frame == repeat(reshape([1.0, 1.25, 1.5, 1.75], :, 1), 1, 4)
    @test global_reset_det.params.timing_model.exposure_mode == GlobalResetExposure()

    interval_source = FunctionExposureFrameSource((start_offset_s, exposure_duration) ->
        fill(start_offset_s <= 1.4 < start_offset_s + exposure_duration ? 10.0 : 0.0, 4, 4))
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

    windowed_rolling = Detector(exposure_duration=1.0, noise=NoiseNone(),
        qe=1.0, response_model=NullFrameResponse(),
        sensor=CMOSSensor(timing_model=RollingShutter(0.25)),
        readout_window=FrameWindow(3:4, 1:2))
    windowed_frame = capture!(windowed_rolling,
        InPlaceFrameSource((out, time) -> fill!(out, time), (4, 4));
        rng=Xoshiro(134))
    @test windowed_frame ==
        repeat(reshape([0.5, 0.75], :, 1), 1, 2)
    windowed_metadata = detector_export_metadata(windowed_rolling)
    @test windowed_metadata.sampling_acquisition_duration == 2.0
    @test windowed_metadata.output_size == (2, 2)

    mtf_detector = Detector(exposure_duration=1.0, noise=NoiseNone(),
        qe=1.0,
        response_model=GaussianPixelResponse(response_width_px=0.7),
        sensor=CMOSSensor(timing_model=RollingShutter(0.25)))
    expected_mtf = detector_mtf(mtf_detector, 0.17, -0.09)
    capture!(mtf_detector, ones(4, 4), Xoshiro(135))
    @test detector_mtf(mtf_detector, 0.17, -0.09) == expected_mtf

    exact_pipeline = Detector(exposure_duration=1.0, noise=NoiseNone(),
        qe=1.0, gain=2.0, full_well=5.0, bits=3, output_type=UInt8,
        response_model=NullFrameResponse(),
        sensor=CMOSSensor(output_model=StaticCMOSOutputPattern(
            2, [0.5], [1.0])),
        defect_model=CompositeDetectorDefectModel(
            PixelResponseNonuniformity(fill(2.0, 2, 2)),
            DarkSignalNonuniformity(fill(1.0, 2, 2))))
    # 1 input × 2 PRNU + 1 DSNU = 3 e⁻; detector gain and output
    # pattern produce 4 DN; 3-bit/full-well scaling gives 5.6 -> 6.
    @test capture!(exact_pipeline, ones(2, 2), Xoshiro(136)) ==
        fill(UInt8(6), 2, 2)

end

@testset "CMOS-family spatial-noise qualification" begin
    frame_size = 32
    frame_count = 1024
    sample_count = frame_count * (frame_size ÷ 2)

    column_sigma = 1.5
    column_variance = column_sigma^2
    column_samples = cmos_covariance_samples(
        Detector(noise=NoiseNone(), qe=1.0,
            sensor=CMOSSensor(column_readout_sigma=column_sigma)),
        Xoshiro(1401), frame_count, frame_size)
    column_bounds = cmos_moment_bounds(sample_count, column_variance)
    @test abs(mean(column_samples.x11)) <= column_bounds.mean_limit
    @test abs(var(column_samples.x11) - column_variance) <=
        column_bounds.variance_limit
    @test abs(cmos_sample_covariance(column_samples.x11,
        column_samples.x21) - column_variance) <=
        column_bounds.covariance_limit
    @test abs(cmos_sample_covariance(column_samples.x11,
        column_samples.x12)) <=
        cmos_moment_bounds(sample_count, column_variance;
            covariance=0.0).covariance_limit

    row_sigma = 1.25
    row_variance = row_sigma^2
    row_samples = cmos_covariance_samples(
        Detector(noise=NoiseNone(), qe=1.0,
            sensor=CMOSSensor(row_readout_sigma=row_sigma)),
        Xoshiro(1402), frame_count, frame_size)
    row_bounds = cmos_moment_bounds(sample_count, row_variance)
    @test abs(mean(row_samples.x11)) <= row_bounds.mean_limit
    @test abs(var(row_samples.x11) - row_variance) <=
        row_bounds.variance_limit
    @test abs(cmos_sample_covariance(row_samples.x11,
        row_samples.x12) - row_variance) <=
        row_bounds.covariance_limit
    @test abs(cmos_sample_covariance(row_samples.x11,
        row_samples.x21)) <=
        cmos_moment_bounds(sample_count, row_variance;
            covariance=0.0).covariance_limit

    pixel_sigma = 0.75
    pixel_variance = pixel_sigma^2
    pixel_samples = cmos_covariance_samples(
        Detector(noise=NoiseNone(), qe=1.0,
            sensor=CMOSSensor(readout_noise_model=CMOSReadNoiseMap(
                fill(pixel_sigma, frame_size, frame_size)))),
        Xoshiro(1403), frame_count, frame_size)
    pixel_bounds = cmos_moment_bounds(sample_count, pixel_variance;
        covariance=0.0)
    @test abs(mean(pixel_samples.x11)) <= pixel_bounds.mean_limit
    @test abs(var(pixel_samples.x11) - pixel_variance) <=
        pixel_bounds.variance_limit
    @test abs(cmos_sample_covariance(pixel_samples.x11,
        pixel_samples.x12)) <= pixel_bounds.covariance_limit
    @test abs(cmos_sample_covariance(pixel_samples.x11,
        pixel_samples.x21)) <= pixel_bounds.covariance_limit

    combined_column_sigma = 0.8
    combined_row_sigma = 1.1
    combined_pixel_sigma = 0.6
    expected_column_covariance = combined_column_sigma^2
    expected_row_covariance = combined_row_sigma^2
    expected_variance = expected_column_covariance +
        expected_row_covariance + combined_pixel_sigma^2
    combined_samples = cmos_covariance_samples(
        Detector(noise=NoiseNone(), qe=1.0,
            sensor=CMOSSensor(
                column_readout_sigma=combined_column_sigma,
                row_readout_sigma=combined_row_sigma,
                readout_noise_model=CMOSReadNoiseMap(fill(
                    combined_pixel_sigma, frame_size, frame_size)))),
        Xoshiro(1404), frame_count, frame_size)
    combined_bounds = cmos_moment_bounds(sample_count, expected_variance)
    @test abs(mean(combined_samples.x11)) <= combined_bounds.mean_limit
    @test abs(var(combined_samples.x11) - expected_variance) <=
        combined_bounds.variance_limit
    @test abs(cmos_sample_covariance(combined_samples.x11,
        combined_samples.x21) - expected_column_covariance) <=
        cmos_moment_bounds(sample_count, expected_variance;
            covariance=expected_column_covariance).covariance_limit
    @test abs(cmos_sample_covariance(combined_samples.x11,
        combined_samples.x12) - expected_row_covariance) <=
        cmos_moment_bounds(sample_count, expected_variance;
            covariance=expected_row_covariance).covariance_limit
    @test abs(cmos_sample_covariance(combined_samples.x11,
        combined_samples.x22)) <=
        cmos_moment_bounds(sample_count, expected_variance;
            covariance=0.0).covariance_limit
end

using AdaptiveOpticsSim
using AdaptiveOpticsSim.Detectors
using Dates
using Random
using Statistics
using TOML

const CMOS_ARTIFACT_DIR =
    joinpath(@__DIR__, "..", "benchmarks", "results", "detectors")
const CMOS_ARTIFACT_PATH = joinpath(CMOS_ARTIFACT_DIR,
    "2026-07-30-cmos-family-qualification.toml")
const CMOS_ARTIFACT_ID = "DET-CMOS-QUAL-2026-07-30"
const CMOS_FRAME_COUNT = 1024
const CMOS_FRAME_SIZE = 32
const CMOS_SAMPLES_PER_CASE =
    CMOS_FRAME_COUNT * (CMOS_FRAME_SIZE ÷ 2)
const CMOS_SIGMA_LIMIT = 6.0

function cmos_artifact_samples(detector, rng)
    samples_per_frame = CMOS_FRAME_SIZE ÷ 2
    x11 = Vector{Float64}(undef, CMOS_SAMPLES_PER_CASE)
    x12 = similar(x11)
    x21 = similar(x11)
    x22 = similar(x11)
    input = zeros(CMOS_FRAME_SIZE, CMOS_FRAME_SIZE)
    offset = 0
    for _ in 1:CMOS_FRAME_COUNT
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

function cmos_artifact_covariance(x, y)
    return sum((x .- mean(x)) .* (y .- mean(y))) / (length(x) - 1)
end

function cmos_artifact_covariance_case(id, detector, seed;
    expected_variance, expected_column_covariance,
    expected_row_covariance, expected_diagonal_covariance=0.0)
    samples = cmos_artifact_samples(detector, Xoshiro(seed))
    sample_count = length(samples.x11)
    observed_mean = mean(samples.x11)
    observed_variance = var(samples.x11)
    observed_column_covariance =
        cmos_artifact_covariance(samples.x11, samples.x21)
    observed_row_covariance =
        cmos_artifact_covariance(samples.x11, samples.x12)
    observed_diagonal_covariance =
        cmos_artifact_covariance(samples.x11, samples.x22)
    mean_limit = CMOS_SIGMA_LIMIT *
        sqrt(expected_variance / sample_count)
    variance_limit = CMOS_SIGMA_LIMIT * expected_variance *
        sqrt(2 / (sample_count - 1))
    covariance_limit(expected) = CMOS_SIGMA_LIMIT *
        sqrt((expected_variance^2 + expected^2) / sample_count)
    gates = Dict{String,Bool}(
        "mean" => abs(observed_mean) <= mean_limit,
        "variance" =>
            abs(observed_variance - expected_variance) <= variance_limit,
        "same_column_covariance" => abs(observed_column_covariance -
            expected_column_covariance) <=
            covariance_limit(expected_column_covariance),
        "same_row_covariance" => abs(observed_row_covariance -
            expected_row_covariance) <=
            covariance_limit(expected_row_covariance),
        "diagonal_covariance" => abs(observed_diagonal_covariance -
            expected_diagonal_covariance) <=
            covariance_limit(expected_diagonal_covariance),
    )
    return Dict{String,Any}(
        "id" => id,
        "seed" => seed,
        "sample_count" => sample_count,
        "expected_mean" => 0.0,
        "expected_variance" => expected_variance,
        "expected_same_column_covariance" =>
            expected_column_covariance,
        "expected_same_row_covariance" => expected_row_covariance,
        "expected_diagonal_covariance" =>
            expected_diagonal_covariance,
        "observed_mean" => observed_mean,
        "observed_variance" => observed_variance,
        "observed_same_column_covariance" =>
            observed_column_covariance,
        "observed_same_row_covariance" => observed_row_covariance,
        "observed_diagonal_covariance" =>
            observed_diagonal_covariance,
        "mean_absolute_limit" => mean_limit,
        "variance_absolute_limit" => variance_limit,
        "same_column_covariance_absolute_limit" =>
            covariance_limit(expected_column_covariance),
        "same_row_covariance_absolute_limit" =>
            covariance_limit(expected_row_covariance),
        "diagonal_covariance_absolute_limit" =>
            covariance_limit(expected_diagonal_covariance),
        "gates" => gates,
    )
end

function cmos_artifact_spatial_noise_cases()
    column_sigma = 1.5
    row_sigma = 1.25
    pixel_sigma = 0.75
    combined_column_sigma = 0.8
    combined_row_sigma = 1.1
    combined_pixel_sigma = 0.6
    return Dict{String,Any}[
        cmos_artifact_covariance_case(
            "column_common",
            Detector(noise=NoiseNone(), qe=1.0,
                sensor=CMOSSensor(
                    column_readout_sigma=column_sigma)),
            4101;
            expected_variance=column_sigma^2,
            expected_column_covariance=column_sigma^2,
            expected_row_covariance=0.0,
        ),
        cmos_artifact_covariance_case(
            "row_common",
            Detector(noise=NoiseNone(), qe=1.0,
                sensor=CMOSSensor(row_readout_sigma=row_sigma)),
            4102;
            expected_variance=row_sigma^2,
            expected_column_covariance=0.0,
            expected_row_covariance=row_sigma^2,
        ),
        cmos_artifact_covariance_case(
            "pixel_independent",
            Detector(noise=NoiseNone(), qe=1.0,
                sensor=CMOSSensor(
                    readout_noise_model=CMOSReadNoiseMap(fill(
                        pixel_sigma, CMOS_FRAME_SIZE,
                        CMOS_FRAME_SIZE)))),
            4103;
            expected_variance=pixel_sigma^2,
            expected_column_covariance=0.0,
            expected_row_covariance=0.0,
        ),
        cmos_artifact_covariance_case(
            "combined",
            Detector(noise=NoiseNone(), qe=1.0,
                sensor=CMOSSensor(
                    column_readout_sigma=combined_column_sigma,
                    row_readout_sigma=combined_row_sigma,
                    readout_noise_model=CMOSReadNoiseMap(fill(
                        combined_pixel_sigma, CMOS_FRAME_SIZE,
                        CMOS_FRAME_SIZE)))),
            4104;
            expected_variance=combined_column_sigma^2 +
                combined_row_sigma^2 + combined_pixel_sigma^2,
            expected_column_covariance=combined_column_sigma^2,
            expected_row_covariance=combined_row_sigma^2,
        ),
    ]
end

function cmos_artifact_deterministic_contract()
    global_detector = Detector(exposure_duration=1.0,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
        sensor=CMOSSensor())
    global_shutter_passed =
        capture!(global_detector, fill(2.0, 4, 4), Xoshiro(4201)) ==
        fill(2.0, 4, 4)

    rolling_detector = Detector(exposure_duration=1.0,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
        sensor=CMOSSensor(timing_model=RollingShutter(0.25;
            row_group_size=2)))
    rolling_source = InPlaceExposureFrameSource(
        (out, start_offset_s, exposure_duration) ->
            fill!(out, start_offset_s + exposure_duration), (4, 4))
    rolling_exposure_passed =
        capture!(rolling_detector, rolling_source, Xoshiro(4202)) ==
        repeat(reshape([1.0, 1.0, 1.25, 1.25], :, 1), 1, 4)

    global_reset_detector = Detector(exposure_duration=1.0,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
        sensor=CMOSSensor(timing_model=RollingShutter(0.25;
            exposure_mode=GlobalResetExposure())))
    global_reset_passed = capture!(global_reset_detector,
        FunctionFrameSource(_ -> ones(4, 4)), Xoshiro(4203)) ==
        repeat(reshape([1.0, 1.25, 1.5, 1.75], :, 1), 1, 4)

    windowed_detector = Detector(exposure_duration=1.0,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
        sensor=CMOSSensor(timing_model=RollingShutter(0.25)),
        readout_window=Detectors.FrameWindow(3:4, 1:2))
    windowed_output = capture!(windowed_detector,
        InPlaceFrameSource((out, time) -> fill!(out, time), (4, 4)),
        Xoshiro(4204))
    windowed_metadata = detector_export_metadata(windowed_detector)
    window_preserves_full_frame_timing =
        windowed_output ==
            repeat(reshape([0.5, 0.75], :, 1), 1, 2) &&
        windowed_metadata.sampling_acquisition_duration == 2.0

    mtf_detector = Detector(exposure_duration=1.0,
        noise=NoiseNone(), qe=1.0,
        response_model=GaussianPixelResponse(response_width_px=0.7),
        sensor=CMOSSensor(timing_model=RollingShutter(0.25)))
    mtf_before = detector_mtf(mtf_detector, 0.17, -0.09)
    capture!(mtf_detector, ones(4, 4), Xoshiro(4205))
    configured_mtf_preserved =
        detector_mtf(mtf_detector, 0.17, -0.09) == mtf_before

    pipeline_detector = Detector(exposure_duration=1.0,
        noise=NoiseNone(), qe=1.0, gain=2.0, full_well=5.0,
        bits=3, output_type=UInt8,
        response_model=NullFrameResponse(),
        sensor=CMOSSensor(output_model=Detectors.StaticCMOSOutputPattern(
            2, [0.5], [1.0])),
        defect_model=CompositeDetectorDefectModel(
            PixelResponseNonuniformity(fill(2.0, 2, 2)),
            DarkSignalNonuniformity(fill(1.0, 2, 2))))
    readout_pipeline_passed =
        capture!(pipeline_detector, ones(2, 2), Xoshiro(4206)) ==
        fill(UInt8(6), 2, 2)

    replay_sensor = CMOSSensor(column_readout_sigma=0.8,
        row_readout_sigma=1.1,
        readout_noise_model=CMOSReadNoiseMap(fill(0.6, 8, 8)))
    replay_a = Detector(noise=NoiseNone(), qe=1.0,
        sensor=replay_sensor)
    replay_b = Detector(noise=NoiseNone(), qe=1.0,
        sensor=replay_sensor)
    replay_input = zeros(8, 8)
    deterministic_replay_passed =
        capture!(replay_a, replay_input, Xoshiro(4207)) ==
        capture!(replay_b, replay_input, Xoshiro(4207))

    allocation_detector = Detector(exposure_duration=1.0,
        noise=NoiseNone(), qe=1.0, response_model=NullFrameResponse(),
        sensor=CMOSSensor(timing_model=RollingShutter(1.0e-3)))
    allocation_source = InPlaceFrameSource(
        (out, time) -> fill!(out, time), (16, 16))
    allocation_rng = Xoshiro(4208)
    capture!(allocation_detector, allocation_source, allocation_rng)
    steady_alloc_bytes = @allocated capture!(allocation_detector,
        allocation_source, allocation_rng)

    return Dict{String,Any}(
        "global_shutter_passed" => global_shutter_passed,
        "rolling_exposure_passed" => rolling_exposure_passed,
        "global_reset_passed" => global_reset_passed,
        "window_preserves_full_frame_timing" =>
            window_preserves_full_frame_timing,
        "configured_mtf_preserved" => configured_mtf_preserved,
        "readout_pipeline_passed" => readout_pipeline_passed,
        "deterministic_replay_passed" =>
            deterministic_replay_passed,
        "steady_alloc_bytes" => steady_alloc_bytes,
        "allocation_gate_passed" => steady_alloc_bytes == 0,
    )
end

function cmos_artifact_git_revision()
    return try
        readchomp(`git rev-parse HEAD`)
    catch
        "unknown"
    end
end

function cmos_artifact_git_dirty()
    return !success(`git diff --quiet HEAD`)
end

function generate_cmos_qualification_artifact()
    cases = cmos_artifact_spatial_noise_cases()
    deterministic = cmos_artifact_deterministic_contract()
    all_noise_gates_passed =
        all(case -> all(values(case["gates"])), cases)
    all_deterministic_gates_passed = all(
        deterministic[key] for key in (
            "global_shutter_passed",
            "rolling_exposure_passed",
            "global_reset_passed",
            "window_preserves_full_frame_timing",
            "configured_mtf_preserved",
            "readout_pipeline_passed",
            "deterministic_replay_passed",
            "allocation_gate_passed",
        ))
    cpu = first(Sys.cpu_info())
    artifact = Dict{String,Any}(
        "schema_version" => 1,
        "artifact_id" => CMOS_ARTIFACT_ID,
        "family" => "parameterized_cmos",
        "all_gates_passed" =>
            all_noise_gates_passed && all_deterministic_gates_passed,
        "qualification" => Dict{String,Any}(
            "frame_count_per_case" => CMOS_FRAME_COUNT,
            "frame_size" => [CMOS_FRAME_SIZE, CMOS_FRAME_SIZE],
            "samples_per_case" => CMOS_SAMPLES_PER_CASE,
            "sigma_limit" => CMOS_SIGMA_LIMIT,
            "spatial_noise_cases" => cases,
            "deterministic" => deterministic,
        ),
        "scope" => Dict{String,Any}(
            "included" => [
                "ordinary, scientific, and quantitative low-noise CMOS through explicit parameters",
                "global, rolling-exposure, and global-reset row-band timing",
                "frame-independent Gaussian row, column, and heteroscedastic pixel read noise",
                "PRNU, DSNU, bad pixels, static output groups, saturation, gain, and quantization",
                "configured presampling response and detector MTF",
            ],
            "excluded" => [
                "temporal 1/f noise",
                "random-telegraph signal noise",
                "calibrated camera maps or profiles",
                "on-camera filtering",
                "vendor-specific readout modes",
            ],
        ),
        "environment" => Dict{String,Any}(
            "timestamp_utc" => string(Dates.now(Dates.UTC)),
            "source_revision" => cmos_artifact_git_revision(),
            "source_dirty" => cmos_artifact_git_dirty(),
            "julia_version" => string(VERSION),
            "adaptive_optics_sim_version" =>
                string(Base.pkgversion(AdaptiveOpticsSim)),
            "kernel" => string(Sys.KERNEL),
            "architecture" => string(Sys.ARCH),
            "cpu_target" => string(Sys.CPU_NAME),
            "cpu_model" => cpu.model,
            "julia_threads" => Threads.nthreads(),
        ),
    )
    mkpath(CMOS_ARTIFACT_DIR)
    open(CMOS_ARTIFACT_PATH, "w") do io
        TOML.print(io, artifact; sorted=true)
    end
    artifact["all_gates_passed"] || error(
        "CMOS qualification artifact contains a failed gate")
    println("wrote ", CMOS_ARTIFACT_PATH)
    return artifact
end

generate_cmos_qualification_artifact()

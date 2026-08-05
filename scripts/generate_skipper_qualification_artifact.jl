using AdaptiveOpticsSim
using AdaptiveOpticsSim.Detectors
using Dates
using Random
using Statistics
using TOML

const SKIPPER_ARTIFACT_PATH = joinpath(
    @__DIR__, "..", "benchmarks", "results", "detectors",
    "2026-07-31-skipper-ccd-qualification.toml")
const SKIPPER_SAMPLE_COUNT = 16_384
const SKIPPER_SIGMA_LIMIT = 6.0

function skipper_noise_case(n_samples::Int, seed::Int)
    single_sample_variance = 16.0
    expected_variance = single_sample_variance / n_samples
    detector = Detector(
        noise=NoiseReadout(4.0), qe=1.0,
        sensor=CCDSensor(
            sampling_mode=SkipperSampling(n_samples),
            sample_duration=2e-6),
        response_model=NullFrameResponse())
    samples = vec(copy(capture!(
        detector, zeros(128, 128), Xoshiro(seed))))
    mean_limit = SKIPPER_SIGMA_LIMIT *
        sqrt(expected_variance / length(samples))
    variance_limit = SKIPPER_SIGMA_LIMIT * expected_variance *
        sqrt(2 / (length(samples) - 1))
    observed_mean = mean(samples)
    observed_variance = var(samples)
    return Dict{String,Any}(
        "id" => "independent_reads_n$(n_samples)",
        "seed" => seed,
        "n_samples" => n_samples,
        "sample_count" => length(samples),
        "single_sample_sigma_electrons" => 4.0,
        "expected_mean" => 0.0,
        "expected_variance" => expected_variance,
        "observed_mean" => observed_mean,
        "observed_variance" => observed_variance,
        "mean_absolute_limit" => mean_limit,
        "variance_absolute_limit" => variance_limit,
        "mean_passed" => abs(observed_mean) <= mean_limit,
        "variance_passed" =>
            abs(observed_variance - expected_variance) <= variance_limit,
    )
end

function skipper_deterministic_contract()
    detector = Detector(
        exposure_duration=1.0, noise=NoiseNone(), qe=1.0,
        gain=2.0, full_well=10.0,
        sensor=CCDSensor(
            sampling_mode=SkipperSampling(8),
            sample_duration=2e-6),
        response_model=NullFrameResponse())
    input = fill(20.0, 16, 16)
    output = capture!(detector, input, Xoshiro(8201))
    products = AdaptiveOpticsSim.Detectors.readout_products(detector)
    metadata = detector_export_metadata(detector)

    cube = fill(3.0, 2, 16, 16)
    original = copy(cube)
    batched_rejected = try
        AdaptiveOpticsSim.Detectors.capture_stack!(
            detector, cube, similar(cube), Xoshiro(8202))
        false
    catch err
        err isa InvalidConfiguration
    end

    allocation_detector = Detector(
        exposure_duration=0.5, noise=NoiseReadout(0.25), qe=1.0,
        sensor=CCDSensor(
            sampling_mode=SkipperSampling(16),
            sample_duration=1e-6),
        response_model=NullFrameResponse())
    allocation_input = fill(2.0, 16, 16)
    allocation_rng = Xoshiro(8203)
    capture!(allocation_detector, allocation_input, allocation_rng)
    steady_alloc_bytes = @allocated capture!(
        allocation_detector, allocation_input, allocation_rng)

    replay_a = Detector(noise=NoiseReadout(1.0), qe=1.0,
        sensor=CCDSensor(sampling_mode=SkipperSampling(4)),
        response_model=NullFrameResponse())
    replay_b = Detector(noise=NoiseReadout(1.0), qe=1.0,
        sensor=CCDSensor(sampling_mode=SkipperSampling(4)),
        response_model=NullFrameResponse())
    replay_input = zeros(8, 8)
    replay_passed =
        capture!(replay_a, replay_input, Xoshiro(8204)) ==
        capture!(replay_b, replay_input, Xoshiro(8204))

    return Dict{String,Any}(
        "exact_mean_and_gain_passed" =>
            output == fill(20.0, 16, 16),
        "retained_charge_packet_passed" =>
            detector.state.response_buffer == fill(10.0, 16, 16),
        "input_referred_full_well_passed" =>
            detector.state.response_buffer == fill(10.0, 16, 16),
        "fixed_frame_storage_passed" =>
            products isa SkipperReadoutProducts &&
            fieldnames(typeof(products)) == (:mean_frame, :sample_count) &&
            size(products.mean_frame) == size(input) &&
            AdaptiveOpticsSim.Detectors.detector_read_cube(detector) ===
                nothing,
        "sample_count_passed" => products.sample_count == 8,
        "sample_duration_seconds" => metadata.sampling_read_duration,
        "sampling_acquisition_duration_s" => metadata.sampling_acquisition_duration,
        "timing_passed" =>
            metadata.sampling_read_duration == 2e-6 &&
            metadata.sampling_acquisition_duration == 1.0 + 16e-6,
        "batched_capture_rejected" => batched_rejected,
        "batched_input_unmodified" => cube == original,
        "deterministic_replay_passed" => replay_passed,
        "steady_alloc_bytes" => steady_alloc_bytes,
        "allocation_gate_passed" => steady_alloc_bytes == 0,
    )
end

git_revision() = try
    readchomp(`git rev-parse HEAD`)
catch
    "unknown"
end

git_dirty() = !success(`git diff --quiet HEAD`)

function generate_skipper_qualification_artifact()
    moment_cases = [
        skipper_noise_case(n, 8210 + index)
        for (index, n) in enumerate((1, 4, 16, 64))
    ]
    deterministic = skipper_deterministic_contract()
    deterministic_keys = (
        "exact_mean_and_gain_passed",
        "retained_charge_packet_passed",
        "input_referred_full_well_passed",
        "fixed_frame_storage_passed",
        "sample_count_passed",
        "timing_passed",
        "batched_capture_rejected",
        "batched_input_unmodified",
        "deterministic_replay_passed",
        "allocation_gate_passed",
    )
    all_gates_passed =
        all(case["mean_passed"] && case["variance_passed"]
            for case in moment_cases) &&
        all(deterministic[key] for key in deterministic_keys)
    artifact = Dict{String,Any}(
        "schema_version" => 1,
        "artifact_id" => "DET-SKIPPER-QUAL-2026-07-31",
        "family" => "skipper_ccd_independent_read",
        "all_gates_passed" => all_gates_passed,
        "environment" => Dict(
            "timestamp_utc" => string(now(UTC)),
            "julia_version" => string(VERSION),
            "julia_threads" => Threads.nthreads(),
            "architecture" => string(Sys.ARCH),
            "kernel" => string(Sys.KERNEL),
            "source_revision" => git_revision(),
            "source_dirty" => git_dirty(),
        ),
        "model" => Dict(
            "estimator" => "streaming arithmetic mean",
            "read_noise_variance" => "sigma_single^2 / n_samples",
            "sample_duration" =>
                "duration of one configured full-frame nondestructive sample",
            "scientific_reference" => "https://arxiv.org/abs/1106.1839",
        ),
        "qualification" => Dict(
            "samples_per_case" => SKIPPER_SAMPLE_COUNT,
            "sigma_limit" => SKIPPER_SIGMA_LIMIT,
            "moment_cases" => moment_cases,
            "deterministic" => deterministic,
        ),
        "scope" => Dict(
            "included" => [
                "fixed-count nondestructive samples of one retained CCD charge packet",
                "independent Gaussian read noise and streaming arithmetic mean",
                "fixed frame-sized readout-product storage",
                "explicit per-sample duration and total sampling duration",
            ],
            "excluded" => [
                "correlated 1/f read noise",
                "dual-slope-integration electronics",
                "adaptive stopping or region-dependent sample counts",
                "raw per-sample cube",
                "charge-transfer inefficiency and camera profiles",
            ],
        ),
    )
    all_gates_passed ||
        error("Skipper CCD qualification artifact contains a failed gate")
    open(SKIPPER_ARTIFACT_PATH, "w") do io
        TOML.print(io, artifact; sorted=true)
    end
    return SKIPPER_ARTIFACT_PATH
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    println(generate_skipper_qualification_artifact())
end

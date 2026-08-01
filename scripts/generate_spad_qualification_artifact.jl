using AdaptiveOpticsSim
using AdaptiveOpticsSim.Detectors
using Dates
using Random
using Statistics
using TOML

const SPAD_ARTIFACT_PATH = joinpath(
    @__DIR__, "..", "benchmarks", "results", "detectors",
    "2026-08-01-spad-qualification.toml")
const SPAD_SAMPLE_COUNT = 16_384
const SPAD_SIGMA_LIMIT = 6.0

function spad_poisson_case(id, detector, input, expected_mean, seed)
    samples = vec(copy(capture!(detector, input, Xoshiro(seed))))
    sample_count = length(samples)
    mean_limit = SPAD_SIGMA_LIMIT * sqrt(expected_mean / sample_count)
    variance_limit = SPAD_SIGMA_LIMIT * sqrt(
        (expected_mean + 2 * expected_mean^2) / (sample_count - 1))
    observed_mean = mean(samples)
    observed_variance = var(samples)
    return Dict{String,Any}(
        "id" => id,
        "seed" => seed,
        "sample_count" => sample_count,
        "expected_mean" => expected_mean,
        "expected_variance" => expected_mean,
        "observed_mean" => observed_mean,
        "observed_variance" => observed_variance,
        "mean_absolute_limit" => mean_limit,
        "variance_absolute_limit" => variance_limit,
        "mean_passed" => abs(observed_mean - expected_mean) <= mean_limit,
        "variance_passed" =>
            abs(observed_variance - expected_mean) <= variance_limit,
    )
end

function spad_poisson_cases()
    dimensions = (128, 128)
    return [
        spad_poisson_case("photon_only",
            SPADArrayDetector(dimensions; noise=NoisePhoton(),
                sensor=SPADArraySensor(
                    active_area_detection_efficiency=1.0)),
            fill(20.0, dimensions), 20.0, 9160),
        spad_poisson_case("dark_only",
            SPADArrayDetector(dimensions; integration_time=2.0,
                noise=NoisePhoton(), sensor=SPADArraySensor(
                    active_area_detection_efficiency=0.0,
                    dark_count_rate=6.0)),
            zeros(dimensions), 12.0, 9161),
        spad_poisson_case("photon_and_dark",
            SPADArrayDetector(dimensions; noise=NoisePhoton(),
                sensor=SPADArraySensor(
                    active_area_detection_efficiency=0.5,
                    dark_count_rate=3.0)),
            fill(10.0, dimensions), 8.0, 9162),
        spad_poisson_case("dead_time_adjusted_surrogate",
            SPADArrayDetector(dimensions; noise=NoisePhoton(),
                sensor=SPADArraySensor(
                    active_area_detection_efficiency=1.0,
                    dead_time_model=NonParalyzableDeadTime(0.05))),
            fill(100.0, dimensions), 100 / 6, 9163),
    ]
end

function spad_dead_time_curves()
    dead_time = 0.1
    dimensionless_rates = (0.0, 1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0)
    curves = Dict{String,Any}[]
    for (id, model, response) in (
        ("nonparalyzable", NonParalyzableDeadTime(dead_time),
            x -> (x / dead_time) / (1 + x)),
        ("paralyzable", ParalyzableDeadTime(dead_time),
            x -> (x / dead_time) * exp(-x)),
    )
        detector = SPADArrayDetector((1, 1); noise=NoiseNone(),
            sensor=SPADArraySensor(active_area_detection_efficiency=1.0,
                dead_time_model=model))
        points = Dict{String,Any}[]
        for (index, x) in enumerate(dimensionless_rates)
            observed = only(capture!(detector,
                fill(x / dead_time, 1, 1), Xoshiro(9170 + index)))
            expected = response(x)
            push!(points, Dict(
                "lambda_tau" => x,
                "expected_count" => expected,
                "observed_count" => observed,
                "passed" => isapprox(observed, expected;
                    rtol=1e-13, atol=1e-13),
            ))
        end
        push!(curves, Dict(
            "id" => id,
            "dead_time_s" => dead_time,
            "points" => points,
            "all_points_passed" => all(point["passed"] for point in points),
        ))
    end
    return curves
end

function spad_deterministic_contract()
    radiometry = SPADArrayDetector((2, 8); integration_time=2.0,
        noise=NoiseNone(), gate_model=DutyCycleGate(0.25),
        sensor=SPADArraySensor(active_area_detection_efficiency=0.5,
            fill_factor=0.8))
    radiometry_output = copy(capture!(radiometry,
        fill(10.0, 2, 8), Xoshiro(9180)))

    afterpulse = SPADArrayDetector((2, 8); noise=NoiseNone(),
        sensor=SPADArraySensor(active_area_detection_efficiency=1.0,
            mean_response_model=FirstOrderAfterpulseMeanResponse(1.5)))
    afterpulse_output = copy(capture!(afterpulse,
        fill(4.0, 2, 8), Xoshiro(9181)))

    redistribution = SPADArrayDetector((3, 3); noise=NoiseNone(),
        sensor=SPADArraySensor(active_area_detection_efficiency=1.0,
            mean_response_model=NearestNeighborCountRedistribution(0.4)))
    impulse = zeros(3, 3)
    impulse[2, 2] = 10.0
    redistribution_output = copy(capture!(redistribution,
        impulse, Xoshiro(9182)))

    fixed_storage = redistribution.state.counts
    before_mismatch = copy(fixed_storage)
    mismatch_rejected = try
        capture!(redistribution, zeros(2, 2), Xoshiro(9183))
        false
    catch error
        error isa DimensionMismatchError
    end
    mismatch_preserved = redistribution.state.counts === fixed_storage &&
        redistribution.state.counts == before_mismatch

    invalid_input_rejected = all((fill(-1.0, 3, 3), fill(Inf, 3, 3),
        fill(NaN, 3, 3))) do input
        try
            capture!(redistribution, input, Xoshiro(9184))
            false
        catch error
            error isa InvalidConfiguration
        end
    end

    function replay_detector()
        return SPADArrayDetector((16, 16); integration_time=0.5,
            noise=NoisePhoton(), gate_model=DutyCycleGate(0.75),
            sensor=SPADArraySensor(active_area_detection_efficiency=0.7,
                dark_count_rate=0.25,
                dead_time_model=NonParalyzableDeadTime(1e-3),
                mean_response_model=FirstOrderAfterpulseMeanResponse(0.1)))
    end
    replay_input = fill(40.0, 16, 16)
    replay_passed =
        capture!(replay_detector(), replay_input, Xoshiro(9185)) ==
        capture!(replay_detector(), replay_input, Xoshiro(9185))

    allocation_detector = SPADArrayDetector((8, 8); noise=NoiseNone(),
        sensor=SPADArraySensor(active_area_detection_efficiency=0.5,
            mean_response_model=NearestNeighborCountRedistribution(0.1)))
    allocation_input = fill(2.0, 8, 8)
    allocation_rng = Xoshiro(9186)
    capture!(allocation_detector, allocation_input, allocation_rng)
    steady_alloc_bytes = @allocated capture!(allocation_detector,
        allocation_input, allocation_rng)

    metadata = detector_export_metadata(afterpulse)
    return Dict{String,Any}(
        "exact_radiometry_and_gate_passed" =>
            radiometry_output == fill(2.0, 2, 8),
        "first_order_afterpulse_mean_passed" =>
            afterpulse_output == fill(10.0, 2, 8),
        "afterpulse_metadata_passed" =>
            metadata.mean_response_model ==
                :first_order_afterpulse_mean_response &&
            metadata.mean_afterpulses_per_detection == 1.5,
        "redistribution_center_passed" => redistribution_output ==
            [0.0 1.0 0.0; 1.0 6.0 1.0; 0.0 1.0 0.0],
        "redistribution_conserves_counts" =>
            sum(redistribution_output) == sum(impulse),
        "fixed_shape_mismatch_rejected" => mismatch_rejected,
        "fixed_shape_storage_preserved" => mismatch_preserved,
        "invalid_input_rejected" => invalid_input_rejected,
        "detector_mtf_not_applicable" =>
            !applicable(detector_mtf, radiometry, 0.1, 0.1),
        "deterministic_replay_passed" => replay_passed,
        "steady_alloc_bytes" => steady_alloc_bytes,
        "allocation_gate_passed" => steady_alloc_bytes == 0,
    )
end

spad_git_revision() = try
    readchomp(`git rev-parse HEAD`)
catch
    "unknown"
end

spad_git_dirty() = !success(`git diff --quiet HEAD`)

function generate_spad_qualification_artifact()
    curves = spad_dead_time_curves()
    poisson_cases = spad_poisson_cases()
    deterministic = spad_deterministic_contract()
    deterministic_gates = (
        "exact_radiometry_and_gate_passed",
        "first_order_afterpulse_mean_passed",
        "afterpulse_metadata_passed",
        "redistribution_center_passed",
        "redistribution_conserves_counts",
        "fixed_shape_mismatch_rejected",
        "fixed_shape_storage_preserved",
        "invalid_input_rejected",
        "detector_mtf_not_applicable",
        "deterministic_replay_passed",
        "allocation_gate_passed",
    )
    all_gates_passed =
        all(curve["all_points_passed"] for curve in curves) &&
        all(case["mean_passed"] && case["variance_passed"]
            for case in poisson_cases) &&
        all(deterministic[key] for key in deterministic_gates)
    artifact = Dict{String,Any}(
        "schema_version" => 1,
        "artifact_id" => "DET-SPAD-QUAL-2026-08-01",
        "family" => "spad_accumulated_count_array",
        "all_gates_passed" => all_gates_passed,
        "environment" => Dict(
            "timestamp_utc" => string(now(UTC)),
            "julia_version" => string(VERSION),
            "julia_threads" => Threads.nthreads(),
            "architecture" => string(Sys.ARCH),
            "kernel" => string(Sys.KERNEL),
            "source_revision" => spad_git_revision(),
            "source_dirty" => spad_git_dirty(),
        ),
        "model" => Dict(
            "input" => "cell-integrated photon-arrival rate per array cell",
            "output" => "accumulated expected-count or sampled-count image",
            "pipeline_order" => "radiometry and dark expectation; dead-time mean law; deterministic mean response; optional Poisson surrogate; output conversion",
            "nonparalyzable_mean_law" => "mu / (1 + mu * tau / live_time)",
            "paralyzable_mean_law" => "mu * exp(-mu * tau / live_time)",
            "statistical_scope" => "Poisson draw from adjusted mean; not the exact dead-time or afterpulse count distribution",
            "scientific_references" => [
                "https://doi.org/10.1109/TCOMM.2018.2822815",
                "https://doi.org/10.1098/rsta.2019.0194",
                "https://doi.org/10.1364/AO.35.001956",
            ],
        ),
        "qualification" => Dict(
            "samples_per_moment_case" => SPAD_SAMPLE_COUNT,
            "sigma_limit" => SPAD_SIGMA_LIMIT,
            "dead_time_curves" => curves,
            "poisson_surrogate_cases" => poisson_cases,
            "deterministic" => deterministic,
        ),
        "scope" => Dict(
            "included" => [
                "fixed array dimensions and preallocated storage",
                "active-area detection efficiency, fill factor, live time, and dark counts",
                "nonparalyzable and paralyzable expected-count dead-time laws",
                "first-order afterpulse mean scaling",
                "count-conserving nearest-neighbor redistribution",
                "Poisson surrogate moments after mean-response adjustment",
                "deterministic replay and warmed CPU allocation",
            ],
            "excluded" => [
                "photon-event timestamps, avalanche waveforms, and energy resolution",
                "exact dead-time and afterpulse count distributions",
                "history-dependent or higher-generation afterpulsing",
                "optical blur, pixel-aperture response, and detector MTF",
                "spectral detection-efficiency curves and calibrated detector maps",
                "named detector modules and camera profiles",
            ],
        ),
    )
    all_gates_passed ||
        error("SPAD qualification artifact contains a failed gate")
    open(SPAD_ARTIFACT_PATH, "w") do io
        TOML.print(io, artifact; sorted=true)
    end
    return SPAD_ARTIFACT_PATH
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    println(generate_spad_qualification_artifact())
end

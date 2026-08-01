using AdaptiveOpticsSim
using AdaptiveOpticsSim.Detectors
using Dates
using Random
using Statistics
using TOML

const LINEAR_APD_ARTIFACT_PATH = joinpath(
    @__DIR__, "..", "benchmarks", "results", "detectors",
    "2026-07-31-linear-apd-qualification.toml")
const LINEAR_APD_SAMPLE_COUNT = 16_384
const LINEAR_APD_SIGMA_LIMIT = 6.0

function linear_apd_moment_case(id, detector, input, expected_mean,
    expected_variance, seed)
    samples = copy(capture!(detector, input, Xoshiro(seed)))
    sample_count = length(samples)
    mean_limit = LINEAR_APD_SIGMA_LIMIT *
        sqrt(expected_variance / sample_count)
    variance_limit = LINEAR_APD_SIGMA_LIMIT * expected_variance *
        sqrt(2 / (sample_count - 1))
    observed_mean = mean(samples)
    observed_variance = var(samples)
    return Dict{String,Any}(
        "id" => id,
        "seed" => seed,
        "sample_count" => sample_count,
        "expected_mean" => expected_mean,
        "expected_variance" => expected_variance,
        "observed_mean" => observed_mean,
        "observed_variance" => observed_variance,
        "mean_absolute_limit" => mean_limit,
        "variance_absolute_limit" => variance_limit,
        "mean_passed" =>
            abs(observed_mean - expected_mean) <= mean_limit,
        "variance_passed" =>
            abs(observed_variance - expected_variance) <= variance_limit,
    )
end

function linear_apd_moment_cases()
    n = LINEAR_APD_SAMPLE_COUNT
    multiplied_shot = LinearAPDDetector(
        topology=LinearAPDChannelBank(n), integration_time=1.0,
        qe=1.0, avalanche_gain=3.0, excess_noise_factor=1.4,
        conversion_gain=2.0, noise=NoisePhoton())
    multiplication_only = LinearAPDDetector(
        topology=LinearAPDChannelBank(n), avalanche_gain=3.0,
        excess_noise_factor=1.4, noise=NoiseNone())
    read_only = LinearAPDDetector(
        topology=LinearAPDChannelBank(n), avalanche_gain=3.0,
        conversion_gain=2.0, noise=NoiseReadout(2.0))
    return [
        linear_apd_moment_case("multiplied_shot", multiplied_shot,
            fill(80.0, n), 480.0, 4032.0, 9040),
        linear_apd_moment_case("multiplication_only",
            multiplication_only, fill(80.0, n), 240.0, 288.0, 9041),
        linear_apd_moment_case("read_only", read_only,
            zeros(n), 0.0, 16.0, 9042),
    ]
end

function linear_apd_deterministic_contract()
    single = LinearAPDDetector(integration_time=0.5, qe=0.5,
        avalanche_gain=4.0, conversion_gain=2.0, noise=NoiseNone())
    single_output = copy(capture!(single, 10.0, Xoshiro(9050)))

    bank = LinearAPDDetector(topology=LinearAPDChannelBank(4),
        integration_time=2.0, qe=0.5, avalanche_gain=4.0,
        dark_current=1.0, conversion_gain=2.0, noise=NoiseNone())
    bank_output = copy(capture!(bank, fill(3.0, 4), Xoshiro(9051)))
    metadata = detector_export_metadata(bank)

    function replay_detector()
        return LinearAPDDetector(topology=LinearAPDChannelBank(64),
            integration_time=0.5, qe=0.75, avalanche_gain=5.0,
            excess_noise_factor=1.2, dark_current=0.5,
            conversion_gain=2.0, noise=NoisePhotonReadout(0.25))
    end
    replay_input = fill(40.0, 64)
    replay_passed =
        capture!(replay_detector(), replay_input, Xoshiro(9052)) ==
        capture!(replay_detector(), replay_input, Xoshiro(9052))

    allocation_detector = LinearAPDDetector(
        topology=LinearAPDChannelBank(64), noise=NoiseNone())
    allocation_input = fill(2.0, 64)
    allocation_rng = Xoshiro(9053)
    capture!(allocation_detector, allocation_input, allocation_rng)
    steady_alloc_bytes = @allocated capture!(allocation_detector,
        allocation_input, allocation_rng)

    return Dict{String,Any}(
        "ambiguous_generic_apd_removed" =>
            !isdefined(Detectors, :APDDetector) &&
            !isdefined(Detectors, :APDSensor),
        "single_element_vector_storage" =>
            single_output == [20.0] && ndims(channel_output(single)) == 1,
        "exact_signal_order_passed" => bank_output == fill(40.0, 4),
        "channel_bank_vector_storage" =>
            ndims(channel_output(bank)) == 1 && length(channel_output(bank)) == 4,
        "topology_metadata_passed" =>
            metadata.topology == :channel_bank && metadata.n_channels == 4,
        "matrix_input_rejected" =>
            !applicable(capture!, bank, fill(3.0, 2, 2)),
        "deterministic_replay_passed" => replay_passed,
        "steady_alloc_bytes" => steady_alloc_bytes,
        "allocation_gate_passed" => steady_alloc_bytes == 0,
    )
end

linear_apd_git_revision() = try
    readchomp(`git rev-parse HEAD`)
catch
    "unknown"
end

linear_apd_git_dirty() = !success(`git diff --quiet HEAD`)

function generate_linear_apd_qualification_artifact()
    moment_cases = linear_apd_moment_cases()
    deterministic = linear_apd_deterministic_contract()
    deterministic_keys = (
        "ambiguous_generic_apd_removed",
        "single_element_vector_storage",
        "exact_signal_order_passed",
        "channel_bank_vector_storage",
        "topology_metadata_passed",
        "matrix_input_rejected",
        "deterministic_replay_passed",
        "allocation_gate_passed",
    )
    all_gates_passed = all(case ->
            case["mean_passed"] && case["variance_passed"], moment_cases) &&
        all(deterministic[key] for key in deterministic_keys)
    artifact = Dict{String,Any}(
        "schema_version" => 1,
        "artifact_id" => "DET-LINEAR-APD-QUAL-2026-07-31",
        "family" => "linear_mode_apd_channels",
        "all_gates_passed" => all_gates_passed,
        "environment" => Dict(
            "timestamp_utc" => string(now(UTC)),
            "julia_version" => string(VERSION),
            "julia_threads" => Threads.nthreads(),
            "architecture" => string(Sys.ARCH),
            "kernel" => string(Sys.KERNEL),
            "source_revision" => linear_apd_git_revision(),
            "source_dirty" => linear_apd_git_dirty(),
        ),
        "model" => Dict(
            "operating_regime" => "linear mode",
            "topology" => "single element or fixed channel bank",
            "input" => "photon flux per channel per second",
            "output" => "one-dimensional analog channel values",
            "excess_noise_factor_definition" => "F = E[M^2] / E[M]^2",
            "multiplication_approximation" =>
                "q' = max(q + sqrt((F - 1) * q) * Z, 0)",
            "qualified_moderate_charge_condition" => "q / (F - 1) >= 25",
            "scientific_references" => [
                "https://doi.org/10.1109/16.47763",
                "https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1065.pdf",
            ],
        ),
        "qualification" => Dict(
            "samples_per_case" => LINEAR_APD_SAMPLE_COUNT,
            "sigma_limit" => LINEAR_APD_SIGMA_LIMIT,
            "moment_cases" => moment_cases,
            "deterministic" => deterministic,
        ),
        "scope" => Dict(
            "included" => [
                "quantum efficiency and integration duration",
                "dark current",
                "linear avalanche gain and moderate-charge excess-noise approximation",
                "additive read noise and conversion gain",
                "single-element and fixed channel-bank topology",
                "deterministic replay and warmed allocation",
            ],
            "excluded" => [
                "Geiger-mode operation and photon timestamps",
                "bandwidth-dependent analog electronics",
                "transimpedance-amplifier response",
                "correlated channel noise",
                "calibrated gain-bandwidth products",
                "named detector modules and camera profiles",
            ],
        ),
    )
    all_gates_passed ||
        error("linear-mode APD qualification artifact contains a failed gate")
    open(LINEAR_APD_ARTIFACT_PATH, "w") do io
        TOML.print(io, artifact; sorted=true)
    end
    return LINEAR_APD_ARTIFACT_PATH
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    println(generate_linear_apd_qualification_artifact())
end

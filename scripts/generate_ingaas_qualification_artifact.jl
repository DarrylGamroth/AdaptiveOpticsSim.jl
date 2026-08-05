using AdaptiveOpticsSim
using AdaptiveOpticsSim.Detectors
using Dates
using Random
using Statistics
using TOML

const INGAAS_ARTIFACT_PATH = joinpath(
    @__DIR__, "..", "benchmarks", "results", "detectors",
    "2026-07-31-ingaas-qualification.toml")
const INGAAS_SAMPLE_COUNT = 16_384
const INGAAS_SIGMA_LIMIT = 6.0

function ingaas_poisson_case()
    expected_mean = 8.0
    detector = Detector(
        exposure_duration=2.0,
        noise=NoiseNone(),
        qe=1.0,
        dark_current=1.5,
        sensor=InGaAsSensor(glow_rate=2.5),
        response_model=NullFrameResponse(),
    )
    samples = vec(copy(capture!(
        detector, zeros(128, 128), Xoshiro(9410))))
    mean_limit = INGAAS_SIGMA_LIMIT *
        sqrt(expected_mean / length(samples))
    variance_limit = INGAAS_SIGMA_LIMIT * expected_mean *
        sqrt(2 / (length(samples) - 1))
    observed_mean = mean(samples)
    observed_variance = var(samples)
    return Dict{String,Any}(
        "id" => "combined_glow_and_dark_current",
        "seed" => 9410,
        "sample_count" => length(samples),
        "exposure_seconds" => 2.0,
        "glow_rate_electrons_per_pixel_second" => 2.5,
        "dark_current_electrons_per_pixel_second" => 1.5,
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

function ingaas_persistence_sequence(gain)
    detector = Detector(
        exposure_duration=1.0,
        noise=NoiseNone(),
        qe=1.0,
        gain=gain,
        sensor=InGaAsSensor(persistence_model=
            ExponentialPersistence(0.25, 0.5)),
        response_model=NullFrameResponse(),
    )
    illumination = fill(8.0, 4, 4)
    dark = zeros(4, 4)
    frame_means = Float64[]
    latent_means = Float64[]
    for (index, input) in enumerate((illumination, dark, dark))
        output = capture!(detector, input, Xoshiro(9420 + index))
        push!(frame_means, mean(output))
        push!(latent_means, mean(detector.state.latent_buffer))
    end
    return frame_means, latent_means
end

function ingaas_deterministic_contract()
    impulse = zeros(7, 7)
    impulse[4, 4] = 1.0
    response_detector = Detector(noise=NoiseNone(), qe=1.0,
        sensor=InGaAsSensor())
    response_output = copy(capture!(
        response_detector, impulse, Xoshiro(9430)))
    response_metadata = detector_export_metadata(response_detector)

    pipeline_detector = Detector(
        exposure_duration=2.0,
        noise=NoiseNone(),
        qe=0.5,
        gain=2.0,
        full_well=5.0,
        sensor=InGaAsSensor(),
        response_model=NullFrameResponse(),
        nonlinearity_model=SaturatingFrameNonlinearity(0.1),
    )
    pipeline_output = capture!(
        pipeline_detector, fill(20.0, 4, 4), Xoshiro(9431))

    unit_frames, unit_latent = ingaas_persistence_sequence(1.0)
    four_frames, four_latent = ingaas_persistence_sequence(4.0)
    expected_charge = [8.0, 2.0, 1.5]
    expected_latent = [2.0, 1.5, 1.125]

    replay_a = Detector(exposure_duration=0.5,
        noise=NoisePhotonReadout(0.25), qe=0.75, dark_current=0.5,
        sensor=InGaAsSensor(glow_rate=0.25,
            persistence_model=ExponentialPersistence(0.1, 0.4)),
        response_model=NullFrameResponse())
    replay_b = Detector(exposure_duration=0.5,
        noise=NoisePhotonReadout(0.25), qe=0.75, dark_current=0.5,
        sensor=InGaAsSensor(glow_rate=0.25,
            persistence_model=ExponentialPersistence(0.1, 0.4)),
        response_model=NullFrameResponse())
    replay_input = fill(4.0, 16, 16)
    replay_passed =
        capture!(replay_a, replay_input, Xoshiro(9432)) ==
        capture!(replay_b, replay_input, Xoshiro(9432)) &&
        replay_a.state.latent_buffer == replay_b.state.latent_buffer

    batched_detector = Detector(noise=NoiseNone(), qe=1.0,
        sensor=InGaAsSensor(persistence_model=
            ExponentialPersistence(0.1, 0.4)),
        response_model=NullFrameResponse())
    cube = fill(3.0, 2, 4, 4)
    original = copy(cube)
    batched_rejected = try
        AdaptiveOpticsSim.Detectors.capture_stack!(
            batched_detector, cube, similar(cube), Xoshiro(9433))
        false
    catch err
        err isa InvalidConfiguration
    end

    allocation_detector = Detector(exposure_duration=0.5,
        noise=NoiseNone(), qe=1.0,
        sensor=InGaAsSensor(persistence_model=
            ExponentialPersistence(0.1, 0.4)),
        response_model=NullFrameResponse())
    allocation_input = fill(2.0, 16, 16)
    allocation_rng = Xoshiro(9434)
    capture!(allocation_detector, allocation_input, allocation_rng)
    steady_alloc_bytes = @allocated capture!(
        allocation_detector, allocation_input, allocation_rng)

    return Dict{String,Any}(
        "default_response_is_null" =>
            response_metadata.frame_response == :none,
        "default_response_preserves_impulse" => response_output == impulse,
        "default_mtf_is_not_claimed" =>
            !AdaptiveOpticsSim.Detectors.supports_detector_mtf(
                response_detector),
        "pipeline_order_passed" => pipeline_output == fill(10.0, 4, 4),
        "persistence_charge_means" => unit_frames,
        "persistence_latent_means" => unit_latent,
        "persistence_recurrence_passed" =>
            unit_frames == expected_charge &&
            unit_latent == expected_latent,
        "persistence_gain_independence_passed" =>
            four_frames == 4 .* unit_frames &&
            four_latent == unit_latent,
        "deterministic_replay_passed" => replay_passed,
        "batched_persistence_rejected" => batched_rejected,
        "batched_input_unmodified" => cube == original,
        "steady_alloc_bytes" => steady_alloc_bytes,
        "allocation_gate_passed" => steady_alloc_bytes == 0,
    )
end

ingaas_git_revision() = try
    readchomp(`git rev-parse HEAD`)
catch
    "unknown"
end

ingaas_git_dirty() = !success(`git diff --quiet HEAD`)

function generate_ingaas_qualification_artifact()
    poisson_case = ingaas_poisson_case()
    deterministic = ingaas_deterministic_contract()
    deterministic_keys = (
        "default_response_is_null",
        "default_response_preserves_impulse",
        "default_mtf_is_not_claimed",
        "pipeline_order_passed",
        "persistence_recurrence_passed",
        "persistence_gain_independence_passed",
        "deterministic_replay_passed",
        "batched_persistence_rejected",
        "batched_input_unmodified",
        "allocation_gate_passed",
    )
    all_gates_passed =
        poisson_case["mean_passed"] && poisson_case["variance_passed"] &&
        all(deterministic[key] for key in deterministic_keys)
    artifact = Dict{String,Any}(
        "schema_version" => 1,
        "artifact_id" => "DET-INGAAS-QUAL-2026-07-31",
        "family" => "ingaas_frame_detector",
        "all_gates_passed" => all_gates_passed,
        "environment" => Dict(
            "timestamp_utc" => string(now(UTC)),
            "julia_version" => string(VERSION),
            "julia_threads" => Threads.nthreads(),
            "architecture" => string(Sys.ARCH),
            "kernel" => string(Sys.KERNEL),
            "source_revision" => ingaas_git_revision(),
            "source_dirty" => ingaas_git_dirty(),
        ),
        "model" => Dict(
            "default_response" => "none",
            "persistence_recurrence" =>
                "latent[k+1] = decay * latent[k] + coupling * charge[k]",
            "persistence_update_boundary" =>
                "after charge transport and before readout electronics",
            "coefficient_domain" => "dimensionless per completed frame",
            "scientific_references" => [
                "https://arxiv.org/abs/1307.1469",
                "https://arxiv.org/abs/1406.2695",
            ],
        ),
        "qualification" => Dict(
            "samples_per_case" => INGAAS_SAMPLE_COUNT,
            "sigma_limit" => INGAAS_SIGMA_LIMIT,
            "poisson_case" => poisson_case,
            "deterministic" => deterministic,
        ),
        "scope" => Dict(
            "included" => [
                "explicit shared frame pipeline effects",
                "independent Poisson glow and dark-current charge",
                "null default response and explicit configured response",
                "discrete charge-domain exponential persistence recurrence",
                "deterministic replay and bounded persistent state",
            ],
            "excluded" => [
                "calibrated wavelength-dependent material response",
                "trap populations or elapsed-time persistence constants",
                "afterglow spectra and exposure-dependent persistence",
                "correlated 1/f or burst noise",
                "universal InGaAs MTF and camera profiles",
            ],
        ),
    )
    all_gates_passed ||
        error("InGaAs qualification artifact contains a failed gate")
    open(INGAAS_ARTIFACT_PATH, "w") do io
        TOML.print(io, artifact; sorted=true)
    end
    return INGAAS_ARTIFACT_PATH
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    println(generate_ingaas_qualification_artifact())
end

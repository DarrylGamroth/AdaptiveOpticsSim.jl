using AdaptiveOpticsSim
using AdaptiveOpticsSim.Detectors
using Dates
using Random
using Statistics
using TOML

const CCD_ARTIFACT_DIR =
    joinpath(@__DIR__, "..", "benchmarks", "results", "detectors")
const CCD_ARTIFACT_PATH = joinpath(CCD_ARTIFACT_DIR,
    "2026-07-30-ccd-single-read-qualification.toml")
const CCD_ARTIFACT_ID = "DET-CCD-QUAL-2026-07-30"
const CCD_FRAME_COUNT = 16
const CCD_FRAME_SIZE = 32
const CCD_SIGMA_LIMIT = 6.0

function ccd_artifact_samples(detector, input, rng)
    samples = Vector{Float64}(undef,
        CCD_FRAME_COUNT * length(input))
    offset = 0
    for _ in 1:CCD_FRAME_COUNT
        frame = capture!(detector, input, rng)
        copyto!(samples, offset + 1, vec(frame), 1, length(frame))
        offset += length(frame)
    end
    return samples
end

function ccd_artifact_moment_case(id, detector, input, seed,
    expected_mean, expected_variance, expected_fourth_central)
    samples = ccd_artifact_samples(detector, input, Xoshiro(seed))
    sample_count = length(samples)
    mean_limit = CCD_SIGMA_LIMIT *
        sqrt(expected_variance / sample_count)
    variance_se = sqrt((expected_fourth_central -
        ((sample_count - 3) / (sample_count - 1)) *
        expected_variance^2) / sample_count)
    variance_limit = CCD_SIGMA_LIMIT * variance_se
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

function ccd_artifact_poisson_case(id, detector, input, seed,
    expected_mean)
    fourth_central = expected_mean + 3 * expected_mean^2
    return ccd_artifact_moment_case(id, detector, input, seed,
        expected_mean, expected_mean, fourth_central)
end

function ccd_artifact_cases()
    zero_input = zeros(CCD_FRAME_SIZE, CCD_FRAME_SIZE)
    signal_input = fill(50.0, CCD_FRAME_SIZE, CCD_FRAME_SIZE)

    cases = Dict{String,Any}[]
    push!(cases, ccd_artifact_poisson_case(
        "shot",
        Detector(exposure_duration=0.5, qe=0.8, noise=NoisePhoton(),
            sensor=CCDSensor()),
        signal_input, 3101, 20.0,
    ))
    push!(cases, ccd_artifact_poisson_case(
        "dark",
        Detector(exposure_duration=2.0, qe=1.0, dark_current=12.0,
            noise=NoiseNone(), sensor=CCDSensor()),
        zero_input, 3102, 24.0,
    ))
    push!(cases, ccd_artifact_poisson_case(
        "clock_induced_charge",
        Detector(exposure_duration=7.0, qe=1.0, noise=NoiseNone(),
            sensor=CCDSensor(clock_induced_charge_per_frame=3.5)),
        zero_input, 3103, 3.5,
    ))

    read_sigma = 3.0
    read_variance = read_sigma^2
    push!(cases, ccd_artifact_moment_case(
        "read_noise",
        Detector(exposure_duration=1.0, qe=1.0,
            noise=NoiseReadout(read_sigma), sensor=CCDSensor()),
        zero_input, 3104, 0.0, read_variance,
        3 * read_variance^2,
    ))

    combined_poisson_mean = 20.0 + 6.0 + 2.0
    combined_sigma = 2.5
    combined_variance = combined_poisson_mean + combined_sigma^2
    push!(cases, ccd_artifact_moment_case(
        "combined",
        Detector(exposure_duration=0.5, qe=0.8, dark_current=12.0,
            noise=NoisePhotonReadout(combined_sigma),
            sensor=CCDSensor(clock_induced_charge_per_frame=2.0)),
        signal_input, 3105, combined_poisson_mean, combined_variance,
        combined_poisson_mean + 3 * combined_variance^2,
    ))
    return cases
end

function ccd_artifact_deterministic_contract()
    zero_input = zeros(16, 16)
    short = Detector(exposure_duration=0.1, noise=NoiseNone(), qe=1.0,
        sensor=CCDSensor(clock_induced_charge_per_frame=3.0))
    long = Detector(exposure_duration=10.0, noise=NoiseNone(), qe=1.0,
        sensor=CCDSensor(clock_induced_charge_per_frame=3.0))
    cic_exposure_invariant =
        capture!(short, zero_input, Xoshiro(3201)) ==
        capture!(long, zero_input, Xoshiro(3201))

    read_duration_rejected = try
        CCDSensor(sample_duration=1e-3)
        false
    catch err
        err isa InvalidConfiguration
    end

    allocation_detector = Detector(exposure_duration=1.0,
        noise=NoiseNone(), qe=1.0,
        sensor=CCDSensor(clock_induced_charge_per_frame=0.75))
    allocation_rng = Xoshiro(3202)
    capture!(allocation_detector, zero_input, allocation_rng)
    steady_alloc_bytes =
        @allocated capture!(allocation_detector, zero_input, allocation_rng)
    metadata = detector_export_metadata(allocation_detector)

    return Dict{String,Any}(
        "cic_exposure_invariant" => cic_exposure_invariant,
        "single_read_read_duration_rejected" => read_duration_rejected,
        "default_response_is_null" => metadata.frame_response == :none,
        "default_supports_mtf" =>
            AdaptiveOpticsSim.Detectors.supports_detector_mtf(
                allocation_detector),
        "sampling_mode" => String(metadata.sampling_mode),
        "sampling_read_duration_is_absent" =>
            isnothing(metadata.sampling_read_duration),
        "steady_alloc_bytes" => steady_alloc_bytes,
        "allocation_gate_passed" => steady_alloc_bytes == 0,
    )
end

function ccd_artifact_git_revision()
    return try
        readchomp(`git rev-parse HEAD`)
    catch
        "unknown"
    end
end

function ccd_artifact_git_dirty()
    return !success(`git diff --quiet HEAD`)
end

function generate_ccd_qualification_artifact()
    cases = ccd_artifact_cases()
    deterministic = ccd_artifact_deterministic_contract()
    all_moments_passed = all(case -> case["mean_passed"] &&
        case["variance_passed"], cases)
    all_deterministic_passed =
        deterministic["cic_exposure_invariant"] &&
        deterministic["single_read_read_duration_rejected"] &&
        deterministic["default_response_is_null"] &&
        !deterministic["default_supports_mtf"] &&
        deterministic["sampling_read_duration_is_absent"] &&
        deterministic["allocation_gate_passed"]

    cpu = first(Sys.cpu_info())
    artifact = Dict{String,Any}(
        "schema_version" => 1,
        "artifact_id" => CCD_ARTIFACT_ID,
        "family" => "conventional_ccd_single_read",
        "all_gates_passed" =>
            all_moments_passed && all_deterministic_passed,
        "qualification" => Dict{String,Any}(
            "frame_count_per_case" => CCD_FRAME_COUNT,
            "frame_size" => [CCD_FRAME_SIZE, CCD_FRAME_SIZE],
            "samples_per_case" =>
                CCD_FRAME_COUNT * CCD_FRAME_SIZE^2,
            "sigma_limit" => CCD_SIGMA_LIMIT,
            "moment_cases" => cases,
            "deterministic" => deterministic,
        ),
        "scope" => Dict{String,Any}(
            "included" => [
                "single-read CCD shared frame pipeline",
                "independent Poisson clock-induced charge per pixel per frame",
                "shot, dark, CIC, and independent Gaussian read-noise moments",
                "Plant-owned readout/readiness timing boundary",
            ],
            "excluded" => [
                "Skipper sampling",
                "charge-transfer inefficiency",
                "serial or parallel register geometry",
                "blooming and transfer smear",
                "radiation damage",
                "camera or vendor profiles",
            ],
        ),
        "environment" => Dict{String,Any}(
            "timestamp_utc" => string(Dates.now(Dates.UTC)),
            "source_revision" => ccd_artifact_git_revision(),
            "source_dirty" => ccd_artifact_git_dirty(),
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

    mkpath(CCD_ARTIFACT_DIR)
    open(CCD_ARTIFACT_PATH, "w") do io
        TOML.print(io, artifact; sorted=true)
    end
    artifact["all_gates_passed"] || error(
        "CCD qualification artifact contains a failed gate")
    println("wrote ", CCD_ARTIFACT_PATH)
    return artifact
end

generate_ccd_qualification_artifact()

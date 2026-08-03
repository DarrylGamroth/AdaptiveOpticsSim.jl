using AdaptiveOpticsSim
using AdaptiveOpticsSim.Backends
using AdaptiveOpticsSim.Detectors
using AdaptiveOpticsSim.Optics
import AdaptiveOpticsSim.Detectors:
    ClippedGaussianMultiplicationApproximation,
    ConditionalGammaMultiplication
using Dates
using Random
using Statistics
using TOML

const EMCCD_ARTIFACT_DIR =
    joinpath(@__DIR__, "..", "benchmarks", "results", "detectors")
const EMCCD_ARTIFACT_PATH = joinpath(EMCCD_ARTIFACT_DIR,
    "2026-07-30-emccd-qualification.toml")
const EMCCD_ARTIFACT_ID = "DET-EMCCD-QUAL-2026-07-30"
const EMCCD_FRAME_COUNT = 16
const EMCCD_FRAME_SIZE = 32
const EMCCD_SIGMA_LIMIT = 6.0

function emccd_artifact_samples(detector, input, seed)
    samples = Vector{Float64}(undef,
        EMCCD_FRAME_COUNT * length(input))
    rng = Xoshiro(seed)
    offset = 0
    for _ in 1:EMCCD_FRAME_COUNT
        frame = capture!(detector, input, rng)
        copyto!(samples, offset + 1, vec(frame), 1, length(frame))
        offset += length(frame)
    end
    return samples
end

function emccd_artifact_moment_record(id, samples, seed,
    expected_mean, expected_variance, expected_fourth_central)
    sample_count = length(samples)
    mean_limit = EMCCD_SIGMA_LIMIT *
        sqrt(expected_variance / sample_count)
    variance_se = sqrt((expected_fourth_central -
        ((sample_count - 3) / (sample_count - 1)) *
        expected_variance^2) / sample_count)
    variance_limit = EMCCD_SIGMA_LIMIT * variance_se
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

function emccd_artifact_cdf_record(id, samples, threshold,
    expected_probability)
    sample_count = length(samples)
    observed_probability =
        count(value -> value <= threshold, samples) / sample_count
    probability_limit = EMCCD_SIGMA_LIMIT * sqrt(
        expected_probability * (1 - expected_probability) / sample_count)
    return Dict{String,Any}(
        "id" => id,
        "sample_count" => sample_count,
        "threshold" => threshold,
        "expected_probability" => expected_probability,
        "observed_probability" => observed_probability,
        "probability_absolute_limit" => probability_limit,
        "passed" => abs(observed_probability - expected_probability) <=
            probability_limit,
    )
end

function emccd_artifact_multiplication_cases()
    n = EMCCD_FRAME_SIZE
    moments = Dict{String,Any}[]
    cdfs = Dict{String,Any}[]

    exponential_gain = 12.0
    exponential_seed = 3401
    exponential = emccd_artifact_samples(
        Detector(integration_time=1.0, qe=1.0, noise=NoiseNone(),
            gain=exponential_gain,
            sensor=EMCCDSensor(excess_noise_factor=sqrt(2.0),
                multiplication_model=ConditionalGammaMultiplication())),
        ones(n, n), exponential_seed)
    push!(moments, emccd_artifact_moment_record(
        "conditional_gamma_exponential", exponential, exponential_seed,
        exponential_gain, exponential_gain^2,
        9 * exponential_gain^4))
    push!(cdfs, emccd_artifact_cdf_record(
        "conditional_gamma_exponential_median", exponential,
        exponential_gain * log(2), 0.5))
    push!(cdfs, emccd_artifact_cdf_record(
        "conditional_gamma_exponential_p90", exponential,
        exponential_gain * log(10), 0.9))

    erlang_charge = 8.0
    erlang_gain = 5.0
    erlang_seed = 3402
    erlang = emccd_artifact_samples(
        Detector(integration_time=1.0, qe=1.0, noise=NoiseNone(),
            gain=erlang_gain,
            sensor=EMCCDSensor(excess_noise_factor=sqrt(2.0),
                multiplication_model=ConditionalGammaMultiplication())),
        fill(erlang_charge, n, n), erlang_seed)
    erlang_variance = erlang_charge * erlang_gain^2
    push!(moments, emccd_artifact_moment_record(
        "conditional_gamma_erlang", erlang, erlang_seed,
        erlang_charge * erlang_gain, erlang_variance,
        3 * erlang_charge * (erlang_charge + 2) * erlang_gain^4))
    erlang_mean_cdf = 1 - exp(-erlang_charge) *
        sum(erlang_charge^order / factorial(order)
            for order in 0:(Int(erlang_charge) - 1))
    push!(cdfs, emccd_artifact_cdf_record(
        "conditional_gamma_erlang_mean", erlang,
        erlang_charge * erlang_gain, erlang_mean_cdf))

    fractional_charge = 10.0
    fractional_gain = 5.0
    fractional_factor2 = 1.4^2 - 1
    fractional_shape = fractional_charge / fractional_factor2
    fractional_variance =
        fractional_charge * fractional_gain^2 * fractional_factor2
    fractional_seed = 3403
    fractional = emccd_artifact_samples(
        Detector(integration_time=1.0, qe=1.0, noise=NoiseNone(),
            gain=fractional_gain,
            sensor=EMCCDSensor(excess_noise_factor=1.4,
                multiplication_model=ConditionalGammaMultiplication())),
        fill(fractional_charge, n, n), fractional_seed)
    push!(moments, emccd_artifact_moment_record(
        "conditional_gamma_fractional_shape", fractional,
        fractional_seed, fractional_charge * fractional_gain,
        fractional_variance,
        (3 + 6 / fractional_shape) * fractional_variance^2))

    approximate_charge = 64.0
    approximate_gain = 5.0
    approximate_variance = approximate_charge * approximate_gain^2
    approximate_seed = 3404
    approximate = emccd_artifact_samples(
        Detector(integration_time=1.0, qe=1.0, noise=NoiseNone(),
            gain=approximate_gain,
            sensor=EMCCDSensor(excess_noise_factor=sqrt(2.0),
                multiplication_model=
                    ClippedGaussianMultiplicationApproximation())),
        fill(approximate_charge, n, n), approximate_seed)
    approximation_record = emccd_artifact_moment_record(
        "clipped_gaussian_moderate_charge", approximate,
        approximate_seed, approximate_charge * approximate_gain,
        approximate_variance, 3 * approximate_variance^2)
    approximation_record["all_nonnegative"] = all(>=(0), approximate)
    push!(moments, approximation_record)

    cic_mean = 8.0
    cic_gain = 5.0
    cic_variance = 2 * cic_mean * cic_gain^2
    cic_seed = 3405
    cic = emccd_artifact_samples(
        Detector(integration_time=0.25, qe=1.0, noise=NoiseNone(),
            gain=cic_gain,
            sensor=EMCCDSensor(excess_noise_factor=sqrt(2.0),
                clock_induced_charge_per_frame=cic_mean,
                multiplication_model=ConditionalGammaMultiplication())),
        zeros(n, n), cic_seed)
    push!(moments, emccd_artifact_moment_record(
        "cic_before_conditional_gamma", cic, cic_seed,
        cic_mean * cic_gain, cic_variance,
        24 * cic_mean * cic_gain^4 + 3 * cic_variance^2))

    return moments, cdfs
end

function emccd_artifact_photon_counting_cases()
    cases = Dict{String,Any}[]
    efficiency = 0.75
    for (case_index, incident_mean) in enumerate((0.05, 0.5, 2.0))
        seed = 3500 + case_index
        samples = emccd_artifact_samples(
            Detector(integration_time=1.0, qe=1.0,
                noise=NoisePhoton(), gain=10.0,
                sensor=EMCCDSensor(excess_noise_factor=1.0,
                    multiplication_model=
                        ClippedGaussianMultiplicationApproximation(),
                    operating_mode=PhotonCountingEMMode(
                        threshold=5.0,
                        detection_efficiency=efficiency))),
            fill(incident_mean, EMCCD_FRAME_SIZE, EMCCD_FRAME_SIZE),
            seed)
        expected = efficiency * (1 - exp(-incident_mean))
        variance = expected * (1 - expected)
        record = emccd_artifact_moment_record(
            "binary_incident_mean_$(incident_mean)", samples, seed,
            expected, variance, variance * (1 - 3 * variance))
        record["incident_mean"] = incident_mean
        record["detection_efficiency"] = efficiency
        record["all_binary"] =
            all(value -> value == 0.0 || value == 1.0, samples)
        record["coincidence_limited"] = mean(samples) <
            efficiency * incident_mean
        push!(cases, record)
    end
    return cases
end

function emccd_artifact_intensity_map(values::AbstractMatrix{T}) where
    {T<:AbstractFloat}
    metadata = OpticalPlaneMetadata(FocalPlane(), values;
        coordinate_domain=AngularCoordinates(),
        sampling=(one(T), one(T)),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition(),
        spectral=MonochromaticChannel(T(0.55e-6)))
    return IntensityMap(metadata, values)
end

function emccd_artifact_deterministic_contract()
    zero_input = zeros(16, 16)
    short = Detector(integration_time=0.1, qe=1.0,
        noise=NoiseNone(), gain=5.0,
        sensor=EMCCDSensor(excess_noise_factor=1.0,
            clock_induced_charge_per_frame=3.0))
    long = Detector(integration_time=10.0, qe=1.0,
        noise=NoiseNone(), gain=5.0,
        sensor=EMCCDSensor(excess_noise_factor=1.0,
            clock_induced_charge_per_frame=3.0))
    cic_exposure_invariant =
        capture!(short, zero_input, Xoshiro(3601)) ==
        capture!(long, zero_input, Xoshiro(3601))

    input_limited = Detector(integration_time=1.0, qe=1.0,
        noise=NoiseNone(), gain=5.0, full_well=10.0,
        sensor=EMCCDSensor(excess_noise_factor=1.0,
            register_full_well=100.0))
    input_limit_passed = capture!(input_limited,
        fill(50.0, 4, 4), Xoshiro(3602)) == fill(50.0, 4, 4)
    register_limited = Detector(integration_time=1.0, qe=1.0,
        noise=NoiseNone(), gain=5.0, full_well=30.0,
        sensor=EMCCDSensor(excess_noise_factor=1.0,
            register_full_well=100.0))
    register_limit_passed = capture!(register_limited,
        fill(50.0, 4, 4), Xoshiro(3603)) == fill(100.0, 4, 4)

    conventional = Detector(integration_time=1.0, qe=1.0,
        noise=NoiseNone(), gain=20.0, full_well=10.0,
        sensor=EMCCDSensor(output_path=ConventionalOutput(),
            em_gain_range=(1.0, 2.0)))
    conventional_passed = capture!(conventional,
        fill(50.0, 4, 4), Xoshiro(3604)) == fill(10.0, 4, 4)

    gain_range_rejected = try
        Detector(sensor=EMCCDSensor(em_gain_range=(2.0, 8.0)),
            gain=9.0)
        false
    catch err
        err isa InvalidConfiguration
    end
    accelerator_gamma_rejected = try
        AdaptiveOpticsSim.Detectors.validate_em_gain_backend(
            EMOutput(), ConditionalGammaMultiplication(),
            CUDABackend())
        false
    catch err
        err isa InvalidConfiguration
    end

    timing_input = fill(2.0, 4, 4)
    sequential = Detector(integration_time=1.0, qe=1.0,
        noise=NoiseNone(), gain=1.0,
        sensor=EMCCDSensor(readout_rate_hz=1000.0,
            acquisition_mode=SequentialAcquisition()))
    frame_transfer = Detector(integration_time=1.0, qe=1.0,
        noise=NoiseNone(), gain=1.0,
        sensor=EMCCDSensor(readout_rate_hz=1000.0,
            acquisition_mode=FrameTransferAcquisition(
                transfer_time=0.002)))
    same_optical_output = capture!(sequential, timing_input,
        Xoshiro(3605)) == capture!(frame_transfer, timing_input,
        Xoshiro(3605))
    sequential_metadata = detector_export_metadata(sequential)
    transfer_metadata = detector_export_metadata(frame_transfer)
    timing_passed =
        sequential_metadata.sampling_wallclock_time == 1.016 &&
        sequential_metadata.steady_state_frame_period == 1.016 &&
        transfer_metadata.sampling_wallclock_time == 1.018 &&
        transfer_metadata.steady_state_frame_period == 1.002

    prepared_values = fill(8.0, 16, 16)
    prepared_map = emccd_artifact_intensity_map(prepared_values)
    allocation_detector = Detector(integration_time=0.25, qe=1.0,
        noise=NoiseNone(), gain=5.0,
        sensor=EMCCDSensor(excess_noise_factor=sqrt(2.0),
            multiplication_model=ConditionalGammaMultiplication()))
    acquisition = prepare_detector_acquisition(allocation_detector,
        prepared_map)
    allocation_rng = Xoshiro(3606)
    capture!(acquisition, allocation_rng)
    steady_alloc_bytes = @allocated capture!(acquisition, allocation_rng)

    signal = 10.0
    dark = 2.0
    cic = 1.0
    read_noise = 30.0
    gain = 100.0
    excess_noise = sqrt(2.0)
    expected_snr = signal / sqrt(excess_noise^2 *
        (signal + dark + cic) + (read_noise / gain)^2)
    snr_passed = isapprox(emccd_snr(signal; dark_electrons=dark,
        cic_electrons=cic, readout_noise=read_noise, gain,
        excess_noise_factor=excess_noise), expected_snr;
        rtol=8eps(Float64), atol=0.0)

    return Dict{String,Any}(
        "cic_exposure_invariant" => cic_exposure_invariant,
        "input_full_well_passed" => input_limit_passed,
        "register_full_well_passed" => register_limit_passed,
        "conventional_output_passed" => conventional_passed,
        "gain_range_rejected" => gain_range_rejected,
        "accelerator_conditional_gamma_rejected" =>
            accelerator_gamma_rejected,
        "frame_transfer_preserves_optical_output" =>
            same_optical_output,
        "timing_contract_passed" => timing_passed,
        "snr_contract_passed" => snr_passed,
        "steady_alloc_bytes" => steady_alloc_bytes,
        "allocation_gate_passed" => steady_alloc_bytes == 0,
    )
end

emccd_artifact_git_revision() = try
    readchomp(`git rev-parse HEAD`)
catch
    "unknown"
end

emccd_artifact_git_dirty() = !success(`git diff --quiet HEAD`)

function generate_emccd_qualification_artifact()
    moments, cdfs = emccd_artifact_multiplication_cases()
    photon_counting = emccd_artifact_photon_counting_cases()
    deterministic = emccd_artifact_deterministic_contract()
    moments_passed = all(record ->
        record["mean_passed"] && record["variance_passed"] &&
        get(record, "all_nonnegative", true), moments)
    cdfs_passed = all(record -> record["passed"], cdfs)
    photon_counting_passed = all(record ->
        record["mean_passed"] && record["variance_passed"] &&
        record["all_binary"] && record["coincidence_limited"],
        photon_counting)
    deterministic_passed = all(key -> deterministic[key], (
        "cic_exposure_invariant",
        "input_full_well_passed",
        "register_full_well_passed",
        "conventional_output_passed",
        "gain_range_rejected",
        "accelerator_conditional_gamma_rejected",
        "frame_transfer_preserves_optical_output",
        "timing_contract_passed",
        "snr_contract_passed",
        "allocation_gate_passed",
    ))

    cpu = first(Sys.cpu_info())
    artifact = Dict{String,Any}(
        "schema_version" => 1,
        "artifact_id" => EMCCD_ARTIFACT_ID,
        "family" => "emccd",
        "all_gates_passed" => moments_passed && cdfs_passed &&
            photon_counting_passed && deterministic_passed,
        "qualification" => Dict{String,Any}(
            "frame_count_per_case" => EMCCD_FRAME_COUNT,
            "frame_size" => [EMCCD_FRAME_SIZE, EMCCD_FRAME_SIZE],
            "samples_per_case" =>
                EMCCD_FRAME_COUNT * EMCCD_FRAME_SIZE^2,
            "sigma_limit" => EMCCD_SIGMA_LIMIT,
            "multiplication_moment_cases" => moments,
            "multiplication_cdf_cases" => cdfs,
            "photon_counting_cases" => photon_counting,
            "deterministic" => deterministic,
        ),
        "scope" => Dict{String,Any}(
            "included" => [
                "conventional and electron-multiplying output paths",
                "CPU conditional-Gamma multiplication",
                "backend-portable clipped-Gaussian moderate/high-charge approximation",
                "CIC before multiplication",
                "input and register-referred saturation",
                "binary coincidence-limited photon-counting output",
                "sequential and frame-transfer timing",
            ],
            "excluded" => [
                "conditional-Gamma accelerator execution",
                "low-charge distribution fidelity for the clipped-Gaussian approximation",
                "frame-transfer smear",
                "charge-transfer inefficiency",
                "stage-by-stage multiplication-register physics",
                "register aging",
                "calibrated gain laws",
                "camera or vendor profiles",
            ],
        ),
        "environment" => Dict{String,Any}(
            "timestamp_utc" => string(Dates.now(Dates.UTC)),
            "source_revision" => emccd_artifact_git_revision(),
            "source_dirty" => emccd_artifact_git_dirty(),
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

    mkpath(EMCCD_ARTIFACT_DIR)
    open(EMCCD_ARTIFACT_PATH, "w") do io
        TOML.print(io, artifact; sorted=true)
    end
    artifact["all_gates_passed"] || error(
        "EMCCD qualification artifact contains a failed gate")
    println("wrote ", EMCCD_ARTIFACT_PATH)
    return artifact
end

generate_emccd_qualification_artifact()

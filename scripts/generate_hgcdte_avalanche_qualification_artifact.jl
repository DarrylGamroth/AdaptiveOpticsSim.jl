using AdaptiveOpticsSim
using AdaptiveOpticsSim.Backends
using AdaptiveOpticsSim.Detectors
import AdaptiveOpticsSim.Detectors:
    ClippedGaussianAvalancheMultiplicationApproximation,
    ConditionalGammaAvalancheMultiplication,
    FrameWindow
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.Plant
import AdaptiveOpticsSim.Plant:
    GlobalShutterAcquisitionState,
    accumulate_exposure_interval!,
    begin_exposure!,
    close_exposure!,
    complete_readout!,
    mark_acquisition_ready!,
    prepare_global_shutter_acquisition,
    take_nondestructive_read!
using Dates
using Random
using Statistics
using TOML

const HGCDTE_AVALANCHE_ARTIFACT_DIR =
    joinpath(@__DIR__, "..", "benchmarks", "results", "detectors")
const HGCDTE_AVALANCHE_ARTIFACT_PATH =
    joinpath(HGCDTE_AVALANCHE_ARTIFACT_DIR,
        "2026-07-31-hgcdte-linear-avalanche-qualification.toml")
const HGCDTE_AVALANCHE_ARTIFACT_ID =
    "DET-HGCDTE-AVALANCHE-QUAL-2026-07-31"
const HGCDTE_AVALANCHE_FRAME_COUNT = 16
const HGCDTE_AVALANCHE_FRAME_SIZE = 32
const HGCDTE_AVALANCHE_SIGMA_LIMIT = 6.0

function hgcdte_avalanche_artifact_samples(
    detector, input, seed)
    samples = Vector{Float64}(undef,
        HGCDTE_AVALANCHE_FRAME_COUNT * length(input))
    rng = Xoshiro(seed)
    offset = 0
    for _ in 1:HGCDTE_AVALANCHE_FRAME_COUNT
        frame = capture!(detector, input, rng)
        copyto!(samples, offset + 1, vec(frame), 1, length(frame))
        offset += length(frame)
    end
    return samples
end

function hgcdte_avalanche_moment_record(
    id, model, charge, gain, excess_noise_factor, seed)
    factor_minus_one = excess_noise_factor - 1
    shape = factor_minus_one == 0 ? Inf :
        charge / factor_minus_one
    expected_mean = charge * gain
    expected_variance =
        charge * gain^2 * factor_minus_one
    expected_fourth = model == "conditional_gamma" ?
        3expected_variance^2 * (1 + 2 / shape) :
        3expected_variance^2
    multiplication_model = model == "conditional_gamma" ?
        ConditionalGammaAvalancheMultiplication() :
        ClippedGaussianAvalancheMultiplicationApproximation()
    detector = Detector(
        integration_time=1.0,
        qe=1.0,
        noise=NoiseNone(),
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(
            avalanche_gain=gain,
            excess_noise_factor=excess_noise_factor,
            multiplication_model=multiplication_model),
    )
    samples = hgcdte_avalanche_artifact_samples(
        detector,
        fill(charge, HGCDTE_AVALANCHE_FRAME_SIZE,
            HGCDTE_AVALANCHE_FRAME_SIZE),
        seed)
    sample_count = length(samples)
    mean_limit = HGCDTE_AVALANCHE_SIGMA_LIMIT *
        sqrt(expected_variance / sample_count)
    variance_se = sqrt((expected_fourth -
        ((sample_count - 3) / (sample_count - 1)) *
        expected_variance^2) / sample_count)
    variance_limit =
        HGCDTE_AVALANCHE_SIGMA_LIMIT * variance_se
    observed_mean = mean(samples)
    observed_variance = var(samples)
    return Dict{String,Any}(
        "id" => id,
        "model" => model,
        "seed" => seed,
        "sample_count" => sample_count,
        "input_charge_electrons" => charge,
        "avalanche_gain" => gain,
        "excess_noise_factor" => excess_noise_factor,
        "conditional_gamma_shape" => shape,
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
        "all_nonnegative" => all(>=(0), samples),
        "moderate_charge_approximation_regime" =>
            model != "clipped_gaussian_approximation" || shape >= 25,
    )
end

function hgcdte_avalanche_moment_cases()
    return Dict{String,Any}[
        hgcdte_avalanche_moment_record(
            "conditional_gamma_single_carrier",
            "conditional_gamma", 1.0, 10.0, 1.5, 6101),
        hgcdte_avalanche_moment_record(
            "conditional_gamma_multiple_carriers",
            "conditional_gamma", 8.0, 5.0, 1.25, 6102),
        hgcdte_avalanche_moment_record(
            "conditional_gamma_low_excess_noise",
            "conditional_gamma", 20.0, 3.0, 1.1, 6103),
        hgcdte_avalanche_moment_record(
            "clipped_gaussian_moderate_charge",
            "clipped_gaussian_approximation",
            64.0, 5.0, 1.5, 6104),
    ]
end

function hgcdte_avalanche_artifact_map(
    values::AbstractMatrix{T}) where {T<:AbstractFloat}
    metadata = OpticalPlaneMetadata(DetectorPlane(), values;
        coordinate_domain=AngularCoordinates(),
        sampling=(one(T), one(T)),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition(),
        spectral=MonochromaticChannel(T(0.8e-6)))
    return IntensityMap(metadata, values)
end

function run_hgcdte_avalanche_scheduled_ramp!(
    prepared, state, values, rng)
    start = PlantTimestamp(0)
    middle = PlantTimestamp(500_000_000)
    close = PlantTimestamp(1_000_000_000)
    begin_exposure!(prepared, state, start)
    take_nondestructive_read!(prepared, state, start, rng)
    accumulate_exposure_interval!(
        prepared, state, start, middle, rng)
    take_nondestructive_read!(prepared, state, middle, rng)
    fill!(values, zero(eltype(values)))
    accumulate_exposure_interval!(
        prepared, state, middle, close, rng)
    take_nondestructive_read!(prepared, state, close, rng)
    cube = detector_ramp_cube(prepared.detector)
    retained = true
    @inbounds for column in axes(cube, 2), row in axes(cube, 1)
        retained &= cube[row, column, 3] == cube[row, column, 2]
    end
    close_exposure!(prepared, state, close)
    complete_readout!(prepared, state, close, rng)
    mark_acquisition_ready!(prepared, state, close)
    return retained
end

function hgcdte_avalanche_scheduled_fixture(seed)
    values = fill(128.0, 8, 8)
    map = hgcdte_avalanche_artifact_map(values)
    detector = Detector(
        integration_time=1.0,
        qe=1.0,
        noise=NoiseNone(),
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(
            avalanche_gain=2.0,
            excess_noise_factor=1.5,
            multiplication_model=
                ClippedGaussianAvalancheMultiplicationApproximation(),
            read_time=0.0,
            sampling_mode=UpTheRampSampling(3)),
    )
    definition = GlobalShutterAcquisitionDefinition(
        PlantDuration(1_000_000_000))
    prepared = prepare_global_shutter_acquisition(
        detector, map, definition)
    state = GlobalShutterAcquisitionState(prepared)
    return (; detector, prepared, state, values, rng=Xoshiro(seed))
end

function hgcdte_avalanche_deterministic_contract()
    exact = Detector(
        integration_time=1.0,
        qe=1.0,
        noise=NoiseNone(),
        gain=3.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(
            avalanche_gain=5.0,
            excess_noise_factor=1.0,
            multiplication_model=
                ConditionalGammaAvalancheMultiplication()),
    )
    exact_gain_passed = capture!(
        exact, fill(7.0, 4, 4), Xoshiro(6201)) ==
        fill(105.0, 4, 4)

    saturated = Detector(
        integration_time=1.0,
        qe=1.0,
        noise=NoiseNone(),
        full_well=100.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(avalanche_gain=5.0),
    )
    input_referred_saturation_passed = capture!(
        saturated, fill(50.0, 4, 4), Xoshiro(6202)) ==
        fill(100.0, 4, 4)

    noisy_signal = Detector(
        integration_time=1.0,
        qe=1.0,
        noise=NoiseReadout(2.0),
        gain=3.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(
            avalanche_gain=4.0,
            excess_noise_factor=1.0),
    )
    noise_only = Detector(
        integration_time=1.0,
        qe=1.0,
        noise=NoiseReadout(2.0),
        gain=3.0,
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(
            avalanche_gain=1.0,
            excess_noise_factor=1.0),
    )
    signal_output = copy(capture!(
        noisy_signal, fill(10.0, 4, 4), Xoshiro(6203)))
    noise_output = copy(capture!(
        noise_only, zeros(4, 4), Xoshiro(6203)))
    read_noise_ordering_passed = isapprox(
        signal_output, fill(120.0, 4, 4) + noise_output;
        rtol=0.0, atol=256eps(Float64))

    generated_charge_ordering_passed =
        all((:dark, :glow)) do generated_charge_kind
        kwargs = generated_charge_kind == :dark ?
            (; dark_current=8.0) : (;)
        sensor_kwargs = generated_charge_kind == :glow ?
            (; glow_rate=8.0) : (;)
        base = Detector(;
            integration_time=1.0, qe=1.0, noise=NoiseNone(),
            response_model=NullFrameResponse(),
            sensor=HgCdTeAvalancheArraySensor(;
                avalanche_gain=1.0, sensor_kwargs...), kwargs...)
        gained = Detector(;
            integration_time=1.0, qe=1.0, noise=NoiseNone(),
            response_model=NullFrameResponse(),
            sensor=HgCdTeAvalancheArraySensor(;
                avalanche_gain=5.0, sensor_kwargs...), kwargs...)
        base_output = copy(capture!(
            base, zeros(4, 4), Xoshiro(6204)))
        gained_output = copy(capture!(
            gained, zeros(4, 4), Xoshiro(6204)))
        gained_output == 5base_output
    end

    sampling_modes = (
        SingleRead(),
        AveragedNonDestructiveReads(4),
        CorrelatedDoubleSampling(),
        FowlerSampling(4),
        UpTheRampSampling(5),
    )
    sampling_modes_passed = all(sampling_modes) do mode
        detector = Detector(
            integration_time=2.0, qe=1.0, noise=NoiseNone(),
            response_model=NullFrameResponse(),
            sensor=HgCdTeAvalancheArraySensor(
                avalanche_gain=3.0,
                sampling_mode=mode))
        capture!(detector, fill(4.0, 4, 4), Xoshiro(6205)) ==
            fill(24.0, 4, 4)
    end

    windowed = Detector(
        integration_time=1.0, qe=1.0, noise=NoiseNone(),
        response_model=NullFrameResponse(),
        readout_window=FrameWindow(2:3, 2:3),
        sensor=HgCdTeAvalancheArraySensor(
            avalanche_gain=2.0, read_time=0.25,
            sampling_mode=CorrelatedDoubleSampling()))
    windowed_output = capture!(
        windowed, fill(4.0, 4, 4), Xoshiro(6206))
    window_timing_passed =
        windowed_output == fill(8.0, 2, 2) &&
        detector_export_metadata(windowed).sampling_wallclock_time ==
            1.5

    response = GaussianPixelResponse(response_width_px=0.7)
    response_detector = Detector(
        qe=1.0, noise=NoiseNone(), response_model=response,
        sensor=HgCdTeAvalancheArraySensor())
    mtf_before = detector_mtf(
        response_detector, 0.17, -0.09)
    capture!(response_detector, ones(4, 4), Xoshiro(6207))
    configured_mtf_preserved = detector_mtf(
        response_detector, 0.17, -0.09) == mtf_before

    ipc_detector = Detector(
        qe=1.0, noise=NoiseNone(),
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(avalanche_gain=2.0),
        charge_coupling_model=InterpixelCapacitance(
            [0.0 0.01 0.0; 0.01 0.96 0.01; 0.0 0.01 0.0]))
    impulse = zeros(5, 5)
    impulse[3, 3] = 50.0
    ipc_output = capture!(ipc_detector, impulse, Xoshiro(6208))
    configured_ipc_passed =
        ipc_output[3, 3] == 96.0 &&
        ipc_output[2, 3] == 1.0 &&
        sum(ipc_output) == 100.0

    architecture_separated =
        Detectors.detector_sensor_symbol(HgCdTeSensor()) == :hgcdte &&
        Detectors.detector_sensor_symbol(HgCdTeAvalancheArraySensor()) ==
            :hgcdte_linear_avalanche_array &&
        !Detectors.supports_avalanche_gain(HgCdTeSensor()) &&
        Detectors.supports_avalanche_gain(
            HgCdTeAvalancheArraySensor()) &&
        !Detectors.supports_photon_counting(
            HgCdTeAvalancheArraySensor())

    gamma_accelerator_rejected = try
        Detectors.validate_hgcdte_avalanche_backend(
            ConditionalGammaAvalancheMultiplication(),
            CUDABackend())
        false
    catch error
        error isa InvalidConfiguration
    end

    replay_a = Detector(
        qe=1.0, noise=NoiseNone(),
        sensor=HgCdTeAvalancheArraySensor(
            avalanche_gain=5.0, excess_noise_factor=1.5,
            multiplication_model=
                ConditionalGammaAvalancheMultiplication()))
    replay_b = Detector(
        qe=1.0, noise=NoiseNone(),
        sensor=HgCdTeAvalancheArraySensor(
            avalanche_gain=5.0, excess_noise_factor=1.5,
            multiplication_model=
                ConditionalGammaAvalancheMultiplication()))
    deterministic_replay_passed =
        capture!(replay_a, fill(8.0, 8, 8), Xoshiro(6209)) ==
        capture!(replay_b, fill(8.0, 8, 8), Xoshiro(6209))

    allocation_input = fill(64.0, 16, 16)
    approximate_allocation_detector = Detector(
        qe=1.0, noise=NoiseNone(),
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(
            avalanche_gain=5.0, excess_noise_factor=1.5))
    approximate_rng = Xoshiro(6210)
    capture!(approximate_allocation_detector,
        allocation_input, approximate_rng)
    approximate_steady_alloc_bytes = @allocated capture!(
        approximate_allocation_detector,
        allocation_input, approximate_rng)

    gamma_allocation_detector = Detector(
        qe=1.0, noise=NoiseNone(),
        response_model=NullFrameResponse(),
        sensor=HgCdTeAvalancheArraySensor(
            avalanche_gain=5.0, excess_noise_factor=1.5,
            multiplication_model=
                ConditionalGammaAvalancheMultiplication()))
    gamma_rng = Xoshiro(6211)
    capture!(gamma_allocation_detector,
        allocation_input, gamma_rng)
    gamma_steady_alloc_bytes = @allocated capture!(
        gamma_allocation_detector, allocation_input, gamma_rng)

    warm_scheduled = hgcdte_avalanche_scheduled_fixture(6212)
    scheduled_retains_prior_multiplication =
        run_hgcdte_avalanche_scheduled_ramp!(
            warm_scheduled.prepared, warm_scheduled.state,
            warm_scheduled.values, warm_scheduled.rng)
    scheduled = hgcdte_avalanche_scheduled_fixture(6213)
    scheduled_steady_alloc_bytes = @allocated(
        run_hgcdte_avalanche_scheduled_ramp!(
            scheduled.prepared, scheduled.state,
            scheduled.values, scheduled.rng))

    return Dict{String,Any}(
        "architecture_separated" => architecture_separated,
        "exact_gain_passed" => exact_gain_passed,
        "input_referred_saturation_passed" =>
            input_referred_saturation_passed,
        "generated_charge_ordering_passed" =>
            generated_charge_ordering_passed,
        "read_noise_and_conversion_gain_ordering_passed" =>
            read_noise_ordering_passed,
        "single_ndr_cds_fowler_ramp_passed" =>
            sampling_modes_passed,
        "window_preserves_full_frame_timing" =>
            window_timing_passed,
        "configured_mtf_preserved" => configured_mtf_preserved,
        "configured_ipc_passed" => configured_ipc_passed,
        "scheduled_retains_prior_multiplication" =>
            scheduled_retains_prior_multiplication,
        "gamma_accelerator_rejected" =>
            gamma_accelerator_rejected,
        "deterministic_replay_passed" =>
            deterministic_replay_passed,
        "approximate_steady_alloc_bytes" =>
            approximate_steady_alloc_bytes,
        "gamma_steady_alloc_bytes" =>
            gamma_steady_alloc_bytes,
        "scheduled_steady_alloc_bytes" =>
            scheduled_steady_alloc_bytes,
        "allocation_gate_passed" =>
            approximate_steady_alloc_bytes == 0 &&
            gamma_steady_alloc_bytes == 0 &&
            scheduled_steady_alloc_bytes == 0,
    )
end

hgcdte_avalanche_git_revision() = try
    readchomp(`git rev-parse HEAD`)
catch
    "unknown"
end

hgcdte_avalanche_git_dirty() =
    !success(`git diff --quiet HEAD`)

function generate_hgcdte_avalanche_qualification_artifact()
    moment_cases = hgcdte_avalanche_moment_cases()
    deterministic = hgcdte_avalanche_deterministic_contract()
    moment_gates_passed = all(case ->
        case["mean_passed"] &&
        case["variance_passed"] &&
        case["all_nonnegative"] &&
        case["moderate_charge_approximation_regime"],
        moment_cases)
    deterministic_gates = (
        "architecture_separated",
        "exact_gain_passed",
        "input_referred_saturation_passed",
        "generated_charge_ordering_passed",
        "read_noise_and_conversion_gain_ordering_passed",
        "single_ndr_cds_fowler_ramp_passed",
        "window_preserves_full_frame_timing",
        "configured_mtf_preserved",
        "configured_ipc_passed",
        "scheduled_retains_prior_multiplication",
        "gamma_accelerator_rejected",
        "deterministic_replay_passed",
        "allocation_gate_passed",
    )
    deterministic_gates_passed =
        all(deterministic[key] for key in deterministic_gates)
    cpu = first(Sys.cpu_info())
    artifact = Dict{String,Any}(
        "schema_version" => 1,
        "artifact_id" => HGCDTE_AVALANCHE_ARTIFACT_ID,
        "family" => "hgcdte_linear_avalanche_photodiode_array",
        "all_gates_passed" =>
            moment_gates_passed && deterministic_gates_passed,
        "qualification" => Dict{String,Any}(
            "frame_count_per_case" =>
                HGCDTE_AVALANCHE_FRAME_COUNT,
            "frame_size" => [
                HGCDTE_AVALANCHE_FRAME_SIZE,
                HGCDTE_AVALANCHE_FRAME_SIZE,
            ],
            "samples_per_case" =>
                HGCDTE_AVALANCHE_FRAME_COUNT *
                HGCDTE_AVALANCHE_FRAME_SIZE^2,
            "sigma_limit" => HGCDTE_AVALANCHE_SIGMA_LIMIT,
            "multiplication_moment_cases" => moment_cases,
            "deterministic" => deterministic,
        ),
        "model" => Dict{String,Any}(
            "excess_noise_factor_definition" =>
                "F = E[M^2] / E[M]^2",
            "conditional_variance" =>
                "avalanche_gain^2 * (F - 1) * input_charge",
            "conditional_gamma" =>
                "shape = input_charge / (F - 1); scale = F - 1 before mean avalanche gain",
            "clipped_gaussian_approximation" =>
                "max(q + sqrt((F - 1)q)Z, 0), qualified only when q / (F - 1) >= 25",
            "scientific_reference" =>
                "https://arxiv.org/abs/1903.08444",
        ),
        "scope" => Dict{String,Any}(
            "included" => [
                "HgCdTe linear-avalanche photodiode area arrays",
                "CPU conditional-Gamma multiplication",
                "backend-portable clipped-Gaussian moderate-charge approximation",
                "input-referred saturation and exact avalanche/read/conversion-gain ordering",
                "dark current and readout glow before avalanche multiplication",
                "single, averaged nondestructive, correlated-double, Fowler, and up-the-ramp sampling",
                "full-frame and windowed output with explicit response/MTF and IPC",
                "scheduled evolving-charge ramps that retain prior multiplication realizations",
            ],
            "excluded" => [
                "Geiger-mode counting or photon timestamp output",
                "conditional-Gamma accelerator execution",
                "low-charge distribution fidelity for the clipped-Gaussian approximation",
                "bias-voltage or gain-temperature calibration",
                "per-pixel avalanche-gain maps",
                "bandwidth-dependent readout electronics",
                "signal-dependent IPC",
                "camera or vendor profiles",
            ],
        ),
        "environment" => Dict{String,Any}(
            "timestamp_utc" => string(Dates.now(Dates.UTC)),
            "source_revision" =>
                hgcdte_avalanche_git_revision(),
            "source_dirty" => hgcdte_avalanche_git_dirty(),
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
    mkpath(HGCDTE_AVALANCHE_ARTIFACT_DIR)
    open(HGCDTE_AVALANCHE_ARTIFACT_PATH, "w") do io
        TOML.print(io, artifact; sorted=true)
    end
    artifact["all_gates_passed"] || error(
        "HgCdTe linear-avalanche qualification artifact contains a failed gate")
    println("wrote ", HGCDTE_AVALANCHE_ARTIFACT_PATH)
    return artifact
end

generate_hgcdte_avalanche_qualification_artifact()

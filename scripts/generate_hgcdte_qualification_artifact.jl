using AdaptiveOpticsSim
using AdaptiveOpticsSim.Detectors
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.Plant
using Dates
using Random
using Statistics
using TOML

const HGCDTE_ARTIFACT_DIR =
    joinpath(@__DIR__, "..", "benchmarks", "results", "detectors")
const HGCDTE_ARTIFACT_PATH = joinpath(HGCDTE_ARTIFACT_DIR,
    "2026-07-30-conventional-hgcdte-qualification.toml")
const HGCDTE_ARTIFACT_ID = "DET-HGCDTE-QUAL-2026-07-30"
const HGCDTE_FRAME_SIZE = 128
const HGCDTE_SAMPLES_PER_CASE = HGCDTE_FRAME_SIZE^2
const HGCDTE_SIGMA_LIMIT = 6.0

function hgcdte_samples(detector, seed)
    input = zeros(HGCDTE_FRAME_SIZE, HGCDTE_FRAME_SIZE)
    return vec(copy(capture!(detector, input, Xoshiro(seed))))
end

function hgcdte_moment_case(id, detector, seed;
    expected_mean, expected_variance, expected_fourth_central)
    samples = hgcdte_samples(detector, seed)
    sample_count = length(samples)
    mean_limit = HGCDTE_SIGMA_LIMIT *
        sqrt(expected_variance / sample_count)
    variance_se = sqrt((expected_fourth_central -
        ((sample_count - 3) / (sample_count - 1)) *
        expected_variance^2) / sample_count)
    variance_limit = HGCDTE_SIGMA_LIMIT * variance_se
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

function hgcdte_gaussian_case(id, mode, seed, expected_sigma)
    variance = expected_sigma^2
    detector = Detector(
        exposure_duration=1.0,
        qe=1.0,
        noise=NoiseReadout(4.0),
        response_model=NullFrameResponse(),
        sensor=HgCdTeSensor(sampling_mode=mode),
    )
    return hgcdte_moment_case(id, detector, seed;
        expected_mean=0.0,
        expected_variance=variance,
        expected_fourth_central=3variance^2)
end

function hgcdte_poisson_case(id, detector, seed, expected_mean)
    return hgcdte_moment_case(id, detector, seed;
        expected_mean=expected_mean,
        expected_variance=expected_mean,
        expected_fourth_central=expected_mean + 3expected_mean^2)
end

function hgcdte_moment_cases()
    sigma = 4.0
    return Dict{String,Any}[
        hgcdte_gaussian_case(
            "single_read_noise", SingleRead(), 5101, sigma),
        hgcdte_gaussian_case(
            "averaged_ndr4_noise", AveragedNonDestructiveReads(4),
            5102, sigma / 2),
        hgcdte_gaussian_case(
            "correlated_double_sampling_noise",
            CorrelatedDoubleSampling(), 5103, sigma * sqrt(2)),
        hgcdte_gaussian_case(
            "fowler8_noise", FowlerSampling(8), 5104, sigma / 2),
        hgcdte_gaussian_case(
            "up_the_ramp16_noise", UpTheRampSampling(16), 5105,
            sigma * sqrt(12 * 15 / (16 * 17))),
        hgcdte_poisson_case(
            "dark_current",
            Detector(exposure_duration=2.0, qe=1.0,
                dark_current=5.0, noise=NoiseNone(),
                response_model=NullFrameResponse(),
                sensor=HgCdTeSensor()),
            5106, 10.0),
        hgcdte_poisson_case(
            "readout_glow",
            Detector(exposure_duration=1.0, qe=1.0,
                noise=NoiseNone(), response_model=NullFrameResponse(),
                sensor=HgCdTeSensor(glow_rate=7.0)),
            5107, 7.0),
    ]
end

function hgcdte_rate_map(values::AbstractMatrix{T}) where {T<:AbstractFloat}
    metadata = OpticalPlaneMetadata(DetectorPlane(), values;
        coordinate_domain=AngularCoordinates(),
        sampling=(one(T), one(T)),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition(),
        spectral=MonochromaticChannel(T(0.55e-6)))
    return IntensityMap(metadata, values)
end

function hgcdte_irregular_scheduled_ramp_contract()
    rate = 1.0e9
    values = fill(rate, 2, 2)
    detector = Detector(
        exposure_duration=1.0e-8,
        qe=1.0,
        noise=NoiseNone(),
        response_model=NullFrameResponse(),
        sensor=HgCdTeSensor(
            read_duration=0.0,
            sampling_mode=UpTheRampSampling(4)),
    )
    definition = GlobalShutterAcquisitionDefinition(PlantDuration(10))
    prepared = Plant.prepare_global_shutter_acquisition(
        detector, hgcdte_rate_map(values), definition)
    state = Plant.GlobalShutterAcquisitionState(prepared)
    rng = Xoshiro(5201)
    start = PlantTimestamp(0)
    Plant.begin_exposure!(prepared, state, start)
    Plant.take_nondestructive_read!(prepared, state, start, rng)
    previous = start
    for read_index in 2:3
        timestamp = start + Plant.nondestructive_read_offset(
            prepared, read_index)
        Plant.accumulate_exposure_interval!(
            prepared, state, previous, timestamp, rng)
        Plant.take_nondestructive_read!(
            prepared, state, timestamp, rng)
        previous = timestamp
    end
    close = start + definition.exposure_duration
    Plant.accumulate_exposure_interval!(
        prepared, state, previous, close, rng)
    Plant.close_exposure!(prepared, state, close)
    Plant.take_nondestructive_read!(prepared, state, close, rng)
    output = Plant.complete_readout!(prepared, state, close, rng)
    return detector_ramp_read_offsets_s(detector) ==
            [0.0, 3.0e-9, 6.0e-9, 1.0e-8] &&
        all(isapprox.(detector_ramp_slope(detector), rate;
            rtol=8eps(Float64), atol=0.0)) &&
        maximum(abs, detector_ramp_intercept(detector)) <=
            8eps(Float64) &&
        all(isapprox.(output, 10.0;
            rtol=8eps(Float64), atol=0.0)) &&
        Detectors.detector_ramp_acquisition(detector) ==
            :scheduled_evolving_charge
end

function hgcdte_deterministic_contract()
    input = fill(3.0, 4, 4)
    modes = (
        SingleRead(),
        AveragedNonDestructiveReads(4),
        CorrelatedDoubleSampling(),
        FowlerSampling(4),
    )
    sampling_modes_passed = all(modes) do mode
        detector = Detector(
            exposure_duration=2.0,
            qe=1.0,
            noise=NoiseNone(),
            response_model=NullFrameResponse(),
            sensor=HgCdTeSensor(sampling_mode=mode),
        )
        capture!(detector, input, Xoshiro(5202)) == fill(6.0, 4, 4)
    end

    direct_ramp = Detector(
        exposure_duration=2.0,
        qe=1.0,
        noise=NoiseNone(),
        response_model=NullFrameResponse(),
        sensor=HgCdTeSensor(
            read_duration=0.1,
            sampling_mode=UpTheRampSampling(5)),
    )
    direct_output = capture!(direct_ramp, input, Xoshiro(5203))
    direct_ramp_passed =
        direct_output == fill(6.0, 4, 4) &&
        detector_ramp_slope(direct_ramp) == fill(3.0, 4, 4) &&
        detector_ramp_intercept(direct_ramp) == zeros(4, 4) &&
        detector_ramp_read_offsets_s(direct_ramp) ==
            [0.0, 0.5, 1.0, 1.5, 2.0] &&
        Detectors.detector_ramp_acquisition(direct_ramp) ==
            :synthesized_final_charge

    windowed = Detector(
        exposure_duration=1.0,
        qe=1.0,
        noise=NoiseNone(),
        response_model=NullFrameResponse(),
        sensor=HgCdTeSensor(
            read_duration=0.25,
            sampling_mode=CorrelatedDoubleSampling()),
        readout_window=Detectors.FrameWindow(2:3, 2:3),
    )
    windowed_output = capture!(windowed, fill(4.0, 4, 4),
        Xoshiro(5204))
    window_preserves_full_frame_timing =
        windowed_output == fill(4.0, 2, 2) &&
        detector_export_metadata(windowed).sampling_read_duration == 0.25 &&
        detector_export_metadata(windowed).sampling_acquisition_duration == 1.5 &&
        detector_ramp_read_offsets_s(windowed) === nothing

    response = GaussianPixelResponse(response_width_px=0.7)
    response_detector = Detector(
        qe=1.0, noise=NoiseNone(), response_model=response,
        sensor=HgCdTeSensor())
    mtf_before = detector_mtf(response_detector, 0.17, -0.09)
    capture!(response_detector, ones(4, 4), Xoshiro(5205))
    configured_mtf_preserved =
        detector_mtf(response_detector, 0.17, -0.09) == mtf_before

    ipc_detector = Detector(
        qe=1.0, noise=NoiseNone(), response_model=NullFrameResponse(),
        sensor=HgCdTeSensor(),
        charge_coupling_model=InterpixelCapacitance(
            [0.0 0.01 0.0; 0.01 0.96 0.01; 0.0 0.01 0.0]))
    impulse = zeros(5, 5)
    impulse[3, 3] = 100.0
    ipc_output = capture!(ipc_detector, impulse, Xoshiro(5206))
    configured_ipc_passed =
        ipc_output[3, 3] == 96.0 &&
        ipc_output[2, 3] == 1.0 &&
        sum(ipc_output) == 100.0

    correction_detector = Detector(
        qe=1.0, noise=NoiseNone(), response_model=NullFrameResponse(),
        sensor=HgCdTeSensor(),
        correction_model=ReferencePixelCommonModeCorrection(1, 1))
    correction_passed = maximum(abs, capture!(correction_detector,
        fill(5.0, 4, 4), Xoshiro(5207))) <= 8eps(Float64)

    saturation_detector = Detector(
        qe=1.0, noise=NoiseNone(), response_model=NullFrameResponse(),
        full_well=100.0, sensor=HgCdTeSensor())
    saturation_passed = capture!(saturation_detector,
        fill(150.0, 2, 2), Xoshiro(5208)) == fill(100.0, 2, 2)

    nonlinear_detector = Detector(
        qe=1.0, noise=NoiseNone(), response_model=NullFrameResponse(),
        nonlinearity_model=SaturatingFrameNonlinearity(0.1),
        sensor=HgCdTeSensor())
    nonlinearity_passed = capture!(nonlinear_detector,
        fill(10.0, 2, 2), Xoshiro(5209)) == fill(5.0, 2, 2)

    persistence_detector = Detector(
        qe=1.0, noise=NoiseNone(), response_model=NullFrameResponse(),
        sensor=HgCdTeSensor(
            persistence_model=ExponentialPersistence(0.5, 0.0)))
    capture!(persistence_detector, fill(4.0, 2, 2), Xoshiro(5210))
    persistence_passed = capture!(persistence_detector,
        zeros(2, 2), Xoshiro(5211)) == fill(2.0, 2, 2)

    shared_readout = Detectors.HgCdTeReadout(
        read_duration=0.1, sampling_mode=FowlerSampling(2))
    architecture_separated =
        detector_export_metadata(Detector(
            sensor=HgCdTeSensor(shared_readout))).sensor == :hgcdte &&
        !Detectors.supports_avalanche_gain(HgCdTeSensor(shared_readout)) &&
        Detectors.supports_avalanche_gain(HgCdTeAvalancheArraySensor(
            shared_readout; avalanche_gain=3.0))

    replay_a = Detector(
        qe=1.0, noise=NoiseReadout(1.0),
        sensor=HgCdTeSensor(sampling_mode=FowlerSampling(2)))
    replay_b = Detector(
        qe=1.0, noise=NoiseReadout(1.0),
        sensor=HgCdTeSensor(sampling_mode=FowlerSampling(2)))
    deterministic_replay_passed =
        capture!(replay_a, zeros(8, 8), Xoshiro(5212)) ==
        capture!(replay_b, zeros(8, 8), Xoshiro(5212))

    allocation_rng = Xoshiro(5213)
    capture!(direct_ramp, input, allocation_rng)
    steady_alloc_bytes =
        @allocated capture!(direct_ramp, input, allocation_rng)

    return Dict{String,Any}(
        "architecture_separated" => architecture_separated,
        "single_ndr_cds_fowler_passed" => sampling_modes_passed,
        "direct_synthesized_ramp_passed" => direct_ramp_passed,
        "irregular_scheduled_ramp_passed" =>
            hgcdte_irregular_scheduled_ramp_contract(),
        "window_preserves_full_frame_timing" =>
            window_preserves_full_frame_timing,
        "configured_mtf_preserved" => configured_mtf_preserved,
        "configured_ipc_passed" => configured_ipc_passed,
        "reference_correction_passed" => correction_passed,
        "saturation_passed" => saturation_passed,
        "nonlinearity_passed" => nonlinearity_passed,
        "persistence_passed" => persistence_passed,
        "deterministic_replay_passed" => deterministic_replay_passed,
        "steady_alloc_bytes" => steady_alloc_bytes,
        "allocation_gate_passed" => steady_alloc_bytes == 0,
    )
end

hgcdte_git_revision() = try
    readchomp(`git rev-parse HEAD`)
catch
    "unknown"
end

hgcdte_git_dirty() = !success(`git diff --quiet HEAD`)

function generate_hgcdte_qualification_artifact()
    moment_cases = hgcdte_moment_cases()
    deterministic = hgcdte_deterministic_contract()
    moment_gates_passed = all(case ->
        case["mean_passed"] && case["variance_passed"], moment_cases)
    deterministic_gates = (
        "architecture_separated",
        "single_ndr_cds_fowler_passed",
        "direct_synthesized_ramp_passed",
        "irregular_scheduled_ramp_passed",
        "window_preserves_full_frame_timing",
        "configured_mtf_preserved",
        "configured_ipc_passed",
        "reference_correction_passed",
        "saturation_passed",
        "nonlinearity_passed",
        "persistence_passed",
        "deterministic_replay_passed",
        "allocation_gate_passed",
    )
    deterministic_gates_passed =
        all(deterministic[key] for key in deterministic_gates)
    cpu = first(Sys.cpu_info())
    artifact = Dict{String,Any}(
        "schema_version" => 1,
        "artifact_id" => HGCDTE_ARTIFACT_ID,
        "family" => "conventional_hgcdte",
        "all_gates_passed" =>
            moment_gates_passed && deterministic_gates_passed,
        "qualification" => Dict{String,Any}(
            "frame_size" => [HGCDTE_FRAME_SIZE, HGCDTE_FRAME_SIZE],
            "samples_per_case" => HGCDTE_SAMPLES_PER_CASE,
            "sigma_limit" => HGCDTE_SIGMA_LIMIT,
            "moment_cases" => moment_cases,
            "deterministic" => deterministic,
        ),
        "scope" => Dict{String,Any}(
            "included" => [
                "conventional HgCdTe area arrays without avalanche multiplication",
                "single, averaged nondestructive, correlated-double, Fowler, and up-the-ramp sampling",
                "direct synthesized-final-charge and Plant scheduled evolving-charge ramps",
                "dark current, readout glow, independent Gaussian read noise, nonlinearity, persistence, saturation, and reference correction",
                "explicit presampling response, MTF, and post-collection IPC",
            ],
            "excluded" => [
                "cosmic-ray jump detection or ramp segmentation",
                "saturation-aware ramp fitting",
                "correlated 1/f read noise",
                "reciprocity failure",
                "calibrated persistence or detector maps",
                "signal-dependent IPC",
                "camera or vendor profiles",
                "linear-avalanche multiplication qualification",
            ],
        ),
        "environment" => Dict{String,Any}(
            "timestamp_utc" => string(Dates.now(Dates.UTC)),
            "source_revision" => hgcdte_git_revision(),
            "source_dirty" => hgcdte_git_dirty(),
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
    mkpath(HGCDTE_ARTIFACT_DIR)
    open(HGCDTE_ARTIFACT_PATH, "w") do io
        TOML.print(io, artifact; sorted=true)
    end
    artifact["all_gates_passed"] || error(
        "HgCdTe qualification artifact contains a failed gate")
    println("wrote ", HGCDTE_ARTIFACT_PATH)
    return artifact
end

generate_hgcdte_qualification_artifact()

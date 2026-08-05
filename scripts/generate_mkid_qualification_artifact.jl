using AdaptiveOpticsSim
using AdaptiveOpticsSim.Detectors
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.WavefrontSensors
using Dates
using Random
using Statistics
using TOML

const MKID_ARTIFACT_PATH = joinpath(
    @__DIR__, "..", "benchmarks", "results", "detectors",
    "2026-08-01-mkid-qualification.toml")
const MKID_SAMPLE_COUNT = 16_384
const MKID_SIGMA_LIMIT = 6.0

function mkid_poisson_case(id, detector, input, expected_mean, seed)
    samples = vec(copy(capture!(detector, input, Xoshiro(seed))))
    sample_count = length(samples)
    mean_limit = MKID_SIGMA_LIMIT * sqrt(expected_mean / sample_count)
    variance_limit = MKID_SIGMA_LIMIT * sqrt(
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

function mkid_poisson_cases()
    dimensions = (128, 128)
    return [
        mkid_poisson_case("photon_only",
            MKIDArrayDetector(noise=NoisePhoton(),
                sensor=MKIDArraySensor(qe=1.0)),
            fill(20.0, dimensions), 20.0, 9260),
        mkid_poisson_case("dark_only",
            MKIDArrayDetector(exposure_duration=2.0, noise=NoisePhoton(),
                sensor=MKIDArraySensor(qe=0.0, dark_count_rate=6.0)),
            zeros(dimensions), 12.0, 9261),
        mkid_poisson_case("photon_and_dark",
            MKIDArrayDetector(noise=NoisePhoton(),
                sensor=MKIDArraySensor(qe=0.5, dark_count_rate=3.0)),
            fill(10.0, dimensions), 8.0, 9262),
        mkid_poisson_case("dead_time_adjusted_surrogate",
            MKIDArrayDetector(noise=NoisePhoton(),
                sensor=MKIDArraySensor(qe=1.0,
                    dead_time_model=NonParalyzableDeadTime(0.05))),
            fill(100.0, dimensions), 100 / 6, 9263),
    ]
end

function mkid_dead_time_curves()
    dead_time = 0.1
    dimensionless_rates = (0.0, 1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0)
    curves = Dict{String,Any}[]
    for (id, model, response) in (
        ("no_dead_time", NoDeadTime(), x -> x / dead_time),
        ("nonparalyzable", NonParalyzableDeadTime(dead_time),
            x -> (x / dead_time) / (1 + x)),
        ("paralyzable", ParalyzableDeadTime(dead_time),
            x -> (x / dead_time) * exp(-x)),
    )
        detector = MKIDArrayDetector(noise=NoiseNone(),
            sensor=MKIDArraySensor(qe=1.0, dead_time_model=model))
        points = Dict{String,Any}[]
        for (index, x) in enumerate(dimensionless_rates)
            observed = only(capture!(detector,
                fill(x / dead_time, 1, 1), Xoshiro(9270 + index)))
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

function mkid_rate_map(values, source)
    T = eltype(values)
    metadata = OpticalPlaneMetadata(DetectorPlane(), values;
        coordinate_domain=AngularCoordinates(),
        sampling=(one(T), one(T)),
        spectral=MonochromaticChannel(T(wavelength(source))),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    return IntensityMap(metadata, values)
end

function mkid_deterministic_contract()
    radiometry = MKIDArrayDetector(exposure_duration=2.0,
        noise=NoiseNone(), gate_model=DutyCycleGate(0.25),
        sensor=MKIDArraySensor(qe=0.5, fill_factor=0.8))
    radiometry_output = copy(capture!(radiometry,
        fill(10.0, 2, 8), Xoshiro(9280)))

    gated_dark = MKIDArrayDetector(exposure_duration=2.0,
        noise=NoiseNone(), gate_model=DutyCycleGate(0.25),
        sensor=MKIDArraySensor(qe=0.0, dark_count_rate=4.0))
    gated_dark_output = copy(capture!(gated_dark,
        zeros(1, 1), Xoshiro(9281)))

    passband = (0.8e-6, 1.4e-6)
    characteristics = MKIDArrayCharacteristics(
        energy_resolving_power=12.0,
        photon_arrival_time_resolution_s=2e-6,
        wavelength_passband_m=passband)
    filtered = MKIDArrayDetector(noise=NoiseNone(),
        sensor=MKIDArraySensor(qe=1.0,
            characteristics=characteristics))
    input = fill(10.0, 2, 2)
    inside = Source(band=:custom, wavelength=1.0e-6)
    lower = Source(band=:custom, wavelength=first(passband))
    upper = Source(band=:custom, wavelength=last(passband))
    outside = Source(band=:custom, wavelength=0.6e-6)
    inside_output = copy(capture!(filtered, input, inside, Xoshiro(9282)))
    lower_output = copy(capture!(filtered, input, lower, Xoshiro(9283)))
    upper_output = copy(capture!(filtered, input, upper, Xoshiro(9284)))
    outside_output = copy(capture!(filtered, input, outside, Xoshiro(9285)))
    matrix_output = copy(capture!(filtered, input, Xoshiro(9286)))
    spectral = with_spectrum(inside,
        SpectralBundle([0.6e-6, 0.8e-6, 1.0e-6, 1.4e-6, 1.6e-6],
            [0.1, 0.2, 0.3, 0.15, 0.25]))
    spectral_output = copy(capture!(filtered, input, spectral,
        Xoshiro(9287)))

    integer = MKIDArrayDetector(noise=NoiseNone(),
        sensor=MKIDArraySensor(qe=1.0), output_type=UInt16)
    rounded_output = copy(capture!(integer,
        fill(2.6, 1, 2), Xoshiro(9288)))
    saturated_output = copy(capture!(integer,
        fill(1e9, 1, 2), Xoshiro(9289)))

    failure_detector = MKIDArrayDetector(noise=NoiseNone(),
        sensor=MKIDArraySensor(qe=1.0), output_type=UInt16)
    capture!(failure_detector, fill(2.0, 2, 3), Xoshiro(9290))
    arrays = (failure_detector.state.counts,
        failure_detector.state.noise_buffer,
        failure_detector.state.host_buffer,
        failure_detector.state.output_buffer,
        failure_detector.state.output_buffer_host)
    snapshots = map(copy, arrays)
    invalid_input_rejected = all((fill(-1.0, 3, 2), fill(Inf, 3, 2),
        fill(NaN, 3, 2))) do invalid_input
        try
            capture!(failure_detector, invalid_input, Xoshiro(9291))
            false
        catch error
            error isa InvalidConfiguration
        end
    end
    invalid_input_preserved = all(begin
        current = (failure_detector.state.counts,
            failure_detector.state.noise_buffer,
            failure_detector.state.host_buffer,
            failure_detector.state.output_buffer,
            failure_detector.state.output_buffer_host)
        current[index] === arrays[index] &&
            current[index] == snapshots[index]
    end for index in eachindex(arrays))

    function replay_detector()
        return MKIDArrayDetector(exposure_duration=0.5,
            noise=NoisePhoton(), gate_model=DutyCycleGate(0.75),
            sensor=MKIDArraySensor(qe=0.7, dark_count_rate=0.25,
                dead_time_model=NonParalyzableDeadTime(1e-3)))
    end
    replay_input = fill(40.0, 16, 16)
    replay_passed =
        capture!(replay_detector(), replay_input, Xoshiro(9292)) ==
        capture!(replay_detector(), replay_input, Xoshiro(9292))

    allocation_detector = MKIDArrayDetector(noise=NoiseNone(),
        sensor=MKIDArraySensor(qe=0.5))
    allocation_input = fill(2.0, 8, 8)
    allocation_rng = Xoshiro(9293)
    capture!(allocation_detector, allocation_input, allocation_rng)
    steady_alloc_bytes = @allocated capture!(allocation_detector,
        allocation_input, allocation_rng)

    rate = mkid_rate_map(copy(input), inside)
    observation = WFSObservation(zeros(2, 2);
        units=:photon_count, layout=:counting_channels)
    prepared_detector = MKIDArrayDetector(noise=NoiseNone(),
        sensor=MKIDArraySensor(qe=1.0,
            characteristics=characteristics))
    plan = prepare_wfs_acquisition(prepared_detector, rate, observation;
        source=inside)
    acquire_wfs_observation!(observation, rate, plan, Xoshiro(9294))

    unfiltered_plan = prepare_wfs_acquisition(
        MKIDArrayDetector(noise=NoiseNone(),
            sensor=MKIDArraySensor(qe=1.0)),
        rate,
        WFSObservation(zeros(2, 2);
            units=:photon_count, layout=:counting_channels))

    metadata = detector_export_metadata(filtered)
    default_metadata = detector_export_metadata(
        MKIDArrayDetector(noise=NoiseNone()))
    return Dict{String,Any}(
        "exact_radiometry_and_gate_passed" =>
            radiometry_output == fill(2.0, 2, 8),
        "gated_dark_live_time_passed" =>
            gated_dark_output == fill(2.0, 1, 1),
        "inclusive_passband_passed" =>
            inside_output == input && lower_output == input &&
            upper_output == input && all(iszero, outside_output),
        "weighted_spectral_bundle_passed" =>
            spectral_output ≈ fill(6.5, 2, 2),
        "matrix_prefilter_contract_passed" => matrix_output == input,
        "integer_rounding_and_saturation_passed" =>
            rounded_output == fill(UInt16(3), 1, 2) &&
            saturated_output == fill(typemax(UInt16), 1, 2),
        "characteristics_separated_from_observable" =>
            metadata.observable == :accumulated_count_image &&
            metadata.characteristics.energy_resolving_power == 12.0 &&
            metadata.characteristics.photon_arrival_time_resolution_s == 2e-6 &&
            metadata.characteristics.wavelength_passband_m == passband &&
            !metadata.provides_energy_estimates &&
            !metadata.provides_photon_arrival_timestamps,
        "optional_characteristics_default_to_nothing" =>
            default_metadata.characteristics.energy_resolving_power === nothing &&
            default_metadata.characteristics.photon_arrival_time_resolution_s === nothing &&
            default_metadata.characteristics.wavelength_passband_m === nothing,
        "spad_mean_response_absent" =>
            metadata.counting.mean_response_model == :none &&
            metadata.counting.mean_afterpulses_per_detection === nothing &&
            metadata.counting.redistribution_fraction === nothing,
        "invalid_input_rejected" => invalid_input_rejected,
        "invalid_input_preserved_storage" => invalid_input_preserved,
        "prepared_source_throughput_bound" =>
            plan.source_throughput == 1.0 && observation.storage == input,
        "source_free_preparation_without_passband" =>
            unfiltered_plan.source === nothing &&
            unfiltered_plan.source_throughput == 1.0,
        "detector_mtf_not_applicable" =>
            !applicable(detector_mtf, filtered, 0.1, 0.1),
        "deterministic_replay_passed" => replay_passed,
        "steady_alloc_bytes" => steady_alloc_bytes,
        "allocation_gate_passed" => steady_alloc_bytes == 0,
    )
end

mkid_git_revision() = try
    readchomp(`git rev-parse HEAD`)
catch
    "unknown"
end

mkid_git_dirty() = !success(`git diff --quiet HEAD`)

function generate_mkid_qualification_artifact()
    curves = mkid_dead_time_curves()
    poisson_cases = mkid_poisson_cases()
    deterministic = mkid_deterministic_contract()
    deterministic_gates = (
        "exact_radiometry_and_gate_passed",
        "gated_dark_live_time_passed",
        "inclusive_passband_passed",
        "weighted_spectral_bundle_passed",
        "matrix_prefilter_contract_passed",
        "integer_rounding_and_saturation_passed",
        "characteristics_separated_from_observable",
        "optional_characteristics_default_to_nothing",
        "spad_mean_response_absent",
        "invalid_input_rejected",
        "invalid_input_preserved_storage",
        "prepared_source_throughput_bound",
        "source_free_preparation_without_passband",
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
        "artifact_id" => "DET-MKID-QUAL-2026-08-01",
        "family" => "mkid_accumulated_count_array",
        "all_gates_passed" => all_gates_passed,
        "environment" => Dict(
            "timestamp_utc" => string(now(UTC)),
            "julia_version" => string(VERSION),
            "julia_threads" => Threads.nthreads(),
            "architecture" => string(Sys.ARCH),
            "kernel" => string(Sys.KERNEL),
            "array_backend" => "CPU",
            "source_revision" => mkid_git_revision(),
            "source_dirty" => mkid_git_dirty(),
        ),
        "model" => Dict(
            "input" => "cell-integrated, spectrally prefiltered photon-arrival rate in photons per second per detector pixel",
            "observable" => "accumulated expected-count or sampled-count image per exposure",
            "pipeline_order" => "QE, fill factor, live time, and source throughput; dark expectation; dead-time mean law; optional Poisson surrogate; output conversion",
            "nonparalyzable_mean_law" =>
                "mu / (1 + mu * tau / live_time)",
            "paralyzable_mean_law" =>
                "mu * exp(-mu * tau / live_time)",
            "characteristics" =>
                "optional energy resolving power, photon-arrival-time resolution, and inclusive wavelength passband",
            "statistical_scope" =>
                "Poisson draw from the adjusted accumulated-count mean; not an event-resolved MKID distribution",
            "scientific_references" => [
                "https://doi.org/10.1038/nature02037",
                "https://arxiv.org/abs/1007.0752",
                "https://doi.org/10.1086/674013",
            ],
        ),
        "qualification" => Dict(
            "samples_per_moment_case" => MKID_SAMPLE_COUNT,
            "sigma_limit" => MKID_SIGMA_LIMIT,
            "dead_time_curves" => curves,
            "poisson_surrogate_cases" => poisson_cases,
            "deterministic" => deterministic,
        ),
        "scope" => Dict(
            "included" => [
                "quantum efficiency, scalar fill factor, live time, and dark counts",
                "inclusive scalar passband and weighted spectral-source throughput",
                "nonparalyzable and paralyzable expected-count dead-time laws",
                "Poisson surrogate moments after dead-time mean adjustment",
                "optional physical-characteristic metadata separated from observables",
                "prepared WFS source-throughput binding",
                "deterministic replay and warmed CPU allocation",
            ],
            "excluded" => [
                "per-photon energy estimates and photon-arrival timestamps",
                "measured energy distributions and exact dead-time count distributions",
                "pulse pile-up, resonator dynamics, and quasiparticle physics",
                "readout-tone multiplexing and calibrated spectral-efficiency curves",
                "detector spatial response, pixel aperture, and detector MTF",
                "calibrated detector maps, named modules, and camera profiles",
            ],
        ),
    )
    all_gates_passed ||
        error("MKID qualification artifact contains a failed gate")
    open(MKID_ARTIFACT_PATH, "w") do io
        TOML.print(io, artifact; sorted=true)
    end
    return MKID_ARTIFACT_PATH
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    println(generate_mkid_qualification_artifact())
end

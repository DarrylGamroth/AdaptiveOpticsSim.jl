using AdaptiveOpticsSim
using Base64
using Dates
using HdrHistogram
using LinearAlgebra
using Random
using Statistics
using TOML

const GATE0_CONTRACT_PATH = get(ENV, "AOS_GATE0_CONTRACT",
    joinpath(@__DIR__, "contracts", "pre_hil_gate0.toml"))
const GATE0_BASELINE_PATH = get(ENV, "AOS_GATE0_BASELINE", "")
const GATE0_OUTPUT_PATH = get(ENV, "AOS_GATE0_OUTPUT", "")

struct Gate0LatencyCard{F}
    id::String
    label::String
    kind::String
    operation::F
    parameters::Dict{String,Any}
    absolute_p99_ns::Int64
    max_alloc_bytes::Int64
end

@inline run_gate0_card!(card::Gate0LatencyCard) = card.operation()

function gate0_opd_ramp!(tel::Telescope)
    n = tel.params.resolution
    @inbounds for j in axes(tel.state.opd, 2), i in axes(tel.state.opd, 1)
        tel.state.opd[i, j] = 8e-9 * (i - 1) / n -
            5e-9 * (j - 1) / n
    end
    return tel
end

function gate0_card_parameters(raw::AbstractDict)
    excluded = Set(("id", "label", "absolute_p99_ns", "max_alloc_bytes"))
    return Dict(String(key) => value for (key, value) in raw
        if !(String(key) in excluded))
end

function make_gate0_card(raw::AbstractDict)
    id = String(raw["id"])
    label = String(raw["label"])
    kind = String(raw["kind"])
    resolution = Int(raw["resolution"])
    absolute_p99_ns = Int64(raw["absolute_p99_ns"])
    max_alloc_bytes = Int64(raw["max_alloc_bytes"])
    parameters = gate0_card_parameters(raw)

    operation = if kind == "electric_field"
        zero_padding = Int(raw["zero_padding"])
        tel = Telescope(resolution=resolution, diameter=8.0,
            central_obstruction=0.2)
        gate0_opd_ramp!(tel)
        src = Source(band=:I, magnitude=1.0)
        wavefront = PupilFunction(tel)
        apply_opd!(wavefront, opd_map(tel))
        field = ElectricField(wavefront, src; zero_padding=zero_padding)
        plan = prepare_pupil_field(tel, wavefront, src, field)
        let field=field, wavefront=wavefront, plan=plan
            () -> fill_electric_field!(field, wavefront, plan)
        end
    elseif kind == "direct_psf"
        # `kind` is a stable historical benchmark-card identifier retained so
        # post-refactor results remain comparable with the frozen pre-HIL run.
        zero_padding = Int(raw["zero_padding"])
        tel = Telescope(resolution=resolution, diameter=8.0,
            central_obstruction=0.2)
        gate0_opd_ramp!(tel)
        src = Source(band=:I, magnitude=1.0)
        wavefront = PupilFunction(tel)
        apply_opd!(wavefront, opd_map(tel))
        prepared = prepare_direct_imaging(tel, wavefront, src;
            zero_padding=zero_padding)
        let prepared=prepared
            () -> form_direct_image!(prepared)
        end
    elseif kind == "spatial_filter"
        zero_padding = Int(raw["zero_padding"])
        tel = Telescope(resolution=resolution, diameter=8.0,
            central_obstruction=0.0)
        gate0_opd_ramp!(tel)
        src = Source(band=:I, magnitude=1.0)
        spatial_filter = SpatialFilter(tel; shape=SquareFilter(),
            diameter=resolution / 3, zero_padding=zero_padding)
        wavefront = PupilFunction(tel)
        apply_opd!(wavefront, opd_map(tel))
        field = ElectricField(wavefront, src; zero_padding=zero_padding,
            normalization=DimensionlessNormalization(),
            spatial_measure=PointSampledMeasure(),
            coherence=CoherentFieldCombination())
        formation = prepare_pupil_field(tel, wavefront, src, field;
            center_even_grid=false, amplitude_scale=1)
        fill_electric_field!(field, wavefront, formation)
        output = PupilFunction(tel)
        plan = prepare_spatial_filter(tel, spatial_filter, field, output)
        workspace = SpatialFilterWorkspace(spatial_filter)
        let output=output, field=field, spatial_filter=spatial_filter,
            plan=plan, workspace=workspace
            () -> filter!(output, field, spatial_filter, plan, workspace)
        end
    elseif kind == "atmosphere_direction"
        tel = Telescope(resolution=resolution, diameter=8.0,
            central_obstruction=0.0)
        atmosphere = MultiLayerAtmosphere(tel; r0=0.2, L0=25.0,
            fractional_cn2=[0.6, 0.4], wind_speed=[7.0, 13.0],
            wind_direction=[0.0, 120.0], altitude=[0.0, 6000.0])
        src = Source(band=:I, magnitude=0.0, coordinates=(3.0, 45.0))
        rng = runtime_rng(Int(raw["rng_seed"]))
        duration = 1e-3
        renderer = prepare_atmosphere_renderer(atmosphere, tel, src)
        output = PupilFunction(tel)
        let atmosphere=atmosphere, renderer=renderer, output=output, rng=rng,
            duration=duration
            () -> begin
                epoch = advance_by!(atmosphere, duration; rng=rng)
                render_atmosphere!(output, renderer, atmosphere, epoch)
            end
        end
    elseif kind == "shack_hartmann"
        tel = Telescope(resolution=resolution, diameter=8.0,
            central_obstruction=0.0)
        gate0_opd_ramp!(tel)
        src = Source(band=:I, magnitude=0.0)
        wfs = ShackHartmannWFS(tel; n_lenslets=Int(raw["n_lenslets"]),
            n_pix_subap=Int(raw["n_pix_subap"]), mode=Diffractive())
        let wfs=wfs, tel=tel, src=src
            () -> measure!(wfs, tel, src)
        end
    elseif kind == "pyramid"
        tel = Telescope(resolution=resolution, diameter=8.0,
            central_obstruction=0.0)
        gate0_opd_ramp!(tel)
        src = Source(band=:I, magnitude=0.0)
        pupil = PupilFunction(tel)
        apply_opd!(pupil, opd_map(tel))
        wfs = PyramidWFS(tel; pupil_samples=Int(raw["pupil_samples"]),
            modulation=3.0,
            modulation_points=Int(raw["modulation_points"]),
            mode=Diffractive())
        front_end = PyramidOpticalFrontEnd(wfs, src)
        rate = pyramid_rate_map(front_end, pupil)
        optical_plan = prepare_wfs_optical_formation(front_end, pupil, rate)
        detector = Detector(noise=NoiseNone(), integration_time=1.0,
            qe=1.0, response_model=NullFrameResponse())
        observation = WFSObservation(similar(rate.values);
            units=:electron_count, layout=:four_pupil_mosaic)
        acquisition_plan = prepare_wfs_acquisition(detector, rate,
            observation)
        set_pyramid_calibration!(wfs,
            zeros(size(wfs.estimator.state.reference_signal_2d));
            wavelength_m=wavelength(src), signature=UInt(0x47305036))
        measurement = WFSMeasurement(similar(slopes(wfs));
            units=:dimensionless, kind=:differential_slopes)
        estimator_plan = prepare_wfs_estimation(wfs, observation,
            measurement)
        rng = runtime_rng(61)
        let rate=rate, pupil=pupil, optical_plan=optical_plan,
            observation=observation, acquisition_plan=acquisition_plan,
            measurement=measurement, estimator_plan=estimator_plan, rng=rng
            () -> begin
                form_wfs_optical_products!(rate, pupil, optical_plan)
                acquire_wfs_observation!(observation, rate,
                    acquisition_plan, rng)
                estimate_wfs_measurement!(measurement, observation,
                    estimator_plan)
            end
        end
    elseif kind == "direct_science_two_detectors"
        zero_padding = Int(raw["zero_padding"])
        tel = Telescope(resolution=resolution, diameter=8.0,
            central_obstruction=0.2)
        gate0_opd_ramp!(tel)
        src = Source(band=:I, magnitude=1.0, coordinates=(0.08, 90.0))
        pupil = PupilFunction(tel)
        apply_opd!(pupil, opd_map(tel))
        imaging = prepare_direct_imaging(tel, pupil, src;
            zero_padding=zero_padding)
        rate_map = direct_imaging_output(imaging)
        detector_a = Detector(noise=NoiseNone(), integration_time=0.003,
            qe=0.8)
        detector_b = Detector(noise=NoiseNone(), integration_time=0.007,
            qe=0.8)
        rng_a = runtime_rng(Int(raw["rng_seed_a"]))
        rng_b = runtime_rng(Int(raw["rng_seed_b"]))
        acquisition_a = prepare_detector_acquisition(detector_a, rate_map)
        acquisition_b = prepare_detector_acquisition(detector_b, rate_map)
        let imaging=imaging, rate_map=rate_map, detector_a=detector_a,
            detector_b=detector_b, acquisition_a=acquisition_a,
            acquisition_b=acquisition_b, rng_a=rng_a, rng_b=rng_b
            () -> begin
                form_direct_image!(imaging)
                capture!(detector_a, rate_map, acquisition_a, rng_a)
                capture!(detector_b, rate_map, acquisition_b, rng_b)
            end
        end
    elseif kind == "closed_loop_step"
        tel = Telescope(resolution=resolution, diameter=8.0,
            central_obstruction=0.0)
        src = Source(band=:I, magnitude=0.0)
        atmosphere = KolmogorovAtmosphere(tel; r0=0.2, L0=25.0)
        dm = DeformableMirror(tel; n_act=Int(raw["n_actuators"]),
            influence_width=0.3)
        wfs = ShackHartmannWFS(tel;
            n_lenslets=Int(raw["n_lenslets"]))
        simulation = AOSimulation(tel, src, atmosphere, dm, wfs)
        interaction = interaction_matrix(dm, wfs, tel; amplitude=0.05)
        reconstructor = ModalReconstructor(interaction; gain=0.5)
        runtime = AdaptiveOpticsSim.ClosedLoopRuntime(simulation, reconstructor;
            atmosphere_step=1e-3,
            rng=runtime_rng(Int(raw["rng_seed"])))
        let runtime=runtime
            () -> step!(runtime)
        end
    elseif kind == "zernike"
        tel = Telescope(resolution=resolution, diameter=8.0,
            central_obstruction=0.0)
        gate0_opd_ramp!(tel)
        src = Source(band=:I, magnitude=0.0)
        pupil = PupilFunction(tel)
        apply_opd!(pupil, opd_map(tel))
        wfs = ZernikeWFS(tel;
            pupil_samples=Int(raw["pupil_samples"]))
        front_end = ZernikeOpticalFrontEnd(wfs, src)
        rate = zernike_rate_map(front_end, pupil)
        optical_plan = prepare_wfs_optical_formation(front_end, pupil,
            rate)
        detector = Detector(noise=NoiseNone(), integration_time=1.0,
            qe=1.0, response_model=NullFrameResponse())
        observation = WFSObservation(similar(rate.values);
            units=:electron_count, layout=:zernike_pupil_image)
        acquisition_plan = prepare_wfs_acquisition(detector, rate,
            observation)
        set_zernike_calibration!(wfs,
            zeros(size(wfs.estimator.state.reference_signal_2d));
            wavelength_m=wavelength(src), signature=UInt(0x47305039))
        measurement = WFSMeasurement(similar(slopes(wfs));
            units=:dimensionless, kind=:normalized_pupil_signal)
        estimator_plan = prepare_wfs_estimation(wfs, observation,
            measurement; source=src)
        rng = runtime_rng(Int(raw["rng_seed"]))
        let rate=rate, pupil=pupil, optical_plan=optical_plan,
            observation=observation, acquisition_plan=acquisition_plan,
            measurement=measurement, estimator_plan=estimator_plan, rng=rng
            () -> begin
                form_wfs_optical_products!(rate, pupil, optical_plan)
                acquire_wfs_observation!(observation, rate,
                    acquisition_plan, rng)
                estimate_wfs_measurement!(measurement, observation,
                    estimator_plan)
            end
        end
    elseif kind == "curvature_two_detectors"
        tel = Telescope(resolution=resolution, diameter=8.0,
            central_obstruction=0.0)
        gate0_opd_ramp!(tel)
        src = Source(band=:I, magnitude=0.0)
        pupil = PupilFunction(tel)
        apply_opd!(pupil, opd_map(tel))
        wfs = CurvatureWFS(tel;
            pupil_samples=Int(raw["pupil_samples"]))
        front_end = CurvatureOpticalFrontEnd(wfs, src)
        rates = curvature_rate_maps(front_end, pupil)
        optical_plan = prepare_wfs_optical_formation(front_end, pupil,
            rates)
        plus_detector = Detector(noise=NoiseNone(), integration_time=0.5,
            qe=1.0, response_model=NullFrameResponse())
        minus_detector = Detector(noise=NoiseNone(), integration_time=1.0,
            qe=1.0, response_model=NullFrameResponse())
        plus_observation = WFSObservation(similar(rates[1].values);
            units=:electron_count, layout=:curvature_branch_image)
        minus_observation = WFSObservation(similar(rates[2].values);
            units=:electron_count, layout=:curvature_branch_image)
        observations = (plus_observation, minus_observation)
        acquisition_plan = prepare_wfs_acquisition(
            (plus_detector, minus_detector), rates, observations)
        set_curvature_calibration!(wfs,
            zeros(size(wfs.estimator.state.reference_signal_2d));
            wavelength_m=wavelength(src), signature=UInt(0x4730503a))
        measurement = WFSMeasurement(similar(slopes(wfs));
            units=:dimensionless, kind=:curvature_signal)
        estimator_plan = prepare_wfs_estimation(wfs, observations,
            measurement; branch_rate_scales=(2.0, 1.0))
        rng = runtime_rng(Int(raw["rng_seed"]))
        let rates=rates, pupil=pupil, optical_plan=optical_plan,
            observations=observations, acquisition_plan=acquisition_plan,
            measurement=measurement, estimator_plan=estimator_plan, rng=rng
            () -> begin
                form_wfs_optical_products!(rates, pupil, optical_plan)
                acquire_wfs_observation!(observations, rates,
                    acquisition_plan, rng)
                estimate_wfs_measurement!(measurement, observations,
                    estimator_plan)
            end
        end
    else
        throw(ArgumentError("unsupported Gate 0 latency-card kind '$kind'"))
    end

    return Gate0LatencyCard(id, label, kind, operation, parameters,
        absolute_p99_ns, max_alloc_bytes)
end

function configure_gate0_benchmark!()
    Threads.nthreads() == 1 || error(
        "Gate 0 latency cards require one Julia thread")
    BLAS.set_num_threads(1)
    AdaptiveOpticsSim.set_fft_provider_threads!(1)
    return nothing
end

function gc_delta(before, after)
    return Dict(String(name) => Int64(
        getfield(after, name) - getfield(before, name))
        for name in fieldnames(typeof(before)))
end

function gate0_histogram_summary(histogram::HdrHistogram.Histogram,
    wall_start_ns::UInt64, wall_end_ns::UInt64, samples::Int, gc_counters)
    encoded = HdrHistogram.encode_into_compressed_byte_buffer(histogram)
    decoded = HdrHistogram.decode_from_compressed_byte_buffer(encoded)
    HdrHistogram.total_count(decoded) == samples || error(
        "encoded Gate 0 histogram lost samples")
    HdrHistogram.value_at_percentile(decoded, 99.0) ==
        HdrHistogram.value_at_percentile(histogram, 99.0) || error(
        "encoded Gate 0 histogram changed p99")
    wall_ns = Int64(wall_end_ns - wall_start_ns)
    return Dict{String,Any}(
        "samples" => samples,
        "completed_operations" => samples,
        "failed_operations" => 0,
        "monotonic_start_ns" => wall_start_ns,
        "monotonic_end_ns" => wall_end_ns,
        "wall_ns" => wall_ns,
        "throughput_hz" => 1.0e9 * samples / wall_ns,
        "min_ns" => min(histogram),
        "mean_ns" => HdrHistogram.mean(histogram),
        "p50_ns" => HdrHistogram.value_at_percentile(histogram, 50.0),
        "p90_ns" => HdrHistogram.value_at_percentile(histogram, 90.0),
        "p99_ns" => HdrHistogram.value_at_percentile(histogram, 99.0),
        "p99_9_ns" => HdrHistogram.value_at_percentile(histogram, 99.9),
        "max_ns" => max(histogram),
        "gc" => gc_counters,
        "histogram_encoding" => "HdrHistogram compressed V2 base64",
        "histogram_base64" => base64encode(encoded),
    )
end

function measure_gate0_run!(card::Gate0LatencyCard, samples::Int,
    lowest_ns::Int64, highest_ns::Int64, significant_figures::Int)
    histogram = HdrHistogram.Histogram(lowest_ns, highest_ns,
        significant_figures)
    GC.gc()
    gc_before = Base.gc_num()
    wall_start = time_ns()
    @inbounds for _ in 1:samples
        start = time_ns()
        run_gate0_card!(card)
        elapsed = Int64(time_ns() - start)
        HdrHistogram.record_value!(histogram, elapsed)
    end
    wall_end = time_ns()
    gc_after = Base.gc_num()
    return gate0_histogram_summary(histogram, wall_start, wall_end, samples,
        gc_delta(gc_before, gc_after))
end

median_integer(values) = round(Int64, median(collect(values)))

function summarize_gate0_runs(runs::Vector{Dict{String,Any}})
    return Dict{String,Any}(
        "median_p50_ns" => median_integer(run["p50_ns"] for run in runs),
        "median_p90_ns" => median_integer(run["p90_ns"] for run in runs),
        "median_p99_ns" => median_integer(run["p99_ns"] for run in runs),
        "median_p99_9_ns" => median_integer(
            run["p99_9_ns"] for run in runs),
        "worst_p99_ns" => maximum(run["p99_ns"] for run in runs),
        "worst_p99_9_ns" => maximum(run["p99_9_ns"] for run in runs),
        "median_throughput_hz" => median(collect(
            run["throughput_hz"] for run in runs)),
    )
end

function load_gate0_baselines(path::AbstractString)
    isempty(path) && return Dict{String,Any}()
    artifact = TOML.parsefile(path)
    return Dict(String(card["id"]) => card for card in artifact["cards"])
end

function gate0_regression(card::Gate0LatencyCard, summary,
    steady_alloc_bytes::Int64, relative_factor::Float64, baselines,
    p99_gate_supported::Bool)
    absolute_observed = Int64(summary["worst_p99_ns"])
    relative_observed = Int64(summary["median_p99_ns"])
    allocation_passed = steady_alloc_bytes <= card.max_alloc_bytes
    absolute_passed = !p99_gate_supported ||
        absolute_observed <= card.absolute_p99_ns
    result = Dict{String,Any}(
        "allocation_gate_passed" => allocation_passed,
        "allocation_limit_bytes" => card.max_alloc_bytes,
        "absolute_p99_gate_evaluated" => p99_gate_supported,
        "absolute_p99_gate_passed" => absolute_passed,
        "absolute_p99_observed_ns" => absolute_observed,
        "absolute_p99_limit_ns" => card.absolute_p99_ns,
        "relative_p99_gate_evaluated" => false,
        "relative_p99_gate_passed" => true,
    )
    p99_gate_supported || return result
    haskey(baselines, card.id) || return result
    baseline = baselines[card.id]
    baseline["kind"] == card.kind || error(
        "baseline kind mismatch for $(card.id)")
    baseline["parameters"] == card.parameters || error(
        "baseline parameters mismatch for $(card.id)")
    baseline_p99 = Int64(baseline["summary"]["median_p99_ns"])
    relative_limit = ceil(Int64, relative_factor * baseline_p99)
    result["relative_p99_gate_evaluated"] = true
    result["relative_p99_gate_passed"] = relative_observed <= relative_limit
    result["relative_p99_observed_ns"] = relative_observed
    result["baseline_median_p99_ns"] = baseline_p99
    result["relative_p99_limit_ns"] = relative_limit
    result["relative_p99_factor"] = relative_factor
    return result
end

function command_output(command)
    try
        return readchomp(command)
    catch
        return "unknown"
    end
end

function optional_file(path::AbstractString)
    try
        return strip(read(path, String))
    catch
        return "unknown"
    end
end

function allowed_cpu_list()
    status = optional_file("/proc/self/status")
    status == "unknown" && return status
    for line in eachline(IOBuffer(status))
        startswith(line, "Cpus_allowed_list:") || continue
        return strip(only(split(line, ':'; limit=2)[2:2]))
    end
    return "unknown"
end

function gate0_environment()
    cpu = first(Sys.cpu_info())
    git_status = command_output(`git status --porcelain=v1`)
    return Dict{String,Any}(
        "timestamp_utc" => string(Dates.now(Dates.UTC)),
        "git_commit" => command_output(`git rev-parse HEAD`),
        "git_dirty" => !isempty(git_status) && git_status != "unknown",
        "git_status_porcelain" => git_status,
        "julia_version" => string(VERSION),
        "adaptive_optics_sim_version" => string(
            Base.pkgversion(AdaptiveOpticsSim)),
        "hdrhistogram_version" => string(Base.pkgversion(HdrHistogram)),
        "active_project" => something(Base.active_project(), "unknown"),
        "kernel" => string(Sys.KERNEL),
        "kernel_release" => command_output(`uname -r`),
        "architecture" => string(Sys.ARCH),
        "cpu_target" => string(Sys.CPU_NAME),
        "cpu_model" => cpu.model,
        "logical_cpu_threads" => Sys.CPU_THREADS,
        "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
        "blas_config" => string(BLAS.get_config()),
        "fft_threads" => 1,
        "julia_cpu_target_env" => get(ENV, "JULIA_CPU_TARGET", "default"),
        "allowed_cpus" => allowed_cpu_list(),
        "scaling_governor_cpu0" => optional_file(
            "/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor"),
        "command" => "julia --project=benchmarks benchmarks/benchmark_gate0_latency.jl",
    )
end

function timer_recording_overhead(samples::Int, lowest_ns::Int64,
    highest_ns::Int64, significant_figures::Int)
    histogram = HdrHistogram.Histogram(lowest_ns, highest_ns,
        significant_figures)
    wall_start = time_ns()
    @inbounds for _ in 1:samples
        start = time_ns()
        elapsed = Int64(time_ns() - start)
        HdrHistogram.record_value!(histogram, elapsed)
    end
    wall_ns = Int64(time_ns() - wall_start)
    return Dict{String,Any}(
        "samples" => samples,
        "p50_clock_delta_ns" => HdrHistogram.value_at_percentile(
            histogram, 50.0),
        "p99_clock_delta_ns" => HdrHistogram.value_at_percentile(
            histogram, 99.0),
        "wall_ns_per_clock_pair_and_record" => wall_ns / samples,
    )
end

function write_gate0_artifact(path::AbstractString, artifact)
    isempty(path) && return nothing
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, artifact; sorted=true)
    end
    println("wrote Gate 0 latency artifact: ", path)
    return nothing
end

function select_gate0_cards(cards::Tuple)
    raw_ids = strip(get(ENV, "AOS_GATE0_CARD_IDS", ""))
    isempty(raw_ids) && return cards, "all"

    requested_ids = strip.(split(raw_ids, ','))
    any(isempty, requested_ids) && error(
        "AOS_GATE0_CARD_IDS must be a comma-separated list of nonempty card IDs")
    length(unique(requested_ids)) == length(requested_ids) || error(
        "AOS_GATE0_CARD_IDS must not contain duplicate card IDs")

    cards_by_id = Dict(card.id => card for card in cards)
    unknown_ids = filter(id -> !haskey(cards_by_id, id), requested_ids)
    isempty(unknown_ids) || error(
        "AOS_GATE0_CARD_IDS contains unknown card IDs: $(join(unknown_ids, ", "))")
    return Tuple(cards_by_id[id] for id in requested_ids), "explicit"
end

function run_gate0_latency_benchmarks()
    configure_gate0_benchmark!()
    contract = TOML.parsefile(GATE0_CONTRACT_PATH)
    samples = parse(Int, get(ENV, "AOS_GATE0_SAMPLES",
        string(contract["samples_per_run"])))
    runs_count = parse(Int, get(ENV, "AOS_GATE0_RUNS",
        string(contract["runs"])))
    warmup = parse(Int, get(ENV, "AOS_GATE0_WARMUP",
        string(contract["warmup_operations"])))
    samples > 0 || error("AOS_GATE0_SAMPLES must be > 0")
    runs_count > 0 || error("AOS_GATE0_RUNS must be > 0")
    warmup >= 0 || error("AOS_GATE0_WARMUP must be >= 0")
    minimum_p99_samples = Int(contract["minimum_samples_for_p99_gate"])
    p99_gate_supported = samples >= minimum_p99_samples
    lowest_ns = Int64(contract["histogram_lowest_ns"])
    highest_ns = Int64(contract["histogram_highest_ns"])
    significant_figures = Int(contract["histogram_significant_figures"])
    relative_factor = Float64(contract["relative_p99_factor"])
    baselines = load_gate0_baselines(GATE0_BASELINE_PATH)
    all_cards = Tuple(make_gate0_card(raw) for raw in contract["cards"])
    cards, card_selection = select_gate0_cards(all_cards)
    results = Vector{Dict{String,Any}}()
    all_gates_passed = true

    println("pre_hil_gate0_latency_contract")
    println("  load_model: warmed serial closed-loop")
    println("  samples_per_run: ", samples)
    println("  runs: ", runs_count)
    println("  warmup_operations: ", warmup)
    println("  p99_relative_factor: ", relative_factor)
    println("  card_selection: ", card_selection)
    println("  card_ids: ", join((card.id for card in cards), ","))
    if !p99_gate_supported
        println("  p99_gates: skipped (requires at least ",
            minimum_p99_samples, " samples per run)")
    end

    for card in cards
        GC.gc()
        first_start = time_ns()
        run_gate0_card!(card)
        first_use_ns = Int64(time_ns() - first_start)
        for _ in 1:warmup
            run_gate0_card!(card)
        end
        GC.gc()
        steady_alloc_bytes = Int64(@allocated run_gate0_card!(card))
        runs = Vector{Dict{String,Any}}(undef, runs_count)
        for run_index in 1:runs_count
            runs[run_index] = measure_gate0_run!(card, samples,
                lowest_ns, highest_ns, significant_figures)
            println(card.id, " run=", run_index,
                " p50_ns=", runs[run_index]["p50_ns"],
                " p99_ns=", runs[run_index]["p99_ns"])
        end
        summary = summarize_gate0_runs(runs)
        regression = gate0_regression(card, summary, steady_alloc_bytes,
            relative_factor, baselines, p99_gate_supported)
        card_passed = regression["allocation_gate_passed"] &&
            regression["absolute_p99_gate_passed"] &&
            regression["relative_p99_gate_passed"]
        all_gates_passed &= card_passed
        push!(results, Dict{String,Any}(
            "id" => card.id,
            "label" => card.label,
            "kind" => card.kind,
            "parameters" => card.parameters,
            "first_use_ns" => first_use_ns,
            "steady_alloc_bytes" => steady_alloc_bytes,
            "runs" => runs,
            "summary" => summary,
            "regression" => regression,
        ))
    end

    environment = gate0_environment()
    if !isempty(GATE0_OUTPUT_PATH) && environment["git_dirty"]
        error("refusing to write durable Gate 0 evidence from a dirty worktree")
    end
    artifact = Dict{String,Any}(
        "schema_version" => 1,
        "benchmark" => String(contract["name"]),
        "evidence_class" => String(contract["evidence_class"]),
        "source_contract" => relpath(abspath(GATE0_CONTRACT_PATH), @__DIR__),
        "characterized_source_revision" => environment["git_commit"],
        "baseline_source_revision" => String(contract["source_revision"]),
        "configured_samples_per_run" => samples,
        "configured_runs" => runs_count,
        "configured_warmup_operations" => warmup,
        "configured_card_selection" => card_selection,
        "configured_card_ids" => [card.id for card in cards],
        "minimum_samples_for_p99_gate" => minimum_p99_samples,
        "p99_gate_supported" => p99_gate_supported,
        "contract" => contract["contract"],
        "scope_exclusions" => contract["scope_exclusions"],
        "histogram" => Dict(
            "lowest_ns" => lowest_ns,
            "highest_ns" => highest_ns,
            "significant_figures" => significant_figures,
        ),
        "timer_recording_overhead" => timer_recording_overhead(
            min(samples, 100_000), lowest_ns, highest_ns,
            significant_figures),
        "environment" => environment,
        "cards" => results,
        "all_gates_passed" => all_gates_passed,
    )
    write_gate0_artifact(GATE0_OUTPUT_PATH, artifact)
    all_gates_passed || error(
        "one or more pre-HIL Gate 0 latency gates failed")
    return artifact
end

run_gate0_latency_benchmarks()

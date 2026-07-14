using AdaptiveOpticsSim
using Dates
using HdrHistogram
using LinearAlgebra
using Random
using Statistics
using TOML

const DETECTOR_HIL_SAMPLES = parse(Int,
    get(ENV, "AOS_DETECTOR_HIL_SAMPLES", "100000"))
const DETECTOR_HIL_WARMUP = parse(Int,
    get(ENV, "AOS_DETECTOR_HIL_WARMUP", "1000"))
const DETECTOR_HIL_RUNS = parse(Int,
    get(ENV, "AOS_DETECTOR_HIL_RUNS", "3"))
const DETECTOR_HIL_SIZE = parse(Int,
    get(ENV, "AOS_DETECTOR_HIL_SIZE", "128"))
const DETECTOR_HIL_BASELINE = get(ENV, "AOS_DETECTOR_HIL_BASELINE", "")
const DETECTOR_HIL_OUTPUT = get(ENV, "AOS_DETECTOR_HIL_OUTPUT", "")
const DETECTOR_HIL_REGRESSION_FACTOR = parse(Float64,
    get(ENV, "AOS_DETECTOR_HIL_REGRESSION_FACTOR", "1.25"))

const DETECTOR_HIL_HIGHEST_NS = Int64(60_000_000_000)
const DETECTOR_HIL_SIGNIFICANT_FIGURES = 3

struct DetectorHILCard{D,A<:AbstractMatrix,R<:AbstractRNG}
    id::String
    label::String
    detector::D
    input::A
    rng::R
end

function configure_detector_hil_benchmark!()
    Threads.nthreads() == 1 || error(
        "detector HIL latency cards require one Julia thread")
    BLAS.set_num_threads(1)
    AdaptiveOpticsSim.set_fft_provider_threads!(1)
    DETECTOR_HIL_SAMPLES > 0 || error("AOS_DETECTOR_HIL_SAMPLES must be > 0")
    DETECTOR_HIL_WARMUP >= 0 || error("AOS_DETECTOR_HIL_WARMUP must be >= 0")
    DETECTOR_HIL_RUNS > 0 || error("AOS_DETECTOR_HIL_RUNS must be > 0")
    DETECTOR_HIL_SIZE > 0 || error("AOS_DETECTOR_HIL_SIZE must be > 0")
    DETECTOR_HIL_REGRESSION_FACTOR >= 1 || error(
        "AOS_DETECTOR_HIL_REGRESSION_FACTOR must be >= 1")
    return nothing
end

function detector_hil_input(::Type{T}, n::Int) where {T<:AbstractFloat}
    input = Matrix{T}(undef, n, n)
    @inbounds for j in axes(input, 2), i in axes(input, 1)
        input[i, j] = T(45_000 + 5_000 * sinpi(T(i + j) / T(n)))
    end
    return input
end

function detector_hil_cards(n::Int=DETECTOR_HIL_SIZE)
    T = Float32
    input = detector_hil_input(T, n)
    ipc = InterpixelCapacitance(T[
        0.00 0.01 0.00
        0.01 0.96 0.01
        0.00 0.01 0.00
    ]; T=T)

    cmos = Detector(
        integration_time=1e-3,
        qe=0.8,
        noise=NoisePhotonReadout(1.2),
        full_well=30_000,
        bits=16,
        sensor=CMOSSensor(row_readout_sigma=0.08,
            column_readout_sigma=0.10, T=T),
        response_model=NullFrameResponse(),
        T=T,
    )
    cmos_mtf_ipc = Detector(
        integration_time=1e-3,
        qe=0.8,
        noise=NoisePhotonReadout(1.2),
        full_well=30_000,
        bits=16,
        sensor=CMOSSensor(row_readout_sigma=0.08,
            column_readout_sigma=0.10, T=T),
        response_model=GaussianPixelResponse(response_width_px=0.35, T=T),
        charge_coupling_model=ipc,
        T=T,
    )
    ccd = Detector(
        integration_time=1e-3,
        qe=0.9,
        noise=NoisePhotonReadout(2.0),
        full_well=100_000,
        bits=16,
        sensor=CCDSensor(clock_induced_charge_per_frame=0.002, T=T),
        response_model=NullFrameResponse(),
        T=T,
    )
    emccd = Detector(
        integration_time=1e-3,
        qe=0.9,
        noise=NoisePhotonReadout(20.0),
        gain=100,
        full_well=800_000,
        bits=16,
        sensor=EMCCDSensor(excess_noise_factor=sqrt(T(2)),
            clock_induced_charge_per_frame=0.002, T=T),
        response_model=NullFrameResponse(),
        T=T,
    )
    hgcdte = Detector(
        integration_time=1e-3,
        qe=0.7,
        noise=NoisePhotonReadout(1.0),
        full_well=100_000,
        bits=16,
        sensor=HgCdTeAvalancheArraySensor(avalanche_gain=20,
            excess_noise_factor=1.2, glow_rate=0.01, read_time=2e-5,
            sampling_mode=CorrelatedDoubleSampling(), T=T),
        response_model=NullFrameResponse(),
        charge_coupling_model=ipc,
        T=T,
    )
    skipper = Detector(
        integration_time=1e-3,
        qe=0.9,
        noise=NoisePhotonReadout(3.0),
        full_well=100_000,
        bits=16,
        sensor=CCDSensor(read_time=2e-6,
            sampling_mode=SkipperSampling(16), T=T),
        response_model=NullFrameResponse(),
        T=T,
    )

    return (
        DetectorHILCard("DET-HIL-01", "CMOS global-shutter capture", cmos,
            input, runtime_rng(101)),
        DetectorHILCard("DET-HIL-02", "CMOS capture with MTF and IPC",
            cmos_mtf_ipc, input, runtime_rng(102)),
        DetectorHILCard("DET-HIL-03", "CCD capture", ccd, input,
            runtime_rng(103)),
        DetectorHILCard("DET-HIL-04", "EMCCD fast linear capture", emccd,
            input, runtime_rng(104)),
        DetectorHILCard("DET-HIL-05", "HgCdTe avalanche CDS capture", hgcdte,
            input, runtime_rng(105)),
        DetectorHILCard("DET-HIL-06", "Skipper CCD 16-sample capture", skipper,
            input, runtime_rng(106)),
    )
end

@inline detector_hil_capture!(card::DetectorHILCard) =
    capture!(card.detector, card.input, card.rng)

function detector_hil_first_capture_ns(card::DetectorHILCard)
    GC.gc()
    start = time_ns()
    detector_hil_capture!(card)
    return Int64(time_ns() - start)
end

function detector_hil_steady_alloc_bytes(card::DetectorHILCard)
    for _ in 1:10
        detector_hil_capture!(card)
    end
    GC.gc()
    return @allocated detector_hil_capture!(card)
end

function gc_counter_delta(before, after)
    return Dict(String(name) => Int64(getfield(after, name) - getfield(before, name))
        for name in fieldnames(typeof(before)))
end

function histogram_summary(histogram::HdrHistogram.Histogram, wall_ns::Int64,
    samples::Int, gc_delta)
    return Dict{String,Any}(
        "samples" => samples,
        "wall_ns" => wall_ns,
        "throughput_hz" => 1.0e9 * samples / wall_ns,
        "min_ns" => min(histogram),
        "mean_ns" => HdrHistogram.mean(histogram),
        "p50_ns" => HdrHistogram.value_at_percentile(histogram, 50.0),
        "p90_ns" => HdrHistogram.value_at_percentile(histogram, 90.0),
        "p99_ns" => HdrHistogram.value_at_percentile(histogram, 99.0),
        "p99_9_ns" => HdrHistogram.value_at_percentile(histogram, 99.9),
        "max_ns" => max(histogram),
        "gc" => gc_delta,
    )
end

function measure_detector_hil_run!(card::DetectorHILCard, samples::Int)
    histogram = HdrHistogram.Histogram(1, DETECTOR_HIL_HIGHEST_NS,
        DETECTOR_HIL_SIGNIFICANT_FIGURES)
    GC.gc()
    gc_before = Base.gc_num()
    wall_start = time_ns()
    @inbounds for _ in 1:samples
        start = time_ns()
        detector_hil_capture!(card)
        elapsed = Int64(time_ns() - start)
        HdrHistogram.record_value!(histogram, elapsed)
    end
    wall_ns = Int64(time_ns() - wall_start)
    gc_after = Base.gc_num()
    return histogram_summary(histogram, wall_ns, samples,
        gc_counter_delta(gc_before, gc_after))
end

function print_detector_hil_run(card::DetectorHILCard, run_index::Int,
    summary::Dict{String,Any})
    println(card.id, " run=", run_index, " label=", card.label)
    println("  p50_ns: ", summary["p50_ns"])
    println("  p90_ns: ", summary["p90_ns"])
    println("  p99_ns: ", summary["p99_ns"])
    println("  p99_9_ns: ", summary["p99_9_ns"])
    println("  max_ns: ", summary["max_ns"])
    println("  throughput_hz: ", summary["throughput_hz"])
    return nothing
end

function median_integer(values)
    return round(Int64, median(collect(values)))
end

function summarize_detector_hil_runs(runs::Vector{Dict{String,Any}})
    return Dict{String,Any}(
        "median_p50_ns" => median_integer(run["p50_ns"] for run in runs),
        "median_p90_ns" => median_integer(run["p90_ns"] for run in runs),
        "median_p99_ns" => median_integer(run["p99_ns"] for run in runs),
        "median_p99_9_ns" => median_integer(run["p99_9_ns"] for run in runs),
        "worst_p99_ns" => maximum(run["p99_ns"] for run in runs),
        "worst_p99_9_ns" => maximum(run["p99_9_ns"] for run in runs),
        "median_throughput_hz" => median(collect(
            run["throughput_hz"] for run in runs)),
    )
end

function detector_card_metadata(card::DetectorHILCard)
    metadata = detector_export_metadata(card.detector)
    return Dict{String,Any}(
        "sensor" => String(metadata.sensor),
        "noise" => String(metadata.noise),
        "frame_size" => collect(metadata.frame_size),
        "output_size" => collect(metadata.output_size),
        "frame_response" => String(metadata.frame_response),
        "charge_coupling" => String(metadata.charge_coupling),
        "sampling_mode" => String(metadata.sampling_mode),
        "sampling_reads" => something(metadata.sampling_reads, 0),
        "output_type" => string(metadata.output_type),
    )
end

function baseline_cards(path::AbstractString)
    isempty(path) && return Dict{String,Any}()
    data = TOML.parsefile(path)
    cards = get(data, "cards", Any[])
    return Dict(String(card["id"]) => card for card in cards)
end

function require_compatible_detector_baseline!(card_id::String,
    current::Dict{String,Any}, baseline::Dict{String,Any})
    baseline_detector = baseline["detector"]
    for key in ("sensor", "noise", "frame_size", "output_size",
        "frame_response", "charge_coupling", "sampling_mode",
        "sampling_reads", "output_type")
        get(baseline_detector, key, nothing) == current[key] || error(
            "baseline detector configuration mismatch for $(card_id) field $(key)")
    end
    return nothing
end

function detector_hil_regression(card_id::String, metadata::Dict{String,Any},
    summary::Dict{String,Any}, allocation_gate_passed::Bool,
    baselines::AbstractDict)
    result = Dict{String,Any}(
        "allocation_gate_passed" => allocation_gate_passed,
        "latency_gate_evaluated" => false,
        "latency_gate_passed" => true,
    )
    haskey(baselines, card_id) || return result
    baseline = baselines[card_id]
    require_compatible_detector_baseline!(card_id, metadata, baseline)
    baseline_summary = baseline["summary"]
    baseline_p99 = Int64(baseline_summary["median_p99_ns"])
    limit = ceil(Int64, DETECTOR_HIL_REGRESSION_FACTOR * baseline_p99)
    observed = Int64(summary["median_p99_ns"])
    result["latency_gate_evaluated"] = true
    result["latency_gate_passed"] = observed <= limit
    result["baseline_median_p99_ns"] = baseline_p99
    result["p99_limit_ns"] = limit
    result["regression_factor"] = DETECTOR_HIL_REGRESSION_FACTOR
    return result
end

function git_commit()
    try
        return readchomp(`git rev-parse HEAD`)
    catch
        return "unknown"
    end
end

function detector_hil_environment()
    cpu = first(Sys.cpu_info())
    return Dict{String,Any}(
        "timestamp_utc" => string(Dates.now(Dates.UTC)),
        "git_commit" => git_commit(),
        "julia_version" => string(VERSION),
        "adaptive_optics_sim_version" => string(Base.pkgversion(AdaptiveOpticsSim)),
        "hdrhistogram_version" => string(Base.pkgversion(HdrHistogram)),
        "kernel" => string(Sys.KERNEL),
        "architecture" => string(Sys.ARCH),
        "cpu_target" => string(Sys.CPU_NAME),
        "cpu_model" => cpu.model,
        "logical_cpu_threads" => Sys.CPU_THREADS,
        "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
        "blas_config" => string(BLAS.get_config()),
        "fft_threads" => 1,
    )
end

function detector_hil_contract()
    return Dict{String,Any}(
        "boundary" => "input detector frame available to converted output frame ready",
        "load_model" => "warmed serial closed-loop",
        "arrival_model" => "next capture starts after the previous capture completes",
        "coordinated_omission_correction" => false,
        "coordinated_omission_reason" => "no independent arrival schedule is modeled",
        "timer" => "time_ns; timer overhead is included and not subtracted",
        "single_writer" => true,
        "julia_threads" => 1,
        "blas_threads" => 1,
        "fft_threads" => 1,
        "first_use_separate" => true,
        "first_use_note" => "first capture of each constructed card in card order; package loading is excluded",
        "samples_per_run" => DETECTOR_HIL_SAMPLES,
        "runs" => DETECTOR_HIL_RUNS,
        "warmup_captures" => DETECTOR_HIL_WARMUP,
        "tail_resolution_note" => DETECTOR_HIL_SAMPLES >= 100_000 ?
            "at least 100 observations are expected beyond p99.9 per run" :
            "fewer than 100000 samples; p99.9 is diagnostic only",
        "scope_exclusions" => [
            "external RTC transport and scheduling",
            "fixed-rate open-loop arrivals and overload",
            "camera-link or frame-grabber I/O",
            "CPU affinity, isolation, and frequency-governor control",
        ],
    )
end

function write_detector_hil_artifact(path::AbstractString, artifact)
    isempty(path) && return nothing
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, artifact; sorted=true)
    end
    println("wrote detector HIL artifact: ", path)
    return nothing
end

function run_detector_hil_latency_benchmarks()
    configure_detector_hil_benchmark!()
    baselines = baseline_cards(DETECTOR_HIL_BASELINE)
    results = Vector{Dict{String,Any}}()
    all_gates_passed = true

    println("detector_hil_latency_contract")
    println("  load_model: warmed serial closed-loop")
    println("  samples_per_run: ", DETECTOR_HIL_SAMPLES)
    println("  runs: ", DETECTOR_HIL_RUNS)
    println("  frame_size: ", DETECTOR_HIL_SIZE, "x", DETECTOR_HIL_SIZE)
    println("  coordinated_omission_correction: false (no independent arrivals)")

    for card in detector_hil_cards()
        first_capture_ns = detector_hil_first_capture_ns(card)
        for _ in 1:DETECTOR_HIL_WARMUP
            detector_hil_capture!(card)
        end
        steady_alloc_bytes = detector_hil_steady_alloc_bytes(card)
        allocation_gate_passed = steady_alloc_bytes == 0
        runs = Vector{Dict{String,Any}}(undef, DETECTOR_HIL_RUNS)
        for run_index in 1:DETECTOR_HIL_RUNS
            runs[run_index] = measure_detector_hil_run!(card,
                DETECTOR_HIL_SAMPLES)
            print_detector_hil_run(card, run_index, runs[run_index])
        end
        summary = summarize_detector_hil_runs(runs)
        metadata = detector_card_metadata(card)
        regression = detector_hil_regression(card.id, metadata, summary,
            allocation_gate_passed, baselines)
        all_gates_passed &= allocation_gate_passed &&
            regression["latency_gate_passed"]
        println("  first_capture_ns: ", first_capture_ns)
        println("  steady_alloc_bytes: ", steady_alloc_bytes)
        println("  allocation_gate_passed: ", allocation_gate_passed)
        println("  median_p99_ns: ", summary["median_p99_ns"])
        println("  median_p99_9_ns: ", summary["median_p99_9_ns"])
        println("  latency_gate_evaluated: ",
            regression["latency_gate_evaluated"])

        push!(results, Dict{String,Any}(
            "id" => card.id,
            "label" => card.label,
            "first_capture_ns" => first_capture_ns,
            "steady_alloc_bytes" => steady_alloc_bytes,
            "detector" => metadata,
            "runs" => runs,
            "summary" => summary,
            "regression" => regression,
        ))
    end

    artifact = Dict{String,Any}(
        "schema_version" => 1,
        "benchmark" => "detector_hil_latency",
        "evidence_class" => "warmed in-process detector capture latency",
        "contract" => detector_hil_contract(),
        "environment" => detector_hil_environment(),
        "cards" => results,
        "all_gates_passed" => all_gates_passed,
    )
    write_detector_hil_artifact(DETECTOR_HIL_OUTPUT, artifact)
    all_gates_passed || error("one or more detector HIL regression gates failed")
    return artifact
end

run_detector_hil_latency_benchmarks()

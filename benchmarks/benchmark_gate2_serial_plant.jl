using AdaptiveOpticsSim
using Base64
using Dates
using HdrHistogram
using LinearAlgebra
using SHA
using Statistics
using TOML

include(joinpath(@__DIR__, "support", "gate2_serial_plant.jl"))
using .Gate2SerialPlantBenchmark

const GATE2_CONTRACT_PATH = get(ENV, "AOS_GATE2_CONTRACT",
    joinpath(@__DIR__, "contracts", "gate2_serial_plant.toml"))
const GATE2_BASELINE_PATH = get(ENV, "AOS_GATE2_BASELINE", "")
const GATE2_OUTPUT_PATH = get(ENV, "AOS_GATE2_OUTPUT", "")

function configure_gate2_benchmark!()
    Threads.nthreads() == 1 || error(
        "Gate 2 serial-plant evidence requires one Julia thread")
    BLAS.set_num_threads(1)
    AdaptiveOpticsSim.set_fft_provider_threads!(1)
    return nothing
end

function gc_delta(before, after)
    return Dict(String(name) => Int64(
        getfield(after, name) - getfield(before, name))
        for name in fieldnames(typeof(before)))
end

function histogram_summary(histogram::HdrHistogram.Histogram,
    wall_start_ns::UInt64, wall_end_ns::UInt64, samples::Int, gc_counters)
    encoded = HdrHistogram.encode_into_compressed_byte_buffer(histogram)
    decoded = HdrHistogram.decode_from_compressed_byte_buffer(encoded)
    HdrHistogram.total_count(decoded) == samples || error(
        "encoded Gate 2 histogram lost samples")
    HdrHistogram.value_at_percentile(decoded, 99.0) ==
        HdrHistogram.value_at_percentile(histogram, 99.0) || error(
        "encoded Gate 2 histogram changed p99")
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

function measure_run!(operation::SerialPlantOperation, samples::Int,
    lowest_ns::Int64, highest_ns::Int64, significant_figures::Int)
    histogram = HdrHistogram.Histogram(lowest_ns, highest_ns,
        significant_figures)
    GC.gc()
    gc_before = Base.gc_num()
    wall_start = time_ns()
    @inbounds for _ in 1:samples
        start = time_ns()
        operation()
        elapsed = Int64(time_ns() - start)
        HdrHistogram.record_value!(histogram, elapsed)
    end
    wall_end = time_ns()
    gc_after = Base.gc_num()
    return histogram_summary(histogram, wall_start, wall_end, samples,
        gc_delta(gc_before, gc_after))
end

median_integer(values) = round(Int64, median(collect(values)))

function summarize_runs(runs::Vector{Dict{String,Any}})
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

function load_baseline(path::AbstractString, benchmark_name::AbstractString,
    contract_definition, workload,
    relative_factor::Float64, observed_p99_ns::Int64,
    p99_supported::Bool)
    result = Dict{String,Any}(
        "relative_p99_gate_evaluated" => false,
        "relative_p99_gate_passed" => true,
    )
    isempty(strip(path)) && return result, Dict{String,Any}()
    baseline = TOML.parsefile(path)
    baseline["benchmark"] == benchmark_name || error(
        "Gate 2 baseline benchmark identity does not match")
    baseline["contract"] == contract_definition || error(
        "Gate 2 baseline timed-boundary contract does not match")
    baseline["workload"] == workload || error(
        "Gate 2 baseline workload does not match the current contract")
    source = Dict{String,Any}(
        "path" => relpath(abspath(path), @__DIR__),
        "sha256" => bytes2hex(SHA.sha256(read(path))),
        "characterized_source_revision" => get(baseline,
            "characterized_source_revision", "unknown"),
    )
    p99_supported || return result, source
    baseline_p99_ns = Int64(baseline["summary"]["median_p99_ns"])
    relative_limit_ns = ceil(Int64, relative_factor * baseline_p99_ns)
    result["relative_p99_gate_evaluated"] = true
    result["relative_p99_gate_passed"] = observed_p99_ns <=
        relative_limit_ns
    result["relative_p99_observed_ns"] = observed_p99_ns
    result["baseline_median_p99_ns"] = baseline_p99_ns
    result["relative_p99_limit_ns"] = relative_limit_ns
    result["relative_p99_factor"] = relative_factor
    return result, source
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

function gate2_environment()
    cpu = first(Sys.cpu_info())
    git_status = command_output(`git status --porcelain=v1`)
    manifest_path = joinpath(@__DIR__, "Manifest.toml")
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
        "manifest_sha256" => isfile(manifest_path) ?
            bytes2hex(SHA.sha256(read(manifest_path))) : "missing",
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
        "julia_depot_path" => get(ENV, "JULIA_DEPOT_PATH", "default"),
        "allowed_cpus" => allowed_cpu_list(),
        "scaling_governor_cpu0" => optional_file(
            "/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor"),
        "command" => "julia --threads=1 --project=benchmarks " *
            "benchmarks/benchmark_gate2_serial_plant.jl",
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

function write_artifact(path::AbstractString, artifact)
    isempty(path) && return nothing
    buffer = IOBuffer()
    TOML.print(buffer, artifact; sorted=true)
    bytes = take!(buffer)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, bytes)
    end
    println("wrote Gate 2 serial-plant artifact: ", path)
    return nothing
end

function run_gate2_serial_plant_benchmark()
    configure_gate2_benchmark!()
    contract = TOML.parsefile(GATE2_CONTRACT_PATH)
    workload = contract["workload"]
    samples = parse(Int, get(ENV, "AOS_GATE2_SAMPLES",
        string(contract["samples_per_run"])))
    runs_count = parse(Int, get(ENV, "AOS_GATE2_RUNS",
        string(contract["runs"])))
    warmup = parse(Int, get(ENV, "AOS_GATE2_WARMUP",
        string(contract["warmup_operations"])))
    samples > 0 || error("AOS_GATE2_SAMPLES must be > 0")
    runs_count > 0 || error("AOS_GATE2_RUNS must be > 0")
    warmup >= 0 || error("AOS_GATE2_WARMUP must be >= 0")
    minimum_p99_samples = Int(contract["minimum_samples_for_p99"])
    p99_supported = samples >= minimum_p99_samples
    lowest_ns = Int64(contract["histogram_lowest_ns"])
    highest_ns = Int64(contract["histogram_highest_ns"])
    significant_figures = Int(contract["histogram_significant_figures"])
    relative_factor = Float64(contract["relative_p99_factor"])
    max_alloc_bytes = Int64(contract["max_alloc_bytes"])

    println("gate2_serial_plant_contract")
    println("  load_model: warmed serial self-paced service time")
    println("  samples_per_run: ", samples)
    println("  runs: ", runs_count)
    println("  warmup_operations: ", warmup)
    println("  p99_claim_supported: ", p99_supported)
    println("  baseline_artifact: ", isempty(GATE2_BASELINE_PATH) ?
        "none" : GATE2_BASELINE_PATH)

    operation = prepare_serial_plant_operation(workload)
    GC.gc()
    first_start = time_ns()
    operation()
    first_use_ns = Int64(time_ns() - first_start)
    correctness = validate_replay_and_reordering(workload)
    for _ in 1:warmup
        operation()
    end
    GC.gc()
    steady_alloc_bytes = Int64(@allocated operation())
    allocation_gate_passed = steady_alloc_bytes <= max_alloc_bytes
    println("  steady_alloc_bytes: ", steady_alloc_bytes)

    runs = Vector{Dict{String,Any}}(undef, runs_count)
    for run_index in 1:runs_count
        runs[run_index] = measure_run!(operation, samples, lowest_ns,
            highest_ns, significant_figures)
        println("  run=", run_index,
            " p50_ns=", runs[run_index]["p50_ns"],
            " p99_ns=", runs[run_index]["p99_ns"])
    end
    summary = summarize_runs(runs)
    relative_regression, baseline_source = load_baseline(
        GATE2_BASELINE_PATH, String(contract["name"]), contract["contract"],
        workload, relative_factor,
        Int64(summary["median_p99_ns"]), p99_supported)
    regression = merge(relative_regression, Dict{String,Any}(
        "allocation_gate_passed" => allocation_gate_passed,
        "allocation_limit_bytes" => max_alloc_bytes,
        "allocation_observed_bytes" => steady_alloc_bytes,
        "absolute_latency_gate_evaluated" => false,
        "absolute_latency_gate_passed" => true,
    ))
    all_gates_passed = allocation_gate_passed &&
        regression["relative_p99_gate_passed"]

    environment = gate2_environment()
    if !isempty(GATE2_OUTPUT_PATH) && environment["git_dirty"]
        error("refusing to write durable Gate 2 evidence from a dirty worktree")
    end
    artifact = Dict{String,Any}(
        "schema_version" => Int(contract["schema_version"]),
        "benchmark" => String(contract["name"]),
        "evidence_class" => String(contract["evidence_class"]),
        "source_contract" => relpath(abspath(GATE2_CONTRACT_PATH), @__DIR__),
        "characterized_source_revision" => environment["git_commit"],
        "configured_samples_per_run" => samples,
        "configured_runs" => runs_count,
        "configured_warmup_operations" => warmup,
        "minimum_samples_for_p99" => minimum_p99_samples,
        "p99_claim_supported" => p99_supported,
        "contract" => contract["contract"],
        "scope_exclusions" => contract["scope_exclusions"],
        "workload" => workload,
        "histogram" => Dict(
            "lowest_ns" => lowest_ns,
            "highest_ns" => highest_ns,
            "significant_figures" => significant_figures,
        ),
        "timer_recording_overhead" => timer_recording_overhead(
            min(samples, 100_000), lowest_ns, highest_ns,
            significant_figures),
        "environment" => environment,
        "correctness" => correctness,
        "first_use_ns" => first_use_ns,
        "steady_alloc_bytes" => steady_alloc_bytes,
        "runs" => runs,
        "summary" => summary,
        "regression" => regression,
        "baseline_artifact" => baseline_source,
        "final_observations" => final_observation_summary(operation),
        "all_gates_passed" => all_gates_passed,
    )
    write_artifact(GATE2_OUTPUT_PATH, artifact)
    all_gates_passed || error(
        "one or more Gate 2 serial-plant benchmark gates failed")
    return artifact
end

run_gate2_serial_plant_benchmark()

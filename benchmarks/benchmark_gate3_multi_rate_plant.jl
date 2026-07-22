using AdaptiveOpticsSim
using Base64
using Dates
using HdrHistogram
using LinearAlgebra
using SHA
using Statistics
using TOML

include(joinpath(@__DIR__, "support", "gate3_multi_rate_plant.jl"))

const AOS = AdaptiveOpticsSim
const Gate3Plant = Gate3MultiRatePlantBenchmark
const GATE3_MULTI_RATE_CONTRACT_PATH = get(ENV,
    "AOS_GATE3_MULTI_RATE_CONTRACT",
    joinpath(@__DIR__, "contracts", "gate3_multi_rate_plant.toml"))
const GATE3_MULTI_RATE_OUTPUT_PATH = get(ENV,
    "AOS_GATE3_MULTI_RATE_OUTPUT", "")
const GATE3_SCHEDULER_MANIFEST_PATH = get(ENV,
    "AOS_GATE3_SCHEDULER_MANIFEST",
    joinpath(@__DIR__, "results", "gate3", "manifest.toml"))

function configure_gate3_multi_rate_benchmark!()
    Threads.nthreads() == 1 || error(
        "Gate 3 multi-rate evidence requires one Julia thread")
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
    wall_start_ns::UInt64, wall_end_ns::UInt64, samples::Int, gc_counters,
    simulated_start_ns::Int64, simulated_end_ns::Int64,
    product_sequence_delta::Vector{UInt64})
    encoded = HdrHistogram.encode_into_compressed_byte_buffer(histogram)
    decoded = HdrHistogram.decode_from_compressed_byte_buffer(encoded)
    HdrHistogram.total_count(decoded) == samples || error(
        "encoded Gate 3 multi-rate histogram lost samples")
    HdrHistogram.value_at_percentile(decoded, 99.0) ==
        HdrHistogram.value_at_percentile(histogram, 99.0) || error(
        "encoded Gate 3 multi-rate histogram changed p99")
    wall_ns = Int64(wall_end_ns - wall_start_ns)
    return Dict{String,Any}(
        "samples" => samples,
        "completed_operations" => samples,
        "failed_operations" => 0,
        "monotonic_start_ns" => wall_start_ns,
        "monotonic_end_ns" => wall_end_ns,
        "wall_ns" => wall_ns,
        "throughput_hz" => 1.0e9 * samples / wall_ns,
        "simulated_start_ns" => simulated_start_ns,
        "simulated_end_ns" => simulated_end_ns,
        "simulated_span_ns" => simulated_end_ns - simulated_start_ns,
        "product_sequence_delta" => Int64.(product_sequence_delta),
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

function measure_run!(operation::Gate3Plant.MultiRatePlantOperation,
    samples::Int, lowest_ns::Int64, highest_ns::Int64,
    significant_figures::Int)
    histogram = HdrHistogram.Histogram(lowest_ns, highest_ns,
        significant_figures)
    sequences_before = Gate3Plant.product_sequence_vector(operation)
    simulated_start = AOS.scheduler_timestamp(operation.state.scheduler)
    simulated_end = simulated_start
    GC.gc()
    gc_before = Base.gc_num()
    wall_start = time_ns()
    @inbounds for _ in 1:samples
        start = time_ns()
        simulated_end = operation()
        elapsed = Int64(time_ns() - start)
        HdrHistogram.record_value!(histogram, elapsed)
    end
    wall_end = time_ns()
    gc_after = Base.gc_num()
    sequences_after = Gate3Plant.product_sequence_vector(operation)
    return histogram_summary(histogram, wall_start, wall_end, samples,
        gc_delta(gc_before, gc_after),
        AOS.plant_nanoseconds(simulated_start),
        AOS.plant_nanoseconds(simulated_end),
        sequences_after .- sequences_before)
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

function gate3_multi_rate_environment()
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
            "benchmarks/benchmark_gate3_multi_rate_plant.jl",
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

function validate_scheduler_evidence(path::AbstractString)
    manifest = TOML.parsefile(path)
    entries = manifest["artifacts"]
    isempty(entries) && error("Gate 3 scheduler evidence manifest is empty")
    entry = last(entries)
    artifact_path = normpath(joinpath(dirname(path), entry["path"]))
    expected_sha = String(entry["sha256"])
    observed_sha = bytes2hex(SHA.sha256(read(artifact_path)))
    observed_sha == expected_sha || error(
        "Gate 3 scheduler artifact does not match its manifest digest")
    artifact = TOML.parsefile(artifact_path)
    artifact["all_gates_passed"] || error(
        "Gate 3 scheduler artifact has a failed gate")
    results = artifact["generator_results"]
    counts = Int[result["generator_count"] for result in results]
    counts == [1, 8, 32, 128, 256] || error(
        "Gate 3 scheduler artifact does not cover the maintained generator matrix")
    all(result -> result["allocation_gate_passed"] &&
            result["sequence_valid"] &&
            result["simultaneous_order_valid"], results) || error(
        "Gate 3 scheduler artifact has incomplete correctness evidence")
    return Dict{String,Any}(
        "manifest" => relpath(abspath(path), @__DIR__),
        "artifact" => relpath(abspath(artifact_path), @__DIR__),
        "artifact_sha256" => observed_sha,
        "artifact_id" => String(entry["id"]),
        "characterized_source_revision" => String(
            artifact["characterized_source_revision"]),
        "generator_counts" => counts,
        "all_gates_passed" => true,
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
    println("wrote Gate 3 multi-rate plant artifact: ", path)
    return nothing
end

function run_gate3_multi_rate_plant_benchmark()
    configure_gate3_multi_rate_benchmark!()
    contract = TOML.parsefile(GATE3_MULTI_RATE_CONTRACT_PATH)
    workload = contract["workload"]
    samples = parse(Int, get(ENV, "AOS_GATE3_MULTI_RATE_SAMPLES",
        string(contract["samples_per_run"])))
    runs_count = parse(Int, get(ENV, "AOS_GATE3_MULTI_RATE_RUNS",
        string(contract["runs"])))
    warmup = parse(Int, get(ENV, "AOS_GATE3_MULTI_RATE_WARMUP",
        string(contract["warmup_timestamps"])))
    allocation_timestamps = parse(Int, get(ENV,
        "AOS_GATE3_MULTI_RATE_ALLOCATION_TIMESTAMPS",
        string(contract["allocation_timestamps"])))
    samples > 0 || error("AOS_GATE3_MULTI_RATE_SAMPLES must be > 0")
    runs_count > 0 || error("AOS_GATE3_MULTI_RATE_RUNS must be > 0")
    warmup >= 0 || error("AOS_GATE3_MULTI_RATE_WARMUP must be >= 0")
    allocation_timestamps > 0 || error(
        "AOS_GATE3_MULTI_RATE_ALLOCATION_TIMESTAMPS must be > 0")
    minimum_p99_samples = Int(contract["minimum_samples_for_p99"])
    p99_supported = samples >= minimum_p99_samples
    lowest_ns = Int64(contract["histogram_lowest_ns"])
    highest_ns = Int64(contract["histogram_highest_ns"])
    significant_figures = Int(contract["histogram_significant_figures"])
    max_alloc_per_timestamp = Int64(
        contract["max_alloc_bytes_per_timestamp"])

    println("gate3_multi_rate_plant_contract")
    println("  load_model: warmed serial self-paced timestamp service time")
    println("  samples_per_run: ", samples)
    println("  runs: ", runs_count)
    println("  warmup_timestamps: ", warmup)
    println("  p99_claim_supported: ", p99_supported)

    _, first_operation = Gate3Plant.prepare_multi_rate_operation(workload)
    GC.gc()
    first_start = time_ns()
    first_operation()
    first_use_ns = Int64(time_ns() - first_start)

    correctness = Dict{String,Any}(
        "representative_replay" =>
            Gate3Plant.validate_replay_and_reordering(workload),
        "faulted_trigger_fanout" =>
            Gate3Plant.validate_faulted_trigger_fanout(workload),
        "fixed_storage" => Gate3Plant.validate_fixed_storage(workload,
            Int64(contract["long_run_horizon_ns"]),
            Int(contract["max_long_run_timestamps"])),
        "direct_jump" => Gate3Plant.validate_direct_jump(workload,
            Int(contract["direct_jump_period_multiplier"])),
        "scheduler_scaling" => validate_scheduler_evidence(
            GATE3_SCHEDULER_MANIFEST_PATH),
    )

    plant, operation = Gate3Plant.prepare_multi_rate_operation(workload)
    Gate3Plant.run_timestamp_window!(operation, warmup)
    Gate3Plant.run_timestamp_window!(operation, allocation_timestamps)
    GC.gc()
    allocated_bytes = Int64(@allocated Gate3Plant.run_timestamp_window!(
        operation, allocation_timestamps))
    alloc_bytes_per_timestamp = allocated_bytes / allocation_timestamps
    allocation_limit_bytes = max_alloc_per_timestamp *
        allocation_timestamps
    allocation_gate_passed = allocated_bytes <= allocation_limit_bytes
    println("  allocation_window_bytes: ", allocated_bytes)
    println("  allocation_bytes_per_timestamp: ",
        alloc_bytes_per_timestamp)

    runs = Vector{Dict{String,Any}}(undef, runs_count)
    for run_index in 1:runs_count
        runs[run_index] = measure_run!(operation, samples, lowest_ns,
            highest_ns, significant_figures)
        println("  run=", run_index,
            " p50_ns=", runs[run_index]["p50_ns"],
            " p99_ns=", runs[run_index]["p99_ns"])
    end
    summary = summarize_runs(runs)
    regression = Dict{String,Any}(
        "allocation_gate_passed" => allocation_gate_passed,
        "allocation_window_timestamps" => allocation_timestamps,
        "allocation_limit_bytes_per_timestamp" =>
            max_alloc_per_timestamp,
        "allocation_limit_window_bytes" => allocation_limit_bytes,
        "allocation_observed_window_bytes" => allocated_bytes,
        "allocation_observed_bytes_per_timestamp" =>
            alloc_bytes_per_timestamp,
        "absolute_latency_gate_evaluated" => false,
        "absolute_latency_gate_passed" => true,
        "relative_latency_gate_evaluated" => false,
        "relative_latency_gate_passed" => true,
    )
    all_gates_passed = allocation_gate_passed

    environment = gate3_multi_rate_environment()
    if !isempty(GATE3_MULTI_RATE_OUTPUT_PATH) && environment["git_dirty"]
        error("refusing to write durable Gate 3 evidence from a dirty worktree")
    end
    artifact = Dict{String,Any}(
        "schema_version" => Int(contract["schema_version"]),
        "benchmark" => String(contract["name"]),
        "evidence_class" => String(contract["evidence_class"]),
        "source_contract" => relpath(abspath(
            GATE3_MULTI_RATE_CONTRACT_PATH), @__DIR__),
        "characterized_source_revision" => environment["git_commit"],
        "configured_samples_per_run" => samples,
        "configured_runs" => runs_count,
        "configured_warmup_timestamps" => warmup,
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
            max(samples, 10_000), lowest_ns, highest_ns,
            significant_figures),
        "environment" => environment,
        "correctness" => correctness,
        "first_use_ns" => first_use_ns,
        "regression" => regression,
        "runs" => runs,
        "summary" => summary,
        "final_products" => Gate3Plant.snapshot_dict(
            Gate3Plant.product_snapshot(plant, operation)),
        "all_gates_passed" => all_gates_passed,
    )
    all_gates_passed || error(
        "Gate 3 multi-rate plant allocation gate failed")
    write_artifact(GATE3_MULTI_RATE_OUTPUT_PATH, artifact)
    return artifact
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_gate3_multi_rate_plant_benchmark()
end

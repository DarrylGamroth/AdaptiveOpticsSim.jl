using AdaptiveOpticsSim
using AdaptiveOpticsSim.Plant
using Dates
using HdrHistogram
using LinearAlgebra
using SHA
using Statistics
using TOML

include(joinpath(@__DIR__, "support", "gate4_command_plant.jl"))
include(joinpath(@__DIR__, "support", "hdr_histogram_artifact.jl"))

const AOSPlant = AdaptiveOpticsSim.Plant
const Gate4Plant = Gate4CommandPlantBenchmark
const HdrArtifact = HdrHistogramArtifact
const GATE4_COMMAND_CONTRACT_PATH = get(ENV,
    "AOS_GATE4_COMMAND_CONTRACT",
    joinpath(@__DIR__, "contracts", "gate4_command_plant.toml"))
const GATE4_COMMAND_OUTPUT_PATH = get(ENV,
    "AOS_GATE4_COMMAND_OUTPUT", "")

function configure_gate4_command_benchmark!()
    Threads.nthreads() == 1 ||
        error("Gate 4 command evidence requires one Julia thread")
    BLAS.set_num_threads(1)
    AdaptiveOpticsSim.set_fft_provider_threads!(1)
    return nothing
end

function require_durable_contract_configuration!(
    output_path::AbstractString,
    contract::AbstractDict,
    samples::Int,
    runs::Int,
    warmup_cycles::Int,
    allocation_cycles::Int,
    endpoint_allocation_cycles::Int,
    replay_cycles::Int,
    long_run_cycles::Int,
)
    isempty(output_path) && return nothing
    observed = (
        samples=samples,
        runs=runs,
        warmup_cycles=warmup_cycles,
        allocation_cycles=allocation_cycles,
        endpoint_allocation_cycles=endpoint_allocation_cycles,
        replay_cycles=replay_cycles,
        long_run_cycles=long_run_cycles,
    )
    expected = (
        samples=Int(contract["samples_per_run"]),
        runs=Int(contract["runs"]),
        warmup_cycles=Int(contract["warmup_cycles"]),
        allocation_cycles=Int(contract["allocation_cycles"]),
        endpoint_allocation_cycles=Int(
            contract["endpoint_allocation_cycles"]),
        replay_cycles=Int(contract["replay_cycles"]),
        long_run_cycles=Int(contract["long_run_cycles"]),
    )
    observed == expected ||
        error("refusing to write durable Gate 4 evidence with " *
              "quick-run configuration overrides")
    samples >= Int(contract["minimum_samples_for_p99"]) ||
        error("durable Gate 4 evidence does not support its p99 claim")
    return nothing
end

function gc_delta(before, after)
    return Dict(String(name) => Int64(
        getfield(after, name) - getfield(before, name))
        for name in fieldnames(typeof(before)))
end

function command_histogram_summary(
    histogram::HdrHistogram.Histogram,
    wall_start_ns::UInt64,
    wall_end_ns::UInt64,
    samples::Int,
    gc_counters,
    simulated_start_ns::Int64,
    simulated_end_ns::Int64,
    presented_commands::UInt64,
    terminal_dispositions::UInt64,
    applied_commands::UInt64,
    lowest_ns::Int64,
    highest_ns::Int64,
    significant_figures::Int,
)
    presented_commands == terminal_dispositions == applied_commands ||
        error("Gate 4 measured run lost exact terminal accounting")
    encoded = HdrArtifact.verified_sparse_histogram(histogram, lowest_ns,
        highest_ns, significant_figures, samples)
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
        "presented_commands" => Int64(presented_commands),
        "terminal_dispositions" => Int64(terminal_dispositions),
        "applied_commands" => Int64(applied_commands),
        "exact_terminal_accounting" => true,
        "min_ns" => min(histogram),
        "mean_ns" => HdrHistogram.mean(histogram),
        "p50_ns" => HdrHistogram.value_at_percentile(histogram, 50.0),
        "p90_ns" => HdrHistogram.value_at_percentile(histogram, 90.0),
        "p99_ns" => HdrHistogram.value_at_percentile(histogram, 99.0),
        "p99_9_ns" =>
            HdrHistogram.value_at_percentile(histogram, 99.9),
        "max_ns" => max(histogram),
        "gc" => gc_counters,
        encoded...,
    )
end

function measure_command_run!(
    operation::Gate4Plant.CommandPlantOperation,
    samples::Int,
    lowest_ns::Int64,
    highest_ns::Int64,
    significant_figures::Int,
)
    histogram = HdrHistogram.Histogram(
        lowest_ns, highest_ns, significant_figures)
    presented_before = operation.presented
    terminal_before = operation.terminal
    applied_before = operation.applied
    simulated_start = AOSPlant.scheduler_timestamp(
        operation.state.scheduler)
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
    return command_histogram_summary(
        histogram,
        wall_start,
        wall_end,
        samples,
        gc_delta(gc_before, gc_after),
        AOSPlant.plant_nanoseconds(simulated_start),
        AOSPlant.plant_nanoseconds(simulated_end),
        operation.presented - presented_before,
        operation.terminal - terminal_before,
        operation.applied - applied_before,
        lowest_ns,
        highest_ns,
        significant_figures,
    )
end

median_integer(values) = round(Int64, median(collect(values)))

function summarize_runs(runs::Vector{Dict{String,Any}})
    return Dict{String,Any}(
        "median_p50_ns" => median_integer(run["p50_ns"] for run in runs),
        "median_p90_ns" => median_integer(run["p90_ns"] for run in runs),
        "median_p99_ns" => median_integer(run["p99_ns"] for run in runs),
        "median_p99_9_ns" =>
            median_integer(run["p99_9_ns"] for run in runs),
        "worst_p99_ns" => maximum(run["p99_ns"] for run in runs),
        "worst_p99_9_ns" =>
            maximum(run["p99_9_ns"] for run in runs),
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

function gate4_command_environment()
    cpu = first(Sys.cpu_info())
    git_status = command_output(`git status --porcelain=v1`)
    manifest_path = joinpath(@__DIR__, "Manifest.toml")
    return Dict{String,Any}(
        "timestamp_utc" => string(Dates.now(Dates.UTC)),
        "git_commit" => command_output(`git rev-parse HEAD`),
        "git_dirty" => !isempty(git_status) && git_status != "unknown",
        "git_status_porcelain" => git_status,
        "julia_version" => string(VERSION),
        "adaptive_optics_sim_version" =>
            string(Base.pkgversion(AdaptiveOpticsSim)),
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
        "julia_cpu_target_env" =>
            get(ENV, "JULIA_CPU_TARGET", "default"),
        "julia_depot_path" =>
            get(ENV, "JULIA_DEPOT_PATH", "default"),
        "allowed_cpus" => allowed_cpu_list(),
        "scaling_governor_cpu0" => optional_file(
            "/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor"),
        "command" => "julia --threads=1 --project=benchmarks " *
            "benchmarks/benchmark_gate4_command_plant.jl",
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
    println("wrote Gate 4 command-plant artifact: ", path)
    return nothing
end

function run_gate4_command_plant_benchmark()
    configure_gate4_command_benchmark!()
    contract = TOML.parsefile(GATE4_COMMAND_CONTRACT_PATH)
    workload = contract["workload"]
    samples = parse(Int, get(ENV, "AOS_GATE4_COMMAND_SAMPLES",
        string(contract["samples_per_run"])))
    runs_count = parse(Int, get(ENV, "AOS_GATE4_COMMAND_RUNS",
        string(contract["runs"])))
    warmup = parse(Int, get(ENV, "AOS_GATE4_COMMAND_WARMUP",
        string(contract["warmup_cycles"])))
    allocation_cycles = parse(Int, get(ENV,
        "AOS_GATE4_COMMAND_ALLOCATION_CYCLES",
        string(contract["allocation_cycles"])))
    endpoint_allocation_cycles = parse(Int, get(ENV,
        "AOS_GATE4_COMMAND_ENDPOINT_ALLOCATION_CYCLES",
        string(contract["endpoint_allocation_cycles"])))
    replay_cycles = parse(Int, get(ENV,
        "AOS_GATE4_COMMAND_REPLAY_CYCLES",
        string(contract["replay_cycles"])))
    long_run_cycles = parse(Int, get(ENV,
        "AOS_GATE4_COMMAND_LONG_RUN_CYCLES",
        string(contract["long_run_cycles"])))
    samples > 0 || error("AOS_GATE4_COMMAND_SAMPLES must be > 0")
    runs_count > 0 || error("AOS_GATE4_COMMAND_RUNS must be > 0")
    warmup >= 0 || error("AOS_GATE4_COMMAND_WARMUP must be >= 0")
    allocation_cycles > 0 ||
        error("AOS_GATE4_COMMAND_ALLOCATION_CYCLES must be > 0")
    endpoint_allocation_cycles > 0 ||
        error("AOS_GATE4_COMMAND_ENDPOINT_ALLOCATION_CYCLES must be > 0")
    replay_cycles > 0 ||
        error("AOS_GATE4_COMMAND_REPLAY_CYCLES must be > 0")
    long_run_cycles > 0 ||
        error("AOS_GATE4_COMMAND_LONG_RUN_CYCLES must be > 0")
    require_durable_contract_configuration!(
        GATE4_COMMAND_OUTPUT_PATH,
        contract,
        samples,
        runs_count,
        warmup,
        allocation_cycles,
        endpoint_allocation_cycles,
        replay_cycles,
        long_run_cycles,
    )
    minimum_p99_samples = Int(contract["minimum_samples_for_p99"])
    p99_supported = samples >= minimum_p99_samples
    lowest_ns = Int64(contract["histogram_lowest_ns"])
    highest_ns = Int64(contract["histogram_highest_ns"])
    significant_figures = Int(contract["histogram_significant_figures"])
    max_alloc_per_cycle = Int64(contract["max_alloc_bytes_per_cycle"])

    println("gate4_command_plant_contract")
    println("  load_model: warmed serial self-paced command-cycle service time")
    println("  samples_per_run: ", samples)
    println("  runs: ", runs_count)
    println("  warmup_cycles: ", warmup)
    println("  p99_claim_supported: ", p99_supported)

    _, first_operation = Gate4Plant.prepare_command_operation(workload)
    GC.gc()
    first_start = time_ns()
    first_operation()
    first_use_ns = Int64(time_ns() - first_start)

    correctness = Dict{String,Any}(
        "independent_exposure" =>
            Gate4Plant.validate_independent_exposure_and_reordering(workload),
        "atomic_application" =>
            Gate4Plant.validate_atomic_application(workload),
        "zero_allocation_endpoint" =>
            Gate4Plant.validate_zero_allocation_endpoint(
                endpoint_allocation_cycles),
        "operation_replay" =>
            Gate4Plant.validate_operation_replay(workload, replay_cycles),
        "fixed_storage" => Gate4Plant.validate_fixed_storage(
            workload, warmup, long_run_cycles),
    )

    _, operation = Gate4Plant.prepare_command_operation(workload)
    Gate4Plant.run_command_cycles!(operation, warmup)
    Gate4Plant.run_command_cycles!(operation, allocation_cycles)
    GC.gc()
    allocated_bytes = Int64(@allocated Gate4Plant.run_command_cycles!(
        operation, allocation_cycles))
    alloc_bytes_per_cycle = allocated_bytes / allocation_cycles
    allocation_limit_bytes = max_alloc_per_cycle * allocation_cycles
    allocation_gate_passed = allocated_bytes <= allocation_limit_bytes
    println("  allocation_window_bytes: ", allocated_bytes)
    println("  allocation_bytes_per_cycle: ", alloc_bytes_per_cycle)

    runs = Vector{Dict{String,Any}}(undef, runs_count)
    for run_index in 1:runs_count
        runs[run_index] = measure_command_run!(operation, samples,
            lowest_ns, highest_ns, significant_figures)
        println("  run=", run_index,
            " p50_ns=", runs[run_index]["p50_ns"],
            " p99_ns=", runs[run_index]["p99_ns"])
    end
    summary = summarize_runs(runs)
    regression = Dict{String,Any}(
        "allocation_gate_passed" => allocation_gate_passed,
        "allocation_window_cycles" => allocation_cycles,
        "allocation_limit_bytes_per_cycle" => max_alloc_per_cycle,
        "allocation_limit_window_bytes" => allocation_limit_bytes,
        "allocation_observed_window_bytes" => allocated_bytes,
        "allocation_observed_bytes_per_cycle" => alloc_bytes_per_cycle,
        "absolute_latency_gate_evaluated" => false,
        "absolute_latency_gate_passed" => true,
        "relative_latency_gate_evaluated" => false,
        "relative_latency_gate_passed" => true,
    )
    all_gates_passed = allocation_gate_passed

    environment = gate4_command_environment()
    if !isempty(GATE4_COMMAND_OUTPUT_PATH) && environment["git_dirty"]
        error("refusing to write durable Gate 4 evidence from a dirty worktree")
    end
    artifact = Dict{String,Any}(
        "schema_version" => Int(contract["schema_version"]),
        "benchmark" => String(contract["name"]),
        "evidence_class" => String(contract["evidence_class"]),
        "source_contract" => relpath(abspath(
            GATE4_COMMAND_CONTRACT_PATH), @__DIR__),
        "characterized_source_revision" => environment["git_commit"],
        "configured_samples_per_run" => samples,
        "configured_runs" => runs_count,
        "configured_warmup_cycles" => warmup,
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
        "final_accounting" => Dict{String,Any}(
            "processed_cycles" => Int64(operation.cycles),
            "presented_commands" => Int64(operation.presented),
            "terminal_dispositions" => Int64(operation.terminal),
            "applied_commands" => Int64(operation.applied),
            "maximum_pending_commands" =>
                Int64(operation.maximum_pending),
            "maximum_terminal_dispositions" =>
                Int64(operation.maximum_dispositions),
            "exact_terminal_accounting" =>
                operation.presented == operation.terminal ==
                operation.applied,
        ),
        "all_gates_passed" => all_gates_passed,
    )
    all_gates_passed ||
        error("Gate 4 command-plant allocation gate failed")
    write_artifact(GATE4_COMMAND_OUTPUT_PATH, artifact)
    return artifact
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_gate4_command_plant_benchmark()
end

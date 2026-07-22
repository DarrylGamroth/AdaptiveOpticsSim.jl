using AdaptiveOpticsSim
using Base64
using Dates
using HdrHistogram
using LinearAlgebra
using SHA
using Statistics
using TOML

include(joinpath(@__DIR__, "support", "hdr_histogram_artifact.jl"))

const AOS = AdaptiveOpticsSim
const HdrArtifact = HdrHistogramArtifact
const GATE3_SCHEDULER_CONTRACT_PATH = get(ENV,
    "AOS_GATE3_SCHEDULER_CONTRACT",
    joinpath(@__DIR__, "contracts", "gate3_event_scheduler.toml"))
const GATE3_SCHEDULER_OUTPUT_PATH = get(ENV,
    "AOS_GATE3_SCHEDULER_OUTPUT", "")

struct SchedulerBenchmarkOperation
    scheduler::AOS.PreparedEventScheduler
    state::AOS.EventSchedulerState
    workspace::AOS.EventSchedulerWorkspace
    recurrence::PlantDuration
end

@inline function (operation::SchedulerBenchmarkOperation)()
    claim = AOS.claim_next_event!(operation.workspace,
        operation.scheduler, operation.state)
    claim === nothing && error("scheduler benchmark exhausted active events")
    next_timestamp = AOS.claimed_event_key(claim).timestamp +
        operation.recurrence
    return AOS.reschedule_event!(operation.scheduler, operation.state,
        claim, next_timestamp)
end

function configure_scheduler_benchmark!()
    Threads.nthreads() == 1 || error(
        "Gate 3 scheduler evidence requires one Julia thread")
    BLAS.set_num_threads(1)
    AdaptiveOpticsSim.set_fft_provider_threads!(1)
    return nothing
end

function prepare_scheduler_operation(generator_count::Int)
    generator_count > 0 || error("generator count must be positive")
    definitions = Vector{AOS.EventGeneratorDefinition}(
        undef, generator_count)
    for index in 1:generator_count
        definitions[index] = AOS.EventGeneratorDefinition(
            PlantTimestamp(index - 1), AOS.OpticalSamplePhase, index)
    end
    scheduler = AOS.prepare_event_scheduler(definitions;
        capacity=generator_count)
    state = AOS.EventSchedulerState(scheduler)
    workspace = AOS.EventSchedulerWorkspace(scheduler)
    return SchedulerBenchmarkOperation(scheduler, state, workspace,
        PlantDuration(generator_count))
end

function validate_scheduler_sequence(generator_count::Int)
    operation = prepare_scheduler_operation(generator_count)
    event_count = 4 * generator_count
    for event_index in 1:event_count
        claim = AOS.claim_next_event!(operation.workspace,
            operation.scheduler, operation.state)
        claim === nothing && error("scheduler sequence ended early")
        key = AOS.claimed_event_key(claim)
        expected_timestamp = PlantTimestamp(event_index - 1)
        expected_ordinal = UInt32(mod(event_index - 1, generator_count) + 1)
        expected_occurrence = UInt64(fld(event_index - 1,
            generator_count) + 1)
        key.timestamp == expected_timestamp || error(
            "scheduler timestamp sequence mismatch")
        key.ordinal == expected_ordinal || error(
            "scheduler ordinal sequence mismatch")
        key.occurrence == expected_occurrence || error(
            "scheduler occurrence sequence mismatch")
        AOS.reschedule_event!(operation.scheduler, operation.state, claim,
            key.timestamp + operation.recurrence)
    end
    return true
end

function validate_simultaneous_order(generator_count::Int)
    definitions = Vector{AOS.EventGeneratorDefinition}(
        undef, generator_count)
    for index in 1:generator_count
        ordinal = generator_count - index + 1
        definitions[index] = AOS.EventGeneratorDefinition(
            PlantTimestamp(0), AOS.OpticalSamplePhase, ordinal)
    end
    scheduler = AOS.prepare_event_scheduler(definitions;
        capacity=generator_count)
    state = AOS.EventSchedulerState(scheduler)
    workspace = AOS.EventSchedulerWorkspace(scheduler)
    AOS.scan_due_events!(workspace, scheduler, state) == generator_count ||
        error("simultaneous due count mismatch")
    for index in 1:generator_count
        key = AOS.due_event_key(workspace, scheduler, state, index)
        key.timestamp == PlantTimestamp(0) || error(
            "simultaneous timestamp mismatch")
        key.ordinal == UInt32(index) || error(
            "simultaneous event ordering depends on declaration order")
    end
    return true
end

function gc_delta(before, after)
    return Dict(String(name) => Int64(
        getfield(after, name) - getfield(before, name))
        for name in fieldnames(typeof(before)))
end

function histogram_summary(histogram::HdrHistogram.Histogram,
    wall_start_ns::UInt64, wall_end_ns::UInt64, samples::Int, gc_counters,
    lowest_ns::Int64, highest_ns::Int64, significant_figures::Int)
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
        "min_ns" => min(histogram),
        "mean_ns" => HdrHistogram.mean(histogram),
        "p50_ns" => HdrHistogram.value_at_percentile(histogram, 50.0),
        "p90_ns" => HdrHistogram.value_at_percentile(histogram, 90.0),
        "p99_ns" => HdrHistogram.value_at_percentile(histogram, 99.0),
        "p99_9_ns" => HdrHistogram.value_at_percentile(histogram, 99.9),
        "max_ns" => max(histogram),
        "gc" => gc_counters,
        encoded...,
    )
end

function measure_run!(operation::SchedulerBenchmarkOperation, samples::Int,
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
        gc_delta(gc_before, gc_after), lowest_ns, highest_ns,
        significant_figures)
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

function scheduler_environment()
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
            "benchmarks/benchmark_gate3_event_scheduler.jl",
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
    println("wrote Gate 3 event-scheduler artifact: ", path)
    return nothing
end

function run_gate3_event_scheduler_benchmark()
    configure_scheduler_benchmark!()
    contract = TOML.parsefile(GATE3_SCHEDULER_CONTRACT_PATH)
    generator_counts = Int.(contract["generator_counts"])
    samples = parse(Int, get(ENV, "AOS_GATE3_SCHEDULER_SAMPLES",
        string(contract["samples_per_run"])))
    runs_count = parse(Int, get(ENV, "AOS_GATE3_SCHEDULER_RUNS",
        string(contract["runs"])))
    warmup = parse(Int, get(ENV, "AOS_GATE3_SCHEDULER_WARMUP",
        string(contract["warmup_operations"])))
    samples > 0 || error("AOS_GATE3_SCHEDULER_SAMPLES must be > 0")
    runs_count > 0 || error("AOS_GATE3_SCHEDULER_RUNS must be > 0")
    warmup >= 0 || error("AOS_GATE3_SCHEDULER_WARMUP must be >= 0")
    lowest_ns = Int64(contract["histogram_lowest_ns"])
    highest_ns = Int64(contract["histogram_highest_ns"])
    significant_figures = Int(contract["histogram_significant_figures"])
    allocation_limit = Int64(contract["max_alloc_bytes"])
    minimum_p99_samples = Int(contract["minimum_samples_for_p99"])
    p99_supported = samples >= minimum_p99_samples

    first_operation = prepare_scheduler_operation(first(generator_counts))
    GC.gc()
    first_start = time_ns()
    first_operation()
    first_use_ns = Int64(time_ns() - first_start)

    results = Vector{Dict{String,Any}}(undef, length(generator_counts))
    all_gates_passed = true
    for (result_index, generator_count) in enumerate(generator_counts)
        sequence_valid = validate_scheduler_sequence(generator_count)
        simultaneous_valid = validate_simultaneous_order(generator_count)
        operation = prepare_scheduler_operation(generator_count)
        for _ in 1:warmup
            operation()
        end
        operation()
        steady_alloc_bytes = Int64(@allocated operation())
        allocation_gate_passed = steady_alloc_bytes <= allocation_limit
        all_gates_passed &= sequence_valid && simultaneous_valid &&
            allocation_gate_passed
        runs = Vector{Dict{String,Any}}(undef, runs_count)
        for run_index in 1:runs_count
            runs[run_index] = measure_run!(operation, samples,
                lowest_ns, highest_ns, significant_figures)
        end
        storage_bytes = Base.summarysize((operation.scheduler,
            operation.state, operation.workspace))
        results[result_index] = Dict{String,Any}(
            "generator_count" => generator_count,
            "capacity" => AOS.event_scheduler_capacity(operation.scheduler),
            "registry_length" => length(getfield(operation.scheduler,
                :definitions)),
            "cursor_length" => length(operation.state.cursors),
            "due_slot_length" => length(operation.workspace.due_slots),
            "prepared_storage_bytes" => storage_bytes,
            "steady_alloc_bytes" => steady_alloc_bytes,
            "allocation_gate_passed" => allocation_gate_passed,
            "sequence_valid" => sequence_valid,
            "simultaneous_order_valid" => simultaneous_valid,
            "summary" => summarize_runs(runs),
            "runs" => runs,
        )
        println("generator_count=", generator_count,
            " storage_bytes=", storage_bytes,
            " alloc_bytes=", steady_alloc_bytes,
            " median_p99_ns=",
            results[result_index]["summary"]["median_p99_ns"])
    end

    artifact = Dict{String,Any}(
        "schema_version" => 1,
        "benchmark" => contract["name"],
        "evidence_class" =>
            "warmed in-process serial CPU scheduler service cost",
        "characterized_source_revision" => command_output(
            `git rev-parse HEAD`),
        "source_contract" => relpath(GATE3_SCHEDULER_CONTRACT_PATH,
            @__DIR__),
        "configured_samples_per_run" => samples,
        "configured_runs" => runs_count,
        "configured_warmup_operations" => warmup,
        "minimum_samples_for_p99" => minimum_p99_samples,
        "p99_claim_supported" => p99_supported,
        "first_use_ns" => first_use_ns,
        "all_gates_passed" => all_gates_passed,
        "contract" => contract["contract"],
        "workload" => contract["workload"],
        "scope_exclusions" => contract["scope_exclusions"],
        "histogram" => Dict(
            "lowest_ns" => lowest_ns,
            "highest_ns" => highest_ns,
            "significant_figures" => significant_figures,
        ),
        "timer_recording_overhead" => timer_recording_overhead(
            max(samples, 10_000), lowest_ns, highest_ns,
            significant_figures),
        "environment" => scheduler_environment(),
        "generator_results" => results,
    )
    all_gates_passed || error(
        "Gate 3 event-scheduler correctness or allocation gate failed")
    write_artifact(GATE3_SCHEDULER_OUTPUT_PATH, artifact)
    return artifact
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_gate3_event_scheduler_benchmark()
end

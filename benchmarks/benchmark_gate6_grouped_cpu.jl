using AdaptiveOpticsSim
using AdaptiveOpticsSim.Atmospheres
using AdaptiveOpticsSim.Plant
using Base64
using Dates
using HdrHistogram
using LinearAlgebra
using SHA
using Statistics
using TOML

include(joinpath(@__DIR__, "support", "gate3_multi_rate_plant.jl"))
include(joinpath(@__DIR__, "support", "gate6_grouped_cpu.jl"))
include(joinpath(@__DIR__, "support", "hdr_histogram_artifact.jl"))

const AOSPlant = AdaptiveOpticsSim.Plant
const Gate3Plant = Gate3MultiRatePlantBenchmark
const Gate6CPU = Gate6GroupedCPUBenchmark
const HdrArtifact = HdrHistogramArtifact
const GATE6_CPU_CONTRACT_PATH = get(
    ENV,
    "AOS_GATE6_CPU_CONTRACT",
    joinpath(@__DIR__, "contracts", "gate6_grouped_cpu.toml"),
)
const GATE6_CPU_OUTPUT_PATH = get(ENV, "AOS_GATE6_CPU_OUTPUT", "")

function configure_gate6_cpu!(contract)
    expected_threads = Int(contract["julia_threads"])
    Threads.nthreads() == expected_threads || error(
        "Gate 6 grouped CPU evidence requires $expected_threads Julia threads")
    BLAS.set_num_threads(Int(contract["blas_threads"]))
    AdaptiveOpticsSim.Backends.set_fft_provider_threads!(
        Int(contract["fft_threads"]))
    budget = AOSPlant.grouped_cpu_execution_budget(
        cpu_context_count=Int(contract["cpu_contexts"]),
        julia_thread_count=expected_threads,
        outer_owner_count=Int(contract["outer_group_owners"]),
        group_julia_thread_count=Int(contract["group_julia_threads"]),
        fft_thread_count=Int(contract["fft_threads"]),
        blas_thread_count=Int(contract["blas_threads"]),
    )
    environment = AOSPlant.CPUExecutionEnvironment(
        available_cpu_context_count=Int(contract["cpu_contexts"]),
        fft_thread_count=Int(contract["fft_threads"]),
    )
    AOSPlant.validate_cpu_execution_budget(budget, environment)
    return budget
end

function gate6_gc_delta(before, after)
    return Dict(
        String(name) =>
            Int64(getfield(after, name) - getfield(before, name))
        for name in fieldnames(typeof(before))
    )
end

function gate6_histogram_summary(
    histogram::HdrHistogram.Histogram,
    wall_start_ns::UInt64,
    wall_end_ns::UInt64,
    samples::Int,
    gc_counters,
    simulated_start::AOSPlant.PlantTimestamp,
    simulated_end::AOSPlant.PlantTimestamp,
    lowest_ns::Int64,
    highest_ns::Int64,
    significant_figures::Int,
)
    encoded = HdrArtifact.verified_sparse_histogram(
        histogram,
        lowest_ns,
        highest_ns,
        significant_figures,
        samples,
    )
    wall_ns = Int64(wall_end_ns - wall_start_ns)
    return Dict{String,Any}(
        "samples" => samples,
        "completed_operations" => samples,
        "failed_operations" => 0,
        "monotonic_start_ns" => wall_start_ns,
        "monotonic_end_ns" => wall_end_ns,
        "wall_ns" => wall_ns,
        "throughput_hz" => 1.0e9 * samples / wall_ns,
        "simulated_start_ns" => AOSPlant.plant_nanoseconds(
            simulated_start),
        "simulated_end_ns" => AOSPlant.plant_nanoseconds(simulated_end),
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

function operation_storage(operation::Gate3Plant.MultiRatePlantOperation)
    return operation
end

function operation_storage(operation::Gate6CPU.GroupedPlantOperation)
    return operation.serial_storage
end

function measure_gate6_run!(
    operation,
    samples::Int,
    lowest_ns::Int64,
    highest_ns::Int64,
    significant_figures::Int,
)
    histogram = HdrHistogram.Histogram(
        lowest_ns, highest_ns, significant_figures)
    storage = operation_storage(operation)
    simulated_start =
        AOSPlant.scheduler_timestamp(storage.state.scheduler)
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
    return gate6_histogram_summary(
        histogram,
        wall_start,
        wall_end,
        samples,
        gate6_gc_delta(gc_before, gc_after),
        simulated_start,
        simulated_end,
        lowest_ns,
        highest_ns,
        significant_figures,
    )
end

gate6_median_integer(values) = round(Int64, median(collect(values)))

function summarize_gate6_runs(runs::Vector{Dict{String,Any}})
    p50_values = Int64[run["p50_ns"] for run in runs]
    p90_values = Int64[run["p90_ns"] for run in runs]
    p99_values = Int64[run["p99_ns"] for run in runs]
    return Dict{String,Any}(
        "median_p50_ns" => gate6_median_integer(p50_values),
        "minimum_p50_ns" => minimum(p50_values),
        "maximum_p50_ns" => maximum(p50_values),
        "median_p90_ns" => gate6_median_integer(p90_values),
        "minimum_p90_ns" => minimum(p90_values),
        "maximum_p90_ns" => maximum(p90_values),
        "median_p99_ns" => gate6_median_integer(p99_values),
        "minimum_p99_ns" => minimum(p99_values),
        "maximum_p99_ns" => maximum(p99_values),
        "median_p99_9_ns" =>
            gate6_median_integer(run["p99_9_ns"] for run in runs),
        "median_throughput_hz" => median(
            collect(run["throughput_hz"] for run in runs),
        ),
    )
end

function copy_path_rng_states(prepared)
    return map(1:AOSPlant.path_execution_group_count(prepared)) do ordinal
        group = AOSPlant.path_execution_group(prepared, ordinal)
        copy(AOSPlant.rng_stream_state(group.rngs, Val(:provider)))
    end
end

function copy_atmosphere_rng_states(prepared)
    return map(prepared.atmosphere_rng.streams) do stream
        copy(stream.state)
    end
end

function copy_acquisition_rng_states(prepared)
    return map(prepared.acquisitions) do acquisition
        copy(acquisition.rng)
    end
end

function validate_serial_grouped(raw, horizon_ns::Integer)
    serial_plant, serial =
        Gate3Plant.prepare_multi_rate_operation(raw)
    grouped_plant, grouped_storage =
        Gate3Plant.prepare_multi_rate_operation(raw)
    grouped = Gate6CPU.GroupedPlantOperation(grouped_storage)
    try
        horizon = AOSPlant.PlantTimestamp(horizon_ns)
        serial_count = AOSPlant.run_plant_events_until!(
            serial.prepared, serial.state, serial.workspace, horizon)
        grouped_count = AOSPlant.run_plant_events_until!(
            grouped_storage.prepared,
            grouped_storage.state,
            grouped_storage.workspace,
            horizon,
            grouped.executor,
        )
        serial_snapshot = Gate3Plant.product_snapshot(serial_plant, serial)
        grouped_snapshot =
            Gate3Plant.product_snapshot(grouped_plant, grouped_storage)
        serial_count == grouped_count ||
            error("serial and grouped timestamp counts differ")
        serial_snapshot == grouped_snapshot ||
            error("serial and grouped acquisition products differ")
        serial.state.scheduler.cursors ==
            grouped_storage.state.scheduler.cursors ||
            error("serial and grouped scheduler cursors differ")
        serial.state.product_sequences ==
            grouped_storage.state.product_sequences ||
            error("serial and grouped product sequences differ")
        serial.state.product_ready_timestamps ==
            grouped_storage.state.product_ready_timestamps ||
            error("serial and grouped readiness timestamps differ")
        copy_path_rng_states(serial.prepared) ==
            copy_path_rng_states(grouped_storage.prepared) ||
            error("serial and grouped path RNG states differ")
        copy_atmosphere_rng_states(serial.prepared) ==
            copy_atmosphere_rng_states(grouped_storage.prepared) ||
            error("serial and grouped atmosphere RNG states differ")
        copy_acquisition_rng_states(serial.prepared) ==
            copy_acquisition_rng_states(grouped_storage.prepared) ||
            error("serial and grouped acquisition RNG states differ")
        serial_epoch = current_epoch(serial.prepared.atmosphere)
        grouped_epoch = current_epoch(grouped_storage.prepared.atmosphere)
        epoch_time(serial_epoch) == epoch_time(grouped_epoch) ||
            error("serial and grouped atmosphere times differ")
        epoch_sequence(serial_epoch) == epoch_sequence(grouped_epoch) ||
            error("serial and grouped atmosphere sequences differ")
        all(id -> !iszero(id), grouped.executor.task_ids) ||
            error("a grouped path owner never completed work")
        length(unique(collect(grouped.executor.task_ids))) ==
            length(grouped.executor.tasks) ||
            error("two path groups shared one benchmark task")
        grouped.executor.materialization_count ==
            grouped.executor.execution_count ||
            error("grouped materialization/execution accounting differs")
        return Dict{String,Any}(
            "exact_serial_grouped_replay" => true,
            "processed_timestamp_count" => serial_count,
            "path_count" =>
                AOSPlant.path_execution_group_count(serial.prepared),
            "acquisition_count" =>
                AOSPlant.plant_event_acquisition_count(serial.prepared),
            "worker_task_count" => length(grouped.executor.tasks),
            "worker_tasks_reused" => grouped.executor.batch_count > 1,
            "materialized_group_calls" =>
                grouped.executor.materialization_count,
            "executed_group_calls" => grouped.executor.execution_count,
            "products" => Gate3Plant.snapshot_dict(serial_snapshot),
        )
    finally
        Gate6CPU.close_executor!(grouped.executor)
    end
end

function gate6_allocation_evidence(raw, warmup::Int, window::Int, contract)
    _, serial = Gate3Plant.prepare_multi_rate_operation(raw)
    _, grouped_storage = Gate3Plant.prepare_multi_rate_operation(raw)
    grouped = Gate6CPU.GroupedPlantOperation(grouped_storage)
    try
        Gate3Plant.run_timestamp_window!(serial, warmup)
        Gate6CPU.run_timestamp_window!(grouped, warmup)
        GC.gc()
        serial_bytes = Int64(@allocated Gate3Plant.run_timestamp_window!(
            serial, window))
        Gate6CPU.reset_allocation_measurement!(grouped.executor)
        GC.gc()
        grouped_bytes = Int64(@allocated Gate6CPU.run_timestamp_window!(
            grouped, window))
        maximum_group_bytes =
            grouped.executor.maximum_call_allocation_bytes
        serial_limit = Int64(contract[
            "max_serial_alloc_bytes_per_timestamp"]) * window
        grouped_limit = Int64(contract[
            "max_grouped_alloc_bytes_per_timestamp"]) * window
        group_call_limit = Int64(contract[
            "max_alloc_bytes_per_group_call"])
        return Dict{String,Any}(
            "window_timestamps" => window,
            "serial_bytes" => serial_bytes,
            "serial_bytes_per_timestamp" => serial_bytes / window,
            "serial_limit_bytes" => serial_limit,
            "serial_gate_passed" => serial_bytes <= serial_limit,
            "grouped_bytes" => grouped_bytes,
            "grouped_bytes_per_timestamp" => grouped_bytes / window,
            "grouped_limit_bytes" => grouped_limit,
            "grouped_gate_passed" => grouped_bytes <= grouped_limit,
            "measured_group_calls" =>
                grouped.executor.measured_call_count,
            "group_call_bytes" =>
                grouped.executor.measured_allocation_bytes,
            "maximum_group_call_bytes" => maximum_group_bytes,
            "group_call_limit_bytes" => group_call_limit,
            "group_call_gate_passed" =>
                maximum_group_bytes <= group_call_limit,
        )
    finally
        Gate6CPU.close_executor!(grouped.executor)
    end
end

function gate6_command_output(command)
    try
        return readchomp(command)
    catch
        return "unknown"
    end
end

function gate6_optional_file(path::AbstractString)
    try
        return strip(read(path, String))
    catch
        return "unknown"
    end
end

function gate6_allowed_cpu_list()
    status = gate6_optional_file("/proc/self/status")
    status == "unknown" && return status
    for line in eachline(IOBuffer(status))
        startswith(line, "Cpus_allowed_list:") || continue
        return strip(only(split(line, ':'; limit=2)[2:2]))
    end
    return "unknown"
end

function gate6_environment()
    cpu = first(Sys.cpu_info())
    git_status = gate6_command_output(`git status --porcelain=v1`)
    manifest_path = joinpath(@__DIR__, "Manifest.toml")
    return Dict{String,Any}(
        "timestamp_utc" => string(Dates.now(Dates.UTC)),
        "git_commit" => gate6_command_output(`git rev-parse HEAD`),
        "git_dirty" =>
            !isempty(git_status) && git_status != "unknown",
        "git_status_porcelain" => git_status,
        "julia_version" => string(VERSION),
        "adaptive_optics_sim_version" =>
            string(Base.pkgversion(AdaptiveOpticsSim)),
        "hdrhistogram_version" =>
            string(Base.pkgversion(HdrHistogram)),
        "active_project" =>
            something(Base.active_project(), "unknown"),
        "manifest_sha256" => isfile(manifest_path) ?
            bytes2hex(SHA.sha256(read(manifest_path))) : "missing",
        "kernel" => string(Sys.KERNEL),
        "kernel_release" => gate6_command_output(`uname -r`),
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
        "openblas_num_threads_env" =>
            get(ENV, "OPENBLAS_NUM_THREADS", "unset"),
        "omp_num_threads_env" =>
            get(ENV, "OMP_NUM_THREADS", "unset"),
        "mkl_num_threads_env" =>
            get(ENV, "MKL_NUM_THREADS", "unset"),
        "allowed_cpus" => gate6_allowed_cpu_list(),
        "scaling_governor_cpu0" => gate6_optional_file(
            "/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor"),
        "command" => "julia --threads=4 --project=benchmarks " *
            "--startup-file=no " *
            "benchmarks/benchmark_gate6_grouped_cpu.jl",
    )
end

function gate6_timer_overhead(
    samples::Int,
    lowest_ns::Int64,
    highest_ns::Int64,
    significant_figures::Int,
)
    histogram = HdrHistogram.Histogram(
        lowest_ns, highest_ns, significant_figures)
    wall_start = time_ns()
    @inbounds for _ in 1:samples
        start = time_ns()
        HdrHistogram.record_value!(
            histogram, Int64(time_ns() - start))
    end
    wall_ns = Int64(time_ns() - wall_start)
    return Dict{String,Any}(
        "samples" => samples,
        "p50_clock_delta_ns" =>
            HdrHistogram.value_at_percentile(histogram, 50.0),
        "p99_clock_delta_ns" =>
            HdrHistogram.value_at_percentile(histogram, 99.0),
        "wall_ns_per_clock_pair_and_record" => wall_ns / samples,
    )
end

function require_durable_gate6_configuration!(
    output_path::AbstractString,
    contract,
    samples::Int,
    runs::Int,
    warmup::Int,
    allocation_window::Int,
)
    isempty(output_path) && return nothing
    expected = (
        Int(contract["samples_per_run"]),
        Int(contract["runs"]),
        Int(contract["warmup_timestamps"]),
        Int(contract["allocation_timestamps"]),
    )
    (samples, runs, warmup, allocation_window) == expected ||
        error("refusing to write durable Gate 6 evidence with quick overrides")
    samples >= Int(contract["minimum_samples_for_p99"]) ||
        error("durable Gate 6 evidence does not support its p99 claim")
    return nothing
end

function write_gate6_artifact(path::AbstractString, artifact)
    isempty(path) && return nothing
    buffer = IOBuffer()
    TOML.print(buffer, artifact; sorted=true)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, take!(buffer))
    end
    println("wrote Gate 6 grouped CPU artifact: ", path)
    return nothing
end

function gate6_first_use_probe(mode::AbstractString)
    contract = TOML.parsefile(GATE6_CPU_CONTRACT_PATH)
    configure_gate6_cpu!(contract)
    workload = contract["workload"]
    if mode == "serial"
        _, operation = Gate3Plant.prepare_multi_rate_operation(workload)
        GC.gc()
        start = time_ns()
        operation()
        return Int64(time_ns() - start)
    elseif mode == "grouped"
        _, storage = Gate3Plant.prepare_multi_rate_operation(workload)
        operation = Gate6CPU.GroupedPlantOperation(storage)
        try
            GC.gc()
            start = time_ns()
            operation()
            return Int64(time_ns() - start)
        finally
            Gate6CPU.close_executor!(operation.executor)
        end
    end
    error("unknown Gate 6 first-use probe mode: $mode")
end

function gate6_fresh_first_use_ns(mode::AbstractString)
    command = `$(Base.julia_cmd()) --threads=$(Threads.nthreads()) --project=$(joinpath(@__DIR__)) --startup-file=no $(@__FILE__)`
    output = readchomp(addenv(
        command,
        "AOS_GATE6_FIRST_USE_MODE" => String(mode),
        "AOS_GATE6_CPU_CONTRACT" => abspath(GATE6_CPU_CONTRACT_PATH),
        "OPENBLAS_NUM_THREADS" => string(BLAS.get_num_threads()),
        "OMP_NUM_THREADS" => "1",
        "MKL_NUM_THREADS" => "1",
    ))
    prefix = "GATE6_FIRST_USE_NS="
    for line in reverse(split(output, '\n'))
        startswith(line, prefix) || continue
        return parse(Int64, only(split(line, '='; limit=2)[2:2]))
    end
    error("fresh Gate 6 $mode first-use probe returned no timing")
end

function run_gate6_grouped_cpu_benchmark()
    contract = TOML.parsefile(GATE6_CPU_CONTRACT_PATH)
    budget = configure_gate6_cpu!(contract)
    workload = contract["workload"]
    samples = parse(Int, get(ENV, "AOS_GATE6_CPU_SAMPLES",
        string(contract["samples_per_run"])))
    run_count = parse(Int, get(ENV, "AOS_GATE6_CPU_RUNS",
        string(contract["runs"])))
    warmup = parse(Int, get(ENV, "AOS_GATE6_CPU_WARMUP",
        string(contract["warmup_timestamps"])))
    allocation_window = parse(Int, get(
        ENV,
        "AOS_GATE6_CPU_ALLOCATION_TIMESTAMPS",
        string(contract["allocation_timestamps"]),
    ))
    samples > 0 || error("AOS_GATE6_CPU_SAMPLES must be positive")
    run_count > 0 || error("AOS_GATE6_CPU_RUNS must be positive")
    warmup >= 0 || error("AOS_GATE6_CPU_WARMUP must be nonnegative")
    allocation_window > 0 ||
        error("AOS_GATE6_CPU_ALLOCATION_TIMESTAMPS must be positive")
    require_durable_gate6_configuration!(
        GATE6_CPU_OUTPUT_PATH,
        contract,
        samples,
        run_count,
        warmup,
        allocation_window,
    )
    lowest_ns = Int64(contract["histogram_lowest_ns"])
    highest_ns = Int64(contract["histogram_highest_ns"])
    significant_figures =
        Int(contract["histogram_significant_figures"])
    p99_supported =
        samples >= Int(contract["minimum_samples_for_p99"])

    println("gate6_grouped_cpu_contract")
    println("  load_model: paired closed-loop self-paced service cost")
    println("  samples_per_mode_per_run: ", samples)
    println("  runs: ", run_count)
    println("  p99_claim_supported: ", p99_supported)

    first_serial_ns = gate6_fresh_first_use_ns("serial")
    first_grouped_ns = gate6_fresh_first_use_ns("grouped")

    correctness = validate_serial_grouped(
        workload, contract["correctness_horizon_ns"])
    allocation = gate6_allocation_evidence(
        workload, warmup, allocation_window, contract)
    println(
        "  serial_allocation_bytes_per_timestamp=",
        allocation["serial_bytes_per_timestamp"],
    )
    println(
        "  grouped_allocation_bytes_per_timestamp=",
        allocation["grouped_bytes_per_timestamp"],
    )
    println(
        "  maximum_group_call_allocation_bytes=",
        allocation["maximum_group_call_bytes"],
    )

    serial_runs = Vector{Dict{String,Any}}(undef, run_count)
    grouped_runs = Vector{Dict{String,Any}}(undef, run_count)
    for run_index in 1:run_count
        serial_plant, serial =
            Gate3Plant.prepare_multi_rate_operation(workload)
        grouped_plant, grouped_storage =
            Gate3Plant.prepare_multi_rate_operation(workload)
        grouped = Gate6CPU.GroupedPlantOperation(grouped_storage)
        try
            Gate3Plant.run_timestamp_window!(serial, warmup)
            Gate6CPU.run_timestamp_window!(grouped, warmup)
            if isodd(run_index)
                serial_runs[run_index] = measure_gate6_run!(
                    serial, samples, lowest_ns, highest_ns,
                    significant_figures)
                grouped_runs[run_index] = measure_gate6_run!(
                    grouped, samples, lowest_ns, highest_ns,
                    significant_figures)
            else
                grouped_runs[run_index] = measure_gate6_run!(
                    grouped, samples, lowest_ns, highest_ns,
                    significant_figures)
                serial_runs[run_index] = measure_gate6_run!(
                    serial, samples, lowest_ns, highest_ns,
                    significant_figures)
            end
            Gate3Plant.product_snapshot(serial_plant, serial) ==
                Gate3Plant.product_snapshot(
                    grouped_plant, grouped_storage) ||
                error("paired timing runs ended with different products")
        finally
            Gate6CPU.close_executor!(grouped.executor)
        end
        println(
            "  run=", run_index,
            " serial_p50_ns=", serial_runs[run_index]["p50_ns"],
            " grouped_p50_ns=", grouped_runs[run_index]["p50_ns"],
        )
    end

    serial_summary = summarize_gate6_runs(serial_runs)
    grouped_summary = summarize_gate6_runs(grouped_runs)
    p50_ratios = [
        grouped_runs[index]["p50_ns"] /
            serial_runs[index]["p50_ns"]
        for index in eachindex(serial_runs)
    ]
    p99_ratios = [
        grouped_runs[index]["p99_ns"] /
            serial_runs[index]["p99_ns"]
        for index in eachindex(serial_runs)
    ]
    comparison = Dict{String,Any}(
        "median_p50_ratio_grouped_over_serial" =>
            median(p50_ratios),
        "minimum_p50_ratio_grouped_over_serial" =>
            minimum(p50_ratios),
        "maximum_p50_ratio_grouped_over_serial" =>
            maximum(p50_ratios),
        "median_p99_ratio_grouped_over_serial" =>
            median(p99_ratios),
        "minimum_p99_ratio_grouped_over_serial" =>
            minimum(p99_ratios),
        "maximum_p99_ratio_grouped_over_serial" =>
            maximum(p99_ratios),
        "run_order" =>
            "serial-first for odd runs; grouped-first for even runs",
        "performance_gate_evaluated" => false,
        "performance_gate_passed" => true,
        "interpretation" =>
            "diagnostic paired service-cost comparison; no speedup claim",
    )
    println(
        "  grouped_over_serial_median_p50_ratio=",
        comparison["median_p50_ratio_grouped_over_serial"],
    )
    all_gates_passed =
        allocation["serial_gate_passed"] &&
        allocation["grouped_gate_passed"] &&
        allocation["group_call_gate_passed"]
    environment = gate6_environment()
    if !isempty(GATE6_CPU_OUTPUT_PATH) && environment["git_dirty"]
        error("refusing to write durable Gate 6 evidence from a dirty worktree")
    end
    artifact = Dict{String,Any}(
        "schema_version" => Int(contract["schema_version"]),
        "benchmark" => String(contract["name"]),
        "evidence_class" => String(contract["evidence_class"]),
        "source_contract" => relpath(
            abspath(GATE6_CPU_CONTRACT_PATH), @__DIR__),
        "characterized_source_revision" => environment["git_commit"],
        "configured_samples_per_mode_per_run" => samples,
        "configured_runs" => run_count,
        "configured_warmup_timestamps" => warmup,
        "minimum_samples_for_p99" =>
            Int(contract["minimum_samples_for_p99"]),
        "p99_claim_supported" => p99_supported,
        "contract" => contract["contract"],
        "scope_exclusions" => contract["scope_exclusions"],
        "workload" => workload,
        "cpu_execution_budget" => Dict(
            "cpu_contexts" => budget.cpu_context_count,
            "julia_threads" => budget.julia_thread_count,
            "outer_group_owners" => budget.outer_owner_count,
            "group_julia_threads" =>
                budget.group_julia_thread_count,
            "fft_threads_per_owner" => budget.fft_thread_count,
            "blas_threads_per_owner" => budget.blas_thread_count,
        ),
        "histogram" => Dict(
            "lowest_ns" => lowest_ns,
            "highest_ns" => highest_ns,
            "significant_figures" => significant_figures,
        ),
        "timer_recording_overhead" => gate6_timer_overhead(
            max(samples, 10_000), lowest_ns, highest_ns,
            significant_figures),
        "environment" => environment,
        "first_use" => Dict(
            "fresh_process_serial_ns" => first_serial_ns,
            "fresh_process_grouped_ns" => first_grouped_ns,
            "method" =>
                "each mode measured after preparation in a separate fresh Julia process",
        ),
        "correctness" => correctness,
        "allocation" => allocation,
        "serial_runs" => serial_runs,
        "grouped_runs" => grouped_runs,
        "serial_summary" => serial_summary,
        "grouped_summary" => grouped_summary,
        "comparison" => comparison,
        "all_gates_passed" => all_gates_passed,
    )
    all_gates_passed ||
        error("Gate 6 grouped CPU allocation gate failed")
    write_gate6_artifact(GATE6_CPU_OUTPUT_PATH, artifact)
    return artifact
end

if abspath(PROGRAM_FILE) == @__FILE__
    first_use_mode = get(ENV, "AOS_GATE6_FIRST_USE_MODE", "")
    if isempty(first_use_mode)
        run_gate6_grouped_cpu_benchmark()
    else
        println(
            "GATE6_FIRST_USE_NS=",
            gate6_first_use_probe(first_use_mode),
        )
    end
end

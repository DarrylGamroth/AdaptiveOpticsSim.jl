using AdaptiveOpticsSim
using Dates
using HdrHistogram
import KernelAbstractions
using LinearAlgebra
using SHA
using Statistics
using TOML

const GATE7_BACKEND_NAME = lowercase(get(
    ENV,
    "AOS_GATE7_BACKEND",
    isempty(ARGS) ? "cpu" : ARGS[1],
))
const GATE7_PLACEMENT = get(
    ENV,
    "AOS_GATE7_PLACEMENT",
    length(ARGS) >= 2 ? ARGS[2] : "",
)
const GATE7_CONTRACT_PATH = get(
    ENV,
    "AOS_GATE7_CONTRACT",
    joinpath(@__DIR__, "contracts", "gate7_single_gpu.toml"),
)
const GATE7_OUTPUT_PATH = get(ENV, "AOS_GATE7_OUTPUT", "")
const GATE7_ROOT = normpath(joinpath(@__DIR__, ".."))

if GATE7_BACKEND_NAME == "cuda"
    @eval import CUDA
elseif GATE7_BACKEND_NAME == "amdgpu"
    @eval import AMDGPU
elseif GATE7_BACKEND_NAME != "cpu"
    error(
        "unsupported AOS_GATE7_BACKEND '$GATE7_BACKEND_NAME'; " *
        "use cpu, cuda, or amdgpu",
    )
end

include(joinpath(@__DIR__, "support", "gate7_single_gpu.jl"))
include(joinpath(@__DIR__, "support", "hdr_histogram_artifact.jl"))

const Gate7Workload = Gate7SingleGPUBenchmark
const HdrArtifact = HdrHistogramArtifact

function gate7_command_output(command)
    try
        return readchomp(command)
    catch
        return "unknown"
    end
end

function gate7_optional_file(path::AbstractString)
    try
        return strip(read(path, String))
    catch
        return "unknown"
    end
end

function gate7_allowed_cpu_list()
    status = gate7_optional_file("/proc/self/status")
    status == "unknown" && return status
    for line in eachline(IOBuffer(status))
        startswith(line, "Cpus_allowed_list:") || continue
        return strip(last(split(line, ':'; limit=2)))
    end
    return "unknown"
end

function gate7_sha256_file(path::AbstractString)
    isfile(path) || return "unavailable"
    return bytes2hex(SHA.sha256(read(path)))
end

function gate7_active_manifest_path(project_path::AbstractString)
    versioned = joinpath(
        dirname(project_path),
        "Manifest-v$(VERSION.major).$(VERSION.minor).toml",
    )
    isfile(versioned) && return versioned
    return joinpath(dirname(project_path), "Manifest.toml")
end

function configure_gate7_benchmark!()
    Threads.nthreads() == 1 ||
        error("Gate 7 benchmark requires exactly one Julia thread")
    BLAS.set_num_threads(1)
    AdaptiveOpticsSim.set_fft_provider_threads!(1)
    return nothing
end

function gate7_backend_selector()
    if GATE7_BACKEND_NAME == "cpu"
        return CPUBackend()
    elseif GATE7_BACKEND_NAME == "cuda"
        CUDA.functional() ||
            error("Gate 7 CUDA benchmark requires a functional device")
        AdaptiveOpticsSim.disable_scalar_backend!(
            AdaptiveOpticsSim.CUDABackendTag,
        )
        return CUDABackend()
    end
    AMDGPU.functional() ||
        error("Gate 7 AMDGPU benchmark requires a functional device")
    AdaptiveOpticsSim.disable_scalar_backend!(
        AdaptiveOpticsSim.AMDGPUBackendTag,
    )
    return AMDGPUBackend()
end

function gate7_backend_module()
    GATE7_BACKEND_NAME == "cuda" && return CUDA
    GATE7_BACKEND_NAME == "amdgpu" && return AMDGPU
    return nothing
end

function gate7_nvidia_smi()
    arguments = (
        "--query-gpu=name,driver_version,pstate,memory.total",
        "--format=csv,noheader",
    )
    for executable in ("nvidia-smi", "/usr/lib/wsl/lib/nvidia-smi")
        output = gate7_command_output(`$executable $arguments`)
        output == "unknown" || return output
    end
    return "unknown"
end

function gate7_hip_driver_version_code()
    version = Ref{Cint}()
    status = AMDGPU.HIP.hipDriverGetVersion(version)
    status == AMDGPU.HIP.hipSuccess ||
        error("hipDriverGetVersion failed with status $status")
    return Int(version[])
end

function gate7_accelerator_environment()
    GATE7_BACKEND_NAME == "cpu" && return Dict{String,Any}(
        "device" => "not applicable",
    )
    if GATE7_BACKEND_NAME == "cuda"
        device = CUDA.device()
        return Dict{String,Any}(
            "device" => CUDA.name(device),
            "compute_capability" => string(CUDA.capability(device)),
            "runtime_version" => string(CUDA.runtime_version()),
            "driver_version" => string(CUDA.driver_version()),
            "compiler_version" => string(CUDA.compiler_version()),
            "nvidia_smi" => gate7_nvidia_smi(),
        )
    end
    device = AMDGPU.device()
    return Dict{String,Any}(
        "device" => AMDGPU.HIP.name(device),
        "gcn_architecture" => AMDGPU.HIP.gcn_arch(device),
        "wavefront_size" => Int(AMDGPU.HIP.wavefrontsize(device)),
        "hip_driver_version_code" => gate7_hip_driver_version_code(),
        "hip_runtime_version" => string(AMDGPU.HIP.runtime_version()),
    )
end

function gate7_source_environment()
    cpu = first(Sys.cpu_info())
    git_status = gate7_command_output(
        `git -C $GATE7_ROOT status --porcelain=v1`,
    )
    project_path = something(Base.active_project(), "unknown")
    manifest_path = project_path == "unknown" ? "unknown" :
        gate7_active_manifest_path(project_path)
    backend_module = gate7_backend_module()
    benchmark_project = GATE7_BACKEND_NAME == "cpu" ?
        "benchmarks" : "benchmarks/$GATE7_BACKEND_NAME"
    return Dict{String,Any}(
        "timestamp_utc" => string(Dates.now(Dates.UTC)),
        "host_name" => gate7_command_output(`hostname`),
        "git_commit" => gate7_command_output(
            `git -C $GATE7_ROOT rev-parse HEAD`,
        ),
        "git_branch" => gate7_command_output(
            `git -C $GATE7_ROOT branch --show-current`,
        ),
        "git_dirty" => !isempty(git_status) && git_status != "unknown",
        "git_status_porcelain" => git_status,
        "julia_version" => string(VERSION),
        "adaptive_optics_sim_version" =>
            string(Base.pkgversion(AdaptiveOpticsSim)),
        "hdrhistogram_version" => string(Base.pkgversion(HdrHistogram)),
        "kernel_abstractions_version" =>
            string(Base.pkgversion(KernelAbstractions)),
        "accelerator_package_version" => backend_module === nothing ?
            "not applicable" : string(Base.pkgversion(backend_module)),
        "accelerator" => gate7_accelerator_environment(),
        "active_project" => project_path,
        "active_project_sha256" => gate7_sha256_file(project_path),
        "active_manifest" => manifest_path,
        "active_manifest_sha256" =>
            gate7_sha256_file(manifest_path),
        "kernel" => string(Sys.KERNEL),
        "kernel_release" => gate7_command_output(`uname -r`),
        "architecture" => string(Sys.ARCH),
        "cpu_target" => string(Sys.CPU_NAME),
        "cpu_model" => cpu.model,
        "logical_cpu_threads" => Sys.CPU_THREADS,
        "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
        "blas_config" => string(BLAS.get_config()),
        "fft_threads" => 1,
        "allowed_cpus" => gate7_allowed_cpu_list(),
        "scaling_governor_cpu0" => gate7_optional_file(
            "/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor",
        ),
        "julia_cpu_target_env" =>
            get(ENV, "JULIA_CPU_TARGET", "default"),
        "julia_depot_path_env" =>
            get(ENV, "JULIA_DEPOT_PATH", "default"),
        "julia_num_precompile_tasks_env" => get(
            ENV,
            "JULIA_NUM_PRECOMPILE_TASKS",
            "default",
        ),
        "backend" => GATE7_BACKEND_NAME,
        "placement" => GATE7_PLACEMENT,
        "command" =>
            "julia --threads=1 --project=$benchmark_project " *
            "--startup-file=no benchmarks/benchmark_gate7_single_gpu.jl " *
            "$GATE7_BACKEND_NAME $GATE7_PLACEMENT",
    )
end

function require_clean_gate7_source_for_output!()
    isempty(GATE7_OUTPUT_PATH) && return nothing
    status = gate7_command_output(
        `git -C $GATE7_ROOT status --porcelain=v1`,
    )
    isempty(status) || error(
        "refusing durable Gate 7 evidence from a dirty worktree",
    )
    return nothing
end

function require_durable_gate7_configuration!()
    isempty(GATE7_OUTPUT_PATH) && return nothing
    for name in (
        "AOS_GATE7_SAMPLES",
        "AOS_GATE7_RUNS",
        "AOS_GATE7_WARMUP",
    )
        haskey(ENV, name) && error(
            "durable Gate 7 evidence must use the contract default $name",
        )
    end
    return nothing
end

function gate7_gc_delta(before, after)
    return Dict(
        String(name) =>
            Int64(getfield(after, name) - getfield(before, name))
        for name in fieldnames(typeof(before))
    )
end

function gate7_histogram_summary(
    histogram::HdrHistogram.Histogram,
    wall_start_ns::UInt64,
    wall_end_ns::UInt64,
    samples::Int,
    gc_counters,
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
        "min_ns" => min(histogram),
        "mean_ns" => HdrHistogram.mean(histogram),
        "p50_ns" =>
            HdrHistogram.value_at_percentile(histogram, 50.0),
        "p90_ns" =>
            HdrHistogram.value_at_percentile(histogram, 90.0),
        "p95_ns" =>
            HdrHistogram.value_at_percentile(histogram, 95.0),
        "max_ns" => max(histogram),
        "gc" => gc_counters,
        encoded...,
    )
end

function measure_gate7_boundary_run!(
    boundary,
    samples::Int,
    lowest_ns::Int64,
    highest_ns::Int64,
    significant_figures::Int,
)
    histogram = HdrHistogram.Histogram(
        lowest_ns,
        highest_ns,
        significant_figures,
    )
    GC.gc()
    gc_before = Base.gc_num()
    wall_start = time_ns()
    @inbounds for _ in 1:samples
        start = time_ns()
        Gate7Workload.run_gate7_boundary!(boundary)
        elapsed = Int64(time_ns() - start)
        HdrHistogram.record_value!(histogram, elapsed)
    end
    wall_end = time_ns()
    gc_after = Base.gc_num()
    return gate7_histogram_summary(
        histogram,
        wall_start,
        wall_end,
        samples,
        gate7_gc_delta(gc_before, gc_after),
        lowest_ns,
        highest_ns,
        significant_figures,
    )
end

gate7_median_integer(values) = round(Int64, median(collect(values)))

function summarize_gate7_runs(runs::Vector{Dict{String,Any}})
    return Dict{String,Any}(
        "median_p50_ns" =>
            gate7_median_integer(run["p50_ns"] for run in runs),
        "minimum_p50_ns" =>
            minimum(run["p50_ns"] for run in runs),
        "maximum_p50_ns" =>
            maximum(run["p50_ns"] for run in runs),
        "median_p90_ns" =>
            gate7_median_integer(run["p90_ns"] for run in runs),
        "minimum_p90_ns" =>
            minimum(run["p90_ns"] for run in runs),
        "maximum_p90_ns" =>
            maximum(run["p90_ns"] for run in runs),
        "median_p95_ns" =>
            gate7_median_integer(run["p95_ns"] for run in runs),
        "minimum_p95_ns" =>
            minimum(run["p95_ns"] for run in runs),
        "maximum_p95_ns" =>
            maximum(run["p95_ns"] for run in runs),
        "maximum_observed_ns" =>
            maximum(run["max_ns"] for run in runs),
        "median_throughput_hz" =>
            median(collect(run["throughput_hz"] for run in runs)),
    )
end

function gate7_timer_recording_overhead(
    samples::Int,
    lowest_ns::Int64,
    highest_ns::Int64,
    significant_figures::Int,
)
    histogram = HdrHistogram.Histogram(
        lowest_ns,
        highest_ns,
        significant_figures,
    )
    wall_start = time_ns()
    @inbounds for _ in 1:samples
        start = time_ns()
        elapsed = Int64(time_ns() - start)
        HdrHistogram.record_value!(histogram, elapsed)
    end
    wall_ns = Int64(time_ns() - wall_start)
    return Dict{String,Any}(
        "samples" => samples,
        "p50_clock_delta_ns" =>
            HdrHistogram.value_at_percentile(histogram, 50.0),
        "p95_clock_delta_ns" =>
            HdrHistogram.value_at_percentile(histogram, 95.0),
        "wall_ns_per_clock_pair_and_record" => wall_ns / samples,
    )
end

function gate7_first_use!(
    boundary_ids::Vector{String},
    selector,
    workload,
)
    boundaries = Dict{String,Any}()
    records = Vector{Dict{String,Any}}(undef, length(boundary_ids))
    @inbounds for index in eachindex(boundary_ids)
        id = boundary_ids[index]
        boundary =
            Gate7Workload.prepare_gate7_boundary(id, selector, workload)
        start = time_ns()
        Gate7Workload.run_gate7_boundary!(boundary)
        elapsed_ns = Int64(time_ns() - start)
        boundaries[id] = boundary
        records[index] = Dict{String,Any}(
            "id" => id,
            "elapsed_ns" => elapsed_ns,
            "declared_order" => index,
        )
    end
    return boundaries, records
end

function gate7_cpu_correctness(operation)
    operation.completed_steps == 1 ||
        error("Gate 7 CPU first-use context did not complete one step")
    AdaptiveOpticsSim.Plant.device_path_batch_owner_count(
        operation.prepared,
    ) == 0 ||
        error("Gate 7 CPU context unexpectedly has a device owner")
    paths = Dict{String,Any}()
    @inbounds for id in Gate7Workload.gate7_path_ids()
        values = Array(
            AdaptiveOpticsSim.Plant.prepared_path(
                operation.plant,
                id,
            ).result.values,
        )
        all(isfinite, values) ||
            error("Gate 7 CPU product $id contains non-finite values")
        paths[String(id)] = Dict{String,Any}(
            "shape" => collect(size(values)),
            "checksum" => sum(Float64, values),
        )
    end
    return Dict{String,Any}(
        "passed" => true,
        "mode" => "independent CPU contextual oracle",
        "completed_steps" => operation.completed_steps,
        "paths" => paths,
    )
end

function gate7_boundary_results(
    boundary_ids::Vector{String},
    selector,
    workload,
    placement,
    samples::Int,
    run_count::Int,
    warmup::Int,
    minimum_samples_for_p95::Int,
    lowest_ns::Int64,
    highest_ns::Int64,
    significant_figures::Int,
)
    results = Vector{Dict{String,Any}}(undef, length(boundary_ids))
    p95_supported = samples >= minimum_samples_for_p95
    @inbounds for boundary_index in eachindex(boundary_ids)
        id = boundary_ids[boundary_index]
        runs = Vector{Dict{String,Any}}(undef, run_count)
        allocations = Vector{Int64}(undef, run_count)
        for run_index in 1:run_count
            boundary = Gate7Workload.prepare_gate7_boundary(
                id,
                selector,
                workload,
            )
            for _ in 1:warmup
                Gate7Workload.run_gate7_boundary!(boundary)
            end
            GC.gc()
            allocated = @allocated Gate7Workload.run_gate7_boundary!(
                boundary,
            )
            allocations[run_index] = Int64(allocated)
            runs[run_index] = measure_gate7_boundary_run!(
                boundary,
                samples,
                lowest_ns,
                highest_ns,
                significant_figures,
            )
            println(
                id,
                " run=",
                run_index,
                " p50_ns=",
                runs[run_index]["p50_ns"],
                " p95_ns=",
                runs[run_index]["p95_ns"],
                " max_ns=",
                runs[run_index]["max_ns"],
                " alloc_bytes=",
                allocations[run_index],
            )
        end
        summary = summarize_gate7_runs(runs)
        config = placement[id]
        allocation_limit = Int64(config["max_alloc_bytes"])
        p95_limit = Int64(config["absolute_p95_ns"])
        maximum_allocation = maximum(allocations)
        absolute_observed = Int64(summary["maximum_p95_ns"])
        gate = Dict{String,Any}(
            "p95_evaluated" => p95_supported,
            "absolute_p95_passed" => !p95_supported ||
                absolute_observed <= p95_limit,
            "absolute_p95_observed_ns" => absolute_observed,
            "absolute_p95_limit_ns" => p95_limit,
            "allocation_passed" =>
                maximum_allocation <= allocation_limit,
            "allocation_observed_bytes" => maximum_allocation,
            "allocation_limit_bytes" => allocation_limit,
        )
        results[boundary_index] = Dict{String,Any}(
            "id" => id,
            "runs" => runs,
            "warmed_allocation_bytes" => allocations,
            "summary" => summary,
            "gate" => gate,
            "passed" => gate["absolute_p95_passed"] &&
                gate["allocation_passed"],
        )
    end
    return results
end

function gate7_boundary_result(results, id::AbstractString)
    index = findfirst(result -> result["id"] == id, results)
    isnothing(index) && error("missing Gate 7 boundary result '$id'")
    return results[index]
end

function gate7_relative_comparison(results, factor::Float64)
    ids = Set(String(result["id"]) for result in results)
    required = Set(("independent_device_ready", "batched_device_ready"))
    if !issubset(required, ids)
        return Dict{String,Any}(
            "evaluated" => false,
            "reason" => "placement has no prepared accelerator owner",
            "passed" => true,
            "limit_factor" => factor,
        )
    end
    p95_supported = all(
        result["gate"]["p95_evaluated"] for result in results
        if result["id"] in required
    )
    if !p95_supported
        return Dict{String,Any}(
            "evaluated" => false,
            "reason" => "sample count does not support the declared p95 gate",
            "passed" => true,
            "limit_factor" => factor,
        )
    end
    independent = gate7_boundary_result(
        results,
        "independent_device_ready",
    )["summary"]["median_p95_ns"]
    batched = gate7_boundary_result(
        results,
        "batched_device_ready",
    )["summary"]["median_p95_ns"]
    ratio = batched / independent
    return Dict{String,Any}(
        "evaluated" => true,
        "independent_median_p95_ns" => independent,
        "batched_median_p95_ns" => batched,
        "batched_to_independent_p95_ratio" => ratio,
        "limit_factor" => factor,
        "passed" => ratio <= factor,
    )
end

function write_gate7_artifact(path::AbstractString, artifact)
    isempty(path) && return nothing
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, artifact; sorted=true)
    end
    println("wrote Gate 7 artifact: ", path)
    return nothing
end

function run_gate7_single_gpu_benchmark()
    configure_gate7_benchmark!()
    require_clean_gate7_source_for_output!()
    require_durable_gate7_configuration!()
    isempty(GATE7_PLACEMENT) && error(
        "set AOS_GATE7_PLACEMENT or pass placement as the second argument",
    )
    contract = TOML.parsefile(GATE7_CONTRACT_PATH)
    placements = contract["placements"]
    haskey(placements, GATE7_PLACEMENT) ||
        error("unknown Gate 7 placement '$GATE7_PLACEMENT'")
    placement = placements[GATE7_PLACEMENT]
    String(placement["backend"]) == GATE7_BACKEND_NAME || error(
        "placement '$GATE7_PLACEMENT' requires backend " *
        "$(placement["backend"]), not $GATE7_BACKEND_NAME",
    )
    samples = parse(
        Int,
        get(
            ENV,
            "AOS_GATE7_SAMPLES",
            string(contract["samples_per_run"]),
        ),
    )
    run_count = parse(
        Int,
        get(ENV, "AOS_GATE7_RUNS", string(contract["runs"])),
    )
    warmup = parse(
        Int,
        get(
            ENV,
            "AOS_GATE7_WARMUP",
            string(contract["warmup_operations"]),
        ),
    )
    samples > 0 || error("AOS_GATE7_SAMPLES must be > 0")
    run_count > 0 || error("AOS_GATE7_RUNS must be > 0")
    warmup >= 0 || error("AOS_GATE7_WARMUP must be >= 0")

    minimum_samples_for_p95 =
        Int(contract["minimum_samples_for_p95"])
    lowest_ns = Int64(contract["histogram_lowest_ns"])
    highest_ns = Int64(contract["histogram_highest_ns"])
    significant_figures =
        Int(contract["histogram_significant_figures"])
    relative_factor =
        Float64(contract["batched_relative_p95_factor"])
    workload = contract["workload"]
    boundary_ids = String.(placement["boundaries"])
    selector = gate7_backend_selector()
    environment = gate7_source_environment()
    if !isempty(GATE7_OUTPUT_PATH) && environment["git_dirty"]
        error("refusing durable Gate 7 evidence from a dirty worktree")
    end

    println("gate7_single_gpu")
    println("  placement: ", GATE7_PLACEMENT)
    println("  backend: ", GATE7_BACKEND_NAME)
    println("  samples_per_run: ", samples)
    println("  runs: ", run_count)
    println("  warmup_operations: ", warmup)
    println("  boundaries: ", join(boundary_ids, ','))

    first_boundaries, first_use = gate7_first_use!(
        boundary_ids,
        selector,
        workload,
    )
    independent_first =
        first_boundaries["independent_device_ready"].operation
    correctness = if GATE7_BACKEND_NAME == "cpu"
        gate7_cpu_correctness(independent_first)
    else
        batched_first =
            first_boundaries["batched_device_ready"].operation
        Gate7Workload.gate7_correctness_evidence(
            independent_first,
            batched_first;
            rtol=Float64(workload["parity_rtol"]),
            atol=Float64(workload["parity_atol"]),
        )
    end
    submission_proxy = if GATE7_BACKEND_NAME == "cpu"
        Dict{String,Any}(
            "passed" => true,
            "derivation" =>
                "prepared CPU path topology; contextual oracle only",
            "observed" => Dict{String,Any}(
                "independent" =>
                    Gate7Workload.gate7_submission_proxy(
                        independent_first,
                        Val(:independent),
                    ),
            ),
        )
    else
        Gate7Workload.gate7_submission_proxy_evidence(
            independent_first,
            first_boundaries["batched_device_ready"].operation,
            contract["submission_proxy"],
        )
    end
    residency = Dict{String,Any}()
    @inbounds for id in boundary_ids
        residency[id] = Gate7Workload.gate7_residency(
            first_boundaries[id].operation,
        )
    end

    results = gate7_boundary_results(
        boundary_ids,
        selector,
        workload,
        placement,
        samples,
        run_count,
        warmup,
        minimum_samples_for_p95,
        lowest_ns,
        highest_ns,
        significant_figures,
    )
    comparison = gate7_relative_comparison(results, relative_factor)
    all_gates_passed = correctness["passed"] &&
        submission_proxy["passed"] &&
        all(result["passed"] for result in results) &&
        comparison["passed"]

    if !isempty(GATE7_OUTPUT_PATH)
        status = gate7_command_output(
            `git -C $GATE7_ROOT status --porcelain=v1`,
        )
        isempty(status) || error(
            "refusing to write durable Gate 7 evidence after source changed",
        )
    end
    artifact = Dict{String,Any}(
        "schema_version" => Int(contract["schema_version"]),
        "benchmark" => String(contract["name"]),
        "evidence_class" => String(contract["evidence_class"]),
        "source_contract" =>
            relpath(abspath(GATE7_CONTRACT_PATH), GATE7_ROOT),
        "characterized_source_revision" => environment["git_commit"],
        "placement" => GATE7_PLACEMENT,
        "backend" => GATE7_BACKEND_NAME,
        "host_class" => String(placement["host_class"]),
        "configured_samples_per_run" => samples,
        "configured_runs" => run_count,
        "configured_warmup_operations" => warmup,
        "minimum_samples_for_p95" => minimum_samples_for_p95,
        "p95_supported" => samples >= minimum_samples_for_p95,
        "contract" => contract["contract"],
        "scope_exclusions" => contract["scope_exclusions"],
        "workload" => workload,
        "declared_submission_proxy" => contract["submission_proxy"],
        "submission_proxy_evidence" => submission_proxy,
        "histogram" => Dict(
            "lowest_ns" => lowest_ns,
            "highest_ns" => highest_ns,
            "significant_figures" => significant_figures,
        ),
        "timer_recording_overhead" =>
            gate7_timer_recording_overhead(
                max(samples, 10_000),
                lowest_ns,
                highest_ns,
                significant_figures,
            ),
        "environment" => environment,
        "first_use_order" => first_use,
        "correctness" => correctness,
        "residency" => residency,
        "boundaries" => results,
        "relative_comparison" => comparison,
        "all_gates_passed" => all_gates_passed,
    )
    write_gate7_artifact(GATE7_OUTPUT_PATH, artifact)
    all_gates_passed ||
        error("one or more Gate 7 benchmark gates failed")
    return artifact
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_gate7_single_gpu_benchmark()
end

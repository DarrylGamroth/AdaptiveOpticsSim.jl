using AdaptiveOpticsSim
using Base64
using Dates
using HdrHistogram
import KernelAbstractions
using LinearAlgebra
using Statistics
using TOML

const PRE_HIL_BACKEND_NAME = lowercase(get(ENV, "AOS_PRE_HIL_BACKEND",
    isempty(ARGS) ? "cpu" : ARGS[1]))
const PRE_HIL_PLACEMENT = get(ENV, "AOS_PRE_HIL_PLACEMENT",
    length(ARGS) >= 2 ? ARGS[2] : "")
const PRE_HIL_CONTRACT_PATH = get(ENV, "AOS_PRE_HIL_BACKEND_CONTRACT",
    joinpath(@__DIR__, "contracts", "pre_hil_backend_latency.toml"))
const PRE_HIL_OUTPUT_PATH = get(ENV, "AOS_PRE_HIL_BACKEND_OUTPUT", "")
const PRE_HIL_ROOT = normpath(joinpath(@__DIR__, ".."))
const PRE_HIL_CONFIG_DIR = joinpath(@__DIR__, "assets", "revolt_like")

if PRE_HIL_BACKEND_NAME == "cuda"
    @eval using CUDA
elseif PRE_HIL_BACKEND_NAME == "amdgpu"
    @eval using AMDGPU
elseif PRE_HIL_BACKEND_NAME != "cpu"
    error("unsupported AOS_PRE_HIL_BACKEND '$PRE_HIL_BACKEND_NAME'; use cpu, cuda, or amdgpu")
end

include(joinpath(@__DIR__, "support", "revolt_like_hil_common.jl"))

struct PreHILBoundary{F}
    id::String
    operation::F
end

@inline run_pre_hil_boundary!(boundary::PreHILBoundary) = boundary.operation()

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
        return strip(last(split(line, ':'; limit=2)))
    end
    return "unknown"
end

function configure_pre_hil_benchmark!()
    Threads.nthreads() == 1 || error(
        "pre-HIL backend latency evidence requires one Julia thread")
    BLAS.set_num_threads(1)
    AdaptiveOpticsSim.set_fft_provider_threads!(1)
    return nothing
end

function backend_module()
    PRE_HIL_BACKEND_NAME == "cuda" && return CUDA
    PRE_HIL_BACKEND_NAME == "amdgpu" && return AMDGPU
    return nothing
end

function backend_version_info()
    module_ref = backend_module()
    module_ref === nothing && return "not applicable"
    versioninfo = getproperty(module_ref, :versioninfo)
    try
        return sprint(io -> Base.invokelatest(versioninfo, io))
    catch err
        return "versioninfo unavailable: $(sprint(showerror, err))"
    end
end

function source_environment()
    cpu = first(Sys.cpu_info())
    git_status = command_output(`git -C $PRE_HIL_ROOT status --porcelain=v1`)
    module_ref = backend_module()
    return Dict{String,Any}(
        "timestamp_utc" => string(Dates.now(Dates.UTC)),
        "host_name" => command_output(`hostname`),
        "git_commit" => command_output(`git -C $PRE_HIL_ROOT rev-parse HEAD`),
        "git_dirty" => !isempty(git_status) && git_status != "unknown",
        "git_status_porcelain" => git_status,
        "julia_version" => string(VERSION),
        "adaptive_optics_sim_version" => string(Base.pkgversion(AdaptiveOpticsSim)),
        "hdrhistogram_version" => string(Base.pkgversion(HdrHistogram)),
        "kernel_abstractions_version" => string(Base.pkgversion(KernelAbstractions)),
        "accelerator_package_version" => module_ref === nothing ?
            "not applicable" : string(Base.pkgversion(module_ref)),
        "accelerator_versioninfo" => backend_version_info(),
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
        "allowed_cpus" => allowed_cpu_list(),
        "scaling_governor_cpu0" => optional_file(
            "/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor"),
        "julia_cpu_target_env" => get(ENV, "JULIA_CPU_TARGET", "default"),
        "backend" => PRE_HIL_BACKEND_NAME,
        "placement" => PRE_HIL_PLACEMENT,
        "command" => "julia --startup-file=no --project=<benchmark-project> benchmarks/benchmark_pre_hil_backend_latency.jl $PRE_HIL_BACKEND_NAME $PRE_HIL_PLACEMENT",
    )
end

function assert_backend_residency(ctx::RevoltLikeHILContext)
    arrays = Dict(
        "telescope_opd" => ctx.tel.state.opd,
        "dm_coefficients" => ctx.dm.state.coefs,
        "wfs_spot_cube" => ctx.wfs.acquisition.spot_cube,
        "output_frame" => ctx.tiled_frame,
    )
    if PRE_HIL_BACKEND_NAME == "cpu"
        all(array isa Array for array in values(arrays)) || error(
            "CPU placement contains a non-Array maintained working buffer")
    else
        any(array isa Array for array in values(arrays)) && error(
            "$PRE_HIL_BACKEND_NAME placement contains a host Array in a maintained backend-resident working buffer")
    end
    return Dict(name => string(typeof(array)) for (name, array) in arrays)
end

function build_benchmark_context()
    return build_revolt_like_hil_context(;
        backend_name=PRE_HIL_BACKEND_NAME,
        config_dir=PRE_HIL_CONFIG_DIR,
        sensor=CMOSSensor(T=Float32),
        T=Float32,
        rng=runtime_rng(20260718),
    )
end

function correctness_evidence(rtol::Float64, atol::Float64)
    cpu_ctx = build_revolt_like_hil_context(;
        backend_name="cpu",
        config_dir=PRE_HIL_CONFIG_DIR,
        sensor=CMOSSensor(T=Float32),
        T=Float32,
        rng=runtime_rng(20260718),
    )
    target_ctx = build_benchmark_context()
    revolt_like_step!(cpu_ctx)
    revolt_like_step!(target_ctx)
    reference = cpu_ctx.tiled_frame
    observed = Array(target_ctx.tiled_frame)
    size(observed) == (352, 352) || error(
        "unexpected REVOLT-like output size $(size(observed))")
    all(isfinite, observed) || error("backend output contains non-finite values")
    isapprox(observed, reference; rtol=rtol, atol=atol) || error(
        "$PRE_HIL_BACKEND_NAME output failed CPU parity")
    difference = Float64.(observed) .- Float64.(reference)
    reference64 = Float64.(reference)
    reference_norm = norm(reference64)
    reference_peak = maximum(abs, reference64)
    return Dict{String,Any}(
        "passed" => true,
        "rtol" => rtol,
        "atol" => atol,
        "cpu_checksum" => sum(reference64),
        "backend_checksum" => sum(Float64, observed),
        "relative_l2_error" => reference_norm == 0 ? 0.0 :
            norm(difference) / reference_norm,
        "max_error_relative_to_cpu_peak" => reference_peak == 0 ? 0.0 :
            maximum(abs, difference) / reference_peak,
        "output_shape" => collect(size(observed)),
    )
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
        "encoded backend histogram lost samples")
    HdrHistogram.value_at_percentile(decoded, 95.0) ==
        HdrHistogram.value_at_percentile(histogram, 95.0) || error(
        "encoded backend histogram changed p95")
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
        "p95_ns" => HdrHistogram.value_at_percentile(histogram, 95.0),
        "max_ns" => max(histogram),
        "gc" => gc_counters,
        "histogram_encoding" => "HdrHistogram compressed V2 base64",
        "histogram_base64" => base64encode(encoded),
    )
end

function measure_boundary_run!(boundary::PreHILBoundary, samples::Int,
    lowest_ns::Int64, highest_ns::Int64, significant_figures::Int)
    histogram = HdrHistogram.Histogram(lowest_ns, highest_ns,
        significant_figures)
    GC.gc()
    gc_before = Base.gc_num()
    wall_start = time_ns()
    @inbounds for _ in 1:samples
        start = time_ns()
        run_pre_hil_boundary!(boundary)
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
        "median_p95_ns" => median_integer(run["p95_ns"] for run in runs),
        "worst_p95_ns" => maximum(run["p95_ns"] for run in runs),
        "maximum_observed_ns" => maximum(run["max_ns"] for run in runs),
        "median_mean_ns" => median(collect(run["mean_ns"] for run in runs)),
        "median_throughput_hz" => median(collect(
            run["throughput_hz"] for run in runs)),
    )
end

function boundary_regression(config::AbstractDict, summary,
    steady_alloc_bytes::Int64, p95_gate_supported::Bool,
    relative_factor::Float64)
    absolute_observed = Int64(summary["worst_p95_ns"])
    relative_observed = Int64(summary["median_p95_ns"])
    absolute_limit = Int64(config["absolute_p95_ns"])
    allocation_limit = Int64(config["max_alloc_bytes"])
    relative_requested = Bool(config["relative_gate"])
    baseline = Int64(config["baseline_median_p95_ns"])
    relative_evaluated = p95_gate_supported && relative_requested
    relative_limit = relative_evaluated ?
        ceil(Int64, relative_factor * baseline) : Int64(0)
    return Dict{String,Any}(
        "allocation_gate_passed" => steady_alloc_bytes <= allocation_limit,
        "allocation_observed_bytes" => steady_alloc_bytes,
        "allocation_limit_bytes" => allocation_limit,
        "absolute_p95_gate_evaluated" => p95_gate_supported,
        "absolute_p95_gate_passed" => !p95_gate_supported ||
            absolute_observed <= absolute_limit,
        "absolute_p95_observed_ns" => absolute_observed,
        "absolute_p95_limit_ns" => absolute_limit,
        "relative_p95_gate_evaluated" => relative_evaluated,
        "relative_p95_gate_passed" => !relative_evaluated ||
            relative_observed <= relative_limit,
        "relative_p95_observed_ns" => relative_observed,
        "baseline_median_p95_ns" => baseline,
        "relative_p95_limit_ns" => relative_limit,
        "relative_p95_factor" => relative_factor,
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
        "p95_clock_delta_ns" => HdrHistogram.value_at_percentile(
            histogram, 95.0),
        "wall_ns_per_clock_pair_and_record" => wall_ns / samples,
    )
end

function make_boundaries(ctx::RevoltLikeHILContext,
    boundary_ids::Vector{String})
    host_frame = Matrix{Float32}(undef, size(ctx.tiled_frame))
    boundaries = Dict{String,PreHILBoundary}()
    for id in boundary_ids
        boundaries[id] = if id == "backend_ready"
            let ctx=ctx
                PreHILBoundary(id, () -> revolt_like_step!(ctx))
            end
        elseif id == "host_ready"
            PRE_HIL_BACKEND_NAME == "cpu" && error(
                "host_ready is redundant for the CPU placement")
            let ctx=ctx, host_frame=host_frame
                PreHILBoundary(id, () -> begin
                    revolt_like_step!(ctx)
                    copyto!(host_frame, ctx.tiled_frame)
                    nothing
                end)
            end
        elseif id == "transfer_only"
            PRE_HIL_BACKEND_NAME == "cpu" && error(
                "transfer_only is redundant for the CPU placement")
            let ctx=ctx, host_frame=host_frame
                PreHILBoundary(id, () -> begin
                    copyto!(host_frame, ctx.tiled_frame)
                    nothing
                end)
            end
        else
            error("unknown backend latency boundary '$id'")
        end
    end
    return boundaries
end

function write_artifact(path::AbstractString, artifact)
    isempty(path) && return nothing
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, artifact; sorted=true)
    end
    println("wrote pre-HIL backend artifact: ", path)
    return nothing
end

function run_pre_hil_backend_latency()
    configure_pre_hil_benchmark!()
    isempty(PRE_HIL_PLACEMENT) && error(
        "set AOS_PRE_HIL_PLACEMENT or pass the placement as the second argument")
    contract = TOML.parsefile(PRE_HIL_CONTRACT_PATH)
    placements = contract["placements"]
    haskey(placements, PRE_HIL_PLACEMENT) || error(
        "unknown placement '$PRE_HIL_PLACEMENT'")
    placement = placements[PRE_HIL_PLACEMENT]
    String(placement["backend"]) == PRE_HIL_BACKEND_NAME || error(
        "placement '$PRE_HIL_PLACEMENT' requires backend $(placement["backend"]), not $PRE_HIL_BACKEND_NAME")
    samples = parse(Int, get(ENV, "AOS_PRE_HIL_SAMPLES",
        string(contract["samples_per_run"])))
    run_count = parse(Int, get(ENV, "AOS_PRE_HIL_RUNS",
        string(contract["runs"])))
    warmup = parse(Int, get(ENV, "AOS_PRE_HIL_WARMUP",
        string(contract["warmup_operations"])))
    samples > 0 || error("AOS_PRE_HIL_SAMPLES must be > 0")
    run_count > 0 || error("AOS_PRE_HIL_RUNS must be > 0")
    warmup >= 0 || error("AOS_PRE_HIL_WARMUP must be >= 0")
    minimum_p95_samples = Int(contract["minimum_samples_for_p95_gate"])
    p95_gate_supported = samples >= minimum_p95_samples
    lowest_ns = Int64(contract["histogram_lowest_ns"])
    highest_ns = Int64(contract["histogram_highest_ns"])
    significant_figures = Int(contract["histogram_significant_figures"])
    relative_factor = Float64(contract["relative_p95_factor"])
    workload = contract["workload"]
    rtol = Float64(workload["parity_rtol"])
    atol = Float64(workload["parity_atol"])

    correctness = correctness_evidence(rtol, atol)
    ctx = build_benchmark_context()
    residency = assert_backend_residency(ctx)
    boundary_ids = String.(placement["boundaries"])
    boundaries = make_boundaries(ctx, boundary_ids)
    results = Vector{Dict{String,Any}}()
    all_gates_passed = true

    println("pre_hil_backend_latency")
    println("  placement: ", PRE_HIL_PLACEMENT)
    println("  backend: ", PRE_HIL_BACKEND_NAME)
    println("  samples_per_run: ", samples)
    println("  runs: ", run_count)
    println("  warmup_operations: ", warmup)
    println("  boundaries: ", join(boundary_ids, ','))
    if !p95_gate_supported
        println("  p95 gates: skipped (requires at least ",
            minimum_p95_samples, " samples per run)")
    end

    for id in boundary_ids
        boundary = boundaries[id]
        first_start = time_ns()
        run_pre_hil_boundary!(boundary)
        first_use_ns = Int64(time_ns() - first_start)
        for _ in 1:warmup
            run_pre_hil_boundary!(boundary)
        end
        GC.gc()
        steady_alloc_bytes = Int64(@allocated run_pre_hil_boundary!(boundary))
        runs = Vector{Dict{String,Any}}(undef, run_count)
        for run_index in 1:run_count
            runs[run_index] = measure_boundary_run!(boundary, samples,
                lowest_ns, highest_ns, significant_figures)
            println(id, " run=", run_index,
                " p50_ns=", runs[run_index]["p50_ns"],
                " p95_ns=", runs[run_index]["p95_ns"],
                " max_ns=", runs[run_index]["max_ns"])
        end
        summary = summarize_runs(runs)
        regression = boundary_regression(placement[id], summary,
            steady_alloc_bytes, p95_gate_supported, relative_factor)
        passed = regression["allocation_gate_passed"] &&
            regression["absolute_p95_gate_passed"] &&
            regression["relative_p95_gate_passed"]
        all_gates_passed &= passed
        push!(results, Dict{String,Any}(
            "id" => id,
            "first_use_ns" => first_use_ns,
            "steady_alloc_bytes" => steady_alloc_bytes,
            "runs" => runs,
            "summary" => summary,
            "regression" => regression,
            "passed" => passed,
        ))
    end

    environment = source_environment()
    if !isempty(PRE_HIL_OUTPUT_PATH) && environment["git_dirty"]
        error("refusing to write durable pre-HIL backend evidence from a dirty worktree")
    end
    artifact = Dict{String,Any}(
        "schema_version" => Int(contract["schema_version"]),
        "benchmark" => String(contract["name"]),
        "evidence_class" => String(contract["evidence_class"]),
        "source_contract" => relpath(abspath(PRE_HIL_CONTRACT_PATH), PRE_HIL_ROOT),
        "baseline_artifact" => String(contract["baseline_artifact"]),
        "characterized_source_revision" => environment["git_commit"],
        "placement" => PRE_HIL_PLACEMENT,
        "backend" => PRE_HIL_BACKEND_NAME,
        "host_class" => String(placement["host_class"]),
        "configured_samples_per_run" => samples,
        "configured_runs" => run_count,
        "configured_warmup_operations" => warmup,
        "minimum_samples_for_p95_gate" => minimum_p95_samples,
        "p95_gate_supported" => p95_gate_supported,
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
        "residency" => residency,
        "boundaries" => results,
        "all_gates_passed" => all_gates_passed,
    )
    write_artifact(PRE_HIL_OUTPUT_PATH, artifact)
    all_gates_passed || error(
        "one or more pre-HIL backend latency gates failed")
    return artifact
end

run_pre_hil_backend_latency()

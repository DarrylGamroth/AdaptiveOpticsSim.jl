const REQUESTED_BACKEND = lowercase(strip(get(
    ENV,
    "AOS_REVOLT_GRAPH_BACKEND",
    "cpu",
)))

if REQUESTED_BACKEND == "cuda"
    @eval using CUDA
elseif REQUESTED_BACKEND == "amdgpu"
    @eval using AMDGPU
elseif REQUESTED_BACKEND != "cpu"
    error("AOS_REVOLT_GRAPH_BACKEND must be cpu, cuda, or amdgpu")
end

using AdaptiveOpticsSim
using AdaptiveOpticsSim.AlgorithmGraphs
using AdaptiveOpticsSim.Backends
using HdrHistogram
using LinearAlgebra
using TOML

include(joinpath(
    dirname(@__DIR__),
    "examples",
    "support",
    "revolt_hsdm277.jl",
))
include(joinpath(
    dirname(@__DIR__),
    "examples",
    "support",
    "revolt_hil_graphs.jl",
))

const AOG = AdaptiveOpticsSim.AlgorithmGraphs
const REPOSITORY_ROOT = dirname(@__DIR__)
const PROFILE_SAMPLES = parse(Int, get(
    ENV,
    "AOS_REVOLT_GRAPH_SAMPLES",
    "50",
))
const PROFILE_WARMUP = parse(Int, get(
    ENV,
    "AOS_REVOLT_GRAPH_WARMUP",
    "5",
))
const PROFILE_GRAPH_PROFILE = Symbol(lowercase(strip(get(
    ENV,
    "AOS_REVOLT_GRAPH_PROFILE",
    "full_optical",
))))
const PROFILE_EXECUTION = Symbol(lowercase(strip(get(
    ENV,
    "AOS_REVOLT_GRAPH_EXECUTION",
    "stream",
))))
const PROFILE_ARCHITECTURES = Tuple(Symbol(strip(value)) for value in split(
    get(ENV, "AOS_REVOLT_GRAPH_ARCHITECTURES", "classic,copper"),
    ',';
    keepempty=false,
))
const HISTOGRAM_HIGHEST_NS = Int64(60_000_000_000)
const HISTOGRAM_SIGNIFICANT_FIGURES = 3

struct RevoltGraphProfile
    name::Symbol
    graph_file::String
    command_index::Int
    frame_output::Symbol
end

mutable struct LatencyRecorder{H}
    histogram::H
    total_ns::Int128
    samples::Int
end

function latency_recorder()
    histogram = HdrHistogram.Histogram(
        Int64(1),
        HISTOGRAM_HIGHEST_NS,
        HISTOGRAM_SIGNIFICANT_FIGURES,
    )
    return LatencyRecorder(histogram, Int128(0), 0)
end

@inline function record_latency!(recorder::LatencyRecorder, elapsed_ns)
    elapsed = max(Int64(1), Int64(elapsed_ns))
    HdrHistogram.record_value!(recorder.histogram, elapsed)
    recorder.total_ns += elapsed
    recorder.samples += 1
    return nothing
end

function latency_summary(recorder::LatencyRecorder)
    recorder.samples > 0 || error("cannot summarize an empty latency recorder")
    histogram = recorder.histogram
    mean_ns = Float64(recorder.total_ns) / recorder.samples
    return Dict{String,Any}(
        "samples" => recorder.samples,
        "minimum_ns" => min(histogram),
        "mean_ns" => mean_ns,
        "p50_ns" => HdrHistogram.value_at_percentile(histogram, 50.0),
        "p90_ns" => HdrHistogram.value_at_percentile(histogram, 90.0),
        "p95_ns" => HdrHistogram.value_at_percentile(histogram, 95.0),
        "p99_ns" => HdrHistogram.value_at_percentile(histogram, 99.0),
        "maximum_ns" => max(histogram),
        "throughput_from_mean_hz" => 1.0e9 / mean_ns,
    )
end

function revolt_graph_profile(name::Symbol)
    if name === :classic
        return RevoltGraphProfile(
            name,
            REVOLTHILGraphs.graph_path(name, PROFILE_GRAPH_PROFILE),
            143,
            :shwfs_frame,
        )
    elseif name === :copper
        return RevoltGraphProfile(
            name,
            REVOLTHILGraphs.graph_path(name, PROFILE_GRAPH_PROFILE),
            139,
            :pwfs_frame,
        )
    end
    error("unknown REVOLT graph architecture '$name'")
end

function configure_profile!()
    Threads.nthreads() == 1 || error(
        "REVOLT graph latency profiling requires one Julia thread",
    )
    BLAS.set_num_threads(1)
    Backends.set_fft_provider_threads!(1)
    PROFILE_SAMPLES > 0 || error("AOS_REVOLT_GRAPH_SAMPLES must be positive")
    PROFILE_WARMUP >= 0 || error(
        "AOS_REVOLT_GRAPH_WARMUP must be nonnegative",
    )
    isempty(PROFILE_ARCHITECTURES) && error(
        "AOS_REVOLT_GRAPH_ARCHITECTURES must select at least one graph",
    )
    PROFILE_GRAPH_PROFILE in REVOLTHILGraphs.supported_profiles() || error(
        "AOS_REVOLT_GRAPH_PROFILE must select one of " *
        join(REVOLTHILGraphs.supported_profiles(), ", "),
    )
    PROFILE_EXECUTION in (:stream, :captured) || error(
        "AOS_REVOLT_GRAPH_EXECUTION must be stream or captured",
    )
    return nothing
end

@inline requested_graph_execution() = PROFILE_EXECUTION === :captured ?
    CapturedGraphExecution() : StreamGraphExecution()

function requested_backend()
    if REQUESTED_BACKEND == "cpu"
        return (
            name="cpu",
            array_type=Array,
            package="none",
            package_version="none",
        )
    end

    tag = if REQUESTED_BACKEND == "cuda"
        Backends.CUDABackendTag
    else
        Backends.AMDGPUBackendTag
    end
    package_name = if REQUESTED_BACKEND == "cuda"
        "CUDA"
    else
        "AMDGPU"
    end
    package_module = getfield(Main, Symbol(package_name))
    Base.invokelatest(getproperty(package_module, :functional)) || error(
        "$package_name is not functional on this host",
    )
    Backends.disable_scalar_backend!(tag)
    array_type = Backends.gpu_backend_array_type(tag)
    isnothing(array_type) && error(
        "$package_name did not register a GPU array type",
    )
    return (
        name=REQUESTED_BACKEND,
        array_type,
        package=package_name,
        package_version=string(Base.pkgversion(package_module)),
    )
end

@inline function backend_array(::Type{A}, values) where {A}
    return A(values)
end

function prepare_profile_graph(
    profile::RevoltGraphProfile,
    backend,
)
    command_host = zeros(Float32, 277)
    command_host[profile.command_index] = 5.0f-8
    coordinates_host = REVOLTHSDM277.actuator_coordinates(Float32)
    grid_indices_host = REVOLTHSDM277.actuator_grid_indices(Int32)
    command = backend_array(backend.array_type, command_host)
    coordinates = backend_array(backend.array_type, coordinates_host)
    grid_indices = backend_array(backend.array_type, grid_indices_host)
    target = compute_device(command)
    definition = load_algorithm_graph(
        profile.graph_file;
        bindings=(;
            pdm_command=command,
            pdm_actuator_coordinates=coordinates,
            pdm_actuator_grid_indices=grid_indices,
        ),
    )
    start = time_ns()
    graph = prepare_algorithm_graph(
        definition;
        target,
        execution=requested_graph_execution(),
    )
    preparation_ns = Int64(time_ns() - start)
    return (; graph, target, preparation_ns)
end

function synchronized_node_step!(graph, node)
    AOG._with_prepared_device_execution_context(graph.context) do
        try
            AOG.step_graph_node!(node.owner)
        finally
            AOG._synchronize_prepared_device_execution_context!(graph.context)
        end
    end
    return nothing
end

@inline warm_node_cycle!(graph, ::Tuple{}) = nothing

@inline function warm_node_cycle!(graph, nodes::Tuple)
    synchronized_node_step!(graph, first(nodes))
    warm_node_cycle!(graph, Base.tail(nodes))
    return nothing
end

@inline record_node_cycle!(graph, ::Tuple{}, ::Tuple{}) = nothing

@inline function record_node_cycle!(
    graph,
    nodes::Tuple,
    recorders::Tuple,
)
    start = time_ns()
    synchronized_node_step!(graph, first(nodes))
    elapsed = time_ns() - start
    record_latency!(first(recorders), elapsed)
    record_node_cycle!(
        graph,
        Base.tail(nodes),
        Base.tail(recorders),
    )
    return nothing
end

function node_results(names::Tuple, recorders::Tuple)
    results = Vector{Dict{String,Any}}(undef, length(names))
    for index in eachindex(names, recorders)
        result = latency_summary(recorders[index])
        result["name"] = String(names[index])
        results[index] = result
    end
    return results
end

function profile_graph_execution(profile::RevoltGraphProfile, backend)
    prepared = prepare_profile_graph(profile, backend)
    graph = prepared.graph
    names = keys(graph.nodes)
    nodes = values(graph.nodes)
    recorders = ntuple(_ -> latency_recorder(), length(nodes))

    for _ in 1:PROFILE_WARMUP
        warm_node_cycle!(graph, nodes)
    end
    for _ in 1:PROFILE_SAMPLES
        record_node_cycle!(graph, nodes, recorders)
    end
    per_node = node_results(names, recorders)

    reset_graph!(graph)
    for _ in 1:PROFILE_WARMUP
        step_graph!(graph)
    end
    atmosphere_before = Array(graph_output(graph, Val(:atmosphere_opd)))
    submission_recorder = latency_recorder()
    completion_wait_recorder = latency_recorder()
    target_ready_recorder = latency_recorder()
    for _ in 1:PROFILE_SAMPLES
        start = time_ns()
        ticket = step_graph_async!(graph)
        submitted = time_ns()
        wait_graph_step!(ticket)
        ready = time_ns()
        record_latency!(submission_recorder, submitted - start)
        record_latency!(completion_wait_recorder, ready - submitted)
        record_latency!(target_ready_recorder, ready - start)
    end
    atmosphere_after = Array(graph_output(graph, Val(:atmosphere_opd)))
    frame = Array(graph_output(graph, profile.frame_output))
    all(isfinite, frame) || error("$(profile.name) graph produced non-finite frame values")
    sum(frame) > 0 || error("$(profile.name) graph produced an empty frame")
    atmosphere_after != atmosphere_before || error(
        "$(profile.name) atmosphere did not evolve during profiling",
    )

    submission_summary = latency_summary(submission_recorder)
    completion_wait_summary = latency_summary(completion_wait_recorder)
    target_ready_summary = latency_summary(target_ready_recorder)
    node_mean_sum = sum(result["mean_ns"] for result in per_node)
    return Dict{String,Any}(
        "preparation_ns" => prepared.preparation_ns,
        "target" => string(prepared.target),
        "execution_policy" => String(PROFILE_EXECUTION),
        "captured_graph_nodes" => captured_graph_node_count(graph),
        "nodes" => per_node,
        "node_mean_sum_ns" => node_mean_sum,
        "graph_submission" => submission_summary,
        "graph_completion_wait" => completion_wait_summary,
        "graph_target_ready" => target_ready_summary,
        "separately_synchronized_node_sum_over_graph_mean" =>
            node_mean_sum / target_ready_summary["mean_ns"],
        "frame_shape" => collect(size(frame)),
        "atmosphere_evolved" => true,
        "frame_finite" => true,
        "frame_nonzero" => true,
    )
end

function profile_hil_execution(profile::RevoltGraphProfile, backend)
    prepared = prepare_profile_graph(profile, backend)
    graph = prepared.graph
    boundary = prepare_graph_hil_boundary(
        graph;
        command_input=:pdm_command,
        frame_output=profile.frame_output,
    )
    for _ in 1:PROFILE_WARMUP
        sequence = step_hil_frame!(boundary)
        adopt_hil_command!(boundary, sequence)
    end

    frame_recorder = latency_recorder()
    command_recorder = latency_recorder()
    cycle_recorder = latency_recorder()
    for _ in 1:PROFILE_SAMPLES
        cycle_start = time_ns()
        frame_start = cycle_start
        sequence = step_hil_frame!(boundary)
        frame_end = time_ns()
        adopt_hil_command!(boundary, sequence)
        cycle_end = time_ns()
        record_latency!(frame_recorder, frame_end - frame_start)
        record_latency!(command_recorder, cycle_end - frame_end)
        record_latency!(cycle_recorder, cycle_end - cycle_start)
    end

    frame = hil_frame_buffer(boundary)
    all(isfinite, frame) || error(
        "$(profile.name) HIL boundary produced non-finite frame values",
    )
    sum(frame) > 0 || error(
        "$(profile.name) HIL boundary produced an empty frame",
    )
    return Dict{String,Any}(
        "preparation_ns" => prepared.preparation_ns,
        "target" => string(prepared.target),
        "execution_policy" => String(PROFILE_EXECUTION),
        "captured_graph_nodes" => captured_graph_node_count(graph),
        "frame_host_ready" => latency_summary(frame_recorder),
        "command_target_ready" => latency_summary(command_recorder),
        "lockstep_cycle" => latency_summary(cycle_recorder),
        "frame_shape" => collect(size(frame)),
        "frame_finite" => true,
        "frame_nonzero" => true,
    )
end

function git_revision()
    try
        return readchomp(Cmd(`git rev-parse HEAD`; dir=REPOSITORY_ROOT))
    catch
        return "unknown"
    end
end

function tracked_source_dirty()
    try
        return !success(Cmd(`git diff --quiet HEAD --`; dir=REPOSITORY_ROOT))
    catch
        return true
    end
end

function profile_environment(backend, target)
    cpu_model = isempty(Sys.cpu_info()) ? "unknown" : Sys.cpu_info()[1].model
    return Dict{String,Any}(
        "julia_version" => string(VERSION),
        "kernel" => string(Sys.KERNEL),
        "architecture" => string(Sys.ARCH),
        "cpu_model" => cpu_model,
        "logical_cpu_threads" => Sys.CPU_THREADS,
        "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
        "fft_threads" => 1,
        "backend" => backend.name,
        "backend_package" => backend.package,
        "backend_package_version" => backend.package_version,
        "compute_device" => string(target),
    )
end

function main()
    configure_profile!()
    backend = requested_backend()
    profiles = Tuple(revolt_graph_profile(name) for name in PROFILE_ARCHITECTURES)
    probe = backend_array(backend.array_type, zeros(Float32, 1))
    target = compute_device(probe)
    results = Vector{Dict{String,Any}}(undef, length(profiles))

    for index in eachindex(profiles)
        profile = profiles[index]
        graph_execution = profile_graph_execution(profile, backend)
        GC.gc(true)
        hil_execution = profile_hil_execution(profile, backend)
        GC.gc(true)
        results[index] = Dict{String,Any}(
            "name" => String(profile.name),
            "graph_profile" => String(PROFILE_GRAPH_PROFILE),
            "execution_policy" => String(PROFILE_EXECUTION),
            "graph_file" => relpath(profile.graph_file, REPOSITORY_ROOT),
            "graph_execution" => graph_execution,
            "hil_execution" => hil_execution,
        )
    end

    artifact = Dict{String,Any}(
        "artifact_id" => "REVOLT-ALGORITHM-GRAPH-NODE-PROFILE",
        "evidence_class" =>
            "warmed submission and synchronized-completion diagnostic profile",
        "source_revision" => git_revision(),
        "source_tracked_dirty" => tracked_source_dirty(),
        "contract" => Dict{String,Any}(
            "graph_profile" => String(PROFILE_GRAPH_PROFILE),
            "samples" => PROFILE_SAMPLES,
            "warmup_cycles" => PROFILE_WARMUP,
            "timer" => "time_ns; overhead included and not subtracted",
            "load_model" => "warmed serial closed loop",
            "arrival_model" => "closed loop; next cycle starts after completion",
            "coordinated_omission_correction" => false,
            "per_node_boundary" =>
                "one node invocation plus exact graph-context synchronization",
            "graph_submission_boundary" =>
                "one bounded-capacity graph submission with no graph-boundary synchronization",
            "graph_completion_wait_boundary" =>
                "completion-ticket wait through exact graph-context synchronization",
            "graph_boundary" =>
                "submission immediately followed by completion-ticket wait",
            "frame_host_ready_boundary" =>
                "graph boundary plus completed target-to-host frame copy",
            "command_target_ready_boundary" =>
                "complete finite host command validation and completed host-to-target copy",
            "lockstep_cycle_boundary" =>
                "frame-host-ready boundary immediately followed by " *
                "same-sequence command adoption",
            "per_node_sum_caveat" =>
                "each node is synchronized separately, so its sum includes " *
                "more synchronization than graph execution",
            "tail_note" =>
                "percentiles are diagnostic; interpret tail levels only when " *
                "sample count supplies enough observations",
        ),
        "environment" => profile_environment(backend, target),
        "architectures" => results,
    )
    TOML.print(stdout, artifact; sorted=true)
    println()
    return artifact
end

main()

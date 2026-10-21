# Run in an isolated environment containing CUDA, local AOS, and Copper:
# julia --project=/tmp/copper-profile-env --startup-file=no \
#   benchmarks/profile_copper_cuda.jl /path/to/REVOLTCopperSim.jl output.toml 128
# AOS_COPPER_PROFILE=1 adds separate, untimed CUDA profiler passes.
# No environment activation, package installation, or source overrides occur here.
using AdaptiveOpticsSim, AdaptiveOpticsSim.AlgorithmGraphs
using REVOLTCopperSim
using LinearAlgebra, Statistics, TOML, SHA
import CUDA

const AG = AdaptiveOpticsSim.AlgorithmGraphs
const BACKENDS = AdaptiveOpticsSim.Backends

function source_record(mod)
    root = dirname(dirname(pathof(mod)))
    return Dict("root" => root, "version" => string(pkgversion(mod)),
        "revision" => readchomp(`git -C $root rev-parse HEAD`),
        "status" => readchomp(`git -C $root status --porcelain`))
end

function distribution(samples)
    return Dict("samples_ms" => samples, "minimum_ms" => minimum(samples),
        "median_ms" => median(samples), "p95_ms" => quantile(samples, 0.95),
        "maximum_ms" => maximum(samples))
end

function complete_samples!(graph, repeats, warmups)
    for _ in 1:warmups
        step_graph!(graph)
    end
    samples = Vector{Float64}(undef, repeats)
    for index in eachindex(samples)
        started = time_ns()
        step_graph!(graph)
        samples[index] = (time_ns() - started) / 1e6
    end
    result = distribution(samples)
    result["host_bytes"] = [@allocated(step_graph!(graph)) for _ in 1:16]
    return result
end

# Diagnostic traversal only: use the actual prepared context and node protocol.
# Events exclude the completion wait; host samples include enqueue and completion.
# Synchronizing each node perturbs scheduling, so these are attribution costs,
# not additive estimates of uninterrupted graph latency. Private ownership access
# is deliberate and confined to this benchmark, not exposed as package API.
function node_samples!(graph, repeats, warmups, captured::Bool)
    isempty(graph.delayed_links) || error("node attribution requires no delayed links")
    for _ in 1:warmups
        step_graph!(graph)
    end
    names = propertynames(graph.nodes)
    gpu_samples = [Vector{Float64}(undef, repeats) for _ in names]
    host_samples = [Vector{Float64}(undef, repeats) for _ in names]
    start_event, stop_event = CUDA.CuEvent(), CUDA.CuEvent()
    for repeat in 1:repeats
        # Match the capture protocol's frame preflight/publication ordering.
        # Preparation itself exercises captured enqueue outside native capture
        # (captured_execution.jl), so this is a supported diagnostic use of
        # those exact owners. No native graph is launched during attribution.
        if captured
            for node in values(graph.nodes)
                AG._preflight_captured_graph_node!(node.owner)
            end
        end
        for (index, name) in enumerate(names)
            owner = AG.prepared_graph_node(graph, name)
            started = time_ns()
            BACKENDS._with_prepared_device_execution_context(graph.context) do
                CUDA.record(start_event)
                if captured
                    AG.enqueue_captured_graph_node!(owner)
                else
                    AG.enqueue_graph_node!(owner)
                end
                CUDA.record(stop_event)
                CUDA.synchronize(stop_event)
            end
            host_samples[index][repeat] = (time_ns() - started) / 1e6
            gpu_samples[index][repeat] = 1000 * CUDA.elapsed(start_event, stop_event)
        end
        if captured
            for node in values(graph.nodes)
                AG._complete_captured_graph_node!(node.owner)
            end
        end
    end
    return Dict(string(name) => Dict(
        "device_event" => distribution(gpu_samples[index]),
        "host_complete" => distribution(host_samples[index]),
    ) for (index, name) in enumerate(names))
end

function main(args)
    2 <= length(args) <= 3 || error(
        "usage: profile_copper_cuda.jl COPPER_ROOT OUTPUT_TOML [REPEATS=128]")
    copper_root = realpath(args[1])
    realpath(dirname(dirname(pathof(REVOLTCopperSim)))) == copper_root ||
        error("loaded Copper source does not match COPPER_ROOT")
    realpath(dirname(dirname(pathof(AdaptiveOpticsSim)))) == realpath(joinpath(@__DIR__, "..")) ||
        error("isolated environment must load this AOS checkout")
    repeats = length(args) == 3 ? parse(Int, args[3]) : 128
    repeats > 0 || error("REPEATS must be positive")
    warmups = 32
    CUDA.functional() || error("functional CUDA device required")
    CUDA.allowscalar(false)
    BLAS.set_num_threads(1)
    BACKENDS.set_fft_provider_threads!(1)
    target = BACKENDS.compute_device(CUDA.zeros(Float32, 1))
    result = Dict{String,Any}(
        "julia" => string(VERSION), "cuda" => string(pkgversion(CUDA)),
        "device" => CUDA.name(CUDA.device()), "repeats" => repeats,
        "warmups" => warmups, "aos" => source_record(AdaptiveOpticsSim),
        "copper" => source_record(REVOLTCopperSim),
        "harness_sha256" => bytes2hex(sha256(read(@__FILE__))),
        "active_project" => something(Base.active_project(), ""),
        "script" => abspath(@__FILE__), "arguments" => copy(args),
        "workload" => "default noisy five-layer Copper grid_gaussian; zero commands; evolving atmosphere/RNG",
        "total_boundary" => "host elapsed time around step_graph! including prepared graph completion; no HIL staging or transport",
        "node_boundary" => "stream or captured-protocol node enqueue in prepared context, CUDA events and per-node completion; serialized attribution, excludes graph bookkeeping and frame preflight/publication",
        "notes" => "Node traversal advances owner state but bypasses graph sequence bookkeeping; separate disposable stream/captured graphs. Captured attribution invokes capture-specific enqueue after regular capture preparation, with frame preflight/publication, but does not measure native replay scheduling. Node costs are not additive full-frame costs. Timings do not establish optical accuracy.",
    )
    for (label, execution) in (("stream", StreamGraphExecution()),
            ("captured", CapturedGraphExecution()))
        system = prepare_hil_system(; target, execution)
        result[label] = complete_samples!(system.graph, repeats, warmups)
        if get(ENV, "AOS_COPPER_PROFILE", "0") == "1"
            println("Untimed CUDA profile: ", label)
            profile = CUDA.@profile begin
                for _ in 1:4
                    step_graph!(system.graph)
                end
            end
            show(stdout, profile)
            println()
        end
        system = nothing
        GC.gc(true)
        CUDA.reclaim()
    end
    for (label, execution, captured) in (
            ("nodes_stream", StreamGraphExecution(), false),
            ("nodes_captured", CapturedGraphExecution(), true))
        system = prepare_hil_system(; target, execution)
        result[label] = node_samples!(system.graph, repeats, warmups, captured)
        system = nothing
        GC.gc(true)
        CUDA.reclaim()
    end
    mkpath(dirname(abspath(args[2])))
    open(args[2], "w") do io
        TOML.print(io, result; sorted=true)
    end
    println("Wrote ", abspath(args[2]))
    for label in ("stream", "captured")
        println(label, " complete median: ", result[label]["median_ms"], " ms")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end

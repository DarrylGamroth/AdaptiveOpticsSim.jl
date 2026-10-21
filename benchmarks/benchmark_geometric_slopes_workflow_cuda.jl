# Complete Gaussian-DM -> pupil-OPD -> geometric wavefront slope workflow.
# Compare the retained serial GPU kernel with the selected implementation using
# device events, and measure the public synchronous workflow separately.
# AOS_BENCHMARK_ROOT permits running this same driver against another checkout.
import CUDA
using TOML, Statistics, LinearAlgebra
const _GEOMETRIC_BENCH_ROOT = get(ENV, "AOS_BENCHMARK_ROOT",
    normpath(joinpath(@__DIR__, "..")))
include(joinpath(_GEOMETRIC_BENCH_ROOT, "benchmarks", "profile_aos_gpu_plant.jl"))
const _GEOMETRIC_HAS_COOPERATIVE = isdefined(WavefrontSensors,
    :cooperative_geometric_slopes_kernel!)

function _geometric_fields(io, fields::NamedTuple)
    for (key, value) in pairs(fields)
        TOML.print(io, Dict{String,typeof(value)}(String(key)=>value))
    end
end

function _geometric_definition(command, opd, coordinates, n, lenslets)
    complete = _plant_definition(command, opd, coordinates;
        name=:geometric_workflow, pupil_resolution=n, lenslets_per_axis=lenslets)
    return algorithm_graph(complete.nodes[1:2]; name=:geometric_workflow,
        inputs=complete.inputs, outputs=complete.outputs[1:2],
        links=(first(complete.links),), parameters=complete.parameters)
end

function _geometric_workflow!(graph, slopes, mask, sampling)
    _complete_step!(graph)
    WavefrontSensors.geometric_wavefront_slopes!(slopes,
        graph_output(graph, Val(:pupil_opd)), mask, sampling)
    return nothing
end

function _geometric_kernel!(style, args, n_sub, ::Val{:serial})
    Backends.launch_kernel_async!(style, WavefrontSensors.geometric_slopes_kernel!,
        args...; ndrange=(n_sub, n_sub))
    return nothing
end
function _geometric_kernel!(style, args, n_sub, ::Val{:selected})
    if _GEOMETRIC_HAS_COOPERATIVE
        WavefrontSensors._enqueue_geometric_slopes!(style, args...)
    else
        _geometric_kernel!(style, args, n_sub, Val(:serial))
    end
    return nothing
end
function _geometric_kernel_graph(style, args, n_sub, variant)
    _geometric_kernel!(style, args, n_sub, variant)
    CUDA.synchronize()
    graph = CUDA.capture() do
        for _ in 1:32
            _geometric_kernel!(style, args, n_sub, variant)
        end
    end
    executable = CUDA.instantiate(graph)
    CUDA.launch(executable)
    CUDA.synchronize()
    return (; graph, executable)
end
function _geometric_event_block(graph, start_event, stop_event)
    times = zeros(Float64, 32)
    for i in eachindex(times)
        CUDA.record(start_event)
        CUDA.launch(graph.executable)
        CUDA.record(stop_event)
        CUDA.synchronize(stop_event)
        times[i] = CUDA.elapsed(start_event, stop_event) / 32
    end
    return times
end

function _geometric_workflow_case(io, n, lenslets)
    command, coordinates = _pdm_command(), _actuator_coordinates()
    opd = zeros(Float32, n, n)
    mask = [mod(i + 2j, 5) != 0 for i in 1:lenslets, j in 1:lenslets]
    sampling = (8f0 / n, 8f0 / n)
    cpu = prepare_algorithm_graph(_geometric_definition(command, opd,
        coordinates, n, lenslets))
    reference = zeros(Float32, 2lenslets^2)
    gpu_command = CUDA.CuArray(command)
    gpu = prepare_algorithm_graph(_geometric_definition(gpu_command,
        CUDA.CuArray(opd), CUDA.CuArray(coordinates), n, lenslets);
        target=Backends.compute_device(gpu_command), execution=CapturedGraphExecution())
    slopes = CUDA.zeros(Float32, length(reference))
    gpu_mask = CUDA.CuArray(mask)
    CUDA.synchronize()
    maximum_error = 0.0
    for scale in (0f0, -0.7f0, 1.3f0, 1f0)
        copyto!(command, scale .* _pdm_command())
        copyto!(gpu_command, command)
        CUDA.synchronize()
        _geometric_workflow!(cpu, reference, mask, sampling)
        _geometric_workflow!(gpu, slopes, gpu_mask, sampling)
        actual = Array(slopes)
        @assert isapprox(actual, reference; rtol=5f-4, atol=2f-12)
        maximum_error = max(maximum_error, maximum(abs, actual-reference))
    end
    for _ in 1:32
        _geometric_workflow!(gpu, slopes, gpu_mask, sampling)
    end
    # Dynamic KA launch allocation is reported, not confused with prepared
    # native graph replay's separate zero-host-allocation contract.
    host_bytes = @allocated _geometric_workflow!(gpu, slopes, gpu_mask, sampling)
    times, gc_bytes, gc_collections = zeros(UInt64, 128), zeros(Int, 128), zeros(Int, 128)
    memory_before = CUDA.used_memory()
    for i in eachindex(times)
        before = Base.gc_num()
        start = time_ns()
        _geometric_workflow!(gpu, slopes, gpu_mask, sampling)
        times[i] = time_ns() - start
        delta = Base.GC_Diff(Base.gc_num(), before)
        gc_bytes[i], gc_collections[i] = delta.allocd, delta.pause
    end
    @assert CUDA.used_memory() == memory_before
    @assert isapprox(Array(slopes), reference; rtol=5f-4, atol=2f-12)
    sub = div(n, lenslets)
    args = (slopes, graph_output(gpu, Val(:pupil_opd)), gpu_mask,
        sub, lenslets, lenslets^2, inv(sampling[1]), inv(sampling[2]))
    style = Backends.execution_style(slopes)
    serial = _geometric_kernel_graph(style, args, lenslets, Val(:serial))
    selected = _geometric_kernel_graph(style, args, lenslets, Val(:selected))
    start_event, stop_event = CUDA.CuEvent(), CUDA.CuEvent()
    blocks = ntuple(16) do index
        variant = mod(index-1, 4) in (0, 3) ? "serial" : "selected"
        graph = variant == "serial" ? serial : selected
        samples = _geometric_event_block(graph, start_event, stop_event)
        return (; round=div(index-1, 4)+1, position=mod(index-1, 4)+1,
            variant, event_seconds_per_kernel=samples)
    end
    @assert isapprox(Array(slopes), reference; rtol=5f-4, atol=2f-12)
    println(io, "\n[[case]]")
    _geometric_fields(io, (; pupil_resolution=n, lenslets_per_axis=lenslets,
        pixels_per_subaperture=sub, maximum_cpu_difference=maximum_error,
        workflow_samples_ns=times, workflow_p50_ns=median(times),
        host_bytes, gc_bytes, gc_collections,
        kernel_boundary="CUDA events around32 repeated kernel launches; no reset or public wait, scheduling amortized per kernel"))
    for block in blocks
        println(io, "\n[[case.blocks]]")
        _geometric_fields(io, block)
    end
    flush(io)
    println("N$n L$lenslets public_p50_ns=$(median(times)) host_bytes=$host_bytes cpu_error=$maximum_error")
    return nothing
end

function benchmark_geometric_workflow()
    CUDA.allowscalar(false)
    BLAS.set_num_threads(1)
    Backends.set_fft_provider_threads!(1)
    output = get(ENV, "AOS_GEOMETRIC_RESULTS", "/tmp/geometric-workflow.toml")
    open(output, "w") do io
        _geometric_fields(io, (; source_revision=readchomp(`git -C $_GEOMETRIC_BENCH_ROOT rev-parse HEAD`),
            selected_cooperative=_GEOMETRIC_HAS_COOPERATIVE,
            device=CUDA.name(CUDA.device()), julia=string(VERSION),
            workflow_boundary="Host time_ns around captured DM+pupil graph completion followed by public geometric_wavefront_slopes! including completion;128 samples,32 warmups"))
        for (n, lenslets) in ((96, 32), (128, 32), (160, 32), (512, 32), (1024, 32))
            _geometric_workflow_case(io, n, lenslets)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    benchmark_geometric_workflow()
end

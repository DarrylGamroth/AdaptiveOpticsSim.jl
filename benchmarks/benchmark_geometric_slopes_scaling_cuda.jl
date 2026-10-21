# Kernel-only geometry sweep using the same serial/selected launch adapters as
# the complete workflow benchmark. This probes output-grid density separately.
include("benchmark_geometric_slopes_workflow_cuda.jl")

function _scaling_graph(style, args, n_sub, variant)
    _geometric_kernel!(style, args, n_sub, variant)
    CUDA.synchronize()
    graph = CUDA.capture() do
        for _ in 1:8
            _geometric_kernel!(style, args, n_sub, variant)
        end
    end
    executable = CUDA.instantiate(graph)
    CUDA.launch(executable)
    CUDA.synchronize()
    return (; graph, executable)
end
function _scaling_events(graph, start_event, stop_event)
    samples = zeros(Float64, 8)
    for i in eachindex(samples)
        CUDA.record(start_event)
        CUDA.launch(graph.executable)
        CUDA.record(stop_event)
        CUDA.synchronize(stop_event)
        samples[i] = CUDA.elapsed(start_event, stop_event) / 8
    end
    return samples
end
function _scaling_case(io, lenslets, sub)
    n = lenslets * sub
    host_opd = [1f-7 * (sin(0.017f0 * i) + cos(0.013f0 * j)) for i in 1:n, j in 1:n]
    mask = [mod(i+2j, 5) != 0 for i in 1:lenslets, j in 1:lenslets]
    reference = zeros(Float32, 2lenslets^2)
    WavefrontSensors.geometric_slopes!(reference, host_opd, mask)
    output, opd, valid = CUDA.zeros(Float32, length(reference)), CUDA.CuArray(host_opd), CUDA.CuArray(mask)
    args = (output, opd, valid, sub, lenslets, lenslets^2, 1f0, 1f0)
    style = Backends.execution_style(output)
    CUDA.synchronize()
    serial = _scaling_graph(style, args, lenslets, Val(:serial))
    selected = _scaling_graph(style, args, lenslets, Val(:selected))
    @assert isapprox(Array(output), reference; rtol=5f-4, atol=2f-12)
    start_event, stop_event = CUDA.CuEvent(), CUDA.CuEvent()
    blocks = ntuple(8) do index
        round = div(index-1, 4)+1
        position = mod(index-1, 4)+1
        isserial = isodd(round) ? position in (1,4) : position in (2,3)
        samples = _scaling_events(isserial ? serial : selected, start_event, stop_event)
        return (; round, position, variant=isserial ? "serial" : "selected",
            event_seconds_per_kernel=samples)
    end
    @assert isapprox(Array(output), reference; rtol=5f-4, atol=2f-12)
    println(io, "\n[[case]]")
    _geometric_fields(io, (; pupil_resolution=n, lenslets_per_axis=lenslets,
        pixels_per_subaperture=sub,
        maximum_cpu_difference=maximum(abs, Array(output)-reference)))
    for block in blocks
        println(io, "\n[[case.blocks]]")
        _geometric_fields(io, block)
    end
    flush(io)
    times(v) = median([t for b in blocks if b.variant == v for t in b.event_seconds_per_kernel])
    serial_us = times("serial") * 1e6
    selected_us = times("selected") * 1e6
    println("L$lenslets sub$sub serial_us=$serial_us selected_us=$selected_us")
    return nothing
end
function benchmark_geometric_scaling()
    CUDA.allowscalar(false)
    BLAS.set_num_threads(1)
    output = get(ENV, "AOS_GEOMETRIC_SCALING_RESULTS", "/tmp/geometric-scaling.toml")
    open(output, "w") do io
        _geometric_fields(io, (; source_revision=readchomp(`git rev-parse HEAD`),
            samples_per_block=8, rounds=2, replays_per_event=8,
            boundary="CUDA events around8 repeated kernels; includes amortized native graph scheduling"))
        for lenslets in (8,16,32,64,128,256), sub in (5,16)
            _scaling_case(io, lenslets, sub)
        end
        for lenslets in (32,128)
            _scaling_case(io, lenslets, 32)
        end
    end
end
if abspath(PROGRAM_FILE) == @__FILE__
    benchmark_geometric_scaling()
end

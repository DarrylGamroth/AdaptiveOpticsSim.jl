Base.find_package("Proper") === nothing &&
    error("Proper.jl is not available in the active environment. Install it with `using Pkg; Pkg.develop(path=\"../proper.jl\")` before running this benchmark.")

using AdaptiveOpticsSim
using Proper
using Statistics: mean

include(joinpath(@__DIR__, "..", "examples", "support", "proper_hil_coronagraph_common.jl"))
using .ProperHILCoronagraphCommon

function timed_stats!(f!::F; warmup::Int=3, samples::Int=10) where {F<:Function}
    for _ in 1:warmup
        f!()
    end
    GC.gc()
    timings = Vector{Int}(undef, samples)
    @inbounds for i in 1:samples
        t0 = time_ns()
        f!()
        timings[i] = time_ns() - t0
    end
    sorted = sort(timings)
    p95_idx = clamp(round(Int, 0.95 * samples), 1, samples)
    return mean(timings), sorted[p95_idx]
end

function allocated_bytes(f!::F; warmup::Int=3) where {F<:Function}
    for _ in 1:warmup
        f!()
    end
    GC.gc()
    return @allocated f!()
end

function main()
    backend_name = length(ARGS) >= 1 ? ARGS[1] : "cpu"
    warmup = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 3
    samples = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 10

    backend = resolve_hil_backend(backend_name)
    ctx = build_proper_hil_context(; backend=backend)

    ao_mean_ns, ao_p95_ns = timed_stats!(() -> ao_step!(ctx); warmup=warmup, samples=samples)
    science_mean_ns, science_p95_ns = timed_stats!(() -> science_step!(ctx); warmup=warmup, samples=samples)
    total_mean_ns, total_p95_ns = timed_stats!(() -> hil_step!(ctx); warmup=warmup, samples=samples)

    ao_alloc = allocated_bytes(() -> ao_step!(ctx); warmup=warmup)
    science_alloc = allocated_bytes(() -> science_step!(ctx); warmup=warmup)
    total_alloc = allocated_bytes(() -> hil_step!(ctx); warmup=warmup)

    psf, sampling = science_step!(ctx)

    println("proper_hil_coronagraph_profile")
    println("  backend: ", backend_name)
    println("  wavelength_um: ", ctx.wavelength_um)
    println("  ao_step_mean_ns: ", ao_mean_ns)
    println("  ao_step_p95_ns: ", ao_p95_ns)
    println("  science_step_mean_ns: ", science_mean_ns)
    println("  science_step_p95_ns: ", science_p95_ns)
    println("  total_mean_ns: ", total_mean_ns)
    println("  total_p95_ns: ", total_p95_ns)
    println("  frame_rate_hz: ", 1e9 / total_mean_ns)
    println("  ao_step_alloc_bytes: ", ao_alloc)
    println("  science_step_alloc_bytes: ", science_alloc)
    println("  total_alloc_bytes: ", total_alloc)
    println("  slopes_length: ", length(slopes(ctx.scenario)))
    println("  psf_shape: ", size(psf))
    println("  science_sampling: ", sampling)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

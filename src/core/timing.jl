struct RuntimeTimingStats{T<:AbstractFloat}
    samples::Int
    min_ns::Int
    max_ns::Int
    mean_ns::T
    std_ns::T
    p50_ns::T
    p95_ns::T
    p99_ns::T
end

"""
    runtime_timing(f!; warmup=10, samples=1000, gc_before=true)

Measure warmed end-to-end latency with `time_ns()`. This small helper is
intended for low-overhead runtime audits; use the maintained benchmark harness
when coordinated-omission handling, environment capture, or regression
evidence is required.
"""
function runtime_timing(f!::F; warmup::Int=10, samples::Int=1000,
    gc_before::Bool=true) where {F<:Function}
    warmup >= 0 || throw(InvalidConfiguration("warmup must be >= 0"))
    samples > 0 || throw(InvalidConfiguration("samples must be positive"))
    for _ in 1:warmup
        f!()
    end
    gc_before && GC.gc()
    timings = Vector{Int}(undef, samples)
    @inbounds for i in eachindex(timings)
        t0 = time_ns()
        f!()
        timings[i] = time_ns() - t0
    end
    sorted = sort(copy(timings))
    p50_idx = clamp(round(Int, 0.50 * samples), 1, samples)
    p95_idx = clamp(round(Int, 0.95 * samples), 1, samples)
    p99_idx = clamp(round(Int, 0.99 * samples), 1, samples)
    return RuntimeTimingStats{Float64}(
        samples,
        first(sorted),
        last(sorted),
        mean(timings),
        samples > 1 ? std(timings; corrected=false) : 0.0,
        sorted[p50_idx],
        sorted[p95_idx],
        sorted[p99_idx],
    )
end

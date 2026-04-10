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

struct RuntimePhaseTimingStats{T<:AbstractFloat}
    samples::Int
    sense_mean_ns::T
    reconstruct_mean_ns::T
    delay_mean_ns::T
    apply_mean_ns::T
    snapshot_mean_ns::T
    total_mean_ns::T
    total_p95_ns::T
end

"""
    runtime_timing(f!; warmup=10, samples=1000, gc_before=true)

Measure end-to-end runtime latency using warmed `time_ns()` samples.

This is intended for maintained low-overhead HIL timing audits rather than
microbenchmarking individual kernels.
"""
function runtime_timing(f!::F; warmup::Int=10, samples::Int=1000, gc_before::Bool=true) where {F<:Function}
    warmup >= 0 || throw(InvalidConfiguration("warmup must be >= 0"))
    samples > 0 || throw(InvalidConfiguration("samples must be positive"))
    for _ in 1:warmup
        f!()
    end
    gc_before && GC.gc()
    timings = Vector{Int}(undef, samples)
    @inbounds for i in 1:samples
        t0 = time_ns()
        f!()
        timings[i] = time_ns() - t0
    end
    sorted = sort(copy(timings))
    T = Float64
    mean_ns = mean(timings)
    std_ns = samples > 1 ? std(timings; corrected=false) : 0.0
    p50_idx = clamp(round(Int, 0.50 * samples), 1, samples)
    p95_idx = clamp(round(Int, 0.95 * samples), 1, samples)
    p99_idx = clamp(round(Int, 0.99 * samples), 1, samples)
    return RuntimeTimingStats{T}(
        samples,
        sorted[1],
        sorted[end],
        mean_ns,
        std_ns,
        sorted[p50_idx],
        sorted[p95_idx],
        sorted[p99_idx],
    )
end

runtime_timing(runtime::ClosedLoopRuntime; kwargs...) = runtime_timing(() -> step!(runtime); kwargs...)
runtime_timing(interface::SimulationInterface; kwargs...) = runtime_timing(() -> step!(interface); kwargs...)
runtime_timing(interface::CompositeSimulationInterface; kwargs...) = runtime_timing(() -> step!(interface); kwargs...)
runtime_timing(sim::AbstractControlSimulation; kwargs...) = runtime_timing(() -> step!(sim); kwargs...)

@inline function _sync_sense_outputs!(runtime::ClosedLoopRuntime)
    synchronize_backend!(execution_style(runtime.slopes))
    if !isnothing(runtime.wfs_detector)
        synchronize_backend!(execution_style(output_frame(runtime.wfs_detector)))
    end
    if requires_runtime_science_pixels(runtime)
        synchronize_backend!(execution_style(output_frame(runtime.science_detector)))
    end
    return nothing
end

"""
    runtime_phase_timing(runtime; warmup=10, samples=1000, gc_before=true)

Measure the mean latency of the major closed-loop runtime phases.

The returned statistics separate sensing, reconstruction, command application,
and interface snapshot costs so backend regressions can be localized more
easily than with one aggregate timing number.
"""
function runtime_phase_timing(runtime::ClosedLoopRuntime; warmup::Int=10, samples::Int=1000, gc_before::Bool=true)
    warmup >= 0 || throw(InvalidConfiguration("warmup must be >= 0"))
    samples > 0 || throw(InvalidConfiguration("samples must be positive"))
    for _ in 1:warmup
        step!(runtime)
    end
    gc_before && GC.gc()
    sense_times = Vector{Int}(undef, samples)
    reconstruct_times = Vector{Int}(undef, samples)
    delay_times = Vector{Int}(undef, samples)
    apply_times = Vector{Int}(undef, samples)
    total_times = Vector{Int}(undef, samples)
    @inbounds for i in 1:samples
        t0 = time_ns()
        sense!(runtime)
        _sync_sense_outputs!(runtime)
        t1 = time_ns()
        reconstruct!(runtime.reconstruct_buffer, runtime.reconstructor, runtime.slopes)
        synchronize_backend!(execution_style(runtime.reconstruct_buffer))
        t2 = time_ns()
        stage_runtime_command!(runtime)
        synchronize_backend!(execution_style(runtime.command))
        t3 = time_ns()
        delayed = shift_delay!(runtime.dm_delay, runtime.command)
        stage_command!(runtime.optic, delayed, runtime.control_sign)
        synchronize_backend!(execution_style(command_storage(runtime.optic)))
        t4 = time_ns()
        sense_times[i] = t1 - t0
        reconstruct_times[i] = t2 - t1
        delay_times[i] = t3 - t2
        apply_times[i] = t4 - t3
        total_times[i] = t4 - t0
    end
    sorted_total = sort(copy(total_times))
    p95_idx = clamp(round(Int, 0.95 * samples), 1, samples)
    return RuntimePhaseTimingStats{Float64}(
        samples,
        mean(sense_times),
        mean(reconstruct_times),
        mean(delay_times),
        mean(apply_times),
        0.0,
        mean(total_times),
        sorted_total[p95_idx],
    )
end

function runtime_phase_timing(interface::SimulationInterface; warmup::Int=10, samples::Int=1000, gc_before::Bool=true)
    warmup >= 0 || throw(InvalidConfiguration("warmup must be >= 0"))
    samples > 0 || throw(InvalidConfiguration("samples must be positive"))
    for _ in 1:warmup
        step!(interface)
    end
    gc_before && GC.gc()
    sense_times = Vector{Int}(undef, samples)
    reconstruct_times = Vector{Int}(undef, samples)
    delay_times = Vector{Int}(undef, samples)
    apply_times = Vector{Int}(undef, samples)
    snapshot_times = Vector{Int}(undef, samples)
    total_times = Vector{Int}(undef, samples)
    runtime = interface.runtime
    @inbounds for i in 1:samples
        t0 = time_ns()
        sense!(runtime)
        _sync_sense_outputs!(runtime)
        t1 = time_ns()
        reconstruct!(runtime.reconstruct_buffer, runtime.reconstructor, runtime.slopes)
        synchronize_backend!(execution_style(runtime.reconstruct_buffer))
        t2 = time_ns()
        stage_runtime_command!(runtime)
        synchronize_backend!(execution_style(runtime.command))
        t3 = time_ns()
        delayed = shift_delay!(runtime.dm_delay, runtime.command)
        stage_command!(runtime.optic, delayed, runtime.control_sign)
        synchronize_backend!(execution_style(command_storage(runtime.optic)))
        t4 = time_ns()
        snapshot_outputs!(interface)
        synchronize_backend!(execution_style(interface.command))
        synchronize_backend!(execution_style(interface.slopes))
        if !isnothing(interface.wfs_frame)
            synchronize_backend!(execution_style(interface.wfs_frame))
        end
        if !isnothing(interface.science_frame)
            synchronize_backend!(execution_style(interface.science_frame))
        end
        t5 = time_ns()
        sense_times[i] = t1 - t0
        reconstruct_times[i] = t2 - t1
        delay_times[i] = t3 - t2
        apply_times[i] = t4 - t3
        snapshot_times[i] = t5 - t4
        total_times[i] = t5 - t0
    end
    sorted_total = sort(copy(total_times))
    p95_idx = clamp(round(Int, 0.95 * samples), 1, samples)
    return RuntimePhaseTimingStats{Float64}(
        samples,
        mean(sense_times),
        mean(reconstruct_times),
        mean(delay_times),
        mean(apply_times),
        mean(snapshot_times),
        mean(total_times),
        sorted_total[p95_idx],
    )
end

function runtime_phase_timing(interface::CompositeSimulationInterface; warmup::Int=10, samples::Int=1000, gc_before::Bool=true)
    warmup >= 0 || throw(InvalidConfiguration("warmup must be >= 0"))
    samples > 0 || throw(InvalidConfiguration("samples must be positive"))
    for _ in 1:warmup
        step!(interface)
    end
    gc_before && GC.gc()
    sense_times = Vector{Int}(undef, samples)
    reconstruct_times = Vector{Int}(undef, samples)
    delay_times = Vector{Int}(undef, samples)
    apply_times = Vector{Int}(undef, samples)
    snapshot_times = Vector{Int}(undef, samples)
    total_times = Vector{Int}(undef, samples)
    @inbounds for i in 1:samples
        t0 = time_ns()
        for child in interface.interfaces
            sense!(child.runtime)
            _sync_sense_outputs!(child.runtime)
        end
        t1 = time_ns()
        for child in interface.interfaces
            reconstruct!(child.runtime.reconstruct_buffer, child.runtime.reconstructor, child.runtime.slopes)
            synchronize_backend!(execution_style(child.runtime.reconstruct_buffer))
        end
        t2 = time_ns()
        for child in interface.interfaces
            stage_runtime_command!(child.runtime)
            synchronize_backend!(execution_style(child.runtime.command))
        end
        t3 = time_ns()
        for child in interface.interfaces
            delayed = shift_delay!(child.runtime.dm_delay, child.runtime.command)
            stage_command!(child.runtime.optic, delayed, child.runtime.control_sign)
            synchronize_backend!(execution_style(command_storage(child.runtime.optic)))
        end
        t4 = time_ns()
        snapshot_outputs!(interface)
        synchronize_backend!(execution_style(interface.command))
        synchronize_backend!(execution_style(interface.slopes))
        t5 = time_ns()
        sense_times[i] = t1 - t0
        reconstruct_times[i] = t2 - t1
        delay_times[i] = t3 - t2
        apply_times[i] = t4 - t3
        snapshot_times[i] = t5 - t4
        total_times[i] = t5 - t0
    end
    sorted_total = sort(copy(total_times))
    p95_idx = clamp(round(Int, 0.95 * samples), 1, samples)
    return RuntimePhaseTimingStats{Float64}(
        samples,
        mean(sense_times),
        mean(reconstruct_times),
        mean(delay_times),
        mean(apply_times),
        mean(snapshot_times),
        mean(total_times),
        sorted_total[p95_idx],
    )
end

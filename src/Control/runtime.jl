@kernel function apply_command_kernel!(coefs, cmd, sign, n::Int)
    i = @index(Global, Linear)
    if i <= n
        @inbounds coefs[i] = sign * cmd[i]
    end
end

mutable struct ClosedLoopRuntime{SIM<:AOSimulation,TEL,A,S,DM,W,R,V,WD,SD,RNG,T<:AbstractFloat}
    simulation::SIM
    tel::TEL
    atm::A
    src::S
    dm::DM
    wfs::W
    reconstructor::R
    command::V
    wfs_detector::WD
    science_detector::SD
    rng::RNG
    control_sign::T
    science_zero_padding::Int
end

mutable struct RTCBoundary{RT,C,S,W,SF}
    runtime::RT
    command::C
    slopes::S
    wfs_frame::W
    science_frame::SF
end

mutable struct MultiRTCBoundary{BT,C,S,WF,SF}
    boundaries::BT
    command::C
    slopes::S
    wfs_frames::WF
    science_frames::SF
end

@inline function apply_command!(::ScalarCPUStyle, coefs::AbstractVector{T}, cmd::AbstractVector{T}, sign::T) where {T<:AbstractFloat}
    @inbounds for i in eachindex(coefs, cmd)
        coefs[i] = sign * cmd[i]
    end
    return coefs
end

@inline function apply_command!(style::AcceleratorStyle, coefs::AbstractVector{T}, cmd::AbstractVector{T}, sign::T) where {T<:AbstractFloat}
    launch_kernel!(style, apply_command_kernel!, coefs, cmd, sign, length(coefs); ndrange=length(coefs))
    return coefs
end

@inline function apply_command!(coefs::AbstractVector{T}, cmd::AbstractVector{T}, sign::T) where {T<:AbstractFloat}
    apply_command!(execution_style(coefs), coefs, cmd, sign)
end

function ClosedLoopRuntime(simulation::AOSimulation, reconstructor;
    wfs_detector=nothing, science_detector=nothing, rng=MersenneTwister(0),
    control_sign::Real=-1.0, science_zero_padding::Int=2)
    T = eltype(simulation.dm.state.coefs)
    command = similar(simulation.dm.state.coefs)
    fill!(command, zero(T))
    return ClosedLoopRuntime{typeof(simulation), typeof(simulation.tel), typeof(simulation.atm),
        typeof(simulation.src), typeof(simulation.dm), typeof(simulation.wfs), typeof(reconstructor),
        typeof(command), typeof(wfs_detector), typeof(science_detector), typeof(rng), T}(
        simulation,
        simulation.tel,
        simulation.atm,
        simulation.src,
        simulation.dm,
        simulation.wfs,
        reconstructor,
        command,
        wfs_detector,
        science_detector,
        rng,
        T(control_sign),
        science_zero_padding,
    )
end

function RTCBoundary(runtime::ClosedLoopRuntime)
    command = similar(runtime.command)
    copyto!(command, runtime.command)
    slopes = similar(runtime.wfs.state.slopes)
    copyto!(slopes, runtime.wfs.state.slopes)
    wfs_frame = isnothing(runtime.wfs_detector) ? nothing : similar(output_frame(runtime.wfs_detector))
    science_frame = isnothing(runtime.science_detector) ? nothing : similar(output_frame(runtime.science_detector))
    boundary = RTCBoundary{typeof(runtime), typeof(command), typeof(slopes), typeof(wfs_frame), typeof(science_frame)}(
        runtime,
        command,
        slopes,
        wfs_frame,
        science_frame,
    )
    return snapshot_outputs!(boundary)
end

function MultiRTCBoundary(boundaries::RTCBoundary...)
    length(boundaries) > 0 || throw(InvalidConfiguration("MultiRTCBoundary requires at least one boundary"))
    first_boundary = boundaries[1]
    total_command = mapreduce(b -> length(rtc_command(b)), +, boundaries)
    total_slopes = mapreduce(b -> length(rtc_slopes(b)), +, boundaries)
    command = similar(rtc_command(first_boundary), total_command)
    slopes = similar(rtc_slopes(first_boundary), total_slopes)
    wfs_frames = ntuple(i -> begin
        frame = rtc_wfs_frame(boundaries[i])
        isnothing(frame) ? nothing : similar(frame)
    end, length(boundaries))
    science_frames = ntuple(i -> begin
        frame = rtc_science_frame(boundaries[i])
        isnothing(frame) ? nothing : similar(frame)
    end, length(boundaries))
    multi = MultiRTCBoundary{typeof(boundaries), typeof(command), typeof(slopes), typeof(wfs_frames), typeof(science_frames)}(
        boundaries,
        command,
        slopes,
        wfs_frames,
        science_frames,
    )
    return snapshot_outputs!(multi)
end

MultiRTCBoundary(runtimes::ClosedLoopRuntime...) = MultiRTCBoundary(map(RTCBoundary, runtimes)...)

function with_reconstructor(runtime::ClosedLoopRuntime, reconstructor)
    refreshed = ClosedLoopRuntime(
        runtime.simulation,
        reconstructor;
        wfs_detector=runtime.wfs_detector,
        science_detector=runtime.science_detector,
        rng=runtime.rng,
        control_sign=runtime.control_sign,
        science_zero_padding=runtime.science_zero_padding,
    )
    copyto!(refreshed.command, runtime.command)
    copyto!(refreshed.dm.state.coefs, runtime.dm.state.coefs)
    return refreshed
end

with_reconstructor(boundary::RTCBoundary, reconstructor) = RTCBoundary(with_reconstructor(boundary.runtime, reconstructor))

function with_reconstructors(boundary::MultiRTCBoundary, reconstructors::Vararg{Any,N}) where {N}
    length(boundary.boundaries) == N ||
        throw(DimensionMismatchError("reconstructor count must match boundary count"))
    refreshed = ntuple(i -> with_reconstructor(boundary.boundaries[i], reconstructors[i]), N)
    return MultiRTCBoundary(refreshed...)
end

@inline function set_command!(runtime::ClosedLoopRuntime, command::AbstractVector)
    copyto!(runtime.command, command)
    copyto!(runtime.dm.state.coefs, command)
    return runtime.dm.state.coefs
end

@inline function set_command!(boundary::RTCBoundary, command::AbstractVector)
    copyto!(boundary.command, command)
    set_command!(boundary.runtime, command)
    return boundary.command
end

@inline function set_command!(boundary::MultiRTCBoundary, command::AbstractVector)
    length(boundary.command) == length(command) ||
        throw(DimensionMismatchError("command length must match aggregated RTC command length"))
    command_offset = 1
    @inbounds for child in boundary.boundaries
        n_command = length(child.command)
        @views set_command!(child, command[command_offset:command_offset + n_command - 1])
        command_offset += n_command
    end
    copyto!(boundary.command, command)
    return boundary.command
end

@inline function snapshot_outputs!(boundary::RTCBoundary)
    copyto!(boundary.command, boundary.runtime.command)
    copyto!(boundary.slopes, boundary.runtime.wfs.state.slopes)
    if !isnothing(boundary.wfs_frame)
        copyto!(boundary.wfs_frame, output_frame(boundary.runtime.wfs_detector))
    end
    if !isnothing(boundary.science_frame)
        copyto!(boundary.science_frame, output_frame(boundary.runtime.science_detector))
    end
    return boundary
end

@inline function snapshot_outputs!(multi::MultiRTCBoundary)
    command_offset = 1
    slope_offset = 1
    @inbounds for i in eachindex(multi.boundaries)
        boundary = snapshot_outputs!(multi.boundaries[i])
        n_command = length(boundary.command)
        n_slopes = length(boundary.slopes)
        @views copyto!(multi.command[command_offset:command_offset + n_command - 1], boundary.command)
        @views copyto!(multi.slopes[slope_offset:slope_offset + n_slopes - 1], boundary.slopes)
        if !isnothing(multi.wfs_frames[i])
            copyto!(multi.wfs_frames[i], boundary.wfs_frame)
        end
        if !isnothing(multi.science_frames[i])
            copyto!(multi.science_frames[i], boundary.science_frame)
        end
        command_offset += n_command
        slope_offset += n_slopes
    end
    return multi
end

@inline rtc_slopes(boundary::RTCBoundary) = boundary.slopes
@inline rtc_command(boundary::RTCBoundary) = boundary.command
@inline rtc_wfs_frame(boundary::RTCBoundary) = boundary.wfs_frame
@inline rtc_science_frame(boundary::RTCBoundary) = boundary.science_frame
@inline rtc_wfs_metadata(boundary::RTCBoundary) = isnothing(boundary.runtime.wfs_detector) ? nothing : detector_export_metadata(boundary.runtime.wfs_detector)
@inline rtc_science_metadata(boundary::RTCBoundary) = isnothing(boundary.runtime.science_detector) ? nothing : detector_export_metadata(boundary.runtime.science_detector)
@inline rtc_slopes(boundary::MultiRTCBoundary) = boundary.slopes
@inline rtc_command(boundary::MultiRTCBoundary) = boundary.command
@inline rtc_wfs_frame(boundary::MultiRTCBoundary) = boundary.wfs_frames
@inline rtc_science_frame(boundary::MultiRTCBoundary) = boundary.science_frames
@inline rtc_wfs_metadata(boundary::MultiRTCBoundary) = map(rtc_wfs_metadata, boundary.boundaries)
@inline rtc_science_metadata(boundary::MultiRTCBoundary) = map(rtc_science_metadata, boundary.boundaries)

@inline function sense_core!(atm::AbstractAtmosphere, tel::Telescope, dm::DeformableMirror,
    wfs::AbstractWFS, src::AbstractSource, rng::AbstractRNG)
    advance!(atm, tel, rng)
    propagate!(atm, tel)
    apply!(dm, tel, DMAdditive())
    measure!(wfs, tel, src)
    return nothing
end

@inline function sense_core!(atm::AbstractAtmosphere, tel::Telescope, dm::DeformableMirror,
    wfs::AbstractWFS, src::AbstractSource, det::AbstractDetector, rng::AbstractRNG)
    advance!(atm, tel, rng)
    propagate!(atm, tel)
    apply!(dm, tel, DMAdditive())
    measure!(wfs, tel, src, det; rng=rng)
    return nothing
end

@inline function capture_science_core!(tel::Telescope, src::AbstractSource,
    det::AbstractDetector, rng::AbstractRNG, zero_padding::Int)
    psf = compute_psf!(tel, src; zero_padding=zero_padding)
    capture!(det, psf; rng=rng)
    return nothing
end

@inline function step_core!(atm::AbstractAtmosphere, tel::Telescope, dm::DeformableMirror,
    wfs::AbstractWFS, src::AbstractSource, reconstructor, command::AbstractVector{T},
    rng::AbstractRNG, control_sign::T) where {T<:AbstractFloat}
    sense_core!(atm, tel, dm, wfs, src, rng)
    reconstruct!(command, reconstructor, wfs.state.slopes)
    apply_command!(dm.state.coefs, command, control_sign)
    return nothing
end

@inline function step_core!(atm::AbstractAtmosphere, tel::Telescope, dm::DeformableMirror,
    wfs::AbstractWFS, src::AbstractSource, wfs_detector::AbstractDetector, reconstructor,
    command::AbstractVector{T}, rng::AbstractRNG, control_sign::T) where {T<:AbstractFloat}
    sense_core!(atm, tel, dm, wfs, src, wfs_detector, rng)
    reconstruct!(command, reconstructor, wfs.state.slopes)
    apply_command!(dm.state.coefs, command, control_sign)
    return nothing
end

@inline function step_core!(atm::AbstractAtmosphere, tel::Telescope, dm::DeformableMirror,
    wfs::AbstractWFS, src::AbstractSource, reconstructor, command::AbstractVector{T},
    science_detector::AbstractDetector, rng::AbstractRNG, control_sign::T, science_zero_padding::Int) where {T<:AbstractFloat}
    sense_core!(atm, tel, dm, wfs, src, rng)
    capture_science_core!(tel, src, science_detector, rng, science_zero_padding)
    reconstruct!(command, reconstructor, wfs.state.slopes)
    apply_command!(dm.state.coefs, command, control_sign)
    return nothing
end

@inline function step_core!(atm::AbstractAtmosphere, tel::Telescope, dm::DeformableMirror,
    wfs::AbstractWFS, src::AbstractSource, wfs_detector::AbstractDetector, reconstructor,
    command::AbstractVector{T}, science_detector::AbstractDetector, rng::AbstractRNG,
    control_sign::T, science_zero_padding::Int) where {T<:AbstractFloat}
    sense_core!(atm, tel, dm, wfs, src, wfs_detector, rng)
    capture_science_core!(tel, src, science_detector, rng, science_zero_padding)
    reconstruct!(command, reconstructor, wfs.state.slopes)
    apply_command!(dm.state.coefs, command, control_sign)
    return nothing
end

function sense!(runtime::ClosedLoopRuntime{SIM,TEL,A,S,DM,W,R,V,Nothing,Nothing,RNG,T}) where {SIM<:AOSimulation,TEL,A,S,DM,W,R,V,RNG,T<:AbstractFloat}
    sense_core!(runtime.atm, runtime.tel, runtime.dm, runtime.wfs, runtime.src, runtime.rng)
    return runtime
end

function sense!(boundary::RTCBoundary)
    sense!(boundary.runtime)
    return snapshot_outputs!(boundary)
end

function sense!(boundary::MultiRTCBoundary)
    @inbounds for child in boundary.boundaries
        sense!(child.runtime)
    end
    return snapshot_outputs!(boundary)
end

function sense!(runtime::ClosedLoopRuntime{SIM,TEL,A,S,DM,W,R,V,WD,Nothing,RNG,T}) where {SIM<:AOSimulation,TEL,A,S,DM,W,R,V,WD<:AbstractDetector,RNG,T<:AbstractFloat}
    sense_core!(runtime.atm, runtime.tel, runtime.dm, runtime.wfs, runtime.src, runtime.wfs_detector, runtime.rng)
    return runtime
end

function sense!(runtime::ClosedLoopRuntime{SIM,TEL,A,S,DM,W,R,V,Nothing,SD,RNG,T}) where {SIM<:AOSimulation,TEL,A,S,DM,W,R,V,SD<:AbstractDetector,RNG,T<:AbstractFloat}
    sense_core!(runtime.atm, runtime.tel, runtime.dm, runtime.wfs, runtime.src, runtime.rng)
    capture_science_core!(runtime.tel, runtime.src, runtime.science_detector, runtime.rng, runtime.science_zero_padding)
    return runtime
end

function sense!(runtime::ClosedLoopRuntime{SIM,TEL,A,S,DM,W,R,V,WD,SD,RNG,T}) where {SIM<:AOSimulation,TEL,A,S,DM,W,R,V,WD<:AbstractDetector,SD<:AbstractDetector,RNG,T<:AbstractFloat}
    sense_core!(runtime.atm, runtime.tel, runtime.dm, runtime.wfs, runtime.src, runtime.wfs_detector, runtime.rng)
    capture_science_core!(runtime.tel, runtime.src, runtime.science_detector, runtime.rng, runtime.science_zero_padding)
    return runtime
end

function step!(runtime::ClosedLoopRuntime{SIM,TEL,A,S,DM,W,R,V,Nothing,Nothing,RNG,T}) where {SIM<:AOSimulation,TEL,A,S,DM,W,R,V,RNG,T<:AbstractFloat}
    step_core!(runtime.atm, runtime.tel, runtime.dm, runtime.wfs, runtime.src, runtime.reconstructor,
        runtime.command, runtime.rng, runtime.control_sign)
    return runtime
end

function step!(runtime::ClosedLoopRuntime{SIM,TEL,A,S,DM,W,R,V,WD,Nothing,RNG,T}) where {SIM<:AOSimulation,TEL,A,S,DM,W,R,V,WD<:AbstractDetector,RNG,T<:AbstractFloat}
    step_core!(runtime.atm, runtime.tel, runtime.dm, runtime.wfs, runtime.src, runtime.wfs_detector,
        runtime.reconstructor, runtime.command, runtime.rng, runtime.control_sign)
    return runtime
end

function step!(runtime::ClosedLoopRuntime{SIM,TEL,A,S,DM,W,R,V,Nothing,SD,RNG,T}) where {SIM<:AOSimulation,TEL,A,S,DM,W,R,V,SD<:AbstractDetector,RNG,T<:AbstractFloat}
    step_core!(runtime.atm, runtime.tel, runtime.dm, runtime.wfs, runtime.src, runtime.reconstructor,
        runtime.command, runtime.science_detector, runtime.rng, runtime.control_sign,
        runtime.science_zero_padding)
    return runtime
end

function step!(runtime::ClosedLoopRuntime{SIM,TEL,A,S,DM,W,R,V,WD,SD,RNG,T}) where {SIM<:AOSimulation,TEL,A,S,DM,W,R,V,WD<:AbstractDetector,SD<:AbstractDetector,RNG,T<:AbstractFloat}
    step_core!(runtime.atm, runtime.tel, runtime.dm, runtime.wfs, runtime.src, runtime.wfs_detector,
        runtime.reconstructor, runtime.command, runtime.science_detector, runtime.rng,
        runtime.control_sign, runtime.science_zero_padding)
    return runtime
end

function step!(boundary::RTCBoundary)
    step!(boundary.runtime)
    return snapshot_outputs!(boundary)
end

function step!(boundary::MultiRTCBoundary)
    @inbounds for child in boundary.boundaries
        step!(child.runtime)
    end
    return snapshot_outputs!(boundary)
end

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
    apply_mean_ns::T
    snapshot_mean_ns::T
    total_mean_ns::T
    total_p95_ns::T
end

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
runtime_timing(boundary::RTCBoundary; kwargs...) = runtime_timing(() -> step!(boundary); kwargs...)
runtime_timing(boundary::MultiRTCBoundary; kwargs...) = runtime_timing(() -> step!(boundary); kwargs...)

@inline function _sync_sense_outputs!(runtime::ClosedLoopRuntime)
    synchronize_backend!(execution_style(runtime.wfs.state.slopes))
    if !isnothing(runtime.wfs_detector)
        synchronize_backend!(execution_style(output_frame(runtime.wfs_detector)))
    end
    if !isnothing(runtime.science_detector)
        synchronize_backend!(execution_style(output_frame(runtime.science_detector)))
    end
    return nothing
end

function runtime_phase_timing(runtime::ClosedLoopRuntime; warmup::Int=10, samples::Int=1000, gc_before::Bool=true)
    warmup >= 0 || throw(InvalidConfiguration("warmup must be >= 0"))
    samples > 0 || throw(InvalidConfiguration("samples must be positive"))
    for _ in 1:warmup
        step!(runtime)
    end
    gc_before && GC.gc()
    sense_times = Vector{Int}(undef, samples)
    reconstruct_times = Vector{Int}(undef, samples)
    apply_times = Vector{Int}(undef, samples)
    total_times = Vector{Int}(undef, samples)
    @inbounds for i in 1:samples
        t0 = time_ns()
        sense!(runtime)
        _sync_sense_outputs!(runtime)
        t1 = time_ns()
        reconstruct!(runtime.command, runtime.reconstructor, runtime.wfs.state.slopes)
        synchronize_backend!(execution_style(runtime.command))
        t2 = time_ns()
        apply_command!(runtime.dm.state.coefs, runtime.command, runtime.control_sign)
        synchronize_backend!(execution_style(runtime.dm.state.coefs))
        t3 = time_ns()
        sense_times[i] = t1 - t0
        reconstruct_times[i] = t2 - t1
        apply_times[i] = t3 - t2
        total_times[i] = t3 - t0
    end
    sorted_total = sort(copy(total_times))
    p95_idx = clamp(round(Int, 0.95 * samples), 1, samples)
    return RuntimePhaseTimingStats{Float64}(
        samples,
        mean(sense_times),
        mean(reconstruct_times),
        mean(apply_times),
        0.0,
        mean(total_times),
        sorted_total[p95_idx],
    )
end

function runtime_phase_timing(boundary::RTCBoundary; warmup::Int=10, samples::Int=1000, gc_before::Bool=true)
    warmup >= 0 || throw(InvalidConfiguration("warmup must be >= 0"))
    samples > 0 || throw(InvalidConfiguration("samples must be positive"))
    for _ in 1:warmup
        step!(boundary)
    end
    gc_before && GC.gc()
    sense_times = Vector{Int}(undef, samples)
    reconstruct_times = Vector{Int}(undef, samples)
    apply_times = Vector{Int}(undef, samples)
    snapshot_times = Vector{Int}(undef, samples)
    total_times = Vector{Int}(undef, samples)
    runtime = boundary.runtime
    @inbounds for i in 1:samples
        t0 = time_ns()
        sense!(runtime)
        _sync_sense_outputs!(runtime)
        t1 = time_ns()
        reconstruct!(runtime.command, runtime.reconstructor, runtime.wfs.state.slopes)
        synchronize_backend!(execution_style(runtime.command))
        t2 = time_ns()
        apply_command!(runtime.dm.state.coefs, runtime.command, runtime.control_sign)
        synchronize_backend!(execution_style(runtime.dm.state.coefs))
        t3 = time_ns()
        snapshot_outputs!(boundary)
        synchronize_backend!(execution_style(boundary.command))
        synchronize_backend!(execution_style(boundary.slopes))
        if !isnothing(boundary.wfs_frame)
            synchronize_backend!(execution_style(boundary.wfs_frame))
        end
        if !isnothing(boundary.science_frame)
            synchronize_backend!(execution_style(boundary.science_frame))
        end
        t4 = time_ns()
        sense_times[i] = t1 - t0
        reconstruct_times[i] = t2 - t1
        apply_times[i] = t3 - t2
        snapshot_times[i] = t4 - t3
        total_times[i] = t4 - t0
    end
    sorted_total = sort(copy(total_times))
    p95_idx = clamp(round(Int, 0.95 * samples), 1, samples)
    return RuntimePhaseTimingStats{Float64}(
        samples,
        mean(sense_times),
        mean(reconstruct_times),
        mean(apply_times),
        mean(snapshot_times),
        mean(total_times),
        sorted_total[p95_idx],
    )
end

function runtime_phase_timing(boundary::MultiRTCBoundary; warmup::Int=10, samples::Int=1000, gc_before::Bool=true)
    warmup >= 0 || throw(InvalidConfiguration("warmup must be >= 0"))
    samples > 0 || throw(InvalidConfiguration("samples must be positive"))
    for _ in 1:warmup
        step!(boundary)
    end
    gc_before && GC.gc()
    sense_times = Vector{Int}(undef, samples)
    reconstruct_times = Vector{Int}(undef, samples)
    apply_times = Vector{Int}(undef, samples)
    snapshot_times = Vector{Int}(undef, samples)
    total_times = Vector{Int}(undef, samples)
    @inbounds for i in 1:samples
        t0 = time_ns()
        for child in boundary.boundaries
            sense!(child.runtime)
            _sync_sense_outputs!(child.runtime)
        end
        t1 = time_ns()
        for child in boundary.boundaries
            reconstruct!(child.runtime.command, child.runtime.reconstructor, child.runtime.wfs.state.slopes)
            synchronize_backend!(execution_style(child.runtime.command))
        end
        t2 = time_ns()
        for child in boundary.boundaries
            apply_command!(child.runtime.dm.state.coefs, child.runtime.command, child.runtime.control_sign)
            synchronize_backend!(execution_style(child.runtime.dm.state.coefs))
        end
        t3 = time_ns()
        snapshot_outputs!(boundary)
        synchronize_backend!(execution_style(boundary.command))
        synchronize_backend!(execution_style(boundary.slopes))
        t4 = time_ns()
        sense_times[i] = t1 - t0
        reconstruct_times[i] = t2 - t1
        apply_times[i] = t3 - t2
        snapshot_times[i] = t4 - t3
        total_times[i] = t4 - t0
    end
    sorted_total = sort(copy(total_times))
    p95_idx = clamp(round(Int, 0.95 * samples), 1, samples)
    return RuntimePhaseTimingStats{Float64}(
        samples,
        mean(sense_times),
        mean(reconstruct_times),
        mean(apply_times),
        mean(snapshot_times),
        mean(total_times),
        sorted_total[p95_idx],
    )
end

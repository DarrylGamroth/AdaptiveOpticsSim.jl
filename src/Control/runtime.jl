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

@inline rtc_slopes(boundary::RTCBoundary) = boundary.slopes
@inline rtc_command(boundary::RTCBoundary) = boundary.command
@inline rtc_wfs_frame(boundary::RTCBoundary) = boundary.wfs_frame
@inline rtc_science_frame(boundary::RTCBoundary) = boundary.science_frame

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

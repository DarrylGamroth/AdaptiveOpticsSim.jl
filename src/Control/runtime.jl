#
# Closed-loop runtime execution and simulation I/O
#
# This file defines the low-latency control loop used by HIL and deterministic
# runtime benchmarks.
#
# The core runtime step is:
# 1. advance and propagate the atmosphere
# 2. apply the current DM command to the telescope OPD
# 3. measure WFS outputs, optionally through a detector
# 4. reconstruct a new command from the sensed slopes
# 5. apply the signed command back to the DM coefficient vector
#
# `SimulationInterface` and `SimulationReadout` expose that state as the
# controller-facing I/O surface without changing the underlying runtime object.
#
@kernel function apply_command_kernel!(coefs, cmd, sign, n::Int)
    i = @index(Global, Linear)
    if i <= n
        @inbounds coefs[i] = sign * cmd[i]
    end
end

abstract type AbstractControlSimulation end
abstract type AbstractExecutionPolicy end

struct SequentialExecution <: AbstractExecutionPolicy end
struct ThreadedExecution <: AbstractExecutionPolicy end
struct BackendStreamExecution <: AbstractExecutionPolicy end

"""
    VectorDelayLine(ref, delay_frames)

Preallocated fixed-length frame delay for vector-valued control signals.

The delay line stores `delay_frames` historical samples and returns the signal
that exits the tail when a new sample is shifted in.
"""
mutable struct VectorDelayLine{A<:AbstractMatrix,V<:AbstractVector}
    buffer::A
    scratch::V
end

function VectorDelayLine(ref::AbstractVector{T}, delay_frames::Int) where {T<:AbstractFloat}
    delay_frames >= 0 || throw(InvalidConfiguration("delay_frames must be >= 0"))
    buffer = similar(ref, T, length(ref), delay_frames)
    fill!(buffer, zero(T))
    scratch = similar(ref, T, length(ref))
    fill!(scratch, zero(T))
    return VectorDelayLine{typeof(buffer),typeof(scratch)}(buffer, scratch)
end

function shift_delay!(line::VectorDelayLine, sample::AbstractVector)
    n_delay = size(line.buffer, 2)
    if n_delay == 0
        copyto!(line.scratch, sample)
        return line.scratch
    end
    copyto!(line.scratch, @view(line.buffer[:, 1]))
    @inbounds for i in 1:n_delay-1
        copyto!(@view(line.buffer[:, i]), @view(line.buffer[:, i + 1]))
    end
    copyto!(@view(line.buffer[:, n_delay]), sample)
    return line.scratch
end

prepare!(sim::AbstractControlSimulation) = sim
supports_prepared_runtime(::Type) = false
supports_prepared_runtime(sim) = supports_prepared_runtime(typeof(sim))
supports_detector_output(::Type) = false
supports_detector_output(sim) = supports_detector_output(typeof(sim))
supports_stacked_sources(::Type) = false
supports_stacked_sources(sim) = supports_stacked_sources(typeof(sim))
supports_grouped_execution(::Type) = false
supports_grouped_execution(sim) = supports_grouped_execution(typeof(sim))

"""
    supports_prepared_runtime(wfs, src)

Return whether a WFS/source pairing exposes a meaningful `prepare_runtime_wfs!`
precomputation hook for repeated runtime stepping.
"""
supports_prepared_runtime(::AbstractWFS, ::Any) = false

"""
    supports_stacked_sources(wfs, src)

Return whether the WFS/source pairing has maintained stacked-source execution
support, typically for `Asterism` inputs.
"""
supports_stacked_sources(::AbstractWFS, ::Any) = false

init_execution_state(::AbstractExecutionPolicy, ref) = nothing

"""
    ClosedLoopRuntime(simulation, reconstructor; ...)

Create the executable closed-loop control runtime for an AO simulation.

The runtime owns the mutable command vector, optional detector outputs, and the
objects needed to execute repeated `sense!` and `step!` calls without
reconstructing the simulation graph.
"""
mutable struct ClosedLoopRuntime{SIM<:AOSimulation,TEL,A,S,DM,W,R,V,WD,SD,RNG,T<:AbstractFloat} <: AbstractControlSimulation
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

"""
    SimulationInterface(runtime)

Allocate a controller-facing I/O view for a runtime.

The interface owns copied command, slope, and frame buffers so external code
can snapshot or stream runtime state without mutating the runtime internals.
"""
mutable struct SimulationInterface{RT,C,S,W,SF}
    runtime::RT
    command::C
    slopes::S
    wfs_frame::W
    science_frame::SF
end

"""
    CompositeSimulationInterface(interfaces...)

Aggregate several compatible simulation interfaces into one combined command,
slope, and frame surface.

This is used for grouped multi-branch execution while preserving per-runtime
state internally.
"""
mutable struct CompositeSimulationInterface{IT,C,S,WF,SF}
    interfaces::IT
    command::C
    slopes::S
    wfs_frames::WF
    science_frames::SF
end

"""
    SimulationReadout

Zero-copy view of the current output side of a simulation runtime or interface.

This bundles the command, slopes, optional frames, and export metadata that an
external controller or logger would observe at the current simulation state.
"""
struct SimulationReadout{C,S,W,SF,WM,SM}
    command::C
    slopes::S
    wfs_frame::W
    science_frame::SF
    wfs_metadata::WM
    science_metadata::SM
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

@inline wfs_output_frame(wfs::AbstractWFS, det::AbstractDetector) = output_frame(det)
@inline wfs_output_frame(wfs::ShackHartmann{<:Diffractive}, det::AbstractDetector) = wfs.state.spot_cube

simulation_interface(runtime::ClosedLoopRuntime) = SimulationInterface(runtime)
simulation_interface(interface::SimulationInterface) = interface
simulation_interface(interface::CompositeSimulationInterface) = interface

@inline supports_prepared_runtime(::ShackHartmann{<:Diffractive}, ::AbstractSource) = true
@inline supports_prepared_runtime(::ShackHartmann{<:Diffractive}, ::Asterism) = true
@inline supports_stacked_sources(::ShackHartmann, ::Asterism) = true
@inline supports_stacked_sources(::PyramidWFS, ::Asterism) = true
@inline supports_stacked_sources(::BioEdgeWFS, ::Asterism) = true

"""
    prepare_runtime_wfs!(wfs, tel, src)

Run any WFS-specific runtime precomputation needed before repeated `step!`
calls. The default implementation is a no-op and returns `wfs`.
"""
@inline function prepare_runtime_wfs!(wfs::AbstractWFS, tel::Telescope, src)
    return wfs
end

@inline function prepare_runtime_wfs!(wfs::ShackHartmann{<:Diffractive}, tel::Telescope, src::AbstractSource)
    prepare_sampling!(wfs, tel, src)
    ensure_sh_calibration!(wfs, tel, src)
    return wfs
end

@inline function prepare_runtime_wfs!(wfs::ShackHartmann{<:Diffractive}, tel::Telescope, ast::Asterism)
    isempty(ast.sources) && throw(InvalidConfiguration("asterism must contain at least one source"))
    prepare_sampling!(wfs, tel, ast.sources[1])
    ensure_sh_calibration!(wfs, tel, ast.sources[1])
    return wfs
end

supports_prepared_runtime(runtime::ClosedLoopRuntime) = supports_prepared_runtime(runtime.wfs, runtime.src)
supports_detector_output(runtime::ClosedLoopRuntime) = !isnothing(runtime.wfs_detector) || !isnothing(runtime.science_detector)
supports_stacked_sources(runtime::ClosedLoopRuntime) = supports_stacked_sources(runtime.wfs, runtime.src)
supports_grouped_execution(::CompositeSimulationInterface) = true

function prepare!(runtime::ClosedLoopRuntime)
    prepare_runtime_wfs!(runtime.wfs, runtime.tel, runtime.src)
    return runtime
end

function prepare!(interface::SimulationInterface)
    prepare!(interface.runtime)
    return snapshot_outputs!(interface)
end

function prepare!(interface::CompositeSimulationInterface)
    @inbounds for child in interface.interfaces
        prepare!(child.runtime)
    end
    return snapshot_outputs!(interface)
end

function SimulationInterface(runtime::ClosedLoopRuntime)
    command = similar(runtime.command)
    copyto!(command, runtime.command)
    slopes = similar(runtime.wfs.state.slopes)
    copyto!(slopes, runtime.wfs.state.slopes)
    wfs_frame = isnothing(runtime.wfs_detector) ? nothing : similar(wfs_output_frame(runtime.wfs, runtime.wfs_detector))
    science_frame = isnothing(runtime.science_detector) ? nothing : similar(output_frame(runtime.science_detector))
    interface = SimulationInterface{typeof(runtime), typeof(command), typeof(slopes), typeof(wfs_frame), typeof(science_frame)}(
        runtime,
        command,
        slopes,
        wfs_frame,
        science_frame,
    )
    return snapshot_outputs!(interface)
end

function CompositeSimulationInterface(interfaces::SimulationInterface...)
    length(interfaces) > 0 || throw(InvalidConfiguration("CompositeSimulationInterface requires at least one interface"))
    first_interface = interfaces[1]
    total_command = mapreduce(i -> length(simulation_command(i)), +, interfaces)
    total_slopes = mapreduce(i -> length(simulation_slopes(i)), +, interfaces)
    command = similar(simulation_command(first_interface), total_command)
    slopes = similar(simulation_slopes(first_interface), total_slopes)
    wfs_frames = ntuple(i -> begin
        frame = simulation_wfs_frame(interfaces[i])
        isnothing(frame) ? nothing : similar(frame)
    end, length(interfaces))
    science_frames = ntuple(i -> begin
        frame = simulation_science_frame(interfaces[i])
        isnothing(frame) ? nothing : similar(frame)
    end, length(interfaces))
    multi = CompositeSimulationInterface{typeof(interfaces), typeof(command), typeof(slopes), typeof(wfs_frames), typeof(science_frames)}(
        interfaces,
        command,
        slopes,
        wfs_frames,
        science_frames,
    )
    return snapshot_outputs!(multi)
end

CompositeSimulationInterface(runtimes::ClosedLoopRuntime...) = CompositeSimulationInterface(map(SimulationInterface, runtimes)...)

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

with_reconstructor(interface::SimulationInterface, reconstructor) = SimulationInterface(with_reconstructor(interface.runtime, reconstructor))

function with_reconstructors(interface::CompositeSimulationInterface, reconstructors::Vararg{Any,N}) where {N}
    length(interface.interfaces) == N ||
        throw(DimensionMismatchError("reconstructor count must match interface count"))
    refreshed = ntuple(i -> with_reconstructor(interface.interfaces[i], reconstructors[i]), N)
    return CompositeSimulationInterface(refreshed...)
end

@inline function set_command!(runtime::ClosedLoopRuntime, command::AbstractVector)
    copyto!(runtime.command, command)
    copyto!(runtime.dm.state.coefs, command)
    return runtime.dm.state.coefs
end

@inline function set_command!(interface::SimulationInterface, command::AbstractVector)
    copyto!(interface.command, command)
    set_command!(interface.runtime, command)
    return interface.command
end

@inline function set_command!(interface::CompositeSimulationInterface, command::AbstractVector)
    length(interface.command) == length(command) ||
        throw(DimensionMismatchError("command length must match aggregated RTC command length"))
    command_offset = 1
    @inbounds for child in interface.interfaces
        n_command = length(child.command)
        @views set_command!(child, command[command_offset:command_offset + n_command - 1])
        command_offset += n_command
    end
    copyto!(interface.command, command)
    return interface.command
end

"""
    snapshot_outputs!(interface)

Refresh the exported command, slopes, and frame buffers from the underlying
runtime state.

This is the explicit copy boundary between the internal mutable runtime objects
and the externally consumed `SimulationInterface` buffers.
"""
@inline function snapshot_outputs!(interface::SimulationInterface)
    copyto!(interface.command, interface.runtime.command)
    copyto!(interface.slopes, interface.runtime.wfs.state.slopes)
    if !isnothing(interface.wfs_frame)
        copyto!(interface.wfs_frame, wfs_output_frame(interface.runtime.wfs, interface.runtime.wfs_detector))
    end
    if !isnothing(interface.science_frame)
        copyto!(interface.science_frame, output_frame(interface.runtime.science_detector))
    end
    return interface
end

@inline function snapshot_outputs!(multi::CompositeSimulationInterface)
    command_offset = 1
    slope_offset = 1
    @inbounds for i in eachindex(multi.interfaces)
        interface = snapshot_outputs!(multi.interfaces[i])
        n_command = length(interface.command)
        n_slopes = length(interface.slopes)
        @views copyto!(multi.command[command_offset:command_offset + n_command - 1], interface.command)
        @views copyto!(multi.slopes[slope_offset:slope_offset + n_slopes - 1], interface.slopes)
        if !isnothing(multi.wfs_frames[i])
            copyto!(multi.wfs_frames[i], interface.wfs_frame)
        end
        if !isnothing(multi.science_frames[i])
            copyto!(multi.science_frames[i], interface.science_frame)
        end
        command_offset += n_command
        slope_offset += n_slopes
    end
    return multi
end

@inline simulation_command(readout::SimulationReadout) = readout.command
@inline simulation_slopes(readout::SimulationReadout) = readout.slopes
@inline simulation_wfs_frame(readout::SimulationReadout) = readout.wfs_frame
@inline simulation_science_frame(readout::SimulationReadout) = readout.science_frame
@inline simulation_wfs_metadata(readout::SimulationReadout) = readout.wfs_metadata
@inline simulation_science_metadata(readout::SimulationReadout) = readout.science_metadata

@inline simulation_command(runtime::ClosedLoopRuntime) = runtime.command
@inline simulation_slopes(runtime::ClosedLoopRuntime) = runtime.wfs.state.slopes
@inline simulation_wfs_frame(runtime::ClosedLoopRuntime) = isnothing(runtime.wfs_detector) ? nothing : wfs_output_frame(runtime.wfs, runtime.wfs_detector)
@inline simulation_science_frame(runtime::ClosedLoopRuntime) = isnothing(runtime.science_detector) ? nothing : output_frame(runtime.science_detector)
@inline simulation_wfs_metadata(runtime::ClosedLoopRuntime) = isnothing(runtime.wfs_detector) ? nothing : detector_export_metadata(runtime.wfs_detector)
@inline simulation_science_metadata(runtime::ClosedLoopRuntime) = isnothing(runtime.science_detector) ? nothing : detector_export_metadata(runtime.science_detector)

@inline simulation_command(sim::AbstractControlSimulation) = simulation_command(simulation_readout(sim))
@inline simulation_slopes(sim::AbstractControlSimulation) = simulation_slopes(simulation_readout(sim))
@inline simulation_wfs_frame(sim::AbstractControlSimulation) = simulation_wfs_frame(simulation_readout(sim))
@inline simulation_science_frame(sim::AbstractControlSimulation) = simulation_science_frame(simulation_readout(sim))
@inline simulation_wfs_metadata(sim::AbstractControlSimulation) = simulation_wfs_metadata(simulation_readout(sim))
@inline simulation_science_metadata(sim::AbstractControlSimulation) = simulation_science_metadata(simulation_readout(sim))

@inline simulation_command(interface::SimulationInterface) = interface.command
@inline simulation_slopes(interface::SimulationInterface) = interface.slopes
@inline simulation_wfs_frame(interface::SimulationInterface) = interface.wfs_frame
@inline simulation_science_frame(interface::SimulationInterface) = interface.science_frame
@inline simulation_wfs_metadata(interface::SimulationInterface) = isnothing(interface.runtime.wfs_detector) ? nothing : detector_export_metadata(interface.runtime.wfs_detector)
@inline simulation_science_metadata(interface::SimulationInterface) = isnothing(interface.runtime.science_detector) ? nothing : detector_export_metadata(interface.runtime.science_detector)

@inline simulation_command(interface::CompositeSimulationInterface) = interface.command
@inline simulation_slopes(interface::CompositeSimulationInterface) = interface.slopes
@inline simulation_wfs_frame(interface::CompositeSimulationInterface) = interface.wfs_frames
@inline simulation_science_frame(interface::CompositeSimulationInterface) = interface.science_frames
@inline simulation_wfs_metadata(interface::CompositeSimulationInterface) = map(simulation_wfs_metadata, interface.interfaces)
@inline simulation_science_metadata(interface::CompositeSimulationInterface) = map(simulation_science_metadata, interface.interfaces)

@inline simulation_readout(readout::SimulationReadout) = readout
@inline function simulation_readout(runtime::ClosedLoopRuntime)
    return SimulationReadout(
        simulation_command(runtime),
        simulation_slopes(runtime),
        simulation_wfs_frame(runtime),
        simulation_science_frame(runtime),
        simulation_wfs_metadata(runtime),
        simulation_science_metadata(runtime),
    )
end
@inline simulation_readout(sim::AbstractControlSimulation) = simulation_readout(simulation_interface(sim))

@inline function simulation_readout(interface::SimulationInterface)
    return SimulationReadout(
        simulation_command(interface),
        simulation_slopes(interface),
        simulation_wfs_frame(interface),
        simulation_science_frame(interface),
        simulation_wfs_metadata(interface),
        simulation_science_metadata(interface),
    )
end

@inline function simulation_readout(interface::CompositeSimulationInterface)
    return SimulationReadout(
        simulation_command(interface),
        simulation_slopes(interface),
        simulation_wfs_frame(interface),
        simulation_science_frame(interface),
        simulation_wfs_metadata(interface),
        simulation_science_metadata(interface),
    )
end

"""
    sense_core!(...)

Execute the sensing half of the closed-loop update.

This advances the atmosphere, propagates the telescope, applies the current DM
surface, and measures the WFS with or without an explicit detector model.
"""
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

"""
    step_core!(...)

Execute one full closed-loop update.

`step_core!` is the hot-path composition of sensing, optional science capture,
reconstruction, and DM command application. The various method overloads differ
only in whether WFS and science detectors participate in the step.
"""
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

function sense!(interface::SimulationInterface)
    sense!(interface.runtime)
    return snapshot_outputs!(interface)
end

function sense!(interface::CompositeSimulationInterface)
    @inbounds for child in interface.interfaces
        sense!(child.runtime)
    end
    return snapshot_outputs!(interface)
end

@inline function reconstruct!(runtime::ClosedLoopRuntime)
    reconstruct!(runtime.command, runtime.reconstructor, runtime.wfs.state.slopes)
    return runtime.command
end

@inline function apply_runtime_command!(runtime::ClosedLoopRuntime)
    apply_command!(runtime.dm.state.coefs, runtime.command, runtime.control_sign)
    return runtime.dm.state.coefs
end

@inline function step_grouped!(interface::CompositeSimulationInterface)
    @inbounds for child in interface.interfaces
        sense!(child.runtime)
    end
    @inbounds for child in interface.interfaces
        reconstruct!(child.runtime)
    end
    @inbounds for child in interface.interfaces
        apply_runtime_command!(child.runtime)
    end
    return snapshot_outputs!(interface)
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

function step!(interface::SimulationInterface)
    step!(interface.runtime)
    return snapshot_outputs!(interface)
end

"""
    step!(interface::CompositeSimulationInterface)

Run grouped multi-branch execution for a composite interface.

The grouped path performs all sensing phases first, then all reconstruction
phases, then all command-application phases, which reduces unnecessary
inter-branch serialization compared with stepping each child end to end.
"""
function step!(interface::CompositeSimulationInterface)
    return step_grouped!(interface)
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
    synchronize_backend!(execution_style(runtime.wfs.state.slopes))
    if !isnothing(runtime.wfs_detector)
        synchronize_backend!(execution_style(output_frame(runtime.wfs_detector)))
    end
    if !isnothing(runtime.science_detector)
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

function runtime_phase_timing(interface::SimulationInterface; warmup::Int=10, samples::Int=1000, gc_before::Bool=true)
    warmup >= 0 || throw(InvalidConfiguration("warmup must be >= 0"))
    samples > 0 || throw(InvalidConfiguration("samples must be positive"))
    for _ in 1:warmup
        step!(interface)
    end
    gc_before && GC.gc()
    sense_times = Vector{Int}(undef, samples)
    reconstruct_times = Vector{Int}(undef, samples)
    apply_times = Vector{Int}(undef, samples)
    snapshot_times = Vector{Int}(undef, samples)
    total_times = Vector{Int}(undef, samples)
    runtime = interface.runtime
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
        snapshot_outputs!(interface)
        synchronize_backend!(execution_style(interface.command))
        synchronize_backend!(execution_style(interface.slopes))
        if !isnothing(interface.wfs_frame)
            synchronize_backend!(execution_style(interface.wfs_frame))
        end
        if !isnothing(interface.science_frame)
            synchronize_backend!(execution_style(interface.science_frame))
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

function runtime_phase_timing(interface::CompositeSimulationInterface; warmup::Int=10, samples::Int=1000, gc_before::Bool=true)
    warmup >= 0 || throw(InvalidConfiguration("warmup must be >= 0"))
    samples > 0 || throw(InvalidConfiguration("samples must be positive"))
    for _ in 1:warmup
        step!(interface)
    end
    gc_before && GC.gc()
    sense_times = Vector{Int}(undef, samples)
    reconstruct_times = Vector{Int}(undef, samples)
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
            reconstruct!(child.runtime.command, child.runtime.reconstructor, child.runtime.wfs.state.slopes)
            synchronize_backend!(execution_style(child.runtime.command))
        end
        t2 = time_ns()
        for child in interface.interfaces
            apply_command!(child.runtime.dm.state.coefs, child.runtime.command, child.runtime.control_sign)
            synchronize_backend!(execution_style(child.runtime.dm.state.coefs))
        end
        t3 = time_ns()
        snapshot_outputs!(interface)
        synchronize_backend!(execution_style(interface.command))
        synchronize_backend!(execution_style(interface.slopes))
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

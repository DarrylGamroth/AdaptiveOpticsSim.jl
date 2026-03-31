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
    copyto!(interface.slopes, interface.runtime.slopes)
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
@inline simulation_slopes(runtime::ClosedLoopRuntime) = requires_runtime_slopes(runtime) ? runtime.slopes : nothing
@inline simulation_wfs_frame(runtime::ClosedLoopRuntime) = requires_runtime_wfs_pixels(runtime) ? wfs_output_frame(runtime.wfs, runtime.wfs_detector) : nothing
@inline simulation_science_frame(runtime::ClosedLoopRuntime) = requires_runtime_science_pixels(runtime) ? output_frame(runtime.science_detector) : nothing
@inline simulation_wfs_metadata(runtime::ClosedLoopRuntime) =
    requires_runtime_wfs_pixels(runtime) ? detector_export_metadata(runtime.wfs_detector) : nothing
@inline simulation_science_metadata(runtime::ClosedLoopRuntime) =
    requires_runtime_science_pixels(runtime) ? detector_export_metadata(runtime.science_detector) : nothing
@inline runtime_profile(runtime::ClosedLoopRuntime) = runtime.profile
@inline runtime_latency(runtime::ClosedLoopRuntime) = runtime.latency

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
@inline simulation_wfs_metadata(interface::SimulationInterface) = simulation_wfs_metadata(interface.runtime)
@inline simulation_science_metadata(interface::SimulationInterface) = simulation_science_metadata(interface.runtime)

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

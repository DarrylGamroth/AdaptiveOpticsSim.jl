abstract type AbstractRuntimeExportPlan end

struct DirectRuntimeExportPlan <: AbstractRuntimeExportPlan end
struct CompositeRuntimeExportPlan <: AbstractRuntimeExportPlan end

@inline runtime_export_plan(::SimulationInterface) = DirectRuntimeExportPlan()
@inline runtime_export_plan(::CompositeSimulationInterface) = CompositeRuntimeExportPlan()

@inline command_layout(runtime::ClosedLoopRuntime) = runtime.command_layout
@inline command_layout(interface::SimulationInterface) = command_layout(interface.runtime)

@inline _runtime_layout_matches_optic(runtime::ClosedLoopRuntime) = runtime.command_layout == command_layout(runtime.optic)

@inline function set_command!(runtime::ClosedLoopRuntime, command::AbstractVector)
    length(runtime.command) == length(command) ||
        throw(DimensionMismatchError("command length must match runtime command length"))
    copyto!(runtime.command, command)
    set_command!(runtime.optic, runtime.command)
    return command_storage(runtime.optic)
end

@inline function set_command!(runtime::ClosedLoopRuntime, command::NamedTuple)
    if _runtime_layout_matches_optic(runtime)
        set_command!(runtime.optic, command)
        copyto!(runtime.command, command_storage(runtime.optic))
    else
        _set_command_segments!(length(runtime.command), command_layout(runtime), command) do seg, segment
            rng = command_segment_range(seg)
            copyto!(@view(runtime.command[rng]), segment)
        end
        set_command!(runtime.optic, runtime.command)
    end
    return runtime.command
end

@inline function set_command!(interface::SimulationInterface, command::AbstractVector)
    length(interface.command) == length(command) ||
        throw(DimensionMismatchError("command length must match exported simulation interface command length"))
    copyto!(interface.command, command)
    set_command!(interface.runtime, command)
    return interface.command
end

@inline function set_command!(interface::SimulationInterface, command::NamedTuple)
    set_command!(interface.runtime, command)
    copyto!(interface.command, interface.runtime.command)
    return interface.command
end

@inline function update_command!(runtime::ClosedLoopRuntime, command::NamedTuple)
    if _runtime_layout_matches_optic(runtime)
        update_command!(runtime.optic, command)
        copyto!(runtime.command, command_storage(runtime.optic))
    else
        _set_command_segments!(length(runtime.command), runtime.command_layout, command; require_all=false) do seg, segment
            rng = command_segment_range(seg)
            copyto!(@view(runtime.command[rng]), segment)
        end
        set_command!(runtime.optic, runtime.command)
    end
    return runtime.command
end

@inline function update_command!(interface::SimulationInterface, command::NamedTuple)
    update_command!(interface.runtime, command)
    copyto!(interface.command, interface.runtime.command)
    return interface.command
end

@inline function update_command!(interface::CompositeSimulationInterface, command::NamedTuple)
    throw(InvalidConfiguration("update_command! on a composite interface should be routed through a RuntimeScenario or child SimulationInterface"))
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
@inline snapshot_outputs!(interface::SimulationInterface) = snapshot_outputs!(runtime_export_plan(interface), interface)
@inline snapshot_outputs!(multi::CompositeSimulationInterface) = snapshot_outputs!(runtime_export_plan(multi), multi)

@inline function snapshot_outputs!(::DirectRuntimeExportPlan, interface::SimulationInterface)
    copyto!(interface.command, interface.runtime.command)
    copyto!(interface.slopes, interface.runtime.slopes)
    if !isnothing(interface.wfs_frame)
        source = wfs_output_frame(interface.runtime.wfs, interface.runtime.wfs_detector)
        if size(interface.wfs_frame) != size(source)
            interface.wfs_frame = similar(source)
        end
        copyto!(interface.wfs_frame, source)
    end
    if !isnothing(interface.science_frame)
        source = output_frame(interface.runtime.science_detector)
        if size(interface.science_frame) != size(source)
            interface.science_frame = similar(source)
        end
        copyto!(interface.science_frame, source)
    end
    return interface
end

@inline function snapshot_outputs!(::CompositeRuntimeExportPlan, multi::CompositeSimulationInterface)
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
        if !isnothing(multi.wfs_stack)
            copyto!(selectdim(multi.wfs_stack, ndims(multi.wfs_stack), i), interface.wfs_frame)
        end
        if !isnothing(multi.science_stack)
            copyto!(selectdim(multi.science_stack, ndims(multi.science_stack), i), interface.science_frame)
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
@inline simulation_grouped_wfs_stack(readout::SimulationReadout) = readout.grouped_wfs_stack
@inline simulation_grouped_science_stack(readout::SimulationReadout) = readout.grouped_science_stack

@inline simulation_command(runtime::ClosedLoopRuntime) = runtime.command
@inline simulation_slopes(runtime::ClosedLoopRuntime) = requires_runtime_slopes(runtime) ? runtime.slopes : nothing
@inline simulation_wfs_frame(runtime::ClosedLoopRuntime) = requires_runtime_wfs_pixels(runtime) ? wfs_output_frame(runtime.wfs, runtime.wfs_detector) : nothing
@inline simulation_science_frame(runtime::ClosedLoopRuntime) = requires_runtime_science_pixels(runtime) ? output_frame(runtime.science_detector) : nothing
@inline simulation_wfs_metadata(runtime::ClosedLoopRuntime) =
    requires_runtime_wfs_pixels(runtime) ? detector_export_metadata(runtime.wfs_detector) : nothing
@inline simulation_science_metadata(runtime::ClosedLoopRuntime) =
    requires_runtime_science_pixels(runtime) ? detector_export_metadata(runtime.science_detector) : nothing
@inline simulation_grouped_wfs_stack(::ClosedLoopRuntime) = nothing
@inline simulation_grouped_science_stack(::ClosedLoopRuntime) = nothing
@inline runtime_profile(runtime::ClosedLoopRuntime) = runtime.profile
@inline runtime_latency(runtime::ClosedLoopRuntime) = runtime.latency

@inline simulation_command(sim::AbstractControlSimulation) = simulation_command(simulation_readout(sim))
@inline simulation_slopes(sim::AbstractControlSimulation) = simulation_slopes(simulation_readout(sim))
@inline simulation_wfs_frame(sim::AbstractControlSimulation) = simulation_wfs_frame(simulation_readout(sim))
@inline simulation_science_frame(sim::AbstractControlSimulation) = simulation_science_frame(simulation_readout(sim))
@inline simulation_wfs_metadata(sim::AbstractControlSimulation) = simulation_wfs_metadata(simulation_readout(sim))
@inline simulation_science_metadata(sim::AbstractControlSimulation) = simulation_science_metadata(simulation_readout(sim))
@inline simulation_grouped_wfs_stack(sim::AbstractControlSimulation) = simulation_grouped_wfs_stack(simulation_readout(sim))
@inline simulation_grouped_science_stack(sim::AbstractControlSimulation) = simulation_grouped_science_stack(simulation_readout(sim))

@inline simulation_command(interface::SimulationInterface) = interface.command
@inline simulation_slopes(interface::SimulationInterface) = interface.slopes
@inline simulation_wfs_frame(interface::SimulationInterface) = interface.wfs_frame
@inline simulation_science_frame(interface::SimulationInterface) = interface.science_frame
@inline simulation_wfs_metadata(interface::SimulationInterface) = simulation_wfs_metadata(interface.runtime)
@inline simulation_science_metadata(interface::SimulationInterface) = simulation_science_metadata(interface.runtime)
@inline simulation_grouped_wfs_stack(::SimulationInterface) = nothing
@inline simulation_grouped_science_stack(::SimulationInterface) = nothing

@inline simulation_command(interface::CompositeSimulationInterface) = interface.command
@inline simulation_slopes(interface::CompositeSimulationInterface) = interface.slopes
@inline simulation_wfs_frame(interface::CompositeSimulationInterface) = requires_grouped_wfs_frames(interface) ? interface.wfs_frames : nothing
@inline simulation_science_frame(interface::CompositeSimulationInterface) = requires_grouped_science_frames(interface) ? interface.science_frames : nothing
@inline simulation_wfs_metadata(interface::CompositeSimulationInterface) = map(simulation_wfs_metadata, interface.interfaces)
@inline simulation_science_metadata(interface::CompositeSimulationInterface) = map(simulation_science_metadata, interface.interfaces)
@inline simulation_grouped_wfs_stack(interface::CompositeSimulationInterface) = requires_grouped_wfs_stack(interface) ? interface.wfs_stack : nothing
@inline simulation_grouped_science_stack(interface::CompositeSimulationInterface) = requires_grouped_science_stack(interface) ? interface.science_stack : nothing

@inline simulation_readout(readout::SimulationReadout) = readout

@inline command(x) = simulation_command(x)
@inline slopes(x) = simulation_slopes(x)
@inline wfs_frame(x) = simulation_wfs_frame(x)
@inline science_frame(x) = simulation_science_frame(x)
@inline wfs_metadata(x) = simulation_wfs_metadata(x)
@inline science_metadata(x) = simulation_science_metadata(x)
@inline grouped_wfs_stack(x) = simulation_grouped_wfs_stack(x)
@inline grouped_science_stack(x) = simulation_grouped_science_stack(x)

@inline function simulation_readout(runtime::ClosedLoopRuntime)
    return SimulationReadout(
        simulation_command(runtime),
        simulation_slopes(runtime),
        simulation_wfs_frame(runtime),
        simulation_science_frame(runtime),
        simulation_wfs_metadata(runtime),
        simulation_science_metadata(runtime),
        simulation_grouped_wfs_stack(runtime),
        simulation_grouped_science_stack(runtime),
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
        simulation_grouped_wfs_stack(interface),
        simulation_grouped_science_stack(interface),
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
        simulation_grouped_wfs_stack(interface),
        simulation_grouped_science_stack(interface),
    )
end

"""
    sense_core!(...)

Execute the sensing half of the closed-loop update.

This advances the atmosphere, propagates the telescope, applies the current DM
surface, and measures the WFS with or without an explicit detector model.
"""

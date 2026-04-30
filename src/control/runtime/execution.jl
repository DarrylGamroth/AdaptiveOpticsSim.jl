@inline function sense_core!(atm::AbstractAtmosphere, tel::Telescope, optic::AbstractControllableOptic,
    wfs::AbstractWFS, src::AbstractSource, rng::AbstractRNG)
    prepropagate_runtime_wfs!(wfs, atm, tel, optic, src, rng)
    measure_runtime_wfs!(wfs, atm, tel, src, rng)
    finish_runtime_wfs_sensing!(wfs, atm, tel, src)
    return nothing
end

@inline function sense_core!(atm::AbstractAtmosphere, tel::Telescope, optic::AbstractControllableOptic,
    wfs::AbstractWFS, src::AbstractSource, det::AbstractDetector, rng::AbstractRNG)
    prepropagate_runtime_wfs!(wfs, atm, tel, optic, src, rng)
    measure_runtime_wfs!(wfs, atm, tel, src, det, rng)
    finish_runtime_wfs_sensing!(wfs, atm, tel, src)
    return nothing
end

@inline function capture_science_core!(tel::Telescope, src::AbstractSource,
    det::AbstractDetector, rng::AbstractRNG, zero_padding::Int)
    psf = compute_psf!(tel, src; zero_padding=zero_padding)
    capture!(det, psf; rng=rng)
    return nothing
end

@inline function stage_sensed_slopes!(runtime::ClosedLoopRuntime)
    synchronize_backend!(execution_style(slopes(runtime.wfs)))
    measured = shift_delay!(runtime.measurement_delay, slopes(runtime.wfs))
    synchronize_backend!(execution_style(measured))
    readout = shift_delay!(runtime.readout_delay, measured)
    synchronize_backend!(execution_style(readout))
    copyto!(runtime.slopes, readout)
    synchronize_backend!(execution_style(runtime.slopes))
    return runtime.slopes
end

@inline function stage_runtime_command!(runtime::ClosedLoopRuntime)
    synchronize_backend!(execution_style(runtime.reconstruct_buffer))
    delayed = shift_delay!(runtime.reconstruction_delay, runtime.reconstruct_buffer)
    synchronize_backend!(execution_style(delayed))
    copyto!(runtime.command, delayed)
    synchronize_backend!(execution_style(runtime.command))
    return runtime.command
end

"""
    step_core!(...)

Execute one full closed-loop update.

`step_core!` is the hot-path composition of sensing, optional science capture,
reconstruction, and DM command application. The various method overloads differ
only in whether WFS and science detectors participate in the step.
"""
@inline function step_core!(atm::AbstractAtmosphere, tel::Telescope, optic::AbstractControllableOptic,
    wfs::AbstractWFS, src::AbstractSource, reconstructor, command::AbstractVector{T},
    rng::AbstractRNG, control_sign::T) where {T<:AbstractFloat}
    sense_core!(atm, tel, optic, wfs, src, rng)
    reconstruct!(command, reconstructor, slopes(wfs))
    _stage_command!(optic, command, control_sign)
    return nothing
end

@inline function step_core!(atm::AbstractAtmosphere, tel::Telescope, optic::AbstractControllableOptic,
    wfs::AbstractWFS, src::AbstractSource, wfs_detector::AbstractDetector, reconstructor,
    command::AbstractVector{T}, rng::AbstractRNG, control_sign::T) where {T<:AbstractFloat}
    sense_core!(atm, tel, optic, wfs, src, wfs_detector, rng)
    reconstruct!(command, reconstructor, slopes(wfs))
    _stage_command!(optic, command, control_sign)
    return nothing
end

@inline function step_core!(atm::AbstractAtmosphere, tel::Telescope, optic::AbstractControllableOptic,
    wfs::AbstractWFS, src::AbstractSource, reconstructor, command::AbstractVector{T},
    science_detector::AbstractDetector, rng::AbstractRNG, control_sign::T, science_zero_padding::Int) where {T<:AbstractFloat}
    sense_core!(atm, tel, optic, wfs, src, rng)
    capture_science_core!(tel, src, science_detector, rng, science_zero_padding)
    reconstruct!(command, reconstructor, slopes(wfs))
    _stage_command!(optic, command, control_sign)
    return nothing
end

@inline function step_core!(atm::AbstractAtmosphere, tel::Telescope, optic::AbstractControllableOptic,
    wfs::AbstractWFS, src::AbstractSource, wfs_detector::AbstractDetector, reconstructor,
    command::AbstractVector{T}, science_detector::AbstractDetector, rng::AbstractRNG,
    control_sign::T, science_zero_padding::Int) where {T<:AbstractFloat}
    sense_core!(atm, tel, optic, wfs, src, wfs_detector, rng)
    capture_science_core!(tel, src, science_detector, rng, science_zero_padding)
    reconstruct!(command, reconstructor, slopes(wfs))
    _stage_command!(optic, command, control_sign)
    return nothing
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
    reconstruct!(runtime.reconstruct_buffer, runtime.reconstructor, runtime.slopes)
    stage_runtime_command!(runtime)
    return runtime.command
end

@inline function apply_runtime_command!(runtime::ClosedLoopRuntime)
    delayed = shift_delay!(runtime.dm_delay, runtime.command)
    _stage_command!(runtime.optic, delayed, runtime.control_sign)
    return command_storage(runtime.optic)
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

function sense!(runtime::ClosedLoopRuntime)
    if isnothing(runtime.wfs_detector)
        sense_core!(runtime.atm, runtime.tel, runtime.optic, runtime.wfs, runtime.src, runtime.rng)
    else
        sense_core!(runtime.atm, runtime.tel, runtime.optic, runtime.wfs, runtime.src, runtime.wfs_detector, runtime.rng)
    end
    if requires_runtime_science_pixels(runtime)
        capture_science_core!(runtime.tel, runtime.src, runtime.science_detector, runtime.rng, runtime.science_zero_padding)
    end
    if requires_runtime_slopes(runtime)
        stage_sensed_slopes!(runtime)
    end
    return runtime
end

function step!(runtime::ClosedLoopRuntime)
    sense!(runtime)
    reconstruct!(runtime)
    apply_runtime_command!(runtime)
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

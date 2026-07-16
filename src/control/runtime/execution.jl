@inline function sense_core!(atm::AbstractAtmosphere, tel::Telescope, optic::AbstractControllableOptic,
    wfs::AbstractWFS, src::AbstractSource, atmosphere_step::Real,
    rng::AbstractRNG)
    prepropagate_runtime_wfs!(wfs, atm, tel, optic, src, atmosphere_step, rng)
    measure_runtime_wfs!(wfs, atm, tel, src, rng)
    finish_runtime_wfs_sensing!(wfs, atm, tel, optic, src)
    return nothing
end

@inline function sense_core_after_advance!(atm::AbstractAtmosphere,
    tel::Telescope, optic::AbstractControllableOptic, wfs::AbstractWFS,
    src::AbstractSource, rng::AbstractRNG)
    render_runtime_wfs_path!(wfs, atm, tel, optic, src)
    measure_runtime_wfs!(wfs, atm, tel, src, rng)
    finish_runtime_wfs_sensing!(wfs, atm, tel, optic, src)
    return nothing
end

@inline function sense_core!(atm::AbstractAtmosphere, tel::Telescope,
    optic::AbstractControllableOptic, wfs::AbstractWFS,
    src::AbstractSource, atmosphere_step::Real, rng::AbstractRNG,
    atmosphere_renderer)
    advance_runtime_atmosphere!(wfs, atm, tel, atmosphere_step, rng)
    render_runtime_wfs_path!(atmosphere_renderer, wfs, atm, tel, optic, src)
    measure_runtime_wfs!(wfs, atm, tel, src, rng)
    finish_runtime_wfs_sensing!(wfs, atm, tel, optic, src)
    return nothing
end

@inline function sense_core_after_advance!(atm::AbstractAtmosphere,
    tel::Telescope, optic::AbstractControllableOptic, wfs::AbstractWFS,
    src::AbstractSource, rng::AbstractRNG, atmosphere_renderer)
    render_runtime_wfs_path!(atmosphere_renderer, wfs, atm, tel, optic, src)
    measure_runtime_wfs!(wfs, atm, tel, src, rng)
    finish_runtime_wfs_sensing!(wfs, atm, tel, optic, src)
    return nothing
end

@inline function sense_core!(atm::AbstractAtmosphere, tel::Telescope, optic::AbstractControllableOptic,
    wfs::AbstractWFS, src::AbstractSource, det::AbstractDetector,
    atmosphere_step::Real, rng::AbstractRNG)
    prepropagate_runtime_wfs!(wfs, atm, tel, optic, src, atmosphere_step, rng)
    measure_runtime_wfs!(wfs, atm, tel, src, det, rng)
    finish_runtime_wfs_sensing!(wfs, atm, tel, optic, src)
    return nothing
end

@inline function sense_core_after_advance!(atm::AbstractAtmosphere,
    tel::Telescope, optic::AbstractControllableOptic, wfs::AbstractWFS,
    src::AbstractSource, det::AbstractDetector, rng::AbstractRNG)
    render_runtime_wfs_path!(wfs, atm, tel, optic, src)
    measure_runtime_wfs!(wfs, atm, tel, src, det, rng)
    finish_runtime_wfs_sensing!(wfs, atm, tel, optic, src)
    return nothing
end

@inline function sense_core!(atm::AbstractAtmosphere, tel::Telescope,
    optic::AbstractControllableOptic, wfs::AbstractWFS,
    src::AbstractSource, det::AbstractDetector, atmosphere_step::Real,
    rng::AbstractRNG,
    atmosphere_renderer)
    advance_runtime_atmosphere!(wfs, atm, tel, atmosphere_step, rng)
    render_runtime_wfs_path!(atmosphere_renderer, wfs, atm, tel, optic, src)
    measure_runtime_wfs!(wfs, atm, tel, src, det, rng)
    finish_runtime_wfs_sensing!(wfs, atm, tel, optic, src)
    return nothing
end

@inline function sense_core_after_advance!(atm::AbstractAtmosphere,
    tel::Telescope, optic::AbstractControllableOptic, wfs::AbstractWFS,
    src::AbstractSource, det::AbstractDetector, rng::AbstractRNG,
    atmosphere_renderer)
    render_runtime_wfs_path!(atmosphere_renderer, wfs, atm, tel, optic, src)
    measure_runtime_wfs!(wfs, atm, tel, src, det, rng)
    finish_runtime_wfs_sensing!(wfs, atm, tel, optic, src)
    return nothing
end

@inline function capture_science_core!(tel::Telescope, src::AbstractSource,
    det::AbstractDetector, rng::AbstractRNG, zero_padding::Int)
    psf = compute_psf!(tel, src; zero_padding=zero_padding)
    capture!(det, psf, src; rng=rng)
    return nothing
end

@inline prepare_science_path!(::ReuseSensedOpticalPath, atm::AbstractAtmosphere,
    tel::Telescope, optic::AbstractControllableOptic, src::AbstractSource) = tel

@inline function prepare_science_path!(::RepropagateScienceOpticalPath,
    atm::AbstractAtmosphere, tel::Telescope, optic::AbstractControllableOptic,
    src::AbstractSource)
    propagate_runtime_atmosphere!(atm, tel, src)
    apply!(optic, tel, DMAdditive())
    return tel
end

@inline prepare_science_path!(::ReuseSensedOpticalPath,
    atmosphere_renderer, atm::AbstractAtmosphere, tel::Telescope,
    optic::AbstractControllableOptic, src::AbstractSource) = tel

@inline function prepare_science_path!(::RepropagateScienceOpticalPath,
    atmosphere_renderer, atm::AbstractAtmosphere, tel::Telescope,
    optic::AbstractControllableOptic, src::AbstractSource)
    render_prepared_atmosphere_path!(atmosphere_renderer, atm, tel, src)
    apply!(optic, tel, DMAdditive())
    return tel
end

@inline function capture_science_core!(plan::AbstractSciencePathPlan,
    atm::AbstractAtmosphere, tel::Telescope, optic::AbstractControllableOptic,
    src::AbstractSource, det::AbstractDetector, rng::AbstractRNG,
    zero_padding::Int)
    prepare_science_path!(plan, atm, tel, optic, src)
    return capture_science_core!(tel, src, det, rng, zero_padding)
end

@inline function capture_science_core!(plan::AbstractSciencePathPlan,
    atmosphere_renderer, atm::AbstractAtmosphere, tel::Telescope,
    optic::AbstractControllableOptic, src::AbstractSource,
    det::AbstractDetector, rng::AbstractRNG, zero_padding::Int)
    prepare_science_path!(plan, atmosphere_renderer, atm, tel, optic, src)
    return capture_science_core!(tel, src, det, rng, zero_padding)
end

@inline _runtime_staging_barrier!(::CPUHILExecutionPlan, array) =
    synchronize_backend!(execution_style(array))
@inline _runtime_staging_barrier!(::DeviceResidentExecutionPlan, array) = nothing

@inline function stage_sensed_slopes!(runtime::ClosedLoopRuntime)
    plan = runtime.execution_plan
    _runtime_staging_barrier!(plan, slopes(runtime.wfs))
    measured = shift_delay!(runtime.measurement_delay, slopes(runtime.wfs))
    _runtime_staging_barrier!(plan, measured)
    readout = shift_delay!(runtime.readout_delay, measured)
    _runtime_staging_barrier!(plan, readout)
    copyto!(runtime.slopes, readout)
    _runtime_staging_barrier!(plan, runtime.slopes)
    return runtime.slopes
end

@inline function stage_runtime_command!(runtime::ClosedLoopRuntime)
    plan = runtime.execution_plan
    _runtime_staging_barrier!(plan, runtime.reconstruct_buffer)
    delayed = shift_delay!(runtime.reconstruction_delay, runtime.reconstruct_buffer)
    _runtime_staging_barrier!(plan, delayed)
    copyto!(runtime.command, delayed)
    _runtime_staging_barrier!(plan, runtime.command)
    return runtime.command
end

"""
    synchronize_runtime!(runtime)

Wait for queued backend work at an explicit observation boundary. Repeated
`step!` calls under `DeviceResidentExecutionPlan` do not require this barrier;
call it before host observation or timing when no exporting interface already
provides the boundary.
"""
function synchronize_runtime!(runtime::ClosedLoopRuntime)
    synchronize_backend!(execution_style(runtime.command))
    return runtime
end

function synchronize_runtime!(interface::SimulationInterface)
    synchronize_runtime!(interface.runtime)
    return interface
end

function synchronize_runtime!(interface::CompositeSimulationInterface)
    @inbounds for child in interface.interfaces
        synchronize_runtime!(child.runtime)
    end
    return interface
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
    atmosphere_step::Real, rng::AbstractRNG,
    control_sign::T) where {T<:AbstractFloat}
    sense_core!(atm, tel, optic, wfs, src, atmosphere_step, rng)
    reconstruct!(command, reconstructor, slopes(wfs))
    _stage_command!(optic, command, control_sign)
    return nothing
end

@inline function step_core!(atm::AbstractAtmosphere, tel::Telescope, optic::AbstractControllableOptic,
    wfs::AbstractWFS, src::AbstractSource, wfs_detector::AbstractDetector, reconstructor,
    command::AbstractVector{T}, atmosphere_step::Real, rng::AbstractRNG,
    control_sign::T) where {T<:AbstractFloat}
    sense_core!(atm, tel, optic, wfs, src, wfs_detector, atmosphere_step, rng)
    reconstruct!(command, reconstructor, slopes(wfs))
    _stage_command!(optic, command, control_sign)
    return nothing
end

@inline function step_core!(atm::AbstractAtmosphere, tel::Telescope, optic::AbstractControllableOptic,
    wfs::AbstractWFS, src::AbstractSource, reconstructor, command::AbstractVector{T},
    science_detector::AbstractDetector, atmosphere_step::Real,
    rng::AbstractRNG, control_sign::T,
    science_zero_padding::Int) where {T<:AbstractFloat}
    sense_core!(atm, tel, optic, wfs, src, atmosphere_step, rng)
    capture_science_core!(tel, src, science_detector, rng, science_zero_padding)
    reconstruct!(command, reconstructor, slopes(wfs))
    _stage_command!(optic, command, control_sign)
    return nothing
end

@inline function step_core!(atm::AbstractAtmosphere, tel::Telescope, optic::AbstractControllableOptic,
    wfs::AbstractWFS, src::AbstractSource, wfs_detector::AbstractDetector, reconstructor,
    command::AbstractVector{T}, science_detector::AbstractDetector,
    atmosphere_step::Real, rng::AbstractRNG,
    control_sign::T, science_zero_padding::Int) where {T<:AbstractFloat}
    sense_core!(atm, tel, optic, wfs, src, wfs_detector, atmosphere_step, rng)
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

abstract type AbstractRuntimeAtmosphereStep end
struct AdvanceRuntimeAtmosphere <: AbstractRuntimeAtmosphereStep end
struct ReuseAdvancedRuntimeAtmosphere <: AbstractRuntimeAtmosphereStep end

@inline function sense_primary_core!(::AdvanceRuntimeAtmosphere,
    runtime::ClosedLoopRuntime)
    if isnothing(runtime.wfs_detector)
        sense_core!(runtime.atm, runtime.tel, runtime.optic, runtime.wfs,
            runtime.src, runtime.atmosphere_step, runtime.rng,
            runtime.wfs_atmosphere_renderer)
    else
        sense_core!(runtime.atm, runtime.tel, runtime.optic, runtime.wfs,
            runtime.src, runtime.wfs_detector, runtime.atmosphere_step,
            runtime.rng,
            runtime.wfs_atmosphere_renderer)
    end
    return nothing
end

@inline function sense_primary_core!(::ReuseAdvancedRuntimeAtmosphere,
    runtime::ClosedLoopRuntime)
    if isnothing(runtime.wfs_detector)
        sense_core_after_advance!(runtime.atm, runtime.tel, runtime.optic,
            runtime.wfs, runtime.src, runtime.rng,
            runtime.wfs_atmosphere_renderer)
    else
        sense_core_after_advance!(runtime.atm, runtime.tel, runtime.optic,
            runtime.wfs, runtime.src, runtime.wfs_detector, runtime.rng,
            runtime.wfs_atmosphere_renderer)
    end
    return nothing
end

@inline function sense_runtime!(advance_policy::AbstractRuntimeAtmosphereStep,
    runtime::ClosedLoopRuntime)
    sense_primary_core!(advance_policy, runtime)
    if requires_runtime_science_pixels(runtime)
        capture_science_core!(runtime.science_path,
            runtime.science_atmosphere_renderer, runtime.atm, runtime.tel,
            runtime.optic, runtime.science_src, runtime.science_detector,
            runtime.rng, runtime.science_zero_padding)
    end
    if requires_runtime_slopes(runtime)
        stage_sensed_slopes!(runtime)
    end
    return runtime
end


function sense!(runtime::ClosedLoopRuntime)
    return sense_runtime!(AdvanceRuntimeAtmosphere(), runtime)
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

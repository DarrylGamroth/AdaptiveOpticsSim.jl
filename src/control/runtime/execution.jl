@inline function sense_core!(atm::AbstractAtmosphere, tel::Telescope,
    pupil::PupilFunction, optic::AbstractControllableOptic,
    wfs::AbstractWFS, src::AbstractSource, atmosphere_step::Real,
    rng::AbstractRNG)
    prepropagate_runtime_wfs!(wfs, atm, tel, pupil, optic, src,
        atmosphere_step, rng)
    measure_runtime_wfs!(wfs, atm, pupil, src, rng)
    finish_runtime_wfs_sensing!(wfs, atm, pupil, optic, src)
    return nothing
end

@inline function sense_core_after_advance!(atm::AbstractAtmosphere,
    pupil::PupilFunction, optic::AbstractControllableOptic, wfs::AbstractWFS,
    src::AbstractSource, rng::AbstractRNG)
    render_runtime_wfs_path!(wfs, atm, pupil, optic, src)
    measure_runtime_wfs!(wfs, atm, pupil, src, rng)
    finish_runtime_wfs_sensing!(wfs, atm, pupil, optic, src)
    return nothing
end

@inline function sense_core!(atm::AbstractAtmosphere, tel::Telescope,
    pupil::PupilFunction, optic::AbstractControllableOptic, wfs::AbstractWFS,
    src::AbstractSource, atmosphere_step::Real, rng::AbstractRNG,
    atmosphere_renderer)
    advance_runtime_atmosphere!(wfs, atm, tel, atmosphere_step, rng)
    render_runtime_wfs_path!(atmosphere_renderer, wfs, atm, pupil, optic, src)
    measure_runtime_wfs!(wfs, atm, pupil, src, rng)
    finish_runtime_wfs_sensing!(wfs, atm, pupil, optic, src)
    return nothing
end

@inline function sense_core_after_advance!(atm::AbstractAtmosphere,
    pupil::PupilFunction, optic::AbstractControllableOptic, wfs::AbstractWFS,
    src::AbstractSource, rng::AbstractRNG, atmosphere_renderer)
    render_runtime_wfs_path!(atmosphere_renderer, wfs, atm, pupil, optic, src)
    measure_runtime_wfs!(wfs, atm, pupil, src, rng)
    finish_runtime_wfs_sensing!(wfs, atm, pupil, optic, src)
    return nothing
end

@inline function sense_core!(atm::AbstractAtmosphere, tel::Telescope,
    pupil::PupilFunction, optic::AbstractControllableOptic,
    wfs::AbstractWFS, src::AbstractSource, det::AbstractDetector,
    atmosphere_step::Real, rng::AbstractRNG)
    prepropagate_runtime_wfs!(wfs, atm, tel, pupil, optic, src,
        atmosphere_step, rng)
    measure_runtime_wfs!(wfs, atm, pupil, src, det, rng)
    finish_runtime_wfs_sensing!(wfs, atm, pupil, optic, src)
    return nothing
end

@inline function sense_core_after_advance!(atm::AbstractAtmosphere,
    pupil::PupilFunction, optic::AbstractControllableOptic, wfs::AbstractWFS,
    src::AbstractSource, det::AbstractDetector, rng::AbstractRNG)
    render_runtime_wfs_path!(wfs, atm, pupil, optic, src)
    measure_runtime_wfs!(wfs, atm, pupil, src, det, rng)
    finish_runtime_wfs_sensing!(wfs, atm, pupil, optic, src)
    return nothing
end

@inline function sense_core!(atm::AbstractAtmosphere, tel::Telescope,
    pupil::PupilFunction, optic::AbstractControllableOptic, wfs::AbstractWFS,
    src::AbstractSource, det::AbstractDetector, atmosphere_step::Real,
    rng::AbstractRNG,
    atmosphere_renderer)
    advance_runtime_atmosphere!(wfs, atm, tel, atmosphere_step, rng)
    render_runtime_wfs_path!(atmosphere_renderer, wfs, atm, pupil, optic, src)
    measure_runtime_wfs!(wfs, atm, pupil, src, det, rng)
    finish_runtime_wfs_sensing!(wfs, atm, pupil, optic, src)
    return nothing
end

@inline function sense_core_after_advance!(atm::AbstractAtmosphere,
    pupil::PupilFunction, optic::AbstractControllableOptic, wfs::AbstractWFS,
    src::AbstractSource, det::AbstractDetector, rng::AbstractRNG,
    atmosphere_renderer)
    render_runtime_wfs_path!(atmosphere_renderer, wfs, atm, pupil, optic, src)
    measure_runtime_wfs!(wfs, atm, pupil, src, det, rng)
    finish_runtime_wfs_sensing!(wfs, atm, pupil, optic, src)
    return nothing
end

@inline function copy_pupil_path!(destination::PupilFunction,
    source::PupilFunction)
    require_same_plane_grid(destination.metadata, source.metadata;
        label="runtime pupil paths")
    destination.aperture_revision == source.aperture_revision || throw(
        InvalidConfiguration(
            "runtime pupil paths use different aperture revisions"))
    copyto!(destination.amplitude, source.amplitude)
    copyto!(destination.opd, source.opd)
    return destination
end

@inline function prepare_science_path!(::ReuseSensedOpticalPath,
    ::AbstractAtmosphere, sensed_pupil::PupilFunction,
    ::AbstractControllableOptic, ::AbstractSource,
    stage::PreparedRuntimeScienceStage)
    return copy_pupil_path!(stage.pupil, sensed_pupil)
end

@inline function prepare_science_path!(::RepropagateScienceOpticalPath,
    atm::AbstractAtmosphere, ::PupilFunction,
    optic::AbstractControllableOptic, src::AbstractSource,
    stage::PreparedRuntimeScienceStage)
    propagate_runtime_atmosphere!(atm, stage.pupil, src)
    update_surface!(optic)
    apply_surface!(stage.pupil, optic, DMAdditive())
    return stage.pupil
end

@inline prepare_science_path!(::ReuseSensedOpticalPath,
    atmosphere_renderer, atm::AbstractAtmosphere,
    sensed_pupil::PupilFunction, optic::AbstractControllableOptic,
    src::AbstractSource, stage::PreparedRuntimeScienceStage) =
    copy_pupil_path!(stage.pupil, sensed_pupil)

@inline function prepare_science_path!(::RepropagateScienceOpticalPath,
    atmosphere_renderer, atm::AbstractAtmosphere,
    ::PupilFunction, optic::AbstractControllableOptic,
    src::AbstractSource, stage::PreparedRuntimeScienceStage)
    render_prepared_atmosphere_path!(atmosphere_renderer, atm, stage.pupil,
        src)
    update_surface!(optic)
    apply_surface!(stage.pupil, optic, DMAdditive())
    return stage.pupil
end

@inline validate_runtime_science_aperture(::Nothing, ::Telescope) = nothing

@inline preflight_runtime_science!(::Nothing, detector) = nothing

@inline function preflight_runtime_science!(
    stage::PreparedRuntimeScienceStage{<:Any,<:Any,<:Any,
        <:DetectorAcquisitionPlan}, detector::Detector)
    _require_prepared_whole_acquisition(detector, stage.output,
        stage.acquisition)
    return nothing
end

@inline function validate_runtime_science_aperture(
    stage::PreparedRuntimeScienceStage, tel::Telescope)
    aperture_revision(tel) == stage.aperture_revision || throw(
        InvalidConfiguration(
            "telescope aperture changed after science-stage preparation; " *
            "reconstruct the runtime"))
    return nothing
end

@inline function capture_science_core!(plan::AbstractSciencePathPlan,
    atm::AbstractAtmosphere, tel::Telescope, sensed_pupil::PupilFunction,
    optic::AbstractControllableOptic,
    src::AbstractSource, stage::PreparedRuntimeScienceStage,
    det::Detector, rng::AbstractRNG)
    validate_runtime_science_aperture(stage, tel)
    prepare_science_path!(plan, atm, sensed_pupil, optic, src, stage)
    form_direct_image!(stage.imaging)
    capture!(det, stage.output, stage.acquisition, rng)
    return nothing
end

@inline function capture_science_core!(plan::AbstractSciencePathPlan,
    atmosphere_renderer, atm::AbstractAtmosphere, tel::Telescope,
    sensed_pupil::PupilFunction, optic::AbstractControllableOptic,
    src::AbstractSource,
    stage::PreparedRuntimeScienceStage, det::Detector, rng::AbstractRNG)
    validate_runtime_science_aperture(stage, tel)
    prepare_science_path!(plan, atmosphere_renderer, atm, sensed_pupil,
        optic, src, stage)
    form_direct_image!(stage.imaging)
    capture!(det, stage.output, stage.acquisition, rng)
    return nothing
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

`step_core!` is the hot-path composition of sensing, reconstruction, and DM
command application. Runtime science capture is a separately prepared stage.
The method overloads differ only in whether a WFS detector participates.
"""
@inline function step_core!(atm::AbstractAtmosphere, tel::Telescope,
    pupil::PupilFunction, optic::AbstractControllableOptic,
    wfs::AbstractWFS, src::AbstractSource, reconstructor, command::AbstractVector{T},
    atmosphere_step::Real, rng::AbstractRNG,
    control_sign::T) where {T<:AbstractFloat}
    sense_core!(atm, tel, pupil, optic, wfs, src, atmosphere_step, rng)
    reconstruct!(command, reconstructor, slopes(wfs))
    _stage_command!(optic, command, control_sign)
    return nothing
end

@inline function step_core!(atm::AbstractAtmosphere, tel::Telescope,
    pupil::PupilFunction, optic::AbstractControllableOptic,
    wfs::AbstractWFS, src::AbstractSource, wfs_detector::AbstractDetector, reconstructor,
    command::AbstractVector{T}, atmosphere_step::Real, rng::AbstractRNG,
    control_sign::T) where {T<:AbstractFloat}
    sense_core!(atm, tel, pupil, optic, wfs, src, wfs_detector,
        atmosphere_step, rng)
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
        sense_core!(runtime.atm, runtime.tel, runtime.wfs_pupil,
            runtime.optic, runtime.wfs, runtime.src,
            runtime.atmosphere_step, runtime.rng,
            runtime.wfs_atmosphere_renderer)
    else
        sense_core!(runtime.atm, runtime.tel, runtime.wfs_pupil,
            runtime.optic, runtime.wfs, runtime.src, runtime.wfs_detector,
            runtime.atmosphere_step, runtime.rng,
            runtime.wfs_atmosphere_renderer)
    end
    return nothing
end

@inline function sense_primary_core!(::ReuseAdvancedRuntimeAtmosphere,
    runtime::ClosedLoopRuntime)
    if isnothing(runtime.wfs_detector)
        sense_core_after_advance!(runtime.atm, runtime.wfs_pupil,
            runtime.optic, runtime.wfs, runtime.src, runtime.rng,
            runtime.wfs_atmosphere_renderer)
    else
        sense_core_after_advance!(runtime.atm, runtime.wfs_pupil,
            runtime.optic, runtime.wfs, runtime.src, runtime.wfs_detector,
            runtime.rng,
            runtime.wfs_atmosphere_renderer)
    end
    return nothing
end

@inline function sense_runtime!(advance_policy::AbstractRuntimeAtmosphereStep,
    runtime::ClosedLoopRuntime)
    validate_runtime_science_aperture(runtime.science_stage, runtime.tel)
    preflight_runtime_science!(runtime.science_stage,
        runtime.science_detector)
    sense_primary_core!(advance_policy, runtime)
    if requires_runtime_science_pixels(runtime)
        capture_science_core!(runtime.science_path,
            runtime.science_atmosphere_renderer, runtime.atm, runtime.tel,
            runtime.wfs_pupil, runtime.optic, runtime.science_src,
            runtime.science_stage,
            runtime.science_detector, runtime.rng)
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

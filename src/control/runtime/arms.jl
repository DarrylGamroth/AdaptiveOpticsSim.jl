"""
    OpticalWFSChannel(wfs; detector=nothing)

One wavefront-sensor consumer attached to a shared optical arm. The optional
detector participates in the WFS-owned detector-coupled measurement path.
"""
struct OpticalWFSChannel{W<:AbstractWFS,D}
    wfs::W
    detector::D
end

function OpticalWFSChannel(wfs::AbstractWFS; detector=nothing)
    require_same_backend(wfs, detector)
    return OpticalWFSChannel{typeof(wfs),typeof(detector)}(wfs, detector)
end

@inline backend(channel::OpticalWFSChannel) = backend(channel.wfs)

@inline _runtime_tuple(::Nothing) = ()
@inline _runtime_tuple(value::Tuple) = value
@inline _runtime_tuple(value) = (value,)

@inline _validate_wfs_channels(::Tuple{}) = nothing
@inline function _validate_wfs_channels(channels::Tuple{<:OpticalWFSChannel,Vararg})
    _validate_wfs_channels(Base.tail(channels))
    return nothing
end
_validate_wfs_channels(::Tuple) =
    throw(InvalidConfiguration("shared optical-arm WFS consumers must be OpticalWFSChannel objects"))

@inline _validate_science_detectors(::Tuple{}) = nothing
@inline function _validate_science_detectors(detectors::Tuple{<:Detector,Vararg})
    _validate_science_detectors(Base.tail(detectors))
    return nothing
end
function _validate_science_detectors(
    detectors::Tuple{<:AbstractDetector,Vararg})
    throw(UnsupportedAlgorithm(
        "shared-arm prepared science acquisition supports Detector frame " *
        "consumers; got $(typeof(first(detectors)))"))
end
_validate_science_detectors(::Tuple) =
    throw(InvalidConfiguration(
        "shared optical-arm science consumers must be Detector objects"))

"""
    SharedOpticalArm(label, source; wfs_channels=(), science_detectors=(),
                     science_zero_padding=1)

A source-specific optical branch with zero or more WFS consumers and zero or
more science detectors. The pupil path is rendered once for the arm, and one
prepared photon-arrival-rate image is shared by every science detector on
that arm.
"""
struct SharedOpticalArm{
    S<:AbstractSource,
    C<:Tuple,
    D<:Tuple,
    B<:AbstractArrayBackend,
}
    label::Symbol
    source::S
    wfs_channels::C
    science_detectors::D
    science_zero_padding::Int
end

function SharedOpticalArm(label::Symbol, source::AbstractSource;
    wfs_channels=(), science_detectors=(), science_zero_padding::Integer=1)
    channels = _runtime_tuple(wfs_channels)
    detectors = _runtime_tuple(science_detectors)
    _validate_wfs_channels(channels)
    _validate_science_detectors(detectors)
    isempty(channels) && isempty(detectors) &&
        throw(InvalidConfiguration("a shared optical arm requires at least one WFS or science consumer"))
    science_zero_padding >= 1 ||
        throw(InvalidConfiguration("shared optical-arm science_zero_padding must be >= 1"))
    selector = require_same_backend(channels..., detectors...)
    return SharedOpticalArm{
        typeof(source),
        typeof(channels),
        typeof(detectors),
        typeof(selector),
    }(
        label,
        source,
        channels,
        detectors,
        Int(science_zero_padding),
    )
end

@inline backend(::SharedOpticalArm{<:Any,<:Any,<:Any,B}) where {B} = B()
@inline wfs_channels(arm::SharedOpticalArm) = arm.wfs_channels
@inline science_detectors(arm::SharedOpticalArm) = arm.science_detectors
@inline science_frames(arm::SharedOpticalArm) = map(output_frame, arm.science_detectors)
@inline wfs_signals(arm::SharedOpticalArm) = map(channel -> slopes(channel.wfs), arm.wfs_channels)

"""
    SharedOpticalRuntime(runtime, arms...)

Execute a primary closed-loop runtime and source-specific auxiliary optical
arms against one atmosphere, telescope, controllable optic, command state, and
RNG. Each `sense!` advances the atmosphere exactly once. Each arm is rendered
into an explicit reusable pupil product; multiple consumers on that arm share
the path when their propagation models permit it.
"""
struct SharedOpticalRuntime{
    RT<:ClosedLoopRuntime,
    A<:Tuple,
    R<:Tuple,
    S<:Tuple,
    P<:PupilFunction,
    B<:AbstractArrayBackend,
} <: AbstractControlSimulation
    runtime::RT
    arms::A
    atmosphere_renderers::R
    science_stages::S
    pupil::P
end

function _shared_arm_with_source(arm::SharedOpticalArm,
    source::AbstractSource)
    return SharedOpticalArm(
        arm.label,
        source;
        wfs_channels=arm.wfs_channels,
        science_detectors=arm.science_detectors,
        science_zero_padding=arm.science_zero_padding,
    )
end

function _freeze_shared_arm_sources(runtime::ClosedLoopRuntime,
    arms::NTuple{N,SharedOpticalArm}) where {N}
    source_map = IdDict{Any,AbstractSource}()
    source_map[wfs_source(runtime.simulation)] = runtime.src
    source_map[science_source(runtime.simulation)] = runtime.science_src
    source_map[runtime.src] = runtime.src
    source_map[runtime.science_src] = runtime.science_src
    return ntuple(N) do i
        arm = arms[i]
        source = get!(source_map, arm.source) do
            freeze_source(arm.source)
        end
        _shared_arm_with_source(arm, source)
    end
end

function SharedOpticalRuntime(runtime::ClosedLoopRuntime,
    arms::Vararg{SharedOpticalArm,N}) where {N}
    N > 0 || throw(InvalidConfiguration("SharedOpticalRuntime requires at least one auxiliary optical arm"))
    labels = ntuple(i -> arms[i].label, N)
    allunique(labels) ||
        throw(InvalidConfiguration("shared optical-arm labels must be unique; got $labels"))
    frozen_arms = _freeze_shared_arm_sources(runtime, arms)
    selector = require_same_backend(runtime, frozen_arms...)
    renderers = map(arm -> prepare_runtime_atmosphere_path(
        runtime.atm, runtime.tel, arm.source), frozen_arms)
    science_stages = map(frozen_arms) do arm
        isempty(arm.science_detectors) && return nothing
        _prepare_runtime_science_stage(runtime.tel, arm.source,
            arm.science_detectors, arm.science_zero_padding)
    end
    pupil = PupilFunction(runtime.tel)
    return SharedOpticalRuntime{typeof(runtime),typeof(frozen_arms),
        typeof(renderers),typeof(science_stages),typeof(pupil),
        typeof(selector)}(
        runtime,
        frozen_arms,
        renderers,
        science_stages,
        pupil,
    )
end

@inline backend(::SharedOpticalRuntime{
    <:Any,<:Any,<:Any,<:Any,<:Any,B}) where {B} = B()
@inline primary_runtime(runtime::SharedOpticalRuntime) = runtime.runtime
@inline optical_arms(runtime::SharedOpticalRuntime) = runtime.arms
@inline wfs_source(runtime::SharedOpticalRuntime) = wfs_source(runtime.runtime)
@inline science_source(runtime::SharedOpticalRuntime) = science_source(runtime.runtime)
@inline runtime_profile(runtime::SharedOpticalRuntime) = runtime_profile(runtime.runtime)
@inline runtime_execution_plan(runtime::SharedOpticalRuntime) = runtime_execution_plan(runtime.runtime)
@inline runtime_latency(runtime::SharedOpticalRuntime) = runtime_latency(runtime.runtime)
@inline runtime_atmosphere_step(runtime::SharedOpticalRuntime) =
    runtime_atmosphere_step(runtime.runtime)
@inline command_layout(runtime::SharedOpticalRuntime) = command_layout(runtime.runtime)
@inline simulation_readout(runtime::SharedOpticalRuntime) = simulation_readout(runtime.runtime)

@inline function prepare_shared_arm_path!(
    runtime::ClosedLoopRuntime, arm::SharedOpticalArm,
    atmosphere_renderer, pupil::PupilFunction)
    render_prepared_atmosphere_path!(atmosphere_renderer, runtime.atm,
        pupil, arm.source)
    update_surface!(runtime.optic)
    apply_surface!(pupil, runtime.optic, DMAdditive())
    return nothing
end

@inline function sense_shared_wfs_channel!(
    channel::OpticalWFSChannel{<:AbstractWFS,Nothing},
    runtime::ClosedLoopRuntime, pupil::PupilFunction,
    source::AbstractSource)
    prepare_shared_runtime_wfs!(channel.wfs, runtime.atm, pupil,
        runtime.optic, source)
    measure_runtime_wfs!(channel.wfs, runtime.atm, pupil, source,
        runtime.rng)
    finish_runtime_wfs_sensing!(channel.wfs, runtime.atm, pupil,
        runtime.optic, source)
    return nothing
end

@inline function sense_shared_wfs_channel!(
    channel::OpticalWFSChannel{<:AbstractWFS,<:AbstractDetector},
    runtime::ClosedLoopRuntime, pupil::PupilFunction,
    source::AbstractSource)
    prepare_shared_runtime_wfs!(channel.wfs, runtime.atm, pupil,
        runtime.optic, source)
    measure_runtime_wfs!(channel.wfs, runtime.atm, pupil, source,
        channel.detector, runtime.rng)
    finish_runtime_wfs_sensing!(channel.wfs, runtime.atm, pupil,
        runtime.optic, source)
    return nothing
end

@inline sense_shared_wfs_channels!(::Tuple{}, runtime::ClosedLoopRuntime,
    pupil::PupilFunction, source::AbstractSource) = nothing

@inline function sense_shared_wfs_channels!(channels::Tuple,
    runtime::ClosedLoopRuntime, pupil::PupilFunction,
    source::AbstractSource)
    sense_shared_wfs_channel!(first(channels), runtime, pupil, source)
    sense_shared_wfs_channels!(Base.tail(channels), runtime, pupil, source)
    return nothing
end

@inline capture_shared_science_detectors!(::Tuple{}, ::Tuple{},
    stage::PreparedRuntimeScienceStage, rng::AbstractRNG) = nothing

@inline function capture_shared_science_detectors!(detectors::Tuple,
    acquisitions::Tuple, stage::PreparedRuntimeScienceStage,
    rng::AbstractRNG)
    capture!(first(detectors), stage.output, first(acquisitions), rng)
    capture_shared_science_detectors!(Base.tail(detectors),
        Base.tail(acquisitions), stage, rng)
    return nothing
end

@inline shared_arm_science_path(::Tuple{}) = ReuseSensedOpticalPath()
@inline shared_arm_science_path(channels::Tuple) =
    shared_arm_science_path(first(channels).wfs, Base.tail(channels))
@inline shared_arm_science_path(::AbstractWFS, channels::Tuple) =
    shared_arm_science_path(channels)
@inline shared_arm_science_path(::CurvatureWFS, ::Tuple) =
    RepropagateScienceOpticalPath()

@inline capture_shared_science_arm!(::Nothing, runtime::ClosedLoopRuntime,
    arm::SharedOpticalArm, atmosphere_renderer,
    pupil::PupilFunction) = nothing

@inline function capture_shared_science_arm!(
    stage::PreparedRuntimeScienceStage,
    runtime::ClosedLoopRuntime, arm::SharedOpticalArm,
    atmosphere_renderer, pupil::PupilFunction)
    prepare_science_path!(shared_arm_science_path(arm.wfs_channels),
        atmosphere_renderer, runtime.atm, pupil, runtime.optic, arm.source,
        stage)
    form_direct_image!(stage.imaging)
    capture_shared_science_detectors!(arm.science_detectors,
        stage.acquisition, stage, runtime.rng)
    return nothing
end

@inline function execute_shared_optical_arm!(runtime::ClosedLoopRuntime,
    arm::SharedOpticalArm, atmosphere_renderer, science_stage,
    pupil::PupilFunction)
    prepare_shared_arm_path!(runtime, arm, atmosphere_renderer, pupil)
    sense_shared_wfs_channels!(arm.wfs_channels, runtime, pupil, arm.source)
    capture_shared_science_arm!(science_stage, runtime, arm,
        atmosphere_renderer, pupil)
    return nothing
end

@inline execute_shared_optical_arms!(::Tuple{}, ::Tuple{}, ::Tuple{},
    runtime::ClosedLoopRuntime, pupil::PupilFunction) = nothing

@inline function execute_shared_optical_arms!(
    arms::Tuple{<:SharedOpticalArm,Vararg}, renderers::Tuple,
    science_stages::Tuple,
    runtime::ClosedLoopRuntime, pupil::PupilFunction)
    execute_shared_optical_arm!(runtime, first(arms), first(renderers),
        first(science_stages), pupil)
    return execute_shared_optical_arms!(Base.tail(arms),
        Base.tail(renderers), Base.tail(science_stages), runtime, pupil)
end

@inline preflight_shared_science_detectors!(::Tuple{}, ::Tuple{},
    ::PreparedRuntimeScienceStage) = nothing

@inline function preflight_shared_science_detectors!(
    detectors::Tuple{<:Detector,Vararg},
    acquisitions::Tuple{<:DetectorAcquisitionPlan,Vararg},
    stage::PreparedRuntimeScienceStage)
    _require_prepared_whole_acquisition(first(detectors), stage.output,
        first(acquisitions))
    preflight_shared_science_detectors!(Base.tail(detectors),
        Base.tail(acquisitions), stage)
    return nothing
end

@inline preflight_shared_science_arm!(::Nothing,
    ::SharedOpticalArm) = nothing

@inline function preflight_shared_science_arm!(
    stage::PreparedRuntimeScienceStage, arm::SharedOpticalArm)
    preflight_shared_science_detectors!(arm.science_detectors,
        stage.acquisition, stage)
    return nothing
end

@inline preflight_shared_science_arms!(::Tuple{}, ::Tuple{}) = nothing

@inline function preflight_shared_science_arms!(
    stages::Tuple, arms::Tuple{<:SharedOpticalArm,Vararg})
    preflight_shared_science_arm!(first(stages), first(arms))
    preflight_shared_science_arms!(Base.tail(stages), Base.tail(arms))
    return nothing
end

@inline validate_shared_science_apertures!(::Tuple{}, ::Telescope) = nothing

@inline function validate_shared_science_apertures!(stages::Tuple,
    tel::Telescope)
    validate_runtime_science_aperture(first(stages), tel)
    validate_shared_science_apertures!(Base.tail(stages), tel)
    return nothing
end

function prepare!(runtime::SharedOpticalRuntime)
    prepare!(runtime.runtime)
    @inbounds for arm in runtime.arms
        @inbounds for channel in arm.wfs_channels
            prepare_runtime_wfs!(channel.wfs, runtime.pupil, arm.source)
        end
    end
    return runtime
end

function sense!(runtime::SharedOpticalRuntime)
    primary = runtime.runtime
    validate_runtime_science_aperture(primary.science_stage, primary.tel)
    validate_shared_science_apertures!(runtime.science_stages, primary.tel)
    preflight_runtime_science!(primary.science_stage,
        primary.science_detector)
    preflight_shared_science_arms!(runtime.science_stages, runtime.arms)
    advance_runtime_atmosphere!(primary.wfs, primary.atm, primary.tel,
        primary.atmosphere_step, primary.rng)
    sense_runtime!(ReuseAdvancedRuntimeAtmosphere(), primary)
    execute_shared_optical_arms!(runtime.arms, runtime.atmosphere_renderers,
        runtime.science_stages, primary, runtime.pupil)
    return runtime
end

function step!(runtime::SharedOpticalRuntime)
    sense!(runtime)
    reconstruct!(runtime.runtime)
    apply_runtime_command!(runtime.runtime)
    return runtime
end

@inline function set_command!(runtime::SharedOpticalRuntime, command)
    return set_command!(runtime.runtime, command)
end

function synchronize_runtime!(runtime::SharedOpticalRuntime)
    synchronize_runtime!(runtime.runtime)
    return runtime
end

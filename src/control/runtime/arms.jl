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
@inline function _validate_science_detectors(detectors::Tuple{<:AbstractDetector,Vararg})
    _validate_science_detectors(Base.tail(detectors))
    return nothing
end
_validate_science_detectors(::Tuple) =
    throw(InvalidConfiguration("shared optical-arm science consumers must be detector objects"))

"""
    SharedOpticalArm(label, source; wfs_channels=(), science_detectors=(),
                     science_zero_padding=1)

A source-specific optical branch with zero or more WFS consumers and zero or
more science detectors. The pupil path is rendered once for the arm, and one
computed PSF is shared by every science detector on that arm.
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
RNG. Each `sense!` advances the atmosphere exactly once. Consecutive arms that
reference the same source object reuse the already-rendered pupil path.
"""
struct SharedOpticalRuntime{
    RT<:ClosedLoopRuntime,
    A<:Tuple,
    B<:AbstractArrayBackend,
} <: AbstractControlSimulation
    runtime::RT
    arms::A
end

function SharedOpticalRuntime(runtime::ClosedLoopRuntime,
    arms::Vararg{SharedOpticalArm,N}) where {N}
    N > 0 || throw(InvalidConfiguration("SharedOpticalRuntime requires at least one auxiliary optical arm"))
    labels = ntuple(i -> arms[i].label, N)
    allunique(labels) ||
        throw(InvalidConfiguration("shared optical-arm labels must be unique; got $labels"))
    selector = require_same_backend(runtime, arms...)
    return SharedOpticalRuntime{typeof(runtime),typeof(arms),typeof(selector)}(
        runtime,
        arms,
    )
end

@inline backend(::SharedOpticalRuntime{<:Any,<:Any,B}) where {B} = B()
@inline primary_runtime(runtime::SharedOpticalRuntime) = runtime.runtime
@inline optical_arms(runtime::SharedOpticalRuntime) = runtime.arms
@inline wfs_source(runtime::SharedOpticalRuntime) = wfs_source(runtime.runtime)
@inline science_source(runtime::SharedOpticalRuntime) = science_source(runtime.runtime)
@inline runtime_profile(runtime::SharedOpticalRuntime) = runtime_profile(runtime.runtime)
@inline runtime_execution_plan(runtime::SharedOpticalRuntime) = runtime_execution_plan(runtime.runtime)
@inline runtime_latency(runtime::SharedOpticalRuntime) = runtime_latency(runtime.runtime)
@inline command_layout(runtime::SharedOpticalRuntime) = command_layout(runtime.runtime)
@inline simulation_readout(runtime::SharedOpticalRuntime) = simulation_readout(runtime.runtime)

abstract type AbstractSharedArmPathPlan end
struct ReuseSharedArmPath <: AbstractSharedArmPathPlan end
struct RenderSharedArmPath <: AbstractSharedArmPathPlan end

@inline shared_arm_path_plan(current_source::AbstractSource,
    target_source::AbstractSource) = current_source === target_source ?
    ReuseSharedArmPath() : RenderSharedArmPath()

@inline prepare_shared_arm_path!(::ReuseSharedArmPath,
    runtime::ClosedLoopRuntime, arm::SharedOpticalArm) = nothing

@inline function prepare_shared_arm_path!(::RenderSharedArmPath,
    runtime::ClosedLoopRuntime, arm::SharedOpticalArm)
    propagate_runtime_atmosphere!(runtime.atm, runtime.tel, arm.source)
    apply!(runtime.optic, runtime.tel, DMAdditive())
    return nothing
end

@inline function sense_shared_wfs_channel!(
    channel::OpticalWFSChannel{<:AbstractWFS,Nothing},
    runtime::ClosedLoopRuntime, source::AbstractSource)
    prepare_shared_runtime_wfs!(channel.wfs, runtime.atm, runtime.tel,
        runtime.optic, source)
    measure_runtime_wfs!(channel.wfs, runtime.atm, runtime.tel, source,
        runtime.rng)
    finish_runtime_wfs_sensing!(channel.wfs, runtime.atm, runtime.tel,
        runtime.optic, source)
    return nothing
end

@inline function sense_shared_wfs_channel!(
    channel::OpticalWFSChannel{<:AbstractWFS,<:AbstractDetector},
    runtime::ClosedLoopRuntime, source::AbstractSource)
    prepare_shared_runtime_wfs!(channel.wfs, runtime.atm, runtime.tel,
        runtime.optic, source)
    measure_runtime_wfs!(channel.wfs, runtime.atm, runtime.tel, source,
        channel.detector, runtime.rng)
    finish_runtime_wfs_sensing!(channel.wfs, runtime.atm, runtime.tel,
        runtime.optic, source)
    return nothing
end

@inline sense_shared_wfs_channels!(::Tuple{}, runtime::ClosedLoopRuntime,
    source::AbstractSource) = nothing

@inline function sense_shared_wfs_channels!(channels::Tuple,
    runtime::ClosedLoopRuntime, source::AbstractSource)
    sense_shared_wfs_channel!(first(channels), runtime, source)
    sense_shared_wfs_channels!(Base.tail(channels), runtime, source)
    return nothing
end

@inline capture_shared_science_detectors!(::Tuple{}, psf,
    source::AbstractSource, rng::AbstractRNG) = nothing

@inline function capture_shared_science_detectors!(detectors::Tuple, psf,
    source::AbstractSource, rng::AbstractRNG)
    capture!(first(detectors), psf, source, rng)
    capture_shared_science_detectors!(Base.tail(detectors), psf, source, rng)
    return nothing
end

@inline capture_shared_science_arm!(::Tuple{}, runtime::ClosedLoopRuntime,
    arm::SharedOpticalArm) = nothing

@inline function capture_shared_science_arm!(detectors::Tuple,
    runtime::ClosedLoopRuntime, arm::SharedOpticalArm)
    psf = compute_psf!(runtime.tel, arm.source;
        zero_padding=arm.science_zero_padding)
    capture_shared_science_detectors!(detectors, psf, arm.source, runtime.rng)
    return nothing
end

@inline function execute_shared_optical_arm!(runtime::ClosedLoopRuntime,
    arm::SharedOpticalArm, current_source::AbstractSource)
    prepare_shared_arm_path!(shared_arm_path_plan(current_source, arm.source),
        runtime, arm)
    sense_shared_wfs_channels!(arm.wfs_channels, runtime, arm.source)
    capture_shared_science_arm!(arm.science_detectors, runtime, arm)
    return arm.source
end

@inline execute_shared_optical_arms!(::Tuple{}, runtime::ClosedLoopRuntime,
    current_source::AbstractSource) = current_source

@inline function execute_shared_optical_arms!(arms::Tuple{<:SharedOpticalArm,Vararg},
    runtime::ClosedLoopRuntime, current_source::AbstractSource)
    next_source = execute_shared_optical_arm!(runtime, first(arms), current_source)
    return execute_shared_optical_arms!(Base.tail(arms), runtime, next_source)
end

@inline primary_terminal_source(runtime::ClosedLoopRuntime) =
    requires_runtime_science_pixels(runtime) ? runtime.science_src : runtime.src

function prepare!(runtime::SharedOpticalRuntime)
    prepare!(runtime.runtime)
    @inbounds for arm in runtime.arms
        @inbounds for channel in arm.wfs_channels
            prepare_runtime_wfs!(channel.wfs, runtime.runtime.tel, arm.source)
        end
    end
    return runtime
end

function sense!(runtime::SharedOpticalRuntime)
    primary = runtime.runtime
    advance_runtime_atmosphere!(primary.wfs, primary.atm, primary.tel,
        primary.rng)
    sense_runtime!(ReuseAdvancedRuntimeAtmosphere(), primary)
    execute_shared_optical_arms!(runtime.arms, primary,
        primary_terminal_source(primary))
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

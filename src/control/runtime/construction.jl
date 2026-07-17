@inline _require_runtime_science_detector(::Detector) = nothing

function _require_runtime_science_detector(detector::AbstractDetector)
    throw(UnsupportedAlgorithm(
        "prepared runtime direct-science acquisition supports Detector " *
        "frame consumers; got $(typeof(detector))"))
end

function _require_runtime_science_detector(detector)
    throw(InvalidConfiguration(
        "runtime science pixel output requires a Detector; got " *
        "$(typeof(detector))"))
end

@inline function _require_runtime_science_configuration(
    outputs::RuntimeOutputRequirements, science_detector)
    outputs.science_pixels || return nothing
    isnothing(science_detector) && throw(InvalidConfiguration(
        "runtime science pixel output requires a Detector; got nothing"))
    _require_runtime_science_detector(science_detector)
    return nothing
end

@inline function _prepare_runtime_acquisition(detector::Detector,
    output::IntensityMap)
    return prepare_detector_acquisition(detector, output)
end

@inline _prepare_runtime_acquisitions(::Tuple{}, output::IntensityMap) = ()

@inline function _prepare_runtime_acquisitions(
    detectors::Tuple{<:Detector,Vararg}, output::IntensityMap)
    return (
        prepare_detector_acquisition(first(detectors), output),
        _prepare_runtime_acquisitions(Base.tail(detectors), output)...,
    )
end

function _prepared_runtime_science_stage(pupil::PupilFunction, imaging,
    output::IntensityMap, detector::Detector, revision::UInt)
    acquisition = _prepare_runtime_acquisition(detector, output)
    return PreparedRuntimeScienceStage(pupil, imaging, output, acquisition,
        revision)
end

function _prepared_runtime_science_stage(pupil::PupilFunction, imaging,
    output::IntensityMap, detectors::Tuple, revision::UInt)
    acquisition = _prepare_runtime_acquisitions(detectors, output)
    return PreparedRuntimeScienceStage(pupil, imaging, output, acquisition,
        revision)
end

function _prepared_runtime_science_stage(::PupilFunction, imaging,
    ::OpticalProductBundle, detector_or_detectors, ::UInt)
    throw(UnsupportedAlgorithm(
        "runtime detector acquisition requires one integrated IntensityMap; " *
        "bundled spectral direct-imaging output must be integrated onto an " *
        "explicit common grid before runtime construction"))
end

function _prepare_runtime_science_stage(tel::Telescope,
    source::AbstractSource, detector_or_detectors, zero_padding::Int)
    zero_padding >= 1 || throw(InvalidConfiguration(
        "runtime science_zero_padding must be >= 1 when science pixels are requested"))
    pupil = PupilFunction(tel)
    imaging = prepare_direct_imaging(tel, pupil, source;
        zero_padding=zero_padding)
    return _prepared_runtime_science_stage(pupil, imaging,
        direct_imaging_output(imaging), detector_or_detectors,
        aperture_revision(tel))
end

@inline function _prepare_primary_runtime_science_stage(
    outputs::RuntimeOutputRequirements, tel::Telescope,
    source::AbstractSource, detector, zero_padding::Int)
    outputs.science_pixels && !isnothing(detector) || return nothing
    return _prepare_runtime_science_stage(tel, source, detector,
        zero_padding)
end

function ClosedLoopRuntime(simulation::AOSimulation, reconstructor;
    atmosphere_step::Real,
    wfs_detector=nothing, science_detector=nothing, rng=runtime_rng(),
    profile::AbstractRuntimeProfile=default_runtime_profile(),
    execution_plan::Union{AbstractRuntimeExecutionPlan,Nothing}=nothing,
    outputs::RuntimeOutputRequirements=default_runtime_outputs(wfs_detector=wfs_detector, science_detector=science_detector),
    latency::RuntimeLatencyModel=default_runtime_latency(profile),
    control_sign::Real=-1.0, science_zero_padding::Union{Int,Nothing}=nothing,
    command_layout::Union{RuntimeCommandLayout,Nothing}=nothing)
    _require_runtime_science_configuration(outputs, science_detector)
    selector = require_same_backend(simulation, wfs_detector, science_detector)
    resolved_execution_plan = something(execution_plan, default_runtime_execution_plan(selector))
    validate_runtime_execution_plan!(resolved_execution_plan, selector, reconstructor)
    outputs.slopes || throw(InvalidConfiguration("ClosedLoopRuntime requires slope outputs for reconstruction"))
    T = eltype(command_storage(simulation.optic))
    resolved_atmosphere_step = _validated_atmosphere_step(T, atmosphere_step)
    command = similar(command_storage(simulation.optic))
    fill!(command, zero(T))
    reconstruct_buffer = similar(command)
    fill!(reconstruct_buffer, zero(T))
    slopes = similar(simulation.wfs.state.slopes)
    fill!(slopes, zero(eltype(slopes)))
    measurement_delay = VectorDelayLine(simulation.wfs.state.slopes, latency.measurement_delay_frames)
    readout_delay = VectorDelayLine(simulation.wfs.state.slopes, latency.readout_delay_frames)
    reconstruction_delay = VectorDelayLine(command, latency.reconstruction_delay_frames)
    dm_delay = VectorDelayLine(command, latency.dm_delay_frames)
    resolved_zero_padding = something(science_zero_padding, runtime_science_zero_padding(profile))
    resolved_command_layout = something(command_layout, AdaptiveOpticsSim.command_layout(simulation.optic))
    resolved_command_layout.total_length == length(command) ||
        throw(InvalidConfiguration("runtime command layout total length must match the runtime command vector length"))
    runtime_source = freeze_source(wfs_source(simulation))
    runtime_science_source = science_source(simulation) === wfs_source(simulation) ?
        runtime_source : freeze_source(science_source(simulation))
    science_stage = _prepare_primary_runtime_science_stage(outputs,
        simulation.tel, runtime_science_source, science_detector,
        resolved_zero_padding)
    resolved_science_path = science_path_plan(runtime_source,
        runtime_science_source)
    wfs_atmosphere_renderer = prepare_runtime_atmosphere_path(
        simulation.atm, simulation.tel, runtime_source)
    science_atmosphere_renderer = runtime_science_source === runtime_source ?
        wfs_atmosphere_renderer : prepare_runtime_atmosphere_path(
            simulation.atm, simulation.tel, runtime_science_source)
    runtime = ClosedLoopRuntime{
        typeof(simulation),
        typeof(simulation.tel),
        typeof(simulation.atm),
        typeof(wfs_atmosphere_renderer),
        typeof(science_atmosphere_renderer),
        typeof(runtime_source),
        typeof(runtime_science_source),
        typeof(simulation.optic),
        typeof(simulation.wfs),
        typeof(reconstructor),
        typeof(command),
        typeof(slopes),
        typeof(wfs_detector),
        typeof(science_detector),
        typeof(science_stage),
        typeof(rng),
        typeof(profile),
        typeof(resolved_execution_plan),
        typeof(outputs),
        typeof(measurement_delay),
        typeof(readout_delay),
        typeof(reconstruction_delay),
        typeof(dm_delay),
        typeof(resolved_command_layout),
        typeof(resolved_science_path),
        T,
        typeof(selector),
    }(
        simulation,
        simulation.tel,
        simulation.atm,
        wfs_atmosphere_renderer,
        science_atmosphere_renderer,
        runtime_source,
        runtime_science_source,
        simulation.optic,
        simulation.wfs,
        reconstructor,
        command,
        reconstruct_buffer,
        slopes,
        wfs_detector,
        science_detector,
        science_stage,
        rng,
        profile,
        resolved_execution_plan,
        outputs,
        latency,
        measurement_delay,
        readout_delay,
        reconstruction_delay,
        dm_delay,
        resolved_command_layout,
        resolved_science_path,
        resolved_atmosphere_step,
        T(control_sign),
        resolved_zero_padding,
        false,
    )
    set_runtime_wfs_output_policy!(runtime.wfs, runtime.outputs)
    return runtime
end

@inline wfs_output_frame(wfs::AbstractWFS, det::AbstractDetector) = output_frame(det)
@inline wfs_output_frame_prototype(wfs::AbstractWFS, det::AbstractDetector) = wfs_output_frame(wfs, det)
@inline wfs_output_metadata(::AbstractWFS) = nothing

simulation_interface(runtime::ClosedLoopRuntime) = SimulationInterface(runtime)
@inline backend(interface::SimulationInterface) = backend(interface.runtime)
@inline backend(interface::CompositeSimulationInterface) = backend(interface.interfaces[1])
@inline wfs_source(interface::SimulationInterface) = wfs_source(interface.runtime)
@inline science_source(interface::SimulationInterface) = science_source(interface.runtime)
simulation_interface(interface::SimulationInterface) = interface
simulation_interface(interface::CompositeSimulationInterface) = interface

"""
    prepare_runtime_wfs!(wfs, tel, src)

Run any WFS-specific runtime precomputation needed before repeated `step!`
calls. The default implementation is a no-op and returns `wfs`.
"""
@inline function prepare_runtime_wfs!(wfs::AbstractWFS, tel::Telescope, src)
    return wfs
end

@inline set_runtime_wfs_output_policy!(wfs::AbstractWFS, outputs) = wfs

supports_prepared_runtime(runtime::ClosedLoopRuntime) = supports_prepared_runtime(runtime.wfs, runtime.src)
supports_detector_output(runtime::ClosedLoopRuntime) = requires_runtime_wfs_pixels(runtime) || requires_runtime_science_pixels(runtime)
supports_stacked_sources(runtime::ClosedLoopRuntime) = supports_stacked_sources(runtime.wfs, runtime.src)
supports_grouped_execution(::CompositeSimulationInterface) = true

function prepare!(runtime::ClosedLoopRuntime)
    prepare_runtime_wfs!(runtime.wfs, runtime.tel, runtime.src)
    set_runtime_wfs_output_policy!(runtime.wfs, runtime.outputs)
    runtime.prepared = true
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
    slopes = similar(runtime.slopes)
    copyto!(slopes, runtime.slopes)
    wfs_frame = requires_runtime_wfs_pixels(runtime) ? similar(wfs_output_frame_prototype(runtime.wfs, runtime.wfs_detector)) : nothing
    science_frame = requires_runtime_science_pixels(runtime) ? similar(output_frame(runtime.science_detector)) : nothing
    interface = SimulationInterface{typeof(runtime), typeof(command), typeof(slopes), typeof(wfs_frame), typeof(science_frame)}(
        runtime,
        command,
        slopes,
        wfs_frame,
        science_frame,
    )
    if !isnothing(interface.wfs_frame)
        source = wfs_output_frame(runtime.wfs, runtime.wfs_detector)
        initial = size(source) == size(interface.wfs_frame) ? source : wfs_output_frame_prototype(runtime.wfs, runtime.wfs_detector)
        copyto!(interface.wfs_frame, initial)
    end
    if !isnothing(interface.science_frame)
        copyto!(interface.science_frame, output_frame(runtime.science_detector))
    end
    return interface
end

function CompositeSimulationInterface(interfaces::SimulationInterface...;
    outputs::GroupedRuntimeOutputRequirements=default_grouped_runtime_outputs(interfaces...))
    length(interfaces) > 0 || throw(InvalidConfiguration("CompositeSimulationInterface requires at least one interface"))
    first_interface = interfaces[1]
    total_command = mapreduce(i -> length(simulation_command(i)), +, interfaces)
    total_slopes = mapreduce(i -> length(simulation_slopes(i)), +, interfaces)
    command = similar(simulation_command(first_interface), total_command)
    slopes = similar(simulation_slopes(first_interface), total_slopes)
    plan = grouped_runtime_output_plan(outputs, interfaces)
    wfs_frames = ntuple(i -> begin
        frame = simulation_wfs_frame(interfaces[i])
        !plan.wfs_frames || isnothing(frame) ? nothing : similar(frame)
    end, length(interfaces))
    science_frames = ntuple(i -> begin
        frame = simulation_science_frame(interfaces[i])
        !plan.science_frames || isnothing(frame) ? nothing : similar(frame)
    end, length(interfaces))
    wfs_stack = plan.wfs_stack ? grouped_frame_stack_buffer(map(simulation_wfs_frame, interfaces)) : nothing
    science_stack = plan.science_stack ? grouped_frame_stack_buffer(map(simulation_science_frame, interfaces)) : nothing
    multi = CompositeSimulationInterface{
        typeof(interfaces),
        typeof(command),
        typeof(slopes),
        typeof(wfs_frames),
        typeof(science_frames),
        typeof(outputs),
        typeof(wfs_stack),
        typeof(science_stack),
    }(
        interfaces,
        command,
        slopes,
        wfs_frames,
        science_frames,
        outputs,
        wfs_stack,
        science_stack,
    )
    return snapshot_outputs!(multi)
end

function CompositeSimulationInterface(; outputs::GroupedRuntimeOutputRequirements=default_grouped_runtime_outputs())
    throw(InvalidConfiguration("CompositeSimulationInterface requires at least one interface"))
end

function CompositeSimulationInterface(runtimes::ClosedLoopRuntime...;
    outputs::GroupedRuntimeOutputRequirements=default_grouped_runtime_outputs(map(SimulationInterface, runtimes)...))
    return CompositeSimulationInterface(map(SimulationInterface, runtimes)...; outputs=outputs)
end

function with_reconstructor(runtime::ClosedLoopRuntime, reconstructor)
    refreshed = ClosedLoopRuntime(
        runtime.simulation,
        reconstructor;
        atmosphere_step=runtime.atmosphere_step,
        wfs_detector=runtime.wfs_detector,
        science_detector=runtime.science_detector,
        rng=runtime.rng,
        profile=runtime.profile,
        execution_plan=runtime.execution_plan,
        outputs=runtime.outputs,
        latency=runtime.latency,
        control_sign=runtime.control_sign,
        science_zero_padding=runtime.science_zero_padding,
        command_layout=runtime.command_layout,
    )
    copyto!(refreshed.command, runtime.command)
    copyto!(refreshed.reconstruct_buffer, runtime.reconstruct_buffer)
    copyto!(refreshed.slopes, runtime.slopes)
    set_command!(refreshed.optic, command_storage(runtime.optic))
    refreshed.prepared = runtime.prepared
    return refreshed
end

with_reconstructor(interface::SimulationInterface, reconstructor) = SimulationInterface(with_reconstructor(interface.runtime, reconstructor))

function with_reconstructors(interface::CompositeSimulationInterface, reconstructors::Vararg{Any,N}) where {N}
    length(interface.interfaces) == N ||
        throw(DimensionMismatchError("reconstructor count must match interface count"))
    refreshed = ntuple(i -> with_reconstructor(interface.interfaces[i], reconstructors[i]), N)
    return CompositeSimulationInterface(refreshed...)
end

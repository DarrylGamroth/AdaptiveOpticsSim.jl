function ClosedLoopRuntime(simulation::AOSimulation, reconstructor;
    wfs_detector=nothing, science_detector=nothing, rng=MersenneTwister(0),
    profile::AbstractRuntimeProfile=default_runtime_profile(),
    products::RuntimeProductRequirements=default_runtime_products(wfs_detector=wfs_detector, science_detector=science_detector),
    latency::RuntimeLatencyModel=default_runtime_latency(profile),
    control_sign::Real=-1.0, science_zero_padding::Union{Int,Nothing}=nothing,
    command_layout::Union{RuntimeCommandLayout,Nothing}=nothing)
    selector = require_same_backend(simulation, wfs_detector, science_detector)
    products.slopes || throw(InvalidConfiguration("ClosedLoopRuntime requires slope products for reconstruction"))
    T = eltype(command_storage(simulation.optic))
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
    runtime = ClosedLoopRuntime{
        typeof(simulation),
        typeof(simulation.tel),
        typeof(simulation.atm),
        typeof(simulation.src),
        typeof(simulation.optic),
        typeof(simulation.wfs),
        typeof(reconstructor),
        typeof(command),
        typeof(slopes),
        typeof(wfs_detector),
        typeof(science_detector),
        typeof(rng),
        typeof(profile),
        typeof(products),
        typeof(measurement_delay),
        typeof(readout_delay),
        typeof(reconstruction_delay),
        typeof(dm_delay),
        typeof(resolved_command_layout),
        T,
        typeof(selector),
    }(
        simulation,
        simulation.tel,
        simulation.atm,
        simulation.src,
        simulation.optic,
        simulation.wfs,
        reconstructor,
        command,
        reconstruct_buffer,
        slopes,
        wfs_detector,
        science_detector,
        rng,
        profile,
        products,
        latency,
        measurement_delay,
        readout_delay,
        reconstruction_delay,
        dm_delay,
        resolved_command_layout,
        T(control_sign),
        resolved_zero_padding,
        false,
    )
    set_runtime_wfs_product_policy!(runtime.wfs, runtime.products)
    return runtime
end

@inline wfs_output_frame(wfs::AbstractWFS, det::AbstractDetector) = output_frame(det)
@inline wfs_output_frame(wfs::ShackHartmann{<:Diffractive}, ::Nothing) = sh_exported_spot_cube(wfs)
@inline wfs_output_frame(wfs::ZernikeWFS, ::Nothing) = camera_frame(wfs)
@inline wfs_output_frame(wfs::CurvatureWFS, ::Nothing) = camera_frame(wfs)
@inline wfs_output_frame_prototype(wfs::PyramidWFS, ::Nothing) = camera_frame(wfs)
@inline wfs_output_frame_prototype(wfs::BioEdgeWFS, ::Nothing) = camera_frame(wfs)
@inline wfs_output_frame_prototype(wfs::ShackHartmann{<:Diffractive}, ::Nothing) = sh_exported_spot_cube(wfs)
@inline wfs_output_frame_prototype(wfs::ZernikeWFS, ::Nothing) = camera_frame(wfs)
@inline wfs_output_frame_prototype(wfs::CurvatureWFS, ::Nothing) = camera_frame(wfs)
@inline wfs_output_frame(wfs::ShackHartmann{<:Diffractive}, det::AbstractDetector) = sh_exported_spot_cube(wfs)
@inline wfs_output_frame_prototype(wfs::AbstractWFS, det::AbstractDetector) = wfs_output_frame(wfs, det)
@inline wfs_output_frame_prototype(wfs::PyramidWFS, det::AbstractDetector) = camera_frame(wfs)
@inline wfs_output_frame_prototype(wfs::BioEdgeWFS, det::AbstractDetector) = camera_frame(wfs)
@inline wfs_output_frame_prototype(wfs::ShackHartmann{<:Diffractive}, det::AbstractDetector) = sh_exported_spot_cube(wfs)
@inline wfs_output_metadata(wfs::ShackHartmann) = (
    n_subap=subaperture_layout(wfs).n_subap,
    n_valid_subap=n_valid_subapertures(subaperture_layout(wfs)),
    subap_pixels=subaperture_layout(wfs).subap_pixels,
    pitch_m=subaperture_layout(wfs).pitch_m,
    slopes_units=subaperture_calibration(wfs).slopes_units,
    calibrated=subaperture_calibration(wfs).calibrated,
)
@inline wfs_output_frame(wfs::ZernikeWFS, det::AbstractDetector) = camera_frame(wfs)
@inline wfs_output_frame(wfs::CurvatureWFS, det::AbstractDetector) = camera_frame(wfs)
@inline wfs_output_frame_prototype(wfs::ZernikeWFS, det::AbstractDetector) = camera_frame(wfs)
@inline wfs_output_frame_prototype(wfs::CurvatureWFS, det::AbstractDetector) = camera_frame(wfs)
@inline wfs_output_metadata(::AbstractWFS) = nothing
@inline wfs_output_metadata(wfs::CurvatureWFS) = wfs_output_metadata(wfs.params.readout_model, wfs)
@inline wfs_output_metadata(::CurvatureFrameReadout, wfs::CurvatureWFS) = nothing
@inline wfs_output_metadata(::CurvatureCountingReadout, wfs::CurvatureWFS) =
    CountingReadoutMetadata(:branch_by_channel, size(wfs.state.camera_frame), length(wfs.state.camera_frame))

simulation_interface(runtime::ClosedLoopRuntime) = SimulationInterface(runtime)
@inline backend(interface::SimulationInterface) = backend(interface.runtime)
@inline backend(interface::CompositeSimulationInterface) = backend(interface.interfaces[1])
simulation_interface(interface::SimulationInterface) = interface
simulation_interface(interface::CompositeSimulationInterface) = interface

@inline supports_prepared_runtime(::ShackHartmann{<:Diffractive}, ::AbstractSource) = true
@inline supports_prepared_runtime(::ShackHartmann{<:Diffractive}, ::Asterism) = true
@inline supports_prepared_runtime(::PyramidWFS, ::AbstractSource) = true
@inline supports_prepared_runtime(::PyramidWFS, ::Asterism) = true
@inline supports_prepared_runtime(::ZernikeWFS, ::AbstractSource) = true
@inline supports_prepared_runtime(::CurvatureWFS, ::AbstractSource) = true
@inline supports_detector_output(::ShackHartmann{<:Diffractive}, ::AbstractDetector) = true
@inline supports_detector_output(::PyramidWFS{<:Diffractive}, ::AbstractDetector) = true
@inline supports_detector_output(::BioEdgeWFS{<:Diffractive}, ::AbstractDetector) = true
@inline supports_detector_output(::ZernikeWFS, ::AbstractDetector) = true
@inline supports_detector_output(wfs::CurvatureWFS, det::AbstractDetector) = supports_detector_output(wfs.params.readout_model, det)
@inline supports_detector_output(::CurvatureFrameReadout, ::AbstractDetector) = true
@inline supports_detector_output(::CurvatureCountingReadout, ::AbstractDetector) = false
@inline supports_detector_output(::CurvatureCountingReadout, ::AbstractCountingDetector) = true
@inline supports_stacked_sources(::ShackHartmann, ::Asterism) = true
@inline supports_stacked_sources(::ShackHartmann, ::SpectralSource) = true
@inline supports_stacked_sources(::ShackHartmann, ::ExtendedSource) = true
@inline supports_stacked_sources(::PyramidWFS, ::Asterism) = true
@inline supports_stacked_sources(::PyramidWFS, ::SpectralSource) = true
@inline supports_stacked_sources(::PyramidWFS, ::ExtendedSource) = true
@inline supports_stacked_sources(::BioEdgeWFS, ::Asterism) = true
@inline supports_grouped_execution(::ShackHartmann{<:Diffractive}, ::Asterism) = true
@inline supports_grouped_execution(::ShackHartmann{<:Diffractive}, ::SpectralSource) = true
@inline supports_grouped_execution(::ShackHartmann{<:Diffractive}, ::ExtendedSource) = true
@inline supports_grouped_execution(::PyramidWFS{<:Diffractive}, ::Asterism) = true
@inline supports_grouped_execution(::PyramidWFS{<:Diffractive}, ::SpectralSource) = true
@inline supports_grouped_execution(::PyramidWFS{<:Diffractive}, ::ExtendedSource) = true
@inline supports_grouped_execution(::BioEdgeWFS{<:Diffractive}, ::Asterism) = true

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

@inline function prepare_runtime_wfs!(wfs::ShackHartmann{<:Diffractive}, tel::Telescope, src::SpectralSource)
    prepare_sampling!(wfs, tel, spectral_reference_source(src))
    ensure_sh_calibration!(wfs, tel, src)
    return wfs
end

@inline function prepare_runtime_wfs!(wfs::ShackHartmann{<:Diffractive}, tel::Telescope, ast::Asterism)
    isempty(ast.sources) && throw(InvalidConfiguration("asterism must contain at least one source"))
    prepare_sampling!(wfs, tel, ast.sources[1])
    ensure_sh_calibration!(wfs, tel, ast.sources[1])
    return wfs
end

@inline function prepare_runtime_wfs!(wfs::ZernikeWFS, tel::Telescope, src::AbstractSource)
    ensure_zernike_calibration!(wfs, tel, src)
    return wfs
end

@inline function prepare_runtime_wfs!(wfs::CurvatureWFS, tel::Telescope, src::AbstractSource)
    ensure_curvature_calibration!(wfs, tel, src)
    return wfs
end

@inline function prepare_runtime_wfs!(wfs::PyramidWFS, tel::Telescope, src::AbstractSource)
    prepare_pyramid_sampling!(wfs, tel)
    ensure_pyramid_calibration!(wfs, tel, src)
    return wfs
end

@inline function prepare_runtime_wfs!(wfs::PyramidWFS, tel::Telescope, src::SpectralSource)
    prepare_pyramid_sampling!(wfs, tel)
    ensure_pyramid_calibration!(wfs, tel, src)
    return wfs
end

@inline function prepare_runtime_wfs!(wfs::PyramidWFS, tel::Telescope, ast::Asterism)
    isempty(ast.sources) && throw(InvalidConfiguration("asterism must contain at least one source"))
    prepare_pyramid_sampling!(wfs, tel)
    ensure_pyramid_calibration!(wfs, tel, ast.sources[1])
    return wfs
end

@inline set_runtime_wfs_product_policy!(wfs::AbstractWFS, products::RuntimeProductRequirements) = wfs
@inline function set_runtime_wfs_product_policy!(wfs::ShackHartmann{<:Diffractive}, products::RuntimeProductRequirements)
    wfs.state.export_pixels_enabled = products.wfs_pixels
    return wfs
end

supports_prepared_runtime(runtime::ClosedLoopRuntime) = supports_prepared_runtime(runtime.wfs, runtime.src)
supports_detector_output(runtime::ClosedLoopRuntime) = requires_runtime_wfs_pixels(runtime) || requires_runtime_science_pixels(runtime)
supports_stacked_sources(runtime::ClosedLoopRuntime) = supports_stacked_sources(runtime.wfs, runtime.src)
supports_grouped_execution(::CompositeSimulationInterface) = true

function prepare!(runtime::ClosedLoopRuntime)
    prepare_runtime_wfs!(runtime.wfs, runtime.tel, runtime.src)
    set_runtime_wfs_product_policy!(runtime.wfs, runtime.products)
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
    products::GroupedRuntimeProductRequirements=default_grouped_runtime_products(interfaces...))
    length(interfaces) > 0 || throw(InvalidConfiguration("CompositeSimulationInterface requires at least one interface"))
    first_interface = interfaces[1]
    total_command = mapreduce(i -> length(simulation_command(i)), +, interfaces)
    total_slopes = mapreduce(i -> length(simulation_slopes(i)), +, interfaces)
    command = similar(simulation_command(first_interface), total_command)
    slopes = similar(simulation_slopes(first_interface), total_slopes)
    plan = grouped_runtime_product_plan(products, interfaces)
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
        typeof(products),
        typeof(wfs_stack),
        typeof(science_stack),
    }(
        interfaces,
        command,
        slopes,
        wfs_frames,
        science_frames,
        products,
        wfs_stack,
        science_stack,
    )
    return snapshot_outputs!(multi)
end

function CompositeSimulationInterface(; products::GroupedRuntimeProductRequirements=default_grouped_runtime_products())
    throw(InvalidConfiguration("CompositeSimulationInterface requires at least one interface"))
end

function CompositeSimulationInterface(runtimes::ClosedLoopRuntime...;
    products::GroupedRuntimeProductRequirements=default_grouped_runtime_products(map(SimulationInterface, runtimes)...))
    return CompositeSimulationInterface(map(SimulationInterface, runtimes)...; products=products)
end

function with_reconstructor(runtime::ClosedLoopRuntime, reconstructor)
    refreshed = ClosedLoopRuntime(
        runtime.simulation,
        reconstructor;
        wfs_detector=runtime.wfs_detector,
        science_detector=runtime.science_detector,
        rng=runtime.rng,
        profile=runtime.profile,
        products=runtime.products,
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

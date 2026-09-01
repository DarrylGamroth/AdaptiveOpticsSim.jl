#
# Prepared Pyramid WFS stages
#

"""
Run-immutable physical and numerical contract for one pyramid detector-plane
photon-arrival-rate map.
"""
struct PyramidOpticsPlan{P,O,C,S,L<:AbstractPreparedFourPupilLGS} <:
        AbstractWFSOpticsPlan
    propagation::P
    operating_modulation::O
    calibration_modulation::C
    source::S
    lgs_model::L
    propagation_revision::UInt
end

"""Exact live owner for one prepared pyramid optics execution."""
struct PreparedPyramidOptics{P,F,W,I,O,R,B,D}
    plan::P
    front_end::F
    workspace::W
    input::I
    output::O
    workspace_binding::R
    backend::B
    device::D
end

"""Run-immutable contract for one fixed Pyramid optics-product bundle."""
struct PyramidOpticsBundlePlan{P<:Tuple} <: AbstractWFSOpticsPlan
    plans::P
end

"""Exact live owner for one prepared Pyramid optics-product bundle."""
struct PreparedPyramidOpticsBundle{P,C<:Tuple,I,O}
    plan::P
    components::C
    input::I
    output::O
end

@inline wfs_optical_products(prepared::PreparedPyramidOptics) =
    prepared.output
@inline wfs_optical_products(prepared::PreparedPyramidOpticsBundle) =
    prepared.output

struct PyramidCalibrationBinding{T<:AbstractFloat,R,A}
    revision::UInt
    wavelength_m::T
    signature::UInt
    reference_signal::R
    valid_support::A
end

"""Run-immutable pyramid differential-estimation contract."""
struct PyramidEstimationPlan{E,P<:AbstractWFSMeasurementPath,C,S,T} <:
        AbstractWFSEstimationPlan
    params::E
    path::P
    calibration_binding::C
    source::S
    normalization_scale::T
end

"""Exact live owner for one prepared pyramid estimator."""
struct PreparedPyramidEstimator{P,W,ST,WS,PR,I,M,SB,WB,PB,B,D}
    plan::P
    sensor::W
    state::ST
    workspace::WS
    products::PR
    input::I
    measurement::M
    state_binding::SB
    workspace_binding::WB
    products_binding::PB
    backend::B
    device::D
end

@inline wfs_measurement_path(prepared::PreparedPyramidEstimator) =
    prepared.plan.path

@inline _pyramid_modulation_batch_workspace_binding(
    ::NoPyramidModulationBatchWorkspace) = ()

@inline function _pyramid_modulation_batch_workspace_binding(
    batch::PyramidModulationBatchWorkspace)
    return (batch.field_stack, batch.operating_weights,
        batch.calibration_weights, batch.fft_plan, batch.bfft_plan,
        batch.batch_size)
end


@inline function _pyramid_modulation_batch_workspace_binding(
    batch::PyramidShiftedMaskModulationWorkspace)
    return (
        batch.field_stack,
        batch.shifted_masks,
        batch.operating_weights,
        batch.axis_1_shifts_rad,
        batch.axis_2_shifts_rad,
        batch.bfft_plan,
        batch.batch_size,
    )
end

@inline function _pyramid_propagation_workspace_binding(workspace)
    return (workspace.field, workspace.focal_field, workspace.pupil_field,
        workspace.pyramid_mask, workspace.phasor, workspace.intensity,
        workspace.temp, workspace.scratch, workspace.asterism_stack,
        workspace.fft_plan, workspace.ifft_plan,
        workspace.elongation_kernel, workspace.lgs_kernel_fft,
        _pyramid_modulation_batch_workspace_binding(
            workspace.modulation_batch))
end

@inline modulated_wfs_propagation_storage(
    front_end::PyramidOpticalFrontEnd) =
    pyramid_propagation_workspace(front_end).field

@inline function _pyramid_estimator_state_binding(state)
    return (state.valid_mask, state.optical_gain, state.valid_i4q,
        state.reference_signal_2d)
end

@inline function _pyramid_estimator_workspace_binding(workspace)
    return (workspace.valid_i4q_host, workspace.valid_signal,
        workspace.valid_signal_indices, workspace.valid_signal_indices_host,
        workspace.valid_signal_count, workspace.valid_flux_sum_buffer,
        workspace.valid_flux_sum_host, workspace.valid_flux_i4q_host,
        workspace.flux_i4q, workspace.signal_2d)
end

@inline _pyramid_estimator_products_binding(products) = (products.slopes,)

@inline function pyramid_output_sampling_factor(
    front_end::PyramidOpticalFrontEnd, pupil_resolution::Int)
    pupil_resolution % front_end.pupil_samples == 0 ||
        throw(WFSPreparationError(:wfs_optics, :shape,
            "pyramid pupil resolution must be divisible by pupil_samples"))
    return div(pupil_resolution, front_end.pupil_samples) *
        front_end.binning
end

function pyramid_output_dimensions(front_end::PyramidOpticalFrontEnd,
    pupil_resolution::Int)
    factor = pyramid_output_sampling_factor(front_end, pupil_resolution)
    side = size(pyramid_propagation_workspace(front_end).intensity, 1)
    side % factor == 0 || throw(WFSPreparationError(:wfs_optics,
        :shape, "pyramid sampling does not evenly divide the detector plane"))
    output_side = div(side, factor)
    return (output_side, output_side)
end

@inline function _pyramid_front_end_wavelength(
    front_end::PyramidOpticalFrontEnd, input::PupilFunction)
    return modulated_input_wavelength(input, front_end.source)
end

@inline function _pyramid_front_end_wavelength(
    ::PyramidOpticalFrontEnd, input::ElectricField)
    return modulated_input_wavelength(input)
end

function _require_pyramid_source(front_end::PyramidOpticalFrontEnd,
    ::PupilFunction)
    source = front_end.source
    source === nothing && throw(WFSPreparationError(:wfs_optics,
        :radiometry, "pyramid WFS optics require a source for PupilFunction input"))
    return _require_single_pyramid_source(source)
end

@inline _require_single_pyramid_source(source) = source

function _require_single_pyramid_source(source::SpectralSource)
    throw(WFSPreparationError(
        :wfs_optics, :plane_count,
        "spectral pyramid optics require an OpticalProductBundle"))
end

function _require_single_pyramid_source(source::Asterism)
    throw(WFSPreparationError(:wfs_optics,
        :plane_count,
        "asterism pyramid optics require path-local pupil inputs"))
end

function _require_single_pyramid_source(source::ExtendedSource)
    throw(WFSPreparationError(
        :wfs_optics, :plane_count,
        "extended pyramid optics require path-local pupil inputs"))
end

function _require_pyramid_source(front_end::PyramidOpticalFrontEnd,
    ::ElectricField)
    front_end.source === nothing || throw(WFSPreparationError(
        :wfs_optics, :radiometry,
        "photon-rate ElectricField input must not also supply a source"))
    return nothing
end

function prepare_wfs_optics(front_end::PyramidOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField}, output::IntensityMap)
    require_modulated_wfs_input(input)
    _require_pyramid_source(front_end, input)
    resolution = input.metadata.dimensions[1]
    input.metadata.dimensions == (resolution, resolution) ||
        throw(WFSPreparationError(:wfs_optics, :shape,
            "pyramid pupil input must be square"))
    size(front_end.modulation.phases, 1) == resolution ||
        throw(WFSPreparationError(:wfs_optics, :shape,
            "pyramid modulation was prepared for another pupil resolution"))
    expected = pyramid_output_dimensions(front_end, resolution)
    wavelength_m = _pyramid_front_end_wavelength(front_end, input)
    require_four_pupil_rate_map(output, expected, wavelength_m)
    require_modulated_wfs_domains(front_end, input, output)
    eltype(pyramid_propagation_workspace(front_end).intensity) ===
        output.metadata.numeric_type || throw(WFSPreparationError(
            :wfs_optics, :numeric_type,
            "pyramid output precision differs from prepared propagation"))
    lgs_model = prepare_four_pupil_lgs(front_end.source, input, front_end)
    propagation = front_end.propagation
    propagation_plan = pyramid_propagation_plan(propagation)
    workspace = pyramid_propagation_workspace(propagation)
    plan = PyramidOpticsPlan(propagation_plan, front_end.modulation,
        front_end.calibration_modulation,
        front_end.source, lgs_model, workspace.revision)
    return PreparedPyramidOptics(plan, front_end, workspace, input, output,
        _pyramid_propagation_workspace_binding(workspace),
        input.metadata.backend, input.metadata.device)
end

function prepare_wfs_optics(front_end::PyramidOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField},
    output::OpticalProductBundle)
    return prepare_pyramid_optical_bundle(front_end, input, output,
        front_end.source)
end

function prepare_pyramid_optical_bundle(front_end::PyramidOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField}, output::OpticalProductBundle,
    source::SpectralSource)
    samples = spectral_bundle(source).samples
    length(output) == length(samples) || throw(WFSPreparationError(
        :wfs_optics, :plane_count,
        "pyramid spectral output count does not match the source"))
    T = eltype(pyramid_propagation_workspace(front_end).intensity)
    plans = ntuple(length(samples)) do index
        sample = samples[index]
        component = FourPupilSpectralComponent(source.source,
            T(sample.wavelength),
            T(photon_irradiance(source)) * T(sample.weight))
        prepare_wfs_optics(
            pyramid_front_end_with_source(front_end, component), input,
            output[index])
    end
    return PreparedPyramidOpticsBundle(
        PyramidOpticsBundlePlan(map(component -> component.plan, plans)),
        plans, input, output)
end

function prepare_wfs_optics(front_end::PyramidOpticalFrontEnd,
    inputs::Union{Tuple,AbstractVector}, output::OpticalProductBundle)
    return prepare_pyramid_optical_bundle(front_end, inputs, output,
        front_end.source)
end

function prepare_pyramid_optical_bundle(front_end::PyramidOpticalFrontEnd,
    inputs::Union{Tuple,AbstractVector}, output::OpticalProductBundle,
    source::Union{Asterism,ExtendedSource})
    sources = four_pupil_path_sources(source)
    length(inputs) == length(sources) || throw(WFSPreparationError(
        :wfs_optics, :plane_count,
        "pyramid path-local pupil count does not match the source count"))
    length(output) == length(sources) || throw(WFSPreparationError(
        :wfs_optics, :plane_count,
        "pyramid path-local output count does not match the source count"))
    isempty(sources) && throw(WFSPreparationError(:wfs_optics,
        :plane_count, "pyramid path-local source collection is empty"))
    plans = ntuple(length(sources)) do index
        prepare_wfs_optics(
            pyramid_front_end_with_source(front_end, sources[index]),
            inputs[index], output[index])
    end
    return PreparedPyramidOpticsBundle(
        PyramidOpticsBundlePlan(map(component -> component.plan, plans)),
        plans, inputs, output)
end

function prepare_pyramid_optical_bundle(front_end, input, output, source)
    throw(WFSPreparationError(
        :wfs_optics, :plane_count,
        "pyramid product bundles require a spectral or path-expanded source"))
end

function _pyramid_native_rate!(front_end::PyramidOpticalFrontEnd,
    input::PupilFunction)
    propagation = pyramid_propagation_workspace(front_end)
    source = front_end.source
    resolution = size(input.opd, 1)
    pad = size(propagation.field, 1)
    offset = div(pad - resolution, 2)
    T = eltype(propagation.intensity)
    amplitude_scale = sqrt(T(photon_irradiance(source)) *
        T(input.metadata.sampling[1] * input.metadata.sampling[2]))
    opd_to_cycles = T(2) / T(wavelength(source))
    fill!(propagation.intensity, zero(T))
    _pyramid_pupil_modulation_batch!(
        propagation.modulation_batch,
        propagation.intensity,
        front_end,
        input.amplitude,
        input.opd,
        front_end.modulation,
        amplitude_scale,
        opd_to_cycles,
        offset,
        resolution,
    ) && return propagation.intensity
    @inbounds for point in 1:modulation_point_count(front_end.modulation)
        fill!(propagation.field, zero(eltype(propagation.field)))
        weight = front_end.modulation.amplitude_weights[point]
        @views @. propagation.field[offset+1:offset+resolution,
            offset+1:offset+resolution] = amplitude_scale * weight *
            input.amplitude * front_end.modulation.phases[:, :, point] *
            cispi(opd_to_cycles * input.opd)
        copyto!(propagation.focal_field, propagation.field)
        accumulate_pyramid_focal_intensity!(propagation.intensity, front_end)
    end
    return propagation.intensity
end

function _pyramid_native_rate!(front_end::PyramidOpticalFrontEnd,
    input::ElectricField)
    propagation = pyramid_propagation_workspace(front_end)
    resolution = size(input.values, 1)
    pad = size(propagation.field, 1)
    offset = div(pad - resolution, 2)
    T = eltype(propagation.intensity)
    fill!(propagation.intensity, zero(T))
    _pyramid_electric_field_modulation_batch!(
        propagation.modulation_batch,
        propagation.intensity,
        front_end,
        input.values,
        front_end.modulation,
        offset,
        resolution,
    ) && return propagation.intensity
    @inbounds for point in 1:modulation_point_count(front_end.modulation)
        fill!(propagation.field, zero(eltype(propagation.field)))
        weight = front_end.modulation.amplitude_weights[point]
        @views @. propagation.field[offset+1:offset+resolution,
            offset+1:offset+resolution] = weight * input.values *
            front_end.modulation.phases[:, :, point]
        copyto!(propagation.focal_field, propagation.field)
        accumulate_pyramid_focal_intensity!(propagation.intensity, front_end)
    end
    return propagation.intensity
end

function _apply_prepared_pyramid_lgs!(
    plan::PreparedPyramidOptics)
    propagation = plan.workspace
    apply_prepared_four_pupil_lgs!(plan.plan.lgs_model, propagation.intensity,
        propagation.scratch, propagation.focal_field,
        propagation.fft_plan, propagation.pupil_field,
        propagation.ifft_plan)
    return propagation.intensity
end

function form_wfs_optical_products!(output::IntensityMap,
    input::Union{PupilFunction,ElectricField},
    plan::PreparedPyramidOptics)
    validate_wfs_optics_binding(output, input, plan)
    native = _pyramid_native_rate!(plan.front_end, input)
    _apply_prepared_pyramid_lgs!(plan)
    factor = pyramid_output_sampling_factor(plan.front_end,
        input.metadata.dimensions[1])
    bin2d!(output.values, native, factor)
    return output
end

@inline _enqueue_prepared_pyramid_lgs!(::NoPreparedFourPupilLGS) = nothing

function enqueue_wfs_optical_products!(
    output::IntensityMap,
    input::PupilFunction,
    plan::PreparedPyramidOptics,
)
    propagation = plan.workspace
    native = _enqueue_pyramid_rate_batch!(
        execution_style(propagation.intensity),
        propagation.intensity,
        plan.front_end,
        input,
        plan.front_end.source,
        propagation.modulation_batch,
    )
    _enqueue_prepared_pyramid_lgs!(plan.plan.lgs_model)
    factor = pyramid_output_sampling_factor(
        plan.front_end,
        input.metadata.dimensions[1],
    )
    style = execution_style(output.values)
    launch_kernel_async!(style, bin2d_kernel!, output.values, native, factor,
        size(output.values, 1), size(output.values, 2);
        ndrange=size(output.values))
    return output
end

function validate_wfs_optics_binding(output::IntensityMap,
    input::Union{PupilFunction,ElectricField},
    plan::PreparedPyramidOptics)
    output === plan.output && input === plan.input ||
        throw(WFSPreparationError(:wfs_optics, :prepared_binding,
            "pyramid optical products do not match prepared storage"))
    workspace = pyramid_propagation_workspace(plan.front_end)
    workspace === plan.workspace &&
        _pyramid_propagation_workspace_binding(workspace) ===
            plan.workspace_binding &&
        workspace.revision == plan.plan.propagation_revision ||
        throw(WFSPreparationError(:wfs_optics, :prepared_binding,
            "pyramid propagation workspace changed after preparation"))
    plan.front_end.propagation.plan === plan.plan.propagation &&
        plan.front_end.phase_mask === plan.plan.propagation.phase_mask &&
        plan.front_end.modulation === plan.plan.operating_modulation &&
        plan.front_end.calibration_modulation ===
            plan.plan.calibration_modulation &&
        plan.front_end.source === plan.plan.source ||
        throw(WFSPreparationError(:wfs_optics, :prepared_binding,
            "pyramid optics definition changed after preparation"))
    return nothing
end

function form_wfs_optical_products!(output::OpticalProductBundle, input,
    plan::PreparedPyramidOpticsBundle)
    validate_wfs_optics_binding(output, input, plan)
    return form_four_pupil_bundle!(output, input, plan.components)
end

function validate_wfs_optics_binding(
    output::OpticalProductBundle, input,
    plan::PreparedPyramidOpticsBundle)
    output === plan.output && input === plan.input ||
        throw(WFSPreparationError(:wfs_optics, :prepared_binding,
            "pyramid spectral products do not match prepared storage"))
    validate_four_pupil_bundle_binding(output, input, plan.components)
    return nothing
end

function pyramid_rate_map(sensor::PyramidWFS{<:Diffractive},
    inputs::Union{Tuple,AbstractVector}, source)
    return pyramid_rate_map(PyramidOpticalFrontEnd(sensor, source), inputs)
end

function pyramid_rate_map(front_end::PyramidOpticalFrontEnd,
    inputs::Union{Tuple,AbstractVector})
    return pyramid_path_rate_bundle(front_end, inputs, front_end.source)
end

function pyramid_path_rate_bundle(front_end::PyramidOpticalFrontEnd,
    inputs::Union{Tuple,AbstractVector},
    source::Union{Asterism,ExtendedSource})
    sources = four_pupil_path_sources(source)
    length(inputs) == length(sources) || throw(WFSPreparationError(
        :wfs_optics, :plane_count,
        "pyramid path-local pupil count does not match the source count"))
    isempty(sources) && throw(WFSPreparationError(:wfs_optics,
        :plane_count, "pyramid path-local source collection is empty"))
    first_map = pyramid_rate_map(
        pyramid_front_end_with_source(front_end, sources[1]), inputs[1])
    maps = Vector{typeof(first_map)}(undef, length(sources))
    maps[1] = first_map
    @inbounds for index in 2:length(sources)
        maps[index] = pyramid_rate_map(
            pyramid_front_end_with_source(front_end, sources[index]),
            inputs[index])
    end
    return OpticalProductBundle(maps)
end

function pyramid_path_rate_bundle(front_end, inputs, source)
    throw(WFSPreparationError(:wfs_optics, :plane_count,
        "path-local pyramid inputs require an Asterism or ExtendedSource"))
end

function pyramid_rate_map(sensor::PyramidWFS{<:Diffractive},
    input::Union{PupilFunction,ElectricField}, source=nothing)
    return pyramid_rate_map(PyramidOpticalFrontEnd(sensor, source), input)
end

function pyramid_rate_map(front_end::PyramidOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField})
    return _pyramid_rate_map(front_end, input, front_end.source)
end

@inline _pyramid_rate_map(front_end::PyramidOpticalFrontEnd, input,
    source::SpectralSource) =
    _pyramid_spectral_rate_bundle(front_end, input, source)

function _pyramid_rate_map(front_end::PyramidOpticalFrontEnd, input,
    source::Union{Asterism,ExtendedSource})
    throw(WFSPreparationError(:wfs_optics, :plane_count,
        "path-expanded pyramid sources require path-local pupil inputs"))
end

function _pyramid_rate_map(
    front_end::PyramidOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField},
    values::AbstractMatrix,
)
    wavelength_m = _pyramid_front_end_wavelength(front_end, input)
    dimensions = pyramid_output_dimensions(front_end,
        input.metadata.dimensions[1])
    T = eltype(pyramid_propagation_workspace(front_end).intensity)
    Base.require_one_based_indexing(values)
    eltype(values) === T || throw(WFSPreparationError(
        :wfs_optics,
        :numeric_type,
        "pyramid rate output numeric type must match its optics workspace",
    ))
    size(values) == dimensions || throw(WFSPreparationError(
        :wfs_optics,
        :shape,
        "pyramid rate output must match the configured four-pupil frame",
    ))
    compute_device(values) == compute_device(_modulated_input_storage(input)) ||
        throw(WFSPreparationError(
            :wfs_optics,
            :device,
            "pyramid input and rate output occupy different compute devices",
        ))
    factor = pyramid_output_sampling_factor(front_end,
        input.metadata.dimensions[1])
    normalized_sampling = T(factor / input.metadata.dimensions[1])
    metadata = OpticalPlaneMetadata(DetectorPlane(), values;
        coordinate_domain=NormalizedPupilCoordinates(),
        sampling=(normalized_sampling, normalized_sampling),
        spectral=MonochromaticChannel(T(wavelength_m)),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    return IntensityMap(metadata, values)
end

"""
    pyramid_rate_map(front_end, input, values)

Bind an already allocated caller-owned, one-based output matrix as the
detector-plane photon-rate frame for a prepared Pyramid optical front end.
The matrix must match the front end's numeric type, shape, and compute device.
Repeated optical execution remains allocation-free.
"""
function pyramid_rate_map(
    front_end::PyramidOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField},
    values::AbstractMatrix,
)
    return _pyramid_rate_map(front_end, input, values)
end

function _pyramid_rate_map(front_end::PyramidOpticalFrontEnd, input, source)
    dimensions = pyramid_output_dimensions(front_end,
        input.metadata.dimensions[1])
    T = eltype(pyramid_propagation_workspace(front_end).intensity)
    values = similar(_modulated_input_storage(input), T, dimensions...)
    fill!(values, zero(T))
    return _pyramid_rate_map(front_end, input, values)
end

@inline _require_pyramid_estimation_source(::WFSNormalization, source) =
    nothing

function _require_pyramid_estimation_source(
    ::IncidenceFluxNormalization, ::Nothing)
    throw(WFSPreparationError(:estimation, :radiometry,
        "incidence-normalized pyramid estimation requires a source"))
end

function _pyramid_spectral_rate_bundle(front_end::PyramidOpticalFrontEnd,
    input, source::SpectralSource)
    samples = spectral_bundle(source).samples
    T = eltype(pyramid_propagation_workspace(front_end).intensity)
    function component_map(sample)
        component = FourPupilSpectralComponent(source.source,
            T(sample.wavelength),
            T(photon_irradiance(source)) * T(sample.weight))
        component_front_end = pyramid_front_end_with_source(front_end,
            component)
        return pyramid_rate_map(component_front_end, input)
    end
    first_map = component_map(first(samples))
    maps = Vector{typeof(first_map)}(undef, length(samples))
    maps[1] = first_map
    @inbounds for index in 2:length(samples)
        maps[index] = component_map(samples[index])
    end
    return OpticalProductBundle(maps)
end

function _pyramid_calibration_binding(sensor::PyramidWFS)
    state = sensor.estimator.state
    state.calibrated || throw(WFSPreparationError(:estimation, :estimator,
        "pyramid estimation requires explicit calibration"))
    return PyramidCalibrationBinding(state.calibration_revision,
        state.calibration_wavelength,
        state.calibration_signature, state.reference_signal_2d,
        state.valid_i4q)
end

function _require_pyramid_calibration(sensor::PyramidWFS,
    binding::PyramidCalibrationBinding)
    state = sensor.estimator.state
    state.calibrated &&
        state.calibration_revision == binding.revision &&
        state.calibration_wavelength == binding.wavelength_m &&
        state.calibration_signature == binding.signature &&
        state.reference_signal_2d === binding.reference_signal &&
        state.valid_i4q === binding.valid_support ||
        throw(WFSPreparationError(:estimation, :prepared_binding,
            "pyramid calibration changed after estimator preparation"))
    return nothing
end

function _prepare_pyramid_estimator_owner(sensor::PyramidWFS, input,
    measurement::WFSMeasurement, path::AbstractWFSMeasurementPath,
    calibration_binding, source, normalization_scale)
    state = pyramid_estimator_state(sensor)
    workspace = pyramid_estimator_workspace(sensor)
    products = pyramid_estimator_products(sensor)
    plan = PyramidEstimationPlan(sensor.estimator.params, path,
        calibration_binding, source, normalization_scale)
    return PreparedPyramidEstimator(plan, sensor, state, workspace, products,
        input, measurement, _pyramid_estimator_state_binding(state),
        _pyramid_estimator_workspace_binding(workspace),
        _pyramid_estimator_products_binding(products),
        measurement.metadata.backend, measurement.metadata.device)
end

function _require_pyramid_estimation_geometry(sensor::PyramidWFS,
    frame_size::Int)
    iseven(frame_size) || throw(WFSPreparationError(:estimation, :shape,
        "pyramid observations require an even detector-frame size"))
    nominal = pyramid_acquisition_workspace(sensor).nominal_detector_resolution
    binning = pyramid_acquisition_plan(sensor).binning
    nominal > 0 || throw(WFSPreparationError(:estimation, :shape,
        "pyramid nominal detector resolution has not been prepared"))
    nominal % binning == 0 || throw(WFSPreparationError(:estimation, :shape,
        "pyramid binning does not divide the nominal detector resolution"))
    sampled_size = div(nominal, binning)
    sampled_size % frame_size == 0 || throw(WFSPreparationError(
        :estimation, :shape,
        "detector sampling does not evenly divide the pyramid frame"))
    total_sampling = binning * div(sampled_size, frame_size)
    n_pixels, half_separation, edge_padding = pyramid_sampled_geometry(
        sensor.estimator.params.pupil_samples,
        sensor.front_end.phase_mask.n_pix_separation,
        sensor.front_end.phase_mask.n_pix_edge, total_sampling)
    n_pixels >= 1 || throw(WFSPreparationError(:estimation, :shape,
        "detector sampling removed every pyramid pupil sample"))
    if sensor.front_end.phase_mask.n_pix_separation === nothing
        frame_size >= 2 * n_pixels || throw(WFSPreparationError(
            :estimation, :shape,
            "pyramid frame does not contain four complete pupil images"))
    else
        frame_size == 2 * (n_pixels + half_separation + edge_padding) ||
            throw(WFSPreparationError(:estimation, :shape,
                "pyramid frame does not preserve the configured pupil geometry"))
    end
    return nothing
end

function prepare_wfs_estimation(sensor::PyramidWFS{<:Diffractive},
    observation::WFSObservation, measurement::WFSMeasurement;
    source=nothing, normalization_scale::Real=1)
    validate_wfs_observation(observation)
    validate_wfs_measurement(measurement)
    isequal(observation.metadata.layout, :four_pupil_mosaic) ||
        throw(WFSPreparationError(:estimation, :detector_mapping,
            "pyramid estimator requires :four_pupil_mosaic layout"))
    isequal(measurement.metadata.kind, :differential_slopes) ||
        throw(WFSPreparationError(:estimation, :estimator,
            "pyramid measurement kind must be :differential_slopes"))
    isequal(measurement.units, :dimensionless) ||
        throw(WFSPreparationError(:estimation, :units,
            "pyramid differential slopes are dimensionless"))
    frame_size = _require_real_square_wfs_observation(observation,
        "pyramid")
    measurement.metadata.numeric_type <: AbstractFloat ||
        throw(WFSPreparationError(:estimation, :numeric_type,
            "pyramid measurement storage must be floating point"))
    _require_wfs_storage_domain(:estimation, observation.metadata,
        pyramid_estimator_workspace(sensor).signal_2d, "pyramid observation")
    _require_wfs_storage_domain(:estimation, measurement.metadata,
        pyramid_estimator_products(sensor).slopes, "pyramid measurement")
    _require_pyramid_estimation_geometry(sensor, frame_size)
    resize_pyramid_signal_buffers!(sensor, frame_size)
    size(measurement.storage) == size(pyramid_estimator_products(sensor).slopes) ||
        throw(WFSPreparationError(:estimation, :shape,
            "pyramid measurement storage has the wrong slope shape"))
    _require_pyramid_estimation_source(
        sensor.estimator.params.normalization, source)
    scale = eltype(pyramid_estimator_products(sensor).slopes)(normalization_scale)
    isfinite(scale) && scale >= zero(scale) || throw(WFSPreparationError(
        :estimation, :radiometry,
        "pyramid normalization scale must be finite and nonnegative"))
    binding = _pyramid_calibration_binding(sensor)
    return _prepare_pyramid_estimator_owner(sensor, observation, measurement,
        AcquiredObservationPath(), binding, source, scale)
end

function estimate_wfs_measurement!(measurement::WFSMeasurement,
    observation::WFSObservation,
    plan::PreparedPyramidEstimator)
    validate_wfs_estimation_binding(measurement, observation, plan)
    sensor = plan.sensor
    _require_pyramid_calibration(sensor, plan.plan.calibration_binding)
    pyramid_signal!(execution_style(observation.storage), sensor,
        observation.storage, plan.plan.source, plan.plan.normalization_scale)
    slopes = pyramid_estimator_products(sensor).slopes
    @. slopes *= sensor.estimator.state.optical_gain
    copyto!(measurement.storage, slopes)
    return measurement
end

function validate_wfs_estimation_binding(measurement::WFSMeasurement, input,
    plan::PreparedPyramidEstimator)
    measurement === plan.measurement && input === plan.input || throw(
        WFSPreparationError(:estimation, :prepared_binding,
            "pyramid estimator storage does not match its plan"))
    state = pyramid_estimator_state(plan.sensor)
    workspace = pyramid_estimator_workspace(plan.sensor)
    products = pyramid_estimator_products(plan.sensor)
    state === plan.state && workspace === plan.workspace &&
        products === plan.products &&
        _pyramid_estimator_state_binding(state) === plan.state_binding &&
        _pyramid_estimator_workspace_binding(workspace) ===
            plan.workspace_binding &&
        _pyramid_estimator_products_binding(products) ===
            plan.products_binding || throw(WFSPreparationError(
                :estimation, :prepared_binding,
                "pyramid estimator state, workspace, or products changed after preparation"))
    return nothing
end

function prepare_wfs_estimation(sensor::PyramidWFS{<:Geometric},
    input::PupilFunction, measurement::WFSMeasurement)
    require_modulated_wfs_input(input)
    validate_wfs_measurement(measurement)
    input.metadata.dimensions == (sensor.estimator.params.pupil_resolution,
        sensor.estimator.params.pupil_resolution) ||
        throw(WFSPreparationError(:estimation, :shape,
            "geometric pyramid input has the wrong pupil dimensions"))
    isequal(measurement.metadata.kind, :geometric_slopes) ||
        throw(WFSPreparationError(:estimation, :estimator,
            "geometric pyramid measurement kind must be :geometric_slopes"))
    isequal(measurement.units, :metre) || throw(WFSPreparationError(
        :estimation, :units,
        "geometric pyramid OPD differences are expressed in metres"))
    measurement.metadata.numeric_type <: AbstractFloat ||
        throw(WFSPreparationError(:estimation, :numeric_type,
            "geometric pyramid measurement storage must be floating point"))
    _require_wfs_storage_domain(:estimation, input.metadata,
        pyramid_estimator_products(sensor).slopes, "geometric pyramid input")
    _require_wfs_storage_domain(:estimation, measurement.metadata,
        pyramid_estimator_products(sensor).slopes, "geometric pyramid measurement")
    size(measurement.storage) == size(pyramid_estimator_products(sensor).slopes) ||
        throw(WFSPreparationError(:estimation, :shape,
            "geometric pyramid measurement has the wrong slope shape"))
    return _prepare_pyramid_estimator_owner(sensor, input, measurement,
        DirectMeasurementPath(), nothing, nothing,
        one(eltype(pyramid_estimator_products(sensor).slopes)))
end

function estimate_wfs_measurement!(measurement::WFSMeasurement,
    input::PupilFunction,
    plan::PreparedPyramidEstimator)
    validate_wfs_estimation_binding(measurement, input, plan)
    sensor = plan.sensor
    state = sensor.estimator.state
    products = pyramid_estimator_products(sensor)
    geometric_slopes!(products.slopes, input.opd, state.valid_mask)
    gain = inv(1 + sensor.estimator.params.geometric_modulation_radius)
    @. products.slopes = gain * products.slopes * state.optical_gain
    copyto!(measurement.storage, products.slopes)
    return measurement
end

function set_pyramid_calibration!(sensor::PyramidWFS,
    reference::AbstractMatrix; wavelength_m::Real,
    signature::UInt=UInt(0), valid_support=nothing)
    state = sensor.estimator.state
    size(reference) == size(state.reference_signal_2d) ||
        throw(DimensionMismatchError(
            "pyramid reference dimensions do not match estimator storage"))
    require_same_backend(state.reference_signal_2d, reference)
    reference_host = Array(reference)
    all(isfinite, reference_host) || throw(InvalidConfiguration(
        "pyramid calibration reference must contain only finite values"))
    wavelength_value = eltype(pyramid_estimator_products(sensor).slopes)(
        wavelength_m)
    isfinite(wavelength_value) && wavelength_value > zero(wavelength_value) ||
        throw(InvalidConfiguration(
            "pyramid calibration wavelength must be finite and positive"))
    support_host = _prepare_pyramid_calibration_support(sensor, valid_support)
    copyto!(state.reference_signal_2d, reference_host)
    if support_host === nothing
        fill!(state.valid_i4q, true)
    else
        copyto!(state.valid_i4q, support_host)
    end
    update_pyramid_valid_signal!(sensor)
    update_pyramid_valid_signal_indices!(sensor)
    resize_pyramid_slope_buffers!(sensor)
    state.calibration_wavelength = wavelength_value
    state.calibration_signature = signature
    state.calibrated = true
    state.calibration_revision += UInt(1)
    return sensor
end

function _prepare_pyramid_calibration_support(sensor::PyramidWFS, ::Nothing)
    if !iszero(sensor.estimator.params.light_ratio)
        throw(InvalidConfiguration(
            "nonzero pyramid light_ratio requires explicit valid_support"))
    end
    return nothing
end

function _prepare_pyramid_calibration_support(sensor::PyramidWFS,
    valid_support::AbstractMatrix{Bool})
    state = sensor.estimator.state
    size(valid_support) == size(state.valid_i4q) ||
        throw(DimensionMismatchError(
            "pyramid calibration support has the wrong dimensions"))
    require_same_backend(state.valid_i4q, valid_support)
    support_host = Array(valid_support)
    any(support_host) || throw(InvalidConfiguration(
        "pyramid calibration support must select at least one sample"))
    return support_host
end

function _prepare_pyramid_calibration_support(sensor::PyramidWFS,
    valid_support)
    throw(InvalidConfiguration(
        "pyramid calibration support must be a Boolean matrix"))
end

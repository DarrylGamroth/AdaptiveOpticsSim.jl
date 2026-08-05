#
# Prepared Bi-O-edge WFS stages
#

"""
Run-immutable physical and numerical contract for one Bi-O-edge detector-plane
photon-arrival-rate map.
"""
struct BiOEdgeOpticsPlan{P,O,C,S,L<:AbstractPreparedFourPupilLGS} <:
        AbstractWFSOpticsPlan
    propagation::P
    operating_modulation::O
    calibration_modulation::C
    source::S
    lgs_model::L
    propagation_revision::UInt
end

"""Exact live owner for one prepared Bi-O-edge optics execution."""
struct PreparedBiOEdgeOptics{P,F,W,I,O,R,B,D}
    plan::P
    front_end::F
    workspace::W
    input::I
    output::O
    workspace_binding::R
    backend::B
    device::D
end

"""Run-immutable contract for one fixed Bi-O-edge optics-product bundle."""
struct BiOEdgeOpticsBundlePlan{P<:Tuple} <: AbstractWFSOpticsPlan
    plans::P
end

"""Exact live owner for one prepared Bi-O-edge optics-product bundle."""
struct PreparedBiOEdgeOpticsBundle{P,C<:Tuple,I,O}
    plan::P
    components::C
    input::I
    output::O
end

@inline wfs_optical_products(prepared::PreparedBiOEdgeOptics) =
    prepared.output
@inline wfs_optical_products(prepared::PreparedBiOEdgeOpticsBundle) =
    prepared.output

struct BiOEdgeCalibrationBinding{T<:AbstractFloat,R,A}
    revision::UInt
    wavelength_m::T
    signature::UInt
    reference_signal::R
    valid_support::A
end

"""Run-immutable Bi-O-edge differential-estimation contract."""
struct BiOEdgeEstimationPlan{E,P<:AbstractWFSMeasurementPath,C,S,T} <:
        AbstractWFSEstimationPlan
    params::E
    path::P
    calibration_binding::C
    source::S
    normalization_scale::T
end

"""Exact live owner for one prepared Bi-O-edge estimator."""
struct PreparedBiOEdgeEstimator{P,W,ST,WS,PR,I,M,SB,WB,PB,B,D}
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

@inline wfs_measurement_path(prepared::PreparedBiOEdgeEstimator) =
    prepared.plan.path

@inline function _bi_o_edge_propagation_workspace_binding(workspace)
    return (workspace.field, workspace.focal_field, workspace.pupil_field,
        workspace.bi_o_edge_masks, workspace.phasor, workspace.intensity,
        workspace.temp, workspace.scratch, workspace.asterism_stack,
        workspace.fft_buffer, workspace.fft_plan, workspace.ifft_plan,
        workspace.elongation_kernel, workspace.lgs_kernel_fft)
end

@inline modulated_wfs_propagation_storage(
    front_end::BiOEdgeOpticalFrontEnd) =
    bi_o_edge_propagation_workspace(front_end).field

@inline function _bi_o_edge_estimator_state_binding(state)
    return (state.valid_mask, state.edge_mask, state.optical_gain,
        state.valid_i4q, state.reference_signal_2d)
end

@inline function _bi_o_edge_estimator_workspace_binding(workspace)
    return (workspace.valid_i4q_host, workspace.valid_signal,
        workspace.valid_signal_indices, workspace.valid_signal_indices_host,
        workspace.valid_signal_count, workspace.valid_flux_sum_buffer,
        workspace.valid_flux_sum_host, workspace.valid_flux_i4q_host,
        workspace.flux_i4q, workspace.signal_2d, workspace.binned_phase,
        workspace.edge_mask_binned, workspace.binned_resolution)
end

@inline _bi_o_edge_estimator_products_binding(products) = (products.slopes,)

@inline function bi_o_edge_output_sampling_factor(
    front_end::BiOEdgeOpticalFrontEnd, pupil_resolution::Int)
    pupil_resolution % front_end.pupil_samples == 0 ||
        throw(WFSPreparationError(:wfs_optics, :shape,
            "Bi-O-edge pupil resolution must be divisible by pupil_samples"))
    pupil_sample = div(pupil_resolution, front_end.pupil_samples)
    return front_end.binning == 1 ? pupil_sample :
        2 * pupil_sample * front_end.binning
end

function bi_o_edge_output_dimensions(front_end::BiOEdgeOpticalFrontEnd,
    pupil_resolution::Int)
    factor = bi_o_edge_output_sampling_factor(front_end, pupil_resolution)
    side = size(bi_o_edge_propagation_workspace(front_end).intensity, 1)
    side % factor == 0 || throw(WFSPreparationError(:wfs_optics,
        :shape, "Bi-O-edge sampling does not evenly divide the detector plane"))
    output_side = div(side, factor)
    return (output_side, output_side)
end

@inline function _bi_o_edge_front_end_wavelength(
    front_end::BiOEdgeOpticalFrontEnd, input::PupilFunction)
    return modulated_input_wavelength(input, front_end.source)
end

@inline function _bi_o_edge_front_end_wavelength(
    ::BiOEdgeOpticalFrontEnd, input::ElectricField)
    return modulated_input_wavelength(input)
end

function _require_bi_o_edge_source(front_end::BiOEdgeOpticalFrontEnd,
    ::PupilFunction)
    source = front_end.source
    source === nothing && throw(WFSPreparationError(:wfs_optics,
        :radiometry, "Bi-O-edge WFS optics require a source for PupilFunction input"))
    return _require_single_bi_o_edge_source(source)
end

@inline _require_single_bi_o_edge_source(source) = source

function _require_single_bi_o_edge_source(source::SpectralSource)
    throw(WFSPreparationError(
        :wfs_optics, :plane_count,
        "spectral Bi-O-edge optics require an OpticalProductBundle"))
end

function _require_single_bi_o_edge_source(source::Asterism)
    throw(WFSPreparationError(:wfs_optics,
        :plane_count,
        "asterism Bi-O-edge optics require path-local pupil inputs"))
end

function _require_single_bi_o_edge_source(source::ExtendedSource)
    throw(WFSPreparationError(
        :wfs_optics, :plane_count,
        "extended Bi-O-edge optics require path-local pupil inputs"))
end

function _require_bi_o_edge_source(front_end::BiOEdgeOpticalFrontEnd,
    ::ElectricField)
    front_end.source === nothing || throw(WFSPreparationError(
        :wfs_optics, :radiometry,
        "photon-rate ElectricField input must not also supply a source"))
    return nothing
end

function prepare_wfs_optics(front_end::BiOEdgeOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField}, output::IntensityMap)
    require_modulated_wfs_input(input)
    _require_bi_o_edge_source(front_end, input)
    resolution = input.metadata.dimensions[1]
    input.metadata.dimensions == (resolution, resolution) ||
        throw(WFSPreparationError(:wfs_optics, :shape,
            "Bi-O-edge pupil input must be square"))
    size(front_end.modulation.phases, 1) == resolution ||
        throw(WFSPreparationError(:wfs_optics, :shape,
            "Bi-O-edge modulation was prepared for another pupil resolution"))
    expected = bi_o_edge_output_dimensions(front_end, resolution)
    wavelength_m = _bi_o_edge_front_end_wavelength(front_end, input)
    require_four_pupil_rate_map(output, expected, wavelength_m)
    require_modulated_wfs_domains(front_end, input, output)
    eltype(bi_o_edge_propagation_workspace(front_end).intensity) ===
        output.metadata.numeric_type || throw(WFSPreparationError(
            :wfs_optics, :numeric_type,
            "Bi-O-edge output precision differs from prepared propagation"))
    lgs_model = prepare_four_pupil_lgs(front_end.source, input, front_end)
    propagation = front_end.propagation
    propagation_plan = bi_o_edge_propagation_plan(propagation)
    workspace = bi_o_edge_propagation_workspace(propagation)
    plan = BiOEdgeOpticsPlan(propagation_plan, front_end.modulation,
        front_end.calibration_modulation,
        front_end.source, lgs_model, workspace.revision)
    return PreparedBiOEdgeOptics(plan, front_end, workspace, input, output,
        _bi_o_edge_propagation_workspace_binding(workspace),
        input.metadata.backend, input.metadata.device)
end

function prepare_wfs_optics(front_end::BiOEdgeOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField},
    output::OpticalProductBundle)
    return prepare_bi_o_edge_optical_bundle(front_end, input, output,
        front_end.source)
end

function prepare_bi_o_edge_optical_bundle(front_end::BiOEdgeOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField}, output::OpticalProductBundle,
    source::SpectralSource)
    samples = spectral_bundle(source).samples
    length(output) == length(samples) || throw(WFSPreparationError(
        :wfs_optics, :plane_count,
        "Bi-O-edge spectral output count does not match the source"))
    T = eltype(bi_o_edge_propagation_workspace(front_end).intensity)
    plans = ntuple(length(samples)) do index
        sample = samples[index]
        component = FourPupilSpectralComponent(source.source,
            T(sample.wavelength),
            T(photon_irradiance(source)) * T(sample.weight))
        prepare_wfs_optics(
            bi_o_edge_front_end_with_source(front_end, component), input,
            output[index])
    end
    return PreparedBiOEdgeOpticsBundle(
        BiOEdgeOpticsBundlePlan(map(component -> component.plan, plans)),
        plans, input, output)
end

function prepare_wfs_optics(front_end::BiOEdgeOpticalFrontEnd,
    inputs::Union{Tuple,AbstractVector}, output::OpticalProductBundle)
    return prepare_bi_o_edge_optical_bundle(front_end, inputs, output,
        front_end.source)
end

function prepare_bi_o_edge_optical_bundle(front_end::BiOEdgeOpticalFrontEnd,
    inputs::Union{Tuple,AbstractVector}, output::OpticalProductBundle,
    source::Union{Asterism,ExtendedSource})
    sources = four_pupil_path_sources(source)
    length(inputs) == length(sources) || throw(WFSPreparationError(
        :wfs_optics, :plane_count,
        "Bi-O-edge path-local pupil count does not match the source count"))
    length(output) == length(sources) || throw(WFSPreparationError(
        :wfs_optics, :plane_count,
        "Bi-O-edge path-local output count does not match the source count"))
    isempty(sources) && throw(WFSPreparationError(:wfs_optics,
        :plane_count, "Bi-O-edge path-local source collection is empty"))
    plans = ntuple(length(sources)) do index
        prepare_wfs_optics(
            bi_o_edge_front_end_with_source(front_end, sources[index]),
            inputs[index], output[index])
    end
    return PreparedBiOEdgeOpticsBundle(
        BiOEdgeOpticsBundlePlan(map(component -> component.plan, plans)),
        plans, inputs, output)
end

function prepare_bi_o_edge_optical_bundle(front_end, input, output, source)
    throw(WFSPreparationError(
        :wfs_optics, :plane_count,
        "Bi-O-edge product bundles require a spectral or path-expanded source"))
end

function _bi_o_edge_native_rate!(front_end::BiOEdgeOpticalFrontEnd,
    input::PupilFunction, lgs_model::AbstractPreparedFourPupilLGS)
    propagation = bi_o_edge_propagation_workspace(front_end)
    source = front_end.source
    resolution = size(input.opd, 1)
    pad = size(propagation.field, 1)
    offset = div(pad - resolution, 2)
    T = eltype(propagation.intensity)
    amplitude_scale = sqrt(T(photon_irradiance(source)) *
        T(input.metadata.sampling[1] * input.metadata.sampling[2]))
    opd_to_cycles = T(2) / T(wavelength(source))
    fill!(propagation.intensity, zero(T))
    @inbounds for point in 1:modulation_point_count(front_end.modulation)
        fill!(propagation.field, zero(eltype(propagation.field)))
        weight = front_end.modulation.amplitude_weights[point]
        @views @. propagation.field[offset+1:offset+resolution,
            offset+1:offset+resolution] = amplitude_scale * weight *
            input.amplitude * front_end.modulation.phases[:, :, point] *
            cispi(opd_to_cycles * input.opd)
        copyto!(propagation.focal_field, propagation.field)
        accumulate_bi_o_edge_masked_pupils!(propagation.intensity, front_end,
            lgs_model)
    end
    return propagation.intensity
end

function _bi_o_edge_native_rate!(front_end::BiOEdgeOpticalFrontEnd,
    input::ElectricField, lgs_model::AbstractPreparedFourPupilLGS)
    propagation = bi_o_edge_propagation_workspace(front_end)
    resolution = size(input.values, 1)
    pad = size(propagation.field, 1)
    offset = div(pad - resolution, 2)
    T = eltype(propagation.intensity)
    fill!(propagation.intensity, zero(T))
    @inbounds for point in 1:modulation_point_count(front_end.modulation)
        fill!(propagation.field, zero(eltype(propagation.field)))
        weight = front_end.modulation.amplitude_weights[point]
        @views @. propagation.field[offset+1:offset+resolution,
            offset+1:offset+resolution] = weight * input.values *
            front_end.modulation.phases[:, :, point]
        copyto!(propagation.focal_field, propagation.field)
        accumulate_bi_o_edge_masked_pupils!(propagation.intensity, front_end,
            lgs_model)
    end
    return propagation.intensity
end

function form_wfs_optical_products!(output::IntensityMap,
    input::Union{PupilFunction,ElectricField},
    plan::PreparedBiOEdgeOptics)
    validate_wfs_optics_binding(output, input, plan)
    native = _bi_o_edge_native_rate!(plan.front_end, input,
        plan.plan.lgs_model)
    factor = bi_o_edge_output_sampling_factor(plan.front_end,
        input.metadata.dimensions[1])
    bin2d!(output.values, native, factor)
    return output
end

function validate_wfs_optics_binding(output::IntensityMap,
    input::Union{PupilFunction,ElectricField},
    plan::PreparedBiOEdgeOptics)
    output === plan.output && input === plan.input ||
        throw(WFSPreparationError(:wfs_optics, :prepared_binding,
            "Bi-O-edge optical products do not match prepared storage"))
    workspace = bi_o_edge_propagation_workspace(plan.front_end)
    workspace === plan.workspace &&
        _bi_o_edge_propagation_workspace_binding(workspace) ===
            plan.workspace_binding &&
        workspace.revision == plan.plan.propagation_revision ||
        throw(WFSPreparationError(:wfs_optics, :prepared_binding,
            "Bi-O-edge propagation workspace changed after preparation"))
    plan.front_end.propagation.plan === plan.plan.propagation &&
        plan.front_end.amplitude_mask ===
            plan.plan.propagation.amplitude_mask &&
        plan.front_end.modulation === plan.plan.operating_modulation &&
        plan.front_end.calibration_modulation ===
            plan.plan.calibration_modulation &&
        plan.front_end.source === plan.plan.source ||
        throw(WFSPreparationError(:wfs_optics, :prepared_binding,
            "Bi-O-edge optics definition changed after preparation"))
    return nothing
end

function form_wfs_optical_products!(output::OpticalProductBundle, input,
    plan::PreparedBiOEdgeOpticsBundle)
    validate_wfs_optics_binding(output, input, plan)
    return form_four_pupil_bundle!(output, input, plan.components)
end

function validate_wfs_optics_binding(
    output::OpticalProductBundle, input,
    plan::PreparedBiOEdgeOpticsBundle)
    output === plan.output && input === plan.input ||
        throw(WFSPreparationError(:wfs_optics, :prepared_binding,
            "Bi-O-edge spectral products do not match prepared storage"))
    validate_four_pupil_bundle_binding(output, input, plan.components)
    return nothing
end

function bi_o_edge_rate_map(sensor::BiOEdgeWFS{<:Diffractive},
    inputs::Union{Tuple,AbstractVector}, source)
    return bi_o_edge_rate_map(BiOEdgeOpticalFrontEnd(sensor, source), inputs)
end

function bi_o_edge_rate_map(front_end::BiOEdgeOpticalFrontEnd,
    inputs::Union{Tuple,AbstractVector})
    return bi_o_edge_path_rate_bundle(front_end, inputs, front_end.source)
end

function bi_o_edge_path_rate_bundle(front_end::BiOEdgeOpticalFrontEnd,
    inputs::Union{Tuple,AbstractVector},
    source::Union{Asterism,ExtendedSource})
    sources = four_pupil_path_sources(source)
    length(inputs) == length(sources) || throw(WFSPreparationError(
        :wfs_optics, :plane_count,
        "Bi-O-edge path-local pupil count does not match the source count"))
    isempty(sources) && throw(WFSPreparationError(:wfs_optics,
        :plane_count, "Bi-O-edge path-local source collection is empty"))
    first_map = bi_o_edge_rate_map(
        bi_o_edge_front_end_with_source(front_end, sources[1]), inputs[1])
    maps = Vector{typeof(first_map)}(undef, length(sources))
    maps[1] = first_map
    @inbounds for index in 2:length(sources)
        maps[index] = bi_o_edge_rate_map(
            bi_o_edge_front_end_with_source(front_end, sources[index]),
            inputs[index])
    end
    return OpticalProductBundle(maps)
end

function bi_o_edge_path_rate_bundle(front_end, inputs, source)
    throw(WFSPreparationError(:wfs_optics, :plane_count,
        "path-local Bi-O-edge inputs require an Asterism or ExtendedSource"))
end

function bi_o_edge_rate_map(sensor::BiOEdgeWFS{<:Diffractive},
    input::Union{PupilFunction,ElectricField}, source=nothing)
    return bi_o_edge_rate_map(BiOEdgeOpticalFrontEnd(sensor, source), input)
end

function bi_o_edge_rate_map(front_end::BiOEdgeOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField})
    return _bi_o_edge_rate_map(front_end, input, front_end.source)
end

@inline _bi_o_edge_rate_map(front_end::BiOEdgeOpticalFrontEnd, input,
    source::SpectralSource) =
    _bi_o_edge_spectral_rate_bundle(front_end, input, source)

function _bi_o_edge_rate_map(front_end::BiOEdgeOpticalFrontEnd, input,
    source::Union{Asterism,ExtendedSource})
    throw(WFSPreparationError(:wfs_optics, :plane_count,
        "path-expanded Bi-O-edge sources require path-local pupil inputs"))
end

function _bi_o_edge_rate_map(front_end::BiOEdgeOpticalFrontEnd, input, source)
    wavelength_m = _bi_o_edge_front_end_wavelength(front_end, input)
    dimensions = bi_o_edge_output_dimensions(front_end,
        input.metadata.dimensions[1])
    T = eltype(bi_o_edge_propagation_workspace(front_end).intensity)
    values = similar(_modulated_input_storage(input), T, dimensions...)
    fill!(values, zero(T))
    factor = bi_o_edge_output_sampling_factor(front_end,
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

@inline _require_bi_o_edge_estimation_source(::WFSNormalization, source) =
    nothing

function _require_bi_o_edge_estimation_source(
    ::IncidenceFluxNormalization, ::Nothing)
    throw(WFSPreparationError(:estimation, :radiometry,
        "incidence-normalized Bi-O-edge estimation requires a source"))
end

function _bi_o_edge_spectral_rate_bundle(front_end::BiOEdgeOpticalFrontEnd,
    input, source::SpectralSource)
    samples = spectral_bundle(source).samples
    T = eltype(bi_o_edge_propagation_workspace(front_end).intensity)
    function component_map(sample)
        component = FourPupilSpectralComponent(source.source,
            T(sample.wavelength),
            T(photon_irradiance(source)) * T(sample.weight))
        component_front_end = bi_o_edge_front_end_with_source(front_end,
            component)
        return bi_o_edge_rate_map(component_front_end, input)
    end
    first_map = component_map(first(samples))
    maps = Vector{typeof(first_map)}(undef, length(samples))
    maps[1] = first_map
    @inbounds for index in 2:length(samples)
        maps[index] = component_map(samples[index])
    end
    return OpticalProductBundle(maps)
end

function _bi_o_edge_calibration_binding(sensor::BiOEdgeWFS)
    state = sensor.estimator.state
    state.calibrated || throw(WFSPreparationError(:estimation, :estimator,
        "Bi-O-edge estimation requires explicit calibration"))
    return BiOEdgeCalibrationBinding(state.calibration_revision,
        state.calibration_wavelength,
        state.calibration_signature, state.reference_signal_2d,
        state.valid_i4q)
end

function _require_bi_o_edge_calibration(sensor::BiOEdgeWFS,
    binding::BiOEdgeCalibrationBinding)
    state = sensor.estimator.state
    state.calibrated &&
        state.calibration_revision == binding.revision &&
        state.calibration_wavelength == binding.wavelength_m &&
        state.calibration_signature == binding.signature &&
        state.reference_signal_2d === binding.reference_signal &&
        state.valid_i4q === binding.valid_support ||
        throw(WFSPreparationError(:estimation, :prepared_binding,
            "Bi-O-edge calibration changed after estimator preparation"))
    return nothing
end

function _prepare_bi_o_edge_estimator_owner(sensor::BiOEdgeWFS, input,
    measurement::WFSMeasurement, path::AbstractWFSMeasurementPath,
    calibration_binding, source, normalization_scale)
    state = bi_o_edge_estimator_state(sensor)
    workspace = bi_o_edge_estimator_workspace(sensor)
    products = bi_o_edge_estimator_products(sensor)
    plan = BiOEdgeEstimationPlan(sensor.estimator.params, path,
        calibration_binding, source, normalization_scale)
    return PreparedBiOEdgeEstimator(plan, sensor, state, workspace, products,
        input, measurement, _bi_o_edge_estimator_state_binding(state),
        _bi_o_edge_estimator_workspace_binding(workspace),
        _bi_o_edge_estimator_products_binding(products),
        measurement.metadata.backend, measurement.metadata.device)
end

function _require_bi_o_edge_estimation_geometry(sensor::BiOEdgeWFS,
    frame_size::Int)
    iseven(frame_size) || throw(WFSPreparationError(:estimation, :shape,
        "Bi-O-edge observations require an even detector-frame size"))
    nominal = bi_o_edge_acquisition_workspace(sensor).nominal_detector_resolution
    binning = bi_o_edge_acquisition_plan(sensor).binning
    nominal > 0 || throw(WFSPreparationError(:estimation, :shape,
        "Bi-O-edge nominal detector resolution has not been prepared"))
    nominal % binning == 0 || throw(WFSPreparationError(:estimation, :shape,
        "Bi-O-edge binning does not divide the nominal detector resolution"))
    sampled_rows = binning == 1 ? 2 * nominal : div(nominal, binning)
    sampled_rows % frame_size == 0 || throw(WFSPreparationError(
        :estimation, :shape,
        "detector sampling does not evenly divide the Bi-O-edge frame"))
    return div(sampled_rows, frame_size)
end

function prepare_wfs_estimation(sensor::BiOEdgeWFS{<:Diffractive},
    observation::WFSObservation, measurement::WFSMeasurement;
    source=nothing, normalization_scale::Real=1)
    validate_wfs_observation(observation)
    validate_wfs_measurement(measurement)
    isequal(observation.metadata.layout, :four_pupil_mosaic) ||
        throw(WFSPreparationError(:estimation, :detector_mapping,
            "Bi-O-edge estimator requires :four_pupil_mosaic layout"))
    isequal(measurement.metadata.kind, :differential_slopes) ||
        throw(WFSPreparationError(:estimation, :estimator,
            "Bi-O-edge measurement kind must be :differential_slopes"))
    isequal(measurement.units, :dimensionless) ||
        throw(WFSPreparationError(:estimation, :units,
            "Bi-O-edge differential slopes are dimensionless"))
    frame_size = _require_real_square_wfs_observation(observation,
        "Bi-O-edge")
    measurement.metadata.numeric_type <: AbstractFloat ||
        throw(WFSPreparationError(:estimation, :numeric_type,
            "Bi-O-edge measurement storage must be floating point"))
    _require_wfs_storage_domain(:estimation, observation.metadata,
        bi_o_edge_estimator_workspace(sensor).signal_2d, "Bi-O-edge observation")
    _require_wfs_storage_domain(:estimation, measurement.metadata,
        bi_o_edge_estimator_products(sensor).slopes, "Bi-O-edge measurement")
    detector_reduction = _require_bi_o_edge_estimation_geometry(sensor,
        frame_size)
    resize_bi_o_edge_signal_buffers!(sensor, frame_size, detector_reduction)
    size(measurement.storage) == size(bi_o_edge_estimator_products(sensor).slopes) ||
        throw(WFSPreparationError(:estimation, :shape,
            "Bi-O-edge measurement storage has the wrong slope shape"))
    _require_bi_o_edge_estimation_source(
        sensor.estimator.params.normalization, source)
    scale = eltype(bi_o_edge_estimator_products(sensor).slopes)(normalization_scale)
    isfinite(scale) && scale >= zero(scale) || throw(WFSPreparationError(
        :estimation, :radiometry,
        "Bi-O-edge normalization scale must be finite and nonnegative"))
    binding = _bi_o_edge_calibration_binding(sensor)
    return _prepare_bi_o_edge_estimator_owner(sensor, observation,
        measurement, AcquiredObservationPath(), binding, source, scale)
end

function estimate_wfs_measurement!(measurement::WFSMeasurement,
    observation::WFSObservation,
    plan::PreparedBiOEdgeEstimator)
    validate_wfs_estimation_binding(measurement, observation, plan)
    sensor = plan.sensor
    _require_bi_o_edge_calibration(sensor, plan.plan.calibration_binding)
    bi_o_edge_signal!(execution_style(observation.storage), sensor,
        observation.storage, plan.plan.source, plan.plan.normalization_scale)
    copyto!(measurement.storage, bi_o_edge_estimator_products(sensor).slopes)
    return measurement
end

function validate_wfs_estimation_binding(measurement::WFSMeasurement, input,
    plan::PreparedBiOEdgeEstimator)
    measurement === plan.measurement && input === plan.input || throw(
        WFSPreparationError(:estimation, :prepared_binding,
            "Bi-O-edge estimator storage does not match its plan"))
    state = bi_o_edge_estimator_state(plan.sensor)
    workspace = bi_o_edge_estimator_workspace(plan.sensor)
    products = bi_o_edge_estimator_products(plan.sensor)
    state === plan.state && workspace === plan.workspace &&
        products === plan.products &&
        _bi_o_edge_estimator_state_binding(state) === plan.state_binding &&
        _bi_o_edge_estimator_workspace_binding(workspace) ===
            plan.workspace_binding &&
        _bi_o_edge_estimator_products_binding(products) ===
            plan.products_binding || throw(WFSPreparationError(
                :estimation, :prepared_binding,
                "Bi-O-edge estimator state, workspace, or products changed after preparation"))
    return nothing
end

function prepare_wfs_estimation(sensor::BiOEdgeWFS{<:Geometric},
    input::PupilFunction, measurement::WFSMeasurement)
    require_modulated_wfs_input(input)
    validate_wfs_measurement(measurement)
    input.metadata.dimensions == (sensor.estimator.params.pupil_resolution,
        sensor.estimator.params.pupil_resolution) ||
        throw(WFSPreparationError(:estimation, :shape,
            "geometric Bi-O-edge input has the wrong pupil dimensions"))
    isequal(measurement.metadata.kind, :geometric_slopes) ||
        throw(WFSPreparationError(:estimation, :estimator,
            "geometric Bi-O-edge measurement kind must be :geometric_slopes"))
    isequal(measurement.units, :metre) || throw(WFSPreparationError(
        :estimation, :units,
        "geometric Bi-O-edge OPD differences are expressed in metres"))
    measurement.metadata.numeric_type <: AbstractFloat ||
        throw(WFSPreparationError(:estimation, :numeric_type,
            "geometric Bi-O-edge measurement storage must be floating point"))
    _require_wfs_storage_domain(:estimation, input.metadata,
        bi_o_edge_estimator_products(sensor).slopes, "geometric Bi-O-edge input")
    _require_wfs_storage_domain(:estimation, measurement.metadata,
        bi_o_edge_estimator_products(sensor).slopes, "geometric Bi-O-edge measurement")
    size(measurement.storage) == size(bi_o_edge_estimator_products(sensor).slopes) ||
        throw(WFSPreparationError(:estimation, :shape,
            "geometric Bi-O-edge measurement has the wrong slope shape"))
    return _prepare_bi_o_edge_estimator_owner(sensor, input, measurement,
        DirectMeasurementPath(), nothing, nothing,
        one(eltype(bi_o_edge_estimator_products(sensor).slopes)))
end

function estimate_wfs_measurement!(measurement::WFSMeasurement,
    input::PupilFunction,
    plan::PreparedBiOEdgeEstimator)
    validate_wfs_estimation_binding(measurement, input, plan)
    sensor = plan.sensor
    state = sensor.estimator.state
    products = bi_o_edge_estimator_products(sensor)
    edge_geometric_slopes!(products.slopes, input.opd, state.valid_mask,
        state.edge_mask)
    @. products.slopes *= state.optical_gain
    copyto!(measurement.storage, products.slopes)
    return measurement
end

function set_bi_o_edge_calibration!(sensor::BiOEdgeWFS,
    reference::AbstractMatrix; wavelength_m::Real,
    signature::UInt=UInt(0), valid_support=nothing)
    state = sensor.estimator.state
    size(reference) == size(state.reference_signal_2d) ||
        throw(DimensionMismatchError(
            "Bi-O-edge reference dimensions do not match estimator storage"))
    require_same_backend(state.reference_signal_2d, reference)
    reference_host = Array(reference)
    all(isfinite, reference_host) || throw(InvalidConfiguration(
        "Bi-O-edge calibration reference must contain only finite values"))
    wavelength_value = eltype(bi_o_edge_estimator_products(sensor).slopes)(
        wavelength_m)
    isfinite(wavelength_value) && wavelength_value > zero(wavelength_value) ||
        throw(InvalidConfiguration(
            "Bi-O-edge calibration wavelength must be finite and positive"))
    support_host = _prepare_bi_o_edge_calibration_support(sensor, valid_support)
    copyto!(state.reference_signal_2d, reference_host)
    if support_host === nothing
        fill!(state.valid_i4q, true)
    else
        copyto!(state.valid_i4q, support_host)
    end
    update_bi_o_edge_valid_signal!(sensor)
    update_bi_o_edge_valid_signal_indices!(sensor)
    resize_bi_o_edge_slope_buffers!(sensor)
    state.calibration_wavelength = wavelength_value
    state.calibration_signature = signature
    state.calibrated = true
    state.calibration_revision += UInt(1)
    return sensor
end

function _prepare_bi_o_edge_calibration_support(sensor::BiOEdgeWFS, ::Nothing)
    if !iszero(sensor.estimator.params.light_ratio)
        throw(InvalidConfiguration(
            "nonzero Bi-O-edge light_ratio requires explicit valid_support"))
    end
    return nothing
end

function _prepare_bi_o_edge_calibration_support(sensor::BiOEdgeWFS,
    valid_support::AbstractMatrix{Bool})
    state = sensor.estimator.state
    size(valid_support) == size(state.valid_i4q) ||
        throw(DimensionMismatchError(
            "Bi-O-edge calibration support has the wrong dimensions"))
    require_same_backend(state.valid_i4q, valid_support)
    support_host = Array(valid_support)
    any(support_host) || throw(InvalidConfiguration(
        "Bi-O-edge calibration support must select at least one sample"))
    return support_host
end

function _prepare_bi_o_edge_calibration_support(sensor::BiOEdgeWFS,
    valid_support)
    throw(InvalidConfiguration(
        "Bi-O-edge calibration support must be a Boolean matrix"))
end

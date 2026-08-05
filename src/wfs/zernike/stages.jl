#
# Prepared Zernike WFS stages
#

"""
Run-immutable physical and numerical contract for one Zernike detector-plane
photon-arrival-rate map.
"""
struct ZernikeOpticsPlan{P,S} <: AbstractWFSOpticsPlan
    propagation::P
    source::S
    propagation_revision::UInt
end

"""Exact live owner for one prepared Zernike optics execution."""
struct PreparedZernikeOptics{P,F,W,I,O,R,B,D}
    plan::P
    front_end::F
    workspace::W
    input::I
    output::O
    workspace_binding::R
    backend::B
    device::D
end

@inline wfs_optical_products(prepared::PreparedZernikeOptics) =
    prepared.output

struct ZernikeCalibrationBinding{T<:AbstractFloat,R,A}
    revision::UInt
    wavelength_m::T
    signature::UInt
    reference_signal::R
    valid_support::A
end

"""Run-immutable normalized-pupil estimation contract."""
struct ZernikeEstimationPlan{E,P<:AbstractWFSMeasurementPath,C,S,T} <:
        AbstractWFSEstimationPlan
    params::E
    path::P
    calibration_binding::C
    source::S
    normalization_scale::T
end

"""Exact live owner for one prepared Zernike estimator."""
struct PreparedZernikeEstimator{P,W,ST,WS,PR,I,M,SB,WB,PB,B,D}
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

@inline wfs_measurement_path(prepared::PreparedZernikeEstimator) =
    prepared.plan.path

@inline function _zernike_propagation_workspace_binding(workspace)
    return (workspace.field, workspace.focal_field, workspace.pupil_field,
        workspace.phasor, workspace.phase_mask, workspace.pupil_intensity,
        workspace.nominal_frame, workspace.fft_plan, workspace.ifft_plan)
end

@inline function _zernike_estimator_state_binding(state)
    return (state.valid_mask, state.reference_signal_2d)
end

@inline function _zernike_estimator_workspace_binding(workspace)
    return (workspace.valid_signal_indices, workspace.signal_2d,
        workspace.normalization_frame, workspace.normalization_partials,
        workspace.normalization_sum, workspace.normalization_sum_host)
end

@inline _zernike_estimator_products_binding(products) = (products.signal,)

@inline _zernike_input_storages(input::PupilFunction) =
    (input.amplitude, input.opd, input.support)
@inline _zernike_input_storages(input::ElectricField) = (input.values,)

@inline _zernike_mightalias_any(value, ::Tuple{}) = false
@inline function _zernike_mightalias_any(value, values::Tuple)
    return _wfs_storage_mightalias(value, first(values)) ||
        _zernike_mightalias_any(value, Base.tail(values))
end

@inline _zernike_any_alias(::Tuple{}) = false
@inline function _zernike_any_alias(values::Tuple)
    remaining = Base.tail(values)
    return _zernike_mightalias_any(first(values), remaining) ||
        _zernike_any_alias(remaining)
end

function _require_zernike_optics_aliases(input, output, workspace)
    storages = (_zernike_input_storages(input)..., output.values,
        _zernike_propagation_workspace_binding(workspace)...)
    _zernike_any_alias(storages) && throw(WFSPreparationError(
        :wfs_optics, :aliasing,
        "Zernike input, rate product, and propagation workspace must not alias"))
    return nothing
end

function _require_zernike_estimation_aliases(sensor::ZernikeWFS,
    observation::WFSObservation, measurement::WFSMeasurement)
    state = zernike_estimator_state(sensor)
    workspace = zernike_estimator_workspace(sensor)
    products = zernike_estimator_products(sensor)
    storages = (observation.storage, measurement.storage,
        _zernike_estimator_state_binding(state)...,
        _zernike_estimator_workspace_binding(workspace)...,
        _zernike_estimator_products_binding(products)...)
    _zernike_any_alias(storages) && throw(WFSPreparationError(
        :estimation, :aliasing,
        "Zernike observation, measurement, calibration, workspace, and product storage must not alias"))
    return nothing
end

@inline modulated_wfs_propagation_storage(
    front_end::ZernikeOpticalFrontEnd) =
    zernike_propagation_workspace(front_end).field

function ZernikeOpticalFrontEnd(sensor::ZernikeWFS, source=nothing)
    front_end = sensor.front_end
    return ZernikeOpticalFrontEnd(front_end.phase_spot,
        front_end.propagation, front_end.binning, source)
end

@inline function zernike_rate_dimensions(front_end::ZernikeOpticalFrontEnd)
    pupil_samples = zernike_propagation_plan(
        front_end.propagation).pupil_samples
    return (div(pupil_samples, front_end.binning),
        div(pupil_samples, front_end.binning))
end

function _require_zernike_front_end_source(
    front_end::ZernikeOpticalFrontEnd, ::PupilFunction)
    source = front_end.source
    source === nothing && throw(WFSPreparationError(:wfs_optics,
        :radiometry, "Zernike WFS optics require a source for PupilFunction input"))
    require_leaf_source(source, "prepared Zernike optics")
    return source
end

function _require_zernike_front_end_source(
    front_end::ZernikeOpticalFrontEnd, ::ElectricField)
    front_end.source === nothing || throw(WFSPreparationError(
        :wfs_optics, :radiometry,
        "photon-rate ElectricField input must not also supply a Zernike source"))
    return nothing
end

@inline _zernike_front_end_wavelength(front_end::ZernikeOpticalFrontEnd,
    input::PupilFunction) = modulated_input_wavelength(input,
        front_end.source)
@inline _zernike_front_end_wavelength(::ZernikeOpticalFrontEnd,
    input::ElectricField) = modulated_input_wavelength(input)

@inline _require_zernike_rate_coordinates(
    ::NormalizedPupilCoordinates) = nothing

function _require_zernike_rate_coordinates(::AbstractPlaneCoordinateDomain)
    throw(WFSPreparationError(:wfs_optics, :plane_metadata,
        "Zernike detector output must use normalized pupil coordinates"))
end

@inline _require_zernike_rate_measure(::CellIntegratedMeasure) = nothing

function _require_zernike_rate_measure(::AbstractSpatialMeasure)
    throw(WFSPreparationError(:wfs_optics, :radiometry,
        "Zernike detector output must carry cell-integrated rate"))
end

function _require_zernike_rate_wavelength(channel::MonochromaticChannel,
    wavelength_m)
    channel.wavelength_m == wavelength_m || throw(
        WFSPreparationError(:wfs_optics, :plane_metadata,
            "Zernike detector output wavelength differs from its input"))
    return nothing
end

function _require_zernike_rate_wavelength(
    ::AbstractSpectralCoordinate, ::Any)
    throw(WFSPreparationError(:wfs_optics, :plane_metadata,
        "Zernike detector output wavelength differs from its input"))
end

function _require_zernike_rate_map(output::IntensityMap,
    expected_dimensions, wavelength_m)
    validate_wfs_optical_products(output)
    _require_zernike_rate_coordinates(output.metadata.coordinate_domain)
    _require_zernike_rate_measure(output.metadata.spatial_measure)
    size(output.values) == expected_dimensions || throw(
        WFSPreparationError(:wfs_optics, :shape,
            "Zernike detector output has the wrong prepared dimensions"))
    _require_zernike_rate_wavelength(output.metadata.spectral, wavelength_m)
    return output
end

function _require_zernike_input_geometry(front_end::ZernikeOpticalFrontEnd,
    input::PupilFunction)
    pupil_resolution = zernike_propagation_plan(
        front_end.propagation).pupil_resolution
    input.metadata.dimensions == (pupil_resolution,
        pupil_resolution) || throw(WFSPreparationError(
        :wfs_optics, :shape,
        "Zernike pupil input dimensions differ from the prepared relay"))
    return nothing
end

function _require_zernike_input_geometry(front_end::ZernikeOpticalFrontEnd,
    input::ElectricField)
    workspace = zernike_propagation_workspace(front_end)
    input.metadata.dimensions == size(workspace.field) || throw(
        WFSPreparationError(:wfs_optics, :shape,
            "Zernike ElectricField dimensions differ from the prepared diffraction grid"))
    return nothing
end

function prepare_wfs_optics(front_end::ZernikeOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField}, output::IntensityMap)
    require_modulated_wfs_input(input)
    _require_zernike_front_end_source(front_end, input)
    _require_zernike_input_geometry(front_end, input)
    wavelength_m = _zernike_front_end_wavelength(front_end, input)
    _require_zernike_rate_map(output, zernike_rate_dimensions(front_end),
        wavelength_m)
    require_modulated_wfs_domains(front_end, input, output)
    workspace = zernike_propagation_workspace(front_end)
    eltype(workspace.pupil_intensity) ===
        output.metadata.numeric_type || throw(WFSPreparationError(
            :wfs_optics, :numeric_type,
            "Zernike output precision differs from prepared propagation"))
    _require_zernike_optics_aliases(input, output, workspace)
    propagation = front_end.propagation
    propagation_plan = zernike_propagation_plan(propagation)
    plan = ZernikeOpticsPlan(propagation_plan, front_end.source,
        workspace.revision)
    return PreparedZernikeOptics(plan, front_end, workspace, input, output,
        _zernike_propagation_workspace_binding(workspace),
        input.metadata.backend, input.metadata.device)
end

function zernike_rate_map(sensor::ZernikeWFS,
    input::Union{PupilFunction,ElectricField}, source=nothing)
    return zernike_rate_map(ZernikeOpticalFrontEnd(sensor, source), input)
end

function zernike_rate_map(front_end::ZernikeOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField})
    wavelength_m = _zernike_front_end_wavelength(front_end, input)
    dimensions = zernike_rate_dimensions(front_end)
    propagation_plan = zernike_propagation_plan(front_end.propagation)
    T = eltype(zernike_propagation_workspace(front_end).pupil_intensity)
    values = similar(_modulated_input_storage(input), T, dimensions...)
    fill!(values, zero(T))
    normalized_sampling = T(front_end.binning / propagation_plan.pupil_samples)
    metadata = OpticalPlaneMetadata(DetectorPlane(), values;
        coordinate_domain=NormalizedPupilCoordinates(),
        sampling=(normalized_sampling, normalized_sampling),
        spectral=MonochromaticChannel(T(wavelength_m)),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    return IntensityMap(metadata, values)
end

function _require_zernike_optical_binding(
    prepared::PreparedZernikeOptics, input, output)
    input === prepared.input && output === prepared.output || throw(
        WFSPreparationError(:wfs_optics, :prepared_binding,
            "Zernike optical products do not match their prepared plan"))
    workspace = prepared.workspace
    prepared.backend === input.metadata.backend &&
        prepared.device == input.metadata.device || throw(
        WFSPreparationError(:wfs_optics, :prepared_binding,
            "Zernike optical input target changed after preparation"))
    prepared.plan.propagation ===
        zernike_propagation_plan(prepared.front_end.propagation) &&
        prepared.front_end.phase_spot ===
            prepared.plan.propagation.phase_spot &&
        prepared.front_end.source === prepared.plan.source &&
        workspace === zernike_propagation_workspace(prepared.front_end) &&
        _zernike_propagation_workspace_binding(workspace) ===
            prepared.workspace_binding &&
        workspace.revision == prepared.plan.propagation_revision || throw(
        WFSPreparationError(:wfs_optics, :prepared_binding,
            "Zernike propagation storage changed after preparation"))
    return nothing
end

@inline validate_wfs_optics_binding(output::IntensityMap, input,
    plan::PreparedZernikeOptics) =
    _require_zernike_optical_binding(plan, input, output)

function _form_zernike_input_field!(front_end::ZernikeOpticalFrontEnd,
    input::PupilFunction)
    propagation = zernike_propagation_workspace(front_end)
    propagation_plan = zernike_propagation_plan(front_end.propagation)
    T = eltype(propagation.pupil_intensity)
    n = propagation_plan.pupil_resolution
    pad = size(propagation.field, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    cell_area = T(input.metadata.sampling[1] * input.metadata.sampling[2])
    amplitude_scale = sqrt(T(photon_irradiance(front_end.source)) *
        cell_area)
    opd_to_cycles = T(2) / T(wavelength(front_end.source))
    fill!(propagation.field, zero(eltype(propagation.field)))
    @views @. propagation.field[ox+1:ox+n, oy+1:oy+n] =
        amplitude_scale * input.amplitude * cispi(opd_to_cycles * input.opd)
    return propagation.field
end

function _form_zernike_input_field!(front_end::ZernikeOpticalFrontEnd,
    input::ElectricField)
    workspace = zernike_propagation_workspace(front_end)
    copyto!(workspace.field, input.values)
    return workspace.field
end

function _form_zernike_rate!(output::AbstractMatrix,
    front_end::ZernikeOpticalFrontEnd, input)
    propagation = zernike_propagation_workspace(front_end)
    propagation_plan = zernike_propagation_plan(front_end.propagation)
    _form_zernike_input_field!(front_end, input)
    copyto!(propagation.focal_field, propagation.field)
    @. propagation.focal_field *= propagation.phasor
    execute_fft_plan!(propagation.focal_field, propagation.fft_plan)
    @. propagation.focal_field *= propagation.phase_mask
    copyto!(propagation.pupil_field, propagation.focal_field)
    execute_fft_plan!(propagation.pupil_field, propagation.ifft_plan)
    n = propagation_plan.pupil_resolution
    pad = size(propagation.pupil_field, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    @views @. propagation.pupil_intensity =
        abs2(propagation.pupil_field[ox+1:ox+n, oy+1:oy+n])
    sampling = div(n, propagation_plan.pupil_samples)
    bin2d!(propagation.nominal_frame, propagation.pupil_intensity, sampling)
    if front_end.binning == 1
        copyto!(output, propagation.nominal_frame)
    else
        bin2d!(output, propagation.nominal_frame, front_end.binning)
    end
    return output
end

function form_wfs_optical_products!(output::IntensityMap,
    input::Union{PupilFunction,ElectricField},
    plan::PreparedZernikeOptics)
    validate_wfs_optics_binding(output, input, plan)
    _form_zernike_rate!(output.values, plan.front_end, input)
    return output
end

function _zernike_calibration_binding(sensor::ZernikeWFS)
    state = zernike_estimator_state(sensor)
    state.calibrated || throw(WFSPreparationError(:estimation, :estimator,
        "Zernike estimation requires explicit calibration"))
    return ZernikeCalibrationBinding(state.calibration_revision,
        state.calibration_wavelength, state.calibration_signature,
        state.reference_signal_2d, state.valid_mask)
end

function _require_zernike_calibration(estimator,
    binding::ZernikeCalibrationBinding)
    state = zernike_estimator_state(estimator)
    state.calibrated &&
        state.calibration_revision == binding.revision &&
        state.calibration_wavelength == binding.wavelength_m &&
        state.calibration_signature == binding.signature &&
        state.reference_signal_2d === binding.reference_signal &&
        state.valid_mask === binding.valid_support || throw(
        WFSPreparationError(:estimation, :prepared_binding,
            "Zernike calibration changed after estimator preparation"))
    return nothing
end

function _prepare_zernike_estimator_owner(sensor::ZernikeWFS, input,
    measurement::WFSMeasurement, path::AbstractWFSMeasurementPath,
    calibration_binding, source, normalization_scale)
    state = zernike_estimator_state(sensor)
    workspace = zernike_estimator_workspace(sensor)
    products = zernike_estimator_products(sensor)
    plan = ZernikeEstimationPlan(zernike_estimator_params(sensor), path,
        calibration_binding, source, normalization_scale)
    return PreparedZernikeEstimator(plan, sensor, state, workspace, products,
        input, measurement, _zernike_estimator_state_binding(state),
        _zernike_estimator_workspace_binding(workspace),
        _zernike_estimator_products_binding(products),
        measurement.metadata.backend, measurement.metadata.device)
end

@inline _require_zernike_estimation_source(
    ::MeanValidFluxNormalization, source) = nothing

function _require_zernike_estimation_source(
    ::IncidenceFluxNormalization, ::Nothing)
    throw(WFSPreparationError(:estimation, :radiometry,
        "Zernike incidence normalization requires its source"))
end

function _require_zernike_estimation_source(
    ::IncidenceFluxNormalization, source::AbstractSource)
    require_leaf_source(source, "Zernike incidence normalization")
    _require_physical_photon_irradiance(source,
        "Zernike incidence normalization")
    return source
end

function prepare_wfs_estimation(sensor::ZernikeWFS,
    observation::WFSObservation, measurement::WFSMeasurement;
    source=nothing, normalization_scale::Real=1)
    validate_wfs_observation(observation)
    validate_wfs_measurement(measurement)
    isequal(observation.metadata.layout, :zernike_pupil_image) || throw(
        WFSPreparationError(:estimation, :detector_mapping,
            "Zernike estimator requires :zernike_pupil_image layout"))
    isequal(measurement.metadata.kind, :normalized_pupil_signal) || throw(
        WFSPreparationError(:estimation, :estimator,
            "Zernike measurement kind must be :normalized_pupil_signal"))
    isequal(measurement.units, :dimensionless) || throw(
        WFSPreparationError(:estimation, :units,
            "Zernike pupil signals are dimensionless"))
    observation.metadata.numeric_type <: Real || throw(
        WFSPreparationError(:estimation, :numeric_type,
            "Zernike observations require real detector samples"))
    measurement.metadata.numeric_type <: AbstractFloat || throw(
        WFSPreparationError(:estimation, :numeric_type,
            "Zernike measurement storage must be floating point"))
    state = zernike_estimator_state(sensor)
    workspace = zernike_estimator_workspace(sensor)
    products = zernike_estimator_products(sensor)
    observation.metadata.dimensions == size(workspace.signal_2d) || throw(
        WFSPreparationError(:estimation, :shape,
            "Zernike observation has the wrong pupil-image dimensions"))
    size(measurement.storage) == size(products.signal) || throw(
        WFSPreparationError(:estimation, :shape,
            "Zernike measurement has the wrong signal-vector dimensions"))
    _require_wfs_storage_domain(:estimation, observation.metadata,
        workspace.signal_2d, "Zernike observation")
    _require_wfs_storage_domain(:estimation, measurement.metadata,
        products.signal, "Zernike measurement")
    measurement.metadata.numeric_type === eltype(products.signal) || throw(
        WFSPreparationError(:estimation, :numeric_type,
            "Zernike measurement precision differs from its estimator"))
    _require_zernike_estimation_aliases(sensor, observation, measurement)
    params = zernike_estimator_params(sensor)
    _require_zernike_estimation_source(params.normalization, source)
    T = eltype(products.signal)
    scale = T(normalization_scale)
    isfinite(scale) && scale >= zero(T) || throw(WFSPreparationError(
        :estimation, :radiometry,
        "Zernike normalization scale must be finite and nonnegative"))
    binding = _zernike_calibration_binding(sensor)
    return _prepare_zernike_estimator_owner(sensor, observation, measurement,
        AcquiredObservationPath(), binding, source, scale)
end

@inline function _zernike_incidence_multiplier(sensor::ZernikeWFS,
    source::AbstractSource, normalization_scale)
    params = zernike_estimator_params(sensor)
    T = eltype(zernike_estimator_products(sensor).signal)
    pupil_sample = T(params.pupil_diameter_m) / T(params.pupil_resolution)
    irradiance = T(_require_physical_photon_irradiance(source,
        "Zernike incidence normalization"))
    return irradiance * abs2(pupil_sample) * normalization_scale /
        T(zernike_normalization_count(sensor))
end

function _zernike_scalar_normalization(::MeanValidFluxNormalization,
    sensor::ZernikeWFS, frame::AbstractMatrix, source,
    normalization_scale::S) where {S<:AbstractFloat}
    state = zernike_estimator_state(sensor)
    summed = zero(S)
    @inbounds for j in axes(frame, 2), i in axes(frame, 1)
        state.valid_mask[i, j] && (summed += S(frame[i, j]))
    end
    return max(summed / S(zernike_normalization_count(sensor)), eps(S))
end

function _zernike_scalar_normalization(::IncidenceFluxNormalization,
    sensor::ZernikeWFS, frame::AbstractMatrix, source,
    normalization_scale::S) where {S<:AbstractFloat}
    state = zernike_estimator_state(sensor)
    workspace = zernike_estimator_workspace(sensor)
    summed = zero(S)
    @inbounds for j in axes(workspace.normalization_frame, 2),
            i in axes(workspace.normalization_frame, 1)
        state.valid_mask[i, j] &&
            (summed += workspace.normalization_frame[i, j])
    end
    return summed * _zernike_incidence_multiplier(sensor, source,
        normalization_scale)
end

function _estimate_zernike_signal!(::ScalarCPUStyle, sensor::ZernikeWFS,
    frame::AbstractMatrix{F}, source, normalization_scale::S) where {
    F<:Real,S<:AbstractFloat,
}
    params = zernike_estimator_params(sensor)
    state = zernike_estimator_state(sensor)
    workspace = zernike_estimator_workspace(sensor)
    products = zernike_estimator_products(sensor)
    count = zernike_normalization_count(sensor)
    fill!(workspace.signal_2d, zero(S))
    count == 0 && (fill!(products.signal, zero(S)); return products.signal)
    normalization = _zernike_scalar_normalization(
        params.normalization, sensor, frame, source,
        normalization_scale)
    usable = isfinite(normalization) && normalization > eps(S)
    if !usable
        fill!(products.signal, zero(S))
        return products.signal
    end
    @inbounds for j in axes(frame, 2), i in axes(frame, 1)
        if state.valid_mask[i, j]
            workspace.signal_2d[i, j] = S(frame[i, j]) / normalization -
                state.reference_signal_2d[i, j]
        end
    end
    @inbounds for index in eachindex(workspace.valid_signal_indices)
        products.signal[index] =
            workspace.signal_2d[workspace.valid_signal_indices[index]]
    end
    return products.signal
end

function _queue_zernike_stage_normalization!(phase::KernelLaunchPhase,
    normalization::MeanValidFluxNormalization, sensor::ZernikeWFS,
    frame, source, normalization_scale)
    queue_zernike_masked_sum!(phase, sensor, frame)
    T = eltype(zernike_estimator_products(sensor).signal)
    return inv(T(zernike_normalization_count(sensor))), true
end

function _queue_zernike_stage_normalization!(phase::KernelLaunchPhase,
    normalization::IncidenceFluxNormalization, sensor::ZernikeWFS,
    frame, source, normalization_scale)
    workspace = zernike_estimator_workspace(sensor)
    queue_zernike_masked_sum!(phase, sensor,
        workspace.normalization_frame)
    return _zernike_incidence_multiplier(sensor, source,
        normalization_scale), false
end

function _estimate_zernike_signal!(style::AcceleratorStyle,
    sensor::ZernikeWFS, frame::AbstractMatrix{F}, source,
    normalization_scale::S) where {F<:Real,S<:AbstractFloat}
    params = zernike_estimator_params(sensor)
    state = zernike_estimator_state(sensor)
    workspace = zernike_estimator_workspace(sensor)
    products = zernike_estimator_products(sensor)
    count = zernike_normalization_count(sensor)
    if count == 0
        fill!(workspace.signal_2d, zero(S))
        fill!(products.signal, zero(S))
        return products.signal
    end
    phase = begin_kernel_phase(style)
    multiplier, clamp_to_epsilon = _queue_zernike_stage_normalization!(
        phase, params.normalization, sensor, frame, source,
        normalization_scale)
    queue_kernel!(phase, zernike_signal_kernel!, workspace.signal_2d, frame,
        state.valid_mask, state.reference_signal_2d,
        workspace.normalization_sum, multiplier, clamp_to_epsilon,
        size(frame, 1), size(frame, 2); ndrange=size(frame))
    queue_kernel!(phase, gather_zernike_signal_kernel!, products.signal,
        workspace.signal_2d, workspace.valid_signal_indices, count;
        ndrange=count)
    finish_kernel_phase!(phase)
    return products.signal
end

function estimate_wfs_measurement!(measurement::WFSMeasurement,
    observation::WFSObservation, plan::PreparedZernikeEstimator)
    measurement === plan.measurement && observation === plan.input || throw(
        WFSPreparationError(:estimation, :prepared_binding,
            "Zernike estimator storage does not match its plan"))
    validate_wfs_estimation_binding(measurement, observation, plan)
    sensor = plan.sensor
    _require_zernike_calibration(sensor, plan.plan.calibration_binding)
    _estimate_zernike_signal!(execution_style(observation.storage), sensor,
        observation.storage, plan.plan.source, plan.plan.normalization_scale)
    copyto!(measurement.storage, zernike_estimator_products(sensor).signal)
    return measurement
end

function validate_wfs_estimation_binding(measurement::WFSMeasurement, input,
    plan::PreparedZernikeEstimator)
    measurement === plan.measurement && input === plan.input || throw(
        WFSPreparationError(:estimation, :prepared_binding,
            "Zernike estimator storage does not match its plan"))
    state = zernike_estimator_state(plan.sensor)
    workspace = zernike_estimator_workspace(plan.sensor)
    products = zernike_estimator_products(plan.sensor)
    state === plan.state && workspace === plan.workspace &&
        products === plan.products &&
        _zernike_estimator_state_binding(state) === plan.state_binding &&
        _zernike_estimator_workspace_binding(workspace) ===
            plan.workspace_binding &&
        _zernike_estimator_products_binding(products) ===
            plan.products_binding &&
        plan.backend === measurement.metadata.backend &&
        plan.device == measurement.metadata.device || throw(
        WFSPreparationError(:estimation, :prepared_binding,
            "Zernike estimator state, workspace, products, or target changed after preparation"))
    return nothing
end

function set_zernike_calibration!(sensor::ZernikeWFS,
    reference::AbstractMatrix; wavelength_m::Real,
    signature::UInt=UInt(0))
    state = zernike_estimator_state(sensor)
    workspace = zernike_estimator_workspace(sensor)
    products = zernike_estimator_products(sensor)
    size(reference) == size(state.reference_signal_2d) || throw(
        InvalidConfiguration(
            "Zernike calibration reference has the wrong dimensions"))
    eltype(reference) <: Real || throw(InvalidConfiguration(
        "Zernike calibration reference must contain real values"))
    typeof(backend(reference)) === typeof(backend(state.reference_signal_2d)) ||
        throw(InvalidConfiguration(
            "Zernike calibration reference backend differs from the estimator"))
    compute_device(reference) == compute_device(state.reference_signal_2d) ||
        throw(InvalidConfiguration(
            "Zernike calibration reference occupies another device"))
    all(isfinite, host_array(reference)) || throw(InvalidConfiguration(
        "Zernike calibration reference must contain finite values"))
    T = eltype(state.reference_signal_2d)
    wavelength_value = T(wavelength_m)
    isfinite(wavelength_value) && wavelength_value > zero(T) || throw(
        InvalidConfiguration(
            "Zernike calibration wavelength must be finite and positive"))
    copyto!(state.reference_signal_2d, reference)
    fill!(workspace.signal_2d, zero(T))
    fill!(products.signal, zero(T))
    state.calibrated = true
    state.calibration_wavelength = wavelength_value
    state.calibration_signature = signature
    state.calibration_revision += UInt(1)
    return sensor
end

#
# Prepared Zernike WFS stages
#

struct ZernikePropagationBinding{F,Q,U,H,M,N,Pf,Pi}
    field::F
    focal_field::Q
    pupil_field::U
    phasor::H
    phase_mask::M
    nominal_frame::N
    fft_plan::Pf
    ifft_plan::Pi
    revision::UInt
end

struct PreparedZernikeOpticalFormation{F,I,O,B}
    front_end::F
    input::I
    output::O
    binding::B
end

struct ZernikeCalibrationBinding{T<:AbstractFloat,R,A}
    revision::UInt
    wavelength_m::T
    signature::UInt
    reference_signal::R
    valid_support::A
end

struct PreparedZernikeEstimator{E,I,M,P<:AbstractWFSMeasurementPath,C,S,T}
    estimator::E
    input::I
    measurement::M
    path::P
    calibration_binding::C
    source::S
    normalization_scale::T
end

@inline wfs_measurement_path(plan::PreparedZernikeEstimator) = plan.path

function ZernikeOpticalFrontEnd(sensor::ZernikeWFS, source=nothing)
    front_end = sensor.front_end
    return ZernikeOpticalFrontEnd(front_end.phase_spot,
        front_end.propagation, front_end.pupil_resolution,
        front_end.pupil_diameter_m, front_end.pupil_samples,
        front_end.binning, source)
end

@inline zernike_rate_dimensions(front_end::ZernikeOpticalFrontEnd) =
    (div(front_end.pupil_samples, front_end.binning),
        div(front_end.pupil_samples, front_end.binning))

function _require_zernike_front_end_source(
    front_end::ZernikeOpticalFrontEnd, ::PupilFunction)
    source = front_end.source
    source === nothing && throw(WFSPreparationError(:optical_formation,
        :radiometry, "Zernike PupilFunction formation requires a source"))
    require_leaf_source(source, "prepared Zernike optical formation")
    return source
end

function _require_zernike_front_end_source(
    front_end::ZernikeOpticalFrontEnd, ::ElectricField)
    front_end.source === nothing || throw(WFSPreparationError(
        :optical_formation, :radiometry,
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
    throw(WFSPreparationError(:optical_formation, :plane_metadata,
        "Zernike detector output must use normalized pupil coordinates"))
end

@inline _require_zernike_rate_measure(::CellIntegratedMeasure) = nothing

function _require_zernike_rate_measure(::AbstractSpatialMeasure)
    throw(WFSPreparationError(:optical_formation, :radiometry,
        "Zernike detector output must carry cell-integrated rate"))
end

function _require_zernike_rate_wavelength(channel::MonochromaticChannel,
    wavelength_m)
    channel.wavelength_m == wavelength_m || throw(
        WFSPreparationError(:optical_formation, :plane_metadata,
            "Zernike detector output wavelength differs from its input"))
    return nothing
end

function _require_zernike_rate_wavelength(
    ::AbstractSpectralCoordinate, ::Any)
    throw(WFSPreparationError(:optical_formation, :plane_metadata,
        "Zernike detector output wavelength differs from its input"))
end

function _require_zernike_rate_map(output::IntensityMap,
    expected_dimensions, wavelength_m)
    validate_wfs_optical_products(output)
    _require_zernike_rate_coordinates(output.metadata.coordinate_domain)
    _require_zernike_rate_measure(output.metadata.spatial_measure)
    size(output.values) == expected_dimensions || throw(
        WFSPreparationError(:optical_formation, :shape,
            "Zernike detector output has the wrong prepared dimensions"))
    _require_zernike_rate_wavelength(output.metadata.spectral, wavelength_m)
    return output
end

function _require_zernike_input_geometry(front_end::ZernikeOpticalFrontEnd,
    input::PupilFunction)
    input.metadata.dimensions == (front_end.pupil_resolution,
        front_end.pupil_resolution) || throw(WFSPreparationError(
        :optical_formation, :shape,
        "Zernike pupil input dimensions differ from the prepared relay"))
    return nothing
end

function _require_zernike_input_geometry(front_end::ZernikeOpticalFrontEnd,
    input::ElectricField)
    input.metadata.dimensions == size(front_end.propagation.field) || throw(
        WFSPreparationError(:optical_formation, :shape,
            "Zernike ElectricField dimensions differ from the prepared diffraction grid"))
    return nothing
end

function prepare_wfs_optical_formation(front_end::ZernikeOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField}, output::IntensityMap)
    require_modulated_wfs_input(input)
    _require_zernike_front_end_source(front_end, input)
    _require_zernike_input_geometry(front_end, input)
    wavelength_m = _zernike_front_end_wavelength(front_end, input)
    _require_zernike_rate_map(output, zernike_rate_dimensions(front_end),
        wavelength_m)
    require_modulated_wfs_domains(front_end, input, output)
    eltype(front_end.propagation.pupil_intensity) ===
        output.metadata.numeric_type || throw(WFSPreparationError(
            :optical_formation, :numeric_type,
            "Zernike output precision differs from prepared propagation"))
    propagation = front_end.propagation
    binding = ZernikePropagationBinding(propagation.field,
        propagation.focal_field, propagation.pupil_field,
        propagation.phasor, propagation.phase_mask,
        propagation.nominal_frame, propagation.fft_plan,
        propagation.ifft_plan, propagation.revision)
    return PreparedZernikeOpticalFormation(front_end, input, output,
        binding)
end

function zernike_rate_map(sensor::ZernikeWFS,
    input::Union{PupilFunction,ElectricField}, source=nothing)
    return zernike_rate_map(ZernikeOpticalFrontEnd(sensor, source), input)
end

function zernike_rate_map(front_end::ZernikeOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField})
    wavelength_m = _zernike_front_end_wavelength(front_end, input)
    dimensions = zernike_rate_dimensions(front_end)
    T = eltype(front_end.propagation.pupil_intensity)
    values = similar(_modulated_input_storage(input), T, dimensions...)
    fill!(values, zero(T))
    normalized_sampling = T(front_end.binning / front_end.pupil_samples)
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
    plan::PreparedZernikeOpticalFormation, input, output)
    input === plan.input && output === plan.output || throw(
        WFSPreparationError(:optical_formation, :prepared_binding,
            "Zernike optical products do not match their prepared plan"))
    propagation = plan.front_end.propagation
    binding = plan.binding
    propagation.field === binding.field &&
        propagation.focal_field === binding.focal_field &&
        propagation.pupil_field === binding.pupil_field &&
        propagation.phasor === binding.phasor &&
        propagation.phase_mask === binding.phase_mask &&
        propagation.nominal_frame === binding.nominal_frame &&
        propagation.fft_plan === binding.fft_plan &&
        propagation.ifft_plan === binding.ifft_plan &&
        propagation.revision == binding.revision || throw(
        WFSPreparationError(:optical_formation, :prepared_binding,
            "Zernike propagation storage changed after preparation"))
    return nothing
end

@inline validate_wfs_optical_formation_binding(output::IntensityMap, input,
    plan::PreparedZernikeOpticalFormation) =
    _require_zernike_optical_binding(plan, input, output)

function _form_zernike_input_field!(front_end::ZernikeOpticalFrontEnd,
    input::PupilFunction)
    propagation = front_end.propagation
    T = eltype(propagation.pupil_intensity)
    n = front_end.pupil_resolution
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
    copyto!(front_end.propagation.field, input.values)
    return front_end.propagation.field
end

function _form_zernike_rate!(output::AbstractMatrix,
    front_end::ZernikeOpticalFrontEnd, input)
    propagation = front_end.propagation
    _form_zernike_input_field!(front_end, input)
    copyto!(propagation.focal_field, propagation.field)
    @. propagation.focal_field *= propagation.phasor
    execute_fft_plan!(propagation.focal_field, propagation.fft_plan)
    @. propagation.focal_field *= propagation.phase_mask
    copyto!(propagation.pupil_field, propagation.focal_field)
    execute_fft_plan!(propagation.pupil_field, propagation.ifft_plan)
    n = front_end.pupil_resolution
    pad = size(propagation.pupil_field, 1)
    ox = div(pad - n, 2)
    oy = div(pad - n, 2)
    @views @. propagation.pupil_intensity =
        abs2(propagation.pupil_field[ox+1:ox+n, oy+1:oy+n])
    sampling = div(n, front_end.pupil_samples)
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
    plan::PreparedZernikeOpticalFormation)
    validate_wfs_optical_formation_binding(output, input, plan)
    _form_zernike_rate!(output.values, plan.front_end, input)
    return output
end

function _zernike_calibration_binding(sensor::ZernikeWFS)
    state = sensor.estimator.state
    state.calibrated || throw(WFSPreparationError(:estimation, :estimator,
        "Zernike estimation requires explicit calibration"))
    return ZernikeCalibrationBinding(state.calibration_revision,
        state.calibration_wavelength, state.calibration_signature,
        state.reference_signal_2d, state.valid_mask)
end

function _require_zernike_calibration(estimator,
    binding::ZernikeCalibrationBinding)
    state = estimator.estimator.state
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
    state = sensor.estimator.state
    observation.metadata.dimensions == size(state.signal_2d) || throw(
        WFSPreparationError(:estimation, :shape,
            "Zernike observation has the wrong pupil-image dimensions"))
    size(measurement.storage) == size(state.slopes) || throw(
        WFSPreparationError(:estimation, :shape,
            "Zernike measurement has the wrong signal-vector dimensions"))
    _require_wfs_storage_domain(:estimation, observation.metadata,
        state.signal_2d, "Zernike observation")
    _require_wfs_storage_domain(:estimation, measurement.metadata,
        state.slopes, "Zernike measurement")
    measurement.metadata.numeric_type === eltype(state.slopes) || throw(
        WFSPreparationError(:estimation, :numeric_type,
            "Zernike measurement precision differs from its estimator"))
    _require_zernike_estimation_source(sensor.params.normalization, source)
    T = eltype(state.slopes)
    scale = T(normalization_scale)
    isfinite(scale) && scale >= zero(T) || throw(WFSPreparationError(
        :estimation, :radiometry,
        "Zernike normalization scale must be finite and nonnegative"))
    binding = _zernike_calibration_binding(sensor)
    return PreparedZernikeEstimator(sensor, observation, measurement,
        AcquiredObservationPath(), binding, source, scale)
end

@inline function _zernike_incidence_multiplier(sensor::ZernikeWFS,
    source::AbstractSource, normalization_scale)
    T = eltype(sensor.estimator.state.slopes)
    pupil_sample = T(sensor.params.pupil_diameter_m) /
        T(sensor.params.pupil_resolution)
    irradiance = T(_require_physical_photon_irradiance(source,
        "Zernike incidence normalization"))
    return irradiance * abs2(pupil_sample) * normalization_scale /
        T(zernike_normalization_count(sensor))
end

function _zernike_scalar_normalization(::MeanValidFluxNormalization,
    sensor::ZernikeWFS, frame::AbstractMatrix, source,
    normalization_scale::S) where {S<:AbstractFloat}
    state = sensor.estimator.state
    summed = zero(S)
    @inbounds for j in axes(frame, 2), i in axes(frame, 1)
        state.valid_mask[i, j] && (summed += S(frame[i, j]))
    end
    return max(summed / S(zernike_normalization_count(sensor)), eps(S))
end

function _zernike_scalar_normalization(::IncidenceFluxNormalization,
    sensor::ZernikeWFS, frame::AbstractMatrix, source,
    normalization_scale::S) where {S<:AbstractFloat}
    state = sensor.estimator.state
    summed = zero(S)
    @inbounds for j in axes(state.normalization_frame, 2),
            i in axes(state.normalization_frame, 1)
        state.valid_mask[i, j] &&
            (summed += state.normalization_frame[i, j])
    end
    return summed * _zernike_incidence_multiplier(sensor, source,
        normalization_scale)
end

function _estimate_zernike_signal!(::ScalarCPUStyle, sensor::ZernikeWFS,
    frame::AbstractMatrix{F}, source, normalization_scale::S) where {
    F<:Real,S<:AbstractFloat,
}
    state = sensor.estimator.state
    count = zernike_normalization_count(sensor)
    fill!(state.signal_2d, zero(S))
    count == 0 && (fill!(state.slopes, zero(S)); return state.slopes)
    normalization = _zernike_scalar_normalization(
        sensor.params.normalization, sensor, frame, source,
        normalization_scale)
    usable = isfinite(normalization) && normalization > eps(S)
    if !usable
        fill!(state.slopes, zero(S))
        return state.slopes
    end
    @inbounds for j in axes(frame, 2), i in axes(frame, 1)
        if state.valid_mask[i, j]
            state.signal_2d[i, j] = S(frame[i, j]) / normalization -
                state.reference_signal_2d[i, j]
        end
    end
    @inbounds for index in eachindex(state.valid_signal_indices)
        state.slopes[index] =
            state.signal_2d[state.valid_signal_indices[index]]
    end
    return state.slopes
end

function _queue_zernike_stage_normalization!(phase::KernelLaunchPhase,
    normalization::MeanValidFluxNormalization, sensor::ZernikeWFS,
    frame, source, normalization_scale)
    queue_zernike_masked_sum!(phase, sensor, frame)
    T = eltype(sensor.estimator.state.slopes)
    return inv(T(zernike_normalization_count(sensor))), true
end

function _queue_zernike_stage_normalization!(phase::KernelLaunchPhase,
    normalization::IncidenceFluxNormalization, sensor::ZernikeWFS,
    frame, source, normalization_scale)
    queue_zernike_masked_sum!(phase, sensor,
        sensor.estimator.state.normalization_frame)
    return _zernike_incidence_multiplier(sensor, source,
        normalization_scale), false
end

function _estimate_zernike_signal!(style::AcceleratorStyle,
    sensor::ZernikeWFS, frame::AbstractMatrix{F}, source,
    normalization_scale::S) where {F<:Real,S<:AbstractFloat}
    state = sensor.estimator.state
    count = zernike_normalization_count(sensor)
    if count == 0
        fill!(state.signal_2d, zero(S))
        fill!(state.slopes, zero(S))
        return state.slopes
    end
    phase = begin_kernel_phase(style)
    multiplier, clamp_to_epsilon = _queue_zernike_stage_normalization!(
        phase, sensor.params.normalization, sensor, frame, source,
        normalization_scale)
    queue_kernel!(phase, zernike_signal_kernel!, state.signal_2d, frame,
        state.valid_mask, state.reference_signal_2d,
        state.normalization_sum, multiplier, clamp_to_epsilon,
        size(frame, 1), size(frame, 2); ndrange=size(frame))
    queue_kernel!(phase, gather_zernike_signal_kernel!, state.slopes,
        state.signal_2d, state.valid_signal_indices, count;
        ndrange=count)
    finish_kernel_phase!(phase)
    return state.slopes
end

function estimate_wfs_measurement!(measurement::WFSMeasurement,
    observation::WFSObservation, plan::PreparedZernikeEstimator)
    measurement === plan.measurement && observation === plan.input || throw(
        WFSPreparationError(:estimation, :prepared_binding,
            "Zernike estimator storage does not match its plan"))
    sensor = plan.estimator
    _require_zernike_calibration(sensor, plan.calibration_binding)
    _estimate_zernike_signal!(execution_style(observation.storage), sensor,
        observation.storage, plan.source, plan.normalization_scale)
    copyto!(measurement.storage, sensor.estimator.state.slopes)
    return measurement
end

function validate_wfs_estimation_binding(measurement::WFSMeasurement, input,
    plan::PreparedZernikeEstimator)
    measurement === plan.measurement && input === plan.input || throw(
        WFSPreparationError(:estimation, :prepared_binding,
            "Zernike estimator storage does not match its plan"))
    return nothing
end

function set_zernike_calibration!(sensor::ZernikeWFS,
    reference::AbstractMatrix; wavelength_m::Real,
    signature::UInt=UInt(0))
    state = sensor.estimator.state
    size(reference) == size(state.reference_signal_2d) || throw(
        InvalidConfiguration(
            "Zernike calibration reference has the wrong dimensions"))
    eltype(reference) <: Real || throw(InvalidConfiguration(
        "Zernike calibration reference must contain real values"))
    typeof(backend(reference)) === typeof(backend(state.reference_signal_2d)) ||
        throw(InvalidConfiguration(
            "Zernike calibration reference backend differs from the estimator"))
    plane_device(reference) == plane_device(state.reference_signal_2d) ||
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
    fill!(state.signal_2d, zero(T))
    fill!(state.slopes, zero(T))
    state.calibrated = true
    state.calibration_wavelength = wavelength_value
    state.calibration_signature = signature
    state.calibration_revision += UInt(1)
    return sensor
end

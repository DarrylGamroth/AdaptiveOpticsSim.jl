#
# Prepared BioEdge WFS stages
#

struct PreparedBioEdgeOpticalFormation{F,I,O,L<:AbstractPreparedFourPupilLGS}
    front_end::F
    input::I
    output::O
    lgs_model::L
end

struct PreparedBioEdgeOpticalBundleFormation{P<:Tuple,I,O}
    plans::P
    input::I
    output::O
end

struct BioEdgeCalibrationBinding{T<:AbstractFloat,R,A}
    revision::UInt
    wavelength_m::T
    signature::UInt
    reference_signal::R
    valid_support::A
end

struct PreparedBioEdgeEstimator{W,I,M,P<:AbstractWFSMeasurementPath,C,S,T}
    sensor::W
    input::I
    measurement::M
    path::P
    calibration_binding::C
    source::S
    normalization_scale::T
end

@inline wfs_measurement_path(plan::PreparedBioEdgeEstimator) = plan.path

@inline function bioedge_output_sampling_factor(
    front_end::BioEdgeOpticalFrontEnd, pupil_resolution::Int)
    pupil_resolution % front_end.pupil_samples == 0 ||
        throw(WFSPreparationError(:optical_formation, :shape,
            "BioEdge pupil resolution must be divisible by pupil_samples"))
    pupil_sample = div(pupil_resolution, front_end.pupil_samples)
    return front_end.binning == 1 ? pupil_sample :
        2 * pupil_sample * front_end.binning
end

function bioedge_output_dimensions(front_end::BioEdgeOpticalFrontEnd,
    pupil_resolution::Int)
    factor = bioedge_output_sampling_factor(front_end, pupil_resolution)
    side = size(front_end.propagation.intensity, 1)
    side % factor == 0 || throw(WFSPreparationError(:optical_formation,
        :shape, "BioEdge sampling does not evenly divide the detector plane"))
    output_side = div(side, factor)
    return (output_side, output_side)
end

@inline function _bioedge_front_end_wavelength(
    front_end::BioEdgeOpticalFrontEnd, input::PupilFunction)
    return modulated_input_wavelength(input, front_end.source)
end

@inline function _bioedge_front_end_wavelength(
    ::BioEdgeOpticalFrontEnd, input::ElectricField)
    return modulated_input_wavelength(input)
end

function _require_bioedge_source(front_end::BioEdgeOpticalFrontEnd,
    ::PupilFunction)
    source = front_end.source
    source === nothing && throw(WFSPreparationError(:optical_formation,
        :radiometry, "BioEdge PupilFunction formation requires a source"))
    return _require_single_bioedge_source(source)
end

@inline _require_single_bioedge_source(source) = source

function _require_single_bioedge_source(source::SpectralSource)
    throw(WFSPreparationError(
        :optical_formation, :plane_count,
        "spectral BioEdge formation requires an OpticalProductBundle"))
end

function _require_single_bioedge_source(source::Asterism)
    throw(WFSPreparationError(:optical_formation,
        :plane_count,
        "asterism BioEdge formation requires path-local pupil inputs"))
end

function _require_single_bioedge_source(source::ExtendedSource)
    throw(WFSPreparationError(
        :optical_formation, :plane_count,
        "extended BioEdge formation requires path-local pupil inputs"))
end

function _require_bioedge_source(front_end::BioEdgeOpticalFrontEnd,
    ::ElectricField)
    front_end.source === nothing || throw(WFSPreparationError(
        :optical_formation, :radiometry,
        "photon-rate ElectricField input must not also supply a source"))
    return nothing
end

function prepare_wfs_optical_formation(front_end::BioEdgeOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField}, output::IntensityMap)
    require_modulated_wfs_input(input)
    _require_bioedge_source(front_end, input)
    resolution = input.metadata.dimensions[1]
    input.metadata.dimensions == (resolution, resolution) ||
        throw(WFSPreparationError(:optical_formation, :shape,
            "BioEdge pupil input must be square"))
    size(front_end.modulation.phases, 1) == resolution ||
        throw(WFSPreparationError(:optical_formation, :shape,
            "BioEdge modulation was prepared for another pupil resolution"))
    expected = bioedge_output_dimensions(front_end, resolution)
    wavelength_m = _bioedge_front_end_wavelength(front_end, input)
    require_four_pupil_rate_map(output, expected, wavelength_m)
    require_modulated_wfs_domains(front_end, input, output)
    eltype(front_end.propagation.intensity) ===
        output.metadata.numeric_type || throw(WFSPreparationError(
            :optical_formation, :numeric_type,
            "BioEdge output precision differs from prepared propagation"))
    lgs_model = prepare_four_pupil_lgs(front_end.source, input, front_end)
    return PreparedBioEdgeOpticalFormation(front_end, input, output,
        lgs_model)
end

function prepare_wfs_optical_formation(front_end::BioEdgeOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField},
    output::OpticalProductBundle)
    return prepare_bioedge_optical_bundle(front_end, input, output,
        front_end.source)
end

function prepare_bioedge_optical_bundle(front_end::BioEdgeOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField}, output::OpticalProductBundle,
    source::SpectralSource)
    samples = spectral_bundle(source).samples
    length(output) == length(samples) || throw(WFSPreparationError(
        :optical_formation, :plane_count,
        "BioEdge spectral output count does not match the source"))
    T = eltype(front_end.propagation.intensity)
    plans = ntuple(length(samples)) do index
        sample = samples[index]
        component = FourPupilSpectralComponent(source.source,
            T(sample.wavelength),
            T(photon_irradiance(source)) * T(sample.weight))
        prepare_wfs_optical_formation(
            bioedge_front_end_with_source(front_end, component), input,
            output[index])
    end
    return PreparedBioEdgeOpticalBundleFormation(plans, input, output)
end

function prepare_wfs_optical_formation(front_end::BioEdgeOpticalFrontEnd,
    inputs::Union{Tuple,AbstractVector}, output::OpticalProductBundle)
    return prepare_bioedge_optical_bundle(front_end, inputs, output,
        front_end.source)
end

function prepare_bioedge_optical_bundle(front_end::BioEdgeOpticalFrontEnd,
    inputs::Union{Tuple,AbstractVector}, output::OpticalProductBundle,
    source::Union{Asterism,ExtendedSource})
    sources = four_pupil_path_sources(source)
    length(inputs) == length(sources) || throw(WFSPreparationError(
        :optical_formation, :plane_count,
        "BioEdge path-local pupil count does not match the source count"))
    length(output) == length(sources) || throw(WFSPreparationError(
        :optical_formation, :plane_count,
        "BioEdge path-local output count does not match the source count"))
    isempty(sources) && throw(WFSPreparationError(:optical_formation,
        :plane_count, "BioEdge path-local source collection is empty"))
    plans = ntuple(length(sources)) do index
        prepare_wfs_optical_formation(
            bioedge_front_end_with_source(front_end, sources[index]),
            inputs[index], output[index])
    end
    return PreparedBioEdgeOpticalBundleFormation(plans, inputs, output)
end

function prepare_bioedge_optical_bundle(front_end, input, output, source)
    throw(WFSPreparationError(
        :optical_formation, :plane_count,
        "BioEdge product bundles require a spectral or path-expanded source"))
end

function _bioedge_native_rate!(front_end::BioEdgeOpticalFrontEnd,
    input::PupilFunction, lgs_model::AbstractPreparedFourPupilLGS)
    propagation = front_end.propagation
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
        accumulate_bioedge_masked_pupils!(propagation.intensity, front_end,
            lgs_model)
    end
    return propagation.intensity
end

function _bioedge_native_rate!(front_end::BioEdgeOpticalFrontEnd,
    input::ElectricField, lgs_model::AbstractPreparedFourPupilLGS)
    propagation = front_end.propagation
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
        accumulate_bioedge_masked_pupils!(propagation.intensity, front_end,
            lgs_model)
    end
    return propagation.intensity
end

function form_wfs_optical_products!(output::IntensityMap,
    input::Union{PupilFunction,ElectricField},
    plan::PreparedBioEdgeOpticalFormation)
    output === plan.output && input === plan.input ||
        throw(WFSPreparationError(:optical_formation, :prepared_binding,
            "BioEdge optical products do not match prepared storage"))
    native = _bioedge_native_rate!(plan.front_end, input, plan.lgs_model)
    factor = bioedge_output_sampling_factor(plan.front_end,
        input.metadata.dimensions[1])
    bin2d!(output.values, native, factor)
    return output
end

function form_wfs_optical_products!(output::OpticalProductBundle, input,
    plan::PreparedBioEdgeOpticalBundleFormation)
    output === plan.output && input === plan.input ||
        throw(WFSPreparationError(:optical_formation, :prepared_binding,
            "BioEdge spectral products do not match prepared storage"))
    return form_four_pupil_bundle!(output, input, plan.plans)
end

function bioedge_rate_map(sensor::BioEdgeWFS{<:Diffractive},
    inputs::Union{Tuple,AbstractVector}, source)
    return bioedge_rate_map(BioEdgeOpticalFrontEnd(sensor, source), inputs)
end

function bioedge_rate_map(front_end::BioEdgeOpticalFrontEnd,
    inputs::Union{Tuple,AbstractVector})
    return bioedge_path_rate_bundle(front_end, inputs, front_end.source)
end

function bioedge_path_rate_bundle(front_end::BioEdgeOpticalFrontEnd,
    inputs::Union{Tuple,AbstractVector},
    source::Union{Asterism,ExtendedSource})
    sources = four_pupil_path_sources(source)
    length(inputs) == length(sources) || throw(WFSPreparationError(
        :optical_formation, :plane_count,
        "BioEdge path-local pupil count does not match the source count"))
    isempty(sources) && throw(WFSPreparationError(:optical_formation,
        :plane_count, "BioEdge path-local source collection is empty"))
    first_map = bioedge_rate_map(
        bioedge_front_end_with_source(front_end, sources[1]), inputs[1])
    maps = Vector{typeof(first_map)}(undef, length(sources))
    maps[1] = first_map
    @inbounds for index in 2:length(sources)
        maps[index] = bioedge_rate_map(
            bioedge_front_end_with_source(front_end, sources[index]),
            inputs[index])
    end
    return OpticalProductBundle(maps)
end

function bioedge_path_rate_bundle(front_end, inputs, source)
    throw(WFSPreparationError(:optical_formation, :plane_count,
        "path-local BioEdge inputs require an Asterism or ExtendedSource"))
end

function bioedge_rate_map(sensor::BioEdgeWFS{<:Diffractive},
    input::Union{PupilFunction,ElectricField}, source=nothing)
    return bioedge_rate_map(BioEdgeOpticalFrontEnd(sensor, source), input)
end

function bioedge_rate_map(front_end::BioEdgeOpticalFrontEnd,
    input::Union{PupilFunction,ElectricField})
    return _bioedge_rate_map(front_end, input, front_end.source)
end

@inline _bioedge_rate_map(front_end::BioEdgeOpticalFrontEnd, input,
    source::SpectralSource) =
    _bioedge_spectral_rate_bundle(front_end, input, source)

function _bioedge_rate_map(front_end::BioEdgeOpticalFrontEnd, input,
    source::Union{Asterism,ExtendedSource})
    throw(WFSPreparationError(:optical_formation, :plane_count,
        "path-expanded BioEdge sources require path-local pupil inputs"))
end

function _bioedge_rate_map(front_end::BioEdgeOpticalFrontEnd, input, source)
    wavelength_m = _bioedge_front_end_wavelength(front_end, input)
    dimensions = bioedge_output_dimensions(front_end,
        input.metadata.dimensions[1])
    T = eltype(front_end.propagation.intensity)
    values = similar(_modulated_input_storage(input), T, dimensions...)
    fill!(values, zero(T))
    factor = bioedge_output_sampling_factor(front_end,
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

@inline _require_bioedge_estimation_source(::WFSNormalization, source) =
    nothing

function _require_bioedge_estimation_source(
    ::IncidenceFluxNormalization, ::Nothing)
    throw(WFSPreparationError(:estimation, :radiometry,
        "incidence-normalized BioEdge estimation requires a source"))
end

function _bioedge_spectral_rate_bundle(front_end::BioEdgeOpticalFrontEnd,
    input, source::SpectralSource)
    samples = spectral_bundle(source).samples
    T = eltype(front_end.propagation.intensity)
    function component_map(sample)
        component = FourPupilSpectralComponent(source.source,
            T(sample.wavelength),
            T(photon_irradiance(source)) * T(sample.weight))
        component_front_end = bioedge_front_end_with_source(front_end,
            component)
        return bioedge_rate_map(component_front_end, input)
    end
    first_map = component_map(first(samples))
    maps = Vector{typeof(first_map)}(undef, length(samples))
    maps[1] = first_map
    @inbounds for index in 2:length(samples)
        maps[index] = component_map(samples[index])
    end
    return OpticalProductBundle(maps)
end

function _bioedge_calibration_binding(sensor::BioEdgeWFS)
    state = sensor.estimator.state
    state.calibrated || throw(WFSPreparationError(:estimation, :estimator,
        "BioEdge estimation requires explicit calibration"))
    return BioEdgeCalibrationBinding(state.calibration_revision,
        state.calibration_wavelength,
        state.calibration_signature, state.reference_signal_2d,
        state.valid_i4q)
end

function _require_bioedge_calibration(sensor::BioEdgeWFS,
    binding::BioEdgeCalibrationBinding)
    state = sensor.estimator.state
    state.calibrated &&
        state.calibration_revision == binding.revision &&
        state.calibration_wavelength == binding.wavelength_m &&
        state.calibration_signature == binding.signature &&
        state.reference_signal_2d === binding.reference_signal &&
        state.valid_i4q === binding.valid_support ||
        throw(WFSPreparationError(:estimation, :prepared_binding,
            "BioEdge calibration changed after estimator preparation"))
    return nothing
end

function prepare_wfs_estimation(sensor::BioEdgeWFS{<:Diffractive},
    observation::WFSObservation, measurement::WFSMeasurement;
    source=nothing, normalization_scale::Real=1)
    validate_wfs_observation(observation)
    validate_wfs_measurement(measurement)
    isequal(observation.metadata.layout, :four_pupil_mosaic) ||
        throw(WFSPreparationError(:estimation, :detector_mapping,
            "BioEdge estimator requires :four_pupil_mosaic layout"))
    isequal(measurement.metadata.kind, :differential_slopes) ||
        throw(WFSPreparationError(:estimation, :estimator,
            "BioEdge measurement kind must be :differential_slopes"))
    isequal(measurement.units, :dimensionless) ||
        throw(WFSPreparationError(:estimation, :units,
            "BioEdge differential slopes are dimensionless"))
    measurement.metadata.numeric_type <: AbstractFloat ||
        throw(WFSPreparationError(:estimation, :numeric_type,
            "BioEdge measurement storage must be floating point"))
    resize_bioedge_signal_buffers!(sensor, size(observation.storage, 1))
    size(measurement.storage) == size(sensor.estimator.state.slopes) ||
        throw(WFSPreparationError(:estimation, :shape,
            "BioEdge measurement storage has the wrong slope shape"))
    _require_bioedge_estimation_source(
        sensor.estimator.params.normalization, source)
    scale = eltype(sensor.estimator.state.slopes)(normalization_scale)
    isfinite(scale) && scale >= zero(scale) || throw(WFSPreparationError(
        :estimation, :radiometry,
        "BioEdge normalization scale must be finite and nonnegative"))
    binding = _bioedge_calibration_binding(sensor)
    return PreparedBioEdgeEstimator(sensor, observation, measurement,
        AcquiredObservationPath(), binding, source, scale)
end

function estimate_wfs_measurement!(measurement::WFSMeasurement,
    observation::WFSObservation,
    plan::PreparedBioEdgeEstimator{
        <:Any,<:Any,<:Any,<:AcquiredObservationPath})
    measurement === plan.measurement && observation === plan.input ||
        throw(WFSPreparationError(:estimation, :prepared_binding,
            "BioEdge estimator storage does not match its plan"))
    sensor = plan.sensor
    _require_bioedge_calibration(sensor, plan.calibration_binding)
    bioedge_signal!(execution_style(observation.storage), sensor,
        observation.storage, plan.source, plan.normalization_scale)
    copyto!(measurement.storage, sensor.estimator.state.slopes)
    return measurement
end

function prepare_wfs_estimation(sensor::BioEdgeWFS{<:Geometric},
    input::PupilFunction, measurement::WFSMeasurement)
    require_modulated_wfs_input(input)
    validate_wfs_measurement(measurement)
    input.metadata.dimensions == (sensor.estimator.params.pupil_resolution,
        sensor.estimator.params.pupil_resolution) ||
        throw(WFSPreparationError(:estimation, :shape,
            "geometric BioEdge input has the wrong pupil dimensions"))
    isequal(measurement.metadata.kind, :geometric_slopes) ||
        throw(WFSPreparationError(:estimation, :estimator,
            "geometric BioEdge measurement kind must be :geometric_slopes"))
    isequal(measurement.units, :metre) || throw(WFSPreparationError(
        :estimation, :units,
        "geometric BioEdge OPD differences are expressed in metres"))
    size(measurement.storage) == size(sensor.estimator.state.slopes) ||
        throw(WFSPreparationError(:estimation, :shape,
            "geometric BioEdge measurement has the wrong slope shape"))
    return PreparedBioEdgeEstimator(sensor, input, measurement,
        DirectMeasurementPath(), nothing, nothing,
        one(eltype(sensor.estimator.state.slopes)))
end

function estimate_wfs_measurement!(measurement::WFSMeasurement,
    input::PupilFunction,
    plan::PreparedBioEdgeEstimator{
        <:Any,<:Any,<:Any,<:DirectMeasurementPath})
    measurement === plan.measurement && input === plan.input ||
        throw(WFSPreparationError(:estimation, :prepared_binding,
            "geometric BioEdge estimator storage does not match its plan"))
    sensor = plan.sensor
    state = sensor.estimator.state
    edge_geometric_slopes!(state.slopes, input.opd, state.valid_mask,
        state.edge_mask)
    @. state.slopes *= state.optical_gain
    copyto!(measurement.storage, state.slopes)
    return measurement
end

function set_bioedge_calibration!(sensor::BioEdgeWFS,
    reference::AbstractMatrix; wavelength_m::Real,
    signature::UInt=UInt(0), valid_support=nothing)
    state = sensor.estimator.state
    size(reference) == size(state.reference_signal_2d) ||
        throw(DimensionMismatchError(
            "BioEdge reference dimensions do not match estimator storage"))
    require_same_backend(state.reference_signal_2d, reference)
    reference_host = Array(reference)
    all(isfinite, reference_host) || throw(InvalidConfiguration(
        "BioEdge calibration reference must contain only finite values"))
    wavelength_value = eltype(state.slopes)(wavelength_m)
    isfinite(wavelength_value) && wavelength_value > zero(wavelength_value) ||
        throw(InvalidConfiguration(
            "BioEdge calibration wavelength must be finite and positive"))
    support_host = _prepare_bioedge_calibration_support(sensor, valid_support)
    copyto!(state.reference_signal_2d, reference_host)
    if support_host === nothing
        fill!(state.valid_i4q, true)
    else
        copyto!(state.valid_i4q, support_host)
    end
    update_bioedge_valid_signal!(sensor)
    update_bioedge_valid_signal_indices!(sensor)
    resize_bioedge_slope_buffers!(sensor)
    state.calibration_wavelength = wavelength_value
    state.calibration_signature = signature
    state.calibrated = true
    state.calibration_revision += UInt(1)
    return sensor
end

function _prepare_bioedge_calibration_support(sensor::BioEdgeWFS, ::Nothing)
    if !iszero(sensor.estimator.params.light_ratio)
        throw(InvalidConfiguration(
            "nonzero BioEdge light_ratio requires explicit valid_support"))
    end
    return nothing
end

function _prepare_bioedge_calibration_support(sensor::BioEdgeWFS,
    valid_support::AbstractMatrix{Bool})
    state = sensor.estimator.state
    size(valid_support) == size(state.valid_i4q) ||
        throw(DimensionMismatchError(
            "BioEdge calibration support has the wrong dimensions"))
    require_same_backend(state.valid_i4q, valid_support)
    support_host = Array(valid_support)
    any(support_host) || throw(InvalidConfiguration(
        "BioEdge calibration support must select at least one sample"))
    return support_host
end

function _prepare_bioedge_calibration_support(sensor::BioEdgeWFS,
    valid_support)
    throw(InvalidConfiguration(
        "BioEdge calibration support must be a Boolean matrix"))
end

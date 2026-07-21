#
# Prepared plant ownership
#
# Preparation is intentionally model-extensible rather than a universal
# optical graph. Cold model definitions dispatch to concrete preparation
# methods; the resulting owners bind exact products and prepared stage plans.
#

"""Semantic contract for when an optical path result samples the plant."""
abstract type AbstractOpticalSamplingContract end

"""The path result is a photon-arrival-rate sample at one plant instant."""
struct InstantaneousOpticalSample <: AbstractOpticalSamplingContract end

struct _PathResultKeyToken end
const _PATH_RESULT_KEY_TOKEN = _PathResultKeyToken()

"""
    PathResultKey

Run-owned compatibility description for one prepared optical result. The key
records physical source geometry and radiometry, spectral sampling, optical
and propagation model identities, the sampling contract, output-plane
semantics, revisions, backend, and physical device. Acquisitions compare all
fields during preparation and retain the exact key and result binding.
"""
struct PathResultKey{G,S,R,M,C,P,O,V,B<:AbstractArrayBackend,
    D<:AbstractPlaneDevice}
    source_geometry::G
    spectral_sampling::S
    radiometry::R
    optical_model::M
    sampling_contract::C
    propagation_model::P
    output_plane::O
    revisions::V
    backend::B
    device::D

    function PathResultKey(::_PathResultKeyToken, source_geometry,
        spectral_sampling, radiometry, optical_model, sampling_contract,
        propagation_model, output_plane, revisions,
        backend::B, device::D) where {
        B<:AbstractArrayBackend,D<:AbstractPlaneDevice,
    }
        geometry_snapshot = deepcopy(source_geometry)
        spectral_snapshot = deepcopy(spectral_sampling)
        radiometry_snapshot = deepcopy(radiometry)
        optical_snapshot = deepcopy(optical_model)
        sampling_snapshot = deepcopy(sampling_contract)
        propagation_snapshot = deepcopy(propagation_model)
        output_snapshot = deepcopy(output_plane)
        revision_snapshot = deepcopy(revisions)
        return new{
            typeof(geometry_snapshot),typeof(spectral_snapshot),
            typeof(radiometry_snapshot),typeof(optical_snapshot),
            typeof(sampling_snapshot),typeof(propagation_snapshot),
            typeof(output_snapshot),typeof(revision_snapshot),B,D,
        }(geometry_snapshot, spectral_snapshot, radiometry_snapshot,
            optical_snapshot, sampling_snapshot, propagation_snapshot,
            output_snapshot, revision_snapshot, backend, device)
    end
end

function PathResultKey(source_geometry, spectral_sampling, radiometry,
    optical_model, sampling_contract, propagation_model, output_plane,
    revisions, backend::AbstractArrayBackend, device::AbstractPlaneDevice)
    return PathResultKey(_PATH_RESULT_KEY_TOKEN, source_geometry,
        spectral_sampling, radiometry, optical_model, sampling_contract,
        propagation_model, output_plane, revisions, backend, device)
end

const _PATH_RESULT_KEY_SNAPSHOT_FIELDS = (
    :source_geometry,
    :spectral_sampling,
    :radiometry,
    :optical_model,
    :sampling_contract,
    :propagation_model,
    :output_plane,
    :revisions,
)

"""Return defensive copies for value descriptions retained by the key."""
function Base.getproperty(key::PathResultKey, name::Symbol)
    name in _PATH_RESULT_KEY_SNAPSHOT_FIELDS &&
        return deepcopy(getfield(key, name))
    return getfield(key, name)
end

function Base.isequal(left::PathResultKey, right::PathResultKey)
    return isequal(getfield(left, :source_geometry),
            getfield(right, :source_geometry)) &&
        isequal(getfield(left, :spectral_sampling),
            getfield(right, :spectral_sampling)) &&
        isequal(getfield(left, :radiometry), getfield(right, :radiometry)) &&
        isequal(getfield(left, :optical_model),
            getfield(right, :optical_model)) &&
        isequal(getfield(left, :sampling_contract),
            getfield(right, :sampling_contract)) &&
        isequal(getfield(left, :propagation_model),
            getfield(right, :propagation_model)) &&
        isequal(getfield(left, :output_plane),
            getfield(right, :output_plane)) &&
        isequal(getfield(left, :revisions), getfield(right, :revisions)) &&
        typeof(getfield(left, :backend)) ===
            typeof(getfield(right, :backend)) &&
        isequal(getfield(left, :device), getfield(right, :device))
end

Base.:(==)(left::PathResultKey, right::PathResultKey) =
    isequal(left, right)

function Base.hash(key::PathResultKey, seed::UInt)
    return hash((getfield(key, :source_geometry),
        getfield(key, :spectral_sampling), getfield(key, :radiometry),
        getfield(key, :optical_model), getfield(key, :sampling_contract),
        getfield(key, :propagation_model), getfield(key, :output_plane),
        getfield(key, :revisions), typeof(getfield(key, :backend)),
        getfield(key, :device)), seed)
end

@inline function _leaf_source_geometry_key(src::Source)
    return (
        kind=Source,
        direction_arcsec=coordinates_xy_arcsec(src),
        height_m=source_height_m(src),
    )
end

@inline function _leaf_source_geometry_key(src::LGSSource)
    params = src.params
    profile = isnothing(params.na_profile) ? nothing : copy(params.na_profile)
    return (
        kind=LGSSource,
        direction_arcsec=coordinates_xy_arcsec(src),
        height_m=source_height_m(src),
        laser_coordinates_m=params.laser_coordinates,
        elongation_factor=params.elongation_factor,
        uplink_fwhm_arcsec=params.fwhm_spot_up,
        sodium_profile=profile,
    )
end

@inline path_source_geometry_key(src::Union{Source,LGSSource}) =
    _leaf_source_geometry_key(src)
@inline path_source_geometry_key(src::SpectralSource) =
    _leaf_source_geometry_key(src.source)

function path_source_geometry_key(src::Asterism)
    return map(path_source_geometry_key, src.sources)
end

function path_source_geometry_key(src::ExtendedSource)
    return map(path_source_geometry_key, src.quadrature.sources)
end

function path_source_geometry_key(src::AbstractSource)
    throw(PlantPreparationError(:path, :source_geometry,
        "source type $(typeof(src)) must implement a prepared path geometry key"))
end

@inline function _leaf_source_spectral_key(src::Union{Source,LGSSource})
    return ((wavelength_m=wavelength(src), weight=one(wavelength(src))),)
end

@inline path_source_spectral_key(src::Union{Source,LGSSource}) =
    _leaf_source_spectral_key(src)

function path_source_spectral_key(src::SpectralSource)
    return map(src.bundle.samples) do sample
        (wavelength_m=sample.wavelength, weight=sample.weight)
    end
end

function path_source_spectral_key(src::Asterism)
    return map(path_source_spectral_key, src.sources)
end

@inline path_source_spectral_key(src::ExtendedSource) =
    _leaf_source_spectral_key(src.source)

function path_source_spectral_key(src::AbstractSource)
    throw(PlantPreparationError(:path, :spectral_sampling,
        "source type $(typeof(src)) must implement a prepared spectral key"))
end

@inline function _leaf_source_radiometry_key(src::Union{Source,LGSSource})
    return (
        policy=typeof(source_radiometry(src)),
        value=source_radiometric_value(src),
    )
end

@inline path_source_radiometry_key(src::Union{Source,LGSSource}) =
    _leaf_source_radiometry_key(src)
@inline path_source_radiometry_key(src::SpectralSource) =
    _leaf_source_radiometry_key(src.source)

function path_source_radiometry_key(src::Asterism)
    return map(path_source_radiometry_key, src.sources)
end

function path_source_radiometry_key(src::ExtendedSource)
    return map(path_source_radiometry_key, src.quadrature.sources)
end

function path_source_radiometry_key(src::AbstractSource)
    throw(PlantPreparationError(:path, :radiometry,
        "source type $(typeof(src)) must implement a prepared radiometry key"))
end

@inline _require_path_result_plane(::FocalPlane) = nothing
@inline _require_path_result_plane(::DetectorPlane) = nothing

function _require_path_result_plane(kind::AbstractOpticalPlaneKind)
    throw(PlantPreparationError(:path, :output_plane,
        "acquisition-facing path results must be on a focal or detector plane; got $(typeof(kind))"))
end

@inline _require_path_result_rate(::PhotonRateNormalization) = nothing

function _require_path_result_rate(normalization::AbstractOpticalNormalization)
    throw(PlantPreparationError(:path, :radiometry,
        "acquisition-facing path results must use photon-rate normalization; got $(typeof(normalization))"))
end

@inline _require_path_result_measure(::SpatialDensityMeasure) = nothing
@inline _require_path_result_measure(::CellIntegratedMeasure) = nothing

function _require_path_result_measure(measure::AbstractSpatialMeasure)
    throw(PlantPreparationError(:path, :radiometry,
        "acquisition-facing path results must use spatial-density or cell-integrated samples; got $(typeof(measure))"))
end

@inline _require_path_result_coherence(::IncoherentIntensityAddition) =
    nothing

function _require_path_result_coherence(coherence::AbstractCombinationPolicy)
    throw(PlantPreparationError(:path, :radiometry,
        "acquisition-facing path results must declare incoherent intensity addition; got $(typeof(coherence))"))
end

@inline _require_path_result_spectral(::MonochromaticChannel) = nothing
@inline _require_path_result_spectral(::IntegratedSpectralChannel) = nothing

function _require_path_result_spectral(spectral::AbstractSpectralCoordinate)
    throw(PlantPreparationError(:path, :spectral_sampling,
        "acquisition-facing path results require a declared spectral channel; got $(typeof(spectral))"))
end

function _path_output_contract(product::IntensityMap)
    metadata = validate_plane_storage(product.metadata, product.values;
        label="prepared optical-path result")
    _require_path_result_plane(metadata.kind)
    _require_path_result_rate(metadata.normalization)
    _require_path_result_measure(metadata.spatial_measure)
    _require_path_result_coherence(metadata.coherence)
    _require_path_result_spectral(metadata.spectral)
    return (
        kind=metadata.kind,
        coordinate_domain=metadata.coordinate_domain,
        dimensions=metadata.dimensions,
        sampling=metadata.sampling,
        origin=metadata.origin,
        centering=metadata.centering,
        orientation=metadata.orientation,
        spectral=metadata.spectral,
        numeric_type=metadata.numeric_type,
        normalization=metadata.normalization,
        spatial_measure=metadata.spatial_measure,
        coherence=metadata.coherence,
    )
end

@inline _path_output_contract(bundle::OpticalProductBundle) =
    _path_output_contract(bundle.products)
@inline _path_output_contract(products::Tuple) =
    map(_path_output_contract, products)

function _path_output_contract(products::_FixedOpticalProductVector)
    isempty(products) && throw(PlantPreparationError(:path, :output_plane,
        "prepared optical-path result bundle must not be empty"))
    return map(_path_output_contract, products)
end

function _path_output_contract(result)
    throw(PlantPreparationError(:path, :output_plane,
        "prepared optical-path result must be an IntensityMap or concrete tuple/bundle of IntensityMap values; got $(typeof(result))"))
end

@inline _first_path_result(result::IntensityMap) = result
@inline _first_path_result(bundle::OpticalProductBundle) =
    _first_path_result(bundle.products)
@inline _first_path_result(products::Tuple{<:Any,Vararg}) =
    _first_path_result(first(products))

function _first_path_result(products::_FixedOpticalProductVector)
    isempty(products) && throw(PlantPreparationError(:path, :output_plane,
        "prepared optical-path result bundle must not be empty"))
    return _first_path_result(first(products))
end

function _first_path_result(::Tuple{})
    throw(PlantPreparationError(:path, :output_plane,
        "prepared optical-path result tuple must not be empty"))
end

@inline _require_path_result_domain(::Tuple{}, ::AbstractArrayBackend,
    ::AbstractPlaneDevice) = nothing

@inline function _require_path_result_domain(product::IntensityMap,
    selector::AbstractArrayBackend, device::AbstractPlaneDevice)
    typeof(backend(product)) === typeof(selector) || throw(
        PlantPreparationError(:path, :backend,
            "prepared path result leaves must use one array backend"))
    plane_device(product.values) == device || throw(PlantPreparationError(
        :path, :device,
        "prepared path result leaves must occupy one physical device"))
    return nothing
end

@inline _require_path_result_domain(bundle::OpticalProductBundle,
    selector::AbstractArrayBackend, device::AbstractPlaneDevice) =
    _require_path_result_domain(bundle.products, selector, device)

@inline function _require_path_result_domain(products::Tuple,
    selector::AbstractArrayBackend, device::AbstractPlaneDevice)
    _require_path_result_domain(first(products), selector, device)
    return _require_path_result_domain(Base.tail(products), selector, device)
end

function _require_path_result_domain(products::_FixedOpticalProductVector,
    selector::AbstractArrayBackend, device::AbstractPlaneDevice)
    @inbounds for product in products
        _require_path_result_domain(product, selector, device)
    end
    return nothing
end

function _validate_path_input(input::PupilFunction)
    typeof(input.metadata.kind) === PupilPlane || throw(
        PlantPreparationError(:path, :input_plane,
            "prepared path PupilFunction input must be on a pupil plane"))
    validate_plane_storage(input.metadata, input.amplitude;
        label="prepared path pupil amplitude")
    validate_plane_storage(input.metadata, input.opd;
        label="prepared path pupil OPD")
    return input
end

function _validate_path_input(input::ElectricField)
    typeof(input.metadata.kind) === PupilPlane || throw(
        PlantPreparationError(:path, :input_plane,
            "prepared path ElectricField input must be on a pupil plane"))
    validate_plane_storage(input.metadata, input.values;
        label="prepared path electric field")
    return input
end

@inline _validate_path_input(::Tuple{}) = throw(PlantPreparationError(
    :path, :input_plane, "prepared path input tuple must not be empty"))

@inline function _validate_path_input(inputs::Tuple)
    _validate_path_input(first(inputs))
    _validate_remaining_path_inputs(Base.tail(inputs))
    return inputs
end

@inline _validate_remaining_path_inputs(::Tuple{}) = nothing

@inline function _validate_remaining_path_inputs(inputs::Tuple)
    _validate_path_input(first(inputs))
    return _validate_remaining_path_inputs(Base.tail(inputs))
end

function _validate_path_input(input)
    throw(PlantPreparationError(:path, :input_plane,
        "prepared path input must be a pupil function, pupil-plane electric field, or concrete tuple of them; got $(typeof(input))"))
end

@inline _path_input_storage(input::PupilFunction) = input.opd
@inline _path_input_storage(input::ElectricField) = input.values

@inline function _require_path_input_domain(
    input::Union{PupilFunction,ElectricField},
    selector::AbstractArrayBackend,
    device::AbstractPlaneDevice,
)
    typeof(backend(input)) === typeof(selector) || throw(
        PlantPreparationError(:path, :backend,
            "prepared path input and result backends differ"))
    plane_device(_path_input_storage(input)) == device || throw(
        PlantPreparationError(:path, :device,
            "prepared path input and result occupy different physical devices"))
    return nothing
end

@inline _require_path_input_domain(::Tuple{}, ::AbstractArrayBackend,
    ::AbstractPlaneDevice) = nothing

@inline function _require_path_input_domain(inputs::Tuple,
    selector::AbstractArrayBackend, device::AbstractPlaneDevice)
    _require_path_input_domain(first(inputs), selector, device)
    return _require_path_input_domain(Base.tail(inputs), selector, device)
end

@inline _require_path_input_revision(::ElectricField, ::UInt) = nothing

function _require_path_input_revision(input::PupilFunction, revision::UInt)
    aperture_revision(input) == revision || throw(PlantPreparationError(
        :path, :revision,
        "prepared path pupil aperture revision does not match its telescope"))
    return nothing
end

@inline _require_path_input_revisions(::Tuple{}, ::UInt) = nothing

@inline function _require_path_input_revisions(inputs::Tuple, revision::UInt)
    _require_path_input_revision(first(inputs), revision)
    return _require_path_input_revisions(Base.tail(inputs), revision)
end

@inline _require_path_input_revisions(input, revision::UInt) =
    _require_path_input_revision(input, revision)

"""Qualified extension seam for cold path execution-binding validation."""
function validate_path_execution_binding(execution, input, result)
    throw(PlantPreparationError(:path, :unsupported_binding_validation,
        "prepared path execution type $(typeof(execution)) does not validate its input/result binding"))
end

struct _PreparedPathExecutorToken end
const _PREPARED_PATH_EXECUTOR_TOKEN = _PreparedPathExecutorToken()

"""
    PreparedPathExecutor

Concrete single-writer owner for one prepared optical path. `input` and
`result` are explicit path-local products, while `execution` owns or references
the concrete prepared propagation/front-end workspace used to update `result`.
"""
struct PreparedPathExecutor{D,S,T,I,R,E,K<:PathResultKey}
    definition::D
    source::S
    telescope::T
    input::I
    result::R
    execution::E
    key::K

    function PreparedPathExecutor(::_PreparedPathExecutorToken,
        definition::D, source::S, telescope::T, input::I, result::R,
        execution::E, key::K) where {D,S,T,I,R,E,K<:PathResultKey}
        return new{D,S,T,I,R,E,K}(definition, source, telescope, input,
            result, execution, key)
    end
end

function PreparedPathExecutor(definition::OpticalPathDefinition,
    source::AbstractSource, telescope::AbstractTelescope, input, result,
    execution;
    optical_model,
    sampling_contract::AbstractOpticalSamplingContract=
        InstantaneousOpticalSample(),
    propagation_model,
    model_revisions=())
    _validate_path_input(input)
    output_plane = _path_output_contract(result)
    first_result = _first_path_result(result)
    selector = backend(first_result)
    device = plane_device(first_result.values)
    _require_path_result_domain(result, selector, device)
    _require_path_input_domain(input, selector, device)
    typeof(backend(telescope)) === typeof(selector) || throw(
        PlantPreparationError(:path, :backend,
            "prepared path and telescope backends differ"))
    plane_device(pupil_reflectivity(telescope)) == device || throw(
        PlantPreparationError(:path, :device,
            "prepared path and telescope occupy different physical devices"))
    revision = aperture_revision(telescope)
    _require_path_input_revisions(input, revision)
    validate_path_execution_binding(execution, input, result)
    revisions = (telescope=revision, model=model_revisions)
    key = PathResultKey(
        path_source_geometry_key(source),
        path_source_spectral_key(source),
        path_source_radiometry_key(source),
        optical_model,
        sampling_contract,
        propagation_model,
        output_plane,
        revisions,
        selector,
        device,
    )
    return PreparedPathExecutor(_PREPARED_PATH_EXECUTOR_TOKEN, definition,
        source, telescope, input, result, execution, key)
end

@inline path_input(path::PreparedPathExecutor) = path.input
@inline path_result(path::PreparedPathExecutor) = path.result
@inline path_result_key(path::PreparedPathExecutor) = path.key

function _require_current_path_revision(path::PreparedPathExecutor)
    revision = getfield(getfield(path.key, :revisions), :telescope)
    aperture_revision(path.telescope) == revision || throw(
        PlantPreparationError(:path, :revision,
            "telescope aperture changed after path preparation"))
    _require_path_input_revisions(path.input, revision)
    return nothing
end

"""Execute one already prepared path without selecting or advancing time."""
function execute_path!(path::PreparedPathExecutor)
    _require_current_path_revision(path)
    return execute_path!(path.result, path.input, path.execution)
end

function execute_path!(result, input, execution)
    throw(PlantPreparationError(:path, :unsupported_execution,
        "prepared path execution type $(typeof(execution)) does not implement execute_path!"))
end

@inline function _require_direct_path_input(
    prepared::PreparedDirectImaging, input)
    prepared.input === input || throw(PlantPreparationError(:path,
        :prepared_binding,
        "direct-imaging input does not match its prepared path"))
    return nothing
end

function _require_direct_path_input(prepared::Union{
    PreparedIncoherentDirectImaging,PreparedBundledDirectImaging}, input)
    @inbounds for component in prepared.components
        component.input === input || throw(PlantPreparationError(:path,
            :prepared_binding,
            "direct-imaging component input does not match its prepared path"))
    end
    return nothing
end

function execute_path!(result, input, prepared::Union{
    PreparedDirectImaging,PreparedIncoherentDirectImaging,
    PreparedBundledDirectImaging})
    validate_path_execution_binding(prepared, input, result)
    return form_direct_image!(prepared)
end

function validate_path_execution_binding(prepared::Union{
    PreparedDirectImaging,PreparedIncoherentDirectImaging,
    PreparedBundledDirectImaging}, input, result)
    direct_imaging_output(prepared) === result || throw(
        PlantPreparationError(:path, :prepared_binding,
            "direct-imaging output does not match its prepared path"))
    _require_direct_path_input(prepared, input)
    return nothing
end

"""Adapter from a Gate 0 prepared WFS optical plan to a plant path."""
struct WFSOpticalPathExecution{P}
    plan::P
end

@inline function validate_path_execution_binding(
    execution::WFSOpticalPathExecution, input, result)
    return validate_wfs_optical_formation_binding(result, input,
        execution.plan)
end

@inline function execute_path!(result, input,
    execution::WFSOpticalPathExecution)
    return form_wfs_optical_products!(result, input, execution.plan)
end

"""Caller-owned products published by one prepared acquisition."""
struct AcquisitionProducts{O,M}
    observation::O
    measurement::M
end

@inline AcquisitionProducts(observation) =
    AcquisitionProducts(observation, nothing)

struct _FrameAcquisitionExecutionToken end
const _FRAME_ACQUISITION_EXECUTION_TOKEN = _FrameAcquisitionExecutionToken()

"""Prepared detector capture and copy into a distinct caller-owned frame."""
struct FrameAcquisitionExecution{D,P,F}
    detector::D
    plan::P
    observation::F

    function FrameAcquisitionExecution(::_FrameAcquisitionExecutionToken,
        detector::D, plan::P, observation::F) where {D,P,F}
        return new{D,P,F}(detector, plan, observation)
    end
end

function FrameAcquisitionExecution(detector::Detector,
    optical_result::IntensityMap, observation::AbstractArray)
    plan = prepare_detector_acquisition(detector, optical_result)
    return _frame_acquisition_execution(detector, plan, observation)
end

function FrameAcquisitionExecution(detector::Detector,
    optical_result::IntensityMap)
    plan = prepare_detector_acquisition(detector, optical_result)
    observation = similar(output_frame(detector))
    fill!(observation, zero(eltype(observation)))
    return _frame_acquisition_execution(detector, plan, observation)
end

function _frame_acquisition_execution(detector::Detector, plan,
    observation::AbstractArray)
    _require_frame_acquisition_observation(detector, observation)
    return FrameAcquisitionExecution(_FRAME_ACQUISITION_EXECUTION_TOKEN,
        detector, plan, observation)
end

function _require_frame_acquisition_observation(detector::Detector,
    observation::AbstractArray)
    frame = output_frame(detector)
    size(observation) == size(frame) || throw(PlantPreparationError(
        :acquisition, :shape,
        "acquisition observation shape must match detector output"))
    eltype(observation) === eltype(frame) || throw(PlantPreparationError(
        :acquisition, :numeric_type,
        "acquisition observation element type must match detector output"))
    typeof(backend(observation)) === typeof(backend(frame)) || throw(
        PlantPreparationError(:acquisition, :backend,
            "acquisition observation and detector output backends differ"))
    plane_device(observation) == plane_device(frame) || throw(
        PlantPreparationError(:acquisition, :device,
            "acquisition observation and detector output occupy different devices"))
    Base.mightalias(observation, frame) && throw(PlantPreparationError(
        :acquisition, :ownership,
        "caller-owned acquisition observation must not alias detector state"))
    return nothing
end

struct _WFSAcquisitionExecutionToken end
const _WFS_ACQUISITION_EXECUTION_TOKEN = _WFSAcquisitionExecutionToken()

"""Prepared Gate 0 WFS acquisition and optional estimator composition."""
struct WFSAcquisitionExecution{A,E,O,M}
    acquisition::A
    estimator::E
    observation::O
    measurement::M

    function WFSAcquisitionExecution(::_WFSAcquisitionExecutionToken,
        acquisition::A, estimator::E, observation::O,
        measurement::M) where {A,E,O,M}
        return new{A,E,O,M}(acquisition, estimator, observation, measurement)
    end
end

@inline WFSAcquisitionExecution(acquisition, observation) =
    WFSAcquisitionExecution(_WFS_ACQUISITION_EXECUTION_TOKEN, acquisition,
        nothing, observation, nothing)

function WFSAcquisitionExecution(acquisition, estimator,
    observation, measurement::WFSMeasurement)
    return WFSAcquisitionExecution(_WFS_ACQUISITION_EXECUTION_TOKEN,
        acquisition, estimator, observation, measurement)
end

"""Qualified extension seam for cold acquisition execution validation."""
function validate_acquisition_execution_binding(execution, path_result,
    products::AcquisitionProducts)
    throw(PlantPreparationError(:acquisition,
        :unsupported_binding_validation,
        "prepared acquisition execution type $(typeof(execution)) does not validate its path/product binding"))
end

function validate_acquisition_execution_binding(
    execution::FrameAcquisitionExecution,
    path_result::IntensityMap,
    products::AcquisitionProducts{<:AbstractArray,Nothing})
    products.observation === execution.observation || throw(
        PlantPreparationError(:acquisition, :prepared_binding,
            "acquisition observation does not match its prepared storage"))
    _require_frame_acquisition_observation(execution.detector,
        execution.observation)
    plan = execution.plan
    path_result.metadata === plan.input_metadata &&
        path_result.values === plan.input_values || throw(
        PlantPreparationError(:acquisition, :prepared_binding,
            "acquisition path result does not match its detector plan"))
    execution.detector.params === plan.detector_params &&
        execution.detector.state === plan.detector_state &&
        execution.detector.state.frame === plan.detector_frame || throw(
        PlantPreparationError(:acquisition, :prepared_binding,
            "acquisition detector storage changed after preparation"))
    return nothing
end

@inline _require_acquired_wfs_estimator(::AcquiredObservationPath) = nothing

function _require_acquired_wfs_estimator(::DirectMeasurementPath)
    throw(PlantPreparationError(:acquisition, :estimator,
        "WFS acquisition composition requires an acquired-observation estimator"))
end

function _require_acquired_wfs_estimator(path)
    throw(PlantPreparationError(:acquisition, :estimator,
        "WFS estimator declared unsupported measurement path $(typeof(path))"))
end

function validate_acquisition_execution_binding(
    execution::WFSAcquisitionExecution{A,Nothing,O,Nothing}, path_result,
    products::AcquisitionProducts{O,Nothing}) where {A,O}
    products.observation === execution.observation || throw(
        PlantPreparationError(:acquisition, :prepared_binding,
            "WFS observation does not match its prepared acquisition"))
    validate_wfs_acquisition_binding(products.observation, path_result,
        execution.acquisition)
    return nothing
end

function validate_acquisition_execution_binding(
    execution::WFSAcquisitionExecution{A,E,O,M}, path_result,
    products::AcquisitionProducts{O,M}) where {A,E,O,M<:WFSMeasurement}
    products.observation === execution.observation &&
        products.measurement === execution.measurement || throw(
        PlantPreparationError(:acquisition, :prepared_binding,
            "WFS products do not match their prepared acquisition"))
    validate_wfs_acquisition_binding(products.observation, path_result,
        execution.acquisition)
    _require_acquired_wfs_estimator(wfs_measurement_path(execution.estimator))
    validate_wfs_estimation_binding(products.measurement,
        products.observation, execution.estimator)
    return nothing
end

struct _PreparedAcquisitionOwnerToken end
const _PREPARED_ACQUISITION_OWNER_TOKEN = _PreparedAcquisitionOwnerToken()

"""
    PreparedAcquisitionOwner

Concrete single-writer owner for detector/readout and optional WFS estimator
state. It borrows one exact prepared path result as read-only input and owns
separate caller-visible observation/measurement products.
"""
struct PreparedAcquisitionOwner{D,K<:PathResultKey,R,E,P}
    definition::D
    path_key::K
    path_result::R
    execution::E
    products::P

    function PreparedAcquisitionOwner(::_PreparedAcquisitionOwnerToken,
        definition::D, path_key::K, path_result::R, execution::E,
        products::P) where {D,K<:PathResultKey,R,E,P}
        return new{D,K,R,E,P}(definition, path_key, path_result, execution,
            products)
    end
end

@inline acquisition_products(owner::PreparedAcquisitionOwner) = owner.products
@inline acquisition_observation(owner::PreparedAcquisitionOwner) =
    owner.products.observation
@inline acquisition_measurement(owner::PreparedAcquisitionOwner) =
    owner.products.measurement

function _path_key_mismatch_reason(actual::PathResultKey,
    expected::PathResultKey)
    !isequal(getfield(actual, :source_geometry),
        getfield(expected, :source_geometry)) &&
        return :source_geometry
    !isequal(getfield(actual, :spectral_sampling),
        getfield(expected, :spectral_sampling)) &&
        return :spectral_sampling
    !isequal(getfield(actual, :radiometry),
        getfield(expected, :radiometry)) && return :radiometry
    !isequal(getfield(actual, :optical_model),
        getfield(expected, :optical_model)) &&
        return :optical_model
    !isequal(getfield(actual, :sampling_contract),
        getfield(expected, :sampling_contract)) &&
        return :sampling_contract
    !isequal(getfield(actual, :propagation_model),
        getfield(expected, :propagation_model)) &&
        return :propagation_model
    !isequal(getfield(actual, :output_plane),
        getfield(expected, :output_plane)) &&
        return :output_plane
    !isequal(getfield(actual, :revisions), getfield(expected, :revisions)) &&
        return :revision
    typeof(getfield(actual, :backend)) ===
        typeof(getfield(expected, :backend)) || return :backend
    !isequal(getfield(actual, :device), getfield(expected, :device)) &&
        return :device
    return :compatibility
end

function _require_path_result_key(actual::PathResultKey,
    expected::PathResultKey)
    isequal(actual, expected) && return actual
    reason = _path_key_mismatch_reason(actual, expected)
    throw(PlantPreparationError(:acquisition, reason,
        "acquisition requires a path result incompatible in $(reason)"))
end

"""
    require_path_result(path; ...)

Validate an acquisition's cold compatibility requirements against a prepared
path before constructing or mutating acquisition destinations. Omitted fields
retain the prepared path's value. This is a preparation-time extension seam,
not a hot-path lookup.
"""
function require_path_result(path::PreparedPathExecutor;
    source_geometry=getfield(path.key, :source_geometry),
    spectral_sampling=getfield(path.key, :spectral_sampling),
    radiometry=getfield(path.key, :radiometry),
    optical_model=getfield(path.key, :optical_model),
    sampling_contract=getfield(path.key, :sampling_contract),
    propagation_model=getfield(path.key, :propagation_model),
    output_plane=getfield(path.key, :output_plane),
    revisions=getfield(path.key, :revisions),
    backend::AbstractArrayBackend=getfield(path.key, :backend),
    device::AbstractPlaneDevice=getfield(path.key, :device))
    required = PathResultKey(source_geometry, spectral_sampling, radiometry,
        optical_model, sampling_contract, propagation_model, output_plane,
        revisions, backend, device)
    _require_path_result_key(path.key, required)
    return path
end

function PreparedAcquisitionOwner(definition::AcquisitionDefinition,
    path::PreparedPathExecutor, execution, products::AcquisitionProducts)
    acquisition_path_id(definition) == path_id(path.definition) || throw(
        PlantPreparationError(:acquisition, :unknown_path,
            "acquisition $(definition.id) does not reference prepared path $(path.definition.id)"))
    validate_acquisition_execution_binding(execution, path.result, products)
    return PreparedAcquisitionOwner(_PREPARED_ACQUISITION_OWNER_TOKEN,
        definition, path.key, path.result, execution, products)
end

"""Execute one acquisition from its already formed path result."""
@inline function execute_acquisition!(owner::PreparedAcquisitionOwner, rng)
    return execute_acquisition!(owner.products, owner.path_result,
        owner.execution, rng)
end

function execute_acquisition!(products, path_result, execution, rng)
    throw(PlantPreparationError(:acquisition, :unsupported_execution,
        "prepared acquisition execution type $(typeof(execution)) does not implement execute_acquisition!"))
end

function execute_acquisition!(products::AcquisitionProducts{<:AbstractArray,
    Nothing}, path_result::IntensityMap,
    execution::FrameAcquisitionExecution, rng::AbstractRNG)
    validate_acquisition_execution_binding(execution, path_result, products)
    frame = capture!(execution.detector, path_result, execution.plan, rng)
    copyto!(products.observation, frame)
    return products
end

function execute_acquisition!(products::AcquisitionProducts{O,Nothing},
    path_result,
    execution::WFSAcquisitionExecution{A,Nothing,O,Nothing}, rng) where {A,O}
    validate_acquisition_execution_binding(execution, path_result, products)
    acquire_wfs_observation!(products.observation, path_result,
        execution.acquisition, rng)
    return products
end

function execute_acquisition!(products::AcquisitionProducts{O,M},
    path_result,
    execution::WFSAcquisitionExecution{A,E,O,M}, rng) where {
    A,E,O,M<:WFSMeasurement,
}
    validate_acquisition_execution_binding(execution, path_result, products)
    acquire_wfs_observation!(products.observation, path_result,
        execution.acquisition, rng)
    estimate_wfs_measurement!(products.measurement, products.observation,
        execution.estimator)
    return products
end

struct _PreparedPlantToken end
const _PREPARED_PLANT_TOKEN = _PreparedPlantToken()

"""Prepared, schedule-free plant with concrete path and acquisition tuples."""
struct PreparedPlant{D,P<:Tuple,A<:Tuple}
    definition::D
    paths::P
    acquisitions::A

    function PreparedPlant(::_PreparedPlantToken, definition::D, paths::P,
        acquisitions::A) where {D,P<:Tuple,A<:Tuple}
        return new{D,P,A}(definition, paths, acquisitions)
    end
end

@inline prepared_paths(plant::PreparedPlant) = plant.paths
@inline prepared_acquisitions(plant::PreparedPlant) = plant.acquisitions

function prepared_path(plant::PreparedPlant, id)
    resolved = _as_optical_path_id(id)
    for path in plant.paths
        path_id(path.definition) == resolved && return path
    end
    throw(PlantPreparationError(:path, :unknown_id,
        "prepared plant has no optical path $resolved"))
end

function prepared_acquisition(plant::PreparedPlant, id)
    resolved = _as_acquisition_id(id)
    for acquisition in plant.acquisitions
        acquisition_id(acquisition.definition) == resolved &&
            return acquisition
    end
    throw(PlantPreparationError(:acquisition, :unknown_id,
        "prepared plant has no acquisition $resolved"))
end

function prepare_path_executor(definition::OpticalPathDefinition,
    telescope::AbstractTelescope, atmosphere::AbstractAtmosphere)
    source = freeze_source(path_source(definition))
    prepared = prepare_path_executor(path_model(definition), definition,
        source, telescope, atmosphere)
    return _require_prepared_path_executor(prepared, definition, source,
        telescope)
end

function _require_prepared_path_executor(prepared::PreparedPathExecutor,
    definition::OpticalPathDefinition, source::AbstractSource,
    telescope::AbstractTelescope)
    prepared.definition === definition || throw(PlantPreparationError(
        :path, :prepared_binding,
        "prepared path does not retain its exact definition"))
    prepared.source === source || throw(PlantPreparationError(:path,
        :prepared_binding,
        "prepared path does not retain its run-owned frozen source"))
    prepared.telescope === telescope || throw(PlantPreparationError(:path,
        :prepared_binding,
        "prepared path does not retain its plant telescope"))
    validate_path_execution_binding(prepared.execution, prepared.input,
        prepared.result)
    return prepared
end


function _require_prepared_path_executor(prepared,
    ::OpticalPathDefinition, ::AbstractSource, ::AbstractTelescope)
    throw(PlantPreparationError(
        :path, :invalid_preparation,
        "path model preparation must return PreparedPathExecutor; got $(typeof(prepared))"))
end

function prepare_path_executor(model, definition::OpticalPathDefinition,
    source::AbstractSource, telescope::AbstractTelescope,
    atmosphere::AbstractAtmosphere)
    throw(PlantPreparationError(:path, :unsupported_model,
        "path model $(typeof(model)) does not implement prepare_path_executor"))
end

function prepare_acquisition_owner(definition::AcquisitionDefinition,
    path::PreparedPathExecutor)
    prepared = prepare_acquisition_owner(acquisition_model(definition),
        definition, path)
    return _require_prepared_acquisition_owner(prepared, definition, path)
end

function _require_prepared_acquisition_owner(
    prepared::PreparedAcquisitionOwner,
    definition::AcquisitionDefinition, path::PreparedPathExecutor)
    prepared.definition === definition || throw(PlantPreparationError(
        :acquisition, :prepared_binding,
        "prepared acquisition does not retain its exact definition"))
    prepared.path_key === path.key || throw(PlantPreparationError(
        :acquisition, :prepared_binding,
        "prepared acquisition does not retain the exact path-result key"))
    prepared.path_result === path.result || throw(PlantPreparationError(
        :acquisition, :prepared_binding,
        "prepared acquisition does not retain the exact path result"))
    validate_acquisition_execution_binding(prepared.execution,
        prepared.path_result, prepared.products)
    return prepared
end

function _require_prepared_acquisition_owner(prepared,
    ::AcquisitionDefinition, ::PreparedPathExecutor)
    throw(PlantPreparationError(
        :acquisition, :invalid_preparation,
        "acquisition model preparation must return PreparedAcquisitionOwner; got $(typeof(prepared))"))
end

function prepare_acquisition_owner(model,
    definition::AcquisitionDefinition, path::PreparedPathExecutor)
    throw(PlantPreparationError(:acquisition, :unsupported_model,
        "acquisition model $(typeof(model)) does not implement prepare_acquisition_owner"))
end

@inline _prepare_path_executors(::Tuple{}, telescope, atmosphere) = ()

@inline function _prepare_path_executors(definitions::Tuple, telescope,
    atmosphere)
    return (
        prepare_path_executor(first(definitions), telescope, atmosphere),
        _prepare_path_executors(Base.tail(definitions), telescope,
            atmosphere)...,
    )
end

@inline _prepare_acquisition_owners(::Tuple{}, paths) = ()

@inline function _prepared_path_for_acquisition(definition, paths::Tuple)
    id = acquisition_path_id(definition)
    for path in paths
        path_id(path.definition) == id && return path
    end
    throw(PlantPreparationError(:acquisition, :unknown_path,
        "acquisition $(definition.id) references an unprepared path $id"))
end

@inline function _prepare_acquisition_owners(definitions::Tuple, paths)
    definition = first(definitions)
    path = _prepared_path_for_acquisition(definition, paths)
    return (
        prepare_acquisition_owner(definition, path),
        _prepare_acquisition_owners(Base.tail(definitions), paths)...,
    )
end

"""
    prepare_plant(definition)

Prepare all declared paths and acquisitions without scheduling or executing
them. Model-specific construction dispatches on the cold model-definition
types. Preparation may allocate and perform fallible backend/device/revision
validation; repeated execution uses the concrete tuples stored in the result.
"""
function prepare_plant(definition::PlantDefinition)
    paths = _prepare_path_executors(path_definitions(definition),
        plant_telescope(definition), plant_atmosphere(definition))
    acquisitions = _prepare_acquisition_owners(
        acquisition_definitions(definition), paths)
    return PreparedPlant(_PREPARED_PLANT_TOKEN, definition, paths,
        acquisitions)
end

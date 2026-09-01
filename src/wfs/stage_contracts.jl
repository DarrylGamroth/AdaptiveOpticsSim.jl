#
# Prepared wavefront-sensor stage contracts
#
# This is deliberately a semantic protocol, not a universal optical graph.
# Concrete front ends, acquisitions, and estimators own their immutable plans
# and single-writer workspaces. Products and destinations remain caller-owned.
#

"""Immutable storage description for one acquired WFS observation."""
struct WFSObservationMetadata{
    N,
    E,
    L,
    B<:AbstractArrayBackend,
    D<:AbstractComputeDevice,
}
    dimensions::NTuple{N,Int}
    numeric_type::Type{E}
    layout::L
    backend::B
    device::D
end

"""Immutable storage description for one typed WFS measurement."""
struct WFSMeasurementMetadata{
    N,
    E,
    K,
    B<:AbstractArrayBackend,
    D<:AbstractComputeDevice,
}
    dimensions::NTuple{N,Int}
    numeric_type::Type{E}
    kind::K
    backend::B
    device::D
end

@inline _wfs_storage_dimensions(storage::AbstractArray) = size(storage)
@inline _wfs_storage_dimensions(::Base.RefValue) = ()
@inline _wfs_storage_numeric_type(storage::AbstractArray) = eltype(storage)
@inline _wfs_storage_numeric_type(::Base.RefValue{T}) where {T} = T
@inline _wfs_storage_backend(storage::AbstractArray) = backend(storage)
@inline _wfs_storage_backend(::Base.RefValue) = CPUBackend()
@inline _wfs_storage_device(storage::AbstractArray) = compute_device(storage)
@inline _wfs_storage_device(::Base.RefValue) = HostComputeDevice()
@inline _wfs_storage_length(storage::AbstractArray) = length(storage)
@inline _wfs_storage_length(::Base.RefValue) = 1

function _require_wfs_storage(storage, stage::Symbol)
    throw(WFSPreparationError(stage, :shape,
        "WFS product storage must be an AbstractArray or host Ref"))
end

function _require_wfs_storage(storage::Union{AbstractArray,Base.RefValue},
    stage::Symbol)
    _wfs_storage_length(storage) > 0 || throw(WFSPreparationError(stage,
        :shape, "WFS product storage must not be empty"))
    numeric_type = _wfs_storage_numeric_type(storage)
    isconcretetype(numeric_type) && numeric_type <: Number ||
        throw(WFSPreparationError(stage, :numeric_type,
            "WFS product storage must have a concrete numeric element type"))
    return storage
end

function WFSObservationMetadata(storage; layout=:dense)
    _require_wfs_storage(storage, :acquisition)
    _require_declared_wfs_descriptor(layout, :acquisition, :detector_mapping,
        "observation layout")
    dimensions = _wfs_storage_dimensions(storage)
    E = _wfs_storage_numeric_type(storage)
    selector = _wfs_storage_backend(storage)
    device = _wfs_storage_device(storage)
    return WFSObservationMetadata{
        length(dimensions),E,typeof(layout),typeof(selector),typeof(device),
    }(dimensions, E, layout, selector, device)
end

function WFSMeasurementMetadata(storage; kind)
    _require_wfs_storage(storage, :estimation)
    _require_declared_wfs_descriptor(kind, :estimation, :estimator,
        "measurement kind")
    dimensions = _wfs_storage_dimensions(storage)
    E = _wfs_storage_numeric_type(storage)
    selector = _wfs_storage_backend(storage)
    device = _wfs_storage_device(storage)
    return WFSMeasurementMetadata{
        length(dimensions),E,typeof(kind),typeof(selector),typeof(device),
    }(dimensions, E, kind, selector, device)
end

"""
    WFSObservation(storage; units, layout=:dense)
    WFSObservation(storage, units, metadata)

Wrap caller-owned, preallocated observation storage. `storage` may be a
host `Ref` for a scalar or an array of any rank. Multiple observations use a
concrete tuple of wrappers. `units` is required and may be a `Symbol` or an
application-defined singleton; no unit package is imposed by core.
"""
struct WFSObservation{S,U,M<:WFSObservationMetadata}
    storage::S
    units::U
    metadata::M

    function WFSObservation(storage::S, units::U, metadata::M) where {
        S,U,M<:WFSObservationMetadata,
    }
        _require_declared_wfs_units(units, :acquisition)
        observation = new{S,U,M}(storage, units, metadata)
        validate_wfs_observation(observation)
        return observation
    end
end

function WFSObservation(storage; units, layout=:dense)
    metadata = WFSObservationMetadata(storage; layout=layout)
    return WFSObservation(storage, units, metadata)
end

"""
    WFSMeasurement(storage; units, kind)
    WFSMeasurement(storage, units, metadata)

Wrap caller-owned, preallocated estimator output. `kind` declares its semantic
quantity (for example `:centroid_slopes`, `:phase`, or an application-owned
singleton); `units` declares the output units independently of its shape.
"""
struct WFSMeasurement{S,U,M<:WFSMeasurementMetadata}
    storage::S
    units::U
    metadata::M

    function WFSMeasurement(storage::S, units::U, metadata::M) where {
        S,U,M<:WFSMeasurementMetadata,
    }
        _require_declared_wfs_units(units, :estimation)
        measurement = new{S,U,M}(storage, units, metadata)
        validate_wfs_measurement(measurement)
        return measurement
    end
end

function WFSMeasurement(storage; units, kind)
    metadata = WFSMeasurementMetadata(storage; kind=kind)
    return WFSMeasurement(storage, units, metadata)
end

@inline observation_storage(observation::WFSObservation) = observation.storage
@inline observation_units(observation::WFSObservation) = observation.units
@inline observation_metadata(observation::WFSObservation) = observation.metadata
@inline measurement_storage(measurement::WFSMeasurement) = measurement.storage
@inline measurement_units(measurement::WFSMeasurement) = measurement.units
@inline measurement_metadata(measurement::WFSMeasurement) = measurement.metadata

@inline _require_declared_wfs_units(units, ::Symbol) = units

@inline _require_declared_wfs_descriptor(value, ::Symbol, ::Symbol,
    ::AbstractString) = value

function _require_declared_wfs_units(::Nothing, stage::Symbol)
    throw(WFSPreparationError(stage, :units,
        "WFS product units must be declared"))
end

function _require_declared_wfs_units(units::Symbol, stage::Symbol)
    isempty(String(units)) && throw(WFSPreparationError(stage, :units,
        "WFS product unit symbol must not be empty"))
    return units
end

function _require_declared_wfs_descriptor(::Nothing, stage::Symbol,
    reason::Symbol, label::AbstractString)
    throw(WFSPreparationError(stage, reason,
        "WFS $label must be declared"))
end

function _require_declared_wfs_descriptor(value::Symbol, stage::Symbol,
    reason::Symbol, label::AbstractString)
    isempty(String(value)) && throw(WFSPreparationError(stage, reason,
        "WFS $label symbol must not be empty"))
    return value
end

function _validate_wfs_storage(metadata, storage, stage::Symbol)
    _require_wfs_storage(storage, stage)
    _wfs_storage_dimensions(storage) == metadata.dimensions ||
        throw(WFSPreparationError(stage, :shape,
            "WFS product storage dimensions do not match prepared metadata"))
    _wfs_storage_numeric_type(storage) === metadata.numeric_type ||
        throw(WFSPreparationError(stage, :shape,
            "WFS product storage element type does not match prepared metadata"))
    typeof(_wfs_storage_backend(storage)) === typeof(metadata.backend) ||
        throw(WFSPreparationError(stage, :backend,
            "WFS product storage backend does not match prepared metadata"))
    _wfs_storage_device(storage) == metadata.device ||
        throw(WFSPreparationError(stage, :device,
            "WFS product storage device does not match prepared metadata"))
    return metadata
end

function _require_wfs_storage_domain(stage::Symbol, metadata, storage,
    label::AbstractString)
    typeof(metadata.backend) === typeof(backend(storage)) ||
        throw(WFSPreparationError(stage, :backend,
            "$label backend does not match the prepared WFS stage"))
    metadata.device == compute_device(storage) ||
        throw(WFSPreparationError(stage, :device,
            "$label device does not match the prepared WFS stage"))
    return nothing
end

function _require_real_square_wfs_observation(observation::WFSObservation,
    label::AbstractString)
    observation.metadata.numeric_type <: Real ||
        throw(WFSPreparationError(:estimation, :numeric_type,
            "$label observations require real detector samples"))
    dimensions = observation.metadata.dimensions
    length(dimensions) == 2 || throw(WFSPreparationError(:estimation, :shape,
        "$label observations require a two-dimensional detector frame"))
    dimensions[1] == dimensions[2] || throw(WFSPreparationError(
        :estimation, :shape,
        "$label observations require a square detector frame"))
    return dimensions[1]
end

function validate_wfs_observation(observation::WFSObservation)
    _require_declared_wfs_units(observation.units, :acquisition)
    _require_declared_wfs_descriptor(observation.metadata.layout,
        :acquisition, :detector_mapping, "observation layout")
    _validate_wfs_storage(observation.metadata, observation.storage,
        :acquisition)
    return observation
end

function validate_wfs_measurement(measurement::WFSMeasurement)
    _require_declared_wfs_units(measurement.units, :estimation)
    _require_declared_wfs_descriptor(measurement.metadata.kind, :estimation,
        :estimator, "measurement kind")
    _validate_wfs_storage(measurement.metadata, measurement.storage,
        :estimation)
    return measurement
end

@inline validate_wfs_observations(::Tuple{}) = throw(WFSPreparationError(
    :acquisition, :plane_count,
    "a prepared WFS acquisition requires at least one observation"))

@inline function validate_wfs_observations(observations::Tuple)
    validate_wfs_observation(first(observations))
    _validate_remaining_wfs_observations(Base.tail(observations))
    return observations
end

@inline _validate_remaining_wfs_observations(::Tuple{}) = nothing

@inline function _validate_remaining_wfs_observations(observations::Tuple)
    validate_wfs_observation(first(observations))
    return _validate_remaining_wfs_observations(Base.tail(observations))
end

@inline validate_wfs_observations(observation::WFSObservation) =
    validate_wfs_observation(observation)

@inline function _require_wfs_input_plane(::PupilPlane)
    return nothing
end

function _require_wfs_input_plane(::AbstractOpticalPlaneKind)
    throw(WFSPreparationError(:wfs_optics, :plane_metadata,
        "WFS optical input must be declared on a pupil plane"))
end

function validate_wfs_optical_input(input::PupilFunction)
    _require_wfs_input_plane(input.metadata.kind)
    validate_plane_storage(input.metadata, input.amplitude;
        label="WFS pupil amplitude")
    validate_plane_storage(input.metadata, input.opd;
        label="WFS pupil OPD")
    return input
end

function validate_wfs_optical_input(input::ElectricField)
    _require_wfs_input_plane(input.metadata.kind)
    validate_plane_storage(input.metadata, input.values;
        label="WFS pupil electric field")
    return input
end

function validate_wfs_optical_input(input)
    throw(WFSPreparationError(:wfs_optics, :plane_metadata,
        "WFS optics require a PupilFunction or pupil-plane ElectricField"))
end

@inline _require_wfs_output_plane(::DetectorPlane) = nothing

function _require_wfs_output_plane(::AbstractOpticalPlaneKind)
    throw(WFSPreparationError(:wfs_optics, :plane_metadata,
        "WFS optical output must be declared on a detector plane"))
end

@inline _require_wfs_rate_normalization(::PhotonRateNormalization) = nothing

function _require_wfs_rate_normalization(::AbstractOpticalNormalization)
    throw(WFSPreparationError(:wfs_optics, :radiometry,
        "WFS detector-facing optical output must use photon-rate normalization"))
end

@inline _require_wfs_rate_measure(::SpatialDensityMeasure) = nothing
@inline _require_wfs_rate_measure(::CellIntegratedMeasure) = nothing

function _require_wfs_rate_measure(::AbstractSpatialMeasure)
    throw(WFSPreparationError(:wfs_optics, :radiometry,
        "WFS detector-facing optical output must declare spatial-density or cell-integrated measure"))
end

@inline _require_wfs_rate_coherence(::IncoherentIntensityAddition) = nothing

function _require_wfs_rate_coherence(::AbstractCombinationPolicy)
    throw(WFSPreparationError(:wfs_optics, :radiometry,
        "WFS detector-facing optical output must declare incoherent-intensity semantics"))
end

@inline _require_wfs_rate_spectral(::MonochromaticChannel) = nothing
@inline _require_wfs_rate_spectral(::IntegratedSpectralChannel) = nothing

function _require_wfs_rate_spectral(::AbstractSpectralCoordinate)
    throw(WFSPreparationError(:wfs_optics, :plane_metadata,
        "WFS detector-facing optical output requires a monochromatic or integrated spectral channel"))
end

function validate_wfs_optical_products(product::IntensityMap)
    validate_plane_storage(product.metadata, product.values;
        label="WFS detector-facing intensity map")
    _require_wfs_output_plane(product.metadata.kind)
    _require_wfs_rate_normalization(product.metadata.normalization)
    _require_wfs_rate_measure(product.metadata.spatial_measure)
    _require_wfs_rate_coherence(product.metadata.coherence)
    _require_wfs_rate_spectral(product.metadata.spectral)
    return product
end

@inline function validate_wfs_optical_products(bundle::OpticalProductBundle)
    validate_wfs_optical_products(bundle.products)
    return bundle
end

function validate_wfs_optical_products(
    products::AbstractVector{<:AbstractOpticalProduct})
    isempty(products) && throw(WFSPreparationError(:wfs_optics,
        :plane_count,
        "WFS optics require at least one detector-facing product"))
    @inbounds for product in products
        validate_wfs_optical_products(product)
    end
    return products
end

@inline validate_wfs_optical_products(::Tuple{}) = throw(WFSPreparationError(
    :wfs_optics, :plane_count,
    "WFS optics require at least one detector-facing product"))

@inline function validate_wfs_optical_products(products::Tuple)
    validate_wfs_optical_products(first(products))
    _validate_remaining_wfs_optical_products(Base.tail(products))
    return products
end

@inline _validate_remaining_wfs_optical_products(::Tuple{}) = nothing

@inline function _validate_remaining_wfs_optical_products(products::Tuple)
    validate_wfs_optical_products(first(products))
    return _validate_remaining_wfs_optical_products(Base.tail(products))
end

function validate_wfs_optical_products(product)
    throw(WFSPreparationError(:wfs_optics, :plane_metadata,
        "WFS optical products must be IntensityMap values or a concrete tuple/bundle of them"))
end

"""
    AbstractWFSOpticsPlan

Nominal interface for a run-immutable WFS optics contract. Concrete
plans retain validated physical and numerical data only; the corresponding
prepared owner binds exact input, workspace, products, backend, and device.
Implementations provide the existing `prepare_wfs_optics`,
`form_wfs_optical_products!`, `wfs_optical_products`, and
`validate_wfs_optics_binding` protocol.
"""
abstract type AbstractWFSOpticsPlan end

"""
    AbstractWFSAcquisitionPlan

Nominal interface for a run-immutable WFS acquisition contract. Concrete
plans contain validated detector/radiometric and observation metadata, while
prepared owners bind exact detector state, workspace, products, optical input,
observation, backend, and device. Obtain the plan from a generic prepared
owner with [`wfs_acquisition_plan`](@ref).
"""
abstract type AbstractWFSAcquisitionPlan end

"""
    AbstractWFSEstimationPlan

Nominal interface for a run-immutable WFS estimation contract. Concrete plans
retain validated calibration and estimator data only; the corresponding
prepared owner binds exact observation or direct optical input, workspace,
measurement product, backend, and device. Implementations provide the existing
`prepare_wfs_estimation`, `estimate_wfs_measurement!`, and
`validate_wfs_estimation_binding` protocol.
"""
abstract type AbstractWFSEstimationPlan end

"""
    wfs_acquisition_plan(prepared)

Return the run-immutable acquisition plan bound by a generic prepared WFS
acquisition owner.
"""
function wfs_acquisition_plan end

"""
    prepare_wfs_optics(model, input, output)
    form_wfs_optical_products!(output, input, prepared)

Preparation and mutating execution protocol for a WFS optical front end.
Implementations consume an explicit `PupilFunction` or pupil-plane
`ElectricField` and write caller-owned detector-plane photon-arrival-rate
products. Preparation performs all fallible structural validation.
"""
function prepare_wfs_optics(model, input, output)
    throw(WFSPreparationError(:wfs_optics, :unsupported,
        "$(typeof(model)) does not implement prepared WFS optics"))
end

function form_wfs_optical_products! end

"""
    enqueue_wfs_optical_products!(output, input, prepared)

Internal prepared-runtime seam for a fixed-address accelerator owner. Methods
enqueue the complete optical-product operation on the current device stream
without synchronizing. Structural admission remains the responsibility of
`prepare_wfs_optics`; callers must synchronize before reading or reusing any
bound storage.
"""
function enqueue_wfs_optical_products! end

"""
    wfs_optical_products(prepared)

Return the caller-visible detector-facing photon-arrival-rate product or
products bound by an exact prepared WFS optics owner. The result is an
`IntensityMap`, `OpticalProductBundle`, or family-specific concrete tuple of
rate products; it is not an acquired [`WFSObservation`](@ref) or estimated
[`WFSMeasurement`](@ref).
"""
function wfs_optical_products end

"""
    validate_wfs_optics_binding(output, input, prepared)

Validate without mutation that a prepared WFS optics plan is bound to the
exact input and output supplied by a containing prepared plant path. WFS
extensions that provide a prepared WFS optics plan must implement this qualified
extension seam alongside `form_wfs_optical_products!`.
"""
function validate_wfs_optics_binding(output, input, prepared)
    throw(WFSPreparationError(:wfs_optics,
        :unsupported_binding_validation,
        "$(typeof(prepared)) does not validate its prepared WFS optics binding"))
end

"""
    prepare_wfs_acquisition(model, optical_products, observation)
    acquire_wfs_observation!(observation, optical_products, prepared, rng)

Preparation and mutating execution protocol for detector acquisition. The
prepared owner binds detector state, workspace, products, and explicit
durations; optical inputs remain photon-arrival rates. Concrete tuples express
fan-out to multiple detectors without a runtime-selected stage graph.
"""
function prepare_wfs_acquisition(model, optical_products, observation)
    throw(WFSPreparationError(:acquisition, :unsupported,
        "$(typeof(model)) does not implement prepared WFS acquisition"))
end

function _prepare_wfs_acquisition_candidate(
    model, optical_products, observation; source=nothing)
    throw(WFSPreparationError(:acquisition, :unsupported,
        "$(typeof(model)) does not implement transactional WFS acquisition preparation"))
end

function acquire_wfs_observation! end

"""
    validate_wfs_acquisition_binding(observation, optical_products, prepared)

Validate without mutation that a prepared WFS acquisition is bound to the
exact optical products and observation supplied by a containing plant owner.
"""
function validate_wfs_acquisition_binding(observation, optical_products,
    prepared)
    throw(WFSPreparationError(:acquisition,
        :unsupported_binding_validation,
        "$(typeof(prepared)) does not validate its prepared acquisition binding"))
end

struct _WFSDetectorAcquisitionPlanToken end
const _WFS_DETECTOR_ACQUISITION_PLAN_TOKEN =
    _WFSDetectorAcquisitionPlanToken()

"""Run-immutable WFS contract layered over a frame-detector plan."""
struct WFSDetectorAcquisitionPlan{
    P<:DetectorAcquisitionPlan,
    M<:WFSObservationMetadata,
} <: AbstractWFSAcquisitionPlan
    detector_plan::P
    observation_metadata::M

    function WFSDetectorAcquisitionPlan(
        ::_WFSDetectorAcquisitionPlanToken, detector_plan::P,
        observation_metadata::M) where {
        P<:DetectorAcquisitionPlan,M<:WFSObservationMetadata,
    }
        return new{P,M}(detector_plan, observation_metadata)
    end
end

struct _PreparedWFSDetectorAcquisitionToken end
const _PREPARED_WFS_DETECTOR_ACQUISITION_TOKEN =
    _PreparedWFSDetectorAcquisitionToken()

"""Exact prepared owner for one frame-detector WFS acquisition."""
struct PreparedWFSDetectorAcquisition{
    P<:WFSDetectorAcquisitionPlan,
    A<:PreparedDetectorAcquisition,
    O<:WFSObservation,
}
    plan::P
    acquisition::A
    observation::O

    function PreparedWFSDetectorAcquisition(
        ::_PreparedWFSDetectorAcquisitionToken, plan::P,
        acquisition::A, observation::O) where {
        P<:WFSDetectorAcquisitionPlan,A<:PreparedDetectorAcquisition,
        O<:WFSObservation,
    }
        return new{P,A,O}(plan, acquisition, observation)
    end
end

@inline wfs_acquisition_plan(
    prepared::PreparedWFSDetectorAcquisition) = prepared.plan

@inline function _wfs_storage_mightalias(
    first::AbstractArray, second::AbstractArray)
    return Base.mightalias(first, second)
end

@inline _wfs_storage_mightalias(first, second) = false

function _require_wfs_detector_observation_ownership(
    detector::Detector, optical_product::IntensityMap,
    observation::WFSObservation)
    storage = observation.storage
    (_detector_storage_mightalias(detector, storage) ||
        _wfs_storage_mightalias(storage, optical_product.values)) && throw(
        WFSPreparationError(:acquisition, :ownership,
            "WFS observation must not alias detector storage or its optical input"))
    return nothing
end

function _validate_wfs_detector_observation_contract(
    detector::Detector, optical_product::IntensityMap,
    observation::WFSObservation)
    output_shape = detector_output_shape(
        detector, size(optical_product.values))
    size(observation.storage) == output_shape ||
        throw(WFSPreparationError(:acquisition, :shape,
            "WFS observation storage must match the prepared detector output"))
    output_type = detector_output_type(detector)
    detector_output = output_frame(detector)
    expected_type = output_type === nothing ?
        eltype(detector_output) : output_type
    observation.metadata.numeric_type === expected_type ||
        throw(WFSPreparationError(:acquisition, :numeric_type,
            "WFS observation element type must match the detector output"))
    typeof(backend(observation.storage)) === typeof(backend(detector)) ||
        throw(WFSPreparationError(:acquisition, :backend,
            "WFS observation and detector backends differ"))
    compute_device(observation.storage) ==
        compute_device(detector_output) ||
        throw(WFSPreparationError(:acquisition, :device,
            "WFS observation and detector output occupy different devices"))
    return nothing
end

function _prepare_wfs_acquisition_candidate(detector::Detector,
    optical_product::IntensityMap, observation::WFSObservation;
    source=nothing)
    validate_wfs_optical_products(optical_product)
    validate_wfs_observation(observation)
    observation.metadata.numeric_type <: Real ||
        throw(WFSPreparationError(:acquisition, :numeric_type,
            "WFS detector observations require real sample storage"))
    _require_wfs_detector_observation_ownership(
        detector, optical_product, observation)
    candidate = _prepare_detached_detector_acquisition(
        detector, optical_product)
    candidate_detector = detector_acquisition_detector(candidate)
    _validate_wfs_detector_observation_contract(
        candidate_detector, optical_product, observation)
    _require_wfs_detector_observation_ownership(
        candidate_detector, optical_product, observation)
    acquisition = _rebind_prepared_detector_acquisition(detector, candidate)
    plan = WFSDetectorAcquisitionPlan(
        _WFS_DETECTOR_ACQUISITION_PLAN_TOKEN,
        detector_acquisition_plan(acquisition), observation.metadata)
    prepared = PreparedWFSDetectorAcquisition(
        _PREPARED_WFS_DETECTOR_ACQUISITION_TOKEN,
        plan, acquisition, observation)
    return prepared
end

@inline function _commit_wfs_acquisition_candidate!(
    prepared::PreparedWFSDetectorAcquisition, ::Detector,
    ::IntensityMap, ::WFSObservation)
    _commit_prepared_detector_acquisition!(prepared.acquisition)
    return prepared
end

function prepare_wfs_acquisition(detector::Detector,
    optical_product::IntensityMap, observation::WFSObservation;
    source=nothing)
    candidate = _prepare_wfs_acquisition_candidate(
        detector, optical_product, observation; source=source)
    return _commit_wfs_acquisition_candidate!(
        candidate, detector, optical_product, observation)
end

function acquire_wfs_observation!(observation::WFSObservation,
    optical_product::IntensityMap, plan::PreparedWFSDetectorAcquisition,
    rng::AbstractRNG)
    validate_wfs_acquisition_binding(observation, optical_product, plan)
    frame = capture!(plan.acquisition, rng)
    copyto!(observation.storage, frame)
    return observation
end

function validate_wfs_acquisition_binding(observation::WFSObservation,
    optical_product::IntensityMap, plan::PreparedWFSDetectorAcquisition)
    observation === plan.observation || throw(WFSPreparationError(
        :acquisition, :prepared_binding,
            "WFS observation does not match prepared storage"))
    observation.metadata === plan.plan.observation_metadata || throw(
        WFSPreparationError(:acquisition, :prepared_binding,
            "WFS observation metadata changed after preparation"))
    detector_acquisition_plan(plan.acquisition) ===
        plan.plan.detector_plan || throw(WFSPreparationError(
            :acquisition, :prepared_binding,
            "WFS detector acquisition plan changed after preparation"))
    prepared_input = detector_acquisition_input(plan.acquisition)
    optical_product === prepared_input || throw(
        WFSPreparationError(:acquisition, :prepared_binding,
            "WFS optical product does not match its prepared detector input"))
    try
        _require_prepared_acquisition(plan.acquisition)
    catch error
        error isa Union{InvalidConfiguration,DimensionMismatchError} ||
            rethrow()
        throw(WFSPreparationError(:acquisition, :prepared_binding,
            "WFS detector ownership changed after preparation"))
    end
    return nothing
end

struct _WFSCountingAcquisitionPlanToken end
const _WFS_COUNTING_ACQUISITION_PLAN_TOKEN =
    _WFSCountingAcquisitionPlanToken()

"""Run-immutable radiometric contract for accumulated-count WFS acquisition."""
struct WFSCountingAcquisitionPlan{
    I<:OpticalPlaneMetadata,
    O<:WFSObservationMetadata,
    T<:AbstractFloat,
} <: AbstractWFSAcquisitionPlan
    input_metadata::I
    observation_metadata::O
    source_throughput::T

    function WFSCountingAcquisitionPlan(
        ::_WFSCountingAcquisitionPlanToken, input_metadata::I,
        observation_metadata::O,
        source_throughput::T) where {
        I<:OpticalPlaneMetadata,O<:WFSObservationMetadata,T<:AbstractFloat,
    }
        return new{I,O,T}(
            input_metadata, observation_metadata, source_throughput)
    end
end

struct _PreparedWFSCountingAcquisitionToken end
const _PREPARED_WFS_COUNTING_ACQUISITION_TOKEN =
    _PreparedWFSCountingAcquisitionToken()

"""Exact prepared owner for one accumulated-count WFS acquisition."""
struct PreparedWFSCountingAcquisition{
    P<:WFSCountingAcquisitionPlan,
    D<:AbstractCountingDetector,
    I<:IntensityMap,
    O<:WFSObservation,
    S,
    W,
    R,
    B,
}
    plan::P
    detector::D
    optical_product::I
    observation::O
    state::S
    workspace::W
    products::R
    detector_binding::B

    function PreparedWFSCountingAcquisition(
        ::_PreparedWFSCountingAcquisitionToken, plan::P,
        detector::D, optical_product::I, observation::O,
        state::S, workspace::W, products::R,
        detector_binding::B) where {
        P<:WFSCountingAcquisitionPlan,D<:AbstractCountingDetector,
        I<:IntensityMap,O<:WFSObservation,S,W,R,B,
    }
        return new{P,D,I,O,S,W,R,B}(plan, detector, optical_product,
            observation, state, workspace, products, detector_binding)
    end
end

@inline wfs_acquisition_plan(
    prepared::PreparedWFSCountingAcquisition) = prepared.plan

@inline function _counting_wfs_detector_binding(
    detector::AbstractCountingDetector)
    return (
        thermal_state=thermal_state(detector),
        noise_buffer=counting_noise_buffer(detector),
        host_buffer=counting_host_buffer(detector),
        output_buffer=counting_output_buffer(detector),
        output_buffer_host=counting_output_host_buffer(detector),
    )
end

@inline function _require_counting_wfs_detector_binding(
    detector::AbstractCountingDetector, binding::NamedTuple)
    thermal_state(detector) === binding.thermal_state &&
        counting_noise_buffer(detector) === binding.noise_buffer &&
        counting_host_buffer(detector) === binding.host_buffer &&
        counting_output_buffer(detector) === binding.output_buffer &&
        counting_output_host_buffer(detector) === binding.output_buffer_host ||
        throw(WFSPreparationError(:acquisition, :prepared_binding,
            "counting detector state, workspace, or products changed after WFS preparation"))
    return nothing
end

@inline function _counting_wfs_storage_mightalias(
    detector::AbstractCountingDetector, value::AbstractArray)
    output = counting_output_buffer(detector)
    output_host = counting_output_host_buffer(detector)
    return _wfs_storage_mightalias(value, counting_array(detector)) ||
        _wfs_storage_mightalias(value, counting_noise_buffer(detector)) ||
        _wfs_storage_mightalias(value, counting_host_buffer(detector)) ||
        _wfs_storage_mightalias(value, output) ||
        _wfs_storage_mightalias(value, output_host)
end

function _require_counting_wfs_storage_ownership(
    detector::AbstractCountingDetector, optical_product::IntensityMap,
    observation::WFSObservation)
    input = optical_product.values
    output = observation.storage
    (_counting_wfs_storage_mightalias(detector, input) ||
        _counting_wfs_storage_mightalias(detector, output) ||
        _wfs_storage_mightalias(input, output)) && throw(
        WFSPreparationError(:acquisition, :ownership,
            "counting WFS input and observation must not alias detector storage or each other"))
    return nothing
end

@inline _counting_wfs_source_required(::AbstractCountingDetector) = false
@inline function _counting_wfs_source_required(detector::MKIDArrayDetector)
    return detector.params.sensor.characteristics.wavelength_passband_m !==
        nothing
end

function _require_counting_wfs_source(detector::AbstractCountingDetector,
    ::Nothing)
    _counting_wfs_source_required(detector) && throw(WFSPreparationError(
        :acquisition, :radiometry,
        "wavelength-selective counting acquisition requires its source"))
    return nothing
end

function _require_counting_wfs_source(::AbstractCountingDetector,
    source::AbstractSource)
    require_leaf_source(source, "counting WFS acquisition")
    return source
end

function _require_counting_wfs_source(::AbstractCountingDetector, source)
    throw(WFSPreparationError(:acquisition, :radiometry,
        "counting WFS acquisition source must be a leaf source or nothing"))
end

function _require_counting_wfs_spectral_match(
    metadata::OpticalPlaneMetadata, ::Nothing)
    return nothing
end

function _require_counting_wfs_spectral_match(
    metadata::OpticalPlaneMetadata, source::AbstractSource)
    return _require_counting_wfs_spectral_match(metadata.spectral, source)
end

function _require_counting_wfs_spectral_match(
    channel::MonochromaticChannel, source::AbstractSource)
    channel.wavelength_m == wavelength(source) || throw(
        WFSPreparationError(:acquisition, :plane_metadata,
            "counting WFS source wavelength differs from its rate product"))
    return nothing
end

function _require_counting_wfs_spectral_match(
    ::AbstractSpectralCoordinate, ::AbstractSource)
    throw(WFSPreparationError(:acquisition, :plane_metadata,
        "counting WFS acquisition requires a monochromatic rate product"))
end

@inline _require_counting_wfs_measure(::CellIntegratedMeasure) = nothing

function _require_counting_wfs_measure(::AbstractSpatialMeasure)
    throw(WFSPreparationError(:acquisition, :radiometry,
        "counting WFS inputs must carry cell-integrated photon rate"))
end

function _prepare_wfs_acquisition_candidate(
    detector::AbstractCountingDetector,
    optical_product::IntensityMap, observation::WFSObservation;
    source=nothing)
    validate_wfs_optical_products(optical_product)
    validate_wfs_observation(observation)
    _require_counting_wfs_measure(optical_product.metadata.spatial_measure)
    observation.metadata.numeric_type <: Real || throw(
        WFSPreparationError(:acquisition, :numeric_type,
            "counting WFS observations require real sample storage"))
    _require_counting_wfs_source(detector, source)
    _require_counting_wfs_spectral_match(optical_product.metadata, source)
    input = optical_product.values
    eltype(input) === eltype(counting_array(detector)) || throw(
        WFSPreparationError(:acquisition, :numeric_type,
            "counting detector and WFS rate product precision differ"))
    typeof(backend(detector)) === typeof(backend(input)) || throw(
        WFSPreparationError(:acquisition, :backend,
            "counting detector and WFS rate product backends differ"))
    compute_device(counting_array(detector)) == compute_device(input) || throw(
        WFSPreparationError(:acquisition, :device,
            "counting detector and WFS rate product occupy different devices"))
    _require_finite_nonnegative_intensity(input)
    size(observation.storage) == size(input) || throw(WFSPreparationError(
        :acquisition, :shape,
        "counting WFS observation storage must match detector output"))
    output_type = detector_output_type(detector)
    expected_type = output_type === nothing ? eltype(counting_array(detector)) :
        output_type
    observation.metadata.numeric_type === expected_type || throw(
        WFSPreparationError(:acquisition, :numeric_type,
            "counting WFS observation element type must match detector output"))
    typeof(backend(observation.storage)) === typeof(backend(detector)) ||
        throw(WFSPreparationError(:acquisition, :backend,
            "counting WFS observation and detector backends differ"))
    compute_device(observation.storage) ==
        compute_device(counting_array(detector)) || throw(
        WFSPreparationError(:acquisition, :device,
            "counting WFS observation and detector output occupy different devices"))
    _require_counting_wfs_storage_ownership(
        detector, optical_product, observation)
    output_type = eltype(counting_array(detector))
    source_throughput = source === nothing ? one(output_type) :
        counting_source_throughput(detector, source, output_type)
    plan = WFSCountingAcquisitionPlan(
        _WFS_COUNTING_ACQUISITION_PLAN_TOKEN,
        optical_product.metadata, observation.metadata,
        source_throughput)
    return plan
end

function _commit_wfs_acquisition_candidate!(
    plan::WFSCountingAcquisitionPlan,
    detector::AbstractCountingDetector,
    optical_product::IntensityMap, observation::WFSObservation)
    input = optical_product.values
    ensure_buffers!(detector, size(input))
    return PreparedWFSCountingAcquisition(
        _PREPARED_WFS_COUNTING_ACQUISITION_TOKEN,
        plan, detector, optical_product, observation,
        counting_detector_state(detector),
        counting_detector_workspace(detector),
        counting_detector_products(detector),
        _counting_wfs_detector_binding(detector))
end

function prepare_wfs_acquisition(detector::AbstractCountingDetector,
    optical_product::IntensityMap, observation::WFSObservation;
    source=nothing)
    candidate = _prepare_wfs_acquisition_candidate(
        detector, optical_product, observation; source=source)
    return _commit_wfs_acquisition_candidate!(
        candidate, detector, optical_product, observation)
end

function acquire_wfs_observation!(observation::WFSObservation,
    optical_product::IntensityMap,
    plan::PreparedWFSCountingAcquisition, rng::AbstractRNG)
    validate_wfs_acquisition_binding(observation, optical_product, plan)
    frame = _capture_prevalidated_counting!(plan.detector,
        optical_product.values, plan.plan.source_throughput, rng)
    copyto!(observation.storage, frame)
    return observation
end

function validate_wfs_acquisition_binding(observation::WFSObservation,
    optical_product::IntensityMap, plan::PreparedWFSCountingAcquisition)
    observation === plan.observation &&
        optical_product === plan.optical_product || throw(
        WFSPreparationError(:acquisition, :prepared_binding,
            "counting WFS storage does not match its prepared acquisition"))
    counting_detector_state(plan.detector) === plan.state &&
        counting_detector_workspace(plan.detector) === plan.workspace &&
        counting_detector_products(plan.detector) === plan.products || throw(
        WFSPreparationError(:acquisition, :prepared_binding,
            "counting detector owner changed after WFS preparation"))
    optical_product.metadata === plan.plan.input_metadata &&
        size(optical_product.values) == plan.plan.input_metadata.dimensions &&
        observation.metadata === plan.plan.observation_metadata || throw(
        WFSPreparationError(:acquisition, :prepared_binding,
            "counting WFS plan metadata changed after preparation"))
    _require_counting_wfs_detector_binding(
        plan.detector, plan.detector_binding)
    _require_counting_wfs_storage_ownership(
        plan.detector, optical_product, observation)
    return nothing
end

struct _WFSMultipleDetectorAcquisitionPlanToken end
const _WFS_MULTIPLE_DETECTOR_ACQUISITION_PLAN_TOKEN =
    _WFSMultipleDetectorAcquisitionPlanToken()

"""Run-immutable static fan-out of concrete WFS acquisition plans."""
struct WFSMultipleDetectorAcquisitionPlan{
    P<:Tuple,
} <: AbstractWFSAcquisitionPlan
    component_plans::P

    function WFSMultipleDetectorAcquisitionPlan(
        ::_WFSMultipleDetectorAcquisitionPlanToken,
        component_plans::P) where {P<:Tuple}
        return new{P}(component_plans)
    end
end

struct _PreparedWFSMultipleDetectorAcquisitionToken end
const _PREPARED_WFS_MULTIPLE_DETECTOR_ACQUISITION_TOKEN =
    _PreparedWFSMultipleDetectorAcquisitionToken()

"""Exact prepared owner for a static tuple of WFS acquisitions."""
struct PreparedWFSMultipleDetectorAcquisition{
    P<:WFSMultipleDetectorAcquisitionPlan,
    A<:Tuple,
    I<:Tuple,
    O<:Tuple,
}
    plan::P
    acquisitions::A
    optical_products::I
    observations::O

    function PreparedWFSMultipleDetectorAcquisition(
        ::_PreparedWFSMultipleDetectorAcquisitionToken,
        plan::P, acquisitions::A, optical_products::I,
        observations::O) where {
        P<:WFSMultipleDetectorAcquisitionPlan,A<:Tuple,
        I<:Tuple,O<:Tuple,
    }
        return new{P,A,I,O}(
            plan, acquisitions, optical_products, observations)
    end
end

@inline wfs_acquisition_plan(
    prepared::PreparedWFSMultipleDetectorAcquisition) = prepared.plan

@inline function _wfs_detector_storage_mightalias(
    detector::Detector, value::AbstractArray)
    return _detector_storage_mightalias(detector, value)
end

@inline _wfs_detector_storage_mightalias(detector, value) = false

@inline function _wfs_detector_storage_mightalias(
    detector::AbstractCountingDetector, value::AbstractArray)
    return _counting_wfs_storage_mightalias(detector, value)
end

function _require_multiple_wfs_acquisition_ownership(
    detectors::Tuple, optical_products::Tuple, observations::Tuple)
    count = length(detectors)
    @inbounds for first_index in 1:count
        first_detector = detectors[first_index]
        first_observation = observations[first_index].storage
        first_input = optical_products[first_index].values
        for second_index in (first_index + 1):count
            first_detector === detectors[second_index] && throw(
                WFSPreparationError(:acquisition, :ownership,
                    "multiple-detector WFS acquisition requires distinct detector owners"))
            _wfs_storage_mightalias(first_observation,
                observations[second_index].storage) && throw(
                WFSPreparationError(:acquisition, :ownership,
                    "multiple-detector WFS observations must not alias"))
        end
        for component_index in 1:count
            component_input = optical_products[component_index].values
            component_observation = observations[component_index].storage
            component_detector = detectors[component_index]
            (_wfs_storage_mightalias(first_observation, component_input) ||
                _wfs_detector_storage_mightalias(
                    component_detector, first_observation) ||
                _wfs_detector_storage_mightalias(
                    component_detector, first_input)) && throw(
                WFSPreparationError(:acquisition, :ownership,
                    "multiple-detector WFS inputs, observations, and detector storage must be disjoint"))
            _wfs_storage_mightalias(
                first_input, component_observation) && throw(
                WFSPreparationError(:acquisition, :ownership,
                    "multiple-detector WFS observations must not alias an optical input"))
        end
    end
    return nothing
end

function prepare_wfs_acquisition(detectors::Tuple,
    optical_products::Tuple, observations::Tuple; source=nothing)
    isempty(detectors) && throw(WFSPreparationError(:acquisition,
        :plane_count, "multiple-detector WFS acquisition must not be empty"))
    length(detectors) == length(optical_products) == length(observations) ||
        throw(WFSPreparationError(:acquisition, :plane_count,
            "detector, optical-product, and observation counts must match"))
    _require_multiple_wfs_acquisition_ownership(
        detectors, optical_products, observations)
    candidates = map(detectors, optical_products, observations) do detector,
            product, observation
        _prepare_wfs_acquisition_candidate(detector, product, observation;
            source=source)
    end
    acquisitions = map(candidates, detectors, optical_products,
        observations) do candidate, detector, product, observation
        _commit_wfs_acquisition_candidate!(
            candidate, detector, product, observation)
    end
    plan = WFSMultipleDetectorAcquisitionPlan(
        _WFS_MULTIPLE_DETECTOR_ACQUISITION_PLAN_TOKEN,
        map(wfs_acquisition_plan, acquisitions))
    return PreparedWFSMultipleDetectorAcquisition(
        _PREPARED_WFS_MULTIPLE_DETECTOR_ACQUISITION_TOKEN,
        plan, acquisitions, optical_products, observations)
end

@inline _acquire_multiple_wfs_observations!(::Tuple{}, ::Tuple{},
    ::Tuple{}, ::AbstractRNG) = nothing

@inline function _acquire_multiple_wfs_observations!(plans::Tuple,
    observations::Tuple, optical_products::Tuple, rng::AbstractRNG)
    acquire_wfs_observation!(first(observations), first(optical_products),
        first(plans), rng)
    return _acquire_multiple_wfs_observations!(Base.tail(plans),
        Base.tail(observations), Base.tail(optical_products), rng)
end

function acquire_wfs_observation!(observations::Tuple,
    optical_products::Tuple,
    plan::PreparedWFSMultipleDetectorAcquisition, rng::AbstractRNG)
    validate_wfs_acquisition_binding(observations, optical_products, plan)
    _acquire_multiple_wfs_observations!(plan.acquisitions, observations,
        optical_products, rng)
    return observations
end

@inline _validate_multiple_wfs_acquisition_bindings(::Tuple{}, ::Tuple{},
    ::Tuple{}) = nothing

@inline _require_multiple_wfs_plan_bindings(::Tuple{}, ::Tuple{}) = nothing

@inline function _require_multiple_wfs_plan_bindings(
    component_plans::Tuple, acquisitions::Tuple)
    first(component_plans) === wfs_acquisition_plan(first(acquisitions)) ||
        throw(WFSPreparationError(:acquisition, :prepared_binding,
            "multiple-detector WFS component plan changed after preparation"))
    return _require_multiple_wfs_plan_bindings(
        Base.tail(component_plans), Base.tail(acquisitions))
end

@inline function _validate_multiple_wfs_acquisition_bindings(plans::Tuple,
    observations::Tuple, optical_products::Tuple)
    validate_wfs_acquisition_binding(first(observations),
        first(optical_products), first(plans))
    return _validate_multiple_wfs_acquisition_bindings(Base.tail(plans),
        Base.tail(observations), Base.tail(optical_products))
end

function validate_wfs_acquisition_binding(observations::Tuple,
    optical_products::Tuple,
    plan::PreparedWFSMultipleDetectorAcquisition)
    observations === plan.observations &&
        optical_products === plan.optical_products || throw(
        WFSPreparationError(:acquisition, :prepared_binding,
            "multiple-detector WFS storage does not match its prepared acquisition"))
    _require_multiple_wfs_plan_bindings(
        plan.plan.component_plans, plan.acquisitions)
    _validate_multiple_wfs_acquisition_bindings(plan.acquisitions, observations,
        optical_products)
    return nothing
end

"""
    prepare_wfs_estimation(model, input, measurement)
    estimate_wfs_measurement!(measurement, input, prepared)

Preparation and mutating execution protocol for a WFS estimator. An acquired
estimator consumes one observation or a concrete tuple of observations. A
declared direct estimator instead consumes explicit pupil/field input and owns
no fictitious optical or detector products.
"""
function prepare_wfs_estimation(model, input, measurement)
    throw(WFSPreparationError(:estimation, :unsupported,
        "$(typeof(model)) does not implement prepared WFS estimation"))
end

function estimate_wfs_measurement! end

"""
    validate_wfs_estimation_binding(measurement, input, prepared)

Validate without mutation that a prepared WFS estimator is bound to the exact
input observation and measurement supplied by a containing plant owner.
"""
function validate_wfs_estimation_binding(measurement, input, prepared)
    throw(WFSPreparationError(:estimation,
        :unsupported_binding_validation,
        "$(typeof(prepared)) does not validate its prepared estimator binding"))
end

abstract type AbstractWFSMeasurementPath end

"""The estimator consumes one or more acquired WFS observations."""
struct AcquiredObservationPath <: AbstractWFSMeasurementPath end

"""The estimator intentionally produces a measurement without optical/acquisition intermediates."""
struct DirectMeasurementPath <: AbstractWFSMeasurementPath end

function wfs_measurement_path(prepared)
    throw(WFSPreparationError(:estimation, :estimator,
        "$(typeof(prepared)) must declare AcquiredObservationPath() or DirectMeasurementPath()"))
end

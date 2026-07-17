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
    D<:AbstractPlaneDevice,
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
    D<:AbstractPlaneDevice,
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
@inline _wfs_storage_device(storage::AbstractArray) = plane_device(storage)
@inline _wfs_storage_device(::Base.RefValue) = HostPlaneDevice()
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
    throw(WFSPreparationError(:optical_formation, :plane_metadata,
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
    throw(WFSPreparationError(:optical_formation, :plane_metadata,
        "WFS optical formation requires a PupilFunction or pupil-plane ElectricField"))
end

@inline _require_wfs_output_plane(::DetectorPlane) = nothing

function _require_wfs_output_plane(::AbstractOpticalPlaneKind)
    throw(WFSPreparationError(:optical_formation, :plane_metadata,
        "WFS optical output must be declared on a detector plane"))
end

@inline _require_wfs_rate_normalization(::PhotonRateNormalization) = nothing

function _require_wfs_rate_normalization(::AbstractOpticalNormalization)
    throw(WFSPreparationError(:optical_formation, :radiometry,
        "WFS detector-facing optical output must use photon-rate normalization"))
end

@inline _require_wfs_rate_measure(::SpatialDensityMeasure) = nothing
@inline _require_wfs_rate_measure(::CellIntegratedMeasure) = nothing

function _require_wfs_rate_measure(::AbstractSpatialMeasure)
    throw(WFSPreparationError(:optical_formation, :radiometry,
        "WFS detector-facing optical output must declare spatial-density or cell-integrated measure"))
end

@inline _require_wfs_rate_coherence(::IncoherentIntensityAddition) = nothing

function _require_wfs_rate_coherence(::AbstractCombinationPolicy)
    throw(WFSPreparationError(:optical_formation, :radiometry,
        "WFS detector-facing optical output must declare incoherent-intensity semantics"))
end

@inline _require_wfs_rate_spectral(::MonochromaticChannel) = nothing
@inline _require_wfs_rate_spectral(::IntegratedSpectralChannel) = nothing

function _require_wfs_rate_spectral(::AbstractSpectralCoordinate)
    throw(WFSPreparationError(:optical_formation, :plane_metadata,
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

@inline validate_wfs_optical_products(::Tuple{}) = throw(WFSPreparationError(
    :optical_formation, :plane_count,
    "WFS optical formation requires at least one detector-facing product"))

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
    throw(WFSPreparationError(:optical_formation, :plane_metadata,
        "WFS optical products must be IntensityMap values or a concrete tuple/bundle of them"))
end

"""
    prepare_wfs_optical_formation(model, input, output)
    form_wfs_optical_products!(output, input, prepared)

Preparation and mutating execution protocol for a WFS optical front end.
Implementations consume an explicit `PupilFunction` or pupil-plane
`ElectricField` and write caller-owned detector-plane photon-arrival-rate
products. Preparation performs all fallible structural validation.
"""
function prepare_wfs_optical_formation(model, input, output)
    throw(WFSPreparationError(:optical_formation, :unsupported,
        "$(typeof(model)) does not implement prepared WFS optical formation"))
end

function form_wfs_optical_products! end

"""
    prepare_wfs_acquisition(model, optical_products, observation)
    acquire_wfs_observation!(observation, optical_products, prepared, rng)

Preparation and mutating execution protocol for detector acquisition. The
prepared model owns detector/acquisition state and explicit durations; optical
inputs remain photon-arrival rates. Concrete tuples express fan-out to multiple
detectors without a runtime-selected stage graph.
"""
function prepare_wfs_acquisition(model, optical_products, observation)
    throw(WFSPreparationError(:acquisition, :unsupported,
        "$(typeof(model)) does not implement prepared WFS acquisition"))
end

function acquire_wfs_observation! end

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

abstract type AbstractWFSMeasurementPath end

"""The estimator consumes one or more acquired WFS observations."""
struct AcquiredObservationPath <: AbstractWFSMeasurementPath end

"""The estimator intentionally produces a measurement without optical/acquisition intermediates."""
struct DirectMeasurementPath <: AbstractWFSMeasurementPath end

function wfs_measurement_path(prepared)
    throw(WFSPreparationError(:estimation, :estimator,
        "$(typeof(prepared)) must declare AcquiredObservationPath() or DirectMeasurementPath()"))
end

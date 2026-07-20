#
# Prepared illumination entry boundaries
#
# Calibration is an application role, not a propagation mode. These values
# describe where an evaluator writes an ordinary optical product and which
# composition semantics the evaluator has already applied.
#

"""Entry whose evaluator writes a caller-owned `PupilFunction`."""
struct PupilFunctionIlluminationEntry end

"""Entry whose evaluator writes a caller-owned `ElectricField`."""
struct ElectricFieldIlluminationEntry end

"""Entry whose evaluator writes a declared `IntensityMap`."""
struct IntensityMapIlluminationEntry end

"""Entry for a typed field or intensity result formed by external optics."""
struct ExternalOpticsResultIlluminationEntry end

"""Entry for an `IntensityMap` immediately upstream of detector acquisition."""
struct DetectorInputIlluminationEntry end

"""One illumination contribution writes the complete entry payload."""
struct SingleIllumination end

"""An evaluator explicitly selects exactly one of its illumination inputs."""
struct ExclusiveIlluminationSelection end

struct _UnsupportedIlluminationCombination end

"""
    illumination_combination(::Type{T})

Trait declaring how a prepared illumination evaluator forms its complete
payload. Supported values are `SingleIllumination()`,
`ExclusiveIlluminationSelection()`, `CoherentFieldCombination()`, and
`IncoherentIntensityAddition()`. Core never infers this choice from a
calibration role or source count.
"""
illumination_combination(::Type) = _UnsupportedIlluminationCombination()

@inline illumination_combination(evaluator) =
    illumination_combination(typeof(evaluator))

struct _IlluminationArrayContract{N,E,B<:AbstractArrayBackend,
    D<:AbstractPlaneDevice}
    dimensions::NTuple{N,Int}
    numeric_type::Type{E}
    backend::B
    device::D
end

struct _PupilIlluminationPayloadContract{M,S,A,O}
    metadata::M
    support::S
    amplitude::A
    opd::O
    aperture_revision::UInt
end

struct _OpticalIlluminationPayloadContract{M,V}
    metadata::M
    values::V
end

@inline function _illumination_array_contract(storage::AbstractArray)
    selector = backend(storage)
    device = plane_device(storage)
    return _IlluminationArrayContract(size(storage), eltype(storage),
        selector, device)
end

@inline function _validate_illumination_pupil_support(
    payload::PupilFunction)
    size(payload.support) == payload.metadata.dimensions || throw(
        PlantPreparationError(:illumination, :shape,
            "pupil illumination support shape does not match its metadata"))
    typeof(backend(payload.support)) === typeof(payload.metadata.backend) ||
        throw(PlantPreparationError(:illumination, :backend,
            "pupil illumination support backend does not match its metadata"))
    plane_device(payload.support) == payload.metadata.device || throw(
        PlantPreparationError(:illumination, :device,
            "pupil illumination support device does not match its metadata"))
    return payload
end

function illumination_payload_contract(payload::PupilFunction)
    _validate_path_input(payload)
    _validate_illumination_pupil_support(payload)
    return _PupilIlluminationPayloadContract(deepcopy(payload.metadata),
        _illumination_array_contract(payload.support),
        _illumination_array_contract(payload.amplitude),
        _illumination_array_contract(payload.opd),
        aperture_revision(payload))
end

function illumination_payload_contract(
    payload::Union{ElectricField,IntensityMap})
    _validate_path_input(payload)
    return _OpticalIlluminationPayloadContract(deepcopy(payload.metadata),
        _illumination_array_contract(payload.values))
end

function illumination_payload_contract(payload)
    throw(PlantPreparationError(:illumination, :unsupported_payload,
        "illumination entry payload type $(typeof(payload)) is unsupported"))
end

@inline function _validate_illumination_metadata(metadata, expected,
    label::AbstractString)
    typeof(metadata) === typeof(expected) || throw(PlantPreparationError(
        :illumination, :metadata,
        "$label metadata type changed after preparation"))
    isequal(metadata, expected) || throw(PlantPreparationError(
        :illumination, :metadata,
        "$label metadata changed after preparation"))
    return metadata
end

@inline function _validate_illumination_array(storage::AbstractArray,
    contract::_IlluminationArrayContract, label::AbstractString)
    size(storage) == contract.dimensions || throw(PlantPreparationError(
        :illumination, :shape,
        "$label shape does not match its prepared contract"))
    eltype(storage) === contract.numeric_type || throw(
        PlantPreparationError(:illumination, :numeric_type,
            "$label element type does not match its prepared contract"))
    typeof(backend(storage)) === typeof(contract.backend) || throw(
        PlantPreparationError(:illumination, :backend,
            "$label backend does not match its prepared contract"))
    plane_device(storage) == contract.device || throw(
        PlantPreparationError(:illumination, :device,
            "$label device does not match its prepared contract"))
    return storage
end

function validate_illumination_payload_contract(payload::PupilFunction,
    contract::_PupilIlluminationPayloadContract)
    _validate_path_input(payload)
    _validate_illumination_metadata(payload.metadata, contract.metadata,
        "pupil illumination")
    _validate_illumination_array(payload.support, contract.support,
        "pupil support")
    _validate_illumination_array(payload.amplitude, contract.amplitude,
        "pupil amplitude")
    _validate_illumination_array(payload.opd, contract.opd,
        "pupil OPD")
    aperture_revision(payload) == contract.aperture_revision || throw(
        PlantPreparationError(:illumination, :revision,
            "pupil illumination aperture revision changed"))
    return payload
end

function validate_illumination_payload_contract(
    payload::Union{ElectricField,IntensityMap},
    contract::_OpticalIlluminationPayloadContract)
    _validate_path_input(payload)
    _validate_illumination_metadata(payload.metadata, contract.metadata,
        "illumination payload")
    _validate_illumination_array(payload.values, contract.values,
        "illumination payload")
    return payload
end

function validate_illumination_payload_contract(payload, contract)
    throw(PlantPreparationError(:illumination, :payload_type,
        "illumination payload type $(typeof(payload)) is incompatible with contract $(typeof(contract))"))
end

@inline _require_declared_illumination_spectral(
    ::AbstractSpectralCoordinate) = nothing

function _require_declared_illumination_spectral(
    ::UnspecifiedSpectralCoordinate)
    throw(PlantPreparationError(:illumination, :spectral_sampling,
        "illumination payload spectral coordinates must be declared"))
end

@inline _require_declared_illumination_normalization(
    ::AbstractOpticalNormalization) = nothing

function _require_declared_illumination_normalization(
    ::UnspecifiedNormalization)
    throw(PlantPreparationError(:illumination, :radiometry,
        "illumination payload normalization must be declared"))
end

@inline _require_declared_illumination_measure(
    ::AbstractSpatialMeasure) = nothing

function _require_declared_illumination_measure(::UnspecifiedSpatialMeasure)
    throw(PlantPreparationError(:illumination, :radiometry,
        "illumination payload spatial measure must be declared"))
end

@inline _require_declared_illumination_coherence(
    ::AbstractCombinationPolicy) = nothing

function _require_declared_illumination_coherence(::UnspecifiedCoherence)
    throw(PlantPreparationError(:illumination, :combination,
        "illumination payload coherence must be declared"))
end

@inline function _require_declared_illumination_metadata(payload)
    metadata = payload.metadata
    _require_declared_illumination_spectral(metadata.spectral)
    _require_declared_illumination_normalization(metadata.normalization)
    _require_declared_illumination_measure(metadata.spatial_measure)
    _require_declared_illumination_coherence(metadata.coherence)
    return payload
end

@inline _require_pupil_illumination_coherence(
    ::CoherentFieldCombination) = nothing

function _require_pupil_illumination_coherence(coherence)
    throw(PlantPreparationError(:illumination, :combination,
        "pupil-function illumination must declare coherent field semantics; got $(typeof(coherence))"))
end

@inline _require_field_illumination_coherence(
    ::CoherentFieldCombination) = nothing
@inline _require_field_illumination_coherence(
    ::NonCombinableProduct) = nothing

function _require_field_illumination_coherence(coherence)
    throw(PlantPreparationError(:illumination, :combination,
        "electric-field illumination must declare coherent or non-combinable semantics; got $(typeof(coherence))"))
end

@inline _require_intensity_illumination_coherence(
    ::IncoherentIntensityAddition) = nothing
@inline _require_intensity_illumination_coherence(
    ::NonCombinableProduct) = nothing

function _require_intensity_illumination_coherence(coherence)
    throw(PlantPreparationError(:illumination, :combination,
        "intensity illumination must declare incoherent or non-combinable semantics; got $(typeof(coherence))"))
end

@inline function _validate_illumination_product_semantics(
    payload::PupilFunction)
    _require_pupil_illumination_coherence(payload.metadata.coherence)
    return payload
end

@inline function _validate_illumination_product_semantics(
    payload::ElectricField)
    _require_field_illumination_coherence(payload.metadata.coherence)
    return payload
end

@inline function _validate_illumination_product_semantics(
    payload::IntensityMap)
    _require_intensity_illumination_coherence(payload.metadata.coherence)
    return payload
end

@inline function validate_illumination_entry_payload(
    ::PupilFunctionIlluminationEntry, payload::PupilFunction)
    _validate_path_input(payload)
    _require_declared_illumination_metadata(payload)
    _validate_illumination_product_semantics(payload)
    return payload
end

@inline function validate_illumination_entry_payload(
    ::ElectricFieldIlluminationEntry, payload::ElectricField)
    _validate_path_input(payload)
    _require_declared_illumination_metadata(payload)
    _validate_illumination_product_semantics(payload)
    return payload
end

@inline function validate_illumination_entry_payload(
    ::IntensityMapIlluminationEntry, payload::IntensityMap)
    _validate_path_input(payload)
    _require_declared_illumination_metadata(payload)
    _validate_illumination_product_semantics(payload)
    return payload
end

@inline function validate_illumination_entry_payload(
    ::ExternalOpticsResultIlluminationEntry,
    payload::Union{ElectricField,IntensityMap})
    _validate_path_input(payload)
    _require_declared_illumination_metadata(payload)
    _validate_illumination_product_semantics(payload)
    return payload
end

@inline _require_detector_entry_plane(::FocalPlane) = nothing
@inline _require_detector_entry_plane(::DetectorPlane) = nothing

function _require_detector_entry_plane(kind::AbstractOpticalPlaneKind)
    throw(PlantPreparationError(:illumination, :entry_plane,
        "detector-input illumination must enter on a focal or detector plane; got $(typeof(kind))"))
end

@inline _require_detector_entry_normalization(
    ::PhotonRateNormalization) = nothing
@inline _require_detector_entry_normalization(
    ::DimensionlessNormalization) = nothing

function _require_detector_entry_normalization(
    normalization::AbstractOpticalNormalization)
    throw(PlantPreparationError(:illumination, :radiometry,
        "detector-input illumination must declare photon-rate or dimensionless normalization; got $(typeof(normalization))"))
end

@inline _require_detector_entry_measure(::SpatialDensityMeasure) = nothing
@inline _require_detector_entry_measure(::CellIntegratedMeasure) = nothing

function _require_detector_entry_measure(measure::AbstractSpatialMeasure)
    throw(PlantPreparationError(:illumination, :radiometry,
        "detector-input illumination must use spatial-density or cell-integrated samples; got $(typeof(measure))"))
end

@inline _require_detector_entry_coherence(
    ::IncoherentIntensityAddition) = nothing

function _require_detector_entry_coherence(
    coherence::AbstractCombinationPolicy)
    throw(PlantPreparationError(:illumination, :combination,
        "detector-input illumination must declare incoherent intensity semantics; got $(typeof(coherence))"))
end

function validate_illumination_entry_payload(
    ::DetectorInputIlluminationEntry, payload::IntensityMap)
    _validate_path_input(payload)
    _require_declared_illumination_spectral(payload.metadata.spectral)
    _require_detector_entry_plane(payload.metadata.kind)
    _require_detector_entry_normalization(payload.metadata.normalization)
    _require_detector_entry_measure(payload.metadata.spatial_measure)
    _require_detector_entry_coherence(payload.metadata.coherence)
    return payload
end

function validate_illumination_entry_payload(boundary, payload)
    throw(PlantPreparationError(:illumination, :entry_payload,
        "illumination entry boundary $(typeof(boundary)) does not accept payload $(typeof(payload))"))
end

@inline _require_illumination_combination(::SingleIllumination,
    payload) = SingleIllumination()

@inline _require_illumination_combination(
    ::ExclusiveIlluminationSelection, payload) =
    ExclusiveIlluminationSelection()

function _require_illumination_combination(
    ::CoherentFieldCombination, payload::ElectricField)
    typeof(payload.metadata.coherence) === CoherentFieldCombination ||
        throw(PlantPreparationError(:illumination, :combination,
            "coherently combined illumination requires a coherent ElectricField payload"))
    return CoherentFieldCombination()
end

function _require_illumination_combination(
    ::IncoherentIntensityAddition, payload::IntensityMap)
    typeof(payload.metadata.coherence) === IncoherentIntensityAddition ||
        throw(PlantPreparationError(:illumination, :combination,
            "incoherently combined illumination requires an incoherent IntensityMap payload"))
    return IncoherentIntensityAddition()
end

function _require_illumination_combination(
    ::_UnsupportedIlluminationCombination, evaluator, payload)
    throw(PlantPreparationError(:illumination, :missing_combination,
        "illumination evaluator type $(typeof(evaluator)) has not declared combination semantics"))
end

function _require_illumination_combination(combination, payload)
    throw(PlantPreparationError(:illumination, :combination,
        "illumination combination $(typeof(combination)) is incompatible with payload $(typeof(payload))"))
end

@inline function _prepared_illumination_combination(
    combination::_UnsupportedIlluminationCombination, evaluator, payload)
    return _require_illumination_combination(combination, evaluator, payload)
end

@inline function _prepared_illumination_combination(combination, evaluator,
    payload)
    return _require_illumination_combination(combination, payload)
end

@inline function _prepared_illumination_combination(evaluator, payload)
    return _prepared_illumination_combination(
        illumination_combination(evaluator), evaluator, payload)
end

"""Qualified cold preparation seam for an illumination evaluator."""
function prepare_illumination_evaluator(definition, destination, boundary)
    throw(PlantPreparationError(:illumination, :unsupported_evaluator,
        "illumination definition type $(typeof(definition)) does not implement prepare_illumination_evaluator"))
end

"""Qualified preparation-time evaluator binding-validation seam."""
function validate_illumination_evaluator_binding(evaluator, destination,
    boundary, contract)
    throw(PlantPreparationError(:illumination,
        :unsupported_evaluator_validation,
        "illumination evaluator type $(typeof(evaluator)) does not validate its prepared binding"))
end

"""Qualified mutating evaluator seam receiving explicit plant time and RNG."""
function evaluate_illumination!(destination, evaluator, model_time,
    rng::AbstractRNG)
    throw(PlantPreparationError(:illumination,
        :unsupported_evaluator_execution,
        "illumination evaluator type $(typeof(evaluator)) does not implement evaluate_illumination!"))
end

"""
    PreparedIlluminationEntry

Run-immutable typed entry boundary, evaluator, destination binding, combination
declaration, downstream-visibility description, and payload contract. The
evaluator may contain a separate mutable single-writer state/workspace object;
the evaluator wrapper itself must be immutable.
"""
struct PreparedIlluminationEntry{B,E,D,C,V,K}
    boundary::B
    evaluator::E
    destination::D
    combination::C
    visibility::V
    contract::K
end

function Base.getproperty(entry::PreparedIlluminationEntry, name::Symbol)
    name === :visibility && return deepcopy(getfield(entry, :visibility))
    return getfield(entry, name)
end

@inline illumination_entry_boundary(entry::PreparedIlluminationEntry) =
    entry.boundary
@inline illumination_evaluator(entry::PreparedIlluminationEntry) =
    entry.evaluator
@inline illumination_destination(entry::PreparedIlluminationEntry) =
    entry.destination
@inline illumination_combination(entry::PreparedIlluminationEntry) =
    entry.combination
@inline illumination_visibility(entry::PreparedIlluminationEntry) =
    deepcopy(getfield(entry, :visibility))

function PreparedIlluminationEntry(boundary, evaluator, destination;
    visibility)
    isnothing(visibility) && throw(PlantPreparationError(
        :illumination, :visibility,
        "illumination downstream visibility must be declared"))
    ismutabletype(typeof(evaluator)) && throw(PlantPreparationError(
        :illumination, :mutable_evaluator,
        "prepared illumination evaluator wrappers must be immutable and retain mutable state separately"))
    validate_illumination_entry_payload(boundary, destination)
    contract = illumination_payload_contract(destination)
    combination = _prepared_illumination_combination(evaluator,
        destination)
    visibility_snapshot = deepcopy(visibility)
    validate_illumination_evaluator_binding(evaluator, destination,
        boundary, contract)
    return PreparedIlluminationEntry(boundary, evaluator, destination,
        combination, visibility_snapshot, contract)
end

"""Prepare and bind one illumination evaluator to an exact typed destination."""
function prepare_illumination_entry(definition, destination, boundary;
    visibility)
    ismutabletype(typeof(definition)) && throw(PlantPreparationError(
        :illumination, :mutable_definition,
        "illumination evaluator definitions must be immutable parameter values"))
    evaluator = prepare_illumination_evaluator(definition, destination,
        boundary)
    return PreparedIlluminationEntry(boundary, evaluator, destination;
        visibility)
end

function validate_illumination_entry_binding(
    entry::PreparedIlluminationEntry, destination)
    getfield(entry, :destination) === destination || throw(
        PlantPreparationError(:illumination, :prepared_binding,
            "illumination entry does not retain the exact destination"))
    validate_illumination_entry_payload(getfield(entry, :boundary),
        destination)
    validate_illumination_payload_contract(destination,
        getfield(entry, :contract))
    current_combination = _prepared_illumination_combination(
        getfield(entry, :evaluator), destination)
    typeof(current_combination) ===
        typeof(getfield(entry, :combination)) || throw(
        PlantPreparationError(:illumination, :combination,
            "illumination evaluator combination semantics changed"))
    validate_illumination_evaluator_binding(getfield(entry, :evaluator),
        destination, getfield(entry, :boundary), getfield(entry, :contract))
    return entry
end

"""Evaluate one prepared entry at explicit plant time into its exact payload."""
function evaluate_illumination!(entry::PreparedIlluminationEntry,
    model_time::Real, rng::AbstractRNG)
    isfinite(model_time) || throw(PlantPreparationError(
        :illumination, :model_time,
        "illumination evaluator model time must be finite"))
    destination = getfield(entry, :destination)
    validate_illumination_entry_binding(entry, destination)
    result = evaluate_illumination!(destination,
        getfield(entry, :evaluator), model_time, rng)
    result === destination || throw(PlantPreparationError(
        :illumination, :evaluator_result,
        "illumination evaluator must return its exact caller-owned destination"))
    validate_illumination_payload_contract(destination,
        getfield(entry, :contract))
    return destination
end

@inline additional_path_materialization_rng_owner_roles(
    ::PreparedIlluminationEntry) = (:illumination,)

function validate_path_materialization_binding(
    entry::PreparedIlluminationEntry, input, atmosphere, source)
    validate_illumination_entry_binding(entry, input)
    return nothing
end

function validate_path_materialization(entry::PreparedIlluminationEntry,
    input, atmosphere, epoch::AtmosphereEpoch)
    validate_illumination_entry_binding(entry, input)
    return input
end

function materialize_path_input!(entry::PreparedIlluminationEntry, input,
    atmosphere, epoch::AtmosphereEpoch)
    throw(PlantPreparationError(:illumination, :rng_owner,
        "prepared illumination materialization requires its path-owned or an explicit RNG"))
end

function materialize_path_input!(entry::PreparedIlluminationEntry, input,
    atmosphere, epoch::AtmosphereEpoch, rng::AbstractRNG)
    validate_illumination_entry_binding(entry, input)
    return evaluate_illumination!(entry, epoch_time(epoch), rng)
end

function materialize_path_input_rngs!(entry::PreparedIlluminationEntry,
    input, atmosphere, epoch::AtmosphereEpoch, rngs::PreparedOwnerRNGs)
    validate_illumination_entry_binding(entry, input)
    rng = rng_stream_state(rngs, Val(:illumination))
    return evaluate_illumination!(entry, epoch_time(epoch), rng)
end

"""
    UniformIntensityIllumination(value; combination)

Native immutable definition for a spatially uniform intensity entry. `value`
is interpreted in the destination `IntensityMap`'s explicitly declared
normalization and spatial measure. Combination semantics are mandatory.
"""
struct UniformIntensityIllumination{T<:Real,C}
    value::T
    combination::C

    function UniformIntensityIllumination(value::T,
        combination::C) where {T<:Real,C}
        isfinite(value) && value >= zero(value) || throw(
            InvalidConfiguration(
                "uniform illumination value must be finite and nonnegative"))
        return new{T,C}(value, combination)
    end
end

UniformIntensityIllumination(value::Real; combination) =
    UniformIntensityIllumination(value, combination)

struct PreparedUniformIntensityIllumination{
    T<:AbstractFloat,C,B<:AbstractArrayBackend,D<:AbstractPlaneDevice}
    value::T
    combination::C
    backend::B
    device::D
end

@inline illumination_combination(
    evaluator::PreparedUniformIntensityIllumination) = evaluator.combination

function prepare_illumination_evaluator(
    definition::UniformIntensityIllumination,
    destination::IntensityMap, boundary)
    T = eltype(destination.values)
    value = T(definition.value)
    isfinite(value) || throw(PlantPreparationError(
        :illumination, :numeric_type,
        "uniform illumination value is not representable in the destination numeric type"))
    return PreparedUniformIntensityIllumination(value,
        definition.combination, backend(destination),
        plane_device(destination.values))
end

function validate_illumination_evaluator_binding(
    evaluator::PreparedUniformIntensityIllumination,
    destination::IntensityMap, boundary, contract)
    validate_illumination_payload_contract(destination, contract)
    typeof(backend(destination)) === typeof(evaluator.backend) || throw(
        PlantPreparationError(:illumination, :backend,
            "uniform illumination evaluator and destination backends differ"))
    plane_device(destination.values) == evaluator.device || throw(
        PlantPreparationError(:illumination, :device,
            "uniform illumination evaluator and destination devices differ"))
    return nothing
end

@inline function evaluate_illumination!(destination::IntensityMap,
    evaluator::PreparedUniformIntensityIllumination, model_time,
    rng::AbstractRNG)
    fill!(destination.values, evaluator.value)
    return destination
end

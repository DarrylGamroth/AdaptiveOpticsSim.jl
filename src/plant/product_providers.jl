#
# Prepared acquisition-product providers
#
# Provider selection is a cold, per-acquisition decision. Every provider writes
# into one caller-owned AcquisitionProducts value and returns that exact value.
# The prepared product contract is shared by full-optical, reduced-order, and
# synthetic/replay implementations.
#

"""Trait value for a provider that evaluates the prepared optical path."""
struct FullOpticalProviderStyle end

"""Trait value for a command-responsive reduced-order AO provider."""
struct CommandResponsiveReducedOrderProviderStyle end

"""Trait value for a nonresponsive synthetic or completed-product replay provider."""
struct SyntheticReplayProviderStyle end

struct _UnsupportedAcquisitionProviderStyle end

"""
    acquisition_provider_style(::Type{T})

Trait for a prepared acquisition-provider implementation. Implementations must
return one of `FullOpticalProviderStyle()`,
`CommandResponsiveReducedOrderProviderStyle()`, or
`SyntheticReplayProviderStyle()`.
"""
acquisition_provider_style(::Type) = _UnsupportedAcquisitionProviderStyle()

@inline acquisition_provider_style(provider) =
    acquisition_provider_style(typeof(provider))

@inline _require_acquisition_provider_style(style::Union{
    FullOpticalProviderStyle,
    CommandResponsiveReducedOrderProviderStyle,
    SyntheticReplayProviderStyle,
}, provider) = style

function _require_acquisition_provider_style(
    ::_UnsupportedAcquisitionProviderStyle, provider)
    throw(PlantPreparationError(:acquisition, :unsupported_provider,
        "acquisition provider type $(typeof(provider)) has not declared an acquisition-provider style"))
end

function _require_acquisition_provider_style(style, provider)
    throw(PlantPreparationError(:acquisition, :invalid_provider_style,
        "acquisition_provider_style($(typeof(provider))) returned unsupported value $(typeof(style))"))
end

@inline function _prepared_acquisition_provider_style(provider)
    style = acquisition_provider_style(typeof(provider))
    return _require_acquisition_provider_style(style, provider)
end

"""
    acquisition_provider_payload_work(::Type{T})

Cold introspection seam declaring the principal payload work performed by a
prepared provider. Extension providers return a nonempty `Symbol`, such as
`:regenerate_reduced_order_product`. The declaration describes work; it is not
performance evidence.
"""
function acquisition_provider_payload_work(::Type{T}) where {T}
    throw(PlantPreparationError(:acquisition, :missing_payload_work,
        "acquisition provider type $T must declare its payload work"))
end

@inline acquisition_provider_payload_work(provider) =
    acquisition_provider_payload_work(typeof(provider))

const _EMPTY_PROVIDER_PAYLOAD_WORK = Symbol("")

function _require_provider_payload_work_value(work::Symbol, ::Any)
    work === _EMPTY_PROVIDER_PAYLOAD_WORK && throw(PlantPreparationError(
        :acquisition, :invalid_payload_work,
        "acquisition provider payload work must not be empty"))
    return work
end

function _require_provider_payload_work_value(work, provider)
    throw(PlantPreparationError(:acquisition,
        :invalid_payload_work,
        "acquisition provider payload work must be a Symbol; got $(typeof(work))"))
end

function _require_provider_payload_work(provider)
    work = acquisition_provider_payload_work(provider)
    return _require_provider_payload_work_value(work, provider)
end

"""
    AcquisitionProducts(observation, measurement; metadata)
    AcquisitionProducts(observation; metadata)

Caller-owned observation and optional measurement storage plus run-immutable
logical product metadata. `metadata` is deliberately application-extensible;
it must declare any geometry, radiometry, units, layout, or semantic labels not
already carried by the typed products.
"""
struct AcquisitionProducts{O,M,D}
    observation::O
    measurement::M
    metadata::D

    function AcquisitionProducts(observation::O, measurement::M,
        metadata::D) where {O,M,D}
        isnothing(metadata) && throw(PlantPreparationError(:acquisition,
            :metadata, "acquisition product metadata must be declared"))
        return new{O,M,D}(observation, measurement, metadata)
    end
end

@inline AcquisitionProducts(observation; metadata) =
    AcquisitionProducts(observation, nothing, metadata)

@inline AcquisitionProducts(observation, measurement; metadata) =
    AcquisitionProducts(observation, measurement, metadata)

@inline acquisition_product_metadata(products::AcquisitionProducts) =
    products.metadata

struct _ArrayProductContract{N,E,B<:AbstractArrayBackend,
    D<:AbstractPlaneDevice}
    dimensions::NTuple{N,Int}
    numeric_type::Type{E}
    backend::B
    device::D
end

struct _RefProductContract{E}
    numeric_type::Type{E}
end

struct _OpticalProductContract{M,C}
    metadata::M
    storage::C
end

struct _WFSObservationProductContract{U,M,C}
    units::U
    metadata::M
    storage::C
end

struct _WFSMeasurementProductContract{U,M,C}
    units::U
    metadata::M
    storage::C
end

struct _AbsentProductContract end

"""
    AcquisitionProductContract

Run-immutable compatibility snapshot for one `AcquisitionProducts` value. It
captures observation and measurement shape, numeric type, memory domain,
typed metadata and units, plus the acquisition-level metadata supplied by the
caller.
"""
struct AcquisitionProductContract{O,M,D}
    observation::O
    measurement::M
    metadata::D
end

@inline function _array_product_contract(storage::AbstractArray)
    selector = backend(storage)
    device = plane_device(storage)
    return _ArrayProductContract(size(storage), eltype(storage), selector,
        device)
end

@inline acquisition_product_contract(product::AbstractArray) =
    _array_product_contract(product)

@inline acquisition_product_contract(product::Base.RefValue{T}) where {T} =
    _RefProductContract(T)

function acquisition_product_contract(product::IntensityMap)
    validate_plane_storage(product.metadata, product.values;
        label="acquisition intensity product")
    return _OpticalProductContract(deepcopy(product.metadata),
        _array_product_contract(product.values))
end

function acquisition_product_contract(product::WFSObservation)
    validate_wfs_observation(product)
    return _WFSObservationProductContract(deepcopy(product.units),
        deepcopy(product.metadata),
        acquisition_product_contract(product.storage))
end

function acquisition_product_contract(product::WFSMeasurement)
    validate_wfs_measurement(product)
    return _WFSMeasurementProductContract(deepcopy(product.units),
        deepcopy(product.metadata),
        acquisition_product_contract(product.storage))
end

@inline acquisition_product_contract(::Nothing) = _AbsentProductContract()

@inline acquisition_product_contract(::Tuple{}) = ()

@inline function acquisition_product_contract(products::Tuple)
    return (acquisition_product_contract(first(products)),
        acquisition_product_contract(Base.tail(products))...)
end

function acquisition_product_contract(product)
    throw(PlantPreparationError(:acquisition, :unsupported_product,
        "acquisition product type $(typeof(product)) must implement an acquisition product contract"))
end

@inline function _validate_product_metadata(value, expected,
    label::AbstractString)
    typeof(value) === typeof(expected) || throw(PlantPreparationError(
        :acquisition, :metadata,
        "$label metadata type does not match its prepared product contract"))
    isequal(value, expected) || throw(PlantPreparationError(
        :acquisition, :metadata,
        "$label metadata does not match its prepared product contract"))
    return value
end

@inline function _validate_product_units(value, expected,
    label::AbstractString)
    typeof(value) === typeof(expected) || throw(PlantPreparationError(
        :acquisition, :units,
        "$label units type does not match its prepared product contract"))
    isequal(value, expected) || throw(PlantPreparationError(
        :acquisition, :units,
        "$label units do not match its prepared product contract"))
    return value
end

"""Return a defensive compatibility snapshot for caller-owned products."""
function acquisition_product_contract(products::AcquisitionProducts)
    return AcquisitionProductContract(
        acquisition_product_contract(products.observation),
        acquisition_product_contract(products.measurement),
        deepcopy(products.metadata),
    )
end

@inline function _validate_array_product_contract(storage::AbstractArray,
    contract::_ArrayProductContract, label::AbstractString)
    size(storage) == contract.dimensions || throw(PlantPreparationError(
        :acquisition, :shape,
        "$label shape does not match its prepared product contract"))
    eltype(storage) === contract.numeric_type || throw(
        PlantPreparationError(:acquisition, :numeric_type,
            "$label element type does not match its prepared product contract"))
    typeof(backend(storage)) === typeof(contract.backend) || throw(
        PlantPreparationError(:acquisition, :backend,
            "$label backend does not match its prepared product contract"))
    plane_device(storage) == contract.device || throw(
        PlantPreparationError(:acquisition, :device,
            "$label device does not match its prepared product contract"))
    return storage
end

@inline function validate_acquisition_product_contract(
    product::AbstractArray, contract::_ArrayProductContract,
    label::AbstractString)
    return _validate_array_product_contract(product, contract, label)
end

@inline function validate_acquisition_product_contract(
    product::Base.RefValue{T}, contract::_RefProductContract{T},
    ::AbstractString) where {T}
    return product
end

function validate_acquisition_product_contract(product::IntensityMap,
    contract::_OpticalProductContract, label::AbstractString)
    validate_plane_storage(product.metadata, product.values; label=label)
    _validate_product_metadata(product.metadata, contract.metadata, label)
    _validate_array_product_contract(product.values, contract.storage,
        label)
    return product
end

function validate_acquisition_product_contract(product::WFSObservation,
    contract::_WFSObservationProductContract, label::AbstractString)
    _validate_product_units(product.units, contract.units, label)
    _validate_product_metadata(product.metadata, contract.metadata, label)
    validate_acquisition_product_contract(product.storage,
        contract.storage, label)
    return product
end

function validate_acquisition_product_contract(product::WFSMeasurement,
    contract::_WFSMeasurementProductContract, label::AbstractString)
    _validate_product_units(product.units, contract.units, label)
    _validate_product_metadata(product.metadata, contract.metadata, label)
    validate_acquisition_product_contract(product.storage,
        contract.storage, label)
    return product
end

@inline validate_acquisition_product_contract(::Nothing,
    ::_AbsentProductContract, ::AbstractString) = nothing

@inline validate_acquisition_product_contract(::Tuple{}, ::Tuple{},
    ::AbstractString) = nothing

function validate_acquisition_product_contract(::Tuple{}, ::Tuple,
    label::AbstractString)
    throw(PlantPreparationError(:acquisition, :shape,
        "$label tuple is shorter than its prepared product contract"))
end

function validate_acquisition_product_contract(::Tuple, ::Tuple{},
    label::AbstractString)
    throw(PlantPreparationError(:acquisition, :shape,
        "$label tuple is longer than its prepared product contract"))
end

@inline function validate_acquisition_product_contract(products::Tuple,
    contracts::Tuple, label::AbstractString)
    validate_acquisition_product_contract(first(products),
        first(contracts), label)
    return validate_acquisition_product_contract(Base.tail(products),
        Base.tail(contracts), label)
end

function validate_acquisition_product_contract(product, contract,
    label::AbstractString)
    throw(PlantPreparationError(:acquisition, :product_type,
        "$label type $(typeof(product)) is incompatible with prepared contract $(typeof(contract))"))
end

function validate_acquisition_product_contract(
    products::AcquisitionProducts,
    contract::AcquisitionProductContract)
    validate_acquisition_product_contract(products.observation,
        contract.observation, "acquisition observation")
    validate_acquisition_product_contract(products.measurement,
        contract.measurement, "acquisition measurement")
    _validate_product_metadata(products.metadata, contract.metadata,
        "acquisition")
    return products
end

@inline function copy_acquisition_product!(destination::AbstractArray,
    source::AbstractArray)
    copyto!(destination, source)
    return destination
end

@inline function copy_acquisition_product!(destination::Base.RefValue,
    source::Base.RefValue)
    destination[] = source[]
    return destination
end

@inline function copy_acquisition_product!(destination::IntensityMap,
    source::IntensityMap)
    copyto!(destination.values, source.values)
    return destination
end

@inline function copy_acquisition_product!(destination::WFSObservation,
    source::WFSObservation)
    copy_acquisition_product!(destination.storage, source.storage)
    return destination
end

@inline function copy_acquisition_product!(destination::WFSMeasurement,
    source::WFSMeasurement)
    copy_acquisition_product!(destination.storage, source.storage)
    return destination
end

@inline copy_acquisition_product!(::Nothing, ::Nothing) = nothing
@inline copy_acquisition_product!(::Tuple{}, ::Tuple{}) = nothing

@inline function copy_acquisition_product!(destination::Tuple,
    source::Tuple)
    copy_acquisition_product!(first(destination), first(source))
    copy_acquisition_product!(Base.tail(destination), Base.tail(source))
    return destination
end

function copy_acquisition_product!(destination, source)
    throw(PlantPreparationError(:acquisition, :unsupported_product_copy,
        "cannot copy acquisition product $(typeof(source)) into $(typeof(destination))"))
end

@inline function _copy_acquisition_products!(
    destination::AcquisitionProducts, source::AcquisitionProducts)
    copy_acquisition_product!(destination.observation, source.observation)
    copy_acquisition_product!(destination.measurement, source.measurement)
    return destination
end

"""
    PreparedAcquisitionProvider(implementation, products)

Prepared provider implementation, exact caller-owned destination products,
and their run-immutable compatibility contract. Provider implementations are
selected during preparation and cannot be replaced in a prepared plant.
"""
struct PreparedAcquisitionProvider{I,P<:AcquisitionProducts,
    C<:AcquisitionProductContract,S}
    implementation::I
    products::P
    contract::C
    style::S
    payload_work::Symbol
end

function PreparedAcquisitionProvider(implementation,
    products::AcquisitionProducts)
    style = _prepared_acquisition_provider_style(implementation)
    payload_work = _require_provider_payload_work(implementation)
    contract = acquisition_product_contract(products)
    validate_acquisition_product_contract(products, contract)
    return PreparedAcquisitionProvider(implementation, products, contract,
        style, payload_work)
end

@inline acquisition_provider_style(provider::PreparedAcquisitionProvider) =
    provider.style
@inline acquisition_provider_payload_work(
    provider::PreparedAcquisitionProvider) = provider.payload_work
@inline acquisition_products(provider::PreparedAcquisitionProvider) =
    provider.products
@inline acquisition_product_contract(provider::PreparedAcquisitionProvider) =
    provider.contract
@inline acquisition_product_metadata(provider::PreparedAcquisitionProvider) =
    provider.products.metadata

"""Prepared acquisition wrapper that requires an ordinary full-optical path."""
struct PreparedFullOpticalProvider{E}
    execution::E
end

@inline acquisition_provider_style(
    ::Type{<:PreparedFullOpticalProvider}) = FullOpticalProviderStyle()
@inline acquisition_provider_payload_work(
    ::Type{<:PreparedFullOpticalProvider}) = :full_optical_evaluation
@inline additional_acquisition_rng_owner_roles(
    provider::PreparedFullOpticalProvider) =
    additional_acquisition_rng_owner_roles(provider.execution)

"""Prepare an ordinary full-optical acquisition execution behind the provider seam."""
function prepare_full_optical_provider(execution,
    products::AcquisitionProducts)
    return PreparedAcquisitionProvider(
        PreparedFullOpticalProvider(execution), products)
end

"""Prepared nonresponsive provider that republishes unchanged destination contents."""
struct PreparedUnchangedSyntheticProvider end

@inline acquisition_provider_style(
    ::Type{PreparedUnchangedSyntheticProvider}) = SyntheticReplayProviderStyle()
@inline acquisition_provider_payload_work(
    ::Type{PreparedUnchangedSyntheticProvider}) = :reuse_unchanged

"""Prepare a nonresponsive provider that performs no payload writes."""
function prepare_unchanged_synthetic_provider(
    products::AcquisitionProducts)
    return PreparedAcquisitionProvider(
        PreparedUnchangedSyntheticProvider(), products)
end

"""Prepared nonresponsive provider that copies one owned product snapshot."""
struct PreparedCopiedSyntheticProvider{P<:AcquisitionProducts}
    source::P
end

@inline acquisition_provider_style(
    ::Type{<:PreparedCopiedSyntheticProvider}) = SyntheticReplayProviderStyle()
@inline acquisition_provider_payload_work(
    ::Type{<:PreparedCopiedSyntheticProvider}) = :copy_prepared_product

"""Prepare a nonresponsive provider that copies one product on every event."""
function prepare_copied_synthetic_provider(products::AcquisitionProducts,
    source::AcquisitionProducts)
    contract = acquisition_product_contract(products)
    validate_acquisition_product_contract(source, contract)
    owned_source = deepcopy(source)
    validate_acquisition_product_contract(owned_source, contract)
    implementation = PreparedCopiedSyntheticProvider(owned_source)
    return PreparedAcquisitionProvider(implementation, products)
end

"""Immutable cyclic replay configuration."""
struct CyclicReplayProviderParams
    corpus_size::Int
end

"""Single-writer cyclic replay cursor."""
mutable struct CyclicReplayProviderState
    next_index::Int
end

"""Prepared nonresponsive bounded cyclic completed-product replay."""
struct PreparedCyclicReplayProvider{P<:AcquisitionProducts,
    M<:Memory{P},S<:CyclicReplayProviderState}
    params::CyclicReplayProviderParams
    corpus::M
    state::S
end

@inline acquisition_provider_style(
    ::Type{<:PreparedCyclicReplayProvider}) = SyntheticReplayProviderStyle()
@inline acquisition_provider_payload_work(
    ::Type{<:PreparedCyclicReplayProvider}) = :bounded_cyclic_replay

function _owned_replay_corpus(corpus::Union{Tuple,AbstractVector})
    isempty(corpus) && throw(PlantPreparationError(:acquisition,
        :empty_replay, "cyclic replay corpus must not be empty"))
    first_product = first(corpus)
    _require_replay_product(first_product)
    P = typeof(first_product)
    owned = Memory{P}(undef, length(corpus))
    for (index, product) in enumerate(corpus)
        typeof(product) === P || throw(PlantPreparationError(:acquisition,
            :product_type,
            "cyclic replay corpus must have one concrete product type"))
        owned[index] = deepcopy(product)
    end
    return owned
end

function _owned_replay_corpus(corpus)
    throw(PlantPreparationError(:acquisition, :invalid_replay,
        "cyclic replay corpus must be a finite Tuple or AbstractVector; got $(typeof(corpus))"))
end

@inline _require_replay_product(product::AcquisitionProducts) = product

function _require_replay_product(product)
    throw(PlantPreparationError(:acquisition, :product_type,
        "cyclic replay corpus must contain AcquisitionProducts values"))
end

"""Prepare a nonresponsive bounded completed-product replay that cycles."""
function prepare_cyclic_replay_provider(products::AcquisitionProducts,
    corpus)
    destination_contract = acquisition_product_contract(products)
    owned = _owned_replay_corpus(corpus)
    @inbounds for source in owned
        validate_acquisition_product_contract(source, destination_contract)
    end
    implementation = PreparedCyclicReplayProvider(
        CyclicReplayProviderParams(length(owned)), owned,
        CyclicReplayProviderState(1))
    return PreparedAcquisitionProvider(implementation, products)
end

"""Qualified preparation-time binding-validation seam for providers."""
function validate_acquisition_provider_binding(implementation, path_result,
    products::AcquisitionProducts, contract::AcquisitionProductContract)
    throw(PlantPreparationError(:acquisition,
        :unsupported_provider_validation,
        "acquisition provider type $(typeof(implementation)) does not validate its path/product binding"))
end

function validate_acquisition_provider_binding(
    provider::PreparedAcquisitionProvider, path_result)
    validate_acquisition_product_contract(provider.products,
        provider.contract)
    style = _prepared_acquisition_provider_style(provider.implementation)
    typeof(style) === typeof(provider.style) || throw(
        PlantPreparationError(:acquisition, :provider_style,
            "prepared provider style changed"))
    acquisition_provider_payload_work(provider.implementation) ===
        provider.payload_work || throw(PlantPreparationError(:acquisition,
            :payload_work,
            "prepared provider payload-work declaration changed"))
    return validate_acquisition_provider_binding(provider.implementation,
        path_result, provider.products, provider.contract)
end

function validate_acquisition_provider_binding(
    implementation::PreparedUnchangedSyntheticProvider, path_result,
    products::AcquisitionProducts, contract::AcquisitionProductContract)
    validate_acquisition_product_contract(products, contract)
    return nothing
end

function validate_acquisition_provider_binding(
    implementation::PreparedCopiedSyntheticProvider, path_result,
    products::AcquisitionProducts, contract::AcquisitionProductContract)
    validate_acquisition_product_contract(products, contract)
    validate_acquisition_product_contract(implementation.source, contract)
    return nothing
end

function validate_acquisition_provider_binding(
    implementation::PreparedCyclicReplayProvider, path_result,
    products::AcquisitionProducts, contract::AcquisitionProductContract)
    validate_acquisition_product_contract(products, contract)
    implementation.params.corpus_size == length(implementation.corpus) ||
        throw(PlantPreparationError(:acquisition, :replay_topology,
            "prepared cyclic replay corpus size changed"))
    index = implementation.state.next_index
    1 <= index <= implementation.params.corpus_size || throw(
        PlantPreparationError(:acquisition, :replay_cursor,
            "prepared cyclic replay cursor is out of bounds"))
    source = @inbounds implementation.corpus[index]
    validate_acquisition_product_contract(source, contract)
    return nothing
end

"""Qualified mutating execution seam for prepared providers."""
function execute_acquisition_provider!(products::AcquisitionProducts,
    path_result, implementation, rngs)
    throw(PlantPreparationError(:acquisition, :unsupported_provider_execution,
        "acquisition provider type $(typeof(implementation)) does not implement execute_acquisition_provider!"))
end

@inline function execute_acquisition_provider!(
    products::AcquisitionProducts, path_result,
    ::PreparedUnchangedSyntheticProvider, rngs)
    return products
end

@inline function execute_acquisition_provider!(
    products::AcquisitionProducts, path_result,
    implementation::PreparedCopiedSyntheticProvider, rngs)
    return _copy_acquisition_products!(products, implementation.source)
end

@inline function execute_acquisition_provider!(
    products::AcquisitionProducts, path_result,
    implementation::PreparedCyclicReplayProvider, rngs)
    index = implementation.state.next_index
    source = @inbounds implementation.corpus[index]
    _copy_acquisition_products!(products, source)
    implementation.state.next_index = ifelse(
        index == implementation.params.corpus_size, 1, index + 1)
    return products
end

function execute_acquisition_provider!(provider::PreparedAcquisitionProvider,
    path_result, rngs)
    validate_acquisition_provider_binding(provider, path_result)
    result = execute_acquisition_provider!(provider.products, path_result,
        provider.implementation, rngs)
    result === provider.products || throw(PlantPreparationError(
        :acquisition, :provider_result,
        "acquisition provider must return its exact caller-owned AcquisitionProducts destination"))
    return result
end

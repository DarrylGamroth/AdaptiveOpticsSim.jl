struct ProviderFullModel{T<:AbstractFloat}
    exposure::T
end

struct ProviderReducedModel{T<:AbstractFloat}
    disturbance::T
end

struct ProviderUnchangedModel{T<:AbstractFloat}
    value::T
end

struct ProviderCopiedModel{T<:AbstractFloat}
    value::T
end

struct ProviderReplayModel{T<:AbstractFloat}
    first_value::T
    second_value::T
end

mutable struct ProviderReducedState{T<:AbstractFloat}
    disturbance::T
    command::T
end

struct TestCommandResponsiveReducedProvider{S<:ProviderReducedState}
    state::S
end

struct UnsupportedProviderStyle end
struct InvalidProviderStyle end
struct MissingPayloadWorkProvider end
struct EmptyPayloadWorkProvider end
struct InvalidPayloadWorkProvider end
struct UnsupportedProviderValidation end
struct UnsupportedProviderExecution end
struct InvalidProviderResult end

for model in (
    ProviderFullModel,
    ProviderReducedModel,
    ProviderUnchangedModel,
    ProviderCopiedModel,
    ProviderReplayModel,
)
    @eval AdaptiveOpticsSim.plant_model_definition_style(
        ::Type{<:$model}) = ColdPlantModelDefinition()
end

AdaptiveOpticsSim.acquisition_provider_style(
    ::Type{<:TestCommandResponsiveReducedProvider}) =
    CommandResponsiveReducedOrderProviderStyle()
AdaptiveOpticsSim.acquisition_provider_payload_work(
    ::Type{<:TestCommandResponsiveReducedProvider}) =
    :regenerate_reduced_order_product

AdaptiveOpticsSim.acquisition_provider_style(
    ::Type{InvalidProviderStyle}) = :invalid
AdaptiveOpticsSim.acquisition_provider_style(
    ::Type{MissingPayloadWorkProvider}) = SyntheticReplayProviderStyle()
AdaptiveOpticsSim.acquisition_provider_style(
    ::Type{EmptyPayloadWorkProvider}) = SyntheticReplayProviderStyle()
AdaptiveOpticsSim.acquisition_provider_payload_work(
    ::Type{EmptyPayloadWorkProvider}) = Symbol("")
AdaptiveOpticsSim.acquisition_provider_style(
    ::Type{InvalidPayloadWorkProvider}) = SyntheticReplayProviderStyle()
AdaptiveOpticsSim.acquisition_provider_payload_work(
    ::Type{InvalidPayloadWorkProvider}) = "not-a-symbol"
AdaptiveOpticsSim.acquisition_provider_style(
    ::Type{UnsupportedProviderValidation}) = SyntheticReplayProviderStyle()
AdaptiveOpticsSim.acquisition_provider_payload_work(
    ::Type{UnsupportedProviderValidation}) = :unsupported_validation
AdaptiveOpticsSim.acquisition_provider_style(
    ::Type{UnsupportedProviderExecution}) = SyntheticReplayProviderStyle()
AdaptiveOpticsSim.acquisition_provider_payload_work(
    ::Type{UnsupportedProviderExecution}) = :unsupported_execution
AdaptiveOpticsSim.acquisition_provider_style(
    ::Type{InvalidProviderResult}) =
    CommandResponsiveReducedOrderProviderStyle()
AdaptiveOpticsSim.acquisition_provider_payload_work(
    ::Type{InvalidProviderResult}) = :regenerate_invalid_result

function provider_test_metadata(path::PreparedPathExecutor)
    return (
        kind=:detector_frame,
        units=:detected_electrons,
        geometry=path.result.metadata,
        radiometry=(exposure_s=0.5, response=:ideal_frame),
        semantics=:complete_acquisition,
    )
end

function provider_test_products(path::PreparedPathExecutor, value)
    storage = similar(path.result.values)
    fill!(storage, convert(eltype(storage), value))
    return AcquisitionProducts(storage;
        metadata=provider_test_metadata(path))
end

function AdaptiveOpticsSim.prepare_acquisition_provider(
    model::ProviderFullModel,
    ::AcquisitionDefinition,
    path::PreparedPathExecutor,
)
    require_path_result(path)
    T = eltype(path.result.values)
    detector = Detector(integration_time=T(model.exposure),
        noise=NoiseNone(), qe=one(T), response_model=NullFrameResponse(),
        T=T, backend=path.key.backend)
    execution = FrameAcquisitionExecution(detector, path.result)
    products = AcquisitionProducts(execution.observation;
        metadata=provider_test_metadata(path))
    return prepare_full_optical_provider(execution, products)
end

function AdaptiveOpticsSim.prepare_acquisition_provider(
    model::ProviderReducedModel,
    ::AcquisitionDefinition,
    path::PreparedPathExecutor,
)
    require_path_result(path)
    products = provider_test_products(path, zero(model.disturbance))
    state = ProviderReducedState(model.disturbance,
        zero(model.disturbance))
    return PreparedAcquisitionProvider(
        TestCommandResponsiveReducedProvider(state), products)
end

function AdaptiveOpticsSim.validate_acquisition_provider_binding(
    ::TestCommandResponsiveReducedProvider,
    path_result,
    products::AcquisitionProducts,
)
    size(products.observation) == size(path_result.values) || throw(
        PlantPreparationError(:acquisition, :shape,
            "test reduced-order product shape must match its declared path"))
    return nothing
end

function AdaptiveOpticsSim.execute_acquisition_provider!(
    products::AcquisitionProducts,
    path_result,
    provider::TestCommandResponsiveReducedProvider,
    rngs,
)
    value = provider.state.disturbance - provider.state.command
    fill!(products.observation, value)
    return products
end

function AdaptiveOpticsSim.prepare_acquisition_provider(
    model::ProviderUnchangedModel,
    ::AcquisitionDefinition,
    path::PreparedPathExecutor,
)
    require_path_result(path)
    return prepare_unchanged_synthetic_provider(
        provider_test_products(path, model.value))
end

function AdaptiveOpticsSim.prepare_acquisition_provider(
    model::ProviderCopiedModel,
    ::AcquisitionDefinition,
    path::PreparedPathExecutor,
)
    require_path_result(path)
    destination = provider_test_products(path, zero(model.value))
    source = provider_test_products(path, model.value)
    return prepare_copied_synthetic_provider(destination, source)
end

function AdaptiveOpticsSim.prepare_acquisition_provider(
    model::ProviderReplayModel,
    ::AcquisitionDefinition,
    path::PreparedPathExecutor,
)
    require_path_result(path)
    destination = provider_test_products(path, zero(model.first_value))
    corpus = [provider_test_products(path, model.first_value),
        provider_test_products(path, model.second_value)]
    return prepare_cyclic_replay_provider(destination, corpus)
end

function AdaptiveOpticsSim.validate_acquisition_provider_binding(
    ::InvalidProviderResult,
    path_result,
    products::AcquisitionProducts,
)
    return nothing
end

AdaptiveOpticsSim.execute_acquisition_provider!(
    products::AcquisitionProducts, path_result,
    ::InvalidProviderResult, rngs) = nothing

function AdaptiveOpticsSim.validate_acquisition_provider_binding(
    ::UnsupportedProviderExecution,
    path_result,
    products::AcquisitionProducts,
)
    return nothing
end

@inline provider_test_consumer(products::AcquisitionProducts) =
    sum(products.observation)

function provider_test_owner_rngs(plant::PreparedPlant,
    owner::PreparedAcquisitionOwner)
    @inbounds for index in eachindex(plant.acquisitions)
        plant.acquisitions[index] === owner &&
            return plant.rngs.acquisitions[index]
    end
    error("test owner is not part of the prepared plant")
end

function provider_execution_allocations(owner, rngs)
    execute_acquisition!(owner, rngs)
    return @allocated execute_acquisition!(owner, rngs)
end

@testset "Invariant acquisition product-provider seam" begin
    T = Float64
    telescope = Telescope(resolution=8, diameter=T(4),
        central_obstruction=zero(T), T=T)
    atmosphere = KolmogorovAtmosphere(telescope; r0=T(0.2), L0=T(25), T=T)
    source = Source(band=:custom, wavelength=T(0.8e-6),
        photon_irradiance=T(3), T=T)
    path_definition = OpticalPathDefinition(:provider_path, source,
        CountedDirectSciencePathModel(2, UInt(41)))
    acquisitions = (
        AcquisitionDefinition(:full, :provider_path,
            ProviderFullModel(T(0.5))),
        AcquisitionDefinition(:reduced, :provider_path,
            ProviderReducedModel(T(4))),
        AcquisitionDefinition(:unchanged, :provider_path,
            ProviderUnchangedModel(T(5))),
        AcquisitionDefinition(:copied, :provider_path,
            ProviderCopiedModel(T(6))),
        AcquisitionDefinition(:replay, :provider_path,
            ProviderReplayModel(T(7), T(8))),
    )
    definition = PlantDefinition(; telescope, atmosphere,
        paths=(path_definition,), acquisitions)
    plant = prepare_plant(definition; run_seed=0x7300)
    path = prepared_path(plant, :provider_path)
    full = prepared_acquisition(plant, :full)
    reduced = prepared_acquisition(plant, :reduced)
    unchanged = prepared_acquisition(plant, :unchanged)
    copied = prepared_acquisition(plant, :copied)
    replay = prepared_acquisition(plant, :replay)

    @test acquisition_provider_style(full) isa FullOpticalProviderStyle
    @test acquisition_provider_style(reduced) isa
        CommandResponsiveReducedOrderProviderStyle
    for owner in (unchanged, copied, replay)
        @test acquisition_provider_style(owner) isa
            SyntheticReplayProviderStyle
    end
    @test acquisition_provider_payload_work(reduced) ===
        :regenerate_reduced_order_product
    @test acquisition_provider_payload_work(unchanged) === :reuse_unchanged
    @test acquisition_provider_payload_work(copied) ===
        :copy_prepared_product
    @test acquisition_provider_payload_work(replay) ===
        :bounded_cyclic_replay

    reduced_implementation = acquisition_provider(reduced).implementation
    provider_values = (
        acquisition_provider(full).implementation,
        reduced_implementation,
        acquisition_provider(unchanged).implementation,
        acquisition_provider(copied).implementation,
        acquisition_provider(replay).implementation,
    )
    @test acquisition_provider_style(provider_values[1]) isa
        FullOpticalProviderStyle
    @test acquisition_provider_style(provider_values[2]) isa
        CommandResponsiveReducedOrderProviderStyle
    for implementation in provider_values[3:5]
        @test acquisition_provider_style(implementation) isa
            SyntheticReplayProviderStyle
    end
    @test map(acquisition_provider_payload_work, provider_values) == (
        :full_optical_evaluation,
        :regenerate_reduced_order_product,
        :reuse_unchanged,
        :copy_prepared_product,
        :bounded_cyclic_replay,
    )
    @test acquisition_products(acquisition_provider(reduced)) ===
        acquisition_products(reduced)
    @test acquisition_product_metadata(acquisition_provider(reduced)) ===
        acquisition_product_metadata(acquisition_products(reduced))

    full_contract = acquisition_product_contract(full)
    for owner in (reduced, unchanged, copied, replay)
        @test validate_acquisition_product_contract(
            acquisition_products(owner), full_contract) ===
            acquisition_products(owner)
    end

    nonoptical_selection = prepare_acquisition_selection(plant,
        (:reduced, :unchanged, :copied, :replay))
    @test prepared_paths(nonoptical_selection) === ()
    @test path.execution.executions[] == 0
    @test execute_acquisition_selection_at!(nonoptical_selection,
        T(0.01)) === nonoptical_selection
    @test path.execution.executions[] == 0

    pixel_count = length(acquisition_observation(reduced))
    @test provider_test_consumer(acquisition_products(reduced)) ==
        T(4) * pixel_count
    @test provider_test_consumer(acquisition_products(unchanged)) ==
        T(5) * pixel_count
    @test provider_test_consumer(acquisition_products(copied)) ==
        T(6) * pixel_count
    @test provider_test_consumer(acquisition_products(replay)) ==
        T(7) * pixel_count

    reduced_implementation.state.command = T(1.5)
    @test @inferred(execute_acquisition!(reduced,
        provider_test_owner_rngs(plant, reduced))) ===
        acquisition_products(reduced)
    @test provider_test_consumer(acquisition_products(reduced)) ==
        T(2.5) * pixel_count

    fill!(path.result.values, T(99))
    for owner in (unchanged, copied)
        execute_acquisition!(owner, provider_test_owner_rngs(plant, owner))
    end
    @test provider_test_consumer(acquisition_products(unchanged)) ==
        T(5) * pixel_count
    @test provider_test_consumer(acquisition_products(copied)) ==
        T(6) * pixel_count

    execute_acquisition!(replay, provider_test_owner_rngs(plant, replay))
    @test provider_test_consumer(acquisition_products(replay)) ==
        T(8) * pixel_count
    execute_acquisition!(replay, provider_test_owner_rngs(plant, replay))
    @test provider_test_consumer(acquisition_products(replay)) ==
        T(7) * pixel_count

    full_selection = prepare_acquisition_selection(plant, :full)
    @test prepared_paths(full_selection) == (path,)
    @test execute_acquisition_selection!(full_selection,
        current_epoch(atmosphere)) === full_selection
    @test path.execution.executions[] == 1
    @test isfinite(provider_test_consumer(acquisition_products(full)))

    @test !ismutabletype(typeof(reduced))
    @test !ismutabletype(typeof(acquisition_provider(reduced)))

    mutable_contract_products = AcquisitionProducts(zeros(T, 2, 2);
        metadata=(labels=Int[1, 2],))
    defensive_provider = prepare_unchanged_synthetic_provider(
        mutable_contract_products)
    contract_snapshot = acquisition_product_contract(defensive_provider)
    contract_snapshot.metadata.labels[1] = 99
    @test acquisition_product_contract(defensive_provider).metadata.labels ==
        [1, 2]
    property_snapshot = defensive_provider.contract
    property_snapshot.metadata.labels[2] = 99
    @test acquisition_product_contract(defensive_provider).metadata.labels ==
        [1, 2]
    @test validate_acquisition_provider_binding(defensive_provider,
        nothing) === nothing
    @test !applicable(validate_acquisition_provider_binding,
        reduced_implementation, path.result,
        acquisition_products(reduced), full_contract)
    @test !applicable(AcquisitionProductContract,
        nothing, nothing, nothing)
    @test !applicable(PreparedAcquisitionProvider,
        PreparedUnchangedSyntheticProvider(), mutable_contract_products,
        contract_snapshot, SyntheticReplayProviderStyle(),
        :reuse_unchanged)

    valid_metadata = (kind=:contract, units=:count,
        geometry=(dimensions=(2, 2), layout=:dense),
        radiometry=(quantity=:detected_electrons, normalization=:per_sample))
    invalid_metadata = AcquisitionProducts(zeros(T, 2, 2);
        metadata=(kind=:contract, units=:count,
            geometry=(dimensions=(1, 2), layout=:dense),
            radiometry=(quantity=:detected_electrons,
                normalization=:per_sample)))
    invalid_metadata_type = AcquisitionProducts(zeros(T, 2, 2);
        metadata=(kind=:contract, units=:count,
            geometry=(dimensions=(2, 2), layout=:dense),
            radiometry=(quantity=:detected_electrons,
                normalization="per_sample")))
    valid_destination = AcquisitionProducts(zeros(T, 2, 2);
        metadata=valid_metadata)
    invalid_shape = AcquisitionProducts(zeros(T, 1, 2);
        metadata=valid_destination.metadata)
    invalid_type = AcquisitionProducts(zeros(Float32, 2, 2);
        metadata=valid_destination.metadata)
    invalid_device = AcquisitionProducts(ContractDeviceArray(
        zeros(T, 2, 2), ContractPlaneDevice(91));
        metadata=valid_destination.metadata)
    for (source_product, reason) in (
        (invalid_metadata, :metadata),
        (invalid_metadata_type, :metadata),
        (invalid_shape, :shape),
        (invalid_type, :numeric_type),
        (invalid_device, :device),
    )
        assert_plant_preparation_error(
            () -> prepare_copied_synthetic_provider(valid_destination,
                source_product),
            :acquisition, reason)
    end
    assert_plant_preparation_error(
        () -> prepare_cyclic_replay_provider(valid_destination, ()),
        :acquisition, :empty_replay)
    assert_plant_preparation_error(
        () -> prepare_cyclic_replay_provider(valid_destination, nothing),
        :acquisition, :invalid_replay)
    assert_plant_preparation_error(
        () -> prepare_cyclic_replay_provider(valid_destination, Any[1]),
        :acquisition, :product_type)

    observation_destination = WFSObservation(zeros(T, 2, 2);
        units=:electron_count, layout=:provider_frame)
    observation_source = WFSObservation(zeros(T, 2, 2);
        units=:adu, layout=:provider_frame)
    wfs_metadata = (kind=:wfs_observation, geometry=:test_grid)
    wfs_destination = AcquisitionProducts(observation_destination;
        metadata=wfs_metadata)
    wfs_source = AcquisitionProducts(observation_source;
        metadata=wfs_metadata)
    assert_plant_preparation_error(
        () -> prepare_copied_synthetic_provider(wfs_destination,
            wfs_source),
        :acquisition, :units)

    intensity_destination = IntensityMap(path.result.metadata,
        fill(zero(T), size(path.result.values)))
    intensity_source = IntensityMap(path.result.metadata,
        fill(T(21), size(path.result.values)))
    composite_observation_destination = WFSObservation(zeros(T, 2, 2);
        units=:electron_count, layout=:provider_frame)
    composite_observation_source = WFSObservation(fill(T(22), 2, 2);
        units=:electron_count, layout=:provider_frame)
    composite_measurement_destination = WFSMeasurement(Ref(zero(T));
        units=:metre, kind=:provider_scalar)
    composite_measurement_source = WFSMeasurement(Ref(T(23));
        units=:metre, kind=:provider_scalar)
    composite_metadata = (kind=:composite_provider_product,)
    composite_destination = AcquisitionProducts(
        (intensity_destination, composite_observation_destination),
        (composite_measurement_destination,);
        metadata=composite_metadata)
    composite_source = AcquisitionProducts(
        (intensity_source, composite_observation_source),
        (composite_measurement_source,);
        metadata=composite_metadata)
    composite_provider = prepare_copied_synthetic_provider(
        composite_destination, composite_source)
    @test acquisition_product_metadata(composite_destination) ===
        composite_metadata
    @test AdaptiveOpticsSim.execute_acquisition_provider!(
        composite_provider, nothing, Xoshiro(0x7303)) ===
        composite_destination
    @test all(==(T(21)), intensity_destination.values)
    @test all(==(T(22)), composite_observation_destination.storage)
    @test composite_measurement_destination.storage[] == T(23)

    tuple_metadata = (kind=:tuple_contract,)
    tuple_contract_products = AcquisitionProducts(
        (zeros(T, 1), zeros(T, 1)); metadata=tuple_metadata)
    tuple_contract = acquisition_product_contract(tuple_contract_products)
    for products in (
        AcquisitionProducts((zeros(T, 1),); metadata=tuple_metadata),
        AcquisitionProducts((zeros(T, 1), zeros(T, 1), zeros(T, 1));
            metadata=tuple_metadata),
    )
        assert_plant_preparation_error(
            () -> validate_acquisition_product_contract(products,
                tuple_contract),
            :acquisition, :shape)
    end
    tuple_type_mismatch = AcquisitionProducts(
        (zeros(T, 1), Ref(zero(T))); metadata=tuple_metadata)
    assert_plant_preparation_error(
        () -> validate_acquisition_product_contract(tuple_type_mismatch,
            tuple_contract),
        :acquisition, :product_type)
    assert_plant_preparation_error(
        () -> acquisition_product_contract(AcquisitionProducts(
            :unsupported; metadata=(kind=:unsupported,))),
        :acquisition, :unsupported_product)
    assert_plant_preparation_error(
        () -> AdaptiveOpticsSim.copy_acquisition_product!(1, 2),
        :acquisition, :unsupported_product_copy)

    owned_source = AcquisitionProducts(fill(T(12), 2, 2);
        metadata=valid_destination.metadata)
    owned_copy = prepare_copied_synthetic_provider(valid_destination,
        owned_source)
    fill!(owned_source.observation, T(13))
    AdaptiveOpticsSim.execute_acquisition_provider!(owned_copy, nothing,
        Xoshiro(0x7301))
    @test all(==(T(12)), valid_destination.observation)

    replay_source_a = AcquisitionProducts(fill(T(14), 2, 2);
        metadata=valid_destination.metadata)
    replay_source_b = AcquisitionProducts(fill(T(15), 2, 2);
        metadata=valid_destination.metadata)
    replay_sources = [replay_source_a, replay_source_b]
    owned_replay = prepare_cyclic_replay_provider(valid_destination,
        replay_sources)
    empty!(replay_sources)
    fill!(replay_source_a.observation, T(16))
    AdaptiveOpticsSim.execute_acquisition_provider!(owned_replay, nothing,
        Xoshiro(0x7302))
    @test all(==(T(14)), valid_destination.observation)
    assert_plant_preparation_error(
        () -> prepare_cyclic_replay_provider(valid_destination,
            Any[replay_source_a, wfs_destination]),
        :acquisition, :product_type)

    assert_plant_preparation_error(
        () -> AcquisitionProducts(zeros(T, 2, 2), nothing, nothing),
        :acquisition, :metadata)
    assert_plant_preparation_error(
        () -> PreparedAcquisitionProvider(UnsupportedProviderStyle(),
            valid_destination),
        :acquisition, :unsupported_provider)
    assert_plant_preparation_error(
        () -> PreparedAcquisitionProvider(InvalidProviderStyle(),
            valid_destination),
        :acquisition, :invalid_provider_style)
    assert_plant_preparation_error(
        () -> PreparedAcquisitionProvider(MissingPayloadWorkProvider(),
            valid_destination),
        :acquisition, :missing_payload_work)
    assert_plant_preparation_error(
        () -> PreparedAcquisitionProvider(EmptyPayloadWorkProvider(),
            valid_destination),
        :acquisition, :invalid_payload_work)
    assert_plant_preparation_error(
        () -> PreparedAcquisitionProvider(InvalidPayloadWorkProvider(),
            valid_destination),
        :acquisition, :invalid_payload_work)

    unsupported_validation = PreparedAcquisitionProvider(
        UnsupportedProviderValidation(), valid_destination)
    assert_plant_preparation_error(
        () -> validate_acquisition_provider_binding(
            unsupported_validation, nothing),
        :acquisition, :unsupported_provider_validation)
    unsupported_execution = PreparedAcquisitionProvider(
        UnsupportedProviderExecution(), valid_destination)
    assert_plant_preparation_error(
        () -> AdaptiveOpticsSim.execute_acquisition_provider!(
            unsupported_execution, nothing, Xoshiro(0x7304)),
        :acquisition, :unsupported_provider_execution)

    invalid_provider = PreparedAcquisitionProvider(InvalidProviderResult(),
        provider_test_products(path, zero(T)))
    invalid_definition = AcquisitionDefinition(:invalid_result,
        :provider_path, ProviderReducedModel(zero(T)))
    invalid_owner = PreparedAcquisitionOwner(invalid_definition, path,
        invalid_provider)
    assert_plant_preparation_error(
        () -> execute_acquisition!(invalid_owner, Xoshiro(0x7310)),
        :acquisition, :provider_result)

    if coverage_instrumented()
        @test_skip "provider allocation assertions are disabled under coverage instrumentation"
    else
        for owner in (full, reduced, unchanged, copied, replay)
            rngs = provider_test_owner_rngs(plant, owner)
            @test provider_execution_allocations(owner, rngs) == 0
        end
    end
end

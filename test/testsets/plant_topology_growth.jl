using InteractiveUtils: code_native

struct TopologyGrowthPathModel end
struct TopologyGrowthAcquisitionModel end
struct TopologyGrowthPathFamilyModel{N} end
struct TopologyGrowthAcquisitionFamilyModel{N} end
struct TopologyGrowthExecutionFamily{N} end

struct TopologyGrowthMaterialization{P}
    destination::P
end

function topology_growth_native_generic_apply_count(function_value,
    argument_types)
    native = sprint() do io
        code_native(io, function_value, argument_types;
            debuginfo=:none, binary=false)
    end
    return count("jl_apply_generic", native)
end

function topology_growth_selection_allocations(selection, epoch)
    execute_acquisition_selection!(selection, epoch)
    return @allocated execute_acquisition_selection!(selection, epoch)
end

function topology_growth_preparation_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

struct TopologyGrowthPathExecution{F,P,R}
    family::F
    input::P
    result::R
end

struct TopologyGrowthAcquisitionExecution{F,R,O}
    family::F
    source::R
    destination::O
end

Plant.plant_model_definition_style(::Type{TopologyGrowthPathModel}) =
    ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{TopologyGrowthAcquisitionModel}) =
    ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{<:TopologyGrowthPathFamilyModel}) =
    ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{<:TopologyGrowthAcquisitionFamilyModel}) =
    ColdPlantModelDefinition()

function Plant.validate_path_materialization_binding(
    materialization::TopologyGrowthMaterialization,
    input::PupilFunction,
    ::AbstractAtmosphere,
    ::AbstractSource,
)
    materialization.destination === input || throw(PlantPreparationError(
        :path, :prepared_binding,
        "topology-growth materialization binding changed"))
    return nothing
end

function Plant.validate_path_materialization_target(
    materialization::TopologyGrowthMaterialization,
    input::PupilFunction,
    ::AbstractAtmosphere,
    target::AdaptiveOpticsSim.Backends.AbstractComputeDevice,
)
    materialization.destination === input || throw(PlantPreparationError(
        :path, :prepared_binding,
        "topology-growth materialization binding changed"))
    Plant._require_exact_plant_product_target(
        materialization.destination, target,
        "topology-growth materialization destination")
    return materialization
end

function Plant.validate_path_materialization(
    materialization::TopologyGrowthMaterialization,
    input::PupilFunction,
    ::AbstractAtmosphere,
    ::AtmosphereEpoch,
)
    materialization.destination === input || throw(PlantPreparationError(
        :path, :prepared_binding,
        "topology-growth materialization binding changed"))
    return nothing
end

function Plant.materialize_path_input!(
    materialization::TopologyGrowthMaterialization,
    input::PupilFunction,
    ::AbstractAtmosphere,
    ::AtmosphereEpoch,
)
    materialization.destination === input || throw(PlantPreparationError(
        :path, :prepared_binding,
        "topology-growth materialization binding changed"))
    fill!(input.opd, zero(eltype(input.opd)))
    return input
end

function Plant.validate_path_execution_binding(
    execution::TopologyGrowthPathExecution,
    input::PupilFunction,
    result::IntensityMap,
)
    execution.input === input && execution.result === result || throw(
        PlantPreparationError(:path, :prepared_binding,
            "topology-growth path execution binding changed"))
    return nothing
end

function Plant.validate_path_execution_target(
    execution::TopologyGrowthPathExecution,
    target::AdaptiveOpticsSim.Backends.AbstractComputeDevice,
)
    Plant._require_exact_plant_product_target(
        execution.input, target, "topology-growth path input")
    Plant._require_exact_plant_product_target(
        execution.result, target, "topology-growth path result")
    return execution
end

function Plant.execute_path!(
    result::IntensityMap,
    input::PupilFunction,
    execution::TopologyGrowthPathExecution{
        TopologyGrowthExecutionFamily{N}},
) where {N}
    Plant.validate_path_execution_binding(execution, input, result)
    family_value = convert(eltype(result.values), N)
    @. result.values = input.opd + family_value
    return result
end

function topology_growth_prepare_path_executor(
    family,
    definition::OpticalPathDefinition,
    source::AbstractSource,
    telescope::Telescope,
    atmosphere::AbstractTimedAtmosphere,
    context,
)
    T = eltype(pupil_reflectivity(telescope))
    pupil = PupilFunction(telescope; T, backend=backend(telescope))
    values = similar(pupil.opd)
    fill!(values, zero(T))
    metadata = OpticalPlaneMetadata(
        FocalPlane(),
        values;
        coordinate_domain=MetricCoordinates(),
        sampling=pupil.metadata.sampling,
        origin=pupil.metadata.origin,
        centering=pupil.metadata.centering,
        orientation=pupil.metadata.orientation,
        spectral=MonochromaticChannel(T(wavelength(source))),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition(),
    )
    result = IntensityMap(metadata, values)
    execution = TopologyGrowthPathExecution(family, pupil, result)
    return PreparedPathExecutor(
        definition,
        source,
        telescope,
        atmosphere,
        pupil,
        result,
        execution;
        context=context,
        materialization=TopologyGrowthMaterialization(pupil),
        optical_model=:topology_growth_handoff,
        propagation_model=:test_typed_handoff,
        model_revisions=UInt(1),
    )
end

function Plant.prepare_path_executor(
    ::TopologyGrowthPathModel,
    definition::OpticalPathDefinition,
    source::AbstractSource,
    telescope::Telescope,
    atmosphere::AbstractTimedAtmosphere,
    context,
)
    return topology_growth_prepare_path_executor(
        TopologyGrowthExecutionFamily{1}(), definition, source, telescope,
        atmosphere, context)
end

function Plant.prepare_path_executor(
    ::TopologyGrowthPathFamilyModel{N},
    definition::OpticalPathDefinition,
    source::AbstractSource,
    telescope::Telescope,
    atmosphere::AbstractTimedAtmosphere,
    context,
) where {N}
    return topology_growth_prepare_path_executor(
        TopologyGrowthExecutionFamily{N}(), definition, source, telescope,
        atmosphere, context)
end

function Plant.validate_acquisition_execution_binding(
    execution::TopologyGrowthAcquisitionExecution,
    path_result::IntensityMap,
    products::AcquisitionProducts,
)
    execution.source === path_result &&
        execution.destination === products.observation || throw(
        PlantPreparationError(:acquisition, :prepared_binding,
            "topology-growth acquisition binding changed"))
    return nothing
end

function Plant.validate_acquisition_execution_target(
    execution::TopologyGrowthAcquisitionExecution,
    target::AdaptiveOpticsSim.Backends.AbstractComputeDevice,
)
    Plant._require_exact_plant_product_target(
        execution.source, target, "topology-growth acquisition source")
    Plant._require_exact_plant_product_target(
        execution.destination, target,
        "topology-growth acquisition destination")
    return execution
end

function Plant.execute_acquisition!(
    products::AcquisitionProducts{<:AbstractMatrix,Nothing},
    path_result::IntensityMap,
    execution::TopologyGrowthAcquisitionExecution{
        TopologyGrowthExecutionFamily{N}},
    ::AbstractRNG,
) where {N}
    Plant.validate_acquisition_execution_binding(
        execution, path_result, products)
    family_value = convert(eltype(products.observation), 100 * N)
    @. products.observation = path_result.values + family_value
    return products
end

function topology_growth_prepare_acquisition_provider(
    family,
    ::AcquisitionDefinition,
    path::PreparedPathExecutor,
)
    require_path_result(path)
    observation = similar(path.result.values)
    fill!(observation, zero(eltype(observation)))
    execution =
        TopologyGrowthAcquisitionExecution(
            family, path.result, observation)
    products = AcquisitionProducts(
        observation;
        metadata=(
            kind=:topology_growth_frame,
            units=:test_samples,
            geometry=path.result.metadata,
        ),
    )
    return prepare_full_optical_provider(execution, products)
end
function Plant.prepare_acquisition_provider(
    ::TopologyGrowthAcquisitionModel,
    definition::AcquisitionDefinition,
    path::PreparedPathExecutor,
)
    return topology_growth_prepare_acquisition_provider(
        TopologyGrowthExecutionFamily{1}(), definition, path)
end

function Plant.prepare_acquisition_provider(
    ::TopologyGrowthAcquisitionFamilyModel{N},
    definition::AcquisitionDefinition,
    path::PreparedPathExecutor,
) where {N}
    return topology_growth_prepare_acquisition_provider(
        TopologyGrowthExecutionFamily{N}(), definition, path)
end

function topology_growth_fixture(path_count::Integer; family_count::Int=1,
    path_families=nothing, acquisition_families=nothing)
    path_count > 0 || error("topology-growth path count must be positive")
    1 <= family_count <= path_count ||
        error("topology-growth family count must be within path count")
    path_family_indices = path_families === nothing ?
        [mod1(index, family_count) for index in 1:path_count] :
        collect(Int, path_families)
    acquisition_family_indices = acquisition_families === nothing ?
        copy(path_family_indices) : collect(Int, acquisition_families)
    length(path_family_indices) == path_count ||
        error("topology-growth path-family count must match path count")
    length(acquisition_family_indices) == path_count ||
        error("topology-growth acquisition-family count must match path count")
    all(index -> 1 <= index <= family_count, path_family_indices) ||
        error("topology-growth path-family index is outside family count")
    all(index -> 1 <= index <= family_count, acquisition_family_indices) ||
        error("topology-growth acquisition-family index is outside family count")
    T = Float64
    telescope = Telescope(
        resolution=5,
        diameter=T(5),
        central_obstruction=zero(T),
        T=T,
    )
    atmosphere = MultiLayerAtmosphere(
        telescope;
        r0=T(0.2),
        L0=T(25),
        fractional_cn2=T[1],
        wind_speed=T[0],
        wind_direction_deg=T[0],
        altitude=T[0],
        layer_ids=(:ground,),
        T=T,
    )
    source = Source(
        band=:custom,
        wavelength=T(0.8e-6),
        photon_irradiance=T(100),
        T=T,
    )
    paths = OpticalPathDefinition[
        OpticalPathDefinition(
            Symbol(:path_, index), source,
            family_count == 1 ? TopologyGrowthPathModel() :
            TopologyGrowthPathFamilyModel{
                path_family_indices[index]}())
        for index in 1:path_count
    ]
    acquisitions = AcquisitionDefinition[
        AcquisitionDefinition(
            Symbol(:acquisition_, index),
            path_id(paths[index]),
            family_count == 1 ? TopologyGrowthAcquisitionModel() :
            TopologyGrowthAcquisitionFamilyModel{
                acquisition_family_indices[index]}(),
        )
        for index in 1:path_count
    ]
    definition = PlantDefinition(;
        telescope=plant_test_telescope_definition(telescope),
        atmosphere=plant_test_atmosphere_definition(atmosphere),
        paths,
        acquisitions,
    )
    plant = prepare_plant(
        definition, PLANT_TEST_HOST_TARGET;
        run_seed=UInt64(0x9800 + path_count))
    selection = prepare_acquisition_selection(
        plant,
        Symbol[acquisition_id(acquisition).name
            for acquisition in acquisitions],
    )
    return (;
        definition,
        plant,
        selection,
        atmosphere=prepared_atmosphere(plant),
        path_family_indices,
        acquisition_family_indices,
    )
end

@inline function topology_growth_is_control_statement(statement)
    return statement isa Core.ReturnNode ||
        statement isa Core.GotoNode ||
        statement isa Core.GotoIfNot ||
        (statement isa Expr && statement.head in (
            :gc_preserve_begin,
            :gc_preserve_end,
            :boundscheck,
            :inbounds,
        ))
end

function topology_growth_any_value_count(function_value, argument_types)
    typed = only(Base.code_typed(
        function_value, argument_types; optimize=true))
    code = typed.first
    slot_count = count(==(Any), code.slottypes)
    ssa_count = 0
    for index in eachindex(code.code)
        code.ssavaluetypes[index] === Any || continue
        topology_growth_is_control_statement(code.code[index]) && continue
        ssa_count += 1
    end
    return slot_count + ssa_count
end

@testset "Count-invariant prepared plant topology" begin
    path_counts = (1, 4, 8, 16)
    fixtures = map(topology_growth_fixture, path_counts)

    for field in (:definition, :plant, :selection)
        types = map(fixture -> typeof(getproperty(fixture, field)), fixtures)
        @test length(unique(types)) == 1
    end
    @test length(unique(map(
        fixture -> typeof(fixture.plant.rngs), fixtures))) == 1

    for (fixture, path_count) in zip(fixtures, path_counts)
        definition = fixture.definition
        plant = fixture.plant
        selection = fixture.selection
        @test path_definitions(definition) isa
            AbstractVector{OpticalPathDefinition}
        @test acquisition_definitions(definition) isa
            AbstractVector{AcquisitionDefinition}
        @test prepared_paths(plant) isa Plant._PreparedPathRegistry
        @test prepared_acquisitions(plant) isa
            Plant._PreparedAcquisitionRegistry
        @test prepared_paths(selection) isa Plant._PreparedPathSelection
        @test prepared_acquisitions(selection) isa
            Plant._PreparedAcquisitionOwnerSelection
        @test all(group -> group.values isa
            FixedSizeVector{<:PreparedPathExecutor}, plant.paths.groups)
        @test all(group -> group.values isa
            FixedSizeVector{<:PreparedAcquisitionOwner},
            plant.acquisitions.groups)
        @test length(plant.paths.groups) == 1
        @test length(plant.acquisitions.groups) == 1
        @test plant.rngs.paths isa Plant._PreparedOwnerRNGRegistry
        @test plant.rngs.acquisitions isa Plant._PreparedOwnerRNGRegistry
        @test plant.rngs.paths.slots === plant.paths.slots
        @test plant.rngs.acquisitions.slots === plant.acquisitions.slots
        @test all(group -> group.values isa
            FixedSizeVector{<:PreparedOwnerRNGs},
            plant.rngs.paths.groups)
        @test all(group -> group.values isa
            FixedSizeVector{<:PreparedOwnerRNGs},
            plant.rngs.acquisitions.groups)
        @test !hasfield(typeof(selection), :path_rng_slots)
        @test !hasfield(typeof(selection), :acquisition_rng_slots)
        @test selection.sampled_aberration_path_plans isa
            Plant._PreparedSampledAberrationPlanRegistry
        @test all(group -> group.values isa FixedSizeVector,
            selection.sampled_aberration_path_plans.groups)
        @test length(prepared_paths(plant)) == path_count
        @test length(prepared_acquisitions(plant)) == path_count
        @test length(prepared_paths(selection)) == path_count
        @test length(prepared_acquisitions(selection)) == path_count
        @test !hasfield(eltype(plant.paths.groups[1].values), :definition)
        @test !hasfield(eltype(plant.acquisitions.groups[1].values),
            :definition)
        for path in prepared_paths(selection)
            @test prepared_path(plant, path_id(path)) === path
            @test Plant._prepared_path_definition(plant.paths, path) ===
                path_definitions(definition)[Int(path.definition_slot)]
        end
        for owner in prepared_acquisitions(selection)
            @test prepared_acquisition(plant, acquisition_id(owner)) ===
                owner
            @test Plant._prepared_acquisition_definition(
                plant.acquisitions, owner) ===
                acquisition_definitions(definition)[
                    Int(owner.definition_slot)]
        end
    end

    largest = last(fixtures)
    @test @inferred(execute_acquisition_selection_at!(
        largest.selection, 0.0)) === largest.selection
    for path in prepared_paths(largest.selection)
        @test all(==(1.0), path_result(path).values)
    end
    for owner in prepared_acquisitions(largest.selection)
        @test all(==(101.0), acquisition_observation(owner))
    end

    representative_path = first(prepared_paths(largest.plant))
    representative_acquisition =
        first(prepared_acquisitions(largest.plant))
    representative_execution =
        representative_acquisition.provider.implementation.execution
    representative_products = representative_acquisition.provider.products
    representative_rng = Xoshiro(0x9800)

    @test topology_growth_any_value_count(
        Plant.execute_path!,
        Tuple{
            typeof(representative_path.result),
            typeof(representative_path.input),
            typeof(representative_path.execution),
        },
    ) == 0
    @test topology_growth_any_value_count(
        Plant.execute_acquisition!,
        Tuple{
            typeof(representative_products),
            typeof(representative_acquisition.path_result),
            typeof(representative_execution),
            typeof(representative_rng),
        },
    ) == 0
    @test @inferred(Plant.execute_path!(
        representative_path.result,
        representative_path.input,
        representative_path.execution,
    )) === representative_path.result
    @test @inferred(Plant.execute_acquisition!(
        representative_products,
        representative_acquisition.path_result,
        representative_execution,
        representative_rng,
    )) === representative_products
end

function assert_topology_growth_descriptor_routing(fixture)
    execute_acquisition_selection_at!(fixture.selection, 0.0)
    plant = fixture.plant
    for index in eachindex(fixture.path_family_indices)
        path = prepared_path(plant, Symbol(:path_, index))
        path_slot = plant.paths.slots[index]
        @test Plant._prepared_path_value(plant.paths, path_slot) === path
        @test path_id(path) == OpticalPathID(Symbol(:path_, index))
        @test all(==(Float64(fixture.path_family_indices[index])),
            path_result(path).values)

        owner = prepared_acquisition(
            plant, Symbol(:acquisition_, index))
        acquisition_slot = plant.acquisitions.slots[index]
        @test Plant._prepared_acquisition_value(
            plant.acquisitions, acquisition_slot) === owner
        @test acquisition_id(owner) ==
            AcquisitionID(Symbol(:acquisition_, index))
        @test acquisition_path_id(owner) == path_id(path)
        expected = Float64(fixture.path_family_indices[index] +
            100 * fixture.acquisition_family_indices[index])
        @test all(==(expected), acquisition_observation(owner))
    end
    return fixture
end

topology_growth_group_family_indices(groups, family) =
    map(groups) do group
        family_type = typeof(family(first(group.values)))
        return only(family_type.parameters)
    end

@testset "Mixed exact execution-family code generation" begin
    homogeneous = topology_growth_fixture(4)
    mixed = topology_growth_fixture(6; family_count=6,
        path_families=(4, 1, 6, 2, 5, 3),
        acquisition_families=(2, 6, 3, 1, 5, 4))
    reordered = topology_growth_fixture(6; family_count=6,
        path_families=(3, 5, 2, 6, 1, 4),
        acquisition_families=(4, 1, 5, 3, 6, 2))
    extended = topology_growth_fixture(7; family_count=6,
        path_families=(4, 1, 6, 2, 5, 3, 4),
        acquisition_families=(2, 6, 3, 1, 5, 4, 2))

    @test length(mixed.plant.paths.groups) == 6
    @test length(mixed.plant.acquisitions.groups) == 6
    @test typeof(mixed.plant.paths) != typeof(homogeneous.plant.paths)
    @test typeof(mixed.plant.acquisitions) !=
        typeof(homogeneous.plant.acquisitions)
    path_group_indices = map(fixture ->
        topology_growth_group_family_indices(
            fixture.plant.paths.groups, path -> path.execution.family),
        (mixed, reordered, extended))
    acquisition_group_indices = map(fixture ->
        topology_growth_group_family_indices(
            fixture.plant.acquisitions.groups,
            owner -> owner.provider.implementation.execution.family),
        (mixed, reordered, extended))
    @test path_group_indices[1] == path_group_indices[2] ==
        path_group_indices[3]
    @test acquisition_group_indices[1] == acquisition_group_indices[2] ==
        acquisition_group_indices[3]
    path_group_types_match = typeof(mixed.plant.paths.groups) ==
        typeof(reordered.plant.paths.groups) ==
        typeof(extended.plant.paths.groups)
    acquisition_group_types_match =
        typeof(mixed.plant.acquisitions.groups) ==
        typeof(reordered.plant.acquisitions.groups) ==
        typeof(extended.plant.acquisitions.groups)
    path_registry_types_match = typeof(mixed.plant.paths) ==
        typeof(reordered.plant.paths) == typeof(extended.plant.paths)
    acquisition_registry_types_match =
        typeof(mixed.plant.acquisitions) ==
        typeof(reordered.plant.acquisitions) ==
        typeof(extended.plant.acquisitions)
    plant_types_match = typeof(mixed.plant) == typeof(reordered.plant) ==
        typeof(extended.plant)
    selection_types_match = typeof(mixed.selection) ==
        typeof(reordered.selection) == typeof(extended.selection)
    @test path_group_types_match
    @test acquisition_group_types_match
    @test path_registry_types_match
    @test acquisition_registry_types_match
    @test plant_types_match
    @test selection_types_match

    assert_topology_growth_descriptor_routing(mixed)
    assert_topology_growth_descriptor_routing(reordered)
    assert_topology_growth_descriptor_routing(extended)

    selection = mixed.selection
    atmosphere = prepared_atmosphere(mixed.plant)
    epoch = current_epoch(atmosphere)
    Plant._validate_selection_epoch!(selection, atmosphere, epoch)
    path_groups = mixed.plant.paths.groups
    path_rng_groups = mixed.plant.rngs.paths.groups
    acquisition_groups = mixed.plant.acquisitions.groups
    acquisition_rng_groups = mixed.plant.rngs.acquisitions.groups

    @test topology_growth_native_generic_apply_count(
        Plant._selected_path_materialize_family!, Tuple{
            typeof(path_groups),
            typeof(path_rng_groups),
            Int,
            Int,
            typeof(atmosphere),
            typeof(epoch),
        }) == 0
    @test topology_growth_native_generic_apply_count(
        Plant._selected_path_execute_family!, Tuple{
            typeof(path_groups),
            typeof(path_rng_groups),
            Int,
            Int,
        }) == 0
    @test topology_growth_native_generic_apply_count(
        Plant._selected_acquisition_execute_family!, Tuple{
            typeof(acquisition_groups),
            typeof(acquisition_rng_groups),
            Int,
            Int,
        }) == 0
    @test topology_growth_native_generic_apply_count(
        Plant._materialize_selected_paths!, Tuple{
            typeof(mixed.plant),
            typeof(selection.paths),
            typeof(atmosphere),
            typeof(epoch),
        }) == 0
    @test topology_growth_native_generic_apply_count(
        Plant._apply_selected_sampled_aberrations!, Tuple{
            typeof(selection.paths),
            typeof(selection.sampled_aberration_path_plans),
        }) == 0
    @test topology_growth_native_generic_apply_count(
        Plant._execute_selected_paths!, Tuple{
            typeof(mixed.plant),
            typeof(selection.paths),
        }) == 0
    @test topology_growth_native_generic_apply_count(
        Plant._execute_selected_acquisitions!, Tuple{
            typeof(mixed.plant),
            typeof(selection.acquisitions),
        }) == 0

    if coverage_instrumented()
        @test_skip "mixed-family allocation gate disabled under coverage instrumentation"
    else
        @test topology_growth_selection_allocations(
            selection, current_epoch(atmosphere)) == 0
    end

    first_path = first(prepared_paths(selection))
    second_path = prepared_paths(selection)[2]
    first_rngs = mixed.plant.rngs.paths[Int(first_path.definition_slot)]
    rng_before = copy(rng_stream_state(first_rngs, Val(:provider)))
    owner_error = topology_growth_preparation_error(
        () -> Plant._require_rng_owner_binding(first_rngs, second_path))
    @test owner_error isa PlantPreparationError
    @test owner_error.component === :rng
    @test owner_error.reason === :prepared_binding
    @test rng_stream_state(first_rngs, Val(:provider)) == rng_before

    copied_slots = copy(mixed.plant.paths.slots)
    foreign_registry = Plant._PreparedOwnerRNGRegistry(
        mixed.plant.rngs.paths.groups, copied_slots)
    registry_error = topology_growth_preparation_error(
        () -> Plant._require_rng_registry_binding(
            foreign_registry, mixed.plant.paths, :path))
    @test registry_error isa PlantPreparationError
    @test registry_error.component === :rng
    @test registry_error.reason === :owner_topology
    path_rng_groups = mixed.plant.rngs.paths.groups
    first_rng_family = first(path_rng_groups)
    R = eltype(first_rng_family.values)
    empty_rng_values = FixedSizeVectorDefault{R}(R[])
    empty_rng_family = Plant._PreparedOwnerRNGFamily{
        R,typeof(empty_rng_values)}(empty_rng_values)
    short_registry = Plant._PreparedOwnerRNGRegistry(
        (empty_rng_family, Base.tail(path_rng_groups)...),
        mixed.plant.paths.slots)
    cardinality_error = topology_growth_preparation_error(
        () -> Plant._require_rng_registry_binding(
            short_registry, mixed.plant.paths, :path))
    @test cardinality_error isa PlantPreparationError
    @test cardinality_error.component === :rng
    @test cardinality_error.reason === :owner_topology
    bounds_error = topology_growth_preparation_error(
        () -> Plant._prepared_owner_rng_family_value(
            mixed.plant.rngs.paths.groups,
            1, length(first_rng_family.values) + 1))
    @test bounds_error isa PlantPreparationError
    @test bounds_error.component === :rng
    @test bounds_error.reason === :owner_topology

end

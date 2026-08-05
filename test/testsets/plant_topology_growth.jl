struct TopologyGrowthPathModel end
struct TopologyGrowthAcquisitionModel end

struct TopologyGrowthMaterialization{P}
    destination::P
end

struct TopologyGrowthPathExecution{P,R}
    input::P
    result::R
end

struct TopologyGrowthAcquisitionExecution{R,O}
    source::R
    destination::O
end

Plant.plant_model_definition_style(::Type{TopologyGrowthPathModel}) =
    ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{TopologyGrowthAcquisitionModel}) =
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
    execution::TopologyGrowthPathExecution,
)
    Plant.validate_path_execution_binding(execution, input, result)
    @. result.values = input.opd + one(eltype(result.values))
    return result
end

function Plant.prepare_path_executor(
    ::TopologyGrowthPathModel,
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
    execution = TopologyGrowthPathExecution(pupil, result)
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
    execution::TopologyGrowthAcquisitionExecution,
    ::AbstractRNG,
)
    Plant.validate_acquisition_execution_binding(
        execution, path_result, products)
    copyto!(products.observation, path_result.values)
    return products
end

function Plant.prepare_acquisition_provider(
    ::TopologyGrowthAcquisitionModel,
    ::AcquisitionDefinition,
    path::PreparedPathExecutor,
)
    require_path_result(path)
    observation = similar(path.result.values)
    fill!(observation, zero(eltype(observation)))
    execution =
        TopologyGrowthAcquisitionExecution(path.result, observation)
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

function topology_growth_fixture(path_count::Integer)
    path_count > 0 || error("topology-growth path count must be positive")
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
            Symbol(:path_, index), source, TopologyGrowthPathModel())
        for index in 1:path_count
    ]
    acquisitions = AcquisitionDefinition[
        AcquisitionDefinition(
            Symbol(:acquisition_, index),
            path_id(paths[index]),
            TopologyGrowthAcquisitionModel(),
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
    path_counts = (4, 8, 16)
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
        @test prepared_paths(plant) isa Memory{PreparedPathExecutor}
        @test prepared_acquisitions(plant) isa
            Memory{PreparedAcquisitionOwner}
        @test prepared_paths(selection) isa Memory{PreparedPathExecutor}
        @test prepared_acquisitions(selection) isa
            Memory{PreparedAcquisitionOwner}
        @test plant.rngs.paths isa Memory{PreparedOwnerRNGs}
        @test plant.rngs.acquisitions isa Memory{PreparedOwnerRNGs}
        @test length(prepared_paths(plant)) == path_count
        @test length(prepared_acquisitions(plant)) == path_count
        @test length(prepared_paths(selection)) == path_count
        @test length(prepared_acquisitions(selection)) == path_count
    end

    largest = last(fixtures)
    @test @inferred(execute_acquisition_selection_at!(
        largest.selection, 0.0)) === largest.selection
    for path in prepared_paths(largest.selection)
        @test all(==(1.0), path_result(path).values)
    end
    for owner in prepared_acquisitions(largest.selection)
        @test all(==(1.0), acquisition_observation(owner))
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

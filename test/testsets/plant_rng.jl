struct RNGTestPathModel{R}
    zero_padding::Int
    revision::R
end

struct RNGTestPathExecution{E}
    imaging::E
end

struct RNGTestAcquisitionModel{T<:AbstractFloat}
    exposure::T
    readout_sigma::T
end

struct RNGTestAcquisitionExecution{E}
    frame::E
end

Plant.plant_model_definition_style(
    ::Type{<:RNGTestPathModel}) = ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{<:RNGTestAcquisitionModel}) = ColdPlantModelDefinition()

Plant.additional_path_rng_owner_roles(
    ::RNGTestPathExecution) = (:path_device,)
Plant.additional_acquisition_rng_owner_roles(
    ::RNGTestAcquisitionExecution) = (:device_model,)

function Plant.validate_path_execution_binding(
    execution::RNGTestPathExecution, input, result)
    return Plant.validate_path_execution_binding(
        execution.imaging, input, result)
end

function Plant.execute_path_rngs!(result, input,
    execution::RNGTestPathExecution,
    rngs::Plant.PreparedOwnerRNGs)
    Plant.validate_path_execution_binding(execution, input,
        result)
    Plant.execute_path!(result, input, execution.imaging)
    provider_rng = Plant.rng_stream_state(rngs,
        Val(:provider))
    device_rng = Plant.rng_stream_state(rngs,
        Val(:path_device))
    T = eltype(result.values)
    scale = one(T) + T(1e-3) * randn(provider_rng, T) +
        T(1e-4) * rand(device_rng, T)
    result.values .*= scale
    return result
end

function Plant.prepare_path_executor(
    model::RNGTestPathModel,
    definition::OpticalPathDefinition,
    source::AbstractSource,
    telescope::Telescope,
    atmosphere::AdaptiveOpticsSim.AbstractTimedAtmosphere,
)
    T = eltype(pupil_reflectivity(telescope))
    pupil = PupilFunction(telescope; T=T, backend=backend(telescope))
    imaging = prepare_direct_imaging(pupil, source;
        zero_padding=model.zero_padding)
    execution = RNGTestPathExecution(imaging)
    return PreparedPathExecutor(
        definition,
        source,
        telescope,
        atmosphere,
        pupil,
        direct_imaging_output(imaging),
        execution;
        materialization=prepare_pupil_opd_materialization(atmosphere,
            telescope, source, pupil),
        optical_model=(kind=:rng_test_direct_imaging,
            zero_padding=model.zero_padding),
        propagation_model=:fraunhofer_fft,
        model_revisions=model.revision,
    )
end

function Plant.validate_acquisition_execution_binding(
    execution::RNGTestAcquisitionExecution, path_result,
    products::AcquisitionProducts)
    return Plant.validate_acquisition_execution_binding(
        execution.frame, path_result, products)
end

function Plant.execute_acquisition_rngs!(products, path_result,
    execution::RNGTestAcquisitionExecution,
    rngs::Plant.PreparedOwnerRNGs)
    Plant.validate_acquisition_execution_binding(execution,
        path_result, products)
    detector_rng = Plant.rng_stream_state(rngs,
        Val(:detector))
    device_rng = Plant.rng_stream_state(rngs,
        Val(:device_model))
    Plant.execute_acquisition!(products, path_result,
        execution.frame, detector_rng)
    T = eltype(products.observation)
    products.observation .+= T(1e-3) * randn(device_rng, T)
    return products
end

function Plant.prepare_acquisition_provider(
    model::RNGTestAcquisitionModel,
    definition::AcquisitionDefinition,
    path::PreparedPathExecutor,
)
    require_path_result(path)
    T = eltype(path.result.values)
    detector = Detector(integration_time=T(model.exposure),
        noise=NoisePhotonReadout(T(model.readout_sigma)), qe=one(T),
        response_model=NullFrameResponse(), T=T,
        backend=path.key.backend)
    frame = FrameAcquisitionExecution(detector, path.result)
    execution = RNGTestAcquisitionExecution(frame)
    metadata = (kind=:detector_frame, units=:detected_electrons,
        geometry=path.result.metadata,
        detector=detector_export_metadata(detector))
    products = AcquisitionProducts(frame.observation; metadata)
    return prepare_full_optical_provider(execution, products)
end

function rng_test_definition(;
    path_order::Tuple=(:random, :quiet),
    acquisition_order::Tuple=(:noisy, :quiet_frame),
    layer_order::Tuple=(1, 2),
    layer_ids=(:ground, :high),
    noisy_acquisition_id::Symbol=:noisy,
    atmosphere_model::Symbol=:multilayer,
)
    T = Float64
    telescope = Telescope(resolution=8, diameter=T(4),
        central_obstruction=zero(T), T=T)
    base_cn2 = (T(0.6), T(0.4))
    base_speed = (T(7), T(19))
    base_direction = (T(15), T(105))
    base_altitude = (T(0), T(8_000))
    resolved_layer_ids = isnothing(layer_ids) ? nothing :
        ntuple(index -> layer_ids[layer_order[index]], length(layer_order))
    atmosphere = if atmosphere_model === :multilayer
        MultiLayerAtmosphere(telescope;
            r0=T(0.2),
            L0=T(25),
            fractional_cn2=T[base_cn2[index] for index in layer_order],
            wind_speed=T[base_speed[index] for index in layer_order],
            wind_direction=T[base_direction[index] for index in layer_order],
            altitude=T[base_altitude[index] for index in layer_order],
            layer_ids=resolved_layer_ids,
            T=T,
        )
    elseif atmosphere_model === :kolmogorov
        KolmogorovAtmosphere(telescope; r0=T(0.2), L0=T(25), T=T)
    else
        error("unsupported test atmosphere model $atmosphere_model")
    end
    random_source = Source(band=:custom, wavelength=T(0.8e-6),
        photon_irradiance=T(80), coordinates=(T(1), T(20)), T=T)
    quiet_source = Source(band=:custom, wavelength=T(0.7e-6),
        photon_irradiance=T(60), coordinates=(T(-2), T(15)), T=T)
    path_map = (
        random=OpticalPathDefinition(:random, random_source,
            RNGTestPathModel(2, UInt(31))),
        quiet=OpticalPathDefinition(:quiet, quiet_source,
            DirectSciencePathModel(2, UInt(32))),
    )
    acquisition_map = (
        noisy=AcquisitionDefinition(noisy_acquisition_id, :random,
            RNGTestAcquisitionModel(T(0.5), T(0.25))),
        quiet_frame=AcquisitionDefinition(:quiet_frame, :quiet,
            FramePlantAcquisitionModel(T(0.5), MatchingPlantPath())),
    )
    paths = ntuple(index -> getproperty(path_map, path_order[index]),
        length(path_order))
    acquisitions = ntuple(index ->
        getproperty(acquisition_map, acquisition_order[index]),
        length(acquisition_order))
    return PlantDefinition(; telescope, atmosphere, paths, acquisitions)
end

function rng_test_owner_seed(metadata, category::Symbol,
    component::Symbol, role::Symbol)
    matches = filter(metadata.owners) do owner
        identity = owner.identity
        identity.category === category &&
            identity.component === component && identity.role === role
    end
    return only(matches).derived_seed
end

function rng_test_observations(plant::PreparedPlant,
    noisy_id::Symbol=:noisy)
    return (
        noisy=copy(acquisition_observation(
            prepared_acquisition(plant, noisy_id))),
        quiet=copy(acquisition_observation(
            prepared_acquisition(plant, :quiet_frame))),
    )
end

function rng_test_group(plant::PreparedPlant, id::Symbol)
    @inbounds for index in eachindex(prepared_acquisitions(plant))
        acquisition_id(prepared_acquisitions(plant)[index].definition) ==
            AcquisitionID(id) && return plant.rngs.acquisitions[index]
    end
    error("missing test acquisition RNG group")
end

@testset "Stable prepared plant RNG ownership" begin
    @test RNGDerivationVersion(1) == RNGDerivationVersion(UInt32(1))
    @test hash(RNGDerivationVersion(1)) ==
        hash(RNGDerivationVersion(UInt32(1)))
    @test_throws InvalidConfiguration RNGDerivationVersion(0)
    @test_throws InvalidConfiguration RNGOwnerIdentity(:rng, :owner,
        Symbol(""))
    @test hash(RNGOwnerIdentity(:rng, :owner, :state)) ==
        hash(RNGOwnerIdentity(:rng, :owner, :state))
    @test hash(AtmosphereLayerID(:ground)) ==
        hash(AtmosphereLayerID(:ground))

    definition = rng_test_definition()
    @test_throws UndefKeywordError prepare_plant(definition)
    assert_plant_preparation_error(
        () -> prepare_plant(definition; run_seed=-1),
        :rng, :invalid_seed)
    assert_plant_preparation_error(
        () -> prepare_plant(definition; run_seed=1.0),
        :rng, :invalid_seed)
    assert_plant_preparation_error(
        () -> prepare_plant(definition; run_seed=0x7100,
            rng_derivation_version=0),
        :rng, :invalid_derivation_version)
    assert_plant_preparation_error(
        () -> prepare_plant(definition; run_seed=0x7100,
            rng_derivation_version=1.0),
        :rng, :invalid_derivation_version)

    missing_layer_ids = rng_test_definition(; layer_ids=nothing)
    assert_plant_preparation_error(
        () -> prepare_plant(missing_layer_ids; run_seed=0x7100),
        :rng, :missing_owner_id)
    duplicate_layer_ids = rng_test_definition(
        layer_ids=(:duplicate, :duplicate))
    assert_plant_preparation_error(
        () -> prepare_plant(duplicate_layer_ids; run_seed=0x7100),
        :rng, :duplicate_owner_id)

    plant = prepare_plant(definition; run_seed=0x7100,
        rng_derivation_version=RNGDerivationVersion(1))
    selection = prepare_acquisition_selection(plant,
        (:quiet_frame, :noisy))
    metadata = rng_replay_metadata(plant)
    @test metadata.run_seed == UInt64(0x7100)
    @test metadata.derivation_version == UInt32(1)
    @test metadata.derivation_algorithm === :fnv1a_splitmix64_v1
    @test metadata.stream_algorithm === :xoshiro
    @test length(metadata.owners) == 8
    @test issorted(metadata.owners; by=Plant._rng_metadata_order)
    @test length(unique(map(owner -> owner.derived_seed,
        metadata.owners))) == length(metadata.owners)

    quiet_group = rng_test_group(plant, :quiet_frame)
    noisy_group = rng_test_group(plant, :noisy)
    quiet_before = copy(Plant.rng_stream_state(quiet_group,
        Val(:detector)))
    noisy_before = copy(Plant.rng_stream_state(noisy_group,
        Val(:detector)))
    @test @inferred(execute_acquisition_selection_at!(selection,
        0.01)) === selection
    observations = rng_test_observations(plant)
    quiet_after = copy(Plant.rng_stream_state(quiet_group,
        Val(:detector)))
    noisy_after = copy(Plant.rng_stream_state(noisy_group,
        Val(:detector)))
    @test quiet_before == quiet_after
    @test noisy_before != noisy_after

    replay = prepare_plant(rng_test_definition(); run_seed=0x7100,
        rng_derivation_version=1)
    replay_selection = prepare_acquisition_selection(replay,
        (:noisy, :quiet_frame))
    execute_acquisition_selection_at!(replay_selection, 0.01)
    @test rng_replay_metadata(replay) == metadata
    @test rng_test_observations(replay) == observations

    reordered = prepare_plant(rng_test_definition(
        path_order=(:quiet, :random),
        acquisition_order=(:quiet_frame, :noisy));
        run_seed=0x7100, rng_derivation_version=1)
    reordered_selection = prepare_acquisition_selection(reordered,
        (:quiet_frame, :noisy))
    execute_acquisition_selection_at!(reordered_selection, 0.01)
    @test rng_replay_metadata(reordered) == metadata
    @test rng_test_observations(reordered) == observations

    changed_seed = prepare_plant(rng_test_definition();
        run_seed=0x7101, rng_derivation_version=1)
    changed_version = prepare_plant(rng_test_definition();
        run_seed=0x7100, rng_derivation_version=2)
    renamed = prepare_plant(rng_test_definition(
        noisy_acquisition_id=:renamed_noisy,
        acquisition_order=(:noisy, :quiet_frame));
        run_seed=0x7100, rng_derivation_version=1)
    detector_seed = rng_test_owner_seed(metadata, :acquisition, :noisy,
        :detector)
    @test detector_seed == UInt64(0x324e0837e0347be4)
    @test rng_test_owner_seed(rng_replay_metadata(changed_seed),
        :acquisition, :noisy, :detector) != detector_seed
    @test rng_test_owner_seed(rng_replay_metadata(changed_version),
        :acquisition, :noisy, :detector) != detector_seed
    @test rng_test_owner_seed(rng_replay_metadata(renamed),
        :acquisition, :renamed_noisy, :detector) != detector_seed

    reordered_layers = prepare_plant(rng_test_definition(
        layer_order=(2, 1)); run_seed=0x7100,
        rng_derivation_version=1)
    reordered_layers_selection = prepare_acquisition_selection(
        reordered_layers, (:noisy, :quiet_frame))
    execute_acquisition_selection_at!(reordered_layers_selection, 0.01)
    @test path_input(prepared_path(reordered_layers, :random)).opd ≈
        path_input(prepared_path(plant, :random)).opd
    @test rng_test_owner_seed(rng_replay_metadata(reordered_layers),
        :atmosphere, :ground, :layer_state) ==
        rng_test_owner_seed(metadata, :atmosphere, :ground, :layer_state)

    single_owner = prepare_plant(rng_test_definition(
        atmosphere_model=:kolmogorov); run_seed=0x7100)
    single_owner_selection = prepare_acquisition_selection(single_owner,
        (:noisy,))
    @test execute_acquisition_selection_at!(single_owner_selection,
        0.01) === single_owner_selection

    if coverage_instrumented()
        @test_skip "prepared-RNG execution allocation assertion is disabled under coverage instrumentation"
    else
        @test prepared_selection_execution_allocations(selection,
            current_epoch(plant_atmosphere(plant.definition))) == 0
    end
end

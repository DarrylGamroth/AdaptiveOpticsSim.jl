struct SampledAberrationHandoffPathModel end

struct SampledAberrationAcquisitionModel{T<:AbstractFloat}
    exposure_s::T
end

struct SampledZeroPupilMaterialization{P}
    destination::P
end

struct SampledAberrationHandoffExecution{P,R}
    input::P
    result::R
end

Plant.plant_model_definition_style(
    ::Type{SampledAberrationHandoffPathModel}) =
    ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{<:SampledAberrationAcquisitionModel}) =
    ColdPlantModelDefinition()

function Plant.validate_path_materialization_binding(
    materialization::SampledZeroPupilMaterialization,
    input::PupilFunction,
    ::AdaptiveOpticsSim.Atmospheres.AbstractAtmosphere,
    ::AdaptiveOpticsSim.Optics.AbstractSource,
)
    materialization.destination === input || throw(
        PlantPreparationError(:path, :prepared_binding,
            "sampled-aberration zero materialization binding changed"))
    return nothing
end

function Plant.validate_path_materialization_target(
    materialization::SampledZeroPupilMaterialization,
    input::PupilFunction,
    ::AdaptiveOpticsSim.Atmospheres.AbstractAtmosphere,
    target::AdaptiveOpticsSim.Backends.AbstractComputeDevice,
)
    materialization.destination === input || throw(
        PlantPreparationError(:path, :prepared_binding,
            "sampled-aberration zero materialization binding changed"))
    Plant._require_exact_plant_product_target(
        materialization.destination, target,
        "sampled-aberration zero materialization destination")
    return materialization
end

function Plant.validate_path_materialization(
    materialization::SampledZeroPupilMaterialization,
    input::PupilFunction,
    ::AdaptiveOpticsSim.Atmospheres.AbstractAtmosphere,
    ::AtmosphereEpoch,
)
    materialization.destination === input || throw(
        PlantPreparationError(:path, :prepared_binding,
            "sampled-aberration zero materialization binding changed"))
    return input
end

function Plant.materialize_path_input!(
    materialization::SampledZeroPupilMaterialization,
    input::PupilFunction,
    ::AdaptiveOpticsSim.Atmospheres.AbstractAtmosphere,
    ::AtmosphereEpoch,
)
    materialization.destination === input || throw(
        PlantPreparationError(:path, :prepared_binding,
            "sampled-aberration zero materialization binding changed"))
    fill!(input.opd, zero(eltype(input.opd)))
    return input
end

function Plant.validate_path_execution_binding(
    execution::SampledAberrationHandoffExecution,
    input::PupilFunction,
    result::IntensityMap,
)
    execution.input === input && execution.result === result || throw(
        PlantPreparationError(:path, :prepared_binding,
            "sampled-aberration handoff binding changed"))
    return nothing
end

function Plant.validate_path_execution_target(
    execution::SampledAberrationHandoffExecution,
    target::AdaptiveOpticsSim.Backends.AbstractComputeDevice,
)
    Plant._require_exact_plant_product_target(
        execution.input, target, "sampled-aberration handoff input")
    Plant._require_exact_plant_product_target(
        execution.result, target, "sampled-aberration handoff result")
    return execution
end

function Plant.execute_path!(
    result::IntensityMap,
    input::PupilFunction,
    execution::SampledAberrationHandoffExecution,
)
    Plant.validate_path_execution_binding(execution, input, result)
    copyto!(result.values, input.opd)
    return result
end

function Plant.prepare_path_executor(
    ::SampledAberrationHandoffPathModel,
    definition::OpticalPathDefinition,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    telescope::Telescope,
    atmosphere::AdaptiveOpticsSim.Atmospheres.AbstractTimedAtmosphere,
    context,
)
    T = eltype(pupil_reflectivity(telescope))
    pupil = PupilFunction(telescope; T=T, backend=backend(telescope))
    values = similar(pupil.opd)
    fill!(values, zero(T))
    metadata = OpticalPlaneMetadata(FocalPlane(), values;
        coordinate_domain=MetricCoordinates(),
        sampling=pupil.metadata.sampling,
        origin=pupil.metadata.origin,
        centering=pupil.metadata.centering,
        orientation=pupil.metadata.orientation,
        spectral=MonochromaticChannel(T(wavelength(source))),
        normalization=PhotonRateNormalization(),
        spatial_measure=CellIntegratedMeasure(),
        coherence=IncoherentIntensityAddition())
    result = IntensityMap(metadata, values)
    execution = SampledAberrationHandoffExecution(pupil, result)
    return PreparedPathExecutor(
        definition,
        source,
        telescope,
        atmosphere,
        pupil,
        result,
        execution;
        context=context,
        materialization=SampledZeroPupilMaterialization(pupil),
        optical_model=:sampled_aberration_handoff,
        propagation_model=:external_typed_handoff,
        model_revisions=UInt(1),
    )
end

function Plant.prepare_acquisition_provider(
    model::SampledAberrationAcquisitionModel,
    ::AcquisitionDefinition,
    path::PreparedPathExecutor,
)
    Plant.require_path_result(path)
    T = eltype(path.result.values)
    detector = Detector(
        exposure_duration=T(model.exposure_s),
        noise=NoiseNone(),
        qe=one(T),
        gain=one(T),
        response_model=NullFrameResponse(),
        sensor=CMOSSensor(timing_model=GlobalShutter(), T=T),
        T=T,
        backend=backend(path.result),
    )
    execution = Plant.FrameAcquisitionExecution(detector, path.result)
    products = Plant.AcquisitionProducts(execution.observation;
        metadata=(
            kind=:sampled_aberration_test,
            units=:detected_electrons,
            geometry=path.result.metadata,
        ))
    return prepare_full_optical_provider(execution, products)
end

function sampled_aberration_test_metadata(
    prototype::PupilFunction,
    surface::AbstractMatrix,
)
    return OpticalPlaneMetadata(PupilPlane(), surface;
        coordinate_domain=MetricCoordinates(),
        sampling=prototype.metadata.sampling,
        orientation=prototype.metadata.orientation,
        spectral=AchromaticSpectralCoordinate(),
        normalization=DimensionlessNormalization(),
        spatial_measure=PointSampledMeasure(),
        coherence=NonCombinableProduct())
end

function sampled_aberration_test_atmosphere(telescope::Telescope)
    T = eltype(pupil_reflectivity(telescope))
    return MultiLayerAtmosphere(
        telescope;
        r0=T(0.2),
        L0=T(25),
        fractional_cn2=T[1],
        wind_speed=T[0],
        wind_direction_deg=T[0],
        altitude=T[0],
        layer_ids=(:ground,),
        T=T,
        backend=backend(telescope),
    )
end

function sampled_aberration_lowlevel_path(source::Source)
    T = typeof(wavelength(source))
    telescope = Telescope(
        resolution=5,
        diameter=T(5),
        central_obstruction=zero(T),
        T=T,
    )
    atmosphere = sampled_aberration_test_atmosphere(telescope)
    definition = PlantDefinition(;
        telescope=plant_test_telescope_definition(telescope),
        atmosphere=plant_test_atmosphere_definition(atmosphere),
        paths=(
            OpticalPathDefinition(
                :geometry, source, SampledAberrationHandoffPathModel()),
        ),
    )
    plant = prepare_plant(
        definition, PLANT_TEST_HOST_TARGET; run_seed=0x8700)
    return plant, prepared_path(plant, :geometry)
end

function sampled_aberration_definition_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

function sampled_aberration_dm_schema(::Type{T}) where {T<:AbstractFloat}
    return PlantCommandSchema(
        T,
        (4,);
        id=:event_dm_schema,
        version=1,
        endpoint=:event_dm_command,
        units=:metre,
        sign_convention=:positive_surface_increases_opd,
        basis=CommandBasis(:actuator, :event_dm),
        basis_revision=1,
        semantics=AbsoluteCommand,
        bounds=UnboundedCommandValues(),
        value_policy=CommandValuePolicy(),
        sequence_policy=CommandSequencePolicy(),
        effective_time_policy=CommandEffectiveTimePolicy(),
        silence_policy=CommandSilencePolicy(),
    )
end

function sampled_aberration_path_fixture()
    T = Float64
    telescope = Telescope(
        resolution=5,
        diameter=T(5),
        central_obstruction=zero(T),
        T=T,
    )
    atmosphere = sampled_aberration_test_atmosphere(telescope)
    source = Source(
        band=:custom,
        wavelength=T(0.8e-6),
        photon_irradiance=T(100),
        T=T,
    )
    prototype = PupilFunction(telescope; T=T)
    common_opd = similar(prototype.opd)
    science_opd = similar(prototype.opd)
    fill!(common_opd, T(1))
    fill!(science_opd, T(2))
    common_metadata =
        sampled_aberration_test_metadata(prototype, common_opd)
    science_metadata =
        sampled_aberration_test_metadata(prototype, science_opd)
    common = SampledAberrationDefinition(
        :common_static,
        OPDMap(common_opd),
        common_metadata;
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
        application=DMReplace(),
    )
    science = SampledAberrationDefinition(
        :science_ncpa,
        NCPA(science_opd),
        science_metadata;
        placement=PupilPlanePlacement(),
        visibility=SelectedPathVisibility(:science),
        application=DMAdditive(),
    )
    paths = (
        OpticalPathDefinition(
            :science, source, SampledAberrationHandoffPathModel()),
        OpticalPathDefinition(
            :wfs, source, SampledAberrationHandoffPathModel()),
        OpticalPathDefinition(
            :other_science, source,
            SampledAberrationHandoffPathModel()),
    )
    acquisitions = (
        AcquisitionDefinition(
            :science_camera, :science,
            SampledAberrationAcquisitionModel(T(0.01))),
        AcquisitionDefinition(
            :wfs_camera, :wfs,
            SampledAberrationAcquisitionModel(T(0.01))),
        AcquisitionDefinition(
            :other_camera, :other_science,
            SampledAberrationAcquisitionModel(T(0.01))),
    )
    definition = PlantDefinition(;
        telescope=plant_test_telescope_definition(telescope),
        atmosphere=plant_test_atmosphere_definition(atmosphere),
        sampled_aberrations=(science, common),
        paths,
        acquisitions,
    )
    plant = prepare_plant(
        definition, PLANT_TEST_HOST_TARGET; run_seed=0x8701)
    selection = prepare_acquisition_selection(
        plant, (:science_camera, :wfs_camera, :other_camera))
    return (;
        telescope,
        atmosphere=prepared_atmosphere(plant),
        definition,
        plant,
        selection,
        common_opd,
        science_opd,
    )
end

@testset "Sampled-aberration declarations and preparation" begin
    @test Base.isexported(Plant, :SampledAberrationID)
    @test Base.isexported(Plant, :SampledAberrationDefinition)
    @test !Base.isexported(AdaptiveOpticsSim, :SampledAberrationID)
    @test !Base.isexported(
        AdaptiveOpticsSim, :SampledAberrationDefinition)

    fixture = sampled_aberration_path_fixture()
    definition = fixture.definition
    plant = fixture.plant
    @test Tuple(map(Plant.sampled_aberration_id,
        Plant.sampled_aberration_definitions(definition))) == (
        SampledAberrationID(:science_ncpa),
        SampledAberrationID(:common_static),
    )
    @test Plant.sampled_aberration_definition(
        definition, :science_ncpa) isa SampledAberrationDefinition
    @test Plant.prepared_sampled_aberration(
        plant, :common_static) isa Plant.PreparedSampledAberration

    bindings = Plant.prepared_sampled_aberration_path_bindings(plant)
    @test Plant.prepared_sampled_aberration_path_count(bindings) == 3
    @test Plant.prepared_sampled_aberration_binding_count(bindings) == 4
    @test map(
        ordinal -> Plant.prepared_sampled_aberration_path_id(
            bindings, ordinal),
        1:3,
    ) == [
        OpticalPathID(:other_science),
        OpticalPathID(:science),
        OpticalPathID(:wfs),
    ]
    science_range =
        Plant.prepared_sampled_aberration_binding_range(bindings, :science)
    @test length(science_range) == 2
    @test map(science_range) do binding
        slot = Plant.prepared_sampled_aberration_slot(bindings, binding)
        Plant.sampled_aberration_id(
            Plant.prepared_sampled_aberrations(plant)[slot])
    end == [
        SampledAberrationID(:common_static),
        SampledAberrationID(:science_ncpa),
    ]
    @test length(Plant.prepared_sampled_aberration_binding_range(
        bindings, :wfs)) == 1

    prepared_common = Plant.sampled_aberration_opd(
        Plant.prepared_sampled_aberration(plant, :common_static))
    prepared_science = Plant.sampled_aberration_opd(
        Plant.prepared_sampled_aberration(plant, :science_ncpa))
    @test prepared_common !== fixture.common_opd
    @test prepared_science !== fixture.science_opd
    @test all(==(1.0), prepared_common)
    @test all(==(2.0), prepared_science)
    fill!(fixture.common_opd, 17.0)
    fill!(fixture.science_opd, 19.0)
    @test all(==(1.0), prepared_common)
    @test all(==(2.0), prepared_science)

    @test @inferred(execute_acquisition_selection_at!(
        fixture.selection, 0.0)) === fixture.selection
    @test all(==(3.0), path_input(
        prepared_path(plant, :science)).opd)
    @test all(==(1.0), path_input(prepared_path(plant, :wfs)).opd)
    @test all(==(1.0), path_input(
        prepared_path(plant, :other_science)).opd)
    @test all(==(3.0), path_result(
        prepared_path(plant, :science)).values)
    @test all(==(1.0), path_result(prepared_path(plant, :wfs)).values)

    unrelated_product = Ref(:unrelated_path_product)
    empty_plan = Plant._PreparedSampledAberrationPathPlan(())
    @test Plant._apply_sampled_aberration_path_plan!(
        unrelated_product, empty_plan) === unrelated_product
    science_pupil = path_input(prepared_path(plant, :science))
    @test Plant._apply_sampled_aberration_path_plan!(
        science_pupil, empty_plan) === science_pupil

    if coverage_instrumented()
        @test_skip "sampled-aberration allocation gate disabled under " *
            "coverage instrumentation"
    else
        epoch = current_epoch(fixture.atmosphere)
        allocation_bytes = prepared_selection_execution_allocations(
            fixture.selection, epoch)
        @test allocation_bytes == 0
    end
end

@testset "Sampled-aberration structured rejection" begin
    T = Float64
    telescope = Telescope(
        resolution=5, diameter=T(5), central_obstruction=zero(T), T=T)
    atmosphere = sampled_aberration_test_atmosphere(telescope)
    source = Source(
        band=:custom,
        wavelength=T(0.8e-6),
        photon_irradiance=T(100),
        T=T,
    )
    prototype = PupilFunction(telescope; T=T)
    values = similar(prototype.opd)
    fill!(values, one(T))
    metadata = sampled_aberration_test_metadata(prototype, values)
    telescope = plant_test_telescope_definition(telescope)
    atmosphere = plant_test_atmosphere_definition(atmosphere)

    registered = SampledAberrationDefinition(
        :registered,
        OPDMap(values),
        metadata;
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
        application=DMAdditive(),
        registration=PupilRelayRegistration(),
    )
    @test Plant.sampled_aberration_registration(registered) isa
        PupilRelayRegistration

    cases = (
        (
            () -> SampledAberrationDefinition(
                :invalid_surface,
                values,
                metadata;
                placement=PupilPlanePlacement(),
                visibility=AllPathVisibility(),
                application=DMAdditive(),
            ),
            :invalid_surface,
        ),
        (
            () -> SampledAberrationDefinition(
                :focal,
                OPDMap(values),
                metadata;
                placement=FocalPlanePlacement(),
                visibility=AllPathVisibility(),
                application=DMAdditive(),
            ),
            :unsupported_placement,
        ),
        (
            () -> SampledAberrationDefinition(
                :invalid_application,
                OPDMap(values),
                metadata;
                placement=PupilPlanePlacement(),
                visibility=AllPathVisibility(),
                application=:add,
            ),
            :invalid_application,
        ),
        (
            () -> SampledAberrationDefinition(
                :invalid_registration,
                OPDMap(values),
                metadata;
                placement=PupilPlanePlacement(),
                visibility=AllPathVisibility(),
                application=DMAdditive(),
                registration=:invalid,
            ),
            :invalid_pupil_relay_registration,
        ),
        (
            () -> SampledAberrationDefinition(
                :invalid_metadata,
                OPDMap(values),
                :invalid;
                placement=PupilPlanePlacement(),
                visibility=AllPathVisibility(),
                application=DMAdditive(),
            ),
            :invalid_surface_metadata,
        ),
        (
            () -> SampledAberrationDefinition(
                "invalid_id",
                OPDMap(values),
                metadata;
                placement=PupilPlanePlacement(),
                visibility=AllPathVisibility(),
                application=DMAdditive(),
            ),
            :invalid_id,
        ),
    )
    for (operation, reason) in cases
        error = sampled_aberration_definition_error(operation)
        @test error isa PlantDefinitionError
        @test error.component === :sampled_aberration
        @test error.reason === reason
    end

    path = OpticalPathDefinition(
        :science, source, SampledAberrationHandoffPathModel())
    first_replace = SampledAberrationDefinition(
        :first_replace,
        OPDMap(values),
        metadata;
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
        application=DMReplace(),
    )
    named_definition = PlantDefinition(;
        telescope,
        atmosphere,
        sampled_aberrations=(first_replace=first_replace,),
        paths=(science=path,),
    )
    @test Tuple(Plant.sampled_aberration_definitions(named_definition)) ==
        (first_replace,)
    vector_definition = PlantDefinition(;
        telescope,
        atmosphere,
        sampled_aberrations=[first_replace],
        paths=[path],
    )
    @test Plant.sampled_aberration_definitions(vector_definition)[1] ===
        first_replace
    sampled_registry =
        Plant.sampled_aberration_definitions(vector_definition)
    @test getfield(sampled_registry, :_storage) isa Tuple
    @test_throws CanonicalIndexError setindex!(
        sampled_registry, first_replace, 1)
    @test_throws MethodError setindex!(
        getfield(sampled_registry, :_storage), first_replace, 1)

    topology_cases = (
        (
            () -> PlantDefinition(;
                telescope,
                atmosphere,
                sampled_aberrations=(nothing,),
                paths=(path,),
            ),
            :invalid_definition,
        ),
        (
            () -> PlantDefinition(;
                telescope,
                atmosphere,
                sampled_aberrations=(wrong=first_replace,),
                paths=(path,),
            ),
            :identity_mismatch,
        ),
    )
    for (operation, reason) in topology_cases
        topology_error = sampled_aberration_definition_error(operation)
        @test topology_error isa PlantDefinitionError
        @test topology_error.component === :sampled_aberration
        @test topology_error.reason === reason
    end

    unknown_definition_error = sampled_aberration_definition_error() do
        Plant.sampled_aberration_definition(named_definition, :missing)
    end
    @test unknown_definition_error isa PlantDefinitionError
    @test unknown_definition_error.component === :sampled_aberration
    @test unknown_definition_error.reason === :unknown_id

    named_plant = prepare_plant(
        named_definition, PLANT_TEST_HOST_TARGET; run_seed=0x8702)
    unknown_prepared_error = sampled_aberration_definition_error() do
        Plant.prepared_sampled_aberration(named_plant, :missing)
    end
    @test unknown_prepared_error isa PlantPreparationError
    @test unknown_prepared_error.component === :sampled_aberration
    @test unknown_prepared_error.reason === :unknown_id

    unknown_binding_error = sampled_aberration_definition_error() do
        Plant.prepared_sampled_aberration_binding_range(
            Plant.prepared_sampled_aberration_path_bindings(named_plant),
            :missing,
        )
    end
    @test unknown_binding_error isa PlantPreparationError
    @test unknown_binding_error.component === :path
    @test unknown_binding_error.reason === :unknown_id

    second_replace = SampledAberrationDefinition(
        :second_replace,
        OPDMap(values),
        metadata;
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
        application=DMReplace(),
    )
    ambiguous = PlantDefinition(;
        telescope,
        atmosphere,
        sampled_aberrations=(second_replace, first_replace),
        paths=(path,),
    )
    assert_plant_preparation_error(
        () -> prepare_plant(
            ambiguous, PLANT_TEST_HOST_TARGET; run_seed=0x8702),
        :sampled_aberration,
        :ambiguous_replacement_order,
    )

    unknown = SampledAberrationDefinition(
        :unknown_path,
        OPDMap(values),
        metadata;
        placement=PupilPlanePlacement(),
        visibility=SelectedPathVisibility(:missing),
        application=DMAdditive(),
    )
    error = sampled_aberration_definition_error() do
        PlantDefinition(;
            telescope,
            atmosphere,
            sampled_aberrations=(unknown,),
            paths=(path,),
        )
    end
    @test error isa PlantDefinitionError
    @test error.component === :sampled_aberration
    @test error.reason === :unknown_visible_path

    duplicate_error = sampled_aberration_definition_error() do
        PlantDefinition(;
            telescope,
            atmosphere,
            sampled_aberrations=(first_replace, first_replace),
            paths=(path,),
        )
    end
    @test duplicate_error isa PlantDefinitionError
    @test duplicate_error.component === :sampled_aberration
    @test duplicate_error.reason === :duplicate_id
end

@testset "Replacement finite support and event execution order" begin
    source = Source(
        band=:custom,
        wavelength=0.8e-6,
        photon_irradiance=100.0,
    )
    _, path = sampled_aberration_lowlevel_path(source)
    pupil = path_input(path)
    surface = similar(pupil.opd, 3, 3)
    metadata = sampled_aberration_test_metadata(pupil, surface)
    fill!(surface, 2.0)
    coupling = prepare_sampled_pupil_footprint_coupling(
        metadata, surface, path, PupilPlanePlacement())
    @test coupling isa PreparedPupilFootprintCoupling

    fill!(pupil.opd, 7.0)
    @test apply_sampled_pupil_surface!(
        pupil, surface, coupling, DMReplace()) === pupil
    expected = zeros(5, 5)
    expected[2:4, 2:4] .= 2.0
    @test pupil.opd == expected

    fill!(pupil.opd, 7.0)
    apply_sampled_pupil_surface!(
        pupil, surface, coupling, DMAdditive())
    additive_expected = fill(7.0, 5, 5)
    additive_expected[2:4, 2:4] .= 9.0
    @test pupil.opd == additive_expected
    if !coverage_instrumented()
        @test @allocated(apply_sampled_pupil_surface!(
            pupil, surface, coupling, DMReplace())) == 0
        @test @allocated(apply_sampled_pupil_surface!(
            pupil, surface, coupling, DMAdditive())) == 0
    end

    T = Float64
    telescope = Telescope(
        resolution=5,
        diameter=T(5),
        central_obstruction=zero(T),
        T=T,
    )
    atmosphere = sampled_aberration_test_atmosphere(telescope)
    event_source = Source(
        band=:custom,
        wavelength=T(0.8e-6),
        photon_irradiance=T(100),
        T=T,
    )
    prototype = PupilFunction(telescope; T=T)
    static_opd = similar(prototype.opd)
    fill!(static_opd, T(2))
    static_metadata =
        sampled_aberration_test_metadata(prototype, static_opd)
    static = SampledAberrationDefinition(
        :event_static,
        OPDMap(static_opd),
        static_metadata;
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
        application=DMReplace(),
    )
    event_path = OpticalPathDefinition(
        :science, event_source, SampledAberrationHandoffPathModel())
    acquisition = AcquisitionDefinition(
        :camera,
        :science,
        SampledAberrationAcquisitionModel(T(0.02)),
    )
    optic = ControllableOpticDefinition(
        :event_dm,
        DeformableMirrorModel(
            n_act=2, influence_width=T(0.3), T=T),
        (sampled_aberration_dm_schema(T),);
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
    )
    plant = prepare_plant(
        PlantDefinition(;
            telescope=plant_test_telescope_definition(telescope),
            atmosphere=plant_test_atmosphere_definition(atmosphere),
            controllable_optics=(optic,),
            sampled_aberrations=(static,),
            paths=(event_path,),
            acquisitions=(acquisition,),
        ), PLANT_TEST_HOST_TARGET;
        run_seed=0x8703,
        command_endpoints=(
            CommandEndpointConfiguration(
                :event_dm_command, T[1e-9, 0, 0, 0]; capacity=2),
        ),
    )
    event_loop = prepare_plant_event_loop(
        plant,
        PlantEventLoopDefinition(
            (
                OpticalSampleDefinition(
                    :science,
                    PeriodicSchedule(period_ns=100_000_000),
                ),
            ),
            (
                AcquisitionEventDefinition(
                    :camera,
                    GlobalShutterAcquisitionDefinition(
                        PlantDuration(20_000_000)),
                    PeriodicAcquisitionStart(
                        PeriodicSchedule(period_ns=1_000_000_000)),
                ),
            ),
        ),
    )
    state = PlantEventLoopState(event_loop)
    workspace = PlantEventLoopWorkspace(event_loop)
    @test Plant.plant_event_sampled_aberration_count(event_loop) == 1
    @test step_plant_events!(event_loop, state, workspace) ==
        PlantTimestamp(0)
    dm_surface = surface_opd(only(state.controllable_optics).active)
    expected_event_opd = T(2) .+ dm_surface
    @test path_input(prepared_path(plant, :science)).opd ≈
        expected_event_opd
    @test path_result(prepared_path(plant, :science)).values ≈
        expected_event_opd
end

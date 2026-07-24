struct PlacedOpticsTestPathModel end
struct PlacedOpticsTestAcquisitionModel{T<:AbstractFloat}
    exposure_s::T
end
struct PlacedOpticsTestOpticModel end

struct PreparedPlacedOpticsTestOptic
    endpoint::CommandEndpointID
end

mutable struct PlacedOpticsTestOpticState{T<:AbstractFloat}
    visible::T
end

mutable struct PlacedOpticsTestOpticWorkspace{T<:AbstractFloat}
    staged::T
end

struct UnsupportedPlacedOpticPlacement <:
       AbstractControllableOpticPlacement end
struct UnsupportedPlacedOpticVisibility <:
       AbstractControllableOpticVisibility end

Plant.plant_model_definition_style(
    ::Type{PlacedOpticsTestPathModel}) = ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{<:PlacedOpticsTestAcquisitionModel}) =
    ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{PlacedOpticsTestOpticModel}) = ColdPlantModelDefinition()

function Plant.prepare_controllable_optic(
    ::PlacedOpticsTestOpticModel,
    definition::ControllableOpticDefinition,
    ::AbstractTelescope,
    ::AbstractAtmosphere,
)
    return PreparedPlacedOpticsTestOptic(
        only(command_endpoint_ids(definition)))
end

function Plant.prepare_controllable_optic_state(
    prepared::PreparedPlacedOpticsTestOptic,
    ::ControllableOpticDefinition,
    endpoint_ids::Tuple,
    initial_commands::Tuple,
)
    only(endpoint_ids) == prepared.endpoint || throw(
        PlantPreparationError(:controllable_optic, :prepared_binding,
            "placed-optics test endpoint binding changed"))
    initial = only(initial_commands)
    return PlacedOpticsTestOpticState(float(initial))
end

function Plant.prepare_controllable_optic_workspace(
    ::PreparedPlacedOpticsTestOptic)
    return PlacedOpticsTestOpticWorkspace(0.0)
end

function Plant.stage_controllable_optic_command!(
    prepared::PreparedPlacedOpticsTestOptic,
    ::PlacedOpticsTestOpticState,
    workspace::PlacedOpticsTestOpticWorkspace,
    endpoint::CommandEndpointID,
    command,
    ::PlantTimestamp,
)
    endpoint == prepared.endpoint || throw(PlantCommandError(
        :physical_application, :endpoint_mismatch,
        "placed-optics test endpoint changed"))
    workspace.staged = command
    return nothing
end

function Plant.commit_controllable_optic_command!(
    ::PreparedPlacedOpticsTestOptic,
    state::PlacedOpticsTestOpticState,
    workspace::PlacedOpticsTestOpticWorkspace,
    ::CommandEndpointID,
    ::PlantTimestamp,
)
    state.visible = workspace.staged
    return nothing
end

function Plant.apply_controllable_optic_surface!(
    input::PupilFunction,
    ::PreparedPlacedOpticsTestOptic,
    state::PlacedOpticsTestOpticState,
)
    @. input.opd += state.visible
    return input
end

function Plant.prepare_path_executor(
    ::PlacedOpticsTestPathModel,
    definition::OpticalPathDefinition,
    source::AbstractSource,
    telescope::Telescope,
    atmosphere::AbstractTimedAtmosphere,
)
    T = eltype(pupil_reflectivity(telescope))
    pupil = PupilFunction(telescope; T, backend=backend(telescope))
    imaging = prepare_direct_imaging(pupil, source; zero_padding=1)
    return PreparedPathExecutor(
        definition,
        source,
        telescope,
        atmosphere,
        pupil,
        direct_imaging_output(imaging),
        imaging;
        materialization=prepare_pupil_opd_materialization(
            atmosphere, telescope, source, pupil),
        optical_model=:placed_optics_test,
        propagation_model=:fraunhofer_fft,
        model_revisions=UInt(1),
    )
end

function Plant.prepare_acquisition_provider(
    model::PlacedOpticsTestAcquisitionModel,
    ::AcquisitionDefinition,
    path::PreparedPathExecutor,
)
    require_path_result(path)
    T = eltype(path.result.values)
    detector = Detector(
        integration_time=T(model.exposure_s),
        noise=NoiseNone(),
        qe=one(T),
        gain=one(T),
        response_model=NullFrameResponse(),
        sensor=CMOSSensor(timing_model=GlobalShutter(), T=T),
        T=T,
        backend=path.key.backend,
    )
    execution = FrameAcquisitionExecution(detector, path.result)
    products = AcquisitionProducts(execution.observation;
        metadata=(kind=:placed_optics_test,
            units=:detected_electrons,
            geometry=path.result.metadata))
    return prepare_full_optical_provider(execution, products)
end

function placed_optics_test_schema(id::Symbol)
    endpoint = Symbol(id, :_command)
    return PlantCommandSchema(
        Float64,
        ();
        id=Symbol(endpoint, :_schema),
        version=1,
        endpoint,
        units=:metre,
        sign_convention=:positive_surface_increases_opd,
        basis=CommandBasis(:modal, id),
        basis_revision=1,
        semantics=AbsoluteCommand,
        bounds=UniformCommandBounds(-10.0, 10.0),
        value_policy=CommandValuePolicy(),
        sequence_policy=CommandSequencePolicy(),
        effective_time_policy=CommandEffectiveTimePolicy(),
        silence_policy=CommandSilencePolicy(),
    )
end

function placed_optics_test_definition(id::Symbol, placement, visibility)
    return ControllableOpticDefinition(
        id,
        PlacedOpticsTestOpticModel(),
        (placed_optics_test_schema(id),);
        placement,
        visibility,
    )
end

function placed_optics_test_components(; reverse::Bool=false)
    T = Float64
    telescope = Telescope(
        resolution=4,
        diameter=T(4),
        central_obstruction=zero(T),
        T=T,
    )
    atmosphere = MultiLayerAtmosphere(
        telescope;
        r0=T(0.2),
        L0=T(25),
        fractional_cn2=T[1],
        wind_speed=T[0],
        wind_direction=T[0],
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
    paths = (
        OpticalPathDefinition(:science, source, PlacedOpticsTestPathModel()),
        OpticalPathDefinition(:lgs, source, PlacedOpticsTestPathModel()),
        OpticalPathDefinition(:ngs, source, PlacedOpticsTestPathModel()),
    )
    optics = (
        placed_optics_test_definition(
            :d_high_all,
            AtmosphericConjugatePlacement(12_000.0),
            AllPathVisibility(),
        ),
        placed_optics_test_definition(
            :a_science_pupil,
            PupilPlanePlacement(),
            SelectedPathVisibility(:science),
        ),
        placed_optics_test_definition(
            :e_mid_science,
            AtmosphericConjugatePlacement(5_000.0),
            SelectedPathVisibility(:science),
        ),
        placed_optics_test_definition(
            :b_pupil_all,
            PupilPlanePlacement(),
            AllPathVisibility(),
        ),
        placed_optics_test_definition(
            :c_high_wfs,
            AtmosphericConjugatePlacement(12_000.0f0),
            reverse ? SelectedPathVisibility(:science, :ngs) :
                SelectedPathVisibility(:ngs, :science),
        ),
    )
    configurations = map(optics) do optic
        CommandEndpointConfiguration(
            only(command_endpoint_ids(optic)),
            0.0;
            capacity=1,
        )
    end
    if reverse
        paths = Base.reverse(paths)
        optics = Base.reverse(optics)
        configurations = Base.reverse(configurations)
    end
    definition = PlantDefinition(;
        telescope,
        atmosphere,
        controllable_optics=optics,
        paths,
    )
    plant = prepare_plant(
        definition;
        run_seed=0x8501,
        command_endpoints=configurations,
    )
    return (; plant, definition)
end

function placed_optics_binding_ids(plant::PreparedPlant, path)
    bindings = prepared_controllable_optic_path_bindings(plant)
    optics = prepared_controllable_optics(plant)
    return map(
        prepared_controllable_optic_binding_range(bindings, path),
    ) do binding
        slot = prepared_controllable_optic_slot(bindings, binding)
        controllable_optic_id(optics[slot].definition)
    end
end

function placed_optics_group_ids(plant::PreparedPlant, path)
    bindings = prepared_controllable_optic_path_bindings(plant)
    optics = prepared_controllable_optics(plant)
    return map(
        prepared_controllable_optic_plane_group_range(bindings, path),
    ) do group_index
        group = prepared_controllable_optic_plane_group(
            bindings, group_index)
        map(prepared_controllable_optic_plane_group_binding_range(group)) do binding
            slot = prepared_controllable_optic_slot(bindings, binding)
            controllable_optic_id(optics[slot].definition)
        end
    end
end

@inline function placed_optics_binding_checksum(
    bindings::PreparedControllableOpticPathBindings,
    path::OpticalPathID,
)
    checksum = 0
    @inbounds for binding in
        prepared_controllable_optic_binding_range(bindings, path)
        checksum += prepared_controllable_optic_slot(bindings, binding)
    end
    return checksum
end

@inline function placed_optics_binding_allocations(bindings, path)
    return @allocated placed_optics_binding_checksum(bindings, path)
end

function placed_optics_captured_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

function placed_optics_event_fixture()
    T = Float64
    telescope = Telescope(
        resolution=4,
        diameter=T(4),
        central_obstruction=zero(T),
        T=T,
    )
    atmosphere = MultiLayerAtmosphere(
        telescope;
        r0=T(0.2),
        L0=T(25),
        fractional_cn2=T[1],
        wind_speed=T[0],
        wind_direction=T[0],
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
    paths = (
        OpticalPathDefinition(:science, source, PlacedOpticsTestPathModel()),
        OpticalPathDefinition(:ngs, source, PlacedOpticsTestPathModel()),
    )
    acquisitions = (
        AcquisitionDefinition(:science_camera, :science,
            PlacedOpticsTestAcquisitionModel(T(0.02))),
        AcquisitionDefinition(:ngs_camera, :ngs,
            PlacedOpticsTestAcquisitionModel(T(0.02))),
    )
    common = placed_optics_test_definition(
        :common,
        PupilPlanePlacement(),
        AllPathVisibility(),
    )
    science = placed_optics_test_definition(
        :science_only,
        PupilPlanePlacement(),
        SelectedPathVisibility(:science),
    )
    plant = prepare_plant(
        PlantDefinition(;
            telescope,
            atmosphere,
            controllable_optics=(science, common),
            paths,
            acquisitions,
        );
        run_seed=0x8502,
        command_endpoints=(
            CommandEndpointConfiguration(:science_only_command, 2.0;
                capacity=1),
            CommandEndpointConfiguration(:common_command, 1.0;
                capacity=1),
        ),
    )
    loop = prepare_plant_event_loop(
        plant,
        PlantEventLoopDefinition(
            (
                OpticalSampleDefinition(:ngs,
                    PeriodicSchedule(period_ns=100_000_000)),
                OpticalSampleDefinition(:science,
                    PeriodicSchedule(period_ns=100_000_000)),
            ),
            (
                AcquisitionEventDefinition(
                    :ngs_camera,
                    GlobalShutterAcquisitionDefinition(
                        PlantDuration(20_000_000)),
                    PeriodicAcquisitionStart(
                        PeriodicSchedule(period_ns=1_000_000_000)),
                ),
                AcquisitionEventDefinition(
                    :science_camera,
                    GlobalShutterAcquisitionDefinition(
                        PlantDuration(20_000_000)),
                    PeriodicAcquisitionStart(
                        PeriodicSchedule(period_ns=1_000_000_000)),
                ),
            ),
        ),
    )
    return plant, loop, PlantEventLoopState(loop),
        PlantEventLoopWorkspace(loop)
end

@testset "Placed-optic declarations and structural errors" begin
    for name in (
        :PupilPlanePlacement,
        :AtmosphericConjugatePlacement,
        :FocalPlanePlacement,
        :AllPathVisibility,
        :SelectedPathVisibility,
    )
        @test Base.isexported(Plant, name)
        @test !Base.isexported(AdaptiveOpticsSim, name)
    end
    for name in (
        :controllable_optic_placement,
        :controllable_optic_visibility,
        :prepared_controllable_optic_path_bindings,
        :prepared_controllable_optic_binding_range,
    )
        @test Base.ispublic(Plant, name)
        @test !Base.isexported(Plant, name)
    end

    selected_input = Symbol[:science, :ngs]
    visibility = SelectedPathVisibility(selected_input)
    @test @inferred(SelectedPathVisibility(:science, :ngs)) ==
        visibility
    selected_input[1] = :changed
    @test selected_path_ids(visibility) ==
        (OpticalPathID(:ngs), OpticalPathID(:science))
    @test selected_path_ids(SelectedPathVisibility(:science, :lgs, :ngs)) ==
        OpticalPathID.((:lgs, :ngs, :science))
    @test !Base.ismutable(visibility)
    @test @inferred(conjugate_altitude_m(
        AtmosphericConjugatePlacement(12_000.0f0))) === 12_000.0f0
    @test AtmosphericConjugatePlacement(12_000.0f0) ==
        AtmosphericConjugatePlacement(12_000.0)
    @test isequal(
        AtmosphericConjugatePlacement(12_000.0f0),
        AtmosphericConjugatePlacement(12_000.0),
    )
    @test hash(AtmosphericConjugatePlacement(12_000.0f0)) ==
        hash(AtmosphericConjugatePlacement(12_000.0))
    @test !signbit(conjugate_altitude_m(
        AtmosphericConjugatePlacement(-0.0)))

    schema = placed_optics_test_schema(:required)
    definition = @inferred ControllableOpticDefinition(
        :required,
        PlacedOpticsTestOpticModel(),
        (schema,);
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
    )
    @test controllable_optic_placement(definition) isa
        PupilPlanePlacement
    @test controllable_optic_visibility(definition) isa
        AllPathVisibility

    missing_placement = placed_optics_captured_error() do
        ControllableOpticDefinition(
            :missing_placement,
            PlacedOpticsTestOpticModel(),
            (schema,);
            visibility=AllPathVisibility(),
        )
    end
    @test missing_placement isa UndefKeywordError

    for (operation, reason) in (
        (() -> AtmosphericConjugatePlacement(-1.0),
            :invalid_conjugate_altitude),
        (() -> AtmosphericConjugatePlacement(Inf),
            :invalid_conjugate_altitude),
        (() -> SelectedPathVisibility(()), :empty_path_visibility),
        (() -> SelectedPathVisibility(:science, :science),
            :duplicate_visible_path),
        (() -> ControllableOpticDefinition(
                :bad_placement,
                PlacedOpticsTestOpticModel(),
                (placed_optics_test_schema(:bad_placement),);
                placement=UnsupportedPlacedOpticPlacement(),
                visibility=AllPathVisibility(),
            ), :invalid_placement),
        (() -> ControllableOpticDefinition(
                :bad_visibility,
                PlacedOpticsTestOpticModel(),
                (placed_optics_test_schema(:bad_visibility),);
                placement=PupilPlanePlacement(),
                visibility=UnsupportedPlacedOpticVisibility(),
            ), :invalid_path_visibility),
    )
        error = placed_optics_captured_error(operation)
        @test error isa PlantDefinitionError
        @test error.component === :controllable_optic
        @test error.reason === reason
    end

    components = placed_optics_test_components()
    unknown = placed_optics_test_definition(
        :unknown_path,
        PupilPlanePlacement(),
        SelectedPathVisibility(:missing),
    )
    unknown_error = placed_optics_captured_error() do
        PlantDefinition(;
            telescope=plant_telescope(components.definition),
            atmosphere=plant_atmosphere(components.definition),
            controllable_optics=(unknown,),
            paths=path_definitions(components.definition),
        )
    end
    @test unknown_error isa PlantDefinitionError
    @test unknown_error.component === :controllable_optic
    @test unknown_error.reason === :unknown_visible_path
end

@testset "Canonical path bindings and co-placed groups" begin
    forward = placed_optics_test_components()
    reversed = placed_optics_test_components(; reverse=true)
    plant = forward.plant
    bindings = @inferred prepared_controllable_optic_path_bindings(plant)

    @test prepared_controllable_optic_path_count(bindings) == 3
    @test prepared_controllable_optic_binding_count(bindings) == 10
    @test prepared_controllable_optic_plane_group_count(bindings) == 7
    @test Tuple(map(
        ordinal -> prepared_controllable_optic_path_id(bindings, ordinal),
        1:prepared_controllable_optic_path_count(bindings),
    )) == OpticalPathID.((:lgs, :ngs, :science))

    expected = (
        lgs=(
            ControllableOpticID(:b_pupil_all),
            ControllableOpticID(:d_high_all),
        ),
        ngs=(
            ControllableOpticID(:b_pupil_all),
            ControllableOpticID(:c_high_wfs),
            ControllableOpticID(:d_high_all),
        ),
        science=(
            ControllableOpticID(:a_science_pupil),
            ControllableOpticID(:b_pupil_all),
            ControllableOpticID(:e_mid_science),
            ControllableOpticID(:c_high_wfs),
            ControllableOpticID(:d_high_all),
        ),
    )
    for path in keys(expected)
        @test Tuple(placed_optics_binding_ids(plant, path)) ==
            getproperty(expected, path)
        @test placed_optics_binding_ids(plant, path) ==
            placed_optics_binding_ids(reversed.plant, path)
        @test placed_optics_group_ids(plant, path) ==
            placed_optics_group_ids(reversed.plant, path)
    end

    @test Tuple(Tuple.(placed_optics_group_ids(plant, :science))) == (
        (
            ControllableOpticID(:a_science_pupil),
            ControllableOpticID(:b_pupil_all),
        ),
        (ControllableOpticID(:e_mid_science),),
        (
            ControllableOpticID(:c_high_wfs),
            ControllableOpticID(:d_high_all),
        ),
    )
    science_slot =
        prepared_controllable_optic_path_slot(bindings, :science)
    for group_index in
        prepared_controllable_optic_plane_group_range(bindings, :science)
        group = prepared_controllable_optic_plane_group(
            bindings, group_index)
        @test prepared_controllable_optic_plane_group_path_slot(group) ==
            science_slot
        @test !isempty(
            prepared_controllable_optic_plane_group_binding_range(group))
    end

    path = OpticalPathID(:science)
    @test @inferred(placed_optics_binding_checksum(bindings, path)) == 15
    placed_optics_binding_checksum(bindings, path)
    @test placed_optics_binding_allocations(bindings, path) == 0
end

@testset "Selected-path execution uses prepared bindings" begin
    plant, loop, state, workspace = placed_optics_event_fixture()
    @test step_plant_events!(loop, state, workspace) == PlantTimestamp(0)
    science = path_input(prepared_path(plant, :science))
    ngs = path_input(prepared_path(plant, :ngs))
    @test all(isapprox.(science.opd .- ngs.opd, 2.0;
        rtol=0.0, atol=8eps(Float64)))
end

@testset "Unimplemented conjugate geometry fails before execution" begin
    components = placed_optics_test_components()
    plant = components.plant
    acquisition_definition = AcquisitionDefinition(
        :science_camera,
        :science,
        PlacedOpticsTestAcquisitionModel(0.02),
    )
    definition = components.definition
    with_acquisition = PlantDefinition(;
        telescope=plant_telescope(definition),
        atmosphere=plant_atmosphere(definition),
        controllable_optics=controllable_optic_definitions(definition),
        paths=path_definitions(definition),
        acquisitions=(acquisition_definition,),
    )
    configurations = map(
        controllable_optic_definitions(with_acquisition),
    ) do optic
        CommandEndpointConfiguration(
            only(command_endpoint_ids(optic)), 0.0; capacity=1)
    end
    prepared = prepare_plant(
        with_acquisition;
        run_seed=0x8503,
        command_endpoints=configurations,
    )
    error = placed_optics_captured_error() do
        prepare_plant_event_loop(
            prepared,
            PlantEventLoopDefinition(
                (OpticalSampleDefinition(
                    :science,
                    PeriodicSchedule(period_ns=100_000_000),
                ),),
                (AcquisitionEventDefinition(
                    :science_camera,
                    GlobalShutterAcquisitionDefinition(
                        PlantDuration(20_000_000)),
                    PeriodicAcquisitionStart(
                        PeriodicSchedule(period_ns=1_000_000_000)),
                ),),
            ),
        )
    end
    @test error isa PlantScheduleError
    @test error.reason === :unsupported_conjugate_geometry
end

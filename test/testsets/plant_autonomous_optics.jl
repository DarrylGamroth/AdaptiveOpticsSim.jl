struct AutonomousPyramidPathModel{T<:AbstractFloat}
    radius::T
    phase_offset::T
    samples::Int
end

struct AutonomousPyramidAcquisitionModel{T<:AbstractFloat}
    exposure::T
end

struct UnsupportedAutonomousPhaseReference <:
    Plant.AbstractWaveformPhaseReference end

Plant.plant_model_definition_style(
    ::Type{<:AutonomousPyramidPathModel}) = ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{<:AutonomousPyramidAcquisitionModel}) =
    ColdPlantModelDefinition()

function Plant.prepare_path_executor(
    model::AutonomousPyramidPathModel,
    definition::OpticalPathDefinition,
    source::AbstractSource,
    telescope::Telescope,
    atmosphere::AdaptiveOpticsSim.AbstractTimedAtmosphere,
)
    T = eltype(pupil_reflectivity(telescope))
    pupil = PupilFunction(telescope; T=T, backend=backend(telescope))
    sensor = PyramidWFS(telescope;
        pupil_samples=2,
        mode=Diffractive(),
        modulation=model.radius,
        modulation_points=model.samples,
        delta_theta=model.phase_offset,
        calib_modulation=model.radius,
        diffraction_padding=2,
        T=T,
        backend=backend(telescope),
    )
    front_end = PyramidOpticalFrontEnd(sensor, source)
    output = pyramid_rate_map(front_end, pupil)
    plan = prepare_wfs_optical_formation(front_end, pupil, output)
    execution = WFSOpticalPathExecution(plan)
    return PreparedPathExecutor(
        definition,
        source,
        telescope,
        atmosphere,
        pupil,
        output,
        execution;
        materialization=AtmosphereIndependentPath(),
        optical_model=(
            kind=:autonomous_pyramid,
            pupil_samples=2,
            modulation_samples=model.samples,
        ),
        propagation_model=:pyramid_focal_plane_mask,
        model_revisions=UInt(1),
    )
end

function Plant.prepare_acquisition_provider(
    model::AutonomousPyramidAcquisitionModel,
    ::AcquisitionDefinition,
    path::PreparedPathExecutor,
)
    require_path_result(path)
    T = eltype(path.result.values)
    detector = Detector(
        integration_time=T(model.exposure),
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
        metadata=(
            kind=:autonomous_pyramid_frame,
            units=:detected_electrons,
            geometry=path.result.metadata,
            detector=detector_export_metadata(detector),
        ))
    return prepare_full_optical_provider(execution, products)
end

function autonomous_pyramid_schema(::Type{T}, endpoint::Symbol;
    dimensions=(),
    units::Symbol,
    sign_convention::Symbol,
    basis=CommandBasis(:waveform_setpoint, :pyramid_modulation),
    semantics=AbsoluteCommand,
    bounds,
    silence_policy=CommandSilencePolicy(),
) where {T<:Real}
    return PlantCommandSchema(
        T,
        dimensions;
        id=Symbol(endpoint, :_schema),
        version=1,
        endpoint,
        units,
        sign_convention,
        basis,
        basis_revision=1,
        semantics,
        bounds,
        value_policy=CommandValuePolicy(),
        sequence_policy=CommandSequencePolicy(),
        effective_time_policy=CommandEffectiveTimePolicy(
            supersession=PreservePendingCommands),
        silence_policy,
    )
end

function autonomous_pyramid_schemas(;
    prefix::Symbol=:mod,
    radius_units::Symbol=:lambda_over_d,
    radius_silence_policy=CommandSilencePolicy(),
)
    T = Float64
    radius_endpoint = Symbol(prefix, :_radius)
    frequency_endpoint = Symbol(prefix, :_frequency)
    phase_endpoint = Symbol(prefix, :_phase)
    enabled_endpoint = Symbol(prefix, :_enabled)
    return (
        radius=autonomous_pyramid_schema(T, radius_endpoint;
            units=radius_units,
            sign_convention=:nonnegative_radius,
            bounds=UniformCommandBounds(T(0), T(8)),
            silence_policy=radius_silence_policy),
        frequency=autonomous_pyramid_schema(T, frequency_endpoint;
            units=:hertz,
            sign_convention=:nonnegative_frequency,
            bounds=UniformCommandBounds(T(0), T(2_000))),
        phase=autonomous_pyramid_schema(T, phase_endpoint;
            units=:radian,
            sign_convention=:counterclockwise_phase,
            bounds=UniformCommandBounds(T(-4pi), T(4pi))),
        enabled=autonomous_pyramid_schema(UInt8, enabled_endpoint;
            units=:binary_state,
            sign_convention=:zero_disabled_one_enabled,
            bounds=UniformCommandBounds(UInt8(0), UInt8(1))),
    )
end

function autonomous_pyramid_schema_like(schema;
    T=command_numeric_type(schema),
    dimensions=command_dimensions(schema),
    units=command_units(schema),
    sign_convention=command_sign_convention(schema),
    basis=command_basis(schema),
    semantics=command_semantics(schema),
    bounds=command_bounds(schema),
)
    return PlantCommandSchema(
        T,
        dimensions;
        id=command_schema_id(schema),
        version=command_schema_version(schema),
        endpoint=command_endpoint_id(schema),
        units,
        sign_convention,
        basis,
        basis_revision=command_basis_revision(schema),
        semantics,
        bounds,
        value_policy=command_value_policy(schema),
        sequence_policy=command_sequence_policy(schema),
        effective_time_policy=command_effective_time_policy(schema),
        silence_policy=command_silence_policy(schema),
    )
end

function autonomous_pyramid_optic(schemas=autonomous_pyramid_schemas();
    id=:pyramid_modulator,
    visibility=SelectedPathVisibility(:pyramid))
    model = CircularPyramidModulator(
        radius_endpoint=command_endpoint_id(schemas.radius),
        frequency_endpoint=command_endpoint_id(schemas.frequency),
        phase_endpoint=command_endpoint_id(schemas.phase),
        enabled_endpoint=command_endpoint_id(schemas.enabled),
    )
    return ControllableOpticDefinition(id, model,
        Tuple(schemas);
        placement=FocalPlanePlacement(),
        visibility)
end

function autonomous_pyramid_event_definition(;
    phase_reference=FreeRunningPhaseReference(),
    trigger_topology=nothing,
    acquisition_start=PeriodicAcquisitionStart(
        PeriodicSchedule(period_ns=1_000_000_000, phase_ns=500_000_000)),
    exposure_ns=20_000_000,
    sample_period_ns=100_000_000,
    optic=:pyramid_modulator,
    path=:pyramid,
)
    return PlantEventLoopDefinition(
        (OpticalSampleDefinition(path,
            PeriodicSchedule(period_ns=sample_period_ns, phase_ns=0)),),
        (AcquisitionEventDefinition(
            :camera,
            GlobalShutterAcquisitionDefinition(
                PlantDuration(exposure_ns)),
            acquisition_start,
        ),);
        trigger_topology,
        autonomous_optics=(
            AutonomousPeriodicOpticDefinition(
                optic,
                path;
                phase_reference,
            ),
        ),
    )
end

function autonomous_pyramid_fixture(;
    radius=1.25,
    frequency_hz=5.0,
    phase_offset=0.25,
    enabled=UInt8(1),
    radius_silence_policy=CommandSilencePolicy(),
    safe_radius=nothing,
    phase_reference=FreeRunningPhaseReference(),
    trigger_topology=nothing,
    acquisition_start=PeriodicAcquisitionStart(
        PeriodicSchedule(period_ns=1_000_000_000, phase_ns=500_000_000)),
    exposure_ns=20_000_000,
    sample_period_ns=100_000_000,
)
    T = Float64
    telescope = Telescope(resolution=4, diameter=T(4),
        central_obstruction=zero(T), T=T)
    atmosphere = MultiLayerAtmosphere(telescope;
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
        wavelength=T(0.75e-6),
        photon_irradiance=T(100),
        T=T,
    )
    schemas = autonomous_pyramid_schemas(; radius_silence_policy)
    optic = autonomous_pyramid_optic(schemas)
    path = OpticalPathDefinition(:pyramid, source,
        AutonomousPyramidPathModel(T(radius), T(phase_offset), 8))
    acquisition = AcquisitionDefinition(:camera, :pyramid,
        AutonomousPyramidAcquisitionModel(T(exposure_ns) / T(1e9)))
    definition = PlantDefinition(;
        telescope,
        atmosphere,
        controllable_optics=(optic,),
        paths=(path,),
        acquisitions=(acquisition,),
    )
    configurations = (
        CommandEndpointConfiguration(:mod_radius, T(radius);
            capacity=4, safe_command=safe_radius),
        CommandEndpointConfiguration(:mod_frequency, T(frequency_hz);
            capacity=4),
        CommandEndpointConfiguration(:mod_phase, T(phase_offset);
            capacity=4),
        CommandEndpointConfiguration(:mod_enabled, enabled;
            capacity=4),
    )
    plant = prepare_plant(definition;
        run_seed=0x7c00,
        command_endpoints=configurations,
    )
    event_definition = autonomous_pyramid_event_definition(;
        phase_reference,
        trigger_topology,
        acquisition_start,
        exposure_ns,
        sample_period_ns,
    )
    prepared = prepare_plant_event_loop(plant, event_definition)
    return (;
        plant,
        prepared,
        state=PlantEventLoopState(prepared),
        workspace=PlantEventLoopWorkspace(prepared),
        schemas,
    )
end

function autonomous_pyramid_submit!(
    fixture, schema, sequence, effective_ns, value;
    admission_ns=0,
)
    return admit_plant_command!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        PlantCommand(
            schema,
            sequence,
            PlantTimestamp(effective_ns),
            value,
        ),
        PlantTimestamp(admission_ns),
    )
end

function autonomous_pyramid_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

function autonomous_pyramid_model_preparation_error(fixture, schemas)
    definition = autonomous_pyramid_optic(schemas)
    plant_definition = fixture.plant.definition
    return autonomous_pyramid_error() do
        Plant.prepare_controllable_optic(
            Plant.controllable_optic_model(definition),
            definition,
            plant_telescope(plant_definition),
            plant_atmosphere(plant_definition),
        )
    end
end

function autonomous_pyramid_event_preparation_error(fixture; kwargs...)
    return autonomous_pyramid_error() do
        prepare_plant_event_loop(
            fixture.plant,
            autonomous_pyramid_event_definition(; kwargs...),
        )
    end
end

@inline function autonomous_pyramid_evaluation_allocations(
    implementation, state, coupling, timestamp)
    return @allocated evaluate_autonomous_periodic_optic!(
        implementation, state, coupling, timestamp)
end

@inline function autonomous_pyramid_phase_allocations(
    implementation, state, timestamp)
    return @allocated autonomous_waveform_phase(
        implementation, state, timestamp)
end

@testset "Autonomous Pyramid construction and preparation contract" begin
    @test Base.isexported(Plant, :AutonomousPeriodicOpticDefinition)
    @test Base.isexported(Plant, :CircularPyramidModulator)
    @test !Base.isexported(AdaptiveOpticsSim,
        :AutonomousPeriodicOpticDefinition)

    duplicate_error = autonomous_pyramid_error() do
        CircularPyramidModulator(
            radius_endpoint=:same,
            frequency_endpoint=:same,
            phase_endpoint=:phase,
            enabled_endpoint=:enabled,
        )
    end
    @test duplicate_error isa PlantDefinitionError
    @test duplicate_error.reason == :duplicate_endpoint

    endpoint = CommandEndpointID(:same)
    positional_duplicate_error = autonomous_pyramid_error() do
        CircularPyramidModulator(
            endpoint,
            endpoint,
            CommandEndpointID(:phase),
            CommandEndpointID(:enabled),
        )
    end
    @test positional_duplicate_error isa PlantDefinitionError
    @test positional_duplicate_error.reason == :duplicate_endpoint

    invalid_schemas = autonomous_pyramid_schemas(
        radius_units=:metre)
    fixture = autonomous_pyramid_fixture()
    definition = PlantDefinition(;
        telescope=plant_telescope(fixture.plant.definition),
        atmosphere=plant_atmosphere(fixture.plant.definition),
        controllable_optics=(autonomous_pyramid_optic(invalid_schemas;
            visibility=AllPathVisibility()),),
    )
    configurations = (
        CommandEndpointConfiguration(:mod_radius, 1.0; capacity=1),
        CommandEndpointConfiguration(:mod_frequency, 5.0; capacity=1),
        CommandEndpointConfiguration(:mod_phase, 0.0; capacity=1),
        CommandEndpointConfiguration(:mod_enabled, UInt8(1); capacity=1),
    )
    schema_error = autonomous_pyramid_error() do
        prepare_plant(definition; run_seed=0x7c01,
            command_endpoints=configurations)
    end
    @test schema_error isa PlantPreparationError
    @test schema_error.reason == :command_units

    schema_cases = (
        (:command_shape, :radius, (; dimensions=(1,))),
        (:command_semantics, :radius,
            (; semantics=IncrementalCommand)),
        (:command_sign_convention, :radius,
            (; sign_convention=CommandSignConvention(:signed_radius))),
        (:command_basis, :radius,
            (; basis=CommandBasis(:actuator_command, :tip_tilt))),
        (:unbounded_setpoint, :radius,
            (; bounds=UnboundedCommandValues())),
        (:negative_setpoint_bound, :radius,
            (; bounds=UniformCommandBounds(-1.0, 8.0))),
        (:enabled_type, :enabled,
            (; T=Float64, bounds=UniformCommandBounds(0.0, 1.0))),
        (:enabled_bounds, :enabled,
            (; bounds=UniformCommandBounds(UInt8(0), UInt8(2)))),
        (:numeric_type, :frequency,
            (; T=Float32,
                bounds=UniformCommandBounds(0.0f0, 2_000.0f0))),
    )
    for (reason, role, overrides) in schema_cases
        replacement = autonomous_pyramid_schema_like(
            getproperty(fixture.schemas, role);
            overrides...,
        )
        schemas = merge(
            fixture.schemas,
            NamedTuple{(role,)}((replacement,)),
        )
        error = autonomous_pyramid_model_preparation_error(
            fixture, schemas)
        @test error isa PlantPreparationError
        @test error.reason == reason
    end

    extra_schema = autonomous_pyramid_schema(
        Float64, :mod_extra;
        units=:radian,
        sign_convention=:counterclockwise_phase,
        bounds=UniformCommandBounds(-1.0, 1.0),
    )
    base_optic = autonomous_pyramid_optic(fixture.schemas)
    extra_definition = ControllableOpticDefinition(
        :pyramid_modulator,
        Plant.controllable_optic_model(base_optic),
        (command_schemas(base_optic)..., extra_schema);
        placement=FocalPlanePlacement(),
        visibility=SelectedPathVisibility(:pyramid),
    )
    count_error = autonomous_pyramid_error() do
        plant_definition = fixture.plant.definition
        Plant.prepare_controllable_optic(
            Plant.controllable_optic_model(extra_definition),
            extra_definition,
            plant_telescope(plant_definition),
            plant_atmosphere(plant_definition),
        )
    end
    @test count_error isa PlantPreparationError
    @test count_error.reason == :command_schema_count

    missing_binding = PlantEventLoopDefinition(
        (OpticalSampleDefinition(:pyramid,
            PeriodicSchedule(period_ns=100_000_000)),),
        (AcquisitionEventDefinition(
            :camera,
            GlobalShutterAcquisitionDefinition(
                PlantDuration(20_000_000)),
            PeriodicAcquisitionStart(
                PeriodicSchedule(period_ns=1_000_000_000)),
        ),),
    )
    binding_error = autonomous_pyramid_error() do
        prepare_plant_event_loop(fixture.plant, missing_binding)
    end
    @test binding_error isa PlantScheduleError
    @test binding_error.reason == :missing_autonomous_binding

    unsupported_reference_error =
        autonomous_pyramid_event_preparation_error(
            fixture;
            phase_reference=UnsupportedAutonomousPhaseReference(),
        )
    @test unsupported_reference_error isa PlantScheduleError
    @test unsupported_reference_error.reason == :unsupported_phase_reference

    for reference in (
        TriggerSourcePhaseReference(:sync),
        TriggerResetPhaseReference(:modulation_reset),
    )
        topology_error = autonomous_pyramid_event_preparation_error(
            fixture;
            phase_reference=reference,
        )
        @test topology_error isa PlantScheduleError
        @test topology_error.reason == :missing_trigger_topology
    end

    source = TriggerSourceDefinition(
        :sync,
        PeriodicSchedule(period_ns=100_000_000),
    )
    camera_consumer = TriggerConsumerDefinition(
        :camera_trigger,
        TriggerSourceID(:sync),
    )
    topology = prepare_trigger_topology(
        (source,),
        (),
        (camera_consumer,);
        in_flight_capacity=2,
    )
    triggered_start = TriggeredAcquisitionStart(:camera_trigger)

    trigger_binding_cases = (
        (TriggerSourcePhaseReference(:missing_source),
            :unknown_trigger_source),
        (TriggerResetPhaseReference(:missing_consumer),
            :unknown_trigger_consumer),
        (TriggerResetPhaseReference(:camera_trigger),
            :duplicate_trigger_consumer),
    )
    for (reference, reason) in trigger_binding_cases
        error = autonomous_pyramid_event_preparation_error(
            fixture;
            phase_reference=reference,
            trigger_topology=topology,
            acquisition_start=triggered_start,
        )
        @test error isa PlantScheduleError
        @test error.reason == reason
    end

    unknown_optic_error = autonomous_pyramid_event_preparation_error(
        fixture;
        optic=:missing_modulator,
    )
    @test unknown_optic_error isa PlantScheduleError
    @test unknown_optic_error.reason == :unknown_controllable_optic

    unknown_inspection_error = autonomous_pyramid_error() do
        autonomous_waveform_radius(
            fixture.prepared,
            fixture.state,
            :missing_modulator,
        )
    end
    @test unknown_inspection_error isa PlantScheduleError
    @test unknown_inspection_error.reason == :unknown_autonomous_optic

    second_schemas = autonomous_pyramid_schemas(prefix=:second)
    second_optic = autonomous_pyramid_optic(second_schemas;
        id=:second_modulator)
    base_definition = fixture.plant.definition
    conflicting_definition = PlantDefinition(;
        telescope=plant_telescope(base_definition),
        atmosphere=plant_atmosphere(base_definition),
        controllable_optics=(
            only(controllable_optic_definitions(base_definition)),
            second_optic,
        ),
        paths=path_definitions(base_definition),
        acquisitions=acquisition_definitions(base_definition),
    )
    conflicting_configurations = (
        CommandEndpointConfiguration(:mod_radius, 1.25; capacity=1),
        CommandEndpointConfiguration(:mod_frequency, 5.0; capacity=1),
        CommandEndpointConfiguration(:mod_phase, 0.25; capacity=1),
        CommandEndpointConfiguration(:mod_enabled, UInt8(1); capacity=1),
        CommandEndpointConfiguration(:second_radius, 1.0; capacity=1),
        CommandEndpointConfiguration(:second_frequency, 4.0; capacity=1),
        CommandEndpointConfiguration(:second_phase, 0.0; capacity=1),
        CommandEndpointConfiguration(:second_enabled, UInt8(1);
            capacity=1),
    )
    conflicting_plant = prepare_plant(conflicting_definition;
        run_seed=0x7c02,
        command_endpoints=conflicting_configurations)
    conflicting_loop = PlantEventLoopDefinition(
        (OpticalSampleDefinition(:pyramid,
            PeriodicSchedule(period_ns=100_000_000)),),
        (AcquisitionEventDefinition(
            :camera,
            GlobalShutterAcquisitionDefinition(
                PlantDuration(20_000_000)),
            PeriodicAcquisitionStart(
                PeriodicSchedule(period_ns=1_000_000_000)),
        ),);
        autonomous_optics=(
            AutonomousPeriodicOpticDefinition(
                :pyramid_modulator, :pyramid;
                phase_reference=FreeRunningPhaseReference()),
            AutonomousPeriodicOpticDefinition(
                :second_modulator, :pyramid;
                phase_reference=FreeRunningPhaseReference()),
        ),
    )
    conflict_error = autonomous_pyramid_error() do
        prepare_plant_event_loop(conflicting_plant, conflicting_loop)
    end
    @test conflict_error isa PlantScheduleError
    @test conflict_error.reason == :conflicting_autonomous_coupling

    @test plant_event_autonomous_optic_count(fixture.prepared) == 1
end

@testset "Cycle-averaged Pyramid optical equivalence and setpoints" begin
    fixture = autonomous_pyramid_fixture()
    path = prepared_path(fixture.plant, :pyramid)
    modulation = path.execution.plan.front_end.modulation
    initial_phases = copy(modulation.phases)
    @test sum(abs2, modulation.amplitude_weights) ≈ 1.0

    execute_path!(path)
    frozen_output = copy(path_result(path).values)
    @test step_plant_events!(
        fixture.prepared, fixture.state, fixture.workspace) ==
        PlantTimestamp(0)
    @test modulation.phases == initial_phases
    @test path_result(path).values ≈ frozen_output atol=1e-12 rtol=1e-12

    autonomous_pyramid_submit!(
        fixture, fixture.schemas.radius, 1, 100_000_000, 2.0;
        admission_ns=1)
    @test step_plant_events!(
        fixture.prepared, fixture.state, fixture.workspace) ==
        PlantTimestamp(100_000_000)
    @test autonomous_waveform_radius(
        fixture.prepared, fixture.state, :pyramid_modulator) == 2.0
    @test modulation.phases != initial_phases
    @test path_result(path).values != frozen_output

    disabled = autonomous_pyramid_fixture()
    disabled_path = prepared_path(disabled.plant, :pyramid)
    disabled_modulation =
        disabled_path.execution.plan.front_end.modulation
    autonomous_pyramid_submit!(
        disabled, disabled.schemas.enabled, 1, 0, UInt8(0))
    @test step_plant_events!(
        disabled.prepared, disabled.state, disabled.workspace) ==
        PlantTimestamp(0)
    @test !autonomous_waveform_enabled(
        disabled.prepared, disabled.state, :pyramid_modulator)
    @test all(==(1.0 + 0.0im), disabled_modulation.phases)
    @test autonomous_waveform_offset(
        disabled.prepared, disabled.state, :pyramid_modulator,
        PlantTimestamp(75_000_000)) == (0.0, 0.0)

    clear_command_dispositions!(disabled.workspace)
    autonomous_pyramid_submit!(
        disabled, disabled.schemas.enabled, 2, 100_000_000, UInt8(1);
        admission_ns=1,
    )
    @test step_plant_events!(
        disabled.prepared, disabled.state, disabled.workspace) ==
        PlantTimestamp(100_000_000)
    @test disabled_modulation.phases == initial_phases
end

@testset "Analytic phase and equal-time command ordering" begin
    float32_implementation =
        PreparedCircularPyramidModulator{Float32,UInt8}(
            CommandEndpointID(:radius),
            CommandEndpointID(:frequency),
            CommandEndpointID(:phase),
            CommandEndpointID(:enabled),
        )
    origin = zero(PlantTimestamp)
    float32_state = CircularPyramidModulatorState(
        1.0f0,
        997.25f0,
        0.25f0,
        true,
        origin,
        0.0f0,
        origin,
        UInt64(0),
        UInt64(0),
        false,
    )
    @test @inferred(autonomous_waveform_phase(
        float32_implementation,
        float32_state,
        PlantTimestamp(1_000_000_000_000),
    )) == 0.25f0

    unsupported_initialization_error = autonomous_pyramid_error() do
        initialize_autonomous_periodic_optic!(
            float32_implementation,
            float32_state,
            UnsupportedAutonomousPhaseReference(),
        )
    end
    @test unsupported_initialization_error isa PlantPreparationError
    @test unsupported_initialization_error.reason ==
        :unsupported_phase_reference

    reference = FreeRunningPhaseReference(
        origin=PlantTimestamp(20_000_000))
    fixture = autonomous_pyramid_fixture(; phase_reference=reference)
    expected = mod(0.25 + 2pi * 5.0 * 0.05, 2pi)
    @test autonomous_waveform_phase(
        fixture.prepared, fixture.state, :pyramid_modulator,
        PlantTimestamp(70_000_000)) ≈ expected
    expected_offset = (1.25 * cos(expected), 1.25 * sin(expected))
    observed_offset = autonomous_waveform_offset(
        fixture.prepared, fixture.state, :pyramid_modulator,
        PlantTimestamp(70_000_000))
    @test observed_offset[1] ≈ expected_offset[1]
    @test observed_offset[2] ≈ expected_offset[2]

    path = prepared_path(fixture.plant, :pyramid)
    @test step_plant_events!(
        fixture.prepared, fixture.state, fixture.workspace) ==
        PlantTimestamp(0)
    baseline = copy(path_result(path).values)
    autonomous_pyramid_submit!(
        fixture, fixture.schemas.frequency, 1, 100_000_000, 10.0;
        admission_ns=1,
    )
    autonomous_pyramid_submit!(
        fixture, fixture.schemas.phase, 1, 100_000_000, 0.5;
        admission_ns=1,
    )
    @test step_plant_events!(
        fixture.prepared, fixture.state, fixture.workspace) ==
        PlantTimestamp(100_000_000)
    expected_at_command = mod(
        2pi * 5.0 * 0.08 + 0.5,
        2pi,
    )
    @test autonomous_waveform_phase(
        fixture.prepared, fixture.state, :pyramid_modulator,
        PlantTimestamp(100_000_000)) ≈ expected_at_command
    @test autonomous_waveform_phase(
        fixture.prepared, fixture.state, :pyramid_modulator,
        PlantTimestamp(150_000_000)) ≈
        mod(expected_at_command + 2pi * 10.0 * 0.05, 2pi)
    @test path_result(path).values ≈ baseline atol=1e-12 rtol=1e-12
end

@testset "Trigger source and delivered-reset phase relationships" begin
    source = TriggerSourceDefinition(:sync,
        PeriodicSchedule(period_ns=100_000_000, phase_ns=0))
    source_camera = TriggerConsumerDefinition(
        :camera_trigger, TriggerSourceID(:sync))
    source_topology = prepare_trigger_topology(
        (source,), (), (source_camera,); in_flight_capacity=2)
    source_fixture = autonomous_pyramid_fixture(
        phase_reference=TriggerSourcePhaseReference(:sync),
        trigger_topology=source_topology,
        acquisition_start=TriggeredAcquisitionStart(:camera_trigger),
        sample_period_ns=10_000_000,
    )
    source_modulation = prepared_path(
        source_fixture.plant, :pyramid).execution.plan.front_end.modulation
    initial_source_phases = copy(source_modulation.phases)
    autonomous_pyramid_submit!(
        source_fixture, source_fixture.schemas.radius,
        1, 100_000_000, 2.0)
    autonomous_pyramid_submit!(
        source_fixture, source_fixture.schemas.frequency,
        1, 100_000_000, 10.0)
    autonomous_pyramid_submit!(
        source_fixture, source_fixture.schemas.phase,
        1, 100_000_000, 0.5)
    run_plant_events_until!(
        source_fixture.prepared,
        source_fixture.state,
        source_fixture.workspace,
        PlantTimestamp(100_000_000),
    )
    @test autonomous_waveform_reference_timestamp(
        source_fixture.prepared, source_fixture.state,
        :pyramid_modulator) == PlantTimestamp(100_000_000)
    @test autonomous_waveform_radius(
        source_fixture.prepared, source_fixture.state,
        :pyramid_modulator) == 2.0
    @test autonomous_waveform_phase(
        source_fixture.prepared, source_fixture.state,
        :pyramid_modulator, PlantTimestamp(100_000_000)) ≈ 0.5
    @test autonomous_waveform_phase(
        source_fixture.prepared, source_fixture.state,
        :pyramid_modulator, PlantTimestamp(150_000_000)) ≈
        mod(0.5 + 2pi * 10.0 * 0.05, 2pi)
    @test source_modulation.phases != initial_source_phases
    run_plant_events_until!(
        source_fixture.prepared,
        source_fixture.state,
        source_fixture.workspace,
        PlantTimestamp(200_000_000),
    )
    @test autonomous_waveform_reference_count(
        source_fixture.prepared, source_fixture.state,
        :pyramid_modulator) == UInt64(3)
    @test autonomous_waveform_reference_sequence(
        source_fixture.prepared, source_fixture.state,
        :pyramid_modulator) == UInt64(3)
    @test autonomous_waveform_reference_timestamp(
        source_fixture.prepared, source_fixture.state,
        :pyramid_modulator) == PlantTimestamp(200_000_000)

    dropped_source = TriggerSourceDefinition(:faulted_sync,
        PeriodicSchedule(period_ns=100_000_000, phase_ns=0);
        faults=TriggerFaultTrace(TriggerFaultTraceEntry(
            2, :source_drop; action=DropTriggerEdge)))
    dropped_camera = TriggerConsumerDefinition(
        :faulted_camera_trigger, TriggerSourceID(:faulted_sync))
    dropped_topology = prepare_trigger_topology(
        (dropped_source,), (), (dropped_camera,); in_flight_capacity=2)
    dropped_fixture = autonomous_pyramid_fixture(
        phase_reference=TriggerSourcePhaseReference(:faulted_sync),
        trigger_topology=dropped_topology,
        acquisition_start=TriggeredAcquisitionStart(
            :faulted_camera_trigger),
        sample_period_ns=10_000_000,
    )
    run_plant_events_until!(
        dropped_fixture.prepared,
        dropped_fixture.state,
        dropped_fixture.workspace,
        PlantTimestamp(225_000_000),
    )
    @test acquisition_product_sequence(
        dropped_fixture.prepared, dropped_fixture.state, :camera) ==
        UInt64(2)
    @test autonomous_waveform_reference_count(
        dropped_fixture.prepared, dropped_fixture.state,
        :pyramid_modulator) == UInt64(3)
    @test autonomous_waveform_reference_sequence(
        dropped_fixture.prepared, dropped_fixture.state,
        :pyramid_modulator) == UInt64(3)
    @test autonomous_waveform_reference_timestamp(
        dropped_fixture.prepared, dropped_fixture.state,
        :pyramid_modulator) == PlantTimestamp(200_000_000)

    link_faults = TriggerFaultTrace(
        TriggerFaultTraceEntry(2, :modulation_drop;
            action=DropTriggerEdge),
        TriggerFaultTraceEntry(3, :modulation_phase_error;
            phase_step=PlantTimeOffset(5_000_000),
            jitter=PlantTimeOffset(3_000_000)),
        TriggerFaultTraceEntry(4, :modulation_duplicate;
            action=DuplicateTriggerEdge,
            duplicate_delay=PlantDuration(5_000_000)),
    )
    link = TriggerLinkDefinition(
        :modulation_link,
        TriggerSourceID(:sync);
        propagation_delay=PlantDuration(10_000_000),
        timing_skew=PlantTimeOffset(2_000_000),
        timestamp_label_offset=PlantTimeOffset(100_000_000),
        faults=link_faults,
    )
    camera = TriggerConsumerDefinition(
        :camera_trigger, TriggerSourceID(:sync))
    modulator = TriggerConsumerDefinition(
        :modulation_reset, TriggerLinkID(:modulation_link))
    reset_topology = prepare_trigger_topology(
        (source,), (link,), (camera, modulator);
        in_flight_capacity=4)
    reset_fixture = autonomous_pyramid_fixture(
        phase_reference=TriggerResetPhaseReference(:modulation_reset),
        trigger_topology=reset_topology,
        acquisition_start=TriggeredAcquisitionStart(:camera_trigger),
        sample_period_ns=10_000_000,
    )
    run_plant_events_until!(
        reset_fixture.prepared,
        reset_fixture.state,
        reset_fixture.workspace,
        PlantTimestamp(225_000_000),
    )
    @test acquisition_product_sequence(
        reset_fixture.prepared, reset_fixture.state, :camera) ==
        UInt64(3)
    @test autonomous_waveform_reference_count(
        reset_fixture.prepared, reset_fixture.state,
        :pyramid_modulator) == UInt64(2)
    @test autonomous_waveform_reference_sequence(
        reset_fixture.prepared, reset_fixture.state,
        :pyramid_modulator) == UInt64(3)
    @test autonomous_waveform_reference_timestamp(
        reset_fixture.prepared, reset_fixture.state,
        :pyramid_modulator) == PlantTimestamp(220_000_000)
    @test autonomous_waveform_phase(
        reset_fixture.prepared, reset_fixture.state,
        :pyramid_modulator, PlantTimestamp(225_000_000)) ≈
        mod(0.25 + 2pi * 5.0 * 0.005, 2pi)
    run_plant_events_until!(
        reset_fixture.prepared,
        reset_fixture.state,
        reset_fixture.workspace,
        PlantTimestamp(325_000_000),
    )
    @test acquisition_product_sequence(
        reset_fixture.prepared, reset_fixture.state, :camera) ==
        UInt64(4)
    @test autonomous_waveform_reference_count(
        reset_fixture.prepared, reset_fixture.state,
        :pyramid_modulator) == UInt64(4)
    @test autonomous_waveform_reference_sequence(
        reset_fixture.prepared, reset_fixture.state,
        :pyramid_modulator) == UInt64(4)
    @test autonomous_waveform_reference_timestamp(
        reset_fixture.prepared, reset_fixture.state,
        :pyramid_modulator) == PlantTimestamp(322_000_000)
end

@testset "Autonomous setpoint command-silence response" begin
    safe_policy = CommandSilencePolicy(
        ApplySafeCommand,
        AgeFromApplication;
        timeout=PlantDuration(150_000_000),
    )
    fixture = autonomous_pyramid_fixture(
        radius_silence_policy=safe_policy,
        safe_radius=0.4,
    )
    path = prepared_path(fixture.plant, :pyramid)
    modulation = path.execution.plan.front_end.modulation
    @test step_plant_events!(
        fixture.prepared, fixture.state, fixture.workspace) ==
        PlantTimestamp(0)
    initial_phases = copy(modulation.phases)
    @test step_plant_events!(
        fixture.prepared, fixture.state, fixture.workspace) ==
        PlantTimestamp(100_000_000)
    @test step_plant_events!(
        fixture.prepared, fixture.state, fixture.workspace) ==
        PlantTimestamp(150_000_000)
    @test autonomous_waveform_radius(
        fixture.prepared, fixture.state, :pyramid_modulator) == 0.4
    @test modulation.phases == initial_phases
    @test step_plant_events!(
        fixture.prepared, fixture.state, fixture.workspace) ==
        PlantTimestamp(200_000_000)
    @test modulation.phases != initial_phases
end

@testset "Autonomous waveform bounded storage and allocation budget" begin
    fixture = autonomous_pyramid_fixture()
    @test length(fixture.prepared.autonomous_optics) == 1
    @test length(fixture.state.controllable_optics) ==
        plant_event_controllable_optic_count(fixture.prepared)
    @test length(fixture.workspace.controllable_optics) ==
        plant_event_controllable_optic_count(fixture.prepared)

    @test step_plant_events!(
        fixture.prepared, fixture.state, fixture.workspace) ==
        PlantTimestamp(0)
    binding = only(fixture.prepared.autonomous_optics)
    optic_state =
        fixture.state.controllable_optics[Int(binding.optic_slot)]
    evaluate_autonomous_periodic_optic!(
        binding.implementation,
        optic_state,
        binding.coupling,
        PlantTimestamp(1),
    )
    optic_state.optical_dirty = true
    phase_bytes = autonomous_pyramid_phase_allocations(
        binding.implementation, optic_state, PlantTimestamp(2))
    evaluation_bytes = autonomous_pyramid_evaluation_allocations(
        binding.implementation,
        optic_state,
        binding.coupling,
        PlantTimestamp(2),
    )
    if coverage_instrumented()
        @test phase_bytes >= 0
        @test evaluation_bytes >= 0
    else
        @test phase_bytes == 0
        @test evaluation_bytes == 0
    end
end

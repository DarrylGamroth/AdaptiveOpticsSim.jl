struct ReducedOrderTestPathModel end
struct ReducedOrderTestOpticModel end

struct ReducedOrderCountedPathExecution{E}
    imaging::E
    executions::Base.RefValue{Int}
end

struct PreparedReducedOrderTestOptic
    endpoint::CommandEndpointID
    command_count::Int
end

mutable struct ReducedOrderTestOpticState{V<:AbstractVector}
    visible::V
end

mutable struct ReducedOrderTestOpticWorkspace{V<:AbstractVector}
    staged::V
end

Plant.plant_model_definition_style(
    ::Type{ReducedOrderTestPathModel}) = ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{ReducedOrderTestOpticModel}) = ColdPlantModelDefinition()

function Plant.validate_path_execution_binding(
    execution::ReducedOrderCountedPathExecution, input, result)
    return Plant.validate_path_execution_binding(
        execution.imaging, input, result)
end

function Plant.execute_path!(result, input,
    execution::ReducedOrderCountedPathExecution)
    execution.executions[] += 1
    return Plant.execute_path!(result, input, execution.imaging)
end

function Plant.prepare_path_executor(
    ::ReducedOrderTestPathModel,
    definition::OpticalPathDefinition,
    source::AbstractSource,
    telescope::Telescope,
    atmosphere::AdaptiveOpticsSim.AbstractTimedAtmosphere)
    T = eltype(pupil_reflectivity(telescope))
    pupil = PupilFunction(telescope; T, backend=backend(telescope))
    imaging = prepare_direct_imaging(pupil, source; zero_padding=1)
    execution = ReducedOrderCountedPathExecution(imaging, Ref(0))
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
        optical_model=:reduced_order_unused_direct_imaging,
        propagation_model=:fraunhofer_fft,
        model_revisions=UInt(1),
    )
end

function Plant.prepare_controllable_optic(
    ::ReducedOrderTestOpticModel,
    definition::ControllableOpticDefinition,
    ::Telescope,
    ::AdaptiveOpticsSim.AbstractAtmosphere)
    schema = only(command_schemas(definition))
    dimensions = command_dimensions(schema)
    length(dimensions) == 1 || throw(PlantPreparationError(
        :controllable_optic, :invalid_dimensions,
        "reduced-order test optic requires a vector command"))
    return PreparedReducedOrderTestOptic(command_endpoint_id(schema),
        only(dimensions))
end

function Plant.prepare_controllable_optic_state(
    prepared::PreparedReducedOrderTestOptic,
    ::ControllableOpticDefinition,
    endpoint_ids::Tuple,
    initial_commands::Tuple)
    only(endpoint_ids) == prepared.endpoint || throw(PlantPreparationError(
        :controllable_optic, :prepared_binding,
        "reduced-order test optic endpoint changed"))
    initial = only(initial_commands)
    length(initial) == prepared.command_count || throw(
        PlantPreparationError(:controllable_optic, :prepared_binding,
            "reduced-order test optic command shape changed"))
    return ReducedOrderTestOpticState(initial)
end

function Plant.prepare_controllable_optic_workspace(
    prepared::PreparedReducedOrderTestOptic)
    staged = Vector{Float64}(undef, prepared.command_count)
    fill!(staged, 0.0)
    return ReducedOrderTestOpticWorkspace(staged)
end

function Plant.stage_controllable_optic_command!(
    prepared::PreparedReducedOrderTestOptic,
    ::ReducedOrderTestOpticState,
    workspace::ReducedOrderTestOpticWorkspace,
    endpoint::CommandEndpointID,
    command::AbstractVector,
    ::PlantTimestamp)
    endpoint == prepared.endpoint || throw(PlantCommandError(
        :physical_application, :endpoint_mismatch,
        "reduced-order test optic received another endpoint"))
    copyto!(workspace.staged, command)
    return nothing
end

function Plant.commit_controllable_optic_command!(
    ::PreparedReducedOrderTestOptic,
    state::ReducedOrderTestOpticState,
    workspace::ReducedOrderTestOpticWorkspace,
    ::CommandEndpointID,
    ::PlantTimestamp)
    copyto!(state.visible, workspace.staged)
    return nothing
end

function reduced_order_test_schema(endpoint::Symbol;
    T::Type{<:AbstractFloat}=Float64,
    dimensions=(2,),
    units=:metre,
    sign_convention=:positive_command_increases_residual,
    basis=CommandBasis(:modal, endpoint),
    basis_revision=1,
    bounds=UniformCommandBounds(T(-20), T(20)))
    return PlantCommandSchema(
        T,
        dimensions;
        id=Symbol(endpoint, :_schema),
        version=1,
        endpoint,
        units,
        sign_convention,
        basis,
        basis_revision,
        semantics=AbsoluteCommand,
        bounds,
        value_policy=CommandValuePolicy(
            range_stage=EnforceOnApplication,
            out_of_range=ClipInvalidCommand),
        sequence_policy=CommandSequencePolicy(),
        effective_time_policy=CommandEffectiveTimePolicy(
            supersession=PreservePendingCommands),
        silence_policy=CommandSilencePolicy(),
    )
end

function reduced_order_test_response(schema::PlantCommandSchema,
    operator;
    units=command_units(schema),
    sign_convention=command_sign_convention(schema),
    basis=command_basis(schema),
    basis_revision=command_basis_revision(schema))
    return Plant.ReducedOrderCommandResponse(
        command_endpoint_id(schema), operator;
        units, sign_convention, basis, basis_revision)
end

function reduced_order_test_fixture(;
    schemas=(reduced_order_test_schema(:dm),),
    responses=nothing,
    disturbance=Plant.HarmonicDisturbanceModel(
        [0.30, -0.22], [19.0, 31.0];
        offsets=[0.05, -0.03], phases_rad=[0.2, -0.4]),
    path_projection=Matrix{Float64}(I, 2, 2),
    sensor_operator=Matrix{Float64}(I, 2, 2),
    calibration_revision=1,
    exposure_ns=1_000_000,
    acquisition_period_ns=2_000_000,
    optical_period_ns=1_000_000,
    command_capacity=128,
    acquisition_start=nothing,
    trigger_topology=nothing)
    response_values = responses === nothing ?
        (reduced_order_test_response(only(schemas),
            Matrix{Float64}(I, 2, 2)),) : responses
    model = Plant.LinearReducedOrderAcquisitionModel(
        disturbance, path_projection, sensor_operator, response_values;
        measurement_kind=:modal_residual,
        measurement_units=:metre,
        residual_kind=:modal_wavefront_error,
        residual_units=:metre,
        calibration_revision,
        operating_envelope=(
            maximum_absolute_residual_m=2.0,
            maximum_disturbance_frequency_hz=40.0,
            sample_period_ns=optical_period_ns,
        ),
        omitted_effects=(
            :diffraction,
            :spatial_aliasing,
            :detector_noise,
            :pyramid_optical_gain,
            :coronagraph_propagation,
        ),
    )
    T = eltype(disturbance.amplitudes)
    telescope = Telescope(resolution=8, diameter=T(8),
        central_obstruction=zero(T), T=T)
    atmosphere = MultiLayerAtmosphere(telescope; r0=T(0.2), L0=T(25),
        fractional_cn2=T[1], wind_speed=T[0],
        wind_direction=T[0], altitude=T[0],
        layer_ids=(:ground,), T=T)
    source = Source(band=:custom, wavelength=T(0.8e-6),
        photon_irradiance=T(1), T=T)
    path = OpticalPathDefinition(:reduced_wfs_path, source,
        ReducedOrderTestPathModel())
    acquisition = AcquisitionDefinition(:reduced_wfs, :reduced_wfs_path,
        model)
    optics = map(schemas) do schema
        endpoint = command_endpoint_id(schema)
        ControllableOpticDefinition(Symbol(endpoint.name, :_optic),
            ReducedOrderTestOpticModel(), (schema,);
            placement=PupilPlanePlacement(),
            visibility=AllPathVisibility())
    end
    configurations = map(schemas) do schema
        T_command = command_numeric_type(schema)
        initial = zeros(T_command, command_dimensions(schema)...)
        CommandEndpointConfiguration(command_endpoint_id(schema), initial;
            capacity=command_capacity)
    end
    definition = PlantDefinition(; telescope, atmosphere,
        controllable_optics=optics, paths=(path,),
        acquisitions=(acquisition,))
    plant = prepare_plant(definition; run_seed=0x7b00,
        command_endpoints=configurations)
    start = acquisition_start === nothing ?
        PeriodicAcquisitionStart(PeriodicSchedule(
            period_ns=acquisition_period_ns, phase_ns=0)) :
        acquisition_start
    loop_definition = PlantEventLoopDefinition(
        (OpticalSampleDefinition(:reduced_wfs_path,
            PeriodicSchedule(period_ns=optical_period_ns, phase_ns=0)),),
        (AcquisitionEventDefinition(:reduced_wfs,
            DirectMeasurementAcquisitionDefinition(
                PlantDuration(exposure_ns)),
            start),);
        trigger_topology,
    )
    prepared = prepare_plant_event_loop(plant, loop_definition)
    return (; plant, prepared, state=PlantEventLoopState(prepared),
        workspace=PlantEventLoopWorkspace(prepared), schemas,
        acquisition_period_ns, exposure_ns)
end

function reduced_order_measurement(fixture)
    measurement = acquisition_measurement(
        prepared_acquisition(fixture.plant, :reduced_wfs))
    return measurement_storage(measurement)
end

function run_to_reduced_order_product!(fixture, sequence::Integer)
    while acquisition_product_sequence(fixture.prepared, fixture.state,
        :reduced_wfs) < sequence
        step_plant_events!(fixture.prepared, fixture.state,
            fixture.workspace)
    end
    return copy(reduced_order_measurement(fixture))
end

function submit_reduced_order_command!(fixture, schema, sequence,
    effective_ns, command, admission_ns)
    clear_command_dispositions!(fixture.workspace)
    return admit_plant_command!(fixture.prepared, fixture.state,
        fixture.workspace,
        PlantCommand(schema, sequence, PlantTimestamp(effective_ns),
            command),
        PlantTimestamp(admission_ns))
end

function reduced_order_closed_loop_trace(; mode=:matched,
    response_gain=1.0, frame_count=48, delay_frames=0,
    stale_frames=0, drop_stride=0)
    schema = reduced_order_test_schema(:dm)
    response = reduced_order_test_response(schema,
        response_gain .* Matrix{Float64}(I, 2, 2))
    fixture = reduced_order_test_fixture(; schemas=(schema,),
        responses=(response,))
    command = zeros(2)
    measurement_history = Vector{Vector{Float64}}()
    metric = Float64[]
    for frame in 1:frame_count
        measurement = run_to_reduced_order_product!(fixture, frame)
        push!(measurement_history, measurement)
        push!(metric, sum(abs2, measurement))
        drop_stride > 0 && frame % drop_stride != 0 && continue
        source_index = max(1, frame - stale_frames)
        feedback = measurement_history[source_index]
        sign = mode === :wrong_sign ? 1.0 : -1.0
        @. command += sign * 0.65 * feedback
        next_start_ns = frame * fixture.acquisition_period_ns
        effective_ns = next_start_ns +
            delay_frames * fixture.acquisition_period_ns
        admission_ns = (frame - 1) * fixture.acquisition_period_ns +
            fixture.exposure_ns + 1
        submit_reduced_order_command!(fixture, schema, frame,
            effective_ns, command, admission_ns)
    end
    return (; fixture, metric, trace=measurement_history)
end

@testset "Direct-measurement acquisition lifecycle" begin
    @test Base.isexported(Plant, :AcquisitionEventDefinition)
    @test Base.isexported(Plant, :DirectMeasurementAcquisitionDefinition)
    @test !isdefined(Plant, :DetectorEventDefinition)
    for name in (
        :AbstractAcquisitionLifecycleDefinition,
        :AbstractPreparedAcquisitionLifecycle,
        :AbstractAcquisitionLifecycleState,
        :HarmonicDisturbanceModel,
        :ReducedOrderCommandResponse,
        :LinearReducedOrderAcquisitionModel,
        :reduced_order_residual,
        :reduced_order_residual_rms,
    )
        @test Base.ispublic(Plant, name)
        @test !Base.isexported(Plant, name)
    end

    measurement = WFSMeasurement(zeros(2);
        units=:metre, kind=:modal_residual)
    definition = DirectMeasurementAcquisitionDefinition(
        PlantDuration(2_000_000_000);
        readout_duration=PlantDuration(5),
        readiness_delay=PlantDuration(7))
    prepared = Plant.prepare_direct_measurement_acquisition(measurement,
        definition)
    state = Plant.DirectMeasurementAcquisitionState(prepared)
    begin_exposure!(prepared, state, PlantTimestamp(10))
    prepared.instantaneous_sample .= [1.0, 3.0]
    Plant.accumulate_direct_measurement_interval!(prepared, state,
        PlantTimestamp(10), PlantTimestamp(1_000_000_010))
    prepared.instantaneous_sample .= [3.0, 5.0]
    Plant.accumulate_direct_measurement_interval!(prepared, state,
        PlantTimestamp(1_000_000_010), PlantTimestamp(2_000_000_010))
    close_exposure!(prepared, state, PlantTimestamp(2_000_000_010))
    output = complete_readout!(prepared, state,
        PlantTimestamp(2_000_000_015), Xoshiro(0x7b01))
    @test output === measurement
    @test measurement_storage(output) == [2.0, 4.0]
    mark_acquisition_ready!(prepared, state,
        PlantTimestamp(2_000_000_022))
    @test Plant.direct_measurement_acquisition_status(state) ==
        Plant.DirectMeasurementReady
end

@testset "Reduced-order provider causality, replay, and path bypass" begin
    first = reduced_order_test_fixture()
    first_product = run_to_reduced_order_product!(first, 1)
    @test first.plant.paths[1].execution.executions[] == 0
    @test first_product != zeros(2)
    @test acquisition_product_ready_timestamp(first.prepared, first.state,
        :reduced_wfs) == PlantTimestamp(1_000_000)
    owner = prepared_acquisition(first.plant, :reduced_wfs)
    @test acquisition_provider_style(owner) isa
        CommandResponsiveReducedOrderProviderStyle
    @test acquisition_provider_payload_work(owner) ===
        :linear_reduced_order_direct_measurement
    @test Plant.reduced_order_sample_timestamp(owner) ==
        PlantTimestamp(1_000_000)
    @test Plant.reduced_order_residual_rms(owner) >= 0
    metadata = acquisition_product_metadata(owner)
    @test metadata.model.calibration_revision == UInt32(1)
    @test metadata.model.residual_metric == :root_mean_square
    @test :diffraction in metadata.model.omitted_effects

    schema = only(first.schemas)
    submit_reduced_order_command!(first, schema, 1, 2_000_000,
        [-0.10, 0.06], 1_000_001)
    second_product = run_to_reduced_order_product!(first, 2)
    @test effective_command(first.prepared, first.state, :dm) ==
        [-0.10, 0.06]
    @test second_product != first_product
    @test first.plant.paths[1].execution.executions[] == 0

    replay_a = reduced_order_test_fixture()
    replay_b = reduced_order_test_fixture()
    trace_a = [run_to_reduced_order_product!(replay_a, index)
        for index in 1:12]
    trace_b = [run_to_reduced_order_product!(replay_b, index)
        for index in 1:12]
    @test trace_a == trace_b
end

@testset "Reduced-order trigger faults and effective command policy" begin
    schema = reduced_order_test_schema(:dm)
    constant_disturbance = Plant.HarmonicDisturbanceModel(
        [0.0, 0.0], [0.0, 0.0]; offsets=[1.0, 1.0])
    clipped = reduced_order_test_fixture(; schemas=(schema,),
        disturbance=constant_disturbance)
    admission = submit_reduced_order_command!(clipped, schema, 1, 0,
        [-30.0, 30.0], 0)
    @test command_admission_status(admission) == CommandAdmittedReady
    @test run_to_reduced_order_product!(clipped, 1) == [-19.0, 21.0]
    @test effective_command(clipped.prepared, clipped.state, :dm) ==
        [-20.0, 20.0]

    trigger_source = TriggerSourceDefinition(:reduced_wfs_trigger,
        PeriodicSchedule(period_ns=2_000_000, phase_ns=0);
        faults=TriggerFaultTrace(
            TriggerFaultTraceEntry(2, :dropped_reduced_wfs_trigger;
                action=DropTriggerEdge),
        ))
    trigger_topology = prepare_trigger_topology((trigger_source,), (),
        (TriggerConsumerDefinition(:reduced_wfs_camera,
            TriggerSourceID(:reduced_wfs_trigger)),);
        in_flight_capacity=2)
    triggered = reduced_order_test_fixture(;
        disturbance=constant_disturbance,
        acquisition_start=TriggeredAcquisitionStart(:reduced_wfs_camera),
        trigger_topology)
    run_plant_events_until!(triggered.prepared, triggered.state,
        triggered.workspace, PlantTimestamp(4_999_999))
    @test acquisition_product_sequence(triggered.prepared, triggered.state,
        :reduced_wfs) == 1
    run_plant_events_until!(triggered.prepared, triggered.state,
        triggered.workspace, PlantTimestamp(5_000_000))
    @test acquisition_product_sequence(triggered.prepared, triggered.state,
        :reduced_wfs) == 2
    @test acquisition_product_ready_timestamp(triggered.prepared,
        triggered.state, :reduced_wfs) == PlantTimestamp(5_000_000)
end

@testset "Independent command timing and endpoint schema validation" begin
    first_schema = reduced_order_test_schema(:woofer)
    second_schema = reduced_order_test_schema(:tweeter)
    first_response = reduced_order_test_response(first_schema,
        [1.0 0.0; 0.0 0.0])
    second_response = reduced_order_test_response(second_schema,
        [0.0 0.0; 0.0 1.0])
    constant_disturbance = Plant.HarmonicDisturbanceModel(
        [0.0, 0.0], [0.0, 0.0]; offsets=[1.0, 1.0])
    fixture = reduced_order_test_fixture(
        schemas=(first_schema, second_schema),
        responses=(first_response, second_response),
        disturbance=constant_disturbance)
    step_plant_events!(fixture.prepared, fixture.state, fixture.workspace)
    @test Plant.reduced_order_residual(
        prepared_acquisition(fixture.plant, :reduced_wfs)) == [1.0, 1.0]
    submit_reduced_order_command!(fixture, first_schema, 1, 1_000_000,
        [-1.0, 0.0], 1)
    clear_command_dispositions!(fixture.workspace)
    admit_plant_command!(fixture.prepared, fixture.state, fixture.workspace,
        PlantCommand(second_schema, 1, PlantTimestamp(2_000_000),
            [0.0, -1.0]), PlantTimestamp(1))
    step_plant_events!(fixture.prepared, fixture.state, fixture.workspace)
    @test Plant.reduced_order_residual(
        prepared_acquisition(fixture.plant, :reduced_wfs)) == [0.0, 1.0]
    step_plant_events!(fixture.prepared, fixture.state, fixture.workspace)
    @test Plant.reduced_order_residual(
        prepared_acquisition(fixture.plant, :reduced_wfs)) == [0.0, 0.0]

    base_schema = reduced_order_test_schema(:dm)
    cases = (
        (:command_dimensions,
            reduced_order_test_response(base_schema, ones(2, 1))),
        (:command_units,
            reduced_order_test_response(base_schema, ones(2, 2);
                units=:radian)),
        (:command_sign_convention,
            reduced_order_test_response(base_schema, ones(2, 2);
                sign_convention=:opposite_sign)),
        (:command_basis,
            reduced_order_test_response(base_schema, ones(2, 2);
                basis=CommandBasis(:modal, :other_basis))),
        (:command_basis_revision,
            reduced_order_test_response(base_schema, ones(2, 2);
                basis_revision=2)),
    )
    for (reason, response) in cases
        error = try
            reduced_order_test_fixture(; schemas=(base_schema,),
                responses=(response,))
            nothing
        catch caught
            caught
        end
        @test error isa PlantPreparationError
        @test error.component == :reduced_order
        @test error.reason == reason
    end

    scalar_schema = reduced_order_test_schema(:focus; dimensions=())
    scalar_endpoint = prepare_command_endpoint(scalar_schema;
        capacity=1, ordinal=1)
    scalar_binding = Plant._PreparedPlantCommandEndpoint(scalar_endpoint,
        UInt32(1), 0.25, nothing)
    scalar_response = reduced_order_test_response(scalar_schema, [2.0, -1.0])
    prepared_scalar_response = Plant._prepare_reduced_order_event_response(
        Plant._prepare_reduced_order_response(scalar_response, CPUBackend()),
        (scalar_binding,))
    scalar_endpoint_state = CommandEndpointState(scalar_endpoint)
    scalar_application = CommandApplicationState(scalar_endpoint,
        scalar_endpoint_state, 0.25)
    scalar_applications = Memory{CommandApplicationState}(undef, 1)
    scalar_applications[1] = scalar_application
    scalar_residual = [1.0, 1.0]
    scalar_workspace = zeros(2)
    @test Plant._apply_reduced_order_response!(scalar_residual,
        scalar_workspace, prepared_scalar_response,
        scalar_applications) === scalar_residual
    @test scalar_residual == [1.5, 0.75]

    matrix_schema = reduced_order_test_schema(:segmented_dm;
        dimensions=(1, 2))
    matrix_endpoint = prepare_command_endpoint(matrix_schema;
        capacity=1, ordinal=1)
    matrix_initial = reshape([0.25, -0.5], 1, 2)
    matrix_binding = Plant._PreparedPlantCommandEndpoint(matrix_endpoint,
        UInt32(1), matrix_initial, nothing)
    matrix_response = reduced_order_test_response(matrix_schema,
        [1.0 2.0; -1.0 0.5])
    prepared_matrix_response = Plant._prepare_reduced_order_event_response(
        Plant._prepare_reduced_order_response(matrix_response, CPUBackend()),
        (matrix_binding,))
    matrix_endpoint_state = CommandEndpointState(matrix_endpoint)
    matrix_application = CommandApplicationState(matrix_endpoint,
        matrix_endpoint_state, matrix_initial)
    matrix_applications = Memory{CommandApplicationState}(undef, 1)
    matrix_applications[1] = matrix_application
    matrix_residual = [1.0, 1.0]
    matrix_workspace = zeros(2)
    @test Plant._apply_reduced_order_response!(matrix_residual,
        matrix_workspace, prepared_matrix_response,
        matrix_applications) === matrix_residual
    @test matrix_residual == [0.25, 0.5]

    @test_throws PlantPreparationError Plant.LinearReducedOrderAcquisitionModel(
        Plant.HarmonicDisturbanceModel([1.0], [1.0]),
        ones(1, 1), ones(1, 1),
        (reduced_order_test_response(base_schema, ones(1, 2)),);
        measurement_kind=:modal_residual,
        measurement_units=:metre,
        residual_kind=:modal_wavefront_error,
        residual_units=:metre,
        calibration_revision=0,
        operating_envelope=(maximum=1.0,),
        omitted_effects=(:diffraction,))
end

@testset "Matched and degraded reduced-order closed loops" begin
    open_loop = reduced_order_closed_loop_trace(; mode=:matched,
        drop_stride=typemax(Int))
    matched = reduced_order_closed_loop_trace()
    wrong_sign = reduced_order_closed_loop_trace(; mode=:wrong_sign)
    delayed = reduced_order_closed_loop_trace(; delay_frames=7)
    stale = reduced_order_closed_loop_trace(; stale_frames=7)
    dropped = reduced_order_closed_loop_trace(; drop_stride=4)
    mismatch = reduced_order_closed_loop_trace(; response_gain=0.25)

    tail = 29:48
    score(result) = sqrt(mean(@view result.metric[tail]))
    open_score = score(open_loop)
    matched_score = score(matched)
    @test all(sample -> all(isfinite, sample), matched.trace)
    @test maximum(sample -> maximum(abs, sample), matched.trace) < 2.0
    @test maximum(abs, effective_command(matched.fixture.prepared,
        matched.fixture.state, :dm)) < 20.0
    @test matched_score < 0.45open_score
    @test score(wrong_sign) > open_score
    @test score(delayed) > 1.5matched_score
    @test score(stale) > 1.5matched_score
    @test score(dropped) > 1.25matched_score
    @test score(mismatch) > 1.25matched_score
end

@testset "Reduced-order event-loop allocation budget" begin
    fixture = reduced_order_test_fixture()
    run_to_reduced_order_product!(fixture, 6)
    acquisition = fixture.prepared.acquisitions[1]
    destination = acquisition.lifecycle.instantaneous_sample
    provider = acquisition.sample_provider
    applications = fixture.state.command_applications
    timestamp = Plant.reduced_order_sample_timestamp(provider.provider)
    Plant.evaluate_linear_reduced_order_sample!(destination, provider,
        timestamp, applications)
    sample_allocated = @allocated Plant.evaluate_linear_reduced_order_sample!(
        destination, provider, timestamp, applications)
    maximum_allocated = 0
    for _ in 1:8
        maximum_allocated = max(maximum_allocated,
            @allocated step_plant_events!(fixture.prepared, fixture.state,
                fixture.workspace))
    end
    if coverage_instrumented()
        @test sample_allocated >= 0
        @test maximum_allocated >= 0
    else
        @test sample_allocated == 0
        @test maximum_allocated <= 2048
    end
end

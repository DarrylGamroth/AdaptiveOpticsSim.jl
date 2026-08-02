struct MixedSerialEventPathModel end
struct MixedSerialEventAcquisitionModel
    exposure_s::Float64
    mode::Symbol
end

Plant.plant_model_definition_style(::Type{MixedSerialEventPathModel}) =
    ColdPlantModelDefinition()
Plant.plant_model_definition_style(
    ::Type{MixedSerialEventAcquisitionModel}) = ColdPlantModelDefinition()

function Plant.prepare_pupil_opd_publication_materialization(
    ::MixedSerialEventPathModel,
    atmosphere::AdaptiveOpticsSim.Atmospheres.AbstractTimedAtmosphere,
    telescope::Telescope,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    pupil::PupilFunction,
)
    return prepare_pupil_opd_materialization(
        atmosphere, telescope, source, pupil)
end

struct MixedSerialEventPathExecution{I,R}
    input::I
    result::R
end

function Plant.validate_path_execution_binding(
    execution::MixedSerialEventPathExecution,
    input::PupilFunction,
    result::IntensityMap,
)
    execution.input === input && execution.result === result || throw(
        PlantPreparationError(
            :path,
            :prepared_binding,
            "mixed serial event test execution lost its exact products",
        ))
    return nothing
end

function Plant.validate_path_execution_target(
    execution::MixedSerialEventPathExecution,
    target::AbstractComputeDevice,
)
    compute_device(execution.input.opd) == target &&
        compute_device(execution.result.values) == target || throw(
            PlantPreparationError(
                :path,
                :wrong_device,
                "mixed serial event test execution occupies another target",
            ))
    return execution
end


function Plant.execute_path!(
    result::IntensityMap,
    input::PupilFunction,
    execution::MixedSerialEventPathExecution,
)
    validate_path_execution_binding(execution, input, result)
    @. result.values = 10.0 + abs(input.opd) * 1.0e8
    return result
end


function Plant.structural_resource_fact(
    execution::MixedSerialEventPathExecution,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    validate_path_execution_target(execution, target)
    resident = structural_array_bytes(execution.input.support, target) +
        structural_array_bytes(execution.input.amplitude, target) +
        structural_array_bytes(execution.input.opd, target) +
        structural_array_bytes(execution.result.values, target)
    return KnownStructuralResourceFact(
        id, target, resident, zero(UInt64))
end

function Plant.prepare_target_local_path_resources(
    ::MixedSerialEventPathModel,
    definition::OpticalPathDefinition,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    telescope::Telescope,
    context,
)
    pupil, result, _ = path_input_publication_test_path_parts(telescope)
    execution = MixedSerialEventPathExecution(pupil, result)
    return PreparedTargetLocalPathResources(
        definition,
        source,
        telescope,
        pupil,
        result,
        execution;
        context,
        optical_model=:mixed_serial_event_test,
        propagation_model=:positive_pupil_opd,
        model_revisions=UInt(1),
    )
end

function Plant.prepare_target_local_acquisition_provider(
    model::MixedSerialEventAcquisitionModel,
    ::AcquisitionDefinition,
    path::PreparedTargetLocalPathResources,
)
    result = path_result(path)
    T = eltype(result.values)
    sensor = if model.mode === :global
        CMOSSensor(timing_model=GlobalShutter(), T=T)
    elseif model.mode === :rolling
        CMOSSensor(
            timing_model=RollingShutter(T(0.002); row_group_size=4),
            T=T,
        )
    elseif model.mode === :frame_transfer
        EMCCDSensor(
            excess_noise_factor=one(T),
            acquisition_mode=FrameTransferAcquisition(
                transfer_time=T(0.004)),
            T=T,
        )
    elseif model.mode === :up_the_ramp
        HgCdTeSensor(
            read_time=zero(T),
            sampling_mode=UpTheRampSampling(3),
            T=T,
        )
    else
        error("unsupported mixed serial event detector mode $(model.mode)")
    end
    detector = Detector(
        integration_time=T(model.exposure_s),
        noise=NoiseNone(),
        qe=one(T),
        gain=one(T),
        response_model=NullFrameResponse(),
        sensor=sensor,
        T=T,
        backend=path_result_key(path).backend,
    )
    execution = FrameAcquisitionExecution(detector, result)
    products = AcquisitionProducts(
        execution.observation;
        metadata=(
            kind=:mixed_serial_event_frame,
            units=:detected_electrons,
            geometry=result.metadata,
            detector=detector_export_metadata(detector),
        ),
    )
    return prepare_full_optical_provider(execution, products)
end

function mixed_serial_event_lifecycle(mode::Symbol)
    exposure = PlantDuration(80_000_000)
    if mode === :rolling
        return RollingShutterAcquisitionDefinition(
            exposure; readiness_delay=PlantDuration(5_000_000))
    elseif mode === :frame_transfer
        return FrameTransferAcquisitionDefinition(
            exposure; readout_duration=PlantDuration(10_000_000))
    end
    return GlobalShutterAcquisitionDefinition(
        exposure;
        readout_duration=PlantDuration(10_000_000),
        readiness_delay=PlantDuration(5_000_000),
    )
end

function mixed_serial_event_fixture(;
    beta_target::AbstractComputeDevice=HANDOFF_TEST_ACCELERATOR,
    with_command::Bool=false,
    command_authority_target::AbstractComputeDevice=HostComputeDevice(),
    silence_policy::CommandSilencePolicy=CommandSilencePolicy(),
    safe_command=nothing,
    detector_mode::Symbol=:global,
    triggered::Bool=false,
    run_seed::Integer=0x9a3c,
)
    camera = AcquisitionDefinition(
        :alpha_camera,
        :alpha,
        MixedSerialEventAcquisitionModel(0.08, detector_mode),
    )
    schema = with_command ?
        effective_command_route_test_schema(
            dimensions=(), silence_policy=silence_policy) : nothing
    optic = with_command ? ControllableOpticDefinition(
        :effective_command_route_optic,
        EffectiveCommandRouteTestOpticModel(),
        (schema,);
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
    ) : nothing
    definition = path_input_publication_test_definition(
        alpha_model=MixedSerialEventPathModel(),
        controllable_optics=with_command ? (optic,) : (),
        acquisitions=(camera,),
    )
    assignment = resolve_plant_partition_assignment(
        definition,
        HostComputeDevice(),
        :alpha => HostComputeDevice(),
        :beta => beta_target,
    )
    command_endpoints = if with_command
        (CommandEndpointConfiguration(
            command_endpoint_id(schema),
            effective_command_route_test_initial(schema);
            capacity=4,
            safe_command,
        ),)
    else
        ()
    end
    partitions = with_path_input_publication_cold_scalar_indexing() do
        prepare_plant_partitions(
            definition,
            assignment;
            run_seed,
            command_authority_target,
            command_endpoints,
        )
    end
    trigger_topology = if triggered
        prepare_trigger_topology(
            (
                TriggerSourceDefinition(
                    :alpha_camera_source,
                    PeriodicSchedule(
                        period_ns=200_000_000,
                        phase_ns=10_000_000,
                    ),
                ),
            ),
            (),
            (
                TriggerConsumerDefinition(
                    :alpha_camera_trigger,
                    TriggerSourceID(:alpha_camera_source),
                ),
            );
            in_flight_capacity=1,
        )
    else
        nothing
    end
    start = triggered ?
        TriggeredAcquisitionStart(:alpha_camera_trigger) :
        PeriodicAcquisitionStart(
            PeriodicSchedule(
                period_ns=200_000_000,
                phase_ns=10_000_000,
            ),
        )
    events = PlantEventLoopDefinition(
        (
            OpticalSampleDefinition(
                :alpha, PeriodicSchedule(period_ns=100_000_000)),
            OpticalSampleDefinition(
                :beta, PeriodicSchedule(period_ns=150_000_000)),
        ),
        (
            AcquisitionEventDefinition(
                :alpha_camera,
                mixed_serial_event_lifecycle(detector_mode),
                start,
            ),
        ),
        trigger_topology=trigger_topology,
    )
    prepared = prepare_plant_event_loop(partitions, events)
    state = Plant.MixedResourcePlantEventLoopState(prepared)
    workspace = Plant.MixedResourcePlantEventLoopWorkspace(
        prepared, state)
    return (; partitions, events, prepared, state, workspace, schema)
end

function mixed_serial_transaction_event_fixture(;
    beta_target::AbstractComputeDevice=HANDOFF_TEST_ACCELERATOR,
)
    alpha_schema = effective_command_route_test_schema(
        dimensions=(),
        id=:mixed_serial_transaction_alpha_schema,
        endpoint=:mixed_serial_transaction_alpha_endpoint,
        basis=:mixed_serial_transaction_alpha_basis,
    )
    beta_schema = effective_command_route_test_schema(
        dimensions=(),
        id=:mixed_serial_transaction_beta_schema,
        endpoint=:mixed_serial_transaction_beta_endpoint,
        basis=:mixed_serial_transaction_beta_basis,
    )
    alpha_optic = ControllableOpticDefinition(
        :mixed_serial_transaction_alpha_optic,
        EffectiveCommandRouteTestOpticModel(),
        (alpha_schema,);
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
    )
    beta_optic = ControllableOpticDefinition(
        :mixed_serial_transaction_beta_optic,
        EffectiveCommandRouteTestOpticModel(),
        (beta_schema,);
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
    )
    camera = AcquisitionDefinition(
        :mixed_serial_transaction_camera,
        :alpha,
        MixedSerialEventAcquisitionModel(0.08, :global),
    )
    definition = path_input_publication_test_definition(
        alpha_model=MixedSerialEventPathModel(),
        controllable_optics=(alpha_optic, beta_optic),
        acquisitions=(camera,),
    )
    assignment = resolve_plant_partition_assignment(
        definition,
        HostComputeDevice(),
        :alpha => HostComputeDevice(),
        :beta => beta_target,
    )
    configurations = (
        CommandEndpointConfiguration(
            command_endpoint_id(alpha_schema), 1.0; capacity=4),
        CommandEndpointConfiguration(
            command_endpoint_id(beta_schema), 1.0; capacity=4),
    )
    partitions = with_path_input_publication_cold_scalar_indexing() do
        prepare_plant_partitions(
            definition,
            assignment;
            run_seed=0x9a3c,
            command_authority_target=HostComputeDevice(),
            command_endpoints=configurations,
        )
    end
    events = PlantEventLoopDefinition(
        (
            OpticalSampleDefinition(
                :alpha, PeriodicSchedule(period_ns=100_000_000)),
            OpticalSampleDefinition(
                :beta, PeriodicSchedule(period_ns=150_000_000)),
        ),
        (
            AcquisitionEventDefinition(
                :mixed_serial_transaction_camera,
                mixed_serial_event_lifecycle(:global),
                PeriodicAcquisitionStart(PeriodicSchedule(
                    period_ns=200_000_000,
                    phase_ns=10_000_000,
                )),
            ),
        ),
    )
    prepared = prepare_plant_event_loop(partitions, events)
    state = Plant.MixedResourcePlantEventLoopState(prepared)
    workspace = Plant.MixedResourcePlantEventLoopWorkspace(prepared, state)
    return (; partitions, prepared, state, workspace, alpha_schema, beta_schema)
end

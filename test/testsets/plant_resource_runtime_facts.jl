@inline resource_runtime_memory_bytes(memory::Memory) =
    UInt64(length(memory) * Base.elsize(typeof(memory)))

function resource_runtime_command_schema(::Type{T}=Float64) where {
    T<:AbstractFloat}
    return PlantCommandSchema(
        T,
        (4,);
        id=:resource_runtime_dm_schema,
        version=1,
        endpoint=:resource_runtime_dm_command,
        units=:metre,
        sign_convention=:positive_surface_increases_opd,
        basis=CommandBasis(:actuator, :resource_runtime_dm),
        basis_revision=1,
        semantics=AbsoluteCommand,
        bounds=UnboundedCommandValues(),
        value_policy=CommandValuePolicy(),
        sequence_policy=CommandSequencePolicy(),
        effective_time_policy=CommandEffectiveTimePolicy(),
        silence_policy=CommandSilencePolicy(),
    )
end

@testset "Command structural resource owners" begin
    T = Float64
    target = HostComputeDevice()
    capacity = 3
    dimensions = (4,)
    payload_bytes = UInt64(prod(dimensions) * sizeof(T))
    schema = resource_runtime_command_schema(T)
    endpoint = prepare_command_endpoint(schema;
        capacity, sequence_window=capacity, ordinal=1)
    endpoint_state = CommandEndpointState(endpoint)
    application_state = CommandApplicationState(
        endpoint, endpoint_state, zeros(T, dimensions))
    disposition_workspace = CommandDispositionWorkspace(endpoint)

    endpoint_id = StructuralResourceOwnerID(
        :command_endpoint_state, :resource_runtime_dm_command)
    endpoint_fact = structural_resource_fact(
        endpoint_state, endpoint_id, target)
    expected_endpoint_resident =
        UInt64(capacity) * payload_bytes +
        UInt64(capacity * Base.elsize(typeof(endpoint_state.payloads.values))) +
        resource_runtime_memory_bytes(endpoint_state.slots) +
        resource_runtime_memory_bytes(endpoint_state.calendar) +
        resource_runtime_memory_bytes(endpoint_state.accepted_sequences)
    expected_endpoint_workspace = payload_bytes + UInt64(
        Base.elsize(typeof(endpoint_state.payloads.values)))
    @test structural_resource_known(endpoint_fact)
    @test structural_resident_bytes(endpoint_fact) ==
        expected_endpoint_resident
    @test structural_workspace_bytes(endpoint_fact) ==
        expected_endpoint_workspace

    application_fact = structural_resource_fact(
        application_state,
        StructuralResourceOwnerID(
            :command_application_state, :resource_runtime_dm_command),
        target,
    )
    @test structural_resident_bytes(application_fact) == payload_bytes
    @test structural_workspace_bytes(application_fact) == payload_bytes

    disposition_fact = structural_resource_fact(
        disposition_workspace,
        StructuralResourceOwnerID(
            :command_disposition_workspace, :resource_runtime_dm_command),
        target,
    )
    @test structural_resident_bytes(disposition_fact) == 0
    @test structural_workspace_bytes(disposition_fact) == UInt64(
        capacity * Base.elsize(typeof(disposition_workspace.dispositions)))

    initial = zeros(T, dimensions)
    binding = Plant._PreparedPlantCommandEndpoint(
        endpoint, UInt32(1), initial, nothing)
    binding_fact = structural_resource_fact(
        binding,
        StructuralResourceOwnerID(
            :command_endpoint_plan, :resource_runtime_dm_command),
        target,
    )
    @test structural_resident_bytes(binding_fact) == payload_bytes
    @test structural_workspace_bytes(binding_fact) == 0

    wrong_target = AcceleratorComputeDevice(
        ResourceFactFakeBackend(), UInt32(1))
    @test !structural_resource_known(structural_resource_fact(
        endpoint_state, endpoint_id, wrong_target))
    @test !structural_resource_known(structural_resource_fact(
        application_state,
        StructuralResourceOwnerID(
            :command_application_state, :resource_runtime_dm_command),
        wrong_target,
    ))
end

@testset "Native deformable-mirror structural resource owners" begin
    T = Float64
    resolution = 4
    n_act = 2
    target = HostComputeDevice()
    telescope = Telescope(
        resolution=resolution,
        diameter=T(4),
        central_obstruction=zero(T),
        T=T,
    )
    atmosphere = KolmogorovAtmosphere(
        telescope; r0=T(0.2), L0=T(25), T=T)
    schema = resource_runtime_command_schema(T)
    model = DeformableMirrorModel(
        n_act=n_act, influence_width=T(0.3), T=T)
    definition = ControllableOpticDefinition(
        :resource_runtime_dm,
        model,
        (schema,);
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
    )
    prepared = prepare_controllable_optic(
        model, definition, telescope, atmosphere)
    endpoint = command_endpoint_id(schema)
    initial = (zeros(T, n_act^2),)
    state = prepare_controllable_optic_state(
        prepared, definition, (endpoint,), initial)
    workspace = prepare_controllable_optic_workspace(prepared)

    prepared_fact = structural_resource_fact(
        prepared,
        StructuralResourceOwnerID(
            :controllable_optic_plan, :resource_runtime_dm),
        target,
    )
    topology_bytes = UInt64(
        2 * (2 * n_act^2 * sizeof(T)) +
        n_act^2 * sizeof(Bool) +
        n_act^2 * sizeof(Int))
    operator_bytes = UInt64(
        resolution^2 * sizeof(Bool) +
        2 * n_act^2 * sizeof(T))
    separable_bytes = UInt64(
        2 * resolution * n_act * sizeof(T))
    @test structural_resident_bytes(prepared_fact) ==
        topology_bytes + operator_bytes + separable_bytes
    @test structural_workspace_bytes(prepared_fact) == 0

    runtime_bytes = UInt64(
        resolution^2 * sizeof(T) +
        2 * n_act^2 * sizeof(T) +
        resolution * n_act * sizeof(T))
    state_fact = structural_resource_fact(
        state,
        StructuralResourceOwnerID(
            :controllable_optic_state, :resource_runtime_dm),
        target,
    )
    workspace_fact = structural_resource_fact(
        workspace,
        StructuralResourceOwnerID(
            :controllable_optic_workspace, :resource_runtime_dm),
        target,
    )
    @test structural_resident_bytes(state_fact) == runtime_bytes
    @test structural_workspace_bytes(state_fact) == 0
    @test structural_resident_bytes(workspace_fact) == 0
    @test structural_workspace_bytes(workspace_fact) == runtime_bytes

    # These are borrowed aliases of plan storage and must not reappear in the
    # state/workspace totals.
    @test state.active.state.modes === prepared.modes
    @test workspace.staged.state.modes === prepared.modes
    @test Base.mightalias(
        state.active.state.opd_vec, state.active.state.opd)
    @test Base.mightalias(
        state.active.state.coefs_grid, state.active.state.actuator_coefs)

    wrong_target = AcceleratorComputeDevice(
        ResourceFactFakeBackend(), UInt32(1))
    @test !structural_resource_known(structural_resource_fact(
        prepared,
        StructuralResourceOwnerID(
            :controllable_optic_plan, :resource_runtime_dm),
        wrong_target,
    ))
    @test !structural_resource_known(structural_resource_fact(
        state,
        StructuralResourceOwnerID(
            :controllable_optic_state, :resource_runtime_dm),
        wrong_target,
    ))

    sampled_modes = zeros(T, resolution^2, n_act^2)
    measured_model = MeasuredInfluenceFunctions(
        sampled_modes; metadata=(calibration=copy(sampled_modes),))
    measured_dm_model = DeformableMirrorModel(
        topology=ActuatorGridTopology(n_act; T),
        influence_model=measured_model,
        T=T,
    )
    measured_definition = ControllableOpticDefinition(
        :resource_runtime_measured_dm,
        measured_dm_model,
        (schema,);
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
    )
    measured_prepared = prepare_controllable_optic(
        measured_dm_model,
        measured_definition,
        telescope,
        atmosphere,
    )
    measured_fact = structural_resource_fact(
        measured_prepared,
        StructuralResourceOwnerID(
            :controllable_optic_plan, :resource_runtime_measured_dm),
        target,
    )
    @test !structural_resource_known(measured_fact)
    @test structural_resource_unknown_reason(measured_fact) ==
        :unsupported_influence_storage
end

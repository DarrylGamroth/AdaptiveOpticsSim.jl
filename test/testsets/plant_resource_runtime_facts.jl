@inline resource_runtime_memory_bytes(memory::Memory) =
    UInt64(length(memory) * Base.elsize(typeof(memory)))
@inline resource_runtime_array_bytes(array::AbstractArray) =
    UInt64(length(array) * sizeof(eltype(array)))
@inline resource_runtime_sum_bytes(arrays) =
    sum(resource_runtime_array_bytes, arrays; init=UInt64(0))

struct ResourceRuntimeFakeBackend <: AbstractArrayBackend end
struct ResourceRuntimeUnsupportedInfluence <: AbstractDMInfluenceModel end

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

function resource_runtime_scalar_command_schema(::Type{T}=Float64) where {
    T<:AbstractFloat}
    return PlantCommandSchema(
        T,
        ();
        id=:resource_runtime_scalar_schema,
        version=1,
        endpoint=:resource_runtime_scalar_command,
        units=:metre,
        sign_convention=:positive_surface_increases_opd,
        basis=CommandBasis(:focus, :resource_runtime_scalar),
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
        ResourceRuntimeFakeBackend(), UInt32(1))
    @test !structural_resource_known(structural_resource_fact(
        endpoint_state, endpoint_id, wrong_target))
    @test !structural_resource_known(structural_resource_fact(
        application_state,
        StructuralResourceOwnerID(
            :command_application_state, :resource_runtime_dm_command),
        wrong_target,
    ))
end

@testset "Scalar command structural resource owners" begin
    target = HostComputeDevice()
    wrong_target = AcceleratorComputeDevice(
        ResourceRuntimeFakeBackend(), UInt32(2))
    schema = resource_runtime_scalar_command_schema()
    endpoint = prepare_command_endpoint(
        schema; capacity=2, sequence_window=2, ordinal=1)
    state = CommandEndpointState(endpoint)
    application = CommandApplicationState(endpoint, state, 0.0)
    workspace = CommandDispositionWorkspace(endpoint)
    id = StructuralResourceOwnerID(
        :command_endpoint, :resource_runtime_scalar_command)

    state_fact = structural_resource_fact(state, id, target)
    @test structural_resource_known(state_fact)
    @test structural_resident_bytes(state_fact) > UInt64(0)
    @test structural_workspace_bytes(state_fact) > UInt64(0)

    application_fact = structural_resource_fact(application, id, target)
    @test structural_resource_known(application_fact)
    @test iszero(structural_resident_bytes(application_fact))
    @test iszero(structural_workspace_bytes(application_fact))

    binding = Plant._PreparedPlantCommandEndpoint(
        endpoint, UInt32(1), 0.0, nothing)
    binding_fact = structural_resource_fact(binding, id, target)
    @test structural_resource_known(binding_fact)
    @test iszero(structural_resident_bytes(binding_fact))

    @test !structural_resource_known(
        structural_resource_fact(state, id, wrong_target))
    @test !structural_resource_known(
        structural_resource_fact(application, id, wrong_target))
    @test !structural_resource_known(
        structural_resource_fact(workspace, id, wrong_target))
end

@testset "Deformable-mirror resource dispatch branches" begin
    target = HostComputeDevice()
    wrong_target = AcceleratorComputeDevice(
        ResourceRuntimeFakeBackend(), UInt32(3))

    @test Plant._dm_topology_metadata_is_exact((;))
    @test !Plant._dm_topology_metadata_is_exact((calibration=:external,))
    @test Plant._dm_influence_storage_is_exact(
        GaussianInfluenceWidth(0.3), nothing)

    dense = DenseInfluenceMatrix(zeros(4, 2))
    @test Plant._dm_influence_storage_is_exact(dense, dense.modes)
    @test !Plant._dm_influence_storage_is_exact(dense, copy(dense.modes))
    @test !Plant._dm_influence_storage_is_exact(
        ResourceRuntimeUnsupportedInfluence(), dense.modes)

    sampled = SampledActuatorTopology(
        [-0.5 0.5; 0.0 0.0])
    sampled_bytes = Plant._dm_topology_host_bytes(sampled, target)
    @test sampled_bytes == resource_runtime_sum_bytes((
        sampled.coords,
        sampled.active_coords,
        sampled.valid_actuators,
        sampled.active_indices,
    ))
    sampled_with_metadata = SampledActuatorTopology(
        [-0.5 0.5; 0.0 0.0]; metadata=(calibration=:external,))
    @test isnothing(Plant._dm_topology_host_bytes(
        sampled_with_metadata, target))
    @test iszero(Plant._dm_topology_host_bytes(sampled, wrong_target))

    modes = zeros(Float64, 4, 2)
    @test Plant._dm_backend_mode_bytes(modes, target) ==
        resource_runtime_array_bytes(modes)
    @test iszero(Plant._dm_backend_mode_bytes(modes, wrong_target))
    @test iszero(Plant._dm_host_mode_bytes(modes, target))

    clipped = ClippedActuators(-1.0, 1.0)
    health = ActuatorHealthMap([1.0, 0.5])
    composite = CompositeDMActuatorModel(clipped, health)
    @test iszero(Plant._dm_actuator_model_bytes(clipped, target))
    @test Plant._dm_actuator_model_bytes(health, target) ==
        resource_runtime_array_bytes(health.gains)
    @test iszero(Plant._dm_actuator_model_bytes(health, wrong_target))
    @test Plant._dm_actuator_model_bytes(composite, target) ==
        resource_runtime_array_bytes(health.gains)
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
        ResourceRuntimeFakeBackend(), UInt32(1))
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

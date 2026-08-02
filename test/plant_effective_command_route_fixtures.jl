# Focused command-fanout fixture built on the explicit fake transfer backend.

const EFFECTIVE_COMMAND_ROUTE_TEST_THROW_STAGE_AFTER_MUTATION = Ref(false)

struct EffectiveCommandRouteTestOpticModel end

struct EffectiveCommandRouteTestPreparedOptic{D,S}
    target::D
    schema::S
end

mutable struct EffectiveCommandRouteTestOpticState{V}
    active::V
    timestamp::PlantTimestamp
end

mutable struct EffectiveCommandRouteTestOpticWorkspace{V}
    staging::V
    timestamp::PlantTimestamp
end

Plant.plant_model_definition_style(
    ::Type{EffectiveCommandRouteTestOpticModel}) =
    ColdPlantModelDefinition()

function Plant.prepare_target_local_controllable_optic(
    ::EffectiveCommandRouteTestOpticModel,
    definition::ControllableOpticDefinition,
    ::Telescope,
    ::MultiLayerAtmosphereDefinition,
    target::AbstractComputeDevice,
)
    return EffectiveCommandRouteTestPreparedOptic(
        target, only(command_schemas(definition)))
end

@inline effective_command_route_test_copy(value::Real) = value
@inline effective_command_route_test_copy(value::AbstractArray) = copy(value)

@inline effective_command_route_test_host_value(value::Array) = value
@inline effective_command_route_test_host_value(value::HandoffTestArray) =
    value.storage
@inline effective_command_route_test_host_value(value::Real) = value

@inline effective_command_route_test_active_device(::HostComputeDevice) =
    zero(UInt32)
@inline effective_command_route_test_active_device(
    target::AcceleratorComputeDevice,
) = UInt32(compute_device_identifier(target))

function Plant.prepare_controllable_optic_state(
    prepared::EffectiveCommandRouteTestPreparedOptic,
    ::ControllableOpticDefinition,
    endpoint_ids::Tuple,
    initial_commands::Tuple,
)
    only(endpoint_ids) == command_endpoint_id(prepared.schema) || throw(
        PlantPreparationError(
            :controllable_optic,
            :endpoint_mismatch,
            "command-route test optic endpoint changed",
        ))
    return EffectiveCommandRouteTestOpticState(
        effective_command_route_test_copy(only(initial_commands)),
        zero(PlantTimestamp),
    )
end

@inline function effective_command_route_test_workspace_value(
    schema::PlantCommandSchema{T,0},
    ::AbstractComputeDevice,
) where {T}
    return zero(T)
end

function effective_command_route_test_workspace_value(
    schema::PlantCommandSchema{T,N},
    target::AbstractComputeDevice,
) where {T,N}
    value = allocate_device_array(target, T, command_dimensions(schema)...)
    fill!(value, zero(T))
    return value
end

function Plant.prepare_controllable_optic_workspace(
    prepared::EffectiveCommandRouteTestPreparedOptic,
)
    return EffectiveCommandRouteTestOpticWorkspace(
        effective_command_route_test_workspace_value(
            prepared.schema, prepared.target),
        zero(PlantTimestamp),
    )
end

@inline function effective_command_route_test_stage!(
    workspace::EffectiveCommandRouteTestOpticWorkspace{<:Real},
    value::Real,
)
    workspace.staging = value
    return nothing
end

@inline function effective_command_route_test_stage!(
    workspace::EffectiveCommandRouteTestOpticWorkspace{<:AbstractArray},
    value::AbstractArray,
)
    copyto!(workspace.staging, value)
    return nothing
end

function Plant.stage_controllable_optic_command!(
    prepared::EffectiveCommandRouteTestPreparedOptic,
    ::EffectiveCommandRouteTestOpticState,
    workspace::EffectiveCommandRouteTestOpticWorkspace,
    endpoint::CommandEndpointID,
    value,
    timestamp::PlantTimestamp,
)
    HANDOFF_TEST_ACTIVE_DEVICE[] ==
        effective_command_route_test_active_device(prepared.target) || error(
        "command-route test optic staged outside its exact target context")
    endpoint == command_endpoint_id(prepared.schema) || throw(
        PlantCommandError(
            :physical_application,
            :foreign_endpoint,
            "command-route test optic received another endpoint",
        ))
    effective_command_route_test_stage!(workspace, value)
    workspace.timestamp = timestamp
    EFFECTIVE_COMMAND_ROUTE_TEST_THROW_STAGE_AFTER_MUTATION[] && error(
        "injected command-route physical-stage failure after mutation")
    return nothing
end

function Plant.commit_controllable_optic_command!(
    prepared::EffectiveCommandRouteTestPreparedOptic,
    state::EffectiveCommandRouteTestOpticState,
    workspace::EffectiveCommandRouteTestOpticWorkspace,
    ::CommandEndpointID,
    timestamp::PlantTimestamp,
)
    HANDOFF_TEST_ACTIVE_DEVICE[] ==
        effective_command_route_test_active_device(prepared.target) || error(
        "command-route test optic committed outside its exact target context")
    workspace.timestamp == timestamp || error(
        "command-route test optic timestamp changed before commit")
    state.active, workspace.staging = workspace.staging, state.active
    state.timestamp = timestamp
    return nothing
end

function Plant.validate_controllable_optic_target(
    prepared::EffectiveCommandRouteTestPreparedOptic,
    target::AbstractComputeDevice,
)
    prepared.target == target || throw(PlantPreparationError(
        :controllable_optic,
        :wrong_target,
        "command-route test optic preparation occupies another target",
    ))
    return prepared
end

@inline function effective_command_route_test_value_target(
    ::Real,
    ::AbstractComputeDevice,
)
    return true
end

@inline function effective_command_route_test_value_target(
    value::AbstractArray,
    target::AbstractComputeDevice,
)
    return compute_device(value) == target
end

function Plant.validate_controllable_optic_state_target(
    ::EffectiveCommandRouteTestPreparedOptic,
    state::EffectiveCommandRouteTestOpticState,
    target::AbstractComputeDevice,
)
    effective_command_route_test_value_target(state.active, target) || throw(
        PlantPreparationError(
            :controllable_optic,
            :wrong_state_target,
            "command-route test optic state occupies another target",
        ))
    return state
end

function Plant.validate_controllable_optic_workspace_target(
    ::EffectiveCommandRouteTestPreparedOptic,
    workspace::EffectiveCommandRouteTestOpticWorkspace,
    target::AbstractComputeDevice,
)
    effective_command_route_test_value_target(
        workspace.staging, target) || throw(PlantPreparationError(
        :controllable_optic,
        :wrong_workspace_target,
        "command-route test optic workspace occupies another target",
    ))
    return workspace
end

@inline function effective_command_route_test_value_bytes(
    ::Real,
    ::AbstractComputeDevice,
)
    return UInt64(0)
end

function effective_command_route_test_value_bytes(
    value::AbstractArray,
    target::AbstractComputeDevice,
)
    compute_device(value) == target || return nothing
    return structural_array_bytes(value, target)
end

function Plant.structural_resource_fact(
    prepared::EffectiveCommandRouteTestPreparedOptic,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    prepared.target == target ||
        return UnknownStructuralResourceFact(id, target, :owner_not_on_device)
    return KnownStructuralResourceFact(id, target, UInt64(0), UInt64(0))
end

function Plant.structural_resource_fact(
    state::EffectiveCommandRouteTestOpticState,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    bytes = effective_command_route_test_value_bytes(state.active, target)
    isnothing(bytes) &&
        return UnknownStructuralResourceFact(id, target, :owner_not_on_device)
    return KnownStructuralResourceFact(id, target, bytes, UInt64(0))
end

function Plant.structural_resource_fact(
    workspace::EffectiveCommandRouteTestOpticWorkspace,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    bytes = effective_command_route_test_value_bytes(
        workspace.staging, target)
    isnothing(bytes) &&
        return UnknownStructuralResourceFact(id, target, :owner_not_on_device)
    return KnownStructuralResourceFact(id, target, UInt64(0), bytes)
end

function effective_command_route_test_schema(;
    dimensions=(2,),
    semantics::CommandValueSemantics=IncrementalCommand,
)
    return PlantCommandSchema(
        Float64,
        dimensions;
        id=:effective_command_route_schema,
        version=1,
        endpoint=:effective_command_route_endpoint,
        units=:metre,
        sign_convention=:positive_surface_increases_opd,
        basis=CommandBasis(:actuator, :effective_command_route_basis),
        basis_revision=1,
        semantics,
        bounds=UnboundedCommandValues(),
        value_policy=CommandValuePolicy(),
        sequence_policy=CommandSequencePolicy(),
        effective_time_policy=CommandEffectiveTimePolicy(),
        silence_policy=CommandSilencePolicy(),
    )
end

function effective_command_route_test_definition(;
    dimensions=(2,),
    semantics::CommandValueSemantics=IncrementalCommand,
)
    schema = effective_command_route_test_schema(; dimensions, semantics)
    optic = ControllableOpticDefinition(
        :effective_command_route_optic,
        EffectiveCommandRouteTestOpticModel(),
        (schema,);
        placement=PupilPlanePlacement(),
        visibility=AllPathVisibility(),
    )
    return path_input_publication_test_definition(
        controllable_optics=(optic,))
end

@inline effective_command_route_test_initial(
    ::PlantCommandSchema{T,0}) where {T} = one(T)
@inline effective_command_route_test_initial(
    schema::PlantCommandSchema{T,N}) where {T,N} =
    fill(one(T), command_dimensions(schema))

function effective_command_route_test_partitions(;
    alpha_target::AbstractComputeDevice=HostComputeDevice(),
    beta_target::AbstractComputeDevice=HANDOFF_TEST_ACCELERATOR,
    atmosphere_target::AbstractComputeDevice=HostComputeDevice(),
    command_authority_target::AbstractComputeDevice=HostComputeDevice(),
    dimensions=(2,),
    semantics::CommandValueSemantics=IncrementalCommand,
)
    definition = effective_command_route_test_definition(;
        dimensions, semantics)
    schema = only(command_schemas(
        only(controllable_optic_definitions(definition))))
    configurations = (
        CommandEndpointConfiguration(
            command_endpoint_id(schema),
            effective_command_route_test_initial(schema);
            capacity=2,
        ),)
    assignment = resolve_plant_partition_assignment(
        definition,
        atmosphere_target,
        :alpha => alpha_target,
        :beta => beta_target,
    )
    partitions = with_path_input_publication_cold_scalar_indexing() do
        prepare_plant_partitions(
            definition,
            assignment;
            run_seed=0x219,
            command_authority_target,
            command_endpoints=configurations,
        )
    end
    return partitions, schema
end

function effective_command_route_test_candidate!(
    authority::PreparedCommandAuthority,
    state::CommandAuthorityState,
    payload,
)
    endpoint_id = CommandEndpointID(:effective_command_route_endpoint)
    application = command_authority_application_state(
        authority, state, endpoint_id)
    staged = Plant._stage_application_candidate!(
        getfield(application, :values),
        command_schema(prepared_command_authority_endpoint(
            authority, endpoint_id)),
        payload,
    )
    staged.decision == Plant._AcceptCommandCandidate || error(
        "command-route test candidate was rejected")
    return Plant._staged_effective_command(application)
end

@inline effective_command_route_test_prepare!(
    ::Tuple{},
    ::Tuple{},
    value,
    ::PlantTimestamp,
    ::UInt64,
) = nothing

@inline function effective_command_route_test_prepare!(
    routes::Tuple,
    states::Tuple,
    value,
    timestamp::PlantTimestamp,
    sequence::UInt64,
)
    route = first(routes)
    publication = Plant._effective_command_route_publication(
        route, timestamp, sequence)
    Plant._prepare_effective_command_publication_route!(
        route, first(states), publication, value)
    return effective_command_route_test_prepare!(
        Base.tail(routes), Base.tail(states), value, timestamp, sequence)
end

@inline effective_command_route_test_stage!(
    ::Tuple{},
    ::Tuple{},
    value,
    ::PlantTimestamp,
    ::UInt64,
) = nothing

@inline function effective_command_route_test_stage!(
    routes::Tuple,
    states::Tuple,
    value,
    timestamp::PlantTimestamp,
    sequence::UInt64,
)
    route = first(routes)
    publication = Plant._effective_command_route_publication(
        route, timestamp, sequence)
    Plant._stage_effective_command_publication_route!(
        route, first(states), publication, value)
    return effective_command_route_test_stage!(
        Base.tail(routes), Base.tail(states), value, timestamp, sequence)
end

@inline effective_command_route_test_reclaim!(::Tuple{}, ::Tuple{}) = nothing

@inline function effective_command_route_test_reclaim!(
    routes::Tuple,
    states::Tuple,
)
    Plant._reclaim_effective_command_publication_route!(
        first(routes), first(states))
    return effective_command_route_test_reclaim!(
        Base.tail(routes), Base.tail(states))
end

@inline effective_command_route_test_commit!(::Tuple{}, ::Tuple{}) = nothing

@inline function effective_command_route_test_commit!(
    routes::Tuple,
    states::Tuple,
)
    Plant._commit_effective_command_publication_route!(
        first(routes), first(states))
    return effective_command_route_test_commit!(
        Base.tail(routes), Base.tail(states))
end

function effective_command_route_test_publish!(
    routes::Plant._PreparedEffectiveCommandPublicationRoutes,
    states::Plant._EffectiveCommandPublicationRoutesState,
    value,
    timestamp::PlantTimestamp,
    sequence::UInt64,
)
    storage = routes.storage
    state_storage = states.storage
    effective_command_route_test_prepare!(
        storage, state_storage, value, timestamp, sequence)
    effective_command_route_test_stage!(
        storage, state_storage, value, timestamp, sequence)
    effective_command_route_test_reclaim!(storage, state_storage)
    effective_command_route_test_commit!(storage, state_storage)
    return nothing
end

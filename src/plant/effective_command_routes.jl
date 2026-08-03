#
# Target-specific effective-command publication routes
#
# The sole authority applies command semantics once. These private prepared
# routes deliver that already-effective value to every role-neutral
# target-local replica. Scalars and same-target arrays are direct; only
# cross-domain arrays acquire one dedicated typed handoff slot.
#

"""
Internal prepared-owner interface for direct and remote effective-command
publication. Implementations define
`_prepare_effective_command_publication_route_state`,
`_prepare_effective_command_publication_route!`,
`_stage_effective_command_publication_route!`,
`_reclaim_effective_command_publication_route!`,
`_commit_effective_command_publication_route!`, and
`_abandon_effective_command_publication_route!`. One route binds an exact
authority, endpoint, destination, optic, and device context; remote routes
additionally bind one handoff contract. Matching state is externally serialized
and non-reentrant. Foreign bindings and reported transfer failures raise
structured Plant errors before publication commit. An uncertain provider
failure preserves fail-stop route ownership and may rethrow the provider error.
"""
abstract type _PreparedEffectiveCommandPublicationRoute end

# Unlike the metadata-oriented `_FixedPlantRegistry`, this execution registry
# retains the concrete type of every prepared route. The warmed serial fan-out
# can therefore specialize direct and remote lanes without runtime dispatch.
struct _PreparedEffectiveCommandPublicationRoutes{R<:Tuple}
    storage::R
end

struct _EffectiveCommandPublicationRoutesState{S<:Tuple}
    storage::S
end

Base.length(routes::_PreparedEffectiveCommandPublicationRoutes) =
    length(routes.storage)
Base.iterate(routes::_PreparedEffectiveCommandPublicationRoutes, state...) =
    iterate(routes.storage, state...)
Base.length(state::_EffectiveCommandPublicationRoutesState) =
    length(state.storage)
Base.iterate(state::_EffectiveCommandPublicationRoutesState, cursor...) =
    iterate(state.storage, cursor...)

mutable struct _EffectiveCommandPublicationRouteIdentity end

@enum _EffectiveCommandPublicationRoutePhase::UInt8 begin
    _EffectiveCommandRouteIdle = 0x01
    _EffectiveCommandRoutePreparing = 0x02
    _EffectiveCommandRouteValidated = 0x03
    _EffectiveCommandRouteSubmitted = 0x04
    _EffectiveCommandRouteCompleted = 0x05
    _EffectiveCommandRouteStaged = 0x06
    _EffectiveCommandRouteReclaimed = 0x07
    _EffectiveCommandRouteFailed = 0x08
    _EffectiveCommandRouteUncertain = 0x09
end

struct _EffectiveCommandPublicationRouteBinding{A,E}
    authority::A
    endpoint::E
end

mutable struct _DirectEffectiveCommandPublicationRouteState{I,P,V}
    route_identity::I
    validated_publication::Base.RefValue{P}
    validated_value::Base.RefValue{V}
    phase::_EffectiveCommandPublicationRoutePhase
    has_validated_publication::Bool
end

struct _PreparedDirectEffectiveCommandPublicationRoute{P,I,B,A,O,X} <:
    _PreparedEffectiveCommandPublicationRoute
    partitions::P
    identity::I
    binding::B
    authority_endpoint::A
    destination::O
    optic::ControllableOpticID
    destination_context::X
end

struct _TransferredEffectiveCommandPublication{R,P}
    route_identity::R
    publication::P
    handoff_slot::UInt32
    handoff_generation::UInt64
end

struct _EffectiveCommandHandoffContract{M,R,I,D,S,A} <:
    AbstractHandoffPayloadContract{M}
    route_identity::R
    authority_identity::I
    destination_target::D
    optic::ControllableOpticID
    schema::S
    payload_axes::A
end

mutable struct _RemoteEffectiveCommandPublicationRouteState{R,H,F,M,D}
    route_identity::R
    handoff::H
    reference::Base.RefValue{F}
    transferred_publication::Base.RefValue{M}
    borrowed_publication::Base.RefValue{M}
    borrowed_payload::Base.RefValue{D}
    phase::_EffectiveCommandPublicationRoutePhase
    has_reference::Bool
end

struct _PreparedRemoteEffectiveCommandPublicationRoute{P,I,B,A,O,C,X} <:
    _PreparedEffectiveCommandPublicationRoute
    partitions::P
    identity::I
    binding::B
    authority_endpoint::A
    destination::O
    optic::ControllableOpticID
    handoff_contract::C
    destination_context::X
end

@noinline function _effective_command_route_preparation_error(
    reason::Symbol,
    message::AbstractString,
)
    throw(PlantPreparationError(
        :effective_command_route, reason, String(message)))
end

@noinline function _effective_command_route_error(
    reason::Symbol,
    message::AbstractString,
)
    throw(PlantCommandError(
        :effective_command_route, reason, String(message)))
end

@inline function _require_effective_command_route_state(route, state)
    state.route_identity === route.identity || _effective_command_route_error(
        :foreign_state,
        "effective-command route state belongs to another prepared route",
    )
    return nothing
end

@inline handoff_payload_eltype(
    contract::_EffectiveCommandHandoffContract) =
    command_numeric_type(contract.schema)
@inline handoff_payload_axes(
    contract::_EffectiveCommandHandoffContract) = contract.payload_axes

@inline function _require_effective_command_handoff_publication(
    condition::Bool,
    reason::Symbol,
    message::AbstractString,
)
    condition || throw(PlantPreparationError(
        :handoff_publication, reason, String(message)))
    return nothing
end

function validate_handoff_publication(
    contract::_EffectiveCommandHandoffContract,
    transferred::_TransferredEffectiveCommandPublication,
)
    publication = transferred.publication
    _require_effective_command_handoff_publication(
        transferred.route_identity === contract.route_identity,
        :foreign_route,
        "effective-command transfer belongs to another prepared route",
    )
    _require_effective_command_handoff_publication(
        command_authority_identity(publication) ===
            contract.authority_identity,
        :foreign_authority,
        "effective-command transfer belongs to another authority",
    )
    _require_effective_command_handoff_publication(
        compute_device(publication) == contract.destination_target,
        :wrong_target,
        "effective-command transfer names another destination target",
    )
    _require_effective_command_handoff_publication(
        controllable_optic_id(publication) == contract.optic,
        :foreign_optic,
        "effective-command transfer names another controllable optic",
    )
    _require_effective_command_handoff_publication(
        _effective_command_publication_schema_matches(
            contract.schema, command_schema(publication)),
        :schema,
        "effective-command transfer schema does not match its route",
    )
    _require_effective_command_handoff_publication(
        !iszero(transferred.handoff_slot),
        :invalid_slot,
        "effective-command transfer handoff slot must be positive",
    )
    _require_effective_command_handoff_publication(
        !iszero(transferred.handoff_generation),
        :invalid_generation,
        "effective-command transfer handoff generation must be positive",
    )
    return nothing
end

@inline _effective_command_route_destination(
    route::_PreparedEffectiveCommandPublicationRoute) = route.destination
@inline compute_device(route::_PreparedEffectiveCommandPublicationRoute) =
    compute_device(_effective_command_route_destination(route))
@inline command_endpoint_id(
    route::_PreparedEffectiveCommandPublicationRoute) =
    command_endpoint_id(route.authority_endpoint)
@inline controllable_optic_id(
    route::_PreparedEffectiveCommandPublicationRoute) = route.optic

function _effective_command_route_publication(
    route::_PreparedEffectiveCommandPublicationRoute,
    timestamp::PlantTimestamp,
    sequence::UInt64,
)
    endpoint = route.authority_endpoint.endpoint
    return _effective_command_publication(
        command_authority_identity(
            prepared_command_authority(route.partitions)),
        compute_device(route),
        route.optic,
        command_schema(endpoint),
        timestamp,
        sequence,
    )
end

function _same_effective_command_publication(
    left::EffectiveCommandPublication,
    right::EffectiveCommandPublication,
)
    return command_authority_identity(left) ===
            command_authority_identity(right) &&
        compute_device(left) == compute_device(right) &&
        controllable_optic_id(left) == controllable_optic_id(right) &&
        _effective_command_publication_schema_matches(
            command_schema(left), command_schema(right)) &&
        effective_command_publication_timestamp(left) ==
            effective_command_publication_timestamp(right) &&
        effective_command_publication_sequence(left) ==
            effective_command_publication_sequence(right)
end

function _require_effective_command_route_binding(
    route::_PreparedEffectiveCommandPublicationRoute,
)
    authority = prepared_command_authority(route.partitions)
    destination = route.destination
    authority.binding === route.binding.authority ||
        _effective_command_route_error(
            :foreign_authority,
            "effective-command route belongs to another prepared authority",
        )
    route.authority_endpoint.endpoint.binding === route.binding.endpoint ||
        _effective_command_route_error(
            :foreign_endpoint,
            "effective-command route belongs to another authority endpoint",
        )
    _prepared_device_execution_compute_device(route.destination_context) ==
        compute_device(route) || _effective_command_route_error(
            :wrong_context,
            "effective-command route context occupies another target",
        )
    command_authority_identity(destination.endpoint) ===
        command_authority_identity(authority) ||
        _effective_command_route_error(
            :foreign_destination,
            "effective-command route destination belongs to another authority",
        )
    controllable_optic_id(destination) == route.optic ||
        _effective_command_route_error(
            :foreign_optic,
            "effective-command route destination names another optic",
        )
    command_endpoint_id(destination) ==
        command_endpoint_id(route.authority_endpoint) ||
        _effective_command_route_error(
            :foreign_endpoint,
            "effective-command route destination names another endpoint",
        )
    _effective_command_publication_schema_matches(
        command_schema(route.authority_endpoint.endpoint),
        command_schema(destination.endpoint),
    ) || _effective_command_route_error(
        :schema,
        "effective-command route endpoint schemas no longer match",
    )
    return nothing
end

function _require_effective_command_route_destination_available(
    route::_PreparedEffectiveCommandPublicationRoute,
)
    destination = route.destination
    endpoint_state = destination.state
    endpoint_state.has_staged_publication &&
        _effective_command_route_error(
            :endpoint_stage_pending,
            "effective-command route destination already has staged state",
        )
    destination.optic.state.has_staged_publication &&
        _effective_command_route_error(
            :optic_stage_pending,
            "effective-command route optic already has staged state",
        )
    return nothing
end

function _target_local_effective_command_optic_owner(
    partition::PreparedTargetPartition,
    optic_id::ControllableOpticID,
)
    found = nothing
    for owner in target_local_controllable_optic_owners(partition)
        controllable_optic_id(owner) == optic_id || continue
        isnothing(found) || _effective_command_route_preparation_error(
            :duplicate_replica,
            "target $(compute_device(partition)) has more than one replica " *
            "of $optic_id",
        )
        found = owner
    end
    return found
end

function _target_local_effective_command_endpoint_owner(
    owner::TargetLocalControllableOpticOwner,
    endpoint_id::CommandEndpointID,
)
    count = 0
    for endpoint in owner.prepared.endpoints
        command_endpoint_id(endpoint) == endpoint_id && (count += 1)
    end
    count == 1 || _effective_command_route_preparation_error(
        iszero(count) ? :missing_replica : :duplicate_replica,
        "target-local optic $(controllable_optic_id(owner)) has $count " *
        "replicas of $endpoint_id",
    )
    return target_local_command_endpoint_owner(owner, endpoint_id)
end

function _require_prepared_effective_command_route_endpoint(
    authority::PreparedCommandAuthority,
    binding::_PreparedPlantCommandEndpoint,
    partition::PreparedTargetPartition,
    destination::TargetLocalCommandEndpointOwner,
    optic_id::ControllableOpticID,
)
    compute_device(destination) == compute_device(partition) ||
        _effective_command_route_preparation_error(
            :wrong_target,
            "target-local endpoint does not occupy its partition target",
        )
    command_authority_identity(destination.endpoint) ===
        command_authority_identity(authority) ||
        _effective_command_route_preparation_error(
            :foreign_authority,
            "target-local endpoint belongs to another command authority",
        )
    controllable_optic_id(destination) == optic_id ||
        _effective_command_route_preparation_error(
            :foreign_optic,
            "target-local endpoint belongs to another controllable optic",
        )
    _effective_command_publication_schema_matches(
        command_schema(binding.endpoint),
        command_schema(destination.endpoint),
    ) || _effective_command_route_preparation_error(
        :schema,
        "authority and target-local endpoint schemas do not match",
    )
    return destination
end

function _prepare_direct_effective_command_publication_route(
    partitions::PreparedPlantPartitions,
    authority::PreparedCommandAuthority,
    partition::PreparedTargetPartition,
    binding::_PreparedPlantCommandEndpoint,
    destination::TargetLocalCommandEndpointOwner,
    optic_id::ControllableOpticID,
)
    return _PreparedDirectEffectiveCommandPublicationRoute(
        partitions,
        _EffectiveCommandPublicationRouteIdentity(),
        _EffectiveCommandPublicationRouteBinding(
            authority.binding, binding.endpoint.binding),
        binding,
        destination,
        optic_id,
        partition.context,
    )
end

function _prepare_remote_effective_command_publication_route(
    partitions::PreparedPlantPartitions,
    authority::PreparedCommandAuthority,
    partition::PreparedTargetPartition,
    binding::_PreparedPlantCommandEndpoint,
    destination::TargetLocalCommandEndpointOwner,
    optic_id::ControllableOpticID,
    schema::PlantCommandSchema{T,N},
) where {T,N}
    route_identity = _EffectiveCommandPublicationRouteIdentity()
    P = EffectiveCommandPublication{
        typeof(command_authority_identity(authority)),
        typeof(compute_device(partition)),
        typeof(schema),
    }
    M = _TransferredEffectiveCommandPublication{
        typeof(route_identity),P,
    }
    payload_axes = map(Base.OneTo, command_dimensions(schema))
    contract = _EffectiveCommandHandoffContract{
        M,
        typeof(route_identity),
        typeof(command_authority_identity(authority)),
        typeof(compute_device(partition)),
        typeof(schema),
        typeof(payload_axes),
    }(
        route_identity,
        command_authority_identity(authority),
        compute_device(partition),
        optic_id,
        schema,
        payload_axes,
    )
    return _PreparedRemoteEffectiveCommandPublicationRoute(
        partitions,
        route_identity,
        _EffectiveCommandPublicationRouteBinding(
            authority.binding, binding.endpoint.binding),
        binding,
        destination,
        optic_id,
        contract,
        partition.context,
    )
end

@inline function _prepare_effective_command_publication_route(
    partitions::PreparedPlantPartitions,
    authority::PreparedCommandAuthority,
    partition::PreparedTargetPartition,
    binding::_PreparedPlantCommandEndpoint{
        <:PreparedCommandEndpoint{<:PlantCommandSchema{<:Any,0}},
    },
    destination::TargetLocalCommandEndpointOwner,
    optic_id::ControllableOpticID,
)
    return _prepare_direct_effective_command_publication_route(
        partitions, authority, partition, binding, destination, optic_id)
end

function _prepare_effective_command_publication_route(
    partitions::PreparedPlantPartitions,
    authority::PreparedCommandAuthority,
    partition::PreparedTargetPartition,
    binding::_PreparedPlantCommandEndpoint,
    destination::TargetLocalCommandEndpointOwner,
    optic_id::ControllableOpticID,
)
    schema = command_schema(binding.endpoint)
    command_authority_target(authority) == compute_device(partition) &&
        return _prepare_direct_effective_command_publication_route(
            partitions, authority, partition, binding, destination, optic_id)
    return _prepare_remote_effective_command_publication_route(
        partitions,
        authority,
        partition,
        binding,
        destination,
        optic_id,
        schema,
    )
end

function _prepare_effective_command_publication_routes(
    partitions::PreparedPlantPartitions,
)
    authority = prepared_command_authority(partitions)
    routes = _PreparedEffectiveCommandPublicationRoute[]
    @inbounds for binding in authority.endpoints
        optic_slot = Int(binding.optic_slot)
        checkbounds(authority.optic_definitions, optic_slot)
        optic_id = controllable_optic_id(
            authority.optic_definitions[optic_slot])
        replica_count = 0
        for partition in prepared_partitions(partitions)
            owner = _target_local_effective_command_optic_owner(
                partition, optic_id)
            isnothing(owner) && continue
            destination = _target_local_effective_command_endpoint_owner(
                owner, command_endpoint_id(binding))
            _require_prepared_effective_command_route_endpoint(
                authority, binding, partition, destination, optic_id)
            push!(routes, _prepare_effective_command_publication_route(
                partitions,
                authority,
                partition,
                binding,
                destination,
                optic_id,
            ))
            replica_count += 1
        end
        iszero(replica_count) && _effective_command_route_preparation_error(
            :missing_replica,
            "command endpoint $(command_endpoint_id(binding)) has no " *
            "target-local replica",
        )
    end
    return _PreparedEffectiveCommandPublicationRoutes(Tuple(routes))
end

function _prepare_effective_command_publication_route_state(
    route::_PreparedDirectEffectiveCommandPublicationRoute,
)
    authority = prepared_command_authority(route.partitions)
    schema = command_schema(route.authority_endpoint.endpoint)
    P = EffectiveCommandPublication{
        typeof(command_authority_identity(authority)),
        typeof(compute_device(route)),
        typeof(schema),
    }
    V = typeof(initial_effective_command(route.authority_endpoint))
    return _DirectEffectiveCommandPublicationRouteState(
        route.identity,
        Ref{P}(),
        Ref{V}(),
        _EffectiveCommandRouteIdle,
        false,
    )
end

function _prepare_effective_command_publication_route_state(
    route::_PreparedRemoteEffectiveCommandPublicationRoute,
)
    return _prepare_effective_command_publication_route_state(
        route, route.handoff_contract)
end

function _prepare_effective_command_publication_route_state(
    route::_PreparedRemoteEffectiveCommandPublicationRoute,
    contract::_EffectiveCommandHandoffContract{M},
) where {M}
    contract === route.handoff_contract ||
        _effective_command_route_preparation_error(
            :foreign_contract,
            "effective-command handoff contract belongs to another route",
        )
    authority = prepared_command_authority(route.partitions)
    schema = command_schema(route.authority_endpoint.endpoint)
    T = command_numeric_type(schema)
    dimensions = command_dimensions(schema)
    source = _with_completed_prepared_device_execution_context(
        authority.context) do
        allocate_device_array(
            command_authority_target(authority), T, dimensions...)
    end
    destination_slot = _with_completed_prepared_device_execution_context(
        route.destination_context) do
        allocate_device_array(compute_device(route), T, dimensions...)
    end
    _require_exact_plant_array_target(
        source,
        command_authority_target(authority),
        "effective-command route source slot",
    )
    _require_exact_plant_array_target(
        destination_slot,
        compute_device(route),
        "effective-command route destination slot",
    )
    handoff = prepare_cross_domain_handoff(
        route.partitions, contract, (source,), (destination_slot,))
    F = HandoffSlotReference{typeof(handoff.identity)}
    D = typeof(destination_slot)
    return _RemoteEffectiveCommandPublicationRouteState(
        route.identity,
        handoff,
        Ref{F}(),
        Ref{M}(),
        Ref{M}(),
        Ref{D}(),
        _EffectiveCommandRouteIdle,
        false,
    )
end

function _prepare_effective_command_publication_routes_state(
    routes::_PreparedEffectiveCommandPublicationRoutes,
)
    states = map(
        _prepare_effective_command_publication_route_state,
        routes.storage,
    )
    return _EffectiveCommandPublicationRoutesState(states)
end

function _prepare_effective_command_publication_route!(
    route::_PreparedDirectEffectiveCommandPublicationRoute,
    state::_DirectEffectiveCommandPublicationRouteState,
    publication::EffectiveCommandPublication,
    value,
)
    _require_effective_command_route_state(route, state)
    state.phase == _EffectiveCommandRouteIdle ||
        _effective_command_route_error(
            :route_busy,
            "direct effective-command route is not idle",
        )
    state.phase = _EffectiveCommandRoutePreparing
    try
        _require_effective_command_route_binding(route)
        _require_effective_command_route_destination_available(route)
        destination = route.destination
        _require_routed_effective_command_publication(
            destination.endpoint,
            destination.state,
            publication,
            value,
        )
        state.validated_publication[] = publication
        state.validated_value[] = value
        state.has_validated_publication = true
    catch
        state.phase = _EffectiveCommandRouteFailed
        rethrow()
    end
    state.phase = _EffectiveCommandRouteValidated
    return nothing
end

function _prepare_effective_command_publication_route!(
    route::_PreparedRemoteEffectiveCommandPublicationRoute,
    state::_RemoteEffectiveCommandPublicationRouteState,
    publication::EffectiveCommandPublication,
    value::AbstractArray,
)
    _require_effective_command_route_state(route, state)
    state.phase == _EffectiveCommandRouteIdle ||
        _effective_command_route_error(
            :route_busy,
            "remote effective-command route is not idle",
        )
    state.phase = _EffectiveCommandRoutePreparing
    try
        _require_effective_command_route_binding(route)
        _require_effective_command_route_destination_available(route)
        schema = command_schema(route.authority_endpoint.endpoint)
        _require_target_local_effective_value_layout(
            schema,
            value,
            command_authority_target(
                prepared_command_authority(route.partitions)),
        )
    catch
        state.phase = _EffectiveCommandRouteFailed
        rethrow()
    end
    status = try_next_free_handoff_slot!(state.reference, state.handoff)
    status == HandoffTransitionSucceeded || begin
        state.phase = _EffectiveCommandRouteUncertain
        _effective_command_route_error(
            :handoff_slot,
            "remote effective-command route could not select its dedicated " *
            "handoff slot ($status)",
        )
    end
    state.has_reference = true
    reference = state.reference[]
    source = producer_handoff_payload(state.handoff, reference)
    try
        authority = prepared_command_authority(route.partitions)
        _with_completed_prepared_device_execution_context(
            authority.context) do
            copyto!(source, value)
        end
    catch
        state.phase = _EffectiveCommandRouteFailed
        rethrow()
    end
    transferred = _TransferredEffectiveCommandPublication(
        route.identity,
        publication,
        reference.slot,
        reference.generation,
    )
    state.transferred_publication[] = transferred
    submission = try
        submit_handoff!(state.handoff, reference, transferred)
    catch
        state.phase = _EffectiveCommandRouteUncertain
        rethrow()
    end
    if submission != HandoffTransitionSucceeded
        state.phase = submission == HandoffSubmissionFailed ?
            _EffectiveCommandRouteFailed :
            _EffectiveCommandRouteUncertain
        _effective_command_route_error(
            submission == HandoffSubmissionFailed ?
                :transfer_submission_failed : :transfer_submission_uncertain,
            "remote effective-command transfer submission failed ($submission)",
        )
    end
    state.phase = _EffectiveCommandRouteSubmitted
    completion = try
        complete_handoff!(state.handoff, reference)
    catch
        state.phase = _EffectiveCommandRouteUncertain
        rethrow()
    end
    if completion != HandoffTransitionSucceeded
        state.phase = completion == HandoffCompletionFailed ?
            _EffectiveCommandRouteFailed :
            _EffectiveCommandRouteUncertain
        _effective_command_route_error(
            completion == HandoffCompletionFailed ?
                :transfer_completion_failed : :transfer_completion_uncertain,
            "remote effective-command transfer completion failed ($completion)",
        )
    end
    state.phase = _EffectiveCommandRouteCompleted
    borrowed = try_borrow_completed_handoff!(
        state.borrowed_payload,
        state.borrowed_publication,
        state.handoff,
        reference,
    )
    borrowed == HandoffTransitionSucceeded || begin
        state.phase = _EffectiveCommandRouteUncertain
        _effective_command_route_error(
            :transfer_borrow_uncertain,
            "completed effective-command transfer could not be borrowed " *
            "($borrowed)",
        )
    end
    borrowed_publication = state.borrowed_publication[]
    if !(borrowed_publication.route_identity === route.identity &&
            borrowed_publication.handoff_slot == reference.slot &&
            borrowed_publication.handoff_generation == reference.generation &&
            _same_effective_command_publication(
                borrowed_publication.publication, publication))
        state.phase = _EffectiveCommandRouteUncertain
        _effective_command_route_error(
            :transfer_publication_uncertain,
            "completed effective-command transfer metadata changed",
        )
    end
    try
        destination = route.destination
        _require_routed_effective_command_publication(
            destination.endpoint,
            destination.state,
            borrowed_publication.publication,
            state.borrowed_payload[],
        )
    catch
        state.phase = _EffectiveCommandRouteFailed
        rethrow()
    end
    return nothing
end

function _prepare_effective_command_publication_route!(
    route::_PreparedRemoteEffectiveCommandPublicationRoute,
    state::_RemoteEffectiveCommandPublicationRouteState,
    ::EffectiveCommandPublication,
    value,
)
    _require_effective_command_route_state(route, state)
    state.phase == _EffectiveCommandRouteIdle ||
        _effective_command_route_error(
            :route_busy,
            "remote effective-command route is not idle",
        )
    state.phase = _EffectiveCommandRouteFailed
    _effective_command_route_error(
        :payload,
        "remote effective-command route requires an array value; got " *
        "$(typeof(value))",
    )
end

function _stage_effective_command_publication_route!(
    route::_PreparedDirectEffectiveCommandPublicationRoute,
    state::_DirectEffectiveCommandPublicationRouteState,
    publication::EffectiveCommandPublication,
    value,
)
    _require_effective_command_route_state(route, state)
    state.phase == _EffectiveCommandRouteValidated ||
        _effective_command_route_error(
            :route_not_validated,
            "direct effective-command route was not validated",
        )
    state.has_validated_publication || begin
        state.phase = _EffectiveCommandRouteFailed
        _effective_command_route_error(
            :route_not_validated,
            "direct effective-command route has no validated publication",
        )
    end
    _same_effective_command_publication(
        state.validated_publication[], publication) || begin
        state.phase = _EffectiveCommandRouteFailed
        _effective_command_route_error(
            :validated_publication_changed,
            "direct effective-command publication changed before staging",
        )
    end
    state.validated_value[] === value || begin
        state.phase = _EffectiveCommandRouteFailed
        _effective_command_route_error(
            :validated_value_changed,
            "direct effective-command value changed before staging",
        )
    end
    try
        _with_completed_prepared_device_execution_context(
            route.destination_context) do
            _stage_prevalidated_effective_command_publication!(
                route.destination, publication, value)
        end
    catch
        state.phase = _EffectiveCommandRouteUncertain
        rethrow()
    end
    state.phase = _EffectiveCommandRouteStaged
    return nothing
end

function _stage_effective_command_publication_route!(
    route::_PreparedRemoteEffectiveCommandPublicationRoute,
    state::_RemoteEffectiveCommandPublicationRouteState,
    publication::EffectiveCommandPublication,
    value,
)
    _require_effective_command_route_state(route, state)
    state.phase == _EffectiveCommandRouteCompleted ||
        _effective_command_route_error(
            :route_not_completed,
            "remote effective-command route has no completed payload",
        )
    borrowed = state.borrowed_publication[]
    _same_effective_command_publication(
        borrowed.publication, publication) || begin
        state.phase = _EffectiveCommandRouteUncertain
        _effective_command_route_error(
            :transfer_publication_uncertain,
            "remote effective-command publication changed before staging",
        )
    end
    try
        _with_completed_prepared_device_execution_context(
            route.destination_context) do
            _stage_prevalidated_effective_command_publication!(
                route.destination,
                borrowed.publication,
                state.borrowed_payload[],
            )
        end
    catch
        state.phase = _EffectiveCommandRouteUncertain
        rethrow()
    end
    state.phase = _EffectiveCommandRouteStaged
    return nothing
end

@inline function _reclaim_effective_command_publication_route!(
    route::_PreparedDirectEffectiveCommandPublicationRoute,
    state::_DirectEffectiveCommandPublicationRouteState,
)
    _require_effective_command_route_state(route, state)
    state.phase == _EffectiveCommandRouteStaged ||
        _effective_command_route_error(
            :route_not_staged,
            "direct effective-command route must be staged before commit",
        )
    state.phase = _EffectiveCommandRouteReclaimed
    return nothing
end

function _reclaim_effective_command_publication_route!(
    route::_PreparedRemoteEffectiveCommandPublicationRoute,
    state::_RemoteEffectiveCommandPublicationRouteState,
)
    _require_effective_command_route_state(route, state)
    state.phase == _EffectiveCommandRouteStaged ||
        _effective_command_route_error(
            :route_not_staged,
            "remote effective-command route must be staged before reclaim",
        )
    status = reclaim_handoff!(state.handoff, state.reference[])
    status == HandoffTransitionSucceeded || begin
        state.phase = _EffectiveCommandRouteUncertain
        _effective_command_route_error(
            :handoff_reclaim_uncertain,
            "remote effective-command handoff could not be reclaimed ($status)",
        )
    end
    state.phase = _EffectiveCommandRouteReclaimed
    state.has_reference = false
    return nothing
end

@inline function _commit_effective_command_publication_route!(
    route::_PreparedDirectEffectiveCommandPublicationRoute,
    state::_DirectEffectiveCommandPublicationRouteState,
)
    _require_effective_command_route_state(route, state)
    state.phase == _EffectiveCommandRouteReclaimed ||
        _effective_command_route_error(
            :route_not_reclaimed,
            "direct effective-command route must be reclaimed before commit",
        )
    try
        _with_completed_prepared_device_execution_context(
            route.destination_context) do
            _commit_effective_command_publication!(route.destination)
        end
    catch
        state.phase = _EffectiveCommandRouteUncertain
        rethrow()
    end
    state.phase = _EffectiveCommandRouteIdle
    state.has_validated_publication = false
    return nothing
end

function _commit_effective_command_publication_route!(
    route::_PreparedRemoteEffectiveCommandPublicationRoute,
    state::_RemoteEffectiveCommandPublicationRouteState,
)
    _require_effective_command_route_state(route, state)
    state.phase == _EffectiveCommandRouteReclaimed ||
        _effective_command_route_error(
            :handoff_not_reclaimed,
            "remote effective-command handoff must be reclaimed before commit",
        )
    try
        _with_completed_prepared_device_execution_context(
            route.destination_context) do
            _commit_effective_command_publication!(route.destination)
        end
    catch
        state.phase = _EffectiveCommandRouteUncertain
        rethrow()
    end
    state.phase = _EffectiveCommandRouteIdle
    return nothing
end

@inline function _abandon_effective_command_publication_route!(
    route::_PreparedDirectEffectiveCommandPublicationRoute,
    state::_DirectEffectiveCommandPublicationRouteState,
)
    _require_effective_command_route_state(route, state)
    if state.phase == _EffectiveCommandRouteIdle ||
            state.phase == _EffectiveCommandRoutePreparing ||
            state.phase == _EffectiveCommandRouteValidated ||
            state.phase == _EffectiveCommandRouteFailed
        state.phase = _EffectiveCommandRouteFailed
        state.has_validated_publication = false
        return true
    end
    state.phase = _EffectiveCommandRouteUncertain
    return false
end

function _abandon_effective_command_publication_route!(
    route::_PreparedRemoteEffectiveCommandPublicationRoute,
    state::_RemoteEffectiveCommandPublicationRouteState,
)
    _require_effective_command_route_state(route, state)
    state.phase == _EffectiveCommandRouteUncertain && return false
    state.has_reference || begin
        state.phase = _EffectiveCommandRouteFailed
        return true
    end
    reference = state.reference[]
    slot_status = handoff_slot_status(state.handoff, reference)
    if slot_status == HandoffSlotFree
        state.phase = _EffectiveCommandRouteFailed
        state.has_reference = false
        return true
    end
    if slot_status == HandoffTransferCompleted ||
            slot_status == HandoffTransferFailed
        reclaimed = reclaim_handoff!(state.handoff, reference)
        reclaimed == HandoffTransitionSucceeded || begin
            state.phase = _EffectiveCommandRouteUncertain
            return false
        end
        state.phase = _EffectiveCommandRouteFailed
        state.has_reference = false
        return true
    end
    state.phase = _EffectiveCommandRouteUncertain
    return false
end

function effective_command_route_test_replica_values(partitions)
    endpoint_id = CommandEndpointID(:effective_command_route_endpoint)
    values = Any[]
    physical = Any[]
    for partition in prepared_partitions(partitions)
        owner = only(target_local_controllable_optic_owners(partition))
        push!(values, effective_command(
            target_local_command_endpoint_state(owner, endpoint_id)))
        push!(physical, getfield(
            target_local_controllable_optic_state(owner), :physical).active)
    end
    return values, physical
end

function captured_effective_command_route_error(f)
    try
        f()
        return nothing
    catch error
        return error
    end
end

@testset "Direct and remote effective-command publication routes" begin
    reset_handoff_test_transfer_controls!()
    partitions, _ = effective_command_route_test_partitions()
    routes = Plant._prepare_effective_command_publication_routes(partitions)
    route_states =
        Plant._prepare_effective_command_publication_routes_state(routes)
    @test length(routes) == 2
    direct = only(route for route in routes if route isa
        Plant._PreparedDirectEffectiveCommandPublicationRoute)
    remote = only(route for route in routes if route isa
        Plant._PreparedRemoteEffectiveCommandPublicationRoute)
    direct_state = route_states.storage[findfirst(
        route -> route isa
            Plant._PreparedDirectEffectiveCommandPublicationRoute,
        routes.storage)]
    remote_state = route_states.storage[findfirst(
        route -> route isa
            Plant._PreparedRemoteEffectiveCommandPublicationRoute,
        routes.storage)]
    @test compute_device(direct) == HostComputeDevice()
    @test compute_device(remote) == HANDOFF_TEST_ACCELERATOR
    @test handoff_source_target(remote_state.handoff) == HostComputeDevice()
    @test handoff_destination_target(remote_state.handoff) ==
        HANDOFF_TEST_ACCELERATOR
    @test handoff_capacity(remote_state.handoff) == 1
    @test handoff_payload_bytes(remote_state.handoff) == 2 * sizeof(Float64)

    authority = prepared_command_authority(partitions)
    authority_state = CommandAuthorityState(authority)
    candidate = effective_command_route_test_candidate!(
        authority, authority_state, [2.0, 3.0])
    @test candidate == [3.0, 4.0]
    timestamp = PlantTimestamp(10)
    run_effective_command_publication_route_conformance!(
        routes, route_states, candidate, timestamp, UInt64(1))

    values, physical = effective_command_route_test_replica_values(partitions)
    @test all(value ->
        effective_command_route_test_host_value(value) == [3.0, 4.0], values)
    @test all(value ->
        effective_command_route_test_host_value(value) == [3.0, 4.0], physical)
    for partition in prepared_partitions(partitions)
        owner = only(target_local_controllable_optic_owners(partition))
        state = target_local_command_endpoint_state(
            owner, CommandEndpointID(:effective_command_route_endpoint))
        @test last_effective_command_publication_timestamp(state) == timestamp
        @test last_effective_command_publication_sequence(state) ==
            EffectiveCommandPublicationSequence(1)
    end
    @test remote_state.phase == Plant._EffectiveCommandRouteIdle
    @test direct_state.phase == Plant._EffectiveCommandRouteIdle
    @test !direct_state.has_validated_publication
    @test HANDOFF_TEST_ACTIVE_DEVICE[] == zero(UInt32)

    accelerator_partitions, _ = effective_command_route_test_partitions(
        alpha_target=HANDOFF_TEST_ACCELERATOR,
        beta_target=HANDOFF_TEST_ACCELERATOR,
        command_authority_target=HANDOFF_TEST_ACCELERATOR,
    )
    accelerator_routes =
        Plant._prepare_effective_command_publication_routes(
            accelerator_partitions)
    @test length(accelerator_routes) == 1
    @test only(accelerator_routes) isa
        Plant._PreparedDirectEffectiveCommandPublicationRoute
    @test compute_device(only(accelerator_routes)) == HANDOFF_TEST_ACCELERATOR

    scalar_partitions, _ = effective_command_route_test_partitions(
        dimensions=(), semantics=AbsoluteCommand)
    scalar_routes = Plant._prepare_effective_command_publication_routes(
        scalar_partitions)
    scalar_route_states =
        Plant._prepare_effective_command_publication_routes_state(scalar_routes)
    @test length(scalar_routes) == 2
    @test all(route -> route isa
        Plant._PreparedDirectEffectiveCommandPublicationRoute, scalar_routes)
    scalar_authority = prepared_command_authority(scalar_partitions)
    scalar_authority_state = CommandAuthorityState(scalar_authority)
    scalar_candidate = effective_command_route_test_candidate!(
        scalar_authority, scalar_authority_state, 2.0)
    run_effective_command_publication_route_conformance!(
        scalar_routes,
        scalar_route_states,
        scalar_candidate,
        PlantTimestamp(11),
        UInt64(1),
    )
    scalar_values, scalar_physical =
        effective_command_route_test_replica_values(scalar_partitions)
    @test all(value -> value == 2.0, scalar_values)
    @test all(value -> value == 2.0, scalar_physical)
    @test HANDOFF_TEST_ACTIVE_DEVICE[] == zero(UInt32)
end

@testset "Effective-command route failure ownership" begin
    reset_handoff_test_transfer_controls!()
    partitions, _ = effective_command_route_test_partitions()
    routes = Plant._prepare_effective_command_publication_routes(partitions)
    route_states =
        Plant._prepare_effective_command_publication_routes_state(routes)
    remote = only(route for route in routes if route isa
        Plant._PreparedRemoteEffectiveCommandPublicationRoute)
    remote_state = route_states.storage[findfirst(
        route -> route isa
            Plant._PreparedRemoteEffectiveCommandPublicationRoute,
        routes.storage)]
    authority = prepared_command_authority(partitions)
    authority_state = CommandAuthorityState(authority)
    candidate = effective_command_route_test_candidate!(
        authority, authority_state, [1.0, 2.0])
    publication = Plant._effective_command_route_publication(
        remote, PlantTimestamp(20), UInt64(1))

    direct = only(route for route in routes if route isa
        Plant._PreparedDirectEffectiveCommandPublicationRoute)
    direct_state = route_states.storage[findfirst(
        route -> route isa
            Plant._PreparedDirectEffectiveCommandPublicationRoute,
        routes.storage)]
    direct_publication = Plant._effective_command_route_publication(
        direct, PlantTimestamp(20), UInt64(1))
    direct_stage_error = captured_effective_command_route_error(() ->
        Plant._stage_effective_command_publication_route!(
            direct, direct_state, direct_publication, candidate))
    @test direct_stage_error isa PlantCommandError
    if direct_stage_error isa PlantCommandError
        @test direct_stage_error.reason === :route_not_validated
    end

    foreign_partitions, _ = effective_command_route_test_partitions()
    foreign_routes =
        Plant._prepare_effective_command_publication_routes(
            foreign_partitions)
    foreign_states =
        Plant._prepare_effective_command_publication_routes_state(
            foreign_routes)
    foreign_direct_state = foreign_states.storage[findfirst(
        route -> route isa
            Plant._PreparedDirectEffectiveCommandPublicationRoute,
        foreign_routes.storage)]
    foreign_state_error = captured_effective_command_route_error(() ->
        Plant._prepare_effective_command_publication_route!(
            direct,
            foreign_direct_state,
            direct_publication,
            candidate,
        ))
    @test foreign_state_error isa PlantCommandError
    if foreign_state_error isa PlantCommandError
        @test foreign_state_error.reason === :foreign_state
    end
    @test foreign_direct_state.phase == Plant._EffectiveCommandRouteIdle
    foreign_remote = only(route for route in foreign_routes if route isa
        Plant._PreparedRemoteEffectiveCommandPublicationRoute)
    foreign_contract_error = captured_effective_command_route_error(() ->
        Plant._prepare_effective_command_publication_route_state(
            remote, foreign_remote.handoff_contract))
    @test foreign_contract_error isa PlantPreparationError
    if foreign_contract_error isa PlantPreparationError
        @test foreign_contract_error.component === :effective_command_route
        @test foreign_contract_error.reason === :foreign_contract
    end

    HANDOFF_TEST_FAIL_COMPLETION[] = true
    failure = captured_effective_command_route_error(() ->
        Plant._prepare_effective_command_publication_route!(
            remote, remote_state, publication, candidate))
    HANDOFF_TEST_FAIL_COMPLETION[] = false
    @test failure isa PlantCommandError
    if failure isa PlantCommandError
        @test failure.stage === :effective_command_route
        @test failure.reason === :transfer_completion_failed
    end
    reference = remote_state.reference[]
    @test handoff_slot_status(remote_state.handoff, reference) ==
        HandoffTransferFailed
    @test Plant._abandon_effective_command_publication_route!(
        remote, remote_state)
    @test remote_state.phase == Plant._EffectiveCommandRouteFailed
    @test !remote_state.has_reference
    @test all(value ->
        effective_command_route_test_host_value(value) == [1.0, 1.0],
        first(effective_command_route_test_replica_values(partitions)))

    uncertain_partitions, _ = effective_command_route_test_partitions()
    uncertain_routes = Plant._prepare_effective_command_publication_routes(
        uncertain_partitions)
    uncertain_route_states =
        Plant._prepare_effective_command_publication_routes_state(
            uncertain_routes)
    uncertain_remote = only(route for route in uncertain_routes if route isa
        Plant._PreparedRemoteEffectiveCommandPublicationRoute)
    uncertain_remote_state = uncertain_route_states.storage[findfirst(
        route -> route isa
            Plant._PreparedRemoteEffectiveCommandPublicationRoute,
        uncertain_routes.storage)]
    uncertain_authority = prepared_command_authority(uncertain_partitions)
    uncertain_state = CommandAuthorityState(uncertain_authority)
    uncertain_candidate = effective_command_route_test_candidate!(
        uncertain_authority, uncertain_state, [1.0, 2.0])
    uncertain_publication = Plant._effective_command_route_publication(
        uncertain_remote, PlantTimestamp(20), UInt64(1))
    HANDOFF_TEST_THROW_COMPLETION[] = true
    @test_throws ErrorException begin
        Plant._prepare_effective_command_publication_route!(
            uncertain_remote,
            uncertain_remote_state,
            uncertain_publication,
            uncertain_candidate,
        )
    end
    HANDOFF_TEST_THROW_COMPLETION[] = false
    uncertain_reference = uncertain_remote_state.reference[]
    @test handoff_slot_status(
        uncertain_remote_state.handoff, uncertain_reference) ==
        HandoffTransferUncertain
    @test !Plant._abandon_effective_command_publication_route!(
        uncertain_remote, uncertain_remote_state)
    @test uncertain_remote_state.phase ==
        Plant._EffectiveCommandRouteUncertain
    @test all(value ->
        effective_command_route_test_host_value(value) == [1.0, 1.0],
        first(effective_command_route_test_replica_values(
            uncertain_partitions)))

    changed_publication_partitions, _ =
        effective_command_route_test_partitions()
    changed_publication_routes =
        Plant._prepare_effective_command_publication_routes(
            changed_publication_partitions)
    changed_publication_states =
        Plant._prepare_effective_command_publication_routes_state(
            changed_publication_routes)
    changed_publication_direct = only(
        route for route in changed_publication_routes if route isa
            Plant._PreparedDirectEffectiveCommandPublicationRoute)
    changed_publication_state =
        changed_publication_states.storage[findfirst(
            route -> route isa
                Plant._PreparedDirectEffectiveCommandPublicationRoute,
            changed_publication_routes.storage)]
    changed_publication_authority =
        prepared_command_authority(changed_publication_partitions)
    changed_publication_authority_state =
        CommandAuthorityState(changed_publication_authority)
    changed_publication_candidate = effective_command_route_test_candidate!(
        changed_publication_authority,
        changed_publication_authority_state,
        [1.0, 2.0],
    )
    validated_publication = Plant._effective_command_route_publication(
        changed_publication_direct, PlantTimestamp(21), UInt64(1))
    Plant._prepare_effective_command_publication_route!(
        changed_publication_direct,
        changed_publication_state,
        validated_publication,
        changed_publication_candidate,
    )
    changed_publication = Plant._effective_command_route_publication(
        changed_publication_direct, PlantTimestamp(21), UInt64(2))
    changed_publication_error = captured_effective_command_route_error(() ->
        Plant._stage_effective_command_publication_route!(
            changed_publication_direct,
            changed_publication_state,
            changed_publication,
            changed_publication_candidate,
        ))
    @test changed_publication_error isa PlantCommandError
    if changed_publication_error isa PlantCommandError
        @test changed_publication_error.reason ===
            :validated_publication_changed
    end
    @test changed_publication_state.phase ==
        Plant._EffectiveCommandRouteFailed

    changed_value_partitions, _ = effective_command_route_test_partitions()
    changed_value_routes =
        Plant._prepare_effective_command_publication_routes(
            changed_value_partitions)
    changed_value_states =
        Plant._prepare_effective_command_publication_routes_state(
            changed_value_routes)
    changed_value_direct = only(
        route for route in changed_value_routes if route isa
            Plant._PreparedDirectEffectiveCommandPublicationRoute)
    changed_value_state = changed_value_states.storage[findfirst(
        route -> route isa
            Plant._PreparedDirectEffectiveCommandPublicationRoute,
        changed_value_routes.storage)]
    changed_value_authority = prepared_command_authority(
        changed_value_partitions)
    changed_value_authority_state =
        CommandAuthorityState(changed_value_authority)
    validated_value = effective_command_route_test_candidate!(
        changed_value_authority,
        changed_value_authority_state,
        [1.0, 2.0],
    )
    changed_value_publication = Plant._effective_command_route_publication(
        changed_value_direct, PlantTimestamp(22), UInt64(1))
    Plant._prepare_effective_command_publication_route!(
        changed_value_direct,
        changed_value_state,
        changed_value_publication,
        validated_value,
    )
    changed_value_error = captured_effective_command_route_error(() ->
        Plant._stage_effective_command_publication_route!(
            changed_value_direct,
            changed_value_state,
            changed_value_publication,
            copy(validated_value),
        ))
    @test changed_value_error isa PlantCommandError
    if changed_value_error isa PlantCommandError
        @test changed_value_error.reason === :validated_value_changed
    end
    @test changed_value_state.phase == Plant._EffectiveCommandRouteFailed

    staged_partitions, _ = effective_command_route_test_partitions()
    staged_routes = Plant._prepare_effective_command_publication_routes(
        staged_partitions)
    staged_states = Plant._prepare_effective_command_publication_routes_state(
        staged_routes)
    staged_direct = only(route for route in staged_routes if route isa
        Plant._PreparedDirectEffectiveCommandPublicationRoute)
    staged_direct_state = staged_states.storage[findfirst(
        route -> route isa
            Plant._PreparedDirectEffectiveCommandPublicationRoute,
        staged_routes.storage)]
    staged_authority = prepared_command_authority(staged_partitions)
    staged_authority_state = CommandAuthorityState(staged_authority)
    staged_candidate = effective_command_route_test_candidate!(
        staged_authority, staged_authority_state, [1.0, 2.0])
    staged_publication = Plant._effective_command_route_publication(
        staged_direct, PlantTimestamp(23), UInt64(1))
    Plant._prepare_effective_command_publication_route!(
        staged_direct,
        staged_direct_state,
        staged_publication,
        staged_candidate,
    )
    Plant._stage_effective_command_publication_route!(
        staged_direct,
        staged_direct_state,
        staged_publication,
        staged_candidate,
    )
    @test !Plant._abandon_effective_command_publication_route!(
        staged_direct, staged_direct_state)
    @test staged_direct_state.phase ==
        Plant._EffectiveCommandRouteUncertain
    @test all(value ->
        effective_command_route_test_host_value(value) == [1.0, 1.0],
        first(effective_command_route_test_replica_values(
            staged_partitions)))

    throwing_partitions, _ = effective_command_route_test_partitions()
    throwing_routes = Plant._prepare_effective_command_publication_routes(
        throwing_partitions)
    throwing_states =
        Plant._prepare_effective_command_publication_routes_state(
            throwing_routes)
    throwing_direct_index = findfirst(route -> route isa
        Plant._PreparedDirectEffectiveCommandPublicationRoute,
        throwing_routes.storage)
    throwing_remote_index = findfirst(route -> route isa
        Plant._PreparedRemoteEffectiveCommandPublicationRoute,
        throwing_routes.storage)
    throwing_direct = throwing_routes.storage[throwing_direct_index]
    throwing_remote = throwing_routes.storage[throwing_remote_index]
    throwing_direct_state = throwing_states.storage[throwing_direct_index]
    throwing_remote_state = throwing_states.storage[throwing_remote_index]
    throwing_authority = prepared_command_authority(throwing_partitions)
    throwing_authority_state = CommandAuthorityState(throwing_authority)
    throwing_candidate = effective_command_route_test_candidate!(
        throwing_authority, throwing_authority_state, [1.0, 2.0])
    throwing_direct_publication = Plant._effective_command_route_publication(
        throwing_direct, PlantTimestamp(24), UInt64(1))
    Plant._prepare_effective_command_publication_route!(
        throwing_direct,
        throwing_direct_state,
        throwing_direct_publication,
        throwing_candidate,
    )
    EFFECTIVE_COMMAND_ROUTE_TEST_THROW_STAGE_AFTER_MUTATION[] = true
    direct_physical_stage_error = captured_effective_command_route_error(() ->
        Plant._stage_effective_command_publication_route!(
            throwing_direct,
            throwing_direct_state,
            throwing_direct_publication,
            throwing_candidate,
        ))
    EFFECTIVE_COMMAND_ROUTE_TEST_THROW_STAGE_AFTER_MUTATION[] = false
    @test direct_physical_stage_error isa ErrorException
    @test throwing_direct_state.phase ==
        Plant._EffectiveCommandRouteUncertain
    @test !Plant._abandon_effective_command_publication_route!(
        throwing_direct, throwing_direct_state)
    direct_reuse_error = captured_effective_command_route_error(() ->
        Plant._prepare_effective_command_publication_route!(
            throwing_direct,
            throwing_direct_state,
            throwing_direct_publication,
            throwing_candidate,
        ))
    @test direct_reuse_error isa PlantCommandError
    if direct_reuse_error isa PlantCommandError
        @test direct_reuse_error.reason === :route_busy
    end

    throwing_remote_publication = Plant._effective_command_route_publication(
        throwing_remote, PlantTimestamp(24), UInt64(1))
    Plant._prepare_effective_command_publication_route!(
        throwing_remote,
        throwing_remote_state,
        throwing_remote_publication,
        throwing_candidate,
    )
    EFFECTIVE_COMMAND_ROUTE_TEST_THROW_STAGE_AFTER_MUTATION[] = true
    remote_physical_stage_error = captured_effective_command_route_error(() ->
        Plant._stage_effective_command_publication_route!(
            throwing_remote,
            throwing_remote_state,
            throwing_remote_publication,
            throwing_candidate,
        ))
    EFFECTIVE_COMMAND_ROUTE_TEST_THROW_STAGE_AFTER_MUTATION[] = false
    @test remote_physical_stage_error isa ErrorException
    @test throwing_remote_state.phase ==
        Plant._EffectiveCommandRouteUncertain
    @test !Plant._abandon_effective_command_publication_route!(
        throwing_remote, throwing_remote_state)
    remote_reuse_error = captured_effective_command_route_error(() ->
        Plant._prepare_effective_command_publication_route!(
            throwing_remote,
            throwing_remote_state,
            throwing_remote_publication,
            throwing_candidate,
        ))
    @test remote_reuse_error isa PlantCommandError
    if remote_reuse_error isa PlantCommandError
        @test remote_reuse_error.reason === :route_busy
    end
    @test all(value ->
        effective_command_route_test_host_value(value) == [1.0, 1.0],
        first(effective_command_route_test_replica_values(
            throwing_partitions)))

    busy_partitions, _ = effective_command_route_test_partitions()
    busy_routes = Plant._prepare_effective_command_publication_routes(
        busy_partitions)
    busy_states = Plant._prepare_effective_command_publication_routes_state(
        busy_routes)
    busy_remote = only(route for route in busy_routes if route isa
        Plant._PreparedRemoteEffectiveCommandPublicationRoute)
    busy_remote_state = busy_states.storage[findfirst(
        route -> route isa
            Plant._PreparedRemoteEffectiveCommandPublicationRoute,
        busy_routes.storage)]
    busy_authority = prepared_command_authority(busy_partitions)
    busy_authority_state = CommandAuthorityState(busy_authority)
    busy_candidate = effective_command_route_test_candidate!(
        busy_authority, busy_authority_state, [1.0, 2.0])
    busy_publication = Plant._effective_command_route_publication(
        busy_remote, PlantTimestamp(24), UInt64(1))
    Plant._prepare_effective_command_publication_route!(
        busy_remote, busy_remote_state, busy_publication, busy_candidate)
    busy_payload_error = captured_effective_command_route_error(() ->
        Plant._prepare_effective_command_publication_route!(
            busy_remote, busy_remote_state, busy_publication, 1.0))
    @test busy_payload_error isa PlantCommandError
    if busy_payload_error isa PlantCommandError
        @test busy_payload_error.reason === :route_busy
    end
    @test busy_remote_state.phase == Plant._EffectiveCommandRouteCompleted
    @test HANDOFF_TEST_ACTIVE_DEVICE[] == zero(UInt32)
end

@testset "Warmed effective-command route allocation" begin
    reset_handoff_test_transfer_controls!()
    partitions, _ = effective_command_route_test_partitions()
    routes = Plant._prepare_effective_command_publication_routes(partitions)
    route_states =
        Plant._prepare_effective_command_publication_routes_state(routes)
    candidate = [4.0, 5.0]
    run_effective_command_publication_route_conformance!(
        routes, route_states, candidate, PlantTimestamp(1), UInt64(1))
    storage = routes.storage
    state_storage = route_states.storage
    direct = first(storage)
    remote = last(storage)
    direct_state = first(state_storage)
    remote_state = last(state_storage)
    direct_publication = Plant._effective_command_route_publication(
        direct, PlantTimestamp(2), UInt64(2))
    remote_publication = Plant._effective_command_route_publication(
        remote, PlantTimestamp(2), UInt64(2))
    direct_prepare_allocations = @allocated begin
        Plant._prepare_effective_command_publication_route!(
            direct, direct_state, direct_publication, candidate)
    end
    remote_prepare_allocations = @allocated begin
        Plant._prepare_effective_command_publication_route!(
            remote, remote_state, remote_publication, candidate)
    end
    stage_allocations = @allocated effective_command_route_test_stage!(
        storage, state_storage, candidate, PlantTimestamp(2), UInt64(2))
    reclaim_allocations = @allocated effective_command_route_test_reclaim!(
        storage, state_storage)
    commit_allocations = @allocated effective_command_route_test_commit!(
        storage, state_storage)
    @test @inferred(run_effective_command_publication_route_conformance!(
        routes,
        route_states,
        candidate,
        PlantTimestamp(3),
        UInt64(3),
    )) === nothing
    if !coverage_instrumented()
        @test direct_prepare_allocations == 0
        @test remote_prepare_allocations == 0
        @test stage_allocations == 0
        @test reclaim_allocations == 0
        @test commit_allocations == 0
    end
end

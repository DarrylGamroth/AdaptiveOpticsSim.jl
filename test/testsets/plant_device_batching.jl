function device_batch_test_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

function device_batch_test_group_ordinal(prepared, id::Symbol)
    target = OpticalPathID(id)
    ordinal = findfirst(
        group -> path_execution_group_path_id(group) == target,
        prepared.path_groups,
    )
    isnothing(ordinal) && error("missing device-batch test path $id")
    return ordinal
end

function device_batch_test_inferred_owner_validation(owner, prepared)
    return device_batch_test_inferred_owner_validation(
        owner.implementation,
        owner,
        prepared,
    )
end

function device_batch_test_inferred_owner_validation(
    implementation::I,
    owner,
    prepared,
) where {I}
    return @inferred Plant._validate_device_path_batch_owner_implementation(
        implementation,
        owner,
        prepared,
    )
end

function compare_device_batch_test_runs(serial, batched)
    @test scheduler_timestamp(batched.state.scheduler) ==
        scheduler_timestamp(serial.state.scheduler)
    @test batched.state.scheduler.revision ==
        serial.state.scheduler.revision
    @test batched.state.scheduler.cursors ==
        serial.state.scheduler.cursors
    @test batched.state.path_sampled == serial.state.path_sampled
    @test batched.state.product_sequences ==
        serial.state.product_sequences
    @test batched.state.product_ready_timestamps ==
        serial.state.product_ready_timestamps

    serial_epoch = current_epoch(serial.prepared.atmosphere)
    batched_epoch = current_epoch(batched.prepared.atmosphere)
    @test epoch_time(batched_epoch) == epoch_time(serial_epoch)
    @test epoch_sequence(batched_epoch) == epoch_sequence(serial_epoch)
    @test length(batched.prepared.atmosphere_rng.streams) ==
        length(serial.prepared.atmosphere_rng.streams)
    for index in eachindex(batched.prepared.atmosphere_rng.streams)
        @test copy(batched.prepared.atmosphere_rng.streams[index].state) ==
            copy(serial.prepared.atmosphere_rng.streams[index].state)
    end

    for id in serial.path_ids
        serial_path = prepared_path(serial.plant, id)
        batched_path = prepared_path(batched.plant, id)
        @test isapprox(
            batched_path.input.opd,
            serial_path.input.opd;
            rtol=4f-5,
            atol=4f-6,
        )
        @test isapprox(
            batched_path.result.values,
            serial_path.result.values;
            rtol=5f-5,
            atol=5f-5,
        )
        serial_group = path_execution_group(
            serial.prepared,
            device_batch_test_group_ordinal(serial.prepared, id),
        )
        batched_group = path_execution_group(
            batched.prepared,
            device_batch_test_group_ordinal(batched.prepared, id),
        )
        @test copy(Plant.rng_stream_state(
            batched_group.rngs,
            Val(:provider),
        )) == copy(Plant.rng_stream_state(
            serial_group.rngs,
            Val(:provider),
        ))
    end

    for id in serial.acquisition_ids
        @test acquisition_product_sequence(
            batched.prepared,
            batched.state,
            id,
        ) == acquisition_product_sequence(
            serial.prepared,
            serial.state,
            id,
        )
        @test acquisition_product_ready_timestamp(
            batched.prepared,
            batched.state,
            id,
        ) == acquisition_product_ready_timestamp(
            serial.prepared,
            serial.state,
            id,
        )
        serial_acquisition = prepared_acquisition(serial.plant, id)
        batched_acquisition = prepared_acquisition(batched.plant, id)
        @test isapprox(
            acquisition_observation(batched_acquisition),
            acquisition_observation(serial_acquisition);
            rtol=5f-5,
            atol=5f-5,
        )
    end
    compare_device_batch_test_command_state(serial, batched)
    return nothing
end

function device_batch_test_allocation_bytes(fixture, owner)
    run_plant_events_until!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        PlantTimestamp(199_000_000),
    )
    claim = Plant.begin_optical_path_batch!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        PlantTimestamp(200_000_000),
    )
    materialization_bytes = @allocated materialize_device_path_batch!(
        owner,
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        claim,
    )
    Plant.seal_optical_path_batch_materialization!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        claim,
    )
    execution_bytes = @allocated execute_device_path_batch!(
        owner,
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        claim,
    )
    Plant.complete_optical_path_batch!(
        fixture.prepared,
        fixture.state,
        fixture.workspace,
        claim,
    )
    return (; materialization_bytes, execution_bytes)
end

@testset "Prepared device path-batch ownership" begin
    qualified_api = (
        :PreparedDevicePathBatchOwner,
        :device_path_batch_owner_count,
        :device_path_batch_owner,
        :device_path_batch_compute_device,
        :device_path_batch_backend,
        :device_path_batch_group_count,
        :device_path_batch_group_ordinal,
        :path_execution_group_device_batch_owner_ordinal,
        :materialize_device_path_batch!,
        :execute_device_path_batch!,
    )
    for name in qualified_api
        @test isdefined(Plant, name)
        @test Base.ispublic(Plant, name)
        @test !Base.isexported(Plant, name)
        @test name ∉ names(AdaptiveOpticsSim)
    end

    public_cpu = device_batch_test_fixture()
    @test @inferred(device_path_batch_owner_count(
        public_cpu.prepared,
    )) == 0
    @test all(
        ordinal -> isnothing(
            path_execution_group_device_batch_owner_ordinal(
                public_cpu.prepared,
                ordinal,
            ),
        ),
        1:path_execution_group_count(public_cpu.prepared),
    )

    singleton = device_batch_test_fixture(
        selection=Val(:all),
        include_beta=false,
        include_unequal_rate=false,
        include_lgs=false,
    )
    @test device_path_batch_owner_count(singleton.prepared) == 0

    batched = device_batch_test_fixture(selection=Val(:all))
    @test @inferred(device_path_batch_owner_count(batched.prepared)) == 1
    owner = @inferred(device_path_batch_owner(batched.prepared, 1))
    @test owner isa PreparedDevicePathBatchOwner
    @test device_batch_test_inferred_owner_validation(
        owner,
        batched.prepared,
    ) == 1
    @test device_path_batch_compute_device(owner) ==
        HostComputeDevice()
    @test device_path_batch_backend(owner) == CPUBackend()
    @test @inferred(device_path_batch_group_count(owner)) == 2
    alpha_ordinal =
        device_batch_test_group_ordinal(batched.prepared, :alpha)
    beta_ordinal =
        device_batch_test_group_ordinal(batched.prepared, :beta)
    unequal_ordinal =
        device_batch_test_group_ordinal(batched.prepared, :unequal)
    phase_offset_ordinal =
        device_batch_test_group_ordinal(batched.prepared, :phase_offset)
    origin_offset_ordinal =
        device_batch_test_group_ordinal(batched.prepared, :origin_offset)
    lgs_ordinal = device_batch_test_group_ordinal(batched.prepared, :lgs)
    @test Tuple(device_path_batch_group_ordinal(owner, index) for
        index in 1:device_path_batch_group_count(owner)) ==
        (alpha_ordinal, beta_ordinal)
    @test path_execution_group_device_batch_owner_ordinal(
        batched.prepared,
        alpha_ordinal,
    ) == 1
    @test path_execution_group_device_batch_owner_ordinal(
        batched.prepared,
        beta_ordinal,
    ) == 1
    @test isnothing(path_execution_group_device_batch_owner_ordinal(
        batched.prepared,
        unequal_ordinal,
    ))
    @test isnothing(path_execution_group_device_batch_owner_ordinal(
        batched.prepared,
        phase_offset_ordinal,
    ))
    @test isnothing(path_execution_group_device_batch_owner_ordinal(
        batched.prepared,
        origin_offset_ordinal,
    ))
    @test isnothing(path_execution_group_device_batch_owner_ordinal(
        batched.prepared,
        lgs_ordinal,
    ))

    claim = Plant.begin_optical_path_batch!(
        batched.prepared,
        batched.state,
        batched.workspace,
        PlantTimestamp(0),
    )
    retained_context = owner.implementation.context
    retained_fft_plan =
        owner.implementation.optical_batch.workspace.fft_plan
    @test Plant.optical_path_batch_due_group_count(
        batched.prepared,
        batched.state,
        batched.workspace,
        claim,
    ) == 4

    independent_materialization_error = device_batch_test_error() do
        Plant.materialize_path_execution_group!(
            batched.prepared,
            batched.state,
            batched.workspace,
            claim,
            alpha_ordinal,
        )
    end
    @test independent_materialization_error isa PlantScheduleError
    @test independent_materialization_error.reason ==
        :device_path_batch_owner_required

    foreign = device_batch_test_fixture(selection=Val(:all))
    foreign_owner = device_path_batch_owner(foreign.prepared, 1)
    foreign_owner_error = device_batch_test_error() do
        materialize_device_path_batch!(
            foreign_owner,
            batched.prepared,
            batched.state,
            batched.workspace,
            claim,
        )
    end
    @test foreign_owner_error isa PlantScheduleError
    @test foreign_owner_error.reason ==
        :foreign_device_path_batch_owner

    @test materialize_device_path_batch!(
        owner,
        batched.prepared,
        batched.state,
        batched.workspace,
        claim,
    ) === nothing
    duplicate_materialization_error = device_batch_test_error() do
        materialize_device_path_batch!(
            owner,
            batched.prepared,
            batched.state,
            batched.workspace,
            claim,
        )
    end
    @test duplicate_materialization_error isa PlantScheduleError
    @test duplicate_materialization_error.reason ==
        :invalid_device_path_batch_group_status
    for ordinal in (unequal_ordinal, lgs_ordinal)
        Plant.materialize_path_execution_group!(
            batched.prepared,
            batched.state,
            batched.workspace,
            claim,
            ordinal,
        )
    end
    Plant.seal_optical_path_batch_materialization!(
        batched.prepared,
        batched.state,
        batched.workspace,
        claim,
    )

    independent_execution_error = device_batch_test_error() do
        Plant.execute_path_execution_group!(
            batched.prepared,
            batched.state,
            batched.workspace,
            claim,
            beta_ordinal,
        )
    end
    @test independent_execution_error isa PlantScheduleError
    @test independent_execution_error.reason ==
        :device_path_batch_owner_required
    @test execute_device_path_batch!(
        owner,
        batched.prepared,
        batched.state,
        batched.workspace,
        claim,
    ) === nothing
    duplicate_execution_error = device_batch_test_error() do
        execute_device_path_batch!(
            owner,
            batched.prepared,
            batched.state,
            batched.workspace,
            claim,
        )
    end
    @test duplicate_execution_error isa PlantScheduleError
    @test duplicate_execution_error.reason ==
        :invalid_device_path_batch_group_status
    for ordinal in (unequal_ordinal, lgs_ordinal)
        Plant.execute_path_execution_group!(
            batched.prepared,
            batched.state,
            batched.workspace,
            claim,
            ordinal,
        )
    end
    @test Plant.complete_optical_path_batch!(
        batched.prepared,
        batched.state,
        batched.workspace,
        claim,
    ) == PlantTimestamp(0)
    @test owner.implementation.context === retained_context
    @test owner.implementation.optical_batch.workspace.fft_plan ===
        retained_fft_plan
    stale_owner_error = device_batch_test_error() do
        execute_device_path_batch!(
            owner,
            batched.prepared,
            batched.state,
            batched.workspace,
            claim,
        )
    end
    @test stale_owner_error isa PlantScheduleError
    @test stale_owner_error.reason == :stale_optical_path_batch_claim

    serial = device_batch_test_fixture()
    comparison = device_batch_test_fixture(selection=Val(:all))
    @test command_admission_status(
        submit_device_batch_test_command!(serial),
    ) == CommandAdmittedPending
    @test command_admission_status(
        submit_device_batch_test_command!(comparison),
    ) == CommandAdmittedPending
    horizon = PlantTimestamp(450_000_000)
    serial_count = run_plant_events_until!(
        serial.prepared,
        serial.state,
        serial.workspace,
        horizon,
    )
    comparison_count = run_plant_events_until!(
        comparison.prepared,
        comparison.state,
        comparison.workspace,
        horizon,
    )
    @test comparison_count == serial_count
    compare_device_batch_test_runs(serial, comparison)

    @test typeof(serial.prepared.path_groups) ===
        typeof(comparison.prepared.path_groups)
    @test typeof(singleton.prepared.device_path_batch_owners) ===
        typeof(comparison.prepared.device_path_batch_owners)

    if coverage_instrumented()
        @test_skip "coverage instrumentation changes allocation counts"
    else
        allocation_fixture = device_batch_test_fixture(
            selection=Val(:all),
            include_unequal_rate=false,
            include_phase_offset=false,
            include_origin_offset=false,
            include_lgs=false,
            include_physical_state=false,
        )
        allocation_owner =
            device_path_batch_owner(allocation_fixture.prepared, 1)
        allocation_bytes = device_batch_test_allocation_bytes(
            allocation_fixture,
            allocation_owner,
        )
        # The test-only host selection crosses the same topology-bounded
        # concrete implementation barrier used by accelerators. Numerical
        # outputs and workspaces are preallocated; these ceilings cover only
        # the bounded dynamic-dispatch boxes at that registry boundary.
        @test allocation_bytes.materialization_bytes <= 1_024
        @test allocation_bytes.execution_bytes <= 1_024
    end
end

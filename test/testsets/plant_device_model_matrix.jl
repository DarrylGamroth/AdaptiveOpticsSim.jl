function compare_device_model_matrix_wfs_runs(oracle, owned)
    @test scheduler_timestamp(owned.state.scheduler) ==
        scheduler_timestamp(oracle.state.scheduler)
    @test owned.state.scheduler.cursors ==
        oracle.state.scheduler.cursors
    @test owned.state.path_sampled == oracle.state.path_sampled
    oracle_epoch = current_epoch(oracle.prepared.atmosphere)
    owned_epoch = current_epoch(owned.prepared.atmosphere)
    @test epoch_time(owned_epoch) == epoch_time(oracle_epoch)
    @test epoch_sequence(owned_epoch) == epoch_sequence(oracle_epoch)
    @test length(owned.prepared.atmosphere_rng.streams) ==
        length(oracle.prepared.atmosphere_rng.streams)
    @inbounds for index in eachindex(
        owned.prepared.atmosphere_rng.streams)
        @test copy(owned.prepared.atmosphere_rng.streams[index].state) ==
            copy(oracle.prepared.atmosphere_rng.streams[index].state)
    end
    for id in oracle.path_ids
        oracle_path = prepared_path(oracle.plant, id)
        owned_path = prepared_path(owned.plant, id)
        @test isapprox(
            Array(owned_path.input.opd),
            Array(oracle_path.input.opd);
            rtol=4f-5,
            atol=4f-6,
        )
        owned_products =
            device_model_matrix_product_host(owned_path.result)
        oracle_products =
            device_model_matrix_product_host(oracle_path.result)
        @test device_model_matrix_products_approx(
            owned_products,
            oracle_products,
        )
    end
    return nothing
end

@inline device_model_matrix_products_approx(
    owned::AbstractArray,
    oracle::AbstractArray,
) = isapprox(owned, oracle; rtol=8f-5, atol=8f-5)

function device_model_matrix_products_approx(
    owned::Tuple,
    oracle::Tuple,
)
    length(owned) == length(oracle) || return false
    return all(
        device_model_matrix_products_approx(owned[index], oracle[index])
        for index in eachindex(owned)
    )
end

@testset "Prepared WFS device path-batch ownership" begin
    for (family, direction, spectral) in device_model_matrix_wfs_rows()
        oracle = device_model_matrix_wfs_fixture(
            family;
            selection=Val(:none),
            direction,
            spectral,
        )
        owned = device_model_matrix_wfs_fixture(
            family;
            selection=Val(:all),
            direction,
            spectral,
        )
        @test device_path_batch_owner_count(oracle.prepared) == 0
        @test device_path_batch_owner_count(owned.prepared) == 1
        for id in owned.path_ids
            group = only(
                group for group in owned.prepared.path_groups
                if group.id == OpticalPathID(id)
            )
            @test path_execution_requires_full_optical(
                path_execution_group_requirements(group),
            )
            @test path_execution_group_acquisition_count(group) == 0
        end
        owner = device_path_batch_owner(owned.prepared, 1)
        @test owner.implementation isa
            Plant._PreparedWFSDevicePathBatch
        @test device_path_batch_group_count(owner) == 2
        retained_context = owner.implementation.context
        retained_plans = ntuple(2) do index
            group_slot = Int(owner.group_slots[index])
            owned.prepared.path_groups[group_slot].path.execution.plan
        end
        horizon = PlantTimestamp(200_000_000)
        @test run_plant_events_until!(
            owned.prepared,
            owned.state,
            owned.workspace,
            horizon,
        ) == run_plant_events_until!(
            oracle.prepared,
            oracle.state,
            oracle.workspace,
            horizon,
        )
        @test owner.implementation.context === retained_context
        current_plans = ntuple(2) do index
            group_slot = Int(owner.group_slots[index])
            owned.prepared.path_groups[group_slot].path.execution.plan
        end
        @test current_plans === retained_plans
        compare_device_model_matrix_wfs_runs(oracle, owned)
    end
end

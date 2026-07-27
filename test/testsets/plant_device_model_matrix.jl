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
    compare_device_batch_test_command_state(oracle, owned)
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

function device_model_matrix_path_group(fixture, id::Symbol)
    return only(
        group for group in fixture.prepared.path_groups
        if group.id == OpticalPathID(id)
    )
end

function device_model_matrix_inferred_owner_validation(owner, prepared)
    return device_model_matrix_inferred_owner_validation(
        owner.implementation,
        owner,
        prepared,
    )
end

function device_model_matrix_inferred_owner_validation(
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
        @test device_model_matrix_inferred_owner_validation(
            owner,
            owned.prepared,
        ) == 1
        @test device_path_batch_group_count(owner) == 2
        retained_context = owner.implementation.context
        retained_plans = ntuple(2) do index
            group_slot = Int(owner.group_slots[index])
            owned.prepared.path_groups[group_slot].path.execution.plan
        end
        initial_command = Array(effective_command(
            owned.prepared,
            owned.state,
            :device_batch_dm_command,
        ))
        @test initial_command ==
            eltype(initial_command)[1e-9, 0, -5e-10, 0]
        @test command_admission_status(
            submit_device_batch_test_command!(oracle),
        ) == CommandAdmittedPending
        @test command_admission_status(
            submit_device_batch_test_command!(owned),
        ) == CommandAdmittedPending
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
        applied_command = Array(effective_command(
            owned.prepared,
            owned.state,
            :device_batch_dm_command,
        ))
        @test applied_command ==
            eltype(applied_command)[-5e-10, 2e-9, 0, 1e-9]
        @test last_command_application_timestamp(
            owned.state.command_applications[1],
        ) == horizon
        compare_device_model_matrix_wfs_runs(oracle, owned)

        detector_row =
            device_model_matrix_wfs_detector_row(family, direction)
        oracle_path = prepared_path(oracle.plant, :wfs_alpha)
        owned_path = prepared_path(owned.plant, :wfs_alpha)
        oracle_input =
            device_model_matrix_detector_facing_product(oracle_path.result)
        owned_input =
            device_model_matrix_detector_facing_product(owned_path.result)
        owned_input_before = copy(owned_input.values)
        oracle_detector = device_model_matrix_execute_detector(
            detector_row;
            input_map=oracle_input,
        )
        owned_detector = device_model_matrix_execute_detector(
            detector_row;
            input_map=owned_input,
        )
        @test owned_detector.prepared.plan.input_values === owned_input.values
        @test owned_detector.status_trace == oracle_detector.status_trace
        @test owned_detector.event_times == oracle_detector.event_times
        @test owned_detector.metadata_after.sensor ==
            device_model_matrix_detector_sensor_symbol(detector_row)
        @test Array(owned_detector.output) ≈
            device_model_matrix_expected_detector_frame(
                detector_row,
                owned_detector.detector,
                Array(owned_input.values),
                Float64,
            ) atol=2e-12 rtol=2e-12
        @test Array(owned_detector.output) ≈
            Array(oracle_detector.output) atol=2e-12 rtol=2e-12
        @test owned_input.values == owned_input_before
    end
end

@testset "WFS device path-batch fallback matrix" begin
    singleton = device_model_matrix_wfs_fixture(
        DeviceModelMatrixShackHartmann();
        include_second=false,
    )
    @test device_path_batch_owner_count(singleton.prepared) == 0
    @test Plant._device_path_batch_candidate(
        Val(:all),
        device_model_matrix_path_group(singleton, :wfs_alpha),
    ) !== nothing

    unequal_schedule = device_model_matrix_wfs_fixture(
        DeviceModelMatrixShackHartmann();
        second_period_ns=150_000_000,
    )
    @test device_path_batch_owner_count(unequal_schedule.prepared) == 0
    schedule_keys = map((:wfs_alpha, :wfs_beta)) do id
        Plant._device_path_batch_candidate(
            Val(:all),
            device_model_matrix_path_group(unequal_schedule, id),
        )
    end
    @test schedule_keys[1].schedule != schedule_keys[2].schedule

    unequal_origin = device_model_matrix_wfs_fixture(
        DeviceModelMatrixShackHartmann();
        second_origin=PlantTimestamp(50_000_000),
    )
    @test device_path_batch_owner_count(unequal_origin.prepared) == 0
    origin_keys = map((:wfs_alpha, :wfs_beta)) do id
        Plant._device_path_batch_candidate(
            Val(:all),
            device_model_matrix_path_group(unequal_origin, id),
        )
    end
    @test origin_keys[1].origin != origin_keys[2].origin

    mixed_family = device_model_matrix_wfs_fixture(
        DeviceModelMatrixShackHartmann();
        second_family=DeviceModelMatrixPyramid(),
    )
    @test device_path_batch_owner_count(mixed_family.prepared) == 0
    family_keys = map((:wfs_alpha, :wfs_beta)) do id
        Plant._device_path_batch_candidate(
            Val(:all),
            device_model_matrix_path_group(mixed_family, id),
        )
    end
    @test family_keys[1].execution_type !== family_keys[2].execution_type

    mixed_signature = device_model_matrix_wfs_fixture(
        DeviceModelMatrixShackHartmann();
        second_variant=1,
    )
    @test device_path_batch_owner_count(mixed_signature.prepared) == 0
    signature_keys = map((:wfs_alpha, :wfs_beta)) do id
        Plant._device_path_batch_candidate(
            Val(:all),
            device_model_matrix_path_group(mixed_signature, id),
        )
    end
    @test signature_keys[1].plan_contract !=
        signature_keys[2].plan_contract

    mixed_product = device_model_matrix_wfs_fixture(
        DeviceModelMatrixPyramid();
        second_spectral=Val(:spectral),
    )
    @test device_path_batch_owner_count(mixed_product.prepared) == 0
    product_keys = map((:wfs_alpha, :wfs_beta)) do id
        Plant._device_path_batch_candidate(
            Val(:all),
            device_model_matrix_path_group(mixed_product, id),
        )
    end
    @test product_keys[1].product_contract !=
        product_keys[2].product_contract

    for unsupported_family in (
        DeviceModelMatrixZernike(),
        DeviceModelMatrixCurvature(),
    )
        unsupported = device_model_matrix_wfs_fixture(unsupported_family)
        @test device_path_batch_owner_count(unsupported.prepared) == 0
        @test all(
            id -> Plant._device_path_batch_candidate(
                Val(:all),
                device_model_matrix_path_group(unsupported, id),
            ) === nothing,
            unsupported.path_ids,
        )
    end

    coexistence = device_model_matrix_wfs_fixture(
        DeviceModelMatrixBioEdge();
        spectral=Val(:spectral),
    )
    @test device_path_batch_owner_count(coexistence.prepared) == 1
    witness_group = device_model_matrix_path_group(
        coexistence,
        :readout_witness,
    )
    witness_ordinal = findfirst(
        group -> group.id == witness_group.id,
        coexistence.prepared.path_groups,
    )
    @test !isnothing(witness_ordinal)
    @test isnothing(path_execution_group_device_batch_owner_ordinal(
        coexistence.prepared,
        something(witness_ordinal),
    ))
end

@testset "Conventional detector device-model matrix" begin
    for row in device_model_matrix_detector_rows()
        result = device_model_matrix_execute_detector(row)
        detector = result.detector
        metadata = result.metadata_after
        expected = device_model_matrix_expected_detector_frame(
            row,
            detector,
            Float64,
        )
        @test Array(result.output) ≈ expected atol=2e-14 rtol=2e-14
        @test detector_acquisition_sequence(result.state) == UInt64(1)
        @test metadata.sensor ==
            device_model_matrix_detector_sensor_symbol(row)
        @test metadata.frame_response ==
            device_model_matrix_detector_response_symbol(row)
        @test metadata.timing_model ==
            device_model_matrix_detector_timing_symbol(row)
        @test metadata.sampling_mode ==
            device_model_matrix_detector_sampling_symbol(row)
        @test metadata.acquisition_mode ==
            device_model_matrix_detector_acquisition_symbol(row)
        @test device_model_matrix_response_metadata_signature(
            result.metadata_before,
        ) == device_model_matrix_response_metadata_signature(
            result.metadata_after,
        )

        supports_mtf =
            device_model_matrix_detector_response_symbol(row) != :none
        @test supports_detector_mtf(detector) == supports_mtf
        frequency_x = 0.25
        frequency_y = 0.125
        @test detector_mtf(detector, frequency_x, frequency_y) ≈
            device_model_matrix_expected_detector_mtf(
                detector,
                frequency_x,
                frequency_y,
                Float64,
            ) atol=2e-14 rtol=2e-14

        if row isa DeviceModelMatrixM2FrameTransferEMCCD
            @test result.status_trace == (
                (true, true, false),
                (false, true, false),
                (false, true, false),
                (true, false, true),
                (true, true, false),
            )
            @test result.event_times == (
                start=PlantTimestamp(100_000_000),
                close=PlantTimestamp(1_100_000_000),
                transfer=PlantTimestamp(1_200_000_000),
                readout=PlantTimestamp(1_400_000_000),
            )
            @test frame_transfer_product_sequence(result.state) == UInt64(1)
            @test acquisition_product_ready_timestamp(result.state) ==
                result.event_times.readout
            @test result.prepared.storage_frame !== detector.state.frame
            @test !Base.mightalias(
                result.prepared.storage_frame,
                detector.state.frame,
            )
        else
            @test result.status_trace == (
                DetectorAcquisitionReady,
                DetectorExposureActive,
                DetectorReadoutPending,
                DetectorReadoutComplete,
                DetectorAcquisitionReady,
            )
            @test detector_acquisition_status(result.state) ==
                DetectorAcquisitionReady
        end

        if row isa DeviceModelMatrixM5RollingCMOS
            @test rolling_band_count(result.prepared) == 3
            @test rolling_opened_band_count(result.state) == 3
            @test rolling_closed_band_count(result.state) == 3
            @test ntuple(
                index -> rolling_band_open_timestamp(
                    result.prepared,
                    result.state,
                    index,
                ),
                3,
            ) == (
                PlantTimestamp(100_000_000),
                PlantTimestamp(200_000_000),
                PlantTimestamp(300_000_000),
            )
            @test ntuple(
                index -> rolling_band_close_timestamp(
                    result.prepared,
                    result.state,
                    index,
                ),
                3,
            ) == (
                PlantTimestamp(1_100_000_000),
                PlantTimestamp(1_200_000_000),
                PlantTimestamp(1_300_000_000),
            )
            @test result.event_times == (
                start=PlantTimestamp(100_000_000),
                readout=PlantTimestamp(1_400_000_000),
                readiness=PlantTimestamp(1_500_000_000),
            )
        elseif row isa DeviceModelMatrixM6UpTheRampHgCdTe
            response = expected
            expected_cube = cat(
                zeros(9, 9),
                response .* 0.25,
                response;
                dims=3,
            )
            @test nondestructive_read_count(result.prepared) == 3
            @test ntuple(
                index -> nondestructive_read_offset(
                    result.prepared,
                    index,
                ),
                3,
            ) == (
                PlantDuration(0),
                PlantDuration(500_000_000),
                PlantDuration(1_000_000_000),
            )
            @test Array(detector_ramp_cube(detector)) ≈ expected_cube
            @test Array(detector_ramp_slope(detector)) ≈ response
            @test Array(detector_ramp_intercept(detector)) ≈
                -response ./ 12
            @test detector_ramp_times(detector) == [0.0, 0.5, 1.0]
            @test result.event_times == (
                start=PlantTimestamp(100_000_000),
                middle=PlantTimestamp(600_000_000),
                close=PlantTimestamp(1_100_000_000),
                readout=PlantTimestamp(1_300_000_000),
                readiness=PlantTimestamp(1_400_000_000),
            )
        elseif !(row isa DeviceModelMatrixM2FrameTransferEMCCD)
            @test exposure_start_timestamp(result.state) ==
                PlantTimestamp(100_000_000)
            @test integrated_through_timestamp(result.state) ==
                PlantTimestamp(1_100_000_000)
            @test result.event_times == (
                start=PlantTimestamp(100_000_000),
                close=PlantTimestamp(1_100_000_000),
                readout=PlantTimestamp(1_300_000_000),
                readiness=PlantTimestamp(1_400_000_000),
            )
        end
    end
end

struct PlantResourceGraphFakeBackend <: AbstractArrayBackend end

@inline function plant_resource_graph_storage_bytes(storage)
    return UInt64(length(storage)) * UInt64(Base.elsize(typeof(storage)))
end

@inline plant_resource_graph_optional_bytes(::Nothing) = UInt64(0)
@inline plant_resource_graph_optional_bytes(storage) =
    plant_resource_graph_storage_bytes(storage)

function plant_resource_graph_storage_bytes(storages::Tuple)
    bytes = UInt64(0)
    for storage in storages
        bytes += plant_resource_graph_optional_bytes(storage)
    end
    return bytes
end

@inline plant_resource_graph_fact_bytes(fact) = (
    resident=structural_resident_bytes(fact),
    workspace=structural_workspace_bytes(fact),
)

function plant_resource_graph_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

@testset "Whole-graph selected-partition ownership" begin
    target = HostComputeDevice()
    id = StructuralResourceOwnerID(:path, :off_partition)
    facts = AbstractStructuralResourceFact[]
    Plant._append_structural_graph_fact!(facts,
        UnknownStructuralResourceFact(id, target, :owner_not_on_device),
        target)
    @test length(facts) == 1
    @test structural_resource_known(only(facts))
    @test structural_resource_owner_id(only(facts)) == id
    @test structural_resident_bytes(only(facts)) == UInt64(0)
    @test structural_workspace_bytes(only(facts)) == UInt64(0)

    unsupported_id = StructuralResourceOwnerID(:path, :unsupported)
    Plant._append_structural_graph_fact!(facts,
        UnknownStructuralResourceFact(
            unsupported_id, target, :unsupported_prepared_path),
        target)
    @test !structural_resource_known(last(facts))
    @test structural_resource_unknown_reason(last(facts)) ==
        :unsupported_prepared_path

    invalid_error = try
        Plant._append_structural_graph_fact!(facts, 1, target)
        nothing
    catch caught
        caught
    end
    @test invalid_error isa StructuralResourceError
    @test invalid_error.reason == :invalid_fact
end

function plant_resource_graph_atmosphere_bytes(atmosphere)
    resident = UInt64(0)
    workspace = UInt64(0)
    for layer in atmosphere.layers
        telescope = layer.generator_telescope
        state = layer.generator.state
        resident += plant_resource_graph_storage_bytes((
            telescope.aperture.pupil,
            telescope.aperture.reflectivity,
            state.opd,
            state.psd,
            state.freqs,
        ))
        workspace += plant_resource_graph_storage_bytes((
            state.spectrum,
            state.noise_re,
            state.noise_im,
            state.noise_re_host,
            state.noise_im_host,
        ))
    end
    return (; resident, workspace)
end

function plant_resource_graph_command_bytes(
    endpoint,
    state,
    application,
    workspace,
)
    binding = endpoint.binding
    payloads = state.payloads
    payload_values = payloads.values
    capacity = length(payload_values) - 1

    resident = plant_resource_graph_storage_bytes((
        binding.initial_command,
        binding.safe_command,
        state.slots,
        state.calendar,
        state.accepted_sequences,
    ))
    resident += UInt64(capacity) * UInt64(Base.elsize(typeof(payload_values)))
    for index in 1:capacity
        resident += plant_resource_graph_storage_bytes(payload_values[index])
    end
    resident += plant_resource_graph_storage_bytes((
        application.values.effective,
        application.values.safe,
    ))

    workspace_bytes = UInt64(Base.elsize(typeof(payload_values)))
    workspace_bytes += plant_resource_graph_storage_bytes(
        payload_values[payloads.staging_slot])
    workspace_bytes += plant_resource_graph_storage_bytes(
        application.values.staging)
    workspace_bytes += plant_resource_graph_storage_bytes(
        workspace.dispositions)
    return (; resident, workspace=workspace_bytes)
end

function plant_resource_graph_optic_bytes(
    prepared,
    state,
    workspace,
)
    implementation = prepared.implementation
    topology = implementation.params.topology
    modes = implementation.modes
    resident = plant_resource_graph_storage_bytes((
        modes.pupil_backend,
        modes.coordinates_backend,
        implementation.separable_x,
        implementation.separable_y_t,
        topology.coords,
        topology.active_coords,
        topology.valid_actuators,
        topology.active_indices,
        state.active.state.opd,
        state.active.state.coefs,
        state.active.state.actuator_coefs,
        state.active.state.separable_tmp,
    ))
    workspace_bytes = plant_resource_graph_storage_bytes((
        workspace.staged.state.opd,
        workspace.staged.state.coefs,
        workspace.staged.state.actuator_coefs,
        workspace.staged.state.separable_tmp,
    ))
    return (; resident, workspace=workspace_bytes)
end

function plant_resource_graph_path_bytes(path)
    execution = path.execution
    renderer = path.materialization.renderer
    resident = plant_resource_graph_storage_bytes((
        path.input.support,
        path.input.amplitude,
        path.input.opd,
        execution.field.values,
        execution.output.values,
        renderer.shift_x,
        renderer.shift_y,
        renderer.footprint_scale,
        renderer.pupil,
    ))
    workspace = plant_resource_graph_storage_bytes((
        execution.plan.propagation.state.scratch,
        execution.plan.unshifted_intensity,
    ))
    return (; resident, workspace)
end

function plant_resource_graph_acquisition_bytes(acquisition)
    lifecycle = acquisition.lifecycle
    detector = lifecycle.detector
    state = detector.state
    resident = plant_resource_graph_storage_bytes((
        state.frame,
        state.accum_buffer,
        state.latent_buffer,
        state.output_buffer,
        lifecycle.read_offsets,
        acquisition.products.observation,
    ))
    workspace = plant_resource_graph_storage_bytes((
        state.presampling_buffer,
        state.presampling_scratch,
        state.response_buffer,
        state.bin_buffer,
        state.temporal_buffer,
        state.noise_buffer,
        state.noise_buffer_host,
        state.batched_buffer_host,
        state.output_buffer_host,
    ))
    return (; resident, workspace)
end

function plant_resource_graph_fact(report, id)
    return only(fact for fact in structural_resource_facts(report) if
        structural_resource_owner_id(fact) == id)
end

@testset "Whole event-loop structural resource graph" begin
    target = HostComputeDevice()
    fixture = device_batch_test_fixture(
        selection=Val(:none),
        include_beta=false,
        include_unequal_rate=false,
        include_phase_offset=false,
        include_origin_offset=false,
        include_lgs=false,
    )
    prepared = fixture.prepared
    state = fixture.state
    workspace = fixture.workspace
    path = only(prepared.path_groups).path
    telescope = path.telescope

    telescope_bytes = (
        resident=plant_resource_graph_storage_bytes((
            telescope.aperture.pupil,
            telescope.aperture.reflectivity,
        )),
        workspace=UInt64(0),
    )
    atmosphere_bytes =
        plant_resource_graph_atmosphere_bytes(prepared.atmosphere)
    endpoint_bytes = plant_resource_graph_command_bytes(
        only(prepared.command_endpoints),
        only(state.command_endpoints),
        only(state.command_applications),
        only(workspace.command_endpoints),
    )
    optic_bytes = plant_resource_graph_optic_bytes(
        only(prepared.optics),
        only(state.controllable_optics),
        only(workspace.controllable_optics),
    )
    path_bytes = plant_resource_graph_path_bytes(path)
    acquisition_bytes = plant_resource_graph_acquisition_bytes(
        only(prepared.acquisitions))
    aberration_bytes = Tuple((
        resident=plant_resource_graph_storage_bytes(aberration.opd),
        workspace=UInt64(0),
    ) for aberration in prepared.sampled_aberrations)

    expected = (
        (id=StructuralResourceOwnerID(:acquisition, :alpha_camera),
            bytes=acquisition_bytes),
        (id=StructuralResourceOwnerID(:atmosphere, :primary),
            bytes=atmosphere_bytes),
        (id=StructuralResourceOwnerID(
                :command_endpoint, :device_batch_dm_command),
            bytes=endpoint_bytes),
        (id=StructuralResourceOwnerID(
                :controllable_optic, :device_batch_dm),
            bytes=optic_bytes),
        (id=StructuralResourceOwnerID(:path, :alpha), bytes=path_bytes),
        (id=StructuralResourceOwnerID(
                :sampled_aberration, :device_batch_alpha_ncpa),
            bytes=aberration_bytes[2]),
        (id=StructuralResourceOwnerID(
                :sampled_aberration, :device_batch_static),
            bytes=aberration_bytes[1]),
        (id=StructuralResourceOwnerID(:telescope, :primary),
            bytes=telescope_bytes),
    )

    report = require_exact_structural_resource_facts(
        prepared, state, workspace, target)
    facts = structural_resource_facts(report)
    @test Tuple(structural_resource_owner_id(fact) for fact in facts) ==
        Tuple(entry.id for entry in expected)
    @test all(structural_resource_known, facts)
    for entry in expected
        @test plant_resource_graph_fact_bytes(
            plant_resource_graph_fact(report, entry.id)) == entry.bytes
    end
    expected_resident = sum(
        entry -> entry.bytes.resident, expected; init=UInt64(0))
    expected_workspace = sum(
        entry -> entry.bytes.workspace, expected; init=UInt64(0))
    @test structural_resident_bytes(report) == expected_resident
    @test structural_workspace_bytes(report) == expected_workspace
    @test opaque_resource_reserve_bytes(report) == UInt64(0)
    @test telescope_bytes.resident == UInt64(576)
    @test sum(entry -> entry.resident, aberration_bytes;
        init=UInt64(0)) == UInt64(1_024)

    repeated = require_exact_structural_resource_facts(
        prepared, state, workspace, target)
    @test structural_resource_facts(repeated) == facts
    @test structural_resident_bytes(repeated) == expected_resident
    @test structural_workspace_bytes(repeated) == expected_workspace

    foreign = device_batch_test_fixture(
        selection=Val(:none),
        include_beta=false,
        include_unequal_rate=false,
        include_phase_offset=false,
        include_origin_offset=false,
        include_lgs=false,
    )
    foreign_state_error = plant_resource_graph_error() do
        require_exact_structural_resource_facts(
            prepared, foreign.state, workspace, target)
    end
    @test typeof(foreign_state_error) === PlantScheduleError
    @test foreign_state_error.reason == :foreign_state
    foreign_workspace_error = plant_resource_graph_error() do
        require_exact_structural_resource_facts(
            prepared, state, foreign.workspace, target)
    end
    @test typeof(foreign_workspace_error) === PlantScheduleError
    @test foreign_workspace_error.reason == :foreign_workspace

    wrong_target = AcceleratorComputeDevice(
        PlantResourceGraphFakeBackend(), UInt32(7))
    wrong_target_error = plant_resource_graph_error() do
        require_exact_structural_resource_facts(
            prepared, state, workspace, wrong_target)
    end
    @test typeof(wrong_target_error) === StructuralResourceError
    @test wrong_target_error.reason == :wrong_device

    reserve = OpaqueResourceReserve(
        StructuralResourceOwnerID(:path, :alpha),
        target,
        :test_fft_provider,
        4_096,
        ResourceEstimateMethod(:declared_test_fft_reserve, 1),
    )
    reserved_report = require_exact_structural_resource_facts(
        prepared, state, workspace, target; opaque_reserves=(reserve,))
    @test structural_resident_bytes(reserved_report) == expected_resident
    @test structural_workspace_bytes(reserved_report) == expected_workspace
    @test opaque_resource_reserve_bytes(reserved_report) == UInt64(4_096)
    @test only(opaque_resource_reserves(reserved_report)) == reserve
    @test structural_resource_owner_id(
        only(opaque_resource_reserves(reserved_report))) ==
        StructuralResourceOwnerID(:path, :alpha)

    unknown_owner_reserve = OpaqueResourceReserve(
        StructuralResourceOwnerID(:path, :not_prepared),
        target,
        :test_fft_provider,
        1,
        ResourceEstimateMethod(:declared_test_fft_reserve, 1),
    )
    reserve_error = plant_resource_graph_error() do
        require_exact_structural_resource_facts(
            prepared, state, workspace, target;
            opaque_reserves=(unknown_owner_reserve,))
    end
    @test typeof(reserve_error) === StructuralResourceError
    @test reserve_error.reason == :unknown_owner
end

@testset "Direct device batch workspace structural ownership" begin
    fixture = device_batch_test_fixture(
        selection=Val(:all),
        include_unequal_rate=false,
        include_phase_offset=false,
        include_origin_offset=false,
        include_lgs=false,
        include_physical_state=false,
    )
    prepared = fixture.prepared
    target = HostComputeDevice()
    owner = only(prepared.device_path_batch_owners)
    implementation = owner.implementation
    atmosphere_workspace = implementation.atmosphere_batch.workspace
    optical_workspace = implementation.optical_batch.workspace
    independent_storage = (
        atmosphere_workspace.shift_x,
        atmosphere_workspace.shift_y,
        atmosphere_workspace.footprint_scale,
        atmosphere_workspace.pupil,
        atmosphere_workspace.output,
        optical_workspace.field_stack,
        optical_workspace.output_stack,
        optical_workspace.shift_axis1,
        optical_workspace.shift_axis2,
    )
    expected_workspace =
        plant_resource_graph_storage_bytes(independent_storage)

    report = require_exact_structural_resource_facts(
        prepared, fixture.state, fixture.workspace, target)
    batch_id = StructuralResourceOwnerID(:direct_batch_workspace, :alpha)
    batch_fact = plant_resource_graph_fact(report, batch_id)
    @test structural_resource_known(batch_fact)
    @test structural_resident_bytes(batch_fact) == UInt64(0)
    @test structural_workspace_bytes(batch_fact) == expected_workspace

    count = length(owner.group_slots)
    @test count == 2
    for member in 1:count
        @test parent(implementation.atmosphere_outputs[member]) ===
            atmosphere_workspace.output
        @test parent(implementation.optical_batch.fields[member].values) ===
            optical_workspace.field_stack
        @test parent(implementation.optical_batch.output[member].values) ===
            optical_workspace.output_stack
        slot = Int(owner.group_slots[member])
        path = prepared.path_groups[slot].path
        @test implementation.path_inputs[member] === path.input
        @test implementation.path_results[member] === path.result
    end
    aliased_view_bytes = plant_resource_graph_storage_bytes(Tuple(
        implementation.atmosphere_outputs,
    )) + plant_resource_graph_storage_bytes(Tuple(
        implementation.optical_batch.fields[index].values for
        index in 1:count
    )) + plant_resource_graph_storage_bytes(Tuple(
        implementation.optical_batch.output[index].values for
        index in 1:count
    ))
    @test aliased_view_bytes > UInt64(0)
    @test structural_workspace_bytes(batch_fact) !=
        expected_workspace + aliased_view_bytes
end

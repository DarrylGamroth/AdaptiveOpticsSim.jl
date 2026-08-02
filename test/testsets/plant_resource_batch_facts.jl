struct WFSBatchResourceFactFakeBackend <: AbstractArrayBackend end

@testset "WFS device path-batch structural resource facts" begin
    fixture = device_model_matrix_wfs_fixture(
        DeviceModelMatrixShackHartmann();
        selection=Val(:all),
        direction=Val(:ngs),
        spectral=Val(:monochromatic),
    )
    @test device_path_batch_owner_count(fixture.prepared) == 1
    owner = device_path_batch_owner(fixture.prepared, 1)
    @test typeof(owner.implementation) <:
        Plant._PreparedWFSDevicePathBatch

    target = HostComputeDevice()
    id = StructuralResourceOwnerID(:workspace, :sh_device_path_batch)
    fact = structural_resource_fact(
        owner, fixture.prepared.atmosphere, id, target)
    atmosphere_workspace = owner.implementation.atmosphere_batch.workspace
    independent_arrays = (
        atmosphere_workspace.shift_x,
        atmosphere_workspace.shift_y,
        atmosphere_workspace.footprint_scale,
        atmosphere_workspace.pupil,
        atmosphere_workspace.output,
    )
    expected_workspace = sum(independent_arrays; init=UInt64(0)) do array
        UInt64(length(array)) * UInt64(sizeof(eltype(array)))
    end

    @test structural_resource_known(fact)
    @test structural_resource_owner_id(fact) == id
    @test compute_device(fact) == target
    @test structural_resident_bytes(fact) == UInt64(0)
    @test structural_workspace_bytes(fact) == expected_workspace
    @test all(array -> compute_device(array) == target, independent_arrays)

    graph_report = require_exact_structural_resource_facts(
        fixture.prepared, fixture.state, fixture.workspace, target)
    graph_id = StructuralResourceOwnerID(
        :wfs_batch_workspace, :wfs_alpha)
    graph_fact = only(candidate for candidate in
        structural_resource_facts(graph_report) if
        structural_resource_owner_id(candidate) == graph_id)
    @test structural_resident_bytes(graph_fact) == UInt64(0)
    @test structural_workspace_bytes(graph_fact) == expected_workspace

    off_partition = AcceleratorComputeDevice(
        WFSBatchResourceFactFakeBackend(), UInt32(1))
    off_partition_fact = structural_resource_fact(
        owner, fixture.prepared.atmosphere, id, off_partition)
    @test !structural_resource_known(off_partition_fact)
    @test structural_resource_owner_id(off_partition_fact) == id
    @test compute_device(off_partition_fact) == off_partition
    @test structural_resource_unknown_reason(off_partition_fact) ==
        :owner_not_on_device

    unsupported = Plant._device_path_batch_structural_resource_fact(
        nothing,
        StructuralResourceOwnerID(:workspace, :unsupported_batch),
        target,
    )
    @test !structural_resource_known(unsupported)
    @test structural_resource_unknown_reason(unsupported) ==
        :unsupported_device_path_batch
end

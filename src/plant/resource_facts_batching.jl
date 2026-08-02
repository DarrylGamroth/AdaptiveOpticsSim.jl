#
# Gate 9A sampled-aberration and direct device path-batch resource facts.
#
# Prepared sampled aberrations own only their copied OPD storage. A direct
# device path-batch owner owns the independent atmosphere and optical batch
# stacks below. Slice views, path inputs/results, binding registries, and FFT
# plans are borrowed or opaque and are not structural storage of these owners.
#

"""Report the target-local OPD owned by one prepared sampled aberration."""
function structural_resource_fact(aberration::PreparedSampledAberration,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    _require_exact_sampled_aberration_target(
        aberration, aberration.metadata.device)
    resident = _structural_array_target_bytes(
        (aberration.opd,), target, :resident_bytes)
    return _targeted_structural_resource_fact(id, target, resident,
        (present=false, bytes=UInt64(0)))
end

@inline function _atmosphere_batch_structural_workspace_bytes(
    implementation::_PreparedDevicePathBatchImplementation,
    target::AbstractComputeDevice,
)
    workspace = implementation.atmosphere_batch.workspace
    return _structural_array_target_bytes((
        workspace.shift_x,
        workspace.shift_y,
        workspace.footprint_scale,
        workspace.pupil,
        workspace.output,
    ), target, :workspace_bytes)
end

function _device_path_batch_structural_resource_fact(
    ::Any,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    return UnknownStructuralResourceFact(
        id, target, :unsupported_device_path_batch)
end

function _device_path_batch_structural_resource_fact(
    implementation::_PreparedDirectImagingDevicePathBatch,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    optical_workspace = implementation.optical_batch.workspace
    workspace = _atmosphere_batch_structural_workspace_bytes(
        implementation, target)
    optical_bytes = _structural_array_target_bytes((
        optical_workspace.field_stack,
        optical_workspace.output_stack,
        optical_workspace.shift_axis1,
        optical_workspace.shift_axis2,
    ), target, :workspace_bytes)
    workspace = _combine_structural_target_bytes(
        workspace, optical_bytes, :workspace_bytes)
    return _targeted_structural_resource_fact(
        id, target, (present=false, bytes=UInt64(0)), workspace)
end

function _device_path_batch_structural_resource_fact(
    implementation::_PreparedWFSDevicePathBatch,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    workspace = _atmosphere_batch_structural_workspace_bytes(
        implementation, target)
    return _targeted_structural_resource_fact(
        id, target, (present=false, bytes=UInt64(0)), workspace)
end

"""
Report storage owned by one prepared device path-batch owner against its
authoritative atmosphere. Unsupported validated batch implementations return
an explicit unknown fact.
"""
function structural_resource_fact(
    owner::PreparedDevicePathBatchOwner,
    atmosphere::AbstractTimedAtmosphere,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    _require_exact_device_path_batch_owner_target(
        owner, atmosphere, owner.device)
    return _device_path_batch_structural_resource_fact(
        owner.implementation, id, target)
end

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

function _device_path_batch_structural_resource_fact(
    ::_PreparedDevicePathBatchImplementation,
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
    atmosphere_workspace = implementation.atmosphere_batch.workspace
    optical_workspace = implementation.optical_batch.workspace
    workspace = _structural_array_target_bytes((
        atmosphere_workspace.shift_x,
        atmosphere_workspace.shift_y,
        atmosphere_workspace.footprint_scale,
        atmosphere_workspace.pupil,
        atmosphere_workspace.output,
        optical_workspace.field_stack,
        optical_workspace.output_stack,
        optical_workspace.shift_axis1,
        optical_workspace.shift_axis2,
    ), target, :workspace_bytes)
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

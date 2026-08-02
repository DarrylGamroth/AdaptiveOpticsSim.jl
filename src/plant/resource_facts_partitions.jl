#
# Exact structural facts for preparation-only target-local resources
#
# A target-local path has no atmosphere renderer. Its resource fact therefore
# counts only the path input, optical products, and propagation workspace. A
# target-local acquisition counts its provider-owned detector/state and public
# products, but no event lifecycle or scheduler storage.
#

function _target_local_path_resource_fact(
    input,
    result,
    execution,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    return structural_resource_fact(execution, id, target)
end

function _target_local_path_resource_fact(
    input::Union{PupilFunction,ElectricField},
    result::IntensityMap,
    execution::PreparedDirectImaging,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    validate_path_execution_binding(execution, input, result)
    return structural_resource_fact(execution, id, target)
end

function _target_local_path_resource_fact(
    input::PupilFunction,
    result::IntensityMap,
    execution::WFSOpticalPathExecution{
        <:PreparedShackHartmannOpticalFormation},
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    plan = execution.plan
    validate_path_execution_binding(execution, input, result)
    resident = _combine_structural_target_bytes(
        _pupil_structural_resident_bytes(input, target),
        _structural_array_target_bytes(
            (result.values,), target, :resident_bytes),
        :resident_bytes,
    )
    path_fact = _targeted_structural_resource_fact(
        id,
        target,
        resident,
        (present=false, bytes=UInt64(0)),
    )
    layout_fact = structural_resource_fact(plan.front_end.layout, id, target)
    propagation_fact = structural_resource_fact(
        plan.front_end.propagation, id, target)
    return _combine_structural_owner_facts(
        id, target, (path_fact, layout_fact, propagation_fact))
end

"""Report one non-executable target-local optical resource owner."""
function structural_resource_fact(
    resources::PreparedTargetLocalPathResources,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    execution_target = getfield(resources.key, :device)
    target == execution_target || return UnknownStructuralResourceFact(
        id, target, :owner_not_on_device)
    _require_prepared_target_local_path_resources(
        resources,
        resources.definition,
        resources.source,
        resources.telescope,
        resources.context,
    )
    return _target_local_path_resource_fact(
        resources.input, resources.result, resources.execution, id, target)
end

function _target_local_acquisition_resource_fact(
    implementation,
    products::AcquisitionProducts,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    return _combine_structural_owner_facts(id, target, (
        structural_resource_fact(implementation, id, target),
        structural_resource_fact(products, id, target),
    ))
end

function structural_resource_fact(
    ::PreparedUnchangedSyntheticProvider,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    return KnownStructuralResourceFact(
        id, target, UInt64(0), UInt64(0))
end

function _target_local_acquisition_resource_fact(
    implementation::PreparedFullOpticalProvider{
        <:FrameAcquisitionExecution},
    products::AcquisitionProducts,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    execution = implementation.execution
    return _combine_structural_owner_facts(id, target, (
        structural_resource_fact(execution.detector, id, target),
        structural_resource_fact(products, id, target),
    ))
end

"""Report one non-executable target-local acquisition resource owner."""
function structural_resource_fact(
    resources::PreparedTargetLocalAcquisitionResources,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    execution_target = getfield(resources.path_key, :device)
    target == execution_target || return UnknownStructuralResourceFact(
        id, target, :owner_not_on_device)
    provider = resources.provider
    _prepared_device_execution_compute_device(resources.context) == target ||
        _structural_resource_error(:acquisition, :wrong_device,
            "target-local acquisition context does not match its resource fact")
    _require_exact_plant_product_target(
        resources.path_result, target, "target-local acquisition path result")
    validate_acquisition_provider_binding(provider, resources.path_result)
    _require_exact_acquisition_products_target(provider.products, target)
    validate_acquisition_provider_target(
        provider.implementation,
        resources.path_result,
        provider.products,
        target,
    )
    return _target_local_acquisition_resource_fact(
        provider.implementation, provider.products, id, target)
end

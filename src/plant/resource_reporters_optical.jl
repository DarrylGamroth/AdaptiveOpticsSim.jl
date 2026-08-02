#
# Dispatch-owned optical structural-resource reporters
#
# Each method below names its arrays directly.  It deliberately does not walk
# object fields: product and workspace aliases are assigned to one explicit
# owner, while FFT plans and provider runtime memory remain opaque.
#

@inline function _owned_structural_array_bytes(
    target::AbstractComputeDevice,
    arrays::AbstractArray...,
)
    total = UInt64(0)
    for array in arrays
        total = _checked_resource_add(total,
            structural_array_bytes(array, target), :resident_bytes)
    end
    return total
end

@inline function _owned_structural_workspace_bytes(
    target::AbstractComputeDevice,
    arrays::AbstractArray...,
)
    total = UInt64(0)
    for array in arrays
        total = _checked_resource_add(total,
            structural_array_bytes(array, target), :workspace_bytes)
    end
    return total
end

@inline function _require_distinct_structural_storage(left::AbstractArray,
    right::AbstractArray, owner::Symbol)
    left === right && _structural_resource_error(owner, :unexpected_alias,
        "explicit structural storage owners must not alias")
    return nothing
end

@inline function _telescope_structural_resident_bytes(telescope::Telescope,
    target::AbstractComputeDevice)
    return _owned_structural_array_bytes(target, telescope.aperture.pupil,
        telescope.aperture.reflectivity)
end

"""Report the explicitly owned aperture mask and reflectivity arrays."""
function structural_resource_fact(telescope::Telescope,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    return KnownStructuralResourceFact(id, target,
        _telescope_structural_resident_bytes(telescope, target), 0)
end

@inline function _pupil_structural_resident_bytes(pupil::PupilFunction,
    target::AbstractComputeDevice)
    return _owned_structural_array_bytes(target, pupil.support,
        pupil.amplitude, pupil.opd)
end

"""Report one path-local pupil product; metadata is scalar-only."""
function structural_resource_fact(pupil::PupilFunction,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    return KnownStructuralResourceFact(id, target,
        _pupil_structural_resident_bytes(pupil, target), 0)
end

@inline function _kolmogorov_structural_resident_bytes(
    atmosphere::KolmogorovAtmosphere,
    target::AbstractComputeDevice,
)
    state = atmosphere.state
    return _owned_structural_array_bytes(target, state.opd, state.psd,
        state.spectrum, state.noise_re, state.noise_im, state.freqs)
end

@inline function _kolmogorov_structural_resident_bytes(
    atmosphere::KolmogorovAtmosphere,
    target::HostComputeDevice,
)
    state = atmosphere.state
    # These staging matrices reside on the CPU graph.  The accelerator method
    # above intentionally excludes them rather than claiming host staging as
    # accelerator storage.
    return _owned_structural_array_bytes(target, state.opd, state.psd,
        state.spectrum, state.noise_re, state.noise_im, state.noise_re_host,
        state.noise_im_host, state.freqs)
end

"""Report a standalone timed Kolmogorov atmosphere, excluding its FFT plan."""
function structural_resource_fact(atmosphere::KolmogorovAtmosphere,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    return KnownStructuralResourceFact(id, target,
        _kolmogorov_structural_resident_bytes(atmosphere, target), 0)
end

"""
Report the finite multilayer atmosphere as one owner.  Its fixed host layer
registry and cold parameter vectors are excluded; each layer's prepared screen
generator and its private screen telescope belong to this atmosphere owner.
"""
function structural_resource_fact(atmosphere::MultiLayerAtmosphere,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    resident = UInt64(0)
    for layer in atmosphere.layers
        resident = _checked_resource_add(resident,
            _telescope_structural_resident_bytes(layer.generator_telescope,
                target), :resident_bytes)
        resident = _checked_resource_add(resident,
            _kolmogorov_structural_resident_bytes(layer.generator, target),
            :resident_bytes)
    end
    return KnownStructuralResourceFact(id, target, resident, 0)
end

@inline function _renderer_structural_resident_bytes(
    renderer::AtmosphereDirectionRenderer,
    target::AbstractComputeDevice,
)
    return _owned_structural_array_bytes(target, renderer.shift_x,
        renderer.shift_y, renderer.footprint_scale, renderer.pupil)
end

function _materialization_structural_resource_fact(
    ::PreparedPupilOPDMaterialization,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    return UnknownStructuralResourceFact(id, target, :unsupported_renderer)
end

function _materialization_structural_resource_fact(
    materialization::PreparedPupilOPDMaterialization{
        <:AtmosphereDirectionRenderer},
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    # The destination is the path's PupilFunction owner and is not counted
    # here.  The frozen renderer arrays are owned by materialization.
    return KnownStructuralResourceFact(id, target,
        _renderer_structural_resident_bytes(materialization.renderer, target),
        0)
end

"""Report one prepared atmosphere-renderer/materialization owner."""
function structural_resource_fact(
    materialization::PreparedPupilOPDMaterialization,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    return _materialization_structural_resource_fact(materialization, id,
        target)
end

@inline function _direct_imaging_resident_bytes(input::PupilFunction,
    field::ElectricField, output::IntensityMap,
    target::AbstractComputeDevice)
    input_bytes = _pupil_structural_resident_bytes(input, target)
    return _checked_resource_add(input_bytes,
        _direct_imaging_noninput_resident_bytes(field, output, target),
        :resident_bytes)
end

function _direct_imaging_resident_bytes(input::ElectricField,
    field::ElectricField, output::IntensityMap,
    target::AbstractComputeDevice)
    input === field || _structural_resource_error(:direct_imaging,
        :invalid_binding,
        "preformed direct-imaging input must be its prepared field")
    return _direct_imaging_noninput_resident_bytes(field, output, target)
end

@inline function _direct_imaging_noninput_resident_bytes(field::ElectricField,
    output::IntensityMap, target::AbstractComputeDevice)
    return _owned_structural_array_bytes(target, field.values, output.values)
end

@inline function _direct_imaging_workspace_bytes(
    propagation::FraunhoferPropagation,
    unshifted_intensity::AbstractMatrix,
    target::AbstractComputeDevice,
)
    _require_distinct_structural_storage(propagation.state.scratch,
        unshifted_intensity, :direct_imaging)
    return _owned_structural_workspace_bytes(target,
        propagation.state.scratch, unshifted_intensity)
end

function _direct_imaging_workspace_bytes(propagation,
    ::AbstractMatrix, ::AbstractComputeDevice)
    _structural_resource_error(:direct_imaging, :unsupported_propagation,
        "no exact structural reporter is defined for $(typeof(propagation))")
end

@inline function _validate_direct_imaging_resource_bindings(
    prepared::PreparedDirectImaging)
    prepared.plan.field_values === prepared.field.values ||
        _structural_resource_error(:direct_imaging, :invalid_binding,
            "prepared direct-imaging field storage changed")
    prepared.plan.output_values === prepared.output.values ||
        _structural_resource_error(:direct_imaging, :invalid_binding,
            "prepared direct-imaging output storage changed")
    prepared.workspace.propagation === prepared.plan.propagation ||
        _structural_resource_error(:direct_imaging, :invalid_binding,
            "prepared direct-imaging propagation changed")
    prepared.workspace.unshifted_intensity === prepared.plan.unshifted_intensity ||
        _structural_resource_error(:direct_imaging, :invalid_binding,
            "prepared direct-imaging scratch storage changed")
    _require_distinct_structural_storage(prepared.field.values,
        prepared.output.values, :direct_imaging)
    return nothing
end

"""
Report a prepared direct-imaging leaf.  Its plan only aliases the explicit
input/field/output products and workspace, so those arrays are counted once
from the leaf fields; the FFT plan itself is intentionally opaque.
"""
function structural_resource_fact(prepared::PreparedDirectImaging,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    _validate_direct_imaging_resource_bindings(prepared)
    resident = _direct_imaging_resident_bytes(prepared.input, prepared.field,
        prepared.output, target)
    workspace = _direct_imaging_workspace_bytes(prepared.plan.propagation,
        prepared.plan.unshifted_intensity, target)
    return KnownStructuralResourceFact(id, target, resident, workspace)
end

@inline function _validate_direct_imaging_batch_resource_bindings(
    prepared::PreparedDirectImagingBatch)
    workspace = prepared.workspace
    bindings = prepared.workspace_bindings
    workspace.field_stack === bindings.field_stack ||
        _structural_resource_error(:direct_imaging_batch, :invalid_binding,
            "prepared direct-imaging field stack changed")
    workspace.output_stack === bindings.output_stack ||
        _structural_resource_error(:direct_imaging_batch, :invalid_binding,
            "prepared direct-imaging output stack changed")
    workspace.shift_axis1 === bindings.shift_axis1 ||
        _structural_resource_error(:direct_imaging_batch, :invalid_binding,
            "prepared direct-imaging axis-1 shifts changed")
    workspace.shift_axis2 === bindings.shift_axis2 ||
        _structural_resource_error(:direct_imaging_batch, :invalid_binding,
            "prepared direct-imaging axis-2 shifts changed")
    workspace.fft_plan === bindings.fft_plan ||
        _structural_resource_error(:direct_imaging_batch, :invalid_binding,
            "prepared direct-imaging FFT plan changed")
    prepared.fields === bindings.fields ||
        _structural_resource_error(:direct_imaging_batch, :invalid_binding,
            "prepared direct-imaging field views changed")
    prepared.output === bindings.output ||
        _structural_resource_error(:direct_imaging_batch, :invalid_binding,
            "prepared direct-imaging output views changed")
    return nothing
end

"""
Report the native stacked direct-science workspace.  Field and product views
alias the two stacks, therefore the stacks are their sole structural owner.
The FFT plan/runtime reserve is excluded.
"""
function structural_resource_fact(prepared::PreparedDirectImagingBatch,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    _validate_direct_imaging_batch_resource_bindings(prepared)
    workspace = prepared.workspace
    bytes = _owned_structural_workspace_bytes(target, workspace.field_stack,
        workspace.output_stack, workspace.shift_axis1, workspace.shift_axis2)
    return KnownStructuralResourceFact(id, target, 0, bytes)
end

function _prepared_direct_path_resource_fact(
    ::Any,
    ::Any,
    ::Any,
    ::Any,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    return UnknownStructuralResourceFact(id, target,
        :unsupported_prepared_path)
end

function _prepared_direct_path_resource_fact(
    input::PupilFunction,
    result::IntensityMap,
    materialization::PreparedPupilOPDMaterialization,
    execution::PreparedDirectImaging,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    execution.input === input || _structural_resource_error(:path,
        :invalid_binding,
        "direct-imaging path input does not match its prepared execution")
    execution.output === result || _structural_resource_error(:path,
        :invalid_binding,
        "direct-imaging path result does not match its prepared execution")
    materialization.destination === input ||
        _structural_resource_error(:path, :invalid_binding,
            "path materialization destination does not match its input")
    _validate_direct_imaging_resource_bindings(execution)
    materialization = _materialization_structural_resource_fact(
        materialization, id, target)
    structural_resource_known(materialization) || return materialization
    resident = _pupil_structural_resident_bytes(input, target)
    resident = _checked_resource_add(resident,
        _direct_imaging_noninput_resident_bytes(execution.field,
            execution.output, target),
        :resident_bytes)
    resident = _checked_resource_add(resident,
        structural_resident_bytes(materialization), :resident_bytes)
    workspace = _direct_imaging_workspace_bytes(
        execution.plan.propagation,
        execution.plan.unshifted_intensity, target)
    return KnownStructuralResourceFact(id, target, resident, workspace)
end

"""
Report a direct-science prepared path as one owner.  The path owns its input,
result, direct-imaging products/workspace, and renderer geometry; internal
plan aliases and the renderer destination are counted only through that owner.
"""
function structural_resource_fact(path::PreparedPathExecutor,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    return _prepared_direct_path_resource_fact(path.input, path.result,
        path.materialization, path.execution, id, target)
end

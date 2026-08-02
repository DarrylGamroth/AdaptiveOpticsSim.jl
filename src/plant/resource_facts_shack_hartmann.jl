#
# Gate 9A deterministic Shack--Hartmann structural resource facts.
#
# Ownership table (all reports retain the caller-supplied owner ID):
#
# - PreparedShackHartmannOpticalFormation, PreparedWFSDetectorAcquisition,
#   DetectorAcquisitionPlan, PreparedShackHartmannEstimator, and
#   AcquisitionProducts only retain borrowed prepared bindings/products.
# - IntensityMap, WFSObservation, and WFSMeasurement own their nominated
#   caller-visible product storage.
# - Shack--Hartmann acquisition/estimator/calibration/layout states own the
#   arrays named below. Each exact-device fact selects only storage resident on
#   that device, including explicit host mirrors in a separate host fact.
# - PreparedMicrolensPropagation owns reusable numerical work arrays. FFT
#   plans remain opaque library storage and are never inferred here.
#

@inline _known_structural_resource_fact(id::StructuralResourceOwnerID,
    target::AbstractComputeDevice, resident::UInt64, workspace::UInt64) =
    KnownStructuralResourceFact(id, target, resident, workspace)

# These prepared wrappers carry exact bindings only. Their referenced
# front-end, detector, plan, input, output, observation, and measurement
# storage is reported by its concrete nominated owner below.
function structural_resource_fact(::PreparedShackHartmannOpticalFormation,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    return _known_structural_resource_fact(id, target, UInt64(0), UInt64(0))
end

function structural_resource_fact(::PreparedWFSDetectorAcquisition,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    return _known_structural_resource_fact(id, target, UInt64(0), UInt64(0))
end

function structural_resource_fact(::DetectorAcquisitionPlan,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    return _known_structural_resource_fact(id, target, UInt64(0), UInt64(0))
end

function structural_resource_fact(::PreparedShackHartmannEstimator,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    return _known_structural_resource_fact(id, target, UInt64(0), UInt64(0))
end

function structural_resource_fact(::AcquisitionProducts,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    return _known_structural_resource_fact(id, target, UInt64(0), UInt64(0))
end

function structural_resource_fact(product::IntensityMap,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    resident = _structural_array_target_bytes(
        (product.values,), target, :resident_bytes)
    return _targeted_structural_resource_fact(id, target, resident,
        (present=false, bytes=UInt64(0)))
end

function structural_resource_fact(observation::WFSObservation{<:AbstractArray},
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    resident = _structural_array_target_bytes(
        (observation.storage,), target, :resident_bytes)
    return _targeted_structural_resource_fact(id, target, resident,
        (present=false, bytes=UInt64(0)))
end

function structural_resource_fact(measurement::WFSMeasurement{<:AbstractArray},
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    resident = _structural_array_target_bytes(
        (measurement.storage,), target, :resident_bytes)
    return _targeted_structural_resource_fact(id, target, resident,
        (present=false, bytes=UInt64(0)))
end

function structural_resource_fact(state::ShackHartmannAcquisitionState,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    resident = _structural_array_target_bytes(
        (state.exported_spot_cube,), target, :resident_bytes)
    workspace = _structural_array_target_bytes(
        (state.spot_cube, state.detector_noise_cube), target,
        :workspace_bytes)
    return _targeted_structural_resource_fact(
        id, target, resident, workspace)
end

function structural_resource_fact(state::ShackHartmannEstimatorState,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    resident = _structural_array_target_bytes(
        (state.slopes,), target, :resident_bytes)
    workspace = _structural_array_target_bytes(
        (state.spot_stats, state.spot_stats_accum, state.slopes_host,
            state.centroid_host), target, :workspace_bytes)
    return _targeted_structural_resource_fact(
        id, target, resident, workspace)
end

function structural_resource_fact(calibration::SubapertureCalibration,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    resident = _structural_array_target_bytes(
        (calibration.reference_signal_2d,
            calibration.reference_signal_host), target, :resident_bytes)
    return _targeted_structural_resource_fact(id, target, resident,
        (present=false, bytes=UInt64(0)))
end

function structural_resource_fact(layout::SubapertureLayout,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    resident = _structural_array_target_bytes(
        (layout.valid_mask, layout.valid_mask_host,
            layout.valid_indices_host), target, :resident_bytes)
    return _targeted_structural_resource_fact(id, target, resident,
        (present=false, bytes=UInt64(0)))
end

function structural_resource_fact(state::DetectorState{
    <:Any,<:AbstractArray,Nothing,Nothing,NoFrameReadoutProducts},
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    resident = _structural_array_target_bytes((
        state.frame,
        state.accum_buffer,
        state.latent_buffer,
    ), target, :resident_bytes)
    workspace = _structural_array_target_bytes((
        state.presampling_buffer,
        state.presampling_scratch,
        state.response_buffer,
        state.bin_buffer,
        state.temporal_buffer,
        state.noise_buffer,
        state.noise_buffer_host,
        state.batched_buffer_host,
    ), target, :workspace_bytes)
    return _targeted_structural_resource_fact(
        id, target, resident, workspace)
end

function structural_resource_fact(workspace::PreparedMicrolensPropagation,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    workspace_bytes = _structural_array_target_bytes((
        workspace.field,
        workspace.phasor,
        workspace.fft_buffer,
        workspace.fft_stack,
        workspace.intensity,
        workspace.intensity_stack,
        workspace.intensity_tmp_stack,
        workspace.temp,
        workspace.bin_buffer,
        workspace.spot,
        workspace.sampled_spot_cube,
        workspace.spot_cube_accum,
        workspace.elongation_kernel,
        workspace.lgs_kernel_fft,
        workspace.fft_asterism_stack,
        workspace.amp_scales,
        workspace.amp_scales_host,
        workspace.opd_to_cycles,
        workspace.opd_to_cycles_host,
    ), target, :workspace_bytes)
    return _targeted_structural_resource_fact(id, target,
        (present=false, bytes=UInt64(0)), workspace_bytes)
end

function _prepared_path_resource_fact(
    input::PupilFunction,
    result::IntensityMap,
    materialization::PreparedPupilOPDMaterialization{
        <:AtmosphereDirectionRenderer},
    execution::WFSOpticalPathExecution{
        <:PreparedShackHartmannOpticalFormation},
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    plan = execution.plan
    plan.input === input || _structural_resource_error(
        :path, :invalid_binding,
        "Shack-Hartmann path input does not match its optical plan")
    plan.output === result || _structural_resource_error(
        :path, :invalid_binding,
        "Shack-Hartmann path result does not match its optical plan")
    materialization.destination === input || _structural_resource_error(
        :path, :invalid_binding,
        "Shack-Hartmann materialization destination does not match its input")

    resident = _combine_structural_target_bytes(
        _pupil_structural_resident_bytes(input, target),
        _structural_array_target_bytes(
            (result.values,), target, :resident_bytes),
        :resident_bytes)
    resident = _combine_structural_target_bytes(resident,
        _renderer_structural_resident_bytes(
            materialization.renderer, target), :resident_bytes)
    path_fact = _targeted_structural_resource_fact(
        id, target, resident,
        (present=false, bytes=UInt64(0)))
    layout_fact = structural_resource_fact(
        plan.front_end.layout, id, target)
    propagation_fact = structural_resource_fact(
        plan.front_end.propagation, id, target)
    return _combine_structural_owner_facts(
        id, target, (path_fact, layout_fact, propagation_fact))
end

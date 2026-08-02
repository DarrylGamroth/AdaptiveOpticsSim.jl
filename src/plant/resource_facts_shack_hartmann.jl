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
#   arrays named below; accelerator facts exclude their explicit host mirrors.
# - PreparedMicrolensPropagation owns reusable numerical work arrays. FFT
#   plans remain opaque library storage and are never inferred here.
#

@inline function _structural_array_total(target::AbstractComputeDevice,
    arrays::Tuple, component::Symbol, total::UInt64=UInt64(0))
    isempty(arrays) && return total
    bytes = structural_array_bytes(first(arrays), target)
    return _structural_array_total(target, Base.tail(arrays),
        component, _checked_resource_add(total, bytes, component))
end

@inline _structural_resident_array_total(target::AbstractComputeDevice,
    arrays...) = _structural_array_total(target, arrays, :resident_bytes)

@inline _structural_workspace_array_total(target::AbstractComputeDevice,
    arrays...) = _structural_array_total(target, arrays, :workspace_bytes)

@inline function _structural_host_staging_total(target::HostComputeDevice,
    arrays::Tuple, component::Symbol)
    return _structural_array_total(target, arrays, component)
end

@inline _structural_host_staging_total(::AcceleratorComputeDevice,
    ::Tuple, ::Symbol) = UInt64(0)

function _structural_host_staging_total(target::AbstractComputeDevice,
    ::Tuple, ::Symbol)
    _structural_resource_error(:array_storage, :unsupported_device,
        "structural host-staging accounting is not defined for $target")
end

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
    return _known_structural_resource_fact(id, target,
        structural_array_bytes(product.values, target), UInt64(0))
end

function structural_resource_fact(observation::WFSObservation{<:AbstractArray},
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    return _known_structural_resource_fact(id, target,
        structural_array_bytes(observation.storage, target), UInt64(0))
end

function structural_resource_fact(measurement::WFSMeasurement{<:AbstractArray},
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    return _known_structural_resource_fact(id, target,
        structural_array_bytes(measurement.storage, target), UInt64(0))
end

function structural_resource_fact(state::ShackHartmannAcquisitionState,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    resident = _structural_resident_array_total(target, state.spot_cube,
        state.exported_spot_cube, state.detector_noise_cube)
    return _known_structural_resource_fact(id, target, resident, UInt64(0))
end

function structural_resource_fact(state::ShackHartmannEstimatorState,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    resident = _structural_resident_array_total(target, state.slopes, state.spot_stats,
        state.spot_stats_accum)
    resident = _checked_resource_add(resident,
        _structural_host_staging_total(target,
            (state.slopes_host, state.centroid_host), :resident_bytes),
        :resident_bytes)
    return _known_structural_resource_fact(id, target, resident, UInt64(0))
end

function structural_resource_fact(calibration::SubapertureCalibration,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    resident = structural_array_bytes(calibration.reference_signal_2d, target)
    resident = _checked_resource_add(resident,
        _structural_host_staging_total(target,
            (calibration.reference_signal_host,), :resident_bytes),
        :resident_bytes)
    return _known_structural_resource_fact(id, target, resident, UInt64(0))
end

function structural_resource_fact(layout::SubapertureLayout,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    resident = structural_array_bytes(layout.valid_mask, target)
    resident = _checked_resource_add(resident,
        _structural_host_staging_total(target,
            (layout.valid_mask_host, layout.valid_indices_host),
            :resident_bytes),
        :resident_bytes)
    return _known_structural_resource_fact(id, target, resident, UInt64(0))
end

function structural_resource_fact(state::DetectorState{
    <:Any,<:AbstractArray,Nothing,Nothing,NoFrameReadoutProducts},
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    resident = _structural_resident_array_total(target,
        state.frame,
        state.presampling_buffer,
        state.presampling_scratch,
        state.response_buffer,
        state.bin_buffer,
        state.temporal_buffer,
        state.noise_buffer,
        state.accum_buffer,
        state.latent_buffer,
    )
    resident = _checked_resource_add(resident,
        _structural_host_staging_total(target,
            (state.noise_buffer_host, state.batched_buffer_host),
            :resident_bytes),
        :resident_bytes)
    return _known_structural_resource_fact(id, target, resident, UInt64(0))
end

function structural_resource_fact(workspace::PreparedMicrolensPropagation,
    id::StructuralResourceOwnerID, target::AbstractComputeDevice)
    workspace_bytes = _structural_workspace_array_total(target,
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
        workspace.opd_to_cycles,
    )
    workspace_bytes = _checked_resource_add(workspace_bytes,
        _structural_host_staging_total(target,
            (workspace.amp_scales_host, workspace.opd_to_cycles_host),
            :workspace_bytes),
        :workspace_bytes)
    return _known_structural_resource_fact(id, target, UInt64(0),
        workspace_bytes)
end

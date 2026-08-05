#
# Gate 9A acquisition and conventional frame-detector resource facts.
#
# Acquisition products own caller-visible storage. Detector state owns
# persistent accumulation, persistence, and thermal state; detector workspaces
# own scratch and host staging. Prepared lifecycles own only storage added by
# the lifecycle itself. Detector plans, path results, lifecycle bindings, and
# the event acquisition's singular `product` field are borrowed and are not
# counted here.
#

import ..Detectors

const _NO_ACQUISITION_STRUCTURAL_BYTES =
    (present=false, bytes=UInt64(0))

@inline _acquisition_optional_array_bytes(
    ::Nothing, ::AbstractComputeDevice, ::Symbol) =
    _NO_ACQUISITION_STRUCTURAL_BYTES

@inline function _acquisition_optional_array_bytes(
    array::AbstractArray,
    target::AbstractComputeDevice,
    component::Symbol,
)
    return _structural_array_target_bytes((array,), target, component)
end

@inline _acquisition_array_tuple_bytes(
    ::Tuple{}, ::AbstractComputeDevice, ::Symbol) =
    _NO_ACQUISITION_STRUCTURAL_BYTES

@inline function _acquisition_array_tuple_bytes(
    arrays::Tuple,
    target::AbstractComputeDevice,
    component::Symbol,
)
    head = _acquisition_optional_array_bytes(
        first(arrays), target, component)
    tail = _acquisition_array_tuple_bytes(
        Base.tail(arrays), target, component)
    return _combine_structural_target_bytes(head, tail, component)
end

# AcquisitionProducts metadata is logical metadata rather than numerical
# product storage. Supported product shapes enumerate their owned arrays
# explicitly; every other product variant fails closed.
function structural_resource_fact(
    products::AcquisitionProducts{<:AbstractArray,Nothing},
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    resident = _structural_array_target_bytes(
        (products.observation,), target, :resident_bytes)
    return _targeted_structural_resource_fact(
        id, target, resident, _NO_ACQUISITION_STRUCTURAL_BYTES)
end

function structural_resource_fact(
    products::AcquisitionProducts{<:WFSObservation,Nothing},
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    resident = _structural_array_target_bytes(
        (products.observation.storage,), target, :resident_bytes)
    return _targeted_structural_resource_fact(
        id, target, resident, _NO_ACQUISITION_STRUCTURAL_BYTES)
end

function structural_resource_fact(
    products::AcquisitionProducts{<:WFSObservation,<:WFSMeasurement},
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    resident = _structural_array_target_bytes((
        products.observation.storage,
        products.measurement.storage,
    ), target, :resident_bytes)
    return _targeted_structural_resource_fact(
        id, target, resident, _NO_ACQUISITION_STRUCTURAL_BYTES)
end

function structural_resource_fact(
    products::AcquisitionProducts{Nothing,<:WFSMeasurement},
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    resident = _structural_array_target_bytes(
        (products.measurement.storage,), target, :resident_bytes)
    return _targeted_structural_resource_fact(
        id, target, resident, _NO_ACQUISITION_STRUCTURAL_BYTES)
end

function structural_resource_fact(
    ::AcquisitionProducts,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    return UnknownStructuralResourceFact(
        id, target, :unsupported_acquisition_products)
end

# Detector parameter ownership is closed by dispatch. A known result either
# proves that the model is scalar/no-array or enumerates every owned array.
@inline _known_detector_parameter_bytes() =
    (known=true, storage=_NO_ACQUISITION_STRUCTURAL_BYTES)

@inline _unknown_detector_parameter_bytes() =
    (known=false, storage=_NO_ACQUISITION_STRUCTURAL_BYTES)

@inline function _known_detector_parameter_arrays(
    arrays::Tuple, target::AbstractComputeDevice)
    return (known=true, storage=_acquisition_array_tuple_bytes(
        arrays, target, :resident_bytes))
end

@inline _detector_parameter_bytes(
    model, ::AbstractComputeDevice) = _unknown_detector_parameter_bytes()

for model in (
    :NoiseNone,
    :NoisePhoton,
    :NoiseReadout,
    :NoisePhotonReadout,
    :SingleRead,
    :AveragedNonDestructiveReads,
    :CorrelatedDoubleSampling,
    :FowlerSampling,
    :UpTheRampSampling,
    :SkipperSampling,
    :GlobalShutter,
    :RollingExposure,
    :GlobalResetExposure,
    :NullFrameResponse,
    :NullChargeCoupling,
    :NullDetectorDefectModel,
    :NullFrameReadoutCorrection,
    :ReferencePixelCommonModeCorrection,
    :ReferenceRowCommonModeCorrection,
    :ReferenceColumnCommonModeCorrection,
    :ReferenceOutputCommonModeCorrection,
    :NullFrameNonlinearity,
    :SaturatingFrameNonlinearity,
    :NullPersistence,
    :ExponentialPersistence,
    :NullDetectorThermalModel,
    :NullTemperatureLaw,
    :ArrheniusRateLaw,
    :LinearTemperatureLaw,
    :ExponentialTemperatureLaw,
    :ScalarQuantumEfficiency,
    :NoBackground,
    :ScalarBackground,
    :NullCMOSReadNoise,
    :NullCMOSOutputModel,
    :ClippedGaussianMultiplicationApproximation,
    :ConditionalGammaMultiplication,
    :LinearEMMode,
    :PhotonCountingEMMode,
    :EMOutput,
    :ConventionalOutput,
    :SequentialAcquisition,
    :FrameTransferAcquisition,
    :ConditionalGammaAvalancheMultiplication,
    :ClippedGaussianAvalancheMultiplicationApproximation,
)
    @eval @inline _detector_parameter_bytes(
        ::Detectors.$model, ::AbstractComputeDevice) =
        _known_detector_parameter_bytes()
end

@inline function _detector_parameter_bytes(
    model::Detectors.RollingShutter,
    target::AbstractComputeDevice,
)
    return _detector_parameter_bytes(model.exposure_mode, target)
end

@inline function _detector_parameter_bytes(
    model::Detectors.FixedTemperature,
    target::AbstractComputeDevice,
)
    return _detector_parameter_tuple_bytes((
        model.dark_current_law,
        model.glow_rate_law,
        model.dark_count_law,
        model.cic_per_frame_law,
    ), target)
end

@inline function _detector_parameter_bytes(
    model::Detectors.FirstOrderThermalModel,
    target::AbstractComputeDevice,
)
    return _detector_parameter_tuple_bytes((
        model.dark_current_law,
        model.glow_rate_law,
        model.dark_count_law,
        model.cic_per_frame_law,
    ), target)
end

@inline function _detector_parameter_bytes(
    model::Detectors.GaussianPixelResponse,
    target::AbstractComputeDevice,
)
    return _known_detector_parameter_arrays((model.kernel,), target)
end

@inline function _detector_parameter_bytes(
    model::Detectors.SampledFrameResponse,
    target::AbstractComputeDevice,
)
    return _known_detector_parameter_arrays((model.kernel,), target)
end

@inline function _detector_parameter_bytes(
    model::Detectors.RectangularPixelAperture,
    target::AbstractComputeDevice,
)
    return _known_detector_parameter_arrays(
        (model.kernel_x, model.kernel_y), target)
end

@inline function _detector_parameter_bytes(
    model::Detectors.InterpixelCapacitance,
    target::AbstractComputeDevice,
)
    return _detector_parameter_bytes(model.response, target)
end

@inline function _detector_parameter_bytes(
    model::Detectors.PixelResponseNonuniformity,
    target::AbstractComputeDevice,
)
    return _known_detector_parameter_arrays((model.gain_map,), target)
end

@inline function _detector_parameter_bytes(
    model::Detectors.DarkSignalNonuniformity,
    target::AbstractComputeDevice,
)
    return _known_detector_parameter_arrays((model.dark_map,), target)
end

@inline function _detector_parameter_bytes(
    model::Detectors.BadPixelMask,
    target::AbstractComputeDevice,
)
    return _known_detector_parameter_arrays((model.mask,), target)
end

@inline function _detector_parameter_tuple_bytes(
    ::Tuple{}, ::AbstractComputeDevice)
    return _known_detector_parameter_bytes()
end

@inline function _combine_detector_parameter_bytes(left, right)
    left.known && right.known || return _unknown_detector_parameter_bytes()
    return (known=true, storage=_combine_structural_target_bytes(
        left.storage, right.storage, :resident_bytes))
end

@inline function _detector_parameter_tuple_bytes(
    models::Tuple, target::AbstractComputeDevice)
    return _combine_detector_parameter_bytes(
        _detector_parameter_bytes(first(models), target),
        _detector_parameter_tuple_bytes(Base.tail(models), target))
end

@inline function _detector_parameter_bytes(
    model::Detectors.CompositeDetectorDefectModel,
    target::AbstractComputeDevice,
)
    return _detector_parameter_tuple_bytes(model.stages, target)
end

@inline function _detector_parameter_bytes(
    model::Detectors.CompositeFrameReadoutCorrection,
    target::AbstractComputeDevice,
)
    return _detector_parameter_tuple_bytes(model.stages, target)
end

@inline function _detector_parameter_bytes(
    model::Detectors.SampledQuantumEfficiency,
    target::AbstractComputeDevice,
)
    return _known_detector_parameter_arrays(
        (model.wavelengths, model.values), target)
end

@inline function _detector_parameter_bytes(
    model::Detectors.BackgroundFrame,
    target::AbstractComputeDevice,
)
    return _known_detector_parameter_arrays((model.map,), target)
end

@inline function _detector_parameter_bytes(
    model::Detectors.CMOSReadNoiseMap,
    target::AbstractComputeDevice,
)
    return _known_detector_parameter_arrays((model.sigma,), target)
end

@inline function _detector_parameter_bytes(
    model::Detectors.StaticCMOSOutputPattern,
    target::AbstractComputeDevice,
)
    return _known_detector_parameter_arrays(
        (model.gains, model.offsets), target)
end

@inline function _detector_parameter_bytes(
    sensor::Detectors.CCDSensor,
    target::AbstractComputeDevice,
)
    return _detector_parameter_bytes(sensor.sampling_mode, target)
end

@inline function _detector_parameter_bytes(
    sensor::Detectors.CMOSSensor,
    target::AbstractComputeDevice,
)
    return _detector_parameter_tuple_bytes((
        sensor.readout_noise_model,
        sensor.output_model,
        sensor.timing_model,
    ), target)
end

@inline function _detector_parameter_bytes(
    sensor::Detectors.EMCCDSensor,
    target::AbstractComputeDevice,
)
    return _detector_parameter_tuple_bytes((
        sensor.multiplication_model,
        sensor.operating_mode,
        sensor.output_path,
        sensor.acquisition_mode,
    ), target)
end

@inline function _detector_parameter_bytes(
    sensor::Detectors.InGaAsSensor,
    target::AbstractComputeDevice,
)
    return _detector_parameter_bytes(sensor.persistence_model, target)
end

@inline function _detector_parameter_bytes(
    readout::Detectors.HgCdTeReadout,
    target::AbstractComputeDevice,
)
    return _detector_parameter_bytes(readout.sampling_mode, target)
end

@inline function _detector_parameter_bytes(
    sensor::Detectors.HgCdTeSensor,
    target::AbstractComputeDevice,
)
    return _detector_parameter_tuple_bytes(
        (sensor.readout, sensor.persistence_model), target)
end

@inline function _detector_parameter_bytes(
    sensor::Detectors.HgCdTeAvalancheArraySensor,
    target::AbstractComputeDevice,
)
    return _detector_parameter_tuple_bytes((
        sensor.multiplication_model,
        sensor.readout,
        sensor.persistence_model,
    ), target)
end

@inline _known_detector_state_storage() = (
    known=true,
    resident=_NO_ACQUISITION_STRUCTURAL_BYTES,
    workspace=_NO_ACQUISITION_STRUCTURAL_BYTES,
)

@inline _unknown_detector_state_storage() = (
    known=false,
    resident=_NO_ACQUISITION_STRUCTURAL_BYTES,
    workspace=_NO_ACQUISITION_STRUCTURAL_BYTES,
)

@inline _detector_thermal_state_storage(
    ::Detectors.NoThermalState) = true
@inline _detector_thermal_state_storage(
    ::Detectors.DetectorThermalState) = true
@inline _detector_thermal_state_storage(state) = false

@inline _detector_readout_product_storage(
    ::Detectors.FrameReadoutProducts,
    ::AbstractComputeDevice,
) = _unknown_detector_state_storage()

@inline _detector_readout_product_storage(
    ::Detectors.NoFrameReadoutProducts,
    ::AbstractComputeDevice,
) = _known_detector_state_storage()

@inline function _detector_readout_product_storage(
    products::Detectors.SkipperReadoutProducts,
    target::AbstractComputeDevice,
)
    resident = _structural_array_target_bytes(
        (products.mean_frame,), target, :resident_bytes)
    return (known=true, resident,
        workspace=_NO_ACQUISITION_STRUCTURAL_BYTES)
end

@inline function _detector_readout_product_storage(
    products::Detectors.SampledFrameReadoutProducts,
    target::AbstractComputeDevice,
)
    resident = _acquisition_array_tuple_bytes((
        products.reference_frame,
        products.signal_frame,
        products.read_cube,
    ), target, :resident_bytes)
    return (known=true, resident,
        workspace=_NO_ACQUISITION_STRUCTURAL_BYTES)
end

@inline function _detector_readout_product_storage(
    products::Detectors.MultiReadFrameReadoutProducts,
    target::AbstractComputeDevice,
)
    resident = _acquisition_array_tuple_bytes((
        products.reference_frame,
        products.signal_frame,
        products.combined_frame,
        products.reference_cube,
        products.signal_cube,
        products.read_cube,
        products.read_offsets_s,
    ), target, :resident_bytes)
    return (; known=true, resident,
        workspace=_NO_ACQUISITION_STRUCTURAL_BYTES)
end

@inline function _detector_readout_product_storage(
    products::Detectors.UpTheRampReadoutProducts,
    target::AbstractComputeDevice,
)
    resident = _structural_array_target_bytes((
        products.slope_frame,
        products.intercept_frame,
        products.integrated_frame,
        products.read_cube,
        products.read_offsets_s,
    ), target, :resident_bytes)
    return (; known=true, resident,
        workspace=_NO_ACQUISITION_STRUCTURAL_BYTES)
end

@inline _detector_readout_workspace_storage(
    ::Detectors.NoFrameReadoutWorkspace,
    ::AbstractComputeDevice,
) = _known_detector_state_storage()

@inline _detector_readout_workspace_storage(
    ::Detectors.FrameReadoutWorkspace,
    ::AbstractComputeDevice,
) = _unknown_detector_state_storage()

@inline function _detector_readout_workspace_storage(
    workspace::Detectors.SkipperReadoutWorkspace,
    target::AbstractComputeDevice,
)
    bytes = _structural_array_target_bytes((workspace.baseline_frame,
        workspace.sample_sum), target, :workspace_bytes)
    return (known=true, resident=_NO_ACQUISITION_STRUCTURAL_BYTES,
        workspace=bytes)
end

@inline function _detector_readout_workspace_storage(
    workspace::Detectors.MultiReadFrameReadoutWorkspace,
    target::AbstractComputeDevice,
)
    bytes = _acquisition_array_tuple_bytes((
        workspace.reference_average,
        workspace.signal_average,
        workspace.reference_cube,
        workspace.signal_cube,
    ), target, :workspace_bytes)
    return (known=true, resident=_NO_ACQUISITION_STRUCTURAL_BYTES,
        workspace=bytes)
end

@inline function _detector_readout_workspace_storage(
    workspace::Detectors.UpTheRampReadoutWorkspace,
    target::AbstractComputeDevice,
)
    bytes = _structural_array_target_bytes((workspace.slope,
        workspace.intercept, workspace.integrated, workspace.cube), target,
        :workspace_bytes)
    return (known=true, resident=_NO_ACQUISITION_STRUCTURAL_BYTES,
        workspace=bytes)
end

function _detector_runtime_resource_fact(
    state::DetectorState,
    workspace::Detectors.DetectorWorkspace,
    products::Detectors.DetectorProducts,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    _detector_thermal_state_storage(state.thermal_state) ||
        return UnknownStructuralResourceFact(
            id, target, :unsupported_detector_thermal_state)
    readout_products = _detector_readout_product_storage(
        products.readout, target)
    readout_products.known || return UnknownStructuralResourceFact(
        id, target, :unsupported_detector_readout_products)
    readout_workspace = _detector_readout_workspace_storage(
        workspace.readout, target)
    readout_workspace.known || return UnknownStructuralResourceFact(
        id, target, :unsupported_detector_readout_workspace)

    resident = _acquisition_array_tuple_bytes((
        products.frame,
        state.accum_buffer,
        state.latent_buffer,
        products.output_buffer,
    ), target, :resident_bytes)
    resident = _combine_structural_target_bytes(
        resident, readout_products.resident, :resident_bytes)
    workspace_bytes = _acquisition_array_tuple_bytes((
        workspace.presampling_buffer,
        workspace.presampling_scratch,
        workspace.response_buffer,
        workspace.bin_buffer,
        workspace.temporal_buffer,
        workspace.noise_buffer,
        workspace.noise_buffer_host,
        workspace.batched_buffer_host,
        workspace.output_buffer_host,
    ), target, :workspace_bytes)
    workspace_bytes = _combine_structural_target_bytes(
        workspace_bytes, readout_workspace.workspace, :workspace_bytes)
    return _targeted_structural_resource_fact(
        id, target, resident, workspace_bytes)
end

function structural_resource_fact(
    state::DetectorState,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    _detector_thermal_state_storage(state.thermal_state) ||
        return UnknownStructuralResourceFact(
            id, target, :unsupported_detector_thermal_state)
    resident = _structural_array_target_bytes((state.accum_buffer,
        state.latent_buffer), target, :resident_bytes)
    return _targeted_structural_resource_fact(id, target, resident,
        _NO_ACQUISITION_STRUCTURAL_BYTES)
end

function structural_resource_fact(
    detector::Detector,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    parameters = _detector_parameter_tuple_bytes((
        detector.noise,
        detector.params.sensor,
        detector.params.quantum_efficiency_model,
        detector.params.response_model,
        detector.params.charge_coupling_model,
        detector.params.defect_model,
        detector.params.timing_model,
        detector.params.correction_model,
        detector.params.nonlinearity_model,
        detector.params.thermal_model,
        detector.background_flux,
        detector.background_map,
    ), target)
    parameters.known || return UnknownStructuralResourceFact(
        id, target, :unsupported_detector_parameter_model)
    state_fact = _detector_runtime_resource_fact(detector.state,
        detector.workspace, detector.products, id, target)
    parameter_fact = _targeted_structural_resource_fact(
        id, target, parameters.storage,
        _NO_ACQUISITION_STRUCTURAL_BYTES)
    return _combine_structural_owner_facts(
        id, target, (state_fact, parameter_fact))
end

@inline function _detector_lifecycle_fact(
    prepared::Union{
        PreparedGlobalShutterAcquisition,
        PreparedRollingShutterAcquisition,
        PreparedFrameTransferAcquisition,
    },
    local_fact::AbstractStructuralResourceFact,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    return _combine_structural_owner_facts(id, target, (
        structural_resource_fact(prepared.detector, id, target),
        local_fact,
    ))
end

function structural_resource_fact(
    prepared::PreparedGlobalShutterAcquisition,
    state::GlobalShutterAcquisitionState,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    _require_detector_event_binding(prepared, state)
    resident = _acquisition_array_tuple_bytes(
        (prepared.read_offsets, prepared.read_offset_binding), target,
        :resident_bytes)
    local_fact = _targeted_structural_resource_fact(
        id, target, resident, _NO_ACQUISITION_STRUCTURAL_BYTES)
    return _detector_lifecycle_fact(
        prepared, local_fact, id, target)
end

function structural_resource_fact(
    prepared::PreparedRollingShutterAcquisition,
    state::RollingShutterAcquisitionState,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    _require_rolling_shutter_event_binding(prepared, state)
    local_fact = UnknownStructuralResourceFact(
        id, target, :owner_not_on_device)
    return _detector_lifecycle_fact(
        prepared, local_fact, id, target)
end

function structural_resource_fact(
    prepared::PreparedFrameTransferAcquisition,
    state::FrameTransferAcquisitionState,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    _require_frame_transfer_event_binding(prepared, state)
    resident = _structural_array_target_bytes(
        (state.storage_frame,), target, :resident_bytes)
    local_fact = _targeted_structural_resource_fact(
        id, target, resident, _NO_ACQUISITION_STRUCTURAL_BYTES)
    return _detector_lifecycle_fact(
        prepared, local_fact, id, target)
end

function structural_resource_fact(
    prepared::PreparedDirectMeasurementAcquisition,
    state::DirectMeasurementAcquisitionState,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    _require_direct_measurement_binding(prepared, state)
    resident = _structural_array_target_bytes(
        (prepared.instantaneous_sample,), target, :resident_bytes)
    workspace = _structural_array_target_bytes(
        (prepared.integrated_sample,), target, :workspace_bytes)
    return _targeted_structural_resource_fact(
        id, target, resident, workspace)
end

function structural_resource_fact(
    ::AbstractPreparedAcquisitionLifecycle,
    ::AbstractAcquisitionLifecycleState,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    return UnknownStructuralResourceFact(
        id, target, :unsupported_acquisition_lifecycle)
end

"""
Report one prepared event acquisition without double-counting its singular
`product` alias or its borrowed optical-path result.
"""
function structural_resource_fact(
    acquisition::_PreparedPlantEventAcquisition,
    lifecycle_state::_AcquisitionEventLifecycleState,
    id::StructuralResourceOwnerID,
    target::AbstractComputeDevice,
)
    lifecycle_fact = structural_resource_fact(
        acquisition.lifecycle, lifecycle_state, id, target)
    products_fact = structural_resource_fact(
        acquisition.products, id, target)
    return _combine_structural_owner_facts(
        id, target, (lifecycle_fact, products_fact))
end

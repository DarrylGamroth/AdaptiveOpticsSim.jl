module AdaptiveOpticsSim

__precompile__(true)

using AbstractFFTs
import FFTW
using KernelAbstractions
using LinearAlgebra
using Logging
using Random
using SparseArrays
using SpecialFunctions
using Statistics

"""
AdaptiveOpticsSim.jl

Julia adaptive optics simulation toolkit (in development).
"""
const PROJECT_STATUS = :in_development

include("core/errors.jl")
include("core/profiles.jl")
include("core/random_services.jl")
include("backends/backends.jl")

# The remaining flat source files are migrated owner by owner. Until their
# modules land, import the canonical Backends identities strictly for
# package-internal composition; the root API does not export or mark them
# public.
import .Backends:
    AMDGPUBackend,
    AMDGPUBackendTag,
    AbstractArrayBackend,
    AbstractComputeDevice,
    AbstractReductionExecutionPlan,
    AcceleratorComputeDevice,
    AcceleratorStyle,
    CPUBackend,
    CUDABackend,
    CUDABackendTag,
    DirectReductionPlan,
    ExecutionStyle,
    GPUBackendTag,
    GPUPrecisionPolicy,
    HostComputeDevice,
    HostMirrorReductionPlan,
    KernelLaunchPhase,
    MetalBackend,
    MetalBackendTag,
    ScalarCPUStyle,
    SplitGPUPrecision,
    UnifiedGPUPrecision,
    _AbstractPreparedDeviceExecutionContext,
    _HostPreparedDeviceExecutionContext,
    _clamp_array!,
    _poisson_noise!,
    _prepare_device_execution_context,
    _prepared_device_execution_compute_device,
    _rand_uniform_backend!,
    _randn_backend!,
    _resolve_array_backend,
    _resolve_backend_selector,
    _synchronize_prepared_device_execution_context!,
    _with_prepared_device_execution_context,
    _write_integer_output!,
    allocate_array,
    array_backend_selector,
    array_backend_type,
    available_gpu_backends,
    backend,
    backend_fill,
    backend_matmul,
    backend_matmul_transpose_right,
    backend_maximum_value,
    backend_rand,
    backend_randn,
    backend_sum_value,
    backend_symmetric_product,
    backend_type,
    backend_zeros,
    begin_kernel_phase,
    clamp_array!,
    compute_device,
    compute_device_backend,
    compute_device_identifier,
    default_gpu_precision_policy,
    disable_scalar_backend!,
    execute_fft_plan!,
    execution_style,
    finish_kernel_phase!,
    gpu_backend_array_type,
    gpu_backend_loaded,
    gpu_backend_name,
    gpu_build_type,
    gpu_runtime_type,
    high_accuracy_gpu_precision_policy,
    host_array,
    launch_kernel!,
    launch_kernel_async!,
    masked_sum2d,
    packed_valid_pair_mean,
    plan_fft_backend!,
    plan_ifft_backend!,
    poisson_noise!,
    poisson_noise_async!,
    queue_kernel!,
    rand_uniform_backend!,
    randn_backend!,
    randn_backend_async!,
    reduction_execution_plan,
    reduction_host_view,
    reduction_parent_source,
    require_same_backend,
    resolve_array_backend,
    same_backend,
    set_fft_provider_threads!,
    synchronize_backend!,
    use_host_build_algebra,
    write_integer_output!

include("core/inverse_policies.jl")
include("core/utils.jl")

# Configuration serialization remains a root-owned cross-domain seam. Declare
# its generic before domain modules add methods, then install the shared
# fallbacks after the first physical owner has loaded.
function config_value end

include("optics/optics.jl")

# The remaining flat source files are migrated owner by owner. Routine Optics
# vocabulary is imported for package-internal composition only; these imports
# are not root exports or qualified-public compatibility bindings.
using .Optics
import .Optics:
    AbstractControllableOptic,
    AbstractDeformableMirror,
    AbstractDMActuatorModel,
    AbstractDMInfluenceModel,
    AbstractDMTopology,
    AbstractOpticalElement,
    AbstractPropagationModel,
    AbstractSource,
    AbstractSourceCompositionStyle,
    AbstractSpectralCoordinate,
    AbstractTelescope,
    DMApplyMode,
    DeformableMirrorParams,
    DeformableMirrorState,
    ExpandedSourceComposition,
    ExtendedSource,
    LeafSourceComposition,
    NCPA,
    LinearStaticActuators,
    PreparedFocalPlaneModulation,
    PreparedMicrolensPropagation,
    PreparedBundledDirectImaging,
    PreparedDirectImaging,
    PreparedIncoherentDirectImaging,
    LGSProfileNaProfile,
    LGSProfileNone,
    _FixedOpticalProductVector,
    _prepare_microlens_propagation,
    _accumulate_field_intensity!,
    _accumulate_field_intensity_async!,
    _converted_nonnegative_finite,
    _converted_positive_finite,
    _direct_imaging_batch_product_contract,
    _pupil_cell_area,
    _pupil_diameter_m,
    _pupil_resolution,
    _require_physical_photon_irradiance,
    aperture_revision,
    actuator_coordinates,
    actuator_model,
    anamorphosis_angle_deg,
    apply_centering_phase!,
    apply_phase_async!,
    axis_centering,
    centered_grid_origin,
    coordinates_xy_arcsec,
    electric_field_wavelength,
    field_embedding_offsets,
    freeze_source,
    fraunhofer_intensity_stack!,
    intensity!,
    is_leaf_source,
    is_lgs_source,
    lgs_elongation_factor,
    lgs_profile,
    photon_irradiance,
    pixel_mask_grid,
    plane_metadata,
    modulation_point_count,
    microlens_array,
    misregistration_component,
    prepare_focal_plane_modulation,
    require_leaf_source,
    require_metric_coordinates,
    require_same_plane_grid,
    require_centered_plane_geometry,
    source_measurement_signature,
    source_composition_style,
    source_height_m,
    source_with_wavelength_and_radiometric_value,
    spectral_bundle,
    sampled_influence_matrix,
    rotation_deg,
    supports_dm_misregistration_identification,
    surface_opd,
    topology,
    topology_axis_count,
    topology_command_count,
    update_cycle_averaged_circular_modulation!,
    validate_dm_actuator_model,
    wavelength,
    validate_plane_storage

include("core/types.jl")
include("core/config.jl")
include("core/kv56.jl")
include("core/workspace.jl")
include("core/parallel.jl")
include("core/telemetry.jl")

include("atmosphere/atmospheres.jl")

# Atmosphere implementations have one authoritative owner. Root composition
# imports the module's routine surface for domains awaiting their own namespace
# migrations; advanced internal seams are imported explicitly below as needed.
using .Atmospheres
import .Atmospheres:
    AbstractAtmosphericFieldModel

include("detectors/detectors.jl")

# Detector implementations have one authoritative type and generic-function
# owner. Root composition imports only the internal seams still required by
# domains that have not yet completed their namespace migrations.
using .Detectors
import .Detectors:
    AbstractChargeCouplingModel,
    AbstractCountingDetector,
    AbstractCountingCorrelationModel,
    AbstractCountingGateModel,
    AbstractDetectorResponse,
    AbstractEMCCDOperatingMode,
    AbstractEMCCDOutputPath,
    AbstractFrameDetector,
    AbstractFrameResponse,
    AbstractFrameTimingModel,
    AbstractEMCCDAcquisitionMode,
    AbstractQuantumEfficiencyModel,
    BackgroundModel,
    CountingReadoutMetadata,
    CountingSensorType,
    FrameSamplingMode,
    FrameSensorType,
    NoBackground,
    NullChargeCoupling,
    NullDetectorDefectModel,
    SampledQuantumEfficiency,
    ScalarQuantumEfficiency,
    _copy_windowed_sampling_plane!,
    _raw_sampling_sigma,
    _require_finite_nonnegative_intensity,
    _require_prepared_acquisition,
    _require_prepared_response_sampling,
    accumulate_incremental_charge_generation!,
    advance_thermal!,
    apply_avalanche_excess_noise!,
    apply_quantization!,
    apply_readout_correction!,
    apply_response!,
    capture_signal_pipeline!,
    capture_stack!,
    capture_stack_with_quantum_efficiency!,
    capture_with_quantum_efficiency!,
    counting_array,
    counting_exposure_time,
    counting_fill_factor,
    counting_integration_time,
    counting_post_gain,
    counting_qe,
    counting_source_throughput,
    detector_noise_symbol,
    detector_output_type,
    detector_output_value,
    detector_readout_sigma,
    effective_qe,
    ensure_buffers!,
    ensure_up_the_ramp_products!,
    finalize_charge_transport!,
    finalize_incremental_capture!,
    finalize_scheduled_up_the_ramp_readout_products!,
    is_global_shutter,
    is_null_cmos_output_model,
    is_null_frame_nonlinearity,
    is_null_persistence,
    line_time,
    persistence_model,
    prepare_signal_frame!,
    qe_at,
    quantum_efficiency_model,
    readout_noise,
    readout_products,
    require_whole_capture_idle,
    response_width_px,
    sample_frame_read!,
    sampling_read_time,
    subtract_background_map!,
    update_sensor_persistence!,
    validate_frame_response_model,
    validate_up_the_ramp_schedule,
    write_output!
include("calibration/modal_basis.jl")
include("calibration/ncpa.jl")
include("wfs/wavefront_sensors.jl")

# Concrete WFS families remain flat only until their ordered namespace gates.
# Import the common owner's identities explicitly so every later method extends
# one authoritative generic rather than creating a root-owned duplicate.
import .WavefrontSensors:
    AbstractGroupedAccumulationPlan,
    AbstractWFS,
    AbstractWFSMeasurementPath,
    AcquiredObservationPath,
    Diffractive,
    DirectMeasurementPath,
    Geometric,
    GroupedDirectAccumulatePlan,
    GroupedStackReducePlan,
    GroupedStaged2DPlan,
    IncidenceFluxNormalization,
    MeanValidFluxNormalization,
    SensingMode,
    WFSMeasurement,
    WFSMeasurementMetadata,
    WFSNormalization,
    WFSObservation,
    WFSObservationMetadata,
    WFSPreparationError,
    _capture_counting_wfs!,
    _require_counting_wfs_source,
    _require_counting_wfs_spectral_match,
    _require_declared_wfs_units,
    _require_real_square_wfs_observation,
    _require_wfs_storage_domain,
    acquire_wfs_observation!,
    accumulate_grouped_sources!,
    apply_shift_wfs!,
    camera_frame,
    deterministic_frame_readout_gain,
    estimate_wfs_measurement!,
    form_wfs_optical_products!,
    grouped_accumulation_plan,
    grouped_stack_view,
    grouped_staging_buffer,
    measure!,
    measurement_metadata,
    measurement_storage,
    measurement_units,
    observation_metadata,
    observation_storage,
    observation_units,
    prepare_runtime_wfs!,
    prepare_wfs_acquisition,
    prepare_wfs_estimation,
    prepare_wfs_optical_formation,
    reference_signal,
    reduce_grouped_blocks_kernel!,
    reduce_grouped_stack!,
    sensing_mode,
    slopes,
    supports_camera_frame,
    supports_detector_output,
    supports_grouped_execution,
    supports_prepared_runtime,
    supports_reference_signal,
    supports_stacked_sources,
    supports_valid_subaperture_mask,
    usable_wfs_normalization,
    validate_wfs_acquisition_binding,
    validate_wfs_estimation_binding,
    validate_wfs_measurement,
    validate_wfs_observation,
    validate_wfs_observations,
    validate_wfs_optical_formation_binding,
    validate_wfs_optical_input,
    validate_wfs_optical_products,
    valid_subaperture_mask,
    wfs_calibration_signature,
    wfs_detector_image,
    wfs_detector_incidence_scale,
    wfs_incident_photon_irradiance,
    wfs_measurement_path,
    wfs_output_frame,
    wfs_output_frame_prototype,
    wfs_output_metadata

include("wfs/calibration.jl")
include("wfs/elongation.jl")
include("wfs/subapertures.jl")
include("wfs/focal_plane_modulation.jl")
include("wfs/shack_hartmann.jl")
include("wfs/pyramid.jl")
include("wfs/bioedge.jl")
include("wfs/zernike.jl")
include("wfs/curvature.jl")
include("wfs/lift.jl")
include("plant/plant.jl")
include("calibration/interaction_matrix.jl")
include("control/controller.jl")
include("control/reconstructors.jl")
include("calibration/control_matrix.jl")
include("calibration/fitting_error.jl")
include("calibration/ao_calibration.jl")
include("calibration/gain_sensing_camera.jl")
include("calibration/misregistration_identification.jl")
include("calibration/ad_sensitivities.jl")
include("control/ensemble.jl")
include("tomography/parameters.jl")
include("tomography/geometry.jl")
include("tomography/fitting.jl")
include("tomography/reconstructors.jl")

include("exports.jl")

end # module AdaptiveOpticsSim

"""
    Plant

Virtual-time adaptive-optics plant construction and execution.

`Plant` owns the HIL-neutral definitions, schedules, command lifecycle,
acquisition events, illumination, preparation, and event composition. Optical,
atmospheric, detector, and wavefront-sensor models remain owned by their
respective package domains and enter through the explicit imports below.
"""
module Plant

using KernelAbstractions
using LinearAlgebra
using Random

import ..Backends:
    AbstractArrayBackend,
    AbstractComputeDevice,
    AcceleratorComputeDevice,
    AcceleratorStyle,
    CPUBackend,
    HostComputeDevice,
    ScalarCPUStyle,
    _prepare_device_execution_context,
    _prepared_device_execution_compute_device,
    _synchronize_prepared_device_execution_context!,
    _with_prepared_device_execution_context,
    allocate_array,
    backend,
    clamp_array!,
    compute_device,
    compute_device_backend,
    execution_style,
    launch_kernel!,
    resolve_array_backend

import ..Optics:
    AbstractCombinationPolicy,
    AbstractDMActuatorModel,
    AbstractDMInfluenceModel,
    AbstractDMTopology,
    AbstractOpticalNormalization,
    AbstractOpticalPlaneKind,
    AbstractOpticalProduct,
    AbstractSource,
    AbstractSourceCompositionStyle,
    AbstractSpatialMeasure,
    AbstractSpectralCoordinate,
    AbstractTelescope,
    AchromaticSpectralCoordinate,
    ActuatorGridTopology,
    Asterism,
    CellIntegratedMeasure,
    CircularModulation,
    CoherentFieldCombination,
    DMAdditive,
    DMApplyMode,
    DMReplace,
    DeformableMirror,
    DeformableMirrorParams,
    DeformableMirrorState,
    DetectorPlane,
    DimensionlessNormalization,
    ElectricField,
    ExpandedSourceComposition,
    ExtendedSource,
    FocalPlane,
    GaussianInfluenceWidth,
    GaussianMechanicalCoupling,
    IncoherentIntensityAddition,
    IntegratedSpectralChannel,
    IntensityMap,
    LGSSource,
    LeafSourceComposition,
    LinearStaticActuators,
    MetricCoordinates,
    Misregistration,
    MonochromaticChannel,
    NCPA,
    NonCombinableProduct,
    OPDMap,
    OpticalPlaneMetadata,
    OpticalProductBundle,
    PhotonRateNormalization,
    PhysicalPhotonIrradianceSource,
    PointSampledMeasure,
    PreparedBundledDirectImaging,
    PreparedDirectImaging,
    PreparedDirectImagingBatch,
    PreparedFocalPlaneModulation,
    PreparedIncoherentDirectImaging,
    PupilFunction,
    PupilPlane,
    Source,
    SpatialDensityMeasure,
    SpectralSource,
    Telescope,
    UnspecifiedCoherence,
    UnspecifiedNormalization,
    UnspecifiedSpatialMeasure,
    UnspecifiedSpectralCoordinate,
    _FixedOpticalProductVector,
    _direct_imaging_batch_product_contract,
    actuator_model,
    aperture_revision,
    coordinates_xy_arcsec,
    direct_imaging_output,
    form_direct_image!,
    freeze_source,
    influence_model,
    n_actuators,
    prepare_direct_imaging_batch,
    pupil_reflectivity,
    set_command!,
    source_composition_style,
    source_height_m,
    source_radiometric_value,
    source_radiometry,
    surface_opd,
    topology,
    topology_command_count,
    update_cycle_averaged_circular_modulation!,
    update_surface!,
    validate_direct_imaging_batch,
    validate_dm_actuator_model,
    validate_plane_storage,
    wavelength

import ..Detectors:
    AbstractEMCCDAcquisitionMode,
    AbstractFrameTimingModel,
    CMOSSensor,
    Detector,
    DetectorAcquisitionPlan,
    EMCCDSensor,
    FrameReadoutProducts,
    FrameSensorType,
    FrameTransferAcquisition,
    GlobalResetExposure,
    HgCdTeAvalancheArraySensor,
    RollingExposure,
    RollingShutter,
    UpTheRampReadoutProducts,
    UpTheRampSampling,
    _copy_windowed_sampling_plane!,
    _raw_sampling_sigma,
    _require_prepared_acquisition,
    accumulate_incremental_charge_generation!,
    advance_thermal!,
    apply_avalanche_excess_noise!,
    apply_quantization!,
    capture!,
    capture_signal_pipeline!,
    ensure_up_the_ramp_products!,
    finalize_charge_transport!,
    finalize_incremental_capture!,
    finalize_scheduled_up_the_ramp_readout_products!,
    is_global_shutter,
    line_time,
    output_frame,
    prepare_detector_acquisition,
    readout_products,
    readout_ready,
    sample_frame_read!,
    sampling_read_time,
    subtract_background_map!,
    update_sensor_persistence!,
    write_output!

import ..Atmospheres:
    AbstractAtmosphere,
    AbstractTimedAtmosphere,
    AtmosphereEpoch,
    AtmosphereLayerID,
    InfiniteAtmosphereLayer,
    InfiniteMultiLayerAtmosphere,
    MovingAtmosphereLayer,
    MultiLayerAtmosphere,
    PreparedAtmosphereDirectionBatch,
    _advance_by_with_rng!,
    _advance_to_with_rng!,
    _validate_atmosphere_destination,
    _validate_atmosphere_renderer_binding,
    _validate_epoch_identity,
    advance_by!,
    advance_to!,
    atmosphere_identity,
    atmosphere_timeline,
    ensure_initialized!,
    epoch_time,
    evolve_atmosphere!,
    evolve_initial_atmosphere!,
    evolve_layer!,
    initialize_atmosphere!,
    prepare_atmosphere_direction_batch,
    prepare_atmosphere_renderer,
    render_atmosphere!,
    render_atmosphere_directions!,
    validate_atmosphere_direction_batch,
    validate_atmosphere_rendering

import ..AdaptiveOpticsSim:
    AdaptiveOpticsSimError,
    InvalidConfiguration,
    PreparedBioEdgeOpticalBundleFormation,
    PreparedBioEdgeOpticalFormation,
    PreparedPyramidOpticalBundleFormation,
    PreparedPyramidOpticalFormation,
    PreparedShackHartmannOpticalBundleFormation,
    PreparedShackHartmannOpticalFormation,
    runtime_rng,
    splitmix64

import ..WavefrontSensors:
    AcquiredObservationPath,
    DirectMeasurementPath,
    WFSMeasurement,
    WFSObservation,
    acquire_wfs_observation!,
    estimate_wfs_measurement!,
    form_wfs_optical_products!,
    measurement_storage,
    validate_wfs_acquisition_binding,
    validate_wfs_estimation_binding,
    validate_wfs_measurement,
    validate_wfs_observation,
    validate_wfs_optical_formation_binding,
    wfs_measurement_path

include("errors.jl")
include("time.jl")
include("scheduling.jl")
include("identities.jl")
include("optical_placement.jl")
include("command_schemas.jl")
include("command_admission.jl")
include("command_application.jl")
include("definitions.jl")
include("controllable_optics.jl")
include("triggers.jl")
include("acquisition_lifecycles.jl")
include("detector_acquisition_events.jl")
include("rolling_frame_transfer_events.jl")
include("direct_measurement_events.jl")
include("randomness.jl")
include("product_providers.jl")
include("preparation.jl")
include("reduced_order.jl")
include("controller_routing.jl")
include("autonomous_periodic_optics.jl")
include("pupil_footprint_coupling.jl")
include("sampled_aberrations.jl")
include("deformable_mirror.jl")
include("cpu_execution.jl")
include("event_composition.jl")
include("device_batching.jl")
include("illumination.jl")
include("api.jl")

end # module Plant

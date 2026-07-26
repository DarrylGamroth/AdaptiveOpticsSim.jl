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

import ..AdaptiveOpticsSim: AbstractArrayBackend, AbstractAtmosphere,
    AbstractCombinationPolicy, AbstractEMCCDAcquisitionMode,
    AbstractDMActuatorModel, AbstractDMInfluenceModel, AbstractDMTopology,
    AbstractFrameTimingModel, AbstractOpticalNormalization,
    AbstractOpticalPlaneKind, AbstractOpticalProduct, AbstractComputeDevice,
    AbstractSource,
    AbstractSourceCompositionStyle,
    AbstractSpatialMeasure, AbstractSpectralCoordinate, AbstractTelescope,
    AbstractTimedAtmosphere, AcquiredObservationPath, AdaptiveOpticsSimError,
    AcceleratorStyle, AchromaticSpectralCoordinate, ActuatorGridTopology,
    Asterism, AtmosphereEpoch, AtmosphereLayerID,
    CMOSSensor, CPUBackend,
    CellIntegratedMeasure, CoherentFieldCombination, Detector,
    DeformableMirror, DeformableMirrorParams, DeformableMirrorState,
    DetectorAcquisitionPlan, DetectorPlane, DimensionlessNormalization,
    CircularModulation, DirectMeasurementPath, DMAdditive, DMApplyMode,
    DMReplace, EMCCDSensor, ElectricField,
    ExpandedSourceComposition, ExtendedSource,
    FocalPlane, FrameReadoutProducts, FrameSensorType,
    FrameTransferAcquisition, GlobalResetExposure,
    GaussianInfluenceWidth, GaussianMechanicalCoupling,
    HgCdTeAvalancheArraySensor, IncoherentIntensityAddition,
    InfiniteAtmosphereLayer, InfiniteMultiLayerAtmosphere,
    IntegratedSpectralChannel, IntensityMap, InvalidConfiguration, LGSSource,
    LeafSourceComposition, LinearStaticActuators, MetricCoordinates,
    Misregistration, MonochromaticChannel,
    MovingAtmosphereLayer,
    MultiLayerAtmosphere,
    NonCombinableProduct, OpticalProductBundle, PhotonRateNormalization,
    PreparedBundledDirectImaging, PreparedDirectImaging,
    PointSampledMeasure, PreparedFocalPlaneModulation,
    PreparedIncoherentDirectImaging,
    PreparedPyramidOpticalBundleFormation, PreparedPyramidOpticalFormation,
    NCPA, OPDMap, OpticalPlaneMetadata, PupilFunction, PupilPlane,
    RollingExposure, RollingShutter, ScalarCPUStyle, Source,
    SpatialDensityMeasure, SpectralSource, Telescope, UnspecifiedCoherence,
    UnspecifiedNormalization, UnspecifiedSpatialMeasure,
    UnspecifiedSpectralCoordinate, UpTheRampReadoutProducts,
    UpTheRampSampling, WFSMeasurement, WFSObservation,
    _FixedOpticalProductVector, _advance_by_with_rng!, _advance_to_with_rng!,
    _copy_windowed_sampling_plane!, _raw_sampling_sigma,
    _require_prepared_acquisition, _validate_atmosphere_destination,
    _validate_atmosphere_renderer_binding, _validate_epoch_identity,
    accumulate_incremental_charge_generation!, acquire_wfs_observation!,
    advance_by!, advance_thermal!, advance_to!, allocate_array,
    actuator_model, aperture_revision, apply_avalanche_excess_noise!,
    apply_quantization!, atmosphere_identity, atmosphere_timeline, backend, capture!,
    capture_signal_pipeline!, clamp_array!, coordinates_xy_arcsec,
    compute_device_backend, direct_imaging_output, ensure_initialized!,
    execution_style,
    ensure_up_the_ramp_products!, epoch_time, estimate_wfs_measurement!,
    evolve_atmosphere!, evolve_initial_atmosphere!, evolve_layer!,
    finalize_charge_transport!, finalize_incremental_capture!,
    finalize_scheduled_up_the_ramp_readout_products!, form_direct_image!,
    form_wfs_optical_products!, freeze_source, initialize_atmosphere!,
    influence_model, is_global_shutter, launch_kernel!, line_time,
    n_actuators, output_frame, compute_device,
    prepare_atmosphere_renderer, prepare_detector_acquisition,
    pupil_reflectivity, readout_products, readout_ready, render_atmosphere!,
    resolve_array_backend, runtime_rng, sample_frame_read!,
    sampling_read_time, source_height_m, source_radiometric_value,
    source_composition_style, source_radiometry, splitmix64,
    set_command!, subtract_background_map!, surface_opd, topology,
    topology_axis_count, topology_command_count,
    update_cycle_averaged_circular_modulation!,
    update_sensor_persistence!, validate_atmosphere_rendering,
    update_surface!, validate_dm_actuator_model, validate_plane_storage,
    validate_wfs_acquisition_binding,
    validate_wfs_estimation_binding, validate_wfs_measurement,
    validate_wfs_observation, validate_wfs_optical_formation_binding,
    measurement_storage, wavelength, wfs_measurement_path, write_output!

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
include("illumination.jl")
include("api.jl")

end # module Plant

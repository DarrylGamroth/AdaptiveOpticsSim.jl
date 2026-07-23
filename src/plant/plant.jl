"""
    Plant

Virtual-time adaptive-optics plant construction and execution.

`Plant` owns the HIL-neutral definitions, schedules, command lifecycle,
acquisition events, illumination, preparation, and event composition. Optical,
atmospheric, detector, and wavefront-sensor models remain owned by their
respective package domains and enter through the explicit imports below.
"""
module Plant

using LinearAlgebra
using Random

import ..AdaptiveOpticsSim: AbstractArrayBackend, AbstractAtmosphere,
    AbstractCombinationPolicy, AbstractEMCCDAcquisitionMode,
    AbstractFrameTimingModel, AbstractOpticalNormalization,
    AbstractOpticalPlaneKind, AbstractPlaneDevice, AbstractSource,
    AbstractSpatialMeasure, AbstractSpectralCoordinate, AbstractTelescope,
    AbstractTimedAtmosphere, AcquiredObservationPath, AdaptiveOpticsSimError,
    Asterism, AtmosphereEpoch, AtmosphereLayerID, CMOSSensor, CPUBackend,
    CellIntegratedMeasure, CoherentFieldCombination, Detector,
    DetectorAcquisitionPlan, DetectorPlane, DimensionlessNormalization,
    DirectMeasurementPath, EMCCDSensor, ElectricField, ExtendedSource,
    FocalPlane, FrameReadoutProducts, FrameSensorType,
    FrameTransferAcquisition, GlobalResetExposure,
    HgCdTeAvalancheArraySensor, IncoherentIntensityAddition,
    InfiniteAtmosphereLayer, InfiniteMultiLayerAtmosphere,
    IntegratedSpectralChannel, IntensityMap, InvalidConfiguration, LGSSource,
    MonochromaticChannel, MovingAtmosphereLayer, MultiLayerAtmosphere,
    NonCombinableProduct, OpticalProductBundle, PhotonRateNormalization,
    PreparedBundledDirectImaging, PreparedDirectImaging,
    PreparedIncoherentDirectImaging, PupilFunction, PupilPlane,
    RollingExposure, RollingShutter, Source,
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
    aperture_revision, apply_avalanche_excess_noise!, apply_quantization!,
    atmosphere_identity, atmosphere_timeline, backend, capture!,
    capture_signal_pipeline!, clamp_array!, coordinates_xy_arcsec,
    direct_imaging_output, ensure_initialized!,
    ensure_up_the_ramp_products!, epoch_time, estimate_wfs_measurement!,
    evolve_atmosphere!, evolve_initial_atmosphere!, evolve_layer!,
    finalize_charge_transport!, finalize_incremental_capture!,
    finalize_scheduled_up_the_ramp_readout_products!, form_direct_image!,
    form_wfs_optical_products!, freeze_source, initialize_atmosphere!,
    is_global_shutter, line_time, output_frame, plane_device,
    prepare_atmosphere_renderer, prepare_detector_acquisition,
    pupil_reflectivity, readout_products, readout_ready, render_atmosphere!,
    resolve_array_backend, runtime_rng, sample_frame_read!,
    sampling_read_time, source_height_m, source_radiometric_value,
    source_radiometry, splitmix64, subtract_background_map!, topology,
    update_sensor_persistence!, validate_atmosphere_rendering,
    validate_plane_storage, validate_wfs_acquisition_binding,
    validate_wfs_estimation_binding, validate_wfs_measurement,
    validate_wfs_observation, validate_wfs_optical_formation_binding,
    measurement_storage, wavelength, wfs_measurement_path, write_output!

include("errors.jl")
include("time.jl")
include("scheduling.jl")
include("identities.jl")
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
include("event_composition.jl")
include("illumination.jl")
include("api.jl")

end # module Plant

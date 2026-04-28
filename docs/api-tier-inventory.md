# API Tier Inventory

Date: 2026-03-31

Status: active

Plan traceability:

- [`PLAN-01`](./package-review-action-plan.md)
- review IDs: `PR-01`, `PR-02`, `PR-04`

## Purpose

This document inventories the current exported API surface and classifies it
into proposed tiers for Phase 1 API curation.

The goal is not to finalize the Phase 1 export set here. The goal is to:

- record the current baseline
- ensure every exported symbol is mapped somewhere
- identify clear de-export candidates
- identify the stable workflow surface that should remain easy to discover

## Baseline

- Phase 0 baseline export count from
  [`src/AdaptiveOpticsSim.jl`](../src/AdaptiveOpticsSim.jl): `566`
- Current pre-`PSP-01` export count after later feature additions:
  - `544`

## Phase 1 implementation status

- Current top-level export count after Phase 1 API curation:
  - `533`
- Change:
  - `-33` exported names

Phase 1 intentionally targeted the first low-risk tranche:

- telemetry, trace, workspace, and config helpers moved to namespaced access
- scenario-builder and misregistration-identification convenience APIs moved to
  namespaced access

Phase 1 intentionally did **not** de-export the backend/build helper surface
yet, because that tooling is still used broadly by maintained benchmark,
backend-validation, and support scripts. That surface remains documented as
advanced/developer-oriented and should be revisited in a later curation pass.

## `PSP-01` implementation status

- Current top-level export count after `PSP-02` curation:
  - `547`
- Change from the current pre-`PSP-02` baseline:
  - `-4` exported names

`PSP-02` intentionally targeted the lowest-risk remaining backend tag tranche:

- `GPUBackendTag` moved to namespaced access
- `CUDABackendTag` moved to namespaced access
- `MetalBackendTag` moved to namespaced access
- `AMDGPUBackendTag` moved to namespaced access

These names remain part of the advanced/developer surface for benchmark, smoke,
and backend-audit scripts, but no longer occupy the default package namespace.

### `PSP-02` de-exported symbols

- `GPUBackendTag`
- `CUDABackendTag`
- `MetalBackendTag`
- `AMDGPUBackendTag`

## `PSP-01` implementation status

- Current top-level export count after `PSP-01` curation:
  - `538`
- Change from the current pre-`PSP-01` baseline:
  - `-6` exported names

`PSP-01` intentionally targeted the lowest-risk backend/build helper tranche:

- build-backend policy types moved to namespaced access
- `default_build_backend` moved to namespaced access
- `set_fft_provider_threads!` moved to namespaced access

`PSP-01` intentionally did **not** de-export the broader GPU backend tag and
allocation-helper surface yet, because that surface is still used heavily by
maintained benchmark, smoke, and backend-audit scripts. That larger tranche
should be revisited only when those developer flows can be updated in one
coherent pass.

### `PSP-01` de-exported symbols

- `BuildBackend`
- `NativeBuildBackend`
- `CPUBuildBackend`
- `GPUArrayBuildBackend`
- `default_build_backend`
- `set_fft_provider_threads!`

### Phase 1 de-exported symbols

Phase 1 moved the following names to documented namespaced access:

- `PROJECT_STATUS`
- `Workspace`, `ensure_psf_buffers!`
- `Telemetry`, `TelemetryRow`, `record!`
- `ClosedLoopTrace`, `ClosedLoopTraceRow`
- `GSCClosedLoopTrace`, `GSCClosedLoopTraceRow`
- `GSCAtmosphereReplayTrace`, `GSCAtmosphereReplayTraceRow`
- `snapshot_config`, `write_config_toml`, `write_config_json`
- `write_telemetry_csv`
- `bin2d`, `poisson_noise!`, `poisson_sample`
- `fast_atmosphere`
- `AOSimulation`, `initialize_ao_shack_hartmann`, `initialize_ao_pyramid`
- `GSCDetectorMetadata`, `detector_metadata`
- `weak_mode_mask`, `attach_detector!`, `detach_detector!`
- `MetaSensitivity`, `compute_meta_sensitivity_matrix`
- `estimate_misregistration`
- `SPRINT`, `estimate!`

## Tier Definitions

- `Stable`
  - end-user workflow API
  - should remain easy to discover and usually exported
- `Advanced`
  - expert-facing simulation, modeling, and analysis API
  - should be documented, but does not necessarily need to stay exported
- `De-export candidate`
  - infrastructure, backend, trait, metadata, or low-level hook surface
  - should likely remain accessible as `AdaptiveOpticsSim.<name>` but not live
    in the default namespace

## Inventory

### Group 1: Project status, exceptions, and core profiles

Proposed tier: `Advanced`

Rationale:

- important for expert usage and diagnostics
- not the main workflow surface for most users

Symbols:

- `PROJECT_STATUS`
- `AdaptiveOpticsSimError`
- `InvalidConfiguration`
- `DimensionMismatchError`
- `UnsupportedAlgorithm`
- `NumericalConditionError`
- `FidelityProfile`
- `ScientificProfile`
- `FastProfile`
- `ProfileBundle`
- `default_fidelity_profile`
- `atmosphere_profile`
- `calibration_profile`
- `detector_profile`
- `lift_profile`
- `tomography_profile`
- `InversePolicy`
- `ExactPseudoInverse`
- `TSVDInverse`
- `TikhonovInverse`
- `InverseStats`
- `inverse_operator`
- `default_modal_inverse_policy`
- `default_calibration_inverse_policy`
- `default_projector_inverse_policy`

### Group 2: Build backends, workspace, telemetry, and config I/O

Proposed tier: `De-export candidate`

Rationale:

- useful infrastructure
- too low-level for the default namespace

Symbols:

- `BuildBackend`
- `NativeBuildBackend`
- `CPUBuildBackend`
- `GPUArrayBuildBackend`
- `default_build_backend`
- `Workspace`
- `ensure_psf_buffers!`
- `Telemetry`
- `TelemetryRow`
- `record!`
- `ClosedLoopTrace`
- `ClosedLoopTraceRow`
- `GSCClosedLoopTrace`
- `GSCClosedLoopTraceRow`
- `GSCAtmosphereReplayTrace`
- `GSCAtmosphereReplayTraceRow`
- `snapshot_config`
- `write_config_toml`
- `write_config_json`
- `write_telemetry_csv`
- `bin2d`
- `poisson_noise!`
- `poisson_sample`

### Group 3: Aperture and mask primitives

Proposed tier: `Stable`

Rationale:

- core optics construction surface
- directly useful and conceptually central

Symbols:

- `AbstractMaskPrimitive`
- `MaskGrid`
- `CircularAperture`
- `AnnularAperture`
- `SpiderMask`
- `RectangularROI`
- `SubapertureGridMask`
- `default_mask_grid`
- `pixel_mask_grid`
- `build_mask!`
- `apply_mask!`

### Group 4: Telescope

Proposed tier: `Stable`

Symbols:

- `Telescope`
- `TelescopeParams`
- `TelescopeState`
- `generate_pupil!`
- `reset_opd!`
- `apply_opd!`
- `set_pupil!`
- `set_pupil_reflectivity!`
- `flux_map`
- `apply_spiders!`

### Group 5: Source and spectrum

Proposed tier: `Stable`

Symbols:

- `Source`
- `SourceParams`
- `LGSSource`
- `LGSSourceParams`
- `wavelength`
- `optical_path`
- `print_optical_path`
- `SpectralSample`
- `SpectralBundle`
- `SpectralSource`
- `with_spectrum`
- `spectral_bundle`
- `spectral_reference_source`
- `has_spectral_bundle`
- `is_polychromatic`
- `weighted_wavelength`

### Group 6: Extended-source models

Proposed tier: `Stable`

Symbols:

- `PointCloudSourceModel`
- `GaussianDiskSourceModel`
- `SampledImageSourceModel`
- `ExtendedSource`
- `with_extended_source`
- `has_extended_source_model`
- `extended_source_model`
- `extended_source_asterism`
- `lgs_elongation_factor`

### Group 7: Electric field and propagation

Proposed tier: `Stable`

Symbols:

- `ElectricField`
- `ElectricFieldParams`
- `ElectricFieldState`
- `fill_from_telescope!`
- `fill_telescope_field!`
- `apply_phase!`
- `apply_amplitude!`
- `intensity!`
- `AbstractPropagationModel`
- `FraunhoferPropagation`
- `FresnelPropagation`
- `propagate_field!`

### Group 8: Atmospheric field propagation

Proposed tier: `Advanced`

Rationale:

- important and maintained
- slightly more expert-facing than the simpler optics workflow surface

Symbols:

- `AbstractAtmosphericFieldModel`
- `GeometricAtmosphericPropagation`
- `LayeredFresnelAtmosphericPropagation`
- `AtmosphericFieldPropagationParams`
- `AtmosphericFieldPropagationState`
- `AtmosphericFieldSlice`
- `AtmosphericFieldPropagation`
- `propagate_atmosphere_field!`
- `atmospheric_intensity!`

### Group 9: PSF, asterism, basis, OPD, NCPA, and spatial filtering

Proposed tier: `Stable`

Symbols:

- `Asterism`
- `coordinates_xy_arcsec`
- `compute_psf!`
- `psf_pixel_scale_arcsec`
- `ensure_psf_state!`
- `ZernikeBasis`
- `compute_zernike!`
- `noll_to_nm`
- `OPDMap`
- `NCPA`
- `NCPABasis`
- `KLBasis`
- `ZernikeModalBasis`
- `M2CBasis`
- `default_ncpa_basis`
- `SpatialFilter`
- `SpatialFilterShape`
- `CircularFilter`
- `SquareFilter`
- `FoucaultFilter`
- `set_spatial_filter!`
- `filter!`

### Group 10: Atmosphere evolution and phase statistics

Proposed tier: `Stable`

Symbols:

- `KolmogorovAtmosphere`
- `KolmogorovParams`
- `KolmogorovState`
- `InfinitePhaseScreenParams`
- `InfinitePhaseScreenState`
- `InfinitePhaseScreen`
- `InfiniteLayerParams`
- `InfiniteLayerState`
- `InfiniteAtmosphereLayer`
- `InfiniteMultiLayerParams`
- `InfiniteMultiLayerState`
- `InfiniteMultiLayerAtmosphere`
- `update_psd!`
- `ensure_psd!`
- `phase_screen_von_karman!`
- `MultiLayerAtmosphere`
- `MultiLayerParams`
- `MultiLayerState`
- `advance!`
- `propagate!`
- `phase_variance`
- `phase_covariance`
- `phase_spectrum`
- `covariance_matrix`
- `ft_phase_screen`
- `ft_sh_phase_screen`
- `PhaseStatsWorkspace`
- `SubharmonicMode`
- `FastSubharmonics`
- `FidelitySubharmonics`
- `default_subharmonic_mode`

### Group 11: Deformable mirror and misregistration

Proposed tier: `Stable`

Symbols:

- `DeformableMirror`
- `DeformableMirrorParams`
- `DeformableMirrorState`
- `build_influence_functions!`
- `apply!`
- `Misregistration`
- `apply_misregistration`
- `rotation_rad`
- `rotation_deg`
- `anamorphosis_angle_rad`
- `anamorphosis_angle_deg`

### Group 12: Detector core families

Proposed tier: `Advanced`

Rationale:

- detector modeling is a major package feature
- but the type surface is broad and more expert-facing than the simple
  `capture!` workflow itself

Symbols:

- `AbstractFrameDetector`
- `AbstractCountingDetector`
- `Detector`
- `DetectorParams`
- `DetectorState`
- `DetectorExportMetadata`
- `AbstractDetectorResponse`
- `AbstractFrameResponse`
- `AbstractFrameMTF`
- `FrameResponseModel`
- `NullFrameResponse`
- `GaussianPixelResponse`
- `SampledFrameResponse`
- `RectangularPixelAperture`
- `SeparablePixelMTF`
- `SeparableGaussianPixelResponse`
- `AbstractDetectorDefectModel`
- `NullDetectorDefectModel`
- `PixelResponseNonuniformity`
- `DarkSignalNonuniformity`
- `BadPixelMask`
- `CompositeDetectorDefectModel`
- `FrameWindow`
- `AbstractFrameTimingModel`
- `GlobalShutter`
- `RollingShutter`
- `FrameSamplingMode`
- `SingleRead`
- `AveragedNonDestructiveReads`
- `CorrelatedDoubleSampling`
- `FowlerSampling`
- `FrameReadoutCorrectionModel`
- `NullFrameReadoutCorrection`
- `ReferencePixelCommonModeCorrection`
- `ReferenceRowCommonModeCorrection`
- `ReferenceColumnCommonModeCorrection`
- `ReferenceOutputCommonModeCorrection`
- `CompositeFrameReadoutCorrection`
- `FrameReadoutProducts`
- `NoFrameReadoutProducts`
- `SampledFrameReadoutProducts`
- `HgCdTeReadoutProducts`
- `AbstractFrameNonlinearityModel`
- `NullFrameNonlinearity`
- `SaturatingFrameNonlinearity`
- `AbstractPersistenceModel`
- `NullPersistence`
- `ExponentialPersistence`
- `AbstractDetectorThermalModel`
- `NullDetectorThermalModel`
- `FixedTemperature`
- `FirstOrderThermalModel`
- `AbstractDetectorThermalState`
- `NoThermalState`
- `DetectorThermalState`
- `AbstractTemperatureLaw`
- `NullTemperatureLaw`
- `ArrheniusRateLaw`
- `LinearTemperatureLaw`
- `ExponentialTemperatureLaw`
- `APDDetector`
- `APDDetectorParams`
- `APDDetectorState`
- `CountingReadoutMetadata`
- `CountingDetectorExportMetadata`
- `CountingDeadTimeModel`
- `NoDeadTime`
- `NonParalyzableDeadTime`
- `ParalyzableDeadTime`
- `AbstractCountingGateModel`
- `NullCountingGate`
- `DutyCycleGate`
- `AbstractCountingCorrelationModel`
- `NullCountingCorrelation`
- `AfterpulsingModel`
- `ChannelCrosstalkModel`
- `CompositeCountingCorrelation`

### Group 13: Detector capture and detector metadata accessors

Proposed tier: `Advanced`

Symbols:

- `capture!`
- `output_frame`
- `readout_products`
- `detector_reference_frame`
- `detector_signal_frame`
- `detector_combined_frame`
- `detector_reference_cube`
- `detector_signal_cube`
- `detector_read_cube`
- `detector_read_times`
- `channel_output`
- `detector_export_metadata`
- `readout_ready`
- `reset_integration!`

### Group 14: Detector capability predicates and readout/thermal helpers

Proposed tier: `De-export candidate`

Rationale:

- useful for power users and extensions
- too broad and low-level for the default namespace

Symbols:

- `supports_detector_mtf`
- `supports_clock_induced_charge`
- `supports_column_readout_noise`
- `supports_avalanche_gain`
- `supports_sensor_glow`
- `supports_nondestructive_reads`
- `supports_reference_read_subtraction`
- `supports_readout_correction`
- `supports_read_cube`
- `supports_detector_defect_maps`
- `supports_detector_persistence`
- `supports_detector_nonlinearity`
- `supports_shutter_timing`
- `supports_counting_noise`
- `supports_dead_time`
- `supports_channel_gain_map`
- `supports_counting_gating`
- `supports_afterpulsing`
- `supports_channel_crosstalk`
- `supports_paralyzable_dead_time`
- `supports_detector_thermal_model`
- `supports_dynamic_thermal_state`
- `supports_temperature_dependent_dark_current`
- `supports_temperature_dependent_glow`
- `supports_temperature_dependent_persistence`
- `supports_temperature_dependent_dark_counts`
- `response_family`
- `response_application_domain`
- `response_support`
- `is_shift_invariant`
- `supports_frequency_domain_application`
- `supports_separable_application`
- `supports_subpixel_geometry`
- `default_response_model`
- `thermal_model`
- `thermal_state`
- `detector_temperature`
- `advance_thermal!`
- `evaluate_temperature_law`
- `effective_dark_current`
- `effective_glow_rate`
- `effective_cic_rate`
- `effective_dark_count_rate`

### Group 15: Sensor family and detector component types

Proposed tier: `Advanced`

Symbols:

- `SensorType`
- `FrameSensorType`
- `CountingSensorType`
- `AvalancheFrameSensorType`
- `HgCdTeAvalancheArraySensorType`
- `CCDSensor`
- `CMOSSensor`
- `EMCCDSensor`
- `InGaAsSensor`
- `HgCdTeAvalancheArraySensor`
- `APDSensor`
- `AbstractCMOSOutputModel`
- `NullCMOSOutputModel`
- `StaticCMOSOutputPattern`
- `AbstractEMGainModel`
- `ExcessNoiseApproximation`
- `StochasticMultiplicationRegister`

### Group 16: Wavefront sensor families

Proposed tier: `Stable`

Symbols:

- `ShackHartmannWFS`
- `ShackHartmannWFSParams`
- `ShackHartmannWFSState`
- `update_valid_mask!`
- `measure!`
- `PyramidWFS`
- `PyramidParams`
- `PyramidState`
- `pyramid_modulation_frame!`
- `BioEdgeWFS`
- `BioEdgeParams`
- `BioEdgeState`
- `ZernikeWFS`
- `ZernikeWFSParams`
- `ZernikeWFSState`
- `CurvatureReadoutModel`
- `CurvatureFrameReadout`
- `CurvatureCountingReadout`
- `CurvatureBranchResponse`
- `CurvatureWFS`
- `CurvatureWFSParams`
- `CurvatureWFSState`
- `ensure_curvature_calibration!`
- `apply_shift_wfs!`
- `set_optical_gain!`

### Group 17: LiFT, calibration, reconstructor, and AO workflow builders

Proposed tier: `Advanced`

Symbols:

- `LiFT`
- `lift_interaction_matrix`
- `lift_interaction_matrix!`
- `LiFTSolveMode`
- `LiFTSolveAuto`
- `LiFTSolveQR`
- `LiFTSolveNormalEquations`
- `LiFTDampingMode`
- `LiFTDampingNone`
- `LiFTLevenbergMarquardt`
- `LiFTAdaptiveLevenbergMarquardt`
- `LiFTDiagnostics`
- `diagnostics`
- `InteractionMatrix`
- `interaction_matrix`
- `ControlMatrix`
- `with_truncation`
- `ModalBasis`
- `KLBasisMethod`
- `KLDMModes`
- `KLHHtPSD`
- `dm_basis`
- `kl_modal_basis`
- `modal_basis`
- `basis_from_m2c`
- `basis_projector`
- `AOCalibration`
- `ao_calibration`
- `fitting_error`
- `fitting_error_dm`
- `fast_atmosphere`
- `AOSimulation`
- `initialize_ao_pyramid`
- `initialize_ao_shack_hartmann`
- `GainSensingCamera`
- `GSCDetectorMetadata`
- `detector_metadata`
- `weak_mode_mask`
- `attach_detector!`
- `detach_detector!`
- `calibrate!`
- `reset_calibration!`
- `compute_optical_gains!`
- `MetaSensitivity`
- `compute_meta_sensitivity_matrix`
- `estimate_misregistration`
- `SPRINT`
- `estimate!`
- `AbstractReconstructorOperator`
- `ModalReconstructor`
- `MappedReconstructor`
- `reconstruct!`
- `reconstruct`

### Group 18: Control and runtime

Proposed tier: `Advanced`

Symbols:

- `AbstractController`
- `DiscreteIntegratorController`
- `update!`
- `AbstractControlSimulation`
- `AbstractExecutionPolicy`
- `SequentialExecution`
- `ThreadedExecution`
- `BackendStreamExecution`
- `AbstractRuntimeProfile`
- `ScientificRuntimeProfile`
- `HILRuntimeProfile`
- `RuntimeOutputRequirements`
- `runtime_outputs`
- `RuntimeLatencyModel`
- `default_runtime_profile`
- `default_runtime_outputs`
- `VectorDelayLine`
- `shift_delay!`
- `prepare!`
- `prepare_runtime_wfs!`
- `init_execution_state`
- `supports_prepared_runtime`
- `supports_detector_output`
- `supports_stacked_sources`
- `supports_grouped_execution`
- `ClosedLoopRuntime`
- `SimulationInterface`
- `CompositeSimulationInterface`
- `SimulationReadout`
- `sense!`
- `step!`
- `set_command!`
- `snapshot_outputs!`
- `simulation_readout`
- `simulation_slopes`
- `simulation_command`
- `simulation_wfs_frame`
- `simulation_science_frame`
- `simulation_wfs_metadata`
- `simulation_science_metadata`
- `simulation_interface`
- `runtime_profile`
- `runtime_latency`
- `RuntimeTimingStats`
- `RuntimePhaseTimingStats`
- `runtime_timing`
- `runtime_phase_timing`
- `with_reconstructor`
- `with_reconstructors`

### Group 19: Abstract interfaces, traits, and normalization models

Proposed tier: `De-export candidate`

Rationale:

- necessary for extension authors
- too low-level and numerous for the default namespace

Symbols:

- `AbstractOpticalElement`
- `AbstractSource`
- `AbstractAtmosphere`
- `AbstractWFS`
- `AbstractDetector`
- `AbstractDeformableMirror`
- `SensingMode`
- `Diffractive`
- `Geometric`
- `AbstractSlopeExtractionModel`
- `CenterOfGravityExtraction`
- `SubapertureLayout`
- `SubapertureCalibration`
- `subaperture_layout`
- `subaperture_calibration`
- `slope_extraction_model`
- `valid_subaperture_indices`
- `n_valid_subapertures`
- `WFSNormalization`
- `MeanValidFluxNormalization`
- `IncidenceFluxNormalization`
- `NoiseModel`
- `NoiseNone`
- `NoisePhoton`
- `NoiseReadout`
- `NoisePhotonReadout`
- `DMApplyMode`
- `DMAdditive`
- `DMReplace`
- `ParallelConfig`
- `with_parallel_config`

### Group 20: Backend and GPU policy surface

Proposed tier: `De-export candidate`

Symbols:

- `set_fft_provider_threads!`
- `gpu_backend_loaded`
- `gpu_backend_array_type`
- `gpu_backend_name`
- `available_gpu_backends`
- `disable_scalar_backend!`
- `backend_rand`
- `backend_randn`
- `backend_zeros`
- `backend_fill`
- `GPUBackendTag`
- `CUDABackendTag`
- `MetalBackendTag`
- `AMDGPUBackendTag`
- `GPUPrecisionPolicy`
- `UnifiedGPUPrecision`
- `SplitGPUPrecision`
- `gpu_runtime_type`
- `gpu_build_type`
- `default_gpu_precision_policy`
- `high_accuracy_gpu_precision_policy`

### Group 21: Tomography and slope-order assembly

Proposed tier: `Advanced`

Symbols:

- `TomographyAtmosphereParams`
- `LGSAsterismParams`
- `LGSWFSParams`
- `TomographyParams`
- `TomographyDMParams`
- `TomographyFitting`
- `AbstractTomographyMethod`
- `ModelBasedTomography`
- `InteractionMatrixTomography`
- `TomographyNoiseModel`
- `RelativeSignalNoise`
- `ScalarMeasurementNoise`
- `DiagonalMeasurementNoise`
- `PhotonReadoutSlopeNoise`
- `AbstractSlopeOrder`
- `SimulationSlopes`
- `InterleavedSlopes`
- `InvertedSlopes`
- `TomographyOperators`
- `TomographicReconstructor`
- `TomographyCommandReconstructor`
- `build_reconstructor`
- `assemble_reconstructor_and_fitting`
- `airmass`
- `zenith_angle_rad`
- `zenith_angle_deg`
- `layer_altitude_m`
- `wind_direction_rad`
- `wind_direction_deg`
- `wind_velocity_components`
- `lgs_height_m`
- `lgs_directions!`
- `lgs_directions`
- `optimization_geometry!`
- `optimization_geometry`
- `direction_vectors!`
- `direction_vectors`
- `valid_lenslet_support`
- `dm_valid_support`
- `support_diameter`
- `sparse_gradient_matrix`
- `auto_correlation`
- `cross_correlation`
- `influence_functions`
- `fit_commands!`
- `fit_commands`
- `reconstruct_wavefront!`
- `reconstruct_wavefront`
- `reconstruct_wavefront_map!`
- `reconstruct_wavefront_map`
- `swap_xy_blocks`
- `interleave_xy_columns`
- `prepare_slope_order`
- `mask_actuators!`
- `dm_commands!`
- `dm_commands`

## Known Current Underdocumented-but-Public Areas

These appear to be effectively public today but are not yet well explained as
public-facing surfaces for ordinary users:

- backend policy helpers and GPU capability surfaces
- detector capability predicates
- many runtime readout and metadata accessors
- expert tomography assembly helpers
- calibration workflow builder internals

These should be revisited in Phase 1 before deciding whether they remain
exported.

## Initial Phase 1 Direction

Conservative direction:

- keep `Stable` groups exported by default
- document `Advanced` groups clearly, but do not assume they all remain
  exported
- target `De-export candidate` groups first when trimming the top-level surface

This should reduce namespace breadth substantially without harming common AO
workflow ergonomics.

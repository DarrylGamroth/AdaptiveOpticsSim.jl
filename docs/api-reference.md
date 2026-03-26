# API Reference

This is a guide to the exported public API. It is organized by subsystem rather
than by source file.

## Core types and utilities

- `Workspace`, `ensure_psf_buffers!`
- `FidelityProfile`, `ScientificProfile`, `FastProfile`, `ProfileBundle`,
  `default_fidelity_profile`
- `atmosphere_profile`, `calibration_profile`, `detector_profile`,
  `lift_profile`, `tomography_profile`
- `Telemetry`, `TelemetryRow`, `record!`
- `ClosedLoopTrace`, `ClosedLoopTraceRow`
- `GSCClosedLoopTrace`, `GSCClosedLoopTraceRow`
- `GSCAtmosphereReplayTrace`, `GSCAtmosphereReplayTraceRow`
- `snapshot_config`, `write_config_toml`, `write_config_json`
- `write_telemetry_csv`
- `bin2d`, `poisson_noise!`, `poisson_sample`
- `ParallelConfig`, `with_parallel_config`
- `InversePolicy`, `ExactPseudoInverse`, `TSVDInverse`, `TikhonovInverse`
- Exceptions: `AdaptiveOpticsSimError`, `InvalidConfiguration`,
  `DimensionMismatchError`

## Optical elements

- `Telescope`, `TelescopeParams`, `TelescopeState`
- `generate_pupil!`, `reset_opd!`, `apply_opd!`, `set_pupil!`, `apply_spiders!`
- `Source`, `SourceParams`, `LGSSource`, `LGSSourceParams`, `wavelength`
- `Asterism`, `coordinates_xy_arcsec`, `compute_psf!`,
  `psf_pixel_scale_arcsec`, `ensure_psf_state!`
- `ZernikeBasis`, `compute_zernike!`, `noll_to_nm`
- `Detector`, `DetectorParams`, `DetectorState`, `DetectorExportMetadata`
- `AbstractFrameDetector`, `AbstractCountingDetector`
- `APDDetector`, `APDDetectorParams`, `APDDetectorState`
- `FrameResponseModel`, `NullFrameResponse`,
  `SeparableGaussianPixelResponse`
- `CountingReadoutMetadata`, `CountingDetectorExportMetadata`
- `CountingDeadTimeModel`, `NoDeadTime`, `NonParalyzableDeadTime`
- `capture!`, `output_frame`, `channel_output`, `detector_export_metadata`,
  `readout_ready`, `reset_integration!`
- `supports_detector_mtf`
- `supports_counting_noise`, `supports_dead_time`,
  `supports_channel_gain_map`
- `SensorType`, `FrameSensorType`, `CountingSensorType`, `CCDSensor`,
  `CMOSSensor`, `EMCCDSensor`, `APDSensor`
- `DeformableMirror`, `DeformableMirrorParams`, `DeformableMirrorState`,
  `build_influence_functions!`, `apply!`
- `Misregistration`, `apply_misregistration`
- `OPDMap`
- `NCPA`, `NCPABasis`, `KLBasis`, `ZernikeModalBasis`, `M2CBasis`,
  `default_ncpa_basis`
- `SpatialFilter`, `SpatialFilterShape`, `CircularFilter`, `SquareFilter`,
  `FoucaultFilter`, `set_spatial_filter!`, `filter!`

## Atmosphere

- `KolmogorovAtmosphere`, `KolmogorovParams`, `KolmogorovState`
- `MultiLayerAtmosphere`, `MultiLayerParams`, `MultiLayerState`
- `advance!`, `propagate!`
- `update_psd!`, `ensure_psd!`, `phase_screen_von_karman!`
- `phase_variance`, `phase_covariance`, `phase_spectrum`, `covariance_matrix`
- `ft_phase_screen`, `ft_sh_phase_screen`, `PhaseStatsWorkspace`
- `SubharmonicMode`, `FastSubharmonics`, `FidelitySubharmonics`,
  `default_subharmonic_mode`
- `fast_atmosphere`

## Wavefront sensing

- `ShackHartmann`, `ShackHartmannParams`, `ShackHartmannState`
- `PyramidWFS`, `PyramidParams`, `PyramidState`
- `BioEdgeWFS`, `BioEdgeParams`, `BioEdgeState`
- `ZernikeWFS`, `ZernikeWFSParams`, `ZernikeWFSState`
- `CurvatureReadoutModel`, `CurvatureFrameReadout`, `CurvatureCountingReadout`
- `CurvatureBranchResponse`
- `update_valid_mask!`, `measure!`, `apply_shift_wfs!`, `set_optical_gain!`
- `LiFT`, `lift_interaction_matrix`, `lift_interaction_matrix!`
- `LiFTDampingMode`, `LiFTDampingNone`, `LiFTLevenbergMarquardt`,
  `LiFTAdaptiveLevenbergMarquardt`
- `lgs_elongation_factor`

## Calibration and reconstruction

- `InteractionMatrix`, `interaction_matrix`
- `CalibrationVault`, `with_truncation`
- `ModalBasis`, `KLBasisMethod`, `KLDMModes`, `KLHHtPSD`
- `dm_basis`, `kl_modal_basis`, `modal_basis`, `basis_from_m2c`,
  `basis_projector`
- `AOCalibration`, `ao_calibration`
- `fitting_error`, `fitting_error_dm`
- `GainSensingCamera`, `calibrate!`, `reset_calibration!`,
  `compute_optical_gains!`
- `MetaSensitivity`, `compute_meta_sensitivity_matrix`,
  `estimate_misregistration`, `SPRINT`, `estimate!`
- `AbstractReconstructorOperator`, `ModalReconstructor`,
  `MappedReconstructor`, `reconstruct!`, `reconstruct`
- `TomographyNoiseModel`, `RelativeSignalNoise`, `ScalarMeasurementNoise`,
  `DiagonalMeasurementNoise`

## Control

- `AbstractController`
- `DiscreteIntegratorController`, `update!`
- `ClosedLoopRuntime`, `SimulationInterface`, `CompositeSimulationInterface`, `SimulationReadout`
- `AbstractControlSimulation`
- `prepare!`, `prepare_runtime_wfs!`, `simulation_interface`
- `simulation_readout`, `simulation_command`, `simulation_slopes`,
  `simulation_wfs_frame`, `simulation_science_frame`
- `simulation_wfs_metadata`, `simulation_science_metadata`
  - these accessors now work on both interface/readout objects and direct
    `AbstractControlSimulation` instances
- `AbstractExecutionPolicy`, `SequentialExecution`, `ThreadedExecution`,
  `BackendStreamExecution`
- `VectorDelayLine`, `shift_delay!`
- `supports_prepared_runtime`, `supports_detector_output`,
  `supports_stacked_sources`, `supports_grouped_execution`

## Traits and interfaces

- Abstract interfaces: `AbstractOpticalElement`, `AbstractSource`,
  `AbstractAtmosphere`, `AbstractWFS`, `AbstractDetector`,
  `AbstractDeformableMirror`
- Sensing-mode traits: `SensingMode`, `Diffractive`, `Geometric`
- Detector/readout traits: `NoiseModel`, `NoiseNone`, `NoisePhoton`,
  `NoiseReadout`, `NoisePhotonReadout`, `SensorType`, `FrameSensorType`,
  `CountingSensorType`
- DM-application traits: `DMApplyMode`, `DMAdditive`, `DMReplace`

## Interface contracts

These contracts summarize the maintained method families that extension code is
expected to implement. The conformance evidence for the maintained surfaces
lives in the `Interface conformance` testset in `test/runtests.jl`.

### `IF-OPT`: optical elements

- `AbstractOpticalElement` implementations that modify telescope phase should
  implement `apply!(element, tel, mode)`.
- Maintained optical-phase elements support both `DMAdditive()` and
  `DMReplace()` semantics.
- Shape mismatches should raise structured errors rather than silently resizing
  or clipping.

### `IF-SRC`: sources

- `AbstractSource` implementations must provide `wavelength(src)`.
- Source-specific flux or spectrum helpers may exist, but `wavelength` is the
  minimal shared contract used across optics and sensing.

### `IF-ATM`: atmospheres

- `AbstractAtmosphere` implementations must provide `advance!(atm, tel; rng)`
  and `propagate!(atm, tel)`.
- `advance!` updates the evolving atmospheric state, while `propagate!`
  applies the resulting OPD to the telescope.

### `IF-WFS`: wavefront sensors

- `AbstractWFS` implementations must provide `update_valid_mask!(wfs, tel)` and
  `measure!(wfs, tel)` or `measure!(wfs, tel, src)`.
- Optional capabilities are surfaced through traits and behavior:
  detector-coupled output, asterism support, prepared runtime support, and
  stacked-source execution.
- `CurvatureWFS` now has two maintained readout families:
  `CurvatureFrameReadout` for frame-style detector coupling and
  `CurvatureCountingReadout` for counting/channel-style readout.
- `CurvatureBranchResponse` models reusable intra-/extra-focal throughput and
  background imbalance independently of any one instrument example.
- `CurvatureWFS` also exposes detector-plane geometry controls through
  `readout_crop_resolution` and `readout_pixels_per_subap` so the exported frame
  sampling can be decoupled from the control-grid sampling.
- Maintained optional-capability queries use
  `supports_prepared_runtime(wfs, src)` and `supports_stacked_sources(wfs, src)`
  in addition to the simulation-level trait surface.
- The maintained runtime expects the measured slope vector to live in
  `wfs.state.slopes`.

### `IF-DET`: detectors

- `AbstractDetector` implementations must provide `capture!(det, psf; rng)`.
- Export-facing code relies on `output_frame(det)` and
  `detector_export_metadata(det)` when detector-coupled outputs are present.
- The maintained `Detector` type is a frame-detector implementation and accepts
  `FrameSensorType` sensors only. Counting sensors such as `APDSensor` are a
  distinct readout family and should be modeled through sensor/readout-specific
  code rather than the generic frame-detector path.
- The maintained frame-detector response family is opt-in and null by default:
  `NullFrameResponse` is the identity model and
  `SeparableGaussianPixelResponse` is the maintained first pixel-response/MTF
  approximation.
- The maintained counting-detector family is currently `APDDetector`, with
  optional capability queries surfaced through `supports_counting_noise`,
  `supports_dead_time`, and `supports_channel_gain_map`.
- Counting-detector dead-time behavior is selected by dispatch through
  `CountingDeadTimeModel`, with `NoDeadTime` as the null model and
  `NonParalyzableDeadTime` as the maintained first physical model.

### `IF-DM`: deformable mirrors

- `AbstractDeformableMirror` implementations must provide
  `build_influence_functions!(dm, tel)` and `apply!(dm, tel, mode)`.
- The maintained DM runtime path assumes `dm.state.coefs` is the active command
  vector and that `apply!` writes phase on the telescope OPD grid.

### `IF-REC`: control reconstructors

- `AbstractReconstructorOperator` is the maintained core slopes-to-command
  operator family, and the required runtime contract is
  `reconstruct!(out, recon, slopes)`.
- `reconstruct(recon, slopes)` is the thin allocation wrapper over that
  mutating operator.
- The maintained control reconstructor families are `ModalReconstructor` and
  `MappedReconstructor`, and both must present the same external
  slopes-to-command contract.

### `IF-CTRL`: controllers

- `AbstractController` implementations should expose a canonical mutating update
  path through `update!(ctrl, input, dt)`.
- The maintained controller family is `DiscreteIntegratorController`, which
  returns the updated command-like controller state.

### `IF-SIM`: control simulations

- `AbstractControlSimulation` implementations must provide `step!(sim)`,
  `simulation_interface(sim)`, and `simulation_readout(sim)`.
- `prepare!(sim)` is the preferred optional hook for runtime precomputation.
- Output-side accessors are `simulation_command`, `simulation_slopes`,
  `simulation_wfs_frame`, `simulation_science_frame`,
  `simulation_wfs_metadata`, and `simulation_science_metadata`.
- Optional behavior is expressed with
  `supports_prepared_runtime`, `supports_detector_output`,
  `supports_stacked_sources`, and `supports_grouped_execution`.

### `IF-CAL`: calibration workflow surfaces

- `interaction_matrix(...)` is the maintained linearized WFS-calibration entry
  point and must return an `InteractionMatrix` with a stored calibration
  `amplitude`.
- `CalibrationVault(D; ...)` is the maintained inverse-storage workflow and
  should retain the forward operator `D`, the selected inverse representation
  when requested, and inversion diagnostics such as singular values, condition
  number, effective rank, and truncation count.
- `modal_basis(...)` is the maintained modal-basis workflow and must return a
  `ModalBasis` with consistent `M2C`, sampled basis vectors, and optional
  projector.
- `ao_calibration(...)` is the maintained packaged AO-calibration workflow and
  should return an `AOCalibration` bundling basis operators with the measured
  `CalibrationVault`.
- `compute_meta_sensitivity_matrix(...)` and `SPRINT(...)` are the maintained
  misregistration-identification workflows and should retain the reference
  calibration, inverted meta-sensitivity operator, finite-difference step
  sizes, and active field ordering.
- `LiFT` reconstruction uses a separate iterative contract:
  `reconstruct!(coeffs, lift, psf_in, mode_ids, ...)` must mutate the supplied
  coefficient buffer in place and retain convergence diagnostics in
  `diagnostics(lift)`.

## Notes on backend-aware APIs

- Constructors typically accept `T=` and `backend=` keywords.
- Methods accept `AbstractArray` inputs where practical, but structs keep
  concrete array fields so the compiler can specialize.
- FFT-heavy paths use `AbstractFFTs.jl`; future GPU-specialized kernels should
  be introduced behind trait-selected backends rather than by branching on
  symbols in hot loops.

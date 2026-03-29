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
- `AbstractRuntimeProfile`, `ScientificRuntimeProfile`, `HILRuntimeProfile`,
  `RuntimeLatencyModel`, `default_runtime_profile`
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
- `FrameWindow`
- `AbstractDetectorThermalModel`, `NullDetectorThermalModel`,
  `FixedTemperature`, `FirstOrderThermalModel`
- `AbstractDetectorThermalState`, `NoThermalState`, `DetectorThermalState`
- `AbstractTemperatureLaw`, `NullTemperatureLaw`, `ArrheniusRateLaw`,
  `LinearTemperatureLaw`, `ExponentialTemperatureLaw`
- `AbstractDetectorResponse`, `AbstractFrameResponse`, `AbstractFrameMTF`
- `FrameResponseModel`, `NullFrameResponse`, `GaussianPixelResponse`,
  `SampledFrameResponse`, `RectangularPixelAperture`, `SeparablePixelMTF`
- `FrameSamplingMode`, `SingleRead`, `AveragedNonDestructiveReads`,
  `CorrelatedDoubleSampling`, `FowlerSampling`
- `CountingReadoutMetadata`, `CountingDetectorExportMetadata`
- `CountingDeadTimeModel`, `NoDeadTime`, `NonParalyzableDeadTime`
- `capture!`, `output_frame`, `channel_output`, `detector_export_metadata`,
  `readout_ready`, `reset_integration!`
- `response_family`, `response_application_domain`, `response_support`
- `default_response_model`
- `thermal_model`, `thermal_state`, `detector_temperature`,
  `advance_thermal!`
- `evaluate_temperature_law`, `effective_dark_current`,
  `effective_glow_rate`, `effective_cic_rate`,
  `effective_dark_count_rate`
- `supports_detector_mtf`, `is_shift_invariant`,
  `supports_frequency_domain_application`,
  `supports_separable_application`, `supports_subpixel_geometry`
- `supports_clock_induced_charge`, `supports_column_readout_noise`
- `supports_avalanche_gain`, `supports_sensor_glow`,
  `supports_nondestructive_reads`, `supports_reference_read_subtraction`
- `supports_counting_noise`, `supports_dead_time`,
  `supports_channel_gain_map`
- `supports_detector_thermal_model`, `supports_dynamic_thermal_state`,
  `supports_temperature_dependent_dark_current`,
  `supports_temperature_dependent_glow`,
  `supports_temperature_dependent_persistence`,
  `supports_temperature_dependent_dark_counts`
- `SensorType`, `FrameSensorType`, `CountingSensorType`, `CCDSensor`,
  `CMOSSensor`, `AvalancheFrameSensorType`,
  `HgCdTeAvalancheArraySensorType`, `EMCCDSensor`, `InGaAsSensor`,
  `HgCdTeAvalancheArraySensor`, `APDSensor`
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
- `runtime_profile`, `runtime_latency`
- `simulation_readout`, `simulation_command`, `simulation_slopes`,
  `simulation_wfs_frame`, `simulation_science_frame`
- `simulation_wfs_metadata`, `simulation_science_metadata`
  - these accessors now work on both interface/readout objects and direct
    `AbstractControlSimulation` instances
- `AbstractExecutionPolicy`, `SequentialExecution`, `ThreadedExecution`,
  `BackendStreamExecution`
- `VectorDelayLine`, `shift_delay!`
- `runtime_timing`, `runtime_phase_timing`
  - `RuntimePhaseTimingStats` now exposes `delay_mean_ns` in addition to the
    existing sense/reconstruct/apply/snapshot timing surface
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
- The maintained multilayer backends also share an internal layer/container
  contract used by the common source-aware rendering path:
  - layer objects implement `sample_layer!`, `sample_layer_accumulate!`,
    `render_layer!`, `render_layer_accumulate!`, and `layer_altitude`
  - container objects provide `layers`, `params.altitude`, `state.opd`, and
    `state.source_geometry`
- Source-aware maintained multilayer atmospheres provide `propagate!(atm, tel,
  src)` and are expected to route through the shared accumulation helpers
  rather than duplicating backend-specific render loops.

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
  `NullFrameResponse` is the identity model, `GaussianPixelResponse` is the
  maintained effective blur-like response, `SampledFrameResponse` is the
  maintained measured/sampled kernel path, `RectangularPixelAperture` is the
  first explicit pixel-geometry model, and `SeparablePixelMTF` is the first
  maintained MTF-specified response family.
- Omitted `response_model` values now resolve by detector family through
  `default_response_model(sensor; ...)`:
  `CCDSensor` and `EMCCDSensor` keep the null response by default,
  `CMOSSensor` and `InGaAsSensor` use a mild Gaussian response,
  and `HgCdTeAvalancheArraySensor` uses a mild sampled-kernel response.
- Detector export metadata now records response family, application domain,
  separability, shift invariance, support, pitch, fill factor, and aperture
  shape, rather than only one response-width scalar.
- Detector thermal behavior is now a separate reusable sublayer:
  `NullDetectorThermalModel` is the default identity path,
  `FixedTemperature` is the first static physical model, and
  `FirstOrderThermalModel` is the first maintained dynamic model.
- Temperature-law infrastructure is shared across detector families through
  `NullTemperatureLaw`, `ArrheniusRateLaw`, `LinearTemperatureLaw`, and
  `ExponentialTemperatureLaw`.
- The maintained thermal accessors are `thermal_model(det)`,
  `thermal_state(det)`, `detector_temperature(det)`, and
  `advance_thermal!(det, dt)`.
- Dynamic thermal evolution is currently modeled as first-order relaxation
  toward a cooling setpoint with explicit initial temperature, ambient
  temperature, time constant, and clamp bounds.
- The maintained temperature-aware detector hooks are
  `effective_dark_current(det)`, `effective_glow_rate(det)`,
  `effective_cic_rate(det)`, and `effective_dark_count_rate(det)`.
- Detector export metadata now also records thermal model family,
  detector temperature, cooling setpoint, and the active temperature-law
  family for dark current, glow, CIC, or dark counts where applicable.
- `FrameWindow(rows, cols)` is a generic frame-readout crop that applies to the
  detector output surface after detector sampling/binning. This is intended for
  subarray/windowed readout and is not detector-family-specific.
- Frame-detector readout correction is now a separate surface from detector
  response. `NullFrameReadoutCorrection` is the null model,
  `ReferencePixelCommonModeCorrection` subtracts one global reference-edge
  bias, `ReferenceRowCommonModeCorrection` and
  `ReferenceColumnCommonModeCorrection` subtract per-row or per-column
  reference-edge bias, `ReferenceOutputCommonModeCorrection` subtracts
  per-output-group bias over column groups, and
  `CompositeFrameReadoutCorrection` composes multiple correction stages.
- Detector export metadata now records correction family, edge support, output
  grouping, and correction stage count.
- Batched detector capture now has two maintained surfaces:
  - the fixed-shape in-place fast path `capture_stack!(det, cube, scratch)`
  - the generalized shape-changing path `capture_stack!(det, out_cube, in_cube)`
- The fixed-shape batched fast path is intended for stacked SH/HIL workloads and
  currently requires:
  - `psf_sampling == 1`
  - `binning == 1`
  - `output_precision === nothing`
  - full-frame readout
  - global-shutter timing
  - null persistence
- The fixed-shape fast path now supports the maintained batched frame-response
  family:
  - `NullFrameResponse`
  - `GaussianPixelResponse`
  - `SampledFrameResponse`
  - `RectangularPixelAperture`
  - `SeparablePixelMTF`
- The fixed-shape fast path also supports the maintained batched readout
  correction family:
  - `NullFrameReadoutCorrection`
  - `ReferencePixelCommonModeCorrection`
  - `ReferenceRowCommonModeCorrection`
  - `ReferenceColumnCommonModeCorrection`
  - `ReferenceOutputCommonModeCorrection`
  - `CompositeFrameReadoutCorrection` when all stages are maintained batched
    corrections
- The generalized batched path preserves detector semantics for
  `psf_sampling`, `binning`, `readout_window`, and `output_precision` by using
  separate input/output stacks. It is more general but not the latency-critical
  fast path.
- Static frame-detector structure is now explicit. The maintained defect-model
  surface includes `NullDetectorDefectModel`,
  `PixelResponseNonuniformity`, `DarkSignalNonuniformity`,
  `BadPixelMask`, and `CompositeDetectorDefectModel`.
- Frame-detector timing is now explicit through `AbstractFrameTimingModel`,
  with `GlobalShutter` as the null/default semantics and `RollingShutter` as
  the maintained timed row-readout model.
- Frame-detector nonlinearity is now explicit through
  `AbstractFrameNonlinearityModel`, with `NullFrameNonlinearity` and
  `SaturatingFrameNonlinearity` as the maintained surface.
- `CCDSensor` supports opt-in clock-induced charge through its constructor.
- `CMOSSensor` supports opt-in column readout noise, grouped-output gain/offset
  patterns through `StaticCMOSOutputPattern`, and explicit shutter timing.
- `EMCCDSensor` supports opt-in excess-noise-factor behavior, EMCCD-specific
  CIC through `cic_rate`, and explicit multiplication-model selection through
  `ExcessNoiseApproximation` or `StochasticMultiplicationRegister`.
- `InGaAsSensor` is a frame-detector family with opt-in glow, persistence
  through `ExponentialPersistence`, and support for explicit detector
  nonlinearity.
- `HgCdTeAvalancheArraySensor` is the maintained avalanche-frame-detector family with
  avalanche gain, excess noise, glow-rate, and optional non-destructive-read
  sampling controls.
- `HgCdTeAvalancheArraySensor` sits under the more explicit
  `HgCdTeAvalancheArraySensorType` family rather than only the generic
  avalanche-frame umbrella.
- The frame-sensor capability traits are now explicit:
  `supports_clock_induced_charge`, `supports_column_readout_noise`,
  `supports_avalanche_gain`, `supports_sensor_glow`, and
  `supports_nondestructive_reads`.
- HgCdTe avalanche-array sensors also expose
  `supports_reference_read_subtraction`, which covers read schemes such as
  correlated double sampling and Fowler sampling.
- `HgCdTeAvalancheArraySensor` also uses an avalanche-aware saturation limit when
  `full_well` is set, so incident charge saturates earlier as avalanche gain
  increases.
- `AveragedNonDestructiveReads(n_reads)` is the maintained first HgCdTe-array
  sampling model and reduces the effective additive readout-noise sigma by
  `1 / sqrt(n_reads)` relative to `SingleRead()`.
- `CorrelatedDoubleSampling()` models one pedestal read and one signal read,
  so additive readout noise increases by `sqrt(2)` relative to `SingleRead()`.
- `FowlerSampling(n_pairs)` models paired pedestal/signal averaging and
  reduces additive readout noise to `sqrt(2 / n_pairs)` times the
  `SingleRead()` sigma.
- `HgCdTeAvalancheArraySensor(read_time=...)` now also threads the read cadence into the
  HgCdTe-array metadata and detector physics, so export metadata records
  reference/signal read counts, per-read time, and total wall-clock time, and
  HgCdTe-array dark/glow accumulation includes the configured read overhead.
- HgCdTe-array readout products are now explicit. The generic accessors
  `detector_reference_frame`, `detector_signal_frame`, `detector_combined_frame`,
  `detector_reference_cube`, `detector_signal_cube`, `detector_read_cube`, and
  `detector_read_times` expose the maintained readout-product surface. The
  dedicated `HgCdTeReadoutProducts` payload makes pedestal/signal averages and
  their corresponding cubes first-class rather than inferring them from only the
  combined detector frame.
- HgCdTe-array subarray timing now scales with active readout rows rather than raw
  cropped area, which better matches row-wise frame-array readout semantics.
- The maintained counting-detector family is currently `APDDetector`, with
  optional capability queries surfaced through `supports_counting_noise`,
  `supports_dead_time`, `supports_channel_gain_map`,
  `supports_counting_gating`, `supports_afterpulsing`,
  `supports_channel_crosstalk`, and `supports_paralyzable_dead_time`.
- Counting-detector dead-time behavior is selected by dispatch through
  `CountingDeadTimeModel`, with `NoDeadTime` as the null model and
  `NonParalyzableDeadTime` and `ParalyzableDeadTime` as maintained physical
  models.
- APD counting control is now explicit through `AbstractCountingGateModel`,
  with `NullCountingGate` and `DutyCycleGate`.
- APD counting correlation is now explicit through
  `AbstractCountingCorrelationModel`, with `NullCountingCorrelation`,
  `AfterpulsingModel`, `ChannelCrosstalkModel`, and
  `CompositeCountingCorrelation`.

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

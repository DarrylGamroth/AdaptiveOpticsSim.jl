# API Reference

Status: active

Related guides:

- [`documentation-map.md`](./documentation-map.md)
- [`user-guide.md`](./user-guide.md)
- [`maintainer-architecture.md`](./maintainer-architecture.md)
- [`runtime-dataflow.md`](./runtime-dataflow.md)

This is a guide to the maintained public API. It is organized by subsystem
rather than by source file.

Project-wide units guidance lives in
[units-policy.md](./units-policy.md).

## Public API policy

The package now uses a tiered API surface:

- `Stable top-level API`
  - exported by default
  - intended for ordinary optics, sensing, calibration, and runtime workflows
- `Advanced documented API`
  - maintained and supported
  - may require qualification as `AdaptiveOpticsSim.<name>`
- `Developer / backend support API`
  - primarily for extensions, benchmark tooling, and backend validation
  - some of this surface remains exported for compatibility, but it should not
    be treated as entry-level workflow API

The practical rule is:

- start with the stable exported workflow surface
- use `AdaptiveOpticsSim.<name>` for advanced helpers such as telemetry,
  scenario builders, build/backend policy utilities, and specialized
  calibration-identification utilities

If you are new to the package, read [`user-guide.md`](./user-guide.md) first.
If you are maintaining the package, pair this document with
[`maintainer-architecture.md`](./maintainer-architecture.md).

## Core types and utilities

- `FidelityProfile`, `ScientificProfile`, `FastProfile`, `ProfileBundle`,
  `default_fidelity_profile`
- `atmosphere_profile`, `calibration_profile`, `detector_profile`,
  `lift_profile`, `tomography_profile`
- `ParallelConfig`, `with_parallel_config`
- `AbstractRuntimeProfile`, `ScientificRuntimeProfile`, `HILRuntimeProfile`,
  `RuntimeLatencyModel`, `default_runtime_profile`
- `InversePolicy`, `ExactPseudoInverse`, `TSVDInverse`, `TikhonovInverse`
- Exceptions: `AdaptiveOpticsSimError`, `InvalidConfiguration`,
  `DimensionMismatchError`
- Aperture/mask primitives: `MaskGrid`, `CircularAperture`,
  `AnnularAperture`, `SpiderMask`, `RectangularROI`,
  `SubapertureGridMask`, `default_mask_grid`, `pixel_mask_grid`,
  `build_mask!`, `apply_mask!`

## Optical elements

- `Telescope`, `TelescopeParams`, `TelescopeState`
- `generate_pupil!`, `reset_opd!`, `apply_opd!`, `set_pupil!`, `apply_spiders!`
- `Source`, `SourceParams`, `LGSSource`, `LGSSourceParams`, `wavelength`
- `ElectricField`, `ElectricFieldParams`, `ElectricFieldState`
- `AbstractPropagationModel`, `FraunhoferPropagation`, `FresnelPropagation`
- `AbstractAtmosphericFieldModel`, `GeometricAtmosphericPropagation`,
  `LayeredFresnelAtmosphericPropagation`
- `AtmosphericFieldPropagation`, `AtmosphericFieldPropagationParams`,
  `AtmosphericFieldPropagationState`, `AtmosphericFieldSlice`
- `fill_from_telescope!`, `fill_telescope_field!`, `apply_phase!`,
  `apply_amplitude!`, `intensity!`, `propagate_field!`,
  `propagate_atmosphere_field!`, `atmospheric_intensity!`
- `Asterism`, `coordinates_xy_arcsec`, `compute_psf!`,
  `psf_pixel_scale_arcsec`, `ensure_psf_state!`
- `ZernikeBasis`, `compute_zernike!`, `noll_to_nm`
- `Detector`, `DetectorParams`, `DetectorState`, `DetectorExportMetadata`
- `AbstractFrameDetector`, `AbstractCountingDetector`
- `APDDetector`, `APDDetectorParams`, `APDDetectorState`
- `SPADArrayDetector`, `SPADArrayDetectorParams`, `SPADArrayDetectorState`
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
- `FrameReadoutProducts`, `NoFrameReadoutProducts`,
  `SampledFrameReadoutProducts`, `MultiReadFrameReadoutProducts`,
  `HgCdTeReadoutProducts`
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
  `HgCdTeAvalancheArraySensor`, `APDSensor`, `SPADArraySensor`
- `AbstractDMTopology`, `ActuatorGridTopology`, `SampledActuatorTopology`
- `AbstractDMInfluenceModel`, `GaussianInfluenceWidth`,
  `GaussianMechanicalCoupling`, `DenseInfluenceMatrix`,
  `MeasuredInfluenceFunctions`
- `AbstractDMActuatorModel`, `LinearStaticActuators`,
  `ClippedActuators`, `ActuatorHealthMap`, `CompositeDMActuatorModel`
- `DeformableMirror`, `DeformableMirrorParams`, `DeformableMirrorState`,
  `build_influence_functions!`, `apply!`, `topology`, `actuator_model`,
  `topology_axis_count`, `topology_command_count`, `actuator_coordinates`,
  `valid_actuator_mask`, `active_actuator_indices`, `topology_metadata`,
  `influence_model`, `influence_width`, `mechanical_coupling`,
  `influence_width_from_mechanical_coupling`
- `Misregistration`, `apply_misregistration`, `rotation_rad`, `rotation_deg`,
  `anamorphosis_angle_rad`, `anamorphosis_angle_deg`
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

## Wavefront sensing

- `ShackHartmann`, `ShackHartmannParams`, `ShackHartmannState`
- `SubapertureLayout`, `SubapertureCalibration`,
  `AbstractSlopeExtractionModel`, `CenterOfGravityExtraction`
- `subaperture_layout`, `subaperture_calibration`,
  `slope_extraction_model`, `valid_subaperture_indices`,
  `n_valid_subapertures`
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
- `AbstractReconstructorOperator`, `ModalReconstructor`,
  `MappedReconstructor`, `reconstruct!`, `reconstruct`
- `TomographyNoiseModel`, `RelativeSignalNoise`, `ScalarMeasurementNoise`,
  `DiagonalMeasurementNoise`
- `airmass`, `zenith_angle_rad`, `zenith_angle_deg`, `layer_altitude_m`,
  `wind_direction_rad`, `wind_direction_deg`, `wind_velocity_components`

## Control

### Primitive control and sensing layer

- `AbstractController`
- `DiscreteIntegratorController`, `update!`
- `NullReconstructor`
- `prepare_runtime_wfs!`

### Runtime execution layer

- `AbstractControlSimulation`
- `ClosedLoopRuntime`, `SimulationReadout`
- `prepare!`, `sense!`, `step!`
- `set_command!`, `update_command!`, `snapshot_outputs!`
  - `set_command!` and `update_command!` accept flat vectors or structured `NamedTuple` commands
- `readout`, `command`, `slopes`, `wfs_frame`, `science_frame`
- `wfs_metadata`, `science_metadata`
- `grouped_wfs_stack`, `grouped_science_stack`
- `runtime_profile`, `runtime_latency`
- `runtime_timing`, `runtime_phase_timing`
  - `RuntimePhaseTimingStats` now exposes `delay_mean_ns` in addition to the existing sense/reconstruct/apply/snapshot timing surface
- `supports_prepared_runtime`, `supports_detector_output`,
  `supports_stacked_sources`, `supports_grouped_execution`

### Orchestration layer

- `RuntimeBranch`, `SingleRuntimeConfig`, `GroupedRuntimeConfig`, `RuntimeScenario`
- `build_runtime_scenario`
- `RuntimeProductRequirements`, `GroupedRuntimeProductRequirements`
- `platform_config`, `platform_boundary`, `platform_name`, `platform_branch_labels`

This is the preferred public runtime assembly surface for normal closed-loop and
HIL usage.

### Advanced single-runtime wrappers

- `SimulationInterface`, `CompositeSimulationInterface`
- `AdaptiveOpticsSim.simulation_interface`

Use these when you are manually assembling or testing a single runtime and do
not need the scenario/config layer.

### Controllable optics used by the runtime layer

- `AbstractControllableOptic`, `CompositeControllableOptic`
- `AbstractModalOpticBasis`
- `FunctionModalBasis`, `MatrixModalBasis`
- `ZernikeOpticBasis`, `CartesianTiltBasis`, `QuadraticFocusBasis`
- `ModalControllableOptic`, `TipTiltMirror`, `FocusStage`, `DeformableMirror`

Preferred construction style:

- `ModalControllableOptic(tel, basis_spec; ...)`
- `TipTiltMirror(...)` and `FocusStage(...)` remain convenience constructors over
  the basis-spec path

### Execution and scheduling helpers

- `AbstractExecutionPolicy`, `SequentialExecution`, `ThreadedExecution`,
  `BackendStreamExecution`
- `VectorDelayLine`, `shift_delay!`

## Advanced documented API

These surfaces are maintained, but Phase 1 moved them out of the default
top-level namespace. Access them as `AdaptiveOpticsSim.<name>`.

### Telemetry, config, and trace containers

- `AdaptiveOpticsSim.PROJECT_STATUS`
- `AdaptiveOpticsSim.Workspace`, `AdaptiveOpticsSim.ensure_psf_buffers!`
- `AdaptiveOpticsSim.Telemetry`, `AdaptiveOpticsSim.TelemetryRow`,
  `AdaptiveOpticsSim.record!`
- `AdaptiveOpticsSim.ClosedLoopTrace`,
  `AdaptiveOpticsSim.ClosedLoopTraceRow`
- `AdaptiveOpticsSim.GSCClosedLoopTrace`,
  `AdaptiveOpticsSim.GSCClosedLoopTraceRow`
- `AdaptiveOpticsSim.GSCAtmosphereReplayTrace`,
  `AdaptiveOpticsSim.GSCAtmosphereReplayTraceRow`
- `AdaptiveOpticsSim.snapshot_config`,
  `AdaptiveOpticsSim.write_config_toml`,
  `AdaptiveOpticsSim.write_config_json`,
  `AdaptiveOpticsSim.write_telemetry_csv`
- `AdaptiveOpticsSim.bin2d`, `AdaptiveOpticsSim.poisson_noise!`,
  `AdaptiveOpticsSim.poisson_sample`

### Scenario builders and calibration workflow helpers

- `AdaptiveOpticsSim.fast_atmosphere`
- `AdaptiveOpticsSim.initialize_ao_shwfs`
- `AdaptiveOpticsSim.initialize_ao_pyramid`
- `AdaptiveOpticsSim.GSCDetectorMetadata`
- `AdaptiveOpticsSim.detector_metadata`
- `AdaptiveOpticsSim.weak_mode_mask`
- `AdaptiveOpticsSim.attach_detector!`
- `AdaptiveOpticsSim.detach_detector!`
- `AdaptiveOpticsSim.MetaSensitivity`
- `AdaptiveOpticsSim.compute_meta_sensitivity_matrix`
- `AdaptiveOpticsSim.estimate_misregistration`
- `AdaptiveOpticsSim.SPRINT`
- `AdaptiveOpticsSim.estimate!`

### Build and backend policy helpers

- `AdaptiveOpticsSim.BuildBackend`
- `AdaptiveOpticsSim.NativeBuildBackend`
- `AdaptiveOpticsSim.CPUBuildBackend`
- `AdaptiveOpticsSim.GPUArrayBuildBackend`
- `AdaptiveOpticsSim.default_build_backend`
- `AdaptiveOpticsSim.set_fft_provider_threads!`

These remain maintained, but they are now explicitly treated as advanced
developer-facing infrastructure rather than normal top-level workflow API.

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
- The maintained wave-optics field surface is `ElectricField`, which is a
  derived propagation object rather than the primary telescope plant state.
- `ElectricField` implementations/builders are expected to provide:
  `fill_from_telescope!`, `apply_phase!`, `apply_amplitude!`, and `intensity!`.
- Maintained propagation models implement `propagate_field!(out, field, model)`
  and `propagate_field!(field, model)`.
- The maintained atmosphere-aware field surface is `AtmosphericFieldPropagation`
  with `propagate_atmosphere_field!` and `atmospheric_intensity!`.
- `FraunhoferPropagation` exports centered pupil-to-focal field propagation with
  `params.output_sampling_rad` describing the angular sample pitch.
- `FresnelPropagation` exports monochromatic transfer-function propagation with
  explicit `params.distance_m` and `params.output_sampling_m`.
- `GeometricAtmosphericPropagation` is the fast phase-only coupled atmosphere
  path, while `LayeredFresnelAtmosphericPropagation` is the maintained
  layer-aware physical propagation path.
- Runtime OPD remains in meters on `Telescope.state.opd`; field-phase
  conversion is explicit and wavelength-aware at the field boundary.

### `IF-SRC`: sources

- `AbstractSource` implementations must provide `wavelength(src)`.
- Source-specific flux or spectrum helpers may exist, but `wavelength` is the
  minimal shared contract used across optics and sensing.
- The maintained spectral wrapper family is `SpectralSample`,
  `SpectralBundle`, and `SpectralSource`, constructed with `with_spectrum`.
- `SpectralBundle` normalizes weights at construction, while `SpectralSource`
  preserves the wrapped source geometry and flux boundary.
- Source-measurement code that supports grouped broad-band execution should use
  `spectral_reference_source(src)` for geometry preparation and
  `source_measurement_signature(src)` for calibration/cache identity.
- The maintained extended-source wrapper family is `PointCloudSourceModel`,
  `GaussianDiskSourceModel`, `SampledImageSourceModel`, and `ExtendedSource`,
  constructed with `with_extended_source`.
- `ExtendedSource` preserves the wrapped source center and total flux while
  expanding the source distribution into weighted measurement samples only in
  the grouped sensing path.

### `IF-MASK`: aperture and support maps

- Aperture and support-map construction is centralized through
  `MaskGrid`, `CircularAperture`, `AnnularAperture`, `SpiderMask`,
  `RectangularROI`, and `SubapertureGridMask`.
- `build_mask!` is the maintained constructor for full-support boolean,
  weighted-real, or weighted-complex masks.
- `apply_mask!` is the maintained mutating path for destructive occluders such
  as spiders.
- Telescope pupil generation, spatial-filter support maps, and WFS valid-mask
  assembly should reuse these builders instead of introducing local geometry
  loops or kernels.

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
- The maintained coupled atmosphere/field surface is expected to consume those
  same source-aware layer helpers rather than introducing a second backend
  split for finite vs infinite atmospheres.

### `IF-WFS`: wavefront sensors

- `AbstractWFS` implementations must provide `update_valid_mask!(wfs, tel)` and
  `measure!(wfs, tel)` or `measure!(wfs, tel, src)`.
- `slopes(wfs)` is the maintained accessor for the exported 1-D WFS signal
  vector.
- Optional exported WFS-state surfaces are accessed through
  `valid_subaperture_mask(wfs)`, `reference_signal(wfs)`, and
  `camera_frame(wfs)`.
- HIL/RTC-facing WFS detector pixels are accessed through
  `wfs_detector_image(wfs)` or `wfs_detector_image(wfs, det)` after
  measurement. For detector-coupled WFSs this returns the configured detector
  output where that is the maintained readout surface. Shack-Hartmann exports
  its detector-coupled lenslet spots as a 2-D mosaic via
  `shack_hartmann_detector_image(...)`.
- Optional WFS-state capabilities are surfaced through
  `supports_valid_subaperture_mask(wfs)`,
  `supports_reference_signal(wfs)`, and
  `supports_camera_frame(wfs)`.
- Optional capabilities are surfaced through traits and behavior:
  detector-coupled output, asterism support, prepared runtime support, and
  stacked-source execution.
- `CurvatureWFS` now has two maintained readout families:
  `CurvatureFrameReadout` for frame-style detector coupling and
  `CurvatureCountingReadout` for counting/channel-style readout.
- `CurvatureWFS` now also has a maintained atmosphere-aware diffractive path
  through `measure!(wfs, tel, src, atm)` and
  `measure!(wfs, tel, ast, atm)` using the shared atmospheric field
  propagation layer.
- `CurvatureBranchResponse` models reusable intra-/extra-focal throughput and
  background imbalance independently of any one instrument example.
- `ShackHartmann` now exposes its maintained subaperture/calibration surface
  through `SubapertureLayout` and `SubapertureCalibration`.
- Grouped-runtime and detector-coupled Shack-Hartmann paths are expected to
  share that layout/calibration contract rather than managing independent
  valid-mask and reference-signal state.
- `CurvatureWFS` also exposes detector-plane geometry controls through
  `readout_crop_resolution` and `readout_pixels_per_subap` so the exported frame
  sampling can be decoupled from the control-grid sampling.
- Maintained optional-capability queries use
  `supports_prepared_runtime(wfs, src)`, `supports_stacked_sources(wfs, src)`,
  `supports_detector_output(wfs, det)`, and
  `supports_grouped_execution(wfs, src)` in addition to the simulation-level
  trait surface.
- The maintained runtime expects `slopes(wfs)` to expose the measured signal
  vector consistently across WFS families.
- Diffractive `ShackHartmann` and `PyramidWFS` now also support grouped
  broad-band execution through `SpectralSource`, including detector-coupled
  readout after wavelength accumulation.
- Diffractive `ShackHartmann` and `PyramidWFS` also support grouped
  extended-source execution through `ExtendedSource`, reusing the same
  maintained grouped accumulation pattern used for asterisms and broad-band
  sensing.

### `IF-DET`: detectors

- `AbstractDetector` implementations must provide `capture!(det, psf; rng)`.
- Export-facing code relies on `output_frame(det)` and
  `detector_export_metadata(det)` when detector-coupled outputs are present.
- The maintained `Detector` type is a frame-detector implementation and accepts
  `FrameSensorType` sensors only. Counting sensors such as `APDSensor` and
  `SPADArraySensor` are a distinct readout family and should be modeled through
  sensor/readout-specific code rather than the generic frame-detector path.
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
- Rolling-shutter timing and detector persistence are intentionally not part of
  either maintained batched-detector path today. They require temporal or
  latent-state semantics that are not a good fit for the current spot-stack HIL
  fast path, so they are deferred until a maintained system surface requires
  them.
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
  shared `MultiReadFrameReadoutProducts` payload is the reusable multi-read
  frame-detector surface, and `HgCdTeReadoutProducts` is maintained as the
  compatibility alias for the current HgCdTe family.
- HgCdTe-array subarray timing now scales with active readout rows rather than raw
  cropped area, which better matches row-wise frame-array readout semantics.
- The maintained counting-detector families are `APDDetector` and
  `SPADArrayDetector`, with optional capability queries surfaced through
  `supports_counting_noise`,
  `supports_dead_time`, `supports_channel_gain_map`,
  `supports_counting_gating`, `supports_afterpulsing`,
  `supports_channel_crosstalk`, and `supports_paralyzable_dead_time`.
- Counting-detector dead-time behavior is selected by dispatch through
  `CountingDeadTimeModel`, with `NoDeadTime` as the null model and
  `NonParalyzableDeadTime` and `ParalyzableDeadTime` as maintained physical
  models.
- Counting-detector gating control is now explicit through
  `AbstractCountingGateModel`,
  with `NullCountingGate` and `DutyCycleGate`.
- Counting-detector correlation is now explicit through
  `AbstractCountingCorrelationModel`, with `NullCountingCorrelation`,
  `AfterpulsingModel`, `ChannelCrosstalkModel`, and
  `CompositeCountingCorrelation`.

### `IF-DM`: deformable mirrors

- `AbstractDeformableMirror` implementations must provide
  `build_influence_functions!(dm, tel)` and `apply!(dm, tel, mode)`.
- The maintained DM runtime path assumes `dm.state.coefs` is the active command
  vector and that `apply!` writes phase on the telescope OPD grid.
- The DM interface now separates:
  - topology through `AbstractDMTopology`
  - static sampled influence bases through `AbstractDMInfluenceModel`
  - actuator behavior through `AbstractDMActuatorModel`
- `ActuatorGridTopology` is the maintained regular grid topology used by the
  concise `n_act=...` constructor path.
- `SampledActuatorTopology` is the maintained path for measured or externally
  converted actuator coordinate sets.
- `MeasuredInfluenceFunctions` is the maintained first-class path for sampled
  manufacturer or lab-measured influence bases. `DenseInfluenceMatrix` remains
  the low-level sampled-basis wrapper.
- `LinearStaticActuators` is the maintained null/default actuator behavior.
  `ClippedActuators`, `ActuatorHealthMap`, and
  `CompositeDMActuatorModel` provide reusable constraint layers without
  modifying the static influence basis.

### `IF-REC`: control reconstructors

- `AbstractReconstructorOperator` is the maintained core slopes-to-command
  operator family, and the required runtime contract is
  `reconstruct!(out, recon, slopes)`.
- `reconstruct(recon, slopes)` is the thin allocation wrapper over that
  mutating operator.
- The maintained control reconstructor families are `ModalReconstructor` and
  `MappedReconstructor`, and both must present the same external
  slopes-to-command contract.
- `inverse_policy(recon)`, `singular_values(recon)`,
  `condition_number(recon)`, and `effective_rank(recon)` are the preferred
  long-form accessors for maintained reconstructor diagnostics.

### `IF-CTRL`: controllers

- `AbstractController` implementations should expose a canonical mutating update
  path through `update!(ctrl, input, dt)`.
- `controller_output(ctrl)` is the maintained accessor for the current
  command-like controller state.
- `reset_controller!(ctrl)` is the preferred optional lifecycle hook for
  controllers that support explicit state reset, and
  `supports_controller_reset(ctrl)` advertises that capability.
- The maintained controller family is `DiscreteIntegratorController`, which
  returns `controller_output(ctrl)` from `update!`.

### `IF-SIM`: control simulations

- `AbstractControlSimulation` implementations must provide `step!(sim)`,
  `AdaptiveOpticsSim.simulation_interface(sim)`, and `readout(sim)`.
- `sense!(sim)` is the plant/sensor-side execution hook when commands are
  supplied externally.
- `prepare!(sim)` is the preferred optional hook for one-time runtime/WFS
  precomputation before repeated `sense!` or `step!` calls.
- `step!(sim)` is the full closed-loop update and must be equivalent to
  `sense!` followed by reconstruction and command application.
- Output-side accessors are `command`, `slopes`,
  `wfs_frame`, `science_frame`,
  `wfs_metadata`, and `science_metadata`.
- Grouped composite simulations may additionally expose compatible-shape grouped
  exports through `grouped_wfs_stack` and
  `grouped_science_stack`.
- Optional behavior is expressed with
  `supports_prepared_runtime`, `supports_detector_output`,
  `supports_stacked_sources`, and `supports_grouped_execution`.

### `IF-CAL`: calibration workflow surfaces

- `interaction_matrix(...)` is the maintained linearized WFS-calibration entry
  point and must return an `InteractionMatrix` with a stored calibration
  `amplitude`. `forward_operator(imat)` and `calibration_amplitude(imat)` are
  the preferred long-form accessors for the stored matrix and drive amplitude.
- `CalibrationVault(D; ...)` is the maintained inverse-storage workflow and
  should retain the forward operator `D`, the selected inverse representation
  when requested, and inversion diagnostics such as singular values, condition
  number, effective rank, and truncation count. `forward_operator(vault)`,
  `inverse_operator_matrix(vault)`, `inverse_policy(vault)`,
  `singular_values(vault)`, `condition_number(vault)`,
  `effective_rank(vault)`, and `truncation_count(vault)` are the preferred
  long-form accessors for the stored operators and diagnostics.
- `modal_basis(...)` is the maintained modal-basis workflow and must return a
  `ModalBasis` with consistent `M2C`, sampled basis vectors, and optional
  projector. `modal_to_command(basis)`, `sampled_basis(basis)`, and
  `modal_projector(basis)` provide the long-form result accessors.
- `ao_calibration(...)` is the maintained packaged AO-calibration workflow and
  should return an `AOCalibration` bundling basis operators with the measured
  `CalibrationVault`. `modal_to_command(calib)`, `sampled_basis(calib)`,
  `modal_projector(calib)`, and `calibration_vault(calib)` provide the
  preferred long-form accessors for the result fields.
- `AdaptiveOpticsSim.compute_meta_sensitivity_matrix(...)` and
  `AdaptiveOpticsSim.SPRINT(...)` are the maintained
  misregistration-identification workflows and should retain the reference
  calibration, inverted meta-sensitivity operator, finite-difference step
  sizes, and active field ordering.
- `LiFT` reconstruction uses a separate iterative contract:
  `reconstruct!(coeffs, lift, psf_in, mode_ids, ...)` must mutate the supplied
  coefficient buffer in place and retain convergence diagnostics in
  `diagnostics(lift)`.

## Notes on backend-aware APIs

- Constructors typically accept `T=` and `backend=` keywords.
- The preferred public backend selectors are `CPUBackend()`, `CUDABackend()`,
  and `AMDGPUBackend()`.
- Use backend selector objects rather than raw array container types in new code.
- Methods accept `AbstractArray` inputs where practical, but structs keep
  concrete array fields so the compiler can specialize.
- FFT-heavy paths use `AbstractFFTs.jl`; future GPU-specialized kernels should
  be introduced behind trait-selected backends rather than by branching on
  symbols in hot loops.

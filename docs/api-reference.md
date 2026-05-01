# API Reference

Status: active

Related guides:

- [`user-guide.md`](user-guide.md)
- [`model-cookbook.md`](model-cookbook.md)
- [`extension-guide.md`](extension-guide.md)
- [`runtime-dataflow.md`](runtime-dataflow.md)

This document describes the maintained public API surface exposed by
`using AdaptiveOpticsSim`. Advanced methods and support utilities may still be
used as `AdaptiveOpticsSim.name`, but they are not part of the ordinary
top-level namespace unless listed here.

## Public API Policy

The exported API is curated in [`../src/exports.jl`](../src/exports.jl).
Export a name only when it is one of these:

- a normal user-facing constructor or workflow function
- a maintained physical model family
- a common mutating operation such as `measure!`, `propagate!`, or `step!`
- a stable extension seam documented in [`extension-guide.md`](extension-guide.md)

Keep these names qualified as `AdaptiveOpticsSim.name`:

- low-level state and workspace containers
- backend launch/allocation helpers
- capability traits and `supports_*` helpers
- telemetry/config internals
- one-off calibration-identification support routines
- benchmark and validation scaffolding

The package intentionally distinguishes three tiers:

- **Stable exported API:** ordinary modeling, calibration, runtime, and HIL use.
- **Advanced documented API:** maintained but usually qualified.
- **Developer support API:** available for package internals, tests, and
  validation tooling, but not promised as the normal user surface.

## Core

- Errors: `AdaptiveOpticsSimError`, `InvalidConfiguration`,
  `DimensionMismatchError`, `UnsupportedAlgorithm`, `NumericalConditionError`
- Profiles and RNG: `FidelityProfile`, `ScientificProfile`, `FastProfile`,
  `default_fidelity_profile`, `runtime_rng`, `deterministic_reference_rng`
- Backend selectors: `CPUBackend`, `CUDABackend`, `AMDGPUBackend`,
  `MetalBackend`, `AbstractArrayBackend`, `backend`
- Inverse policies: `InversePolicy`, `ExactPseudoInverse`, `TSVDInverse`,
  `TikhonovInverse`, `default_modal_inverse_policy`

## Masks And Apertures

- `CircularAperture`
- `AnnularAperture`
- `SpiderMask`
- `RectangularROI`
- `SubapertureGridMask`
- `build_mask!`
- `apply_mask!`

## Optical Models

- Telescope/source: `Telescope`, `Source`, `LGSSource`, `Asterism`
- Source accessors: `wavelength`, `optical_path`
- Telescope mutation: `reset_opd!`, `apply_opd!`, `set_pupil!`,
  `set_pupil_reflectivity!`
- Pupil helpers: `pupil_mask`, `apply_spiders!`
- PSF: `compute_psf!`, `psf_pixel_scale_arcsec`
- Spectral sources: `SpectralSample`, `SpectralBundle`, `SpectralSource`,
  `with_spectrum`
- Extended sources: `GaussianDiskSourceModel`, `PointCloudSourceModel`,
  `SampledImageSourceModel`, `with_extended_source`,
  `extended_source_asterism`
- Fields/propagation: `ElectricField`, `FraunhoferPropagation`,
  `FresnelPropagation`, `GeometricAtmosphericPropagation`,
  `LayeredFresnelAtmosphericPropagation`, `AtmosphericFieldPropagation`
- Zernike/OPD/NCPA: `ZernikeBasis`, `compute_zernike!`, `OPDMap`,
  `Misregistration`, `apply_misregistration`, `NCPA`, `KLBasis`,
  `ZernikeModalBasis`
- Spatial filtering: `SpatialFilter`, `CircularFilter`, `SquareFilter`,
  `FoucaultFilter`, `filter!`

## Atmosphere

- `AbstractAtmosphere`
- `KolmogorovAtmosphere`
- `MultiLayerAtmosphere`
- `InfinitePhaseScreen`
- `InfiniteMultiLayerAtmosphere`
- `advance!`
- `propagate!`

Atmosphere implementations are expected to mutate preallocated state in
`advance!` and apply the resulting OPD in `propagate!`.

## Deformable Mirrors And Controllable Optics

- Topology: `ActuatorGridTopology`, `SampledActuatorTopology`
- Influence models: `GaussianInfluenceWidth`, `GaussianMechanicalCoupling`,
  `DenseInfluenceMatrix`, `MeasuredInfluenceFunctions`
- Actuator behavior: `ClippedActuators`, `ActuatorHealthMap`,
  `CompositeDMActuatorModel`
- Main DM type: `DeformableMirror`
- DM accessors: `influence_model`, `influence_width`,
  `mechanical_coupling`, `n_actuators`
- Mutating application: `apply!`, `update_command!`
- Modal optics: `FunctionModalBasis`, `MatrixModalBasis`,
  `ZernikeOpticBasis`, `CartesianTiltBasis`, `ModalControllableOptic`,
  `TipTiltMirror`, `FocusStage`, `CompositeControllableOptic`
- Application modes: `DMAdditive`, `DMReplace`

The normal DM constructor supports concise Gaussian keywords. Use explicit
topology, influence, and actuator-model objects when modeling measured
influence functions, non-square actuator layouts, clipping, or actuator health.

## Detectors

- Detector types: `Detector`, `APDDetector`, `SPADArrayDetector`
- Noise: `NoiseModel`, `NoiseNone`, `NoisePhoton`, `NoiseReadout`,
  `NoisePhotonReadout`
- Sensor families: `SensorType`, `CCDSensor`, `CMOSSensor`, `EMCCDSensor`,
  `InGaAsSensor`, `HgCdTeAvalancheArraySensor`, `APDSensor`,
  `SPADArraySensor`
- Frame response: `FrameResponseModel`, `NullFrameResponse`,
  `GaussianPixelResponse`, `SampledFrameResponse`,
  `RectangularPixelAperture`, `SeparablePixelMTF`
- Defects: `PixelResponseNonuniformity`, `DarkSignalNonuniformity`,
  `BadPixelMask`, `CompositeDetectorDefectModel`
- Readout timing and correction: `RollingShutter`,
  `CorrelatedDoubleSampling`, `FowlerSampling`,
  `FrameReadoutCorrectionModel`, `NullFrameReadoutCorrection`,
  `ReferencePixelCommonModeCorrection`, `ReferenceRowCommonModeCorrection`,
  `ReferenceColumnCommonModeCorrection`,
  `ReferenceOutputCommonModeCorrection`, `CompositeFrameReadoutCorrection`
- Readout products: `FrameReadoutProducts`, `NoFrameReadoutProducts`,
  `MultiReadFrameReadoutProducts`, `HgCdTeReadoutProducts`
- Nonlinearity and persistence: `SaturatingFrameNonlinearity`,
  `ExponentialPersistence`
- Thermal models: `AbstractDetectorThermalModel`,
  `NullDetectorThermalModel`, `FixedTemperature`, `FirstOrderThermalModel`,
  `ArrheniusRateLaw`, `LinearTemperatureLaw`, `ExponentialTemperatureLaw`
- Counting models: `CountingDeadTimeModel`, `NoDeadTime`,
  `NonParalyzableDeadTime`, `ParalyzableDeadTime`, `DutyCycleGate`,
  `AfterpulsingModel`, `ChannelCrosstalkModel`,
  `CompositeCountingCorrelation`
- Runtime functions: `capture!`, `output_frame`, `detector_export_metadata`,
  `readout_ready`, `reset_integration!`, `thermal_model`

Use `bits` for detector quantization depth and `output_type` for the Julia
element type exported to an RTC/HIL boundary.

## Wavefront Sensors

- Sensing modes: `Diffractive`, `Geometric`
- WFS families: `ShackHartmannWFS`, `PyramidWFS`, `BioEdgeWFS`,
  `ZernikeWFS`, `CurvatureWFS`
- Curvature readout: `CurvatureReadoutModel`, `CurvatureCountingReadout`,
  `CurvatureBranchResponse`
- Shack-Hartmann calibration and extraction: `FluxThresholdValidSubapertures`,
  `CenterOfGravityExtraction`, `SubapertureLayout`,
  `SubapertureCalibration`, `subaperture_layout`,
  `subaperture_calibration`, `slope_extraction_model`,
  `valid_subaperture_indices`, `n_valid_subapertures`
- Flux normalization: `MeanValidFluxNormalization`,
  `IncidenceFluxNormalization`
- Measurement and WFS images: `measure!`, `pyramid_modulation_frame!`,
  `valid_subaperture_mask`, `camera_frame`, `wfs_detector_image`,
  `shack_hartmann_detector_image`, `shack_hartmann_detector_image!`
- LiFT: `LiFT`, `LiFTSolveAuto`, `LiFTSolveQR`,
  `LiFTSolveNormalEquations`, `LiFTLevenbergMarquardt`,
  `LiFTAdaptiveLevenbergMarquardt`

The maintained HIL image boundary is `wfs_detector_image(...)`. For
Shack-Hartmann sensors this returns a detector-like lenslet mosaic assembled
from the spot cube; frame-style WFS families return their maintained camera or
detector frame.

## Calibration And Reconstruction

- Interaction/control matrices: `InteractionMatrix`, `interaction_matrix`,
  `ControlMatrix`
- Modal bases: `ModalBasis`, `KLDMModes`, `KLHHtPSD`,
  `kl_modal_basis`, `modal_basis`, `basis_from_m2c`
- AO calibration: `AOCalibration`, `ao_calibration`, `control_matrix`
- Error and optical-gain calibration: `fitting_error`, `GainSensingCamera`,
  `calibrate!`, `compute_optical_gains!`
- Reconstructors: `NullReconstructor`, `ModalReconstructor`,
  `MappedReconstructor`, `reconstruct!`, `reconstruct`
- Controller: `DiscreteIntegratorController`

## Runtime And HIL

- Simulation/runtime types: `AbstractControlSimulation`,
  `AbstractExecutionPolicy`, `SimulationReadout`, `AOSimulation`,
  `SequentialExecution`, `ThreadedExecution`, `BackendStreamExecution`
- Runtime profiles and outputs: `ScientificRuntimeProfile`,
  `HILRuntimeProfile`, `RuntimeOutputRequirements`,
  `GroupedRuntimeOutputRequirements`
- Delay lines: `VectorDelayLine`, `shift_delay!`
- Runtime setup and commands: `prepare!`, `prepare_runtime_wfs!`,
  `set_command!`, `command_segments`, `command_segment_range`
- Orchestration: `ControlLoopBranch`, `SingleControlLoopConfig`,
  `GroupedControlLoopConfig`, `ControlLoopScenario`,
  `build_control_loop_scenario`, `control_loop_name`,
  `control_loop_branch_labels`
- Runtime execution: `sense!`, `step!`
- Readout accessors: `readout`, `command`, `slopes`, `wfs_frame`,
  `science_frame`, `grouped_wfs_stack`, `runtime_timing`

`ControlLoopScenario` is the preferred public assembly surface for normal
closed-loop and HIL simulations. Lower-level runtime construction remains
available as qualified advanced API for tests and specialized tooling.

## Tomography

- Parameter containers: `TomographyAtmosphereParams`, `LGSAsterismParams`,
  `LGSWFSParams`, `TomographyParams`, `TomographyDMParams`
- Reconstructors: `ModelBasedTomography`, `InteractionMatrixTomography`,
  `build_reconstructor`, `assemble_reconstructor_and_fitting`
- Signal layouts: `SimulationSlopes`, `InterleavedSlopes`, `InvertedSlopes`
- Helpers: `zenith_angle_deg`, `wind_direction_deg`,
  `reconstruct_wavefront_map`, `dm_commands`

## Extension Contracts

Use [`extension-guide.md`](extension-guide.md) for detailed instructions on
adding new detectors, WFS families, DMs, controllable optics, controllers, and
reconstructors. The short version is:

- define small concrete parameter/state types
- implement family-owned methods near the family source file
- use multiple dispatch or traits instead of central type switches
- preallocate workspaces and expose hot-path mutation with `!` methods
- keep optional plotting, file formats, and heavyweight dependencies outside
  the core package unless they are required by the maintained runtime surface

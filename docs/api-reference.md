# API Reference

Status: active

Related guides:

- [`user-guide.md`](user-guide.md)
- [`model-cookbook.md`](model-cookbook.md)
- [`glossary.md`](glossary.md)
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
  `DimensionMismatchError`, `UnsupportedAlgorithm`, `NumericalConditionError`,
  `AtmosphereTimeError`, `AtmosphereEpochError`
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

- Telescope/source: `Telescope`, `Source`, `LGSSource`, `Asterism`;
  source radiometry is declared with `PhysicalPhotonIrradianceSource` or
  `NormalizedTestSource`
- Source accessors: `wavelength`, `photon_irradiance`, `source_radiometry`,
  `source_radiometric_value`, `optical_path`
- Telescope mutation: `reset_opd!`, `apply_opd!`, `set_pupil!`,
  `set_pupil_reflectivity!`
- Pupil helpers: `pupil_mask`, `apply_spiders!`

`pupil_mask(telescope)` and `pupil_reflectivity(telescope)` are zero-copy
read views of backend-resident aperture storage. Treat them as read-only.
Change the aperture through `set_pupil!`, `set_pupil_reflectivity!`, or
`apply_spiders!`; these APIs advance the telescope aperture revision so every
WFS reference calibration is invalidated coherently. Direct array or field
mutation is unsupported because it bypasses that revision boundary.

- Prepared direct imaging: `prepare_direct_imaging`, `form_direct_image!`,
  `direct_imaging_output`, `direct_imaging_components`, and
  `focal_plane_pixel_scale_arcsec`. `DirectImagingPlan` and
  `DirectImagingWorkspace` bind exact caller-owned pupil/field/output storage
  and single-writer scratch. The output is a source-scaled, cell-integrated
  photon-arrival-rate `IntensityMap` on focal-plane angular coordinates before
  exposure; it is not implicitly a normalized PSF. Off-axis placement uses a
  preparation-time-validated, rounded-sample periodic shift on the prepared
  focal grid and respects its declared `:x`/`:y` axis order and signs
- Spectral sources: `SpectralSample`, `SpectralBundle`, `SpectralSource`,
  `with_spectrum`; construct a `SpectralSource` through `with_spectrum` from a
  `Source` or `LGSSource` leaf rather than nesting source expansions. Bundle
  weights are normalized photon-number fractions of the source's total photon
  irradiance, not radiant-energy fractions
- Extended sources: `GaussianDiskSourceModel`, `PointCloudSourceModel`,
  `SampledImageSourceModel`, `with_extended_source`,
  `extended_source_asterism`. For direct imaging, explicitly expand an
  `ExtendedSource` with `extended_source_asterism` and prepare that `Asterism`
- Optical products: `PupilFunction`, `ElectricField`, `IntensityMap`,
  `OpticalProductBundle`, `OpticalPlaneMetadata`; coordinates are declared
  with `MetricCoordinates` or `AngularCoordinates`
- Spectral coordinates: `AchromaticSpectralCoordinate`,
  `MonochromaticChannel`, or `IntegratedSpectralChannel`;
  `UnspecifiedSpectralCoordinate` is rejected by prepared intensity
  accumulation. Detector acquisition is channel-specific and requires a
  monochromatic or integrated channel; sampled QE requires the monochromatic
  case
- Prepared diffractive Shack–Hartmann formation writes same-grid sources to
  one rate mosaic and retains distinct wavelength grids as native-sampling
  leaves in an `OpticalProductBundle`. Passing a distinct-grid source to one
  output is rejected; the legacy single-product `measure!` path remains
  restricted to a common wavelength grid
- Product semantics: `PhotonRateNormalization` or
  `DimensionlessNormalization`; `PointSampledMeasure`,
  `SpatialDensityMeasure`, or `CellIntegratedMeasure`; and
  `CoherentFieldCombination`, `IncoherentIntensityAddition`, or
  `NonCombinableProduct`
- Direct imaging of a same-wavelength `Asterism` produces one compatible
  incoherent output and exposes its leaf products with
  `direct_imaging_components`. Direct imaging of a `SpectralSource` retains
  wavelength-dependent grids in an `OpticalProductBundle`; there is no
  implicit resampling or spectral sum
- Compatible intensity accumulation: `PreparedIncoherentSum`,
  `prepare_incoherent_sum`, `accumulate_intensity!`
- Fields/propagation: `FraunhoferPropagation`,
  `FresnelPropagation`, `GeometricAtmosphericPropagation`,
  `LayeredFresnelAtmosphericPropagation`, `AtmosphericFieldPropagation`
- Zernike/OPD/NCPA: `ZernikeBasis`, `compute_zernike!`, `OPDMap`,
  `Misregistration`, `apply_misregistration`, `NCPA`, `KLBasis`,
  `ZernikeModalBasis`
- Spatial filtering: `SpatialFilter`, `CircularFilter`, `SquareFilter`,
  `FoucaultFilter`, `filter!`

Known photometric bands use physical photon-irradiance radiometry by default.
A custom-band source is a normalized test source unless it supplies
`photon_irradiance` or an explicit radiometry policy. `LGSSource` likewise
requires explicit photon irradiance to claim a physical rate; its default is a
normalized source. Calling `photon_irradiance` on a normalized source is an
error rather than an implicit unit conversion.

## Atmosphere

- `AbstractAtmosphere`
- `KolmogorovAtmosphere`
- `MultiLayerAtmosphere`
- `InfinitePhaseScreen`
- `InfiniteMultiLayerAtmosphere`
- Epochs: `AtmosphereEpoch`, `current_epoch`, `epoch_time`, `epoch_sequence`
- Explicit evolution: `advance_by!`, `advance_to!`
- Direction preparation: `prepare_atmosphere_renderer`,
  `prepare_atmosphere_renderers`, `direction_renderers`
- Caller-owned rendering: `render_atmosphere!`
- Field execution: `propagate_atmosphere_field!`, `atmospheric_intensity!`
- Transitional/static extension verbs: `advance!`, `propagate!`

Timed atmosphere implementations mutate physical layer state only during an
explicit advance and publish a stable epoch. A prepared single-direction
renderer consumes the current epoch, writes a compatible caller-owned
`PupilFunction`, and neither advances the atmosphere nor consumes RNG. The
plural preparation API expands an `Asterism` or `ExtendedSource`; the singular
API rejects multi-direction sources.

`propagate!(atmosphere, telescope, source)` remains a transitional convenience
that prepares a renderer on each call. Repeated and HIL-sensitive paths should
prepare once and call `render_atmosphere!`. `advance!` remains the extension
verb for source-independent untimed/static atmosphere models; maintained timed
models use `advance_by!` or `advance_to!`.

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
Analytic Gaussian models use lazy operator storage with separable or fused
matrix-free application; dense sampled matrices are materialized only by
setup/calibration consumers that require them.
Actuator print-through is not a separate analytic DM model today; represent it
through `DenseInfluenceMatrix` or `MeasuredInfluenceFunctions` when the sampled
influence basis already includes the print-through structure.

## Detectors

- Detector types: `Detector`, `LinearAPDDetector`, `APDDetector`,
  `SPADArrayDetector`, `MKIDArrayDetector`
- Noise: `NoiseModel`, `NoiseNone`, `NoisePhoton`, `NoiseReadout`,
  `NoisePhotonReadout`
- Sensor families: `SensorType`, `CCDSensor`, `CMOSSensor`, `EMCCDSensor`,
  `InGaAsSensor`, `HgCdTeAvalancheArraySensor`, `APDSensor`, `SPADArraySensor`,
  `MKIDArraySensor`
- EMCCD modes and helpers: `LinearEMMode`, `PhotonCountingEMMode`,
  `EMOutput`, `ConventionalOutput`, `emccd_snr`
- Linear APD topology: `SingleElementAPD`, `APDChannelBank`; analog APD output
  uses `channel_output`, while `APDDetector` remains the counting-channel path.
- CMOS readout structure: `CMOSReadNoiseMap`; row and column noise, output
  groups, shutter timing, and detector defect maps compose through
  `CMOSSensor`. Vendor camera profiles are intentionally outside core.
- Frame response: `FrameResponseModel`, `NullFrameResponse`,
  `GaussianPixelResponse`, `SampledFrameResponse`,
  and `RectangularPixelAperture`. These are spatial-domain presampling response
  models. Evaluate the normalized interior, infinite-grid transfer magnitude
  of the realized discrete acquisition kernel with
  `AdaptiveOpticsSim.detector_mtf(model, fx, fy)`, where frequency is in cycles
  per detector pixel. Finite frames use zero extension and therefore have
  boundary-dependent response. This diagnostic is not a continuous
  subpixel-aperture MTF; such a model requires an explicitly prepared
  oversampled optical grid.
- Post-collection coupling: `AdaptiveOpticsSim.NullChargeCoupling` and
  `AdaptiveOpticsSim.InterpixelCapacitance`. Configure these with
  `Detector(...; charge_coupling_model=...)`; this stage runs after photon and
  generated-charge statistics rather than as a pre-shot image blur.
- Defects: `PixelResponseNonuniformity`, `DarkSignalNonuniformity`,
  `BadPixelMask`, `CompositeDetectorDefectModel`
- Readout timing and correction: `GlobalShutter`, `RollingShutter`,
  `RollingExposure`, `GlobalResetExposure`, `SequentialAcquisition`,
  `FrameTransferAcquisition`, `SingleRead`,
  `AveragedNonDestructiveReads`,
  `FunctionFrameSource`, `InPlaceFrameSource`,
  `FunctionExposureFrameSource`, `InPlaceExposureFrameSource`,
  `CorrelatedDoubleSampling`, `FowlerSampling`, `UpTheRampSampling`,
  `SkipperSampling`,
  `FrameReadoutCorrectionModel`, `NullFrameReadoutCorrection`,
  `ReferencePixelCommonModeCorrection`, `ReferenceRowCommonModeCorrection`,
  `ReferenceColumnCommonModeCorrection`,
  `ReferenceOutputCommonModeCorrection`, `CompositeFrameReadoutCorrection`
- Readout products: `FrameReadoutProducts`, `NoFrameReadoutProducts`,
  `MultiReadFrameReadoutProducts`, `UpTheRampReadoutProducts`,
  `SkipperReadoutProducts`, `HgCdTeReadoutProducts`
- Nonlinearity and persistence: `SaturatingFrameNonlinearity`,
  `ExponentialPersistence`
- Thermal models: `AbstractDetectorThermalModel`,
  `NullDetectorThermalModel`, `FixedTemperature`, `FirstOrderThermalModel`,
  `ArrheniusRateLaw`, `LinearTemperatureLaw`, `ExponentialTemperatureLaw`
- Counting models: `CountingDeadTimeModel`, `NoDeadTime`,
  `NonParalyzableDeadTime`, `ParalyzableDeadTime`, `DutyCycleGate`,
  `AfterpulsingModel`, `ChannelCrosstalkModel`,
  `CompositeCountingCorrelation`
- Runtime functions: `capture!`, `output_frame`, `channel_output`,
  `detector_export_metadata`, `readout_ready`, `reset_integration!`,
  `thermal_model`, `detector_ramp_slope`, `detector_ramp_intercept`,
  `detector_ramp_cube`, `detector_ramp_times`
- Prepared intensity-map acquisition: `DetectorAcquisitionPlan`,
  `prepare_detector_acquisition`

Use `bits` for detector quantization depth and `output_type` for the Julia
element type exported to an RTC/HIL boundary. A detector with `bits` must also
provide a fixed positive `full_well`; per-frame peak normalization is not an
ADC model and is rejected.

`Detector(...; qe=...)` accepts either a scalar quantum efficiency or a
qualified QE model such as
`AdaptiveOpticsSim.SampledQuantumEfficiency(wavelengths, values)`. Matrix-only
capture uses the scalar `params.qe` value, which is the supplied scalar or the
peak sampled QE. Source-aware capture, `capture!(det, image, src; rng=...)`,
evaluates the QE model at `wavelength(src)`. For `SpectralSource`, it uses the
flux-weighted effective QE over the spectral bundle. Pyramid frame-detector
paths specialize this boundary by applying sampled QE per wavelength before
incoherent optical-rate accumulation. Prepared diffractive Shack–Hartmann
formation instead retains distinct wavelength-rate products in a bundle so
acquisition can apply channel-specific QE without an implicit grid conversion;
its legacy single-product path accumulates only contributions on one common
grid.
Other source-aware detector paths retain the effective-QE contract unless
explicitly documented otherwise.

Response kernels, sampled QE vectors, and detector defect maps are copied at
their public construction boundaries, and `Detector` takes another run-owned
copy. Treat the resulting parameter arrays as immutable and rebuild the model
or detector to change them; direct mutation of detector parameter fields is
unsupported. WFS calibration keys use the identity of this frozen storage, so
warmed CPU and GPU checks are constant-time and do not copy device arrays to
the host. Detector-aware Pyramid, BioEdge, and Zernike reference frames apply
the same deterministic presampling, sampling, QE/exposure, binning, PRNU, and
bad-pixel-throughput path as ordinary signal acquisition, followed by a
configured built-in homogeneous reference-pixel correction. Noiseless
single-read, averaged nondestructive, correlated-double, Fowler, and valid
up-the-ramp HgCdTe readout all reduce to that transform. Calibration applies
HgCdTe avalanche gain, detector gain, and homogeneous correction in acquisition
order. Up-the-ramp schedules are validated, but calibration does not create or
mutate acquisition readout products.
They fail closed when deterministic stages outside that reference path are
configured, including saturation/quantization, nonlinearity, IPC, DSNU,
persistence, background-map subtraction, or output grouping. Stochastic
detector noise remains an acquisition effect, not part of the reference frame;
custom correction models are rejected unless their calibration behavior is
explicitly implemented.

For a metadata-validated repeated path, call
`prepare_detector_acquisition(detector, intensity_map)` once and pass the
returned plan to `capture!`. Photon-rate maps cannot be rescaled;
dimensionless maps require an explicit `normalized_to_photon_rate` conversion.
Spatial-density maps use their declared cell measure, while cell-integrated
maps are already rates per represented cell. The prepared frame path applies a
non-null presampling response before physical-pixel integration, then applies
QE and the explicit whole or incremental exposure once. A sampled QE model
requires a declared monochromatic channel on this path. Until an explicit
optical-grid mapping is prepared, a non-null response requires
`psf_sampling == 1`. Preparation rejects empty, negative, NaN, or infinite
intensity values. Repeated capture trusts later writes to the prepared storage,
so its producer is responsible for preserving finite nonnegative samples.

`MKIDArrayDetector` is the maintained MKID surface for accumulated counting-array
HIL use. It models photon-counting output with quantum efficiency, fill factor,
dark count rate, optional counting dead time/correlation models, and exported
energy-resolution and timing-jitter metadata. `energy_resolution` is the
dimensionless resolving power `E/ΔE`, and `timing_jitter_s` is in seconds.
Configure its optional inclusive passband in meters with
`wavelength_range_m=(minimum, maximum)`. Source-aware capture applies that
passband, including weighted `SpectralSource` bundles; matrix-only capture
assumes spectrally prefiltered input. The current model does not emit per-photon
timestamp/energy event lists.

`CMOSSensor` covers CMOS, sCMOS, and quantitative low-noise CMOS architectures
through composition rather than camera classes. `NoiseReadout` supplies a
uniform independent component; `row_readout_sigma`, `column_readout_sigma`,
and `CMOSReadNoiseMap` add structured components at the readout stage.
`PixelResponseNonuniformity`, `DarkSignalNonuniformity`, `BadPixelMask`, and
`StaticCMOSOutputPattern` carry measured calibration structure. Core provides
no vendor defaults or named cameras.

`SkipperSampling(n)` configures repeated nondestructive CCD sampling. The
implementation accumulates a mean online and exposes `SkipperReadoutProducts`
without retaining an `n`-plane read cube. This bounds memory independently of
sample count and keeps the warmed repeated-capture path allocation-free. The
current model assumes independent read samples. Configure CCD clock-induced
charge with `clock_induced_charge_per_frame`; unlike dark current, it is not
scaled by integration time.

`UpTheRampSampling(n)` is available on `HgCdTeAvalancheArraySensor`, including
the conventional gain-one configuration. It retains `n` evenly spaced
nondestructive reads, fits an intercept and slope, and returns
`slope * integration_time` so ordinary `capture!` output units do not change.
Use `detector_ramp_slope(det)`, `detector_ramp_intercept(det)`,
`detector_ramp_cube(det)`, and `detector_ramp_times(det)` for diagnostics.
Reads start at exposure time zero and end at the configured integration time;
`read_time` must fit within the resulting cadence. Full-frame and windowed
repeated capture reuse
their products after warmup.

The current ramp model assumes linear accumulation and independent per-read
Gaussian read noise. It shares the exposure's photon/dark realization across
the ramp and does not yet provide cosmic-ray segmentation, saturation-aware
fitting, or correlated 1/f-noise estimation.

`EMCCDSensor(...; acquisition_mode=FrameTransferAcquisition(...))` models frame
transfer as timing only. With `readout_rate_hz` configured, metadata reports
the pixel-read duration, one-frame output latency in `sampling_wallclock_time`,
and the overlapped `steady_state_frame_period`. `SequentialAcquisition()` is
the default. Neither acquisition policy changes the presampling response or its
derived MTF, QE, charge multiplication, or detector noise.

## Wavefront Sensors

- Sensing modes: `Diffractive`, `Geometric`
- Prepared products: `WFSObservationMetadata`, `WFSMeasurementMetadata`,
  `WFSObservation`, `WFSMeasurement`
- Product accessors: `observation_storage`, `observation_units`,
  `observation_metadata`, `measurement_storage`, `measurement_units`,
  `measurement_metadata`
- Prepared stage protocols: `prepare_wfs_optical_formation` /
  `form_wfs_optical_products!`, `prepare_wfs_acquisition` /
  `acquire_wfs_observation!`, and `prepare_wfs_estimation` /
  `estimate_wfs_measurement!`
- Estimation paths: `AbstractWFSMeasurementPath`,
  `AcquiredObservationPath`, `DirectMeasurementPath`,
  `wfs_measurement_path`
- Contract failure: `WFSPreparationError`, whose `stage` and open
  extension-defined `reason` fields identify rejected preparation contracts or
  execution-time prepared-binding violations before mutation
- WFS families: `ShackHartmannWFS`, `PyramidWFS`, `BioEdgeWFS`,
  `ZernikeWFS`, `CurvatureWFS`
- Shack-Hartmann optical composition: `MicrolensArrayParams`,
  `MicrolensArray`, `prepare_microlens_propagation`, `microlens_array`,
  `ShackHartmannOpticalFrontEnd`, and `shack_hartmann_rate_map`. The
  microlens array is the immutable regular-array model and its numerical
  sampling policy; prepared propagation holds only the backend/grid-bound FFT
  plans and reusable optical scratch. A front end can be assembled directly
  from that model, propagation state, and a
  `SubapertureLayout` without constructing or retaining a `ShackHartmannWFS`.
  The concrete `PreparedMicrolensPropagation` implementation type is
  intentionally qualified rather than exported; callers obtain it through the
  preparation function.
- Curvature readout: `CurvatureReadoutModel`, `CurvatureCountingReadout`,
  `CurvatureBranchResponse`
- Shack-Hartmann calibration and extraction: `FluxThresholdValidSubapertures`,
  `AbstractSlopeExtractionModel`, `CenterOfGravityExtraction`,
  `SubapertureLayout`,
  `SubapertureCalibration`, `subaperture_layout`,
  `subaperture_calibration`, `slope_extraction_model`,
  `set_subaperture_calibration!`, `valid_subaperture_indices`,
  `n_valid_subapertures`
- WFS normalization policies (with transitional type names):
  `MeanValidFluxNormalization`, `IncidenceFluxNormalization`
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

The prepared stage API is a static composition protocol, not an optical graph.
An optical front end accepts a caller-owned `PupilFunction` or pupil-plane
`ElectricField` and writes one detector-plane photon-arrival-rate
`IntensityMap`, a statically structured concrete tuple of them, or an
`OpticalProductBundle`.
Acquisition consumes those unchanged rates, owns its explicit duration and
detector state, and writes one `WFSObservation` or a concrete tuple of them.
Estimation writes a `WFSMeasurement`. Observation and measurement storage may
be a host `Ref` or an array of any rank; multiple products use concrete tuples
rather than an abstract container. Units are required and may be symbols or
application-defined singleton values.

Every concrete preparer validates all applicable shape, plane/radiometry,
backend/device, duration, mapping, and estimator contracts before repeated
execution. Mutating execution receives explicit caller-owned products and
destinations and an explicit RNG at acquisition. A direct geometric or
reduced-order estimator declares `DirectMeasurementPath()` and allocates no
fictitious rate plane, observation, or detector workspace. Shack-Hartmann,
Pyramid, and BioEdge implement the generic contract; the remaining WFS
families stay internally coupled until their ordered migration PRs.

The diffractive Shack-Hartmann staged path uses a real-valued
`:lenslet_mosaic` observation whose element type exactly matches the prepared
detector output. The family-neutral detector-acquisition stage validates the
detector-plane product and output storage; the Shack-Hartmann estimator then
validates the `:lenslet_mosaic` layout. Its declared units describe that output
and are not restricted to electrons: `:electron_count`, `:adu`, or an
application-defined singleton are valid examples. The centroid estimator requires a
`:centroid_slopes` measurement whose units match the installed calibration's
`output_units`. Geometric Shack-Hartmann estimation requires a
`:radian`/`:geometric_slopes` measurement. Physical estimator preparation also
requires an explicit finite calibration installed with
`set_subaperture_calibration!`; the prepared estimator binds its revision,
centroid response, output units, wavelength, signature, and reference storage
and rejects later recalibration before mutating its output.

`PyramidOpticalFrontEnd` and `BioEdgeOpticalFrontEnd` bind their distinct
physical masks, a prepared `NoModulation`, `CircularModulation`, or
`SampledModulation` policy, reusable propagation storage, and a source. Use
`pyramid_rate_map` or `bioedge_rate_map` to allocate the caller-owned
four-pupil rate product, then prepare and execute it with the generic optical
stage functions. Spectral sources return an `OpticalProductBundle`. An
`Asterism` or `ExtendedSource` requires one explicit pupil input per rendered
direction and also returns a bundle; no direction-dependent pupils are added
by array index. Install acquired-estimator references with
`set_pyramid_calibration!` or `set_bioedge_calibration!`. The setters validate
finite storage and advance a revision, so a prepared estimator must be rebuilt
after recalibration. Geometric variants instead consume an explicit
`PupilFunction` and report a `:metre`/`:geometric_slopes` direct measurement;
they have no focal-plane or detector workspace.

The legacy constructor keyword `calib_modulation` prepares the broader optical
quadrature used only to select valid estimator support. A zero-aberration
reference is formed with the operating modulation so it subtracts the static
response of the sensor that will actually acquire data. A user-sampled
modulation path is retained for both operations rather than replaced by a
generated circle.

Both geometric and centroid measurements use `[axis 1; axis 2]` block order.
Each `n_lenslets`-by-`n_lenslets` block follows Julia column-major order, so
`reshape(view(measurement, block), n_lenslets, n_lenslets)` has the same
`(i, j)` indexing as the subaperture mask. Invalid subapertures are reported as
zero. The calibration reference is an `n_lenslets^2`-by-2 table whose columns
are the two centroid components in that same lenslet order. The maintained
OOPAO reference harness performs the explicit axis-block and row-major adapter;
the package API does not expose OOPAO's storage convention.

Microlens sampling, the valid-subaperture layout, and centroid calibration are
cold configuration. Updating a layout through the maintained mutation API
advances its revision; prepared optical and estimator plans then reject reuse
until they are prepared again. A layout update made through a complete
`ShackHartmannWFS` also invalidates its calibration. Component-level users must
install a calibration matching the new layout before preparing another
estimator. Direct mutation of layout masks, host mirrors, calibration fields,
or calibration reference storage is unsupported because it bypasses revision
tracking. Caller-owned pupil values, rate destinations, observations, and
measurements remain mutable between repeated executions of a compatible plan.

`CenterOfGravityExtraction` currently implements thresholded, unwindowed
centroiding. Supplying a window is rejected rather than silently ignored;
windowed, correlation, and matched-filter estimators remain future policies.

## Calibration And Reconstruction

- Interaction/control matrices: `InteractionMatrix`, `interaction_matrix`,
  `ControlMatrix`. Caller-owned calibration storage is available through the
  qualified `AdaptiveOpticsSim.interaction_matrix!` API.
- Modal bases: `ModalBasis`, `KLDMModes`, `KLHHtPSD`,
  `kl_modal_basis`, `modal_basis`, `basis_from_m2c`
- AO calibration: `AOCalibration`, `ao_calibration`, `control_matrix`
- Error and optical-gain calibration: `fitting_error`, `GainSensingCamera`,
  `calibrate!`, `compute_optical_gains!`
- Reconstructors: `NullReconstructor`, `ModalReconstructor`,
  `FactorizedReconstructor`, `MappedReconstructor`,
  `ControlledReconstructor`, `reconstruct!`, `reconstruct`
- Controller: `DiscreteIntegratorController`. `ControlledReconstructor`
  composes a reconstructor and stateful controller without adding a runtime
  branch.

## Runtime And HIL

- Simulation/runtime types: `AbstractControlSimulation`,
  `AbstractExecutionPolicy`, `SimulationReadout`, `AOSimulation`,
  `SequentialExecution`, `ThreadedExecution`, `BackendStreamExecution`,
  `DeterministicExecution`, `AcceleratedKernelsExecution`, `DaggerExecution`
- Optical-path sources: `wfs_source`, `science_source`. `AOSimulation` uses its
  WFS source for atmosphere/WFS sensing and may carry a distinct science source
  for science-camera propagation.
- Runtime profiles and outputs: `ScientificRuntimeProfile`,
  `HILRuntimeProfile`, `RuntimeOutputRequirements`,
  `GroupedRuntimeOutputRequirements`
- Runtime residency: `CPUHILExecutionPlan`,
  `DeviceResidentExecutionPlan`, `runtime_execution_plan`,
  `synchronize_runtime!`
- Shared optical arms: `OpticalWFSChannel`, `SharedOpticalArm`,
  `SharedOpticalRuntime`, `primary_runtime`, `optical_arms`,
  `science_frames`, `wfs_signals`
- Delay lines: `VectorDelayLine`, `shift_delay!`
- Runtime setup and commands: `prepare!`, `prepare_runtime_wfs!`,
  `set_command!`, `command_segments`, `command_segment_range`
- Orchestration: `ControlLoopBranch`, `SingleControlLoopConfig`,
  `GroupedControlLoopConfig`, `ControlLoopScenario`,
  `build_control_loop_scenario`, `control_loop_name`,
  `control_loop_branch_labels`
- Coarse ensembles: `SimulationEnsemble`; `prepare!`, `sense!`, and `step!`
  apply to every member. Qualified advanced access is available through
  `AdaptiveOpticsSim.run_ensemble!`, `ensemble_members`, `execution_policy`,
  and `ensemble_readouts`.
- Runtime execution: `sense!`, `step!`
- Readout accessors: `readout`, `command`, `slopes`, `wfs_frame`,
  `science_frame`, `grouped_wfs_stack`, `runtime_timing`,
  `runtime_atmosphere_step`

`ControlLoopScenario` is the preferred public assembly surface for normal
single-plant closed-loop and HIL simulations. `SharedOpticalRuntime` is the
typed surface when auxiliary source/WFS/science arms share that plant.
Single and grouped control-loop configurations require an explicit positive
`atmosphere_step`; it advances model time for one sensing update and is
independent of detector exposure.
`SimulationEnsemble` owns independent plants for offline sweeps and ensemble
simulation. Direct `step!` calls remain the CPU HIL path; Dagger and
AcceleratedKernels are optional coarse schedulers, not inner-loop runtime
plans.
Lower-level runtime construction remains
available as qualified advanced API for tests and specialized tooling.
For external Proper science-arm integration, use
[`proper-integration-guide.md`](./proper-integration-guide.md).

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

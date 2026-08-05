# API Reference

Status: active

Related guides:

- [`user-guide.md`](user-guide.md)
- [`model-cookbook.md`](model-cookbook.md)
- [`glossary.md`](glossary.md)
- [`extension-guide.md`](extension-guide.md)
- [`runtime-dataflow.md`](runtime-dataflow.md)

This document describes the currently implemented ordinary package API
imported by `using AdaptiveOpticsSim`, optical-foundation vocabulary imported
by `using AdaptiveOpticsSim.Optics`, backend vocabulary imported by
`using AdaptiveOpticsSim.Backends`, conventional detector vocabulary imported
by `using AdaptiveOpticsSim.Detectors`, atmosphere vocabulary imported by
`using AdaptiveOpticsSim.Atmospheres`, calibration vocabulary imported by
`using AdaptiveOpticsSim.Calibration`, control vocabulary imported by
`using AdaptiveOpticsSim.Control`, tomography vocabulary imported by
`using AdaptiveOpticsSim.Tomography`, coarse offline execution vocabulary
imported by `using AdaptiveOpticsSim.Ensembles`, the routine plant workflow
imported by `using AdaptiveOpticsSim.Plant`, and stable qualified APIs addressed through
their canonical modules. Qualified public names are maintained but do not
enter the caller's ordinary namespace. The root package exports the canonical
modules, not compatibility exports for their moved contents.

The breaking namespace migration is complete. `Backends`, `Optics`,
`Atmospheres`, `Detectors`, `WavefrontSensors`, `Plant`, `Calibration`,
`Control`, `Tomography`, and `Ensembles` are canonical owners. `Atmospheres`
owns atmosphere
models and state, source-direction rendering and batching, and
atmosphere-coupled propagation. `Detectors` owns both conventional frame/area
and counting/channel sensor APIs. `Optics` owns
apertures, telescopes, sources, optical products and
metadata, pupil/field formation, backend-portable Fraunhofer and Fresnel
propagation, direct imaging, sampled OPD, physical NCPA, controllable optics,
deformable mirrors, spatial filtering, and reusable physical WFS components.
`WavefrontSensors` owns composed sensing, estimation, and its geometric slope
kernels. `Calibration` owns inverse policies, interaction and control
matrices, modal bases, model-derived NCPA synthesis, fitting, and calibration
workflows. The root provides no compatibility bindings for domain-owned
contents. The exact final root and domain allowlists are frozen in
[`../test/contracts/namespace_authority.toml`](../test/contracts/namespace_authority.toml).
The maintained implementation stage is recorded in
[`../test/contracts/namespace_migration_state.toml`](../test/contracts/namespace_migration_state.toml).

## Public API Policy

Root exported and qualified-public names are curated in
[`../src/exports.jl`](../src/exports.jl); the canonical plant surface is
curated separately in [`../src/plant/api.jl`](../src/plant/api.jl), the
conventional detector surface in
[`../src/detectors/api.jl`](../src/detectors/api.jl), and the
canonical backend surface in
[`../src/backends/api.jl`](../src/backends/api.jl). The calibration surface is
curated in [`../src/calibration/api.jl`](../src/calibration/api.jl). The currently implemented
optical surface is curated in
[`../src/optics/api.jl`](../src/optics/api.jl), and the atmosphere surface in
[`../src/atmosphere/api.jl`](../src/atmosphere/api.jl). Each owner PR must converge on
the exact domain allowlists without root forwarding aliases or compatibility
adapters.
Export a name only when it is one of these:

- a normal user-facing constructor or workflow function
- a maintained physical model family
- a common mutating operation such as `measure!`, `propagate!`, or `step!`
- a stable extension seam documented in [`extension-guide.md`](extension-guide.md)

Declare a stable name `public` but keep it qualified under its owning module
when it is one of these:

- low-level state and workspace containers
- backend launch/allocation helpers
- capability traits and `supports_*` helpers
- telemetry/config internals
- one-off calibration-identification support routines
- benchmark and validation scaffolding

The package intentionally distinguishes three tiers:

- **Stable exported API:** ordinary modeling, calibration, runtime, and plant
  construction/execution use.
- **Advanced documented API:** declared public, maintained, and usually
  qualified.
- **Developer support API:** available for package internals, tests, and
  validation tooling, but not promised as the normal user surface.

The final root allowlist is deliberately small:

- canonical modules: `Backends`, `Optics`, `Atmospheres`, `Detectors`,
  `WavefrontSensors`, `Calibration`, `Control`, `Tomography`, `Ensembles`, and
  `Plant`
- shared errors: `AdaptiveOpticsSimError`, `InvalidConfiguration`,
  `DimensionMismatchError`, `UnsupportedAlgorithm`, and
  `NumericalConditionError`
- fidelity and RNG services: `FidelityProfile`, `ScientificProfile`,
  `FastProfile`, `default_fidelity_profile`, `runtime_rng`, and
  `deterministic_reference_rng`
- shared timing vocabulary: `runtime_timing`

All other supported bindings are owned by a domain module or removed. A root
binding is never retained merely for migration compatibility.

## Core

- Errors: `AdaptiveOpticsSimError`, `InvalidConfiguration`,
  `DimensionMismatchError`, `UnsupportedAlgorithm`, `NumericalConditionError`,
  `AtmosphereTimeError`, `AtmosphereEpochError`, and `WFSPreparationError`
- Profiles and RNG: `FidelityProfile`, `ScientificProfile`, `FastProfile`,
  `default_fidelity_profile`, `runtime_rng`, `deterministic_reference_rng`
- Inverse policies: `InversePolicy`, `ExactPseudoInverse`, `TSVDInverse`,
  `TikhonovInverse`, `default_modal_inverse_policy`

## Backends

Import routine backend vocabulary explicitly:

```julia
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Backends

backend = CPUBackend()
```

`Backends` exports `CPUBackend`, `CUDABackend`, `AMDGPUBackend`,
`MetalBackend`, `AbstractArrayBackend`, `backend`, and `compute_device`.
`AbstractComputeDevice`, `HostComputeDevice`, `AcceleratorComputeDevice`,
`compute_device_backend`, `compute_device_identifier`, the exact-device
availability result types/accessors, `ComputeDeviceError`, and
`allocate_device_array(device, T, dims...)` are stable qualified-public names,
such as `Backends.HostComputeDevice()`. Exact-device availability is a cold
runtime fact, not a model-support, capacity, or latency claim. Backend-family
allocation, launch, FFT, reduction, and registration helpers remain developer
or extension seams unless the exact owner allowlist promotes them.

## `AdaptiveOpticsSim.Optics` Foundations

Import routine optical-foundation vocabulary explicitly:

```julia
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Optics

tel = Telescope(resolution=32, diameter=8.0)
src = Source(band=:I, magnitude=8.0)
```

`AbstractSource` is a stable qualified-public extension seam and is addressed
as `Optics.AbstractSource`; it is not exported into an ordinary caller
namespace.

### Masks And Apertures

- `CircularAperture`
- `AnnularAperture`
- `SpiderMask`
- `RectangularROI`
- `SubapertureGridMask`
- `build_mask!`
- `apply_mask!`

### Optical Models

- Telescope/source: `Telescope`, `Source`, `LGSSource`, `Asterism`;
  source radiometry is declared with `PhysicalPhotonIrradianceSource` or
  `NormalizedTestSource`
- Qualified cold telescope preparation: `Optics.AbstractTelescopeDefinition`,
  `Optics.TelescopeDefinition`, `Optics.prepare_telescope(definition, target)`,
  and `Optics.validate_telescope_target(telescope, target)`. A definition
  contains aperture configuration and no prepared arrays, backend, compute
  device, or mutable state; preparation materializes one independent
  numerical `Optics.AbstractTelescope` on the exact target. Custom prepared
  telescope types validate their owned data-plane arrays through the
  fail-closed validation seam
- Source accessors: `wavelength`, `photon_irradiance`, `source_radiometry`,
  `source_radiometric_value`, `optical_path`
- Telescope aperture mutation: `set_pupil!`, `set_pupil_reflectivity!`,
  `apply_spiders!`
- Path-owned pupil products: `PupilFunction`, `pupil_support`,
  `pupil_amplitude`, `opd_map`, `reset_opd!`, `apply_opd!`
- Aperture helpers: `pupil_mask`, `pupil_reflectivity`

`pupil_mask(telescope)` and `pupil_reflectivity(telescope)` are zero-copy
read views of backend-resident aperture storage. Treat them as read-only.
Change the aperture through `set_pupil!`, `set_pupil_reflectivity!`, or
`apply_spiders!`; these APIs advance the telescope aperture revision so every
WFS reference calibration is invalidated coherently. Direct array or field
mutation is unsupported because it bypasses that revision boundary.

`PupilFunction(telescope)` snapshots the current aperture support and converts
intensity reflectivity to field amplitude with its square root. Each optical
path owns its mutable `PupilFunction`; the telescope contains no mutable OPD.
Prepared propagation and sensing APIs accept that explicit product and reject
incompatible geometry revisions, backends, or devices.

- Prepared direct imaging: `prepare_direct_imaging`,
  `prepare_direct_imaging_batch`, `form_direct_image!`,
  `direct_imaging_output`, `direct_imaging_components`, and
  `focal_plane_pixel_scale_arcsec`. The qualified-public
  `Optics.DirectImagingPlan` contains the reusable sampling, field-formation,
  propagation, metadata, and integer-shift contract without retaining exact
  caller arrays. `Optics.DirectImagingWorkspace` contains replaceable
  single-writer propagation and shift scratch. `Optics.PreparedDirectImaging`
  binds the exact caller-owned input, work field, output product, plan, and
  workspace; execute only `form_direct_image!(prepared)`. The output is a
  source-scaled, cell-integrated
  photon-arrival-rate `IntensityMap` on focal-plane angular coordinates before
  exposure; it is not implicitly a normalized PSF. Off-axis placement uses a
  preparation-time-validated, rounded-sample periodic shift on the prepared
  focal grid and respects its declared `:x`/`:y` axis order and signs
- Prepared native Fraunhofer batches accept a physical `Source`,
  same-wavelength ordered `Asterism`, or `SpectralSource{Source}` with one
  shared `PupilFunction`, or one exact path-local pupil per ordered physical
  source leaf. Compatible leaves share one concrete source-storage contract,
  one stacked field allocation, one stacked photon-rate allocation, and one
  FFT plan over the optical axes. Their public result remains an
  `OpticalProductBundle` of ordinary wavelength- and direction-specific
  `IntensityMap` views. Fixed-size host memberships keep runtime cardinality
  out of the type and retain independent exact-binding snapshots so element
  replacement fails before numerical mutation. Batching does not coherently
  combine sources, integrate wavelengths, apply detector exposure, copy to
  host, or regroup paths during warmed execution. The compatibility trait,
  signature, prepared type, and validation accessors are qualified public API
  rather than root exports
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
  with `MetricCoordinates` or `AngularCoordinates`.
  `Backends.compute_device(array)` reports the concrete host or accelerator
  identity independently of its semantic array-backend family
- Spectral coordinates: `AchromaticSpectralCoordinate`,
  `MonochromaticChannel`, or `IntegratedSpectralChannel`;
  `UnspecifiedSpectralCoordinate` is rejected by prepared intensity
  accumulation. Detector acquisition is channel-specific and requires a
  monochromatic or integrated channel; sampled QE requires the monochromatic
  case
- Prepared diffractive Shack–Hartmann optics write same-grid sources to
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
  implicit resampling or spectral sum. Use `prepare_direct_imaging_batch` when
  each leaf must remain an independently consumable product while compatible
  native Fraunhofer work is submitted as one prepared stack
- Compatible intensity accumulation: `PreparedIncoherentSum`,
  `prepare_incoherent_sum`, `accumulate_intensity!`
- Fields and general propagation: `FraunhoferPropagation` and
  `FresnelPropagation`. Each is a prepared owner containing a qualified-public
  `Optics.AbstractPropagationPlan` implementation and a concrete workspace.
  Fraunhofer plans contain sampling and plane metadata; Fresnel plans also own
  their logically immutable backend-resident transfer coefficients. FFT
  handles and transform scratch remain in replaceable single-writer
  workspaces. Use `Optics.propagation_plan`,
  `Optics.propagation_workspace`, and the metadata accessors for advanced
  integration. Atmosphere-coupled propagation is owned by
  `AdaptiveOpticsSim.Atmospheres`
- Sampled OPD and physical NCPA: `OPDMap` and `NCPA`. `NCPA(opd)` stores only
  the explicit backend-resident optical-path-difference map in metres.
  KL/Zernike/M2C basis selection, coefficient generation, and atmosphere- or
  DM-derived synthesis are calibration policy; synthesis returns an
  `Optics.NCPA` but the physical optic does not retain calibration provenance
- Optical bases and registration: `ZernikeBasis`, `compute_zernike!`,
  `Misregistration`, and `apply_misregistration`. Import calibration-owned
  `KLBasis` and `ZernikeModalBasis` from `AdaptiveOpticsSim.Calibration`
- Spatial filtering: `SpatialFilter`, `CircularFilter`, `SquareFilter`,
  `FoucaultFilter`, `prepare_spatial_filter`, and `filter!`.
  `prepare_spatial_filter` returns an exact `Optics.PreparedSpatialFilter`;
  its reusable plan borrows no telescope aperture array, its output
  `PupilFunction` owns support, and its FFT resources remain in a replaceable
  single-writer workspace. Because Base also
  exports an unrelated collection operation named `filter!`, call
  `Optics.filter!` or explicitly
  `import AdaptiveOpticsSim.Optics: filter!`
- Controllable optics: `DeformableMirror`, `ModalControllableOptic`,
  `TipTiltMirror`, `FocusStage`, `set_command!`, and `update_surface!`
- Reusable WFS optics: focal-plane modulation models, `MicrolensArray`,
  `PyramidPhaseMask`, `BiOEdgeAmplitudeMask`, `ZernikePhaseSpot`, and
  `CurvatureDefocusPair`. Composed WFS front ends, detector acquisition, and
  estimators remain outside `Optics`

Known photometric bands use physical photon-irradiance radiometry by default.
A custom-band source is a normalized test source unless it supplies
`photon_irradiance` or an explicit radiometry policy. `LGSSource` likewise
requires explicit photon irradiance to claim a physical rate; its default is a
normalized source. Calling `photon_irradiance` on a normalized source is an
error rather than an implicit unit conversion.

## `AdaptiveOpticsSim.Plant`

`AdaptiveOpticsSim.Plant` is the canonical owner of virtual plant time,
scheduling and trigger distribution, detector-acquisition events, plant
definitions and preparation, illumination and providers, command lifecycle,
and serial event composition. Use:

```julia
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Plant
```

Routine construction and execution vocabulary is then available unqualified.
Policies, mutable state/workspace owners, low-level transitions, traits, and
inspection accessors are stable qualified API such as
`Plant.CommandValuePolicy` and `Plant.EventSchedulerState`. Plant-specific
errors are owned and exported here: `PlantTimeError`, `PlantScheduleError`,
`PlantDefinitionError`, `PlantCommandError`, `PlantPreparationError`, and
`DetectorAcquisitionError`. `StructuralResourceError` remains qualified with
the resource-reporting API below.

The completed Gate 3 API defines canonical plant instants, elapsed durations, nominal periodic
recurrence, and a fixed-capacity serial event calendar without adding wall-clock
pacing. The completed Gate 2 topology API declares one shared telescope and
atmosphere, reusable optical paths, and independent acquisitions. Preparation
adds concrete single-writer owners without implicit atmosphere advancement:

- Canonical plant-time values: `PlantTimestamp`, `PlantDuration`,
  `PlantTimeOffset`, `plant_nanoseconds`, `plant_time_seconds`, and
  `plant_duration_seconds`
- Nominal recurrence: `PeriodicSchedule`, `schedule_period`, `schedule_phase`,
  and `schedule_timestamp`

- Qualified scheduler definitions and ownership: `EventGeneratorDefinition`,
  `PreparedEventScheduler`, `EventSchedulerState`,
  `EventSchedulerWorkspace`, `EventGeneratorHandle`, and `EventClaim`
- Qualified scheduler preparation and inspection: `prepare_event_scheduler`,
  `event_generator_handle`, `event_generator_count`,
  `event_scheduler_capacity`, `scheduler_timestamp`,
  `scan_due_events!`, `due_event_count`, `due_event_timestamp`, and
  `due_event_key`
- Qualified scheduler transitions: `claim_next_event!`,
  `claimed_event_key`, `reschedule_event!`,
  `reschedule_periodic_event!`, `activate_event_generator!`, and
  `deactivate_event_generator!`
- Trigger declarations: `TriggerSourceID`, `TriggerLinkID`,
  `TriggerConsumerID`, `TriggerFaultID`, `TriggerSourceDefinition`,
  `TriggerLinkDefinition`, `TriggerConsumerDefinition`,
  `TriggerFaultTraceEntry`, and `TriggerFaultTrace`; `TriggerEdgeAction` and
  its values remain qualified
- Trigger preparation and qualified state: `prepare_trigger_topology`,
  `PreparedTriggerTopology`, `TriggerTopologyState`,
  `TriggerTopologyWorkspace`, `trigger_source_handle`,
  `trigger_source_count`, `trigger_link_count`, `trigger_consumer_count`,
  `trigger_in_flight_capacity`, and
  `required_trigger_in_flight_capacity`
- Qualified trigger execution: `next_trigger_source`,
  `realize_next_trigger_source!`, `next_trigger_delivery`,
  `next_trigger_delivery_timestamp`, `pop_next_trigger_delivery!`,
  `pending_trigger_delivery_count`, `contains_trigger_fault`, and the
  `trigger_fault_observation_count`, `trigger_fault_observation`,
  `trigger_fault_id`, and `trigger_fault_location` accessors. Delivery records
  retain distinct
  `NominalTriggerEdge`, `DeliveredTriggerEdge`, and
  `ReportedTriggerTimestamp` values
- Global-shutter definition and qualified event lifecycle:
  `GlobalShutterAcquisitionDefinition`,
  `PreparedGlobalShutterAcquisition`, `GlobalShutterAcquisitionState`,
  `DetectorAcquisitionStatus`, `detector_acquisition_status`,
  `detector_acquisition_sequence`, `exposure_start_timestamp`,
  `exposure_close_timestamp`, `integrated_through_timestamp`,
  `readout_complete_timestamp`, `acquisition_readiness_timestamp`,
  `nondestructive_read_count`, `nondestructive_read_offset`,
  `next_nondestructive_read_timestamp`,
  `prepare_global_shutter_acquisition`, `begin_exposure!`,
  `accumulate_exposure_interval!`, `take_nondestructive_read!`,
  `close_exposure!`, `complete_readout!`, and
  `mark_acquisition_ready!`. Integer plant timestamps own all transitions;
  the separately owned state rejects busy retriggers and intervals that cross
  exposure close or a pending scheduled ramp read
- Rolling-shutter definition and qualified event lifecycle:
  `RollingShutterAcquisitionDefinition`,
  `PreparedRollingShutterAcquisition`, `RollingShutterAcquisitionState`,
  `rolling_band_count`, `rolling_band_rows`,
  `rolling_band_open_timestamp`, `rolling_band_close_timestamp`,
  `rolling_opened_band_count`, `rolling_closed_band_count`,
  `next_rolling_band_open_timestamp`,
  `next_rolling_band_close_timestamp`, and
  `prepare_rolling_shutter_acquisition`. The detector's `RollingShutter`
  timing model supplies line time, row-group size, and rolling-exposure or
  global-reset semantics; the ordinary event transition functions above are
  specialized for the prepared rolling lifecycle
- Frame-transfer definition and qualified event lifecycle:
  `FrameTransferAcquisitionDefinition`,
  `PreparedFrameTransferAcquisition`, `FrameTransferAcquisitionState`,
  `frame_transfer_storage_capacity`, `frame_transfer_image_sequence`,
  `frame_transfer_storage_sequence`, `frame_transfer_product_sequence`,
  `frame_transfer_complete_timestamp`, `frame_transfer_image_ready`,
  `frame_transfer_storage_empty`, `frame_transfer_readout_pending`,
  `prepare_frame_transfer_acquisition`, and `complete_frame_transfer!`.
  One prepared device-resident storage frame permits image-area integration
  to overlap storage-area readout when the declared timing fits that capacity
- Direct-measurement definition and qualified lifecycle:
  `DirectMeasurementAcquisitionDefinition`,
  `PreparedDirectMeasurementAcquisition`,
  `DirectMeasurementAcquisitionState`,
  `prepare_direct_measurement_acquisition`, and
  `accumulate_direct_measurement_interval!`. This intentional
  non-detector lifecycle time-averages the held instantaneous
  `WFSMeasurement` over a half-open exposure, then follows explicit readout
  completion and readiness delays. It creates no observation or implicit
  detector physics
- Plant-event composition: `PeriodicAcquisitionStart`,
  `TriggeredAcquisitionStart`, `OpticalSampleDefinition`,
  `AcquisitionEventDefinition`, `PlantEventLoopDefinition`,
  qualified `PreparedPlantEventLoop`, `PlantEventLoopState`,
  `PlantEventLoopWorkspace`, `prepare_plant_event_loop`,
  `plant_event_path_count`, `plant_event_acquisition_count`,
  `plant_event_generator_count`, `next_plant_event_timestamp`,
  `step_plant_events!`, `run_plant_events_until!`,
  `acquisition_product_sequence`, and
  `acquisition_product_ready_timestamp`. This HIL-neutral serial oracle
  composes exact periodic or delivered-trigger acquisition starts with
  independently periodic optical paths and complete acquisition products; it
  retains the prepared plant's exact prepared device execution context. Its
  public step and finite-horizon calls enter that context and restore the
  caller's previous device and stream selections only after establishing the
  retained stream's backend completion boundary. It owns no wall clock, task,
  queue, port, transport, or RTC protocol
- Mixed-resource plant-event composition: the
  `prepare_plant_event_loop(::PreparedPlantPartitions, ...)` method returns a
  qualified `PreparedMixedResourcePlantEventLoop`; callers construct
  `MixedResourcePlantEventLoopState` and
  `MixedResourcePlantEventLoopWorkspace`, then use the same
  `admit_plant_command!`, `admit_plant_command_transaction!`,
  `command_disposition_count`, `command_disposition`,
  `clear_command_dispositions!`, `effective_command`,
  `acquisition_products`, `acquisition_product_sequence`,
  `acquisition_product_ready_timestamp`, `next_plant_event_timestamp`,
  `step_plant_events!`, and `run_plant_events_until!` operations as the
  single-resource loop. This deterministic serial oracle retains one
  atmosphere and command authority while entering every complete
  path/acquisition group's exact prepared target. The initial method supports
  full-optical frame-detector lifecycles and rejects autonomous periodic
  optics during preparation. After preparation and exact-path warmup, its
  successful CPU command-admission and complete event-run paths allocate zero
  Julia heap bytes; preparation, compilation, diagnostics, exceptions, and
  accelerator or KernelAbstractions launch behavior are separate contracts.
  It does not add a worker, queue, pacing, placement, product-egress, or
  transport policy
- Qualified single-device path batching:
  `PreparedDevicePathBatchOwner`, `device_path_batch_owner_count`,
  `device_path_batch_owner`, `device_path_batch_compute_device`,
  `device_path_batch_backend`, `device_path_batch_group_count`,
  `device_path_batch_group_ordinal`,
  `path_execution_group_device_batch_owner_ordinal`,
  `materialize_device_path_batch!`, and `execute_device_path_batch!`.
  Event-loop preparation automatically binds compatible, co-resident,
  equal-schedule accelerator paths into exact capability groups. Native
  Fraunhofer direct science shares atmosphere and optical batches. Maintained
  Shack-Hartmann, Pyramid, or Bi-O-edge groups share the retained device context
  and atmosphere-direction batch while invoking the exact family's existing
  lenslet/modulation/spectral pipeline in canonical member order. Original
  path-local products remain device resident, and return establishes backend
  completion before any member is marked complete. CPU paths, singletons,
  unequal schedule/origin paths, mixed or unsupported WFS families, and
  incompatible signatures/products retain the independent Gate 6 lifecycle.
  Each public independent-group worker call enters the retained context in its
  own task and completes that stream before publishing the group's ready or
  complete phase
- Qualified structural-resource reporting:
  `StructuralResourceOwnerID`, `ResourceEstimateMethod`,
  `KnownStructuralResourceFact`, `UnknownStructuralResourceFact`,
  `OpaqueResourceReserve`, `StructuralResourceReport`,
  `structural_resource_fact`, `structural_array_bytes`,
  `aggregate_structural_resource_facts`, and
  `require_exact_structural_resource_facts`. The event-loop overload derives
  stable owner IDs and reports one non-overlapping numerical-storage graph for
  an exact device; accelerator-prepared loops may also report their explicit
  host mirrors and staging as a separate `HostComputeDevice()` partition.
  Logical owners with no array storage in the selected partition remain
  visible as known zero-byte facts. Whole-graph reporters currently cover the
  promoted deterministic direct-science, conventional-detector, and
  monochromatic diffractive Shack–Hartmann profiles; other path or acquisition
  implementations fail closed until they define exact ownership.
  Resident, reusable workspace, and caller-declared opaque provider reserves
  remain separate. Unknown ownership, unsupported storage, duplicate owners,
  wrong devices, inconsistent estimate versions, and byte-count overflow fail
  structurally. Reporting is cold inspection and does not walk object fields,
  query allocator state, copy accelerator arrays to host, or assign execution
  resources

- Stable identities: `AtmosphereLayerID`, `ControllableOpticID`,
  `SampledAberrationID`, `CommandEndpointID`, `PlantCommandSchemaID`,
  `OpticalPathID`, `AcquisitionID`, `RNGOwnerIdentity`
- RNG derivation and replay: `RNGDerivationVersion`, `rng_replay_metadata`
- Definitions: `ControllableOpticDefinition`,
  `SampledAberrationDefinition`, `PlantCommandSchema`,
  `OpticalPathDefinition`, `AcquisitionDefinition`, `PlantDefinition`
- Native Plant deformable mirror: `DeformableMirrorModel`, an exported
  `AdaptiveOpticsSim.Plant` cold definition that binds actuator topology,
  influence and actuator-response models, device-internal `Misregistration`,
  and `PupilRelayRegistration` without owning live command/surface storage.
  Plant preparation validates one vector actuator endpoint and constructs
  separate active/staged native `DeformableMirror` runtimes
- Optical-component placement and visibility:
  `PupilPlanePlacement`, `AtmosphericConjugatePlacement`,
  `FocalPlanePlacement`, `AllPathVisibility`, `SelectedPathVisibility`;
  controllable-optic and sampled-aberration definitions require explicit
  placement and visibility
- Plant-command value types: `PlantCommandSchemaVersion`,
  `CommandBasisRevision`, `CommandUnit`, `CommandSignConvention`,
  `CommandBasis`, `UnboundedCommandValues`, `UniformCommandBounds`,
  `PlantCommandSequence`, `CommandPresentationID`, `PlantCommand`,
  `PlantCommandOrderKey`, `PlantCommandAdmission`,
  `PlantCommandDisposition`, and `CommandDispositionReason`
- Qualified plant-command policies: `CommandValueSemantics`, `InvalidCommandAction`,
  `CommandRangeStage`, `CommandSequenceAction`, `FutureCommandPolicy`,
  `LateCommandPolicy`, `CommandSupersessionPolicy`, `CommandSilenceAction`,
  `CommandAgeOrigin`, `CommandValuePolicy`, `CommandSequencePolicy`,
  `CommandEffectiveTimePolicy`, and `CommandSilencePolicy`; their corresponding
  enum values are also stable qualified configuration vocabulary
- Cold-model trait: `plant_model_definition_style`,
  `ColdPlantModelDefinition`
- Identity and model accessors: `controllable_optic_id`,
  `sampled_aberration_id`,
  `command_schemas`, `command_endpoint_ids`, `command_schema_id`,
  `command_schema_version`, `command_endpoint_id`, `path_id`, `acquisition_id`,
  `acquisition_path_id`, `controllable_optic_model`, `path_source`,
  `path_model`, `acquisition_model`
- Command-schema accessors and validation: `command_numeric_type`,
  `command_dimensions`, `command_units`, `command_sign_convention`,
  `command_basis`, `command_basis_revision`, `command_semantics`,
  `command_bounds`, `command_value_policy`, `command_sequence_policy`,
  `command_effective_time_policy`, `command_silence_policy`,
  `validate_plant_command_payload`
- Prepared controller-output routing: exported `ControllerOutputRoute` and
  `prepare_controller_output_routing`; qualified
  `PreparedControllerOutputRoute`, `PreparedControllerOutputRouting`,
  `prepared_controller_output_routes`, `controller_output_route`,
  `controller_output_product`, `controller_output_endpoint`,
  `controller_output_schema`, and `controller_output_payload`. Preparation
  binds every named borrowed controller product to one distinct prepared
  endpoint and validates exact numeric type, shape, backend, and compute
  device. Routing adds no command timing, admission, transaction, queue, or
  transport semantics
- Standalone bounded command admission: `PreparedCommandEndpoint`,
  `prepare_command_endpoint`, `validate_plant_command`,
  `admit_plant_command!`, `claim_next_application_ready_command!`,
  `claimed_command_payload`, `mark_plant_command_applied!`,
  `fail_plant_command_application!`, and `fail_pending_plant_commands!`.
  `CommandSequenceClass`, `CommandAdmissionStatus`, `CommandTerminalKind`,
  and their values provide qualified result vocabulary; qualified `command_*`
  accessors inspect admissions, claims, and dispositions.
  `command_requested_effective_timestamp`, `command_scheduled_timestamp`,
  `command_admission_timestamp`, `command_ready_timestamp`, and
  `command_terminal_timestamp` keep the distinct command lifecycle instants
  explicit
- Qualified public standalone effective-command application:
  `effective_command`,
  `last_command_application_timestamp`, `apply_claimed_plant_command!`,
  `next_command_silence_timestamp`, and
  `apply_command_silence_transition!`. `PlantCommandSilenceTransition` and
  its `command_silence_*` accessors describe the exact replayable safe/fail
  transition; `command_endpoint_failed` and
  `last_command_admission_timestamp` expose endpoint lifecycle state
- Qualified public command-endpoint mutable storage:
  `AdaptiveOpticsSim.Plant.CommandEndpointState` and
  `AdaptiveOpticsSim.Plant.CommandApplicationState`, plus the caller-owned
  `AdaptiveOpticsSim.Plant.CommandDispositionWorkspace`. These remain explicit
  state/workspace containers rather than exported model types
- Plant accessors: `telescope_definition`, `atmosphere_definition`,
  `plant_definition`, `prepared_telescope`, `prepared_atmosphere`,
  `controllable_optic_definitions`, `sampled_aberration_definitions`,
  `path_definitions`,
  `acquisition_definitions`, `controllable_optic_definition`,
  `sampled_aberration_definition`,
  `command_endpoint_owner`, `command_schema`, `plant_command_schema`,
  `path_definition`, `acquisition_definition`
- Ordinary prepared boundary: `PreparedPlant`, `prepare_plant`,
  `prepare_acquisition_selection`, `execute_acquisition_selection!`, and
  `execute_acquisition_selection_at!`. Preparation requires one positional
  exact compute-device target; query a prepared plant through
  `Backends.compute_device(plant)`. Qualified extension/execution seams are
  `prepare_pupil_opd_materialization`, `materialize_path_input!`,
  `execute_path!`, and `execute_acquisition!`
- Qualified preparation-only partition boundary:
  `ResolvedPlantPartitionAssignment`, `prepare_plant_partitions`,
  `PreparedPlantPartitions`, `PreparedTargetPartition`,
  `PreparedAtmosphereAuthority`, `prepared_partitions`,
  `prepared_partition`, `prepared_atmosphere_authority`,
  `target_local_controllable_optic_owners`, and
  `partition_resource_report`. Every visible logical optic has at most one
  role-neutral owner per exact target, shared by all co-located paths. A complete
  `CommandEndpointConfiguration` set supplies initial effective values whenever
  the declaration contains endpoints
- Qualified typed cross-device handoff boundary:
  `AbstractHandoffPayloadContract`, `PreparedCrossDomainHandoff`,
  `HandoffPreparationIdentity`, `HandoffSlotReference`,
  `prepare_cross_domain_handoff`, `handoff_source_target`,
  `handoff_destination_target`, `handoff_capacity`, and
  `handoff_payload_bytes`. A contract extends `handoff_payload_eltype`,
  `handoff_payload_axes`, and `validate_handoff_publication`. The caller-driven
  lifecycle uses `try_next_free_handoff_slot!`, `producer_handoff_payload`,
  `submit_handoff!`, `try_complete_handoff!`,
  `try_borrow_completed_handoff!`, and `reclaim_handoff!`; inspect it through
  `handoff_slot_status` and `handoff_slot_failure_reason`. `HandoffStatus` and
  `HandoffSlotStatus` provide allocation-free structured outcomes. This
  object requires one external owner to serialize all lifecycle methods and
  is not thread-safe. An explicit backend failure is reclaimable only when the
  backend guarantees quiescence; an unexpected exception leaves the slot in
  `HandoffTransferUncertain`. The boundary owns neither a queue nor a HIL
  payload lease, wait strategy, placement policy, or transport
- Qualified authority-owned pupil-OPD publication boundary:
  `PreparedPupilOPDPublicationRoute`, prepared through
  `prepare_pupil_opd_publication_route`. Its direct and remote implementations
  are private; inspect the selected case through
  `pupil_opd_publication_authority_target`,
  `pupil_opd_publication_route_target`, and
  `pupil_opd_publication_handoff_capacity`. A path model opts in explicitly by
  extending `prepare_pupil_opd_publication_materialization` and returning one
  native, exactly bound `PreparedPupilOPDMaterialization`. Allocate
  caller-owned metadata output once with
  `prepare_pupil_opd_publication_output`, then drive
  `materialize_pupil_opd_publication!`, `submit_pupil_opd_publication!`,
  `try_complete_pupil_opd_publication!`,
  `apply_pupil_opd_publication!`, and
  `reclaim_pupil_opd_publication!`. `MaterializedPupilOPDPublication` carries
  the exact opaque route identity, authority, epoch, timestamp, path, and
  generation. The route identity binds the frozen source and prepared
  materialization without retaining either in publication metadata;
  `PupilOPDPublicationStatus` reports expected lifecycle outcomes. A direct
  route writes its exact target-local `PupilFunction`. A remote route owns one
  fixed handoff slot and transfers only OPD. One external owner must serialize
  the lifecycle and must not advance the authority until every due route has
  applied. Core provides no polling or all-routes batch coordinator here
- Qualified target-local command inspection:
  `CommandAuthorityIdentity`, `EffectiveCommandPublicationSequence`,
  `EffectiveCommandPublication`, `PreparedTargetLocalCommandEndpoint`,
  `TargetLocalCommandEndpointState`,
  `PreparedTargetLocalControllableOptic`,
  `TargetLocalControllableOpticState`,
  `TargetLocalControllableOpticWorkspace`, and their `command_authority_identity`,
  `effective_command*`, `target_local_*`, and
  `last_effective_command_publication_*` accessors. Publication metadata does
  not own an array payload. The committed target-local active copy is sealed
  from caller, staging, endpoint, and owner aliases. Public mutation remains
  withheld until the HIL composition layer can enforce the sole-publisher and
  bounded-payload-lifetime contract
- Qualified exact-target extension seams:
  `prepare_target_local_controllable_optic`,
  `validate_controllable_optic_target`,
  `validate_controllable_optic_state_target`,
  `validate_controllable_optic_workspace_target`,
  `validate_pupil_surface_coupling_target`,
  `validate_autonomous_optic_coupling_target`,
  `validate_path_execution_target`, `validate_path_materialization_target`,
  `validate_acquisition_execution_target`,
  `validate_acquisition_provider_target`, and
  `validate_illumination_evaluator_target`. Their unknown-type defaults fail
  closed. Controllable-optic preparation, event-loop state, and event-loop
  workspace have distinct seams because those owners are constructed at
  different preparation boundaries. Built-in subsystem validators remain
  implementation details. Custom prepared WFS plans use the WFS-owned qualified
  `WavefrontSensors.validate_wfs_target` seam
- Qualified sampled-aberration preparation:
  `PreparedSampledAberration`,
  `PreparedSampledAberrationPathBindings`,
  `prepared_sampled_aberrations`,
  `prepared_sampled_aberration_path_bindings`, and their `sampled_aberration_*`
  and `prepared_sampled_aberration_*` accessors. Preparation makes a
  exact-target defensive OPD copy and resolves bounded, canonical
  path-local couplings before execution
- Qualified controllable-optic path coupling:
  `AbstractPupilSurfacePathCoupling`,
  `PreparedDirectPupilSurfaceCoupling`,
  `PupilRelayRegistration`,
  `PreparedIdentityPupilFootprintCoupling`,
  `PreparedPupilFootprintCoupling`,
  `prepare_controllable_optic_path_coupling`,
  `prepare_sampled_pupil_footprint_coupling`,
  `apply_controllable_optic_surface!`, and
  `apply_sampled_pupil_surface!`. The latter accepts `DMAdditive()` or
  `DMReplace()`; transformed replacement writes zero outside finite sampled
  support. These values bind an exact prepared path and surface-grid contract;
  they own no mutable surface, command state, cadence, or transport
- Calibration-illumination entry boundary:
  `PupilFunctionIlluminationEntry`, `ElectricFieldIlluminationEntry`,
  `IntensityMapIlluminationEntry`, `ExternalOpticsResultIlluminationEntry`,
  `DetectorInputIlluminationEntry`, `UniformIntensityIllumination`,
  `prepare_illumination_entry` and `evaluate_illumination!`. The combination
  policies `Plant.SingleIllumination` and
  `Plant.ExclusiveIlluminationSelection`, like the `illumination_*`
  accessors, remain qualified. Coherent-field and incoherent-intensity
  combination reuse the root optical policies `CoherentFieldCombination` and
  `IncoherentIntensityAddition`
- Prepared accessors: `prepared_paths`, `prepared_acquisitions`,
  `prepared_path`, `prepared_acquisition`, `path_input`, `path_result`,
  `path_result_key`, `acquisition_provider`, `acquisition_products`,
  `acquisition_observation`, `acquisition_measurement`, and
  `acquisition_product_metadata`
- Prepared execution-target contract: qualified
  `Plant.PathExecutionRequirements`, `path_execution_group_path_id`,
  `path_execution_target_support`, `path_execution_target_supported`, and
  `path_execution_target_rejection_reason`. A single-device prepared group
  supports only its exact device; a same-family alternate device reports
  `:requires_repreparation` rather than moving at runtime.
- Product-provider contract: qualified traits
  `Plant.FullOpticalProviderStyle`,
  `Plant.CommandResponsiveReducedOrderProviderStyle`, and
  `Plant.SyntheticReplayProviderStyle`; qualified
  `acquisition_provider_style`,
  `acquisition_provider_payload_work`, `acquisition_product_contract`,
  `validate_acquisition_product_contract`,
  `prepare_full_optical_provider`,
  `prepare_unchanged_synthetic_provider`,
  `prepare_copied_synthetic_provider`, and
  `prepare_cyclic_replay_provider`
- Qualified model-extension boundary: `PreparedPathExecutor`,
  `PreparedAcquisitionOwner`, `PreparedAcquisitionProvider`,
  `AcquisitionProducts`, `AcquisitionProductContract`, `PathResultKey`,
  `PreparedAcquisitionSelection`, `PreparedPupilOPDMaterialization`,
  `PreparedIlluminationEntry`, `AtmosphereIndependentPath`,
  `AbstractOpticalSamplingContract`, `InstantaneousOpticalSample`,
  `require_path_result`, `prepare_path_executor`,
  `prepare_acquisition_provider`, `validate_path_execution_binding`,
  `validate_path_materialization_binding`,
  `validate_path_materialization`, `validate_atmosphere_rendering`, and
  `validate_acquisition_execution_binding`; custom providers also extend
  qualified `validate_acquisition_provider_binding`,
  `execute_acquisition_provider!`, and, for custom copied products,
  `copy_acquisition_product!`. RNG extensions use qualified
  `additional_path_rng_owner_roles`,
  `additional_path_materialization_rng_owner_roles`,
  `additional_acquisition_rng_owner_roles`, `execute_path_rngs!`,
  `execute_acquisition_rngs!`, and `rng_stream_state`

`PlantTimestamp`, `PlantDuration`, and qualified `PlantTimeOffset` are distinct,
checked integer-nanosecond values: they represent a nonnegative run-local
instant, a nonnegative elapsed interval, and a signed displacement between
nominal and realized instants. Applying an offset must still produce a
representable nonnegative timestamp. `PeriodicSchedule` describes only a
positive nominal period and a nonnegative phase. These values do not execute
events, own mutable cursor state, read wall time, or imply detector timing.

The scheduler surface remains qualified developer API after Gate 3 closure.
Preparation copies definitions into a flat,
canonical, fixed-length registry and allocates no run-length event list.
`EventSchedulerState` is the single writer for compact cursors;
`EventSchedulerWorkspace` owns fixed due slots. One `EventClaim` may be
outstanding at a time and must be resolved by rescheduling or deactivation.
Equal-timestamp rescheduling remains legal only when the incremented occurrence
makes the complete key strictly later. The scheduler has no callback, product,
detector, transport, task, execution-clock, or pacing responsibility.

The trigger surface is likewise qualified. Preparation canonicalizes explicit source, link,
consumer, trace, and fault identities into a flat finite fan-out; checks exact
finite phase-step, jitter, drop, duplicate, label-offset, non-overtaking, and
capacity behavior; and allocates fixed propagation, realization, observation,
and pending-delivery storage. `TriggerTopologyState` owns one source sequence
and phase offset per source, one phase offset per link, and no run-length edge
list. Source
faults remain correlated through every surviving branch, while a link fault
affects only that link and its descendants. Equal-time source realization
precedes delivery removal, and neither operation may silently backdate the
other. One prepared topology supports at most 64 distinct fault identities in
its compact delivery fault set. The record-returning `next_trigger_delivery`
and two-argument
`pop_next_trigger_delivery!` are inspection conveniences; latency-sensitive
code uses `next_trigger_delivery_timestamp` and the caller-owned
`Ref{TriggerDelivery}` overload of `pop_next_trigger_delivery!` to avoid boxing
identity-bearing records. This layer assigns no detector meaning to an edge and
owns no RTC command, wall clock, task, ring, transport, or GPU work.

An illumination entry is a prepared path-input materializer, not a calibration
mode or a second acquisition API. It binds one exact caller-owned optical
product, an explicit downstream-visibility description, immutable evaluator
parameters, separate single-writer state/workspace, combination semantics,
backend/device, plant time, and a path-owned `:illumination` RNG stream. The
ordinary path execution and acquisition provider consume the resulting product.
`UniformIntensityIllumination` is the only native source definition in this
gate; unusual source physics extend the qualified evaluator seams described in
the extension guide. Entry tags do not infer a lamp, relay, instrument, control
authority, or upstream propagation bypass.

Every plant stores one immutable `Optics.AbstractTelescopeDefinition` and one
immutable `Atmospheres.AbstractTimedAtmosphereDefinition`; it does not accept
already prepared or mutable telescope and atmosphere objects. Every
controllable optic, sampled aberration, command endpoint, path, and
acquisition carries an explicit typed identity. Tuples, named tuples, and
vectors organize cold declarations but do not define identity;
`PlantDefinition` copies them into fixed-size topology registries, and named
keys must match the IDs they contain. `PlantDefinition` rejects duplicate
optic/aberration/path/acquisition identities,
duplicate active schema identities, command endpoints with more than one optic
owner, and unknown path references with `PlantDefinitionError`.
Controllable-optic, optical-path, and acquisition model
types are rejected by default and must opt in to the cold-definition contract
by returning `ColdPlantModelDefinition()` from
`plant_model_definition_style(::Type{MyDefinition})`. That opt-in asserts that
instances contain configuration only. Each controllable-optic definition owns
one or more immutable `PlantCommandSchema` values, while preparation workspaces,
mutable optic/simulation/acquisition and command state, schedules, RNG streams,
queues, transport, and HIL descriptors are intentionally absent.

The first eight Gate 4 slices record controllable-optic/endpoint ownership,
versioned semantic payload contracts, bounded endpoint state, physical-optic
preparation, deterministic command/event composition, controller-output
routing, sampled device feedback, and autonomous periodic optics. A
`CommandEndpointConfiguration` supplies each declared endpoint's bounded
calendar/history capacities and copied initial and optional safe values.
It does not select a backend. `prepare_plant` derives array-command storage
from its exact target, requires exactly one configuration for every declared
endpoint, and prepares every declared optic; neither is silently omitted.
Scalar command values remain host-resident. Stable endpoint ordinals and optic
execution order derive from typed identities rather than declaration position.

`prepare_controllable_optic` returns immutable model-specific preparation
data. `prepare_controllable_optic_state` and
`prepare_controllable_optic_workspace` separately construct mutable physical
state and scratch. Array-valued initial commands passed to each state
constructor are fresh state-owned copies rather than caller or prepared-plan
storage. During command application,
`stage_controllable_optic_command!` validates and stages one complete
effective endpoint value without changing the visible surface. Both staging
and commit receive the plant-effective timestamp so a time-dependent physical
model can preserve continuity. After a successful stage,
`commit_controllable_optic_command!` must be a bounded, nonthrowing publication
operation. `prepare_controllable_optic_path_coupling` freezes the exact
path-local optical mapping. `apply_controllable_optic_surface!` receives that
mapping as its fourth argument and applies the visible physical response to an
already materialized path input. Sampled OPD models may delegate preparation
and application to `prepare_sampled_pupil_footprint_coupling` and
`apply_sampled_pupil_surface!`.

`AutonomousPeriodicOpticDefinition` binds one path-local autonomous optic to a
full-optical path, immutable fidelity, and `FreeRunningPhaseReference`,
`TriggerSourcePhaseReference`, or `TriggerResetPhaseReference`.
`CircularPyramidModulator` is the native model and declares separate bounded
radius, frequency, phase, and enabled endpoint roles. The corresponding
`CycleAveragedModulationFidelity` updates the existing prepared circular
quadrature in place from radius and enabled state; frequency and
trigger-relative phase remain analytic inspection state because they do not
change a complete-cycle average. The qualified
`autonomous_waveform_phase`, `autonomous_waveform_offset`, reference
timestamp/sequence/count, enabled, and radius accessors expose that state.
Time-resolved modulation and steering-servo dynamics are not part of this
baseline.

Successful admission copies caller payload storage. Late apply-now commands
retain their requested time but schedule at the current plant time; the
endpoint never backdates. Only successful admission enters sequence history,
so capacity/time rejection is retryable as a new presentation. Incremental
endpoints cannot use the lossy older-pending supersession policy.

Exactly one separately owned `CommandApplicationState` binds an explicit
initial held value and, when configured, a copied safe value to one exact
endpoint-state owner before its first successful admission.
`apply_claimed_plant_command!` transactionally applies absolute replacement or
incremental addition, checks the resulting finite value, enforces declared
application-stage bounds, requires the immutable scheduled timestamp to equal
the claim/application timestamp, and records one terminal disposition.
Rejected or failed application leaves the held value unchanged. Array
application uses preallocated backend-resident staging and buffer exchange.
Before the first qualifying command event, the endpoint's configured initial
timestamp is the silence-age baseline. Thereafter successful admission rebases
`AgeFromAdmission`, while successful effective application rebases
`AgeFromApplication`. A non-hold silence policy produces one exact safe/fail
transition per unchanged age origin, and a due command at the same timestamp
must be resolved first.

The routed `admit_plant_command!` overload selects the exact event-loop-owned
endpoint and arms its command-phase generator. Its admission timestamp may be
the initial timestamp or an unprocessed instant after the current scheduler
timestamp, but it must not follow the next unprocessed plant event. Once any
phase at a timestamp has run, routed admission cannot reuse that timestamp.
This frontier prevents ingress from changing endpoint state across an older
pending plant event while permitting command ingress between scheduled events.
`PlantCommandTransaction` explicitly groups two or more endpoints on distinct
physical optics at one common requested timestamp; every member endpoint must
use `PreservePendingCommands`. `admit_plant_command_transaction!` has the same
frontier precondition and preflights and stages every payload before changing
any endpoint calendar.
At application, every effective command and physical response is staged before
any member is published. A rejected or failed member leaves every transaction
member's held command and visible physical state unchanged. Equal timestamps,
shared placement, and packed payloads never imply atomicity. A transaction
terminated during admission has neither a transaction identity nor a scheduled
timestamp.

Plant state is right-continuous: a command effective at `t` affects an optical
sample at `t`. Detector exposure intervals are half-open, so accumulation
ending at `t` retains the state used strictly before `t`. The event loop copies
terminal endpoint outcomes into a fixed-capacity
`PlantEventLoopWorkspace` ledger; callers inspect them with
`command_disposition_count` and `command_disposition`, then acknowledge them
with `clear_command_dispositions!`.

Plant preparation resolves all-path or selected-path visibility into bounded
canonical `PreparedControllableOpticPathBindings` and co-placed plane groups.
Each due full-optical path applies only its direct visible binding range; a
prepared path-local autonomous role retains its exact focal-plane coupling.
Sampled-surface models can prepare identity or affine finite-support
pupil-footprint couplings for NGS and finite-height LGS paths.
Device-specific stroke, slew, settling, hysteresis, and feedback remain model
extensions or later physical-model work. HIL session, external-clock, lease,
ring, completion-credit, and transport metadata remain outside core command
values. Execution-clock ingress liveness is distinct from replayable
plant-time command silence.

Package-emitted disposition reasons are stable nonempty symbols. Admission may
emit `:endpoint_mismatch`, `:schema_mismatch`,
`:schema_version_mismatch`, payload-validation reason symbols,
`:payload_validation_failure`, `:duplicate_sequence`, `:stale_sequence`,
`:reordered_sequence`, `:skipped_sequence`, `:future_command`,
`:late_command`, `:equal_time_endpoint_conflict`, `:calendar_capacity`,
`:payload_storage_failure`, or `:superseded_by_newer_sequence`; successful
reporting emits `:applied`.
Effective-command application may additionally emit `:applied_clipped`,
`:nonfinite_rejected`, `:nonfinite_failure`, `:out_of_range_rejected`,
`:out_of_range_failure`, `:missed_application_timestamp`, or
`:application_storage_failure`. A fail-on-silence drain emits
`:command_silence`; other application-failure and pending-drain callers may
supply another nonempty `CommandDispositionReason`. Explicit transactions use
`:atomic_transaction_aborted` for otherwise acceptable members when another
member rejects during application, preserve the rejecting member's own reason,
and prefix admission-preflight abort reasons with
`:atomic_transaction_aborted_`. Unexpected physical staging failures use
`:physical_application_failure`; a model may provide a more specific
`PlantCommandError` reason.

`prepare_plant(definition, target; ...)` requires one exact compute-device
target and one explicit `run_seed`, accepts a versioned
`rng_derivation_version`, freezes each path source, and dispatches on the
concrete cold model types to build target-, shape-, and revision-bound owners.
It materializes independent numerical telescope and timed-atmosphere owners
from the reusable cold definitions and retains a target-specific execution
context. Every `PreparedPathExecutor` retains that exact prepared device
execution context, and each corresponding `PreparedAcquisitionOwner` retains
the same context. Public `materialize_path_input!`, `execute_path!`, and
`execute_acquisition!` calls enter the retained context and restore the
caller's previous device and stream selections after success or an exception.
Every prepared data-plane array in the owned graph must report the exact target;
scalar command values, host configuration and registries, and
deliberate host RNG staging are exempt. Target selection restores the caller's
prior accelerator context after success or failure. A `PathResultKey` records
source geometry, spectral sampling, radiometry, optical and propagation model
keys, instantaneous-sample semantics, output-plane contract, revisions,
backend, and device. Its
descriptive values are defensively snapshotted, and its value equality/hash
contract is intended for cold compatibility lookup rather than warmed
execution. Prepared owner constructors validate that concrete execution plans
retain their exact input, result, atmosphere materialization, provider,
detector, observation, estimator, and product storage. Each acquisition has
one run-immutable full-optical, command-responsive reduced-order, or
nonresponsive synthetic/replay provider. All providers write and return the
same caller-owned `AcquisitionProducts` logical contract for that acquisition;
required metadata captures any geometry, radiometry, units, layout, or
semantics not already present in its typed products. Several acquisitions may
borrow the exact same read-only path result while retaining distinct provider,
detector, readout, WFS estimator, observation, and measurement state.
Preparation also derives
separate stateful streams for each declared atmosphere layer, path/provider,
acquisition/detector, and extension-declared device role. Multilayer
atmospheres used by a prepared plant therefore require explicit `layer_ids`;
duplicate owner identities fail before execution.
`rng_replay_metadata(plant)` returns the run seed, derivation
version, derivation and stream algorithms, and canonical
owner-to-derived-seed records without exposing mutable RNG state. RNG
derivation does not use Julia's `hash`.

`prepare_acquisition_selection` resolves caller-supplied acquisition IDs once,
deduplicates only the full-optical paths those providers require, and stores
the selected path and acquisition registries in stable-ID order independent of
declaration and caller-selection order. Reduced-order and synthetic/replay
providers bypass otherwise unused full-optical path execution.
In scheduled execution, an explicit `OpticalSampleDefinition` is itself a
demand for the path's device-ready optical product even when that path has no
acquisition consumer. A scheduled path bypasses full-optical work only when
every acquisition consumer bound to it is reduced-order or synthetic/replay;
an empty consumer set is not an implicit reduced-order declaration.
`execute_acquisition_selection!` preflights every owner and current
`AtmosphereEpoch` before mutation, materializes all unique path inputs, applies
each prepared sampled-aberration plan, invokes each typed path executor once,
and then runs the selected acquisitions. The `_at!` variant first advances the
shared atmosphere to an explicit absolute model time.
Both methods enter the prepared plant's retained device execution context and
restore the caller's previous device and stream selections. They use exact
owner-bound RNG streams and do not accept tuple-position-dependent caller RNGs.
Low-level stochastic model
APIs continue to receive an explicit `AbstractRNG`, while prepared execution
supplies it directly rather than performing a registry lookup. Neither method
introduces cadence, triggers, a scheduler, ports, or a retained atmosphere
snapshot. A provider that requires effective-command and exposure state may
reject this schedule-free execution boundary. In particular, the built-in
linear reduced-order provider must run through `prepare_plant_event_loop` with
a `DirectMeasurementAcquisitionDefinition`; it does not invent initial
commands for `execute_acquisition_selection!`. The warmed successful serial
call is allocation-free on the
maintained Julia version; cold selection/preparation, metadata construction,
compilation, and exceptional paths are outside that budget.

In serial event execution the corresponding full-optical order is atmosphere
materialization, sampled aberrations, controllable surfaces, autonomous
path-local optics, and then the typed native or external path executor.
Replacement sampled aberrations precede additives; preparation rejects more
than one visible replacement on a path because replacement operations do not
commute.

The qualified reduced-order construction vocabulary is:

- `HarmonicDisturbanceModel`, a deterministic explicit-plant-time modal
  disturbance
- `ReducedOrderCommandResponse`, one calibrated endpoint-to-residual operator
  with exact command units, sign convention, basis, and basis revision
- `LinearReducedOrderAcquisitionModel`, which composes a path projection,
  every currently visible endpoint response, and a sensor operator into one
  direct `WFSMeasurement`
- `reduced_order_disturbance`, `reduced_order_residual`,
  `reduced_order_sample_timestamp`, and `reduced_order_residual_rms` for
  validation inspection

Preparation copies numerical operators to the path backend, validates command
schemas and memory domains, records a positive calibration revision, mandatory
operating envelope, residual metric, and omitted effects, and resolves
endpoint identities to fixed slots. At each path sample the event loop first
integrates the prior held measurement through the sample timestamp, then
evaluates the disturbance and current effective commands. A command therefore
affects only later integration intervals. Reduced-order-only paths do not
advance or execute their otherwise unused full-optical path. The maintained
implementation publishes exposure-averaged direct measurements; approximate
raw pixels are not yet a supported reduced-order surface.

The maintained [Gate 2 serial plant
artifact](../benchmarks/results/gate2/2026-07-21-serial-plant.toml) composes
science, off-axis NGS Shack-Hartmann, and finite-height LGS pyramid paths with
four detector acquisitions, including two unequal-exposure acquisitions sharing
one science result. It records deterministic declaration-order replay, zero
warmed allocation, raw HdrHistogram distributions, and exact environment
metadata. This is self-paced in-process CPU service time; it does not establish
an event-scheduler, external-RTC, transport, or fixed-arrival latency claim.

The maintained Gate 3 [scheduler](../benchmarks/results/gate3/2026-07-21-event-scheduler-gate3-closure.toml)
and [composed multi-rate plant](../benchmarks/results/gate3/2026-07-21-multi-rate-plant.toml)
artifacts add current-revision generator scaling, exact science/NGS/LGS replay,
trigger-fault fan-out, conventional detector lifecycle, fixed-storage,
direct-jump, allocation, and service-cost evidence. These qualified APIs remain
the deterministic serial virtual-time oracle; they do not expose a wall clock,
RTC transport, port, task, command endpoint, or parallel placement policy.

## Atmosphere

Import this vocabulary with `using AdaptiveOpticsSim.Atmospheres`.
`AdaptiveOpticsSim` exports the `Atmospheres` module itself and provides no
root forwarding bindings for atmosphere-owned names.

- `AbstractAtmosphere`
- `KolmogorovAtmosphere`
- `MultiLayerAtmosphere`
- `InfinitePhaseScreen`
- `InfiniteMultiLayerAtmosphere`
- Qualified cold timed-atmosphere preparation:
  `Atmospheres.AbstractTimedAtmosphereDefinition`,
  `Atmospheres.KolmogorovAtmosphereDefinition`,
  `Atmospheres.MultiLayerAtmosphereDefinition`,
  `Atmospheres.InfiniteMultiLayerAtmosphereDefinition`, and
  `Atmospheres.prepare_timed_atmosphere(definition, telescope, target)`;
  custom prepared owners implement the fail-closed
  `Atmospheres.validate_timed_atmosphere_target(atmosphere, target)` seam
- Epochs: `AtmosphereEpoch`, `current_epoch`, `epoch_time`, `epoch_sequence`;
  custom timed owners retain a qualified-public `AtmosphereIdentity` and
  `AtmosphereTimelineState`, construct the latter with
  `new_atmosphere_timeline(T)`, and implement `atmosphere_identity` and
  `atmosphere_timeline`
- Explicit evolution: `advance_by!`, `advance_to!`
- Direction preparation: `prepare_atmosphere_renderer`,
  `prepare_atmosphere_renderers`, `direction_renderers`
- Homogeneous direction batching: `prepare_atmosphere_direction_batch`,
  `render_atmosphere_directions!`, `atmosphere_direction_output`
- Caller-owned single-direction rendering: `render_atmosphere!`
- Field execution: `propagate_atmosphere_field!`, `atmospheric_intensity!`
- Static extension verbs: `advance!`, `propagate!`

Timed atmosphere implementations mutate physical layer state only during an
explicit advance and publish a stable current-state epoch token. The token does
not retain layer storage and becomes stale after the next advance. A prepared
single-direction renderer consumes the current epoch, writes a compatible
caller-owned `PupilFunction`, and neither advances the atmosphere nor consumes RNG. The
plural preparation API expands an `Asterism` or `ExtendedSource`; the singular
API rejects multi-direction sources.

`prepare_atmosphere_direction_batch` freezes the direction order and binds a
caller-owned `(pupil axis 1, pupil axis 2, direction)` atmospheric OPD stack.
Finite and infinite multilayer models support this capability. The prepared
contract fixes atmosphere and layer identity, source geometry, pupil grid,
numeric type, backend, concrete compute device, and capacity. Complete
validation precedes output mutation. CPU execution remains the serial
single-direction oracle; accelerator execution uses device-resident geometry
and establishes an explicit device-ready completion boundary.

Maintained timed atmospheres require prepared renderers and
`render_atmosphere!`. `propagate!(atmosphere, pupil[, source])` and `advance!`
remain extension verbs for source-independent untimed/static atmosphere
models; maintained timed models use `advance_by!` or `advance_to!`.

## Deformable Mirrors And Controllable Optics

- Topology: `ActuatorGridTopology`, `SampledActuatorTopology`
- Influence models: `GaussianInfluenceWidth`, `GaussianMechanicalCoupling`,
  `DenseInfluenceMatrix`, `MeasuredInfluenceFunctions`
- Actuator behavior: `ClippedActuators`, `ActuatorHealthMap`,
  `CompositeDMActuatorModel`
- Main DM type: `DeformableMirror`
- DM accessors: `influence_model`, `influence_width`,
  `mechanical_coupling`, `n_actuators`
- Command and surface application: `set_command!`, `update_surface!`,
  `apply_surface!`
- Modal optics: `FunctionModalBasis`, `MatrixModalBasis`,
  `ZernikeOpticBasis`, `CartesianTiltBasis`, `ModalControllableOptic`,
  `TipTiltMirror`, `FocusStage`
- Application modes: `DMAdditive`, `DMReplace`

`update_surface!(optic)` materializes the optic's current surface into
optic-owned storage. `apply_surface!(pupil, optic, mode)` then adds or replaces
the OPD of one explicit path. This separation lets one surface feed multiple
independently owned paths without using telescope state as shared scratch.

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

Import conventional frame and area-detector vocabulary with
`using AdaptiveOpticsSim.Detectors`. The root exports the `Detectors` module,
not the moved detector bindings.

- Conventional frame/area detector: `Detectors.Detector`
- Channel and counting-array detectors: `LinearAPDDetector`,
  `SPADArrayDetector`, and `MKIDArrayDetector`
- Noise: `NoiseModel`, `NoiseNone`, `NoisePhoton`, `NoiseReadout`,
  `NoisePhotonReadout`
- Conventional sensor families: `SensorType`, `CCDSensor`, `CMOSSensor`,
  `EMCCDSensor`, `InGaAsSensor`, and `HgCdTeSensor`
- Linear-avalanche area sensor: `HgCdTeAvalancheArraySensor`; its avalanche
  parameters and multiplication policy remain distinct from the qualified-public
  `Detectors.HgCdTeReadout` sampling configuration shared with
  `HgCdTeSensor`. Select the qualified-public
  `Detectors.ConditionalGammaAvalancheMultiplication` CPU reference model or
  `Detectors.ClippedGaussianAvalancheMultiplicationApproximation` portable
  moderate-charge approximation explicitly when the default policy is not the
  intended statistical contract.
- Counting sensor families: `SPADArraySensor` and `MKIDArraySensor`
- EMCCD modes and helpers: `LinearEMMode`, `PhotonCountingEMMode`,
  `EMOutput`, `ConventionalOutput`, `emccd_snr`
- Linear-mode APD topology: `SingleElementLinearAPD`,
  `LinearAPDChannelBank`; analog APD output uses `channel_output`.
- CMOS readout structure: `CMOSReadNoiseMap` and
  qualified-public `Detectors.StaticCMOSOutputPattern`; row and column noise,
  output groups, shutter timing, and detector defect maps compose through
  `CMOSSensor`. Vendor camera profiles are intentionally outside core.
- Frame response: `FrameResponseModel`, `NullFrameResponse`,
  `GaussianPixelResponse`, `SampledFrameResponse`,
  and `RectangularPixelAperture`. These are spatial-domain presampling response
  models. Evaluate the normalized interior, infinite-grid transfer magnitude
  of the realized discrete acquisition kernel with
  `AdaptiveOpticsSim.Detectors.detector_mtf(model, fx, fy)`, where frequency is in cycles
  per detector pixel. Finite frames use zero extension and therefore have
  boundary-dependent response. This diagnostic is not a continuous
  subpixel-aperture MTF; such a model requires an explicitly prepared
  oversampled optical grid.
- Post-collection coupling: `AdaptiveOpticsSim.Detectors.NullChargeCoupling` and
  `AdaptiveOpticsSim.Detectors.InterpixelCapacitance`. Configure these with
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
  `FirstOrderAfterpulseMeanResponse`, `NearestNeighborCountRedistribution`,
  `CompositeCountingMeanResponse`
- Runtime functions: `capture!`, `output_frame`, `channel_output`,
  `detector_export_metadata`, `readout_ready`, `reset_integration!`,
  `thermal_model`, `detector_ramp_slope`, `detector_ramp_intercept`,
  `detector_ramp_cube`, `detector_ramp_times`
- Prepared intensity-map acquisition: exported
  `prepare_detector_acquisition`; qualified-public
  `Detectors.DetectorAcquisitionPlan`,
  `Detectors.PreparedDetectorAcquisition`,
  `Detectors.detector_acquisition_detector`,
  `Detectors.detector_acquisition_input`,
  `Detectors.detector_acquisition_plan`,
  `Detectors.detector_acquisition_state`,
  `Detectors.detector_acquisition_workspace`, and
  `Detectors.detector_acquisition_products`

`capture!(prepared_acquisition; integration_duration=seconds)` is the
frame-step incremental convenience surface. `integration_duration` is a
positive integration duration, not an absolute timestamp or sample period.
Event-driven global-
shutter, rolling-shutter, and frame-transfer acquisition instead use the
qualified plant API above, with exact scheduler-owned exposure/read timestamps
rather than floating accumulated duration as completion authority. The
composed event loop publishes complete-product sequence/readiness state only;
leases, ports, progressive transport delivery, and wall-clock pacing remain
outside this core surface.

Use `bits` for detector quantization depth and `output_type` for the Julia
element type exported to an RTC/HIL boundary. A detector with `bits` must also
provide a fixed positive `full_well`; per-frame peak normalization is not an
ADC model and is rejected.

`Detector(...; qe=...)` accepts either a scalar quantum efficiency or a
qualified QE model such as
`AdaptiveOpticsSim.Detectors.SampledQuantumEfficiency(wavelengths, values)`. Matrix-only
capture uses the scalar `params.qe` value, which is the supplied scalar or the
peak sampled QE. Source-aware capture, `capture!(det, image, src; rng=...)`,
evaluates the QE model at `wavelength(src)`. For `SpectralSource`, it uses the
flux-weighted effective QE over the spectral bundle. Pyramid frame-detector
paths specialize this boundary by applying sampled QE per wavelength before
incoherent optical-rate accumulation. Prepared diffractive Shack–Hartmann
optics instead retain distinct wavelength-rate products in a bundle so
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
the host. Detector-aware Pyramid, Bi-O-edge, and Zernike reference frames apply
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
returned exact prepared owner to `capture!`:

```julia
acquisition = prepare_detector_acquisition(detector, intensity_map)
frame = capture!(acquisition, rng)
```

The plan owns only run-immutable metadata, shapes, radiometric scaling, and QE.
The prepared owner binds the exact detector and input storage to persistent
state, replaceable workspace, and caller-visible products. Replacing a bound
owner or violating a forbidden alias requires preparation again; the removed
detector/input/plan call shape has no compatibility adapter. Photon-rate maps
cannot be rescaled;
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
The incremental form is
`capture!(acquisition; rng, integration_duration=seconds)`; the duration is not
an absolute timestamp or a sampling period.

For whole-frame acquisition, the common ordered stages are: prepared
radiometric scaling; presampling response; physical-pixel integration;
QE/exposure; configured binning; signal defects and persistence; photon and
background-flux statistics; dark current, DSNU, and sensor-generated charge;
nonlinearity, saturation, charge transfer or pre-read multiplication; IPC;
read noise; post-read gain and sensor output effects; readout correction;
quantization and configured background-map subtraction; then windowing and
output conversion. Unsupported or unconfigured sensor-specific stages are
typed no-ops. IPC therefore does not contribute to `detector_mtf`, and output
conversion does not redefine the charge-domain full-well limit.

For the lowest-cost deterministic RTC load, configure `NoiseNone()` and
`NullFrameResponse()` and omit optional defect, coupling, readout-product, and
output-conversion effects. The prepared path remains type-inferred and
allocation-free after warmup. This is still the ordinary detector pipeline,
not a separate reduced-order detector type or a camera profile.

`SPADArrayDetector((rows, columns); ...)` is the maintained fixed-shape
Geiger-mode area-counting surface. Its deterministic order is cell-integrated
photon-arrival rate times active-area detection efficiency, fill factor, live
integration time, and source throughput; plus integrated dark counts; followed
by the selected dead-time mean law and deterministic mean-response stages. An
optional `NoisePhoton()` stage then draws a Poisson count from that adjusted
expectation. This last step is a bounded accumulated-count surrogate, not the
true non-Poisson distribution produced by dead time or afterpulsing.
Mean-response composites may contain at most one first-order afterpulse stage
and one nearest-neighbor redistribution stage, so their scalar export metadata
remains exact. Direct capture validates finite nonnegative input and exact fixed
dimensions before
mutation; prepared WFS acquisition binds those conditions once and remains
device resident. Preparation validates the initial values; repeated acquisition
trusts the bound producer to preserve finite nonnegative samples. Fill factor is
scalar radiometry and does not create a pixel aperture or detector MTF. The
model does not emit photon timestamps or avalanche events.

`MKIDArrayDetector` is the maintained MKID surface for accumulated-count HIL
use. Its observable is one expected-count or sampled-count image. The shared
pipeline applies quantum efficiency, fill factor, live integration time, source
throughput, dark counts, an optional dead-time mean law, optional Poisson
sampling, and output conversion. It does not apply SPAD afterpulse or
nearest-neighbor redistribution stages.

Optional physical characteristics are grouped in
`MKIDArrayCharacteristics`. `energy_resolving_power` is the dimensionless ratio
`E/ΔE`; `photon_arrival_time_resolution_s` is in seconds; and
`wavelength_passband_m=(minimum, maximum)` is an inclusive interval in meters.
`MKIDArrayExportMetadata` keeps these characteristics separate from its
`observable=:accumulated_count_image` declaration and explicitly reports that
energy estimates and photon-arrival timestamps are absent. Source-aware capture
applies the passband, including weighted `SpectralSource` bundles; matrix-only
capture assumes spectrally prefiltered input. Prepared WFS paths compute source
throughput once during preparation. No detector spatial response or MTF is
modeled or inferred: inputs are already cell-integrated per detector pixel, and
scalar fill factor does not define a pixel aperture or presampling response.

`CMOSSensor` covers CMOS, sCMOS, and quantitative low-noise CMOS architectures
through composition rather than camera classes. `NoiseReadout` supplies a
uniform independent component; `row_readout_sigma`, `column_readout_sigma`,
and `CMOSReadNoiseMap` add structured components at the readout stage.
`PixelResponseNonuniformity`, `DarkSignalNonuniformity`, `BadPixelMask`, and
qualified-public `Detectors.StaticCMOSOutputPattern` carry measured calibration
structure. Core provides no vendor defaults or named cameras.

`InGaAsSensor` is a product-neutral InGaAs area-sensor model. Its `glow_rate`
is an independent Poisson expectation rate in electrons per pixel per second;
the model applies it over integration duration and does not infer a read-cadence
glow law. Ordinary `Detector.dark_current` remains a separate rate. InGaAs
technology does not select a presampling response: the default is
`NullFrameResponse`, and an aperture, sampled response, or other detector MTF
must be configured explicitly. [Astronomical InGaAs characterization](https://arxiv.org/abs/1307.1469)
shows that dark current, readout glow, read noise, nonuniformity, and other
effects depend on the particular focal-plane array and operating condition,
so core supplies no family calibration or camera profile.

`ExponentialPersistence(coupling, decay)` is an optional deterministic
frame-to-frame reduced model. For `InGaAsSensor`, the latent state is updated
as `l[k+1] = decay*l[k] + coupling*q[k]`, where `q[k]` is the charge-domain
frame after nonlinearity, saturation, and charge coupling but before read
noise, conversion gain, correction, quantization, and background subtraction.
The latent state is added to the next frame. Both coefficients are
dimensionless and applied once per completed frame; the model has no elapsed-
time constant, trap population, exposure-dependent decay, or calibrated
persistence claim. Configure `NullPersistence()` when this approximation is
not intended.

Qualified-public `Detectors.FrameWindow` is an output-product selection
applied after completion of the full modeled frame. It does not change
row-band exposure or acquisition timing; use physically cropped detector
dimensions when the active sensor geometry itself is smaller.

`CCDSensor()` is the conventional single-read CCD model. Its
`clock_induced_charge_per_frame` parameter is an independent Poisson
expectation in electrons per pixel per completed frame, so it is not scaled by
exposure duration. A configured presampling response supplies detector MTF;
`CCDSensor` does not infer one from CCD technology. For `SingleRead()`,
`sample_duration` must remain zero: the parameter belongs only to
`SkipperSampling`, while Plant acquisition definitions own conventional
readout and readiness durations. Core does not infer a serial-register timing
model.
The optional thermal field is named `cic_per_frame_law` for the same
dimensional reason.

`SkipperSampling(n)` configures repeated nondestructive CCD sampling. The
implementation accumulates a mean online and exposes `SkipperReadoutProducts`
without retaining an `n`-plane read cube. This bounds memory independently of
sample count and keeps the warmed repeated-capture path allocation-free. The
current model assumes independent read samples. Configure CCD clock-induced
charge with `clock_induced_charge_per_frame`; unlike dark current, it is not
scaled by integration time. `CCDSensor.sample_duration` is the duration of one
configured full-frame nondestructive sample. It is not a sampling period or an
electronics integration-window model.

`UpTheRampSampling(n)` is available on `HgCdTeSensor` and on the
linear-avalanche HgCdTe sensor. It retains `n` evenly spaced
nondestructive reads, fits an intercept and slope, and returns
`slope * integration_time` so ordinary `capture!` output units do not change.
Use `detector_ramp_slope(det)`, `detector_ramp_intercept(det)`,
`detector_ramp_cube(det)`, and `detector_ramp_times(det)` for diagnostics.
`Detectors.detector_ramp_acquisition(det)` returns
`:synthesized_final_charge` for this direct convenience or
`:scheduled_evolving_charge` for a Plant event acquisition.
Reads start at exposure time zero and end at the configured integration time;
`read_time` must fit within the resulting cadence. Full-frame and windowed
repeated capture reuse
their products after warmup.

The current ramp model assumes linear accumulation and independent per-read
Gaussian read noise. It shares the exposure's photon/dark realization across
the ramp by synthesizing fractional reads from the completed frame. It is a
post-exposure lower-fidelity convenience, not a time-resolved simulation of an
evolving charge ramp. The scheduled detector path instead records actual
nondestructive-read events at its prepared, possibly nanosecond-quantized
timestamps. The current model also does not provide cosmic-ray segmentation,
saturation-aware fitting, or correlated 1/f-noise estimation.

`EMCCDSensor(...; acquisition_mode=FrameTransferAcquisition(...))` models frame
transfer as timing only. With `readout_rate_hz` configured, metadata reports
the pixel-read duration, one-frame output latency in `sampling_wallclock_time`,
and the overlapped `steady_state_frame_period`. `SequentialAcquisition()` is
the default. Neither acquisition policy changes the presampling response or its
derived MTF, QE, charge multiplication, or detector noise.

The qualified-public EMCCD multiplication policies live in `Detectors`.
`ClippedGaussianMultiplicationApproximation` is the default nonnegative
moderate/high-charge approximation and is available on maintained CPU and GPU
backends. `ConditionalGammaMultiplication` implements the conditional-Gamma
law on CPU and is rejected during accelerator detector construction rather
than changing algorithms by backend. `em_gain_range` is enforced for
`EMOutput`. If both detector `full_well` and `register_full_well` are set, the
smaller input-referred limit applies. `ConventionalOutput` bypasses EM gain and
its gain-range constraint.

`PhotonCountingEMMode` thresholds the post-EM, post-read-noise frame and then
applies its Bernoulli `detection_efficiency`. Its output is binary per pixel per
frame and therefore coincidence-limited; it is not photon-number resolving
within one frame. The internal capability traits report photon counting and
photon-number resolution separately.

## Wavefront Sensors

Reusable physical components in this section—modulation models, microlens
arrays, masks, phase spots, and defocus pairs—are imported from
`AdaptiveOpticsSim.Optics`. Import the common sensing, product, prepared-stage,
and capability vocabulary from its canonical owner:

```julia
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Optics
using AdaptiveOpticsSim.WavefrontSensors
```

`WavefrontSensors` owns the complete Shack–Hartmann, Pyramid, Bi-O-edge,
Zernike, and Curvature families as well as the common contracts.
`MicrolensArray`, `PyramidPhaseMask`, `BiOEdgeAmplitudeMask`,
`ZernikePhaseSpot`, and `CurvatureDefocusPair` remain independent `Optics`
components. Migrated bindings are not forwarded through the root.

- Sensing modes: `Diffractive`, `Geometric`
- Prepared products: `WFSObservationMetadata`, `WFSMeasurementMetadata`,
  `WFSObservation`, `WFSMeasurement`
- Product accessors: `observation_storage`, `observation_units`,
  `observation_metadata`, `measurement_storage`, `measurement_units`,
  `measurement_metadata`
- Prepared stage protocols: `prepare_wfs_optics` /
  `form_wfs_optical_products!`, `prepare_wfs_acquisition` /
  `acquire_wfs_observation!`, and `prepare_wfs_estimation` /
  `estimate_wfs_measurement!`
- Qualified plan and owner API:
  `WavefrontSensors.AbstractWFSOpticsPlan`,
  `WavefrontSensors.AbstractWFSAcquisitionPlan`,
  `WavefrontSensors.AbstractWFSEstimationPlan`, the generic detector
  acquisition plan/prepared-owner types, and
  `WavefrontSensors.wfs_acquisition_plan`. These names are public but not
  exported; ordinary workflows use the prepared stage protocols above.
- Estimation paths: `AbstractWFSMeasurementPath`,
  `AcquiredObservationPath`, `DirectMeasurementPath`,
  `wfs_measurement_path`
- Contract failure: `WFSPreparationError`, whose `stage` and open
  extension-defined `reason` fields identify rejected preparation contracts or
  execution-time prepared-binding violations before mutation
- Module-owned WFS families: `ShackHartmannWFS`, `PyramidWFS`, `BiOEdgeWFS`,
  `ZernikeWFS`, `CurvatureWFS`
- Zernike optical composition: `Optics.ZernikePhaseSpot`,
  `ZernikeOpticalFrontEnd`, `zernike_rate_map`, and
  `set_zernike_calibration!`
- Shack-Hartmann optical composition: `Optics` owns `MicrolensArrayParams`,
  `MicrolensArray`, `prepare_microlens_propagation`, and `microlens_array`;
  `WavefrontSensors` owns `ShackHartmannDirectFrontEnd`,
  `ShackHartmannOpticalFrontEnd`, `shack_hartmann_optics`, and
  `shack_hartmann_rate_map`. The
  microlens array is the immutable regular-array model and its numerical
  sampling policy. `MicrolensPropagationPlan` is run-immutable;
  `MicrolensPropagationWorkspace` owns backend/grid-bound FFT handles and
  reusable optical scratch; `PreparedMicrolensPropagation` pairs them for one
  single-writer execution. A `ShackHartmannOptics` composes
  that prepared owner with a propagation-free diffractive front-end
  definition, and can be assembled without constructing or retaining a
  `ShackHartmannWFS`. `ShackHartmannWFS.front_end` is always the physical
  definition; `ShackHartmannWFS.optics` is `nothing` for geometric sensing
  and the explicit execution composition for diffractive sensing. The
  superseded
  top-level `microlens_array`, `optical_workspace`, and `layout` fields are not
  emulated, and `microlens_array` accepts the owning front end rather than the
  whole sensor.
  The plan, workspace, prepared-owner, and accessor vocabulary is
  qualified-public from `Optics`; callers normally obtain it through the
  preparation function.
- Curvature optical composition and readout: `Optics.CurvatureDefocusPair`,
  `CurvatureOpticalFrontEnd`, `curvature_rate_maps`,
  `CurvatureReadoutModel`, `CurvatureFrameReadout`,
  `CurvatureChannelReadout`, `CurvaturePackedAcquisition`,
  `CurvatureBranchResponse`, and `set_curvature_calibration!`
- Shack-Hartmann calibration and extraction: `FluxThresholdValidSubapertures`,
  `AbstractSlopeExtractionModel`, `CenterOfGravityExtraction`,
  `SubapertureLayout`,
  `SubapertureCalibration`, `subaperture_layout`,
  `subaperture_calibration`, `slope_extraction_model`,
  `set_subaperture_calibration!`, `valid_subaperture_indices`,
  `n_valid_subapertures`
  (`subaperture_layout` likewise accepts the owning Shack-Hartmann front end,
  not `ShackHartmannWFS`)
- WFS normalization policies:
  `MeanValidFluxNormalization`, `IncidenceFluxNormalization`
- Measurement and detector products: `measure!`,
  `pyramid_modulation_frame!`, `valid_subaperture_mask`,
  `wfs_detector_image`, and `shack_hartmann_spot_cube`
- Prepared WFS optical-product access: the qualified-public
  `WavefrontSensors.wfs_optical_products`
- LiFT forward and observation contracts: `PreparedLiFTForwardModel`,
  `prepare_lift_forward_model`, `lift_forward_output`,
  `evaluate_lift_forward!`, `predict_lift_observation!`,
  `LiFTObservation`, `lift_observation_contract`,
  `LiFTIdentityMapping`, `LiFTFrameMapping`, `LiFTPhotonRate`,
  `LiFTExpectedCounts`, and `LiFTNormalizedIntensity`
- LiFT estimation: `LiFT`, the qualified-only
  `WavefrontSensors.reconstruct!` and `WavefrontSensors.reconstruct`
  generics, `diagnostics`,
  `LiFTSolveAuto`, `LiFTSolveQR`, `LiFTSolveNormalEquations`,
  `LiFTLevenbergMarquardt`, and `LiFTAdaptiveLevenbergMarquardt`

`WavefrontSensors.wfs_optical_products(prepared)` returns the exact typed
photon-arrival-rate product or products bound by a prepared WFS optics owner.
It never implies detector acquisition. After detector-coupled acquisition,
`wfs_detector_image(wfs, detector)` returns the detector's actual output only
when that output is two-dimensional. Use `Detectors.output_frame(detector)`
for non-image detector products and `observation_storage(observation)` for a
prepared `WFSObservation`. No detector-free camera-frame accessor is provided.

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
Pyramid, Bi-O-edge, Zernike, and Curvature implement the generic contract. LiFT
intentionally remains outside the ordinary `AbstractWFS` hierarchy: prepare
its focal-plane model independently, bind caller-owned acquired values with
`LiFTObservation`, and then run `WavefrontSensors.reconstruct!`. This generic
is deliberately distinct from the control-reconstruction generic exported at
the package root. Modal selection is a cold-path
`LiFT(...; mode_ids=...)` choice rather than a per-frame argument. The forward
output is a cell-integrated photon-arrival rate; count or normalized
observations require an explicit observation-domain conversion. Neither
forward evaluation nor estimation reads telescope cadence or invokes a
detector.

LiFT modal arrays are dimensionless OPD shapes. Reconstructed coefficients and
the prepared diversity OPD are in metres; their assembled modal sum is the OPD
map supplied to the focal-plane model.

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

`PyramidOpticalFrontEnd` and `BiOEdgeOpticalFrontEnd` bind their distinct
physical masks, a prepared `NoModulation`, `CircularModulation`, or
`SampledModulation` policy, reusable propagation storage, and a source. Use
`pyramid_rate_map` or `bi_o_edge_rate_map` to allocate the caller-owned
four-pupil rate product, then prepare and execute it with the generic optical
stage functions. Spectral sources return an `OpticalProductBundle`. An
`Asterism` or `ExtendedSource` requires one explicit pupil input per rendered
direction and also returns a bundle; no direction-dependent pupils are added
by array index. Install acquired-estimator references with
`set_pyramid_calibration!` or `set_bi_o_edge_calibration!`. The setters validate
finite storage and advance a revision, so a prepared estimator must be rebuilt
after recalibration. Acquired differential estimation requires a real, square
`:four_pupil_mosaic` on the same backend and compute device as the estimator.
The detector samples may be floating-point electron counts or integer ADU/count
values; arithmetic is performed in the floating-point estimator precision so
unsigned subtraction cannot wrap. Detector sampling and binning may reduce the
mosaic while preserving complete pupil images. Prepared optical plans also bind
the propagation-sampling revision and reject execution after a legacy sampling
resize, before mutating the rate output. Geometric variants instead consume an
explicit `PupilFunction` and report a floating-point
`:metre`/`:geometric_slopes` direct measurement; they have no focal-plane or
detector workspace.

The legacy constructor keyword `calib_modulation` prepares the broader optical
quadrature used only to select valid estimator support. A zero-aberration
reference is formed with the operating modulation so it subtracts the static
response of the sensor that will actually acquire data. A user-sampled
modulation path is retained for both operations rather than replaced by a
generated circle.

`ZernikeOpticalFrontEnd` separates the immutable phase-shifting
`ZernikePhaseSpot` and prepared re-imaged-pupil propagation from detector
acquisition and the referenced pupil estimator. `zernike_rate_map` allocates a
caller-owned photon-arrival-rate plane; its acquired observation uses the
`:zernike_pupil_image` layout. Estimation writes a dimensionless
`:normalized_pupil_signal` and accepts real floating-point or integer detector
samples. `set_zernike_calibration!` installs the reference, wavelength, and
optical signature atomically; prepared estimators reject a later calibration
revision before output mutation.

`CurvatureOpticalFrontEnd` forms a fixed two-element tuple ordered as positive
then negative defocus. Each product is a separate normalized-pupil-coordinate,
cell-integrated photon-arrival-rate plane. A concrete pair of ordinary detector
plans permits independent response models, QE, exposure durations, stochastic
effects, and readout for the two branches. `CurvaturePackedAcquisition` instead
maps both compatible branches to one detector: `CurvatureFrameReadout` uses a
`:curvature_branch_regions` observation and `CurvatureChannelReadout` uses
`:curvature_branch_channels`. Packed branches must have identical geometry,
radiometry, backend, device, numeric type, and detector-owned exposure duration;
different branch exposures require separate detectors. Estimator preparation
accepts either representation and requires explicit `branch_rate_scales` when
acquisition durations or deterministic gains must be removed before the
calibrated branch difference. `CurvatureBranchResponse` is a legacy optical
relay throughput/background approximation, not detector MTF or acquisition
response; detector response remains downstream and branch specific.

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

## Calibration

Import the maintained calibration surface explicitly:

```julia
using AdaptiveOpticsSim.Calibration
```

- Interaction/control matrices: `InteractionMatrix`, `interaction_matrix`,
  `ControlMatrix`. The caller-owned in-place seam is currently internal and
  available as `AdaptiveOpticsSim.Calibration.interaction_matrix!`.
- Modal bases: `ModalBasis`, `KLDMModes`, `KLHHtPSD`,
  `kl_modal_basis`, `modal_basis`, `basis_from_m2c`
- AO calibration: `AOCalibration`, `ao_calibration`, `control_matrix`
- Error and optical-gain calibration: `fitting_error`, `GainSensingCamera`,
  `calibrate!`, `compute_optical_gains!`
- Misregistration identification is an experimental qualified workflow under
  `AdaptiveOpticsSim.Calibration`, including `MetaSensitivity`,
  `compute_meta_sensitivity_matrix`, `estimate_misregistration`, `SPRINT`,
  and `estimate!`. These names are not yet marked as stable public API.
  `MetaSensitivity` is the structured
  result: it retains the reference interaction matrix, sensitivity operator,
  finite-difference validation steps, and ordered parameter names
- Structured configuration snapshots use qualified
  `AdaptiveOpticsSim.config_dict` and `AdaptiveOpticsSim.snapshot_config`

## Control

Import the maintained control surface explicitly:

```julia
using AdaptiveOpticsSim.Control
```

- Reconstructors: `NullReconstructor`, `ModalReconstructor`,
  `FactorizedReconstructor`, `MappedReconstructor`,
  `ControlledReconstructor`, `reconstruct!`, `reconstruct`
- Controller: `DiscreteIntegratorController`. `ControlledReconstructor`
  composes a reconstructor and stateful controller without adding a runtime
  branch.

Core calibration and configuration APIs perform no implicit filesystem I/O.
They accept no cache path or serialization policy and do not select a file
format. Persist a returned `MetaSensitivity`, estimate, or configuration
dictionary explicitly in caller code or through an optional format extension.
The optional JSON3 extension is one such policy; it is not a core dependency or
the canonical representation of these results.

## Runtime And HIL

- HIL-neutral orchestration: the `AdaptiveOpticsSim.Plant` definitions,
  prepared owners, command lifecycle, triggers, detector lifecycles, and event
  loop documented above
- Independent `AdaptiveOpticsSim.Control` primitives: `VectorDelayLine`,
  `shift_delay!`, `DiscreteIntegratorController`, and `reconstruct!`; command
  application remains the `AdaptiveOpticsSim.Optics.set_command!` operation
- WFS preparation helper: `prepare_runtime_wfs!` prepares the retained WFS
  family-specific scratch required by explicit model loops
- Generic timing helper: `runtime_timing`
- Coarse offline execution is imported explicitly with
  `using AdaptiveOpticsSim.Ensembles`
- Execution policies: `AbstractExecutionPolicy`, `SequentialExecution`,
  `ThreadedExecution`, `BackendStreamExecution`, `DeterministicExecution`,
  `AcceleratedKernelsExecution`, and `DaggerExecution`
- Coarse independent work: `SimulationEnsemble`; qualified access through
  `Ensembles.run_ensemble!`, `ensemble_members`,
  `execution_policy`, `ensemble_ownership_roots`,
  `init_ensemble_scheduler`, and `execute_ensemble!`

Model-specific packages may define their own `prepare!`, `step!`, and `readout`
methods for a fixed offline composition. Core does not provide a generic
single-optic or packed-command closed-loop runtime. Independent physical optics
remain independent Plant command endpoints; use prepared controller-output
routing when one RTC buffer supplies several endpoints.

`SimulationEnsemble` is for independent offline sweeps or model instances.
Dagger, AcceleratedKernels, Julia threads, and backend streams are explicit
coarse policies, not the HIL event-loop scheduler.

For external Proper science-arm integration, use
[`proper-integration-guide.md`](./proper-integration-guide.md).

## Tomography

Import the maintained tomography surface explicitly:

```julia
using AdaptiveOpticsSim.Tomography
```

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

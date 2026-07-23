# Maintainer Architecture Guide

Status: active

## Purpose

This document describes the current implemented architecture of
AdaptiveOpticsSim.jl.

It is the maintainer-facing synthesis layer for:

- the public workflow docs in [`user-guide.md`](user-guide.md)
- the symbol inventory in [`api-reference.md`](api-reference.md)
- the canonical terminology in [`glossary.md`](glossary.md)
- the extension seams in [`extension-guide.md`](extension-guide.md)
- the runtime flow in [`runtime-dataflow.md`](runtime-dataflow.md)

## Architectural Principles

The implemented system follows these package-level rules:

- use multiple dispatch and traits instead of OO-style inheritance
- separate params from mutable runtime state
- make hot paths explicit with mutating `!` functions
- keep core algorithms backend-generic where practical
- isolate backend policy away from model code when possible
- keep deterministic validation and frozen-reference evidence first-class
- keep optional integrations outside core or in extension modules
- traverse column-major arrays with the first index in the innermost CPU loop;
  use `eachindex` for linear elementwise passes
- use native `@simd` only on measured, independent iterations with correctness
  and allocation evidence; do not add explicit SIMD packages speculatively

## Package Shape

The codebase is organized into subsystem directories:

- Directory names should be lower-case. Use lower-case subsystem names such as
  `src/wfs`, `src/detectors`, and `src/control`; reserve CamelCase for Julia
  types and constructors, not filesystem directories.
- `src/core`
  - traits, errors, backend services, reduction/random helpers, shared low-level
    infrastructure
- `src/optics`
  - telescope, sources, explicit optical products, propagation, masks, DMs,
    and prepared direct imaging
- `src/atmosphere`
  - Von Karman/Kolmogorov screens, multilayer atmospheres, infinite-screen
    evolution, atmosphere-field propagation support
- `src/wfs`
  - Shack-Hartmann, Pyramid, BioEdge, Curvature, Zernike, grouped execution,
    calibration scaffolding
- `src/plant`
  - the real `AdaptiveOpticsSim.Plant` submodule, with explicit imports from
    already-loaded core, optics, atmosphere, detector, and WFS domains;
    immutable topology, prepared execution, virtual time/events, commands,
    providers, and mutable single-writer owners remain separate layers
- `src/detectors`
  - frame/counting detectors, response models, thermal models, readout pipeline
- `src/calibration`
  - interaction matrices, modal bases, control matrices, AO calibration bundles,
    gain sensing, misregistration sensitivity, and caller-owned calibration
    storage seams
- `src/control`
  - reconstructors, runtime construction, execution, output planning, timing,
    controllers, and dense/factorized composable control operators
- `src/simulation`
  - compact simulation assembly types and maintained scenario-builder helpers
- `ext`
  - optional backend or ecosystem integrations

The root package exports `Plant` but no forwarding aliases for plant-owned
bindings. `src/plant/api.jl` exports routine construction/execution vocabulary,
marks stable advanced records, policies, state/workspace owners, traits, and
accessors `public`, and leaves underscore-prefixed implementation machinery
unmarked. Extensions add methods to the canonical `AdaptiveOpticsSim.Plant`
function rather than introducing a second root function identity.

## Main Data Model

The dominant pattern is:

- immutable params struct
- mutable state struct
- top-level domain object holding params plus state

Examples:

- `PlantDefinition` with explicit `ControllableOpticID`, `CommandEndpointID`,
  `PlantCommandSchemaID`, `OpticalPathID`, and `AcquisitionID` values, immutable
  optic/path/acquisition topology bindings, and one immutable semantic command
  schema per endpoint. It owns no mutable command state, prepared workspace,
  schedule, RNG stream, queue, transport, or HIL descriptor. Controllable-optic,
  path, and acquisition model types must explicitly opt in to the
  configuration-only `plant_model_definition_style` contract; live
  optic/detector/runtime owners fail closed. Its separately owned telescope and
  atmosphere retain their established state semantics, so structural
  immutability of the plant is not a deep-freeze claim for those owners
- `PreparedCommandEndpoint` as one run-immutable schema/capacity/window/
  ordinal/payload-storage-backend binding, with separate single-writer
  `CommandEndpointState`
  and caller-owned `CommandDispositionWorkspace`. Endpoint-owned payload slots,
  one array staging slot, a flat sorted future calendar, and modulo accepted-
  sequence history are fixed at preparation. Admission copies caller payloads,
  never backdates, and assigns one presentation identity whose eventual
  terminal disposition is rejected, superseded, applied, or failed. Scalar
  payload slots are host-resident; array slots use their prepared backend. One
  opaque isbits application-ready claim may be outstanding and is revalidated
  against endpoint-owned metadata before payload access or completion
- qualified `CommandApplicationState` as the unique separate single-writer
  owner bound to one exact endpoint-state instance before first successful
  admission, with an explicit initial/effective command, preallocated array
  staging, optional copied safe command, application timestamp, and replayable
  silence latch. `apply_claimed_plant_command!` transactionally implements
  absolute or incremental semantics, finite/range policy, and exactly one
  terminal disposition. Safe/fail silence transitions use exact plant-time
  deadlines; equal-time commands take precedence. Direct use owns no physical
  optic mutation, while the prepared event loop stages this state together
  with the model-specific physical response
- qualified `PreparedControllableOptic` as the immutable model-specific
  preparation plus independent endpoint-slot binding for one physical device.
  Its separately constructed state owns visible and staged physical response;
  its workspace owns scratch. Staging is fallible, while commit is a bounded,
  nonthrowing publication step so explicit multi-optic transactions cannot
  become partially visible
- `PreparedPlant` as a schedule-free concrete tuple of
  `PreparedControllableOptic`, prepared command-endpoint,
  `PreparedPathExecutor`, and `PreparedAcquisitionOwner` values. Every declared
  endpoint has one explicit `CommandEndpointConfiguration`; canonical identity
  determines endpoint ordinals and optic order. Each path owns one explicit
  input/result pair and prepared optical workspace; each acquisition borrows
  that exact result read-only and binds one immutable
  `PreparedAcquisitionProvider` choice. Full-optical providers own independent
  detector/WFS execution state, reduced-order providers own their declared
  parameter/state split, and synthetic/replay providers own fixed payload
  snapshots or bounded replay state. Every style writes one caller-owned
  `AcquisitionProducts` contract. The plant also owns exact
  stateful RNG groups derived from its required run seed, derivation version,
  and stable owner identities; selected execution references those groups
  directly. `PathResultKey` performs cold source/optics/output/revision/backend/
  device compatibility checks without putting IDs, shapes, rates, or device
  ordinals in type parameters
- `PreparedPlantEventLoop` as the HIL-neutral serial virtual-time composition
  of independent endpoint command generators, detector or direct-measurement
  acquisition lifecycles, optical sample schedules, one atmosphere timeline
  for due full-optical work, and canonical physical-optic application. Plant
  state is right-continuous and exposure intervals are half-open. Explicit
  `PlantCommandTransaction` values provide all-or-none admission and physical
  publication across distinct optics; equal time or placement alone does not.
  The current Gate 4 surface operation treats every prepared optic as one
  common co-conjugated group on every path, pending explicit placement and
  visibility
- `Telescope` with immutable `TelescopeParams` and a revisioned prepared
  `TelescopeAperture`; it owns spatial geometry and intensity reflectivity but
  no mutable OPD, cadence, or exposure duration
- caller-owned `PupilFunction`, `ElectricField`, and `IntensityMap`
  products with immutable `OpticalPlaneMetadata`, including normalization,
  spatial-measure, and coherent/incoherent combination policy. A
  `PupilFunction` snapshots aperture support and field amplitude and owns the
  mutable OPD for exactly one optical path
- `ShackHartmannWFS` composed through one explicit `front_end`, separate
  calibration, acquisition, and estimator owners. The front end is a
  propagation-free `ShackHartmannDirectFrontEnd` for geometric sensing or a
  `ShackHartmannOpticalFrontEnd` containing the immutable `MicrolensArray`,
  backend/grid-bound `PreparedMicrolensPropagation`, and layout for
  diffractive sensing. There are no top-level optical field aliases or a
  whole-WFS optical-owner union
- `Detector` with `DetectorParams` and `DetectorState`
- `DetectorAcquisitionPlan` as the cold-path compatibility and buffer contract
  between one frame detector and one immutable intensity-map description
- `DirectImagingPlan` and `DirectImagingWorkspace` as the fixed-storage,
  single-writer native image-formation contract; composition returns concrete
  prepared values accessed through `direct_imaging_output` and
  `direct_imaging_components` rather than exposing telescope-owned focal state
- caller-owned `WFSObservation` and `WFSMeasurement` products with explicit
  units, layout/kind, shape, numeric type, backend, and physical-device metadata
- concrete prepared WFS optical-formation, acquisition, and estimation plans
  connected through six dispatch functions rather than a universal stage graph
- direct `execute_path!` and `execute_acquisition!` dispatch over concrete
  prepared owners; there is no abstract executor collection, closure field,
  queue, task, or scheduler at this boundary
- `AdaptiveOpticsSim.Plant` with cold definitions, prepared owners,
  independently timed command/acquisition endpoints, and one bounded serial
  virtual-time event loop
- prepared controller-output routing that borrows named controller products
  and binds them to exact endpoints without packing or timing semantics

This gives:

- explicit configuration
- stable memory ownership
- hot-path mutation without repeated allocation

The generic WFS stage protocol exists independently of the `AbstractWFS`
object layout. Shack-Hartmann separates microlens formation, acquisition, and
estimation over caller-owned products. Its geometric and diffractive signals
share one explicit
`[axis 1; axis 2]`, Julia-column-major lenslet convention; OOPAO row-major
reference adaptation remains in the test harness. Its geometric mode declares
`DirectMeasurementPath()` and allocates no placeholder optical or acquisition
workspace. Microlens sampling, synchronized subaperture layout, and
calibration are cold configuration: maintained mutation advances a revision,
and prepared plans reject stale bindings while caller-owned product contents
remain mutable.

This is a breaking refactor. Superseded public and internal representations are
removed and callers are migrated directly; synthetic property forwarding,
state views, deprecated aliases, and permanent compatibility adapters are not
part of the maintained architecture.

Pyramid and BioEdge use separate `PyramidPhaseMask` and
`BioEdgeAmplitudeMask` physical front ends. They share prepared focal-plane
modulation only where the optical quadrature is identical; its normalized
weights average intensity and never integrate detector time. Each front end
writes a normalized-pupil-coordinate, cell-integrated photon-arrival-rate
four-pupil mosaic or a typed spectral/path-local bundle. Generic detector
acquisition applies response, QE, and duration afterward. Their differential
estimators own valid support, normalization, reference subtraction, optical
gain, and a calibration revision that invalidates stale prepared plans.
Geometric Pyramid and BioEdge declare `DirectMeasurementPath()` and construct
neither propagation nor acquisition workspace.

Zernike and Curvature now follow the same ownership boundary. Their optical
front ends own only physical descriptions and single-writer prepared
propagation state; detector acquisition owns observation state; and estimators
own valid support, reference state, normalization, and calibration revisions.
Curvature exposes a fixed positive-/negative-defocus rate tuple that can feed
two independent detectors or an explicitly packed single-detector mapping.
Convenience execution coordinates these same explicit component owners; no
synthetic state view or monolithic optical-owner adapter is retained.
LiFT separately owns a prepared focal-plane forward model, caller-provided
`LiFTObservation`, and iterative estimator state. Its prepared modal subset and
observation contract are cold-path bindings; repeated estimation neither owns
nor triggers detector acquisition.

## Execution Layers

There are three main execution layers:

### 1. Model-local computation

Examples:

- `advance_by!(atm, elapsed_seconds)` or `advance_to!(atm, model_time)`
- `render_atmosphere!(pupil, renderer, atm, epoch)`
- `measure!(wfs, pupil, src, det)`
- `update_surface!(optic)` followed by
  `apply_surface!(pupil, optic, DMAdditive())`

These functions should own local physics and algorithm behavior, but not broad
backend policy.

Timed atmosphere models have one evolution writer. Evolution publishes an
immutable `AtmosphereEpoch` token for the current mutable layer state;
prepared, path-local renderers consume only that current epoch and write
caller-owned products. The token is not retained layer storage. A scheduled
executor must materialize due path inputs before the writer advances or bind a
model-specific retained state. Atmosphere state therefore owns physical layers
and timeline state, not a shared pupil-sized render target or a mutable last-
source geometry cache. The timed atmosphere API reads no telescope timing
value. Explicit model loops pass an elapsed duration directly; the Plant event
loop converts its canonical timestamp to model time at the physical-model
boundary. Atmosphere advancement remains independent of detector exposure and
sample period.

### 2. Shared subsystem services

Examples:

- detector pipeline helpers
- prepared compatible-intensity accumulation and typed detector acquisition
- propagation contexts
- backend reductions and random/noise services

These shared services reduce duplicated orchestration logic across SH, Pyramid,
BioEdge, detectors, explicit model loops, and Plant implementations.

### 3. Plant and model orchestration

Examples:

- `PlantDefinition` and `prepare_plant`
- `prepare_plant_event_loop`, `step_plant_events!`, and
  `run_plant_events_until!`
- prepared controller-output routing
- model-specific `prepare!`, `step!`, and `readout` methods for a fixed offline
  composition
- `SimulationEnsemble` for coarse independent work

Core Plant coordinates virtual-time prepared owners without owning wall-clock
pacing or RTC transport. A model-specific loop may coordinate the same
numerical primitives but is not a generic package runtime.

## Runtime Ownership Model

The maintained ownership model is:

- cold definitions are separate from prepared plans and mutable state
- every mutable endpoint, optic, path workspace, acquisition, and RNG owner has
  one writer
- exported products are distinct from scratch buffers
- command endpoints own independent schema, sequence, effective-time, silence,
  bounded calendar, and effective-command state
- controller-output routing borrows exact named products; successful admission
  is the boundary that copies a payload into endpoint-owned storage
- WFS, science, calibration, and coronagraph paths own distinct
  `PupilFunction` or field products; path reuse requires an explicit prepared
  compatibility contract
- WFS and detector pipelines own their sampled/readout/intermediate products
  explicitly rather than relying on in-place aliasing
- optical formation produces photon-arrival rates or explicitly dimensionless
  products; prepared detector acquisition validates metadata, applies
  presampling response before physical-pixel integration, and integrates its
  explicit exposure exactly once
- CPU and accelerator placement is fixed during preparation; synchronization
  occurs only at explicit dependencies or observation/transport boundaries
- the serial Plant event loop is the deterministic HIL-neutral oracle;
  `SimulationEnsemble` and task-graph policies stay outside its deadline path

For the step-by-step view, see [`runtime-dataflow.md`](runtime-dataflow.md).

## Backend Strategy

Backend support is a first-class concern, but backend policy should be
centralized.

The current intended split is:

- `Core` owns backend traits, reduction helpers, random services, and launch
  abstractions
- domain models should call shared helpers rather than embedding backend
  policy in model code
- `ext` modules should add backend-specific adaptations without broadening the
  public workflow API

This package currently maintains:

- CPU
- optional CUDA
- optional AMDGPU

Validation and benchmark evidence for these backends lives in:

- [`model-validity-matrix.md`](model-validity-matrix.md)
- benchmark artifacts under `benchmarks/results/`

## Validation Structure

The package uses four distinct evidence classes:

- analytic and structural checks
- frozen reference-bundle regression against OOPAO and targeted SPECULA cases
- backend smoke and parity checks
- benchmark evidence, including cross-package benchmark archives

The maintained synthesis doc for this is
[`model-validity-matrix.md`](model-validity-matrix.md).

## Optional Boundaries

The package is intentionally conservative about what belongs in core.

Core owns:

- optical model primitives
- atmosphere and sensing models
- detector/readout physics
- calibration and runtime orchestration

Optional or boundary-managed surfaces include:

- external science-path integrations
- HIL orchestration, telemetry transport, artifact lifecycle, and RTC/testbench
  operation; see [`hil-package-boundary.md`](hil-package-boundary.md)
- heavyweight format or ecosystem adapters
- extra plotting or notebook-facing helpers

See the optional-boundary policy in this guide.

## Recommended Reading Order for Maintainers

1. [`documentation-map.md`](documentation-map.md)
2. [`user-guide.md`](user-guide.md)
3. [`maintainer-architecture.md`](maintainer-architecture.md)
4. [`runtime-dataflow.md`](runtime-dataflow.md)
5. [`extension-guide.md`](extension-guide.md)
6. [`model-validity-matrix.md`](model-validity-matrix.md)

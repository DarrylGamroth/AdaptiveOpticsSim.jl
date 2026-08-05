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
  - Shack-Hartmann, Pyramid, Bi-O-edge, Curvature, Zernike, grouped execution,
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

### Namespace Authority

The package completed its transition from the original flat root surface to
real domain modules alongside `Plant`. `Backends`, `Optics`, `Atmospheres`,
`Detectors`, `WavefrontSensors`, `Calibration`, `Control`, `Tomography`, and
`Ensembles` are canonical owners with exact API allowlists. `Atmospheres` owns
atmosphere models and state,
source-direction rendering and batching, and
atmosphere-coupled propagation. `Detectors` is the canonical owner of both
conventional frame/area and counting/channel detector models and acquisition,
with one shared detector type graph. `Optics` owns apertures, telescopes, sources, optical
location/product metadata, pupil and field formation, general Fraunhofer and
Fresnel propagation, direct imaging, sampled OPD, physical NCPA, controllable
optics, and reusable physical WFS components. `Calibration` owns inverse
policy, interaction/control matrices, modal basis construction, model-derived
NCPA synthesis, fitting, optical-gain calibration, and registration
identification. `Control` owns reconstruction and controller composition;
`Tomography` owns guide-star geometry, atmospheric reconstruction, fitting,
and DM-command projection. `Ensembles` owns coarse offline execution policies,
sweeps, and optional scheduler integrations. The
authoritative pre-migration inventory, final exact API allowlists,
extension-hook owners, and cross-module imports are frozen in
[`../test/contracts/namespace_authority.toml`](../test/contracts/namespace_authority.toml).
The completed implementation state is recorded separately in
[`../test/contracts/namespace_migration_state.toml`](../test/contracts/namespace_migration_state.toml).
These files are the equality contract for the migration; documentation
summaries do not override them.

The owner dependency graph is acyclic:

- `Optics` depends on `Backends`
- `Atmospheres` and `Detectors` depend on `Backends` and `Optics`
- `WavefrontSensors` depends on `Backends`, `Optics`, and `Detectors`
- `Plant` depends on those five physical/runtime foundations
- `Calibration` may additionally compose atmosphere, detector, and WFS models
- `Control` depends on `Backends` and `Calibration`
- `Tomography` depends on the physical domains and `Calibration`
- `Ensembles` depends only on `Backends`

Physical NCPA is owned by `Optics`; its runtime value stores an explicit
sampled OPD map rather than a calibration basis or coefficient history.
KL-, Zernike-, and M2C-based synthesis belongs to `Calibration`, while path
visibility belongs to `Plant`. Microlens arrays, focal-plane modulation,
phase/amplitude masks, phase spots, and defocus optics are implemented physical
`Optics` components. Composed WFS front ends,
observations, measurements, detector bindings, and estimators belong to
`WavefrontSensors`.

The final root surface contains only the canonical modules, shared generic
errors, fidelity profiles, RNG services, and `runtime_timing`. It provides no
compatibility aliases. `AbstractSource` and `AbstractTimedAtmosphere` become
qualified-public seams under `Optics` and `Atmospheres`, respectively, because
the companion HIL package implements test models against those interfaces.

## Main Data Model

Cold model configuration continues to use immutable parameter or definition
types. Prepared numerical execution uses a more precise ownership split:

- an immutable definition describes reusable configuration
- a run-immutable plan owns validated numerical or physical execution data
- mutable state owns values that affect later scientific results
- a replaceable workspace owns scratch and execution resources
- caller-visible products remain distinct from scratch
- a prepared execution owner binds the exact plan, state, workspace, products,
  and backend/device/context for one single writer

Discarding state may change the next result. Recreating a workspace must not.
This distinction is stronger than whether a field is mutable. A plan may retain
a logically immutable and reentrant FFT or backend operator; a scratch-owning,
stream-bound, task-bound, or non-reentrant handle belongs to the workspace or
prepared owner instead. Exact caller-array identity normally belongs to the
prepared owner rather than the reusable plan.

Nominal abstract plan types follow the `AbstractFFTs` model only where one
canonical domain owns a real shared protocol. The owning module documents the
required methods, invariants, traits, failure behavior, and conformance tests.
Concrete prepared owners parameterize their plan, state, and workspace fields;
hot execution does not store an abstract plan root. No implementation field or
collection element type is `Any`; cold input is normalized to concrete
ownership before preparation returns. Unconstrained method-signature wildcards
are not stored representations and remain available for dispatch fallbacks.
There is no universal
`AbstractAdaptiveOpticsPlan`, global `Interfaces` namespace, or generic
`process!` verb. Domain operations retain their accepted names.

The target architecture does not directly use Julia's `Memory` as an
implementation or API type. Homogeneous bounded host storage uses a
FixedSizeArrays.jl `FixedSizeArray` whose element type is concrete and whose
final shape is allocated during preparation. A `Vector` may be a cold builder
but is sealed before execution. Because fixed-size arrays still permit element
replacement, exact prepared owners retain an independent ordered membership
snapshot and reject replacement before numerical mutation. Small fixed
heterogeneous composition uses
concrete tuples or unions. Larger topologies group concrete homogeneous owners
by execution family and use compact concrete descriptors or purpose-built
owners without encoding path count in the type. Backend numerical storage
retains its concrete backend array type. FixedSizeArrays.jl may use `Memory` as
an internal backing; package code does not name or depend on that backing.
Boundedness is an owner/lifecycle invariant, and replacing `Memory{Any}` with
`Vector{Any}` does not repair type erasure. The breaking migration and evidence
gates are tracked in
[issue #225](https://github.com/DarrylGamroth/AdaptiveOpticsSim.jl/issues/225).

The following examples include current mixed owners that issue #225 will
decompose; they describe the implementation baseline rather than exceptions to
the target contract.

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
  for due full-optical work, canonical common-pupil surface application, and
  exact path-local autonomous-optic bindings. Plant state is right-continuous
  and exposure intervals are half-open. Explicit `PlantCommandTransaction`
  values provide all-or-none admission and physical publication across
  distinct optics; equal time or placement alone does not. Required
  placement/visibility declarations resolve into bounded canonical per-path
  ranges and co-placed groups; a due path applies only its visible pupil-plane
  members. `AutonomousPathExecutionRole` retains one exact focal-plane
  `AutonomousPeriodicOpticDefinition` coupling. Sampled atmospheric-conjugate
  models prepare an exact NGS/LGS pupil-footprint coupling.
  `Plant.DeformableMirrorModel` prepares the native DM with distinct
  active/staged runtimes, and serial common-plane MCAO plus target-local MOAO
  scenarios exercise the same exact path bindings in full-optical and
  reduced-order execution. Native `NCPA` and `OPDMap` declarations likewise
  prepare run-owned sampled OPD and exact common or selected-path bindings
  before atmosphere-derived pupil products enter controllable and autonomous
  optical execution
- `Telescope` with immutable `TelescopeParams` and a revisioned prepared
  `TelescopeAperture`; it owns spatial geometry and intensity reflectivity but
  no mutable OPD, cadence, or exposure duration
- caller-owned `PupilFunction`, `ElectricField`, and `IntensityMap`
  products with immutable `OpticalPlaneMetadata`, including normalization,
  spatial-measure, and coherent/incoherent combination policy. A
  `PupilFunction` snapshots aperture support and field amplitude and owns the
  mutable OPD for exactly one optical path
- `WavefrontSensors.ShackHartmannWFS`, owned with its calibration and
  estimation implementation by the canonical WFS module. Its `front_end` is
  a propagation-free physical definition: `ShackHartmannDirectFrontEnd` for
  geometric sensing or `ShackHartmannOpticalFrontEnd` for diffractive sensing.
  Diffractive execution composes that definition with a separately prepared
  microlens plan/workspace owner in `ShackHartmannOptics`.
  Replaceable `ShackHartmannWorkspace`, caller-visible
  `ShackHartmannProducts`, and persistent calibration are distinct owners;
  geometric sensing allocates no spot cube, detector-noise cube, or exported
  spot product. There are no top-level optical aliases, property forwarding,
  or whole-WFS optical-owner union
- frame `Detector` owners with immutable `DetectorParams`, persistent
  `DetectorState`, replaceable `DetectorWorkspace`, and caller-visible
  `DetectorProducts`; frame-readout products and workspace follow the same
  separation, including multi-read, up-the-ramp, and Skipper acquisition
- `Detectors.DetectorAcquisitionPlan` as the run-immutable numerical and
  radiometric contract, and `Detectors.PreparedDetectorAcquisition` as the
  exact single-writer detector/input/plan/state/workspace/product binding.
  These advanced ownership values and their `detector_acquisition_*` accessors
  are qualified-public; routine callers prepare once and execute the returned
  owner with `capture!`
- `DirectImagingPlan` as reusable run-immutable image-formation metadata and
  mapping, `DirectImagingWorkspace` as replaceable single-writer scratch, and
  `PreparedDirectImaging` as the exact input/field/output/plan/workspace
  binding; composition seals homogeneous memberships in concrete
  `FixedSizeVector` storage and returns prepared values accessed through
  `direct_imaging_output` and `direct_imaging_components` rather than exposing
  telescope-owned focal state
- caller-owned `WFSObservation` and `WFSMeasurement` products with explicit
  units, layout/kind, shape, numeric type, backend, and physical-device metadata
- concrete prepared WFS optics, acquisition, and estimation plans
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
object layout. Its qualified nominal plan roots are
`AbstractWFSOpticsPlan`, `AbstractWFSAcquisitionPlan`, and
`AbstractWFSEstimationPlan`; each concrete plan is run-immutable and remains a
separate field of an exact single-writer prepared owner. Generic frame,
accumulated-count, and static fan-out acquisition already follow that split
and expose their plan through `wfs_acquisition_plan`. Family optical and
estimator owners migrate under the same contract in PE-05B through PE-05F;
there is no universal stage graph or compatibility owner. Shack-Hartmann and
the Pyramid/Bi-O-edge families now implement that exact split. Shack-Hartmann
separates microlens optics, acquisition, and estimation over caller-owned
products. Its geometric and diffractive signals share one explicit
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

`WavefrontSensors.PyramidWFS` and `WavefrontSensors.BiOEdgeWFS` own their
composition, acquisition, calibration, and estimation implementations.
`Optics` independently owns their `PyramidPhaseMask` and
`BiOEdgeAmplitudeMask` physical components. The sensor families share prepared
focal-plane modulation only where the optical quadrature is identical; its
normalized weights average intensity and never integrate detector time. Each
front end writes a normalized-pupil-coordinate, cell-integrated
photon-arrival-rate four-pupil mosaic or a typed spectral/path-local bundle.
Generic detector acquisition applies response, QE, and duration afterward.
Their differential
estimators separate persistent calibration/support state, replaceable
single-writer scratch, and caller-visible slopes. Their family convenience
acquisitions likewise separate immutable binning, derived native detector
sampling workspace, and the caller-visible frame. Propagation has an immutable
physical and numerical plan, replaceable backend-bound FFT workspace, and an
exact prepared owner. Prepared optics and estimation owners bind these parts by
identity and reject workspace, state, product, or calibration replacement
before mutating a destination. Recreating equivalent independent workspaces
therefore preserves the numerical result without allowing a stale prepared
owner to follow replacement storage. The estimator state owns valid support,
normalization, reference subtraction,
optical gain, and a calibration revision that invalidates stale prepared plans.
Geometric Pyramid and Bi-O-edge declare `DirectMeasurementPath()` and construct
neither propagation nor acquisition workspace.

Zernike keeps its physical focal-plane phase spot under `Optics`. Its composed
front end binds that component to an immutable propagation plan and a
replaceable backend-bound FFT workspace. Convenience acquisition separately
owns an immutable sampling plan and caller-visible frame. It requires no
placeholder workspace because sampling writes directly into that product. The
normalized-pupil estimator separates immutable parameters,
persistent support/reference calibration state, replaceable reduction scratch,
and its caller-visible signal product. Exact prepared owners reject replacement
of any bound role or calibration revision before destination mutation.
Curvature exposes a fixed positive-/negative-defocus rate tuple that can feed
two independent detectors or an explicitly packed single-detector mapping; its
remaining propagation, acquisition, and estimator owner split is the next
family gate. No synthetic state view or monolithic optical-owner adapter is
retained.
LiFT separately owns a shareable run-immutable `LiFTForwardPlan`, replaceable
per-owner `LiFTForwardWorkspace`, exact prepared forward input and output,
caller-provided `LiFTObservation`, run-immutable `LiFTEstimationPlan`,
replaceable iteration workspace, caller-owned coefficient result, and immutable
diagnostic snapshots. `LiFT` is only the cold estimator definition. Each
`PreparedLiFTEstimator` receives independent propagation and linear-algebra
scratch even when several estimators share the same forward plan. Modal subset
and observation contract are cold-path bindings; repeated estimation neither
owns nor triggers detector acquisition. LiFT reconstruction is the qualified-only
`WavefrontSensors.reconstruct!`/`WavefrontSensors.reconstruct` generic, not
`Control.reconstruct!`/`Control.reconstruct`.

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
Bi-O-edge, detectors, explicit model loops, and Plant implementations.

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

- cold definitions, run-immutable plans, persistent state, replaceable
  workspaces, caller-visible products, and exact prepared bindings are separate
- every mutable endpoint, optic, path workspace, acquisition, and RNG owner has
  one writer
- exported products are distinct from scratch buffers
- boundedness is enforced by prepared-owner lifecycle and capacity checks;
  detector read-offset schedules use concrete fixed-cardinality
  `FixedSizeVector` storage, while remaining direct `Memory` migration outside
  detector acquisition continues under issue #225
- implementation fields and collection element types do not store `Any`
- hot registries do not store abstract element types or uninstantiated
  parametric families; bounded heterogeneous work crosses concrete family
  groups, a concrete prepared owner, or an explicit function barrier
- command endpoints own independent schema, sequence, effective-time, silence,
  bounded calendar, and effective-command state
- controller-output routing borrows exact named products; successful admission
  is the boundary that copies a payload into endpoint-owned storage
- WFS, science, calibration, and coronagraph paths own distinct
  `PupilFunction` or field products; path reuse requires an explicit prepared
  compatibility contract
- WFS and detector pipelines separate persistent state, replaceable scratch,
  and caller-visible sampled/readout products; direct writes into a declared
  product remain valid and do not make that product workspace
- the WFS optics stage produces photon-arrival rates or explicitly dimensionless
  products; the exact prepared detector owner validates its bound input and
  runtime owners, applies presampling response before physical-pixel
  integration, and integrates its explicit exposure exactly once
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

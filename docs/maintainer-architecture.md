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

## Main Data Model

The dominant pattern is:

- immutable params struct
- mutable state struct
- top-level domain object holding params plus state

Examples:

- `Telescope` with immutable `TelescopeParams`, prepared `TelescopeAperture`,
  and transitional `LegacyTelescopePathState`; telescope parameters contain
  spatial geometry but no cadence or exposure duration
- caller-owned `PupilFunction`, `ElectricField`, and `IntensityMap`
  products with immutable `OpticalPlaneMetadata`, including normalization,
  spatial-measure, and coherent/incoherent combination policy
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
- `ClosedLoopRuntime` with runtime profile, output plan, and prepared state

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
- `propagate!(atm, tel, src)` as a transitional convenience adapter
- `measure!(wfs, tel, src, det)`
- `apply!(dm, tel, DMAdditive())`

These functions should own local physics and algorithm behavior, but not broad
backend policy.

Timed atmosphere models have one evolution writer. Evolution publishes an
immutable `AtmosphereEpoch`; prepared, path-local renderers consume only the
current epoch and write caller-owned products. Atmosphere state therefore owns
physical layers and timeline state, not a shared pupil-sized render target or a
mutable last-source geometry cache. The timed atmosphere API reads no telescope
timing value. `ClosedLoopRuntime` and the single/grouped scenario
configurations instead require an explicit positive `atmosphere_step` for each
sensing update; that duration is independent of detector exposure.

### 2. Shared subsystem services

Examples:

- detector pipeline helpers
- grouped WFS execution helpers
- prepared compatible-intensity accumulation and typed detector acquisition
- runtime output planning
- propagation contexts
- backend reductions and random/noise services

These shared services reduce duplicated orchestration logic across SH, Pyramid,
BioEdge, detectors, and runtime.

### 3. Simulation/runtime orchestration

Examples:

- `ClosedLoopRuntime`
- `SimulationInterface`
- `prepare!`
- `runtime_profile`
- `simulation_readout`

This layer coordinates prepared objects, latency staging, exported outputs,
and detector/wfs/science readout ownership.

## Runtime Ownership Model

The current runtime model is:

- product requirements are explicit
- exported outputs are distinct from scratch buffers
- prepared runtime state is separated from per-step mutation
- runtime profiles model delays/output fidelity independently of execution
  residency
- CPU HIL is a direct host-resident plan, while accelerator runtimes use a
  device-resident plan with explicit observation barriers
- shared multi-arm runtimes own one atmosphere advance and command state while
  source-specific arms own their WFS consumers, prepared direct-imaging
  products/workspace, and independent detector acquisitions
- WFS and detector pipelines own their sampled/readout/intermediate products
  explicitly rather than relying on in-place aliasing
- optical formation produces photon-arrival rates or explicitly dimensionless
  products; prepared detector acquisition validates the metadata, applies
  presampling response before physical-pixel integration, and integrates its
  explicit exposure exactly once

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

# Maintainer Architecture Guide

Status: active

Plan traceability:

- [`PLAN-34`](./package-review-action-plan.md)
- review ID: `PR-28`

## Purpose

This document describes the current implemented architecture of
AdaptiveOpticsSim.jl.

It is the maintainer-facing synthesis layer that sits between:

- the platform synthesis guide in
  [`platform-architecture.md`](./platform-architecture.md)
- the public workflow docs in [`user-guide.md`](./user-guide.md)
- the symbol inventory in [`api-reference.md`](./api-reference.md)
- the subsystem plans and audits under [`docs/`](./)

## Architectural Principles

The implemented system follows these package-level rules:

- use multiple dispatch and traits instead of OO-style inheritance
- separate params from mutable runtime state
- make hot paths explicit with mutating `!` functions
- keep core algorithms backend-generic where practical
- isolate backend policy away from model code when possible
- keep deterministic validation and frozen-reference evidence first-class
- keep optional integrations outside core or in extension modules

## Package Shape

The codebase is organized into subsystem directories:

- `src/core`
  - traits, errors, backend services, reduction/random helpers, shared low-level
    infrastructure
- `src/optics`
  - telescope, sources, electric fields, propagation, masks, DMs, PSFs
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
    gain sensing, and misregistration sensitivity
- `src/control`
  - reconstructors, runtime construction, execution, output planning, timing,
    controllers
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

- `Telescope` with `TelescopeParams` and `TelescopeState`
- `ShackHartmannWFS` with `ShackHartmannWFSParams` and `ShackHartmannWFSState`
- `Detector` with `DetectorParams` and `DetectorState`
- `ClosedLoopRuntime` with runtime profile, output plan, and prepared state

This gives:

- explicit configuration
- stable memory ownership
- hot-path mutation without repeated allocation

## Execution Layers

There are three main execution layers:

### 1. Model-local computation

Examples:

- `advance!(atm, tel)`
- `propagate!(atm, tel, src)`
- `measure!(wfs, tel, src, det)`
- `apply!(dm, tel, DMAdditive())`

These functions should own local physics and algorithm behavior, but not broad
backend policy.

### 2. Shared subsystem services

Examples:

- detector pipeline helpers
- grouped WFS execution helpers
- runtime output planning
- propagation contexts
- backend reductions and random/noise services

These were expanded during the package-review cleanup to reduce duplicated
orchestration logic across SH, Pyramid, BioEdge, detectors, and runtime.

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
- WFS and detector pipelines own their sampled/readout/intermediate products
  explicitly rather than relying on in-place aliasing

For the step-by-step view, see [`runtime-dataflow.md`](./runtime-dataflow.md).

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

- [`model-validity-matrix.md`](./model-validity-matrix.md)
- [`benchmark-matrix-plan.md`](./benchmark-matrix-plan.md)
- [`cross-package-benchmark-harness.md`](./cross-package-benchmark-harness.md)

## Validation Structure

The package uses four distinct evidence classes:

- analytic and structural checks
- frozen reference-bundle regression against OOPAO and targeted SPECULA cases
- backend smoke and parity checks
- benchmark evidence, including cross-package benchmark archives

The maintained synthesis doc for this is
[`model-validity-matrix.md`](./model-validity-matrix.md).

## Optional Boundaries

The package is intentionally conservative about what belongs in core.

Core owns:

- optical model primitives
- atmosphere and sensing models
- detector/readout physics
- calibration and runtime orchestration

Optional or boundary-managed surfaces include:

- external science-path integrations
- heavyweight format or ecosystem adapters
- extra plotting or notebook-facing helpers

See [`optional-integration-boundaries.md`](./optional-integration-boundaries.md).

## Recommended Reading Order for Maintainers

1. [`documentation-map.md`](./documentation-map.md)
2. [`platform-architecture.md`](./platform-architecture.md)
3. [`maintainer-architecture.md`](./maintainer-architecture.md)
4. [`runtime-dataflow.md`](./runtime-dataflow.md)
5. [`model-validity-matrix.md`](./model-validity-matrix.md)
6. subsystem plan docs only as needed

# Platform Architecture Guide

Status: active

Plan traceability:

- [`PSP-02`](./platform-strengthening-plan.md)
- direction IDs: `PSR-01`, `PSR-03`, `PSR-06`, `PSR-07`

## Purpose

This guide explains AdaptiveOpticsSim.jl as a simulation platform rather than
as a collection of independent subsystem docs.

Use it when you already understand the core domain concepts and need to answer
questions like:

- how the main subsystem families fit together
- what the runtime actually owns
- where grouped or platform-scale composition belongs
- how validation and benchmark evidence support package-level claims
- which surfaces belong in core versus optional boundaries

This is a stable synthesis guide. It complements:

- [`user-guide.md`](./user-guide.md) for workflow-oriented entry
- [`platform-workflows.md`](./platform-workflows.md) for script-first workflow
  composition
- [`api-reference.md`](./api-reference.md) for maintained symbols
- [`maintainer-architecture.md`](./maintainer-architecture.md) for subsystem
  ownership
- [`runtime-dataflow.md`](./runtime-dataflow.md) for step-by-step runtime flow

## Platform Model

The package should be understood as five connected layers:

1. physical model primitives
2. sensing and detector pipelines
3. calibration and reconstruction surfaces
4. runtime/orchestration surfaces
5. validation and benchmark evidence

Those layers are intentionally Julia-native:

- typed constructors and params/state objects
- explicit mutating hot-path methods
- dispatch and traits for algorithm and backend choice
- scripts as the primary composition interface

The package is not intended to be config-file-first, and it is not intended to
hide the model behind opaque process graphs.

## Major Subsystem Families

### Optics and atmosphere

These subsystems define the plant:

- telescope, masks, sources, asterisms, DMs, NCPA
- atmospheric phase screens and multilayer atmospheres
- electric-field and propagation surfaces

This layer answers:

- what optical state exists
- how phase and field evolve
- what reaches each sensing or science surface

Key docs:

- [`runtime-dataflow.md`](./runtime-dataflow.md)
- [`units-policy.md`](./units-policy.md)
- [`deterministic-simulation.md`](./deterministic-simulation.md)

### Wavefront sensing and detectors

This is the maintained sensor/readout layer:

- Shack-Hartmann, Pyramid, BioEdge, Curvature, Zernike WFS
- frame and counting detectors
- grouped and detector-coupled sensing paths

This layer is responsible for:

- sampled optical signal generation
- detector/readout application
- slope or branch-signal extraction
- explicit ownership of scratch versus exported outputs

### Calibration and reconstruction

This layer turns prepared sensing models into reusable operators:

- interaction matrices
- modal bases
- control matrices
- reconstructors
- LiFT and gain-sensing workflows
- tomography builders and command assembly

This layer is where most “builder-side” heavy algebra belongs. It should stay
typed and explicit rather than being hidden inside runtime loops.

### Runtime and orchestration

This is the platform coordination layer:

- `ClosedLoopRuntime`
- `SimulationInterface`
- `CompositeSimulationInterface`
- runtime output planning
- grouped export policy
- execution policy and prepared runtime state

This is the layer that makes the package a platform rather than only a library
of models.

It is responsible for:

- deciding what products are needed
- staging prepared state outside hot loops
- coordinating detector/wfs/science ownership
- exposing reusable runtime entry points for simulation, validation, and
  benchmarking

### Evidence and validation

The package treats model evidence as a first-class platform layer, not as an
afterthought.

Important evidence classes are:

- analytic and structural checks
- frozen OOPAO and SPECULA reference bundles
- backend smoke and parity checks
- benchmark archives and cross-package contracts

Key docs:

- [`model-validity-matrix.md`](./model-validity-matrix.md)
- [`benchmark-matrix-plan.md`](./benchmark-matrix-plan.md)
- [`cross-package-benchmark-harness.md`](./cross-package-benchmark-harness.md)

## Ownership Model

The dominant architectural rule is:

- immutable params define configuration
- mutable state owns evolving runtime buffers
- top-level domain objects hold params plus state

This pattern supports:

- explicit configuration surfaces
- predictable memory ownership
- preallocation of hot-path workspaces
- backend-generic composition with typed specialization

The same rule carries upward into runtime/orchestration:

- prepared runtime state should be distinct from per-step mutation
- exported outputs should be distinct from scratch buffers
- grouped runtime outputs should be explicit rather than implied by shared
  temporary stacks

## How Platform Composition Works Today

The maintained composition model is script-first.

Typical platform assembly is:

1. construct telescope and source geometry
2. build atmosphere and long-lived screen state
3. build one or more WFS models and optional detectors
4. build DMs plus calibration/reconstruction objects
5. create runtime and product requirements
6. `prepare!` once
7. step, benchmark, or validate through the runtime surface

This composition style is intentional:

- it keeps scenario structure visible in Julia code
- it works naturally with dispatch and typed builders
- it makes benchmark and validation scripts use the same core surfaces as
  normal platform assembly

The next platform-strengthening phases may add richer typed orchestration
objects, but they should still preserve the script-first model.

## Backend Strategy

Backend support is part of the platform architecture, not an optional afterthought.

The intended rule is:

- semantic contracts stay shared
- backend-sensitive execution policy is selected by typed dispatch
- hot-loop selection should be static and type-driven, not symbol- or
  string-driven

This architecture now appears in:

- grouped runtime execution plans
- detector execution plans
- atmospheric-field execution plans
- reduction/export execution plans

The goal is not identical implementation across CPU, CUDA, AMDGPU, and future
Metal paths. The goal is semantic equivalence with backend-appropriate memory
and execution plans.

## Validation And Benchmark Support

Package-level claims should be read together with their evidence.

Use:

- [`model-validity-matrix.md`](./model-validity-matrix.md)
  - to understand whether a model claim is strong, medium-strong, or still
    limited
- [`benchmark-matrix-plan.md`](./benchmark-matrix-plan.md)
  - to understand maintained performance surfaces
- [`cross-package-benchmark-harness.md`](./cross-package-benchmark-harness.md)
  - to understand normalized cross-package comparisons

This means the platform is not only “what can be constructed,” but also “what
is defensibly validated and benchmarked.”

## Core Boundary And Extension Rule

Core should own:

- optics and atmosphere primitives
- sensing and detector physics
- calibration and runtime/orchestration surfaces
- validation and benchmark infrastructure required to support those claims

Core should avoid taking on:

- science-path or focal-plane post-processing stacks
- heavyweight ecosystem adapters as first-class runtime requirements
- config-file-first orchestration as the main interface

Those surfaces should remain optional unless a later explicit boundary review
changes the rule.

See:

- [`optional-integration-boundaries.md`](./optional-integration-boundaries.md)
- [`future-platform-direction.md`](./future-platform-direction.md)

## Reading Path

If you are:

- learning the platform shape:
  - read this guide, then [`platform-workflows.md`](./platform-workflows.md)
- maintaining code structure:
  - read this guide, then [`maintainer-architecture.md`](./maintainer-architecture.md)
- debugging runtime ownership or export behavior:
  - read this guide, then [`runtime-dataflow.md`](./runtime-dataflow.md)
- evaluating whether a claim is supported:
  - read this guide, then [`model-validity-matrix.md`](./model-validity-matrix.md)

Use subsystem plans only after the stable guides above are no longer enough.

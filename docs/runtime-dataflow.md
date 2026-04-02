# Runtime Dataflow Guide

Status: active

Plan traceability:

- [`PLAN-36`](./package-review-action-plan.md)
- review ID: `PR-30`

## Purpose

This document explains the end-to-end runtime and dataflow of the maintained
simulation stack.

It is not a symbol reference. It is the operational picture of:

- what gets built
- what gets prepared
- what mutates each step
- which products are exported
- where validation and benchmark evidence attach

## Main Runtime Objects

The most important runtime-facing objects are:

- `Telescope`
- atmosphere model such as `MultiLayerAtmosphere` or
  `InfiniteMultiLayerAtmosphere`
- one or more WFS objects
- optional `Detector` objects
- `DeformableMirror`
- calibration/reconstructor objects
- `ClosedLoopRuntime`

The orchestration layer is in `src/Control`.

## Build Phase

The build phase constructs long-lived simulation objects:

1. create telescope and source geometry
2. create atmosphere model and persistent screen state
3. create WFS models and detectors
4. create DM and calibration/reconstruction objects
5. build `ClosedLoopRuntime` with:
   - runtime profile
   - latency model
   - product requirements
   - execution policy

Outputs of this phase should be stable objects with owned buffers, not
ad hoc temporary arrays.

## Prepare Phase

`prepare!` and subsystem-specific preparation do the work that should not
happen inside the hot loop:

- allocate exported readout and metadata surfaces
- precompute calibration products
- prepare grouped execution stacks
- prepare delay lines and runtime staging
- ensure detector and WFS pipeline buffers exist

The design goal is:

- construction is allowed to be heavier
- the step path should only mutate prepared state

## Per-Step Flow

At a high level, one closed-loop step is:

1. advance atmosphere
2. reset or update telescope OPD state
3. apply DM command into telescope phase
4. propagate source/field/atmosphere state to the sensing surface
5. run WFS measurement and optional detector readout
6. extract required products:
   - slopes
   - WFS pixels
   - science pixels
7. run reconstruction/controller update
8. stage delayed commands or delayed readout products as needed

This split is visible in runtime timing and benchmark surfaces.

## Product Ownership

The package now makes product ownership more explicit than before.

Important distinctions:

- scratch buffers are not exported products
- sampled detector inputs are not the same thing as detector outputs
- grouped WFS intermediate stacks are not the same thing as archived or
  exported readouts

For example, in diffractive Shack-Hartmann execution there is now a clearer
boundary between:

- sampled pre-detector spot stacks
- post-detector signal stacks
- exported runtime pixel products

That separation matters for:

- correctness
- GPU layout ownership
- benchmarked allocation behavior

## Detector and WFS Dataflow

The maintained detector/WFS pipeline now follows a more explicit shape:

1. produce sampled optical signal
2. apply detector/readout pipeline if configured
3. reduce or extract slopes/signals
4. snapshot only the runtime products that were requested

The runtime product plan decides whether a given simulation step must produce:

- slopes only
- slopes plus WFS pixels
- science pixels
- metadata surfaces

For grouped composite execution, exported products can additionally include:

- per-branch grouped WFS frames
- per-branch grouped science frames
- compatible-shape grouped WFS stacks
- compatible-shape grouped science stacks

This keeps slopes-only runs from paying unnecessary export costs.

## Backend Flow

The intended backend execution model is:

- model code calls shared backend services
- shared services select CPU / CUDA / AMDGPU behavior
- optional backend extensions refine behavior without rewriting the public API

Backend validation is attached at three levels:

- optional backend test coverage
- smoke scripts
- benchmark evidence

Representative performance evidence is intentionally kept out of `Pkg.test()`.

## Benchmarks and Evidence

There are two relevant evidence layers:

### Functional/validation evidence

- [`model-validity-matrix.md`](./model-validity-matrix.md)
- frozen OOPAO and SPECULA reference bundles

### Runtime/engineering evidence

- [`benchmark-matrix-plan.md`](./benchmark-matrix-plan.md)
- [`cross-package-benchmark-harness.md`](./cross-package-benchmark-harness.md)
- archived result files under `benchmarks/results/`

## How to Use This Guide

When debugging or extending runtime behavior:

1. use this doc to identify the layer that owns the behavior
2. use [`maintainer-architecture.md`](./maintainer-architecture.md) for the
   subsystem map
3. use [`api-reference.md`](./api-reference.md) for specific symbols
4. use subsystem plans only after the stable guides above are no longer enough

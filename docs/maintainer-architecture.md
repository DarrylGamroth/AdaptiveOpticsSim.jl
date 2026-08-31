# Maintainer Architecture Guide

Status: active

## Purpose

This guide describes the maintained package architecture. The canonical
terminology is in [`glossary.md`](glossary.md), graph execution is detailed in
[`runtime-dataflow.md`](runtime-dataflow.md), and extension contracts are in
[`extension-guide.md`](extension-guide.md).

## Design Rules

- Use multiple dispatch and traits instead of inheritance-shaped object graphs.
- Keep immutable params separate from mutable state.
- Use `Plan` only for run-immutable numerical or physical contracts.
- Keep persistent scientific state, replaceable workspace, caller-visible
  products, and exact input/output/context bindings in separate owners.
- Give every mutable scientific state one writer.
- Prepare allocation and validation work before repeated `!` operations.
- Parameterize numerical owners by element type, array backend, and exact
  compute device where required.
- Do not store `Any`, abstract element types, or uninstantiated parametric
  families in execution registries.
- Use `FixedSizeArray` with a concrete element type for sealed homogeneous host
  registries. A `Vector` may be a cold builder but is not the armed registry.
- Do not name Julia `Memory` in package or extension implementation.
- Keep file formats and optional ecosystem integrations outside core.
- Fail with structured exceptions; do not print and return.
- Preserve the vocabulary in [`glossary.md`](glossary.md).

## Package Shape

The root module loads independent domain owners and exports routine cross-domain
vocabulary plus the modules themselves:

| Module | Responsibility |
|---|---|
| `Backends` | array allocation, backend/device identity, reductions, random helpers, kernels, synchronization |
| `Optics` | telescope/source geometry, optical products and planes, propagation, direct imaging, NCPA, deformable mirrors and other physical optics |
| `Atmospheres` | turbulence definitions/state, evolution, source-direction rendering and batches |
| `Detectors` | detector response, sensor families, acquisition, readout, products |
| `WavefrontSensors` | composed WFS optics, observations, measurements, estimators |
| `Calibration` | interaction/control matrices, bases, inverse policies, fitting and identification |
| `Control` | reconstruction, controller state, delay lines, closed-loop operations |
| `Tomography` | guide-star geometry, atmospheric reconstruction, fitting, DM projection |
| `Ensembles` | coarse independent runs and optional offline scheduling |
| `AlgorithmGraphs` | static complete-frame graph composition and lockstep HIL exchange |

Dense domain APIs belong to their owner module. Root exports are deliberately
curated. Stable advanced seams are marked `public` in the owning module rather
than forwarded through compatibility aliases.

## Ownership Model

A prepared numerical path separates six roles:

1. Immutable configuration describes user intent.
2. A run-immutable plan owns validated coefficients, mappings, and numerical
   contracts.
3. Mutable state owns values that affect later results.
4. A replaceable workspace owns scratch and execution resources.
5. Products expose caller-visible scientific results.
6. A prepared execution owner binds the exact concrete values, arrays, backend,
   device/context, and single writer.

Discarding state may change the next result. Recreating workspace must not.
Exact external array identity usually belongs to the prepared owner, not the
reusable plan.

Persistent RNGs are state. FFT plans or backend handles belong in the plan only
when they are logically immutable and reentrant; scratch-owning, stream-bound,
or non-reentrant handles belong in workspace or the prepared owner.

## Scientific Execution

Domain operations keep specific verbs such as `render_atmosphere!`,
`form_wfs_optical_products!`, `acquire_wfs_observation!`,
`estimate_wfs_measurement!`, `reconstruct!`, and `update_surface!`. There is no
universal package-wide `process!` interface.

Direct Julia composition is the unrestricted modeling surface. It owns explicit
ordering when a scenario has generated topology, multiple cadences, conditional
execution, or sub-frame optical sampling.

## Algorithm Graphs

`AlgorithmGraphs` is a small, deterministic complete-frame scheduler around
native scientific operations. A definition contains nodes, graph inputs and
outputs, direct links, and explicit one-frame delayed links. Preparation:

- validates node labels and port contracts
- rejects missing, duplicate, incompatible, or cyclic direct connections
- computes one fixed topological order
- allocates and binds exact storage
- prepares each concrete node owner
- fixes the backend/device execution context

`step_graph!` executes that order with one writer and publishes graph outputs
only for a completed step. `reset_graph!` restores node state, delayed links,
RNG seeds, and sequence state according to each node contract.

TOML graph files are a versioned static authoring surface. Julia definitions
remain available for generated or application-specific composition. The graph
runtime owns no wall-clock pacer, transport, dynamic topology, or automatic
placement.

`PreparedGraphHILBoundary` adds one lockstep host exchange: execute a complete
frame, publish it to a host buffer, require a matching complete command, then
allow the next frame. PipeWireAO or another application owns transport and
pacing.

## Backend Strategy

Core algorithms accept generic array storage only when their implementation is
actually generic. GPU-ready paths avoid scalar indexing, use AbstractFFTs and
KernelAbstractions where appropriate, and keep transfers explicit.

A prepared graph has one exact target context. It does not silently copy between
CPU and GPU or fall back to CPU. HIL exchange is an explicit host boundary, so a
GPU simulation may copy one completed detector frame to a CPU RTC and copy the
next complete command back.

Parallelism is coarse grained across independent sources, time steps, or
sweeps. Avoid nested Julia, BLAS, FFT, backend, and graph parallelism.
Deterministic mode uses one thread.

## Optional Integrations

Optional backends and packages live in extensions or companion packages.
Proper.jl prescriptions are application-owned and may be called directly or
adapted as explicit graph nodes. PipeWire embedding is outside this package;
AdaptiveOpticsSim supplies complete arrays and model timestamps at a transport-
neutral boundary.

## Validation

Tests are registered in `test/test_selection.jl`. Select the smallest relevant
suite during development:

~~~bash
julia --project=. test/ci/run_selected_tests.jl algorithm-graphs
julia --project=. test/ci/run_selected_tests.jl core
julia --project=. test/ci/run_selected_tests.jl sensors
~~~

Use bare `Pkg.test()` for a cross-cutting or release gate. Validate optional
AMDGPU locally, CUDA on the WSL host, and Linux aarch64 on the Raspberry Pi only
when the changed surface requires those targets. Record exact commands and
versions with release evidence.

## Maintainer Reading Order

1. [`glossary.md`](glossary.md)
2. [`runtime-dataflow.md`](runtime-dataflow.md)
3. [`extension-guide.md`](extension-guide.md)
4. [`supported-production-surfaces.md`](supported-production-surfaces.md)
5. [`release-validation-runbook.md`](release-validation-runbook.md)

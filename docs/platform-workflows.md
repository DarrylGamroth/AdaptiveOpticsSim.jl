# Platform Workflow Guide

Status: active

Plan traceability:

- [`PSP-03`](./platform-strengthening-plan.md)
- direction IDs: `PSR-01`, `PSR-02`, `PSR-04`

## Purpose

This guide explains how to use AdaptiveOpticsSim.jl in practice once you move
beyond isolated tutorials.

It focuses on the maintained script-first workflow for:

- constructing a realistic simulation
- choosing the right runtime surface
- validating model behavior
- running maintained benchmarks
- deciding when grouped or platform-scale composition is appropriate

This is not a symbol reference and not a plan doc. It is the stable workflow
bridge between:

- [`user-guide.md`](./user-guide.md)
- [`platform-architecture.md`](./platform-architecture.md)
- [`runtime-dataflow.md`](./runtime-dataflow.md)
- [`model-validity-matrix.md`](./model-validity-matrix.md)
- [`benchmark-matrix-plan.md`](./benchmark-matrix-plan.md)
- [`scenario-builder-style.md`](./scenario-builder-style.md)

## Workflow Model

The maintained workflow model is:

1. write a Julia script
2. construct typed model objects directly
3. prepare long-lived runtime state outside hot loops
4. run simulation, validation, or benchmark entry points intentionally
5. record outputs or evidence using the maintained runtime/export surfaces

The package should not be approached as:

- a config-file-first application
- a hidden process graph
- or a benchmark-only harness

The same typed objects should be reusable across:

- ordinary simulations
- validation scripts
- benchmark scripts
- grouped/runtime composition work

## Choose The Right Workflow

### Workflow 1: Normal platform simulation

Use this when you want to build or study a realistic AO system in Julia.

Typical shape:

1. construct telescope and source geometry
2. add atmosphere
3. build WFS and optional detectors
4. build DM and calibration/reconstruction surfaces
5. create runtime or higher-level simulation wrapper
6. `prepare!` if the runtime supports preparation
7. step the simulation and inspect exported products

Use:

- exported workflow APIs from [`api-reference.md`](./api-reference.md)
- tutorial/example support code where it helps
- script-local helper functions for scenario setup

Do not start with:

- benchmark contracts
- frozen reference generators
- cross-package harnesses

### Workflow 2: Validation and reference comparison

Use this when your goal is to support or check a scientific or engineering
claim.

Typical shape:

1. construct a deterministic scenario
2. disable or control noise where appropriate
3. run the maintained reference or analytic comparison surface
4. compare outputs within a documented tolerance
5. record the scope and limitation of the claim

Use:

- [`model-validity-matrix.md`](./model-validity-matrix.md)
- OOPAO and SPECULA frozen bundles
- committed artifact generators when a maintained artifact exists

Typical entry points:

- `Pkg.test()`
- reference harness scripts/generators
- deterministic comparison scripts under `scripts/`

This workflow is for defending claims, not only for “does the code run?”

### Workflow 3: Benchmark and engineering evidence

Use this when your goal is runtime, memory, backend, or scaling evidence.

Typical shape:

1. choose a maintained benchmark family
2. run the benchmark intentionally outside `Pkg.test()`
3. record timing, allocations, and scenario shape
4. compare against archived evidence where appropriate

Use:

- [`benchmark-matrix-plan.md`](./benchmark-matrix-plan.md)
- [`cross-package-benchmark-harness.md`](./cross-package-benchmark-harness.md)
- maintained profile scripts under `scripts/`

Do not treat benchmark scripts as a substitute for normal workflow examples.

### Workflow 4: Grouped or platform-scale composition

Use this when you need:

- multi-source sensing
- grouped WFS products
- composite runtime interfaces
- multi-branch exported readouts

Start with:

- [`grouped-runtime-contract.md`](./grouped-runtime-contract.md)
- [`grouped-runtime-validation.md`](./grouped-runtime-validation.md)
- the grouped runtime surfaces in the runtime API

This is still script-first. The grouped/platform-scale path should be built
from typed Julia objects and explicit runtime products, not from ad hoc hidden
configuration.

## Canonical Script Shape

The recommended script shape is:

```julia
using AdaptiveOpticsSim
using Random

function build_scenario(; T=Float32, backend=Array)
    tel = Telescope(...)
    src = Source(...)
    atm = MultiLayerAtmosphere(tel; ..., T=T, backend=backend)
    wfs = ShackHartmann(tel; ..., T=T, backend=backend)
    det = Detector(...; T=T, backend=backend)
    dm = DeformableMirror(tel; ..., T=T, backend=backend)
    imat = interaction_matrix(dm, wfs, tel, src; ...)
    recon = ModalReconstructor(imat; ...)
    runtime = ClosedLoopRuntime(...)
    return (; tel, src, atm, wfs, det, dm, recon, runtime)
end

function main()
    rng = MersenneTwister(0)
    sim = build_scenario()
    prepare!(sim.runtime)
    for _ in 1:10
        step!(sim.runtime)
    end
    return simulation_readout(sim.runtime)
end
```

The important structural rules are:

- separate scenario construction from execution
- keep deterministic RNG setup visible
- keep preparation explicit
- use exported runtime/readout accessors rather than digging through internal
  buffers

For the maintained style rules behind that script shape, use
[`scenario-builder-style.md`](./scenario-builder-style.md).

## Building A Reproducible Platform Script

Use these rules:

- choose `T` and `backend` at construction time
- keep long-lived scenario parameters in one visible place
- use deterministic seeds in validation and benchmark scripts
- prefer typed helper functions over dynamically assembled dictionaries or
  loose config blobs
- keep any advanced helper usage namespaced as `AdaptiveOpticsSim.<name>`
- follow the maintained script-construction rules in
  [`scenario-builder-style.md`](./scenario-builder-style.md)

Good script responsibilities:

- define the physical scenario
- define runtime/export needs
- invoke maintained simulation or benchmark surfaces

Bad script responsibilities:

- reimplement runtime ownership rules
- stash hidden mutable globals
- mix benchmark contracts into ordinary simulation code

## When To Use Tutorials, Stable APIs, Or Support Scripts

### Tutorials

Use tutorials when:

- you are learning one subsystem family
- you want a runnable minimal example
- you want a starting point for a small script

### Stable APIs

Use stable APIs when:

- you are writing the actual simulation or orchestration script
- you need a maintained public surface
- you want code that should survive refactors cleanly

### Example support modules

Use support modules under `examples/support/` when:

- you want a maintained, realistic scenario skeleton
- you are building AO188/AO3k-style scenarios
- you need a more realistic starting point than a single tutorial

### Benchmark and contract scripts

Use benchmark or contract scripts when:

- you are collecting evidence
- you are re-baselining performance
- you are comparing against OOPAO, SPECULA, or REVOLT-like scenarios

These are engineering surfaces, not the default user path.

## Validation Workflow In Practice

Use this sequence:

1. decide what claim you are making
2. find the model family in [`model-validity-matrix.md`](./model-validity-matrix.md)
3. pick the matching evidence surface:
   - analytic test
   - frozen reference bundle
   - backend parity/smoke
   - maintained benchmark artifact
4. keep deterministic inputs where the evidence requires them
5. record limitations instead of over-claiming equivalence

Examples:

- use frozen OOPAO bundles for PSF/SH/Pyramid/BioEdge parity questions
- use SPECULA-targeted bundles for atmospheric-field, Curvature, and Zernike
  contract questions
- use grouped runtime artifacts for grouped export/runtime ownership questions

## Benchmark Workflow In Practice

Use this sequence:

1. choose whether the question is regression-sized or representative
2. choose the maintained script from
   [`benchmark-matrix-plan.md`](./benchmark-matrix-plan.md)
3. run it outside `Pkg.test()`
4. record:
   - runtime
   - allocations
   - backend
   - scenario scale
5. compare to archived evidence or recent maintained numbers

Examples of maintained local benchmark surfaces:

- `scripts/profile_ao3k_runtime.jl`
- `scripts/profile_multi_source_multi_wfs_runtime.jl`
- `scripts/profile_atmospheric_field_runtime.jl`
- `scripts/profile_revolt_pwfs_runtime.jl`

Cross-package comparison should go through the harness, not through ad hoc
one-off scripts.

## Grouped And Platform-Scale Workflows

Use grouped/platform workflows when:

- multiple sources contribute to the same sensing product
- you need grouped WFS/science export stacks
- you need `CompositeSimulationInterface`
- you are studying multi-branch runtime behavior

The practical rules are:

- request only the runtime products you need
- prefer `RuntimeScenario` plus `SingleRuntimeConfig` or
  `GroupedRuntimeConfig` when a script is assembling reusable runtime
  composition rather than a one-off local demo
- keep grouped export policy explicit
- rely on `simulation_grouped_wfs_stack` and
  `simulation_grouped_science_stack` only when the runtime contract says the
  layout is compatible
- for mixed layouts, use explicit per-branch readout surfaces instead of
  assuming a fake common stack

## How To Move From Small To Large Scenarios

Recommended progression:

1. one source, one WFS, no detector
2. one source, one detector-backed WFS
3. one closed-loop runtime
4. one realistic AO188/AO3k-style support scenario
5. grouped or composite runtime scenario
6. validation or benchmark script on top of the same typed model

This preserves continuity between:

- tutorial learning
- real simulation scripts
- benchmark evidence
- cross-package comparison

## Related Reading

- [`platform-architecture.md`](./platform-architecture.md)
  - package-level platform picture
- [`platform-orchestration.md`](./platform-orchestration.md)
  - typed platform scenario/config layer
- [`runtime-dataflow.md`](./runtime-dataflow.md)
  - step-by-step runtime ownership
- [`api-reference.md`](./api-reference.md)
  - maintained public symbols
- [`model-validity-matrix.md`](./model-validity-matrix.md)
  - current claim/evidence status
- [`benchmark-matrix-plan.md`](./benchmark-matrix-plan.md)
  - maintained benchmark families

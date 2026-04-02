# Scenario Builder Style Guide

Status: active

Plan traceability:

- [`PSP-04`](./platform-strengthening-plan.md)
- direction IDs: `PSR-01`, `PSR-02`

## Purpose

This guide defines the maintained style for script-first scenario construction
in AdaptiveOpticsSim.jl.

Use it when writing:

- normal simulation scripts
- validation scripts
- benchmark scripts
- support modules that build realistic platform scenarios

Its job is to keep scenario authoring:

- Julia-native
- typed
- reproducible
- readable
- and aligned with the package architecture

This guide intentionally does not introduce a manifest or config-file-first
workflow. Scenario manifests remain deferred and secondary.

## Core Rule

The maintained scenario-authoring model is:

- typed Julia objects
- explicit helper functions
- explicit runtime preparation
- explicit validation or benchmark entry points

The package should prefer scripts that are clear Julia programs, not opaque
configuration wrappers.

## Canonical Structure

The preferred script structure is:

1. imports
2. small scenario parameter/config struct or named helper defaults
3. one or more `build_*` functions
4. one `run_*` or `main` function
5. explicit result/readout return

Recommended skeleton:

```julia
using AdaptiveOpticsSim
using Random

Base.@kwdef struct MyScenarioConfig{T<:AbstractFloat}
    resolution::Int = 64
    diameter::T = T(8.0)
    sampling_time::T = T(1e-3)
    magnitude::T = T(8.0)
end

function build_scenario(cfg::MyScenarioConfig{T}; backend=Array) where {T}
    tel = Telescope(
        resolution=cfg.resolution,
        diameter=cfg.diameter,
        sampling_time=cfg.sampling_time,
        T=T,
        backend=backend,
    )
    src = Source(band=:I, magnitude=cfg.magnitude, T=T)
    return (; tel, src)
end

function main(; seed=0)
    rng = MersenneTwister(seed)
    cfg = MyScenarioConfig{Float32}()
    sim = build_scenario(cfg)
    return sim
end
```

The exact scenario contents will vary, but the structure should stay close to
that shape.

## Params And State In Scripts

Follow the package’s own architectural rule:

- keep configuration visible and mostly immutable
- let domain objects own mutable runtime state

In practice:

- do put physical constants, dimensions, gains, and mode choices in an
  explicit config surface
- do not mirror internal mutable state in ad hoc script globals
- do not duplicate package-owned runtime buffers in parallel script-side
  structures

Good:

- `cfg::MyScenarioConfig`
- `tel = Telescope(...)`
- `runtime = ClosedLoopRuntime(...)`

Bad:

- loose dictionaries of mixed types for core scenario definition
- separate script-owned mutable arrays that shadow detector or runtime state

## Deterministic RNG Ownership

Deterministic control should be explicit.

Rules:

- create RNGs in `main` or in a dedicated runner function
- pass RNGs into atmosphere and detector operations where reproducibility
  matters
- keep benchmark and validation seeds visible at the top level
- do not hide seed creation inside deep helper layers unless the helper is
  explicitly a deterministic fixture generator

Recommended pattern:

```julia
function main(; seed=0)
    rng = MersenneTwister(seed)
    sim = build_scenario(...)
    advance!(sim.atm, sim.tel; rng=rng)
    return sim
end
```

For validation and benchmark scripts:

- seed values should be recorded in the script or artifact metadata
- random behavior should be deliberate, not ambient

## Helper Function Boundaries

Use small helper functions with clear ownership.

Preferred helper kinds:

- `build_*`
  - construct scenario objects
- `prepare_*`
  - do scenario-specific preparation that belongs outside the hot loop
- `run_*`
  - execute one intended workflow
- `collect_*` or `summarize_*`
  - turn outputs into a stable result payload

Avoid helper functions that:

- both construct and execute long workflows implicitly
- mutate global state behind the caller’s back
- hide major backend or precision choices

## Backend And Precision Choices

Backend and precision should be explicit at construction time.

Rules:

- choose `T` and `backend` near the top-level script entry
- thread them through builder functions explicitly
- use namespaced advanced helpers where needed, for example:
  - `AdaptiveOpticsSim.CPUBuildBackend()`
  - `AdaptiveOpticsSim.GPUArrayBuildBackend(...)`
- do not bury backend choice inside random utility layers

This keeps:

- CPU/GPU comparisons clearer
- benchmark scripts auditable
- validation scripts reproducible

## Runtime And Product Planning

Scenario scripts should be explicit about what outputs they need.

Rules:

- request runtime products intentionally
- use `prepare!` when the runtime supports a prepared path
- access outputs through maintained readout accessors
- keep grouped export policy explicit in grouped/composite scenarios

Good:

- `simulation_readout(runtime)`
- `simulation_slopes(readout)`
- `simulation_wfs_frame(readout)`
- `simulation_grouped_wfs_stack(readout)` only when the grouped contract says
  the layout is compatible

Bad:

- reaching directly into internal scratch buffers as a normal script pattern
- assuming grouped stack compatibility without checking the grouped runtime
  contract

## Script Categories

### Normal simulation scripts

Purpose:

- build and run one realistic scenario

Should emphasize:

- readability
- clear configuration
- explicit outputs

### Validation scripts

Purpose:

- support a scientific or engineering claim

Must emphasize:

- determinism
- tolerance-aware comparisons
- explicit baseline/evidence references

Validation scripts should link naturally to:

- [`model-validity-matrix.md`](./model-validity-matrix.md)
- frozen bundle or artifact generators where applicable

### Benchmark scripts

Purpose:

- collect engineering evidence

Must emphasize:

- scenario shape
- backend choice
- timing and allocation reporting
- archived or comparable results

Benchmark scripts should stay outside `Pkg.test()` unless they are deliberately
small smoke/regression surfaces.

## What Not To Put In Scripts

Avoid these patterns:

- primary config-file parsing as the main scenario definition mechanism
- hidden mutable singletons
- implicit global backend switching in the middle of scenario execution
- direct manipulation of internals that already have maintained accessors
- benchmark-contract normalization logic mixed into ordinary simulation scripts
- scientific post-processing pipelines that belong behind an optional boundary

If a script keeps growing new responsibilities, prefer:

- factoring shared logic into a support module
- or promoting the pattern into a maintained builder/runtime surface

## Scenario Naming And Decomposition

Use names that reveal intent.

Good examples:

- `build_closed_loop_sh_scenario`
- `run_reference_regression`
- `profile_grouped_runtime_medium`
- `collect_runtime_summary`

Avoid names that hide workflow intent:

- `run_all`
- `helper1`
- `test_case`

Scenario decomposition should track workflow intent:

- construction
- preparation
- execution
- result capture

## Recommended Reading Path

Use these docs together:

- [`platform-architecture.md`](./platform-architecture.md)
  - package-level picture
- [`platform-workflows.md`](./platform-workflows.md)
  - which workflow you are actually writing
- [`runtime-dataflow.md`](./runtime-dataflow.md)
  - runtime ownership and export model
- [`api-reference.md`](./api-reference.md)
  - maintained public surfaces

Use this guide whenever a new scenario script starts to drift toward hidden
configuration, hidden state, or mixed-purpose workflow logic.

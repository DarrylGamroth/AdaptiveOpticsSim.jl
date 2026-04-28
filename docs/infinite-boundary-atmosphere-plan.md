# Infinite Boundary-Injection Atmosphere Plan

This document is the concrete implementation plan for replacing the current
persistent finite periodic moving-screen atmosphere with a higher-fidelity
infinite boundary-injection model.

It refines Milestone 1 in
[algorithmic-implementation-roadmap.md](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/docs/algorithmic-implementation-roadmap.md)
and the `ATM-*` requirements in
[atmosphere-runtime-spec.md](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/docs/atmosphere-runtime-spec.md).

## Purpose

The current [multilayer.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/atmosphere/multilayer.jl)
model is scientifically and operationally better than the older redraw-based
path, but it still uses a finite periodic canvas. That is acceptable for fast
HIL-oriented runs and short sequences, but it is not yet the final frozen-flow
model implied by `ATM-021`.

The next atmosphere objective is:

- preserve the current GPU-capable fast path,
- add an explicit infinite-screen backend for science fidelity,
- keep the API and tests clear about which atmosphere model is being used.

## Top-Level Decision

This work SHOULD introduce a new atmosphere backend rather than mutating the
current `MultiLayerAtmosphere` in place.

Planned model split:

- `MultiLayerAtmosphere`: retained as the finite periodic moving-screen fast path.
- `InfiniteMultiLayerAtmosphere`: new boundary-injection atmosphere for
  long-run frozen-flow fidelity.

Reason:

- the finite model remains useful for HIL and performance studies,
- the infinite model has materially different state and update operators,
- keeping both models explicit avoids hiding a major fidelity/performance tradeoff.

## Scope

In scope:

- per-layer infinite phase-screen state,
- covariance-consistent boundary injection,
- wind-driven transport with subpixel extraction,
- CPU reference implementation,
- GPU-capable state layout and extraction path,
- multilayer composition,
- deterministic replay and regression coverage,
- explicit runtime/profile selection between finite and infinite models.

Out of scope:

- replacing the current fast finite model,
- inter-process transport or hardware interfaces,
- control redesign,
- science propagation through `Proper.jl`,
- atmospheric boiling or non-frozen-flow forcing.

## Traceability Targets

This plan primarily closes the remaining gaps in:

- `ATM-020` persistent per-layer runtime state
- `ATM-021` wind transport and boundary extension instead of redraw
- `ATM-023` subpixel accumulation
- `ATM-024` zero-wind stationarity
- `ATM-025` explicit source-aware footprint extraction
- `ATM-030` explicit RNG ownership
- `ATM-031` deterministic replay

New implementation-specific tracking IDs for this plan:

- `INF-001` Infinite-screen layers MUST not rely on periodic wraparound to
  synthesize newly exposed turbulence.
- `INF-002` Boundary injection MUST be covariance-consistent with the configured
  turbulence model within documented tolerances.
- `INF-003` Long-run atmosphere evolution MUST remain statistically stationary.
- `INF-004` Infinite atmosphere MUST preserve the existing OPD-in-meters runtime contract.
- `INF-005` Infinite atmosphere MUST retain backend-generic state layout and
  keep scalar indexing out of GPU hot paths.

## Target Architecture

### New types

Planned immutable parameter structs:

- `InfinitePhaseScreenParams`
- `InfiniteLayerParams`
- `InfiniteMultiLayerParams`

Planned mutable state structs:

- `InfinitePhaseScreenState`
- `InfiniteLayerState`
- `InfiniteMultiLayerState`

Planned container/runtime structs:

- `InfinitePhaseScreen`
- `InfiniteAtmosphereLayer`
- `InfiniteMultiLayerAtmosphere`

### Per-layer state

Each infinite layer should own:

- current screen buffer
- pupil extraction workspace
- boundary injection stencil definition
- covariance-derived predictor operator
- covariance-derived residual factor
- integer transport accumulator
- subpixel transport accumulator
- wind direction metadata
- altitude metadata

RNG ownership remains external per `ATM-030` and `ATM-031`.

### Reuse from current code

The following should be reused where possible:

- backend-generic FFT/PSD utilities from
  [kolmogorov.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/atmosphere/kolmogorov.jl)
- backend-generic extraction flow from
  [multilayer.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/atmosphere/multilayer.jl)
- shared `K_{5/6}` helper from
  [kv56.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/core/kv56.jl)
- existing `fractional_cn2` semantics and multilayer validation tests

## Algorithm Plan

### Phase A: CPU reference model

Goal:

- build a numerically clear CPU implementation before optimizing or porting.

Algorithm:

1. Initialize a larger-than-pupil phase screen using the existing Von Karman
   generator.
2. Define a one-sided boundary injection stencil for each wind direction class.
3. Precompute the conditional Gaussian boundary model:
   - predictor matrix `A`
   - residual factor `B`
4. Advance transport continuously via accumulated physical wind offsets.
5. When integer transport crosses a pixel boundary:
   - shift the effective screen support,
   - inject one or more new boundary rows or columns,
   - keep the existing interior state.
6. Continue subpixel extraction using the existing interpolation path.

Initial directional simplification:

- first implement axis-aligned injection for the dominant boundary direction,
- then extend to arbitrary wind vectors by decomposing row and column injections.

### Phase B: Boundary injection math

Goal:

- isolate the statistically hard part from the atmosphere container logic.

Deliverables:

- a small internal module or file for boundary-operator construction
- clear separation between:
  - stencil geometry
  - covariance evaluation
  - factorization
  - stochastic boundary sampling

Implementation note:

- do not compute large dense solves inside `advance!`
- all covariance/operator work must be builder-time/precompute work

### Phase C: Infinite multilayer manager

Goal:

- make the infinite model drop into the same simulation stack cleanly.

Deliverables:

- `advance!(::InfiniteMultiLayerAtmosphere, tel, rng)`
- `propagate!(::InfiniteMultiLayerAtmosphere, tel)`
- layer accumulation using the existing `fractional_cn2` semantics
- explicit selection between finite and infinite atmosphere constructors

### Phase D: GPU-capable port

Goal:

- preserve the GPU-ready direction established by the current finite model.

Required rules:

- no scalar GPU indexing in runtime kernels
- no host fallback in the hot path
- factorization/precompute may remain CPU-first initially if runtime state and
  stepping stay on-device afterward

Expected split:

- builder-time covariance/factorization work may start on CPU
- step-time extraction, accumulation, and boundary application must be
  backend-generic

### Phase E: Runtime/profile integration

Goal:

- make the model choice explicit for HIL vs science runs.

Planned interface direction:

- `FiniteMovingAtmosphereProfile` or existing `MultiLayerAtmosphere`
- `InfiniteAtmosphereProfile` or `InfiniteMultiLayerAtmosphere`

Selection must be explicit in constructors or runtime profiles rather than
hidden behind a heuristic.

## Implementation Sequence

### Work package 1: Spec and type skeleton

Files:

- `docs/atmosphere-runtime-spec.md`
- `src/atmosphere/`

Tasks:

- add explicit note that the finite periodic path remains a maintained fast path
- define the public type names and constructor surface for the infinite model
- define which parts of the current extraction path are shared vs model-specific

Acceptance:

- the planned public API exists as skeletons or documented signatures

### Work package 2: Boundary math utilities

Files:

- `src/atmosphere/infinite_screen_math.jl`
- `src/atmosphere/phase_stats.jl`

Tasks:

- implement covariance stencil assembly
- implement conditional predictor/residual construction
- add numerical conditioning checks and structured errors
- document allowed parameter regimes and failure modes

Acceptance:

- CPU tests validate the conditional sampler against analytic covariance

Builder-time constraints and current failure modes:

- `stencil_size` must be odd and larger than the maintained screen size.
- The initial builder path constructs one-sided row/column stencils and
  conditional Gaussian operators on CPU.
- Ill-conditioned stencil covariance or non-PSD residual covariance must fail
  with structured numerical-condition errors instead of silently producing bad
  operators.

### Work package 3: CPU infinite layer stepping

Files:

- `src/atmosphere/infinite_screen.jl`
- `src/atmosphere/multilayer.jl`
- `test/runtests.jl`

Tasks:

- implement per-layer transport and boundary injection
- support row and column insertion paths
- support zero-wind and nonzero-wind evolution
- preserve subpixel accumulation

Acceptance:

- long-run states do not show periodic wrap artifacts
- zero-wind layers remain stationary
- small-step evolution remains shift-correlated

Implementation status:

- `InfiniteMultiLayerAtmosphere` now performs maintained CPU/GPU stepping with
  precomputed row/column boundary operators, integer boundary injection, and
  residual subpixel extraction.
- GPU construction now warms the boundary-injection paths up front so the first
  integer-shift sample does not absorb one-time kernel setup in the runtime
  allocation surface.

### Work package 4: Multilayer and off-axis extraction

Files:

- `src/atmosphere/multilayer.jl`
- `src/optics/source.jl`
- `test/runtests.jl`

Tasks:

- integrate infinite layers into multilayer accumulation
- verify off-axis footprints and source-aware extraction are preserved
- validate `fractional_cn2` variance behavior against the new backend

Acceptance:

- finite and infinite models agree on short-run local statistics
- infinite model avoids long-run periodic repetition

Implementation status:

- finite and infinite atmospheres now expose source-aware propagation from the
  current evolved layer state via per-layer footprint shifts and finite-height
  footprint scaling.
- the infinite backend now matches the finite backend on short-run local
  variance within regression tolerance and preserves `fractional_cn2` variance
  behavior.

### Work package 5: GPU runtime port

Files:

- `src/atmosphere/infinite_screen.jl`
- `src/core/backends.jl`
- `scripts/gpu_smoke_contract.jl`

Tasks:

- port runtime stepping kernels to `KernelAbstractions`
- preallocate boundary workspaces on the chosen backend
- add GPU smoke coverage for infinite atmosphere construction and stepping

Acceptance:

- AMDGPU smoke coverage passes
- CUDA path remains backend-correct where available
- no scalar indexing regressions appear in the new hot path

Implementation status:

- the infinite atmosphere stepping path now runs on both CPU and GPU targets
  through direction-specialized `KernelAbstractions` gather/injection kernels.
- GPU smoke coverage now includes infinite multilayer construction, stepping,
  propagation, and CPU/GPU statistical agreement checks.
- the one-time GPU kernel warmup cost is now paid during construction rather
  than the first integer boundary-injection step.

### Work package 6: Regression and benchmark gates

Files:

- `test/runtests.jl`
- `scripts/profile_*`
- `docs/benchmark-matrix-plan.md`

Tasks:

- add stationarity tests
- add long-run non-periodicity checks
- add deterministic replay for infinite layers
- add benchmark surfaces comparing:
  - finite CPU
  - finite GPU
  - infinite CPU
  - infinite GPU

Acceptance:

- the project has an explicit fidelity/performance comparison between finite and infinite models

Implementation status:

- long-run regression coverage now includes deterministic replay, diagonal-wind
  subpixel accumulation, extended non-periodicity checks, and windowed
  stationarity checks for the infinite model.
- `scripts/profile_atmosphere_runtime.jl` now reports the maintained finite vs
  infinite timing/allocation comparison on CPU or GPU backends.

## Required Tests

### Numerical tests

- covariance reconstruction of injected boundaries vs analytic reference
- short-run PSD agreement with the configured Von Karman model
- long-run covariance stationarity
- `fractional_cn2` multilayer variance preservation

### Behavioral tests

- zero-wind stationarity
- subpixel accumulation
- non-axis-aligned wind evolution
- off-axis footprint extraction
- deterministic replay under fixed RNG

### Regression tests

- infinite model does not repeat with the short periodicity of the finite canvas
- finite model remains available and unchanged for HIL-oriented runs
- CPU and GPU infinite paths agree statistically within documented tolerances

## Performance Gates

The infinite model is expected to be slower than the finite periodic model.
That is acceptable, but the slowdown must be measured and explicit.

Required benchmark outputs:

- build/precompute time
- per-step mean and p95
- per-step allocations
- backend synchronization count in profiling runs

Expected policy:

- finite periodic model remains the default fast/HIL path
- infinite model becomes the explicit science-fidelity path

## Risks And Mitigations

### Risk: unstable covariance conditioning

Mitigation:

- start with small stencils
- add conditioning diagnostics
- fail with structured errors rather than silent bad updates

### Risk: GPU implementation becomes branch-heavy

Mitigation:

- keep builder-time work separate from step-time work
- use direction-specialized kernels if needed
- benchmark axis-aligned and arbitrary-wind cases separately

### Risk: infinite model silently replaces the fast path

Mitigation:

- keep separate constructors/types
- require explicit model choice in runtime/profile setup

## Exit Criteria

This plan is complete when:

- `InfiniteMultiLayerAtmosphere` exists as a maintained model,
- the model satisfies `ATM-020` through `ATM-025` with boundary injection rather
  than periodic wraparound,
- CPU and GPU validation exist for the new runtime path,
- the finite model remains available as the explicit fast path,
- the benchmark and regression story makes the fidelity/performance tradeoff explicit.

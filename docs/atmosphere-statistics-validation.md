# Atmosphere Statistics Validation Note

Status: active

Plan traceability:

- [post-review-platform-plan.md](./post-review-platform-plan.md) `PVP-02`
- [model-validity-matrix.md](./model-validity-matrix.md) `MV-01`

Archived evidence:

- [2026-04-01-phase1-pvp02.toml](../benchmarks/results/atmosphere/2026-04-01-phase1-pvp02.toml)
- [manifest.toml](../benchmarks/results/atmosphere/manifest.toml)
- generator:
  [generate_atmosphere_statistics_artifact.jl](../scripts/generate_atmosphere_statistics_artifact.jl)

## Purpose

This note records the maintained finite/infinite atmosphere statistics artifact
for the current atmosphere stack.

It is not a claim of external frozen parity to OOPAO or SPECULA for the
infinite atmosphere model. It is a committed internal evidence artifact showing
that the maintained finite and infinite implementations satisfy the expected
statistical and behavioral checks under fixed deterministic scenarios.

## Scope

Artifact scope:

- finite model: `MultiLayerAtmosphere`
- infinite model: `InfiniteMultiLayerAtmosphere`
- pupil resolution: `16`
- diameter: `8.0 m`
- sampling time: `1e-3 s`
- `r0 = 0.2 m`
- `L0 = 25.0 m`

The artifact freezes three classes of evidence:

1. finite vs infinite ensemble standard-deviation agreement
2. infinite-model windowed stationarity
3. finite periodic wrap vs infinite non-periodicity

## Archived Results

From [2026-04-01-phase1-pvp02.toml](../benchmarks/results/atmosphere/2026-04-01-phase1-pvp02.toml):

### Ensemble standard deviation

- finite single-layer mean std: `0.1186199240799813`
- infinite single-layer mean std: `0.11074736556307828`
- single-layer ratio: `0.9336320725378746`
- infinite equal two-layer mean std: `0.10575280686420888`
- equal-layer ratio vs infinite single-layer: `0.9549013317519985`
- infinite uneven two-layer mean std: `0.1073086409730575`
- uneven-layer ratio vs infinite single-layer: `0.9689498294380449`

Interpretation:

- the maintained finite and infinite single-layer statistics remain within the
  intended `20%` agreement band used by the regression tests
- the infinite two-layer cases remain close to the infinite single-layer
  baseline rather than collapsing or drifting strongly

### Windowed stationarity

- early-window mean std: `0.09222158816585495`
- late-window mean std: `0.09274818020588944`
- late/early ratio: `1.0057100734275737`

Interpretation:

- the infinite model remains statistically stationary over the maintained
  windowed trajectory check

### Periodicity and non-periodicity

- finite screen period: `48` steps
- finite periodic replay exact: `true`
- infinite correlation after one finite screen period: `-0.32552788909479285`
- infinite long-run correlation: `-0.11522698479565403`

Interpretation:

- the finite periodic model still behaves like the maintained periodic canvas
- the infinite model remains clearly non-periodic against that same replay
  horizon

## Why There Is Still No External Frozen Infinite-Atmosphere Baseline

That remains intentional for now.

Reasons:

- the maintained infinite model is a Julia-native boundary-injection
  implementation rather than a direct frozen parity port
- the current external references do not already exist in this repository as a
  deterministic committed finite/infinite atmosphere statistics bundle
- freezing a cross-package infinite-screen bundle would require stronger
  scenario normalization than currently exists for:
  - stencil geometry
  - boundary-operator conditioning choices
  - RNG and boundary-sampling policy
  - source-aware extraction semantics

So the current maintained claim is:

- strong internal statistical and behavioral evidence
- backend parity and runtime profiling
- no external frozen infinite-atmosphere parity claim yet

## Maintained Acceptance Shape

This artifact is considered acceptable evidence for `MV-01` because it gives a
committed, reproducible summary for the exact properties the infinite model is
expected to preserve:

- variance scale
- stationarity
- deterministic replay under fixed RNG
- non-periodicity relative to the finite periodic canvas

If a later cross-package infinite-atmosphere normalization becomes practical,
this note should be updated to reference that frozen external bundle instead of
claiming only internal evidence.

Status: active

Plan traceability:

- [`PVP-06`](./post-review-platform-plan.md)
- [`MV-10`](./model-validity-matrix.md)

# Tomography Benchmark Scope

## Purpose

This note records the current maintained scope for tomography benchmark
evidence.

The package has strong tomography correctness evidence through:

- frozen pyTomoAO-aligned reference bundles
- dedicated tomography functional tests
- maintained GPU profiling scripts for compact and medium builder phases

What remains weaker is a truly representative, routinely refreshable
tomography benchmark artifact comparable to the closed-loop runtime surfaces.

## Current Evidence

The maintained evidence for tomography and reconstruction currently includes:

- frozen reference comparisons in [reference_harness.jl](../test/reference_harness.jl)
  and [reference_data](../test/reference_data)
- dedicated functional coverage in [tomography.jl](../test/tomography.jl)
- GPU-oriented profiling contracts in:
  - [gpu_profile_model_tomography_contract.jl](../scripts/gpu_profile_model_tomography_contract.jl)
  - [gpu_profile_model_tomography_phases_contract.jl](../scripts/gpu_profile_model_tomography_phases_contract.jl)

This is enough to support the current `medium-strong` validity rating for
`MV-10`, but it is not yet equivalent to a maintained representative benchmark
artifact of the type used for AO3k or the major WFS runtime families.

## Deferred Representative Artifact

A dedicated Phase 1 representative tomography artifact was attempted as part of
`PVP-06`, but it was deferred rather than committed.

Reason:

- even after trimming the candidate case and reducing sampling, the candidate
  model-based tomography builder path remained too expensive to serve as a
  routine refresh artifact
- the resulting refresh cost was not yet acceptable as a maintained evidence
  surface for normal package validation

So `PVP-06` is currently closed as a scoped defer, not as a completed
representative artifact.

Review date:

- 2026-06-30, or earlier if tomography becomes an active performance target

## What Counts As "Representative" Later

A future tomography benchmark artifact should:

- be clearly larger than the current compact builder/profiler contracts
- include both reconstructor build cost and steady-state solve cost
- record the key dimensions:
  - number of LGS
  - number of fit sources
  - lenslet count
  - DM count and command length
  - grid-mask dimensions
- be reproducible in a normal maintenance workflow without impractical wall time

## Review Trigger

This defer should be revisited when one of the following becomes true:

1. the tomography builder path is measurably faster after future algorithm or
   linear-algebra improvements
2. a smaller but still defensible representative-refresh contract is agreed and
   documented in [benchmark-matrix-plan.md](./benchmark-matrix-plan.md)
3. tomography becomes the next major performance focus and warrants a dedicated
   maintained artifact despite a longer refresh time

Until then, tomography remains validated primarily through frozen reference
data and functional coverage, with representative benchmark evidence explicitly
deferred.

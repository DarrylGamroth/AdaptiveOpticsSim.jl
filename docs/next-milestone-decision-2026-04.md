# Next Milestone Decision 2026-04

Status: active

Plan traceability:

- [`PVP-11`](./post-review-platform-plan.md)
- [`PVP-12`](./post-review-platform-plan.md)
- [`PVP-13`](./post-review-platform-plan.md)
- direction IDs: `DIR-03`, `DIR-04`, `DIR-05`

## Purpose

This memo selects the next core-platform breadth milestone after the
post-review validation and benchmark hardening phases.

It compares the three candidates defined in
[post-review-platform-plan.md](./post-review-platform-plan.md):

- `NB-01`: additional atmospheric field propagation scenarios
- `NB-02`: richer grouped sensing/runtime orchestration
- `NB-03`: a narrowly scoped controller/process-family addition

## Evidence Base

The decision is grounded in:

- [future-platform-direction.md](./future-platform-direction.md)
- [model-validity-matrix.md](./model-validity-matrix.md)
- [cross-package-benchmark-harness.md](./cross-package-benchmark-harness.md)
- [cross-package-benchmark-inventory.md](./cross-package-benchmark-inventory.md)
- [gpu-runtime-structural-refactor-plan.md](./gpu-runtime-structural-refactor-plan.md)

## Candidate Comparison

| Candidate | SPECULA gap closed | Current evidence footing | Core-boundary fit | Risk of premature framework expansion | Decision |
| --- | --- | --- | --- | --- | --- |
| `NB-01` additional atmospheric field propagation scenarios | medium | `MV-03` is already strong, and `CP-05` now has an explicit contract/defer record but no runnable external benchmark participant yet | good | low | defer |
| `NB-02` richer grouped sensing/runtime orchestration | high | grouped execution exists internally, runtime ownership was recently cleaned up, and `CP-04` remains a clear Julia-first breadth gap with SPECULA as the stronger architectural reference | strong | low-medium | choose |
| `NB-03` controller/process-family addition | uncertain | runtime and validation are stronger than before, but there is still no benchmark-backed reason to broaden controller/process families yet | mixed | high | defer |

## Decision

Choose `NB-02`: richer grouped sensing/runtime orchestration.

## Why `NB-02` Wins

`NB-02` best matches the current direction guardrails:

- it closes a real platform-breadth gap where SPECULA is the stronger
  reference
- it stays inside the core package boundary
- it builds on recently improved shared runtime/product/scheduling surfaces
- it does not require optional science tooling
- it does not force early controller/process-family expansion

`NB-01` remains important, but the atmospheric-field surface is already in a
better state than the grouped orchestration surface:

- the model-validity footing is already strong for `MV-03`
- frozen SPECULA-targeted contracts already exist
- the main missing benchmark surface is now explicitly documented by the
  `CP-05` defer note rather than being an ambiguous hole

`NB-03` is not selected because the current evidence still does not justify
expanding the controller/process-family surface:

- no realistic cross-package benchmark currently demands it
- it would widen core scope faster than the validation surface
- it conflicts with the explicit direction to defer broad controller growth
  until it is benchmark-justified

## Chosen Scope

The chosen milestone should focus on:

- grouped sensing execution semantics
- grouped runtime product ownership
- grouped detector/WFS/science export planning
- grouped benchmark and validation surfaces

The chosen milestone should not include:

- new broad controller families
- optional science-path integrations
- open-ended atmospheric-field parity chasing

## Deferred Candidates

### Deferred `NB-01`

Revisit only when at least one of the following becomes true:

- a runnable external `CP-05` participant is normalized into the maintained
  harness
- a concrete missing atmospheric-field scenario blocks grouped-runtime work
- a new validation gap appears in `MV-03`

### Deferred `NB-03`

Revisit only when all of the following become true:

- grouped runtime ownership is stable after the chosen `NB-02` work
- a controller/process-family use case has a maintained benchmark target
- validation expectations are explicit before implementation begins

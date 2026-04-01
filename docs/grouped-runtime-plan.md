# Grouped Runtime Plan

Status: active

Plan traceability:

- chosen by [`PVP-12`](./post-review-platform-plan.md)
- decision memo: [next-milestone-decision-2026-04.md](./next-milestone-decision-2026-04.md)
- direction IDs: `DIR-03`, `DIR-04`, `DIR-05`

## Purpose

This plan defines the next selected platform-breadth milestone:

- `NB-02`: richer grouped sensing/runtime orchestration

The goal is to make grouped execution a first-class, benchmarked core runtime
surface rather than a collection of partially shared execution paths.

## Scope

In scope:

- grouped source and grouped WFS execution contracts
- grouped runtime product planning and export ownership
- grouped detector/WFS/science pipeline reuse
- grouped benchmark surfaces and validation notes

Out of scope:

- new controller/process families
- optional science-path adapters
- unrelated atmospheric-field feature expansion unless required by grouped
  runtime execution

## Milestones

### GR-1: Grouped Contract Audit

Deliverables:

- maintained grouped-contract note documenting:
  - grouped source semantics
  - grouped WFS product ownership
  - grouped detector export rules
  - grouped science export rules
- conformance checklist against current grouped helpers

Acceptance:

- one maintained contract note exists
- current grouped runtime entry points are mapped to that contract

### GR-2: Grouped Product Ownership Cleanup

Deliverables:

- shared grouped runtime product plan for:
  - slopes
  - WFS frames
  - grouped detector cubes or mosaics
  - science frames
- explicit separation between grouped scratch state and exported grouped
  products

Acceptance:

- grouped execution no longer relies on ambiguous in-place product ownership
- grouped export semantics are test-covered

### GR-3: Grouped Runtime Pipeline Reuse

Deliverables:

- grouped SH / Pyramid / BioEdge paths use the same grouped execution skeleton
  where behavior is truly shared
- grouped detector-backed execution reuses the shared detector pipeline where
  possible

Acceptance:

- duplication in grouped runtime hot paths is reduced
- no regression in existing grouped runtime tests

### GR-4: Grouped Validation Surface

Deliverables:

- grouped validation note added to the validity matrix
- one committed grouped benchmark artifact for a maintained Julia-first grouped
  scenario family
- optional backend grouped smoke remains explicit and separate from benchmark
  evidence

Acceptance:

- the grouped runtime family has a maintained validation/evidence surface
- benchmark evidence is archived and reproducible

## Validation Expectations

Required evidence for milestone completion:

- functional tests for grouped export/product ownership
- benchmark evidence for at least one maintained grouped runtime family
- explicit limitations note where grouped cross-package equivalence is not yet
  claimed

Preferred scripts and surfaces:

- [profile_multi_source_multi_wfs_runtime.jl](../scripts/profile_multi_source_multi_wfs_runtime.jl)
- [gpu_smoke_contract.jl](../scripts/gpu_smoke_contract.jl)
- grouped runtime tests under [test/testsets](../test/testsets)

## Non-Goals And Guardrails

- Do not widen this plan into controller/process-family work.
- Do not add optional science tooling to core.
- Do not claim SPECULA equivalence where only Julia-side grouped evidence
  exists.
- Do not broaden grouped runtime features without updating the validation and
  benchmark surfaces.

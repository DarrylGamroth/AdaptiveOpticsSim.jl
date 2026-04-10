# Production Hardening Plan: 2026-04

Status: active

Use together with:

- [package-review-status-2026-04.md](./package-review-status-2026-04.md)
- [production-readiness-checklist.md](./production-readiness-checklist.md)
- [supported-production-surfaces.md](./supported-production-surfaces.md)
- [api-tier-inventory.md](./api-tier-inventory.md)
- [modularization-inventory.md](./modularization-inventory.md)
- [model-validity-matrix.md](./model-validity-matrix.md)

## Purpose

This plan turns the current package-review status into an execution order for
getting `AdaptiveOpticsSim.jl` from "credible and well-validated on maintained
surfaces" to "production-ready on a clearly scoped AO platform surface".

This is not a new roadmap for all future features. It is a short-to-medium-term
hardening plan focused on:

- curating the public API
- reducing structural maintenance risk
- widening validation in the most defensible way
- making supported-surface release claims easy to verify

## Current Starting Point

The package is already in a strong position on its supported surfaces:

- maintained CPU test coverage is broad
- CUDA and AMDGPU parity/runtime coverage is materially stronger than it was in
  March 2026
- OOPAO-aligned HEART equivalence is strong on the maintained null/noise-free
  boundary
- a second frozen OOPAO equivalence artifact exists for Pyramid and BioEdge
- a scientist-owned HEART boundary truth artifact now exists

The biggest remaining risks are structural and operational rather than
language-level:

1. the public API is still too broad
2. some source and test surfaces are still too large or too loosely curated
3. production claims depend on maintained operational GPU validation, not just
   checked-in workflow files
4. the package still needs one more layer of maintainer-facing synthesis around
   what is stable, expert-facing, or experimental

## Non-Goals

This plan does not aim to:

- make every exported symbol production-supported
- resolve SPECULA equivalence before shipping supported surfaces
- broaden controller/process breadth to match SPECULA in this pass
- split the package into multiple packages
- redesign the scientific core around another language or runtime

## Success Criteria

This plan is complete when all of the following are true:

1. there is a smaller, explicitly tiered top-level public API
2. the highest-risk oversized files are either split further or explicitly
   justified as stable entry points over modular internals
3. production-supported surfaces and release-gating evidence are easy to audit
   from a small set of docs and scripts
4. routine real-hardware GPU validation is operational through either CI or a
   maintained release-validation host cadence
5. at least one more realistic maintained fidelity surface exists beyond the
   current null/noise-free HEART baseline, or the current scope is explicitly
   frozen as the production boundary for this release train

## Execution Order

### Phase PH-1: Freeze the production boundary and release gate

Goal:

- make the supported claim auditable from a small number of maintained files

Tasks:

- confirm the production boundary docs remain aligned:
  - [supported-production-surfaces.md](./supported-production-surfaces.md)
  - [production-readiness-checklist.md](./production-readiness-checklist.md)
  - [release-validation-runbook.md](./release-validation-runbook.md)
- ensure the release procedure references all required evidence:
  - CPU suite
  - CUDA parity/runtime surfaces
  - AMDGPU parity/runtime surfaces
  - OOPAO equivalence artifacts
  - HEART truth artifact
- define an operational cadence for real-hardware validation:
  - CI, or
  - release-validation host

Exit criteria:

- one person can determine release readiness by reading only:
  - [supported-production-surfaces.md](./supported-production-surfaces.md)
  - [production-readiness-checklist.md](./production-readiness-checklist.md)
  - [release-validation-runbook.md](./release-validation-runbook.md)
- the required runtime/equivalence/truth artifacts are all linked there

Verification:

- dry-run the release validation procedure
- verify doc links and commands are current

### Phase PH-2: Public API curation pass 2

Goal:

- reduce the top-level namespace to a clearer stable workflow surface without
  breaking maintained production flows

Tasks:

- use [api-tier-inventory.md](./api-tier-inventory.md) as the source of truth
- identify the next low-risk tranche of de-export candidates, prioritizing:
  - backend/build helper types not needed for normal users
  - telemetry/config helper surfaces that already work via namespaced access
  - expert-only infrastructure hooks that are documented but not workflow-level
- keep user-facing workflow ergonomics intact for:
  - telescope/source/atmosphere construction
  - maintained detector/WFS constructors
  - runtime and closed-loop entry points
  - validation-facing helper surfaces that are part of supported workflows
- update API docs and migration notes for any de-exports

Exit criteria:

- top-level export surface is measurably smaller
- stable vs advanced vs internal surfaces are clearer in docs
- maintained package scripts and supported workflows still run

Verification:

- export count recorded before/after
- `Pkg.test()` on maintained CPU surface
- optional backend coverage for any moved GPU/backend helpers touched by the pass

### Phase PH-3: Structural modularization pass 2

Goal:

- reduce maintenance risk in the largest remaining files and test surfaces

Priority targets:

1. [src/AdaptiveOpticsSim.jl](../src/AdaptiveOpticsSim.jl)
   - keep as a curated entry point only
   - move grouped export/comment policy into clearer include groupings if needed
2. [test/testsets/detectors_and_wfs.jl](../test/testsets/detectors_and_wfs.jl)
   - split by maintained family or validation role:
     - detectors
     - Shack-Hartmann
     - Pyramid/BioEdge
     - backend parity adjuncts
3. any remaining oversized subsystem file that still mixes:
   - types
   - setup
   - execution
   - calibration
   - backend-specific policy

Tasks:

- use [modularization-inventory.md](./modularization-inventory.md)
- prefer internal file splits over public package-boundary changes
- preserve package load order and include clarity
- keep new file boundaries aligned with subsystem responsibilities, not micro-files

Exit criteria:

- no single maintained test file dominates subsystem coverage unnecessarily
- major subsystem entry files remain short include front doors where practical
- no new public API churn is introduced by the split

Verification:

- `Pkg.test()` full maintained CPU suite
- optional backend suite if touched files participate in GPU code paths

### Phase PH-4: Validation breadth pass 2

Goal:

- strengthen the production claim with one additional maintained realism or
  fidelity surface beyond the current strongest HEART null/noise-free case

Preferred target order:

1. a realistic HEART fidelity surface with turbulence/noise kept deterministic
   enough to compare within defensible tolerances
2. a maintained PWFS runtime/fidelity surface tied to the frozen OOPAO bundle
3. a grouped/runtime platform artifact if it is already stable enough to claim

Tasks:

- choose one surface only for this pass
- define exact contract and tolerances
- record CPU baseline artifact first
- add CUDA/AMDGPU parity or runtime evidence only if the surface is part of the
  supported accelerator claim

Exit criteria:

- one more maintained artifact exists with:
  - explicit contract
  - archived result
  - linked verification guidance
- production docs clearly state whether it is part of the supported scope or a
  supporting confidence artifact

Verification:

- artifact generation command committed
- archived result committed
- docs linked from the release/readiness docs

### Phase PH-5: Operational validation cadence

Goal:

- make the supported GPU claim operational rather than aspirational

Tasks:

- choose the real validation path:
  - self-hosted CI runners, or
  - a maintained release-validation host procedure
- define the minimum cadence:
  - per release
  - per tagged candidate
  - or another explicit release gate
- record where validation logs/artifacts live
- add a short operational checklist for who runs it and when

Exit criteria:

- the team has one real, repeatable GPU validation path
- production claims do not depend on ad hoc memory or private local procedure

Verification:

- run the chosen cadence once and record the result location

## Recommended Implementation Sequence

Do the next work in this order:

1. `PH-1` production boundary/release gate cleanup
2. `PH-2` API curation pass 2
3. `PH-3` structural modularization pass 2
4. `PH-5` operational validation cadence
5. `PH-4` additional realism/fidelity surface

Rationale:

- `PH-1` and `PH-2` reduce confusion immediately
- `PH-3` lowers maintenance risk before broader new validation work
- `PH-5` is required for a defensible production claim on GPU-supported surfaces
- `PH-4` should be added once the package boundary and release gate are already
  clearer

## Immediate Next Slice

The next implementation slice should be:

1. `PH-2` API curation pass 2
2. focused on one coherent tranche from [api-tier-inventory.md](./api-tier-inventory.md)
3. with explicit migration notes and zero ambiguity about supported workflow
   breakage

That is the highest-value codebase change because it improves usability,
maintainability, and public-surface discipline without requiring new scientific
scope.

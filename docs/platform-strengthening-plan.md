# Platform Strengthening Plan

Status: active

Traceability sources:

- [future-platform-direction.md](./future-platform-direction.md)
- [package-review-2026-03.md](./package-review-2026-03.md)
- [model-validity-matrix.md](./model-validity-matrix.md)
- [execution-plan-closeout.md](./execution-plan-closeout.md)

## Purpose

This plan governs the next package phase after the execution-plan rollout.

Its purpose is to address the remaining high-value platform gaps without
drifting into config-first design, parity chasing, or premature science-path
integration in core.

This plan is intentionally focused on:

- API shape and synthesis-oriented usability
- stronger validation and benchmark evidence
- further ROCm/AMDGPU runtime hardening
- richer platform-scale orchestration and composition

This plan intentionally does not make declarative manifest/config files the
primary interface.

## Direction Requirements

These IDs are the stable requirements for this plan.

| ID | Requirement | Source |
| --- | --- | --- |
| `PSR-01` | Scripts and typed Julia builders remain the primary interface. | [future-platform-direction.md](./future-platform-direction.md) |
| `PSR-02` | Config-file-first architecture is out of scope for this phase. | user direction, [future-platform-direction.md](./future-platform-direction.md) |
| `PSR-03` | Address remaining review gaps before broad new capability expansion. | [package-review-2026-03.md](./package-review-2026-03.md) |
| `PSR-04` | Strengthen validation and representative benchmark evidence where current claims are still narrower than ideal. | [model-validity-matrix.md](./model-validity-matrix.md) |
| `PSR-05` | Reduce backend-specific workaround debt where it still affects maintained GPU paths. | [execution-plan-closeout.md](./execution-plan-closeout.md) |
| `PSR-06` | Use SPECULA as the stronger reference for platform-scale breadth, but keep the implementation Julia-native. | [future-platform-direction.md](./future-platform-direction.md) |
| `PSR-07` | Science-path / focal-plane integration remains outside core unless a later explicit boundary review changes that rule. | [future-platform-direction.md](./future-platform-direction.md) |

## Primary Interface Policy

The maintained interface policy for this phase is:

- primary:
  - Julia scripts
  - typed constructors
  - typed scenario/build orchestration objects
  - benchmark and validation runners
- deferred:
  - optional scenario manifests
  - config-file-driven primary workflows

The package should prefer:

- programmable Julia-native composition,
- reproducible script entry points,
- and typed scenario objects

over:

- opaque declarative configuration as the main user experience.

If a manifest/config format is revisited later, it should sit on top of the
typed platform model, not replace it.

## Current Gaps

The remaining gaps this plan is meant to address are:

1. the public API is still broader and flatter than ideal
2. documentation still needs more synthesis-oriented platform guidance
3. tomography representative benchmark evidence is still explicitly deferred
4. SPECULA-aligned validation is strong, but still largely contract-oriented
5. some maintained host-mirror fallback behavior still exists on AMDGPU
6. platform-scale orchestration/composition is still weaker than the strongest
   SPECULA-style platform surfaces
7. science-path / focal-plane integration should be considered, but only as an
   optional boundary, not as new core scope

## Non-Goals

This plan does not directly implement:

- config-file-first scenario control
- broad controller/process-family expansion
- core-package science-camera or coronagraph integration
- parity work justified only because OOPAO or SPECULA have a feature

## Phase Overview

| Phase | Goal | Exit gate |
| --- | --- | --- |
| `PH-1` | tighten platform usability and synthesis docs | API layering and synthesis docs are materially improved without breaking maintained workflows |
| `PH-2` | strengthen evidence where current claims are still narrow | tomography and SPECULA-aligned evidence gaps are either closed or explicitly re-scoped |
| `PH-3` | reduce remaining maintained AMDGPU workaround debt | benchmark-backed ROCm fallback inventory is reduced or explicitly justified |
| `PH-4` | introduce a Julia-native platform orchestration layer | typed scenario/composition surfaces exist for richer multi-branch runtime assembly |
| `PH-5` | validate platform orchestration on realistic maintained scenarios | new orchestration surfaces have benchmark and validation evidence on maintained backends |
| `PH-6` | close out and record optional deferred tracks | next-step decision is explicit and science-path/config tracks remain intentionally deferred |

## Phase 1: Platform Usability And Synthesis

### Goal

Make the platform easier to understand and use without changing the primary
Julia-native workflow model.

### Milestones

| Milestone | Task | Outputs |
| --- | --- | --- |
| `PSP-01` | Curate the next tranche of low-value exports into namespaced access. | export reduction patch, updated API tier inventory |
| `PSP-02` | Add a synthesis-oriented platform architecture guide aimed at users moving beyond single tutorials. | new stable guide |
| `PSP-03` | Add a platform workflow guide explaining how atmosphere, runtime, WFS, detector, and benchmark entry points fit together. | new stable guide or expanded user guide section |
| `PSP-04` | Add a scenario-builder style guide that defines how Julia scripts should express reproducible scenarios without config files. | style guide / conventions note |

### Acceptance

- the maintained top-level API is smaller or at least better tiered
- new platform docs explain composition more directly than plan docs do
- the recommended user path is still scripts + typed objects, not config files

## Phase 2: Evidence Strengthening

### Goal

Close the highest-value remaining validation and benchmark gaps.

### Milestones

| Milestone | Task | Outputs |
| --- | --- | --- |
| `PSP-05` | Convert the tomography representative benchmark defer into a maintained benchmark artifact, or explicitly re-scope it with benchmark-backed reasoning if it is still not practical. | benchmark artifact or updated defer note |
| `PSP-06` | Expand SPECULA-aligned validation beyond current contract-only cases where the platform claim would benefit from broader evidence. | new frozen bundle cases and updated validity matrix |
| `PSP-07` | Add at least one platform-scale SPECULA-informed scenario comparing orchestration/runtime behavior rather than only isolated optical kernels. | new benchmark/validation contract and archived result |

### Acceptance

- `MV-10` is improved or re-scoped with stronger evidence
- SPECULA-aligned evidence covers at least one broader runtime/platform scenario
- new claims are reflected in [model-validity-matrix.md](./model-validity-matrix.md)

## Phase 3: AMDGPU Workaround Reduction

### Goal

Reduce the remaining maintained host-mirror fallback surfaces where that yields
real value.

### Milestones

| Milestone | Task | Outputs |
| --- | --- | --- |
| `PSP-08` | Write and commit a current ROCm fallback inventory after the execution-plan rollout, with each fallback labeled as temporary, acceptable, or next-target-for-removal. | inventory note |
| `PSP-09` | Remove or shrink the highest-value remaining ROCm fallback in a benchmark-backed hot path. | code change + benchmark evidence |
| `PSP-10` | Re-baseline CPU/AMDGPU/CUDA realistic runtime surfaces after the ROCm cleanup. | updated benchmark note/results |

### Acceptance

- at least one maintained ROCm host-mirror fallback is removed or materially reduced
- grouped and AO3k benchmark surfaces remain in-family
- residual fallback behavior is explicit rather than ambient

## Phase 4: Julia-Native Platform Orchestration

### Goal

Build a stronger platform-scale orchestration layer without switching to
config-first design.

### Design Rule

The primary interface remains:

- typed scenario/config structs
- Julia builders
- scripts that construct and run those scenarios

This phase should introduce richer orchestration objects, not external config
files.

### Milestones

| Milestone | Task | Outputs |
| --- | --- | --- |
| `PSP-11` | Define typed platform scenario/config families for representative runtime compositions. | scenario/config types and guide |
| `PSP-12` | Add builders for representative multi-branch platform compositions. | builder surfaces for grouped/multi-WFS/multi-detector cases |
| `PSP-13` | Add runtime entry points that make these scenario objects easy to execute, benchmark, and validate from scripts. | maintained runner/API surfaces |
| `PSP-14` | Add one or two canonical main-platform example scripts that use the new orchestration layer. | examples and docs |

### Acceptance

- platform composition is more explicit and reusable than today’s ad hoc script assembly
- scripts remain the main entry point
- no declarative manifest layer is required to use the new platform model

## Phase 5: Platform-Scale Validation And Benchmarking

### Goal

Back the new orchestration layer with realistic evidence.

### Milestones

| Milestone | Task | Outputs |
| --- | --- | --- |
| `PSP-15` | Add maintained realistic benchmark surfaces for the new scenario/composition layer on CPU and maintained GPUs. | benchmark runners/results |
| `PSP-16` | Add one cross-package platform scenario contract where normalization is strong enough to make comparison meaningful. | contract doc + archived baseline |
| `PSP-17` | Add backend parity/functional coverage for the new orchestration surfaces. | tests/smoke additions |

### Acceptance

- new orchestration surfaces have maintained benchmark evidence
- at least one realistic cross-package scenario is normalized and archived
- CPU/AMDGPU/CUDA claims are explicit and tested where maintained

## Phase 6: Closeout And Deferred Tracks

### Goal

Finish the phase with an explicit decision record rather than drifting into the
next topic.

### Milestones

| Milestone | Task | Outputs |
| --- | --- | --- |
| `PSP-18` | Record what remains behind SPECULA after the orchestration pass. | decision/closeout note |
| `PSP-19` | Explicitly defer or schedule science-path / focal-plane work under the optional boundary rule. | roadmap/boundary update |
| `PSP-20` | Explicitly defer or schedule a future optional manifest/config layer, if still desired, as a secondary interface on top of the typed platform model. | defer note or future-plan link |

### Acceptance

- the next platform step is explicit
- science-path remains optional unless deliberately changed
- config manifests remain deferred and non-primary unless a later plan changes that rule

## Execution Order

1. `PSP-01`
2. `PSP-02`
3. `PSP-03`
4. `PSP-04`
5. `PSP-05`
6. `PSP-06`
7. `PSP-07`
8. `PSP-08`
9. `PSP-09`
10. `PSP-10`
11. `PSP-11`
12. `PSP-12`
13. `PSP-13`
14. `PSP-14`
15. `PSP-15`
16. `PSP-16`
17. `PSP-17`
18. `PSP-18`
19. `PSP-19`
20. `PSP-20`

## Traceability Matrix

| Requirement | Covered by milestones | Completion rule |
| --- | --- | --- |
| `PSR-01` | `PSP-02`, `PSP-03`, `PSP-04`, `PSP-11` to `PSP-14` | scripts and typed builders remain the primary documented workflow |
| `PSR-02` | `PSP-04`, `PSP-20` | config-first design remains explicitly deferred |
| `PSR-03` | `PSP-01` to `PSP-10` | the remaining review-quality gaps are materially reduced |
| `PSR-04` | `PSP-05` to `PSP-07`, `PSP-15` to `PSP-17` | validation and benchmark evidence are stronger where current claims are still narrow |
| `PSR-05` | `PSP-08` to `PSP-10` | ROCm workaround debt is reduced or explicitly justified |
| `PSR-06` | `PSP-06`, `PSP-07`, `PSP-16`, `PSP-18` | SPECULA-informed platform work is evidence-backed and Julia-native |
| `PSR-07` | `PSP-19` | science-path work remains outside core unless deliberately revisited |

## Checklist

### Phase 1

- [ ] `PSP-01`
- [ ] `PSP-02`
- [ ] `PSP-03`
- [ ] `PSP-04`

### Phase 2

- [ ] `PSP-05`
- [ ] `PSP-06`
- [ ] `PSP-07`

### Phase 3

- [ ] `PSP-08`
- [ ] `PSP-09`
- [ ] `PSP-10`

### Phase 4

- [ ] `PSP-11`
- [ ] `PSP-12`
- [ ] `PSP-13`
- [ ] `PSP-14`

### Phase 5

- [ ] `PSP-15`
- [ ] `PSP-16`
- [ ] `PSP-17`

### Phase 6

- [ ] `PSP-18`
- [ ] `PSP-19`
- [ ] `PSP-20`

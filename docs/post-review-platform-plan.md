# Post-Review Platform Plan

Status: active

Traceability sources:

- [future-platform-direction.md](./future-platform-direction.md)
- [model-validity-matrix.md](./model-validity-matrix.md)
- [cross-package-benchmark-inventory.md](./cross-package-benchmark-inventory.md)
- [cross-package-benchmark-harness.md](./cross-package-benchmark-harness.md)

## Purpose

This plan governs the next package phase after the completed review-cleanup and
interface-spec work.

Its goal is to make the next steps explicit and traceable so that:

- validation/benchmark strengthening is done before more platform breadth
- SPECULA-informed expansion happens from evidence rather than intuition
- controller/process growth is deferred until runtime and validation surfaces
  are strong enough
- optional science-path integrations remain outside the core package boundary

## Direction Requirements

These IDs are the stable directional requirements for this plan.

| ID | Requirement | Source |
| --- | --- | --- |
| `DIR-01` | Strengthen validation evidence where the maintained model-validity matrix is still medium or medium-strong rather than strong. | [future-platform-direction.md](./future-platform-direction.md), [model-validity-matrix.md](./model-validity-matrix.md) |
| `DIR-02` | Strengthen benchmark evidence where realistic runtime/HIL claims exist but cross-package or representative evidence is still thin. | [future-platform-direction.md](./future-platform-direction.md), [cross-package-benchmark-inventory.md](./cross-package-benchmark-inventory.md) |
| `DIR-03` | Use SPECULA as the stronger reference for the next broader capability work rather than chasing more OOPAO parity. | [future-platform-direction.md](./future-platform-direction.md) |
| `DIR-04` | Do not start broader controller/process-family expansion until runtime ownership, benchmark surfaces, and validation expectations are explicit. | [future-platform-direction.md](./future-platform-direction.md) |
| `DIR-05` | Keep science-path integrations optional and outside core unless a later explicit boundary review changes that rule. | [future-platform-direction.md](./future-platform-direction.md) |
| `DIR-06` | Use frozen bundles, maintained benchmark contracts, and committed evidence rather than ad hoc comparisons or parity claims. | [model-validity-matrix.md](./model-validity-matrix.md), [cross-package-benchmark-harness.md](./cross-package-benchmark-harness.md) |

## Plan Scope

This plan covers three work classes:

1. closing high-value validation and benchmark gaps already documented
2. extending cross-package evidence on realistic scenarios
3. selecting the next SPECULA-informed breadth milestone from an evidence-backed shortlist

This plan does not directly implement:

- `Proper.jl` or other science-path integrations in core
- broad new controller/process families
- large exported-surface cleanup work

## Current Gaps To Address

These are the concrete surfaces that justify this plan.

### Validation Gaps

- `MV-01`: no frozen external infinite-atmosphere statistics bundle yet
- `MV-02`: no standalone maintained report for `K_{5/6}` helper accuracy
- `MV-03`: no frozen external atmospheric-field bundle yet
- `MV-04`: detector-family frozen references are still limited
- `MV-08`: LiFT/gain-sensing runtime evidence is lighter than core WFS/runtime surfaces
- `MV-10`: tomography representative benchmark evidence remains thin

### Cross-Package Benchmark Gaps

- `CP-02`: REVOLT-like SH runtime exists, but normalization and archived evidence
  can be stronger
- `CP-03`: REVOLT-like PWFS representative comparison is not yet first-class in
  `main`
- `CP-05`: atmospheric field propagation / curvature-through-atmosphere has no
  frozen external SPECULA-aligned benchmark contract yet

### Capability-Selection Gaps

- there is no committed decision artifact yet for which SPECULA-informed
  breadth area comes next after evidence strengthening
- candidate next areas are currently described only at the directional level:
  - additional atmospheric field propagation scenarios
  - richer grouped sensing/runtime orchestration
  - carefully selected controller/process families later

## Phase Overview

| Phase | Goal | Exit gate |
| --- | --- | --- |
| `PH-1` | Strengthen maintained validation evidence | selected medium/medium-strong `MV-*` gaps have committed evidence or an explicit defer decision |
| `PH-2` | Strengthen realistic cross-package benchmark evidence | maintained benchmark contracts and archived baselines exist for selected `CP-*` scenarios |
| `PH-3` | Choose the next SPECULA-informed breadth milestone | a committed decision doc names one next milestone with acceptance criteria and non-goals |

## Phase 1: Validation Evidence Hardening

### Goal

Raise the quality of maintained model evidence before expanding platform
breadth.

### Tasks

| Task ID | Task | Primary outputs |
| --- | --- | --- |
| `PVP-01` | Add an explicit phase-statistics accuracy note for the shared `K_{5/6}` helper and its validation range. | accuracy note doc or section, linked from the validity matrix |
| `PVP-02` | Add a frozen or archived finite/infinite atmosphere statistics artifact, or explicitly document why an external frozen baseline is not practical yet. | bundle or maintained note with acceptance bounds |
| `PVP-03` | Add a frozen external or contract-oriented atmospheric-field reference surface aligned to SPECULA-style scenarios. | committed bundle and generator, or an explicit justified defer note |
| `PVP-04` | Add at least one stronger detector-family validation artifact beyond integrated runtime smoke. | detector-oriented reference or maintained analytic/fixture report |
| `PVP-05` | Add stronger runtime/profile evidence for LiFT and gain-sensing workflows. | maintained benchmark/report artifact |
| `PVP-06` | Add representative tomography benchmark evidence or an explicit scope-limitation note if that benchmark remains deferred. | benchmark script/result note or defer decision |

### Work Steps

1. For each targeted `MV-*` family, write down the exact missing evidence type:
   `A`, `R`, `G`, `P`, or `M`.
2. Prefer the smallest committed artifact that closes the evidence gap:
   - a frozen bundle
   - a maintained benchmark result file
   - a concise explicit non-equivalence or defer note
3. Update [model-validity-matrix.md](./model-validity-matrix.md) immediately as
   each artifact lands.
4. Record generation/refresh surfaces for every new artifact.

### Acceptance Criteria

- `PVP-01` through `PVP-06` are each either:
  - implemented with linked evidence, or
  - explicitly deferred with a reason and review date
- [model-validity-matrix.md](./model-validity-matrix.md) reflects the new
  status without stale wording
- no new parity claim is introduced without either frozen evidence or a scoped
  limitation note

## Phase 2: Cross-Package Benchmark Expansion

### Goal

Turn the current benchmark harness into a stronger decision surface for realistic
cross-package comparisons.

### Tasks

| Task ID | Task | Primary outputs |
| --- | --- | --- |
| `PVP-07` | Normalize and archive the REVOLT-like SH comparison contract (`CP-02`) more strongly. | updated harness contract, archived results, comparability note |
| `PVP-08` | Promote the REVOLT-like PWFS comparison (`CP-03`) to a first-class maintained benchmark family in `main`. | contract entry, runner support, archived baseline |
| `PVP-09` | Add a SPECULA-aligned atmospheric-field benchmark family (`CP-05`) with an explicit scenario contract. | contract entry, runner support, archived baseline or scoped defer |
| `PVP-10` | Add a maintained benchmark-results manifest summarizing which scenarios are mandatory, optional, or environment-gated. | manifest update and benchmark guide update |

### Work Steps

1. For each candidate family, normalize:
   - atmosphere assumptions
   - source/wavelength assumptions
   - detector/readout assumptions
   - calibration/reconstructor assumptions
2. Encode those assumptions in the contract file before archiving new results.
3. Run the maintained harness and archive results in a dated directory or
   versioned baseline file.
4. Update benchmark docs with:
   - how to reproduce
   - what is considered equivalent
   - what is intentionally not normalized yet

### Acceptance Criteria

- `CP-02` has a stronger archived baseline than the current minimal Phase 5
  result
- `CP-03` exists as a maintained first-class family in `main` or is explicitly
  deferred with a blocking dependency note
- `CP-05` has either:
  - a maintained benchmark contract and archived result, or
  - an explicit defer note tied to a missing external scenario asset
- benchmark docs clearly separate:
  - fidelity comparisons
  - runtime comparisons
  - optional/skipped environment-dependent runs

## Phase 3: Next SPECULA-Informed Breadth Selection

### Goal

Choose one next platform-breadth milestone from evidence rather than opinion.

### Candidate Next Milestones

- `NB-01`: additional atmospheric field propagation scenarios
- `NB-02`: richer grouped sensing/runtime orchestration
- `NB-03`: a narrowly scoped controller/process-family addition, only if
  evidence now justifies it

### Tasks

| Task ID | Task | Primary outputs |
| --- | --- | --- |
| `PVP-11` | Write a short decision memo comparing `NB-01`, `NB-02`, and `NB-03` against the strengthened evidence. | decision memo |
| `PVP-12` | Commit one chosen next milestone plan with scope, non-goals, acceptance criteria, and validation expectations. | next-milestone plan doc |
| `PVP-13` | Explicitly defer the non-chosen candidates so they do not become ambient scope creep. | roadmap or decision doc updates |

### Decision Criteria

- Which candidate closes the highest-value capability gap relative to SPECULA?
- Which candidate now has the strongest validation and benchmark footing?
- Which candidate can stay within the core-package boundary without dragging in
  optional science tooling?
- Which candidate does not prematurely force a broad controller/process
  framework?

### Acceptance Criteria

- there is one committed chosen milestone
- the chosen milestone has clear validation and benchmark expectations
- the non-chosen candidates are documented as deferred, not forgotten

## Execution Order

1. `PVP-01`
2. `PVP-02`
3. `PVP-03`
4. `PVP-04`
5. `PVP-05`
6. `PVP-06`
7. `PVP-07`
8. `PVP-08`
9. `PVP-09`
10. `PVP-10`
11. `PVP-11`
12. `PVP-12`
13. `PVP-13`

## Verification And Evidence

Every completed task should record at least one of:

- committed frozen bundle path
- benchmark result artifact
- generator/runner script path
- test path
- maintained note or limitation record

Preferred evidence order:

1. test-backed or frozen-bundle evidence
2. maintained benchmark evidence
3. explicit scoped limitation note

## Traceability Matrix

| Direction ID | Satisfied by tasks | Completion rule |
| --- | --- | --- |
| `DIR-01` | `PVP-01` to `PVP-06` | targeted `MV-*` gaps are closed or explicitly deferred |
| `DIR-02` | `PVP-07` to `PVP-10` | realistic benchmark contracts and archived baselines exist |
| `DIR-03` | `PVP-03`, `PVP-09`, `PVP-11`, `PVP-12` | SPECULA-guided breadth selection is evidence-backed |
| `DIR-04` | `PVP-11`, `PVP-12`, `PVP-13` | controller/process growth remains deferred unless explicitly justified |
| `DIR-05` | `PVP-11`, `PVP-12` | chosen next milestone stays within the core/optional boundary rules |
| `DIR-06` | `PVP-02`, `PVP-03`, `PVP-07`, `PVP-08`, `PVP-09`, `PVP-10` | frozen bundles/contracts/results replace informal parity claims |

## Checklist

### Phase 1

- `[x]` `PVP-01` phase-statistics accuracy evidence
- `[ ]` `PVP-02` atmosphere statistics artifact or defer note
- `[ ]` `PVP-03` atmospheric-field external reference artifact or defer note
- `[ ]` `PVP-04` detector-family validation artifact
- `[ ]` `PVP-05` LiFT/gain-sensing runtime evidence
- `[ ]` `PVP-06` tomography representative evidence or defer note

### Phase 2

- `[ ]` `PVP-07` stronger `CP-02` archived baseline
- `[ ]` `PVP-08` first-class `CP-03` benchmark family
- `[ ]` `PVP-09` `CP-05` benchmark family or scoped defer
- `[ ]` `PVP-10` benchmark-results manifest and guide update

### Phase 3

- `[ ]` `PVP-11` next-milestone comparison memo
- `[ ]` `PVP-12` committed next-milestone plan
- `[ ]` `PVP-13` explicit defer records for non-chosen candidates

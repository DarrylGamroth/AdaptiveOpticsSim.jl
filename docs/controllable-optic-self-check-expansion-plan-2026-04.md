# Controllable Optic Self-Check Expansion Plan 2026-04

Status: active

Plan traceability:

- [`CO-6`](./controllable-optics-plan-2026-04.md)
- validity families: `MV-05`, `MV-06`, `MV-07`, `MV-14`

## Purpose

This plan extends the implemented controllable-optic self-check surface beyond
the current static `ShackHartmann` composite-optic contract.

The current maintained evidence is strong for:

- interface correctness
- low-order directional behavior
- composition and ordering
- small-signal OPD linearity
- deterministic CPU artifacts
- CUDA and AMDGPU parity

But it is still narrow in three important ways:

1. it uses one primary sensing family for the maintained self-check surface
2. it validates one command amplitude regime more strongly than a broader
   response envelope
3. it is still light on long-sequence statefulness and explicit failure-path
   contracts

This plan closes those gaps without blurring the distinction between:

- internal self-check evidence
- backend parity evidence
- external OOPAO / SPECULA equivalence

## Scope

This expansion plan covers six follow-on validation surfaces:

1. sensitivity-envelope checks
2. detector-in-the-loop invariants
3. multi-step statefulness checks
4. backend replay determinism
5. negative and failure-path contracts
6. multi-WFS and richer composite-family expansion

The goal is to keep `MV-14` scientifically useful even when an exact external
reference does not yet exist for a specific composite plant.

## Maintained WFS Matrix

The self-check work should no longer be treated as `ShackHartmann`-only.

The maintained rollout order is:

1. `ShackHartmann`
   - current reference surface
   - remains the first contract for every new composite optic
2. `Pyramid`
   - first non-SH sensing family for maintained low-order self-checks
3. `BioEdge`
   - second non-SH sensing family where comparable low-order observables exist
4. `Zernike` / `Curvature`
   - only where the observable contract is physically meaningful and stable

The intent is not to force every optic family through every WFS immediately.
The intent is to ensure that supported low-order controllable-optic claims are
backed by more than one sensing family.

## Validation Requirements

### OSX-1 Sensitivity Envelope

Each maintained low-order controllable optic must have a documented amplitude
envelope over which its response remains well-behaved on supported WFS
surfaces.

This includes:

- at least three amplitudes per optic on the maintained sensing surface
- amplitude-response monotonicity where that is physically expected
- clear separation between:
  - small-signal regime
  - valid but nonlinear regime
  - unsupported regime

### OSX-2 Detector-In-The-Loop Invariants

Each maintained runtime self-check surface must be revalidated with detector
transforms enabled, not just with the simplest exported pixel path.

This includes:

- non-null detector response with deterministic settings
- invariants that survive detector transforms:
  - command isolation
  - sign symmetry where meaningful
  - bounded energy / frame-sum changes
- explicit note where detector effects are expected to distort raw symmetry

### OSX-3 Multi-Step Statefulness

Composite runtime plants must remain correct over repeated command staging and
sensing steps.

This includes:

- repeated `set_command!` / `update_command!` / `sense!` sequences
- no stale state accumulation across steps
- no unintended mutation of untouched command segments
- deterministic replay of short multi-step sequences on CPU

### OSX-4 Backend Replay Determinism

Where deterministic behavior is claimed, each backend should be able to replay
the maintained self-check surface consistently across repeated runs.

This includes:

- repeated CPU replay of the artifact surface
- repeated CUDA replay where CUDA support is claimed
- repeated AMDGPU replay where AMDGPU support is claimed
- explicit tolerances for backend repeatability where exact replay is not
  realistic

### OSX-5 Negative And Failure-Path Contracts

The runtime control boundary must fail structurally and predictably on invalid
inputs.

This includes:

- wrong command length
- missing segment labels
- partial updates to nonexistent segments
- unsupported grouped/export combinations
- clear exception behavior for invalid composite-command routing

### OSX-6 Multi-WFS And Richer Composite Families

The maintained self-check surface must expand beyond the current
`tiptilt/steering/focus + dm` `ShackHartmann` baseline.

This includes:

- low-order composite self-checks on:
  - `Pyramid`
  - `BioEdge`
- richer composite families such as:
  - `tiptilt + focus + dm`
  - `steering + focus + dm`
  - segmented high-order layouts when they become production-supported

## Action Plan

### OSX-A1 Add amplitude-sweep sensitivity checks

Requirement coverage:

- `OSX-1`

Deliverables:

- deterministic amplitude-sweep helpers for:
  - `tiptilt`
  - `steering`
  - `focus`
- per-optic amplitude-sweep assertions on the maintained WFS surfaces
- archived summary metrics for envelope boundaries where useful

Target locations:

- [test/testsets/control_and_runtime.jl](../test/testsets/control_and_runtime.jl)
- [scripts/generate_multi_optic_runtime_artifact.jl](../scripts/generate_multi_optic_runtime_artifact.jl)

Verification evidence:

- focused runtime test pass
- CPU artifact records the selected amplitude points and response summaries

Status:

- pending

### OSX-A2 Add detector-in-the-loop self-check variants

Requirement coverage:

- `OSX-2`

Deliverables:

- deterministic detector-backed variants of the current low-order runtime cases
- per-case invariant assertions that remain meaningful after detector response
- explicit contract note for invariants that are detector-sensitive

Target locations:

- [test/testsets/control_and_runtime.jl](../test/testsets/control_and_runtime.jl)
- [test/testsets/detectors_and_wfs.jl](../test/testsets/detectors_and_wfs.jl)
- [docs/platform-orchestration-validation.md](./platform-orchestration-validation.md)

Verification evidence:

- detector-backed CPU runtime test pass
- updated platform artifact or companion artifact records detector-backed
  observables

Status:

- pending

### OSX-A3 Add short-sequence statefulness contracts

Requirement coverage:

- `OSX-3`

Deliverables:

- deterministic multi-step command sequences for composite optics
- assertions for:
  - no stale state
  - no hidden accumulation
  - correct partial-update isolation across steps
- grouped-runtime variant where grouped composite plants are maintained

Target locations:

- [test/testsets/control_and_runtime.jl](../test/testsets/control_and_runtime.jl)
- [test/backend_optional_common.jl](../test/backend_optional_common.jl)

Verification evidence:

- CPU replay of short multi-step sequences
- optional GPU parity on the same sequence surface

Status:

- pending

### OSX-A4 Add backend repeatability checks

Requirement coverage:

- `OSX-4`

Deliverables:

- repeated-run wrappers for the maintained composite-optic GPU contract
- backend repeatability thresholds for:
  - command vectors
  - slopes
  - exported WFS frames
- explicit policy for exact vs tolerance-based replay per backend

Target locations:

- [scripts/gpu_runtime_equivalence_contract.jl](../scripts/gpu_runtime_equivalence_contract.jl)
- [docs/backend-validation-guide.md](./backend-validation-guide.md)

Verification evidence:

- repeated CUDA replay logs
- repeated AMDGPU replay logs
- documented tolerance rationale in backend docs

Status:

- pending

### OSX-A5 Add invalid-input and failure-path coverage

Requirement coverage:

- `OSX-5`

Deliverables:

- tests for invalid command layout interactions
- tests for invalid segment names and wrong command lengths
- explicit exception expectations for grouped/composite runtime misuse

Target locations:

- [test/testsets/control_and_runtime.jl](../test/testsets/control_and_runtime.jl)
- [src/Control/runtime/outputs.jl](../src/Control/runtime/outputs.jl)
- [src/Optics/controllable_optics.jl](../src/Optics/controllable_optics.jl)

Verification evidence:

- focused negative-path test pass
- no ambiguous silent acceptance of invalid commands

Status:

- pending

### OSX-A6 Expand self-checks across WFS families and richer composites

Requirement coverage:

- `OSX-6`

Deliverables:

- a maintained WFS-family matrix for low-order composite self-checks:
  - `ShackHartmann`
  - `Pyramid`
  - `BioEdge`
- a staged richer-composite rollout:
  - `tiptilt + focus + dm`
  - `steering + focus + dm`
- explicit defer notes for any WFS family where the low-order observable
  contract is not yet physically meaningful enough to freeze

Target locations:

- [test/testsets/control_and_runtime.jl](../test/testsets/control_and_runtime.jl)
- [scripts/generate_multi_optic_runtime_artifact.jl](../scripts/generate_multi_optic_runtime_artifact.jl)
- [docs/model-validity-matrix.md](./model-validity-matrix.md)
- [docs/platform-orchestration-validation.md](./platform-orchestration-validation.md)

Verification evidence:

- at least one non-SH maintained self-check artifact or explicit deterministic
  replay surface
- optional CUDA/AMDGPU parity for the newly maintained non-SH surface

Status:

- pending

## Recommended Execution Order

1. `OSX-A5`
2. `OSX-A3`
3. `OSX-A1`
4. `OSX-A2`
5. `OSX-A6`
6. `OSX-A4`

Rationale:

- failure-path and statefulness gaps are the most likely to catch real runtime
  bugs quickly
- amplitude sweeps should be added before broadening detector and WFS-family
  surfaces
- multi-WFS expansion should happen only after the base contracts are stable
- backend repeatability should be layered on top of stable per-surface
  semantics rather than used to discover them

## Traceability Matrix

| ID | Requirement | Planned actions | Primary evidence | Status |
| --- | --- | --- | --- | --- |
| `OSX-1` | Sensitivity envelope | `OSX-A1` | amplitude-sweep tests and artifact metrics | pending |
| `OSX-2` | Detector-in-the-loop invariants | `OSX-A2` | detector-backed runtime tests and artifact observables | pending |
| `OSX-3` | Multi-step statefulness | `OSX-A3` | deterministic multi-step replay tests | pending |
| `OSX-4` | Backend replay determinism | `OSX-A4` | repeated-run parity logs and tolerances | pending |
| `OSX-5` | Negative and failure-path contracts | `OSX-A5` | invalid-input tests and explicit exceptions | pending |
| `OSX-6` | Multi-WFS and richer composites | `OSX-A6` | non-SH self-check surface and richer composite coverage | pending |

## Done Definition

This plan is complete when:

- maintained low-order optics have amplitude-envelope evidence, not only
  single-point checks
- the self-check surface has at least one detector-backed maintained variant
- short multi-step composite command sequences replay correctly and without
  stale state
- backend repeatability is documented and checked for maintained GPU-claimed
  surfaces
- invalid composite/runtime command inputs fail predictably
- at least one non-`ShackHartmann` WFS family is part of the maintained
  controllable-optic self-check surface
- richer composite plants are either covered or explicitly deferred with reason

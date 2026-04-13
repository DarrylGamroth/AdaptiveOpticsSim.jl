# Controllable Optic Self-Check Plan 2026-04

Status: active

Plan traceability:

- [`CO-5`](./controllable-optics-plan-2026-04.md)
- [`CO-6`](./controllable-optics-plan-2026-04.md)
- validity family: `MV-14`

## Purpose

This plan defines the maintained self-checking validation surface for
controllable optics and composite runtime plants.

The goal is to make the new optic family scientifically checkable even when an
exact external OOPAO or SPECULA analogue does not yet exist for a given plant.

This plan covers:

- `TipTiltMirror`
- `FocusStage`
- `SteeringMirror`
- `DeformableMirror`
- `CompositeControllableOptic`
- runtime surfaces that expose one or more controllable optics through packed
  command layouts

This is not a replacement for external truth. It is the internal physics and
behavior contract that must hold before backend parity or cross-package
equivalence claims are trusted.

## Scope

The maintained self-check surface is organized into six validation layers:

1. interface tests
2. behavioral invariants
3. optical-effect tests
4. linearity and small-signal tests
5. order and composition tests
6. deterministic CPU truth artifacts plus GPU parity

Every new controllable optic or composite plant that becomes part of the
supported runtime surface should map to this structure.

## Validation Requirements

### OSV-1 Interface Correctness

The runtime-facing command interface for each controllable optic must be
explicit, deterministic, and segment-safe.

This includes:

- correct `command_layout(...)`
- correct segment labels and lengths
- `set_command!` writes all intended segments
- `update_command!` writes only targeted segments
- packed runtime command state matches optic command storage after staging

### OSV-2 Behavioral Invariants

Each low-order optic must demonstrate the expected directional response on a
deterministic static surface.

This includes:

- `+tip` and `-tip` produce opposite-sign response
- `+tilt` and `-tilt` produce opposite-sign response
- tip and tilt excite different dominant slope axes
- zero command returns the baseline response
- partial updates do not mutate untargeted segments

### OSV-3 Optical-Effect Semantics

Each optic must demonstrate the expected physical effect on the measured output.

This includes:

- tip/tilt behaves like a rigid pointing term on the sensing surface
- focus produces a symmetric focus-like perturbation
- steering behaves like a pointing term rather than a high-order distortion
- total signal/energy remains bounded as expected for rigid shift-like
  perturbations

### OSV-4 Small-Signal Linearity

For supported small-command regimes, optic responses must be approximately
linear around zero.

This includes:

- `response(a + b) ≈ response(a) + response(b)` within tolerance
- scaling a small command scales the response proportionally
- linearity checks are only required on explicitly documented small-signal
  surfaces

### OSV-5 Composition And Ordering

Composite optic plants must demonstrate deterministic ordering and expected
composition semantics.

This includes:

- `CompositeControllableOptic` matches the documented child ordering
- composite packed command routing matches child layouts
- sequential application and composite application agree on maintained surfaces
- non-targeted child optics remain unchanged under partial updates

### OSV-6 Deterministic Truth And Backend Parity

Supported controllable-optic validation surfaces must have deterministic CPU
truth and optional maintained GPU parity where the backend surface is claimed.

This includes:

- archived CPU truth artifact
- replay test for the truth artifact
- optional CUDA parity against CPU truth
- optional AMDGPU parity against CPU truth
- tolerances documented per surface rather than by global guesswork

## Action Plan

### OSV-A1 Build a maintained interface-contract test matrix

Requirement coverage:

- `OSV-1`

Deliverables:

- runtime/unit tests for:
  - `command_layout(...)`
  - segment labels
  - segment lengths
  - `set_command!`
  - `update_command!`
  - command-isolation on partial updates
- one shared helper for asserting packed-command integrity on:
  - single-optic surfaces
  - composite-optic surfaces
  - grouped runtime surfaces

Target locations:

- [test/testsets/control_and_runtime.jl](../test/testsets/control_and_runtime.jl)
- [src/Control/runtime/outputs.jl](../src/Control/runtime/outputs.jl)
- [src/Optics/controllable_optics.jl](../src/Optics/controllable_optics.jl)

Verification evidence:

- focused runtime test pass
- no stale command-segment mutation on untargeted optics

Status:

- partial

### OSV-A2 Add deterministic directional behavior tests for low-order optics

Requirement coverage:

- `OSV-2`

Deliverables:

- deterministic static-atmosphere behavior checks for:
  - pure tip
  - pure tilt
  - pure focus
  - pure steering
- explicit assertions for:
  - dominant response axis
  - opposite-sign antisymmetry
  - zero-command baseline
  - segment isolation

Target locations:

- [test/testsets/control_and_runtime.jl](../test/testsets/control_and_runtime.jl)
- optional artifact generator under [scripts/](../scripts/)

Verification evidence:

- exact or tolerance-based deterministic replay on CPU
- archived behavior summary for maintained surfaces

Status:

- partial

Notes:

- `tiptilt + dm` directional checks are already present and should become the
  template for the remaining low-order optics.

### OSV-A3 Add optical-effect observables for each maintained optic family

Requirement coverage:

- `OSV-3`

Deliverables:

- per-optic effect observables such as:
  - slope-axis dominance
  - centroid or frame-shift statistics
  - signal-energy preservation for rigid shifts
  - radial symmetry metrics for focus-like surfaces
- helper functions that summarize the observable contract instead of forcing
  every test to compare full frames only

Target locations:

- [scripts/generate_multi_optic_runtime_artifact.jl](../scripts/generate_multi_optic_runtime_artifact.jl)
- new optic-behavior generators if needed
- [docs/platform-orchestration-validation.md](./platform-orchestration-validation.md)

Verification evidence:

- archived CPU artifact records the observables
- tests assert the observable contract directly

Status:

- pending

### OSV-A4 Add small-signal linearity checks

Requirement coverage:

- `OSV-4`

Deliverables:

- small-signal checks for:
  - tip/tilt
  - focus
  - steering
- documented tolerance envelopes for:
  - additivity
  - scalar response scaling

Target locations:

- [test/testsets/control_and_runtime.jl](../test/testsets/control_and_runtime.jl)
- [docs/model-validity-matrix.md](./model-validity-matrix.md)

Verification evidence:

- deterministic CPU checks at one or more documented amplitudes
- validity note describing where linearity is expected and where it is not

Status:

- pending

### OSV-A5 Add explicit composition and ordering checks

Requirement coverage:

- `OSV-5`

Deliverables:

- tests that compare:
  - sequential application of child optics
  - `CompositeControllableOptic` application
- explicit order-sensitive checks for:
  - `tiptilt + dm`
  - `steering + dm`
  - `focus + dm`
- grouped runtime checks for nested command layouts over composite optics

Target locations:

- [test/testsets/control_and_runtime.jl](../test/testsets/control_and_runtime.jl)
- [docs/controllable-optics-plan-2026-04.md](./controllable-optics-plan-2026-04.md)

Verification evidence:

- sequential-vs-composite equivalence on documented surfaces
- documented ordering semantics in the active architecture docs

Status:

- pending

### OSV-A6 Freeze deterministic truth artifacts and GPU parity contracts

Requirement coverage:

- `OSV-6`

Deliverables:

- one deterministic CPU artifact per maintained composite-optic validation
  surface
- corresponding replay tests
- optional CUDA and AMDGPU parity contracts with explicit tolerances
- manifest updates for archived artifacts

Target locations:

- [benchmarks/results/platform/](../benchmarks/results/platform/)
- [scripts/gpu_runtime_equivalence_contract.jl](../scripts/gpu_runtime_equivalence_contract.jl)
- [test/backend_optional_common.jl](../test/backend_optional_common.jl)
- [docs/backend-validation-guide.md](./backend-validation-guide.md)

Verification evidence:

- archived CPU truth artifact
- GPU parity runs on real hardware
- documented tolerance rationale per surface

Status:

- partial

## Recommended Execution Order

1. `OSV-A1`
2. `OSV-A2`
3. `OSV-A5`
4. `OSV-A3`
5. `OSV-A4`
6. `OSV-A6`

Rationale:

- interface correctness and directional behavior should be locked first
- composition semantics should be proven before broadening observables
- linearity is valuable, but it should be layered on top of a stable
  deterministic behavior surface
- artifact freezing and GPU parity should happen after the CPU meaning checks
  are stable

## Traceability Matrix

| ID | Requirement | Planned actions | Primary evidence | Status |
| --- | --- | --- | --- | --- |
| `OSV-1` | Interface correctness | `OSV-A1` | runtime tests | partial |
| `OSV-2` | Behavioral invariants | `OSV-A2` | deterministic behavior tests and artifact observables | partial |
| `OSV-3` | Optical-effect semantics | `OSV-A3` | optic-specific observables | pending |
| `OSV-4` | Small-signal linearity | `OSV-A4` | linearity tests | pending |
| `OSV-5` | Composition and ordering | `OSV-A5` | sequential-vs-composite checks | pending |
| `OSV-6` | Truth artifacts and GPU parity | `OSV-A6` | CPU artifact plus optional GPU contracts | partial |

## Done Definition

This plan is complete when:

- every maintained controllable-optic family has explicit interface coverage
- every maintained low-order optic has deterministic behavior checks
- composite ordering semantics are documented and tested
- at least one deterministic CPU artifact exists for each supported composite
  runtime surface
- optional CUDA and AMDGPU parity exist for every maintained GPU-claimed
  controllable-optic surface
- the validity matrix and backend-validation docs reflect the resulting
  supported scope

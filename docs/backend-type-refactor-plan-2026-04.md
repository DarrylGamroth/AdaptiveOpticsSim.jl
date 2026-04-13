# Backend Type Refactor Plan 2026-04

Status: implemented

## Purpose

This plan introduces explicit semantic backend parameters on the package's major public stateful types while preserving the existing concrete
array-type parameterization used for performance.

The goal is to make backend intent visible at the type level, propagate it
cleanly through model composition, and enforce a simple production rule:

- one simulation/runtime path uses one backend by default
- mixed CPU/GPU composition requires an explicit transfer boundary

This is a structural refactor. It should not be attempted as an opportunistic
cleanup inside unrelated feature work.

Implementation outcome:

- major public stateful types now carry explicit semantic backend parameters
- homogeneous backend composition is enforced for standard simulation/runtime paths
- public constructors accept semantic backend selectors, not raw array storage types
- `backend(x)`, `backend_type(x)`, `same_backend(...)`, and `require_same_backend(...)` are the maintained identity/query helpers
- maintained OOPAO equivalence remained within tolerance after the refactor
- the full package test suite remained green except for the pre-existing SPECULA reference regression

## Why This Is Needed

The current backend design is improved compared with the earlier raw-array API:

- public construction now uses semantic selectors such as `CPUBackend()` and
  `CUDABackend()`
- internal state remains array-parametric and GPU-capable
- runtime and HIL paths already validate CPU/CUDA/AMDGPU on maintained surfaces

But backend intent is still implicit in storage fields rather than explicit on
major public types like `Telescope`, `DeformableMirror`, `ShackHartmann`, and
`AOSimulation`.

This leaves a few weaknesses:

- backend compatibility is inferred late instead of declared early
- composed model paths can accidentally mix CPU and GPU objects
- dispatch on backend intent is harder than it needs to be
- diagnostics for mixed-backend misuse are weaker than they should be
- the semantic user API is better than the internal type story

## Design Principles

1. Keep the semantic selector API as the user-facing entry point.
- `CPUBackend()`, `CUDABackend()`, `AMDGPUBackend()`, and `MetalBackend()` stay
  the construction surface.
- Users should not be asked to think in terms of raw array container types.

2. Add explicit backend parameters to major public stateful types.
- Backend intent should be visible in type signatures where it materially
  affects model composition and runtime dispatch.

3. Retain concrete array-type parameterization for performance.
- The refactor should add backend identity, not replace concrete array typing.
- Concrete field types remain the primary mechanism for storage specialization.

4. Enforce homogeneous backend composition by default.
- A single simulation/runtime path should use one backend unless the user
  inserts an explicit transfer boundary.

5. Mixed-backend support must be explicit, not accidental.
- If CPU/GPU bridges are needed later, they should be represented by explicit
  staging/copy adapters rather than silent host-device movement.

6. Preserve hot-path behavior.
- No extra allocations in sensing/step hot paths.
- No new scalar indexing or hidden synchronization on GPU surfaces.

## Target Architecture

### Backend type identity

Major public types should gain an explicit backend parameter `B<:AbstractArrayBackend`
(or a backend-tag equivalent if that proves cleaner for dispatch).

Representative targets:

- `Telescope{B,...}`
- `DeformableMirror{B,...}`
- `LowOrderMirror{B,...}`
- `TipTiltMirror(...)` convenience constructor
- `FocusStage(...)` convenience constructor
- `CompositeControllableOptic{B,...}`
- `ShackHartmann{B,...}`
- `PyramidWFS{B,...}`
- `BioEdgeWFS{B,...}`
- `ZernikeWFS{B,...}`
- `CurvatureWFS{B,...}`
- `Detector{B,...}` where practical and non-disruptive
- `AOSimulation{B,...}`
- `ClosedLoopRuntime{B,...}`
- `RuntimeBranch{B,...}` and grouped/runtime boundary types where useful

Not every small helper type needs this immediately. The focus is the major
stateful objects that users compose into one path.

### Backend query helpers

Add simple public helpers:

- `backend(x)`
- `backend_type(x)` if needed
- `same_backend(x, y)`
- `require_same_backend(xs...)`

These should report semantic backend identity, not raw array container types.

### Homogeneous composition rule

The default rule for composed model paths should be:

- a `Telescope`, atmosphere, controllable optic, WFS, and detector inside one
  `AOSimulation` must share the same backend
- grouped runtime branches may differ across branches only if explicitly allowed
  by the scenario/config layer
- within one branch/path, mixed backends are rejected unless an explicit bridge
  object is used

### Explicit transfer boundary

Mixed CPU/GPU composition, if supported later, should use an explicit staging
object or copy adapter such as:

- `BackendTransferBoundary`
- `materialize_on(backend, x)`
- `copy_backend(backend, x)`

This is out of scope for the initial refactor, but the plan should leave room
for it.

## Phase Plan

### BT-1 Backend identity helpers and invariants

Deliverables:

- add `backend(x)` / `same_backend(...)` / `require_same_backend(...)`
- define the semantic backend identity rules for current objects
- add backend conformance tests for existing major types

Success criteria:

- backend identity is queryable without poking into array fields
- tests prove the helper API works on CPU/CUDA/AMDGPU maintained surfaces

### BT-2 Add backend parameter to core optical stateful types

Deliverables:

- add explicit backend parameter to:
  - `Telescope`
  - `DeformableMirror`
  - `LowOrderMirror`
  - `TipTiltMirror`
  - `FocusStage`
  - `CompositeControllableOptic`
- constructors enforce that allocated arrays match the declared backend

Success criteria:

- type display now shows backend identity
- no runtime behavior change on maintained CPU/CUDA/AMDGPU surfaces

### BT-3 Add backend parameter to sensing and detector path types

Deliverables:

- add explicit backend parameter to maintained WFS families
- add explicit backend parameter to detector types where practical
- thread backend identity through readout/export metadata when useful

Success criteria:

- WFS and detector composition no longer infer backend solely from field arrays
- maintained GPU parity surfaces remain green

### BT-4 Generalize `AOSimulation` and runtime types to explicit backend identity

Deliverables:

- `AOSimulation{B,...}`
- `ClosedLoopRuntime{B,...}`
- runtime/scenario helpers propagate backend identity explicitly

Success criteria:

- backend-homogeneous paths are enforced at build time or with clear structured
  errors during construction
- HIL/runtime code can dispatch on backend identity cleanly

### BT-5 Enforce homogeneous single-path composition

Deliverables:

- constructor checks for:
  - `AOSimulation(...)`
  - `RuntimeBranch(...)`
  - grouped runtime composition where backend homogeneity is required
- clear `InvalidConfiguration` errors for accidental CPU/GPU mixing

Success criteria:

- mixed backend paths fail fast with actionable errors
- no silent host/device copy behavior remains in standard path assembly

### BT-6 Compatibility and migration layer

Deliverables:

- compatibility constructors or aliases only where needed to avoid a hard break
- migration notes in docs and examples
- deprecation warnings for any helper patterns made obsolete by the refactor

Success criteria:

- active docs use the new type-explicit backend story
- migration burden is explicit and finite

### BT-7 Optional explicit transfer boundaries

Deliverables:

- only if needed after the main refactor
- design and prototype an explicit CPU/GPU bridge object or conversion helper

Success criteria:

- mixed backend use, if allowed at all, is explicit and auditable
- no ambiguity remains about ownership or copy semantics

## Validation Plan

### Required functional validation

- constructor smoke for all major stateful types on `CPUBackend()`
- optional backend parity on maintained CUDA/AMDGPU surfaces
- runtime/control path validation using the current maintained HIL/RTC surfaces
- grouped runtime validation where branch composition is touched

### Required external validation

- regenerate the maintained OOPAO equivalence artifact
- confirm the current production equivalence cases remain within tolerance
- rerun the maintained HEART runtime artifacts if runtime types are touched

### Required regression gates

- `Pkg.test(test_args=["control_and_runtime"])`
- optional backend test coverage in `test/backend_optional_common.jl`
- GPU smoke contract if constructor/runtime/backend identity surfaces are changed

## Migration Notes

Expected user-facing direction after the refactor:

```julia
tel = Telescope(...; backend=CPUBackend())
dm = DeformableMirror(tel; backend=CPUBackend())
wfs = ShackHartmann(tel; backend=CPUBackend())
```

And the composed path should enforce:

```julia
AOSimulation(tel, src, atm, dm, wfs)
```

only if all components share the same semantic backend.

The refactor should not reintroduce raw array types as public constructor
inputs. Those remain an internal implementation detail.

## Recommended Execution Order

1. `BT-1` backend helpers and invariants
2. `BT-2` core optical/public stateful types
3. `BT-4` `AOSimulation` and runtime types
4. `BT-5` homogeneous composition enforcement
5. `BT-3` sensing/detector type parameter propagation
6. `BT-6` migration cleanup
7. `BT-7` explicit transfer boundaries only if justified

This order keeps the highest-value semantics work first: backend identity and
path validity. It avoids overextending the initial slice into optional mixed
backend bridging before the package has a clean homogeneous model.

# Controllable Optics Plan 2026-04

Status: implemented

## Purpose

## Status Summary

This plan is implemented on the supported runtime surface. The runtime and controllable-optic API now uses the preferred names directly.

Completed phases:

- `CO-1` `AbstractControllableOptic` and shared command-layout interface
- `CO-2` `AOSimulation`/runtime storage generalized from `dm` to `optic`
- `CO-3` multi-optic apply/staging path through `sense!`, `step!`, and timing
- `CO-4` `CompositeControllableOptic` with ordered packed-command routing
- `CO-5` first extra controllable optics: `LowOrderMirror`, `TipTiltMirror`, `FocusStage`
- `CO-6` backend-native command routing covered on the maintained GPU surfaces
- `CO-7` public API/docs migrated to the preferred runtime vocabulary

This plan extends the runtime and plant model from a single deformable mirror to
an explicit family of controllable optics that can be driven through one packed
RTC command boundary while still applying one or more optical surfaces natively
on CPU, CUDA, and AMDGPU.

The immediate goal is to support cases such as:

- a deformable mirror mounted on a tip/tilt platform
- packed woofer/tweeter control surfaces
- steering-mirror plus DM paths
- linear-stage or focus-stage surfaces that perturb focus/defocus

This is a control/runtime and plant-application expansion. It is not yet a
broad passive-propagation expansion for dichroics, beam splitters, OAP/OAE, or
relay optics.

## Why This Is Needed

The current runtime cleanup completed the boundary-side work:

- explicit external-control support via `NullReconstructor`
- clearer runtime/scenario naming
- Julia-style accessors
- structured command layouts and nested branch commands

But the plant still has one hard-coded controllable optic:

- `AOSimulation` stores one `DeformableMirror`
- `ClosedLoopRuntime` stores one `dm`
- `sense!` and `step!` only apply one surface

So the package can now represent a packed command contract like
`(:tiptilt, :dm)` at the runtime boundary, but it cannot yet model that plant as
separate physical surfaces with separate application order.

## Design Principles

1. Separate controllable optics from passive propagation optics.
- Introduce `AbstractControllableOptic <: AbstractOpticalElement`.
- Keep passive optics such as dichroics, beam splitters, OAP/OAE, and relay
  optics under the broader propagation/optical-element family for later work.

2. Preserve the packed RTC command boundary.
- Controllers and HIL systems often use one flat command vector.
- The new optic layer should enrich metadata and application semantics without
  breaking flat-vector interoperability.

3. Keep hot paths allocation-free and GPU-native.
- Command staging must work on backend-native arrays.
- Applying one or more optics in a step must not require host copies.
- Multi-surface support must compose with existing CUDA/AMDGPU validation.

4. Prefer explicit interfaces plus traits over ad hoc field inspection.
- Follow the detector/runtime pattern: capability queries and clear surface
  contracts.

5. Keep compatibility during rollout.
- `DeformableMirror` remains supported as the default/single-optic path while
  the composite optic layer is introduced.

## Target Architecture

### Core abstractions

Add:

- `AbstractControllableOptic <: AbstractOpticalElement`
- `AbstractControllableSurface` only if later needed for finer-grained
  composition; do not introduce it in the first slice unless required

Define a small runtime-facing interface:

- `command_layout(optic)`
- `set_command!(optic, command)`
- `apply!(optic, tel, mode)`
- `n_control_dofs(optic)`
- `controllable_surface_labels(optic)`
- `supports_segmented_command(optic)`

The minimum requirement is that a controllable optic owns or views one packed
command vector and can apply its current state to the telescope without
allocations.

### Concrete optics

Required first:

- `DeformableMirror <: AbstractControllableOptic`
- `CompositeControllableOptic`

Planned next:

- `LowOrderMirror`
- `TipTiltMirror`
- `FocusStage` or `LinearStageFocus` for focus/defocus-like stage perturbation

The first goal is not a perfect mechanical model of every stage. It is to
create a consistent plant/runtime abstraction that can express and apply those
surfaces cleanly.

### Composite optic semantics

`CompositeControllableOptic` should:

- own an ordered tuple of child controllable optics
- expose one packed `command_layout`
- route command segments to child optics without copies where possible
- apply child optics in deterministic order during sensing and stepping

This gives a direct path to:

- `(:tiptilt, :dm)`
- `(:woofer, :tweeter)`
- `(:steering, :focus, :dm)`

## Phase Plan

### CO-1 Introduce `AbstractControllableOptic`

Deliverables:

- new abstract type
- baseline capability/interface functions
- `DeformableMirror` migrated to subtype it
- no behavior change yet

Success criteria:

- current single-DM simulations remain source-compatible
- runtime/control code can dispatch on controllable optics instead of concrete
  `DeformableMirror`

### CO-2 Generalize `AOSimulation` and runtime storage

Deliverables:

- `AOSimulation` parameterized on `AbstractControllableOptic`
- `ClosedLoopRuntime` stores `optic` instead of concrete `dm`
- compatibility accessors for the old `dm` field during transition if needed

Success criteria:

- existing DM-only paths still work
- no regression in current CPU/CUDA/AMDGPU maintained surfaces

### CO-3 Runtime apply pipeline for one or more optics

Deliverables:

- `apply_runtime_command!` generalized from one DM vector to one controllable
  optic command surface
- support for applying one or more optics in order in `sense!`, `step!`, and
  grouped runtime timing paths
- optional `apply_selected!` or equivalent internal helper so a step can update
  only some surfaces when required by future HIL use cases

Success criteria:

- a runtime can update multiple controllable optics in one step without manual
  slice math
- step timing paths still report phase costs correctly

### CO-4 `CompositeControllableOptic`

Deliverables:

- ordered composite optic type
- packed command-layout routing to children
- zero-allocation segmented `set_command!`
- deterministic ordered `apply!`

Success criteria:

- packed woofer/tweeter and tip/tilt+DM plants work through one runtime
- grouped scenarios can nest branch command layouts over composite optics

### CO-5 First additional concrete optics

Deliverables:

- `LowOrderMirror`
- `TipTiltMirror`
- one stage-like focus optic:
  - preferred first name: `FocusStage`
- `TipTiltMirror` should remain a convenience constructor over the general
  low-order modal optic rather than a distinct implementation

Notes:

- `FocusStage` should model a commanded focus-like optical perturbation, not a
  full opto-mechanical rail simulation
- the first implementation can use a compact modal or OPD basis so long as the
  command contract is explicit and validated

Success criteria:

- one maintained runtime test demonstrates `(:tiptilt, :dm)`
- one maintained runtime test demonstrates `(:focus, :dm)` or
  `(:steering, :dm)`

### CO-6 GPU-native command routing and application

Deliverables:

- no host copies in segmented `set_command!` for accelerator arrays
- child optic command views or direct copies remain backend-native
- timing/profile helpers include the multi-optic apply path
- optional backend tests cover CPU/CUDA/AMDGPU parity for composite optics on a
  maintained surface

Success criteria:

- packed composite-optic command injection is parity-tested on CPU vs CUDA and
  CPU vs AMDGPU
- no scalar indexing regressions on GPU

### CO-7 Public API migration and cleanup

Deliverables:

- cookbook/user-guide/runtime examples updated to use controllable-optic
  vocabulary where relevant
- `api-reference.md` updated
- old DM-specific runtime assumptions removed from the active runtime surface

Success criteria:

- HIL examples show multi-surface command layouts directly
- active runtime docs and examples use the direct controllable-optic vocabulary

## Validation Plan

### Functional validation

Add runtime tests for:

- single `DeformableMirror` still works unchanged
- `CompositeControllableOptic((tiptilt, dm))`
- `CompositeControllableOptic((woofer, tweeter))`
- grouped runtime with nested branch commands and composite optics
- `NullReconstructor` plus explicit segmented commands on composite optics

### Backend validation

Add optional backend parity tests for:

- CPU vs CUDA command injection and sensing outputs for composite optics
- CPU vs AMDGPU command injection and sensing outputs for composite optics

Prefer maintained surfaces that reuse current SH runtime/export validation
rather than inventing a new fragile benchmark surface.

### Performance validation

Track:

- segmented command injection cost
- multi-optic `apply!` cost
- total `sense!` / `step!` cost vs the current single-DM baseline

The new path must avoid unnecessary host transfers and should not materially
regress the current single-DM HIL surface.

## Non-Goals For This Plan

These are explicitly out of scope for the first implementation pass:

- full passive relay-optics modeling
- dichroic/beam-splitter path branching
- broad optical train graph execution
- full mechanical stage metrology models
- removing the compatibility layer immediately

Those may become later propagation/optical-graph work, but they should not be
mixed into the controllable-optic runtime expansion.

## Recommended Execution Order

1. `CO-1` abstract controllable-optic layer
2. `CO-2` generalize `AOSimulation` and runtime storage
3. `CO-3` generalized apply pipeline
4. `CO-4` composite controllable optic
5. `CO-6` GPU-native routing and parity validation
6. `CO-5` first extra optics (`TipTiltMirror`, `FocusStage`)
7. `CO-7` doc/API migration and deprecation prep

This order keeps the runtime and GPU semantics stable before adding more optic
families.

## Immediate Next Slice

Implement `CO-1` and `CO-2` first:

- introduce `AbstractControllableOptic`
- make `DeformableMirror` subtype it
- generalize `AOSimulation` and `ClosedLoopRuntime` away from concrete-DM-only
  storage
- keep the external API compatibility-preserving during that migration

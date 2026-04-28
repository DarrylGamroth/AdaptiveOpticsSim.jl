# Layered Interface Plan 2026-04

Status: active

## Purpose

This plan defines the next API/documentation cleanup pass for execution-layer
interfaces.

The goal is not to collapse all simulation verbs into one generic interface.
The goal is to make the maintained layers explicit, standardize the canonical
entry points within each layer, and stop mixing low-level and high-level
surfaces in the main user path.

## Problem Statement

The package currently has the right raw ingredients, but the public story is
still too mixed:

- primitive physics verbs:
  - `advance!`
  - `propagate!`
  - `apply!`
  - `measure!`
- runtime verbs:
  - `prepare!`
  - `sense!`
  - `step!`
  - `set_command!`
  - `readout`
- orchestration builders:
  - `build_control_loop_scenario`
  - `simulation_interface`

This is not a problem in the implementation itself. It is a problem in public
layering:

- the cookbook still mixes `AdaptiveOpticsSim.simulation_interface(runtime)` with
  `build_control_loop_scenario(...)`
- the user path does not yet treat one runtime/orchestration surface as the
  clear default
- low-level and high-level APIs are both visible without enough guidance about
  when each should be used

## Design Decision

Adopt a standard interface per layer, not one universal verb family across all
layers.

### Layer 1: Primitive Physics Layer

This layer exists for:

- algorithm work
- calibration internals
- custom research workflows
- direct optical/sensor experimentation

Canonical verbs:

- `advance!`
- `propagate!`
- `apply!`
- `measure!`
- `capture!`

Policy:

- keep these verbs explicit
- do not rename them to match runtime verbs
- document them as lower-level surfaces

### Layer 2: Runtime Execution Layer

This layer exists for:

- repeated AO simulation execution
- HIL plant execution
- closed-loop execution with package-owned reconstructors
- runtime timing and export surfaces

Canonical interface:

- `prepare!`
- `sense!`
- `step!`
- `set_command!`
- `update_command!`
- `readout`
- `command`
- `slopes`
- `wfs_frame`
- `science_frame`
- `wfs_metadata`
- `science_metadata`

Policy:

- this is the maintained execution surface for normal AO/HIL users
- runtime docs should prefer this layer over direct primitive-verb sequences
- HIL examples should use `set_command!` + `sense!`
- internal closed-loop examples should use `step!`

### Layer 3: Orchestration Layer

This layer exists for:

- named runtime composition
- grouped or multi-branch execution
- explicit export/product policy
- stable script-first scenario construction

Canonical interface:

- `ControlLoopBranch`
- `SingleControlLoopConfig`
- `GroupedControlLoopConfig`
- `ControlLoopScenario`
- `build_control_loop_scenario`
- then the same execution verbs as Layer 2:
  - `prepare!`
  - `sense!`
  - `step!`
  - `readout`

Policy:

- this is the preferred public surface for cookbook/runtime assembly examples
- single-runtime cookbook examples should prefer `build_control_loop_scenario(...)`
  unless the example is explicitly teaching low-level runtime assembly

## Scope Boundary

This pass is about public layering, naming, and documentation.

It is not about:

- replacing primitive verbs with runtime verbs
- redesigning numerical algorithms
- changing the runtime dataflow semantics
- introducing a ModelingToolkit-style DSL

## Target Public Story

The package should present the execution stack like this:

1. Component construction
   - telescope
   - source
   - atmosphere
   - controllable optic
   - sensor / detector
2. Plant assembly
   - `AOSimulation(...)`
3. Runtime assembly
   - `ControlLoopBranch(...)`
   - `SingleControlLoopConfig(...)` or `GroupedControlLoopConfig(...)`
   - `build_control_loop_scenario(...)`
4. Execution
   - `prepare!(...)`
   - `sense!(...)` for external control
   - `step!(...)` for internal control
5. Boundary readout
   - `readout(...)`
   - `command(...)`
   - `slopes(...)`
   - `wfs_frame(...)`
   - `science_frame(...)`

The low-level `ClosedLoopRuntime` + `AdaptiveOpticsSim.simulation_interface(...)` path remains
supported, but it should be documented as a lower-level/power-user surface.

## Concrete Work Items

### LI-1 Documentation Layer Split

Update stable docs so they clearly map to one layer at a time.

Required changes:

- `README.md`
  - keep the user-facing runtime story on the orchestration/runtime surface
- `docs/user-guide.md`
  - explicitly describe Layers 1-3
- `docs/model-cookbook.md`
  - keep normal runtime recipes on `build_control_loop_scenario(...)`
  - label any `AdaptiveOpticsSim.simulation_interface(...)` example as advanced/low-level
- `docs/api-reference.md`
  - group exported execution APIs by layer
- `docs/control-loop-orchestration.md`
  - explicitly state that this is the preferred public runtime assembly layer

Success criteria:

- a normal user can stay on the orchestration/runtime path without seeing the
  low-level runtime wrapper first
- the docs do not imply that `AdaptiveOpticsSim.simulation_interface(...)` and
  `build_control_loop_scenario(...)` are equally preferred entry points

### LI-2 Cookbook Normalization

Normalize the cookbook so it follows one public execution story.

Required changes:

- convert the current closed-loop runtime recipe to the scenario builder path,
  or relabel it explicitly as a low-level runtime recipe
- keep one dedicated advanced recipe for direct `ClosedLoopRuntime` /
  `AdaptiveOpticsSim.simulation_interface(...)`
- make the HIL recipes use:
  - `prepare!`
  - `set_command!`
  - `sense!`
  - `readout`

Success criteria:

- cookbook recipes do not jump layers without saying so
- HIL and internal closed-loop examples are obviously distinct

### LI-3 API Reference Layer Index

Make the reference page reflect the layered API model.

Required changes:

- add explicit subsections:
  - primitive physics layer
  - runtime execution layer
  - orchestration layer
- document when `AdaptiveOpticsSim.simulation_interface(...)` is appropriate
- document the exact `sense!` vs `step!` distinction once and reuse it

Success criteria:

- the reference reads as a coherent interface hierarchy, not just a symbol list

### LI-4 Advanced Surface Labeling

Treat `AdaptiveOpticsSim.simulation_interface(...)` as an advanced surface without removing it.

Required changes:

- document it as:
  - low-level
  - single-runtime
  - useful when manually assembling or testing runtimes
- keep it in API reference and maintainer docs
- avoid making it the first runtime example in user-facing docs

Success criteria:

- users encountering `AdaptiveOpticsSim.simulation_interface(...)` understand why it exists
- users do not need it for the normal HIL/runtime path

### LI-5 Example Classification

Tag or group examples by layer.

Required changes:

- user-facing examples:
  - scenario/orchestration path
- advanced examples:
  - direct primitive layer
  - direct runtime-interface layer

Success criteria:

- the examples tree matches the documented public story

## Optional Follow-On

If the layered story works well after the doc/examples pass, consider a small
API cleanup:

- keep `AdaptiveOpticsSim.simulation_interface(...)` exported if it remains a legitimate advanced
  surface
- otherwise move it to the advanced docs only and consider reducing its
  prominence in exported examples

This follow-on should only happen after the documentation split is complete and
validated.

## Execution Order

1. `LI-1` documentation layer split
2. `LI-2` cookbook normalization
3. `LI-3` API reference layer index
4. `LI-4` advanced surface labeling
5. `LI-5` example classification

## Validation

This plan is complete when:

- the stable user path (`README.md`, `user-guide.md`, `model-cookbook.md`)
  stays on the orchestration/runtime surface by default
- `api-reference.md` explicitly groups the interfaces by layer
- `control-loop-orchestration.md` is the clear first stop for public runtime
  assembly
- at least one advanced doc/example still covers direct
  `AdaptiveOpticsSim.simulation_interface(...)` usage for power users

## Relationship To Existing Plans

This plan complements, but does not replace:

- [`interface-style-spec-plan.md`](./interface-style-spec-plan.md)
  - broad interface/naming contract work across abstract families
- [`archive/2026-04/api-cleanup-plan-2026-04.md`](./archive/2026-04/api-cleanup-plan-2026-04.md)
  - completed runtime/HIL naming and accessors cleanup record
- [`control-loop-orchestration.md`](./control-loop-orchestration.md)
  - maintained orchestration-layer guide

This document is the narrower execution-layer cleanup plan that turns those
general decisions into a coherent public layering story.

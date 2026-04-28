# API Cleanup Plan 2026-04

Status: completed naming cleanup

## Purpose

This plan records the control/runtime API naming cleanup. The package now
prefers names that describe the computational role of each object rather than
historical "platform" or "product" terminology.

The breaking cleanup replaced the historical names:

- `RuntimeBranch`
- `SingleRuntimeConfig`
- `GroupedRuntimeConfig`
- `RuntimeScenario`
- `RuntimeProductRequirements`
- `GroupedRuntimeProductRequirements`
- `simulation_*` accessors

with clearer user-facing vocabulary. No compatibility aliases are kept.

## Problems To Address

1. HIL assembly currently requires a fake reconstructor.
2. Orchestration names are too internal and do not read naturally to new users.
3. `simulation_*` accessors are verbose and un-Julian.
4. Multi-DM control is still represented as one opaque flat command vector.
5. The current docs explain the old surface more strongly than the intended one.

## Compatibility Strategy

The package intentionally does not keep compatibility aliases for the old
runtime orchestration names. New code should use the names below.

## Phase Plan

### AC-1 External-Control Runtime Support

Add an explicit no-op control operator so HIL runtimes do not need a fake
reconstructor.

Deliverables:

- `NullReconstructor`
- clear `step!` / `reconstruct!` failure mode for that operator
- HIL docs using `sense!` + `set_command!`
- tests pinning external-command runtime behavior

Status:

- implemented in the current pass

### AC-2 Control-Loop Orchestration Naming

Current maintained names:

- `ControlLoopBranch`
- `SingleControlLoopConfig`
- `GroupedControlLoopConfig`
- `ControlLoopScenario`
- `build_control_loop_scenario`
- `RuntimeOutputRequirements`
- `GroupedRuntimeOutputRequirements`

Status:

- implemented as the only public surface

### AC-3 Preferred Accessor Layer

Add shorter Julia-style accessors as the preferred public surface.

Target preferred names:

- `readout(x)`
- `command(x)`
- `slopes(x)`
- `wfs_frame(x)`
- `science_frame(x)`
- `wfs_metadata(x)`
- `science_metadata(x)`
- `grouped_wfs_stack(x)`
- `grouped_science_stack(x)`

Status:

- implemented in the current pass

### AC-4 Structured Command Layouts For Multi-DM Control

Replace the current “opaque flat vector only” story with an explicit command
layout contract.

Target additions:

- command-layout metadata for runtime boundaries
- named command segments for grouped and multi-DM paths
- `set_command!` overloads for structured command containers such as
  `NamedTuple`
- retained flat-vector support for RTC compatibility

Success criteria:

- woofer/tweeter configurations can be addressed without manual slice math in
  user code
- the exported boundary still records the packed flat command order

Status:

- implemented in the current pass
- runtime boundaries can now carry explicit segmented command layouts
- grouped runtime scenarios accept nested structured commands by branch
- current structured-layout support is label-and-slice based, which is enough for packed woofer/tweeter, tip/tilt, or steering-mirror surfaces that share one command vector
- future work can add richer surface-specific metadata without changing the packed RTC boundary

### AC-5 Documentation And Example Migration

Switch user-facing docs and maintained examples to the preferred surface.

Scope:

- `README.md`
- `docs/user-guide.md`
- `docs/model-cookbook.md`
- runtime/platform examples

Status:

- started in the current pass

### AC-6 Deprecation Layer

After the preferred surface is stable and validated:

- add deprecations for the old orchestration names
- add deprecations for `simulation_*` accessors in user-facing docs and API
  guidance
- keep the compatibility layer long enough to survive one validated release

Status:

- completed

### AC-7 Compatibility Removal

Remove the legacy layer only after:

- docs and examples use the preferred surface
- release validation passes on CPU, CUDA, and AMDGPU
- cross-package maintained artifacts do not depend on the old names

Status:

- completed

## Current Preferred Surface

For new code, prefer:

```julia
branch = ControlLoopBranch(:main, sim, NullReconstructor(); science_detector=det)
cfg = SingleControlLoopConfig(name=:hil_demo, branch_label=:main,
    outputs=RuntimeOutputRequirements(slopes=true, wfs_pixels=true, science_pixels=true))
scenario = build_control_loop_scenario(cfg, branch)
prepare!(scenario)
set_command!(scenario, cmd)
sense!(scenario)
img = science_frame(scenario)
meta = science_metadata(scenario)
```

## Next Implementation Target

The next highest-value slice after `AC-4` is the plant/runtime expansion in
[`controllable-optics-plan-2026-04.md`](./controllable-optics-plan-2026-04.md),
and the migration/removal slice is now complete for the legacy `dm` and old orchestration surfaces.

# Optional Integration Boundaries

Date: 2026-04-03

Status: active

Plan traceability:

- [`PLAN-16`](./package-review-action-plan.md)
- [`PLAN-40`](./package-review-action-plan.md)
- review IDs: `PR-07`, `PR-08`, `PR-38`

## Purpose

This document records the current package-boundary decision after the Phase 2
modularization pass.

The goal is to keep the core package focused on maintained AO simulation while
keeping optional backend and I/O integrations outside the default core code
path.

## Decision Summary

No new package split is required in Phase 2.

The correct boundary today is:

- keep the scientific simulation core in `AdaptiveOpticsSim.jl`
- keep backend-specific and format-specific integrations in `ext/`
- keep future science-path integrations in optional adapter packages or
  extensions rather than moving them into the core package

This remains true after the package-review cleanup phases. Phase 7 strengthens
the policy rather than changing it.

## Keep In Core

These surfaces remain part of the core package because they are part of the
maintained simulation model and runtime:

- optics and field propagation
- atmosphere models
- detectors and detector response models
- WFS implementations
- calibration and reconstruction
- closed-loop runtime and timing surfaces
- tomography and atmospheric field propagation

## Keep In Extensions

These integrations are correctly placed in `ext/` and should remain optional:

- [`ext/AdaptiveOpticsSimCUDAExt.jl`](../ext/AdaptiveOpticsSimCUDAExt.jl)
- [`ext/AdaptiveOpticsSimAMDGPUExt.jl`](../ext/AdaptiveOpticsSimAMDGPUExt.jl)
- [`ext/AdaptiveOpticsSimMetalExt.jl`](../ext/AdaptiveOpticsSimMetalExt.jl)
- [`ext/AdaptiveOpticsSimTablesExt.jl`](../ext/AdaptiveOpticsSimTablesExt.jl)
- [`ext/AdaptiveOpticsSimCSVExt.jl`](../ext/AdaptiveOpticsSimCSVExt.jl)
- [`ext/AdaptiveOpticsSimJSON3Ext.jl`](../ext/AdaptiveOpticsSimJSON3Ext.jl)

Rationale:

- GPU policy and backend-native implementations should stay outside core model
  code where possible.
- Table and serialization support are useful, but they are not part of the AO
  model itself.

## Keep Out Of Core For Now

These areas should remain outside the core package boundary unless a later plan
justifies a stronger integration:

- `Proper.jl` science-path bridging
- downstream science-camera and coronagraph adapters
- cross-package benchmark harnesses and REVOLT-style scenario runners
- additional file-format/export adapters beyond the current lightweight
  extension surface

Recommended form:

- separate adapter package, or
- optional extension with a narrow maintained handoff API

## Phase 2 Outcome

Phase 2 does not move any existing optional integration out of the package.

Instead, it confirms:

- the current extension split is directionally correct
- further package splitting should wait until after more internal
  modularization and validation work
- future science-path integrations must remain optional by policy

## Phase 7 Outcome

Phase 7 keeps the same package boundary and makes the future-direction policy
explicit:

- OOPAO parity is not, by itself, a reason to move more features into core
- SPECULA-informed breadth work still has to respect the same package boundary
- science-path integration should arrive through a narrow handoff surface plus
  optional extensions or adapter packages

## Platform-Strengthening Phase 6 Outcome

Phase 6 keeps the same science-path boundary in force after the completed
platform-orchestration work.

The explicit rule remains:

- do not move focal-plane or downstream science-path work into core as the next
  default step
- keep future science-path work optional unless a later boundary review changes
  that decision explicitly

The current closeout decision is recorded in:

- [platform-strengthening-closeout.md](./platform-strengthening-closeout.md)

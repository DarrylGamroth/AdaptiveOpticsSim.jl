# Roadmap

Status: active

This roadmap is intentionally short. Historical implementation plans were
removed from the live docs tree during the April 2026 documentation cleanup; use
git history if old plan details are needed.

## Current State

AdaptiveOpticsSim.jl is usable for maintained CPU workflows and selected
hardware-validated GPU workflows. The core package now has:

- typed AO model objects for optics, atmosphere, WFS, detectors, DMs, control,
  and runtime orchestration
- CPU tests across Linux, macOS, and Windows
- optional CUDA and AMDGPU validation entry points for hosts with hardware
- a maintained zero-allocation hot-path focus for CPU HIL-style runtime loops
- committed reference data and regression tests for core WFS and detector paths
- a compact docs set with one extension guide instead of subsystem plan sprawl

## Near-Term Priorities

1. Keep CI and hardware validation healthy.
2. Keep model validity claims tied to tests, committed reference data, or
   benchmark artifacts.
3. Avoid expanding public API surface without updating
   [`api-reference.md`](api-reference.md) and [`extension-guide.md`](extension-guide.md).
4. Preserve allocation-free CPU HIL hot paths when changing runtime,
   detector/WFS, controller, or DM code.
5. Treat CUDA and AMDGPU support as maintained only on surfaces covered by the
   dedicated backend validation targets.

## Active Cleanup Themes

- Keep the README and `user-guide.md` focused on one recommended user path.
- Keep `api-reference.md` aligned with the exported API instead of documenting
  every internal or qualified helper.
- Keep full visual examples in `../AdaptiveOpticsSimPlots.jl`; the core
  package examples remain plotting-free and runnable through
  [`run_core_examples.sh`](../scripts/run_core_examples.sh).
- Prefer a few high-value OOPAO/SPECULA/REVOLT-like equivalence artifacts over
  broad claims that are not release-gated.
- Consolidate validation around maintained entry points:
  `Pkg.test()`, backend-specific hardware targets, the core example runner, and
  the release-validation script.

## Deferred Areas

The following remain valid future directions, but should not drive ad hoc API
growth:

- broader manufacturer-specific DM technology models
- richer detector physical models where reusable readout/thermal/counting layers
  are insufficient
- wider cross-package numerical equivalence beyond the maintained reference data
  and artifacts
- companion visualization and analysis packages outside the core package
- science-path integrations that belong in optional extensions or sibling
  packages

## Documentation Rule

Do not add new one-off roadmap fragments. Update this file, the
[`documentation-map.md`](documentation-map.md), or the relevant maintained guide
instead.

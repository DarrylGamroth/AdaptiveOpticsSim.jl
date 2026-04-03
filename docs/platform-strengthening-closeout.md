# Platform Strengthening Closeout

Date: 2026-04-03

Status: active

Plan traceability:

- [`PSP-18`](./platform-strengthening-plan.md)
- validity direction: `PSR-06`

## Purpose

This note closes the platform-strengthening plan after the orchestration pass.

It answers two questions explicitly:

1. what still remains behind SPECULA after the completed platform work
2. what the next main-platform step should be

## What Improved In This Plan

The completed plan materially improved:

- synthesis-oriented user/platform documentation
- validation and benchmark evidence
- ROCm/AMDGPU execution transparency
- typed platform orchestration
- benchmark-backed platform runtime evidence on maintained backends
- one normalized cross-package platform contract against the neighboring Julia
  legacy tree

The package is now in a stronger position as a Julia-native AO platform, not
just as a collection of subsystem implementations.

## What Still Trails SPECULA

The remaining gaps are no longer basic runtime ownership or grouped export
semantics. They are broader platform-breadth topics.

### 1. Controller And Process Breadth

SPECULA still has the stronger surface for:

- process-object style orchestration breadth
- controller/process family variety
- larger platform workflows built around those process families

AdaptiveOpticsSim should not chase breadth here casually. The remaining gap is
real, but it should only be closed through a focused, benchmark-backed plan.

### 2. Platform-Scale External Comparability

The package now has:

- SPECULA-informed platform/runtime evidence
- a normalized Julia-to-Julia grouped platform contract

It does not yet have:

- a stronger executable external SPECULA platform benchmark ladder
- a broadly normalized external platform runtime comparison surface

That gap is now smaller and better documented, but it still exists.

### 3. Optional Science-Path Breadth

SPECULA-style downstream science-path breadth remains beyond the current core
scope.

That is an intentional boundary, not a hidden backlog item:

- core should expose stable handoff surfaces
- science-path / focal-plane work remains optional

## Explicit Non-Goals After This Plan

The completed plan does **not** justify moving immediately to:

- config-file-first orchestration
- broad science-path integration in core
- controller/process-family expansion without a benchmarked target use case
- parity work justified only because SPECULA has a broader platform surface

## Next Main-Platform Step

The next main-platform step should be:

- a focused controller/process-breadth decision and implementation plan

That next plan should only start once it names:

- a concrete target use case
- the process/controller family to add
- the benchmark and validation surfaces that make the work worth doing

In other words:

- the next step is not more platform cleanup
- the next step is not config manifests
- the next step is not science-path expansion in core
- the next step is a deliberately scoped platform-breadth milestone

## Decision Summary

After this closeout:

- OOPAO remains a behavior/reference baseline
- SPECULA remains the stronger external reference for future platform breadth
- science-path work remains optional
- manifest/config-first orchestration remains deferred
- the next serious expansion should be a focused controller/process-breadth
  plan, not another generic cleanup pass

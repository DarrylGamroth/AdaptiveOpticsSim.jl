# Control-Loop Orchestration Validation

Status: active

Plan traceability:

- [`PSP-15`](./platform-strengthening-plan.md)
- [`PSP-17`](./platform-strengthening-plan.md)
- validity family: `MV-13`

## Purpose

This note records the maintained evidence for the typed control-loop orchestration
layer.

The objective is to make three things explicit:

- the direct `ControlLoopScenario` runtime surface has realistic timing evidence
- CPU and maintained GPU backends exercise that surface intentionally
- the orchestration layer is benchmarked as a first-class API, not only through
  older runtime wrappers

## Maintained Runner

- [profile_control_loop_runtime.jl](../scripts/profile_control_loop_runtime.jl)

Default commands:

```bash
julia --project=. --startup-file=no scripts/profile_control_loop_runtime.jl cpu medium
julia --project=. --startup-file=no scripts/profile_control_loop_runtime.jl amdgpu medium
```

CUDA is maintained on the `spiders` host:

```bash
julia --project=. --startup-file=no scripts/profile_control_loop_runtime.jl cuda medium
```

## Archived Artifact

- [2026-04-03-phase5-psp15.toml](../benchmarks/results/platform/2026-04-03-phase5-psp15.toml)
- [2026-04-13-multi-optic-hil.toml](../benchmarks/results/platform/2026-04-13-multi-optic-hil.toml)
- manifest:
  [manifest.toml](../benchmarks/results/platform/manifest.toml)

Generators:

- [generate_control_loop_orchestration_artifact.jl](../scripts/generate_control_loop_orchestration_artifact.jl)
- [generate_multi_optic_runtime_artifact.jl](../scripts/generate_multi_optic_runtime_artifact.jl)

## Backend Functional Coverage

Optional backend smoke now instantiates the orchestration layer directly:

- [optional_amdgpu_backends.jl](../test/optional_amdgpu_backends.jl)
- [optional_cuda_backends.jl](../test/optional_cuda_backends.jl)
- shared logic in
  [backend_optional_common.jl](../test/backend_optional_common.jl)

The maintained smoke verifies:

- single-branch `ControlLoopScenario` execution
- grouped `ControlLoopScenario` execution
- grouped stack export availability on compatible grouped backends
- backend-specific execution-plan defaults still align with the orchestration
  surface

## Evidence Shape

The maintained control-loop artifacts record:

- single-branch control-loop runtime
  - detector-backed Pyramid branch with science output
- compatible grouped control-loop runtime
  - grouped WFS stack export enabled
- mixed grouped control-loop runtime
  - grouped WFS stack export explicitly unavailable
- deterministic multi-optic HIL runtime
  - `tiptilt + dm`, `steering + dm`, and `focus + dm` composite optics with
    structured command update
  - richer stateful composite sequences:
    - `tiptilt + focus + dm`
    - `steering + focus + dm`
  - exact replay of packed command, slopes, and WFS frame on the maintained CPU
    surface
  - maintained non-SH low-order replay surfaces:
    - `Pyramid`
    - `BioEdge`
  - behavioral observables:
    - directional slope-axis dominance for low-order pointing terms
    - opposite-sign antisymmetry for signed low-order commands
    - bounded frame-sum invariants for rigid-shift-like perturbations
    - detector-backed frame-response invariants under deterministic
      non-null detector response
  - composition and small-signal observables:
    - sequential-vs-composite OPD equivalence
    - small-signal OPD additivity and scaling for low-order optics
    - amplitude-sweep envelopes for maintained low-order response surfaces

For each case, the artifact records:

- build time
- step mean and p95
- average host-side allocation count observed from the runner
- export availability and shape

## Scope Limits

- This is a Julia-native orchestration/runtime artifact, not by itself a
  cross-package equivalence claim.
- AMDGPU still uses explicit host-mirror fallback for the maintained
  Poisson/readout-correction detector surfaces; that debt is tracked in
  [rocm-fallback-inventory.md](./rocm-fallback-inventory.md).
- Cross-package platform comparison is recorded separately in
  [revolt-platform-benchmark-contract.md](./revolt-platform-benchmark-contract.md).

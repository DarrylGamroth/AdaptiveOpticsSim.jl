Status: active

Plan traceability:

- [`PSP-07`](./platform-strengthening-plan.md)
- validity family: `MV-13`

# SPECULA-Informed Platform Runtime Validation

## Purpose

This note records the maintained platform-scale runtime artifact added to
strengthen SPECULA-informed evidence beyond isolated optical contracts.

AdaptiveOpticsSim remains Julia-native and does not claim full platform-level
numerical equivalence to SPECULA. The goal here is narrower and more useful:

- exercise grouped orchestration and heterogeneous multi-WFS runtime behavior
- archive the maintained exported-product behavior for those scenarios
- make the platform-scale claim explicit instead of leaving it implied by
  subsystem results

## Artifact

- [2026-04-02-phase2-psp07.toml](../benchmarks/results/platform/2026-04-02-phase2-psp07.toml)
- manifest:
  [manifest.toml](../benchmarks/results/platform/manifest.toml)

## Generator

- [generate_specula_platform_runtime_artifact.jl](../scripts/generate_specula_platform_runtime_artifact.jl)

Default command:

```bash
julia --project=. --startup-file=no scripts/generate_specula_platform_runtime_artifact.jl
```

## Why This Is SPECULA-Informed

The scenario is aligned to the class of platform compositions where SPECULA is
the stronger reference:

- multiple wfs/runtime branches
- integrated WFS + detector style pipelines
- orchestration behavior rather than only isolated kernel correctness

The current rationale is anchored to SPECULA processing-object tests such as:

- `test_mmse_reconstructor.py::TestMMSEReconstructor::test_mmse_multiple_wfs`
- `test_sprint.py`
- `test_sprint_pyr.py`

AdaptiveOpticsSim does not mirror those Python processing objects directly.
Instead, it measures the corresponding Julia-native grouped/composite runtime
surface and records the result as maintained evidence.

## Evidence Shape

The artifact records `medium` and `representative` CPU grouped-runtime cases
covering:

- per-family grouped asterism runtime:
  - Shack-Hartmann
  - Pyramid
  - BioEdge
- compatible grouped composite runtime
- mixed grouped composite runtime
- grouped export semantics:
  - compatible grouped WFS stack availability
  - mixed grouped stack refusal

## Scope Limits

- This is a Julia-native SPECULA-informed artifact, not a cross-package runtime
  equivalence benchmark.
- Optional AMDGPU/CUDA grouped runtime evidence remains in the maintained
  backend smoke/benchmark surfaces instead of this artifact.
- External SPECULA execution is still stronger for some process-family and
  orchestration breadth questions; this artifact narrows that gap but does not
  eliminate it.

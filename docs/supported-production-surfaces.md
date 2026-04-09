# Supported Production Surfaces

Status: active

## Purpose

This note defines the current production-supported scope for
AdaptiveOpticsSim.jl.

The intent is to be explicit about which surfaces are currently defended by
maintained evidence, backend parity checks, and realistic runtime artifacts,
and which surfaces should still be treated as experimental or research-grade.

This is a support-policy note, not a claim that every exported symbol or model
family has identical maturity.

## Support Rule

A surface should only be treated as production-supported when all of the
following are true:

1. the modeling/workflow surface is documented and maintained
2. functional correctness is covered in the normal test suite
3. backend parity is covered where accelerator support is claimed
4. realistic runtime evidence exists where performance or HIL claims are made
5. external equivalence is either demonstrated or explicitly scoped

Use these supporting docs together:

- [model-validity-matrix.md](./model-validity-matrix.md)
- [backend-validation-guide.md](./backend-validation-guide.md)
- [cross-package-benchmark-harness.md](./cross-package-benchmark-harness.md)
- [production-readiness-checklist.md](./production-readiness-checklist.md)

## Production-Supported Surfaces

### CPU baseline

The primary production baseline is the typed Julia CPU path on maintained
runtime and validation surfaces.

Current CPU-supported families:

- finite and infinite multilayer atmosphere
- detector-family execution on maintained detector surfaces
- Shack-Hartmann, Pyramid, and BioEdge WFS on maintained validated surfaces
- runtime and closed-loop execution on maintained local/runtime artifacts
- grouped/runtime orchestration surfaces with committed validation artifacts

Primary evidence:

- [model-validity-matrix.md](./model-validity-matrix.md)
- [platform-workflows.md](./platform-workflows.md)
- [benchmark-matrix-plan.md](./benchmark-matrix-plan.md)

### CUDA backend

CUDA is the strongest production-supported accelerator backend today.

Current CUDA-supported scope:

- maintained GPU smoke matrix surfaces
- maintained GPU runtime-equivalence surfaces
- maintained Shack-Hartmann exported-pixel parity surfaces
- matched HEART RTC HIL runtime surfaces

Primary evidence:

- [backend-validation-guide.md](./backend-validation-guide.md)
- [gpu-smoke-contract.jl](../scripts/gpu_smoke_contract.jl)
- [gpu_runtime_equivalence_contract.jl](../scripts/gpu_runtime_equivalence_contract.jl)
- [cross-package-benchmark-harness.md](./cross-package-benchmark-harness.md)

Current expectation:

- if a maintained CUDA surface regresses numerically against CPU, that is a
  release-blocking defect for the CUDA-supported scope

### AMDGPU backend

AMDGPU is now production-supported on the maintained surfaces that are covered
by the post-cleanup smoke, parity, and HEART runtime artifacts.

Current AMDGPU-supported scope:

- maintained GPU smoke matrix surfaces
- maintained GPU runtime-equivalence surfaces
- maintained Shack-Hartmann exported-pixel parity surfaces
- matched HEART RTC HIL runtime surfaces

Primary evidence:

- [backend-validation-guide.md](./backend-validation-guide.md)
- [rocm-failure-catalog.md](./rocm-failure-catalog.md)
- [rocm-fallback-inventory.md](./rocm-fallback-inventory.md)
- [cross-package-benchmark-harness.md](./cross-package-benchmark-harness.md)

Important qualifier:

- production support for AMDGPU is defined by the maintained validated
  surfaces, not by every possible ROCm kernel path
- unsupported or not-yet-recovered ROCm surfaces should remain explicit
  fallbacks or experimental paths

### OOPAO-aligned external equivalence

The strongest current external production-alignment story is the OOPAO-matched
HEART Shack-Hartmann surface.

Current production-supported external equivalence surface:

- HEART RTC HIL boundary:
  - `277` DM commands in
  - `352x352` Shack-Hartmann frame out
- matched OOPAO/Julia baseline on the comparison-owned flux-threshold validity
  convention

Primary evidence:

- [oopao-reference-datasets.md](./oopao-reference-datasets.md)
- [/home/dgamroth/workspaces/codex/AdaptiveOpticsComparisons/docs/cross-package-benchmark-harness.md](/home/dgamroth/workspaces/codex/AdaptiveOpticsComparisons/docs/cross-package-benchmark-harness.md)
- [/home/dgamroth/workspaces/codex/AdaptiveOpticsComparisons/results/archived/2026-04-08-heart-hil-oopao-baseline.toml](/home/dgamroth/workspaces/codex/AdaptiveOpticsComparisons/results/archived/2026-04-08-heart-hil-oopao-baseline.toml)

## Explicitly Not Yet Production-Supported

The following should still be treated as experimental, scoped, or unresolved:

- SPECULA pixel-level equivalence on the HEART Shack-Hartmann surface
- Metal backend support
- broad claims about every detector/WFS/backend combination outside the
  maintained evidence surfaces
- cross-package grouped/platform equivalence beyond the currently normalized
  contracts
- instrument-truth alignment beyond the currently maintained comparison
  artifacts

## Release Interpretation

For current release/readiness decisions:

- regressions on the CPU baseline are release-blocking
- regressions on the maintained CUDA surfaces are release-blocking
- regressions on the maintained AMDGPU surfaces are release-blocking
- unresolved SPECULA differences are not release-blocking unless the package
  starts claiming SPECULA equivalence for that surface
- new features outside this supported surface should ship as experimental until
  they gain the same evidence shape

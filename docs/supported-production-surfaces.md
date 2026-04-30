# Supported Production Surfaces

Status: active

## Purpose

This note defines the current production-supported scope for
AdaptiveOpticsSim.jl.

The intent is to be explicit about which surfaces are defended by maintained
evidence, backend parity checks, and real-hardware validation, and which
surfaces are outside the current support claim.

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

- [model-validity-matrix.md](model-validity-matrix.md)
- [backend-validation-guide.md](backend-validation-guide.md)
- [release-validation-runbook.md](release-validation-runbook.md)

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

- [model-validity-matrix.md](model-validity-matrix.md)
- the workflows in [`user-guide.md`](user-guide.md) and [`release-validation-runbook.md`](release-validation-runbook.md)
- benchmark artifacts under `benchmarks/results/`
- [../benchmarks/results/validation_runs/2026-04-10-rtc-devel-cpu.toml](../benchmarks/results/validation_runs/2026-04-10-rtc-devel-cpu.toml)

### CUDA backend

CUDA is a production-supported accelerator backend on the maintained surfaces
covered by the dedicated hardware validation target:

- [../test/runtests_cuda.jl](../test/runtests_cuda.jl)

Current CUDA-supported scope:

- maintained optional backend functional/parity checks
- maintained runtime-equivalence surfaces
- maintained high-accuracy post-command runtime-equivalence surfaces
- maintained Shack-Hartmann exported-pixel parity surfaces
- maintained composite low-order runtime surfaces
- matched HEART RTC HIL runtime surfaces

Primary evidence:

- [backend-validation-guide.md](backend-validation-guide.md)
- [release-validation-runbook.md](release-validation-runbook.md)
- benchmark artifacts under `benchmarks/results/`
- [../benchmarks/results/validation_runs/2026-04-10-spiders-cuda.toml](../benchmarks/results/validation_runs/2026-04-10-spiders-cuda.toml)

Current expectation:

- if a maintained CUDA surface regresses numerically against CPU, that is a
  release-blocking defect for the CUDA-supported scope

### AMDGPU backend

AMDGPU is a production-supported accelerator backend on the maintained surfaces
covered by the dedicated hardware validation target:

- [../test/runtests_amdgpu.jl](../test/runtests_amdgpu.jl)

Current AMDGPU-supported scope:

- maintained optional backend functional/parity checks
- maintained runtime-equivalence surfaces
- maintained high-accuracy post-command runtime-equivalence surfaces
- maintained Shack-Hartmann exported-pixel parity surfaces
- maintained composite low-order runtime surfaces
- matched HEART RTC HIL runtime surfaces

Primary evidence:

- [backend-validation-guide.md](backend-validation-guide.md)
- [release-validation-runbook.md](release-validation-runbook.md)
- benchmark artifacts under `benchmarks/results/`
- [../benchmarks/results/validation_runs/2026-04-10-rtc-devel-amdgpu.toml](../benchmarks/results/validation_runs/2026-04-10-rtc-devel-amdgpu.toml)

Current expectation:

- if a maintained AMDGPU surface regresses numerically against CPU, that is a
  release-blocking defect for the AMDGPU-supported scope

### GPU support-boundary rule

GPU support is defined by the maintained dedicated hardware targets and the
release-validation path.

The package does still contain broader backend-audit and subsystem-investigation
scripts, but those do not by themselves define supported GPU scope. A GPU-touched
surface is only support-claimed when it is also promoted into the maintained
hardware targets and release-validation cadence.

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

- committed OOPAO reference data under `test/reference_data`
- `../AdaptiveOpticsComparisons/docs/cross-package-benchmark-harness.md`
- `../AdaptiveOpticsComparisons/results/archived/2026-04-08-heart-hil-oopao-baseline.toml`

Additional production-supported frozen OOPAO equivalence surfaces:

- diffractive Pyramid ramp from the committed OOPAO reference bundle
- diffractive BioEdge ramp from the committed OOPAO reference bundle

Primary evidence:

- committed OOPAO reference data under `test/reference_data`
- [../scripts/generate_oopao_equivalence_artifact.jl](../scripts/generate_oopao_equivalence_artifact.jl)
- [../benchmarks/results/equivalence/2026-04-09-oopao-production-equivalence.toml](../benchmarks/results/equivalence/2026-04-09-oopao-production-equivalence.toml)

Scientist-owned HEART boundary truth artifact:

- [../scripts/generate_heart_boundary_truth_artifact.py](../scripts/generate_heart_boundary_truth_artifact.py)
- [../benchmarks/results/truth/2026-04-09-heart-boundary-truth.toml](../benchmarks/results/truth/2026-04-09-heart-boundary-truth.toml)

## Explicitly Not Yet Production-Supported

The following are outside the current support claim:

- SPECULA pixel-level equivalence on the HEART Shack-Hartmann surface
- Metal backend support
- backend-audit surfaces that are not part of the maintained hardware targets
  and release-validation cadence
- broad claims about every detector/wfs/backend combination outside the
  maintained evidence surfaces
- cross-package grouped/platform equivalence beyond the currently normalized
  contracts
- full optical or on-sky instrument-truth alignment beyond the maintained boundary artifact

## Release Interpretation

For current release/readiness decisions:

- regressions on the CPU baseline are release-blocking
- regressions on the maintained CUDA surfaces are release-blocking
- regressions on the maintained AMDGPU surfaces are release-blocking
- unresolved SPECULA differences are not release-blocking unless the package
  starts claiming SPECULA equivalence for that surface
- new features outside this supported surface should not be described as
  supported until they gain the same evidence shape
- supported accelerator claims require routine validation on real hardware,
  either through CI or through the maintained release-validation host path

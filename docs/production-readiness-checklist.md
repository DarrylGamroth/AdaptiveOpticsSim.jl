# Production Readiness Checklist

Status: active

## Purpose

This checklist tracks the concrete items required to treat
AdaptiveOpticsSim.jl as production-ready on its currently supported surface,
rather than only as a fast research port.

It is intentionally narrower than the full roadmap. The goal is to track the
remaining gaps between the current validated package and a mature, defensible
AO toolchain release surface.

Use together with:

- [supported-production-surfaces.md](supported-production-surfaces.md)
- [model-validity-matrix.md](model-validity-matrix.md)
- [backend-validation-guide.md](backend-validation-guide.md)
- benchmark artifacts under `benchmarks/results/`

## Current Status

### Validation breadth

- [x] Maintained CPU functional test suite
- [x] Maintained optional CUDA backend smoke/parity coverage
- [x] Maintained optional AMDGPU backend smoke/parity coverage
- [x] Maintained HEART OOPAO/Julia null/noise-free equivalence baseline
- [x] Maintained Julia CPU/CUDA/AMDGPU HEART runtime artifacts
- [x] Maintained second HEART realism runtime artifact (`null` vs `default`)
- [x] At least one additional trusted external equivalence surface beyond HEART SH
- [x] Instrument-truth artifact beyond package-to-package comparison

### Automated enforcement

- [x] Checked-in CPU CI workflow
- [x] Checked-in CUDA backend CI workflow
- [x] Checked-in AMDGPU backend CI workflow
- [x] CUDA validation routinely executed on real hardware (CI or release-validation host)
- [x] AMDGPU validation routinely executed on real hardware (CI or release-validation host)

### Supported-scope clarity

- [x] Production-supported scope explicitly documented
- [x] Explicit non-production/experimental scope documented
- [x] Public release notes or README-level summary of supported scope

### External truth and equivalence

- [x] Strong OOPAO-aligned HEART SH comparison surface
- [x] SPECULA mismatch resolved or explicitly closed out as non-blocking scope
- [x] Scientist-owned or instrument-owned truth artifact for HEART boundary

### Operational hardening

- [x] Maintained cross-package contracts and archived manifests
- [x] Same-host CPU-vs-CUDA runtime artifact on `spiders`
- [x] Versioned bootstrap/runner instructions for production validation hosts
- [x] One-command release validation procedure documented

## Required Before Claiming Release Readiness

There are no known code blockers on the currently supported production surface.
Before claiming a specific release or handoff commit as production-ready, rerun
and archive the release-validation evidence for that exact commit.

Historical operational evidence was recorded on commit `2192140`:

- [../benchmarks/results/validation_runs/2026-04-10-rtc-devel-cpu.toml](../benchmarks/results/validation_runs/2026-04-10-rtc-devel-cpu.toml)
- [../benchmarks/results/validation_runs/2026-04-10-rtc-devel-amdgpu.toml](../benchmarks/results/validation_runs/2026-04-10-rtc-devel-amdgpu.toml)
- [../benchmarks/results/validation_runs/2026-04-10-spiders-cuda.toml](../benchmarks/results/validation_runs/2026-04-10-spiders-cuda.toml)

Current automated CI evidence is provided by the CPU validation and coverage
workflows on `main`. CUDA and AMDGPU validation are hardware-gated manual
workflows and should be treated as release-validation evidence once they have
completed successfully on the target hardware.

## Next Recommended Order

1. Keep the release-validation runbook, cadence note, and bootstrap path exercised on production hosts.
2. Extend the same evidence shape to the next production-candidate surfaces rather than broadening support informally.

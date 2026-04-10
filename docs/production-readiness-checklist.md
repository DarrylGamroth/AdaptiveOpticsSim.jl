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

- [supported-production-surfaces.md](./supported-production-surfaces.md)
- [model-validity-matrix.md](./model-validity-matrix.md)
- [backend-validation-guide.md](./backend-validation-guide.md)
- [cross-package-benchmark-harness.md](./cross-package-benchmark-harness.md)

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
- [ ] CUDA validation routinely executed on real hardware (CI or release-validation host)
- [ ] AMDGPU validation routinely executed on real hardware (CI or release-validation host)

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

## Required Before Claiming Production Readiness

The remaining blocker is:

1. Routine GPU validation on real hardware must be operational, not just described.
   - That can be satisfied either by green self-hosted CI runners or by a maintained release-validation host path that is exercised for supported releases.

## Next Recommended Order

1. Put CUDA and AMDGPU validation on a real operational cadence, either through CI or through a maintained release-validation host.
2. Keep the release-validation runbook, cadence note, and bootstrap path exercised on production hosts.

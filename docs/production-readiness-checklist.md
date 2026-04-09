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
- [ ] At least one additional trusted external equivalence surface beyond HEART SH
- [ ] Instrument-truth artifact beyond package-to-package comparison

### Automated enforcement

- [x] Checked-in CPU CI workflow
- [x] Checked-in CUDA backend CI workflow
- [x] Checked-in AMDGPU backend CI workflow
- [ ] Self-hosted CUDA runner registered and green in routine CI
- [ ] Self-hosted AMDGPU runner registered and green in routine CI

### Supported-scope clarity

- [x] Production-supported scope explicitly documented
- [x] Explicit non-production/experimental scope documented
- [ ] Public release notes or README-level summary of supported scope

### External truth and equivalence

- [x] Strong OOPAO-aligned HEART SH comparison surface
- [ ] SPECULA mismatch resolved or explicitly closed out as non-blocking scope
- [ ] Scientist-owned or instrument-owned truth artifact for HEART boundary

### Operational hardening

- [x] Maintained cross-package contracts and archived manifests
- [x] Same-host CPU-vs-CUDA runtime artifact on `spiders`
- [ ] Versioned bootstrap/runner instructions for production validation hosts
- [x] One-command release validation procedure documented

## Required Before Claiming Production Readiness

The remaining blockers are:

1. GPU CI must be real, not just checked-in.
   - The workflows need actual self-hosted runners with the expected labels.
2. One more trusted equivalence surface is needed.
   - The current HEART SH surface is strong, but one surface is still too narrow.
3. Instrument-truth validation is still missing.
   - Package-to-package agreement is useful, but not sufficient as final truth.
4. SPECULA must be either resolved or explicitly scoped out.
   - Unresolved exploratory comparison should not stay ambiguous.
5. Release operations need a single reproducible validation recipe.

## Next Recommended Order

1. Register and exercise self-hosted CUDA and AMDGPU runners in CI.
2. Add one more trusted equivalence/runtime surface beyond HEART SH.
3. Secure a HEART instrument-truth artifact or scientist-owned truth bundle.
4. Write the release validation runbook.
5. Decide whether SPECULA stays exploratory or becomes a formal target.

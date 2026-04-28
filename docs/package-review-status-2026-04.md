# Package Review Status Update: 2026-04

Status: active

This note updates the March 2026 package review after the production-hardening,
external-equivalence, and GPU backend cleanup work completed in early April.

Use together with:

- [package-review-2026-03.md](./package-review-2026-03.md)
- [production-readiness-checklist.md](./production-readiness-checklist.md)
- [supported-production-surfaces.md](./supported-production-surfaces.md)
- [documentation-map.md](./documentation-map.md)

## Executive Update

The March review is still directionally correct. The package still looks like a
strong Julia codebase that needs continued architectural refactoring rather
than a rewrite in another language.

The biggest changes since that review are:

- production-readiness evidence is materially stronger
- OOPAO-aligned equivalence is now strong on maintained surfaces
- CUDA and AMDGPU parity and runtime validation are much stronger
- documentation/navigation is better curated than it was in late March

So the review remains relevant, but some findings should now be read as
partially mitigated rather than current top blockers.

## Findings Still Correct

### Structural diagnosis

The central diagnosis still holds: the main remaining weaknesses are structural,
not language-level. The package is already aligned with Julia's strengths in:

- backend-parametric arrays and dispatch in [src/core/backends.jl](../src/core/backends.jl)
- a broad AO surface in [src/AdaptiveOpticsSim.jl](../src/AdaptiveOpticsSim.jl)
- mutating hot paths and zero-allocation runtime evidence in [performance-audit.md](./performance-audit.md)

That still argues for continued refactoring of the current architecture rather
than a rewrite in a systems language.

### API breadth and oversized files

These are still real issues. The public surface is still broad, and several
source/test files are still larger than they should be for a mature package.
Those remain good refactoring targets.

### SPECULA breadth gap

The March review was also right that SPECULA still leads in some controller and
process breadth. That is now explicitly scoped as non-blocking for current
production claims, but the breadth gap itself has not disappeared.

## Findings Now Partially Outdated

### Backend leakage

This is still relevant, but no longer as severe as it was in March. Since the
review, the package has gained:

- stronger execution-plan seams
- much better CUDA and AMDGPU parity coverage
- maintained HEART runtime artifacts across CPU, CUDA, and AMDGPU
- substantial ROCm path cleanup on the maintained Shack-Hartmann surfaces

So backend-policy leakage is still worth reducing, but it is no longer the same
level of platform risk described in the March snapshot.

### Documentation quality

The criticism was fair at the time, but the docs are now better curated. The
repo now has clearer stable entry points, especially:

- [documentation-map.md](./documentation-map.md)
- [supported-production-surfaces.md](./supported-production-surfaces.md)
- [production-readiness-checklist.md](./production-readiness-checklist.md)
- [release-validation-runbook.md](./release-validation-runbook.md)

The docs are still rich and plan-heavy, but the information architecture is now
meaningfully better than the March review snapshot.

### External equivalence and operational readiness

The March review predated the current production-hardening pass. Since then the
package has gained:

- a matched OOPAO/Julia HEART baseline
- a second frozen OOPAO equivalence artifact for Pyramid/BioEdge
- scientist-owned HEART boundary truth provenance
- checked-in CPU/CUDA/AMDGPU validation workflows and runbooks

That means the package is now closer to a defendable production AO toolchain on
its supported surfaces than the March review could claim.

## Current Best Summary

A current high-level reading is:

1. The audits remain relevant.
2. The March package review is still correct in its core diagnosis.
3. The package is now materially further along on validation, GPU support, and
   production hardening than that review captured.
4. The next work should stay focused on curation, modularization, supported
   scope, and instrument-truth alignment rather than language-level rewrites.

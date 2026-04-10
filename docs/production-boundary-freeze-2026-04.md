# Production Boundary Freeze: 2026-04

Status: active

## Decision

For the current release train, the production-supported boundary is frozen to
the surfaces already listed in [supported-production-surfaces.md](./supported-production-surfaces.md).

This hardening pass does not expand the scientific production claim beyond those
surfaces.

## Rationale

The package already has enough evidence to support a scoped production claim on
its maintained surfaces:

- HEART OOPAO/Julia null/noise-free equivalence
- HEART CPU/CUDA/AMDGPU runtime artifacts
- HEART detector-response realism runtime artifact (`null` vs `default`)
- frozen OOPAO diffractive Pyramid/BioEdge equivalence artifact
- scientist-owned HEART boundary truth artifact

That is sufficient for the current release train. Additional realism surfaces
remain valuable follow-on work, but they are not required to close the current
production-hardening pass.

## Consequence

For this release train:

- no further fidelity-surface expansion is required before making the scoped
  production claim
- new surfaces should still be added later as supporting confidence artifacts
  or future production-boundary expansions

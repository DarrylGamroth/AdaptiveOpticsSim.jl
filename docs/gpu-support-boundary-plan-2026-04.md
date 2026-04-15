# GPU Support Boundary Plan 2026-04

Status: implemented

## Purpose

This plan makes the GPU support contract explicit and operational.

The package should not have an ambiguous class of "GPU-touched but maybe
supported" paths. A backend surface is either:

- supported and release-gated
- auxiliary validation or audit coverage
- not supported on that backend

## Actions

### GSB-1 Define the supported GPU surface

The supported CUDA and AMDGPU surfaces are the ones covered by the dedicated
hardware validation targets:

- `test/runtests_cuda.jl`
- `test/runtests_amdgpu.jl`

These targets must cover:

- optional backend functional/parity checks
- maintained runtime-equivalence contracts
- maintained high-accuracy post-command runtime-equivalence contracts

Status:

- implemented

### GSB-2 Separate support claims from auxiliary audit coverage

Broader scripts such as `scripts/gpu_smoke_contract.jl` remain useful for
backend auditing and subsystem triage, but they do not by themselves define the
supported GPU surface.

A surface is not support-claimed just because it is touched by a broad GPU
script.

Status:

- implemented

### GSB-3 Encode the policy in maintained docs

The support boundary must be stated consistently in:

- `docs/supported-production-surfaces.md`
- `docs/backend-validation-guide.md`
- `docs/documentation-map.md`

Status:

- implemented

### GSB-4 Prune superseded one-off plan docs

Implemented plan documents that no longer carry live traceability outside the
documentation map should be removed so the docs set stays navigable.

Status:

- implemented

## Outcome

Current release-gated GPU support is defined by the maintained dedicated
hardware targets on real CUDA and AMDGPU hosts.

Broader backend scripts remain valuable for audit and investigation, but they
are not treated as a support claim unless their surfaces are also promoted into
the maintained hardware targets and release-validation path.

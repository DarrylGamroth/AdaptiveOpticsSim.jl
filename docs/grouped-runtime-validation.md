# Grouped Runtime Validation

Date: 2026-04-01

Status: active

Plan traceability:

- [`GR-4`](./grouped-runtime-plan.md)
- validity family: `MV-13`

## Purpose

This note records the maintained validation/evidence surface for grouped
runtime execution.

Grouped execution is a Julia-first runtime capability. The current goal is to
make its product ownership and runtime behavior explicit and reproducible,
without claiming cross-package equivalence where no maintained external
scenario exists yet.

## Artifact

- [2026-04-01-gr.toml](../benchmarks/results/grouped/2026-04-01-gr.toml)
- manifest:
  [manifest.toml](../benchmarks/results/grouped/manifest.toml)

The artifact records:

- CPU grouped runtime timing on the maintained `multi_source_multi_wfs` family
- whether compatible grouped branches expose a grouped WFS stack export
- whether mixed grouped branches correctly refuse a grouped WFS stack export

## Generator

- [generate_grouped_runtime_artifact.jl](../scripts/generate_grouped_runtime_artifact.jl)

Default command:

```bash
julia --project=. --startup-file=no scripts/generate_grouped_runtime_artifact.jl
```

## Evidence Shape

This is a Julia-first grouped-runtime artifact.

It is intended to answer:

- does grouped execution produce a maintained exported-product surface?
- are compatible grouped WFS branches stackable?
- do mixed grouped branches remain explicitly unstacked?
- what is the current grouped runtime cost on maintained CPU scenarios?

## Scope Limits

- The artifact is currently CPU-only.
- Optional AMDGPU/CUDA grouped smoke remains separate in the backend-validation
  surfaces.
- This is grouped-runtime evidence, not a claim of SPECULA equivalence.

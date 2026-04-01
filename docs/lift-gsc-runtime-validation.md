# LiFT and Gain-Sensing Runtime Validation

Date: 2026-04-01

Status: active

Plan traceability:

- [`PVP-05`](./post-review-platform-plan.md)
- validity family: `MV-08`

## Purpose

This note records the maintained workflow-profile artifact for LiFT and the
gain-sensing camera.

The package already had:

- frozen OOPAO-oriented reference coverage
- unit and tutorial execution coverage

What was missing was a maintained runtime/profile surface comparable to the
other optics and runtime families. This artifact closes that gap with a small,
deterministic CPU baseline.

## Artifact

- [2026-04-01-phase1-pvp05.toml](../benchmarks/results/workflows/2026-04-01-phase1-pvp05.toml)
- manifest:
  [manifest.toml](../benchmarks/results/workflows/manifest.toml)

Recorded workflow cases:

- `lift_reconstruct`
  - tutorial-like modal PSF fitting with analytic LiFT
- `gain_sensing_camera`
  - tutorial-like Pyramid gain-sensing calibration and per-frame optical-gain
    estimation

## Generator

- [generate_lift_gsc_profile_artifact.jl](../scripts/generate_lift_gsc_profile_artifact.jl)

Default command:

```bash
julia --project=. --startup-file=no scripts/generate_lift_gsc_profile_artifact.jl
```

## Evidence Shape

This is a maintained workflow-profile artifact, not an external parity bundle.

It records:

- scenario identity and dimensions
- build or first-setup timing where useful
- steady-state mean and p95 timing
- per-call allocation bytes
- one lightweight outcome metric per workflow

This is enough to make LiFT and gain-sensing runtime claims visible in the
validity matrix without turning them into a large benchmark family.

## Scope Limits

- The artifact is currently CPU-only.
- It is workflow evidence, not a cross-package runtime comparison.
- It does not replace the OOPAO-oriented reference bundle for numerical
  agreement.
- Gain-sensing closed-loop atmosphere replay remains covered primarily by
  tutorial/runtime traces rather than by a separate dedicated benchmark family.

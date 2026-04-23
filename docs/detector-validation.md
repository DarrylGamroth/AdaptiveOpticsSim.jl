# Detector Validation

Date: 2026-04-23

Status: active

Plan traceability:

- [`PVP-04`](./post-review-platform-plan.md)
- validity family: `MV-04`

## Purpose

This note describes the committed detector-family validation artifact used to
strengthen detector evidence beyond integrated runtime smoke.

The goal is not to replace the detector testset in
[detectors_and_wfs.jl](../test/testsets/detectors_and_wfs.jl). The goal is to
archive a small fixed-seed fixture report that makes detector-family behavior
visible as maintained evidence in the validity matrix.

## Artifact

- [2026-04-23-phase2-pvp04.toml](../benchmarks/results/detectors/2026-04-23-phase2-pvp04.toml)
- manifest:
  [manifest.toml](../benchmarks/results/detectors/manifest.toml)

Recorded detector families:

- CCD / EMCCD gain-chain comparison
- CMOS structured output with PRNU/DSNU/bad-pixel composition
- InGaAs persistence
- HgCdTe avalanche gain and excess-noise behavior
- HgCdTe multi-read readout products and sampling-mode interactions
- HgCdTe timing, windowing, read-overhead, and readout-correction interactions
- APD counting chain with gating, dead time, and afterpulsing

## Generator

- [generate_detector_validation_artifact.jl](../scripts/generate_detector_validation_artifact.jl)

Default command:

```bash
julia --project=. --startup-file=no scripts/generate_detector_validation_artifact.jl
```

## Evidence Shape

The artifact is a deterministic fixture report, not an external parity bundle.

It records:

- fixed-seed outputs or summary statistics for representative family behaviors
- explicit `contract_holds` booleans for each detector-family case
- an interpretation block summarizing whether the archived contract still holds

This is intentionally narrower than a full detector reference platform. It is a
maintained evidence layer that complements:

- detailed detector unit/functional tests
- optional GPU backend smoke
- integrated AO runtime profiles

## Current Scope Limits

- This artifact is Julia-side maintained evidence, not an external CCD/CMOS/HgCdTe
  parity dataset.
- It does not yet archive full output frames for every detector family.
- It does not replace detector realism evidence in integrated AO scenarios; it
  makes the family-specific transformations easier to audit independently.

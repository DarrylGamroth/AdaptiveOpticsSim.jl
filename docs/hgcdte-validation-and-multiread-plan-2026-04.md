## HgCdTe Validation And Multi-Read Refactor Plan

Status: completed

## Purpose

This plan closes two specific gaps in the detector subsystem:

1. HgCdTe validation is strong at the unit-test level but still too narrow in
   the maintained detector-family artifact.
2. The reusable multi-read frame-readout layer is currently embedded inside the
   HgCdTe detector implementation instead of living as a detector sublayer that
   future sensors can reuse.

## Goals

- expand maintained HgCdTe validation beyond isolated feature checks
- factor reusable multi-read sampling and readout-product assembly out of
  `hgcdte_avalanche_array.jl`
- keep avalanche-specific physics in the HgCdTe family file
- avoid broad detector API churn

## Implementation Steps

### HMR-1: Shared Multi-Read Detector Layer

Move the following into a detector-generic sampling/readout layer:

- `SingleRead`
- `AveragedNonDestructiveReads`
- `CorrelatedDoubleSampling`
- `FowlerSampling`
- read-count and read-time helpers
- read-cube averaging helpers
- windowed frame/cube packing helpers
- reusable explicit multi-read readout-product type

Deliverables:

- shared `MultiReadFrameReadoutProducts`
- shared sampling-mode semantics on `FrameSamplingMode`
- shared helper functions for assembling reference/signal/read cubes and read
  times

### HMR-2: HgCdTe On Shared Multi-Read Layer

Retain HgCdTe-specific behavior in `hgcdte_avalanche_array.jl`:

- avalanche gain
- excess noise
- glow accumulation
- sensor-specific read cadence scaling
- avalanche-aware saturation

But route explicit multi-read readout-product assembly through the shared
multi-read layer.

Deliverables:

- `HgCdTeReadoutProducts` compatibility preserved through a type alias or
  equivalent maintained API surface
- no duplicated sensor-local window/cube/read-time assembly logic

### HMR-3: HgCdTe Validation Expansion

Expand the detector validation artifact and tests so richer HgCdTe behavior is
visible as maintained evidence.

Add artifact cases for:

- NDR read-cube contract
- CDS combined-frame contract
- Fowler read-count/read-time contract
- read-time plus window scaling
- readout correction interaction on explicit multi-read products
- dark/glow accumulation with read overhead

Deliverables:

- expanded detector validation artifact
- updated detector validation docs
- stronger HgCdTe interaction coverage in `test/testsets/detectors.jl`

## Success Criteria

- reusable multi-read machinery no longer lives only inside the HgCdTe file
- HgCdTe keeps only family-specific physics and policy
- maintained validation evidence explicitly covers richer HgCdTe interactions
- future multi-read sensors can reuse the same readout-product layer

## Completion

Completed in the first maintained slice:

- `MultiReadFrameReadoutProducts` is now the reusable multi-read frame-detector
  readout payload, with `HgCdTeReadoutProducts` preserved as a compatibility
  alias
- shared sampling-mode and read-cube assembly logic now lives in
  `src/Detectors/frame_sampling.jl`
- `HgCdTeAvalancheArraySensor` now supplies family-specific physics and
  cadence policy while delegating reusable multi-read assembly to the shared
  layer
- the maintained detector validation artifact now includes richer HgCdTe cases
  for avalanche behavior, multi-read products, timing/window scaling, and
  readout-correction interaction
- validation passed through:
  - `scripts/generate_detector_validation_artifact.jl`
  - `Pkg.test(test_args=["detectors"])`
  - `test/runtests_amdgpu.jl`

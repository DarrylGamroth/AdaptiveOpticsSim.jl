# Reusable Infrastructure Inventory

Date: 2026-03-31

Status: implemented in part

Plan traceability:

- [`PLAN-03`](./package-review-action-plan.md)
- review IDs: `PR-09`, `PR-10`, `PR-11`, `PR-12`, `PR-13`, `PR-14`

## Purpose

This document inventories repeated orchestration and backend-policy patterns
that should be turned into shared infrastructure instead of remaining
model-local code.

The goal is not to refactor the code in this document. The goal is to:

- identify concrete duplication patterns
- point to the current occurrences
- propose a shared abstraction and likely owning module
- record the likely payoff and abstraction risk

Phase 3 implementation status:

- `RI-01` implemented via shared detector pipeline helpers in
  [`src/Detectors/pipeline.jl`](../src/Detectors/pipeline.jl)
- `RI-02` implemented via grouped WFS helpers in
  [`src/WFS/grouped.jl`](../src/WFS/grouped.jl)
- `RI-03` implemented via shared runtime product planning in
  [`src/Control/products.jl`](../src/Control/products.jl)
- `RI-04` implemented via shared propagation context helpers in
  [`src/Optics/propagation_context.jl`](../src/Optics/propagation_context.jl)
- `RI-05` implemented in part through shared reduction infrastructure in
  [`src/Core/reductions.jl`](../src/Core/reductions.jl) and grouped WFS helpers
- `RI-06` implemented in part through shared random/noise services in
  [`src/Core/random_services.jl`](../src/Core/random_services.jl)
- `RI-07` implemented via shared calibration scaffolding in
  [`src/WFS/calibration.jl`](../src/WFS/calibration.jl)

## Inventory

### `RI-01`: Detector frame-capture pipeline scaffolding

Duplicated pattern:

- apply response
- apply photon noise
- apply dark current
- apply pre-readout gain
- apply readout noise
- apply post-readout gain
- finalize readout products
- apply readout correction

Current occurrences:

- [`src/Detectors/frame_capture.jl`](../src/Detectors/frame_capture.jl)
- sensor specializations in:
  - [`src/Detectors/ccd.jl`](../src/Detectors/ccd.jl)
  - [`src/Detectors/cmos.jl`](../src/Detectors/cmos.jl)
  - [`src/Detectors/emccd.jl`](../src/Detectors/emccd.jl)
  - [`src/Detectors/ingaas.jl`](../src/Detectors/ingaas.jl)
  - [`src/Detectors/hgcdte_avalanche_array.jl`](../src/Detectors/hgcdte_avalanche_array.jl)
- batched analog in [`src/Detectors/frame_batched.jl`](../src/Detectors/frame_batched.jl)

Proposed shared abstraction:

- a staged detector-pipeline service layer with reusable phase helpers and a
  clear scalar vs batched split

Likely owner:

- `src/Detectors/`

Risk of abstraction:

- medium

Expected payoff:

- high

Implementation note:

- shared detector pipeline helpers now own reusable scalar capture/readout
  staging in [`src/Detectors/pipeline.jl`](../src/Detectors/pipeline.jl)

### `RI-02`: Grouped WFS accumulation skeleton

Duplicated pattern:

- build grouped per-source or per-sample stack
- run per-slice propagation/sampling
- reduce stack into intensity
- continue with signal extraction

Current occurrences:

- Shack-Hartmann grouped stack paths in
  [`src/WFS/shack_hartmann.jl`](../src/WFS/shack_hartmann.jl)
- Pyramid grouped accumulation in
  [`src/WFS/pyramid.jl`](../src/WFS/pyramid.jl)
- BioEdge grouped accumulation in
  [`src/WFS/bioedge.jl`](../src/WFS/bioedge.jl)

Proposed shared abstraction:

- grouped sensing execution helper layer with explicit stack-capacity,
  prepare/reduce stages, and backend hooks

Likely owner:

- `src/WFS/`

Risk of abstraction:

- medium-high

Expected payoff:

- high

Implementation note:

- shared grouped reduction/view helpers now live in
  [`src/WFS/grouped.jl`](../src/WFS/grouped.jl)

### `RI-03`: Runtime product planning and export policy

Duplicated pattern:

- decide which products are needed
- allocate or skip frame buffers
- snapshot WFS/science products
- expose readout metadata consistently

Current occurrences:

- product policy in [`src/Control/runtime.jl`](../src/Control/runtime.jl)
- SH pixel export behavior in
  [`src/WFS/shack_hartmann.jl`](../src/WFS/shack_hartmann.jl)
- detector readout-product surfaces in
  [`src/Detectors/generic.jl`](../src/Detectors/generic.jl)

Proposed shared abstraction:

- a runtime product-planning layer lowered from high-level requirements into a
  concrete per-runtime product plan

Likely owner:

- `src/Control/`

Risk of abstraction:

- medium

Expected payoff:

- high

Implementation note:

- runtime product requirements and lowering now live in
  [`src/Control/products.jl`](../src/Control/products.jl)

### `RI-04`: Source / field / atmosphere accumulation pipeline

Duplicated pattern:

- derive source-aware geometry
- build field/sampled layer
- accumulate into OPD or intensity output
- provide on-axis fast path where possible

Current occurrences:

- source-aware atmosphere propagation in:
  - [`src/Atmosphere/source_geometry.jl`](../src/Atmosphere/source_geometry.jl)
  - [`src/Atmosphere/multilayer.jl`](../src/Atmosphere/multilayer.jl)
  - [`src/Atmosphere/infinite_screen.jl`](../src/Atmosphere/infinite_screen.jl)
- atmospheric field propagation in
  [`src/Optics/atmospheric_field_propagation.jl`](../src/Optics/atmospheric_field_propagation.jl)
- spectral field slicing in [`src/Optics/spectrum.jl`](../src/Optics/spectrum.jl)

Proposed shared abstraction:

- a maintained propagation-context helper layer for source geometry, spectral
  slicing, and accumulation ownership

Likely owner:

- `src/Optics/` plus shared atmosphere helpers

Risk of abstraction:

- medium

Expected payoff:

- medium-high

Implementation note:

- source-aware layer render context and shared render entry points now live in
  [`src/Optics/propagation_context.jl`](../src/Optics/propagation_context.jl)

### `RI-05`: Backend reduction services

Duplicated pattern:

- masked sum / valid-pixel reductions
- grouped stack reductions
- backend-specific dispatch for accelerator reductions

Current occurrences:

- shared masked reduction in [`src/Core/reductions.jl`](../src/Core/reductions.jl)
- CUDA extension hook in [`ext/AdaptiveOpticsSimCUDAExt.jl`](../ext/AdaptiveOpticsSimCUDAExt.jl)
- model-local grouped reduction kernels in:
  - [`src/WFS/shack_hartmann.jl`](../src/WFS/shack_hartmann.jl)
  - [`src/WFS/pyramid.jl`](../src/WFS/pyramid.jl)
  - [`src/WFS/bioedge.jl`](../src/WFS/bioedge.jl)

Proposed shared abstraction:

- backend-aware reduction service layer with a small maintained set of common
  reduction primitives

Likely owner:

- `src/Core/`

Risk of abstraction:

- medium

Expected payoff:

- high

Implementation note:

- grouped WFS reduction kernels were moved onto shared helpers; masked reduction
  services remain centered in [`src/Core/reductions.jl`](../src/Core/reductions.jl)

### `RI-06`: Backend random-fill and noise services

Duplicated pattern:

- Poisson noise application
- Gaussian random fill
- ROCm/CUDA specialization and host-mirror fallback handling

Current occurrences:

- generic helpers in [`src/Core/utils.jl`](../src/Core/utils.jl)
- AMDGPU specialization in [`ext/AdaptiveOpticsSimAMDGPUExt.jl`](../ext/AdaptiveOpticsSimAMDGPUExt.jl)
- detector-owned wrappers in [`src/Detectors/frame_capture.jl`](../src/Detectors/frame_capture.jl)
- atmosphere Gaussian fill in [`src/Atmosphere/kolmogorov.jl`](../src/Atmosphere/kolmogorov.jl)

Proposed shared abstraction:

- backend service layer separating:
  - device-native random fill
  - detector-owned host-mirror fallback
  - atmosphere-owned host-mirror fallback

Likely owner:

- `src/Core/` with narrow opt-in hooks from detectors/atmosphere

Risk of abstraction:

- medium-high

Expected payoff:

- high

Implementation note:

- random fill and Poisson services were moved out of `utils.jl` into
  [`src/Core/random_services.jl`](../src/Core/random_services.jl), with backend
  extensions remaining responsible for specializations

### `RI-07`: WFS calibration orchestration skeleton

Duplicated pattern:

- prepare sampling
- ensure calibration state
- build or copy reference signal
- cache wavelength/signature

Current occurrences:

- Shack-Hartmann calibration in
  [`src/WFS/shack_hartmann.jl`](../src/WFS/shack_hartmann.jl)
- Pyramid calibration in [`src/WFS/pyramid.jl`](../src/WFS/pyramid.jl)
- BioEdge calibration in [`src/WFS/bioedge.jl`](../src/WFS/bioedge.jl)
- Curvature calibration in [`src/WFS/curvature.jl`](../src/WFS/curvature.jl)
- Zernike calibration in [`src/WFS/zernike.jl`](../src/WFS/zernike.jl)
- runtime preparation hooks in [`src/Control/runtime.jl`](../src/Control/runtime.jl)

Proposed shared abstraction:

- common calibration-state contract plus family-specific prepare/build hooks

Likely owner:

- `src/WFS/` with runtime integration in `src/Control/`

Risk of abstraction:

- medium

Expected payoff:

- high

### `RI-08`: Modal calibration and reconstructor assembly

Duplicated pattern:

- build modal basis
- measure interaction matrix
- construct inverse/operator object
- wrap it into runtime-friendly reconstructor form

Current occurrences:

- [`src/Calibration/modal_basis.jl`](../src/Calibration/modal_basis.jl)
- [`src/Calibration/interaction_matrix.jl`](../src/Calibration/interaction_matrix.jl)
- [`src/Calibration/calibration_vault.jl`](../src/Calibration/calibration_vault.jl)
- [`src/Calibration/reconstructor.jl`](../src/Calibration/reconstructor.jl)
- [`src/Calibration/ao_calibration.jl`](../src/Calibration/ao_calibration.jl)

Proposed shared abstraction:

- explicit calibration-pipeline builder surface, with reusable artifact types
  for basis, IM, inverse, and runtime reconstructor packaging

Likely owner:

- `src/Calibration/`

Risk of abstraction:

- low-medium

Expected payoff:

- medium-high

Implementation note:

- shared calibration signature/wavelength/reference helpers now live in
  [`src/WFS/calibration.jl`](../src/WFS/calibration.jl)

### `RI-09`: Runtime grouped execution orchestration

Duplicated pattern:

- prepare child runtimes
- sense all branches
- reconstruct all branches
- snapshot outputs
- expose grouped readout surfaces

Current occurrences:

- [`src/Control/runtime.jl`](../src/Control/runtime.jl) around
  `CompositeSimulationInterface` and grouped `step!`

Proposed shared abstraction:

- separate grouped-runtime orchestration helper layer rather than leaving all
  grouped logic inside one runtime file

Likely owner:

- `src/Control/`

Risk of abstraction:

- medium

Expected payoff:

- medium

## Priority Order

Recommended first extraction targets:

1. `RI-03` runtime product planning and export policy
2. `RI-07` WFS calibration orchestration skeleton
3. `RI-01` detector frame-capture pipeline scaffolding
4. `RI-05` backend reduction services
5. `RI-02` grouped WFS accumulation skeleton

Reasoning:

- these give the best maintainability payoff with the least architectural risk
- they also align with already-visible structural seams from the recent GPU
  runtime work

## Areas to Avoid Abstracting Too Early

- exact FFT plan ownership inside individual WFS families
- sensor-family physics that only look similar at a superficial level
- grouped WFS execution beyond a small common skeleton
- atmospheric statistical kernels whose shared shape is still evolving

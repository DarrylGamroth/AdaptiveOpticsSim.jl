# Modularization Inventory

Date: 2026-03-31

Status: active

Plan traceability:

- [`PLAN-02`](./package-review-action-plan.md)
- review IDs: `PR-05`, `PR-06`, `PR-07`

## Purpose

This document inventories the oversized source and test files identified in the
package review and proposes target split boundaries for later implementation.

The goal is not to refactor the code in this document. The goal is to:

- record the current file-size baseline
- identify responsibility clusters inside each oversized file
- propose concrete split targets
- identify risk points and likely dependency seams
- rank the splits by payoff and difficulty

## Baseline File Sizes

Current line counts:

- [`src/WFS/shack_hartmann.jl`](../src/WFS/shack_hartmann.jl): `2167`
- [`src/WFS/pyramid.jl`](../src/WFS/pyramid.jl): `1463`
- [`src/WFS/bioedge.jl`](../src/WFS/bioedge.jl): `1241`
- [`src/Control/runtime.jl`](../src/Control/runtime.jl): `1123`
- [`test/runtests.jl`](../test/runtests.jl): `3300`

## Phase 2 implementation status

Phase 2 completed the first internal modularization pass without changing the
package boundary.

Current entry-point line counts:

- [`src/WFS/shack_hartmann.jl`](../src/WFS/shack_hartmann.jl): `12`
- [`src/WFS/pyramid.jl`](../src/WFS/pyramid.jl): `11`
- [`src/WFS/bioedge.jl`](../src/WFS/bioedge.jl): `10`
- [`src/Control/runtime.jl`](../src/Control/runtime.jl): `12`
- [`test/runtests.jl`](../test/runtests.jl): `8`

Current split layout:

- Shack-Hartmann:
  - [`src/WFS/shack_hartmann/setup.jl`](../src/WFS/shack_hartmann/setup.jl)
  - [`src/WFS/shack_hartmann/measure.jl`](../src/WFS/shack_hartmann/measure.jl)
  - [`src/WFS/shack_hartmann/stacks.jl`](../src/WFS/shack_hartmann/stacks.jl)
  - [`src/WFS/shack_hartmann/signals.jl`](../src/WFS/shack_hartmann/signals.jl)
  - [`src/WFS/shack_hartmann/lgs.jl`](../src/WFS/shack_hartmann/lgs.jl)
- Pyramid:
  - [`src/WFS/pyramid/setup.jl`](../src/WFS/pyramid/setup.jl)
  - [`src/WFS/pyramid/measure.jl`](../src/WFS/pyramid/measure.jl)
  - [`src/WFS/pyramid/optics.jl`](../src/WFS/pyramid/optics.jl)
  - [`src/WFS/pyramid/signals.jl`](../src/WFS/pyramid/signals.jl)
- BioEdge:
  - [`src/WFS/bioedge/setup.jl`](../src/WFS/bioedge/setup.jl)
  - [`src/WFS/bioedge/measure.jl`](../src/WFS/bioedge/measure.jl)
  - [`src/WFS/bioedge/signals.jl`](../src/WFS/bioedge/signals.jl)
- Runtime:
  - [`src/Control/runtime/types.jl`](../src/Control/runtime/types.jl)
  - [`src/Control/runtime/construction.jl`](../src/Control/runtime/construction.jl)
  - [`src/Control/runtime/outputs.jl`](../src/Control/runtime/outputs.jl)
  - [`src/Control/runtime/execution.jl`](../src/Control/runtime/execution.jl)
  - [`src/Control/runtime/timing.jl`](../src/Control/runtime/timing.jl)
- Tests:
  - [`test/runtests_head.jl`](../test/runtests_head.jl)
  - [`test/testsets/core_optics.jl`](../test/testsets/core_optics.jl)
  - [`test/testsets/atmosphere.jl`](../test/testsets/atmosphere.jl)
  - [`test/testsets/control_and_runtime.jl`](../test/testsets/control_and_runtime.jl)
  - [`test/testsets/detectors_and_wfs.jl`](../test/testsets/detectors_and_wfs.jl)
  - [`test/testsets/calibration_and_analysis.jl`](../test/testsets/calibration_and_analysis.jl)
  - [`test/testsets/reference_and_tutorials.jl`](../test/testsets/reference_and_tutorials.jl)

## Inventory

### 1. `src/WFS/shack_hartmann.jl`

Current responsibility clusters:

1. type definitions and constructor
   - params/state/type
   - constructor and setup
2. static preparation and workspace management
   - valid mask
   - phasor/sampling prep
   - buffer sizing
3. geometric and diffractive measurement entry points
4. grouped diffractive execution
   - asterism stack paths
   - spectral stack paths
   - batched GPU kernels
5. spot sampling and centroid/statistics logic
6. detector-coupled measurement paths
7. LGS-specific elongation/kernel handling
8. calibration/reference helpers
9. low-level GPU kernels

Recommended split target:

- `src/WFS/shack_hartmann_types.jl`
  - params/state/type
  - constructor
- `src/WFS/shack_hartmann_sampling.jl`
  - valid mask
  - phasor
  - buffer prep
  - sample-spot helpers
- `src/WFS/shack_hartmann_measure.jl`
  - main geometric/diffractive `measure!` dispatch surface
- `src/WFS/shack_hartmann_grouped.jl`
  - spectral
  - asterism
  - batched grouped execution
- `src/WFS/shack_hartmann_detector.jl`
  - detector-coupled paths
  - pixel-product ownership
- `src/WFS/shack_hartmann_lgs.jl`
  - elongation
  - kernel generation
- `src/WFS/shack_hartmann_calibration.jl`
  - reference signal
  - calibration ramp
  - centroid sums
- `src/WFS/shack_hartmann_kernels.jl`
  - KA kernels only

Highest-risk dependencies:

- grouped execution shares state fields with detector and centroid paths
- calibration logic currently reaches into low-level buffers directly
- LGS support is interleaved with ordinary diffractive measurement

Recommended split order:

1. types/constructor
2. kernels
3. calibration/LGS helpers
4. grouped execution
5. detector coupling
6. main `measure!` dispatch surface cleanup

Expected payoff:

- very high

Implementation risk:

- high

### 2. `src/WFS/pyramid.jl`

Current responsibility clusters:

1. low-level GPU kernels
2. params/state/type and constructor
3. valid-mask and signal indexing maintenance
4. grouped accumulation
   - asterism
   - spectral
   - extended-source
5. sampling and signal-buffer setup
6. measurement entry points
7. mask/phasor/modulation construction
8. intensity and signal normalization
9. calibration and optical-gain handling
10. LGS handling

Recommended split target:

- `src/WFS/pyramid_types.jl`
- `src/WFS/pyramid_masks.jl`
  - phasor
  - mask construction
  - modulation phases
- `src/WFS/pyramid_grouped.jl`
  - grouped asterism/spectral/extended-source accumulation
- `src/WFS/pyramid_signal.jl`
  - `pyramid_signal!`
  - normalization
  - valid-i4q selection
- `src/WFS/pyramid_measure.jl`
  - main `measure!` dispatch
- `src/WFS/pyramid_lgs.jl`
- `src/WFS/pyramid_calibration.jl`
- `src/WFS/pyramid_kernels.jl`

Highest-risk dependencies:

- calibration and valid-signal indexing are tightly coupled to signal layout
- grouped execution and detector-coupled execution share setup assumptions

Recommended split order:

1. types
2. kernels
3. masks/modulation
4. calibration
5. grouped execution
6. signal/measure split

Expected payoff:

- high

Implementation risk:

- medium-high

### 3. `src/WFS/bioedge.jl`

Current responsibility clusters:

1. low-level GPU kernels
2. params/state/type and constructor
3. valid-mask and edge-mask logic
4. phasor/mask/modulation setup
5. grouped asterism accumulation
6. sampling and signal-buffer preparation
7. measurement entry points
8. signal normalization
9. valid-i4q and valid-signal selection
10. calibration and optical-gain handling
11. LGS handling

Recommended split target:

- `src/WFS/bioedge_types.jl`
- `src/WFS/bioedge_masks.jl`
  - edge mask
  - phasor
  - mask construction
- `src/WFS/bioedge_grouped.jl`
- `src/WFS/bioedge_signal.jl`
  - signal and normalization
  - valid-i4q handling
- `src/WFS/bioedge_measure.jl`
- `src/WFS/bioedge_lgs.jl`
- `src/WFS/bioedge_calibration.jl`
- `src/WFS/bioedge_kernels.jl`

Highest-risk dependencies:

- edge-mask handling is a core setup surface used throughout the file
- signal-layout assumptions are shared by grouped and detector-facing paths

Recommended split order:

1. types
2. kernels
3. masks
4. calibration/LGS helpers
5. grouped execution
6. signal/measure split

Expected payoff:

- high

Implementation risk:

- medium

### 4. `src/Control/runtime.jl`

Current responsibility clusters:

1. execution-policy and runtime-profile definitions
2. latency model and delay-line infrastructure
3. runtime state structs
   - `ClosedLoopRuntime`
   - `SimulationInterface`
   - `CompositeSimulationInterface`
   - `SimulationReadout`
4. constructor/build logic
5. preparation logic
6. interface composition
7. main sense/step execution
8. timing helpers and profiling
9. runtime-product policy

Recommended split target:

- `src/Control/runtime_types.jl`
  - policies
  - profiles
  - data structs
- `src/Control/runtime_latency.jl`
  - delay lines
  - latency model
- `src/Control/runtime_prepare.jl`
  - build/preparation
  - prepared-state helpers
- `src/Control/runtime_interfaces.jl`
  - interface wrappers
  - readout views
  - composition
- `src/Control/runtime_execute.jl`
  - `sense!`
  - `step!`
  - command application and staged execution
- `src/Control/runtime_timing.jl`
  - timing and profiling helpers

Highest-risk dependencies:

- preparation and execution both reach into many runtime fields directly
- product/export policy and execution are still intertwined
- interfaces are thin wrappers around runtime internals, so type boundaries
  must be preserved carefully

Recommended split order:

1. types
2. timing
3. latency
4. interfaces/readout
5. prepare
6. execute

Expected payoff:

- very high

Implementation risk:

- high

### 5. `test/runtests.jl`

Current responsibility clusters:

1. interface-conformance helpers
2. core optics
3. atmosphere
4. detector
5. WFS families
6. runtime/control
7. tomography and calibration
8. reference harness and parity checks
9. tutorials/examples

Recommended split target:

- `test/core_optics.jl`
- `test/atmosphere.jl`
- `test/detectors.jl`
- `test/wfs_shack_hartmann.jl`
- `test/wfs_pyramid_bioedge_lgs.jl`
- `test/wfs_curvature_zernike_lift.jl`
- `test/runtime_control.jl`
- `test/calibration.jl`
- `test/tomography.jl`
- `test/reference_regression.jl`
- `test/tutorial_examples.jl`
- `test/interface_conformance.jl`
- keep `test/runtests.jl` as a small include driver

Highest-risk dependencies:

- a number of helper functions at the top of the current file are shared by
  many later testsets
- ordering assumptions may exist where one testset seeds reusable artifacts
  used by later sections

Recommended split order:

1. move helper functions into `test/test_helpers.jl`
2. split interface/reference/tutorial sections first
3. split subsystem tests second
4. reduce `runtests.jl` to includes and high-level orchestration

Expected payoff:

- very high

Implementation risk:

- medium

## Files That Should Not Be Split Yet

These are not priority split targets yet:

- [`src/AdaptiveOpticsSim.jl`](../src/AdaptiveOpticsSim.jl)
  - large export surface, but not a monolithic algorithm file
- [`src/Atmosphere/multilayer.jl`](../src/Atmosphere/multilayer.jl)
  - still manageable in size compared with the WFS/runtime files
- [`src/Detectors/frame_capture.jl`](../src/Detectors/frame_capture.jl)
  - worth watching, but not the first modularization target

## Ranked Split Order Across The Package

1. `test/runtests.jl`
2. `src/Control/runtime.jl`
3. `src/WFS/shack_hartmann.jl`
4. `src/WFS/pyramid.jl`
5. `src/WFS/bioedge.jl`

Rationale:

- splitting the test file gives earlier maintainability wins and lowers risk
  for the later source refactors
- runtime and Shack-Hartmann carry the most structural concentration
- Pyramid and BioEdge should follow once the split patterns are established

## Initial Recommendation For Phase 2

Use this inventory as the contract for `PLAN-11` through `PLAN-15`.

The preferred implementation pattern is:

- introduce small include files without changing behavior first
- move helper blocks and structs before moving top-level dispatch functions
- keep stable tests green after each split wave
- avoid package-level fragmentation until the internal file boundaries are
  healthier

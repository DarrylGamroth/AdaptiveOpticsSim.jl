# HIL Package Boundary

Status: active

## Purpose

This note records the recommended split between `AdaptiveOpticsSim.jl` and a
future HIL companion package.

`AdaptiveOpticsSim.jl` should remain the adaptive-optics simulation kernel:
typed optical models, atmosphere, WFS, detectors, DMs, calibration primitives,
control primitives, and deterministic runtime stepping.

A HIL companion package should own operational testbench behavior: fixed-rate
execution, telemetry transport, artifact lifecycle, replay, fault injection,
and real-time validation workflows.

## Rule Of Thumb

Put a feature in `AdaptiveOpticsSim.jl` when it describes physics, a reusable
AO computation, or a stable simulation data contract.

Put a feature in the HIL package when it describes when things run, how data
is moved, how artifacts are selected, or how an RTC/testbench is operated.

## Belongs In AdaptiveOpticsSim.jl

- Physical and computational model objects:
  - `Telescope`, sources, atmosphere, WFS, detectors, DMs, controllable optics
  - detector response/readout/thermal/counting models
  - DM topology, influence functions, actuator behavior, and optic application
- Calibration primitives:
  - `InteractionMatrix`
  - `ControlMatrix`
  - `AOCalibration`
  - subaperture calibration/layout data
  - detector export metadata
  - DM topology and command-layout metadata
  - misregistration models and sensitivity computation
- Runtime stepping seams:
  - `ControlLoopScenario`
  - `prepare!`
  - `set_command!`
  - `sense!`
  - `step!`
  - `readout`, `command`, `slopes`, `wfs_frame`, and `science_frame`
- Stable backend behavior:
  - CPU/CUDA/AMDGPU array ownership
  - no-allocation hot-path simulation surfaces
  - backend parity tests and release validation targets
- Simple user-facing calibration and simulation recipes that are useful without
  external hardware or RTC services.

## Belongs In The HIL Companion Package

- Fixed-rate HIL loop orchestration:
  - wall-clock pacing
  - deadline and jitter tracking
  - dropped-frame accounting
  - warmup and soak-test execution
  - command/readout callback scheduling
- Transport and external I/O:
  - shared memory
  - memory-mapped files
  - sockets or message queues
  - pinned host buffers
  - camera/DM/RTC adapters
  - explicit host/device transfer boundaries
- Telemetry lifecycle:
  - ring buffers
  - record/replay files
  - run manifests
  - timestamp policies
  - long-run allocation, GC, backend-sync, and latency summaries
- Calibration artifact management:
  - artifact manifests and schema versions
  - artifact storage layout
  - artifact compatibility checks
  - promotion of calibration runs to maintained artifacts
  - selection of the active calibration for a HIL campaign
- Multi-rate operational timing:
  - detector exposure/integration windows
  - camera start delay
  - RTC delay
  - DM latency
  - WFS/science-camera cadence differences
  - delayed command queues
- Fault injection:
  - dropped or stale frames
  - stale commands
  - timestamp jitter
  - command clipping
  - dead or drifting actuators
  - bad/hot pixels
  - saturation and quantization stress cases
- Operational workflows:
  - calibrate, archive, load, run, replay
  - bench-specific setup and teardown
  - hardware availability checks
  - HIL-specific benchmark reporting

## Split Features

Some concepts need small stable definitions in `AdaptiveOpticsSim.jl` and
operational lifecycle handling in the HIL package.

- Calibration artifacts:
  - `AdaptiveOpticsSim.jl` owns the meaning of matrices, masks, detector
    metadata, DM topology, subaperture layouts, and validity conventions.
  - The HIL package owns serialization, schema versioning, compatibility
    checks, storage, and promotion workflows.
- Misregistration:
  - `AdaptiveOpticsSim.jl` owns the model and estimator math.
  - The HIL package owns calibration campaigns, artifact comparison, and bench
    procedures.
- Controllers and filters:
  - `AdaptiveOpticsSim.jl` may keep simple reusable controller primitives.
  - The HIL package or a control companion should own RTC-like gain scheduling,
    multi-rate filters, dynamic integrators, and telemetry-driven gain
    optimization.
- Visualization:
  - `AdaptiveOpticsSimPlots.jl` owns plotting.
  - The HIL package should emit clean telemetry products and avoid plotting
    dependencies in its hot path.

## Package Comparison Guidance

Borrow concepts, not package structure.

- From SPECULA, borrow the separation between loop control, calibration data,
  buffers, stores, and display. Put the operational pieces in the HIL package.
- From OOPAO, borrow calibration and tutorial ergonomics. Keep simple model
  recipes in `AdaptiveOpticsSim.jl`; put bench procedures in the HIL package.
- From REVOLT, borrow practical command/telemetry archival and hardware
  assumptions. Keep them out of the simulation kernel.
- From FALCO, borrow focal-plane wavefront-control ideas only in a separate
  coronagraph-control companion package. Pair-wise probing, EFC, and dark-hole
  optimization are not core HIL runtime responsibilities.
- From HCIPy, borrow science-arm breadth only through optional packages such as
  `Proper.jl` or future high-contrast-imaging companions.

## Expected HIL Package Shape

A future HIL package should depend on `AdaptiveOpticsSim.jl`, not the other way
around.

The first useful public surface should be small:

- a fixed-step runner around `ControlLoopScenario`
- a zero-allocation telemetry ring buffer
- a record/replay artifact format
- a calibration artifact manifest
- a small set of transport adapters kept outside the hot simulation kernel

If that package exposes missing seams in `AdaptiveOpticsSim.jl`, add the
smallest stable simulation-facing interface here instead of moving operational
logic into core.

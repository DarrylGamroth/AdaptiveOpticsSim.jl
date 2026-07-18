# HIL Product And Package Architecture

Status: active

## Purpose

This guide is the maintained entry point for using `AdaptiveOpticsSim.jl` as
a scalable, low-latency adaptive-optics plant for development and validation
of an external real-time controller (RTC). Prepared product providers span
fast synthetic boundary load, command-responsive reduced-order simulation, and
high-fidelity optical propagation. This guide defines the product boundary and
points to the subsystem specifications that carry the detailed design and
compliance contracts.

The product evolves the atmosphere and telescope, renders direction-dependent
optical paths, models detectors and controllable optics, presents sensor
products on controlled acquisition schedules or triggers, and accepts commands
from an external RTC. It does not replace a production RTC or a general
optical-propagation package.

These documents describe a target architecture. Current support claims remain
defined by [`supported-production-surfaces.md`](supported-production-surfaces.md)
and require maintained functional and hardware evidence before expansion.

## Subsystem Specifications

| Specification | Scope |
|---|---|
| [`hil/plant-and-optics.md`](hil/plant-and-optics.md) | Pre-HIL optical-plane/radiometric ownership, explicit atmosphere time, direct-science acquisition, and WFS decomposition; shared plant, reusable optical paths, independent acquisition endpoints, prepared fidelity providers, calibration illumination, NCPA, co-conjugated optics, MCAO, and MOAO |
| [`hil/time-and-scheduling.md`](hil/time-and-scheduling.md) | Performance contract, execution clocks, detector-trigger distribution, event causality, detector pacing, command timing, lifecycle, overload, and replay |
| [`hil/rtc-ports.md`](hil/rtc-ports.md) | Ownership, bounded SPSC ports, complete-product leases, sequence domains, publication, and resource-specific capacity semantics |
| [`hil/execution-and-placement.md`](hil/execution-and-placement.md) | Prepared CPU execution, affinity, static mixed CPU/GPU and multi-GPU placement, the multi-host boundary, and offline parallelism |
| [`hil/package-boundaries.md`](hil/package-boundaries.md) | Responsibilities of core, the HIL companion, PROPER, and user RTC integration |
| [`hil/validation.md`](hil/validation.md) | Functional, timing, backend, low-fidelity load, and instrument-profile validation, including NFIRAOS synchronization and MORFEO extreme scale |
| [`hil/compliance-matrix.md`](hil/compliance-matrix.md) | Stable requirement IDs, implementation state, evidence state, and capability gates |

Temporary PR checklists and implementation notes belong in PR descriptions or
issues. The subsystem specifications and compliance matrix remain valid after
individual PRs merge.

## Product Position

The primary product boundary is:

```text
 simulated AO plant                          external RTC

 atmosphere + telescope                     reconstruction
 direction-specific optical paths    <->    control law
 detectors and acquisition triggers         command generation
 controllable-optic response                 RTC telemetry
 timestamps and fault injection
```

The simulator owns the physical plant and the presentation of detector data.
The external RTC owns reconstruction, controller state, and the production of
commands. Simple internal reconstructors and controllers remain useful as
correctness oracles, calibration tools, examples, and stand-alone smoke tests;
they do not define the product's RTC scope.

Primary use cases include:

- closed-loop HIL testing of an external RTC
- fast production-shaped RTC latency, throughput, overload, and recovery tests
  that do not require full optical propagation
- reduced-order command-responsive HIL in which the RTC performs real pixel or
  slope processing, reconstruction, tomography, control, and fault recovery
- independently scheduled or externally triggered WFS and science-camera
  acquisitions over heterogeneous optical paths
- direct-imaging and PROPER-backed coronagraph science paths
- NGS, LGS, and mixed-guide-star tomography through one atmosphere
- SCAO, co-conjugated multi-DM, MCAO, and MOAO plant configurations
- independent command timing for DMs, tip/tilt mirrors, focus stages, and other
  controllable optics
- detector, boundary-delivery, stale-data, and missed-deadline fault injection
- common-source and per-detector trigger delay, skew, jitter, dropped/duplicate
  edges, and timestamp-label fault injection
- record/replay and repeatable regression of RTC behavior
- static prepared CPU, GPU, and mixed CPU/GPU execution
- optional single-host process-boundary adapters, including iceoryx2, without
  replacing direct/SPSC ownership inside one HIL process
- device-resident offline operation using the same physical models

Non-goals include:

- implementing a full production RTC inside `AdaptiveOpticsSim.jl`
- replacing `Proper.jl` for detailed relay, coronagraph, or instrument
  propagation
- claiming hard real-time behavior from an ordinary Julia and general-purpose
  operating-system deployment
- using Dagger or another dynamic task-graph scheduler in the HIL critical path
- baking transport, persistence, or artifact codecs into the canonical HIL data
  plane, or embedding bench-specific transports and camera profiles in core

## Current Foundation And Target Gaps

The detailed state in [`hil/compliance-matrix.md`](hil/compliance-matrix.md) is
authoritative. This table is only an orientation to the principal transitions:

| Area | Current foundation | Target contract |
|---|---|---|
| Optical planes and workspaces | `TelescopeAperture` separates prepared geometry from transitional WFS OPD state; caller-owned `PupilFunction`, `ElectricField`, and `IntensityMap` products use separate propagation, direct-imaging, and spatial-filter workspaces, and the telescope has no cadence or exposure field | Complete migration of WFS paths to explicit products, then remove transitional telescope OPD state |
| Time and radiometry | Timed atmospheres receive explicit elapsed or absolute model time; runtime and control-loop configuration require a positive `atmosphere_step`. Optical formation produces declared photon-rate or dimensionless products without elapsed time, and a prepared detector acquisition applies its explicit exposure once | Add independent trigger/exposure scheduling while preserving the separate atmosphere-step, optical-rate, and detector-exposure contracts |
| Atmosphere and source rendering | One writer publishes an explicit stable atmosphere epoch; path-local prepared renderers freeze NGS/LGS/spectral/extended-source descriptions and write caller-owned products without shared render/cache state | Carry this boundary through the remaining WFS/science path decomposition and multi-rate scheduler |
| WFS composition | Shack–Hartmann separates microlens optics, layout/calibration, acquisition, and estimation. Pyramid and BioEdge expose physically distinct focal-plane front ends over shared prepared modulation. Zernike separates its phase spot and pupil relay; Curvature separates its two branch planes and supports independent or packed detector acquisition. All five use explicit photon-rate products, family-neutral acquisition where applicable, and revision-bound estimators | Complete the same observation, forward-model, and estimator separation for LiFT |
| Shared plant and paths | One telescope/atmosphere with sequential shared arms | Immutable event snapshots, path-local propagation workspace, and reusable optical paths separated from acquisition state |
| Detector acquisition | Detector families, incremental accumulation, rolling shutter, frame transfer, and supported nondestructive reads; `DetectorAcquisitionPlan` validates a typed rate product and prepares response/pixel-integration buffers before allocation-free repeated capture | Independent scheduled acquisition and readout/publication events with explicit equal-time semantics, preserving the prepared rate boundary and detector ordering |
| Guide stars and AO modes | NGS, finite-height LGS, elongation, tomography, co-conjugated additive optics | Mixed NGS/LGS scheduling, independent co-conjugated devices, altitude-conjugated MCAO, and path-specific MOAO |
| Controllable optics | One runtime optic or packed composite at the pupil | Named placed optics with independent command endpoints, effective times, visibility, and bounded late/future policies |
| Science optics and NCPA | Prepared direct imaging owns explicit pupil/field/output/workspace products, supports same-grid incoherent composition and spectral bundles, and reuses one rate image across independent detector acquisitions; sampled NCPA primitives and an asynchronous prepared PROPER handoff remain | Caller-owned native or explicitly converted external photon-arrival-rate products feed independently scheduled detector acquisition; later path composition adds explicit NCPA visibility and external-executor scheduling |
| Product fidelity | Full optical models plus separate geometric WFS implementations and deterministic benchmark inputs | One prepared acquisition seam for full, reduced-order, and synthetic/replay providers with invariant RTC boundary semantics and explicit claim limits |
| Calibration illumination | Optical source models plus static and temporal detector-input sources | A narrow typed path-entry seam where users declare injection, visibility, physics, timing, and composition without instrument assumptions in core |
| CPU and accelerators | Direct CPU execution and one implicit backend/device per plant | Prepared path owners plus fully explicit or constrained deterministic CPU/GPU and multi-GPU placement; fully automatic cost optimization is deferred |
| HIL boundary | Application-specific proving-ground integration | Transport-neutral command and complete-product acquisition ports, bounded pools, lifecycle, and resource-specific overload behavior |
| Trigger, clock, and scheduling | Frame-step runtime without injected operational clock or trigger fan-out | Integer event time, modeled trigger generation/distribution, and separate versioned `Clocks.jl` execution pacing/domain mapping |
| Multi-host | No general HIL deployment contract | Separate future surface with explicit clock, transport, ownership, partition, and recovery evidence |
| Offline orchestration | Independent ensembles with optional AcceleratedKernels and Dagger policies | Preserve those policies outside a detector-to-RTC deadline path |

## Design Review Rules

When refining this architecture:

- prefer dispatch and traits over `isa` branches
- keep immutable configuration separate from mutable state
- keep all repeated-path storage prepared and explicitly owned
- parallelize coarse source/path work and avoid nested parallelism
- retain a deterministic serial oracle at every stage
- keep telescope/aperture definitions separate from caller-owned path products
  and keep optical products separate from propagation plans and scratch
- keep telescope spatial sampling separate from model time; advance shared
  atmosphere state explicitly once per epoch and render directions through
  prepared path-local source/geometry workspaces
- attach run-immutable plane geometry, wavelength/channel, radiometric,
  density-versus-cell-integrated measure, combination, backend, and device
  metadata and reject incompatible handoffs during preparation
- form native, external, or WFS photon-arrival-rate products before
  detector-owned presampling response, physical-pixel integration, QE and
  elapsed-time integration, readout, and estimation; support
  multiple detector planes without assuming one area frame
- combine spectral/source planes by index only when their physical grids and
  incoherent-addition semantics are compatible; otherwise retain a bundle or
  use an explicit prepared mapping
- use direct calls within one owner and bounded SPSC ports only at ownership
  boundaries
- keep iceoryx2 and other middleware at explicit process or external-RTC
  boundaries; do not replace immutable in-process epoch sharing or path-local
  workspace ownership with publish/subscribe
- require lock-free release/acquire publication, cache-line-isolated cursors,
  and explicit full-capacity behavior for HIL data-plane rings
- do not use `Base.Channel` as a HIL data-plane port implementation
- select one immutable fidelity provider per prepared acquisition and run,
  preserve one RTC boundary across providers, and record every approximation
- treat calibration as ordinary source/path composition; do not infer physics,
  bypasses, topology, or control ownership from a calibration label
- do not claim raw-pixel throughput from descriptor-only or metadata-only load
  tests
- distinguish execution-clock time, nominal and delivered detector triggers,
  physical exposure boundaries, and reported detector timestamps
- distinguish ring enqueue, semantic admission, physical effective time, and
  terminal outcome
- define full and recovery behavior per resource rather than one generic
  overload policy
- prepare immutable HIL placement; do not migrate work opportunistically
- require predeclared absolute and relative evidence gates for performance
  promotion
- update the relevant subsystem specification and compliance rows rather than
  adding a one-off plan document

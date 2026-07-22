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

## Specification Authority And Normative Language

The subsystem specifications define behavior; the compliance matrix assigns
stable requirement IDs and records implementation and evidence state. In this
specification set, **MUST** denotes a requirement whose absence is a compliance
gap, **SHOULD** denotes the required default unless an alternative has recorded
rationale and equivalent evidence, and **MAY** denotes optional behavior that
remains inside the contract when implemented. Informative examples and
candidate Julia names do not create public API commitments.

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

### Target API position

The target core API models an adaptive-optics plant; it is not an OOPAO object
compatibility layer and it is not HIL-only. HIL, deterministic virtual-time,
offline, and device-resident execution consume the same scientific plant
definitions and prepared model operations. Core plant types therefore SHOULD
use plant, path, acquisition, product, and preparation terminology rather than
an `HIL` prefix or the class layout of the original Python port.

The API separates four semantic roles:

- run-immutable definitions describing topology and physical parameters
- prepared plans binding shapes, algorithms, capacities, backends, and devices
- mutable single-writer model state and workspaces
- caller-owned products passed through explicit mutating operations

The existing `AOSimulation`, `ClosedLoopRuntime`,
`CompositeControllableOptic`, `RuntimeCommandLayout`, and generic `slopes`
surface are transitional characterization/oracle APIs. New HIL architecture
MUST NOT extend those abstractions merely to preserve source compatibility.
Their replacements are introduced gate by gate, and superseded surfaces are
removed once their numerical oracles have migrated.

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
| Optical planes and workspaces | `TelescopeAperture` owns revisioned geometry and intensity reflectivity; every maintained WFS, science, calibration, atmosphere, and controllable-optic path uses caller- or runtime-owned `PupilFunction`, `ElectricField`, and `IntensityMap` products, and the telescope owns no OPD, cadence, exposure, or focal result. The serial event loop independently schedules prepared paths and acquisitions over these products | Preserve this ownership boundary through command-responsive, placed, and paced HIL execution |
| Time and radiometry | Timed atmospheres receive explicit absolute model time. Optical formation produces declared photon-rate or dimensionless products without elapsed time, and exact event lifecycles apply each detector interval once. The serial event loop uses integer-nanosecond schedules and advances one shared atmosphere once for every due optical timestamp | Add command-effective timing and injected execution-clock pacing while preserving the separate plant-time, optical-rate, detector-integration, and execution-lateness contracts |
| Atmosphere and source rendering | One writer publishes an explicit current-state atmosphere epoch token; path-local prepared renderers freeze NGS/LGS/spectral/extended-source descriptions and write caller-owned products without shared render/cache state | Materialize every due path input before advancing mutable atmosphere state, or bind a model-specific retained state snapshot when cross-timestamp rendering is required |
| WFS composition | Shack–Hartmann separates microlens optics, layout/calibration, acquisition, and estimation. Pyramid and BioEdge expose physically distinct focal-plane front ends over shared prepared modulation. Zernike separates its phase spot and pupil relay; Curvature separates its two branch planes and supports independent or packed detector acquisition. LiFT separately prepares its focal-plane forward model, explicit observation domain, and iterative estimator without detector ownership or `AbstractWFS` inheritance | Preserve the completed Gate 0 staged contract while composing independently scheduled acquisitions |
| Shared plant and paths | One prepared telescope/atmosphere feeds reusable path owners and independent acquisition owners. The serial event loop composes unequal science/NGS/LGS schedules, shared-result fan-out, exact atmosphere-epoch materialization, and stable owner RNGs without retaining run-length event lists | Add command-responsive paths, prepared execution groups, and explicit CPU/GPU placement without weakening the serial oracle |
| Detector acquisition | Prepared exact-time global-shutter, rolling-exposure/global-reset row-band, frame-transfer EMCCD, and evolving-charge nondestructive-read lifecycles preserve the common detector pipeline and complete-product readiness. Direct `capture!` retains a separately labeled frame-step/post-exposure convenience | Carry the same complete-product boundary into leases and HIL acquisition ports; keep progressive delivery in user transport rather than creating fragment events in core |
| Guide stars and AO modes | NGS, finite-height LGS, elongation, tomography, co-conjugated additive optics | Mixed NGS/LGS scheduling, independent co-conjugated devices, altitude-conjugated MCAO, and path-specific MOAO |
| Controllable optics | One runtime optic or packed composite at the pupil | Named placed optics with independent command endpoints, effective times, visibility, and bounded late/future policies |
| Science optics and NCPA | Prepared direct imaging owns explicit pupil/field/output/workspace products, supports same-grid incoherent composition and spectral bundles, and reuses one rate image across independent detector acquisitions; sampled NCPA primitives and an asynchronous prepared PROPER handoff remain | Caller-owned native or explicitly converted external photon-arrival-rate products feed independently scheduled detector acquisition; later path composition adds explicit NCPA visibility and external-executor scheduling |
| Product fidelity | One schedule-free prepared acquisition seam now binds full-optical, command-responsive reduced-order extension, or nonresponsive unchanged/copy/bounded-replay providers to an invariant caller-owned product contract. Core validates metadata and result identity and bypasses otherwise unused full-optical path execution; all declared topology is still prepared, and only the boundary plus a test reduced-order causality fixture are claimed | Carry the same product contract into scheduled descriptors, leases, ports, and overload behavior; add validated production reduced-order envelopes and production-shaped load evidence at their assigned gates |
| Calibration illumination | A schedule-free typed path-entry seam binds native or user evaluators, declared visibility/combination, explicit epoch time, stable RNG ownership, and caller-owned pupil/field/intensity/external/detector-input products to ordinary path and acquisition execution without instrument assumptions | Carry the same seam into independent trigger scheduling and optional user-owned setpoint commands without adding a core calibration mode, instrument profiles, or transport policy |
| CPU and accelerators | Direct CPU execution and one implicit backend/device per plant; the deterministic multi-rate event composition is currently a serial CPU oracle, while direct detector lifecycles retain optional-device residency checks | Prepared path owners plus fully explicit or constrained deterministic CPU/GPU and multi-GPU placement; fully automatic cost optimization is deferred |
| HIL boundary | Application-specific proving-ground integration | Transport-neutral command and complete-product acquisition ports, bounded pools, lifecycle, and resource-specific overload behavior |
| Trigger, clock, and scheduling | Integer plant time, fixed-capacity next-event scheduling, modeled finite trigger fan-out/faults, and mapping of delivered edges to detector starts are implemented without wall-clock ownership | Add separate versioned `Clocks.jl` execution pacing/domain mapping and command-ingress timing without moving modeled trigger faults into the execution clock |
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
  atmosphere state explicitly once per epoch, hold its mutable layer state
  through due-path materialization, and render directions through prepared
  path-local source/geometry workspaces
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
  boundaries; do not replace direct current-epoch materialization,
  caller-owned path products, or path-local workspace ownership with
  publish/subscribe
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
- distinguish ring enqueue, core validation/admission, physical effective time,
  core model disposition, and the correlated HIL terminal outcome
- define full and recovery behavior per resource rather than one generic
  overload policy
- prepare immutable HIL placement; do not migrate work opportunistically
- require predeclared absolute and relative evidence gates for performance
  promotion
- update the relevant subsystem specification and compliance rows rather than
  adding a one-off plan document

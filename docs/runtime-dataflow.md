# Runtime Dataflow Guide

Status: active

## Purpose

This guide describes the maintained execution and ownership model. The
canonical HIL-neutral runtime is `AdaptiveOpticsSim.Plant`; application and
model-specific code composes the independent numerical optics and control
primitives around that boundary.

For package structure, see
[`maintainer-architecture.md`](maintainer-architecture.md). For normative HIL
requirements, see
[`hil-package-boundary.md`](hil-package-boundary.md) and
[`hil/compliance-matrix.md`](hil/compliance-matrix.md).

## Ownership Layers

The runtime distinguishes six roles:

1. Cold definitions describe stable identities, optical topology, acquisition
   topology, command schemas, timing policies, and model choices.
2. Run-immutable plans describe validated numerical or physical execution,
   including mappings, coefficients, compatibility, and backend requirements.
3. Persistent mutable state has a single writer and affects later scientific
   results. It includes the virtual-time scheduler, command and detector
   lifecycles, controllable-optic state, and product-readiness state.
4. Replaceable workspaces own scratch and execution resources. Recreating one
   cannot change the deterministic physical trajectory.
5. Products are explicit caller-visible values such as a `PupilFunction`,
   `IntensityMap`, `WFSObservation`, `WFSMeasurement`, or
   `AcquisitionProducts`; they are not scratch.
6. A prepared execution owner validates and binds the exact concrete plan,
   state, workspace, products, and backend/device/context for one writer.

Preparation may construct all six run-owned roles, but construction time does
not make every result a plan. Likewise, mutability alone does not distinguish
state from workspace: discarding state can change the next result, while a
workspace is replaceable without changing the trajectory. The package does
not retain direct Julia `Memory` use in the target architecture; issue #225
migrates the current representation. Prepared owners enforce final capacity
and armed-state non-growth over FixedSizeArrays.jl fixed-size arrays with
concrete element types, tuples, purpose-built owners, or backend arrays. A
resizable vector is limited to cold construction before sealing. Preparation
does not return a field or collection with stored type `Any`; external input is
normalized into concrete family ownership before execution.

The telescope owns revisioned aperture geometry, intensity reflectivity,
diameter, and spatial sampling. It owns neither a path's mutable OPD or electric
field nor detector cadence, exposure, FFT scratch, or atmosphere model time.

## Definition And Preparation

A plant starts from `PlantDefinition`. It may declare:

- one telescope and atmosphere
- named optical paths
- independent acquisition endpoints
- independent controllable optics
- native sampled `NCPA` or `OPDMap` aberrations with explicit path visibility
- one or more command schemas per optic

`prepare_plant` requires an explicit run seed and, for every command endpoint,
a `CommandEndpointConfiguration`. Preparation:

- canonicalizes stable identities
- makes run-owned exact-target sampled-OPD copies from supported caller input
- resolves bounded path-local sampled-aberration and controllable-optic
  couplings
- prepares path and acquisition implementations
- retains one exact prepared device execution context in the plant, its path
  executors, and their acquisition owners
- allocates backend- and device-compatible destinations
- copies initial and optional safe commands into endpoint-owned storage
- constructs independent RNG streams from stable owner identities
- validates all run-immutable bindings

Successful public preparation and execution returns complete the retained
backend stream before restoring the caller's device and stream selection.
Worker-visible independent path-group calls establish that completion before
publishing their ready or complete phase; coordinator-local serial execution
uses one enclosing completion boundary.

Preparation may allocate and fail. Repeated execution mutates only prepared
state and caller-owned products. Public schedule-free
`materialize_path_input!`, `execute_path!`, and `execute_acquisition!` calls
enter the retained context, including its accelerator stream where applicable,
and restore the caller's previous device and stream selections after success or
an exception. The schedule-free `execute_acquisition_selection!` and
`execute_acquisition_selection_at!` entry points enter the plant context once
and execute their lower-level operations inside it.
Model-specific execution dispatches run inside that boundary and must not
select an alternate device or stream.

## Controller Output Routing

An RTC or controller often computes one flat vector even though the plant owns
several independently timed command endpoints. The maintained bridge is
prepared controller-output routing:

```julia
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Plant

controller_buffer = zeros(Float32, 5)
products = (
    woofer=@view(controller_buffer[1:2]),
    tweeter=@view(controller_buffer[3:5]),
)

routing = prepare_controller_output_routing(
    plant,
    products,
    ControllerOutputRoute(:woofer, :woofer_command),
    ControllerOutputRoute(:tweeter, :tweeter_command),
)
```

Each named product is borrowed without packing or copying. Preparation requires
an exact match with its prepared endpoint's numeric type, shape, array backend,
and compute device. Scalar products use assigned `Ref` storage.

Routing owns no sequence number, effective timestamp, admission, transaction,
queue, transport, or optical-grouping semantics. Integration code obtains a
payload with `controller_output_payload` and constructs a separate
`PlantCommand` for each due endpoint. Plant-command admission then copies the
presented payload into that endpoint's bounded storage. Co-conjugated optics,
flat-buffer adjacency, and equal timestamps do not imply atomic application;
explicit `PlantCommandTransaction` membership does.

## Virtual-Time Event Flow

`prepare_plant_event_loop` composes a finite prepared plant with:

- command endpoint generators
- periodic optical samples
- periodic or delivered-trigger acquisition starts
- detector lifecycle transitions
- one shared atmosphere timeline
- sampled-aberration application
- controllable-optic application

`step_plant_events!` processes all work at the next canonical
`PlantTimestamp`; `run_plant_events_until!` advances through all due timestamps
up to a limit. Equal-time ordering follows stable causal phase, prepared owner
ordinal, and occurrence—not task completion or tuple iteration.
`PreparedPlantEventLoop` retains the exact prepared plant device execution
context. Both public execution calls enter that context and restore the caller's
previous device and stream selections; their lower-level optical work reuses
the already active context.

At a due command timestamp the loop:

1. claims the next application-ready command
2. applies endpoint value policy and stages physical-optic state
3. commits the visible physical state
4. records a bounded disposition

Commands due at the same timestamp are applied before optical samples. Each
endpoint retains independent sequence, effective-time, silence, and capacity
policy. The loop owns plant time only; wall-clock pacing, external timestamp
mapping, payload leases, transport, and RTC protocol belong to the HIL
integration layer.

For every due full-optical path, the serial order is atmosphere
materialization, sampled aberrations, controllable surfaces, autonomous
path-local optics, and then the typed path executor. Replacement sampled
aberrations precede canonical placement/identity-ordered additives.
Preparation rejects more than one replacement visible on a path.

## Optical Sample Flow

For each due full-optical path, execution:

1. advances the shared atmosphere to explicit model time
2. renders the published epoch through a path-local prepared direction
   renderer into caller-owned path storage
3. applies the visible static and controllable surfaces
4. executes the prepared optical model into its exact path result
5. executes each due acquisition provider into its own
   `AcquisitionProducts`

Preparation resolves required optic placement and all-path or selected-path
visibility into canonical bounded ranges and compatible path-coupling groups.
Full-optical execution applies only a path's visible surfaces additively.
Sampled surface models may bind the same-grid fast path or a finite-support
NGS/LGS atmospheric-conjugate transform; unsupported model/placement
combinations fail during preparation. No packed or composite optic object
couples device cadence.

Native direct science uses `prepare_direct_imaging` and
`form_direct_image!`. Prepared PROPER integration meets native optics at the
same declared photon-arrival-rate/acquisition boundary. NGS and LGS WFS paths,
science cameras, calibration paths, NCPA branches, and coronagraph paths remain
separate declared paths when their propagation differs.

## WFS And Detector Flow

The maintained WFS stage contract is:

1. `form_wfs_optical_products!` produces one or more detector-facing
   photon-arrival-rate products.
2. `acquire_wfs_observation!` applies detector acquisition into a typed
   `WFSObservation`.
3. `estimate_wfs_measurement!` writes a typed `WFSMeasurement`.

The optical front end, detector acquisition, and estimator are independent
components. For example, a Shack–Hartmann optical front end composes a
`MicrolensArray`; it does not make the microlens optic part of detector or
centroid state.

Prepared frame acquisition applies operations in this order:

1. validate the photon-arrival-rate product and its spatial measure
2. apply presampling detector response on the optical grid
3. integrate represented cells into physical pixels
4. apply wavelength-dependent QE and explicit exposure duration exactly once
5. apply charge-domain effects, stochastic response, binning, and readout
6. publish a complete product at the detector lifecycle's readiness event

Global shutter, rolling shutter, frame transfer, nondestructive reads, and
up-the-ramp sampling retain separate virtual-time lifecycle definitions.
Exposure, readout completion, product readiness, and re-arming are distinct
events.

## Acquisition Fidelity

Every prepared acquisition binds one run-immutable provider:

- full optical
- command-responsive reduced order
- synthetic or bounded replay

All provider styles write the same prepared logical product contract. This
allows one plant to combine high-fidelity paths with fast causal surrogates for
RTC throughput and latency testing. Fidelity does not change during a run; a
different provider requires another preparation.

Calibration illumination uses the ordinary path and acquisition machinery.
The user or companion model declares where it enters, which paths see it, and
how its state evolves. Core does not infer a special calibration mode from a
label.

## Explicit Model Loops

Package examples such as the Subaru AO188/AO3k model may expose their own
`prepare!`, `step!`, and `readout` functions. These are explicit model
compositions used for numerical examples and benchmarks; they are not a second
generic runtime API. A simple closed loop can likewise compose:

- `advance_by!` and `render_atmosphere!`
- `set_command!` and `apply!`
- the three WFS stages
- `AdaptiveOpticsSim.Control.reconstruct!`
- `AdaptiveOpticsSim.Control.VectorDelayLine` and
  `AdaptiveOpticsSim.Control.DiscreteIntegratorController`

This keeps reusable numerical/control primitives independent of HIL scheduling.

## Backend And Parallel Execution

Prepared storage is parameterized by numeric type, backend, and compute
device. Optical and detector kernels reject implicit host/device mixing and
avoid GPU scalar indexing. Same-owner stages call directly; future ownership
boundaries may exchange bounded descriptors.

`SimulationEnsemble` remains a coarse-grained facility for independent model
runs, sources, time-series sweeps, or offline studies. It is not the
deadline-path HIL scheduler. Sequential execution is the default; Julia
threads, AcceleratedKernels, Dagger, and backend streams are explicit policies
whose use must avoid nested oversubscription.

## Validation And Performance Evidence

Functional support is recorded in
[`model-validity-matrix.md`](model-validity-matrix.md). Production-surface
qualification is recorded in
[`supported-production-surfaces.md`](supported-production-surfaces.md).

Benchmark workloads under `benchmarks/` own their explicit loop composition;
they do not introduce a package runtime abstraction. Representative performance
artifacts live under `benchmarks/results/` and remain outside ordinary
`Pkg.test()`.

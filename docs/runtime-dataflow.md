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

The runtime is intentionally split into four layers:

1. Cold definitions describe stable identities, optical topology, acquisition
   topology, command schemas, timing policies, and model choices.
2. Preparation validates the complete binding and constructs typed plans,
   backend-resident storage, bounded capacities, workspaces, and deterministic
   RNG owners.
3. Mutable state has a single writer. The prepared plant event loop owns the
   virtual-time scheduler, command lifecycle, detector lifecycle, controllable-
   optic state, and product-readiness state.
4. Products are explicit caller-visible values such as a `PupilFunction`,
   `IntensityMap`, `WFSObservation`, `WFSMeasurement`, or
   `AcquisitionProducts`.

The telescope owns revisioned aperture geometry, intensity reflectivity,
diameter, and spatial sampling. It owns neither a path's mutable OPD or electric
field nor detector cadence, exposure, FFT scratch, or atmosphere model time.

## Definition And Preparation

A plant starts from `PlantDefinition`. It may declare:

- one telescope and atmosphere
- named optical paths
- independent acquisition endpoints
- independent controllable optics
- one or more command schemas per optic

`prepare_plant` requires an explicit run seed and, for every command endpoint,
a `CommandEndpointConfiguration`. Preparation:

- canonicalizes stable identities
- prepares path and acquisition implementations
- allocates backend- and device-compatible destinations
- copies initial and optional safe commands into endpoint-owned storage
- constructs independent RNG streams from stable owner identities
- validates all run-immutable bindings

Preparation may allocate and fail. Repeated execution mutates only prepared
state and caller-owned products.

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
and physical device. Scalar products use assigned `Ref` storage.

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
- controllable-optic application

`step_plant_events!` processes all work at the next canonical
`PlantTimestamp`; `run_plant_events_until!` advances through all due timestamps
up to a limit. Equal-time ordering follows stable causal phase, prepared owner
ordinal, and occurrence—not task completion or tuple iteration.

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

## Optical Sample Flow

For each due full-optical path, execution:

1. advances the shared atmosphere to explicit model time
2. renders the published epoch through a path-local prepared direction
   renderer into caller-owned path storage
3. applies the visible static and controllable surfaces
4. executes the prepared optical model into its exact path result
5. executes each due acquisition provider into its own
   `AcquisitionProducts`

The current common controllable-optic execution applies independently
commanded co-conjugated surfaces additively. Explicit altitude placement and
path visibility are later HIL gates; no packed or composite optic object
couples their cadence.

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
4. apply wavelength-dependent QE and explicit exposure time exactly once
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
- `reconstruct!`
- `VectorDelayLine` and `DiscreteIntegratorController`

This keeps reusable numerical/control primitives independent of HIL scheduling.

## Backend And Parallel Execution

Prepared storage is parameterized by numeric type, backend, and physical
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

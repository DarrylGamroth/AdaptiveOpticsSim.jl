# Runtime Dataflow Guide

Status: active

## Purpose

This document explains the end-to-end runtime and dataflow of the maintained
simulation stack.

It is not a symbol reference. It is the operational picture of:

- what gets built
- what gets prepared
- what mutates each step
- which products are exported
- where validation and benchmark evidence attach

For the broader package/platform picture, start with
[`maintainer-architecture.md`](maintainer-architecture.md).

## Main Runtime Objects

The most important runtime-facing objects are:

- `Telescope`
- atmosphere model such as `MultiLayerAtmosphere` or
  `InfiniteMultiLayerAtmosphere`
- one or more WFS objects
- optional `Detector` objects
- `DeformableMirror`
- calibration/reconstructor objects
- `ClosedLoopRuntime`

The orchestration layer is in `src/control`.

## Build Phase

The build phase constructs long-lived simulation objects:

1. create telescope and source geometry
2. create atmosphere model and persistent screen state
3. create WFS models and detectors
4. create DM and calibration/reconstruction objects
5. build `ClosedLoopRuntime` with:
   - runtime profile
   - latency model
   - product requirements
   - execution policy

Outputs of this phase should be stable objects with owned buffers, not
ad hoc temporary arrays.

## Prepare Phase

`prepare!` and subsystem-specific preparation do the work that should not
happen inside the hot loop:

- allocate exported readout and metadata surfaces
- precompute calibration products
- prepare grouped execution stacks
- bind each acquisition to a full-optical, reduced-order, or synthetic/replay
  product provider with preallocated state
- prepare any user-declared illumination evaluator at its supported typed path
  entry or detector input
- freeze mutable source/profile inputs while preserving intentional source
  identity across shared optical arms
- prepare delay lines and runtime staging
- ensure detector and WFS pipeline buffers exist

The design goal is:

- construction is allowed to be heavier
- the step path should only mutate prepared state

## Per-Step Flow

At a high level, one closed-loop step is:

1. advance atmosphere
2. reset or update telescope OPD state
3. apply DM command into telescope phase
4. propagate source/field/atmosphere state to the sensing surface
5. run WFS measurement and optional detector readout
6. extract required products:
   - slopes
   - WFS pixels
   - science pixels
7. run reconstruction/controller update
8. stage delayed commands or delayed readout products as needed

This split is visible in runtime timing and benchmark surfaces.

This is the current implemented single-step flow. Atmosphere evolution and
direction rendering already use the first two explicit ownership boundaries
below; the remaining Gate 0/2 work replaces shared telescope path scratch and
finishes the full prepared flow:

1. advance one shared atmosphere epoch to explicit model time
2. render each due direction through a path-local prepared renderer into a
   caller-owned pupil function or electric field
3. apply visible static and controllable surfaces to that path product
4. form a WFS or direct-science photon-arrival-rate `IntensityMap` through a
   prepared optical front end
5. integrate each compatible detector's explicit sample/exposure duration and
   acquire independently from the immutable arrival-rate product
6. estimate a WFS measurement where the endpoint requires one

The telescope owns aperture, reflectivity, spatial sampling, and geometry; it
does not own cadence/exposure duration, the path OPD/field, PSF result, FFT
plans, or propagation scratch. Atmosphere advancement receives explicit model
time rather than telescope sampling time. Source geometry and extended-source
expansion are frozen or prepared per path rather than mutated in the shared
atmosphere or source definition.

Every prepared plane declares coordinate domain and sampling,
centering/orientation,
wavelength/channel, units/normalization, coherence/combination policy, backend,
device compatibility, and whether values are spatial densities or integrated
over represented cells. A physical detector-facing product is either photon
irradiance (photons·s⁻¹·m⁻²) or a cell-integrated photon rate
(photons·s⁻¹); detector acquisition applies elapsed time exactly once and
preserves the declared spatial measure through response and pixel integration.
Incompatible spectral grids remain a bundle or require an explicit prepared
mapping. Native direct science and prepared PROPER output meet at this same
caller-owned photon-arrival-rate/acquisition boundary. These are concrete
products and functions, not a universal optical graph or resampling framework.

The target then separates reusable optical paths from independently scheduled
or triggered acquisition state, and schedules exposure, optical sampling,
readout, and publication as separate events. Commands arrive from an external
RTC and become effective independently of detector cadence.
Acquisitions may follow independent schedules or delivered edges from a common
trigger source with per-detector delay, skew, jitter, and explicit dropped or
duplicate-edge faults. Physical exposure timing, reported detector timestamps,
and HIL execution-clock lateness remain distinct.

That target admits heterogeneous NGS/LGS WFS paths, direct science cameras,
PROPER-backed coronagraph paths, common MCAO planes, and path-specific MOAO
planes while preserving independently commanded co-conjugated optics. Sampled
NCPA belongs to the selected native branch; detailed instrument and
coronagraph propagation remains at the prepared PROPER seam.

The target also decouples acquisition scheduling and publication from product
generation fidelity. A prepared acquisition may use the full optical path, a
command-responsive geometric/linear/reduced-resolution provider, or a
preallocated synthetic/replay source. All providers produce the same declared
shape, type, geometry/radiometry, metadata, sequence, complete-product lease,
and port behavior, so
an RTC adapter can be tested at production rate without paying for optics that
are outside the test boundary. Fidelity is mixed per acquisition and never
changes during a run; another fidelity tier requires another prepare/arm cycle.

The reduced-order provider remains a causal AO plant: it advances a seeded or
replayed time-correlated disturbance, projects it into each sensing direction,
subtracts the response of commands that are physically effective at the sample
time, and emits calibrated slopes or approximate pixels with selected noise and
detector timing. The external RTC still performs its normal reconstruction,
tomography, controller, command-splitting, and recovery work; only the expensive
optical evaluation is replaced by a validated surrogate.

Calibration illumination follows the same path machinery rather than entering
a special runtime mode. User or companion code declares the supported entry
boundary, downstream visibility, source evaluator, state timing, and any
combination rule. Core executes that prepared evaluator and the ordinary
acquisition pipeline; a calibration label alone never changes propagation.

The current `AOSimulation` stores one controllable optic, so existing
multi-surface plants use `CompositeControllableOptic` as a packed-command and
additive-application adapter. The breaking HIL refactor removes that type and
registers each optic and independently timed command endpoint explicitly, then
derives co-located optical execution groups during preparation. Packed transport
schemas belong to user integration outside the general HIL package. They map to
canonical command transactions submitted through the HIL ports and carry no
physical grouping, clock, or atomicity semantics of their own.

Each canonical command endpoint prepares its payload type/shape, units, basis
and calibration revision, absolute or incremental semantics, limits, session
epoch, sequence behavior, and silence/watchdog policy. Sampled actuator or
device feedback returns through an ordinary acquisition endpoint rather than
being conflated with the terminal outcome of a command.

The maintained target architecture is indexed by
[`hil-package-boundary.md`](hil-package-boundary.md); durable capability IDs,
states, and acceptance gates are tracked in
[`hil/compliance-matrix.md`](hil/compliance-matrix.md).

## Product Ownership

The package now makes product ownership more explicit than before.

Important distinctions:

- scratch buffers are not exported outputs
- sampled detector inputs are not the same thing as detector outputs
- grouped WFS intermediate stacks are not the same thing as archived or
  exported readouts

For example, in diffractive Shack-Hartmann execution there is now a clearer
boundary between:

- sampled pre-detector spot stacks
- post-detector signal stacks
- exported runtime pixel outputs

That separation matters for:

- correctness
- GPU layout ownership
- benchmarked allocation behavior

## Detector and WFS Dataflow

The maintained detector/WFS pipeline now follows a more explicit shape:

1. produce sampled optical signal
2. apply detector/readout pipeline if configured
3. reduce or extract slopes/signals
4. snapshot only the runtime outputs that were requested

The runtime output plan decides whether a given simulation step must produce:

- slopes only
- slopes plus WFS pixels
- science pixels
- metadata surfaces

For grouped composite execution, exported outputs can additionally include:

- per-branch grouped WFS frames
- per-branch grouped science frames
- compatible-shape grouped WFS stacks
- compatible-shape grouped science stacks

This keeps slopes-only runs from paying unnecessary export costs.

For synthetic load tests, the output plan also declares whether payloads are
reused, touched, generated, copied, replayed, and consumed. Descriptor-only
tests measure descriptor throughput; production-shaped pixel tests retain the
representative payload and memory traffic needed for an RTC pixel-processing
claim.

## Backend Flow

The intended backend execution model is:

- model code calls shared backend services
- shared services select CPU / CUDA / AMDGPU behavior
- optional backend extensions refine behavior without rewriting the public API

Backend validation is attached at three levels:

- optional backend test coverage
- smoke scripts
- benchmark evidence

Independent runtime plants may be collected into `SimulationEnsemble` for
coarse sweeps. Sequential execution remains the default. Julia threads,
AcceleratedKernels task partitioning, and Dagger task graphs are explicit
policies above each runtime; they do not replace the runtime's CPU-HIL or
device-resident execution plan. The current CPU HIL baseline calls one plant
directly. The target HIL companion drives one immutable prepared plant plan
through its owned CPU/GPU agents and canonical ports; it does not schedule that
deadline path as a `SimulationEnsemble`.

Within one HIL process, same-owner stages call directly and owner changes use
bounded SPSC descriptors. Iceoryx2 or another middleware is reserved for a
deliberate process or external-RTC boundary; it is not the fan-out mechanism
for ordinary in-process optical arms.

Representative performance evidence is intentionally kept out of `Pkg.test()`.

## Benchmarks and Evidence

There are two relevant evidence layers:

### Functional/validation evidence

- [`model-validity-matrix.md`](model-validity-matrix.md)
- frozen OOPAO and SPECULA reference bundles

### Runtime/engineering evidence

- benchmark artifacts under `benchmarks/results/`
- profile scripts under `scripts/`

## How to Use This Guide

When debugging or extending runtime behavior:

1. use this doc to identify the layer that owns the behavior
2. use [`maintainer-architecture.md`](maintainer-architecture.md) for the
   package-level composition picture
3. use [`extension-guide.md`](extension-guide.md) for new subsystem families
4. use [`api-reference.md`](api-reference.md) for specific symbols

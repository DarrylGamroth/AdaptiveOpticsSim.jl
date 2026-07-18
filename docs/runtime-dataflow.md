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
- optional prepared `DetectorAcquisitionPlan` objects for typed intensity maps
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
   - explicit positive `atmosphere_step`
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
- prepare compatible incoherent intensity accumulation or retain incompatible
  products in an `OpticalProductBundle`
- validate detector-facing map metadata and size acquisition buffers with
  `prepare_detector_acquisition`

The design goal is:

- construction is allowed to be heavier
- the step path should only mutate prepared state

## Per-Step Flow

At a high level, one closed-loop step is:

1. advance atmosphere by the runtime's explicit `atmosphere_step`
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

This is the current implemented single-step flow. Atmosphere evolution,
direction rendering, rate formation, and the prepared frame-detector boundary
provide the corresponding ownership foundations below; the remaining Gate 0/2
work replaces shared telescope/WFS path state and finishes the fully composed
prepared flow:

1. advance one shared atmosphere epoch to explicit model time
2. render each due direction through a path-local prepared renderer into a
   caller-owned pupil function or electric field
3. apply visible static and controllable surfaces to that path product
4. call `form_wfs_optical_products!` to form one WFS photon-arrival-rate
   `IntensityMap`, a concrete tuple of them, or an `OpticalProductBundle`
   through a prepared optical front end
5. call `acquire_wfs_observation!` to integrate each compatible detector's
   explicit sample/exposure duration independently without mutating the
   prepared arrival-rate storage or its immutable metadata/binding
6. call `estimate_wfs_measurement!` to write a typed `WFSMeasurement` where the
   endpoint requires one

The telescope owns aperture, reflectivity, spatial sampling, and geometry; it
does not own cadence/exposure duration, a direct-science path's OPD/field or
focal-plane result, FFT plans, or propagation scratch. Transitional telescope
OPD remains for WFS families awaiting decomposition. Atmosphere advancement receives explicit model
time; the telescope has no timing property. Source geometry and extended-source
expansion are frozen or prepared per path rather than mutated in the shared
atmosphere or source definition.

Every prepared plane declares coordinate domain and sampling,
centering/orientation,
wavelength/channel, units/normalization, coherence/combination policy, backend,
device compatibility, and whether values are spatial densities or integrated
over represented cells. A physical detector-facing product is either photon
irradiance (photons·s⁻¹·m⁻²) or a cell-integrated photon rate
(photons·s⁻¹). Prepared detector acquisition applies presampling response on
the optical grid, integrates spatial-density samples over represented cells
and then physical pixels, and applies QE and elapsed time exactly once.
Incompatible spectral grids remain a bundle or require an explicit prepared
mapping. Native direct science uses `prepare_direct_imaging` and
`form_direct_image!`: same-wavelength asterism components share a compatible
incoherent output, differing spectral grids remain an `OpticalProductBundle`,
and extended sources are explicitly expanded through
`extended_source_asterism`. Its off-axis placement is resolved to a finite
integer shift during preparation from the grid's declared `:x`/`:y` axis order
and signs; it remains periodic, not subpixel interpolation or finite-field
loss. Runtime-length asterism, extended-source, and spectral expansions use
vector-backed prepared components, so quadrature length does not become part
of a recursively specialized method type. Each native leaf still owns its FFT
workspace; a future shared-propagation or convolution optimization must first
declare the source-dependent pupil-field equivalence it relies on. Prepared PROPER output
meets native imaging at the same caller-owned photon-arrival-rate/acquisition boundary. These are concrete
products and functions, not a universal optical graph or resampling framework.

The transitional `SharedOpticalRuntime` forms each arm's native rate image once
and captures its detector tuple serially. Detector state and exposure are
independent, but stochastic draws currently come from one runtime RNG in tuple
order. The later event runtime may assign stable per-endpoint streams when
order-invariant detector noise is required. Before advancing the atmosphere,
the current runtime preflights every science detector's exact prepared binding
and whole-exposure idle state. Detector acquisition preparation has already
sized conventional multi-read products and rejected predictable shape and
up-the-ramp schedule errors before changing detector storage. Unexpected
backend, kernel, RNG, or concurrent external failures remain fail-stop rather
than transactional rollback.

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

For example, the prepared diffractive Shack–Hartmann front end accumulates a
compatible same-wavelength asterism into one optical-rate mosaic and retains
distinct wavelength grids as native-sampling leaves in an
`OpticalProductBundle`. It does not index-add or implicitly resample those
leaves; acquisition consumes each through its declared detector mapping.
Extended-source expansion and the legacy single-product spectral convenience
path retain their documented common-grid behavior. Frozen references for the
retired reference-wavelength index-grid approximation are test-only
characterization data, not a production capability. There is now a clearer
boundary between:

- pre-detector optical-rate spot stacks
- post-detector signal stacks
- exported runtime pixel outputs

That separation matters for:

- correctness
- GPU layout ownership
- benchmarked allocation behavior

## Detector and WFS Dataflow

The prepared frame-detector path now follows this order:

1. produce a photon-arrival-rate intensity map or an explicitly normalized map
   with a prepared physical scale
2. apply the presampling detector response on the optical grid
3. integrate represented cells into physical pixels
4. apply wavelength-channel QE and explicit whole or incremental exposure time
   exactly once
5. apply charge-domain effects, stochastic response, binning/readout, and
   publication
6. reduce or estimate WFS signals where required
7. snapshot only the runtime outputs that were requested

`WFSObservation`, `WFSMeasurement`, and the generic prepared
formation/acquisition/estimation protocols now establish this static boundary.
Most maintained WFS state types still compose several stages internally; their
adoption of the contract remains later Gate 0 family work. Raw matrix detector
entry remains a documented legacy cell-integrated-rate path, while
`DetectorAcquisitionPlan` is the metadata-validated prepared frame boundary.

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

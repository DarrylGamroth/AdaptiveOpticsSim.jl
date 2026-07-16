# Extension Guide

Status: active

This guide describes the maintained extension seams for adding new detector
families, wavefront sensors, deformable mirrors, controllable optics,
controllers, and reconstructors.

Prefer small concrete types plus multiple dispatch. Do not add central
type-switching registries for new families.

## Source Layout

Use lower-case directory names for new source-tree locations. Julia type names
may use CamelCase, but package directories should remain lower-case subsystem
names such as `src/wfs`, `src/detectors`, `src/optics`, or `src/control`.

## Detectors

Physical detector families live in `src/detectors/`.

Use this split:

- immutable parameter/configuration structs for physical constants and readout
  configuration
- mutable state structs for accumulated charge, readout buffers, and thermal or
  persistence state
- `!` methods for capture, integration, correction, and readout hot paths
- shared layers for reusable behavior such as quantization, frame response,
  multi-read readout, thermal state, and counting statistics

New detector families should implement the smallest family-owned methods needed
to connect to the generic detector pipeline. If a feature can apply to multiple
families, put it in a shared layer instead of copying it into each sensor.

A custom frame-response model is not automatically eligible for
detector-coupled WFS reference calibration. Its extension must provide an
instance-complete
`AdaptiveOpticsSim.detector_response_calibration_signature(model, seed)`
method covering every parameter that can change the deterministic response.
Without that method, WFS calibration fails closed rather than caching a
reference by response type alone.

Use `bits` for quantization depth and `output_type` for the Julia element type
returned at a HIL/RTC boundary.

## Wavefront Sensors

WFS families live in `src/wfs/`. Family-specific implementation should stay
near the family, usually in the corresponding directory under `src/wfs/`.

New physical WFS work should first implement the prepared semantic stages:

- `prepare_wfs_optical_formation(model, input, output)` validates an explicit
  caller-owned pupil function or electric field and one or more
  detector-plane photon-arrival-rate outputs; the returned concrete plan is
  executed by `form_wfs_optical_products!`
- `prepare_wfs_acquisition(model, optical_products, observation)` validates
  detector mappings, explicit duration, backend/device placement, and one or
  more caller-owned `WFSObservation` destinations; the returned concrete plan
  is executed by `acquire_wfs_observation!` with an explicit caller-owned RNG
- `prepare_wfs_estimation(model, observation, measurement)` validates the
  estimator and caller-owned `WFSMeasurement`; the returned concrete plan is
  executed by `estimate_wfs_measurement!` and declares
  `AcquiredObservationPath()` through `wfs_measurement_path`

A geometric or reduced-order provider may instead prepare estimation directly
from a `PupilFunction` or pupil-plane `ElectricField`. It must declare
`DirectMeasurementPath()` and must not create unused rate, detector, or
observation storage. `WFSObservation` supports scalar `Ref` storage and arrays
of any rank, as does `WFSMeasurement`; use concrete tuples for multiple
observations. Preserve incompatible spectral or branch rate products in
`OpticalProductBundle`.

Extensions should call the qualified validation seams
`AdaptiveOpticsSim.validate_wfs_optical_input`,
`AdaptiveOpticsSim.validate_wfs_optical_products`,
`AdaptiveOpticsSim.validate_wfs_observation` or
`AdaptiveOpticsSim.validate_wfs_observations`, and
`AdaptiveOpticsSim.validate_wfs_measurement` as applicable before returning a
prepared plan. These remain qualified extension APIs rather than ordinary
exported workflow names.

Prepared types should contain concrete immutable plans/params and separately
typed single-writer workspace, detector, calibration, and RNG state. Bind exact
array/state identities, validate physical device as well as semantic backend,
and create detector-output aliases or packed views only after detector buffers
are prepared. Repeated execution must not resize, rebuild metadata, query a
device, copy to the host, or select stages through an abstract container.
Raise `WFSPreparationError(stage, reason, msg)` for preparation
incompatibility and for execution-time prepared-binding rejection before any
destination mutation. The protocol stages are `:optical_formation`,
`:acquisition`, and `:estimation`; `reason` is an open extension identifier.

The current maintained WFS families are being migrated incrementally. During
that transition, a legacy `AbstractWFS` family also provides:

- a concrete sensor type subtype of the package WFS abstraction
- setup/precomputation methods owned by the WFS family
- measurement methods that write into preallocated state
- signal extraction methods for slopes, intensities, or family-specific
  products
- detector image formation if the sensor has a maintained detector-facing path

The generic runtime and future family adapters call prepared seams or WFS-owned
legacy seams through dispatch rather than branching on sensor families.

## Runtime Source Roles

`AOSimulation` keeps WFS and science sources as separate typed roles. Extension
code should use `wfs_source(simulation)` and `science_source(simulation)` rather
than assuming one source serves both paths.

Maintained timed atmospheres separate evolution from direction rendering. Their
public execution boundary is:

- `advance_by!(atmosphere, elapsed_seconds; rng)` or
  `advance_to!(atmosphere, model_time; rng)` for the single evolution writer
- `prepare_atmosphere_renderer(atmosphere, telescope, source)` during cold
  preparation for one source direction
- `render_atmosphere!(destination, renderer, atmosphere, epoch)` for warmed,
  caller-owned output

The timed atmosphere must not read a telescope timing value, detector cadence,
wall time, or a renderer to determine elapsed time. Rendering must not evolve
the atmosphere or consume RNG. A stale epoch, renderer from another atmosphere,
or incompatible output must fail before output mutation.

`AdaptiveOpticsSim.AbstractTimedAtmosphere` and its model hooks are an advanced
qualified extension seam. A concrete implementation owns an atmosphere
identity and timeline, implements `initialize_atmosphere!` and
`evolve_atmosphere!`, declares its layer-storage element type through
`atmosphere_numeric_type`, and exposes physical layers whose render methods
accept precomputed direction geometry. Wrappers around a timed atmosphere must
delegate its identity, timeline, layer, numeric-type, and rendering seams.
Source-independent static extensions may stay subtypes of `AbstractAtmosphere`
and implement the ordinary two-argument `propagate!`; the source-sensitive
fallback delegates to it.

Use `prepare_atmosphere_renderers` for `Asterism` and `ExtendedSource`. The
singular preparation call rejects them so a multi-direction source cannot be
silently treated as its reference direction. A custom source that owns mutable
profile data should specialize the qualified `AdaptiveOpticsSim.freeze_source`
seam and return a run-owned copy.

Shared multi-arm execution advances once and stores one prepared atmosphere
renderer per source arm. Runtime construction freezes mutable source inputs and
preserves source identity across arms so identical paths can still be reused. A
custom WFS that needs special pupil preparation when consuming an
already-rendered shared arm should extend
`prepare_shared_runtime_wfs!`; the default is a no-op. The Curvature WFS uses
this seam to establish its measurement path and then restores the normal
source path in `finish_runtime_wfs_sensing!`.

## Deformable Mirrors

DM implementation lives in `src/optics/deformable_mirror.jl`.

Use the composed model:

- topology describes actuator layout and coordinates
- influence basis describes how commands map to sampled surface modes
- actuator behavior describes clipping, dead actuators, slaving, or future
  dynamic behavior
- `DeformableMirror` composes those pieces into the runtime optic

Analytic Gaussian influence models are operator-backed: regular grids use a
factored separable operator and other Gaussian topologies evaluate a fused
matrix-free operator. Do not require `state.modes` to be mutable dense storage.
Code that truly needs the full sampled basis should materialize it during setup,
then keep that result out of repeated HIL application. Explicit dense and
measured influence models continue to own backend-native matrix storage.

Measured manufacturer influence functions should enter as sampled influence
bases with explicit metadata. Do not force measured data through an anonymous
dense matrix unless it is truly only an escape hatch.

## Controllable Optics

Use controllable optics for modal or low-order optical surfaces that are not
better represented as a DM technology family.

Common examples:

- tip/tilt surfaces
- low-order modal optics
- Zernike-mode control surfaces
- calibration or NCPA compensation surfaces

Keep command application explicit and ordered in the plant. Runtime code should
operate on the generic controllable-optic interface.

## Controllers And Reconstructors

Controllers own temporal behavior. Reconstructors own slopes-to-command or
measurement-to-command operators.

Use this split:

- reconstructors expose the command operator and any diagnostics needed by
  calibration or validation
- controllers own integrator state, gains, leakage, saturation, or future
  temporal filters
- runtime code calls generic accessors and `!` methods

Prefer concrete state and workspace fields. Avoid hidden allocations in the
per-step control path.

For large calibration/control surfaces:

- use `AdaptiveOpticsSim.interaction_matrix!` with caller-owned storage when
  the calibration matrix should live in a memory map, shared buffer, or
  backend-native allocation
- use `FactorizedReconstructor` when a validated truncated SVD rank can bound
  runtime storage and compute without materializing the dense inverse
- use `ControlledReconstructor` to compose temporal state with any maintained
  reconstructor while preserving the existing runtime contract

Rank selection is an optical/control validation decision. A benchmark win from
truncation is not evidence that discarded modes are acceptable for a real
plant. The dense `ModalReconstructor` remains the compatibility baseline.

Runtime profiles and execution plans are orthogonal. A profile such as
`HILRuntimeProfile` selects modeled delays and output fidelity. An execution
plan selects storage and synchronization behavior:

- `CPUHILExecutionPlan` is the direct, host-resident low-latency path
- `DeviceResidentExecutionPlan` keeps repeated simulation and control work on
  one accelerator and synchronizes at explicit observation/export boundaries

The runtime chooses the plan from the simulation backend unless one is
provided explicitly. Accelerator reconstructors must implement
`AdaptiveOpticsSim.runtime_reconstructor_storage(reconstructor)` and return a tuple containing
their hot-path control matrices and workspaces. This lets construction reject
mixed host/device control paths before the first step. Return `()` only for a
truly backend-agnostic operator. Use `synchronize_runtime!` before direct host
observation; `SimulationInterface` snapshots already provide an export
boundary.

Stateful reconstructors should additionally implement
`AdaptiveOpticsSim.runtime_reconstructor_ownership_roots(reconstructor)` and
return their mutable workspaces. This prevents two threaded ensemble members
from sharing operator scratch state. Stateful controllers used inside
`ControlledReconstructor` implement `runtime_controller_storage` for backend
residency; controller objects are conservatively treated as mutable ownership
roots by default.

## Ensemble Scheduling

`SimulationEnsemble` schedules independent simulation boundaries at coarse
granularity. Keep it outside an external-RTC HIL loop: the direct scalar CPU
runtime has tighter and more predictable latency than task creation or task
graph scheduling.

Available policies have deliberately different scopes:

- `SequentialExecution` is the default and preserves fixed member order
- `DeterministicExecution` additionally requires a one-thread Julia process
  and sets BLAS and FFT providers to one thread
- `ThreadedExecution` uses Julia tasks for independent local members
- `AcceleratedKernelsExecution` uses a reusable task partitioner when the
  AcceleratedKernels weak dependency is loaded
- `DaggerExecution` uses Dagger task graphs and an optional processor scope
  when the Dagger weak dependency is loaded

Use Dagger for large sweeps, multi-process or multi-node work, and locality
constraints. Its callable and simulation state must be serializable when tasks
may leave the current process. The current integration fetches each updated
member at the end of `run_ensemble!`; for distributed use, make the supplied
operation encompass a complete trajectory or sweep unit rather than calling
one remote `step!` at a time. Use AK only after measuring a representative
member size and count on the target many-core host; its scheduling overhead can
outweigh gains on small local workloads.

Ensemble construction rejects shared mutable plant ownership. A custom wrapper
around another runtime should implement
`AdaptiveOpticsSim.ensemble_ownership_roots(wrapper)` and return the mutable
plant/state objects that cannot safely be updated by another member at the
same time. The scheduler operation passed to `run_ensemble!` must mutate and
return normally only after the member is safe for the next coarse operation.

Avoid nested parallelism. When using a threaded ensemble policy, normally keep
BLAS and FFT-provider thread counts at one and let the ensemble own the coarse
parallelism.

## Backend And Allocation Rules

Extension code should follow the package-wide backend rules:

- parametrize by numeric type and array backend where practical
- avoid scalar indexing on GPU arrays
- preallocate workspaces for hot paths
- use explicit `!` methods when mutating state
- centralize RNG ownership in the caller or workspace
- use multiple dispatch or traits instead of `isa` chains

If an extension needs a new reusable behavior, add a small generic seam and one
family implementation first. Then add the second implementation when another
family actually needs it.

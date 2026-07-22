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

## Plant Model Definitions

`ControllableOpticDefinition`, `OpticalPathDefinition`, and
`AcquisitionDefinition` accept only explicitly declared cold model-definition
types. A third-party definition opts in through dispatch:

```julia
AdaptiveOpticsSim.plant_model_definition_style(::Type{MyModelDefinition}) =
    ColdPlantModelDefinition()
```

This method is an ownership assertion, not a recursive mutability test. An
opted-in value must contain configuration only. It must not retain a prepared
plan or workspace, mutable optic, simulation, or detector state, a command
schema or endpoint state, schedule, RNG stream, queue, transport, or HIL
descriptor. Keep those values in separately owned prepared runtime objects.
Types that do not opt in fail closed with `PlantDefinitionError`; do not opt
live controllable optics, detectors, WFSs, atmospheres, or runtime owners into
this trait.

One cold controllable-optic declaration names one physical device and owns one
or more explicit `PlantCommandSchema` values through `command_schemas`. Each
schema targets one `CommandEndpointID` and fixes the exact scalar or
backend-neutral array element type and dimensions, units/sign convention,
basis/revision, absolute or incremental meaning, bounds, and policy vocabulary.
`command_endpoint_ids` is a derived identity view. A named schema key must match
its endpoint identity, and schema/endpoint identities must remain unique in the
plant. These declarations still do not own mutable sequencing, an admission
calendar, applied values, packed controller layout, atomic application,
placement, visibility, or an optical execution group.

Use `validate_plant_command_payload` for non-mutating presentation
compatibility. It deliberately does not clip, admit, sequence, schedule, or
apply a command. Extend neither the schema nor an optic-model definition with
session, external-clock, payload-lease, transport, or HIL descriptor metadata;
that belongs to the later boundary contract.

`prepare_command_endpoint` now binds one exact schema to fixed payload-slot,
accepted-sequence-window, future-calendar, ordinal, and backend capacity. Its
separately owned, qualified `CommandEndpointState` and
`CommandDispositionWorkspace` support warmed `admit_plant_command!`, one
outstanding application-ready claim, and explicit applied/failed/pending-drain
completion without callbacks or run-length storage. Admission copies caller
payloads; array endpoints reserve one additional staging payload so a failed
copy or presentation-time clip cannot corrupt a pending command. Consume and
clear every disposition before reusing its workspace. Give endpoints stable,
unique ordinals when their order keys will be composed.

This is a standalone core endpoint, not yet an optic-model extension hook.
Application-stage state-dependent bounds, effective optic mutation, silence,
safe values, atomic multi-optic latch, and `PreparedPlant` event composition
remain later Gate 4 slices. `prepare_plant` therefore still rejects nonempty
controllable-optic topology rather than silently ignoring it. An incremental
schema must preserve pending deltas; only absolute commands may select
`SupersedeOlderPendingCommands`.

Preparation then dispatches on those same concrete model types. A path method
receives the exact definition, its run-owned frozen source, the plant telescope,
and atmosphere, and returns a `PreparedPathExecutor`:

```julia
function AdaptiveOpticsSim.prepare_path_executor(
    model::MyOpticalModelDefinition,
    definition::AdaptiveOpticsSim.OpticalPathDefinition,
    source::AdaptiveOpticsSim.AbstractSource,
    telescope::AdaptiveOpticsSim.AbstractTelescope,
    atmosphere::AdaptiveOpticsSim.AbstractAtmosphere,
)
    input, result, execution = prepare_my_optics(
        model, source, telescope, atmosphere)
    materialization = AdaptiveOpticsSim.prepare_pupil_opd_materialization(
        atmosphere, telescope, source, input)
    return AdaptiveOpticsSim.PreparedPathExecutor(
        definition, source, telescope, atmosphere, input, result, execution;
        materialization,
        optical_model=my_exact_model_key(model),
        propagation_model=my_exact_propagation_key(model),
        model_revisions=my_revision_key(model),
    )
end
```

The example uses the maintained phase-only path operation, which writes the
current atmosphere OPD into the exact path-local `PupilFunction`. A genuinely
atmosphere-independent model instead passes the qualified
`AdaptiveOpticsSim.AtmosphereIndependentPath()` marker. Do not use that marker
as a fallback for an unsupported atmospheric field or layer-aware model.
Those models provide a concrete materialization owner and extend the qualified
`validate_path_materialization_binding`, `validate_path_materialization`, and
four-argument `materialize_path_input!` dispatches. The validator must check
the exact atmosphere, destination, source, backend, device, shape, and any
model-specific revision without mutating output. The mutating method may then
write only its bound caller-owned path input. This two-phase contract lets a
selection reject every invalid path before materializing the first one.

A `MultiLayerAtmosphere` or `InfiniteMultiLayerAtmosphere` used by
`prepare_plant` declares one stable `AtmosphereLayerID` per layer through its
`layer_ids` keyword. The ordinary atmosphere constructors still permit omitted
IDs for non-plant numerical work, but plant preparation rejects missing or
duplicate stochastic-owner identities. A custom single-owner timed atmosphere
continues to implement its ordinary `initialize_atmosphere!` and
`evolve_atmosphere!` methods against `AbstractRNG`; the prepared plant supplies
the exact owner-bound RNG to those methods.

`input` is a path-local `PupilFunction`, a declared-plane `ElectricField` or
`IntensityMap`, or a concrete tuple of them. The prepared execution determines
which input planes and products it can consume; there is no implicit
resampling or propagation between entry and execution. `result` is an
acquisition-facing photon-rate `IntensityMap` or a concrete
tuple/`OpticalProductBundle` of such maps. The custom telescope must implement
the aperture revision, reflectivity, backend, and physical-device interfaces
consumed by those products. Every leaf in a multi-input or multi-result path
must share that backend and physical device.
The model and propagation keys must be value-comparable and cover every setting
that can alter the result. Do not encode IDs, dimensions, rates, timestamps, or
device ordinals as type parameters. `InstantaneousOpticalSample()` is the
current default and states that the rate product samples one plant instant; it
does not imply a cadence or exposure duration.

`PreparedPathExecutor` snapshots those value descriptions when it builds its
`PathResultKey`. Keep descriptions compact and configuration-only; do not put
live workspaces, device arrays, detector state, or other mutable execution
owners in a key. Key equality and hashing are cold compatibility operations,
not part of warmed optical execution.

A custom `AbstractSource` used in a prepared plant also implements the
qualified `AdaptiveOpticsSim.path_source_geometry_key`,
`AdaptiveOpticsSim.path_source_spectral_key`, and
`AdaptiveOpticsSim.path_source_radiometry_key` methods. Return run-owned,
value-comparable descriptions covering every source property that can change a
path result. If the source contains mutable profile/image storage,
`freeze_source` must copy it before these keys are built.

An acquisition model implements the qualified `prepare_acquisition_provider`
seam. It first calls `require_path_result` for any stricter cold requirements,
then constructs independent detector/WFS state and caller-owned products. The
product metadata is required and run-immutable; include every geometry,
radiometry, unit, layout, or semantic declaration not already carried by the
typed observation or measurement:

```julia
function AdaptiveOpticsSim.prepare_acquisition_provider(
    model::MyAcquisitionDefinition,
    definition::AdaptiveOpticsSim.AcquisitionDefinition,
    path::AdaptiveOpticsSim.PreparedPathExecutor,
)
    AdaptiveOpticsSim.require_path_result(
        path; optical_model=model.required_optical_model)
    execution, observation, measurement = prepare_my_acquisition(model, path)
    products = AdaptiveOpticsSim.AcquisitionProducts(
        observation, measurement;
        metadata=my_product_metadata(model, path, observation, measurement))
    return prepare_full_optical_provider(execution, products)
end
```

Core constructs the `PreparedAcquisitionOwner` after validating the returned
provider against the exact path result. An acquisition extension does not
construct or return the owner itself.

Use qualified `AdaptiveOpticsSim.WFSOpticalPathExecution` to adapt an existing
Gate 0 WFS optical plan, `AdaptiveOpticsSim.FrameAcquisitionExecution` for a
frame detector plus a distinct caller-owned observation, and
`AdaptiveOpticsSim.WFSAcquisitionExecution` to compose already prepared WFS
acquisition and estimator plans. A deterministic concrete path execution type
must extend both the three-argument `execute_path!` dispatch and the qualified
`AdaptiveOpticsSim.validate_path_execution_binding(execution, input, result)`
seam. A stateful stochastic path can instead extend the four-argument form that
receives its prepared provider `AbstractRNG`. If it needs additional
independent device streams, extend qualified
`additional_path_rng_owner_roles` and `execute_path_rngs!`, then obtain each
declared stream directly with `rng_stream_state(group, Val(:role))`. A
different full-optical acquisition execution type similarly extends the
four-argument `execute_acquisition!` dispatch and
`AdaptiveOpticsSim.validate_acquisition_execution_binding(execution,
path_result, products)`. Extra acquisition/device streams use qualified
`additional_acquisition_rng_owner_roles` and `execute_acquisition_rngs!`.
The full-optical provider wrapper delegates these seams to the acquisition
execution object.
Role tuples and their `Val` lookups are prepared once; models never consult a
global RNG registry in the hot path. Each validator must reject mismatched exact
storage or state before mutation. Do not store a `Function`, abstract executor
vector, schedule, RNG registry, queue, or transport in model execution owners.
Preparation may allocate; warmed execution must retain the allocation contract
of its underlying stages.

A custom reduced-order provider instead returns qualified
`AdaptiveOpticsSim.PreparedAcquisitionProvider(implementation, products)` and
implements four qualified methods:

- `acquisition_provider_style(::Type{MyProvider})` returns
  `CommandResponsiveReducedOrderProviderStyle()`
- `acquisition_provider_payload_work(::Type{MyProvider})` returns a nonempty
  symbol describing its principal payload work
- `validate_acquisition_provider_binding(implementation, path_result,
  products)` rejects provider-specific incompatible bindings; core validates
  the products against its private prepared contract first
- `execute_acquisition_provider!(products, path_result, implementation,
  rngs)` mutates and returns the exact `products` value

Keep immutable reduced-order parameters separate from a mutable single-writer
state object. The provider must remain causally responsive to effective optic
commands; replaying completed products is not a reduced-order AO plant.
`acquisition_product_contract` and
`validate_acquisition_product_contract` cover standard arrays,
`IntensityMap`, `WFSObservation`, `WFSMeasurement`, `Ref`, `nothing`, and
concrete tuples. Custom wrappers extend the qualified leaf contract,
validation, and `copy_acquisition_product!` methods.

Core supplies three nonresponsive synthetic implementations:
`prepare_unchanged_synthetic_provider` republishes unchanged destination
contents, `prepare_copied_synthetic_provider` copies one preparation-owned
snapshot, and `prepare_cyclic_replay_provider` cycles through a fixed-size
preparation-owned corpus. These methods validate the full logical contract
before returning. Their payload-work declarations describe execution only;
they do not establish RTC latency, capacity, cache-residency, or optical
evidence.

Callers prepare a fixed acquisition subset with
`prepare_acquisition_selection(plant, ids)`. Its required full-optical paths
and acquisition owners are available through `prepared_paths` and
`prepared_acquisitions` accessors. Repeated execution supplies either one
explicit current `AtmosphereEpoch` to `execute_acquisition_selection!` or one
absolute atmosphere model time to `execute_acquisition_selection_at!`.
`prepare_plant(definition; run_seed, rng_derivation_version)` owns all stateful
streams, so neither selected-execution method accepts an RNG argument.
`rng_replay_metadata(plant)` provides structured replay identity and seed data
without granting another writer access to those streams.
Reduced-order and synthetic/replay selections do not form an otherwise unused
full-optical path. Provider style is fixed by preparation and cannot change in
a prepared plant.

## Calibration Illumination Evaluators

Calibration is a scenario role; do not add a calibration source superclass,
global mode, instrument selector, or propagation-bypass flag. Bind a prepared
evaluator to the ordinary path input with `prepare_illumination_entry`, then
pass that entry as the `materialization` of `PreparedPathExecutor`. The path's
normal execution and acquisition provider remain unchanged.

The containing optical model's `optical_model` and revision key must still
cover the entry boundary and every evaluator parameter that can change the
path result. Do not put evolving evaluator state, model time, RNG state, or
backend arrays in that key.

The maintained entry tags accept these products:

| Entry tag | Accepted caller-owned payload |
|---|---|
| `PupilFunctionIlluminationEntry()` | `PupilFunction` with fully declared metadata |
| `ElectricFieldIlluminationEntry()` | `ElectricField` on its declared optical plane |
| `IntensityMapIlluminationEntry()` | `IntensityMap` on its declared optical plane |
| `ExternalOpticsResultIlluminationEntry()` | `ElectricField` or `IntensityMap` already formed by a prepared external executor |
| `DetectorInputIlluminationEntry()` | focal- or detector-plane `IntensityMap` with declared spectrum, photon-rate or dimensionless normalization, spatial-density or cell-integrated measure, and incoherent semantics |

An immutable user definition prepares an immutable evaluator wrapper. Mutable
single-writer state and backend-native workspace may be referenced by that
wrapper, but remain separate objects. Implement these qualified methods:

```julia
struct MyIlluminationDefinition{T}
    level::T
end

mutable struct MyIlluminationState{T}
    evaluations::Int
    last_time_s::T
end

struct PreparedMyIllumination{P,S,B,D}
    params::P
    state::S
    backend::B
    device::D
end

AdaptiveOpticsSim.illumination_combination(
    ::Type{<:PreparedMyIllumination}) = SingleIllumination()

function AdaptiveOpticsSim.prepare_illumination_evaluator(
    definition::MyIlluminationDefinition,
    destination::IntensityMap,
    ::DetectorInputIlluminationEntry,
)
    state = MyIlluminationState(0, zero(eltype(destination.values)))
    return PreparedMyIllumination(definition, state,
        backend(destination), plane_device(destination.values))
end

function AdaptiveOpticsSim.validate_illumination_evaluator_binding(
    evaluator::PreparedMyIllumination,
    destination::IntensityMap,
    ::DetectorInputIlluminationEntry,
)
    typeof(backend(destination)) === typeof(evaluator.backend) ||
        throw(PlantPreparationError(:illumination, :backend,
            "prepared illumination backend changed"))
    plane_device(destination.values) == evaluator.device ||
        throw(PlantPreparationError(:illumination, :device,
            "prepared illumination device changed"))
    return nothing
end

function AdaptiveOpticsSim.evaluate_illumination!(
    destination::IntensityMap,
    evaluator::PreparedMyIllumination,
    model_time,
    rng::AbstractRNG,
)
    evaluator.state.evaluations += 1
    evaluator.state.last_time_s = model_time
    fill!(destination.values, evaluator.params.level)
    return destination
end
```

Core validates the destination against its private prepared payload contract
before invoking the evaluator-specific binding validator; extension methods
validate only their additional parameter, state, workspace, backend, device,
and boundary bindings. The evaluator receives explicit plant time and its
path-owned `:illumination` RNG stream and must mutate and return the exact
destination without allocating after warmup. It declares one of
`SingleIllumination()`,
`ExclusiveIlluminationSelection()`, `CoherentFieldCombination()`, or
`IncoherentIntensityAddition()` through dispatch; core never guesses how
contributions combine. The `visibility` value supplied to
`prepare_illumination_entry` is a required, defensively snapshotted application
description. Keep it configuration-only rather than storing live state,
workspace, a closure, or transport ownership. Core records but does not
interpret its topology.

The warmed wrapper revalidates product identity, metadata, shape, numeric type,
backend, device, and combination semantics. It deliberately does not scan a
device array for finite or nonnegative samples on every evaluation. A custom
evaluator therefore owns that numerical invariant; downstream prepared stages
may rely on it. The native uniform model enforces a finite nonnegative level at
construction.

The native `UniformIntensityIllumination(value; combination)` is useful for a
spatially uniform intensity entry. Its value is expressed in the destination
map's already declared normalization and spatial measure; it does not perform a
radiometric conversion or invent a spectrum. Detector-input entries still pass
through ordinary detector response, quantum efficiency, exposure, readout, and
product-provider semantics.

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
`OpticalProductBundle`. Bundle membership is fixed at construction; the arrays
owned by its product leaves remain mutable destinations for prepared optical
formation.

Extensions should call the qualified validation seams
`AdaptiveOpticsSim.validate_wfs_optical_input`,
`AdaptiveOpticsSim.validate_wfs_optical_products`,
`AdaptiveOpticsSim.validate_wfs_observation` or
`AdaptiveOpticsSim.validate_wfs_observations`, and
`AdaptiveOpticsSim.validate_wfs_measurement` as applicable before returning a
prepared plan. Every new prepared plan also implements its corresponding exact
binding validator:

- `AdaptiveOpticsSim.validate_wfs_optical_formation_binding`
- `AdaptiveOpticsSim.validate_wfs_acquisition_binding`
- `AdaptiveOpticsSim.validate_wfs_estimation_binding`

The containing plant owner calls these validators during construction, and the
stage calls them again before mutation. They remain qualified extension APIs
rather than ordinary exported workflow names.

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

Concrete WFS family types may provide a composed `measure!` workflow over the
same prepared stages. That workflow is an ordinary high-level API, not a second
ownership model: optical inputs remain caller-owned, stage workspaces remain
sensor-owned and preallocated, and the implementation dispatches through the
prepared seams rather than branching on sensor families.

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
plant. The dense `ModalReconstructor` remains the full-rank accuracy baseline.

Runtime profiles and execution plans are orthogonal. A profile such as
`HILRuntimeProfile` selects modeled delays and output fidelity. An execution
plan selects storage and synchronization behavior. The HIL profile uses the
minimum valid direct-science zero-padding factor of `1` (no focal-grid
oversampling) when science pixels are requested; omitting science pixels still
prepares no science stage. The execution-plan choices are:

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

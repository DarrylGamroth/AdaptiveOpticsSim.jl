# Extension Guide

Status: active

This guide describes the maintained extension seams for adding new detector
families, wavefront sensors, deformable mirrors, controllable optics,
controllers, and reconstructors.

Prefer small concrete types plus multiple dispatch. Do not add central
type-switching registries for new families.

Plant extension seams are owned by `AdaptiveOpticsSim.Plant`. Extend the
qualified function in that module; do not define a root-package forwarding
method. Routine plant constructors may be imported with
`using AdaptiveOpticsSim.Plant`, while the advanced seams below remain
qualified deliberately.

## Canonical Generic Ownership

Every extension hook has one canonical module owner. The exact method inventory
and owner mapping is frozen in
[`../test/contracts/namespace_authority.toml`](../test/contracts/namespace_authority.toml).
The owner categories are:

- `Backends` for array allocation, random generation, FFT planning/execution,
  reduction, compute-device identity, and accelerator registration
- `Atmospheres` for phase-PSD and phase-noise generation
- `Detectors` for detector execution plans, frame noise, counting output, and
  device readout correction
- `WavefrontSensors` for grouped/geometric sensing and LiFT solver fallbacks
- `Calibration` for inverse-operator construction
- `Tomography` for tomography-specific Hermitian solves
- `Ensembles` for optional coarse scheduler execution
- the root package only for cross-domain config and telemetry serialization

Implemented domain owners, including `Backends`, `Optics`, `Tomography`, and
`Ensembles`, receive extension methods at `AdaptiveOpticsSim.Owner.name`; no
root-package forwarding generic remains part of the supported API. In particular,
tomography-specific Hermitian solvers extend
`AdaptiveOpticsSim.Tomography.stable_hermitian_right_division`. Scheduler
extensions likewise extend
`AdaptiveOpticsSim.Ensembles.init_ensemble_scheduler` and
`AdaptiveOpticsSim.Ensembles.execute_ensemble!`. The maintained stage is
recorded in
[`../test/contracts/namespace_migration_state.toml`](../test/contracts/namespace_migration_state.toml).

## Prepared Numerical Execution Interfaces

Prepared algorithm interfaces follow an `AbstractFFTs`-style nominal model,
but each protocol remains with its canonical scientific owner. Add an abstract
plan type only when multiple implementations share documented operations and
invariants. Do not add a universal package plan root, an `Interfaces` module,
or a generic `process!` operation merely to make unrelated algorithms look
uniform.

Use the following ownership split:

| Value | Extension contract |
|---|---|
| Definition | Cold reusable configuration with no prepared storage or live run ownership. |
| Plan | Run-immutable dimensions, mappings, coefficients, compatibility rules, and logically immutable reentrant operators. |
| State | Persistent mutable values that can affect a later scientific result; one writer. |
| Workspace | Replaceable scratch and execution resources; recreation cannot change the deterministic trajectory. |
| Product | Caller-visible output with explicit metadata and ownership. |
| Prepared execution owner | Exact validated binding of the concrete plan, state, workspace, products, and backend/device/context. |

An abstract plan subtype signals nominal participation. It is not sufficient
by itself. The owning module must document and test:

- required domain-specific mutating methods
- input, output, state, and workspace invariants
- supported numeric types, array backends, and exact compute devices
- structured preparation and execution failures
- stale or foreign binding rejection before mutation
- type inference and warmed CPU allocation behavior
- applicable accelerator residency and synchronization behavior

Do not store `Any` in an implementation field, array element type, reference,
or other prepared owner. Preparation must normalize user or extension input to
a concrete internal representation before returning. An unconstrained
`::Any`, `<:Any`, or `Vararg{Any}` in a fallback or matching method signature
does not itself impose boxed storage; review such signatures for dispatch
intent rather than rejecting them with a textual search.

Keep concrete plan, state, and workspace types as type parameters of a prepared
owner. Do not store the abstract plan root in a hot field or collection. Exact
input/output array identity belongs to the prepared owner unless the underlying
operator itself requires that identity. A reentrant logically immutable FFT or
backend handle may be part of a plan; a scratch-owning, stream-bound,
task-bound, or non-reentrant handle belongs to a workspace or prepared owner.

`Optics.AbstractPropagationPlan` is the nominal plan interface for native
Fraunhofer and Fresnel propagation. Implementations provide
`Optics.propagation_input_metadata` and
`Optics.propagation_output_metadata`; domain-specific `propagate_field!`
methods accept a concrete plan and matching concrete workspace. A prepared
propagation model retains that plan/workspace pair. Native direct imaging and
spatial filtering use a stricter exact-owner boundary:
`Optics.PreparedDirectImaging` binds its input, work field, output, plan, and
workspace, while `Optics.PreparedSpatialFilter` additionally binds its physical
filter. Foreign identities, stale metadata, incompatible backend/device
storage, dimensions, numeric types, or forbidden workspace aliases must raise
a structured error before product mutation. Plans are reentrant when their
logically immutable coefficient arrays are not mutated; workspaces and prepared
owners are single-writer. Recreating a compatible workspace must not change a
deterministic optical result. Warmed CPU execution remains allocation-free;
accelerator launch-allocation claims require backend-specific evidence.

Frame-detector acquisition uses the concrete contract owned by
`AdaptiveOpticsSim.Detectors`. `DetectorAcquisitionPlan` retains only
run-immutable detector parameters, optical metadata, dimensions, radiometric
scaling, and channel QE. `PreparedDetectorAcquisition` binds the exact detector
and mutable intensity storage to that plan, persistent `DetectorState`,
replaceable `DetectorWorkspace`, and caller-visible `DetectorProducts`.
Qualified accessors named `detector_acquisition_detector`,
`detector_acquisition_input`, `detector_acquisition_plan`,
`detector_acquisition_state`, `detector_acquisition_workspace`, and
`detector_acquisition_products` expose those ownership roles without exporting
them from the package root. Routine code calls `prepare_detector_acquisition`
once and executes `capture!(prepared, rng)` or the incremental prepared form;
there is no compatibility call that separately passes the detector, input, and
plan.

A frame-detector extension must keep deterministic persistent values out of
workspace, keep readout fitting scratch out of products, grant one exclusive
lifecycle-state owner per prepared detector lifecycle, and reject replacement
or forbidden aliasing of exact bound owners before mutation. Its prepared CPU
path must remain inferred and allocation-free after warmup. Accelerator support
requires same-device input, state, workspace, and products, scalar indexing
disabled, explicit synchronization at the documented completion boundary, and
numerical parity evidence. Kernel or backend launch allocation is measured
separately and is not covered by the CPU zero-allocation contract.

New and migrated package or extension implementation does not directly use
Julia's `Memory` type. Use
[FixedSizeArrays.jl](https://github.com/JuliaArrays/FixedSizeArrays.jl)'s
`FixedSizeArray{ConcreteT}` for homogeneous armed host storage; construction
must reject a non-concrete element type. A `Vector` may be used as a cold
builder, but preparation must seal it into fixed-size storage before execution.
Fixed-size storage prevents cardinality changes, not element replacement;
exact prepared owners must retain an independent ordered membership snapshot
and reject replacement before numerical mutation.
Use a concrete tuple or union for small bounded heterogeneous composition. For
larger topologies, group owners into concrete homogeneous fixed-size arrays by
execution family and select those groups with compact concrete descriptors or
a purpose-built owner whose type does not encode path count. Use the concrete
backend array type for numerical storage.

On Julia versions where FixedSizeArrays.jl uses `Memory` as its internal
backing, that backing remains an implementation detail of the dependency. Do
not name, dispatch on, assert, or expose the backing type in package code or
tests. The fixed-size semantic container does not repair element-type erasure:
changing `Memory{Any}` to `FixedSizeArray{Any}`, `Vector{Any}`, an abstract
vector, or a vector of an uninstantiated parametric family does not satisfy the
execution contract. Use a concrete owner or an explicit prepared function
barrier with inference and allocation evidence.

The implementation migration and its required evidence are tracked in
[issue #225](https://github.com/DarrylGamroth/AdaptiveOpticsSim.jl/issues/225).

## Algorithm Graph Adapters

Optional packages implement the qualified-public
[`Algorithm graph adapter protocol`](./glossary.md) owned by
`AdaptiveOpticsSim.AlgorithmGraphs`. Calculon remains the authority for
scientific declarations, configuration, ports, formats, properties, sparse
parameters, preparation, processing, and reset. The graph adapter translates
that interface into exact graph storage; it must not define a parallel
scientific plan or `process!` API.

The protocol methods are:

| Method | Adapter obligation |
|---|---|
| `prepare_algorithm_instance(Declaration, configuration)` | Prepare a concrete prototype before final graph buffers exist. Do not retain provisional input or output storage. |
| `algorithm_port_contracts(Declaration, prototype)` | Return one fixed tuple of `graph_port_contract(...)` values in declaration order. Names, direction, role, element type, shape, schema, and layout must match the Calculon declaration. |
| `bind_algorithm_instance(Declaration, prototype, inputs, outputs)` | Optionally return a specialized prepared owner bound to the final admitted named buffers. The default returns the prototype unchanged. |
| `replace_algorithm_parameter!(prototype, Val(name), values)` | Apply an admitted startup sparse parameter during preparation. Live parameter replacement is not yet a graph capability. |
| `process_algorithm!(prepared, outputs, inputs)` | Execute bounded repeated-path work against the exact named buffers without replacing their storage. |
| `reset_algorithm!(prepared)` | Reset persistent algorithm state at a serialized graph boundary. |

Use these exact names through the `AlgorithmGraphs` module. They are public but
not exported so extension code makes ownership visible. There are no
underscore-prefixed aliases. An adapter whose underlying Calculon instance can
operate on graph-supplied buffers directly may return that instance unchanged.
An adapter that must prebind FFT plans, wavefronts, or other exact owners uses
`bind_algorithm_instance` only after all formats, topology, startup parameters,
and storage have been admitted.

Preparation may allocate and may fail. `process_algorithm!` and
`reset_algorithm!` execute inside the graph's retained exact-device context;
the adapter must not select a different device or stream. A CPU repeated path
requires inference and an explicit warmed allocation budget. Accelerator claims
require native storage, exact-device residency, scalar-index prohibition, and
real-hardware evidence. Input/output aliasing, replacement, wrong schema or
shape, and foreign or wrapped storage must fail before repeated execution.

## Exact Compute-Device Selection

Gate 9A preparation distinguishes a semantic backend family from one exact
runtime device. `Backends.compute_device_availability(device)` is the cold
availability query. It returns `ComputeDeviceAvailable` or
`ComputeDeviceUnavailable`; inspect the latter through
`compute_device_unavailable_reason`. Availability is not model support or a
memory/deadline admission result.

`Backends.allocate_device_array(device, T, dims...)` is the qualified
preparation-time allocation seam. An accelerator extension that supports exact
selection must extend `Backends.compute_device_availability`,
`Backends._with_compute_device`, and
`Backends._prepare_device_execution_context` for its family-qualified
`AcceleratorComputeDevice`. Device identifiers use the owning runtime's native
same-process convention: CUDA ordinals are zero-based and AMDGPU identifiers
are one-based in the maintained extensions. The extension must:

- validate the identifier and return a structured unavailable reason;
- select exactly that runtime device and restore the caller's previous device
  context even when allocation or context construction throws;
- return storage whose `Backends.compute_device` equals the request; and
- create prepared stream/context state on that same device.

Do not implement exact selection by allocating on the current default device
and inferring its identity afterward. `ComputeDeviceError` reports structured
`operation`, `reason`, and `device` fields when an exact selection or allocation
cannot be honored. Exact-target preparation is cold; repeated kernels retain
their prepared context and do not call this selection path.

Whole-plant preparation uses the same contract through the mandatory
`prepare_plant(definition, target; ...)` positional target. It materializes the
cold telescope and timed-atmosphere definitions, prepares all numerical model
owners, and derives array-command storage inside one retained target context.
Every target-authoritative data-plane array owned by the resulting prepared
execution graph must report exactly `target` through
`Backends.compute_device`. The retained cold definition and its caller-owned
sample inputs, scalar command values, host configuration and registries,
documented host inspection/staging mirrors, and deliberate host RNG staging
are not accelerator data-plane arrays. `Backends.compute_device(plant)` reports
the prepared target. A custom model must allocate its owned execution arrays
through the target selected by the supplied prepared telescope, atmosphere,
inputs, and destinations. Do not copy directly from one accelerator to a
different accelerator during preparation; accept caller-owned host input or
require the caller to stage the transfer explicitly.

Prepared-owner extension points validate exact residency through the qualified
`Optics.validate_telescope_target`,
`Atmospheres.validate_timed_atmosphere_target`,
`Plant.validate_controllable_optic_target`,
`validate_controllable_optic_state_target`,
`validate_controllable_optic_workspace_target`,
`validate_pupil_surface_coupling_target`,
`validate_autonomous_optic_coupling_target`,
`validate_path_execution_target`,
`validate_path_materialization_target`,
`validate_acquisition_execution_target`,
`validate_acquisition_provider_target`, and
`validate_illumination_evaluator_target` seams. The controllable-optic
definition, state, and workspace checks are separate because event-loop state
and workspace are allocated after immutable optic preparation. Each default
fails closed for an unknown prepared implementation. Extend every seam owned by
the preparation object rather than weakening whole-plant traversal or
inspecting fields by reflection. A custom prepared WFS plan extends the
WFS-owned qualified `WavefrontSensors.validate_wfs_target` seam. Lower-level
built-in validators remain implementation details of their owning subsystems.

The maintained structured exact-device reason vocabulary is:

- `:exact_device_selection_unavailable` when no owning extension implements
  exact selection;
- `:invalid_device_identifier` when the identifier does not follow the
  runtime's identifier convention;
- `:backend_runtime_unavailable` when the owning runtime is not functional;
- `:device_unavailable` when the runtime cannot address that otherwise valid
  identifier; and
- `:wrong_device` when allocation violates the exact-residency invariant.

`ComputeDeviceError.operation` is `:select`, `:prepare_context`, or `:allocate`
according to the failed preparation boundary. These values are programmatic
contract data; user-facing detail remains in the exception message.

## Source Layout

Use lower-case directory names for new source-tree locations. Julia type names
may use CamelCase, but package directories should remain lower-case subsystem
names such as `src/wfs`, `src/detectors`, `src/optics`, `src/control`, or
`src/ensembles`.

## Plant Model Definitions

`PlantDefinition` accepts one immutable `Optics.AbstractTelescopeDefinition`
and one immutable `Atmospheres.AbstractTimedAtmosphereDefinition`, not a
prepared numerical telescope or atmosphere. The built-in
`Optics.TelescopeDefinition` and timed-atmosphere definitions are reusable
configuration. Where retained caller inputs satisfy the target preparation
contract, preparing the same plant definition for another target creates
independent numerical owners and evolution state. Inspect a declaration with
`Plant.telescope_definition` and `Plant.atmosphere_definition`; inspect a
prepared owner with `Plant.plant_definition`, `Plant.prepared_telescope`, and
`Plant.prepared_atmosphere`.

A custom `Optics.AbstractTelescopeDefinition` implements the qualified-public
`Optics.prepare_telescope(definition, target)` method and returns a concrete
`Optics.AbstractTelescope`. Its prepared owner implements
`Optics.validate_telescope_target(telescope, target)` and the aperture
revision, reflectivity, backend, and compute-device interfaces consumed by its
paths. A custom `Atmospheres.AbstractTimedAtmosphereDefinition` similarly
implements `Atmospheres.prepare_timed_atmosphere(definition, telescope,
target)` and its prepared `Atmospheres.AbstractTimedAtmosphere` implements
`Atmospheres.validate_timed_atmosphere_target(atmosphere, target)`. Both
validation seams must inspect every owner-controlled data-plane array and fail
closed on the wrong target without mutating state.
Timed-atmosphere owners retain a unique `Atmospheres.AtmosphereIdentity` and an
`Atmospheres.AtmosphereTimelineState` made by
`Atmospheres.new_atmosphere_timeline(T)`, and implement the qualified-public
`atmosphere_identity` and `atmosphere_timeline` accessors. These values let the
shared epoch and owner-bound RNG machinery remain independent of a concrete
atmosphere layout.

`ControllableOpticDefinition`, `OpticalPathDefinition`, and
`AcquisitionDefinition` accept only explicitly declared cold model-definition
types. A third-party definition opts in through dispatch:

```julia
AdaptiveOpticsSim.Plant.plant_model_definition_style(::Type{MyModelDefinition}) =
    AdaptiveOpticsSim.Plant.ColdPlantModelDefinition()
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
or an optical execution group. Each controllable-optic definition does require
one explicit optical `placement` and one `visibility` declaration. Use
`PupilPlanePlacement`, `AtmosphericConjugatePlacement`, or
`FocalPlanePlacement` and `AllPathVisibility` or `SelectedPathVisibility`;
do not encode these contracts in endpoint identity or model labels.

Use `validate_plant_command_payload` for non-mutating presentation
compatibility. It deliberately does not clip, admit, sequence, schedule, or
apply a command. Extend neither the schema nor an optic-model definition with
session, external-clock, payload-lease, transport, or HIL descriptor metadata;
that belongs to the later boundary contract.

`prepare_command_endpoint` now binds one exact schema to fixed payload-slot,
accepted-sequence-window, future-calendar, ordinal, and payload-storage
capacity. Scalar payload slots require `Backends.CPUBackend()` and remain
host-resident; fixed-shape array slots use the selected array backend. Its
separately owned, qualified public `CommandEndpointState` and
`CommandDispositionWorkspace` support warmed `admit_plant_command!`, one
outstanding application-ready claim, and explicit applied/failed/pending-drain
completion without callbacks or run-length storage. Admission copies caller
payloads; array endpoints reserve one additional staging payload so a failed
copy or presentation-time clip cannot corrupt a pending command. Consume and
clear every disposition before reusing its workspace. Give endpoints stable,
unique ordinals when their order keys will be composed; equal scheduled times
order by endpoint ordinal before the endpoint-local sequence.

The separately owned qualified public `CommandApplicationState` binds one explicit
initial effective command and any required copied safe command to one exact
endpoint-state owner before its first successful admission.
`apply_claimed_plant_command!` transactionally implements absolute
replacement or incremental addition, finite-result checking, declared
application-stage bound behavior, and terminal disposition. Array staging and
effective storage remain on the endpoint backend. Its silence operations
schedule exact safe/fail transitions from admission or application age, with a
due equal-time command resolved first.

`CommandEndpointConfiguration` supplies each declared endpoint's run-specific
capacities, initial effective value, and optional safe value. It does not
select storage placement: array-command storage follows the exact plant target
and scalar command values remain host-resident.
Pass exactly one configuration per declared endpoint through
`prepare_plant(definition, target; command_endpoints=configurations, ...)`.
Plant preparation
canonicalizes endpoint and optic order by stable identity, prepares every
declared device, and rejects missing or extra configurations.

A controllable-optic extension separates immutable preparation, mutable
single-writer physical state, and scratch:

```julia
const Plant = AdaptiveOpticsSim.Plant

function Plant.prepare_controllable_optic(
    model::MyOpticDefinition,
    definition::Plant.ControllableOpticDefinition,
    telescope::AdaptiveOpticsSim.Optics.AbstractTelescope,
    atmosphere::AdaptiveOpticsSim.Atmospheres.AbstractAtmosphere,
)
    return prepare_my_immutable_optic_plan(
        model, definition, telescope, atmosphere)
end

function Plant.prepare_controllable_optic_state(
    plan::MyPreparedOptic,
    definition::Plant.ControllableOpticDefinition,
    endpoint_ids::Tuple,
    initial_commands::Tuple,
)
    return MyOpticState(plan, endpoint_ids, initial_commands)
end

Plant.prepare_controllable_optic_workspace(plan::MyPreparedOptic) =
    MyOpticWorkspace(plan)

function Plant.stage_controllable_optic_command!(
    plan::MyPreparedOptic,
    state::MyOpticState,
    workspace::MyOpticWorkspace,
    endpoint::Plant.CommandEndpointID,
    effective_command,
    timestamp::Plant.PlantTimestamp,
)
    stage_complete_physical_response!(
        workspace, plan, state, endpoint, effective_command, timestamp)
    return nothing
end

function Plant.commit_controllable_optic_command!(
    plan::MyPreparedOptic,
    state::MyOpticState,
    workspace::MyOpticWorkspace,
    endpoint::Plant.CommandEndpointID,
    timestamp::Plant.PlantTimestamp,
)
    publish_staged_physical_response!(state, workspace, endpoint, timestamp)
    return nothing
end

function Plant.prepare_controllable_optic_path_coupling(
    plan::MyPreparedOptic,
    definition::Plant.ControllableOpticDefinition,
    path::Plant.PreparedPathExecutor,
)
    return Plant.prepare_sampled_pupil_footprint_coupling(
        plan.surface_metadata,
        plan.surface_prototype,
        path,
        Plant.controllable_optic_placement(definition);
        registration=plan.pupil_relay_registration,
    )
end

function Plant.apply_controllable_optic_surface!(
    input::AdaptiveOpticsSim.Optics.PupilFunction,
    plan::MyPreparedOptic,
    state::MyOpticState,
    coupling::Union{
        Plant.PreparedIdentityPupilFootprintCoupling,
        Plant.PreparedPupilFootprintCoupling,
    },
)
    return Plant.apply_sampled_pupil_surface!(
        input, state.visible_surface, coupling)
end
```

`prepare_controllable_optic` must return immutable data; keep all evolving
response state in the separately constructed state. Every array in
`initial_commands` is a fresh state-owned copy, so the state constructor may
retain and mutate it; it must not substitute caller or prepared-plan storage.
Staging may reject or fail before publication. The supplied timestamp is the
exact plant instant at which the effective physical response becomes visible;
use it to preserve continuity in a time-dependent model. After staging
succeeds, commit must be bounded and nonthrowing so an explicit
`PlantCommandTransaction` cannot expose only part of a multi-optic update. Do
not infer transaction membership from placement, equal timestamps, or packed
command storage. Every endpoint participating in a transaction must use
`PreservePendingCommands`. Outside transactions, an incremental schema must
preserve pending deltas; only absolute commands may select
`SupersedeOlderPendingCommands`.

The default `controllable_optic_execution_role` is
`PupilSurfaceExecutionRole`. Preparation resolves the declared path visibility
into bounded per-path ranges and co-placed groups. The event loop applies only
the visible members of a due materialized path, in canonical optic-identity
order within their compatible coupling group. Extend only path-input types for
which the device has a well-defined surface operation.

The default `prepare_controllable_optic_path_coupling` supports a model-owned
same-grid operation at `PupilPlanePlacement`. A model exposing a sampled OPD
surface uses `prepare_sampled_pupil_footprint_coupling` to bind that surface's
metric `PupilPlane` metadata to the exact prepared path input. Preparation
combines the path source direction, finite-height LGS cone scale, conjugate
altitude, sampled-grid metadata, and an optional `PupilRelayRegistration`.
Repeated application then uses the same-grid fast path or finite-support
bilinear interpolation without host scalar indexing on accelerator arrays.
Expanded multi-direction sources must be represented by separate prepared
paths at an atmospheric conjugate; a `SpectralSource` is accepted because its
wavelength samples share one direction.

`PupilRelayRegistration` is relay/path geometry: magnification, rotation,
parity, and metric decenter from a geometric source footprint to the sampled
surface grid. It is not a DM actuator `Misregistration`. Apply actuator
misregistration while forming the sampled DM surface, then apply the relay
registration exactly once while coupling that completed surface to a path.

A locally generated path-specific waveform instead returns
`AutonomousPathExecutionRole` and provides the qualified autonomous-device
seams:

```julia
struct MyPeriodicOpticFidelity <:
    Plant.AbstractAutonomousOpticFidelity end

Plant.controllable_optic_execution_role(::MyPreparedPeriodicOptic) =
    Plant.AutonomousPathExecutionRole()

Plant.prepare_autonomous_periodic_optic(
    plan::MyPreparedPeriodicOptic,
    path::Plant.PreparedPathExecutor,
    fidelity::MyPeriodicOpticFidelity,
) = prepare_my_path_coupling(plan, path, fidelity)

function Plant.validate_autonomous_optic_coupling_target(
    coupling::MyPreparedPathCoupling,
    target::Backends.AbstractComputeDevice,
)
    validate_my_path_coupling_target(coupling, target)
    return coupling
end

function Plant.initialize_autonomous_periodic_optic!(
    plan::MyPreparedPeriodicOptic,
    state::MyPeriodicOpticState,
    reference::Plant.AbstractWaveformPhaseReference,
)
    initialize_my_phase_reference!(state, reference)
    return nothing
end

function Plant.reset_autonomous_periodic_optic_phase!(
    plan::MyPreparedPeriodicOptic,
    state::MyPeriodicOpticState,
    timestamp::Plant.PlantTimestamp,
    sequence::UInt64,
)
    reset_my_phase_reference!(state, timestamp, sequence)
    return nothing
end

function Plant.evaluate_autonomous_periodic_optic!(
    plan::MyPreparedPeriodicOptic,
    state::MyPeriodicOpticState,
    coupling::MyPreparedPathCoupling,
    timestamp::Plant.PlantTimestamp,
)
    evaluate_my_periodic_optic!(plan, state, coupling, timestamp)
    return nothing
end
```

The cold model still owns ordinary bounded `PlantCommandSchema` setpoints, and
an `AutonomousPeriodicOpticDefinition` must bind the prepared device to one
scheduled full-optical path and one immutable fidelity/reference contract.
Such a device uses `FocalPlanePlacement()` and a
`SelectedPathVisibility` containing exactly that coupled path.
The current event composer accepts the declared `FreeRunningPhaseReference`,
`TriggerSourcePhaseReference`, and `TriggerResetPhaseReference` relationships.
A new phase-reference family also needs explicit event-composition ownership
and notification semantics; subtype construction alone does not register one.
Preparation must reject incompatible paths and conflicting exclusive
couplings. Evaluation mutates only storage owned by that prepared path, stays
bounded, and allocates no steady-state heap storage when it is on the optical
hot path. Do not use the autonomous role merely to encode another placement
kind, and do not hide detector exposure integration in a waveform
quadrature.

Preparation then dispatches on those same concrete model types. A path method
receives the exact definition, its run-owned frozen source, the plant telescope
and atmosphere, and the plant's prepared device execution context, and returns
a `PreparedPathExecutor`:

```julia
function AdaptiveOpticsSim.Plant.prepare_path_executor(
    model::MyOpticalModelDefinition,
    definition::AdaptiveOpticsSim.Plant.OpticalPathDefinition,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    telescope::AdaptiveOpticsSim.Optics.AbstractTelescope,
    atmosphere::AdaptiveOpticsSim.Atmospheres.AbstractAtmosphere,
    context,
)
    input, result, execution = prepare_my_optics(
        model, source, telescope, atmosphere)
    materialization = AdaptiveOpticsSim.Plant.prepare_pupil_opd_materialization(
        atmosphere, telescope, source, input)
    return AdaptiveOpticsSim.Plant.PreparedPathExecutor(
        definition, source, telescope, atmosphere, input, result, execution;
        context=context,
        materialization,
        optical_model=my_exact_model_key(model),
        propagation_model=my_exact_propagation_key(model),
        model_revisions=my_revision_key(model),
    )
end
```

Pass the supplied `context` unchanged. It selects the exact plant compute
device and, for an accelerator, retains the prepared backend stream; a model
extension must not select or prepare an independent stream. The constructed
`PreparedPathExecutor` retains that context, and core gives the same context to
the corresponding `PreparedAcquisitionOwner`. Public
`materialize_path_input!`, `execute_path!`, and `execute_acquisition!` calls
enter the retained context and restore the caller's previous device and stream
selections after success or an exception. A successful return is also a
backend completion boundary for that retained stream, so an extension must not
launch work on an untracked stream or publish a product before its work is
ordered into the supplied context. The lower-level model dispatches run inside
that boundary.

The example uses the maintained phase-only path operation, which writes the
current atmosphere OPD into the exact path-local `PupilFunction`. A genuinely
atmosphere-independent model instead passes the qualified
`AdaptiveOpticsSim.Plant.AtmosphereIndependentPath()` marker. Do not use that marker
as a fallback for an unsupported atmospheric field or layer-aware model.
Those models provide a concrete materialization owner and extend the qualified
`validate_path_materialization_binding`, `validate_path_materialization`, and
four-argument `materialize_path_input!` dispatches. The validator must check
the exact atmosphere, destination, source, backend, device, shape, and any
model-specific revision without mutating output. The mutating method may then
write only its bound caller-owned path input. This two-phase contract lets a
selection reject every invalid path before materializing the first one.

A `MultiLayerAtmosphereDefinition` or
`InfiniteMultiLayerAtmosphereDefinition` stored by a `PlantDefinition`
declares one stable `AtmosphereLayerID` per layer through its `layer_ids`
keyword. Ordinary numerical atmosphere constructors still permit omitted IDs
for non-plant work, but the cold plant declaration rejects missing or duplicate
stochastic-owner identities. A custom single-owner timed atmosphere also
implements its ordinary `initialize_atmosphere!` and
`evolve_atmosphere!` methods against `AbstractRNG`; the prepared plant supplies
the exact owner-bound RNG to those methods. To participate in target-local
partition preparation it also implements qualified
`Plant.partition_atmosphere_layer_ids(definition)` and returns a `Tuple` of
stable `AtmosphereLayerID` values. Return `()` for a single stochastic owner
with no independently identified layers. Core rejects non-tuples, non-layer
identities, and duplicates.

`input` is a path-local `PupilFunction`, a declared-plane `ElectricField` or
`IntensityMap`, or a concrete tuple of them. The prepared execution determines
which input planes and products it can consume; there is no implicit
resampling or propagation between entry and execution. `result` is an
acquisition-facing photon-rate `IntensityMap` or a concrete
tuple/`OpticalProductBundle` of such maps. The custom telescope must implement
the aperture revision, reflectivity, backend, and compute-device interfaces
consumed by those products. Every leaf in a multi-input or multi-result path
must share that backend and compute device.
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
qualified `AdaptiveOpticsSim.Plant.path_source_geometry_key`,
`AdaptiveOpticsSim.Plant.path_source_spectral_key`, and
`AdaptiveOpticsSim.Plant.path_source_radiometry_key` methods. Return run-owned,
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
function AdaptiveOpticsSim.Plant.prepare_acquisition_provider(
    model::MyAcquisitionDefinition,
    definition::AdaptiveOpticsSim.Plant.AcquisitionDefinition,
    path::AdaptiveOpticsSim.Plant.PreparedPathExecutor,
)
    AdaptiveOpticsSim.Plant.require_path_result(
        path; optical_model=model.required_optical_model)
    execution, observation, measurement = prepare_my_acquisition(model, path)
    products = AdaptiveOpticsSim.Plant.AcquisitionProducts(
        observation, measurement;
        metadata=my_product_metadata(model, path, observation, measurement))
    return prepare_full_optical_provider(execution, products)
end
```

Core constructs the `PreparedAcquisitionOwner` after validating the returned
provider against the exact path result. An acquisition extension does not
construct or return the owner itself.

Use qualified `AdaptiveOpticsSim.Plant.WFSOpticalPathExecution` to adapt an existing
Gate 0 WFS optical plan, `AdaptiveOpticsSim.Plant.FrameAcquisitionExecution` for a
frame detector plus a distinct caller-owned observation, and
`AdaptiveOpticsSim.Plant.WFSAcquisitionExecution` to compose already prepared WFS
acquisition and estimator plans. A deterministic concrete path execution type
must extend both the three-argument `execute_path!` dispatch and the qualified
`AdaptiveOpticsSim.Plant.validate_path_execution_binding(execution, input, result)`
seam. A stateful stochastic path can instead extend the four-argument form that
receives its prepared provider `AbstractRNG`. If it needs additional
independent device streams, extend qualified
`additional_path_rng_owner_roles` and `execute_path_rngs!`, then obtain each
declared stream directly with `rng_stream_state(group, Val(:role))`. A
different full-optical acquisition execution type similarly extends the
four-argument `execute_acquisition!` dispatch and
`AdaptiveOpticsSim.Plant.validate_acquisition_execution_binding(execution,
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

### Target-local partition preparation

A model that participates in preparation-only CPU/accelerator partitioning
implements the applicable qualified target-local seams. The path seam prepares
only co-located static data and reusable workspace; it receives no atmosphere
and must not invent a renderer, timeline, transfer, or executor:

```julia
function AdaptiveOpticsSim.Plant.prepare_target_local_path_resources(
    model::MyOpticalModelDefinition,
    definition::AdaptiveOpticsSim.Plant.OpticalPathDefinition,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    telescope::AdaptiveOpticsSim.Optics.AbstractTelescope,
    context,
)
    input, result, execution = prepare_my_target_local_optics(
        model, source, telescope)
    return AdaptiveOpticsSim.Plant.PreparedTargetLocalPathResources(
        definition, source, telescope, input, result, execution;
        context,
        optical_model=my_exact_model_key(model),
        propagation_model=my_exact_propagation_key(model),
        model_revisions=my_revision_key(model),
    )
end

function AdaptiveOpticsSim.Plant.prepare_target_local_acquisition_provider(
    model::MyAcquisitionDefinition,
    definition::AdaptiveOpticsSim.Plant.AcquisitionDefinition,
    path::AdaptiveOpticsSim.Plant.PreparedTargetLocalPathResources,
)
    AdaptiveOpticsSim.Plant.require_path_result(
        path; optical_model=model.required_optical_model)
    execution, observation, measurement =
        prepare_my_acquisition(model, path)
    products = AdaptiveOpticsSim.Plant.AcquisitionProducts(
        observation, measurement;
        metadata=my_product_metadata(model, path, observation, measurement))
    return AdaptiveOpticsSim.Plant.prepare_full_optical_provider(
        execution, products)
end
```

A controllable-optic model separately implements the fail-closed target-local
physical-preparation seam:

```julia
function AdaptiveOpticsSim.Plant.prepare_target_local_controllable_optic(
    model::MyControllableOpticModel,
    definition::AdaptiveOpticsSim.Plant.ControllableOpticDefinition,
    telescope::AdaptiveOpticsSim.Optics.AbstractTelescope,
    atmosphere_definition::AdaptiveOpticsSim.Atmospheres.AbstractTimedAtmosphereDefinition,
    target::AdaptiveOpticsSim.Backends.AbstractComputeDevice,
)
    return prepare_my_controllable_optic(model, definition, telescope, target)
end
```

It receives the cold atmosphere definition only for geometry or capability
decisions; it must not prepare a second timed atmosphere. It returns immutable
physical preparation data on the exact target. Core then invokes the existing
`prepare_controllable_optic_state` and
`prepare_controllable_optic_workspace` seams, validates all three physical
owners, and constructs one shared target-local owner per visible logical optic.
Provide exact `validate_controllable_optic_target`,
`validate_controllable_optic_state_target`,
`validate_controllable_optic_workspace_target`, and
`structural_resource_fact` methods. Unknown model/target combinations fail
closed.

Reuse a shared model-specific helper when ordinary `prepare_path_executor` and
target-local preparation construct the same optical data. Do not implement the
target-local seam by preparing and discarding a second timed atmosphere. The
returned execution/provider types must implement their normal exact-binding
and exact-target validators. They must also provide exact
`Plant.structural_resource_fact` dispatch for every owned resident array and
workspace; an unknown fact remains an admission failure rather than an
estimated byte count.

The caller, not core, resolves placement:

```julia
assignment = AdaptiveOpticsSim.Plant.resolve_plant_partition_assignment(
    definition,
    host_target,
    :wfs => accelerator_target,
    :science => host_target,
)
prepared = AdaptiveOpticsSim.Plant.prepare_plant_partitions(
    assignment;
    run_seed=0x1234,
    command_authority_target=host_target,
    command_endpoints=endpoint_configurations,
)
authority =
    AdaptiveOpticsSim.Plant.prepared_atmosphere_authority(prepared)
gpu_partition = AdaptiveOpticsSim.Plant.prepared_partition(
    prepared, accelerator_target)
```

Every declared path must appear exactly once. Acquisition placement is derived
from its path, so an acquisition cannot be split onto another target. The
assignment admits host resources and at most one exact accelerator target in
Gate 9A. `PreparedPlantPartitions` is an inspectable, non-executable preparation
result. Cold preparation may defensively copy declared static data, such as a
sampled OPD, into each exact target that uses it. The result performs no runtime
atmosphere publication, command admission, authoritative command publication,
product movement or automatic handoff, scheduling, or placement planning. If the
definition declares command endpoints, `endpoint_configurations` must contain
one complete `CommandEndpointConfiguration` per endpoint. Initial and optional
safe array values are cold host-resident configuration: Core seals a detached
host copy before deriving independent exact-target authority and replica
storage, and rejects device-resident configuration arrays. Capacity, history,
safe-command, silence, and admission policy remain authority-side concerns.

### Cross-device handoff transfers

A semantic payload contract subtypes qualified
`AdaptiveOpticsSim.Plant.AbstractHandoffPayloadContract{M}` and extends
`handoff_payload_eltype`, `handoff_payload_axes`, and
`validate_handoff_publication`. The element type must be concrete and isbits;
the axes and immutable inline-stored publication type are fixed at preparation.
The caller and contract must not mutate any semantic state referenced by a
publication between successful submission and reclamation. Construct a
fixed-capacity transfer with `Plant.prepare_cross_domain_handoff` using paired,
nonaliasing caller-provided source and destination arrays on distinct exact
host/accelerator targets.

An accelerator extension opts a concrete storage pair into this boundary by
implementing `Backends._prepare_array_transfer(source, destination,
source_context, destination_context)`,
`Backends._submit_prepared_array_transfer!`, and one nonblocking
`Backends._observe_prepared_array_transfer_completion!` observation. The
deterministic serial oracle additionally requires
`Backends._complete_prepared_array_transfer!`, which blocks inside the backend
until the submitted transfer is quiescent and returns only completed or
failed. Prepared backend state retains the exact context or event state needed
to launch and complete the transfer while restoring the caller's previous
device context.
Submission returns `_PreparedArrayTransferSubmitted`, or returns
`_PreparedArrayTransferSubmissionFailed` only after guaranteeing that no
transfer can access either slot. One completion observation returns
`_PreparedArrayTransferPending`, `_PreparedArrayTransferCompleted`, or
`_PreparedArrayTransferCompletionFailed`; the failed result likewise
guarantees that the transfer is quiescent. An exception or unsupported result
leaves the slot in fail-stop `HandoffTransferUncertain` state and it cannot be
reclaimed. The caller must disarm and dispose of that prepared handoff unless a
later backend-specific recovery contract explicitly establishes quiescence.
There is no generic `copyto!` fallback, global device synchronization, polling
loop, or implicit host staging in Core.

Before submission the single producer exclusively owns the selected source
slot. Successful submission transfers source ownership to the handoff until
terminal reclamation. Preparation gives the handoff exclusive ownership of
destination slots; a consumer accesses one only after successful completion
through `try_borrow_completed_handoff!`. HIL separately owns slot leasing,
descriptor rings, return credit, backpressure, scheduling, and wait policy.
One HIL transfer Agent should serialize every Core lifecycle method; producer
and consumer Agents exchange descriptors or leases with that owner instead of
mutating the handoff concurrently.
The current maintained implementation validates the seam with a deterministic
fake accelerator; it does not qualify CUDA or AMDGPU transfers.

### Authority-owned pupil-OPD routes

Core specializes the handoff boundary for the native atmospheric pupil-OPD
product. `Plant.prepare_pupil_opd_publication_route(partitions, path_id)`
accepts only a target-local path whose input is a native `PupilFunction` and
whose atmosphere materializer is a native `PreparedPupilOPDMaterialization`.
It binds the exact prepared partition set, atmosphere authority, target-local
path, frozen source, element type, axes, and compute targets. Same-target
routes materialize directly into the path input. Remote routes prepare one
authority-local source pupil and one target-local destination slot, then
transfer only the OPD array through the generic handoff seam above.

A path model must opt in deliberately by defining the qualified method

```julia
AdaptiveOpticsSim.Plant.prepare_pupil_opd_publication_materialization(
    ::MyPathModel,
    atmosphere::AdaptiveOpticsSim.Atmospheres.AbstractTimedAtmosphere,
    telescope::AdaptiveOpticsSim.Optics.Telescope,
    source::AdaptiveOpticsSim.Optics.AbstractSource,
    pupil::AdaptiveOpticsSim.Optics.PupilFunction,
) = AdaptiveOpticsSim.Plant.prepare_pupil_opd_materialization(
    atmosphere, telescope, source, pupil)
```

Core requires the returned value to be a native
`PreparedPupilOPDMaterialization` bound to those exact supplied objects.
Its constructor is intentionally closed: use
`prepare_pupil_opd_materialization`, which binds the exact telescope and pupil
and prepares a renderer for the requested frozen source geometry. Models
without this method, values of another materialization type, and stale,
foreign, or disguised optical bindings fail with structured preparation
errors. Do not opt in a model merely because its path input happens to be a
`PupilFunction`; the model must actually define atmosphere-derived pupil OPD
as this publication product.

Drive one publication through materialize, submit, nonblocking completion,
apply, and reclaim before reusing the route. One external owner must serialize
those calls and must not advance atmosphere state until every due route has
applied. An unexpected backend result makes the route uncertain and
non-reclaimable. Core does not poll, retain arbitrary atmosphere history, or
compose the all-routes epoch barrier.

The immutable publication carries an opaque route identity rather than the
source object or prepared resource graph. That identity binds the exact frozen
source and materialization for the route; callers compare it only through
`pupil_opd_publication_route_identity` and do not reconstruct it.

Do not use this seam to infer transfer for an arbitrary path input. Pupil
support and amplitude, sampled static aberrations, controllable-optic surfaces,
autonomous path optics, propagation workspaces, and path results remain
target-local. A future field, intensity, extended-source, asterism,
multi-component, or illumination product requires its own semantic publication
and apply contract; it must not opt in through a generic `copyto!` fallback.

Plant invokes a full-optical path executor only after atmosphere
materialization, prepared native sampled aberrations, controllable surfaces,
and autonomous path-local optics have formed its caller-owned input. A typed
external executor may therefore hand that residual pupil or field to
`Proper.jl` or another optical owner without importing it into core. Do not
reapply the same native aberration or controllable surface inside that
executor.

A custom reduced-order provider instead returns qualified
`AdaptiveOpticsSim.Plant.PreparedAcquisitionProvider(implementation, products)` and
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
Those four methods establish the schedule-free provider boundary only. They do
not grant a custom provider access to composed effective-command state.

For the maintained scheduled linear model, construct qualified
`HarmonicDisturbanceModel`, one `ReducedOrderCommandResponse` per currently
visible prepared endpoint, and `LinearReducedOrderAcquisitionModel`, then bind
the acquisition with `DirectMeasurementAcquisitionDefinition`. Event-loop
preparation validates exact command type, dimensions, units, sign convention,
basis, basis revision, backend, and device before resolving fixed endpoint
slots. The lifecycle integrates the held direct measurement over exposure and
publishes it after the declared readout delay. A reduced-order-only path does
not execute its unused full-optical path. Custom scheduled reduced-order
providers remain unsupported until a bounded event-sampling extension contract
is specified; extensions must not reach into event-loop internals to imitate
one.
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
`prepare_plant(definition, target; run_seed, rng_derivation_version)` owns all
stateful streams, so neither selected-execution method accepts an RNG argument.
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

AdaptiveOpticsSim.Plant.illumination_combination(
    ::Type{<:PreparedMyIllumination}) =
    AdaptiveOpticsSim.Plant.SingleIllumination()

function AdaptiveOpticsSim.Plant.prepare_illumination_evaluator(
    definition::MyIlluminationDefinition,
    destination::IntensityMap,
    ::DetectorInputIlluminationEntry,
)
    state = MyIlluminationState(0, zero(eltype(destination.values)))
    return PreparedMyIllumination(definition, state,
        AdaptiveOpticsSim.Backends.backend(destination),
        AdaptiveOpticsSim.Backends.compute_device(destination.values))
end

function AdaptiveOpticsSim.Plant.validate_illumination_evaluator_binding(
    evaluator::PreparedMyIllumination,
    destination::IntensityMap,
    ::DetectorInputIlluminationEntry,
)
    typeof(AdaptiveOpticsSim.Backends.backend(destination)) ===
            typeof(evaluator.backend) ||
        throw(PlantPreparationError(:illumination, :backend,
            "prepared illumination backend changed"))
    AdaptiveOpticsSim.Backends.compute_device(destination.values) ==
            evaluator.device ||
        throw(PlantPreparationError(:illumination, :device,
            "prepared illumination device changed"))
    return nothing
end

function AdaptiveOpticsSim.Plant.evaluate_illumination!(
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
`AdaptiveOpticsSim.Plant.SingleIllumination()`,
`AdaptiveOpticsSim.Plant.ExclusiveIlluminationSelection()`,
`CoherentFieldCombination()`, or
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
`AdaptiveOpticsSim.WavefrontSensors.detector_response_calibration_signature(model, seed)`
method covering every parameter that can change the deterministic response.
Without that method, WFS calibration fails closed rather than caching a
reference by response type alone.

Use `bits` for quantization depth and `output_type` for the Julia element type
returned at a HIL/RTC boundary.

## Wavefront Sensors

`WavefrontSensors` owns the common WFS contracts. Family-specific
implementation stays near the family, usually in the corresponding directory
under `src/wfs/`, and extends generics imported explicitly from that owner.

New physical WFS work should first implement the prepared semantic stages:

- `prepare_wfs_optics(model, input, output)` validates an explicit
  caller-owned pupil function or electric field and one or more
  detector-plane photon-arrival-rate outputs; the returned concrete prepared
  owner is executed by `form_wfs_optical_products!`
- `prepare_wfs_acquisition(model, optical_products, observation)` validates
  detector mappings, explicit duration, backend/device placement, and one or
  more caller-owned `WFSObservation` destinations; the returned concrete
  prepared owner is executed by `acquire_wfs_observation!` with an explicit
  caller-owned RNG
- `prepare_wfs_estimation(model, observation, measurement)` validates the
  estimator and caller-owned `WFSMeasurement`; the returned concrete prepared
  owner is executed by `estimate_wfs_measurement!` and declares
  `AcquiredObservationPath()` through `wfs_measurement_path`

A geometric or reduced-order provider may instead prepare estimation directly
from a `PupilFunction` or pupil-plane `ElectricField`. It must declare
`DirectMeasurementPath()` and must not create unused rate, detector, or
observation storage. `WFSObservation` supports scalar `Ref` storage and arrays
of any rank, as does `WFSMeasurement`; use concrete tuples for multiple
observations. Preserve incompatible spectral or branch rate products in
`OpticalProductBundle`. Bundle membership is fixed at construction; the arrays
owned by its product leaves remain mutable destinations for prepared WFS
optics.

The qualified nominal plan interfaces are
`WavefrontSensors.AbstractWFSOpticsPlan`,
`WavefrontSensors.AbstractWFSAcquisitionPlan`, and
`WavefrontSensors.AbstractWFSEstimationPlan`. Subtype the interface for the
stage contract owned by the implementation, but keep that run-immutable plan
as one field of a distinct exact prepared owner. A plan may retain validated
physical coefficients, dimensions, mappings, calibration revisions, and
backend/device requirements. It must not retain replaceable scratch,
persistent scientific state, caller-visible products, or incidental exact
array identity. Generic detector acquisition owners expose their plan through
the qualified `WavefrontSensors.wfs_acquisition_plan` accessor. Family optical
and estimator accessors are introduced with their owner migrations rather
than through a universal stage-plan API. Static detector fan-out validates all
component contracts and cross-component ownership before it commits any
detector owner; a later component failure must leave every prior prepared
owner valid.

Extensions should call the qualified validation seams
`AdaptiveOpticsSim.WavefrontSensors.validate_wfs_optical_input`,
`AdaptiveOpticsSim.WavefrontSensors.validate_wfs_optical_products`,
`AdaptiveOpticsSim.WavefrontSensors.validate_wfs_observation` or
`AdaptiveOpticsSim.WavefrontSensors.validate_wfs_observations`, and
`AdaptiveOpticsSim.WavefrontSensors.validate_wfs_measurement` as applicable
before returning a prepared owner. Every new prepared owner also implements its
corresponding exact binding validator:

- `AdaptiveOpticsSim.WavefrontSensors.validate_wfs_optics_binding`
- `AdaptiveOpticsSim.WavefrontSensors.validate_wfs_acquisition_binding`
- `AdaptiveOpticsSim.WavefrontSensors.validate_wfs_estimation_binding`

The containing plant owner calls these validators during construction, and the
stage calls them again before mutation. They remain qualified extension APIs
rather than ordinary exported workflow names.

Prepared types should contain concrete immutable plans/params and separately
typed single-writer workspace, detector, calibration, and RNG state. Bind exact
array/state identities, validate compute device as well as semantic backend,
and create detector-output aliases or packed views only after detector buffers
are prepared. Repeated execution must not resize, rebuild metadata, query a
device, copy to the host, or select stages through an abstract container.
Raise `WFSPreparationError(stage, reason, msg)` for preparation
incompatibility and for execution-time prepared-binding rejection before any
destination mutation. The protocol stages are `:wfs_optics`,
`:acquisition`, and `:estimation`; `reason` is an open extension identifier.

Concrete WFS family types may provide a composed `measure!` workflow over the
same prepared stages. That workflow is an ordinary high-level API, not a second
ownership model: optical inputs remain caller-owned, stage workspaces remain
sensor-owned and preallocated, and the implementation dispatches through the
prepared seams rather than branching on sensor families.

## Optical-Path Source Roles

Each `OpticalPathDefinition` owns its explicit frozen source. Do not assume one
source serves every WFS and science path, and do not infer source roles from
tuple position. A prepared plant may contain several NGS/LGS sensing
directions, science sources, calibration sources, and coronagraph paths.

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

`AdaptiveOpticsSim.Atmospheres.AbstractTimedAtmosphere` and its model hooks are an advanced
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
profile data should specialize the qualified
`AdaptiveOpticsSim.Optics.freeze_source` seam and return a run-owned copy.

Plant preparation freezes mutable source inputs, retains one prepared
atmosphere renderer per full-optical direction, and preserves explicit path
identity. A WFS extension consumes the path's declared `PupilFunction` or
`ElectricField` through the ordinary prepared WFS stage contracts; it must not
depend on a hidden shared-runtime hook.

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

Keep command application explicit and ordered in the plant. Every physical
optic retains independent endpoint state and timing, even when several
co-conjugated surfaces add on the same path.

## Controllers And Reconstructors

Controllers own temporal behavior. Reconstructors own measurement-to-command
operators.

Import the maintained surface with `using AdaptiveOpticsSim.Control`.

Use this split:

- reconstructors expose the command operator and any diagnostics needed by
  calibration or validation
- controllers own integrator state, gains, leakage, saturation, or future
  temporal filters
- model-specific loops or HIL integration call generic accessors and `!`
  methods

Prefer concrete state and workspace fields. Avoid hidden allocations in the
per-step control path.

`DiscreteIntegratorController` is a prepared execution owner rather than an
algorithm plan. Its qualified `discrete_integrator_plan`,
`discrete_integrator_state`, and `discrete_integrator_workspace` accessors
expose the run-immutable coefficients, persistent mathematical history, and
replaceable scratch as separate concrete owners. A graph adapter may package
the state and workspace together for an external interface that provides one
mutable execution slot, but it must not reclassify persistent state as
replaceable workspace.

For large calibration/control surfaces:

- use `AdaptiveOpticsSim.Calibration.interaction_matrix!` with caller-owned
  storage when
  the calibration matrix should live in a memory map, shared buffer, or
  backend-native allocation
- use `FactorizedReconstructor` when a validated truncated SVD rank can bound
  runtime storage and compute without materializing the dense inverse
- use `ControlledReconstructor` to compose temporal state with any maintained
  reconstructor while preserving the control-storage contract

Rank selection is an optical/control validation decision. A benchmark win from
truncation is not evidence that discarded modes are acceptable for a real
plant. The dense `ModalReconstructor` remains the full-rank accuracy baseline.

Accelerator reconstructors must implement
`AdaptiveOpticsSim.Control.runtime_reconstructor_storage(reconstructor)` and return a tuple containing
their hot-path control matrices and workspaces. A model or Plant extension uses
this contract to reject mixed host/device control paths during preparation.
Return `()` only for a truly backend-agnostic operator. Synchronize at an
explicit observation, transport, or measurement boundary rather than inside an
otherwise device-resident chain.

Stateful reconstructors should additionally implement
`AdaptiveOpticsSim.Control.runtime_reconstructor_ownership_roots(reconstructor)` and
return their mutable workspaces. This prevents two threaded ensemble members
from sharing operator scratch state. Stateful controllers used inside
`ControlledReconstructor` implement
`AdaptiveOpticsSim.Control.runtime_controller_storage` for backend
residency; controller objects are conservatively treated as mutable ownership
roots by default.

## Ensemble Scheduling

Import this domain with `using AdaptiveOpticsSim.Ensembles`.

`SimulationEnsemble` schedules independent model boundaries at coarse
granularity. Keep it outside an external-RTC HIL loop: the Plant event loop and
its single-writer owners must not acquire task-creation or task-graph jitter.

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

Ensemble construction rejects shared mutable ownership. A mutable member owns
itself by default. An immutable or mutable wrapper around shared state should
implement
`AdaptiveOpticsSim.Ensembles.ensemble_ownership_roots(wrapper)` and return the
mutable plant/state objects that cannot safely be updated by another member at
the same time. Immutable values with no mutable owners, such as scalar sweep
points, require no extension. The scheduler operation passed to
`run_ensemble!` must mutate and return normally only after the member is safe
for the next coarse operation.

Avoid nested parallelism. When using a threaded ensemble policy, normally keep
BLAS and FFT-provider thread counts at one and let the ensemble own the coarse
parallelism.

## Backend And Allocation Rules

Extension code should follow the package-wide backend rules:

- parametrize by numeric type and array backend where practical
- avoid scalar indexing on GPU arrays
- preallocate workspaces for hot paths
- use explicit `!` methods when mutating state
- centralize evolving RNG state in the caller or an explicit prepared state
  owner; workspace may own only replaceable random-generation scratch
- use multiple dispatch or traits instead of `isa` chains

An accelerator-array extension registers both a semantic backend family and a
concrete compute-device identifier. Implement
`AdaptiveOpticsSim.Backends.array_backend_selector(::Type{<:MyArray})` with a
zero-state `AdaptiveOpticsSim.Backends.AbstractArrayBackend` selector and
`AdaptiveOpticsSim.Backends.compute_device_identifier(::MyArray)` with a
non-`nothing`, allocation-free isbits identifier that distinguishes every
simultaneously addressable device in that backend runtime.
`AdaptiveOpticsSim.Backends.compute_device(array)` combines those two
contracts into an `AdaptiveOpticsSim.Backends.AcceleratorComputeDevice`; its
methods for supported array views preserve the parent identity. Do not use a
backend family alone as device zero, and do not encode task, stream, or future
placement ownership in the identifier.

An extension that defines its own
`AdaptiveOpticsSim.Backends.AbstractComputeDevice` subtype instead of using
`AdaptiveOpticsSim.Backends.AcceleratorComputeDevice` must also implement
`AdaptiveOpticsSim.Backends.compute_device_backend(::MyComputeDevice)`. That
accessor returns the corresponding semantic
`AdaptiveOpticsSim.Backends.AbstractArrayBackend` selector and must agree with
`AdaptiveOpticsSim.Backends.backend(array)` wherever the device is used in
prepared Plant contracts.

If an extension needs a new reusable behavior, add a small generic seam and one
family implementation first. Then add the second implementation when another
family actually needs it.

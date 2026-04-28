# Adding A Deformable Mirror

Status: active

## Purpose

This guide is the maintained authoring path for adding a new deformable-mirror
family, influence model, topology, or actuator-behavior layer to
`AdaptiveOpticsSim.jl`.

Use it when you need to:

- add a new reusable DM topology
- add a new static influence-basis model
- add a new actuator-behavior model
- decide whether a manufacturer artifact belongs in core or in an external
  conversion tool

Use together with:

- [api-reference.md](./api-reference.md)
- [model-cookbook.md](./model-cookbook.md)
- [dm-interface-evolution-plan-2026-04.md](./dm-interface-evolution-plan-2026-04.md)

## Design Rule

The DM subsystem is organized around this split:

- topology owns actuator layout and active-actuator structure
- influence model owns the static surface basis
- actuator model owns command-space behavior before surface application

In practice:

- if a feature changes where actuators are or which are valid, make it a
  topology
- if a feature changes the static basis from actuator command to OPD, make it
  an influence model
- if a feature changes command behavior before the basis is applied, make it an
  actuator model
- do not add new DM behavior as ad hoc branches in the `DeformableMirror`
  constructor

## First Decision: Which Layer Changes?

Start by deciding which DM layer the new feature actually belongs to.

### Topology

Use `AbstractDMTopology` when the feature defines actuator layout, indexing, or
active-mask structure.

Current maintained examples:

- [ActuatorGridTopology](../src/optics/deformable_mirror.jl)
- [SampledActuatorTopology](../src/optics/deformable_mirror.jl)

Topology should answer:

- where actuators live
- which actuators are valid
- how active command channels map to the full actuator set

Topology should not own:

- surface basis data
- actuator clipping, health, hysteresis, or dynamics

### Static Influence Basis

Use `AbstractDMInfluenceModel` when the feature defines the static mapping from
actuator command to sampled OPD basis.

Current maintained examples:

- `GaussianInfluenceWidth`
- `GaussianMechanicalCoupling`
- `DenseInfluenceMatrix`
- `MeasuredInfluenceFunctions`

Influence models should answer:

- what sampled OPD basis each actuator produces
- whether the model supports the separable Gaussian fast path
- whether it supports misregistration-identification workflows

Influence models should not own:

- actuator layout metadata that belongs in topology
- actuator clipping, health, or other command-space behavior

### Actuator Behavior

Use `AbstractDMActuatorModel` when the feature changes the actuator command
vector before the static basis is applied.

Current maintained examples:

- `LinearStaticActuators`
- `ClippedActuators`
- `ActuatorHealthMap`
- `CompositeDMActuatorModel`

Actuator models should answer:

- how command values are transformed before OPD application
- whether certain actuators are weakened, clipped, or effectively dead

Actuator models should not own:

- basis matrices
- sampled OPD surfaces
- actuator coordinates

## File Placement

Most DM extensions should be implemented in
[deformable_mirror.jl](../src/optics/deformable_mirror.jl), because the DM
composition layer is intentionally centralized there.

Typical steps:

1. add the new topology, influence model, or actuator model in
   [deformable_mirror.jl](../src/optics/deformable_mirror.jl)
2. export the new public type or accessor from
   [AdaptiveOpticsSim.jl](../src/AdaptiveOpticsSim.jl) if it is meant to be
   public
3. add tests in
   [control_and_runtime.jl](../test/testsets/control_and_runtime.jl)
4. add calibration restrictions or compatibility hooks only if they are
   structurally required
5. update user docs if the new surface is public

Do not add vendor file parsing or format-specific I/O to core. The core package
should consume normalized Julia data structures, not raw FITS or proprietary
manufacturer formats.

## Adding A New Topology

### 1. Define The Topology Type

Add a concrete subtype of `AbstractDMTopology`.

Example shape:

```julia
struct MyDMTopology{T,M1,M2,V1,V2,MD} <: AbstractDMTopology
    coords::M1
    active_coords::M2
    valid_actuators::V1
    active_indices::V2
    metadata::MD
end
```

Keep topology focused on geometry and indexing.

### 2. Provide The Topology Accessors

At minimum, define the accessors used by the composed DM layer:

- `actuator_coordinates(topology)`
- `valid_actuator_mask(topology)`
- `active_actuator_indices(topology)`
- `topology_metadata(topology)`

And, only when appropriate:

- `topology_axis_count(topology)`
- `supports_separable_topology(topology)`

If the topology is not a regular fully active grid, it should usually not claim
separable support.

### 3. Validate Early

Topology constructors should validate:

- coordinate dimensions
- valid-mask length
- at least one active actuator

Throw structured errors such as `InvalidConfiguration` or
`DimensionMismatchError`.

## Adding A New Influence Model

### 1. Choose Whether It Is Analytic Or Sampled

Analytic models should own only the minimal parameterization required to build
the sampled basis.

Sampled models should store the normalized sampled basis directly.

Examples:

- analytic: `GaussianInfluenceWidth`
- sampled: `MeasuredInfluenceFunctions`

### 2. Implement Conversion And Validation

New influence models must participate in the existing DM normalization path.
That usually means adding:

- conversion into the resolved numeric/backend form
- dimension validation against the chosen topology and telescope resolution

The core rule is:

- author-facing types may be convenient
- runtime state must end up in the package’s normalized dense sampled basis

### 3. Implement The Basis Build Path

Every influence model must support basis construction through the DM build
layer.

For analytic models, build the sampled basis from topology and pupil geometry.

For sampled models, validate and copy into the runtime backend.

### 4. Declare Structural Limits Through Dispatch

If a new influence model does not support:

- separable application
- misregistration-identification workflows
- mechanical-coupling conversion

say so through dispatch helpers rather than special-case branches.

Current examples to mirror:

- `supports_separable_influence(...)`
- `supports_dm_misregistration_identification(...)`

## Adding A New Actuator Model

### 1. Keep It Command-Space Only

Actuator models should transform the actuator command vector, not the OPD
surface.

If the model needs to alter the basis itself, it probably belongs in the
influence layer instead.

### 2. Implement In-Place Command Preparation

New actuator models should participate through the existing command-preparation
path, so the hot loop remains preallocated and explicit.

Examples:

- clipping
- gain/health maps
- future hysteresis or dynamics layers

### 3. Compose When Possible

If a behavior can be applied as a stage in sequence, prefer composing it with
`CompositeDMActuatorModel` instead of inventing a family-specific monolith.

## Manufacturer Influence Data

Manufacturer artifacts such as ALPAO or Boston influence functions are valid
targets for the DM interface, but they should enter core only after conversion
to package-native structures.

The maintained core path is:

1. external tool reads the vendor artifact
2. external tool converts it into:
   - actuator coordinates / validity mask
   - sampled influence basis
   - optional metadata
3. core consumes those via:
   - `SampledActuatorTopology(...)`
   - `MeasuredInfluenceFunctions(...; metadata=...)`

This means:

- `MeasuredInfluenceFunctions` is the public landing zone for converted
  manufacturer-style sampled bases
- FITS parsing, vendor naming, and file-normalization policy should stay
  outside the core package

## Reuse Rule

This is the main architectural test for DM work:

- if a feature affects only one DM family today, it may remain local
- if a second DM family needs it, factor it into topology, influence, or
  actuator-model reuse immediately unless there is a strong reason not to

Good reuse examples:

- sampled topologies
- measured influence bases
- clipping and health maps

Bad pattern:

- adding vendor-specific booleans or branches to the generic `DeformableMirror`
  constructor

## Validation Checklist

Every new DM extension should add evidence in
[control_and_runtime.jl](../test/testsets/control_and_runtime.jl), and when
relevant in
[calibration_and_analysis.jl](../test/testsets/calibration_and_analysis.jl).

Minimum expectations:

1. construction works and validates bad inputs correctly
2. command length and topology metadata are consistent
3. `apply_opd!` produces finite, correctly shaped output
4. backend behavior works on CPU and maintained GPU backends when supported
5. calibration restrictions are explicit when the new DM family is not
   compatible with a workflow

Add cookbook or user-guide coverage if the new surface is public.

## What Not To Do

- do not add raw FITS or vendor file parsing into core
- do not add new DM behavior as `isa` branches in generic code
- do not bypass the topology / influence / actuator split for convenience
- do not hide structural incompatibilities; reject them explicitly
- do not duplicate existing Gaussian or sampled-basis logic under a new name

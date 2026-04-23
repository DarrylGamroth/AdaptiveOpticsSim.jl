Status: completed

## Purpose

This plan defines the next interface evolution for `DeformableMirror`.

The current DM surface is good for:

- Gaussian continuous-face-sheet approximations
- mechanical-coupling parameterization of that Gaussian family
- explicit sampled dense mode matrices for advanced users

It is not yet the best long-term authoring surface for:

- measured manufacturer influence functions
- explicit actuator layout metadata
- valid/slaved actuator masks
- segmented or modal mirror families
- actuator technology constraints such as saturation, dead actuators,
  hysteresis, creep, or dynamics

This plan improves the interface before adding vendor-specific or
technology-specific DM families.

## Motivation

The package should be able to accept manufacturer-style influence-function data
cleanly.

One concrete motivating example is an ALPAO-style influence-function artifact
outside the repository:

- `/mnt/datadrive/DATA/dm/AX307_Influences.fits`

That kind of input should not be forced through a raw `DenseInfluenceMatrix`
alone, because a bare sampled matrix loses important metadata such as:

- actuator coordinates and pitch
- actuator indexing and topology
- valid/slaved/dead actuator masks
- family and technology semantics
- provenance for the sampled basis

`DenseInfluenceMatrix` remains a valid low-level escape hatch, but it should
not be the only expressive path for measured DM data.

## Current Surface

The maintained DM influence-model surface in
[`src/Optics/deformable_mirror.jl`](../src/Optics/deformable_mirror.jl) is:

- `AbstractDMInfluenceModel`
- `GaussianInfluenceWidth`
- `GaussianMechanicalCoupling`
- `DenseInfluenceMatrix`

This conflates three distinct concerns:

1. actuator topology
2. static influence basis
3. actuator technology behavior

That makes the current design harder to extend cleanly.

## Target Interface Shape

Evolve the DM stack toward explicit composition of:

1. `AbstractDMTopology`
   - actuator coordinates
   - actuator indexing
   - pupil-valid mask
   - optional slaving map
   - optional metadata such as pitch and nominal pupil span

2. `AbstractDMStaticInfluenceModel`
   - Gaussian continuous-face-sheet family
   - measured sampled influence-function family
   - segmented/modal surface family later

3. `AbstractDMActuatorModel`
   - linear static actuator behavior first
   - saturation/clipping
   - dead/weak actuators
   - hysteresis/creep later
   - dynamic response later

Then `DeformableMirror` becomes a composition of those pieces rather than a
single influence-family choice.

## Public API Direction

The concise existing constructor must remain supported:

```julia
DeformableMirror(tel; n_act=16, influence_width=0.3)
DeformableMirror(tel; n_act=16, mechanical_coupling=0.08)
DeformableMirror(tel; n_act=16, influence_model=DenseInfluenceMatrix(modes))
```

The new surface should add, not replace, higher-fidelity construction paths:

```julia
DeformableMirror(tel;
    topology=ActuatorGridTopology(...),
    influence_model=GaussianInfluenceWidth(0.3),
    actuator_model=LinearStaticActuators(),
)

DeformableMirror(tel;
    topology=SampledActuatorTopology(coords; valid_actuators=mask, metadata=meta),
    influence_model=MeasuredInfluenceFunctions(sampled_modes; metadata=meta),
    actuator_model=LinearStaticActuators(),
)
```

These names are placeholders for the design direction, not frozen API.

## Phases

### DMIE-1: Introduce topology objects

Add a topology layer without changing maintained behavior yet.

Candidate maintained first types:

- `ActuatorGridTopology`
- `SampledActuatorTopology`

Requirements:

- explicit actuator count
- explicit coordinates
- optional valid-actuator mask
- optional slaving metadata
- conversion to the current dense command ordering

Success criteria:

- existing `n_act` grid-based constructor still works unchanged
- DM command-layout code can obtain actuator count and ordering from topology
- Gaussian DM path can be expressed through `ActuatorGridTopology`

### DMIE-2: Split influence-model role

Rename or refactor the current influence layer toward static sampled surface
bases.

Candidate maintained first types:

- `GaussianInfluenceWidth`
- `GaussianMechanicalCoupling`
- `MeasuredInfluenceFunctions`
- `DenseInfluenceMatrix` retained as low-level basis wrapper

Requirements:

- measured influence-function objects carry topology-aligned metadata
- conversion to the internal dense sampled basis remains centralized
- misregistration rules are explicit for measured bases

Success criteria:

- Gaussian and dense existing paths still pass unchanged
- measured influence-function inputs can be represented without metadata loss

### DMIE-3: Add actuator-model layer

Introduce a separate actuator behavior seam.

Candidate maintained first types:

- `LinearStaticActuators`
- `ClippedActuators`
- `MaskedActuators` or `ActuatorHealthMap`

Requirements:

- no change to current behavior when the default actuator model is used
- command preprocessing stays allocation-free in the hot path
- actuator constraints are applied before OPD synthesis

Success criteria:

- dead or clipped actuators can be modeled without modifying the influence
  basis
- the same static influence basis can be reused with different actuator models

### DMIE-4: Add measured influence-function import boundary

Do not bake FITS or vendor formats into core.

Instead:

- define a package-native measured-basis object
- add a documented external conversion path from manufacturer artifacts
- keep any format-specific loaders optional or external

For the motivating ALPAO-style case, the supported workflow should be:

1. external tool reads `AX307_Influences.fits`
2. external tool writes a package-native sampled influence basis artifact
3. `AdaptiveOpticsSim.jl` consumes the native sampled artifact through the
   measured influence-model API

Success criteria:

- no hard dependency on FITS in core
- measured DM data can be used reproducibly in tests and examples

### DMIE-5: Validation and conformance

Add explicit contract tests for:

- topology ordering and command-layout correctness
- Gaussian vs measured-basis equivalence where expected
- dead/slaved actuator handling
- backend parity on the maintained CPU/CUDA/AMDGPU DM surfaces
- calibration compatibility:
  - interaction matrix
  - modal basis
  - misregistration identification where supported

Success criteria:

- the new interface does not regress current controllable-optic and GPU
  validation surfaces
- unsupported algorithms fail structurally and explicitly

## Non-Goals

This plan does not yet include:

- vendor-specific constructors like `ALPAODM(...)` or `BostonMicromachinesDM(...)`
- dynamic actuator models in the first slice
- segmented mirror families in the first slice
- core FITS I/O support

Those should sit on top of the improved interface, not be used to define it.

## Recommended Order

1. `DMIE-1` topology layer
2. `DMIE-2` measured static influence layer
3. `DMIE-3` actuator-model layer
4. `DMIE-5` validation closeout
5. `DMIE-4` optional external conversion docs/examples in parallel where useful

That order is deliberate:

- topology is the missing foundation
- measured influence functions are the immediate practical need
- actuator behavior should be layered on a clean topology/basis split
- validation must keep pace with each new seam

## Immediate Next Slice

The first implementation slice should be:

- topology objects
- measured influence-function model object
- internal conversion to the current dense runtime basis
- no change yet to the maintained Gaussian constructor path

That gives the package a clean landing zone for manufacturer data without
forcing a full DM rewrite in one step.

## Outcome

Completed in the maintained package surface:

- explicit topology layer:
  - `ActuatorGridTopology`
  - `SampledActuatorTopology`
- explicit sampled influence-basis layer:
  - existing Gaussian families
  - `DenseInfluenceMatrix`
  - `MeasuredInfluenceFunctions`
- explicit actuator-behavior layer:
  - `LinearStaticActuators`
  - `ClippedActuators`
  - `ActuatorHealthMap`
  - `CompositeDMActuatorModel`
- `DeformableMirror` constructor composition through:
  - `topology=...`
  - `influence_model=...`
  - `actuator_model=...`

The concise Gaussian `n_act=...` constructor remains supported unchanged.

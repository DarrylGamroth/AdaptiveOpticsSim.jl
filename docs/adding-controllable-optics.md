# Adding A Controllable Optic

Status: active

## Purpose

This guide is the maintained authoring path for adding a new controllable-optic
family to `AdaptiveOpticsSim.jl`.

Use it when you need to:

- add a new low-order or modal controllable optic
- add a new convenience wrapper over an existing controllable-optic basis
- decide whether a surface belongs in `ModalControllableOptic`,
  `DeformableMirror`, or a new optic family

Use together with:

- [api-reference.md](./api-reference.md)
- [model-cookbook.md](./model-cookbook.md)
- [adding-deformable-mirrors.md](./adding-deformable-mirrors.md)

## Design Rule

The controllable-optic subsystem is organized around this split:

- modal or low-order surfaces should reuse `ModalControllableOptic` when the
  surface is just a basis plus coefficients
- DM-like actuator-grid surfaces should use `DeformableMirror`
- composite application should be expressed through `CompositeControllableOptic`
  rather than ad hoc ordered application logic

In practice:

- if the surface is a small modal basis, prefer `AbstractModalOpticBasis`
- if the surface is a sampled actuator basis with actuator geometry and command
  behavior, prefer `DeformableMirror`
- if the surface is just a named convenience constructor, keep it as a wrapper
  over an existing maintained basis type

## First Decision: Modal Basis Or Actuator Mirror?

### Modal / Low-Order Surface

Use `ModalControllableOptic` when the surface is naturally expressed as a small
set of basis modes.

Current maintained examples:

- `TipTiltMirror`
- `FocusStage`
- `ModalControllableOptic(..., ZernikeOpticBasis(...))`
- `ModalControllableOptic(..., MatrixModalBasis(...))`

This is the right surface when the optic is conceptually:

- a modal corrector
- a low-order steering stage
- a small basis of known shapes

### Actuator-Grid Mirror

Use `DeformableMirror` when the surface is naturally expressed as actuator
commands applied through influence functions.

Current maintained examples:

- Gaussian continuous-sheet DM models
- measured sampled influence-function models
- clipped or health-mapped actuator behavior

If the surface needs:

- actuator topology
- influence-function modeling
- actuator clipping or health logic

then it belongs in the DM layer, not in `ModalControllableOptic`.

## File Placement

Most modal controllable-optic extensions should be implemented in
[controllable_optics.jl](../src/Optics/controllable_optics.jl).

Typical steps:

1. add the new basis or convenience wrapper in
   [controllable_optics.jl](../src/Optics/controllable_optics.jl)
2. export public types or constructors from
   [AdaptiveOpticsSim.jl](../src/AdaptiveOpticsSim.jl) if the surface is meant
   to be public
3. add tests in
   [control_and_runtime.jl](../test/testsets/control_and_runtime.jl)
4. update user docs if the surface is public

If the new optic is really a DM variant, do not add it here. Use the DM guide
instead.

## Adding A New Modal Basis

### 1. Define The Basis Type

Add a concrete subtype of `AbstractModalOpticBasis`.

Current maintained basis types:

- `FunctionModalBasis`
- `MatrixModalBasis`
- `ZernikeOpticBasis`
- `CartesianTiltBasis`
- `QuadraticFocusBasis`

The basis type should own only the parameters needed to construct the sampled
mode matrix.

### 2. Provide Basis Metadata

At minimum, define:

- `_modal_basis_normalize(basis)`
- `_default_modal_labels(basis)`

These keep the family-local policy with the basis rather than in the generic
constructor path.

### 3. Provide The Basis Matrix Builder

Define:

- `_modal_basis_matrix(tel, basis, T, selector)`

That method should produce the host-space sampled basis matrix with the correct
resolution and mode count. The generic constructor will handle the normalized
backend copy and state assembly.

### 4. Add A Convenience Wrapper Only If It Buys Clarity

Named wrappers such as `TipTiltMirror` and `FocusStage` are appropriate when
they improve semantics.

Good wrapper:

- the name describes a meaningful optic family or hardware-facing surface

Bad wrapper:

- the name is just a synonym for an existing basis with no clearer meaning

## Command Layout Rule

Every controllable optic should expose a stable command surface through:

- `command_storage(optic)`
- `command_layout(optic)`

If a new optic has segmented or structured command semantics, encode them in
the `RuntimeCommandLayout` rather than inventing a one-off setter convention.

## Reuse Rule

This is the main architectural test:

- if a new optic is just a new basis, do not add a new optic family
- if a second optic wants the same basis machinery, keep that machinery in the
  shared modal path
- if the feature requires actuator topology and influence models, move it to
  the DM layer instead of stretching the modal API

## Validation Checklist

Every new controllable optic should add evidence in
[control_and_runtime.jl](../test/testsets/control_and_runtime.jl).

Minimum expectations:

1. construction works and validates bad inputs correctly
2. command length and labels are correct
3. `apply!` produces finite, correctly shaped OPD
4. `DMAdditive` and `DMReplace` semantics are correct when relevant
5. CPU and maintained GPU backends work when the surface is supported there

Add reference or external-truth validation if the optic becomes part of a
maintained scientific comparison surface.

## What Not To Do

- do not add a new optic family when a basis type is enough
- do not hide command segmentation outside `RuntimeCommandLayout`
- do not duplicate DM behavior in `ModalControllableOptic`
- do not add new modal families through `isa` branches in generic code

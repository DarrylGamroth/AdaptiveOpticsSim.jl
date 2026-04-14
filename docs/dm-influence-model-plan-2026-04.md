# DM Influence Model Plan 2026-04

Status: active

## Purpose

This plan refactors deformable-mirror influence handling from a pair of scalar
keywords into an explicit influence-model interface.

The goal is to keep the public `DeformableMirror` API clean and idiomatic while
making the internal model surface more extensible for:

- Gaussian influence parameterized by `influence_width`
- Gaussian influence parameterized by nearest-neighbor `mechanical_coupling`
- future dense/measured influence matrices
- future DM-family specialization without keyword proliferation

This is an influence-model refactor, not a change to the maintained DM physics
or runtime semantics.

## Why This Is Needed

The current DM constructor now accepts both:

- `influence_width=...`
- `mechanical_coupling=...`

That is useful, but it does not scale well if the DM family grows. Each new
influence representation would otherwise become another top-level constructor
keyword and another set of special cases inside the DM build path.

The package already uses small abstract model families in other subsystems:

- detector timing, nonlinearity, thermal, and defect models
- modal optic basis specs

The DM should follow the same pattern.

## Design Principles

1. Keep the public API semantic.
- Users should continue to write `DeformableMirror(tel; ...)`.
- The common Gaussian-width and Gaussian-coupling inputs should remain concise.

2. Use small concrete model values plus multiple dispatch.
- Avoid an OO-style hierarchy with mutable strategy objects.
- Keep the interface narrow and computation-oriented.

3. Preserve hot-path performance.
- Influence-model choice should be resolved at construction/build time.
- DM application kernels should remain allocation-free and statically
  specialized where possible.

4. Separate representation from canonical runtime storage.
- Different influence models may build the same stored actuator-to-OPD surface.
- The runtime DM state should continue to expose the same efficient command and
  apply paths.

5. Do not broaden the abstraction prematurely.
- This plan is for DM influence models, not every optical basis concept in the
  package.

## Target API Shape

### New interface

Add:

- `AbstractDMInfluenceModel`

First concrete models:

- `GaussianInfluenceWidth`
- `GaussianMechanicalCoupling`

Planned later when needed:

- `DenseInfluenceMatrix`
- `MeasuredInfluenceMap`

### Public constructor surface

The preferred public DM API remains:

```julia
DeformableMirror(tel; n_act=..., influence_width=...)
DeformableMirror(tel; n_act=..., mechanical_coupling=...)
DeformableMirror(tel; n_act=..., influence_model=GaussianInfluenceWidth(...))
DeformableMirror(tel; n_act=..., influence_model=GaussianMechanicalCoupling(...))
```

Public behavior goals:

- `influence_width=...` remains the default/common path
- `mechanical_coupling=...` remains a supported convenience path
- `influence_model=...` becomes the explicit advanced path
- specifying multiple influence parameterizations at once should keep throwing a
  structured configuration error

### Internal parameter storage

Refactor `DeformableMirrorParams` so it stores an influence-model object rather
than a raw `influence_width` scalar.

Conceptually:

```julia
struct DeformableMirrorParams{T<:AbstractFloat,I<:AbstractDMInfluenceModel}
    n_act::Int
    influence_model::I
    misregistration::Misregistration{T}
end
```

Then expose semantic accessors:

- `influence_model(dm)`
- `influence_width(dm)` where meaningful
- `mechanical_coupling(dm)` where meaningful

The existing `mechanical_coupling(dm)` helper should become a real model
accessor, not just a derived convenience over one scalar field.

## Required Interface Functions

The DM influence-model interface should stay narrow.

Required first:

- `canonical_influence_width(model, n_act)`
- `supports_separable_influence(model, misregistration)`
- `build_dense_influence!(dm, tel, model)`
- `build_separable_influence!(dm, tel, model)` where applicable

Useful derived helpers:

- `mechanical_coupling(model, n_act)` where meaningful
- `influence_model_kind(model)`

Rules:

- Gaussian models should share the same canonical dense/separable build logic
- models without a separable formulation must cleanly opt out
- the runtime `apply!` path should not branch on loosely typed conditionals in
  hot loops

## Phase Plan

### DMI-1 Introduce the DM influence-model interface

Deliverables:

- add `AbstractDMInfluenceModel`
- add `GaussianInfluenceWidth`
- add `GaussianMechanicalCoupling`
- add conversion/validation helpers between width and coupling

Success criteria:

- no behavior change for the maintained DM surface
- the scalar conversion helpers continue to work

### DMI-2 Store the model in `DeformableMirrorParams`

Deliverables:

- `DeformableMirrorParams` stores `influence_model`
- add public accessors:
  - `influence_model(dm)`
  - `influence_width(dm)` where meaningful
  - `mechanical_coupling(dm)` where meaningful
- keep `params.influence_width` out of the active public story

Success criteria:

- no user-visible regression in existing constructors
- existing runtime and calibration code compiles through the new params shape

### DMI-3 Route constructors through semantic model objects

Deliverables:

- `influence_width=...` builds `GaussianInfluenceWidth`
- `mechanical_coupling=...` builds `GaussianMechanicalCoupling`
- add explicit `influence_model=...` constructor support
- retain structured errors for conflicting inputs

Success criteria:

- all currently documented constructors still work
- advanced users can now pass a model object directly

### DMI-4 Dispatch dense and separable influence generation on the model

Deliverables:

- `build_influence_functions!` routes through the influence model
- separable fast path is enabled through model capability, not by ad hoc
  scalar-field assumptions
- no runtime hot-path regressions

Success criteria:

- CPU, CUDA, and AMDGPU maintained DM surfaces remain green
- separable Gaussian path still uses the fast specialized implementation

### DMI-5 Add at least one non-Gaussian explicit model surface

Deliverables:

- add `DenseInfluenceMatrix` or equivalent explicit influence-model object
- allow construction from a precomputed mode matrix

Success criteria:

- `DeformableMirror(...; influence_model=DenseInfluenceMatrix(...))` works
- the model participates in runtime/control flows
- unsupported fast paths fail structurally rather than silently

### DMI-6 Update docs, examples, and validation surfaces

Deliverables:

- update API reference
- update user guide and cookbook where DM parameterization is discussed
- add targeted tests for:
  - width/coupling equivalence
  - direct `influence_model=...` construction
  - dense explicit model path
  - separable-capability behavior
- rerun maintained GPU validation and OOPAO reference suites

Success criteria:

- stable docs present the clean semantic surface
- validation evidence remains green on maintained CPU/CUDA/AMDGPU paths

## Testing And Validation

Required evidence for completion:

1. Focused runtime/control tests pass.
- `Pkg.test(test_args=["control_and_runtime"])`

2. Reference surfaces stay green.
- `Pkg.test(test_args=["reference_and_tutorials"])`

3. Optional backend parity stays green.
- AMDGPU runtime-equivalence script
- CUDA runtime-equivalence script

4. New explicit-model tests exist for:
- Gaussian width vs Gaussian coupling equivalence
- direct model-object construction
- dense explicit influence-model construction
- structured failure on unsupported combinations

## Non-Goals

- changing the maintained Gaussian DM physics
- changing RTC/runtime command semantics
- introducing a broad package-wide `AbstractBasis` unification
- adding measured influence I/O/file-format support in the core package

## Completion Criteria

This plan is complete when:

- `AbstractDMInfluenceModel` is the internal DM influence abstraction
- the public `DeformableMirror` API remains concise and semantic
- the common Gaussian width/coupling paths are preserved
- at least one explicit non-Gaussian influence-model path exists
- maintained CPU/CUDA/AMDGPU and reference validation stay green

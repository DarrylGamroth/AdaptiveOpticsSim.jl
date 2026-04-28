# Adding A Wavefront Sensor

Status: active

## Purpose

This guide is the maintained authoring path for adding a new wavefront-sensor
family to `AdaptiveOpticsSim.jl`.

Use it when you need to:

- add a new `AbstractWFS` family
- add a new detector-coupled or detector-free sensing model
- decide whether a behavior belongs in a shared WFS seam or one sensor family

Use together with:

- [api-reference.md](./api-reference.md)
- [model-cookbook.md](./model-cookbook.md)
- [interface-contract-closeout-plan-2026-04.md](./interface-contract-closeout-plan-2026-04.md)

## Design Rule

The WFS subsystem is organized around this split:

- WFS families own their sensing physics and state layout
- shared runtime and calibration layers consume a small maintained accessor
  contract instead of reaching directly into family fields

In practice:

- if a quantity is part of the maintained WFS contract, expose it through an
  accessor in [interface.jl](../src/wfs/interface.jl)
- if a behavior is reused across multiple WFS families, factor it into a shared
  seam instead of duplicating it
- do not make runtime or calibration code depend on undocumented family-local
  state fields

## First Decision: What Kind Of WFS Is It?

Start by deciding how the sensor fundamentally produces its exported signal.

Common maintained patterns:

- subaperture/centroid style
  - [shack_hartmann](../src/wfs/shack_hartmann/)
- focal-plane modulation and differential pupil-plane style
  - [pyramid](../src/wfs/pyramid/)
  - [bioedge](../src/wfs/bioedge/)
- modal or normalized-signal style
  - [zernike.jl](../src/wfs/zernike.jl)
  - [curvature.jl](../src/wfs/curvature.jl)

The maintained subsystem contract still uses `slopes(wfs)` as the generic
exported 1-D signal vector even when the values are not geometric slopes.

## File Placement

New WFS families should normally add one new file or one new family directory
under [`src/wfs/`](../src/wfs/).

Typical steps:

1. add the WFS family implementation under `src/wfs/`
2. include it from the relevant WFS entry file
3. export public types and constructors from
   [AdaptiveOpticsSim.jl](../src/AdaptiveOpticsSim.jl) if the family is meant
   to be public
4. extend the maintained WFS accessors in
   [interface.jl](../src/wfs/interface.jl) when the family exports a stable
   interface quantity
5. add tests in the family testset plus
   [control_and_runtime.jl](../test/testsets/control_and_runtime.jl) when the
   family participates in maintained runtime flows

## Minimum WFS Shape

Every WFS family should have:

- a params struct for immutable configuration
- a mutable state struct for runtime buffers
- explicit mutating measurement functions for the hot path

The hot-path contract should stay explicit:

- build once
- reuse buffers
- mutate state through `measure!`, `sense!`, or other family-specific `!`
  functions

## Maintained Accessors

These are the maintained subsystem-level accessors in
[interface.jl](../src/wfs/interface.jl):

- `slopes(wfs)`
- `valid_subaperture_mask(wfs)`
- `reference_signal(wfs)`
- `camera_frame(wfs)`

And the corresponding capability helpers:

- `supports_valid_subaperture_mask(wfs)`
- `supports_reference_signal(wfs)`
- `supports_camera_frame(wfs)`

The rule is simple:

- if the quantity is a stable part of the maintained family contract, add the
  accessor
- if it is just a family-local implementation detail, keep it local

Do not extend the interface accessors for ephemeral internal buffers.

## Detector-Coupled Vs Detector-Free Sensors

If the WFS uses a detector-like intermediate image:

- expose a stable `camera_frame(wfs)` only if that frame is part of the
  maintained external contract
- keep detector assembly and capture sequencing explicit

If the WFS is detector-free or only uses an internal intermediate surface:

- do not invent a camera-frame accessor just for symmetry

The accessors are meant to describe stable exported state, not every internal
buffer.

## Prepared Runtime And Compatibility

When adding a new WFS family, decide explicitly whether it supports:

- prepared runtime construction
- grouped/multi-source runtime paths
- detector coupling
- valid-mask export
- reference-signal export

If the answer is no for a maintained workflow, reject it structurally rather
than letting runtime code fail later.

## Reuse Rule

This is the main architectural test for WFS work:

- if one family owns a behavior, it may stay family-local
- if a second family needs the same behavior, factor it into a shared seam
  unless there is a strong reason not to

Good reuse targets:

- interface accessors
- shared detector coupling seams
- shared valid-mask or reference-signal conventions
- shared prepared-runtime helpers

Bad pattern:

- copying another WFS family’s runtime/export logic and then editing field
  names

## Validation Checklist

Every new WFS family should add evidence in the family-specific testset and in
maintained runtime tests when the family is part of the supported runtime
surface.

Minimum expectations:

1. construction works and validates bad inputs correctly
2. `slopes(wfs)` returns a stable, correctly shaped signal vector
3. `valid_subaperture_mask`, `reference_signal`, and `camera_frame` only claim
   support when the family really exports them
4. runtime integration works through `sense!` and maintained scenario builders
   when supported
5. backend behavior works on CPU and maintained GPU backends when supported

Add reference or external-truth validation if the family is part of the
maintained scientific validation surface.

## What Not To Do

- do not make runtime code depend on undocumented WFS state fields
- do not expose every internal buffer as part of the WFS interface
- do not add new WFS family behavior through `isa` branches in shared runtime
  code
- do not bypass the maintained accessors when the quantity is part of the
  supported WFS contract

# Adding A Detector

Status: active

## Purpose

This guide is the maintained authoring path for adding a new detector family to
 `AdaptiveOpticsSim.jl`.

Use it when you need to:

- add a new detector family
- add a reusable detector feature that should apply to multiple families
- decide whether a behavior belongs in a sensor file or a shared detector layer

Use together with:

- [api-reference.md](./api-reference.md)
- [detector-interface-extension.md](./detector-interface-extension.md)
- [detector-validation.md](./detector-validation.md)

## Design Rule

The detector subsystem is organized around this split:

- shared detector layers own reusable models and orchestration
- detector families own family-specific physics and policy

In practice:

- if a behavior is generic across multiple sensors, make it a reusable model
  family or shared detector sublayer
- if a behavior is intrinsic to one detector family, keep it in that family file
- do not add new detector behavior as ad hoc branches in `generic.jl`

## First Decision: Frame Or Counting

Start by choosing which detector category the new family belongs to.

### Frame Detector

Use a `FrameSensorType` when the detector fundamentally produces a frame image.

Current examples:

- [ccd.jl](../src/Detectors/ccd.jl)
- [cmos.jl](../src/Detectors/cmos.jl)
- [ingaas.jl](../src/Detectors/ingaas.jl)
- [hgcdte_avalanche_array.jl](../src/Detectors/hgcdte_avalanche_array.jl)

Frame detectors usually reuse the generic `Detector(...)` constructor and the
shared frame-capture pipeline in [frame_capture.jl](../src/Detectors/frame_capture.jl).

### Counting Detector

Use a `CountingSensorType` when the detector fundamentally produces counted
channels or counts rather than a linear frame sensor pipeline.

Current examples:

- [apd.jl](../src/Detectors/apd.jl)
- [spad_array.jl](../src/Detectors/spad_array.jl)

Counting detectors currently own more of their capture path directly, because
their physics and export shape differ more substantially from the frame
pipeline. But they now share a maintained counting-detector sublayer in
[counting_common.jl](../src/Detectors/counting_common.jl) for common counting
capture semantics and metadata assembly.

## File Placement

A new detector family should normally add one new file under
[`src/Detectors/`](../src/Detectors/) and then be included from
[detector.jl](../src/Detectors/detector.jl).

Typical steps:

1. add the detector family file under `src/Detectors/`
2. include it from [detector.jl](../src/Detectors/detector.jl)
3. export public types and constructors from [AdaptiveOpticsSim.jl](../src/AdaptiveOpticsSim.jl)
4. add tests in [test/testsets/detectors.jl](../test/testsets/detectors.jl)
5. update user docs if the detector is intended to be public

## Adding A Frame Detector

### 1. Define the Sensor Type

Add a concrete sensor type that subtypes `FrameSensorType`.

Example shape:

```julia
struct MySensor{T<:AbstractFloat} <: FrameSensorType
    parameter_a::T
    parameter_b::T
end
```

Keep this type focused on the detector-family parameters, not runtime frame
 state.

### 2. Provide Family Identity And Capabilities

At minimum, define:

- `detector_sensor_symbol(sensor)`

Then add capability predicates only when they are true for the new family, for
 example:

- `supports_detector_response`
- `supports_readout_correction`
- `supports_nondestructive_reads`
- `supports_read_cube`
- `supports_detector_persistence`
- `supports_detector_nonlinearity`
- `supports_shutter_timing`

These methods live with the sensor family, not in
[generic.jl](../src/Detectors/generic.jl).

### 3. Provide Family Defaults

If the family has a natural default reusable model, define it locally:

- `default_response_model(sensor; T=..., backend=...)`
- and, only if needed later, analogous default helpers for other model classes

Do not hard-code those defaults in central detector assembly.

### 4. Reuse Shared Models Before Adding New Ones

Check whether the new family can already use existing shared detector model
 seams:

- frame response / MTF
- defect models
- timing models
- readout correction
- nonlinearity
- thermal models
- frame sampling / multi-read products

If the feature is already represented there, use it. Do not clone the logic in
 the family file.

### 5. Add Family-Specific Physics Hooks

Most frame detectors customize the generic capture pipeline by extending a small
 set of hooks in [frame_capture.jl](../src/Detectors/frame_capture.jl):

- `apply_sensor_statistics!(sensor, det, rng)`
- `apply_pre_readout_gain!(sensor, det, rng)`
- `apply_post_readout_gain!(sensor, det)`

Only add a detector-specific `finalize_capture!` when the family genuinely
 needs a different sequencing contract.

Current example:

- [hgcdte_avalanche_array.jl](../src/Detectors/hgcdte_avalanche_array.jl)

That family keeps avalanche-specific gain/glow/saturation behavior local while
 reusing the shared multi-read layer.

### 6. Use Readout Products Only If Needed

If the detector needs richer outputs than the main frame, use the
 `FrameReadoutProducts` seam.

Options:

- reuse `SampledFrameReadoutProducts`
- reuse `MultiReadFrameReadoutProducts`
- add a new `FrameReadoutProducts` subtype only if the existing payloads do not
  fit

If you add a new payload type, implement only the accessors you actually
 provide:

- `detector_reference_frame`
- `detector_signal_frame`
- `detector_combined_frame`
- `detector_reference_cube`
- `detector_signal_cube`
- `detector_read_cube`
- `detector_read_times`

These accessors default to `nothing`, so a new payload should not require edits
 to central generic code just to expose one extra product.

## Adding A Counting Detector

Counting detectors still expose more family-level policy than frame detectors,
but they no longer need to duplicate the full capture skeleton.

Use [apd.jl](../src/Detectors/apd.jl) and
[spad_array.jl](../src/Detectors/spad_array.jl) as the reference shapes.

Typical steps:

1. define a sensor subtype of `CountingSensorType`
2. define the detector type and its params/state
3. add the small counting accessors needed by the shared counting layer
4. implement only the family-specific counting hooks that remain
5. reuse the shared counting export metadata path unless the family genuinely
   needs a different output contract

Shared counting seams include:

- noise models
- dead-time models
- gate models
- correlation models such as afterpulsing and channel crosstalk
- shared counting capture sequencing in
  [counting_common.jl](../src/Detectors/counting_common.jl)
- shared counting export metadata assembly

If a counting behavior is reusable across multiple families, it should become a
shared counting model or shared counting-layer hook rather than being embedded
in one detector file.

## Reusable Feature Rule

This is the main architectural test when adding new physical behavior:

- if one detector needs it, it may remain family-local
- if a second detector needs it, factor it into a shared seam unless there is a
  strong reason not to

Current examples of good reuse:

- readout correction models
- thermal models
- response models
- counting-correlation models
- multi-read frame-readout assembly

Bad pattern:

- copying a family-local algorithm into another detector file with only small
  parameter changes

## Export Metadata

Every detector family should expose enough metadata to qualify generated output.

For frame detectors, much of this is already assembled in
[generic.jl](../src/Detectors/generic.jl) through `detector_export_metadata(det)`.

For counting detectors, or for families with materially different export
 surfaces, define a detector-specific `detector_export_metadata(...)`.

The rule is simple:

- output arrays alone are not enough
- the metadata should record the configuration that materially explains those
  outputs

## Validation Checklist

Minimum validation for a new detector family:

1. constructor and parameter validation tests
2. family-specific physics tests
3. export metadata tests
4. readout-product tests if applicable
5. backend tests if the maintained backend surface should support that family

Use these existing surfaces:

- [test/testsets/detectors.jl](../test/testsets/detectors.jl)
- [detector-validation.md](./detector-validation.md)
- [scripts/generate_detector_validation_artifact.jl](../scripts/generate_detector_validation_artifact.jl)
- [test/runtests_amdgpu.jl](../test/runtests_amdgpu.jl)
- [test/runtests_cuda.jl](../test/runtests_cuda.jl)

If the detector is intended to be part of the maintained model-validity story,
 add a committed detector artifact case rather than relying only on unit tests.

## When To Add A Shared Layer

Create or extend a shared detector layer when all of the following are true:

- the behavior is not specific to one detector family
- the behavior has clear semantics independent of one family name
- a second detector can already use it or is likely to use it soon
- factoring it out reduces duplicated policy or duplicated runtime/export logic

Recent example:

- the reusable multi-read layer in
  [frame_sampling.jl](../src/Detectors/frame_sampling.jl)

That was the right refactor because explicit multi-read readout products and
 read-time assembly were not really HgCdTe-specific concepts.

## What Not To Do

- do not add detector-family feature tables to `generic.jl`
- do not rely on direct mutation of internal detector state as an extension API
- do not use `isa` branches in package code when a dispatch seam is the real
  abstraction
- do not duplicate reusable physical models in multiple sensor files
- do not make a new detector family carry runtime buffers in its immutable
  sensor params type

## Practical Recipe

If you are adding a new detector family, the shortest correct workflow is:

1. decide frame vs counting
2. define the sensor type
3. add `detector_sensor_symbol` and true capability methods
4. reuse existing shared model families first
5. add only the family-specific physics hooks that remain
6. add readout products only if the detector needs richer export surfaces
7. add tests and, if maintained, a detector artifact case
8. if a second family needs the same behavior, factor it into a shared layer

If that workflow feels awkward for a candidate detector, that is the signal
 that the interface likely still needs work.

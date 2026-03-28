# Detector Response Interface Plan

Status: implemented in core detector response APIs and maintained models,
including a sampled-kernel response path.

## Goal

Define a reusable detector-response interface that cleanly separates:

- sensor optics
- readout/export layout
- detector response physics

The immediate focus is frame-detector spatial response, including pixel
aperture, fill factor, pitch, and MTF-like response models, while keeping the
core design open to non-frame detector families.

This plan is intentionally about interface and layering, not just one concrete
response model.

## Design Rules

- Use abstract types, traits, and multiple dispatch.
- Avoid `isa`-driven branching in hot paths.
- Keep AO188 and other instrument defaults out of core.
- Provide null/identity models explicitly.
- Keep counting-detector behavior out of frame-response abstractions unless a
  capability is genuinely shared.
- Allow both convolution-space and Fourier-space implementations under the same
  interface when they are mathematically the same detector response.

## Scope

This plan covers:

- detector-response interfaces for frame detectors
- detector-response metadata/export
- first physically named pixel-geometry models
- MTF-capability traits
- integration with the existing `Detector` pipeline

This plan does not yet cover:

- APD/counting dead-time behavior
- SPAD/SPAD-array event logic
- AO188-specific detector tuning
- measured instrument-specific calibration files

## Interface Shape

### Core abstractions

Add or formalize:

- `AbstractDetectorResponse`
- `AbstractFrameResponse <: AbstractDetectorResponse`
- `AbstractFrameMTF <: AbstractFrameResponse`

Current `FrameResponseModel` should either become the public abstract family or
be replaced by the hierarchy above. The important requirement is that the
abstraction names describe response physics rather than one implementation.

### Required operations

Frame-response models should support a common interface such as:

- `response_metadata(model)`
- `supports_detector_mtf(model)`
- `response_support(model, dims...)`
- `apply_response!(model, det)`

Optional traits/capabilities:

- `is_shift_invariant(model)`
- `supports_frequency_domain_application(model)`
- `supports_separable_application(model)`
- `supports_subpixel_geometry(model)`

These should be trait queries, not booleans stored in one large response type.

### Null model

Keep an explicit identity model:

- `NullFrameResponse <: AbstractFrameResponse`

This is required for:

- deterministic regression surfaces
- fast paths
- simpler tutorial/default configurations

## First concrete response family

The first physically meaningful maintained family should include:

- `GaussianPixelResponse`
  - effective blur-like response in pixel units
  - successor to the current `SeparableGaussianPixelResponse`
- `RectangularPixelAperture`
  - explicit aperture width in pitch units
  - models fill factor directly
- `SeparablePixelMTF`
  - frequency-domain separable response for cases where MTF is the natural
    specification

The key improvement is that fill factor and aperture shape become first-class
physical parameters rather than being hidden inside one blur width.

## Pixel geometry model

Detector response should be able to represent:

- `pitch_x`, `pitch_y`
- `fill_factor_x`, `fill_factor_y`
- aperture shape
  - initially rectangular
  - possibly later circular or measured masks

Pitch should not duplicate optical-plane pixel scale. Instead:

- optical models continue to determine image sampling in physical/image-plane
  units
- detector response interprets that sampled image through detector-pixel
  geometry

That means pitch/fill factor belong to the detector-response layer, not the WFS
optics layer.

## Metadata

Detector export metadata should expose response information clearly enough for
HIL and downstream consumers:

- response family symbol
- geometric parameters
- whether response was applied in image space or frequency space
- whether the response is separable

Do not hide fill factor or aperture geometry inside one generic width field
once more physical models exist.

## Integration Plan

### Phase 1: Interface cleanup

- Introduce the detector-response abstract family.
- Rename `SeparableGaussianPixelResponse` if needed so the public name reflects
  its role in the new family.
- Keep `NullFrameResponse` as the explicit identity path.
- Update docs and API references.

### Phase 2: Pixel-geometry models

- Implement `RectangularPixelAperture`.
- Add explicit fill-factor parameters.
- Support both CPU and accelerator execution.
- Extend metadata and tests accordingly.

### Phase 3: MTF path

- Add `AbstractFrameMTF` and the first maintained MTF-backed model.
- Allow response application either:
  - by direct convolution in image space, or
  - by transfer multiplication in frequency space
- Keep the same high-level `apply_response!` interface.

### Phase 4: Pipeline integration

- Integrate the new response family cleanly into `Detector`.
- Preserve zero-allocation hot paths where possible.
- Keep batched capture support explicit:
  - supported for null and maintained separable paths
  - rejected clearly for unsupported response models

### Phase 5: Broader detector-family adoption

- Decide which response models are valid for:
  - CCD
  - CMOS
  - EMCCD
  - InGaAs
  - HgCdTe avalanche arrays
- Keep this as capability-based policy, not one universal detector flag set.

## Testing Requirements

Add or maintain tests for:

- null-model identity behavior
- response metadata correctness
- fill-factor and aperture-shape parameter validation
- CPU/GPU numerical agreement for maintained response models
- separable-vs-direct equivalence where both exist
- batched-capture acceptance/rejection behavior

Representative tests should include:

- a response model that changes encircled energy but preserves total flux
- a fill-factor change that measurably alters detector sampling
- a detector-window path after response application

## Naming Guidance

- Prefer names that describe detector physics, not implementation details.
- Avoid detector-family names in function names when dispatch already carries
  the specialization.
- Prefer nouns for models:
  - `RectangularPixelAperture`
  - `GaussianPixelResponse`
  - `SeparablePixelMTF`
- Prefer generic verbs for operations:
  - `apply_response!`
  - `response_metadata`
  - `supports_detector_mtf`

## Open Questions

- Should `FrameResponseModel` be retained as the main public abstract type or
  replaced by `AbstractFrameResponse`?
- Should measured detector responses be introduced as sampled kernels first, or
  only after the analytic families are stable?
- Should detector response eventually be allowed to vary per sensor axis or
  subarray region, or should that stay out of the first interface?

## Recommended next step

Implement Phase 1 and Phase 2 together:

- formalize the detector-response abstract family
- add `RectangularPixelAperture`
- extend metadata to include explicit geometry

That is the smallest useful step that materially improves fidelity and keeps
the interface clean.

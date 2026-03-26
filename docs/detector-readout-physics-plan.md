# Detector And Readout Physics Plan

This plan defines the next detector/readout-physics pass for
`AdaptiveOpticsSim.jl`.

The immediate driver is AO188 curvature sensing, but the goal is broader:

- keep sensor optics reusable and detector-agnostic where possible
- model frame detectors and counting detectors as distinct families
- add higher-fidelity readout physics without collapsing everything into one
  generic `Detector`
- make the detector/readout layer scalable to future sensors such as APD arrays
  and SPAD imaging arrays
- define interfaces with types, traits, and multiple dispatch rather than
  `isa`-driven branching or flag-heavy runtime logic

## Current State

The package now has a cleaner split than before:

- sensor optics live in sensor modules such as
  `src/WFS/shack_hartmann.jl`,
  `src/WFS/pyramid.jl`,
  `src/WFS/bioedge.jl`,
  `src/WFS/zernike.jl`,
  `src/WFS/curvature.jl`
- the generic frame detector lives in `src/Optics/detector.jl`
- detector families are now explicit:
  - `FrameSensorType`
  - `CountingSensorType`
  - `APDSensor` no longer pretends to be a generic frame detector
- `CurvatureWFS` now has:
  - frame readout
  - counting readout
  - branch throughput/background response
  - configurable detector-plane crop/sampling geometry

This is a better foundation, but it is still not a full detector/readout model.

## Main Questions

## Interface Rules

This plan should be implemented with the same interface discipline used
elsewhere in the package:

- use abstract types to define the main detector/readout families
- use traits for optional capabilities
- use multiple dispatch for behavior selection
- avoid `isa` checks in hot paths when dispatch can express the interface
- avoid a single generic detector object with many booleans that emulate type
  differences

Concretely, detector/readout work in this plan should prefer shapes like:

- `AbstractFrameDetector`
- `AbstractCountingDetector`
- `AbstractFrameReadout`
- `AbstractCountingReadout`
- trait queries such as `supports_detector_mtf(x)` or
  `supports_counting_statistics(x)`
- dispatch like `capture!(det::AbstractFrameDetector, frame; rng)` and
  `capture!(det::AbstractCountingDetector, channels; rng)`

The goal is not abstraction for its own sake. The goal is to keep the physics
layer explicit and extensible without hidden concrete-type checks.

### 1. Should detector/readout physics be overhauled per detector family?

Yes, but not by duplicating everything per sensor.

The right structure is:

- shared optical sensor models
- shared readout reduction models where the math is genuinely common
- detector-family-specific response/noise/statistics layers

That means:

- do not create one giant generic detector path that tries to model CCD,
  CMOS, EMCCD, APD, and SPAD through a single set of flags
- do not hardwire detector behavior into each WFS either

Instead:

- frame detectors should share a reusable frame-detector pipeline
- counting detectors should share a reusable counting-detector pipeline
- each sensor should choose which family its readout uses

### 2. Should we implement an MTF path for detector layout?

Yes, for frame detectors.

A detector MTF / pixel response path is physically meaningful for:

- CCD
- CMOS
- EMCCD
- frame-based curvature readout
- Shack-Hartmann/Pyramid/BioEdge/Zernike frame exports

It is not the first priority for counting detectors such as scalar APD-channel
readout, because those are not naturally modeled as a pixel-array convolution
problem.

So:

- add MTF/pixel response to the frame-detector family
- do not force the same model onto APD-like counting readout

### 3. Should we consider SPAD and SPAD imaging arrays?

Yes.

Not as an immediate AO188 blocker, but absolutely as part of the architecture.

SPADs split into two useful categories:

- channel-style SPAD counters
  - conceptually close to APD counting paths
- SPAD imaging arrays
  - conceptually closer to frame detectors, but with counting statistics,
    dead time, fill factor, and sometimes binary/event-style readout behavior

That means SPAD should not be squeezed into either "just APD" or "just CMOS".
It deserves a future branch in the detector/readout taxonomy.

## Design Direction

The detector/readout stack should be organized into three layers.

### Layer 1: Optical Sensor

This layer produces optical intensity or channelized optical readout before
detector-specific electronics/statistics.

Examples:

- `ShackHartmann`
- `PyramidWFS`
- `BioEdgeWFS`
- `ZernikeWFS`
- `CurvatureWFS`

This layer should handle:

- propagation
- masks/modulation
- branch geometry
- detector-plane sampling/crop
- channel packing or frame packing

This layer should not decide:

- EM gain
- read noise
- APD avalanche statistics
- SPAD dead time
- ADC quantization

unless the quantity is inseparable from the sensor optics itself.

### Layer 2: Readout Model

This layer maps optical outputs into detector-family-compatible measurement
channels.

Examples:

- frame readout
- counting readout
- quadrant readout
- per-subap scalar reduction
- APD-channel packing

This layer should stay generic where possible.

### Layer 3: Detector Response

This layer handles the detector/electronics/statistics model.

Frame-detector examples:

- QE
- dark current
- read noise
- EM gain / excess noise
- saturation
- ADC quantization
- pixel response / MTF

Counting-detector examples:

- QE / throughput
- dark counts
- avalanche gain or effective counting response
- dead time
- afterpulsing if ever needed
- channel non-uniformity

## Immediate APD-Like Counting Plan

This is the next practical plan for AO188-like curvature work.

### Phase 1: Counting Readout Contract

Define a small reusable counting-readout surface.

Add core types along these lines:

- `AbstractCountingReadout`
- `CountingReadoutMetadata`
- `channel_output(readout)`

Optional capability traits should be used where behavior is not universal, for
example:

- `supports_frame_export(readout)`
- `supports_channel_metadata(readout)`

The goal is to avoid treating counting output as a fake image.

Deliverables:

- explicit channel-output conventions
- metadata describing channel layout
- runtime/export support analogous to frame metadata

### Phase 2: Counting Detector Family

Add a reusable counting-detector model rather than embedding APD behavior in
`CurvatureWFS`.

Candidate type family:

- `AbstractCountingDetector`
- `APDDetector`

Candidate trait surface:

- `supports_counting_noise(det)`
- `supports_dead_time(det)`
- `supports_channel_gain_map(det)`

Candidate physics:

- QE
- dark count rate
- per-channel gain or sensitivity non-uniformity
- optional Poisson counting noise
- optional additive count-floor / bias model

This phase should remain simple enough for deterministic regression and HIL use.

### Phase 3: Curvature Counting Integration

Wire `CurvatureCountingReadout` to the counting-detector family.

That means:

- `CurvatureWFS` still owns the optics and channel packing
- the counting detector owns the channel response/statistics
- AO188 curvature can use that path without pretending an APD array is a frame
  detector

Deliverables:

- detector-coupled counting measurement path
- deterministic counting mode
- noisy counting mode
- runtime/readout metadata for HIL export

### Phase 4: AO188 Defaults

Once the counting-detector family exists, choose maintained AO188 defaults for:

- counting detector type
- QE
- dark count level
- deterministic vs noisy default profile
- any per-branch imbalance defaults if justified

This phase should stay example-level, not core-level.

## Frame-Detector Physics Plan

Frame detectors deserve a separate improvement track.

### Phase 1: Detector MTF / Pixel Response

Add an optional frame-detector MTF or pixel-response stage in
`src/Optics/detector.jl`.

Requirements:

- opt-in
- precomputable
- backend-friendly
- works with current frame sensors
- selected through detector type/trait dispatch, not detector-kind flags

The first maintained version should be simple:

- separable pixel-response blur or transfer function
- parameterized by fill factor / response width

This should be reusable for:

- SH
- Pyramid
- BioEdge
- Zernike
- frame-style CurvatureWFS

### Phase 2: Sensor-Specific Refinements

Then refine the frame families where needed:

- CCD
- CMOS
- EMCCD

That should be represented with concrete detector types and dispatch, not
deeply nested conditional logic inside one frame-detector implementation.

Examples:

- EM gain / excess-noise handling
- rolling/global-shutter style timing only if it matters later
- sensor-specific noise parameterization

This should not duplicate the optical sensor models.

## SPAD Plan

SPAD support should be planned now, even if implemented later.

### SPAD Channel Path

For small channelized readout, SPAD can reuse much of the counting-detector
architecture.

### SPAD Imaging Array Path

For array readout, SPAD should likely be a distinct detector family:

- `SPADArraySensor` or `SPADImager`

with its own type/trait-dispatch path rather than as a mode flag on `APD` or
`CMOS`.

Likely physics:

- fill factor
- photon-counting statistics
- dead time
- saturation / pile-up behavior
- event-frame or accumulated-count frame modes

This is not the same thing as either:

- `APDDetector`
- `CMOS` with Poisson noise

So the architecture should leave room for it now.

## Recommended Implementation Order

1. counting-readout metadata and interface
2. reusable `APDDetector` / counting-detector path
3. curvature counting integration
4. frame-detector MTF / pixel response
5. AO188 curvature detector default tuning
6. SPAD family design/implementation

## Non-Goals For This Pass

- do not try to model every detector phenomenon immediately
- do not fold SPAD into APD just because both count photons
- do not hardcode AO188 detector behavior into `CurvatureWFS`
- do not require detector MTF for counting-channel readout

## Success Criteria

- APD-like curvature readout no longer depends on the generic frame detector
- frame detectors have an explicit optional MTF/pixel-response path
- detector/readout docs reflect the detector-family split clearly
- AO188 curvature defaults become more instrument-like without making the core
  sensor Subaru-specific
- the design leaves a clean path for SPAD imaging arrays later

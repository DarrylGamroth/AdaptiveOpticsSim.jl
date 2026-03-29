# Detector Family Feature Plan

This plan defines the next detector-specific feature work for the detector
families that are most useful to current modeling work:

- `CMOSSensor`
- `EMCCDSensor`
- `InGaAsSensor`
- `APDDetector`

`HgCdTeAvalancheArraySensor` is intentionally not the focus here. Its current
surface is already broad enough that the next work should be driven by concrete
system needs rather than more speculative detector expansion.

## Goal

Expand the detector families that matter most for practical simulation while
preserving the detector architecture already established in core:

- optics, readout, and detector response remain separate layers
- detector families are expressed with types, traits, and multiple dispatch
- null/default models remain explicit
- AO188, AO3k, and other instrument examples remain example-layer
  configurations, not core API drivers

## Scope

This plan covers:

- detector-family-specific response/statistics/readout refinements
- metadata and interface additions needed to expose those refinements cleanly
- benchmarking and allocation checks on maintained system surfaces
- CPU/GPU compatibility requirements for maintained hot paths

This plan does not cover:

- further `HgCdTeAvalancheArraySensor` expansion unless a concrete gap appears
- SPAD / SPAD-array implementation
- WFS-specific detector fusion work except where profiling proves it is needed

## Design Rules

- Use detector-family types, not detector-name-prefixed functions, to select
  behavior.
- Use traits only for optional capabilities that may be shared across
  families.
- Avoid `isa`-driven branches where dispatch can express the interface.
- Keep family-specific physics in the family file rather than in generic
  detector plumbing.
- Preserve deterministic null paths for regression and HIL baselines.
- Keep hot paths allocation-aware from the start.

## Current State

The current detector families already provide a useful baseline:

- `CMOSSensor`
  - column readout noise
  - mild default detector response
- `EMCCDSensor`
  - excess-noise-factor approximation
- `InGaAsSensor`
  - glow-rate model
  - mild default detector response
- `APDDetector`
  - counting QE / gain / dark counts
  - non-paralyzable dead time
  - optional counting noise
  - optional channel gain map

This is a good foundation, but each family still has clear detector-specific
features worth adding.

## Family Priorities

### CMOS

Highest-value missing features:

- `PixelResponseNonuniformity`
  - static multiplicative PRNU map
- `DarkSignalNonuniformity`
  - static additive DSNU map
- grouped readout structure
  - per-column or per-output gain/offset maps
- shutter/read timing
  - explicit `GlobalShutter`
  - explicit `RollingShutter`

Why these matter:

- CMOS systems often show fixed-pattern structure that is not captured by one
  scalar column noise term
- readout topology can affect realism for frame sequencing and HIL

Recommended order:

1. PRNU / DSNU maps
2. grouped output gain/offset structure
3. shutter/read timing semantics

### EMCCD

Highest-value missing features:

- multiplication-register stochastic model
  - explicit gain-register sampling rather than only an effective excess-noise
    factor
- CIC refinement
  - clock-induced charge as a detector-family effect rather than only a CCD
    concern
- saturation / clipping semantics after EM gain
- optional gain-register aging / nonuniformity later if needed

Why these matter:

- EMCCD behavior is defined by the multiplication chain more than by generic
  frame read noise
- the current excess-noise-factor path is useful, but still a coarse
  approximation

Recommended order:

1. stochastic multiplication-register model
2. EMCCD-specific CIC placement
3. post-gain saturation semantics

### InGaAs

Highest-value missing features:

- persistence / latent-image model
- bad-pixel and hot-pixel maps
- nonlinearity near high signal
- ROIC-dependent glow / offset structure

Why these matter:

- InGaAs cameras are often attractive because of wavelength coverage, but their
  realism is usually limited by persistence, glow, and pixel defects rather
  than by abstract frame noise alone

Recommended order:

1. persistence
2. bad-pixel / hot-pixel maps
3. nonlinear response

### APDDetector

Highest-value missing features:

- `ParalyzableDeadTime`
- afterpulsing
- inter-channel crosstalk
- detector gating / duty-cycle control
- saturation / count-rate compression metadata

Why these matter:

- APD behavior is dominated by counting physics and recovery behavior
- these features affect curvature and other channel-readout systems more
  directly than frame-detector refinements do

Recommended order:

1. `ParalyzableDeadTime`
2. afterpulsing
3. detector gating
4. inter-channel crosstalk

## Shared Abstractions To Add

Some detector-specific work should still reuse a shared interface.

### Static Detector Defect Maps

Add reusable defect-map families for frame detectors:

- `AbstractDetectorDefectModel`
- `NullDetectorDefectModel`
- `PixelResponseNonuniformity`
- `DarkSignalNonuniformity`
- `BadPixelMask`

These should be detector-family-agnostic at the interface level, while each
family decides whether it uses them by default.

### Timing Models

Add explicit timing families where detector timing is part of the physics:

- `AbstractFrameTimingModel`
- `GlobalShutter`
- `RollingShutter`
- `CountingGate`

This should stay separate from optics and separate from detector noise.

### Counting Correlation Models

For `APDDetector`, add a shared counting-correlation surface such as:

- `AbstractCountingCorrelationModel`
- `NullCountingCorrelation`
- `AfterpulsingModel`
- `ChannelCrosstalkModel`

These should be layered on top of ideal counting readout and dead-time logic,
not mixed into the readout layout itself.

## Capability Traits

Add traits only where they express real optional behavior:

- `supports_detector_defect_maps(x)`
- `supports_detector_persistence(x)`
- `supports_detector_nonlinearity(x)`
- `supports_shutter_timing(x)`
- `supports_counting_gating(x)`
- `supports_afterpulsing(x)`
- `supports_channel_crosstalk(x)`
- `supports_paralyzable_dead_time(x)`

Do not add traits that merely restate the concrete detector family.

## Metadata Requirements

Detector metadata should expose family-specific behavior clearly enough for HIL,
profiling, and downstream consumers.

Examples:

- CMOS
  - shutter mode
  - output-group count
  - PRNU / DSNU enabled
- EMCCD
  - EM gain mode
  - excess-noise mode
  - CIC enabled
- InGaAs
  - persistence enabled
  - bad-pixel map enabled
  - glow mode
- APDDetector
  - dead-time family
  - gating enabled
  - afterpulsing enabled
  - crosstalk enabled

## Benchmark And Validation Requirements

Each phase should include:

- deterministic regression coverage
- allocation checks on maintained hot paths
- CPU / AMDGPU / CUDA coverage for maintained paths where supported
- representative system-level benchmarks where the detector family is actually
  used

Recommended maintained surfaces:

- `CMOSSensor`
  - REVOLT-like SH HIL
  - any maintained visible-frame WFS tutorial/runtime surface
- `EMCCDSensor`
  - pixel-output WFS profiling surface
- `InGaAsSensor`
  - NIR frame-detector tutorial/runtime surface
- `APDDetector`
  - AO188-like curvature counting path

## Phased Implementation

### Phase 1: Static Frame-Detector Structure

- add defect-map interfaces
- implement PRNU / DSNU for `CMOSSensor`
- implement bad-pixel / hot-pixel surfaces for `InGaAsSensor`
- expose metadata and null/default models

### Phase 2: Counting Detector Realism

- add `ParalyzableDeadTime`
- add counting gating
- add afterpulsing
- update `APDDetector` metadata and tests

### Phase 3: Frame Timing And Multiplication

- add `GlobalShutter` / `RollingShutter`
- add EMCCD stochastic multiplication-register model
- add EMCCD-specific CIC placement

### Phase 4: Persistence And Nonlinearity

- add `InGaAsPersistenceModel`
- add frame-detector nonlinearity surface
- expose detector-family-specific defaults where justified

### Phase 5: System Adoption And Benchmarking

- adopt the new detector features on maintained example/profiler surfaces
- measure:
  - timing
  - allocation impact
  - output differences
- keep only the detector-specific features that provide real modeling value

## Recommended Execution Order

If the work is done incrementally, the best order is:

1. `APDDetector`
   - smallest family surface
   - immediate value for curvature/counting systems
2. `CMOSSensor`
   - broad usefulness on HIL and visible-frame paths
3. `InGaAsSensor`
   - practical value for NIR frame systems
4. `EMCCDSensor`
   - important, but its next-step realism is more specialized and should be
     driven by actual use cases

If the work is done by modeling priority instead of implementation simplicity,
the best order is:

1. `CMOSSensor`
2. `APDDetector`
3. `InGaAsSensor`
4. `EMCCDSensor`

## Exit Criteria

This plan is in a good stopping state when:

- each target detector family has at least one detector-specific physics
  refinement beyond the current baseline
- null/default behavior remains explicit and tested
- the maintained detector metadata surfaces are clear and stable
- representative profiler surfaces show the performance cost of the new
  features
- no family-specific logic has leaked back into optics or generic readout
  layers

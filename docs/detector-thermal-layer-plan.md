# Detector Thermal Layer Plan

This plan defines a reusable detector thermal layer for
`AdaptiveOpticsSim.jl`.

The goal is to model detector temperature and cooling behavior as a generic
detector concern that can be composed with existing detector families such as:

- `CCDSensor`
- `CMOSSensor`
- `EMCCDSensor`
- `InGaAsSensor`
- `HgCdTeAvalancheArraySensor`
- `APDDetector`

This should stay separate from:

- sensor optics
- readout layout
- detector response / MTF
- detector family examples such as AO188 or AO3k

## Goal

Add a reusable thermal layer that can drive temperature-sensitive detector
physics without baking instrument-specific assumptions into core.

Examples of temperature-sensitive effects that should eventually be expressible
through this layer:

- frame-detector dark current
- counting-detector dark count rate
- glow rate
- persistence strength / decay
- CIC or related charge-generation terms where appropriate
- readout-dependent self-heating terms later if needed

The thermal layer should support both:

- deterministic fixed-temperature modeling
- future dynamic thermal evolution

## Current State

The detector subsystem already models several effects that are thermally
relevant, but only as direct scalar parameters:

- generic frame-detector `dark_current`
- `APDDetector.dark_count_rate`
- `EMCCDSensor.cic_rate`
- `InGaAsSensor.glow_rate`
- `InGaAsSensor.persistence_model`
- `HgCdTeAvalancheArraySensor.glow_rate`
- `HgCdTeAvalancheArraySensor.read_time`

These are useful, but the package does not currently have:

- detector temperature
- a thermal state
- cooling setpoints
- temperature-dependent parameter laws
- thermal time constants
- detector-family thermal metadata

So today the package has thermally relevant parameters, but not a thermal
model.

## Design Rules

- Keep the thermal layer separate from optics, readout, and detector response.
- Use abstract types, traits, and multiple dispatch.
- Avoid `isa`-driven hot-path branching where dispatch can express the model.
- Keep null/default models explicit.
- Keep AO188 and AO3k as example-layer consumers, not drivers of the thermal
  interface.
- Preserve deterministic and allocation-aware hot paths.

## Layering

The intended detector architecture becomes:

- sensor optics
- readout model
- detector response/statistics
- detector thermal model

The thermal layer does not replace detector response. It parameterizes parts of
it.

Examples:

- a frame detector still owns dark-current application
- the thermal model provides the temperature or effective rate law
- a counting detector still owns dead time and afterpulsing
- the thermal model can influence dark-count and afterpulse parameters if
  needed

## Abstract Interface

Add a reusable thermal family such as:

- `AbstractDetectorThermalModel`
- `NullDetectorThermalModel`
- `FixedTemperature`
- `FirstOrderThermalModel`

Optional state types:

- `AbstractDetectorThermalState`
- `NoThermalState`
- `DetectorThermalState`

The null/default path should require no dynamic state.

## Null Model

The null model should remain the default:

- `NullDetectorThermalModel`

That means:

- no explicit temperature is tracked
- existing detector parameters continue to work as-is
- regression and HIL baselines stay deterministic

This is important so the thermal layer does not force every maintained example
or benchmark into a richer detector model.

## Thermal Model Families

### Phase 1: Fixed Temperature

Add a simple static model:

- `FixedTemperature(temperature_K)`

This gives the package a first-class temperature concept without forcing
dynamic thermal evolution.

This is the right first implementation because it already enables:

- cooled EMCCD models
- cooled InGaAs models
- cooled HgCdTe models
- explicit APD operating temperature

### Phase 2: First-Order Thermal Evolution

Add a dynamic model such as:

- `FirstOrderThermalModel`

Candidate fields:

- `ambient_temperature_K`
- `setpoint_temperature_K`
- `time_constant_s`
- `power_coefficient`
- `min_temperature_K`
- `max_temperature_K`

This should support simple evolution toward a setpoint and later allow modest
self-heating approximations driven by read cadence or avalanche operation.

### Phase 3: Temperature Law Families

Temperature itself should not directly encode all detector physics.

Add reusable parameter-law families such as:

- `AbstractTemperatureLaw`
- `NullTemperatureLaw`
- `ArrheniusDarkCurrentLaw`
- `ExponentialGlowLaw`
- `PersistenceTemperatureLaw`
- `LinearTemperatureLaw`

These should be used to map thermal state to detector-family parameters.

## Detector Integration Strategy

Thermal behavior should be layered into detector families by dispatch.

### Generic Frame Detector Surface

Add thermal accessors along these lines:

- `thermal_model(det)`
- `thermal_state(det)`
- `detector_temperature(det)`
- `advance_thermal!(det, dt)`

And temperature-aware hooks such as:

- `effective_dark_current(det)`
- `effective_glow_rate(det)`

The null implementation should simply return the configured detector parameter.

### `EMCCDSensor`

Initial thermal targets:

- dark current through the generic detector path
- CIC temperature dependence later if justified

The thermal layer should not force a detailed EM gain temperature model
initially.

### `InGaAsSensor`

Initial thermal targets:

- glow-rate dependence
- persistence coupling/decay dependence
- generic dark current through the detector core

This is likely the first family where thermal modeling gives a visible realism
benefit.

### `HgCdTeAvalancheArraySensor`

Initial thermal targets:

- glow-rate dependence
- dark-current dependence
- optional read-cadence/self-heating coupling later

This family is already broad enough that the first thermal step should stay
conservative.

### `APDDetector`

Initial thermal targets:

- dark count rate dependence
- optional afterpulsing dependence later

This should stay generic and not assume one instrument-specific APD package.

### `CCDSensor` and `CMOSSensor`

These should also be compatible with the interface, even if early use is
limited to:

- dark current
- possibly fixed-pattern thermal bias later

## Capability Traits

Add traits only for optional behavior:

- `supports_detector_thermal_model(x)`
- `supports_dynamic_thermal_state(x)`
- `supports_temperature_dependent_dark_current(x)`
- `supports_temperature_dependent_glow(x)`
- `supports_temperature_dependent_persistence(x)`
- `supports_temperature_dependent_dark_counts(x)`

Do not add traits that simply restate the detector family.

## Metadata Requirements

Detector metadata should expose thermal configuration clearly enough for HIL,
profiling, and exported runtime surfaces.

Candidate metadata fields:

- thermal model family
- detector temperature
- ambient temperature
- cooling setpoint
- thermal time constant
- active temperature-law families

The null thermal path should report `:none` or equivalent.

## Implementation Phases

### Phase 1: Core Thermal Interface

Add:

- `AbstractDetectorThermalModel`
- `NullDetectorThermalModel`
- `FixedTemperature`
- metadata and trait surface
- detector constructor support for `thermal_model=...`

Deliverables:

- thermal types
- exports
- metadata
- docs and conformance tests

### Phase 2: Temperature-Aware Core Hooks

Add generic hooks in the detector layer:

- `effective_dark_current(...)`
- `effective_dark_count_rate(...)`
- `effective_glow_rate(...)`
- `effective_persistence_model(...)`

Initial implementation can be identity for null models and simple law-driven
mapping for fixed-temperature models.

### Phase 3: Family Wiring

Wire the thermal layer into:

- generic frame-detector dark current
- `APDDetector.dark_count_rate`
- `InGaAsSensor.glow_rate` and persistence
- `HgCdTeAvalancheArraySensor.glow_rate`

Optional early `EMCCDSensor` thermal coupling can stay limited to dark current.

### Phase 4: Dynamic Thermal Evolution

Add:

- `FirstOrderThermalModel`
- `advance_thermal!(det, dt)`
- explicit thermal state storage where needed

This should be opt-in and must not disturb deterministic default runs.

### Phase 5: Profiling and Example Adoption

Use the new thermal layer on maintained surfaces where it matters:

- AO3k NIR detector configurations
- detector family regression cases
- representative HIL profiles where dark current/glow realism matters

AO188 or AO3k temperature defaults should remain example-layer choices.

## Validation Requirements

Each phase should include:

- deterministic null-path regression
- constructor validation
- metadata checks
- allocation checks on maintained hot paths
- CPU / AMDGPU / CUDA coverage for maintained detector paths where relevant

For dynamic thermal models, also add:

- steady-state convergence tests
- deterministic fixed-step evolution tests

## Open Questions

- Should thermal state live in `DetectorState` / `APDDetectorState`, or be
  layered through a smaller shared thermal-state object?
- Should temperature laws be shared across detector families, or should each
  family own its own law surface with only light common helpers?
- Is self-heating worth modeling in Phase 1, or should it wait for dynamic
  thermal evolution?
- Do we want Celsius aliases in public APIs, or should the core remain
  Kelvin-only?

## Recommended Next Step

Implement Phases 1 and 2 first:

- core thermal interface
- fixed-temperature model
- temperature-aware rate hooks

That gives immediate modeling value for cooled detector families without adding
unnecessary dynamic complexity too early.

## SPAD Array Sensor Plan

Status: completed (first maintained pass)

## Purpose

This plan defines a first-class `SPADArraySensor` family that fits the current
detector interface without overloading either:

- `APDDetector` channel-counting semantics
- frame-detector linear-response semantics

The target is a maintained detector family for SPAD imaging arrays useful for
AO/HIL work, not a full event-camera framework on the first pass.

Representative hardware motivation includes modern SPAD imaging cameras such as
Hamamatsu-style SPAD imagers and PI Imaging-style array products like SPAD512
or SPAD Alpha. These names are motivation only, not a promise of exact
vendor-level emulation in the first implementation.

## Scope

The implemented first maintained SPAD array family supports:

- per-pixel photon detection efficiency
- dark count rate
- dead time
- afterpulsing
- nearest-neighbor or local crosstalk
- fill factor / effective active-area scaling
- accumulated-count frame output over a configurable integration interval

The implemented first maintained SPAD array family does not require:

- exact vendor timing microarchitecture
- full event-stream export
- TDC histogram readout
- pixel-by-pixel binary latch timing realism
- every proprietary quenching/readout mode

## Design Position

`SPADArraySensor` should be its own detector family.

It should not be modeled as:

- `APDDetector` with a bigger layout
- `CMOSSensor` with Poisson noise
- `HgCdTeAvalancheArraySensor` with a counting flag

Reason:

- SPAD arrays are photon-counting imagers
- dead time and pile-up are core physics, not optional decorations
- accumulated output is a count image, but the underlying response is not
  frame-linear

## Proposed Type Surface

### Sensor Family

Add a dedicated sensor family:

```julia
abstract type SPADArraySensorType <: CountingSensorType end

struct SPADArraySensor{T<:AbstractFloat,DT,AP,CT} <: SPADArraySensorType
    pde::T
    dark_count_rate::T
    fill_factor::T
    dead_time_model::DT
    afterpulsing_model::AP
    crosstalk_model::CT
end
```

That keeps the family sensor-local and lets detector assembly stay generic.

### Detector Family

Add a dedicated counting-imager detector type:

```julia
struct SPADArrayDetectorParams{T,S,G,TM}
    integration_time::T
    gate_model::G
    thermal_model::TM
    sensor::S
    output_precision::Union{Nothing,DataType}
    frame_mode::Symbol
end

mutable struct SPADArrayDetectorState{T,A,O,TS}
    counts::A
    noise_buffer::A
    output_buffer::O
    thermal_state::TS
end

struct SPADArrayDetector{N,P,S,B} <: AbstractCountingDetector
    noise::N
    params::P
    state::S
end
```

The important distinction from `APDDetector` is that the maintained output
surface is a 2D count array rather than a small channel vector.

## Output Contract

The first maintained output contract should be:

- accumulated counts per pixel over the integration interval
- optional count-rate normalization helper, but not the primary detector output

Public output accessors should follow counting-detector conventions first:

- `channel_output(det)` may flatten or alias the count image where needed
- `output_frame(det)` should return the native 2D accumulated count frame

The initial implementation should avoid pretending this is a standard frame
detector readout with Gaussian read noise and linear gain.

## Physics Contract

### Phase 1: Deterministic Counting-Imager Baseline

Implement a deterministic baseline with:

- PDE scaling
- fill-factor scaling
- dark-count accumulation
- optional dead-time correction

This gives a stable HIL/testing surface before introducing stochastic details.

### Phase 2: Stochastic SPAD Effects

Add opt-in stochastic behavior:

- Poisson photon counting
- dead-time loss / pile-up
- afterpulsing
- crosstalk

These should be small model objects and dispatch seams, not flags buried in one
giant detector update function.

### Phase 3: Gating And Timing

Support gated operation where useful:

- inherit the current counting gate model seam where possible
- add SPAD-specific timing helpers only if the current counting-gate interface
  proves too narrow

## Interface Placement

This family should reuse existing counting-detector abstractions where they fit:

- `CountingSensorType`
- `CountingDeadTimeModel`
- `AbstractCountingGateModel`
- `AbstractCountingCorrelationModel`
- thermal-model integration

But it should introduce family-specific metadata helpers:

- `detector_sensor_symbol(::SPADArraySensor) = :spad_array`
- `supports_dead_time(::SPADArrayDetector) = true`
- `supports_afterpulsing(::SPADArrayDetector) = ...`
- `supports_channel_crosstalk(::SPADArrayDetector) = ...`

If later needed, add:

- `supports_binary_latch_mode`
- `supports_event_export`
- `supports_histogram_readout`

These should not be added until a real model needs them.

## Validation Plan

The initial SPAD array family should be validated internally, not through OOPAO
or SPECULA parity, because the current local OOPAO and SPECULA trees do not
appear to expose SPAD detector models.

Validation should include:

1. Deterministic detector contracts
   - zero-flux counts equal dark-count expectation in deterministic mode
   - PDE and fill-factor scaling are monotone and bounded
   - dead-time loss reduces counts under high flux

2. Stochastic detector contracts
   - mean count curves match expected low-flux and moderate-flux behavior
   - afterpulsing raises counts in the expected direction
   - crosstalk leaks counts spatially in a controlled way

3. Runtime/export contracts
   - maintained 2D accumulated-count output shape
   - deterministic HIL export metadata
   - backend parity on the maintained deterministic surfaces

4. GPU validation
   - CUDA and AMDGPU only after the CPU deterministic contract is stable
   - start with deterministic counting frame export before stochastic kernels

## Implemented Scope

Completed in the first maintained pass:

1. `SPADArraySensor` plus detector metadata/capability methods.
2. Deterministic accumulated-count array output on CPU.
3. Detector export metadata and runtime integration.
4. Dead-time and dark-count modeling.
5. Afterpulsing and crosstalk model seams.

Deferred:

1. CUDA/AMDGPU deterministic parity.
2. Optional stochastic GPU support.

## Non-Goals

- exact vendor emulation of Hamamatsu or PI Imaging devices on the first pass
- full asynchronous event-camera infrastructure
- forcing SPAD arrays through the frame-detector readout stack
- forcing SPAD arrays through the scalar APD channel layout

## Success Criteria

- SPAD arrays have a first-class detector family
- the implementation fits the current detector extension seams
- accumulated-count imaging works in deterministic CPU mode
- HIL/runtime export can consume the maintained output cleanly
- future vendor-specific presets can be layered on top without changing the
  family architecture

# Zernike WFS Plan

Status:
- [x] Phase 1 scalar CPU MVP
- [~] Phase 2 GPU path
- [x] Phase 3 HIL/RTC integration

## Goal

Add a Zernike wavefront sensor to `AdaptiveOpticsSim.jl` in a way that fits the
existing WFS architecture, remains GPU-ready, and is useful for HIL/RTC-facing
workflows.

The preferred Julia-side name is `ZernikeWFS`. The shorter acronym `ZWFS` is
fine in discussion, but the concrete type should read consistently next to
`PyramidWFS` and `BioEdgeWFS`.

## Reference Surface

OOPAO does not currently provide a Zernike/ZWFS implementation in the checked-in
local surface, so there is no direct OOPAO parity target for this sensor.

The primary external reference is SPECULA:

- [zernike_sensor.py](/home/dgamroth/workspaces/codex/SPECULA/specula/processing_objects/zernike_sensor.py)
  - optical sensor model
  - implemented as a no-modulation focal-plane phase-shifting spot derived from
    SPECULA's modulated pyramid machinery
- [zwfs_slopec.py](/home/dgamroth/workspaces/codex/SPECULA/specula/processing_objects/zwfs_slopec.py)
  - pixel-to-signal normalization/output stage
  - useful for signal convention, but not the main optical reference
- [test_zernike_sensor.py](/home/dgamroth/workspaces/codex/SPECULA/test/test_zernike_sensor.py)
  - basic flat-wavefront and focus-behavior checks
- [test_zwfs_slopec.py](/home/dgamroth/workspaces/codex/SPECULA/test/test_zwfs_slopec.py)
  - signal normalization checks

## Architectural Fit

`ZernikeWFS` should be implemented as another `AbstractWFS` under
[`src/wfs`](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/wfs),
following the established split:

- immutable params
- mutable state
- `measure!`-based API
- detector-backed and direct-frame measurement paths
- CPU scalar path first
- accelerator path added without changing the public API

The closest internal structural reference is
[`pyramid.jl`](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/wfs/pyramid.jl),
because the Zernike sensor is also a focal-plane-mask, pupil-reimaging,
diffractive WFS.

## Proposed MVP Scope

Phase 1 should implement a usable diffractive sensor, not a placeholder.

Included in MVP:

- `ZernikeWFSParams`
- `ZernikeWFSState`
- `ZernikeWFS <: AbstractWFS`
- focal-plane phase-mask spot with configurable:
  - phase shift
  - spot radius in `lambda / D`
- diffractive propagation path:
  - pupil field
  - focal plane
  - phase mask application
  - return to pupil/image plane
- output intensity frame
- direct `measure!(wfs, tel, src)` path
- detector-backed `measure!(wfs, tel, src, det; rng=...)` path
- pupil-valid-signal extraction and normalized signal vector
- reference/null measurement calibration

Deferred from MVP:

- custom modulation variants
- extended-source support
- LGS-specific variants
- advanced estimator/reconstructor conventions specific to a given RTC
- physics-complete parity against literature beyond the chosen reference surface

## Public API Shape

Initial constructor sketch:

```julia
ZernikeWFS(
    tel::Telescope;
    pupil_samples::Int,
    phase_shift_pi::Real = 0.5,
    spot_radius_lambda_over_d::Real = 1.0,
    threshold::Real = 0.0,
    normalization::WFSNormalization = MeanValidFluxNormalization(),
    diffraction_padding::Int = 2,
    binning::Int = 1,
    T::Type{<:AbstractFloat} = Float64,
    backend = CPUBackend(),
)
```

Public methods should mirror the other WFS types:

- `measure!(wfs, tel, src)`
- `measure!(wfs, tel, src, det; rng=...)`
- `update_valid_mask!(wfs, tel)`

Possible follow-up public helper:

- `zernike_phase_mask_frame!(out, wfs, tel)`

## Signal Convention

The first version should keep the signal convention simple and explicit:

- retain the camera/intensity frame in state
- extract pupil-valid pixels or pupil-valid subaperture bins
- normalize by reference or mean-valid flux
- store the 1-D signal vector in `state.slopes` for compatibility with the
  current control/runtime surfaces

This matches how SPECULA's `zwfs_slopec.py` treats the ZWFS signal: normalized
intensity over the valid pupil support, rather than SH-style geometric slopes.

That means the naming in the internals may still use `slopes` for consistency
with the current `AbstractWFS` runtime surfaces, but the docs should be explicit
that for `ZernikeWFS` this vector is an intensity-derived signal, not a
centroid-derived slope field.

## Calibration Plan

The MVP should include an explicit reference/null calibration path:

1. build valid pupil support
2. measure a flat/reference source
3. save the reference signal and reference frame
4. subtract/reference-normalize in subsequent `measure!` calls

This is closer to the existing Pyramid/BioEdge pattern than to Shack-Hartmann
reference centroids.

## Implementation Phases

### Phase 1: Scalar CPU MVP

- add `src/wfs/zernike.jl`
- include/export from
  [`src/AdaptiveOpticsSim.jl`](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/AdaptiveOpticsSim.jl)
- implement params/state/constructor
- implement valid-support update
- implement focal-plane phase-mask generation
- implement direct diffractive `measure!`
- implement detector-backed `measure!`
- implement reference calibration and normalized signal extraction
- add unit tests

### Phase 2: GPU Path

- add accelerator-compatible mask generation
- avoid scalar indexing in signal extraction
- add CUDA/AMDGPU smoke coverage analogous to the maintained Pyramid/BioEdge
  surface

### Phase 3: HIL/RTC Integration

- verify `ClosedLoopRuntime` compatibility
- add reconstructor/control example using the Zernike signal vector
- decide whether a dedicated `ZernikeSignalOrder` or other sensor-specific
  convention is needed

Current maintained status:
- `ClosedLoopRuntime` compatibility is covered by tests
- a compact runtime-based example lives in
  [`examples/tutorials/closed_loop_zernike.jl`](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/examples/tutorials/closed_loop_zernike.jl)
- no dedicated `ZernikeSignalOrder` has been introduced yet; the current
  `state.slopes` compatibility surface is sufficient for the maintained MVP
- the compact runtime/profile surface is covered by
  [`scripts/profile_zernike_runtime.jl`](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/scripts/profile_zernike_runtime.jl)
  on CPU, AMDGPU, and CUDA

## Test Plan

### Unit Tests

- constructor/state sizing
- flat-wavefront output frame shape and non-negativity
- reference calibration populates expected state
- normalized valid-pupil signal on a simple synthetic frame

### Behavioral Tests

- focus aberration produces the expected signed intensity pattern, following the
  spirit of SPECULA's
  [test_zernike_sensor.py](/home/dgamroth/workspaces/codex/SPECULA/test/test_zernike_sensor.py)
- signal normalization follows the expected mean-valid-flux convention, in the
  spirit of
  [test_zwfs_slopec.py](/home/dgamroth/workspaces/codex/SPECULA/test/test_zwfs_slopec.py)

### GPU Tests

- CUDA smoke path
- AMDGPU smoke path
- CPU/GPU numerical agreement on a compact deterministic case

## Risks and Open Choices

### 1. Naming of the output vector

The runtime/control surface expects a 1-D `slopes` vector today. For a Zernike
sensor this is conceptually a normalized pupil-intensity signal, not a slope
vector. The first implementation should keep compatibility and document the
semantic difference rather than forcing a broad control-surface rename.

### 2. Detector plane vs re-imaged pupil plane conventions

The optical chain and signal extraction need to be explicit about where the
measured intensity lives. The MVP should pick one clear convention and document
it rather than trying to abstract over multiple literature variants on day one.

### 3. Reference target for fidelity

Because OOPAO does not provide a direct ZWFS target here, fidelity should be
checked against:

- SPECULA behavior where comparable
- qualitative optical behavior tests
- CPU/GPU agreement
- internal deterministic regression fixtures

## Recommendation

Implement `ZernikeWFS` now as a focused Phase 1 scalar CPU MVP, using SPECULA's
`zernike_sensor.py` as the optical reference and `zwfs_slopec.py` as the signal
normalization reference. Keep the first version narrow and correct, then add GPU
support once the scalar path and tests are in place.

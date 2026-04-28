# Atmosphere Runtime Specification

This document freezes the Milestone 0 contract for atmosphere-related runtime
behavior in AdaptiveOpticsSim.jl. It is the normative reference for the
atmosphere/runtime refactor that follows.

General project-wide units guidance lives in
[units-policy.md](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/docs/units-policy.md).

## Scope

This specification covers:

- atmosphere parameter semantics,
- per-layer runtime state ownership,
- wind transport semantics,
- units at the telescope/runtime boundary,
- stochastic ownership for deterministic execution,
- the boundary between runtime OPD fields and phase-statistics helpers.

This specification does not yet define:

- the exact persistent-screen algorithm,
- cross-hardware bitwise determinism,
- inter-process or hardware transport interfaces,
- coronagraph or science-camera propagation details.

## Reference Policy

- OOPAO is a behavioral reference for atmosphere outcomes and tutorial parity.
- SPECULA is a structural reference for persistent-screen decomposition.
- AdaptiveOpticsSim.jl is not required to preserve either architecture.
- When OOPAO and SPECULA disagree with a stronger Julia design or a clearer
  physical contract, AdaptiveOpticsSim.jl may diverge deliberately.

## Definitions

- `OPD`: optical path difference on the pupil grid, in meters.
- `phase`: wavelength-dependent optical phase, in radians.
- `fractional_cn2`: per-layer turbulence-strength fraction. It is a normalized
  layer weight and not a per-layer Fried parameter.
- `persistent screen`: a layer state that evolves from previous state by
  transport and boundary injection rather than by full redraw each step.
- `runtime atmosphere`: the atmosphere object used in `advance!` and
  `propagate!` during simulation.

## Normative Requirements

### Units And Data Semantics

- `ATM-001` Runtime pupil-plane atmosphere fields MUST be represented as OPD in
  meters.
- `ATM-002` `Telescope.state.opd` MUST be interpreted as OPD in meters.
- `ATM-003` Wavelength-dependent phase conversion MUST occur only at
  wavelength-aware optics, PSF, or WFS boundaries.
- `ATM-004` Any helper API that returns or consumes phase in radians MUST say so
  explicitly in its name or documentation.

### Layer Parameterization

- `ATM-010` Multi-layer turbulence weights MUST be parameterized by
  `fractional_cn2`.
- `ATM-011` `fractional_cn2` MUST be non-negative and MUST sum to one within
  documented tolerance.
- `ATM-012` `fractional_cn2` MUST be treated as layer turbulence-strength
  weights, not as direct per-layer `r0` scale factors.
- `ATM-013` No deprecated `fractional_r0` compatibility alias will be kept in
  the active API.

### Runtime Evolution

- `ATM-020` A moving atmosphere runtime MUST own persistent per-layer state.
- `ATM-021` A moving layer MUST evolve by wind transport and boundary extension,
  not by redrawing a fresh whole-layer screen every step.
- `ATM-022` Wind metadata MUST be represented in physical units:
  - speed in meters per second,
  - direction in degrees at the API boundary,
  - Cartesian velocity components in meters per second internally.
- `ATM-023` Subpixel wind motion MUST accumulate across steps.
- `ATM-024` Zero-wind layers MUST remain stationary absent other stochastic
  forcing.
- `ATM-025` Off-axis source footprint extraction MUST remain explicit and
  source-aware.

Current implementation note:

- Milestone 1 now uses a persistent finite periodic canvas per layer with
  subpixel interpolation.
- The finite periodic model remains a maintained fast/HIL-oriented atmosphere path.
- `InfiniteMultiLayerAtmosphere` now exists as a maintained CPU/GPU model with
  builder-time conditional boundary operators, runtime row/column injection,
  and residual subpixel extraction.
- finite and infinite atmospheres now both support source-aware propagation for
  off-axis and finite-height footprints from the current evolved layer state.
- GPU construction now warms the boundary-injection kernels up front so the
  first integer-shift step does not pay a large one-time runtime allocation.
- The staged plan for the next atmosphere backend is in
  [infinite-boundary-atmosphere-plan.md](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/docs/infinite-boundary-atmosphere-plan.md).

### Determinism And Stochastic Ownership

- `ATM-030` Stochastic atmosphere updates MUST accept an explicit `rng`.
- `ATM-031` Deterministic execution MUST be driven by caller- or
  workspace-owned RNG state rather than hidden global draws.
- `ATM-032` The atmosphere runtime MUST remain compatible with the project
  deterministic-mode policy in
  [deterministic-simulation.md](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/docs/deterministic-simulation.md).

## Phase-Statistics Boundary

The runtime OPD contract and the existing phase-statistics helpers are not yet a
fully unified surface.

- `ATM-040` The runtime atmosphere contract is OPD-space.
- `ATM-041` The current `phase_stats.jl` helpers are phase-statistics utilities
  by intent and MUST NOT be assumed interchangeable with runtime OPD arrays
  without an explicit documented conversion.
- `ATM-042` Milestone 2 will either unify runtime and helper normalization or
  split them into explicitly different APIs with conversion helpers.

## Traceability Matrix

| ID | Requirement | Current mapping | Verification | Status |
| --- | --- | --- | --- | --- |
| `ATM-001` | Runtime atmosphere stores OPD meters | [src/optics/telescope.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/optics/telescope.jl), [src/optics/psf.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/optics/psf.jl), [src/wfs/shack_hartmann.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/wfs/shack_hartmann.jl) | Manual code review | Partially mapped |
| `ATM-002` | Telescope OPD is OPD meters | [src/optics/telescope.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/optics/telescope.jl), [src/optics/opd_map.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/optics/opd_map.jl) | Manual code review | Mapped |
| `ATM-003` | Phase conversion occurs at wavelength-aware boundaries | [src/optics/psf.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/optics/psf.jl), [src/wfs/shack_hartmann.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/wfs/shack_hartmann.jl), [src/wfs/lift.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/wfs/lift.jl) | Manual code review | Mapped |
| `ATM-004` | Phase-space helpers are explicit | [src/atmosphere/phase_stats.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/atmosphere/phase_stats.jl) | Documentation review | Partially mapped |
| `ATM-010` | Layer weights use `fractional_cn2` | [src/atmosphere/multilayer.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/atmosphere/multilayer.jl), [src/tomography/parameters.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/tomography/parameters.jl) | API review plus test pass | Mapped |
| `ATM-011` | `fractional_cn2` is normalized | [src/tomography/parameters.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/tomography/parameters.jl) | Unit tests in [test/tomography.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/test/tomography.jl) | Partially covered |
| `ATM-012` | Layer weights are not direct `r0` scale factors | [src/atmosphere/multilayer.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/atmosphere/multilayer.jl), [test/runtests.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/test/runtests.jl) | Code review plus dedicated variance regression | Mapped |
| `ATM-013` | No deprecated alias | Active API surface | Build and full test pass after rename | Mapped |
| `ATM-020` | Persistent per-layer runtime state | [src/atmosphere/multilayer.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/atmosphere/multilayer.jl) | Code review plus test pass | Mapped |
| `ATM-021` | No redraw-every-step moving atmosphere | [src/atmosphere/multilayer.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/atmosphere/multilayer.jl) | Code review plus test pass; infinite boundary injection still pending | Partially mapped |
| `ATM-022` | Wind metadata uses physical units | [src/atmosphere/multilayer.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/atmosphere/multilayer.jl), [src/tomography/geometry.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/tomography/geometry.jl) | Code review | Mapped |
| `ATM-023` | Subpixel accumulation is preserved | [src/atmosphere/multilayer.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/atmosphere/multilayer.jl), [test/runtests.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/test/runtests.jl) | Test pass | Mapped |
| `ATM-024` | Zero-wind layers remain stationary | [src/atmosphere/multilayer.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/atmosphere/multilayer.jl), [test/runtests.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/test/runtests.jl) | Test pass | Mapped |
| `ATM-025` | Source-aware footprint extraction remains explicit | [src/optics/source.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/optics/source.jl), [src/wfs/](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/wfs) | API review | Partially mapped |
| `ATM-030` | Atmosphere updates accept explicit `rng` | [src/atmosphere/kolmogorov.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/atmosphere/kolmogorov.jl), [src/atmosphere/multilayer.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/atmosphere/multilayer.jl) | Existing API review | Mapped |
| `ATM-031` | RNG ownership stays external | [docs/deterministic-simulation.md](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/docs/deterministic-simulation.md) | Documentation review | Mapped |
| `ATM-032` | Atmosphere remains compatible with deterministic mode | [docs/deterministic-simulation.md](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/docs/deterministic-simulation.md) | Manual review | Mapped |
| `ATM-040` | Runtime contract is OPD-space | This specification, [src/optics/psf.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/optics/psf.jl) | Manual review | Mapped |
| `ATM-041` | `phase_stats.jl` is not assumed runtime-equivalent | [src/atmosphere/phase_stats.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/atmosphere/phase_stats.jl) | Documentation review | Partially mapped |
| `ATM-042` | Runtime/helper normalization must be unified or split explicitly | [src/atmosphere/kolmogorov.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/atmosphere/kolmogorov.jl), [src/atmosphere/phase_stats.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/atmosphere/phase_stats.jl), [test/runtests.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/test/runtests.jl) | Milestone 2 implementation plus regression test pass | Mapped |

## Immediate Consequences For Implementation

- New atmosphere-facing APIs should use `fractional_cn2`.
- New runtime work should target persistent moving layers, not redraw-and-shift.
- Runtime and `phase_stats.jl` normalization are now aligned; phase-space
  helpers should still remain explicit at wavelength-aware boundaries.
- Milestone 4 regression coverage now guards:
  - temporal coherence,
  - subpixel transport,
  - `fractional_cn2` weighting,
  - deterministic replay,
  - moving-atmosphere WFS and closed-loop traces.
- The infinite boundary-injection backend now adds:
  - deterministic replay under fixed RNG,
  - long-run stationarity windows,
  - long-run non-periodicity against the finite periodic canvas,
  - GPU smoke coverage plus CPU/GPU statistical agreement checks.

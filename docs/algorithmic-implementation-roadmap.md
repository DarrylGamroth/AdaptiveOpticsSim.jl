# Algorithmic Implementation Roadmap

This roadmap supersedes the earlier "port OOPAO as-is" mindset.
OOPAO remains a behavioral reference where useful, but it is not an
architectural constraint. SPECULA is also a useful reference, especially for
decomposition and persistent-screen atmosphere modeling. The implementation
target for AdaptiveOpticsSim.jl is:

- scientifically stronger than the current OOPAO-derived runtime where needed,
- idiomatic Julia in architecture and execution model,
- optimized for high-fidelity and high-speed HIL simulation first,
- extensible toward richer control and science optics without bloating core.

## Priorities

1. Fix scientific correctness in the atmosphere/runtime path.
2. Build a clean high-fidelity and high-speed HIL execution core.
3. Integrate science focal-plane propagation through `Proper.jl` as an optional extension.
4. Defer advanced controller design/synthesis work until the runtime and plant model are stable.

## Non-goals For This Roadmap

- No deprecated aliases. Breaking renames are acceptable while the package is still in development.
- No `Dagger.jl` integration for now.
- No coronagraph implementation in core AdaptiveOpticsSim.
- No major control-framework investment before atmosphere/runtime fidelity is repaired.
- No sockets, shared memory, RPC, or hardware-specific transport layers in core.

## Reference Policy

- OOPAO: use for behavior checks, tutorial parity, and regression baselines.
- SPECULA: use for algorithm decomposition and persistent-screen implementation ideas.
- Native Julia design: final authority for API, data layout, execution, and extensions.

## Milestone 0: Freeze Targets And Rename Semantics

Goal: lock the scientific and API direction before refactoring.

Deliverables:

- Replace `fractional_r0` terminology with `fractional_cn2` everywhere new code is touched.
- Define one atmosphere-runtime contract for:
  - persistent per-layer state,
  - wind-driven transport,
  - subpixel accumulation,
  - source footprint extraction,
  - deterministic RNG ownership.
- Decide and document atmospheric units explicitly:
  - whether internal arrays represent phase in radians or OPD in meters,
  - where wavelength scaling occurs,
  - which helper APIs are phase-space vs OPD-space.
- Write a short reference note stating:
  - OOPAO behavior is a check,
  - SPECULA informs structure,
  - AdaptiveOpticsSim is free to exceed both where justified.

Acceptance gate:

- A written contract document exists and all new implementation work follows it.

Suggested files:

- `docs/atmosphere-runtime-spec.md`
- `src/Atmosphere/`
- `src/Calibration/initialization.jl`
- `examples/`

## Milestone 1: Replace The Current Moving Atmosphere

Goal: remove the redraw-then-shift atmosphere model and replace it with a
persistent moving-screen model.

Current blocker:

- [multilayer.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/Atmosphere/multilayer.jl) redraws each layer every step and then shifts it.

Deliverables:

- Introduce a persistent per-layer screen type, for example:
  - `InfinitePhaseScreen`
  - `MovingAtmosphereLayer`
  - `MovingAtmosphere`
- Each layer owns:
  - current screen buffer,
  - covariance/update operators for boundary injection,
  - RNG,
  - integer/subpixel motion accumulators,
  - altitude and wind metadata.
- The runtime updates layers by transport and boundary extension, not by full redraw.
- Off-axis footprint extraction remains explicit and source-aware.
- The old `MultiLayerAtmosphere` implementation is removed or rewritten in place to use the new model.

Implementation guidance:

- Use OOPAO as the behavioral reference for moving-layer semantics.
- Use SPECULA's separation between persistent screen state and evolution driver as the structural template.
- Keep CPU first; make the state layout GPU-ready but do not force the first implementation onto GPU.

Acceptance gate:

- Consecutive layer states are shift-correlated under small `dt`.
- Zero-wind layers remain stationary.
- Subpixel wind accumulates correctly across steps.

Suggested files:

- `src/Atmosphere/`
- `test/`
- `docs/deterministic-simulation.md`

## Milestone 2: Unify Phase-Screen Statistics And Runtime Scaling

Goal: eliminate normalization and units ambiguity between static phase-statistics
helpers and the runtime atmosphere.

Deliverables:

- Make [phase_stats.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/Atmosphere/phase_stats.jl) and the runtime atmosphere share the same PSD/normalization conventions, or split them into explicitly different APIs with clear names.
- Add explicit conversion helpers if both phase-space and OPD-space representations are kept.
- Audit:
  - `ft_phase_screen`
  - `ft_sh_phase_screen`
  - runtime Kolmogorov PSD generation
  - layer weighting
  - wavelength scaling
- Remove undocumented hidden scale factors.

Acceptance gate:

- Statistical tests show agreement between:
  - runtime-generated large-screen samples,
  - helper-generated large-screen samples,
  - analytic covariance/variance expectations
  within documented tolerances.

Suggested files:

- `src/Atmosphere/kolmogorov.jl`
- `src/Atmosphere/phase_stats.jl`
- `test/runtests.jl`

## Milestone 3: Build A HIL-Focused Runtime Core

Goal: make the runtime suitable for both high-fidelity simulation and
high-speed HIL-oriented execution without conflating the two.

Scope boundary:

- This milestone is strictly about in-process runtime behavior.
- "HIL-focused" here means latency-aware, deterministic, low-allocation
  execution that can later sit inside a larger HIL stack.
- It does not include sockets, shared memory, RPC, device drivers, or a
  generic transport abstraction layer.

Deliverables:

- Split fidelity policy from execution policy cleanly.
- Define two maintained runtime targets:
  - `ScientificRuntimeProfile`
  - `HILRuntimeProfile`
- The HIL profile may simplify expensive science-only paths, but must do so explicitly.
- Make latency staging first-class for:
  - sensor integration,
  - detector/readout,
  - reconstruction,
  - DM command application.
- Audit the hot loop for:
  - allocations,
  - type stability,
  - unnecessary copies,
  - CPU/GPU transfer points,
  - oversized builder-time work leaking into runtime.
- Build prepared runtime objects that are initialized once and stepped many times.
- Expose stable in-memory runtime readouts for:
  - WFS frames,
  - slopes,
  - DM commands,
  - science-path inputs,
  - per-step metadata.
- Keep coarse-grained parallelism only:
  - multi-source,
  - multi-WFS,
  - per-step grouping,
  - sweep-level batching.

Acceptance gate:

- Representative HIL benchmark scenarios are defined and tracked.
- Runtime performance is measured on:
  - single-WFS closed loop,
  - multi-WFS grouped execution,
  - AO188/3k-style scenario.
- Fixed-seed runs are replayable with stable numerical results and stable timing/allocation envelopes.

Suggested files:

- `src/Control/runtime.jl`
- `src/Core/parallel.jl`
- `docs/benchmark-matrix-plan.md`
- `docs/control-simulation-architecture-plan.md`
- `scripts/profile_*`

## Milestone 4: Lock Algorithmic Regression Coverage

Goal: guard the corrected algorithms with meaningful numerical tests.

Deliverables:

- Add atmosphere regression tests for:
  - temporal coherence,
  - `fractional_cn2` weighting,
  - subpixel transport,
  - deterministic replay.
- Add WFS regression tests tied to corrected atmosphere motion, not just static ramps.
- Add closed-loop traces that exercise:
  - moving turbulence,
  - latency,
  - detector coupling,
  - reconstructor stability.
- Keep OOPAO reference bundles where they are useful, but add Julia-native regression cases for behaviors OOPAO does not cover well.

Acceptance gate:

- A failing atmosphere/scaling regression breaks CI.
- Static smoke tests are no longer the only protection around turbulence evolution.

Suggested files:

- `test/runtests.jl`
- `test/reference_harness.jl`
- `test/reference_data/`

## Milestone 5: Integrate `Proper.jl` For Science Optics

Goal: move coronagraph and advanced focal-plane propagation out of core and into
an optional science-optics extension.

Deliverables:

- Add a `Proper.jl` extension package/module, for example:
  - `ext/AdaptiveOpticsSimProperExt.jl`
- Define a bridge from AdaptiveOpticsSim state to `Proper.WaveFront`:
  - pupil,
  - OPD,
  - reflectivity/apodization,
  - DM surfaces,
  - wavelength,
  - sampling.
- Support at least:
  - diffraction-limited science PSF handoff,
  - aberrated science PSF handoff,
  - DM-corrected science path.
- Keep the extension boundary explicit so core does not depend on `Proper.jl`.

Acceptance gate:

- A maintained example demonstrates science propagation through `Proper.jl`.
- The interface is documented and covered by smoke tests.

Suggested files:

- `ext/AdaptiveOpticsSimProperExt.jl`
- `examples/`
- `docs/`

## Milestone 6: Revisit Control After The Plant Is Stable

Goal: avoid over-investing in controller infrastructure before the atmosphere
and runtime plant model are trustworthy.

Decision:

- Do not prioritize broad control work right now.
- Prioritize only the minimum control support needed for:
  - realistic latency staging,
  - HIL runtime integration,
  - representative closed-loop testing.

Deliverables:

- Keep lightweight native runtime controllers in core.
- Add a design-time pathway based on `ControlSystems.jl` for:
  - discretization,
  - transfer-function analysis,
  - state-space design,
  - IIR/FIR coefficient generation.
- Defer `ModelingToolkit.jl` unless and until plant linearization and symbolic model composition become real bottlenecks or repeated pain points.
- Revisit richer controllers only after Milestones 1-4 are complete:
  - IIR families,
  - modal/controller-bank variants,
  - POLC-style paths,
  - predictive control.

Acceptance gate:

- Controller changes are driven by verified plant/runtime needs, not by framework ambition.

Suggested files:

- `src/Control/`
- `docs/control-simulation-architecture-plan.md`

## Milestone 7: Selective Post-OOPAO Expansion

Goal: move beyond OOPAO where it is scientifically or operationally limiting.

Candidate directions after the core runtime is fixed:

- richer control families informed by SPECULA,
- stronger detector/readout physics for runtime realism,
- science-camera and coronagraph workflows through `Proper.jl`,
- improved tomography builder/runtime scaling,
- multi-rate and multi-stage control for woofer/tweeter or RTC-style systems.

Explicit rule:

- OOPAO parity is a validation layer, not a feature ceiling.

## Recommended Execution Order

1. Milestone 0
2. Milestone 1
3. Milestone 2
4. Milestone 4
5. Milestone 3
6. Milestone 5
7. Milestone 6
8. Milestone 7

Rationale:

- The atmosphere/runtime plant must be fixed before control sophistication matters.
- Regression coverage should land before broad runtime optimization, so performance work does not stabilize wrong behavior.
- `Proper.jl` integration should happen after the wavefront/OPD contract is trustworthy.

## Immediate Next Sprint

1. Define the `ScientificRuntimeProfile` and `HILRuntimeProfile` boundary in code.
2. Audit the closed-loop hot path for per-step allocations and type-instability.
3. Make latency staging first-class across sensing, reconstruction, and command application.
4. Build prepared runtime objects that separate setup work from `step!`.
5. Add benchmark scenarios for single-WFS, grouped multi-WFS, and AO188/3k-style runs.

## Exit Criteria For The Roadmap

This roadmap is complete when:

- atmosphere evolution is persistent and regression-backed,
- phase-statistics helpers and runtime scaling are consistent,
- HIL runtime profiles are explicit and benchmarked,
- `Proper.jl` provides optional science propagation,
- control remains lean in core and only grows after the plant is stable.

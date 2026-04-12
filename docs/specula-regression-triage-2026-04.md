# SPECULA Regression Triage 2026-04

## Status

The remaining `SPECULA reference regression` is **not** currently a frozen-data mismatch.
The SPECULA bundle validates cleanly in isolation. The failure is an **order-dependent regression** that only appears after broader test execution.

## Failing Surface

Current failing case:
- `shack_hartmann_polychromatic_frame`
- manifest entry: [test/reference_data_specula/manifest.toml](../test/reference_data_specula/manifest.toml)
- case kind: `shack_hartmann_frame`
- SPECULA contract:
  - `test_poly_chrom_sh.py::TestPolyChromSH::test_polychrom_sh_basic`
  - "polychromatic Shack-Hartmann sensing yields a stable sampled detector-plane frame on the maintained Julia surface"

## Observed Behavior

### 1. Full `Pkg.test()` / broad test execution

The broad suite still fails at:
- [test/testsets/reference_and_tutorials.jl](../test/testsets/reference_and_tutorials.jl)
- testset: `SPECULA reference regression`

Observed failure shape:
- `result.ok == false`
- failure occurs inside the SPECULA loop after earlier testsets have already run

### 2. SPECULA bundle in isolation

Direct isolated validation of the SPECULA bundle succeeds:
- `n_cases = 9`
- all 9 cases passed
- each case produced `maxabs = 0.0`

This includes:
- `shack_hartmann_polychromatic_frame ok=true maxabs=0.0`

### 3. Narrowed reproducer window

The failing SPECULA case remains clean after:
- `core_optics`
- `atmosphere`
- `control_and_runtime`
- `detectors_and_wfs`

The failing SPECULA case also remains clean after:
- `calibration_and_analysis` alone

The regression appears only after the broader combined sequence that includes both:
- the earlier WFS-heavy suites
- the later calibration / analysis suite

So the current diagnosis is:
- **not** a bad reference bundle
- **not** a stable numerical disagreement
- **yes** an order-dependent interaction between earlier WFS tests and later calibration / analysis work

## Concrete Evidence

### Clean isolated SPECULA run

Observed direct verdicts:
- `atmospheric_chromatic_intensity ok=true maxabs=0.0`
- `atmospheric_fresnel_intensity ok=true maxabs=0.0`
- `atmospheric_geometric_intensity ok=true maxabs=0.0`
- `curvature_flat_signal ok=true maxabs=0.0`
- `curvature_quadratic_signal ok=true maxabs=0.0`
- `pyramid_polychromatic_frame ok=true maxabs=0.0`
- `shack_hartmann_polychromatic_frame ok=true maxabs=0.0`
- `zernike_flat_signal ok=true maxabs=0.0`
- `zernike_quadratic_signal ok=true maxabs=0.0`

### Narrowed SH frame case in clean conditions

For `shack_hartmann_polychromatic_frame` in clean conditions:
- `size_actual = (120, 120)`
- `size_expected = (120, 120)`
- `actual_nan = 14400`
- `expected_nan = 14400`
- `ok = true`
- `maxabs = 0.0`

### Broad failing behavior

In the failing broad run:
- `validate_reference_case(case)` returned `ok = false`
- `maxabs = Inf`

Since the clean run shows identical shapes and identical all-`NaN` patterns, the broad failure is most consistent with:
- a `NaN`-pattern mismatch on the same `120 x 120` surface
- not a small floating-point drift

## Likely Root Cause

Most likely cause:
- order-dependent mutable state contamination in the polychromatic Shack-Hartmann detector-frame path

Most plausible technical classes:
- cached Shack-Hartmann sampling or calibration state leaking across tests
- detector / export policy state leaking into the SPECULA frame contract path
- a broad test modifying the conditions under which the `shack_hartmann_frame` reference case returns an all-`NaN` frame

The most important conclusion is:
- the SPECULA bundle itself is currently trustworthy
- the regression is in **test-order / runtime-state interaction**, not in the frozen reference data

## Relevant Code Paths

Primary reference harness path:
- [test/reference_harness.jl](../test/reference_harness.jl)
- `compute_reference_actual(case)`
- `validate_reference_case(case)`

Relevant SH frame branch:
- [test/reference_harness.jl](../test/reference_harness.jl)
- `case.kind === :shack_hartmann_frame`
- current implementation calls:
  - `AdaptiveOpticsSim.sampled_spots_peak!(wfs, tel, src)`
  - then returns `copy(wfs.state.spot_cube[1, :, :])`

Relevant Shack-Hartmann runtime path:
- [src/WFS/shack_hartmann/signals.jl](../src/WFS/shack_hartmann/signals.jl)
- [src/WFS/shack_hartmann/setup.jl](../src/WFS/shack_hartmann/setup.jl)

## Recommended Next Steps

1. Add a dedicated reproducer script for this single SPECULA case.
- Run the case:
  - in isolation
  - after `detectors_and_wfs`
  - after `calibration_and_analysis`
  - after both
- Print:
  - `size`
  - `count(isnan, ...)`
  - `valid_mask`
  - `peak`
  - `minimum` / `maximum` of finite values

2. Instrument the `:shack_hartmann_frame` reference branch.
- Record before returning:
  - `count(wfs.state.valid_mask)`
  - `count(isnan, wfs.state.spot_cube)`
  - whether `prepare_sampling!` / `ensure_sh_calibration!` materially changes the result

3. Bisect the interaction inside `calibration_and_analysis` after `detectors_and_wfs`.
- Most likely suspects are the tests that touch:
  - gain-sensing camera calibration
  - modal / interaction-matrix calibration
  - WFS-related calibration helpers

4. Decide whether the frozen SPECULA SH frame contract itself is the right target.
- A contract that is "all NaNs" is fragile as a long-term validation surface.
- If that is truly the intended SPECULA-aligned surface, document why.
- Otherwise replace it with a more stable, physically informative output surface.

## Current Recommendation

Do **not** remove or silence the SPECULA regression yet.

The right next action is:
- treat it as an order-dependent SH frame regression
- reproduce it with a dedicated focused script
- then either:
  - fix the state leak, or
  - replace the fragile all-`NaN` contract with a better-maintained SPECULA reference surface

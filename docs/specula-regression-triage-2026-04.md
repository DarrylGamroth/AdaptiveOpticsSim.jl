# SPECULA Regression Triage 2026-04

## Status

The remaining `SPECULA reference regression` has been **resolved**.
The root cause was a missing Shack-Hartmann sampling preparation step in the SPECULA
reference harness, which left the frozen `shack_hartmann_polychromatic_frame`
contract as an unprepared all-`NaN` surface.

## Resolved Surface

Resolved case:
- `shack_hartmann_polychromatic_frame`
- manifest entry: [test/reference_data_specula/manifest.toml](../test/reference_data_specula/manifest.toml)
- case kind: `shack_hartmann_frame`
- SPECULA contract:
  - `test_poly_chrom_sh.py::TestPolyChromSH::test_polychrom_sh_basic`
  - "polychromatic Shack-Hartmann sensing yields a stable sampled detector-plane frame on the maintained Julia surface"

## Observed Behavior

### 1. Root-cause harness path

The failing branch in [test/reference_harness.jl](../test/reference_harness.jl) did:
- optional `apply_reference_opd!`
- `AdaptiveOpticsSim.sampled_spots_peak!(wfs, tel, src)`
- return `copy(wfs.state.spot_cube[1, :, :])`

It did **not** call:
- `AdaptiveOpticsSim.prepare_sampling!(wfs, tel, src)`

That meant the reference surface was captured from an unprepared Shack-Hartmann state.

### 2. Old frozen SH frame contract

Before the fix, the frozen SPECULA SH frame artifact was:
- `shape = [120, 120]`
- `actual_nan = expected_nan = 14400`
- effectively an all-`NaN` contract

That surface was fragile and not scientifically useful as a maintained reference.

### 3. Prepared SH frame surface

After explicitly preparing the Shack-Hartmann sampling state:
- `size = (20, 20)`
- `nan = 0`
- `finite = 400`
- `min = 3.217715424077204e8`
- `max = 1.2234702222504738e13`

## Concrete Evidence

### Harness fix

Resolved by adding:
- `AdaptiveOpticsSim.prepare_sampling!(wfs, tel, src)`

immediately before:
- `AdaptiveOpticsSim.sampled_spots_peak!(wfs, tel, src)`

in the `case.kind === :shack_hartmann_frame` branch.

### Refreshed SPECULA artifact

Regenerated files:
- [test/reference_data_specula/manifest.toml](../test/reference_data_specula/manifest.toml)
- [test/reference_data_specula/shack_hartmann_polychromatic_frame.txt](../test/reference_data_specula/shack_hartmann_polychromatic_frame.txt)

Updated contract:
- `shape = [20, 20]`
- finite prepared detector-plane frame

## Root Cause

The SPECULA Shack-Hartmann frame reference branch was incomplete:
- it sampled spots without preparing the Shack-Hartmann sampling/calibration state first

The resulting frozen reference artifact was therefore not a reliable prepared sensing
surface. The broad test failure exposed that weakness, but the underlying bug was in the
reference harness, not in the runtime physics path.

## Relevant Code Paths

Primary reference harness path:
- [test/reference_harness.jl](../test/reference_harness.jl)
- `compute_reference_actual(case)`
- `validate_reference_case(case)`

Relevant SH frame branch:
- [test/reference_harness.jl](../test/reference_harness.jl)
- `case.kind === :shack_hartmann_frame`
- current implementation now calls:
  - `AdaptiveOpticsSim.prepare_sampling!(wfs, tel, src)`
  - `AdaptiveOpticsSim.sampled_spots_peak!(wfs, tel, src)`
  - then returns `copy(wfs.state.spot_cube[1, :, :])`

Relevant Shack-Hartmann runtime path:
- [src/wfs/shack_hartmann/signals.jl](../src/wfs/shack_hartmann/signals.jl)
- [src/wfs/shack_hartmann/setup.jl](../src/wfs/shack_hartmann/setup.jl)

## Verification

After the harness fix and SPECULA bundle refresh:
- the SPECULA reference regression passes
- the OOPAO reference regression still passes
- the broad `Pkg.test(test_args=["control_and_runtime"])` run completes successfully

Most relevant broad-suite results:
- `SPECULA reference regression | 19 / 19`
- `OOPAO reference regression | 35 / 35`
- `Tutorial examples | 52 / 52`

## Recommendation

Keep the refreshed SPECULA bundle.

The old all-`NaN` `120 x 120` Shack-Hartmann frame should not be restored. The prepared
finite `20 x 20` surface is the correct maintained reference for this SPECULA contract.

# SPECULA Reference Datasets

Date: 2026-03-31

Status: active

Plan traceability:

- [`PLAN-26`](./package-review-action-plan.md)
- review IDs: `PR-19`, `PR-20`, `PR-36`

## Purpose

This document describes the committed SPECULA-targeted frozen reference bundle
used by AdaptiveOpticsSim.jl.

Unlike the OOPAO bundle, this bundle is intentionally narrow and contract-led.
It exists for model families where SPECULA is the stronger external reference
for what should be validated, but where full platform-level numerical
equivalence is not yet the right maintenance burden.

## Bundle Root

- [test/reference_data_specula](../test/reference_data_specula)

The bundle is loaded through the same maintained reference harness used for the
OOPAO and pyTomoAO cases in [reference_harness.jl](../test/reference_harness.jl).

Environment override:

- `ADAPTIVEOPTICS_SPECULA_REFERENCE_ROOT`

Default root:

- [test/reference_data_specula](../test/reference_data_specula)

## Bundle Kind

This is a `contract_bundle`, not a full external-equivalence dump.

That means:

- the scenarios are derived from explicit contracts in SPECULA tests
- the committed expected outputs are frozen Julia-side arrays for those
  scenarios
- the bundle is used to keep these contract surfaces stable and visible in CI
- it is not a claim that AdaptiveOpticsSim matches SPECULA everywhere

## Current Committed Cases

1. `zernike_flat_signal`
   - provenance:
     `test_zernike_sensor.py::TestZernikeSensor::test_flat_wavefront_output_size`
   - contract:
     flat-wavefront Zernike sensing remains stable and shape-correct on the
     maintained exported signal surface
2. `zernike_quadratic_signal`
   - provenance:
     `test_zernike_sensor.py::TestZernikeSensor::test_focus`
   - contract:
     a low-order quadratic phase produces a structured non-flat Zernike signal
3. `curvature_flat_signal`
   - provenance:
     `test_curvature_sensor.py::TestCurvatureSensor::test_flat_wavefront_flux`
   - contract:
     after reference subtraction, flat-wavefront curvature output stays null on
     the maintained exported signal surface
4. `curvature_quadratic_signal`
   - provenance:
     `test_curvature_sensor.py::TestCurvatureSensor::test_full_chain_focus_response`
   - contract:
     a low-order quadratic phase produces a non-trivial curvature signal

## Why These Cases Are Narrow

AdaptiveOpticsSim and SPECULA do not expose identical public sensor surfaces for
all optics chains. In particular:

- parameterization differs
- some internal representations differ
- not every SPECULA processing object has a one-to-one maintained Julia API

So the committed bundle chooses deterministic scenarios that preserve the
validation intent of the SPECULA tests while staying aligned with the current
AdaptiveOpticsSim maintained interfaces.

## Generation Workflow

The maintained generator is:

- [generate_specula_reference_bundle.jl](../scripts/generate_specula_reference_bundle.jl)

Default command:

```bash
julia --project=. --startup-file=no scripts/generate_specula_reference_bundle.jl
```

Custom output root:

```bash
julia --project=. --startup-file=no scripts/generate_specula_reference_bundle.jl /tmp/specula-bundle
```

The generator:

- constructs deterministic Julia-side scenarios aligned to selected SPECULA
  test contracts
- computes the maintained output arrays
- writes the frozen text arrays and `manifest.toml`

## Validation Policy

These cases belong in the normal maintained test suite because they are small
and deterministic.

They should not be expanded into a large platform benchmark harness here.
Broader cross-package runtime and fidelity comparisons belong in the separate
benchmark work tracked after Phase 4.

## Known Limits

- This bundle does not yet cover atmospheric field propagation directly.
- It does not yet cover grouped/polychromatic SPECULA comparisons for Pyramid
  or Shack-Hartmann.
- It should be treated as targeted external validation, not as a full SPECULA
  parity claim.

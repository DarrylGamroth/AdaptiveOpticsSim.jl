Status: active

Plan traceability:

- [`PVP-09`](./post-review-platform-plan.md)
- [`CP-05`](./cross-package-benchmark-inventory.md)

# SPECULA Atmospheric-Field Benchmark Scope

## Purpose

This note records the current maintained scope for the `CP-05`
SPECULA-aligned atmospheric-field benchmark family.

The package already has strong contract-oriented validation for:

- atmospheric geometric propagation
- layered Fresnel propagation
- chromatic atmospheric propagation
- curvature-through-atmosphere behavior

through the committed SPECULA-targeted reference bundle and the local runtime
profiling surface.

What is still missing is a maintained external benchmark runner that can
participate in the cross-package harness the same way `CP-02` and `CP-03` now
do for the Julia-to-Julia REVOLT families.

## Available Assets

The following external assets exist today:

- SPECULA `AtmoPropagation` implementation in
  [../SPECULA/specula/processing_objects/atmo_propagation.py](../../SPECULA/specula/processing_objects/atmo_propagation.py)
- SPECULA atmospheric propagation tests under
  [../SPECULA/test/test_atmo_propagation.py](../../SPECULA/test/test_atmo_propagation.py)
- REVOLT/SPECULA configs that route atmospheric electric fields through
  `AtmoPropagation`, for example:
  - [../REVOLT/Python/specula/params_revolt_modal.yml](../../REVOLT/Python/specula/params_revolt_modal.yml)
  - [../REVOLT/Python/specula/PWFS/params_revolt_modal_PWFS.yml](../../REVOLT/Python/specula/PWFS/params_revolt_modal_PWFS.yml)

The following maintained local asset exists today:

- [profile_atmospheric_field_runtime.jl](../scripts/profile_atmospheric_field_runtime.jl)

## Why `CP-05` Is Deferred

`CP-05` is currently recorded as a scoped defer rather than a runnable harness
family.

Reason:

- the external trees do not currently provide a maintained executable runner
  that emits comparable runtime-profile fields for:
  - atmospheric geometric propagation
  - layered Fresnel propagation
  - curvature-through-atmosphere
- the available external assets are tests and YAML-configured pipelines, not a
  normalized benchmark entrypoint that fits the maintained harness contract

So `CP-05` is tracked in the harness as a contract-only scenario with an
explicit skip reason, not as an implicit missing feature.

## What Would Unblock It

`CP-05` should be promoted from deferred to runnable when at least one of the
following exists:

1. a maintained SPECULA-side runner that reports comparable runtime fields for
   an atmospheric-field propagation scenario
2. a small wrapper in the neighboring SPECULA/REVOLT tree that turns an
   existing `AtmoPropagation` test/config surface into a stable benchmark
   command
3. a frozen benchmark-style result generator that can be refreshed in a normal
   maintainer workflow and archived under the cross-package harness policy

## Current Interpretation

Until that happens:

- validation for atmospheric-field behavior remains primarily in:
  - [specula-reference-datasets.md](./specula-reference-datasets.md)
  - [model-validity-matrix.md](./model-validity-matrix.md)
- runtime evidence remains package-local in:
  - [profile_atmospheric_field_runtime.jl](../scripts/profile_atmospheric_field_runtime.jl)
- `CP-05` remains an explicit deferred benchmark family rather than an
  unstated gap

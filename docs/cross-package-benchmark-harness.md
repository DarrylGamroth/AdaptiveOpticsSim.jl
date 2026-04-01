# Cross-Package Benchmark Harness

Status: active

Plan traceability:

- [`PLAN-28`](./package-review-action-plan.md)
- [`PLAN-29`](./package-review-action-plan.md)
- [`PLAN-30`](./package-review-action-plan.md)
- [`PLAN-31`](./package-review-action-plan.md)
- [`PLAN-32`](./package-review-action-plan.md)
- review IDs: `PR-21`, `PR-22`, `PR-23`, `PR-24`, `PR-25`

## Purpose

This document defines the maintained cross-package benchmark harness for:

- AdaptiveOpticsSim main
- the REVOLT-aligned Julia branch in
  [`../AdaptiveOpticsSim.jl-revolt-real`](../AdaptiveOpticsSim.jl-revolt-real)
- external OOPAO/SPECULA comparison surfaces where frozen bundles or executable
  scenario assets exist

The harness is engineering evidence, not unit-test coverage.

## Contract Source

Scenario contracts are maintained in:

- [`benchmarks/contracts/cross_package.toml`](../benchmarks/contracts/cross_package.toml)

The current maintained ladder is:

- `cp01_compact_reference`
  - compact fidelity-first OOPAO/SPECULA frozen bundle comparison
- `cp02_revolt_sh_medium`
  - medium REVOLT-like Shack-Hartmann HIL runtime comparison
  - normalization rules and accepted differences are documented in
    [`revolt-sh-benchmark-contract.md`](./revolt-sh-benchmark-contract.md)
- `cp03_revolt_pwfs_representative`
  - representative REVOLT-like PWFS comparison between `main` and
    `revolt-real`
  - normalization rules and accepted differences are documented in
    [`revolt-pwfs-benchmark-contract.md`](./revolt-pwfs-benchmark-contract.md)
- `cp05_specula_atmo_field_medium`
  - SPECULA-aligned atmospheric-field benchmark family
  - currently contract-only/deferred with scope note in
    [`specula-atmo-field-benchmark-scope.md`](./specula-atmo-field-benchmark-scope.md)

## Runner

Run the maintained harness with:

```bash
julia --project=. --startup-file=no scripts/run_cross_package_benchmarks.jl
```

Optional arguments:

```bash
julia --project=. --startup-file=no scripts/run_cross_package_benchmarks.jl \
  benchmarks/contracts/cross_package.toml \
  benchmarks/results/cross_package/custom-run.toml
```

## Output Policy

Results are written under:

- [`benchmarks/results/cross_package`](../benchmarks/results/cross_package)

Maintained metadata is recorded in:

- [`benchmarks/results/cross_package/manifest.toml`](../benchmarks/results/cross_package/manifest.toml)

The result schema separates:

- scenario identity and class
- runtime metrics
- fidelity metrics
- deferred or skipped scenarios

## Manifest Policy

The archived cross-package manifest records two policy layers:

- `scenario_policies`
  - classifies each benchmark family as `mandatory`, `optional`, or
    `environment_gated`
- `implementation_policies`
  - records per-implementation exceptions when a scenario is maintained but a
    specific participant depends on neighboring-tree assets or optional runtime
    environment setup

Current interpretation:

- `mandatory`
  - expected to appear in the active archived ladder and pass in a normal
    maintainer environment
- `optional`
  - useful comparison participant, but not required for the harness run to be
    considered valid
- `environment_gated`
  - maintained by contract, but expected to remain deferred until the required
    external runnable asset or environment is available

This keeps the archived benchmark record explicit about what is required,
best-effort, and intentionally deferred.

## Maintained Metrics

### Runtime metrics

Recorded where the scenario provides a maintained profile surface:

- `build_time_ns`
- `total_mean_ns`
- `total_p95_ns`
- `frame_rate_hz`
- `total_alloc_bytes`

### Fidelity metrics

Recorded where a frozen reference bundle exists:

- `case_count`
- `failed_case_count`
- `max_abs_error`
- `max_rel_error`
- `l2_rel_error`

## Test Separation

This harness is intentionally outside [`Pkg.test()`](../test/runtests.jl).

Reasons:

- cross-package runs can depend on neighboring trees and non-core assets
- benchmark runtimes are materially slower than unit/regression tests
- representative comparisons are engineering evidence, not correctness gates

Unit tests remain responsible for:

- local functional correctness
- frozen bundle regression
- backend smoke and parity checks

The cross-package harness is responsible for:

- realistic comparative runtime evidence
- cross-package fidelity evidence where equivalent bundles exist
- archived benchmark records for maintainer review

# REVOLT Platform Benchmark Contract

Status: active

Plan traceability:

- [`PSP-16`](./platform-strengthening-plan.md)

## Purpose

This note records the maintained cross-package platform benchmark contract for
the typed orchestration layer.

The contract is intentionally narrow:

- compare the maintained `multi_source_multi_wfs` medium rung
- normalize the shared runtime-shape fields
- keep the historical comparison fields explicit while the external
  comparison-workspace grouped runner remains deferred

It does not claim full platform-output equivalence.

## Contract

- scenario id: `cp06_revolt_platform_medium`
- contract file:
  [cross_package.toml](../benchmarks/contracts/cross_package.toml)

Normalized comparison fields:

- `backend`
- `scale`
- `runtime_resolution`
- `runtime_n_subap`
- `runtime_branch_count`

Accepted known differences:

- `typed_platform_scenario_surface`
- `grouped_export_stack_contract`
- `composite_runtime_entrypoint_shape`

These differences remain recorded because `main` now routes this rung through
the typed `PlatformScenario` layer, while the former REVOLT-aligned fork was
never normalized into a distinct external grouped-platform runner before the
comparison-workspace migration.

## Maintained Command

The core-package contract is now intentionally `contract_only` for `CP-06`.
The historical archived result remains:

## Archived Result

- [2026-04-03-phase5-psp16.toml](../benchmarks/results/cross_package/2026-04-03-phase5-psp16.toml)
- manifest:
  [manifest.toml](../benchmarks/results/cross_package/manifest.toml)

## Interpretation

This contract answers a narrower question than the earlier REVOLT-like SH and
PWFS contracts:

- which grouped-runtime fields were normalized before the external
  comparison-workspace migration deferred this family?

It does not answer:

- whether grouped export products are identical across trees
- whether the internal composition APIs are equivalent
- whether broader SPECULA-style platform behavior is matched

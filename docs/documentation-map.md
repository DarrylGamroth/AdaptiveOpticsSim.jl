# Documentation Map

Status: active

Plan traceability:

- [`PLAN-33`](./package-review-action-plan.md)
- [`PLAN-37`](./package-review-action-plan.md)
- review IDs: `PR-26`, `PR-27`, `PR-29`

## Purpose

This is the maintained navigation surface for the documentation set.

Use this page first instead of scanning the full `docs/` directory. The goal is
to keep a small set of stable entry points while preserving the more detailed
planning and audit records for subsystem work.

## Start Here

If you are:

- learning the package:
  - [`user-guide.md`](./user-guide.md)
- looking for maintained public APIs:
  - [`api-reference.md`](./api-reference.md)
- maintaining or refactoring the codebase:
  - [`maintainer-architecture.md`](./maintainer-architecture.md)
- trying to understand runtime execution:
  - [`runtime-dataflow.md`](./runtime-dataflow.md)
- checking whether a model is validated:
  - [`model-validity-matrix.md`](./model-validity-matrix.md)
- looking for performance evidence:
  - [`benchmark-matrix-plan.md`](./benchmark-matrix-plan.md)
  - [`cross-package-benchmark-harness.md`](./cross-package-benchmark-harness.md)

## Active Stable Guides

- [`user-guide.md`](./user-guide.md)
  - workflow-oriented user entry point
- [`api-reference.md`](./api-reference.md)
  - maintained public API surface
- [`maintainer-architecture.md`](./maintainer-architecture.md)
  - current system structure and subsystem ownership
- [`runtime-dataflow.md`](./runtime-dataflow.md)
  - end-to-end build, step, and export dataflow
- [`model-validity-matrix.md`](./model-validity-matrix.md)
  - model claims, evidence, and limitations
- [`benchmark-matrix-plan.md`](./benchmark-matrix-plan.md)
  - maintained performance surfaces and benchmark classes
- [`cross-package-benchmark-harness.md`](./cross-package-benchmark-harness.md)
  - cross-package benchmark execution and archived evidence

## Active Supporting Reference

- [`roadmap.md`](./roadmap.md)
  - package-level implementation status and backlog
- [`units-policy.md`](./units-policy.md)
  - unit normalization rules
- [`deterministic-simulation.md`](./deterministic-simulation.md)
  - reproducibility guidance
- [`oopao-reference-datasets.md`](./oopao-reference-datasets.md)
  - OOPAO frozen bundle policy and provenance
- [`specula-reference-datasets.md`](./specula-reference-datasets.md)
  - SPECULA-targeted frozen bundle policy and provenance
- [`optional-integration-boundaries.md`](./optional-integration-boundaries.md)
  - optional integration and package-boundary rules

## Active Implementation Plans

These documents are still useful, but they are subsystem plans rather than
entry-point guides.

- [`algorithmic-implementation-roadmap.md`](./algorithmic-implementation-roadmap.md)
- [`atmospheric-field-propagation-roadmap.md`](./atmospheric-field-propagation-roadmap.md)
- [`core-optics-expansion-roadmap.md`](./core-optics-expansion-roadmap.md)
- [`gpu-runtime-structural-refactor-plan.md`](./gpu-runtime-structural-refactor-plan.md)
- [`infinite-boundary-atmosphere-plan.md`](./infinite-boundary-atmosphere-plan.md)

Use these when actively working in that subsystem. Do not treat them as the
primary orientation docs for new contributors.

## Completed Reviews and Inventories

These remain useful as implementation evidence, but they are no longer the
recommended entry point for normal use.

- [`package-review-2026-03.md`](./package-review-2026-03.md)
- [`package-review-action-plan.md`](./package-review-action-plan.md)
- [`api-tier-inventory.md`](./api-tier-inventory.md)
- [`modularization-inventory.md`](./modularization-inventory.md)
- [`reusable-infrastructure-inventory.md`](./reusable-infrastructure-inventory.md)
- [`model-validation-inventory.md`](./model-validation-inventory.md)
- [`cross-package-benchmark-inventory.md`](./cross-package-benchmark-inventory.md)

## Status Conventions

- `active`
  - maintained as a current guide or currently relevant implementation plan
- `completed`
  - retained as an implemented review, inventory, or plan artifact
- `archived/superseded`
  - no longer current; kept only for historical context

When adding new docs, prefer one of the stable guides above unless the document
is specifically a review, inventory, or subsystem plan.

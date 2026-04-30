# Documentation Map

Status: active

This is the maintained navigation surface for the documentation set. The docs
directory is intentionally small. Historical plans, audits, inventories, and
triage notes are not kept as live documentation; use git history if old
implementation records are needed.

## User Docs

- [`user-guide.md`](./user-guide.md)
- [`model-cookbook.md`](./model-cookbook.md)
- [`api-reference.md`](./api-reference.md)
- [`julia-tutorial-mappings.md`](./julia-tutorial-mappings.md)

## Maintainer Docs

- [`maintainer-architecture.md`](./maintainer-architecture.md)
- [`runtime-dataflow.md`](./runtime-dataflow.md)
- [`extension-guide.md`](./extension-guide.md)
- [`roadmap.md`](./roadmap.md)

## Validation Docs

- [`supported-production-surfaces.md`](./supported-production-surfaces.md)
- [`release-validation-runbook.md`](./release-validation-runbook.md)
- [`backend-validation-guide.md`](./backend-validation-guide.md)
- [`model-validity-matrix.md`](./model-validity-matrix.md)

## Policy Docs

- [`deterministic-simulation.md`](./deterministic-simulation.md)
- [`units-policy.md`](./units-policy.md)

## Archive Policy

The April 2026 cleanup removed completed plans and one-off audit records from
the live docs tree. [`archive/2026-04/README.md`](./archive/2026-04/README.md)
records that decision. Use git history for the deleted documents.

## Adding New Docs

- Prefer updating an existing document.
- Add a new top-level document only when it is a durable guide or runbook.
- Do not add one-off plans, audit notes, or temporary triage documents under
  `docs/`; use issues, PR descriptions, or git history for those records.
- Every top-level Markdown document must include `Status: active` near the top.

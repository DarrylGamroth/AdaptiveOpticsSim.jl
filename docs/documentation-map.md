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

### User Path

Most users should stop after these three:

- [`user-guide.md`](./user-guide.md)
- [`model-cookbook.md`](./model-cookbook.md)
- [`api-reference.md`](./api-reference.md)
- `examples/tutorials/`

### Validation And Production Path

Use these when the task is evidence, release support, or backend validation:

- [`model-validity-matrix.md`](./model-validity-matrix.md)
- [`supported-production-surfaces.md`](./supported-production-surfaces.md)
- [`production-readiness-checklist.md`](./production-readiness-checklist.md)
- [`release-validation-runbook.md`](./release-validation-runbook.md)
- [`self-hosted-gpu-runner-setup.md`](./self-hosted-gpu-runner-setup.md)
- [`benchmark-matrix-plan.md`](./benchmark-matrix-plan.md)
- [`cross-package-benchmark-harness.md`](./cross-package-benchmark-harness.md)

### Developer And Maintainer Path

Use these when changing package structure or internals:

- [`maintainer-architecture.md`](./maintainer-architecture.md)
- [`platform-architecture.md`](./platform-architecture.md)
- [`platform-workflows.md`](./platform-workflows.md)
- [`scenario-builder-style.md`](./scenario-builder-style.md)
- [`runtime-dataflow.md`](./runtime-dataflow.md)

### External Truth And Equivalence

Use these when the task is cross-package or instrument-truth work:

- [`oopao-reference-datasets.md`](./oopao-reference-datasets.md)
- [`../scripts/generate_oopao_equivalence_artifact.jl`](../scripts/generate_oopao_equivalence_artifact.jl)
- [`../scripts/generate_heart_boundary_truth_artifact.py`](../scripts/generate_heart_boundary_truth_artifact.py)
- [`../benchmarks/results/truth/2026-04-09-heart-boundary-truth.toml`](../benchmarks/results/truth/2026-04-09-heart-boundary-truth.toml)

## Active Stable Guides

- [`user-guide.md`](./user-guide.md)
  - workflow-oriented user entry point
- [`model-cookbook.md`](./model-cookbook.md)
  - compact recipe-first entry point for common model types, including HIL-style runtime surfaces
- [`platform-architecture.md`](./platform-architecture.md)
  - stable synthesis guide for platform-scale structure and boundaries
- [`platform-workflows.md`](./platform-workflows.md)
  - stable workflow guide for script-first simulation, validation, and benchmarking
- [`platform-orchestration.md`](./platform-orchestration.md)
  - typed scenario/config layer for Julia-native runtime composition
- [`platform-orchestration-validation.md`](./platform-orchestration-validation.md)
  - maintained benchmark and backend evidence for the `PlatformScenario` layer
- [`platform-strengthening-closeout.md`](./platform-strengthening-closeout.md)
  - explicit post-plan decision on what still trails SPECULA and what comes next
- [`platform-manifest-defer.md`](./platform-manifest-defer.md)
  - explicit defer note for optional scenario manifests/config-style orchestration
- [`external-comparison-workspace-plan.md`](./external-comparison-workspace-plan.md)
  - migration plan for replacing the `revolt-real` fork with a dedicated
    external comparison workspace
- [`scenario-builder-style.md`](./scenario-builder-style.md)
  - maintained conventions for script-first scenario construction
- [`api-cleanup-plan-2026-04.md`](./api-cleanup-plan-2026-04.md)
  - phased cleanup plan for runtime/HIL naming, accessors, and external-control semantics
- [`api-reference.md`](./api-reference.md)
  - maintained public API surface
- [`maintainer-architecture.md`](./maintainer-architecture.md)
  - current system structure and subsystem ownership
- [`runtime-dataflow.md`](./runtime-dataflow.md)
  - end-to-end build, step, and export dataflow
- [`model-validity-matrix.md`](./model-validity-matrix.md)
  - model claims, evidence, and limitations
- [`supported-production-surfaces.md`](./supported-production-surfaces.md)
  - explicit production-supported vs experimental scope
- [`production-readiness-checklist.md`](./production-readiness-checklist.md)
  - remaining blockers before production-readiness claims
- [`production-hardening-plan-2026-04.md`](./production-hardening-plan-2026-04.md)
  - ordered execution plan for finishing the current production-hardening pass
- [`operational-gpu-validation-cadence.md`](./operational-gpu-validation-cadence.md)
  - required real-hardware release-validation cadence for CUDA and AMDGPU
- [`production-boundary-freeze-2026-04.md`](./production-boundary-freeze-2026-04.md)
  - explicit decision to freeze the current supported production boundary for this release train
- [`release-validation-runbook.md`](./release-validation-runbook.md)
- [`self-hosted-gpu-runner-setup.md`](./self-hosted-gpu-runner-setup.md)
  - one-command release validation procedure and optional GPU/comparison tracks
- [`benchmark-matrix-plan.md`](./benchmark-matrix-plan.md)
  - maintained performance surfaces and benchmark classes
- [`cross-package-benchmark-harness.md`](./cross-package-benchmark-harness.md)
  - cross-package benchmark execution, archived evidence, and scenario policy
- [`revolt-sh-benchmark-contract.md`](./revolt-sh-benchmark-contract.md)
  - normalization rules and accepted differences for the maintained `CP-02`
    REVOLT-like Shack-Hartmann comparison
- [`revolt-pwfs-benchmark-contract.md`](./revolt-pwfs-benchmark-contract.md)
  - normalization rules and accepted differences for the maintained `CP-03`
    REVOLT-like Pyramid comparison
- [`revolt-platform-benchmark-contract.md`](./revolt-platform-benchmark-contract.md)
  - normalization rules and accepted differences for the maintained `CP-06`
    grouped platform-runtime comparison
- [`specula-atmo-field-benchmark-scope.md`](./specula-atmo-field-benchmark-scope.md)
  - current scope and defer conditions for the `CP-05`
    SPECULA-aligned atmospheric-field benchmark family
- [`future-platform-direction.md`](./future-platform-direction.md)
  - post-review platform direction and scope guardrails
- [`platform-strengthening-closeout.md`](./platform-strengthening-closeout.md)
  - Phase 6 closeout decision for the main platform-strengthening pass
- [`post-review-platform-plan.md`](./post-review-platform-plan.md)
  - phased next-step plan for validation hardening, benchmark expansion, and
    next capability selection
- [`execution-plan-milestones.md`](./execution-plan-milestones.md)
  - phased rollout for backend-specific execution-plan adoption
- [`execution-plan-closeout.md`](./execution-plan-closeout.md)
  - rollout outcome, builder decision, and follow-up triggers
- [`platform-strengthening-plan.md`](./platform-strengthening-plan.md)
  - next main-platform implementation plan after the execution-plan rollout
- [`next-milestone-decision-2026-04.md`](./next-milestone-decision-2026-04.md)
  - evidence-backed selection memo for the next breadth milestone
- [`grouped-runtime-plan.md`](./grouped-runtime-plan.md)
  - chosen next breadth milestone for grouped sensing/runtime orchestration
- [`grouped-runtime-contract.md`](./grouped-runtime-contract.md)
  - maintained grouped source/runtime/export contract
- [`grouped-runtime-validation.md`](./grouped-runtime-validation.md)
  - grouped runtime validation artifact and scope limits
- [`specula-platform-runtime-validation.md`](./specula-platform-runtime-validation.md)
  - Julia-native SPECULA-informed platform/runtime artifact and scope limits
- [`rocm-fallback-inventory.md`](./rocm-fallback-inventory.md)
  - explicit ROCm fallback classification after the execution-plan rollout and
    Phase 3 cleanup
- [`rocm-phase3-rebaseline.md`](./rocm-phase3-rebaseline.md)
  - CPU/AMDGPU/CUDA realistic runtime rebaseline after the Phase 3 ROCm cleanup

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
- [`phase-statistics-accuracy.md`](./phase-statistics-accuracy.md)
  - maintained accuracy note for the shared `K_{5/6}` covariance helper
- [`atmosphere-statistics-validation.md`](./atmosphere-statistics-validation.md)
  - maintained finite/infinite atmosphere statistics artifact and limitations
- [`detector-validation.md`](./detector-validation.md)
  - maintained detector-family fixture artifact and scope limits
- [`rocm-fallback-inventory.md`](./rocm-fallback-inventory.md)
  - explicit ROCm fallback classification and residual workaround status
- [`rocm-phase3-rebaseline.md`](./rocm-phase3-rebaseline.md)
  - maintained Phase 3 CPU/AMDGPU/CUDA runtime rebaseline
- [`gpu-sh-centroid-redesign-plan.md`](./gpu-sh-centroid-redesign-plan.md)
  - redesign plan for a genuinely GPU-friendly Shack-Hartmann centroid/export
    path shared across CUDA and AMDGPU
- [`gpu-readout-correction-plan.md`](./gpu-readout-correction-plan.md)
  - completed/active plan for stack-aware GPU readout-correction on maintained
    detector surfaces
- [`gpu-detector-followup-plan.md`](./gpu-detector-followup-plan.md)
  - follow-on plan for detector finalization, generalized batching, and
    GPU-native output surfaces after readout-correction recovery
- [`lift-gsc-runtime-validation.md`](./lift-gsc-runtime-validation.md)
  - maintained workflow-profile artifact for LiFT and gain-sensing camera
- [`tomography-benchmark-scope.md`](./tomography-benchmark-scope.md)
  - current maintained scope and archived re-scoping record for representative
    tomography benchmark evidence
- [`optional-integration-boundaries.md`](./optional-integration-boundaries.md)
  - optional integration and package-boundary rules
- [`future-platform-direction.md`](./future-platform-direction.md)
  - post-cleanup expansion policy
- [`backend-validation-guide.md`](./backend-validation-guide.md)
  - functional vs backend smoke vs benchmark separation
- [`execution-plan-architecture.md`](./execution-plan-architecture.md)
  - shared semantic contracts with backend-specific execution-plan policy
- [`execution-plan-closeout.md`](./execution-plan-closeout.md)
  - completed rollout summary and maintained decision record

## Active Implementation Plans

These documents are still useful, but they are subsystem plans rather than
entry-point guides.

- [`algorithmic-implementation-roadmap.md`](./algorithmic-implementation-roadmap.md)
- [`atmospheric-field-propagation-roadmap.md`](./atmospheric-field-propagation-roadmap.md)
- [`core-optics-expansion-roadmap.md`](./core-optics-expansion-roadmap.md)
- [`gpu-runtime-structural-refactor-plan.md`](./gpu-runtime-structural-refactor-plan.md)
- [`execution-plan-architecture.md`](./execution-plan-architecture.md)
- [`infinite-boundary-atmosphere-plan.md`](./infinite-boundary-atmosphere-plan.md)
- [`execution-plan-milestones.md`](./execution-plan-milestones.md)
- [`execution-plan-closeout.md`](./execution-plan-closeout.md)
- [`platform-strengthening-plan.md`](./platform-strengthening-plan.md)
- [`post-review-platform-plan.md`](./post-review-platform-plan.md)
- [`amdgpu-sh-convergence-plan.md`](./amdgpu-sh-convergence-plan.md)
- [`gpu-sh-centroid-redesign-plan.md`](./gpu-sh-centroid-redesign-plan.md)

Use these when actively working in that subsystem. Do not treat them as the
primary orientation docs for new contributors.

## Completed Reviews and Inventories

These remain useful as implementation evidence, but they are no longer the
recommended entry point for normal use.

- [`package-review-2026-03.md`](./package-review-2026-03.md)
- [`package-review-status-2026-04.md`](./package-review-status-2026-04.md)
  - current read of which March review findings still stand vs have been mitigated
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

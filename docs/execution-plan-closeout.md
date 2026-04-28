# Execution Plan Closeout

Status: completed

Plan traceability:

- [execution-plan-architecture.md](./execution-plan-architecture.md)
- [execution-plan-milestones.md](./execution-plan-milestones.md)

## Purpose

This note closes the first execution-plan rollout.

It records:

- what was implemented,
- what was intentionally left alone,
- what benchmark evidence backs the maintained defaults,
- and what should trigger a future follow-up.

## Milestone Outcome

Completed milestones:

- `EP-1` grouped WFS plan formalization
- `EP-2` AMDGPU grouped WFS recovery
- `EP-3` detector execution plans
- `EP-4` atmospheric-field execution plans
- `EP-5` reduction and runtime-export plans
- `EP-6` builder adoption and closeout

Implemented plan families now exist for:

- grouped WFS accumulation
- detector capture and readout staging
- atmospheric-field propagation and spectral accumulation
- reduction policy
- runtime export policy

## Builder Decision

Tomography and calibration builders were reviewed during `EP-6`.

The outcome is:

- no new execution-plan layer was added there in this rollout

Reason:

- the builder surfaces already have an explicit typed policy layer through
  [`BuildBackend`](../src/core/inverse_policies.jl),
- tomography already routes backend-sensitive matrix/materialization decisions
  through [`GPUArrayBuildBackend`](../src/core/inverse_policies.jl) and
  [`prepare_build_matrix`](../src/core/inverse_policies.jl),
- adding a second execution-plan abstraction there now would duplicate an
  existing mechanism without a demonstrated benchmark problem.

So the maintained rule is:

- keep using `BuildBackend` for builder-side host/device algebra policy
- only introduce a dedicated execution-plan layer for builders if a new,
  benchmark-backed backend-sensitive bottleneck appears that `BuildBackend`
  cannot express cleanly

## Benchmarked Defaults

The rollout is backed by maintained benchmark surfaces.

Grouped runtime:

- CPU `multi_source_multi_wfs medium`
- AMDGPU `multi_source_multi_wfs medium`
- CUDA `multi_source_multi_wfs medium`

Atmospheric field:

- CPU `profile_atmospheric_field_runtime.jl geometric|fresnel finite medium`
- CUDA `profile_atmospheric_field_runtime.jl geometric|fresnel finite medium`

Detector-backed runtime:

- CPU `profile_ao3k_runtime.jl cpu medium default`
- CUDA `profile_ao3k_runtime.jl cuda medium default`

These are the surfaces that currently justify:

- `GroupedStaged2DPlan` for AMDGPU grouped Pyramid/BioEdge
- direct detector execution plans on CPU/CUDA
- host-mirror detector execution plans on AMDGPU
- synchronous atmospheric-field plans on CPU
- async atmospheric-field plans on accelerators
- host-mirror reduction plans on AMDGPU
- direct reduction plans on CPU/CUDA

## Validation Summary

The rollout kept the following evidence green:

- full local `Pkg.test()`
- optional backend plan checks
- grouped runtime artifact generation
- CUDA smoke on `spiders`
- OOPAO and SPECULA frozen reference regressions

## Residual Risks

The remaining risks are narrower than when this rollout started:

- AMDGPU performance still depends on a few maintained host-mirror fallbacks
- runtime export plans are explicit now, but still intentionally simple
- builder-side execution plans are not a first-class family yet

These are acceptable for now because:

- they are explicit,
- they are benchmarked where maintained,
- and they no longer leak through ad hoc backend-name branching in hot paths.

## Follow-Up Triggers

Re-open execution-plan work only if one of these happens:

1. a maintained backend regresses on a benchmark-backed surface
2. a new backend such as Metal needs a different default memory plan
3. tomography/calibration builders show a backend-sensitive bottleneck that the
   current `BuildBackend` policy cannot express
4. runtime export needs more than the current direct/composite split

If none of those are true, new work should prefer feature or validation
progress over further execution-plan abstraction.

# Execution Plan Milestones

Status: active

Plan traceability:

- architecture: [execution-plan-architecture.md](./execution-plan-architecture.md)
- grouped follow-up: [grouped-runtime-plan.md](./grouped-runtime-plan.md)
- platform direction: [future-platform-direction.md](./future-platform-direction.md)

## Purpose

This document turns the execution-plan architecture into a short,
implementation-facing milestone sequence.

The immediate motivation is clear:

- grouped runtime semantics improved,
- but AMDGPU grouped diffractive performance regressed,
- and the package needs a traceable way to recover that performance without
  abandoning the shared grouped-runtime contract.

## Core Requirements

- `EPR-01`
  - execution-plan selection must be static and type-driven
- `EPR-02`
  - semantic contracts must remain shared across plans
- `EPR-03`
  - backend-specific defaults must be benchmark-backed
- `EPR-04`
  - semantic equivalence across plans must be tested
- `EPR-05`
  - plan selection must happen outside hot loops
- `EPR-06`
  - new backend-sensitive plan layers must start internal, not public

## Milestones

### EP-1: Grouped WFS Plan Formalization

Status: implemented

Goal:

- introduce an explicit typed execution-plan layer for grouped WFS
  accumulation without changing the current default behavior

In scope:

- [grouped.jl](../src/wfs/grouped.jl)
- grouped Pyramid/BioEdge call sites
- grouped runtime tests

Deliverables:

- typed grouped accumulation plan objects
- typed grouped plan-selection entry point
- grouped accumulation routed through plan execution methods
- tests covering current default-plan selection and grouped behavior

Acceptance:

- current grouped CPU semantics remain unchanged
- the default grouped path still passes the full test suite
- the grouped runtime artifact remains reproducible
- the new plan-selection surface is internal and type-driven

Verification:

- `julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'`
- `julia --project=. --startup-file=no scripts/generate_grouped_runtime_artifact.jl`

Traceability:

- `EPR-01`, `EPR-02`, `EPR-05`, `EPR-06`

### EP-2: AMDGPU Grouped Recovery

Status: implemented

Goal:

- recover AMDGPU grouped Pyramid/BioEdge performance without regressing CUDA or
  CPU grouped semantics

In scope:

- grouped Pyramid/BioEdge accumulation defaults
- benchmark evidence for CPU, AMDGPU, and CUDA grouped surfaces

Deliverables:

- AMDGPU-specific grouped accumulation default for the affected WFS families
- benchmark evidence showing recovered AMDGPU grouped runtime
- no grouped-runtime correctness regression

Acceptance:

- AMDGPU `multi_source_multi_wfs medium` Pyramid/BioEdge grouped surfaces are
  no worse than the pre-regression maintained baseline within normal warmed
  benchmark variance
- CUDA grouped surfaces remain in-family
- CPU grouped surfaces remain in-family

Verification:

- `julia --project=. --startup-file=no scripts/profile_multi_source_multi_wfs_runtime.jl cpu medium`
- `julia --project=. --startup-file=no scripts/profile_multi_source_multi_wfs_runtime.jl amdgpu medium`
- on `spiders`:
  - `julia --project=. --startup-file=no scripts/profile_multi_source_multi_wfs_runtime.jl cuda medium`

Traceability:

- `EPR-02`, `EPR-03`, `EPR-04`

### EP-3: Detector Execution Plans

Status: implemented

Goal:

- separate detector semantic contracts from backend-specific staging and
  readout plans

Deliverables:

- typed detector execution plans for staged capture/readout
- ROCm/CUDA/CPU policy localized behind plan dispatch

Verification:

- `julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'`

Traceability:

- `EPR-01` through `EPR-06`

### EP-4: Atmospheric Field Execution Plans

Status: implemented

Goal:

- formalize execution-plan selection for atmospheric field propagation and
  spectral accumulation

Deliverables:

- typed atmospheric-field execution plans
- benchmark-backed default-plan choices for maintained backends

Acceptance:

- CPU atmospheric-field semantics remain unchanged
- maintained accelerator semantics remain unchanged
- the atmospheric-field benchmark surfaces remain in-family for CPU and CUDA

Verification:

- `julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'`
- `julia --project=. --startup-file=no scripts/profile_atmospheric_field_runtime.jl cpu geometric finite medium`
- `julia --project=. --startup-file=no scripts/profile_atmospheric_field_runtime.jl cpu fresnel finite medium`
- on `spiders`:
  - `julia --project=. --startup-file=no scripts/profile_atmospheric_field_runtime.jl cuda geometric finite medium`
  - `julia --project=. --startup-file=no scripts/profile_atmospheric_field_runtime.jl cuda fresnel finite medium`

Traceability:

- `EPR-01` through `EPR-06`

### EP-5: Reduction And Runtime Export Plans

Status: implemented

Goal:

- move remaining backend-sensitive reduction/export policy behind typed plan
  dispatch

Deliverables:

- reduction-plan and runtime-export-plan surfaces
- benchmark-backed defaults where a maintained backend difference exists

Acceptance:

- reduction semantics remain unchanged across maintained backends
- runtime export semantics remain unchanged for single and composite interfaces
- grouped-runtime benchmark surfaces remain in-family for CPU, AMDGPU, and CUDA

Verification:

- `julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'`
- `julia --project=. --startup-file=no scripts/profile_multi_source_multi_wfs_runtime.jl cpu medium`
- `julia --project=. --startup-file=no scripts/profile_multi_source_multi_wfs_runtime.jl amdgpu medium`
- on `spiders`:
  - `julia --project=. --startup-file=no scripts/profile_multi_source_multi_wfs_runtime.jl cuda medium`

Traceability:

- `EPR-01` through `EPR-06`

### EP-6: Builder Adoption And Closeout

Status: implemented

Goal:

- adopt the same architecture in tomography/calibration builders where it
  materially helps and close the benchmark/validation loop

Deliverables:

- builder-plan note or implementation where warranted
- updated benchmark/validation documentation
- explicit closeout note for the execution-plan rollout

Outcome:

- builder-side adoption is closed as a documented decision, not a new code path
- tomography/calibration builders continue to use the existing typed
  `BuildBackend` policy
- rollout closeout is recorded in
  [execution-plan-closeout.md](./execution-plan-closeout.md)

Verification:

- reviewed current builder policy in
  [inverse_policies.jl](../src/core/inverse_policies.jl)
- confirmed rollout benchmark and validation evidence in the milestone notes
  above
  

Traceability:

- `EPR-01` through `EPR-06`

## Checklist

- [x] `EP-1` grouped WFS plan formalization
- [x] `EP-2` AMDGPU grouped recovery
- [x] `EP-3` detector execution plans
- [x] `EP-4` atmospheric field execution plans
- [x] `EP-5` reduction and runtime export plans
- [x] `EP-6` builder adoption and closeout

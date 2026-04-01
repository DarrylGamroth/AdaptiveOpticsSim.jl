# Model Validity Matrix

Date: 2026-03-31

Status: active

Plan traceability:

- [`PLAN-23`](./package-review-action-plan.md)
- [`PLAN-24`](./package-review-action-plan.md)
- [`PLAN-25`](./package-review-action-plan.md)
- [`PLAN-26`](./package-review-action-plan.md)
- [`PLAN-27`](./package-review-action-plan.md)
- review IDs: `PR-15`, `PR-16`, `PR-17`, `PR-18`, `PR-19`, `PR-20`, `PR-36`

## Purpose

This document is the maintained validation surface for major model families in
AdaptiveOpticsSim.jl.

It answers four questions explicitly:

1. what model family is being validated
2. what class of evidence exists for it
3. what the current primary external baseline is
4. what assumptions, limits, or known non-equivalences still apply

This document supersedes the inventory-only role of
[model-validation-inventory.md](./model-validation-inventory.md). The inventory
remains useful as a planning/support artifact; this file is the maintained
matrix.

## Validation Classes

- `A`: analytic / structural check
  - invariants, conservation checks, deterministic replay, sanity identities,
    shape/range contracts
- `R`: frozen reference-bundle check
  - versioned committed bundles compared within tolerance
- `G`: backend parity / GPU contract
  - CPU, AMDGPU, and CUDA functional parity or smoke surfaces
- `P`: benchmark / runtime evidence
  - maintained runtime profiles or benchmark scripts on realistic cases
- `M`: maintained model note
  - spec, limitations doc, or explicit non-equivalence record

Validation classes are intentionally distinct:

- `A` answers “is the model internally coherent?”
- `R` answers “does it match a frozen external or committed baseline?”
- `G` answers “does the backend implementation behave consistently?”
- `P` answers “does it behave plausibly on realistic maintained workloads?”
- `M` answers “what is intentionally not claimed?”

## Frozen Bundle Roots

- OOPAO and pyTomoAO bundles:
  - [test/reference_data](../test/reference_data)
  - bundle doc: [oopao-reference-datasets.md](./oopao-reference-datasets.md)
- SPECULA-targeted contract bundle:
  - [test/reference_data_specula](../test/reference_data_specula)
  - bundle doc: [specula-reference-datasets.md](./specula-reference-datasets.md)

## Generation / Refresh Surfaces

- OOPAO:
  - [generate_oopao_reference_bundle.py](../scripts/generate_oopao_reference_bundle.py)
- pyTomoAO:
  - [generate_pytomoao_reference_bundle.py](../scripts/generate_pytomoao_reference_bundle.py)
- SPECULA-targeted contract bundle:
  - [generate_specula_reference_bundle.jl](../scripts/generate_specula_reference_bundle.jl)

## Matrix

| Family ID | Model family | Evidence | Primary baseline | Evidence links | Assumptions / limits | Current status |
| --- | --- | --- | --- | --- | --- | --- |
| `MV-01` | Atmosphere: finite and infinite multilayer propagation | `A`, `G`, `P`, `M` | analytic + backend parity | [runtests.jl](../test/runtests.jl), [gpu_smoke_contract.jl](../scripts/gpu_smoke_contract.jl), [profile_atmosphere_runtime.jl](../scripts/profile_atmosphere_runtime.jl), [atmosphere-runtime-spec.md](./atmosphere-runtime-spec.md) | No frozen external infinite-atmosphere statistics bundle yet; finite and infinite are validated mainly by internal statistics and backend parity | strong |
| `MV-02` | Phase statistics and covariance helpers | `A`, `G`, `M` | analytic | [runtests.jl](../test/runtests.jl), [gpu_smoke_contract.jl](../scripts/gpu_smoke_contract.jl), [atmosphere-runtime-spec.md](./atmosphere-runtime-spec.md) | `K_{5/6}` helper accuracy is exercised in tests but not yet summarized as a standalone maintained report | medium |
| `MV-03` | Core optics: electric field, Fraunhofer, Fresnel, atmospheric field propagation | `A`, `G`, `P`, `M` | analytic + SPECULA-informed design | [runtests.jl](../test/runtests.jl), [gpu_smoke_contract.jl](../scripts/gpu_smoke_contract.jl), [profile_atmospheric_field_runtime.jl](../scripts/profile_atmospheric_field_runtime.jl), [atmospheric-field-propagation-roadmap.md](./atmospheric-field-propagation-roadmap.md) | No frozen external atmospheric-field bundle yet; current external SPECULA alignment is via contract and design provenance rather than full-array equivalence | medium-strong |
| `MV-04` | Detectors and detector-family execution | `A`, `G`, `P`, `M` | analytic + runtime behavior | [runtests.jl](../test/runtests.jl), [optional_amdgpu_backends.jl](../test/optional_amdgpu_backends.jl), [optional_cuda_backends.jl](../test/optional_cuda_backends.jl), [profile_ao3k_runtime.jl](../scripts/profile_ao3k_runtime.jl) | Detector realism is strongest in integrated runtime scenarios; family-specific frozen detector references remain limited | medium |
| `MV-05` | Shack-Hartmann WFS | `A`, `R`, `G`, `P`, `M` | OOPAO | [reference_harness.jl](../test/reference_harness.jl), [reference_data](../test/reference_data), [gpu_smoke_contract.jl](../scripts/gpu_smoke_contract.jl), [profile_multi_source_multi_wfs_runtime.jl](../scripts/profile_multi_source_multi_wfs_runtime.jl), [profile_ao3k_runtime.jl](../scripts/profile_ao3k_runtime.jl) | OOPAO remains the primary frozen parity baseline; no SPECULA-targeted frozen SH bundle yet | strong |
| `MV-06` | Pyramid and BioEdge WFS | `A`, `R`, `G`, `P`, `M` | OOPAO | [reference_harness.jl](../test/reference_harness.jl), [reference_data](../test/reference_data), [gpu_smoke_contract.jl](../scripts/gpu_smoke_contract.jl), [profile_multi_source_multi_wfs_runtime.jl](../scripts/profile_multi_source_multi_wfs_runtime.jl) | OOPAO remains the primary parity baseline; grouped/polychromatic SPECULA-targeted bundles are still future work | strong |
| `MV-07` | Curvature and Zernike WFS | `A`, `R`, `G`, `P`, `M` | SPECULA-targeted contract bundle | [reference_harness.jl](../test/reference_harness.jl), [reference_data_specula](../test/reference_data_specula), [specula-reference-datasets.md](./specula-reference-datasets.md), [gpu_smoke_contract.jl](../scripts/gpu_smoke_contract.jl), [profile_atmospheric_field_runtime.jl](../scripts/profile_atmospheric_field_runtime.jl) | Current frozen coverage is contract-oriented and scenario-aligned to SPECULA tests; it is not a claim of complete numerical equivalence to the full SPECULA platform | strong |
| `MV-08` | LiFT and gain-sensing camera | `A`, `R`, `G`, `M` | OOPAO | [reference_harness.jl](../test/reference_harness.jl), [reference_data](../test/reference_data), [runtests.jl](../test/runtests.jl) | Runtime/benchmark evidence is lighter than for core WFS/runtime surfaces | medium-strong |
| `MV-09` | Runtime and closed-loop execution | `A`, `R`, `G`, `P`, `M` | OOPAO + maintained runtime profiles | [reference_harness.jl](../test/reference_harness.jl), [reference_data](../test/reference_data), [profile_ao3k_runtime.jl](../scripts/profile_ao3k_runtime.jl), [profile_revolt_hil_runtime.jl](../scripts/profile_revolt_hil_runtime.jl), [gpu_runtime_equivalence_contract.jl](../scripts/gpu_runtime_equivalence_contract.jl) | Cross-package benchmark harness is still Phase 5 work; current evidence is package-local plus OOPAO traces | strong |
| `MV-10` | Tomography and reconstruction | `A`, `R`, `G`, `M` | pyTomoAO / OOPAO-adjacent frozen references | [tomography.jl](../test/tomography.jl), [reference_harness.jl](../test/reference_harness.jl), [reference_data](../test/reference_data), tomography GPU profile scripts under [scripts](../scripts) | Representative benchmark evidence is still thin relative to closed-loop runtime surfaces | medium-strong |
| `MV-11` | GPU backend execution policy | `G`, `P`, `M` | backend parity | [optional_amdgpu_backends.jl](../test/optional_amdgpu_backends.jl), [optional_cuda_backends.jl](../test/optional_cuda_backends.jl), [backend-validation-guide.md](./backend-validation-guide.md), [gpu_smoke_contract.jl](../scripts/gpu_smoke_contract.jl), [profile_ao3k_runtime.jl](../scripts/profile_ao3k_runtime.jl), [profile_multi_source_multi_wfs_runtime.jl](../scripts/profile_multi_source_multi_wfs_runtime.jl), [gpu-runtime-structural-refactor-plan.md](./gpu-runtime-structural-refactor-plan.md) | Backend evidence is now split more clearly across functional tests, optional smoke, and benchmarks, but benchmark breadth still remains stronger on some runtime families than others | strong |
| `MV-12` | Tutorials and workflow examples | `A`, `M` | local workflow correctness | [examples/tutorials](../examples/tutorials), [reference_and_tutorials.jl](../test/testsets/reference_and_tutorials.jl), [user-guide.md](./user-guide.md) | Tutorials are executed, but they are not external parity evidence by themselves | medium |

## Known Baseline Rules

- OOPAO is the primary frozen parity baseline for:
  - PSF
  - Shack-Hartmann
  - Pyramid
  - BioEdge
  - LiFT
  - gain-sensing camera
  - compact closed-loop traces
- pyTomoAO remains the frozen baseline for tomography model/reconstructor data.
- SPECULA is the targeted external baseline where it is the stronger reference
  for behavior breadth rather than legacy parity, currently captured through the
  committed contract bundle for:
  - Zernike WFS
  - Curvature WFS

## Known Non-Equivalences and Scope Limits

- Infinite-atmosphere fidelity is not currently claimed against a frozen OOPAO
  or SPECULA statistics bundle.
- Atmospheric field propagation is validated by internal regression, backend
  parity, and benchmark evidence, but not yet by a frozen external field bundle.
- SPECULA-targeted contract bundles are intentionally narrower than full
  platform equivalence:
  - they freeze deterministic contract scenarios
  - they do not claim identical internal representation or orchestration
- Tutorial execution is evidence that workflows stay runnable, not evidence of
  scientific parity by itself.

## Acceptance Rule

For a model family to be treated as maintained and strong, the preferred target
is:

- `A` plus `R` or an explicitly justified absence of `R`
- `G` where accelerator support is claimed
- `P` where runtime/HIL or realistic execution claims are made
- `M` whenever equivalence is intentionally partial or scoped

Families that do not yet meet that shape should remain explicitly marked
medium/partial instead of being implied complete.

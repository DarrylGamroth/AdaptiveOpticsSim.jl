# Model Validation Inventory

Date: 2026-03-31

Status: supporting inventory

Plan traceability:

- [`PLAN-04`](./package-review-action-plan.md)
- review IDs: `PR-15`, `PR-16`, `PR-17`

## Purpose

This document inventories the maintained model families in
AdaptiveOpticsSim.jl and records their current validation evidence.

This document now supports the maintained
[model-validity-matrix.md](./model-validity-matrix.md).

The goal is not to redefine the final validation matrix here. The goal is to:

- normalize the evidence classes used across the package
- map each maintained model family to its current evidence
- identify which baseline is primary today:
  - analytic
  - OOPAO
  - SPECULA
  - backend parity
  - benchmark evidence
- record the most important gaps before Phase 4 implementation begins

## Evidence Classes

- `A`: analytic / structural check
  - invariants, shape/range rules, deterministic replay, sanity identities
- `R`: frozen reference-bundle check
  - committed or generated reference bundle compared within tolerance
- `G`: backend parity / GPU contract
  - CPU vs AMDGPU vs CUDA smoke or parity coverage
- `P`: benchmark / runtime evidence
  - maintained profile or benchmark script showing behavior on realistic cases
- `M`: manual or document-backed evidence
  - roadmap/spec statement, audit procedure, or explicit limitation note

## Maintained Model Families

| Family ID | Model family | Primary code surface | Current evidence | Primary baseline today | Current status | Main gaps |
| --- | --- | --- | --- | --- | --- | --- |
| `MV-01` | Atmosphere: finite and infinite multilayer propagation | [`src/Atmosphere/multilayer.jl`](../src/Atmosphere/multilayer.jl), [`src/Atmosphere/infinite_screen.jl`](../src/Atmosphere/infinite_screen.jl) | `A`, `G`, `P`, `M` | analytic + backend parity | strong | no frozen OOPAO/SPECULA bundle for infinite atmosphere statistics |
| `MV-02` | Phase statistics and covariance helpers | [`src/Atmosphere/phase_stats.jl`](../src/Atmosphere/phase_stats.jl), [`src/Core/kv56.jl`](../src/Core/kv56.jl) | `A`, `G`, `M` | analytic | medium | no dedicated benchmark/accuracy report beyond tests and smoke |
| `MV-03` | Core optics: electric field, Fraunhofer, Fresnel, atmospheric field propagation | [`src/Optics/electric_field.jl`](../src/Optics/electric_field.jl), [`src/Optics/propagation.jl`](../src/Optics/propagation.jl), [`src/Optics/atmospheric_field_propagation.jl`](../src/Optics/atmospheric_field_propagation.jl) | `A`, `G`, `P`, `M` | analytic + SPECULA-informed design | medium-strong | no frozen SPECULA comparison bundle yet |
| `MV-04` | Detectors and detector-family execution | [`src/Detectors/`](../src/Detectors) | `A`, `G`, `P`, `M` | analytic + runtime behavior | medium | limited frozen detector reference datasets; realism mostly exercised through runtime scenarios |
| `MV-05` | Shack-Hartmann WFS | [`src/WFS/shack_hartmann.jl`](../src/WFS/shack_hartmann.jl), [`src/WFS/subapertures.jl`](../src/WFS/subapertures.jl) | `A`, `R`, `G`, `P`, `M` | OOPAO | strong | no separate maintained cross-package benchmark report yet |
| `MV-06` | Pyramid and BioEdge WFS | [`src/WFS/pyramid.jl`](../src/WFS/pyramid.jl), [`src/WFS/bioedge.jl`](../src/WFS/bioedge.jl) | `A`, `R`, `G`, `P`, `M` | OOPAO | strong | SPECULA-targeted parity missing for grouped/polychromatic cases |
| `MV-07` | Curvature and Zernike WFS | [`src/WFS/curvature.jl`](../src/WFS/curvature.jl), [`src/WFS/zernike.jl`](../src/WFS/zernike.jl) | `A`, `G`, `P`, `M` | analytic + local regression | medium | no frozen external reference bundles yet |
| `MV-08` | LiFT and gain-sensing camera | [`src/Calibration/lift.jl`](../src/Calibration/lift.jl), gain-sensing surfaces in [`src/`](../src) | `A`, `R`, `G`, `M` | OOPAO | medium-strong | limited runtime/benchmark evidence outside compact regression cases |
| `MV-09` | Runtime and closed-loop execution | [`src/Control/runtime.jl`](../src/Control/runtime.jl) | `A`, `R`, `G`, `P`, `M` | OOPAO + runtime benchmarks | strong | cross-package benchmark harness not yet formalized |
| `MV-10` | Tomography and reconstruction | [`src/Tomography/`](../src/Tomography), [`test/tomography.jl`](../test/tomography.jl) | `A`, `R`, `G`, `M` | pyTomoAO / OOPAO-adjacent frozen references | medium-strong | benchmark evidence and representative runtime coverage are thin |
| `MV-11` | GPU backend execution policy | [`test/optional_gpu_backends.jl`](../test/optional_gpu_backends.jl), [`scripts/gpu_smoke_contract.jl`](../scripts/gpu_smoke_contract.jl), backend extensions in [`ext/`](../ext) | `G`, `P`, `M` | backend parity | strong | still lacks a single synthesized backend validation guide |
| `MV-12` | Tutorials and workflow examples | [`examples/tutorials/`](../examples/tutorials), [`test/runtests.jl`](../test/runtests.jl) | `A`, `M` | local workflow correctness | medium | examples are exercised, but not organized as a user-facing validity story |

## Evidence Inventory by Family

### `MV-01`: Atmosphere

Primary evidence:

- analytic/statistical regression in [`test/runtests.jl`](../test/runtests.jl)
  under `"Atmosphere propagation"`
- maintained atmosphere benchmark in
  [`scripts/profile_atmosphere_runtime.jl`](../scripts/profile_atmosphere_runtime.jl)
- GPU contract coverage in
  [`scripts/gpu_smoke_contract.jl`](../scripts/gpu_smoke_contract.jl)
- normative behavior/spec in
  [`docs/atmosphere-runtime-spec.md`](./atmosphere-runtime-spec.md)

Current gaps:

- no frozen OOPAO or SPECULA bundle for infinite-screen statistical behavior
- benchmark evidence exists, but not yet folded into a formal validity matrix

### `MV-02`: Phase statistics and covariance helpers

Primary evidence:

- direct regression coverage in [`test/runtests.jl`](../test/runtests.jl)
- GPU smoke coverage through atmosphere/tomography consumers in
  [`scripts/gpu_smoke_contract.jl`](../scripts/gpu_smoke_contract.jl)

Current gaps:

- helper-level accuracy reporting is implicit, not documented as a maintained
  validation artifact
- no standalone benchmark or error-summary doc for the `K_{5/6}` approximation

### `MV-03`: Core optics and atmospheric field propagation

Primary evidence:

- regression coverage in [`test/runtests.jl`](../test/runtests.jl)
- maintained field benchmark in
  [`scripts/profile_atmospheric_field_runtime.jl`](../scripts/profile_atmospheric_field_runtime.jl)
- GPU smoke coverage in
  [`scripts/gpu_smoke_contract.jl`](../scripts/gpu_smoke_contract.jl)
- design/acceptance statements in
  [`docs/atmospheric-field-propagation-roadmap.md`](./atmospheric-field-propagation-roadmap.md)

Current gaps:

- SPECULA is the intended primary external reference, but there is no frozen
  SPECULA bundle yet
- chromatic/extended-source field propagation still relies mainly on internal
  regression and smoke evidence

### `MV-04`: Detectors

Primary evidence:

- detector regression set in [`test/runtests.jl`](../test/runtests.jl)
- realistic runtime exercise via
  [`scripts/profile_ao3k_runtime.jl`](../scripts/profile_ao3k_runtime.jl)
- backend exercise in [`test/optional_gpu_backends.jl`](../test/optional_gpu_backends.jl)

Current gaps:

- detector-family realism is validated mostly through integrated runtime paths,
  not family-specific frozen references
- readout/noise model fidelity vs external references is not yet documented as
  a matrix

### `MV-05`: Shack-Hartmann

Primary evidence:

- frozen OOPAO reference harness in
  [`test/reference_harness.jl`](../test/reference_harness.jl)
- bundle provenance and case inventory in
  [`docs/oopao-reference-datasets.md`](./oopao-reference-datasets.md)
- runtime and GPU coverage in:
  - [`scripts/profile_multi_source_multi_wfs_runtime.jl`](../scripts/profile_multi_source_multi_wfs_runtime.jl)
  - [`scripts/profile_ao3k_runtime.jl`](../scripts/profile_ao3k_runtime.jl)
  - [`scripts/gpu_smoke_contract.jl`](../scripts/gpu_smoke_contract.jl)

Current gaps:

- no separate report comparing SH fidelity/runtime against OOPAO and SPECULA on
  the same realistic scenario ladder

### `MV-06`: Pyramid and BioEdge

Primary evidence:

- frozen OOPAO diffractive reference cases via
  [`test/reference_harness.jl`](../test/reference_harness.jl)
- grouped/GPU/runtime coverage in:
  - [`test/runtests.jl`](../test/runtests.jl)
  - [`scripts/profile_multi_source_multi_wfs_runtime.jl`](../scripts/profile_multi_source_multi_wfs_runtime.jl)
  - [`scripts/profile_ao3k_runtime.jl`](../scripts/profile_ao3k_runtime.jl)
  - [`scripts/gpu_smoke_contract.jl`](../scripts/gpu_smoke_contract.jl)

Current gaps:

- no SPECULA-oriented frozen comparisons for polychromatic and extended-source
  diffractive paths

### `MV-07`: Curvature and Zernike WFS

Primary evidence:

- regression coverage in [`test/runtests.jl`](../test/runtests.jl)
- atmospheric-field integrated runtime evidence in
  [`scripts/profile_atmospheric_field_runtime.jl`](../scripts/profile_atmospheric_field_runtime.jl)
- GPU smoke coverage in [`scripts/gpu_smoke_contract.jl`](../scripts/gpu_smoke_contract.jl)

Current gaps:

- no external frozen reference bundle
- current evidence is mostly internal-consistency and backend-parity based

### `MV-08`: LiFT and gain-sensing camera

Primary evidence:

- frozen OOPAO reference cases listed in
  [`docs/oopao-reference-datasets.md`](./oopao-reference-datasets.md)
- regression coverage in [`test/runtests.jl`](../test/runtests.jl)

Current gaps:

- performance evidence is thin outside compact regression surfaces
- long-horizon nonlinear GSC behavior is documented, but only partially reduced
  to frozen branch-point cases

### `MV-09`: Runtime and closed-loop execution

Primary evidence:

- frozen OOPAO closed-loop traces via
  [`test/reference_harness.jl`](../test/reference_harness.jl)
- maintained runtime benchmarks in:
  - [`scripts/profile_ao3k_runtime.jl`](../scripts/profile_ao3k_runtime.jl)
  - [`scripts/profile_revolt_hil_runtime.jl`](../scripts/profile_revolt_hil_runtime.jl)
  - [`scripts/profile_external_optics_hil.jl`](../scripts/profile_external_optics_hil.jl)
- GPU smoke/runtime equivalence scripts:
  - [`scripts/gpu_smoke_contract.jl`](../scripts/gpu_smoke_contract.jl)
  - [`scripts/gpu_runtime_equivalence_contract.jl`](../scripts/gpu_runtime_equivalence_contract.jl)

Current gaps:

- no single maintained doc yet mapping runtime models to exact parity and
  benchmark expectations
- cross-package benchmark methodology is not yet frozen

### `MV-10`: Tomography and reconstruction

Primary evidence:

- dedicated tomography regression in
  [`test/tomography.jl`](../test/tomography.jl)
- frozen pyTomoAO-style bundle coverage described in
  [`docs/oopao-reference-datasets.md`](./oopao-reference-datasets.md)
- GPU profiling/contracts in the tomography GPU scripts under
  [`scripts/`](../scripts)

Current gaps:

- representative benchmark evidence is not yet consolidated with the main
  benchmark ladder
- tomography validity story is split across tests and scripts rather than one
  matrix

### `MV-11`: GPU backend execution policy

Primary evidence:

- optional AMDGPU smoke in [`test/optional_gpu_backends.jl`](../test/optional_gpu_backends.jl)
- CUDA/AMDGPU smoke matrix in
  [`scripts/gpu_smoke_contract.jl`](../scripts/gpu_smoke_contract.jl)
- realistic backend benchmarks in:
  - [`scripts/profile_ao3k_runtime.jl`](../scripts/profile_ao3k_runtime.jl)
  - [`scripts/profile_multi_source_multi_wfs_runtime.jl`](../scripts/profile_multi_source_multi_wfs_runtime.jl)

Current gaps:

- backend validation is real, but the evidence is spread across tests, scripts,
  and planning docs
- there is no single “backend support and evidence” guide yet

### `MV-12`: Tutorials and workflow examples

Primary evidence:

- tutorial execution coverage in [`test/runtests.jl`](../test/runtests.jl)
- tutorial scripts in [`examples/tutorials/`](../examples/tutorials)

Current gaps:

- examples are not yet tied into the same validity taxonomy as model families
- user-facing guidance on “which examples are scientifically validated vs
  illustrative” is weak

## External Baseline Direction

Current recommendation:

- OOPAO remains the primary frozen parity baseline for:
  - SH
  - Pyramid
  - BioEdge
  - LiFT
  - compact closed-loop traces
- SPECULA should be the targeted primary external baseline for:
  - atmospheric field propagation
  - selected polychromatic and extended-source sensing cases
  - future broader platform/process comparisons

## Highest-Value Gaps Before Phase 4

1. Create a formal model-validity matrix doc derived from this inventory rather
   than leaving evidence spread across plans and tests.
2. Add frozen SPECULA comparison bundles for atmospheric field propagation and
   selected diffractive sensing cases.
3. Consolidate backend validation into one maintainer-facing guide.
4. Consolidate realistic runtime and cross-package benchmark methodology into a
   maintained comparison artifact.

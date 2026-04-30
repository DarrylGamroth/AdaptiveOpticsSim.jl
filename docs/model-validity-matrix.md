# Model Validity Matrix

Status: active

## Purpose

This document is the maintained validation surface for major model families in
AdaptiveOpticsSim.jl.

It answers four questions explicitly:

1. what model family is being validated
2. what class of evidence exists for it
3. what the current primary external baseline is
4. what assumptions, limits, or known non-equivalences still apply

This file is the maintained matrix. Historical inventories and planning notes
are intentionally not part of the live docs set.

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
  - bundle doc: committed OOPAO reference data under `test/reference_data`
- SPECULA-targeted contract bundle:
  - [test/reference_data_specula](../test/reference_data_specula)
  - bundle doc: committed SPECULA reference data under `test/reference_data_specula`

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
| `MV-01` | Atmosphere: finite and infinite multilayer propagation | `A`, `G`, `P`, `M` | analytic + backend parity | [runtests.jl](../test/runtests.jl), [gpu_smoke_contract.jl](../scripts/gpu_smoke_contract.jl), [profile_atmosphere_runtime.jl](../scripts/profile_atmosphere_runtime.jl), atmosphere runtime tests, atmosphere statistics artifact under `benchmarks/results/atmosphere/`, [2026-04-01-phase1-pvp02.toml](../benchmarks/results/atmosphere/2026-04-01-phase1-pvp02.toml) | No frozen external infinite-atmosphere parity bundle yet, but a committed fixed-seed finite/infinite statistics artifact now records variance agreement, stationarity, and non-periodicity evidence | strong |
| `MV-02` | Phase statistics and covariance helpers | `A`, `G`, `M` | analytic | [runtests.jl](../test/runtests.jl), [gpu_smoke_contract.jl](../scripts/gpu_smoke_contract.jl), atmosphere runtime tests, phase-statistics tests and tolerances in this matrix | `K_{5/6}` helper accuracy is now summarized in a maintained note over `x ∈ [1e-6, 140]`; there is still no frozen external phase-statistics bundle | medium-strong |
| `MV-03` | Core optics: electric field, Fraunhofer, Fresnel, atmospheric field propagation | `A`, `R`, `G`, `P`, `M` | analytic + SPECULA-targeted contract bundle | [runtests.jl](../test/runtests.jl), [reference_harness.jl](../test/reference_harness.jl), [reference_data_specula](../test/reference_data_specula), committed SPECULA reference data under `test/reference_data_specula`, [gpu_smoke_contract.jl](../scripts/gpu_smoke_contract.jl), [profile_atmospheric_field_runtime.jl](../scripts/profile_atmospheric_field_runtime.jl), atmospheric-field tests and scripts | External atmospheric-field evidence is now contract-oriented through deterministic SPECULA-aligned scenarios; this is still narrower than full platform-level numerical equivalence | strong |
| `MV-04` | Detectors and detector-family execution | `A`, `G`, `P`, `M` | analytic + runtime behavior + committed detector fixture artifact | [runtests.jl](../test/runtests.jl), [optional_amdgpu_backends.jl](../test/optional_amdgpu_backends.jl), [optional_cuda_backends.jl](../test/optional_cuda_backends.jl), [profile_ao3k_runtime.jl](../scripts/profile_ao3k_runtime.jl), detector fixture artifact under `benchmarks/results/detectors/`, [2026-04-23-phase2-pvp04.toml](../benchmarks/results/detectors/2026-04-23-phase2-pvp04.toml) | Detector realism is still strongest in integrated runtime scenarios, but a committed fixed-seed detector-family fixture artifact now records family-specific transformations independently, including richer HgCdTe multi-read interactions | medium-strong |
| `MV-05` | Shack-Hartmann WFS | `A`, `R`, `G`, `P`, `M` | OOPAO + narrow SPECULA detector-frame contracts | [reference_harness.jl](../test/reference_harness.jl), [reference_data](../test/reference_data), [reference_data_specula](../test/reference_data_specula), committed SPECULA reference data under `test/reference_data_specula`, [gpu_smoke_contract.jl](../scripts/gpu_smoke_contract.jl), [profile_multi_source_multi_wfs_runtime.jl](../scripts/profile_multi_source_multi_wfs_runtime.jl), [profile_ao3k_runtime.jl](../scripts/profile_ao3k_runtime.jl) | OOPAO remains the primary frozen parity baseline; SPECULA now also anchors a narrow maintained polychromatic detector-frame contract rather than full SH parity | strong |
| `MV-06` | Pyramid and BioEdge WFS | `A`, `R`, `G`, `P`, `M` | OOPAO + narrow SPECULA detector-frame contracts | [reference_harness.jl](../test/reference_harness.jl), [reference_data](../test/reference_data), [reference_data_specula](../test/reference_data_specula), committed SPECULA reference data under `test/reference_data_specula`, [gpu_smoke_contract.jl](../scripts/gpu_smoke_contract.jl), [profile_multi_source_multi_wfs_runtime.jl](../scripts/profile_multi_source_multi_wfs_runtime.jl) | OOPAO remains the primary parity baseline; SPECULA now also anchors a narrow maintained polychromatic pyramid detector-frame contract rather than grouped/platform equivalence | strong |
| `MV-07` | Curvature and Zernike WFS | `A`, `R`, `G`, `P`, `M` | SPECULA-targeted contract bundle | [reference_harness.jl](../test/reference_harness.jl), [reference_data_specula](../test/reference_data_specula), committed SPECULA reference data under `test/reference_data_specula`, [gpu_smoke_contract.jl](../scripts/gpu_smoke_contract.jl), [profile_atmospheric_field_runtime.jl](../scripts/profile_atmospheric_field_runtime.jl) | Current frozen coverage is contract-oriented and scenario-aligned to SPECULA tests; it is not a claim of complete numerical equivalence to the full SPECULA platform | strong |
| `MV-08` | LiFT and gain-sensing camera | `A`, `R`, `G`, `P`, `M` | OOPAO + maintained workflow profile artifact | [reference_harness.jl](../test/reference_harness.jl), [reference_data](../test/reference_data), [runtests.jl](../test/runtests.jl), LiFT/GSC workflow artifact under `benchmarks/results/workflows/`, [2026-04-01-phase1-pvp05.toml](../benchmarks/results/workflows/2026-04-01-phase1-pvp05.toml) | Runtime/profile evidence is now captured through a committed workflow artifact, but it remains narrower than the heavier closed-loop runtime families | strong |
| `MV-09` | Runtime and closed-loop execution | `A`, `R`, `G`, `P`, `M` | OOPAO + maintained runtime profiles | [reference_harness.jl](../test/reference_harness.jl), [reference_data](../test/reference_data), [profile_ao3k_runtime.jl](../scripts/profile_ao3k_runtime.jl), [profile_revolt_hil_runtime.jl](../scripts/profile_revolt_hil_runtime.jl), [gpu_runtime_equivalence_contract.jl](../scripts/gpu_runtime_equivalence_contract.jl) | Current evidence is package-local plus OOPAO traces; broader cross-package benchmarking remains intentionally scoped | strong |
| `MV-10` | Tomography and reconstruction | `A`, `R`, `G`, `M` | pyTomoAO / OOPAO-adjacent frozen references | [tomography.jl](../test/tomography.jl), [reference_harness.jl](../test/reference_harness.jl), [reference_data](../test/reference_data), tomography GPU profile scripts under [scripts](../scripts), tomography benchmark artifacts under `benchmarks/results/tomography/`, [2026-04-02-phase2-psp05.toml](../benchmarks/results/tomography/2026-04-02-phase2-psp05.toml) | Frozen reference and functional coverage remain strong, but the refreshed representative benchmark review still closed as an explicit scoped defer tied to a routine-maintenance timeout budget | medium-strong |
| `MV-11` | GPU backend execution policy | `G`, `P`, `M` | backend parity | [optional_amdgpu_backends.jl](../test/optional_amdgpu_backends.jl), [optional_cuda_backends.jl](../test/optional_cuda_backends.jl), [backend-validation-guide.md](backend-validation-guide.md), [gpu_smoke_contract.jl](../scripts/gpu_smoke_contract.jl), [profile_ao3k_runtime.jl](../scripts/profile_ao3k_runtime.jl), [profile_multi_source_multi_wfs_runtime.jl](../scripts/profile_multi_source_multi_wfs_runtime.jl), [profile_control_loop_runtime.jl](../scripts/profile_control_loop_runtime.jl), backend execution policy summarized in [`backend-validation-guide.md`](backend-validation-guide.md), AMDGPU fallback status summarized in [`backend-validation-guide.md`](backend-validation-guide.md), platform rebaseline artifact under `benchmarks/results/platform/`, [2026-04-02-phase3-psp10.toml](../benchmarks/results/platform/2026-04-02-phase3-psp10.toml) | Backend evidence is split across functional tests, optional smoke, and benchmarks. The remaining ROCm fallback surfaces are explicit, the realistic CPU/AMDGPU/CUDA runtime ladder has a committed post-cleanup rebaseline artifact, and typed `ControlLoopScenario` backend smoke/benchmark evidence is maintained | strong |
| `MV-12` | Tutorials and workflow examples | `A`, `M` | local workflow correctness | [examples/tutorials](../examples/tutorials), [reference_and_tutorials.jl](../test/testsets/reference_and_tutorials.jl), [user-guide.md](user-guide.md) | Tutorials are executed, but they are not external parity evidence by themselves | medium |
| `MV-13` | Grouped runtime and multi-source / multi-WFS orchestration | `A`, `G`, `P`, `M` | Julia-first grouped runtime artifact + SPECULA-informed platform runtime artifact + typed orchestration runtime evidence | [control_and_runtime.jl](../test/testsets/control_and_runtime.jl), [backend_optional_common.jl](../test/backend_optional_common.jl), [profile_multi_source_multi_wfs_runtime.jl](../scripts/profile_multi_source_multi_wfs_runtime.jl), [profile_control_loop_runtime.jl](../scripts/profile_control_loop_runtime.jl), grouped runtime tests and artifacts, grouped runtime artifacts under `benchmarks/results/grouped/`, SPECULA-informed platform artifacts under `benchmarks/results/platform/`, control-loop runtime tests and artifacts, REVOLT-like platform benchmark artifacts, [2026-04-01-gr.toml](../benchmarks/results/grouped/2026-04-01-gr.toml), [2026-04-02-phase2-psp07.toml](../benchmarks/results/platform/2026-04-02-phase2-psp07.toml), [2026-04-03-phase5-psp15.toml](../benchmarks/results/platform/2026-04-03-phase5-psp15.toml), [2026-04-03-phase5-psp16.toml](../benchmarks/results/cross_package/2026-04-03-phase5-psp16.toml) | This now includes direct typed-orchestration runtime evidence and one normalized cross-package platform contract against the neighboring legacy tree. It still does not claim full cross-package platform equivalence | strong |
| `MV-14` | Controllable optics and composite low-order runtime plants | `A`, `R`, `G`, `P`, `M` | deterministic CPU self-check artifacts + backend parity + narrow OOPAO modal/composite baseline | [control_and_runtime.jl](../test/testsets/control_and_runtime.jl), [backend_optional_common.jl](../test/backend_optional_common.jl), [gpu_runtime_equivalence_contract.jl](../scripts/gpu_runtime_equivalence_contract.jl), [generate_multi_optic_runtime_artifact.jl](../scripts/generate_multi_optic_runtime_artifact.jl), [reference_harness.jl](../test/reference_harness.jl), [reference_data](../test/reference_data), controllable-optic self-check tests, expanded controllable-optic self-check tests, control-loop runtime tests and artifacts, [backend-validation-guide.md](backend-validation-guide.md), committed OOPAO reference data under `test/reference_data`, [2026-04-13-multi-optic-hil.toml](../benchmarks/results/platform/2026-04-13-multi-optic-hil.toml) | Current maintained evidence now includes invalid-input contracts, multi-step statefulness, amplitude sweeps, detector-backed invariants, richer composites (`tiptilt + focus + dm`, `steering + focus + dm`), and non-SH self-check surfaces (`Pyramid`, `BioEdge`) with CPU/CUDA/AMDGPU parity. External evidence is still narrow rather than broad, but it now includes both Cartesian tip/tilt modal responses and one representative `tiptilt + dm` composite plant on diffractive `ShackHartmannWFS`, `Pyramid`, and `BioEdge`. Broader composite families are still internal-artifact and backend-parity validated rather than externally cross-simulator validated | strong |

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
  - atmospheric field propagation
  - Zernike WFS
  - Curvature WFS

## Known Non-Equivalences and Scope Limits

- Infinite-atmosphere fidelity is not currently claimed against a frozen OOPAO
  or SPECULA statistics bundle.
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

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
  - CPU-to-maintained-accelerator functional parity or smoke surfaces;
    currently AMDGPU and CUDA
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
| `MV-01` | Atmosphere: finite and infinite multilayer propagation | `A`, `G`, `P`, `M` | analytic + backend parity | [atmosphere tests](../test/testsets/atmosphere.jl), [gpu_smoke_contract.jl](../scripts/gpu_smoke_contract.jl), [profile_atmosphere_runtime.jl](../scripts/profile_atmosphere_runtime.jl), atmosphere statistics artifact under `benchmarks/results/atmosphere/`, [2026-04-01-phase1-pvp02.toml](../benchmarks/results/atmosphere/2026-04-01-phase1-pvp02.toml) | Explicit elapsed/absolute time, stable epochs, analytic wind motion, fixed-seed replay, NGS/LGS direction order invariance, frozen source inputs, and allocation-free CPU rendering are validated. No frozen external infinite-atmosphere parity bundle exists yet, but a committed fixed-seed finite/infinite statistics artifact records variance agreement, stationarity, and non-periodicity evidence. | strong |
| `MV-02` | Phase statistics and covariance helpers | `A`, `G`, `M` | analytic | [runtests.jl](../test/runtests.jl), [gpu_smoke_contract.jl](../scripts/gpu_smoke_contract.jl), atmosphere runtime tests, phase-statistics tests and tolerances in this matrix | `K_{5/6}` helper accuracy is now summarized in a maintained note over `x ∈ [1e-6, 140]`; there is still no frozen external phase-statistics bundle | medium-strong |
| `MV-03` | Core optics: electric field, Fraunhofer, Fresnel, direct imaging, atmospheric field propagation | `A`, `R`, `G`, `P`, `M` | analytic + SPECULA-targeted contract bundle | [runtests.jl](../test/runtests.jl), [direct_science.jl](../test/testsets/direct_science.jl), [backend_optional_common.jl](../test/backend_optional_common.jl), [reference_harness.jl](../test/reference_harness.jl), [reference_data_specula](../test/reference_data_specula), committed SPECULA reference data under `test/reference_data_specula`, [gpu_smoke_contract.jl](../scripts/gpu_smoke_contract.jl), [Gate 0 science-stage artifact](../benchmarks/results/gate0/2026-07-16-pre-hil-06-science-stage.toml), [profile_atmospheric_field_runtime.jl](../scripts/profile_atmospheric_field_runtime.jl), atmospheric-field tests and scripts | Prepared direct imaging writes a caller-owned focal-plane angular, cell-integrated photon-arrival-rate product; resolves a finite integer off-axis shift during preparation according to the declared axis order/signs; preserves same-grid incoherent asterism composition; and retains differing spectral grids as a bundle. Placement remains nearest-sample and periodic rather than subpixel with finite-field loss. Prepared direct imaging has focused CPU correctness and allocation evidence plus current maintained AMDGPU/CUDA parity and residency evidence with scalar indexing disabled. The Gate 0 artifact establishes only warmed serial CPU service-time regression bounds; it is not external-RTC latency or production-scale capacity evidence. External atmospheric-field evidence remains narrower than full platform-level numerical equivalence. | strong |
| `MV-04` | Detectors and detector-family execution | `A`, `G`, `P`, `M` | analytic + runtime behavior + committed detector fixture and latency artifacts | [runtests.jl](../test/runtests.jl), [optional_amdgpu_backends.jl](../test/optional_amdgpu_backends.jl), [optional_cuda_backends.jl](../test/optional_cuda_backends.jl), [profile_ao3k_runtime.jl](../scripts/profile_ao3k_runtime.jl), detector fixture artifact under `benchmarks/results/detectors/`, [2026-07-12-detector-mkid-validation.toml](../benchmarks/results/detectors/2026-07-12-detector-mkid-validation.toml), [2026-07-14-detector-hil-latency.toml](../benchmarks/results/detectors/2026-07-14-detector-hil-latency.toml) | Detector realism is still strongest in integrated runtime scenarios, but a committed fixed-seed detector-family fixture artifact records family-specific transformations independently, including HgCdTe multi-read interactions plus APD, SPAD, and MKID counting paths. Prepared `IntensityMap` acquisition validates the plane, rate normalization, spatial measure, incoherent policy, spectral channel, numeric type, backend, and device before repeated capture; it also sizes and zero-initializes maintained multi-read products and rejects predictable defect/background/CMOS/readout-schedule shape errors before buffer mutation. Detector exposure is applied exactly once, and both the first prepared conventional capture and warmed CPU path are allocation-free. A non-null presampling response is applied on the optical grid before physical-pixel integration; response with `psf_sampling > 1` is rejected until an explicit optical-grid mapping is prepared. Scalar QE and sampled QE interpolation are supported for a declared monochromatic channel. Generic spectral source-aware capture uses a flux-weighted effective QE. Pyramid frame-detector paths apply sampled QE to each wavelength contribution before incoherent accumulation. Prepared diffractive Shack–Hartmann formation retains distinct wavelength grids as separate rate products so each can use a channel-specific acquisition mapping; its legacy single-product path remains restricted to a common wavelength grid. Presampling detector response and post-collection IPC are explicit, separate models; `detector_mtf` validates the interior, infinite-grid transfer magnitude of the realized discrete acquisition kernel, not a finite-frame global or continuous subpixel-aperture MTF. Finite response application uses non-amplifying zero extension, so signal outside detector support is lost rather than clamped back into edge pixels. No HgCdTe blur or coupling profile is implicit. CMOS, sCMOS, and quantitative low-noise CMOS share a generic sensor architecture with row, column, per-pixel, output-group, timing, and defect-map composition; core intentionally carries no named camera profiles. Skipper CCD sampling uses an online independent-read mean with bounded frame-sized storage and does not yet model correlated 1/f noise or adaptive stopping. EMCCD CIC is per frame, CPU stochastic multiplication uses a conditional Gamma model, accelerator stochastic multiplication remains a moment approximation, and photon-counting efficiency is Bernoulli acceptance after thresholding. Linear APDs use scalar/vector channel topology rather than area frames; the model does not yet include bandwidth-dependent analog electronics or correlated channel noise. The HdrHistogram latency artifact is a 64×64, one-thread, serial in-process detector-boundary baseline with three 100,000-sample repetitions and zero steady-state allocations; it is not external-RTC or open-loop evidence. MKID support is an accumulated counting-array HIL surface with an inclusive meter-valued source passband plus energy-resolution and timing metadata; matrix-only capture assumes prefiltered input, and the model does not simulate per-photon timestamp/energy event streams | medium-strong |
| `MV-05` | Shack-Hartmann WFS | `A`, `R`, `G`, `P`, `M` | OOPAO parity + SPECULA legacy characterization | [reference_harness.jl](../test/reference_harness.jl), [reference_data](../test/reference_data), [reference_data_specula](../test/reference_data_specula), committed SPECULA reference data under `test/reference_data_specula`, [wfs_stage_contracts.jl](../test/testsets/wfs_stage_contracts.jl), [backend_optional_common.jl](../test/backend_optional_common.jl), [Gate 0 Shack-Hartmann artifact](../benchmarks/results/gate0/2026-07-16-pre-hil-07-shack-hartmann-stage.toml), [profile_multi_source_multi_wfs_runtime.jl](../scripts/profile_multi_source_multi_wfs_runtime.jl), [profile_ao3k_runtime.jl](../scripts/profile_ao3k_runtime.jl) | OOPAO remains the primary frozen parity baseline. The narrow SPECULA polychromatic detector-frame data characterizes the retired reference-wavelength index-grid approximation through a test-only adapter; it is not a supported production contract. The maintained prepared path decomposes microlens optics, layout/calibration, detector acquisition, and centroid estimation. It forms cell-integrated photon-rate mosaics from explicit pupil/field inputs, preserves LGS elongation/sodium profiles and same-wavelength asterisms, and retains distinct wavelength grids as separate native-sampling products in an `OpticalProductBundle`; the legacy single-product `measure!` path remains common-grid only. `Asterism` remains a flat, common-wavelength leaf-source list and one acquisition requires a common optical calibration signature, so mixed NGS/LGS leaves or differing LGS geometry/profile leaves use independent WFS paths. Extended-source quadrature and rate conservation are supported by the legacy convenience path, but direct WFS propagation does not yet apply component angular offsets to focal-plane morphology; direction-dependent atmosphere rendering remains path-local. The Gate 0 artifact establishes warmed serial CPU service-time regression bounds only; it is not external-RTC latency or open-loop capacity evidence. | strong |
| `MV-06` | Pyramid and BioEdge WFS | `A`, `R`, `G`, `P`, `M` | OOPAO + narrow SPECULA detector-frame contracts | [reference_harness.jl](../test/reference_harness.jl), [reference_data](../test/reference_data), [reference_data_specula](../test/reference_data_specula), committed SPECULA reference data under `test/reference_data_specula`, [gpu_smoke_contract.jl](../scripts/gpu_smoke_contract.jl), [profile_multi_source_multi_wfs_runtime.jl](../scripts/profile_multi_source_multi_wfs_runtime.jl) | OOPAO remains the primary parity baseline; SPECULA now also anchors a narrow maintained polychromatic pyramid detector-frame contract rather than grouped/platform equivalence. Pyramid has explicit spectral and extended-source paths. Extended-source quadrature conserves optical rate, but direct Pyramid propagation does not yet apply component angular offsets to focal-plane morphology; direction-dependent atmosphere rendering remains path-local. Pyramid and BioEdge asterism acquisitions require a common optical calibration signature; mixed NGS/LGS leaves or differing LGS geometry/profile leaves use independent WFS paths. BioEdge supports leaf sources and a flat common-wavelength asterism; unsupported spectral or extended wrappers fail structurally instead of being treated as their base source. | strong |
| `MV-07` | Curvature and Zernike WFS | `A`, `R`, `G`, `P`, `M` | SPECULA-targeted contract bundle | [reference_harness.jl](../test/reference_harness.jl), [reference_data_specula](../test/reference_data_specula), committed SPECULA reference data under `test/reference_data_specula`, [gpu_smoke_contract.jl](../scripts/gpu_smoke_contract.jl), [profile_atmospheric_field_runtime.jl](../scripts/profile_atmospheric_field_runtime.jl) | Current frozen coverage is contract-oriented and scenario-aligned to SPECULA tests; it is not a claim of complete numerical equivalence to the full SPECULA platform. Zernike and direct Curvature sensing currently accept leaf sources only; unsupported composite wrappers fail structurally. Atmosphere-aware Curvature asterisms are flat directional lists whose leaves share one optical calibration signature; heterogeneous NGS/LGS or LGS-profile paths require independent sensors. | strong |
| `MV-08` | LiFT and gain-sensing camera | `A`, `R`, `G`, `P`, `M` | OOPAO + maintained workflow profile artifact | [reference_harness.jl](../test/reference_harness.jl), [reference_data](../test/reference_data), [runtests.jl](../test/runtests.jl), LiFT/GSC workflow artifact under `benchmarks/results/workflows/`, [2026-04-01-phase1-pvp05.toml](../benchmarks/results/workflows/2026-04-01-phase1-pvp05.toml) | Runtime/profile evidence is now captured through a committed workflow artifact, but it remains narrower than the heavier closed-loop runtime families | strong |
| `MV-09` | Runtime and closed-loop execution | `A`, `R`, `G`, `P`, `M` | OOPAO + maintained runtime profiles | [reference_harness.jl](../test/reference_harness.jl), [reference_data](../test/reference_data), [control_and_runtime.jl](../test/testsets/control_and_runtime.jl), [backend_optional_common.jl](../test/backend_optional_common.jl), [profile_ao3k_runtime.jl](../scripts/profile_ao3k_runtime.jl), [profile_revolt_hil_runtime.jl](../scripts/profile_revolt_hil_runtime.jl), [gpu_runtime_equivalence_contract.jl](../scripts/gpu_runtime_equivalence_contract.jl), [2026-07-13-control-operators.toml](../benchmarks/results/platform/2026-07-13-control-operators.toml) | Current evidence includes allocation-free same-source CPU runtime, explicit CPU-HIL versus device-resident execution plans, validated dense/factorized/controlled reconstructor residency, caller-owned calibration storage, and distinct WFS/science-source atmosphere rendering on CPU and maintained GPU smoke. The factorized performance artifact uses a synthetic truncated operator; it does not validate a production optical rank. Device-resident stepping still contains synchronization owned by individual optical kernels, so the plan is a residency and boundary contract rather than a claim that the complete optical graph is one asynchronous launch chain. It remains package-local plus OOPAO traces; broader cross-package benchmarking is intentionally scoped | strong |
| `MV-10` | Tomography and reconstruction | `A`, `R`, `G`, `M` | pyTomoAO / OOPAO-adjacent frozen references | [tomography.jl](../test/tomography.jl), [reference_harness.jl](../test/reference_harness.jl), [reference_data](../test/reference_data), tomography GPU profile scripts under [scripts](../scripts), tomography benchmark artifacts under `benchmarks/results/tomography/`, [2026-04-02-phase2-psp05.toml](../benchmarks/results/tomography/2026-04-02-phase2-psp05.toml) | Frozen reference and functional coverage remain strong, but the refreshed representative benchmark review still closed as an explicit scoped defer tied to a routine-maintenance timeout budget | medium-strong |
| `MV-11` | GPU backend execution policy | `G`, `P`, `M` | backend parity | [optional_amdgpu_backends.jl](../test/optional_amdgpu_backends.jl), [optional_cuda_backends.jl](../test/optional_cuda_backends.jl), [backend-validation-guide.md](backend-validation-guide.md), [gpu_smoke_contract.jl](../scripts/gpu_smoke_contract.jl), [profile_ao3k_runtime.jl](../scripts/profile_ao3k_runtime.jl), [profile_multi_source_multi_wfs_runtime.jl](../scripts/profile_multi_source_multi_wfs_runtime.jl), [profile_control_loop_runtime.jl](../scripts/profile_control_loop_runtime.jl), backend execution policy summarized in [`backend-validation-guide.md`](backend-validation-guide.md), platform rebaseline artifacts under `benchmarks/results/platform/`, [2026-04-02-phase3-psp10.toml](../benchmarks/results/platform/2026-04-02-phase3-psp10.toml), [2026-07-14-wsl-cuda-local-amdgpu.toml](../benchmarks/results/platform/2026-07-14-wsl-cuda-local-amdgpu.toml) | Backend evidence is split across functional tests, maintained hardware targets, and benchmarks. AMDGPU and CUDA both have current Julia 1.12.6 hardware validation on the benchmarked checkout. The cross-host artifact preserves physical command parity, synchronized device-ready latency, CUDA host-readout latency, allocations, and run dispersion; it remains a two-host diagnostic rather than a platform-independent capacity claim | strong |
| `MV-12` | Tutorials and workflow examples | `A`, `M` | local workflow correctness | [examples/tutorials](../examples/tutorials), [reference_and_tutorials.jl](../test/testsets/reference_and_tutorials.jl), [user-guide.md](user-guide.md) | Tutorials are executed, but they are not external parity evidence by themselves | medium |
| `MV-13` | Grouped runtime and multi-source / multi-WFS orchestration | `A`, `G`, `P`, `M` | Julia-first grouped runtime artifact + SPECULA-informed platform runtime artifact + typed orchestration runtime evidence | [control_and_runtime.jl](../test/testsets/control_and_runtime.jl), [backend_optional_common.jl](../test/backend_optional_common.jl), [profile_multi_source_multi_wfs_runtime.jl](../scripts/profile_multi_source_multi_wfs_runtime.jl), [profile_control_loop_runtime.jl](../scripts/profile_control_loop_runtime.jl), grouped runtime tests and artifacts, grouped runtime artifacts under `benchmarks/results/grouped/`, SPECULA-informed platform artifacts under `benchmarks/results/platform/`, control-loop runtime tests and artifacts, REVOLT-like platform benchmark artifacts, [2026-04-01-gr.toml](../benchmarks/results/grouped/2026-04-01-gr.toml), [2026-04-02-phase2-psp07.toml](../benchmarks/results/platform/2026-04-02-phase2-psp07.toml), [2026-04-03-phase5-psp15.toml](../benchmarks/results/platform/2026-04-03-phase5-psp15.toml), [2026-04-03-phase5-psp16.toml](../benchmarks/results/cross_package/2026-04-03-phase5-psp16.toml) | This includes direct typed-orchestration runtime evidence, allocation-free shared optical-arm execution with one atmosphere advance, reuse of one same-arm photon-arrival-rate image by multiple detectors, and one normalized cross-package platform contract against the neighboring legacy tree. Each detector applies its own response and exposure after optical formation. The replacement prepared direct-imaging path is validated by focused CPU tests and the maintained AMDGPU and CUDA hardware targets, including shared-arm detector fan-out and device-resident science-frame output. Auxiliary WFS signals remain separate control inputs rather than being implicitly concatenated into one reconstructor. It still does not claim full cross-package platform equivalence. | strong |
| `MV-14` | Controllable optics and composite low-order runtime plants | `A`, `R`, `G`, `P`, `M` | deterministic CPU self-check artifacts + backend parity + narrow OOPAO modal/composite baseline | [control_and_runtime.jl](../test/testsets/control_and_runtime.jl), [backend_optional_common.jl](../test/backend_optional_common.jl), [gpu_runtime_equivalence_contract.jl](../scripts/gpu_runtime_equivalence_contract.jl), [generate_multi_optic_runtime_artifact.jl](../scripts/generate_multi_optic_runtime_artifact.jl), [reference_harness.jl](../test/reference_harness.jl), [reference_data](../test/reference_data), controllable-optic self-check tests, expanded controllable-optic self-check tests, control-loop runtime tests and artifacts, [backend-validation-guide.md](backend-validation-guide.md), committed OOPAO reference data under `test/reference_data`, [2026-04-13-multi-optic-hil.toml](../benchmarks/results/platform/2026-04-13-multi-optic-hil.toml) | Current maintained evidence includes invalid-input contracts, multi-step statefulness, amplitude sweeps, detector-backed invariants, richer composites (`tiptilt + focus + dm`, `steering + focus + dm`), and non-SH self-check surfaces (`Pyramid`, `BioEdge`) with current CPU/CUDA/AMDGPU parity. External evidence is still narrow rather than broad, but it includes both Cartesian tip/tilt modal responses and one representative `tiptilt + dm` composite plant on diffractive `ShackHartmannWFS`, `Pyramid`, and `BioEdge`. Broader composite families are still internal-artifact and backend-parity validated rather than externally cross-simulator validated | strong |

Across `MV-03` through `MV-09` and `MV-13`, the telescope now owns spatial
geometry but no cadence. Runtime construction requires a positive explicit
`atmosphere_step`, while optical formation produces a rate independent of both
that step and detector exposure. Focused rate, unequal-exposure, density versus
cell-integrated, compatible-sum, typed-bundle, and prepared WFS-stage tests
provide CPU and representative CUDA/AMDGPU residency evidence for the ownership
boundary. The direct-science stage now has focused CPU and maintained
CUDA/AMDGPU evidence for caller-owned formation and detector fan-out. This
remains a synchronized package-local execution contract rather than independent
multi-rate HIL scheduling or external-RTC timing evidence. Shared runtime
detector state and exposure are independent, but noisy detector tuple order
still determines draws from its single runtime RNG. The stage fixtures
establish composition semantics, not another WFS physics model or a fidelity
promotion. The later HIL event scheduler is still
required for independent multi-rate triggers and exposure events, and most WFS
families still combine optical, detector, and estimator state internally.

Elementwise spectral or source accumulation is valid only for compatible
physical grids with declared incoherent-addition semantics. Incompatible
prepared products remain in an `OpticalProductBundle` or require an explicit
prepared mapping; array shape alone is not evidence of compatibility.

For `MV-04`, HgCdTe up-the-ramp sampling currently fits a linear ramp with
independent read noise and retains the read cube and timestamps. It shares the
exposure-level photon/dark realization across reads and does not yet model
cosmic-ray segmentation, saturation-aware fitting, or correlated 1/f noise.
EMCCD frame transfer is an acquisition-timing model: it changes reported
one-frame latency and steady-state period, but not the optical or charge
response.

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

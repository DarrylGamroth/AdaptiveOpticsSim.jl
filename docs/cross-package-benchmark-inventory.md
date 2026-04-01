# Cross-Package Benchmark Inventory

Date: 2026-03-31

Status: completed inventory backing active harness

Plan traceability:

- [`PLAN-05`](./package-review-action-plan.md)
- review IDs: `PR-21`, `PR-22`, `PR-23`, `PR-24`, `PR-25`

## Purpose

This document inventories the currently available assets for maintained
cross-package benchmarking against OOPAO- and SPECULA-based models.

This document began as the Phase 0 inventory for cross-package benchmarking.
Phase 5 is now implemented, so this file remains as the source inventory behind
the maintained harness rather than as the primary execution guide.

The inventory goals are to:

- record which scenario sources already exist
- identify scenario families that are actually comparable
- assign candidate `compact`, `medium`, and `representative` seeds
- record the main comparability limits before Phase 5 benchmark work starts

For the active harness and archived evidence, see:

- [`cross-package-benchmark-harness.md`](./cross-package-benchmark-harness.md)
- [`../benchmarks/contracts/cross_package.toml`](../benchmarks/contracts/cross_package.toml)
- [`../benchmarks/results/cross_package`](../benchmarks/results/cross_package)

## Scenario Sources

### 1. Current AdaptiveOpticsSim benchmark surfaces

Maintained local scripts:

- [`scripts/profile_ao3k_runtime.jl`](../scripts/profile_ao3k_runtime.jl)
- [`scripts/profile_atmosphere_runtime.jl`](../scripts/profile_atmosphere_runtime.jl)
- [`scripts/profile_atmospheric_field_runtime.jl`](../scripts/profile_atmospheric_field_runtime.jl)
- [`scripts/profile_external_optics_hil.jl`](../scripts/profile_external_optics_hil.jl)
- [`scripts/profile_lgs_sh_runtime.jl`](../scripts/profile_lgs_sh_runtime.jl)
- [`scripts/profile_mixed_sh_asterism_runtime.jl`](../scripts/profile_mixed_sh_asterism_runtime.jl)
- [`scripts/profile_multi_source_multi_wfs_runtime.jl`](../scripts/profile_multi_source_multi_wfs_runtime.jl)
- [`scripts/profile_pixel_output_runtime.jl`](../scripts/profile_pixel_output_runtime.jl)
- [`scripts/profile_revolt_hil_runtime.jl`](../scripts/profile_revolt_hil_runtime.jl)
- [`scripts/profile_zernike_runtime.jl`](../scripts/profile_zernike_runtime.jl)

Reference and parity assets:

- frozen OOPAO bundle described in
  [`docs/oopao-reference-datasets.md`](./oopao-reference-datasets.md)
- harness in [`test/reference_harness.jl`](../test/reference_harness.jl)

### 2. REVOLT / SPECULA assets

Relevant available files in [`../REVOLT`](../REVOLT):

- SPECULA REVOLT sequencer:
  [`../REVOLT/Python/specula/revolt.py`](../REVOLT/Python/specula/revolt.py)
- SPECULA SH modal SCAO config:
  [`../REVOLT/Python/specula/params_revolt_modal.yml`](../REVOLT/Python/specula/params_revolt_modal.yml)
- SPECULA PWFS modal SCAO config:
  [`../REVOLT/Python/specula/PWFS/params_revolt_modal_PWFS.yml`](../REVOLT/Python/specula/PWFS/params_revolt_modal_PWFS.yml)
- associated calibration assets in the same tree:
  - `calib_subaps_revolt.yml`
  - `calib_im_rec_modal_revolt.yml`
  - `generate_IF_and_KL_basis.py`
  - `prepare_pushpull_amplitudes.py`

### 3. AdaptiveOpticsSim REVOLT-aligned branch assets

Relevant files in [`../AdaptiveOpticsSim.jl-revolt-real`](../AdaptiveOpticsSim.jl-revolt-real):

- REVOLT-aligned scenario builders:
  - [`../AdaptiveOpticsSim.jl-revolt-real/scripts/revolt/common.jl`](../AdaptiveOpticsSim.jl-revolt-real/scripts/revolt/common.jl)
  - [`../AdaptiveOpticsSim.jl-revolt-real/scripts/revolt/shwfs.jl`](../AdaptiveOpticsSim.jl-revolt-real/scripts/revolt/shwfs.jl)
  - [`../AdaptiveOpticsSim.jl-revolt-real/scripts/revolt/pwfs.jl`](../AdaptiveOpticsSim.jl-revolt-real/scripts/revolt/pwfs.jl)
  - [`../AdaptiveOpticsSim.jl-revolt-real/scripts/revolt/pwfs_unmod.jl`](../AdaptiveOpticsSim.jl-revolt-real/scripts/revolt/pwfs_unmod.jl)
- TOML parameter assets:
  - `scripts/revolt/parameter_files/common.toml`
  - `scripts/revolt/parameter_files/shwfs.toml`
  - `scripts/revolt/parameter_files/pwfs.toml`
  - `scripts/revolt/parameter_files/pwfs_unmod.toml`
  - `scripts/revolt/parameter_files/cameras.toml`
- benchmark assets:
  - [`../AdaptiveOpticsSim.jl-revolt-real/benchmarks/assets/revolt_like/`](../AdaptiveOpticsSim.jl-revolt-real/benchmarks/assets/revolt_like)
  - [`../AdaptiveOpticsSim.jl-revolt-real/scripts/profile_revolt_hil_runtime.jl`](../AdaptiveOpticsSim.jl-revolt-real/scripts/profile_revolt_hil_runtime.jl)

## Candidate Scenario Families

### Family `CP-01`: Compact frozen parity cases

Purpose:

- fidelity-first comparison
- low-cost regression and cross-code sanity checking

Available in:

- AdaptiveOpticsSim main via frozen reference harness
- OOPAO via the committed bundle generator/provenance

Candidate cases:

- `psf_baseline`
- `shack_hartmann_geometric_ramp_xy`
- `shack_hartmann_diffractive_ramp`
- `pyramid_diffractive_ramp`
- `bioedge_diffractive_ramp`
- `lift_interaction_matrix`
- compact closed-loop traces from the committed bundle

Recommended class:

- `compact`

Comparability quality:

- strong between AdaptiveOpticsSim and OOPAO
- not a REVOLT/SPECULA benchmark family

### Family `CP-02`: REVOLT-like SH pixel-output HIL runtime

Purpose:

- realistic sensor/runtime pipeline comparison
- detector-backed SH execution with active-actuator map and camera model

Available in:

- AdaptiveOpticsSim main:
  [`scripts/profile_revolt_hil_runtime.jl`](../scripts/profile_revolt_hil_runtime.jl)
- AdaptiveOpticsSim REVOLT-aligned branch:
  [`../AdaptiveOpticsSim.jl-revolt-real/scripts/revolt/shwfs.jl`](../AdaptiveOpticsSim.jl-revolt-real/scripts/revolt/shwfs.jl)
- REVOLT/SPECULA SH config:
  [`../REVOLT/Python/specula/params_revolt_modal.yml`](../REVOLT/Python/specula/params_revolt_modal.yml)

Recommended class:

- `medium`

Current comparability limits:

- the main-repo HIL script uses synthetic REVOLT-like assets rather than the
  exact REVOLT scenario builder
- calibration and reconstructor generation pipelines are not yet normalized
  across Julia and SPECULA
- detector and payload details are now normalized only at the coarse contract
  level documented in [revolt-sh-benchmark-contract.md](./revolt-sh-benchmark-contract.md),
  not yet as a full parity surface

### Family `CP-03`: REVOLT-like PWFS modal SCAO

Purpose:

- compare a realistic modal PWFS control/runtime pipeline across Julia and
  SPECULA-style stacks

Available in:

- AdaptiveOpticsSim REVOLT-aligned branch:
  [`../AdaptiveOpticsSim.jl-revolt-real/scripts/revolt/pwfs.jl`](../AdaptiveOpticsSim.jl-revolt-real/scripts/revolt/pwfs.jl)
- REVOLT/SPECULA PWFS config:
  [`../REVOLT/Python/specula/PWFS/params_revolt_modal_PWFS.yml`](../REVOLT/Python/specula/PWFS/params_revolt_modal_PWFS.yml)
- related maintained main-repo runtime surfaces:
  - [`scripts/profile_ao3k_runtime.jl`](../scripts/profile_ao3k_runtime.jl)
  - [`scripts/profile_multi_source_multi_wfs_runtime.jl`](../scripts/profile_multi_source_multi_wfs_runtime.jl)

Recommended class:

- `representative`

Current comparability limits:

- this family is now first-class in `main` for the Julia-to-Julia comparison,
  with the current normalization rules documented in
  [revolt-pwfs-benchmark-contract.md](./revolt-pwfs-benchmark-contract.md)
- controller/modal basis/calibration alignment still needs a stronger external
  frozen or executable benchmark surface

### Family `CP-04`: Multi-source / multi-WFS grouped execution

Purpose:

- compare grouped execution and batching efficiency

Available in:

- AdaptiveOpticsSim main:
  [`scripts/profile_multi_source_multi_wfs_runtime.jl`](../scripts/profile_multi_source_multi_wfs_runtime.jl)
- OOPAO concepts exist for the underlying WFS families
- SPECULA has the stronger platform-level reference for orchestration breadth

Recommended class:

- `medium`

Current comparability limits:

- there is no direct ready-made OOPAO or SPECULA scenario bundle matching the
  current grouped Julia benchmark
- this is currently a Julia-first benchmark family, not a cross-package one

### Family `CP-05`: Atmospheric field propagation and curvature-through-atmosphere

Purpose:

- compare atmospheric field propagation fidelity/runtime on a SPECULA-like
  surface

Available in:

- AdaptiveOpticsSim main:
  [`scripts/profile_atmospheric_field_runtime.jl`](../scripts/profile_atmospheric_field_runtime.jl)
- SPECULA conceptual baseline in:
  [`../REVOLT/Python/specula/revolt.py`](../REVOLT/Python/specula/revolt.py)
  and associated `AtmoPropagation`-style configs

Recommended class:

- `medium`

Current comparability limits:

- SPECULA-targeted validation bundles exist, but there is still no maintained
  executable benchmark runner for this family in the neighboring trees
- the current scope/defer record is documented in
  [specula-atmo-field-benchmark-scope.md](./specula-atmo-field-benchmark-scope.md)

## Candidate Ladder

### Compact candidate

- Family: `CP-01`
- Primary comparison:
  - AdaptiveOpticsSim vs OOPAO
- Why it is viable now:
  - frozen bundle already exists
  - deterministic harness already exists
  - tolerances and provenance are already documented

### Medium candidate

- Family: `CP-02`
- Primary comparison:
  - AdaptiveOpticsSim main / `revolt-real` vs REVOLT/SPECULA SH
- Why it is viable now:
  - all three source trees already contain relevant scenario assets
  - detector-backed SH runtime is already benchmarked locally
  - the gap is normalization, not missing scenario material

Secondary medium candidate:

- Family: `CP-05`
- Primary comparison:
  - AdaptiveOpticsSim vs SPECULA-style atmospheric propagation
- Why it is not the first medium target:
  - less frozen today than the SH REVOLT family

### Representative candidate

- Family: `CP-03`
- Primary comparison:
  - AdaptiveOpticsSim REVOLT-aligned PWFS vs REVOLT/SPECULA PWFS
- Why it is the representative seed:
  - modal SCAO loop
  - multi-layer atmosphere
  - detector-backed sensing
  - calibration and reconstructor path

## Comparability Constraints

These must be recorded explicitly in any later benchmark harness.

### Atmosphere

- same `r0`, `L0`, altitude ladder, fractional `Cn2`, wind ladder
- same seed policy or frozen atmosphere replay
- same finite vs infinite-screen choice

### Sources and wavelengths

- same guide/science wavelengths
- same magnitude or flux-scaling convention
- same source geometry

### WFS and detector

- same WFS family and modulation settings
- same detector geometry and sampling
- same detector noise/readout assumptions
- same reference-pixel/common-mode policy where relevant

### Controller / reconstructor

- same modal basis or actuator basis
- same interaction-matrix build assumptions
- same delay model and gain law

### Benchmark metrics

Every cross-package benchmark should record:

- build/precompute time
- mean and p95 step time
- steady-state allocation or memory footprint
- fidelity metric:
  - slopes
  - pixel products
  - residual OPD
  - Strehl or equivalent science metric

## Recommended Next Benchmark-Harness Targets

1. Formalize `CP-02` as the first maintained cross-package benchmark family.
2. Freeze a PWFS representative contract for `CP-03` using the REVOLT/SPECULA
   PWFS assets and `revolt-real` Julia builder.
3. Add a SPECULA-focused atmospheric field propagation comparison family for
   `CP-05`.
4. Keep `CP-01` as the compact fidelity baseline and do not overload it with
   runtime benchmarking duties.

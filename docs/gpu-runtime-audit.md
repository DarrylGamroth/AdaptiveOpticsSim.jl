# GPU Runtime Audit

This document tracks the state of the device-resident execution path after the
runtime sprint and the current CUDA/AMDGPU validation workflows.

The timing-oriented scripts in this document use warmed `time_ns()`-based
measurement for operational HIL checks. They are intentionally separate from the
`BenchmarkTools` benchmark entry points under `benchmarks/`, which exist for
canonical CPU/CUDA/AMDGPU benchmarking with tighter statistics.

The maintained benchmark entry points are:

```bash
julia --project=benchmarks benchmarks/benchmark_cpu.jl
julia --project=benchmarks benchmarks/benchmark_cuda.jl
julia --project=benchmarks benchmarks/benchmark_amdgpu.jl
```

The GPU benchmark entry points fail fast with a clear message if the requested
backend is not functional on the current host.

## Standard CUDA Validation Workflow

The maintained CUDA validation entry points are:

- `scripts/gpu_smoke_cuda.jl`
  - broad runtime/device-resident smoke coverage
- `scripts/gpu_builder_cuda.jl`
  - reconstructor/calibration builder coverage for modal and tomography paths
- `scripts/gpu_hil_cuda.jl`
  - combined runtime + builder HIL-oriented smoke coverage
- `scripts/gpu_sync_audit_cuda.jl`
  - timing-oriented audit for runtime and builder-heavy CUDA paths
- `scripts/cpu_crossover_sweep.jl`
  - warmed CPU crossover sweep for runtime and builder cases
- `scripts/gpu_crossover_cuda.jl`
  - warmed CUDA crossover sweep for runtime and builder cases
- `scripts/gpu_profile_model_tomography_cuda.jl`
  - sampling-profile the medium CUDA model-based tomography builder
- `scripts/profile_multi_source_multi_wfs_runtime.jl`
  - warmed stacked asterism and composite multi-WFS runtime profile
  - supports `compact`, `medium`, and `representative` scales
- `scripts/profile_mixed_sh_asterism_runtime.jl`
  - warmed mixed NGS/LGS diffractive Shack-Hartmann runtime profile
- `scripts/profile_zernike_runtime.jl`
  - warmed compact `ZernikeWFS` closed-loop runtime profile
  - supports `compact`, `medium`, and `representative` scales
- `scripts/profile_revolt_hil_runtime.jl`
  - warmed REVOLT-like synthetic SH HIL benchmark with `277` active commands
    and full `352 x 352` pixel output
  - supports detector-family and `default` vs `null` response sweeps
  - now reports warmed allocation bytes per phase and for the full step
- `scripts/profile_external_optics_hil.jl`
  - warmed external-optics HIL proxy with `468` active commands and exported
    `640 x 512` phase output
  - now reports warmed allocation bytes per phase and for the full step
- `scripts/profile_ao3k_runtime.jl`
  - warmed AO3k pyramid runtime profile with `default` vs `null`
    high-detector response sweeps

On a CUDA host, the standard workflow is:

```bash
julia --project=. scripts/gpu_smoke_cuda.jl
julia --project=. scripts/gpu_builder_cuda.jl
julia --project=. scripts/gpu_hil_cuda.jl
julia --project=. scripts/gpu_sync_audit_cuda.jl
julia --project=. scripts/cpu_crossover_sweep.jl
julia --project=. scripts/gpu_crossover_cuda.jl
julia --project=. scripts/gpu_profile_model_tomography_cuda.jl
julia --project=. scripts/profile_multi_source_multi_wfs_runtime.jl cuda
julia --project=. scripts/profile_mixed_sh_asterism_runtime.jl cuda
julia --project=. scripts/profile_zernike_runtime.jl cuda
julia --project=. scripts/profile_revolt_hil_runtime.jl cuda
julia --project=. scripts/profile_external_optics_hil.jl cuda
julia --project=. scripts/profile_ao3k_runtime.jl cuda
```

The `spiders` workstation is the current real-hardware validation host for this
workflow.

## Standard AMDGPU Validation Workflow

The maintained AMDGPU validation entry points are:

- `scripts/gpu_smoke_amdgpu.jl`
  - broad runtime/device-resident smoke coverage on `ROCArray`
- `scripts/gpu_builder_amdgpu.jl`
  - reconstructor/calibration builder coverage on `ROCArray`
- `scripts/gpu_hil_amdgpu.jl`
  - combined runtime + builder HIL-oriented smoke coverage on `ROCArray`
- `scripts/gpu_sync_audit_amdgpu.jl`
  - timing-oriented audit for runtime and builder-heavy AMDGPU paths
- `scripts/gpu_crossover_amdgpu.jl`
  - warmed AMDGPU crossover sweep for runtime and builder cases
- `scripts/gpu_profile_model_tomography_amdgpu.jl`
  - sampling-profile the medium AMDGPU model-based tomography builder
- `scripts/gpu_profile_model_tomography_phases_amdgpu.jl`
  - explicit phase timing for the medium AMDGPU model-based tomography builder
- `scripts/profile_multi_source_multi_wfs_runtime.jl`
  - warmed stacked asterism and composite multi-WFS runtime profile
  - supports `compact`, `medium`, and `representative` scales
- `scripts/profile_mixed_sh_asterism_runtime.jl`
  - warmed mixed NGS/LGS diffractive Shack-Hartmann runtime profile
- `scripts/profile_zernike_runtime.jl`
  - warmed compact `ZernikeWFS` closed-loop runtime profile
  - supports `compact`, `medium`, and `representative` scales
- `scripts/profile_revolt_hil_runtime.jl`
  - warmed REVOLT-like synthetic SH HIL benchmark with `277` active commands
    and full `352 x 352` pixel output
  - supports detector-family and `default` vs `null` response sweeps
  - now reports warmed allocation bytes per phase and for the full step
- `scripts/profile_external_optics_hil.jl`
  - warmed external-optics HIL proxy with `468` active commands and exported
    `640 x 512` phase output
  - now reports warmed allocation bytes per phase and for the full step
- `scripts/profile_ao3k_runtime.jl`
  - warmed AO3k pyramid runtime profile with `default` vs `null`
    high-detector response sweeps

On an AMDGPU host, the standard workflow is:

```bash
julia --project=. scripts/gpu_smoke_amdgpu.jl
julia --project=. scripts/gpu_builder_amdgpu.jl
julia --project=. scripts/gpu_hil_amdgpu.jl
julia --project=. scripts/gpu_sync_audit_amdgpu.jl
julia --project=. scripts/gpu_crossover_amdgpu.jl
julia --project=. scripts/gpu_profile_model_tomography_amdgpu.jl
julia --project=. scripts/gpu_profile_model_tomography_phases_amdgpu.jl
julia --project=. scripts/profile_multi_source_multi_wfs_runtime.jl amdgpu
julia --project=. scripts/profile_mixed_sh_asterism_runtime.jl amdgpu
julia --project=. scripts/profile_zernike_runtime.jl amdgpu
julia --project=. scripts/profile_revolt_hil_runtime.jl amdgpu
julia --project=. scripts/profile_external_optics_hil.jl amdgpu
julia --project=. scripts/profile_ao3k_runtime.jl amdgpu
```

Current AMDGPU caveat:

- FFT-backed runtime paths are native on `ROCArray`.
- Modal/calibration inverse operators now use native rocSOLVER SVD on
  `ROCArray`.
- LiFT normal-equation solves now use native rocSOLVER Cholesky solves with a
  concrete `n×1` ROC RHS buffer.
- Tomography Hermitian right-division now uses native rocSOLVER Cholesky on
  `ROCArray`, with native LU fallback for the non-Cholesky case.
- LiFT fallback now also uses native rocSOLVER SVD on `ROCArray`.
- For the maintained smoke, builder, HIL, and audit surfaces, the dense linear
  algebra path is now AMD-native.

Focused `LinearAlgebra` audit result for AMDGPU:

- maintained ROCArray builder and LiFT paths are now covered by AMD-specific
  extension methods rather than generic `svd(ROCArray)` or generic
  `lu(ROCArray)` dispatch,
- the remaining obvious generic CPU-oriented dense-LA call in the codebase is
  [fitting.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/Tomography/fitting.jl),
  which still uses `pinv(Matrix(...))` and is not part of the maintained AMD
  smoke/builder/HIL surface,
- the core generic fallbacks remain in place intentionally for CPU and
  unsupported GPU workflows, but the maintained AMD path now bypasses them.

Current warmed AMDGPU sync-audit snapshot on this host:

- runtime step mean: about `1.48e6 ns`
- modal build: about `1.88e6 ns`
- interaction-matrix tomography build: about `1.54e6 ns`
- model tomography build: about `6.36e7 ns`
- high-accuracy model tomography build: about `7.24e7 ns`

Current warmed stacked multi-source / multi-WFS runtime snapshot:

- CPU on this host
  - SH asterism: about `3.75e4 ns`
  - Pyramid asterism: about `4.37e5 ns`
  - BioEdge asterism: about `9.98e5 ns`
  - compatible composite step: about `2.40e4 ns`
  - mixed composite step: about `7.18e4 ns`
- AMDGPU on this host
  - SH asterism: about `3.61e5 ns`
  - Pyramid asterism: about `2.13e6 ns`
  - BioEdge asterism: about `3.77e6 ns`
  - compatible composite step: about `8.33e5 ns`
  - mixed composite step: about `2.02e6 ns`
- CUDA on `spiders`
  - SH asterism: about `2.05e5 ns`
  - Pyramid asterism: about `1.72e6 ns`
  - BioEdge asterism: about `3.49e6 ns`
  - compatible composite step: about `6.47e5 ns`
  - mixed composite step: about `1.61e6 ns`

Interpretation:

- stacked multi-source and grouped multi-WFS execution now work across CPU,
  AMDGPU, and CUDA on the maintained runtime surface,
- CPU remains the lowest-latency path for these compact cases,
- CUDA is currently ahead of AMDGPU on the same stacked runtime surface,
- mixed heterogeneous composite execution is still materially slower than
  compatible grouped execution, which matches the current grouped-phase design.

Current warmed compact `ZernikeWFS` runtime snapshot:

- CPU on this host
  - build time: about `1.42e7 ns`
  - runtime step mean: about `1.11e4 ns`
  - about `89.9 kHz`
- AMDGPU on this host
  - build time: about `1.11e7 ns`
  - runtime step mean: about `3.86e5 ns`
  - about `2.59 kHz`
- CUDA on `spiders`
  - build time: about `2.82e7 ns`
  - runtime step mean: about `3.64e5 ns`
  - about `2.75 kHz`

Interpretation:

- the compact detector-backed `ZernikeWFS` runtime surface is now covered by a
  maintained warmed profiler on CPU, AMDGPU, and CUDA,
- CPU remains the lowest-latency path for this tiny `4 x 4` `ZernikeWFS`
  closed-loop case,
- CUDA is slightly ahead of AMDGPU on the maintained GPU profile surface,
- the current profile shape is useful as a regression guard before adding the
  next sensor family.

Current warmed AO3k medium runtime snapshot on this host:

- CPU on this host
  - `default` SAPHIRA high-detector path
    - response: sampled default
    - sampling: CDS
    - correction: reference-pixel
    - about `48.7 Hz`
    - `runtime_alloc_bytes`: about `263536`
    - high-order products:
      - reference frame `(64, 64)`
      - signal frame `(64, 64)`
      - combined frame `(64, 64)`
      - reference cube `(64, 64, 1)`
      - signal cube `(64, 64, 1)`
      - read cube `(64, 64, 2)`
  - `null` + `Fowler(8)` + output correction
    - about `27.0 Hz`
    - `runtime_alloc_bytes`: about `1181104`
    - high-order products:
      - reference cube `(64, 64, 8)`
      - signal cube `(64, 64, 8)`
      - read cube `(64, 64, 16)`
  - dimensions:
    - pupil `160`
    - high-order subapertures `32`
    - high-order frame `(64, 64)`
    - high-order slopes `2048`
    - control modes `1024`

Interpretation:

- the new AO3k profiler is now in place on a maintained system surface that
  now actually uses the SAPHIRA/HgCdTe detector family rather than a generic
  frame-detector proxy,
- the maintained AO3k medium rung now exposes the explicit HgCdTe readout
  products and readout-correction surface at the profiler level,
- heavier SAPHIRA sampling and correction modes materially increase both
  runtime and allocations, so those surfaces are now visible to the benchmark
  matrix instead of being hidden inside detector internals.

## Initial Runtime Ladder Snapshot

In the current benchmark matrix, `medium` means "larger than the compact guard
case and useful for showing scaling/crossover trends." It is not automatically
the same thing as "fully representative of a deployed instrument." The
`representative` rung is the first one intended to approximate deployment-scale
workload shape for that family.

For HIL-oriented systems with unusually large detector payloads, a fixed named
representative case is often more informative than another generic scale rung.
The maintained REVOLT-like synthetic and external-optics HIL profilers below
are meant to cover that gap without depending on external instrument repos.

These HIL profilers now also report warmed allocation bytes. On GPU backends,
those values should be read as host/runtime launch-side allocation overhead,
not total device memory traffic.

### ZernikeWFS Closed-Loop Runtime

Current warmed `runtime_step_mean_ns` frame-rate snapshot:

| scale | CPU | AMDGPU | CUDA |
| --- | ---: | ---: | ---: |
| `compact` | `87.9 kHz` | `2.56 kHz` | `2.45 kHz` |
| `medium` | `4.62 kHz` | `2.57 kHz` | `2.79 kHz` |
| `representative` | `305 Hz` | `1.80 kHz` | `2.69 kHz` |

Dimensions:

- `compact`: pupil `16`, frame `(4, 4)`, signal length `16`, command length `9`
- `medium`: pupil `64`, frame `(16, 16)`, signal length `224`, command length `81`
- `representative`: pupil `128`, frame `(32, 32)`, signal length `856`, command
  length `289`

Interpretation:

- CPU dominates the tiny compact case,
- by the first representative rung, both GPU backends are ahead of CPU,
- CUDA is currently ahead of AMDGPU on this maintained `ZernikeWFS` runtime
  surface.

### Pixel-Output Closed-Loop Runtime

Current warmed `runtime_step_mean_ns` frame-rate snapshot for the AO188-style
pixel-output profile (`sequential` / `direct`):

| scale | CPU | AMDGPU | CUDA |
| --- | ---: | ---: | ---: |
| `compact` | `3.25 kHz` | `1.15 kHz` | `1.32 kHz` |
| `medium` | `865 Hz` | `1.02 kHz` | `1.25 kHz` |
| `representative` | `522 Hz` | `737 Hz` | `1.24 kHz` |

Dimensions:

- `compact`: pupil `64`, high-order subap `8`, low-order subap `2`,
  command modes `64`, frame `(64, 8, 8)`
- `medium`: pupil `112`, high-order subap `14`, low-order subap `2`,
  command modes `188`, frame `(196, 8, 8)`
- `representative`: pupil `160`, high-order subap `20`, low-order subap `4`,
  command modes `512`, frame `(400, 8, 8)`

Interpretation:

- the compact AO188-style rung is still CPU-first,
- GPU is already competitive by the medium rung,
- CUDA is clearly ahead at the current representative rung,
- the current medium CPU profile remains effectively allocation-flat with
  `runtime_alloc_bytes = 64`.

### Fixed Representative HIL Cases

Current warmed fixed-case representative snapshot:

| case | CPU | AMDGPU | CUDA |
| --- | ---: | ---: | ---: |
| `revolt_like_hil` | `177 Hz` | `487 Hz` | `1502 Hz` |
| `external_optics_hil` | `2.36 kHz` | `2.00 kHz` | `6.72 kHz` |

REVOLT-like synthetic SH HIL dimensions:

- active commands `277`
- extrapolated command length `277`
- DM grid commands `361` (`19 x 19`)
- synthetic DM277-like active-actuator layout from
  `benchmarks/assets/revolt_like/revolt_like_dmActuatorMap_277.csv`
- synthetic sparse extrapolation from
  `benchmarks/assets/revolt_like/revolt_like_dmExtrapolation.csv`
- pupil resolution `352`
- `16 x 16` SH subapertures
- ROI `22 x 22`
- full pixel-output mosaic `(352, 352)`

Current warmed phase timing snapshot for `revolt_like_hil`:

| backend | command map | DM apply | sense | mosaic | total |
| --- | ---: | ---: | ---: | ---: | ---: |
| CPU | `16.6 µs` | `243 µs` | `5.45 ms` | `5.71 ms` | `5.66 ms` |
| AMDGPU | `204 µs` | `417 µs` | `1.87 ms` | `2.01 ms` | `2.05 ms` |
| CUDA | `60.2 µs` | `128 µs` | `629 µs` | `661 µs` | `666 µs` |

Interpretation:

- this synthetic case is far heavier in detector pixels than the generic AO188-style
  `representative` rung,
- the large `352 x 352` pixel-output contract is enough to make GPU execution
  clearly worthwhile,
- CUDA is currently much stronger than AMDGPU on this fixed SH HIL case,
- CPU remains respectable but is no longer the best latency path once the
  detector payload is large enough.

Focused detector-response overhead snapshot on the same REVOLT-like HIL case:

| backend | sensor | response | rate | total alloc bytes |
| --- | --- | --- | ---: | ---: |
| CPU | `CMOSSensor` | `default` (`gaussian`) | `107 Hz` | `32` |
| CPU | `CMOSSensor` | `null` | `115 Hz` | `32` |
| CPU | `CCDSensor` | `default` (`none`) | `128 Hz` | `32` |
| AMDGPU | `CMOSSensor` | `default` (`gaussian`) | `317 Hz` | `245040` |
| AMDGPU | `CMOSSensor` | `null` | `514 Hz` | `83264` |

Interpretation:

- the first maintained Gaussian pixel-response path has a visible cost on the
  large HIL surface,
- on CPU the effect is mainly latency rather than allocation growth,
- on AMDGPU the current response path adds both latency and host/runtime
  allocation overhead,
- detector-response benchmarking is therefore worth keeping alongside the base
  HIL timing surfaces.

External-optics HIL dimensions:

- active commands `468`
- DM grid commands `576` (`24 x 24`)
- generated DM phase grid `(640, 640)`
- exported phase crop `(640, 512)`
- downstream image contract `(640, 512)`

Current warmed phase timing snapshot for `external_optics_hil`:

| backend | command map | DM phase | export | total |
| --- | ---: | ---: | ---: | ---: |
| CPU | `297 ns` | `245 µs` | `376 µs` | `424 µs` |
| AMDGPU | `198 µs` | `484 µs` | `518 µs` | `499 µs` |
| CUDA | `27.8 µs` | `119 µs` | `150 µs` | `149 µs` |

Interpretation:

- this case benchmarks only the part AdaptiveOpticsSim owns in an external
  optics pipeline: command mapping, DM phase generation, and phase export,
- CPU is still competitive because the work is dominated by one DM apply plus a
  crop rather than a full in-package detector model,
- CUDA is already clearly worthwhile for this surface,
- AMDGPU is closer to CPU than to CUDA on the current implementation,
- the maintained CPU path is allocation-free after warmup on this surface:
  `command_map_alloc_bytes = 0`, `dm_phase_alloc_bytes = 0`,
  `export_alloc_bytes = 0`, `total_alloc_bytes = 0`.

### Multi-Source / Multi-WFS Runtime

Current warmed `composite_compatible_mean_ns` frame-rate snapshot:

| scale | CPU | AMDGPU | CUDA |
| --- | ---: | ---: | ---: |
| `compact` | `38.1 kHz` | `1.31 kHz` | `1.56 kHz` |
| `medium` | `6.91 kHz` | `645 Hz` | `977 Hz` |
| `representative` | `1.14 kHz` | `441 Hz` | `696 Hz` |

Dimensions:

- `compact`: `2` sources, asterism pupil `32`, runtime pupil `16`, `2`
  compatible branches
- `medium`: `4` sources, asterism pupil `64`, runtime pupil `32`, `3`
  compatible branches
- `representative`: `6` sources, asterism pupil `96`, runtime pupil `72`, `4`
  compatible branches

Additional current warmed asterism sensor snapshots:

| metric | compact CPU / AMDGPU / CUDA | medium CPU / AMDGPU / CUDA | representative CPU / AMDGPU / CUDA |
| --- | --- | --- | --- |
| SH asterism | `9.44 kHz / 3.09 kHz / 4.78 kHz` | `3.24 kHz / 2.35 kHz / 4.35 kHz` | `821 Hz / 1.23 kHz / 3.36 kHz` |
| Pyramid asterism | `2.14 kHz / 584 Hz / 660 Hz` | `202 Hz / 328 Hz / 371 Hz` | `26.9 Hz / 108.8 Hz / 247.7 Hz` |
| BioEdge asterism | `962 Hz / 293 Hz / 267 Hz` | `86.7 Hz / 156 Hz / 166 Hz` | `11.4 Hz / 59.3 Hz / 110.3 Hz` |

Interpretation:

- the compact grouped-runtime case remains CPU-first,
- larger multi-source and diffractive asterism workloads become much more
  favorable to GPU backends,
- CUDA is currently ahead of AMDGPU on the maintained stacked-runtime surface.

Current warmed mixed NGS/LGS diffractive SH runtime snapshot:

- CPU on this host
  - detector-backed mixed asterism: about `9.20e5 ns`
  - about `1087 Hz`
- AMDGPU on this host
  - detector-backed mixed asterism: about `8.62e5 ns`
  - about `1160 Hz`

Current mixed NGS/LGS SH runtime equivalence snapshot on AMDGPU:

- `spot_cube` max abs: about `3.91e-2`
- `spot_cube` max rel: about `1.46e-5`
- `slopes` max abs: about `1.62e-6`
- `slopes` max rel: about `1.58e-5`

Interpretation:

- the mixed detector-backed SH path now preserves exported `spot_cube`
  semantics instead of leaving only the last-source frame in the runtime
  readout,
- the maintained AMDGPU mixed SH path remains faster than CPU on this host for
  the profiled `14 x 14` detector-backed case,
- the focused CPU-vs-AMDGPU mixed SH equivalence is now tight for both
  exported pixels and slopes.

## Canonical BenchmarkTools Snapshot

Current `BenchmarkTools` snapshots on this host are:

- CPU
  - `runtime_step`
    - median: about `1.11e4 ns`
    - mean: about `1.17e4 ns`
    - memory: `32 bytes`
    - allocs: `1`
  - `model_tomography_build`
    - median: about `1.08e10 ns`
    - mean: about `1.08e10 ns`
    - memory: about `7.87e8 bytes`
    - allocs: `600`
  - `model_tomography_high_accuracy_build`
    - median: about `1.32e10 ns`
    - mean: about `1.32e10 ns`
    - memory: about `1.57e9 bytes`
    - allocs: `620`
- AMDGPU
  - `runtime_step`
    - median: about `1.56e6 ns`
    - mean: about `1.61e6 ns`
    - memory: about `3.67e5 bytes`
    - allocs: `8511`
  - `interaction_tomography_build`
    - median: about `3.58e6 ns`
    - mean: about `3.69e6 ns`
    - memory: about `3.56e5 bytes`
    - allocs: `9321`
- CUDA on `spiders`
  - `runtime_step`
    - median: about `1.41e6 ns`
    - mean: about `1.40e6 ns`
    - memory: about `1.80e5 bytes`
    - allocs: `4195`
  - `model_tomography_build`
    - median: about `1.49e9 ns`
    - mean: about `1.49e9 ns`
    - memory: about `2.38e5 bytes`
    - allocs: `5371`
  - `model_tomography_high_accuracy_build`
    - median: about `1.53e9 ns`
    - mean: about `1.53e9 ns`
    - memory: about `2.90e5 bytes`
    - allocs: `6714`

The current `BenchmarkTools` AMD suite intentionally benchmarks the maintained
interaction-matrix tomography builder surface rather than the heavier
model-based tomography builder, because the detailed model-based tomography
profiling is already covered by the dedicated AMD profiling scripts above.

## Initial Crossover Sweep on AMDGPU

The maintained AMD crossover sweep was run with:

```bash
julia --project=. scripts/gpu_crossover_amdgpu.jl
```

The current warmed results on this host are:

- runtime `compact`
  - AMDGPU: about `1.48e6 ns`
- runtime `small`
  - AMDGPU: about `5.19e6 ns`
- runtime `medium`
  - AMDGPU: about `2.04e7 ns`
- modal builder `compact`
  - AMDGPU: about `1.24e6 ns`
- modal builder `small`
  - AMDGPU: about `9.84e5 ns`
- modal builder `medium`
  - AMDGPU: about `1.93e6 ns`
- interaction tomography builder `compact`
  - AMDGPU: about `2.63e6 ns`
- interaction tomography builder `small`
  - AMDGPU: about `2.17e6 ns`
- interaction tomography builder `medium`
  - AMDGPU: about `4.53e6 ns`
- model tomography builder `compact`
  - AMDGPU: about `1.21e8 ns`
- model tomography builder `small`
  - AMDGPU: about `1.22e9 ns`
- model tomography builder `medium`
  - AMDGPU: about `1.54e9 ns`
- model tomography builder `compact` high accuracy
  - AMDGPU: about `1.34e8 ns`
- model tomography builder `small` high accuracy
  - AMDGPU: about `1.52e9 ns`
- model tomography builder `medium` high accuracy
  - AMDGPU: about `1.76e9 ns`

The interpretation matches the CUDA story:

- the tiny single-loop runtime is still CPU-first,
- modal and interaction-matrix builder cases are still small enough that GPU
  overhead dominates,
- model-based tomography is the first builder family where GPU acceleration is
  already attractive.

## AMDGPU Model-Based Tomography Profile

The profiling scripts:

```bash
julia --project=. scripts/gpu_profile_model_tomography_amdgpu.jl
julia --project=. scripts/gpu_profile_model_tomography_phases_amdgpu.jl
```

were run on the `medium` AMDGPU model-based tomography case.

The key sampling-profile output was:

- default build time: about `1.40e9 ns`
- high-accuracy build time: about `1.49e9 ns`
- default build allocations: about `4.80e5 bytes`
- high-accuracy build allocations: about `4.74e5 bytes`

The phase-timed breakdown was:

- `auto_correlation`: about `1.90e10 ns`
- `cross_correlation`: about `3.48e9 ns`
- `cnz`: about `1.59e9 ns`
- `recstat`: about `1.60e9 ns`
- `css_signal`: about `5.19e8 ns`
- `fit_source_average`: about `4.37e8 ns`
- `extract_submatrix`: about `3.42e8 ns`

Two current conclusions matter for AMD specifically:

- `auto_correlation` is still the dominant tomography-builder hotspot by a wide
  margin, just as it is on CUDA.
- the AMD sampling profile shows real allocator pressure in `_fit_source_average`
  on `ROCArray`, so that path is worth revisiting if AMD builder latency becomes
  a priority.

## Current AMD Precision Guidance

The maintained AMD precision policies are:

- default: `SplitGPUPrecision(Float32, Float32)`
- high accuracy: `SplitGPUPrecision(Float32, Float64)`

Current recommendation:

- keep the default `Float32` builder policy for runtime-oriented and
  throughput-oriented HIL workflows,
- use the high-accuracy builder policy only when tomography conditioning
  requires it.
- for runtime-equivalence expectations, treat `Float32` and `Float64`
  differently:
  - `scripts/gpu_runtime_equivalence_amdgpu.jl` is the maintained fast-runtime
    `Float32` surface and checks AO188 pixels/slopes/commands plus SH/LGS
    detector outputs.
  - `scripts/gpu_runtime_equivalence_high_accuracy_amdgpu.jl` is the maintained
    scientific `Float64` surface and adds a stricter post-command AO188
    equivalence check (`tel_opd`, post-command pixels, and post-command
    slopes).
  - The maintained AO188 simulation itself now lives in
    `examples/support/subaru_ao188_simulation.jl`; the scripts under `scripts/`
    are its audit and equivalence entry points.
  - Keeping the AO188 runtime reconstructor in mapped two-stage form
    (`modal_reconstructor` followed by `M2C`) is still worth it for the fast
    path: it improves the maintained AO188 runtime throughput on CPU and
    AMDGPU, but it does not remove the stricter post-command `Float32`
    equivalence gap by itself.
  - The current AO188 split uses dense DM application for calibration
    (`interaction_matrix`) and a structured separable DM application in the
    runtime path when the DM misregistration keeps the Gaussian influence basis
    separable.
  - That split materially improves the AO188 runtime path on all maintained
    backends, and the current CUDA-specific runtime DM apply path pushes the
    structured separable operator through dedicated CUDA kernels instead of
    relying on the generic CuArray matmul path.
  - The remaining AO188 `Float32` fast-runtime command drift was not in the
    sensing path. It came from calibrating/building the AO188 interaction
    matrices and reconstructors directly on the GPU backend. The maintained
    AO188 GPU path now defaults to CPU-built calibration/operators, then
    uploads the resulting reconstructors to the runtime backend.
  - The same CPU-build/GPU-run rule now applies to the generic runtime
    calibration surfaces too: `ModalReconstructor`, `CalibrationVault`, and
    `ao_calibration` default to CPU build when their source matrices live on a
    GPU backend, while still materializing the resulting operators back onto
    the runtime backend.
  - Current warmed AO188 simulation rates on the structured runtime path are
    about `1.02 kHz` on local CPU, `1.02 kHz` on local AMDGPU, and about
    `1.25 kHz` on CUDA on `spiders`.
  - On AMDGPU, the stricter post-command `Float32` AO188 `tel_opd` max-abs
    error improved from about `2.68e-7` to about `8.94e-8`.
  - On CUDA, the same stricter post-command `Float32` AO188 `tel_opd` max-abs
    error improved from about `2.68e-7` to about `1.49e-7`.
  - With CPU-built AO188 reconstructors uploaded to GPU, the maintained
    fast-runtime AO188 command surface is now tight on both GPU backends:
    AMDGPU command max abs is about `2.37e-8`, and CUDA command max abs is
    about `3.39e-8`.
  - The maintained generic GPU smoke surfaces also remain green with this
    default. On both AMDGPU and CUDA, the closed-loop runtime and
    interaction-matrix reconstructor smoke paths now run with CPU-built
    calibration operators uploaded to the GPU runtime path.
  - The scientific `Float64` AO188 post-command surface remains very tight on
    both GPUs against CPU:
    AMDGPU `tel_opd` max abs about `3.89e-16`, CUDA `tel_opd` max abs about
    `3.33e-16`.

On this host, the warmed sync-audit medium model tomography builder rises from
about `6.58e7 ns` to about `7.27e7 ns` when switching to the high-accuracy
policy, and the standalone medium build profile rises from about `1.40e9 ns` to
about `1.49e9 ns`. So the high-accuracy policy is viable, but it is still a
real builder-cost trade.

## Warmed CPU vs CUDA Snapshot on `spiders`

With the warmed audit scripts on `spiders`, the compact reference case currently
looks like this:

- CPU runtime step mean: about `1.45e4 ns`
- CUDA runtime step mean: about `1.30e6 ns`
- CPU modal build: about `5.16e4 ns`
- CUDA modal build: about `7.86e5 ns`
- CPU interaction-matrix tomography build: about `5.83e3 ns`
- CUDA interaction-matrix tomography build: about `1.46e6 ns`
- CPU model-based tomography build: about `7.72e8 ns`
- CUDA model-based tomography build: about `9.03e7 ns`

For this compact audit case, the interpretation is:

- CPU is clearly better for the tiny steady-state runtime loop.
- CPU is also better for tiny modal and interaction-matrix builder cases.
- CUDA is already materially better for the heavier model-based tomography
  builder.

So the package is currently in the expected regime:

- use CPU for small, latency-dominated compact cases,
- use GPU when the builder/reconstructor workload is large enough to amortize
  device overhead.

The next benchmarking task should therefore be a size sweep to find the
practical crossover points, not more guessing from a single compact case.

## Initial Crossover Sweep on `spiders`

The maintained crossover sweep was run on `spiders` with:

```bash
julia --project=. scripts/cpu_crossover_sweep.jl
julia --project=. scripts/gpu_crossover_cuda.jl
```

The current warmed results are:

- runtime `compact`
  - CPU: about `1.47e4 ns`
  - CUDA: about `1.24e6 ns`
- runtime `small`
  - CPU: about `8.84e4 ns`
  - CUDA: about `4.12e6 ns`
- runtime `medium`
  - CPU: about `3.76e5 ns`
  - CUDA: about `1.56e7 ns`
- modal builder `compact`
  - CPU: about `3.06e4 ns`
  - CUDA: about `5.53e5 ns`
- modal builder `small`
  - CPU: about `1.43e4 ns`
  - CUDA: about `4.39e5 ns`
- modal builder `medium`
  - CPU: about `4.64e4 ns`
  - CUDA: about `6.69e5 ns`
- interaction tomography builder `compact`
  - CPU: about `1.41e4 ns`
  - CUDA: about `1.62e6 ns`
- interaction tomography builder `small`
  - CPU: about `6.20e5 ns`
  - CUDA: about `5.14e6 ns`
- interaction tomography builder `medium`
  - CPU: about `4.78e6 ns`
  - CUDA: about `1.33e7 ns`
- model tomography builder `compact`
  - CPU: about `1.48e9 ns`
  - CUDA: about `1.81e8 ns`
- model tomography builder `small`
  - CPU: about `1.30e10 ns`
  - CUDA: about `1.26e9 ns`
- model tomography builder `medium`
  - CPU: about `1.41e10 ns`
  - CUDA: about `1.46e9 ns`

The current interpretation is:

- CPU wins decisively for the tested single-loop runtime sizes.
- CPU also wins for the tested modal and interaction-matrix builder sizes.
- CUDA already wins decisively for model-based tomography builder workloads,
  even at the compact case.

So the practical crossover is not “GPU for everything.” It is currently:

- compact and medium-sized single-loop HIL runtime: CPU-first
- heavy model-based tomography/reconstructor generation: GPU-first
- modal and interaction-matrix builders: keep on CPU unless larger problem
  sizes show a crossover later

## Model-Based Tomography CUDA Profile on `spiders`

The profiling script:

```bash
julia --project=. scripts/gpu_profile_model_tomography_cuda.jl
```

was run on the `medium` model-based tomography case on `spiders`.

The key output was:

- default build time: about `1.48e9 ns`
- high-accuracy build time: about `1.51e9 ns`
- default build allocations: about `1.96 MB`
- high-accuracy build allocations: about `2.71 MB`

The dominant hotspot is not generic Julia dispatch. It is repeated stream
synchronization during GPU covariance assembly:

- [launch_kernel!](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/Core/backends.jl#L64)
  currently synchronizes after every accelerator kernel launch
- [_covariance_matrix](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/Tomography/reconstructors.jl#L384)
  launches a kernel for every covariance block
- [cross_correlation](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/Tomography/reconstructors.jl#L710)
  calls `_covariance_matrix` inside nested `fit_idx`, `gs`, and `layer` loops
- [build_reconstructor](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/Tomography/reconstructors.jl#L1160)
  spends most of its time in that `cross_correlation` path

The profile showed thousands of samples in CUDA stream synchronization under
`launch_kernel!` and `cross_correlation`, which means the current GPU builder
path is dominated by per-kernel synchronization overhead rather than arithmetic.

So the next optimization target is now clear:

1. add a non-synchronizing internal kernel-launch path for builder assembly
2. synchronize once at the outer builder boundary instead of once per covariance
   block
3. then evaluate whether covariance blocks should be batched or fused further

A first experiment with an internal non-synchronizing launch helper did not
materially change the medium-case timings. That indicates the effective
synchronization cost is not isolated to the explicit launch helper alone; it is
likely also being reintroduced by the surrounding covariance extraction and
accumulation path. So the next useful optimization is a deeper refactor of
covariance block handling, not just another launch-wrapper tweak.

To separate sampling distortion from real wall-clock cost, an explicit phase
timing script was also added:

```bash
julia --project=. scripts/gpu_profile_model_tomography_phases_cuda.jl
```

On the same `medium` CUDA case on `spiders`, the phase-timed breakdown showed:

- `auto_correlation`: about `1.70e10 ns`
- `cross_correlation`: about `2.09e9 ns`
- `cnz`: about `1.32e9 ns`
- `recstat`: about `5.04e8 ns`
- `css_signal`: about `4.87e8 ns`
- `fit_source_average`: about `3.10e8 ns`
- `extract_submatrix`: about `2.85e8 ns`

That changes the optimization priority materially:

1. `auto_correlation` is the real dominant CUDA builder phase.
2. `cross_correlation` and noise-covariance assembly are second-tier costs.
3. `fit_source_average` and submatrix extraction are no longer the first place
   to spend optimization effort.

So the next tomography CUDA optimization pass should target the internals of
`auto_correlation` first, not keep chasing later-stage averaging/gather code.

To make that target concrete, an additional sub-phase script was added:

```bash
julia --project=. scripts/gpu_profile_auto_correlation_cuda.jl
```

On the `medium` CUDA case on `spiders`, the auto-correlation breakdown showed
roughly:

- guide-grid generation: `1.31e9 ns`
- shifted-coordinate generation: `2.44e9 ns`
- covariance assembly: `1.72e9 ns`
- selected-block accumulation: `2.89e8 ns`
- block scatter into the final matrix: `5.88e8 ns`

That makes the remaining priority inside `auto_correlation` clearer:

1. coordinate generation and covariance assembly dominate,
2. accumulation/scatter are secondary costs,
3. builder-side caching helps code structure, but did not materially change the
   end-to-end `model_tomography_build_ns` in the compact warmed sync audit,
   which remains about `8.95e7 ns` on `spiders`.

A narrower follow-up change then precomputed the LGS shifted-coordinate stack in
GPU [cross_correlation](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/Tomography/reconstructors.jl#L781)
instead of rebuilding those coordinates inside every `fit_idx × gs × layer`
iteration. On the same `medium` CUDA case on `spiders`, that moved
`cross_correlation` from about `4.69e9 ns` down to about `2.37e9 ns` in the
phase profiler.

That is a real win for the secondary hotspot, but it still does not materially
change the warmed end-to-end `model_tomography_build_ns`, which stayed in the
same `~8.9e7 ns` band. The reason is unchanged: `auto_correlation` is still the
dominant builder phase by a wide margin, so future optimization effort should
stay focused there.

The next follow-up then tightened `_covariance_matrix` for GPU builder paths:

- precompute the covariance constants once per outer builder phase instead of
  once per covariance block
- skip `materialize_build` when the covariance inputs are already backend-native
  vectors from the shifted-coordinate stacks

On `spiders`, that moved the dedicated `auto_correlation` phase profile from
about `6.25e9 ns` down to about `6.04e9 ns`, with the largest visible drops in:

- shifted-coordinate work: about `2.42e9 ns` to `2.36e9 ns`
- covariance assembly: about `1.68e9 ns` to `1.54e9 ns`

The warmed end-to-end CUDA builder also moved slightly in the right direction:

- `model_tomography_build_ns`: about `8.95e7 ns` to about `8.93e7 ns`

That is still only a modest end-to-end win, but it is the first change in this
phase that improved both the hotspot profile and the warmed builder timing
without introducing a larger structural rewrite.

The next reuse pass then moved GPU tomography covariance assembly to explicit
workspace reuse:

- `auto_correlation` now reuses one `block` buffer across guide-star pairs
- `auto_correlation` and `cross_correlation` now reuse one covariance workspace
  instead of allocating a new covariance matrix per layer
- the GPU covariance helper gained an in-place path so those workspaces can be
  filled directly

On `spiders`, the resulting profile/sample points were:

- dedicated `auto_correlation` timed phase: about `6.04e9 ns` to about
  `6.09e9 ns`
- full model-tomography phase profile:
  - `auto_correlation`: about `1.69e10 ns` to about `1.54e10 ns`
  - `cross_correlation`: about `2.38e9 ns` to about `2.31e9 ns`

Most importantly, the warmed end-to-end CUDA builder moved materially:

- `model_tomography_build_ns`: about `8.93e7 ns` to about `8.67e7 ns`
- `model_tomography_high_accuracy_build_ns`: about `9.13e7 ns` to about
  `8.26e7 ns`

So this reuse-oriented pass is worth keeping. It is the first tomography CUDA
optimization in this phase that clearly reduced the warmed builder wall time by
more than noise.

The next kept pass replaced the GPU `auto_correlation` inner loop's
`covariance_matrix!` plus selected-block accumulation sequence with a direct
selected-block kernel that accumulates the valid submatrix across all layers in
one launch per guide-star pair.

On `spiders`, the older phase profiler became less representative of the new
implementation because it still reports the legacy helper-stage labels, and the
sampled phase totals were noisier:

- dedicated `auto_correlation` timed phase: about `6.09e9 ns` to about
  `6.54e9 ns`
- full model-tomography phase profile:
  - `auto_correlation`: about `1.54e10 ns` to about `1.70e10 ns`
  - `cross_correlation`: about `2.31e9 ns` to about `3.18e9 ns`

But the warmed end-to-end builder improved again while the CUDA builder
fidelity checks stayed green:

- `model_tomography_build_ns`: about `8.67e7 ns` to about `8.26e7 ns`
- `model_tomography_high_accuracy_build_ns`: stayed near `9.13e7 ns`

So this direct selected-block path is also worth keeping. For this builder, the
warmed end-to-end audit is now a more trustworthy decision metric than the
older helper-phase profiler.

## Validated GPU-Resident Surface

The following paths are currently validated on CUDA with
`CUDA.allowscalar(false)` via `scripts/gpu_smoke_matrix.jl`:

- PSF generation
- asterism PSF generation
- detector capture with `NoiseNone`
- detector capture with photon + readout noise
- atmosphere step
- deformable-mirror apply
- geometric Shack-Hartmann
- diffractive Shack-Hartmann
- Shack-Hartmann LGS
- geometric Pyramid
- diffractive Pyramid
- Pyramid detector path
- geometric BioEdge
- diffractive BioEdge
- closed-loop step
- `ClosedLoopRuntime`
- `ClosedLoopRuntime` with science detector
- `CompositeSimulationInterface` aggregated runtime stepping
- runtime reconstructor refresh via `with_reconstructor`
- interaction-matrix build + reconstruction
- gain sensing camera
- LiFT

These paths now execute without scalar indexing failures on GPU.

## Validated GPU Builder Surface

The following builder paths are now validated on CUDA with
`CUDA.allowscalar(false)` via `scripts/gpu_builder_cuda.jl`:

- `CalibrationVault(...; build_backend=GPUArrayBuildBackend(CUDABackendTag))`
- `ModalReconstructor(...; build_backend=GPUArrayBuildBackend(CUDABackendTag))`
- `build_reconstructor(InteractionMatrixTomography(), ...; build_backend=GPUArrayBuildBackend(CUDABackendTag))`
- `build_reconstructor(ModelBasedTomography(), ...; build_backend=GPUArrayBuildBackend(CUDABackendTag))`

The builder contract now also checks numerical agreement against CPU builds for
the exercised cases:

- modal reconstruction output
- interaction-matrix tomography wavefront reconstruction
- interaction-matrix tomography command reconstruction
- model-based tomography wavefront reconstruction

Those checks currently pass on `spiders` with:

- `rtol=1e-5`
- `atol=1e-6`

## Runtime Fallbacks Removed

The following runtime host fallbacks have been removed:

- `ClosedLoopRuntime` DM command update
  - previously used a host `eachindex` loop over `AbstractVector`
  - now uses backend-dispatched scalar CPU / KA accelerator paths
- backend Gaussian noise generation in `randn_backend!`
  - previously staged through a host `Array`
  - now uses a backend-native stateless KA kernel
- photon noise generation in `poisson_noise!`
  - previously staged through a host matrix
  - now uses a backend-native stateless KA kernel
- diffractive Shack-Hartmann centroid extraction
  - previously copied `valid_mask`, `spot_cube`, and `slopes` to host
  - now computes slopes directly on device with a KA kernel
- diffractive Shack-Hartmann calibration/reference helpers
  - previously used host-side `valid_mask` access and host reductions
  - now keep sampled-spot peak and slope-unit reductions on device in the
    validated calibration path

## Remaining Host Fallbacks and Sync Risks

The main remaining host/device round-trips are:

- some setup-time Pyramid/BioEdge helpers
  - host-built masks / modulation phases copied to device
- some setup-time telescope / detector helpers that materialize host buffers
- some tomography setup/preparation stages
  - sparse/operator assembly is still CPU-originating even though the final
    builder outputs now stay on the selected backend

The main remaining synchronization risks are:

- planned FFT calls, where backend execution is correct but still worth
  profiling for implicit synchronization cost
- calibration/build paths that mix CPU-originating metadata with device arrays
- any future expansion of tomography builders beyond the currently validated
  compact cases

## Priority Assessment

The most important remaining blockers for a truly device-resident runtime loop
are:

1. distinguishing setup-time host copies from true runtime fallbacks
2. optional cleanup of setup-time host-built masks / phases
3. profiling and auditing hidden synchronization in builder-heavy RTC workflows

The remaining setup-time host copies are lower priority than the runtime
fallbacks that have now been removed.

## Recommended Next Steps

1. Run a size-sweep benchmark to find the CPU/GPU crossover points for runtime
   stepping and tomography/reconstructor builds.
2. Re-audit setup-time Pyramid/BioEdge mask and modulation builders only if
   startup cost becomes a practical problem.
3. Keep the maintained CUDA validation set (`gpu_smoke_cuda.jl`,
   `gpu_builder_cuda.jl`, and `gpu_hil_cuda.jl`) green on `spiders`.
4. Use `gpu_sync_audit_cuda.jl` to track the timing surface of runtime and
   builder-heavy RTC paths as the CUDA implementation evolves.

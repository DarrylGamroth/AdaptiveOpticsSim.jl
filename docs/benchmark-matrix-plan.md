# Benchmark Matrix Plan

Status:
- [~] Phase 2 partially implemented
  - `scripts/profile_zernike_runtime.jl` now supports `compact`, `medium`, and
    `representative`
  - `scripts/profile_pixel_output_runtime.jl` now supports `compact`,
    `medium`, and `representative`
  - `scripts/profile_multi_source_multi_wfs_runtime.jl` now supports
    `compact`, `medium`, and `representative`
  - `scripts/profile_revolt_hil_runtime.jl` now covers a REVOLT-like
    synthetic pixel-output SH HIL case with a fixed representative
    configuration and detector-family / detector-response / detector-thermal
    sweeps
  - `scripts/profile_external_optics_hil.jl` now covers a representative
    external-optics HIL case where AdaptiveOpticsSim owns DM phase generation
    and phase export but not the downstream image formation
  - `scripts/profile_ao3k_runtime.jl` now covers the maintained AO3k pyramid
    runtime surface with detector-response, HgCdTe-array sampling-mode,
    readout-correction, and detector-thermal sweeps
- [ ] representative closed-loop and builder ladders are still pending beyond
  these first runtime surfaces
- [x] Phase 5 cross-package harness implemented
  - scenario contracts now live in
    [`../benchmarks/contracts/cross_package.toml`](../benchmarks/contracts/cross_package.toml)
  - the maintained runner is
    [`../scripts/run_cross_package_benchmarks.jl`](../scripts/run_cross_package_benchmarks.jl)
  - archived evidence now lives under
    [`../benchmarks/results/cross_package`](../benchmarks/results/cross_package)
- [x] platform-strengthening Phase 2 evidence added
  - tomography representative benchmark scope now has an archived decision
    record under
    [`../benchmarks/results/tomography`](../benchmarks/results/tomography)
  - a Julia-native SPECULA-informed platform/runtime artifact now lives under
    [`../benchmarks/results/platform`](../benchmarks/results/platform)

## Cross-Package Harness

Maintained cross-package benchmarking is now a first-class engineering surface.

Primary files:

- contract:
  [`../benchmarks/contracts/cross_package.toml`](../benchmarks/contracts/cross_package.toml)
- runner:
  [`../scripts/run_cross_package_benchmarks.jl`](../scripts/run_cross_package_benchmarks.jl)
- archive:
  [`../benchmarks/results/cross_package`](../benchmarks/results/cross_package)
- harness guide:
  [`./cross-package-benchmark-harness.md`](./cross-package-benchmark-harness.md)

The current maintained ladder is:

- `compact`
  - frozen OOPAO and SPECULA reference fidelity bundles
- `medium`
  - REVOLT-like SH HIL runtime comparison between `main` and
    `../AdaptiveOpticsSim.jl-revolt-real`
- `representative`
  - PWFS REVOLT-aligned contract recorded and intentionally deferred until the
    main-repo scenario runner is normalized

These runs are deliberately outside `Pkg.test()` and should be treated as
maintained benchmark evidence rather than unit-test gates.

## Goal

Separate the package's performance surfaces into two explicit classes:

- **Regression benchmarks**
  - small, deterministic, cheap to run
  - protect hot-path latency/allocation behavior and backend correctness
  - suitable for frequent local use and maintained audit workflows
- **Representative benchmarks**
  - medium and large configurations intended to approximate real AO workloads
  - answer CPU vs GPU crossover, scaling, and algorithm-appropriateness
    questions
  - suitable for engineering decisions, not just regression detection

The package already has strong regression-sized timing coverage. What is missing
is a consistent representative benchmark matrix that lets us say when a given
algorithm or backend becomes appropriate for realistic systems.

## Problem Statement

Many current maintained runtime profiles are intentionally compact:

- compact `ClosedLoopRuntime`
- compact AO188 simulation slices
- compact stacked multi-source / multi-WFS profiles
- compact `ZernikeWFS` runtime profile

Those are valuable because they are stable and fast, but they bias the observed
story toward:

- CPU looking much stronger than GPU
- per-step latency dominating throughput considerations
- tiny detector/FFT/gather workloads dominating conclusions

That is the right story for regression guards, but it is not enough to judge:

- whether a GPU-oriented implementation is appropriate at realistic scale
- whether a runtime algorithm scales well with pupil/detector size
- whether a sensor model is structurally suitable for HIL workloads
- where the real CPU/GPU crossover occurs for maintained workflows

## Benchmark Classes

### Class A: Regression Benchmarks

Purpose:

- catch performance regressions
- protect zero-allocation or low-allocation hot paths
- keep backend smoke/audit scripts cheap enough for routine use

Properties:

- deterministic
- warm-started
- small enough to run quickly on developer machines
- stable across hosts and backend versions

Examples already in place:

- `scripts/profile_zernike_runtime.jl`
- `scripts/profile_mixed_sh_asterism_runtime.jl`
- `scripts/profile_multi_source_multi_wfs_runtime.jl`
- `scripts/profile_pixel_output_runtime.jl`
- `scripts/profile_atmosphere_runtime.jl`
- `scripts/profile_atmospheric_field_runtime.jl`
- `scripts/ao188_3k_hil_audit.jl`
- `benchmarks/benchmark_{cpu,cuda,amdgpu}.jl`
- GPU smoke now also guards the maintained monochromatic propagation surface:
  `ElectricField`, `FraunhoferPropagation`, and `FresnelPropagation`
- GPU smoke now also guards the maintained coupled atmosphere-to-field
  propagation surface and curvature-through-atmosphere sensing path
- GPU smoke now also guards maintained broad-band diffractive WFS execution for
  `ShackHartmann` and `PyramidWFS` through `SpectralSource`
- GPU smoke now also guards maintained extended-source diffractive WFS
  execution for `ShackHartmann` and `PyramidWFS` through `ExtendedSource`

### Class B: Representative Benchmarks

Purpose:

- approximate realistic optical/control workloads
- locate CPU/GPU crossover points
- judge whether an algorithm is appropriate for deployment-scale simulations
- provide a better basis for optimization priorities

Properties:

- use physically plausible medium/large configurations
- may be slower and more hardware-sensitive
- should be run intentionally, not on every quick validation pass
- should record both throughput and shape/context, not just a single time

These benchmarks should be treated as engineering evidence rather than
regression gates.

## Representative Benchmark Families

The representative matrix should cover the main runtime and builder families,
not every public function.

### 1. Pixel-Output WFS Runtime

Families:

- Shack-Hartmann
- Pyramid
- BioEdge

### 0. Atmosphere Evolution Runtime

Families:

- finite periodic moving-screen atmosphere
- infinite boundary-injection atmosphere

Purpose:

- make the fidelity/performance tradeoff explicit
- compare CPU and GPU execution for the same multilayer transport workload
- keep the fast/HIL finite path visible even after the infinite model lands

Maintained surface:

- `scripts/profile_atmosphere_runtime.jl`
- `scripts/profile_atmospheric_field_runtime.jl`

Required outputs:

- finite and infinite build/precompute time
- finite and infinite per-step mean and p95
- finite and infinite steady-state allocations
- finite and infinite synchronization count per sample
- explicit infinite/finite slowdown ratio
- geometric and layered-Fresnel atmosphere-aware field propagation timing
- curvature-through-atmosphere timing on the same backend ladder
- ZernikeWFS

Questions:

- how do frame size and pupil size change CPU/GPU crossover?
- when does GPU stop being dominated by launch/sync overhead?
- which sensing families scale most cleanly on CUDA and AMDGPU?

Scale ladder:

- `compact`
  - current regression-sized surface
- `medium`
  - larger pupil and detector sizes where GPU may start to amortize overhead
  - intended as a scaling/crossover rung, not automatically a deployment-scale
    configuration
- `representative`
  - realistic single-WFS closed-loop size meant to approximate deployment use

### 2. Multi-Source / Multi-WFS Runtime

Families:

- stacked asterism sensing
- grouped `CompositeSimulationInterface` execution
- mixed heterogeneous branch execution

Questions:

- how much does batching improve backend utilization?
- when does multi-branch GPU execution become worthwhile?
- when does grouped execution beat branch-by-branch stepping?

Scale ladder:

- `compact`
  - current maintained stacked profiles
- `medium`
  - more sources and larger detector tiles
  - intended as a scaling/crossover rung, not automatically a final
    instrument-like scene
- `representative`
  - realistic multi-WFS HIL or tomography-facing runtime scene

### 3. Closed-Loop HIL Runtime

Families:

- compact maintained runtime
- Subaru AO188 simulation example
- fixed representative HIL cases with large pixel outputs
- later representative SCAO / woofer-tweeter / MCAO examples

Questions:

- what is the realistic per-frame rate for CPU vs GPU?
- how much does latency staging or multi-rate structure matter?
- when do runtime interface choices become the limiting factor?

Scale ladder:

- `compact`
  - existing minimal runtime surfaces
- `medium`
  - AO188-class example
  - larger and more informative than compact, but not necessarily the final
    representative endpoint for every system family
- `representative`
  - system-level examples with realistic detector sizes and control paths

For this family, some representative systems are better modeled as fixed named
cases instead of generic scale rungs. Two maintained examples now exist:

- `revolt_like_hil`
  - REVOLT-like synthetic SH HIL benchmark
  - `277` active commands mapped through a synthetic sparse extrapolation
    operator
  - applied onto a synthetic DM277-like `19 x 19` active-actuator layout
  - benchmark assets live under `benchmarks/assets/revolt_like/`
  - `16 x 16` subapertures
  - `22 x 22` ROI
  - full `352 x 352` pixel output mosaic
- `external_optics_hil`
  - external-propagator / PROPER-style proxy
  - `468` active commands mapped onto a `24 x 24` DM grid
  - full DM OPD generation on a `640 x 640` phase grid
  - exported `640 x 512` phase crop for downstream image formation

### 4. Calibration / Builder Workloads

Families:

- interaction matrices
- modal basis construction
- LiFT
- model tomography

Questions:

- where does GPU already have clear value?
- which builders should default to CPU-build / GPU-run?
- which calibration paths remain too host-oriented?

Scale ladder:

- `compact`
  - current audit/builder surfaces
- `medium`
  - current tomography profiling surfaces
  - intended to expose scaling trends and crossover onset before the fully
    representative builder cases
- `representative`
  - KAPA-like or instrument-scale assembly workloads

## Required Reporting

Every representative benchmark entry should report:

- backend
- configuration class (`compact`, `medium`, `representative`)
- key problem dimensions
  - pupil resolution
  - detector frame shape
  - number of sources
  - number of WFS
  - command length / slope length where relevant
- warmed mean and p95 timing
- effective frame rate or throughput
- phase timing when available
- steady-state allocation bytes where relevant
- response/readout metadata where the benchmark intentionally compares
  detector-response or detector-family variants

Representative reports should avoid naked timings without configuration context.

## Initial Representative Configurations

The first pass does not need instrument-perfect models for every family. It
does need consistent, defensible scale points.

### ZernikeWFS

- `compact`
  - current `16 x 16` / `4 x 4` frame maintained runtime profile
- `medium`
  - increase pupil resolution and sampled frame enough to make GPU execution
    nontrivial
- `representative`
  - a realistic detector-backed closed-loop size intended for RTC-facing use

### Shack-Hartmann

- `compact`
  - current maintained diffractive runtime and mixed-asterism profiles
- `medium`
  - larger detector-backed NGS and LGS cases
- `representative`
  - AO188-class or larger realistic pixel-output closed-loop profile
  - may also be a fixed named HIL case such as REVOLT when the defining feature
    is the detector payload rather than only the command length

### Pyramid / BioEdge

- `compact`
  - current maintained stacked/runtime surfaces
- `medium`
  - larger pupil/reimaged-pupil grids
- `representative`
  - realistic closed-loop detector-backed profiles

### Tomography / Builders

- `compact`
  - current maintained builder surfaces
- `medium`
  - current model-tomography profiling sizes
- `representative`
  - KAPA-like or instrument-scale cases already used for reference validation
  - currently deferred as a routine maintained artifact; see
    [tomography-benchmark-scope.md](./tomography-benchmark-scope.md)

## Backend Questions This Matrix Should Answer

The matrix exists to answer specific engineering questions:

1. At what sizes does GPU become worthwhile for each sensing family?
2. Are CUDA and AMDGPU scaling similarly, or do they diverge by family?
3. Which algorithms are only acceptable as compact regression paths, and which
   remain appropriate at representative scale?
4. Which performance conclusions are genuinely structural, and which are merely
   artifacts of tiny workloads?

If a benchmark does not help answer one of those questions, it should probably
remain a regression benchmark instead of entering the representative matrix.

## Rollout Plan

### Phase 1: Classification

- tag current maintained scripts as regression or representative
- document which current runtime profiles are intentionally compact
- stop using compact profiles as implicit evidence for deployment-scale
  conclusions

### Phase 2: First Representative Runtime Ladder

- add `compact` / `medium` / `representative` modes for:
  - `ZernikeWFS`
  - Shack-Hartmann pixel-output runtime
  - stacked multi-source / multi-WFS runtime

This should be done by extending current profiling scripts rather than creating
entirely new ad hoc entry points where possible.

Current implementation status:

- done for the first three maintained runtime profilers listed above
- still pending for broader closed-loop and builder/runtime families

### Phase 3: Representative Closed-Loop and Builder Cases

- add representative AO188-class and larger closed-loop scenarios
- add representative tomography/calibration benchmark cases
- record CPU, AMDGPU, and CUDA snapshots where available

### Phase 4: Policy Integration

Use the resulting evidence to guide:

- which paths are CPU-first by design
- which paths justify GPU-specific optimization
- which defaults should prefer CPU-build / GPU-run
- which algorithms need redesign before they are suitable at realistic scale

## Success Criteria

This plan is complete when:

1. The package clearly distinguishes regression vs representative benchmarks.
2. Each major runtime family has at least one representative scale point.
3. CPU/GPU crossover conclusions are based on size ladders rather than a single
   compact case.
4. Performance discussions in docs and reviews can cite representative evidence
   instead of inferring too much from tiny workloads.

## Recommendation

Do not replace the compact regression surfaces. Keep them and label them
explicitly. Then add a small, curated representative matrix on top, starting
with:

1. `ZernikeWFS`
2. Shack-Hartmann pixel-output runtime
3. stacked multi-source / multi-WFS runtime
4. representative tomography/build cases

That is the smallest useful matrix that will materially improve performance
decision-making for the package.

## Current Baselines

Snapshot date:

- 2026-03-30

Validated surfaces:

- local CPU representative ladders
- local AMDGPU realistic `medium` ladders
- remote CUDA validation on `spiders`

### AO3k Runtime

Commands:

- `julia --project=. scripts/profile_ao3k_runtime.jl cpu medium default`
- `julia --project=. scripts/profile_ao3k_runtime.jl cpu representative default`
- `julia --project=. scripts/profile_ao3k_runtime.jl amdgpu medium default`
- on `spiders`: `julia --project=. scripts/profile_ao3k_runtime.jl cuda medium default`

Recorded results:

- CPU `medium`
  - `runtime_step_mean_ns`: `1.76839216e7`
  - `frame_rate_hz`: `56.55`
  - `runtime_alloc_bytes`: `4688`
- CPU `representative`
  - `runtime_step_mean_ns`: `4.21010778e7`
  - `frame_rate_hz`: `23.75`
  - `runtime_alloc_bytes`: `4688`
- AMDGPU `medium`
  - `runtime_step_mean_ns`: `7.2786754e6`
  - `frame_rate_hz`: `137.38763511833486`
  - `runtime_alloc_bytes`: `752920`
- CUDA `medium` on `spiders`
  - `runtime_step_mean_ns`: `2.3901578e6`
  - `frame_rate_hz`: `418.38241809808545`
  - `runtime_alloc_bytes`: `242304`

Configuration context:

- pupil resolution: `160`
- `n_subap`: `32`
- high sensor: `hgcdte_avalanche_array`
- high sampling mode: `correlated_double_sampling`
- high readout correction: `reference_pixel_common_mode`

### Multi-Source / Multi-WFS Runtime

Commands:

- `julia --project=. scripts/profile_multi_source_multi_wfs_runtime.jl cpu medium`
- `julia --project=. scripts/profile_multi_source_multi_wfs_runtime.jl cpu representative`
- `julia --project=. scripts/profile_multi_source_multi_wfs_runtime.jl amdgpu medium`

Recorded results:

- CPU `medium`
  - `sh_asterism_mean_ns`: `317637.5`
  - `pyramid_asterism_mean_ns`: `5.005580166666667e6`
  - `bioedge_asterism_mean_ns`: `1.2097988e7`
  - `composite_compatible_mean_ns`: `172208.33`
  - `composite_mixed_mean_ns`: `384359.58`
- CPU `representative`
  - `sh_asterism_mean_ns`: `958680.5`
  - `pyramid_asterism_mean_ns`: `3.5360958666666664e7`
  - `bioedge_asterism_mean_ns`: `8.621651366666667e7`
  - `composite_compatible_mean_ns`: `969146.33`
  - `composite_mixed_mean_ns`: `3.3986123333333335e6`
- AMDGPU `medium`
  - `sh_asterism_mean_ns`: `3.9118879416666664e7`
  - `pyramid_asterism_mean_ns`: `3.0214100833333335e6`
  - `bioedge_asterism_mean_ns`: `6.20568325e6`
  - `composite_compatible_mean_ns`: `5.478529775e7`
  - `composite_mixed_mean_ns`: `3.3314810583333332e7`

Interpretation:

- AO3k now has a maintained realistic-size benchmark on CPU, AMDGPU, and CUDA.
- AMDGPU and CUDA both beat the CPU AO3k `medium` rung, but CUDA is currently
  much farther ahead.
- The multi-WFS family is mixed rather than uniformly GPU-favorable: Pyramid
  and BioEdge scale well on AMDGPU at the `medium` rung, while the current
  diffractive SH asterism path is still expensive there.

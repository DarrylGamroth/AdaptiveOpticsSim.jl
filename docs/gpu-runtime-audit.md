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

On a CUDA host, the standard workflow is:

```bash
julia --project=. scripts/gpu_smoke_cuda.jl
julia --project=. scripts/gpu_builder_cuda.jl
julia --project=. scripts/gpu_hil_cuda.jl
julia --project=. scripts/gpu_sync_audit_cuda.jl
julia --project=. scripts/cpu_crossover_sweep.jl
julia --project=. scripts/gpu_crossover_cuda.jl
julia --project=. scripts/gpu_profile_model_tomography_cuda.jl
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

On an AMDGPU host, the standard workflow is:

```bash
julia --project=. scripts/gpu_smoke_amdgpu.jl
julia --project=. scripts/gpu_builder_amdgpu.jl
julia --project=. scripts/gpu_hil_amdgpu.jl
julia --project=. scripts/gpu_sync_audit_amdgpu.jl
julia --project=. scripts/gpu_crossover_amdgpu.jl
julia --project=. scripts/gpu_profile_model_tomography_amdgpu.jl
julia --project=. scripts/gpu_profile_model_tomography_phases_amdgpu.jl
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
    `Float32` surface and checks pre-command AO188 pixels/slopes plus SH/LGS
    detector outputs.
  - `scripts/gpu_runtime_equivalence_high_accuracy_amdgpu.jl` is the maintained
    scientific `Float64` surface and adds a stricter post-command AO188
    equivalence check (`tel_opd`, post-command pixels, and post-command
    slopes).
  - On this host, the post-command AO188 surface is tightly equivalent in
    `Float64`, while `Float32` still has a visible DM-application accumulation
    gap on GPU.
  - Keeping the AO188 runtime reconstructor in mapped two-stage form
    (`modal_reconstructor` followed by `M2C`) is still worth it for the fast
    path: it improves the maintained AO188 runtime throughput on CPU and
    AMDGPU, but it does not remove the stricter post-command `Float32`
    equivalence gap by itself.
  - Current warmed AO188 surrogate rates on the mapped path are about
    `183 Hz` on local CPU, `215 Hz` on local AMDGPU, and `633 Hz` on CUDA on
    `spiders`, while the maintained fast-runtime equivalence scripts still
    pass on AMDGPU and CUDA.

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
- `MultiRTCBoundary` aggregated runtime stepping
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

# GPU Runtime Audit

This document tracks the state of the device-resident execution path after the
runtime sprint and the current CUDA validation workflow.

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

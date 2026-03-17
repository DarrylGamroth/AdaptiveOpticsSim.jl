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

On a CUDA host, the standard workflow is:

```bash
julia --project=. scripts/gpu_smoke_cuda.jl
julia --project=. scripts/gpu_builder_cuda.jl
julia --project=. scripts/gpu_hil_cuda.jl
julia --project=. scripts/gpu_sync_audit_cuda.jl
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

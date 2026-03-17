# GPU Runtime Audit

This document tracks the state of the device-resident execution path after the
runtime sprint and the current CUDA validation workflow.

## Standard CUDA Validation Workflow

The maintained CUDA validation entry points are:

- `scripts/gpu_smoke_cuda.jl`
  - broad runtime/device-resident smoke coverage
- `scripts/gpu_builder_cuda.jl`
  - reconstructor/calibration builder coverage for modal and tomography paths

On a CUDA host, the standard workflow is:

```bash
julia --project=. scripts/gpu_smoke_cuda.jl
julia --project=. scripts/gpu_builder_cuda.jl
```

The `spiders` workstation is the current real-hardware validation host for this
workflow.

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
  - guide-star coordinate generation still originates on CPU arrays before
    backend materialization
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

1. Re-audit setup-time Pyramid/BioEdge mask and modulation builders only if
   startup cost becomes a practical problem.
2. Keep the maintained CUDA validation pair (`gpu_smoke_cuda.jl` and
   `gpu_builder_cuda.jl`) green on `spiders`.
3. Add profiling-based sync/transfer audits for RTC-facing builder paths.

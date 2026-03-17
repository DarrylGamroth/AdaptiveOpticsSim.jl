# GPU Runtime Audit

This document tracks the state of the device-resident execution path after the
runtime sprint.

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

## Runtime Fallbacks Removed

The following runtime host fallbacks have been removed:

- `ClosedLoopRuntime` DM command update
  - previously used a host `eachindex` loop over `AbstractVector`
  - now uses backend-dispatched scalar CPU / KA accelerator paths
- backend Gaussian noise generation in `randn_backend!`
  - previously staged through a host `Array`
  - now uses a backend-native stateless KA kernel
- diffractive Shack-Hartmann centroid extraction
  - previously copied `valid_mask`, `spot_cube`, and `slopes` to host
  - now computes slopes directly on device with a KA kernel

## Remaining Host Fallbacks

The main remaining host/device round-trips are:

- photon noise in `poisson_noise!`
  - accelerator path still copies through a host matrix
- some Shack-Hartmann calibration/reference helpers
  - `sampled_spots_peak!`
  - `host_mask_view`
  - `mean_valid_signal`
- some setup-time Pyramid/BioEdge helpers
  - host-built masks / modulation phases copied to device
- some setup-time telescope / detector helpers that materialize host buffers

## Priority Assessment

The most important remaining blockers for a truly device-resident runtime loop
are:

1. backend-native Poisson noise
2. detector-noise paths that depend on it
3. calibration/reference helpers if they need to run on device at runtime

The setup-time host copies are lower priority than runtime noise generation.

## Recommended Next Steps

1. Implement a backend-native Poisson-noise kernel.
2. Re-audit detector capture paths with noise enabled after that change.
3. Decide whether calibration/reference helpers must be fully device-resident or
   can remain setup-only host-assisted paths.

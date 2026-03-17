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

## Remaining Host Fallbacks

The main remaining host/device round-trips are:

- some setup-time Pyramid/BioEdge helpers
  - host-built masks / modulation phases copied to device
- some setup-time telescope / detector helpers that materialize host buffers

## Priority Assessment

The most important remaining blockers for a truly device-resident runtime loop
are:

1. distinguishing setup-time host copies from true runtime fallbacks
2. optional cleanup of setup-time host-built masks / phases
3. broader GPU coverage for less-common workflows not yet in the smoke matrix

The remaining setup-time host copies are lower priority than the runtime
fallbacks that have now been removed.

## Recommended Next Steps

1. Re-audit setup-time Pyramid/BioEdge mask and modulation builders only if
   startup cost becomes a practical problem.
2. Add GPU smoke coverage for less-common workflows, especially diffractive
   Shack-Hartmann asterism variants if they matter in practice.
3. Keep separating true runtime fallbacks from acceptable setup-time host work.

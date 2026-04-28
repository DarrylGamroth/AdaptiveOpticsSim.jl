# AMDGPU SH Convergence Plan

Status: active, shared batched SH default recovered

## Purpose

This plan narrows the AMDGPU performance gap on maintained diffractive
Shack-Hartmann runtime surfaces by converging AMDGPU toward the CUDA sensing
path.

The target architecture is:

- shared CPU/CUDA/AMDGPU Shack-Hartmann orchestration,
- shared batched sensing/dataflow by default,
- and narrow ROCm-specific kernels or fallbacks only where ROCm still needs
  them.

This is specifically not a plan to keep the current blanket ROCm-safe path
indefinitely.

For the detailed redesign of the GPU centroid/export algorithm shape itself,
including launch-overhead guardrails and shared CUDA/AMDGPU work, see
[gpu-sh-centroid-redesign-plan.md](./gpu-sh-centroid-redesign-plan.md).

## Current Evidence

On the maintained HEART Julia backend runtime surface:

- CPU `sense_mean_ns`: about `16.0 ms`
- AMDGPU `sense_mean_ns`: about `37.2 ms`
- CUDA `sense_mean_ns`: about `1.6 ms`

The slowdown is concentrated in sensing, not atmosphere or DM application.

The remaining structural causes are:

1. AMDGPU still selects `DetectorHostMirrorPlan`.
2. AMDGPU still uses a CUDA-tolerable but not yet ROCm-optimized centroid
   kernel shape inside the shared batched SH path.
3. Detector capture/finalization still carries ROCm-specific host-mirror
   behavior on maintained SH surfaces.
4. Some non-SH ROCm reductions still require host-mirror fallback.

## Milestones

### SHP-1: Explicit SH Sensing Plan Seam

Goal:

- replace backend-name conditionals with an explicit Shack-Hartmann sensing
  execution plan

Outputs:

- `ShackHartmannWFSScalarPlan`
- `ShackHartmannWFSBatchedPlan`
- `ShackHartmannWFSRocmSafePlan`
- AMDGPU extension override for the ROCm-safe plan
- backend plan tests pinning CUDA vs AMDGPU selection

Acceptance:

- no behavior change on validated CPU/CUDA surfaces
- AMDGPU-specific fork is explicit in code

Status:

- started in this pass

### SHP-2: Device-Native AMDGPU Centroid Extraction

Goal:

- remove per-subap host copies for centroid extraction on the maintained SH
  sensing path

Work:

- add a ROCm-safe centroid kernel operating on the sampled/exported spot cube
- keep existing CPU centroid logic as the reference implementation
- preserve current tolerances against CPU/CUDA

Acceptance:

- no host `centroid_host` copy in the maintained AMDGPU SH sensing path
- reduced `sense_mean_ns` and `sense_alloc_bytes`

Current status:

- first direct recovery attempt was rejected before the ROCm toolchain update
- after the ROCm toolchain update, the shared batched centroid path became
  stable again on the maintained HEART surface
- the remaining issue is performance tuning, not compiler failure

### SHP-3: Batched AMDGPU SH Measurement Path

Goal:

- move maintained AMDGPU diffractive SH back to the batched measurement path

Work:

- route AMDGPU through `ShackHartmannWFSBatchedPlan` for maintained surfaces
- add targeted ROCm-specific kernels only for failing/regressing operations
- keep the ROCm-safe plan as an escape hatch until parity is proven

Acceptance:

- maintained HEART null/noise-free CPU vs AMDGPU equivalence still holds
- AMDGPU `sense_mean_ns` materially improves

Current status:

- completed for the maintained HEART SH surface
- AMDGPU now defaults back to `ShackHartmannWFSBatchedPlan`
- maintained HEART runtime improved from about `90.8 ms` total / `11.0 Hz` to
  about `39.1 ms` total / `25.6 Hz`
- this is still slower than CPU and much slower than CUDA, so convergence work
  remains active

### SHP-4: AMDGPU Device-Side HEART Mosaic Export

Goal:

- remove host tiling from the maintained HEART `352x352` export path

Work:

- replace AMDGPU host tiling in `REVOLT/Julia` with a ROCm-safe device kernel
- keep the current host-tiled path as fallback until parity is proven

Acceptance:

- maintained HEART CPU vs AMDGPU export-equivalence artifact still passes
- AMDGPU `sense_mean_ns` and `sense_alloc_bytes` improve further

### SHP-5: Detector Fallback Narrowing On Maintained SH Surfaces

Goal:

- remove detector host-mirror dependence where it is no longer necessary for
  this maintained SH path

Work:

- profile the remaining detector stages individually on AMDGPU
- move device-safe stages back to direct execution
- keep only narrowly justified host-mirror operations

Acceptance:

- AMDGPU detector execution plan is no longer a blanket host-mirror default for
  the maintained SH surface

## Execution Order

1. `SHP-1` explicit sensing-plan seam
2. `SHP-2` centroid extraction
3. `SHP-3` batched SH sensing
4. `SHP-5` detector narrowing
5. `SHP-4` HEART export tiling

## First Slice Started

This plan has completed the first recovery arc:

- explicit SH sensing-plan seam exists
- AMDGPU default SH sensing is back on the shared batched path
- detector-path narrowing has begun and already recovered direct batched
  device Poisson noise on AMDGPU

The next work is no longer broad SH recovery. It is targeted narrowing of the
remaining detector and reduction fallbacks.

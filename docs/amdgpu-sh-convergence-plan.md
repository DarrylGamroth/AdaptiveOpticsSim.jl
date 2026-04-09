# AMDGPU SH Convergence Plan

Status: active

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
- AMDGPU `sense_mean_ns`: about `91.0 ms`
- CUDA `sense_mean_ns`: about `1.8 ms`

The slowdown is concentrated in sensing, not atmosphere or DM application.

The main structural causes are:

1. AMDGPU still selects `DetectorHostMirrorPlan`.
2. AMDGPU Shack-Hartmann still routes through a ROCm-safe sensing fork instead
   of the batched GPU sensing path.
3. Accelerator centroid extraction still copies spot images to host.
4. The HEART runner still uses host tiling for AMDGPU final `352x352` export.

## Milestones

### SHP-1: Explicit SH Sensing Plan Seam

Goal:

- replace backend-name conditionals with an explicit Shack-Hartmann sensing
  execution plan

Outputs:

- `ShackHartmannScalarPlan`
- `ShackHartmannBatchedPlan`
- `ShackHartmannRocmSafePlan`
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

- first direct recovery attempt was rejected
- replacing the maintained ROCm-safe host loop with the existing batched
  `sh_spot_centroid_kernel!` caused a GPUCompiler / ROCm segfault on the
  maintained HEART HIL path
- a narrower per-spot stats kernel compiled but regressed runtime due to
  per-subap kernel-launch overhead
- the kept baseline remains the host-staged ROCm-safe centroid path until a
  narrower ROCm-specific centroid kernel is proven safe and faster

### SHP-3: Batched AMDGPU SH Measurement Path

Goal:

- move maintained AMDGPU diffractive SH back to the batched measurement path

Work:

- route AMDGPU through `ShackHartmannBatchedPlan` for maintained surfaces
- add targeted ROCm-specific kernels only for failing/regressing operations
- keep the ROCm-safe plan as an escape hatch until parity is proven

Acceptance:

- maintained HEART null/noise-free CPU vs AMDGPU equivalence still holds
- AMDGPU `sense_mean_ns` materially improves

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
4. `SHP-4` HEART export tiling
5. `SHP-5` detector narrowing

## First Slice Started

This pass begins `SHP-1` by introducing an explicit Shack-Hartmann sensing plan
interface and routing AMDGPU through an extension-defined ROCm-safe plan rather
than a backend-name check in core code.

# GPU SH Centroid Redesign Plan

Status: active

Plan traceability:

- [amdgpu-sh-convergence-plan.md](./amdgpu-sh-convergence-plan.md)
- [rocm-failure-catalog.md](./rocm-failure-catalog.md)
- [backend-validation-guide.md](./backend-validation-guide.md)

## Purpose

This plan redesigns the diffractive Shack-Hartmann centroid and export path so
it is genuinely GPU-friendly on both CUDA and AMDGPU.

The current "batched" path is only batched at the API level. Its centroid stage
still behaves like:

- one logical thread per spot,
- one long serial loop over all spot pixels,
- thresholding, reduction, and slope finalization fused into one kernel.

CUDA tolerates that shape. AMDGPU does not perform well on it, even though the
updated ROCm toolchain now appears to compile and run that path correctly.

The redesign goal is:

- keep all maintained SH sensing/export buffers on device,
- increase pixel-level parallelism within each lenslet spot,
- minimize unnecessary kernel launches,
- and share as much of the redesigned GPU path as possible across CUDA and
  AMDGPU.

## Current Evidence

On the maintained HEART null/noise-free SH surface:

- kept AMDGPU ROCm-safe path:
  - `sense_mean_ns`: about `8.54e7`
  - `sense_alloc_bytes`: about `1.05e7`
- forced shared batched SH path on AMDGPU:
  - numerically equivalent to CPU
  - `sense_mean_ns`: about `1.07e8`
  - `sense_alloc_bytes`: about `5.70e6`
- CUDA on the same maintained surface remains much faster

Interpretation:

- the ROCm compiler blocker for the shared centroid kernel is no longer the
  primary issue
- the current batched centroid algorithm is simply not a good GPU kernel shape
- AMDGPU is more sensitive to that kernel shape than CUDA
- reducing allocations alone is not enough if the kernel remains
  under-parallelized

## Design Rules

1. Optimize the algorithmic kernel shape first, not just the backend fork.
2. Keep the maintained SH sensing/export surface device-resident.
3. Prefer shared CUDA/AMDGPU orchestration with backend-tunable launch
   parameters.
4. Split kernels only when the split is justified by parallelism or correctness.
5. Treat kernel-launch overhead as a first-class constraint.

This means:

- do not keep the current one-thread-per-spot serial inner-loop design
- do not blindly split into many tiny kernels if launch overhead erases gains
- do not reintroduce host staging unless it is the only stable fallback

## Target Architecture

### Shared Device Pipeline

For maintained diffractive SH sensing:

1. sampled spot cube already exists on device
2. per-spot threshold/statistics are computed on device
3. slopes are finalized on device
4. exported SH cube remains on device
5. HEART `352x352` mosaic is tiled on device

### Centroid Strategy

The intended centroid shape is:

- one workgroup per spot
- threads cooperate over the `n_pix_subap x n_pix_subap` pixels
- local/shared memory reductions produce:
  - `peak`
  - `total`
  - `sx`
  - `sy`

This is fundamentally different from the current path:

- current path: serial pixel loop inside each thread
- target path: parallel pixel work within each spot

## Kernel Options

### Option A: Two-Pass Workgroup Reduction

Per spot:

1. kernel A computes `peak`
2. kernel B computes thresholded `total`, `sx`, `sy`
3. kernel C finalizes slopes

Pros:

- simplest reduction logic
- easiest to validate

Cons:

- more launch overhead

### Option B: Single Workgroup Kernel With Intra-Group Barrier

Per spot:

1. reduce `peak` into local memory
2. synchronize workgroup
3. reduce thresholded `total`, `sx`, `sy`
4. write one stats record
5. optional separate finalize kernel for slopes

Pros:

- fewer launches
- best chance to beat the current ROCm-safe path

Cons:

- more complex kernel
- may need backend-specific tuning for workgroup size or local memory use

### Preferred Direction

Start with `Option B` if it is stable enough, but keep `Option A` as the
fallback if launch coordination or compiler/runtime behavior becomes fragile.

## Kernel Launch Guardrails

Kernel launch overhead may dominate if this is over-split.

So each implementation attempt must answer:

1. does the split increase useful parallel work enough to justify its launches?
2. does a fused workgroup kernel outperform the split version?
3. is the launch count still acceptable relative to:
   - the current ROCm-safe path
   - the current CUDA path

Required measurements for each attempt:

- `sense_mean_ns`
- `sense_alloc_bytes`
- launch count per SH sensing step
- time spent in:
  - centroid/stat pass
  - slope finalize pass
  - export tiling pass

Stop condition:

- reject any redesign that is numerically correct but slower than the kept
  ROCm-safe path without a clear follow-on path to recover the loss

## Proposed Milestones

### GCP-1: Isolate GPU SH Stats/Finalize Phases

Goal:

- separate the current centroid implementation into explicit stats and finalize
  responsibilities without changing results

Outputs:

- explicit internal interfaces for:
  - spot stats accumulation
  - slope finalization

Acceptance:

- CPU/CUDA behavior unchanged
- no regression in maintained parity tests

### GCP-2: Workgroup-Parallel GPU Stats Kernel

Goal:

- replace the serial per-spot centroid loop with a workgroup-parallel spot
  stats kernel

Outputs:

- one workgroup-per-spot stats kernel
- preallocated stats buffers reused across calls

Acceptance:

- CPU parity on maintained SH export surface
- no host centroid staging on the GPU path under test

### GCP-3: Shared CUDA/AMDGPU GPU Centroid Path

Goal:

- run both CUDA and AMDGPU through the redesigned stats/finalize path

Work:

- start from shared kernel logic
- tune workgroup size by backend if needed
- keep backend-specific special cases narrow and explicit

Acceptance:

- CUDA parity still holds
- AMDGPU parity holds
- at least one GPU backend benefits materially versus the old batched kernel

### GCP-4: Device-Side HEART Mosaic Export

Goal:

- remove remaining host tiling from the maintained HEART export path

Acceptance:

- maintained HEART CPU-vs-GPU export-equivalence still passes
- launch count and runtime do not regress

### GCP-5: AMDGPU Default-Path Decision

Goal:

- decide whether AMDGPU can move from `ShackHartmannWFSRocmSafePlan` to the shared
  GPU SH path by default

Decision rule:

- keep the shared path only if all are true:
  - stable
  - numerically equivalent
  - faster than the kept ROCm-safe path on the maintained HEART surface

If not:

- keep the ROCm-safe default
- retain the redesigned path as an experimental or backend-specific option

## Validation Matrix

Each milestone should be validated on:

- CPU reference
- CUDA maintained SH parity surface
- AMDGPU maintained SH parity surface
- HEART `277 -> 352x352` null/noise-free export
- HEART runtime profile

Required checks:

- `sh_exported_spot_cube(...)`
- `wfs_output_frame(...)`
- HEART tiled frame equivalence
- backend runtime profile

## Immediate Next Step

Implement `GCP-1` and `GCP-2` together:

- introduce explicit stats/finalize seams in the SH centroid path
- prototype a workgroup-parallel stats kernel behind an internal execution-plan
  branch
- compare against:
  - current CUDA batched path
  - current AMDGPU ROCm-safe path
  - current forced batched AMDGPU path

The first success criterion is not "beat CUDA."
It is:

- match CPU numerically,
- beat the kept AMDGPU ROCm-safe sensing runtime,
- and do so without exploding launch overhead.

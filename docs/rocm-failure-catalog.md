# ROCm Failure Catalog

Status: active

Plan traceability:

- [rocm-fallback-inventory.md](./rocm-fallback-inventory.md)
- [amdgpu-sh-convergence-plan.md](./amdgpu-sh-convergence-plan.md)
- [backend-validation-guide.md](./backend-validation-guide.md)

## Purpose

This note records the specific ROCm / AMDGPU failure surfaces that motivated
the current AMDGPU-specific fallbacks and execution-plan overrides.

This is not a list of theoretical incompatibilities.

It is a list of maintained surfaces where one or more of the following
happened on ROCm-backed execution:

- compiler failure,
- runtime crash or segfault,
- numerically incorrect result,
- or severe runtime/allocation regression large enough to reject the direct
  path for maintained use.

The goal is to give future refactor work, `SKILL.md`, or `AGENTS.md`
instructions a concrete source of truth for what is known to fail, what is only
suspected, and what has already been recovered.

## Interpretation Rule

Read each entry as:

- "this exact surface failed or regressed on maintained ROCm paths"
- not "the whole CUDA-style implementation is impossible on ROCm"

The intended end state is still:

- shared accelerator path with CUDA,
- plus narrow ROCm-specific kernels or overrides only where needed.

## Catalog

### 1. Detector direct execution for maintained frame finalization surfaces

Files:

- [AdaptiveOpticsSimAMDGPUExt.jl](../ext/AdaptiveOpticsSimAMDGPUExt.jl)
- [pipeline.jl](../src/detectors/pipeline.jl)
- [frame_capture.jl](../src/detectors/frame_capture.jl)

Current fallback:

- AMDGPU selects `DetectorHostMirrorPlan`

Known failing or rejected direct surfaces:

- Poisson frame noise on maintained ROCm detector paths
- readout-correction execution on maintained ROCm detector paths
- earlier detector finalization surfaces involving saturation / quantization on
  GPU-backed frame buffers

Observed symptom classes:

- runtime crashes / segfaults in maintained detector capture/finalization paths
- ROCm / GPUCompiler IR-validation failures on detector-path kernels

Current status:

- Gaussian frame noise has already been recovered to device-native execution on
  AMDGPU
- batched Poisson noise has also been recovered to device-native execution on
  AMDGPU for maintained SH detector capture
- readout-correction and some scalar frame-finalization operations remain on
  the host-mirror path

### 2. Direct ROCm reductions used by maintained SH and grouped runtime paths

Files:

- [reductions.jl](../src/core/reductions.jl)
- [AdaptiveOpticsSimAMDGPUExt.jl](../ext/AdaptiveOpticsSimAMDGPUExt.jl)

Current fallback:

- AMDGPU selects `HostMirrorReductionPlan`

Known failing or rejected direct surfaces:

- direct ROCm `maximum` in the maintained Shack-Hartmann path
- shared masked-sum / packed-pair reduction surfaces in maintained grouped
  Pyramid/BioEdge calibration/runtime paths

Observed symptom classes:

- segfaults on direct ROCm reduction surfaces
- instability on grouped calibration/runtime reductions

Current status:

- no maintained device-native AMDGPU replacement yet

### 3. Blanket CUDA-style batched Shack-Hartmann sensing on ROCm

Files:

- [setup.jl](../src/wfs/shack_hartmann/setup.jl)
- [signals.jl](../src/wfs/shack_hartmann/signals.jl)
- [measure.jl](../src/wfs/shack_hartmann/measure.jl)
- [stacks.jl](../src/wfs/shack_hartmann/stacks.jl)

Current maintained path:

- AMDGPU now selects `ShackHartmannWFSBatchedPlan`
- `ShackHartmannWFSRocmSafePlan` remains available as an escape hatch, not the
  default

Why this exists:

- the fully batched accelerator sensing path that works on CUDA was not kept as
  the maintained AMDGPU path because multiple ROCm-sensitive sub-surfaces were
  not reliable enough

Known failing or rejected sub-surfaces:

- direct ROCm-safe replacement for all SH reductions/centroid extraction was
  not in place
- several maintained SH sensing surfaces depended on host-side mirrors or safe
  loops to avoid ROCm fragility

Observed symptom classes:

- compiler/runtime fragility on portions of the fully batched SH path
- performance collapse large enough to reject naive direct-path experiments

Important nuance:

- this entry does not mean the CUDA-style SH path is impossible on ROCm
- it means the package currently lacks a sufficiently narrow, maintained
  ROCm-specific kernel set to support it safely

Current status:

- an explicit SH sensing execution-plan seam now exists
- after the ROCm toolchain update and detector-path cleanup, the shared
  batched SH path is again stable enough to be the maintained AMDGPU default
- the remaining gap is performance, not basic viability

### 4. Host centroid extraction in ROCm-safe SH sensing recovery work

Files:

- [signals.jl](../src/wfs/shack_hartmann/signals.jl)

Current behavior:

- the old ROCm-safe SH path stages spot information back to host-owned buffers
  or host-owned slope/reference structures
- this is no longer the maintained default, but it remains relevant as the
  fallback baseline and as a reference for future regressions

Observed symptom classes:

- large runtime overhead
- large allocation overhead on maintained HEART HIL SH surfaces

Evidence:

- [2026-04-08-heart-hil-julia-backend-runtime.toml](../../AdaptiveOpticsComparisons/results/archived/2026-04-08-heart-hil-julia-backend-runtime.toml)
  identified `sense_mean_ns` as the dominant AMDGPU cost on the maintained
  HEART surface

Current status:

- this is no longer the maintained default path
- it remains useful as a debug / recovery baseline

Additional observed failure while attempting recovery:

- a direct AMDGPU device-centroid replacement using the existing batched
  `sh_spot_centroid_kernel!` on the maintained ROCm-safe SH signal path
  triggered a real GPUCompiler / ROCm segfault during the maintained HEART HIL
  runtime (`signal 11` during IR validation in `sh_signal_from_spots!`)
- a narrower per-spot stats kernel did compile, but it launched once per
  subaperture and made the maintained HEART path slower than the kept host
  staging path

Implication:

- "device centroid extraction" is still the right target
- but the CUDA batched centroid kernel cannot currently just be dropped into
  the ROCm-safe maintained SH path and assumed safe
- the next attempt should use a narrower ROCm-specific kernel and must be
  validated first on the maintained HEART surface

Focused repro:

- [amdgpu_sh_localmem_repro.jl](../scripts/amdgpu_sh_localmem_repro.jl)
- `serial` mode:
  - one-thread-per-spot cutoff-stats kernel
  - expected to complete successfully on AMDGPU
- `localmem` mode:
  - one-workgroup-per-spot cutoff-stats kernel using `@localmem` and
    `@synchronize`
  - currently expected to trip a ROCm runtime `GPU Kernel Exception`

Additional findings:

- `KernelAbstractions.@atomic` float accumulators also fail on ROCm in this
  environment with `InvalidIRError` during compilation.
- native `@roc` float atomics fail on the same maintained SH surface with the
  same kind of compile-time `InvalidIRError`.
- native `@roc` workgroup-local reductions using `@ROCStaticLocalArray` and
  `sync_workgroup()` do compile and run, so the current ROCm-friendly route is
  backend-native cooperative reduction rather than KA local-memory or atomic
  accumulation.
- on the maintained HEART SH surface, that native cooperative reduction is
  correct but still slower than the current shared batched path, so it is not
  the default path.

### 5. AMDGPU HEART SH export tiling in `REVOLT/Julia`

Files:

- [profile_runtime.jl](../../REVOLT/Julia/support/profile_runtime.jl)

Current behavior:

- AMDGPU uses host staging for final `352x352` Shack-Hartmann mosaic export in
  the maintained HEART runner

Why:

- earlier direct ROCm export / tiling path was not stable enough to keep as the
  maintained implementation

Observed symptom classes:

- correctness recovery required staged host tiling
- remaining performance is not competitive with CPU or CUDA

Current status:

- correctness is maintained
- performance remains poor until this path is replaced with a device-native,
  ROCm-safe export kernel

### 6. Atmosphere phase-noise staging

Files:

- [kolmogorov.jl](../src/atmosphere/kolmogorov.jl)
- [AdaptiveOpticsSimAMDGPUExt.jl](../ext/AdaptiveOpticsSimAMDGPUExt.jl)

Current behavior:

- Gaussian phase-screen fills still route through `randn_phase_noise!` with a
  host staging buffer

Observed symptom classes:

- not currently the dominant runtime cost on maintained surfaces
- retained because the fully direct ROCm fill path has not been justified
  enough to replace the explicit staged path

Current status:

- acceptable temporary fallback

## What Has Already Been Recovered

These surfaces were previously problematic but are now maintained:

- AMDGPU Gaussian detector frame noise is device-native again
- CUDA SH detector/export parity is now explicitly pinned in tests and full GPU
  smoke
- matched HEART CPU vs AMDGPU frame equivalence is maintained on the Julia-side
  comparison export surface

## Guidance For Future Refactor Work

When working on ROCm paths:

1. Do not assume "CUDA path" is categorically unsafe on ROCm.
2. Replace blanket ROCm-safe forks with explicit execution-plan seams first.
3. Attack one failing sub-surface at a time:
   - reductions
   - centroid extraction
   - detector finalization
   - export tiling
4. Require all of:
   - no compiler/runtime failure
   - CPU numerical equivalence on the maintained surface
   - no catastrophic regression versus the kept fallback
5. Keep explicit notes here when a surface is:
   - recovered,
   - newly observed to fail,
   - or still only suspected.

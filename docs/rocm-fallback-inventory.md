# ROCm Fallback Inventory

Status: active

Plan traceability:

- [platform-strengthening-plan.md](./platform-strengthening-plan.md) `PSP-08`
- [execution-plan-closeout.md](./execution-plan-closeout.md)
- [amdgpu-sh-convergence-plan.md](./amdgpu-sh-convergence-plan.md)
- [rocm-failure-catalog.md](./rocm-failure-catalog.md)

## Purpose

This note records the remaining explicit ROCm-specific fallback surfaces after
the execution-plan rollout and the Phase 3 cleanup pass.

For the underlying "what exactly failed on ROCm and why this fallback exists"
detail, see [rocm-failure-catalog.md](./rocm-failure-catalog.md).

The goal is to keep the remaining workaround debt:

- explicit,
- benchmark-backed,
- and classified by whether it is acceptable, temporary, or the next best
  removal target.

## Current Inventory

### 1. Detector host-mirror path for Poisson and readout correction

Files:

- [frame_capture.jl](../src/detectors/frame_capture.jl)
- [pipeline.jl](../src/detectors/pipeline.jl)
- [AdaptiveOpticsSimAMDGPUExt.jl](../ext/AdaptiveOpticsSimAMDGPUExt.jl)

Current behavior:

- AMDGPU still selects `DetectorHostMirrorPlan` for detector execution.
- Poisson frame noise and readout-correction models still copy through the
  detector-owned `noise_buffer_host`.
- Phase 3 reduced this fallback surface by moving Gaussian frame noise back to
  device-native execution on AMDGPU.

Status:

- acceptable temporary maintained fallback

Reason:

- the remaining host-mirror stages sit on known ROCm-sensitive detector paths
  and are still justified by the maintained AO3k benchmark surface

### 2. Host-mirror reduction plan for ROCm peak and masked reductions

Files:

- [reductions.jl](../src/core/reductions.jl)
- [AdaptiveOpticsSimAMDGPUExt.jl](../ext/AdaptiveOpticsSimAMDGPUExt.jl)

Current behavior:

- AMDGPU selects `HostMirrorReductionPlan`.
- `backend_maximum_value`, `masked_sum2d`, and `packed_valid_pair_mean` stay on
  the host-mirror path.

Status:

- acceptable temporary maintained fallback

Reason:

- direct ROCm `maximum` still segfaults in the maintained Shack-Hartmann path
- the shared GPU masked-sum kernel still segfaults in maintained grouped
  Pyramid/BioEdge calibration

### 3. ROCm-specific Shack-Hartmann safe path

Files:

- [setup.jl](../src/wfs/shack_hartmann/setup.jl)
- [signals.jl](../src/wfs/shack_hartmann/signals.jl)
- [measure.jl](../src/wfs/shack_hartmann/measure.jl)

Current behavior:

- diffractive Shack-Hartmann on ROCm still avoids parts of the stacked GPU
  measurement path
- cached host vectors and mask mirrors owned by SH state are used where ROCm
  remains fragile

Status:

- next target for removal

Reason:

- this is the largest remaining ROCm-specific behavioral fork in a maintained
  hot path
- grouped and detector execution plans are now explicit enough that this path
  can be attacked in isolation later

### 4. Atmosphere phase-noise staging

Files:

- [kolmogorov.jl](../src/atmosphere/kolmogorov.jl)
- [AdaptiveOpticsSimAMDGPUExt.jl](../ext/AdaptiveOpticsSimAMDGPUExt.jl)

Current behavior:

- AMDGPU phase-screen Gaussian fills still route through
  `randn_phase_noise!` with a preallocated host matrix

Status:

- acceptable temporary maintained fallback

Reason:

- the fallback is explicit, preallocated, and not currently the dominant cost
  on maintained realistic runtime surfaces

## Phase 3 Outcome

Phase 3 reduced one maintained ROCm host-mirror fallback surface without
changing global backend policy:

- detector Gaussian frame noise no longer uses the detector host mirror on
  AMDGPU

The remaining fallbacks are now narrower and should be treated as explicit
technical debt rather than ambient backend behavior.

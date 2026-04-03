# ROCm Phase 3 Rebaseline

Status: active

Plan traceability:

- [platform-strengthening-plan.md](./platform-strengthening-plan.md) `PSP-09`
- [platform-strengthening-plan.md](./platform-strengthening-plan.md) `PSP-10`

## Purpose

This note records the maintained runtime surfaces re-baselined after the Phase
3 ROCm cleanup.

Phase 3 changed one hot-path behavior only:

- AMDGPU detector Gaussian frame noise moved from the detector-owned host mirror
  back to device-native execution

The remaining detector Poisson/readout-correction host-mirror path and the
remaining ROCm reduction/SH-safe-path workarounds are unchanged and explicitly
tracked in [rocm-fallback-inventory.md](./rocm-fallback-inventory.md).

## Recorded surfaces

Artifact:

- [2026-04-02-phase3-psp10.toml](../benchmarks/results/platform/2026-04-02-phase3-psp10.toml)

Maintained commands:

- local CPU:
  - `julia --project=. --startup-file=no scripts/profile_ao3k_runtime.jl cpu medium default`
  - `julia --project=. --startup-file=no scripts/profile_multi_source_multi_wfs_runtime.jl cpu medium`
- local AMDGPU:
  - `julia --project=. --startup-file=no scripts/profile_ao3k_runtime.jl amdgpu medium default`
  - `julia --project=. --startup-file=no scripts/profile_multi_source_multi_wfs_runtime.jl amdgpu medium`
- remote CUDA on `spiders`:
  - `julia --project=. --startup-file=no scripts/gpu_smoke_cuda.jl`
  - `julia --project=. --startup-file=no scripts/profile_ao3k_runtime.jl cuda medium default`
  - `julia --project=. --startup-file=no scripts/profile_multi_source_multi_wfs_runtime.jl cuda medium`

## Results

### AO3k medium

| Backend | Step mean ns | Alloc bytes |
| --- | ---: | ---: |
| CPU | `1.73412479e7` | `6016` |
| AMDGPU | `7.1657279e6` | `548328` |
| CUDA | `2.5233259e6` | `268440` |

### Multi-source / multi-WFS medium

| Backend | SH ns | Pyramid ns | BioEdge ns | Compatible composite ns | Mixed composite ns |
| --- | ---: | ---: | ---: | ---: | ---: |
| CPU | `304309.6666666667` | `4.886819583333333e6` | `1.144210275e7` | `168559.66666666666` | `365629.5` |
| AMDGPU | `3.906693175e7` | `2.87203475e6` | `6.168384583333333e6` | `5.3198266416666664e7` | `3.4180570583333336e7` |
| CUDA | `223283.41666666666` | `2.7633235e6` | `6.206203083333333e6` | `1.1169828333333333e6` | `1.857524e6` |

## Interpretation

- CPU remained in-family on both realistic surfaces.
- AMDGPU grouped runtime remained in-family after the grouped execution-plan
  recovery work.
- AMDGPU AO3k remained in-family after the Phase 3 detector cleanup, but the
  runtime allocation surface is still materially larger than CPU and CUDA.
- CUDA remained healthy on `spiders`, and the maintained smoke matrix stayed
  green.

## Conclusion

Phase 3 reduced one benchmark-backed ROCm host-mirror fallback without causing
a maintained CPU, AMDGPU grouped-runtime, or CUDA regression.

The next ROCm reduction target should be chosen from:

1. the Shack-Hartmann safe path
2. the remaining detector Poisson/readout-correction host-mirror path
3. atmosphere phase-noise staging

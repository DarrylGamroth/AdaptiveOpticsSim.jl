# Backend Validation Guide

Status: active

Plan traceability:

- [`PLAN-15`](./package-review-action-plan.md)
- [`PLAN-32`](./package-review-action-plan.md)
- review IDs: `PR-31`, `PR-32`, `PR-33`

## Purpose

This document explains how backend validation is organized.

For backend-specific failure history that motivated current fallback plans,
especially on ROCm/AMDGPU, see
[rocm-failure-catalog.md](./rocm-failure-catalog.md).

The goal is to keep four distinct classes of evidence separate:

- functional tests in `Pkg.test()`
- optional backend smoke/parity checks in `Pkg.test()`
- dedicated hardware-backed backend validation targets
- benchmark and profile evidence outside `Pkg.test()`

For current release/support scope, use:

- [supported-production-surfaces.md](./supported-production-surfaces.md)
- [production-readiness-checklist.md](./production-readiness-checklist.md)
- [self-hosted-gpu-runner-setup.md](./self-hosted-gpu-runner-setup.md)

## Test Layout

### Functional and subsystem tests

These run unconditionally in [`runtests.jl`](../test/runtests.jl) through the
grouped testsets under [`test/testsets`](../test/testsets):

- core optics
- atmosphere
- control and runtime
- detectors and WFS
- calibration and analysis
- reference and tutorial regression

These are normal correctness tests, not backend throughput checks.

### Optional backend smoke in `Pkg.test()`

These run after the functional testsets and skip cleanly if the backend package
or device runtime is unavailable:

- [`optional_amdgpu_backends.jl`](../test/optional_amdgpu_backends.jl)
- [`optional_cuda_backends.jl`](../test/optional_cuda_backends.jl)

Shared smoke scaffolding lives in:

- [`backend_optional_common.jl`](../test/backend_optional_common.jl)

The reduced maintained smoke covers:

- multilayer atmosphere
- infinite atmosphere
- atmospheric field propagation
- polychromatic diffractive SH
- deterministic diffractive SH detector/export equivalence against CPU
- deterministic composite-optic low-order runtime parity against CPU:
  - `tiptilt + dm`
    - `ShackHartmannWFS`
    - `Pyramid`
    - `BioEdge`
  - `steering + dm`
  - `focus + dm`
- curvature-through-atmosphere

For the maintained low-order composite surfaces, optional backend smoke now also
checks:

- short command-sequence correctness after `set_command!` / `update_command!`
- command-isolation failure paths
- backend replay determinism through repeated GPU execution on the runtime
  equivalence contracts

For broader backend-audit coverage, use:

- `ADAPTIVEOPTICS_TEST_FULL_AMDGPU=1`
- `ADAPTIVEOPTICS_TEST_FULL_CUDA=1`

which route through [`gpu_smoke_contract.jl`](../scripts/gpu_smoke_contract.jl).

That full matrix is useful for subsystem audit and backend triage, but it does
not by itself define the supported GPU surface. The supported GPU surface is
defined by the dedicated hardware-backed validation targets below.

For explicit hardware-backed backend validation targets that combine the
optional backend smoke with the maintained runtime-equivalence contracts,
including the high-accuracy post-command equivalence pass, run:

- `julia --project=. --startup-file=no test/runtests_amdgpu.jl`
- `julia --project=. --startup-file=no test/runtests_cuda.jl`

These targets are intended for hosts where the corresponding GPU package and
runtime are actually available. They fail fast instead of skipping when the
backend is unavailable. They are the release-gated support boundary for CUDA
and AMDGPU.

The full GPU smoke matrix now also pins the exact batched Shack-Hartmann
detector/export surface that previously regressed on CUDA:

- null-noise diffractive SH with detector capture
- CPU vs GPU comparison of:
  - [`sh_exported_spot_cube`](../src/wfs/shack_hartmann/setup.jl)
  - [`wfs_output_frame`](../src/control/runtime/construction.jl)

This keeps the public exported-pixel surface under backend parity coverage, not
just the slope output.

## Benchmark Separation

Benchmarks do not belong in `Pkg.test()`.

Representative runtime evidence should be gathered with maintained profile or
benchmark scripts such as:

- [`profile_ao3k_runtime.jl`](../scripts/profile_ao3k_runtime.jl)
- [`profile_control_loop_runtime.jl`](../scripts/profile_control_loop_runtime.jl)
- [`profile_multi_source_multi_wfs_runtime.jl`](../scripts/profile_multi_source_multi_wfs_runtime.jl)
- [`run_cross_package_benchmarks.jl`](../scripts/run_cross_package_benchmarks.jl)
- [`benchmarks/`](../benchmarks)

This separation exists so:

- unit tests stay reasonably fast
- optional GPU availability does not block normal development
- benchmark evidence remains intentional and archived rather than implicit

## Automation

Checked-in CI automation now exists in:

- [../.github/workflows/cpu-validation.yml](../.github/workflows/cpu-validation.yml)
- [../.github/workflows/cuda-backend-validation.yml](../.github/workflows/cuda-backend-validation.yml)
- [../.github/workflows/amdgpu-backend-validation.yml](../.github/workflows/amdgpu-backend-validation.yml)

Current intent:

- CPU workflow:
  - runs the normal `Pkg.test()` suite on a hosted runner
- CUDA workflow:
  - targets a self-hosted runner labeled `self-hosted`, `linux`, `cuda`
  - runs the maintained CUDA smoke and runtime-equivalence scripts
  - runtime equivalence includes the composite-optic low-order HIL surfaces:
    - `tiptilt + dm`
      - `ShackHartmannWFS`
      - `Pyramid`
      - `BioEdge`
    - `steering + dm`
    - `focus + dm`
- AMDGPU workflow:
  - targets a self-hosted runner labeled `self-hosted`, `linux`, `amdgpu`
  - runs the maintained AMDGPU smoke and runtime-equivalence scripts
  - runtime equivalence includes the composite-optic low-order HIL surfaces:
    - `tiptilt + dm`
      - `ShackHartmannWFS`
      - `Pyramid`
      - `BioEdge`
    - `steering + dm`
    - `focus + dm`

These workflows are part of production hardening, but they are only fully
effective once the expected self-hosted GPU runners are actually registered and
kept healthy.

## Validation Rule

Claiming backend support should imply:

- functional correctness in the main testsets
- optional backend smoke coverage in `Pkg.test()`
- maintained benchmark or profile evidence where runtime claims are made

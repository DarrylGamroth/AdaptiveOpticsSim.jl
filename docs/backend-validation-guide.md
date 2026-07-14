# Backend Validation Guide

Status: active

## Purpose

This document explains how backend validation is organized.

For backend-specific failure history that motivated current fallback plans,
especially on ROCm/AMDGPU, see
AMDGPU validation history in git history.

The goal is to keep four distinct classes of evidence separate:

- functional tests in `Pkg.test()`
- optional backend smoke/parity checks in `Pkg.test()`
- dedicated hardware-backed backend validation targets
- benchmark and profile evidence outside `Pkg.test()`

For current release/support scope, use:

- [supported-production-surfaces.md](supported-production-surfaces.md)
- [release-validation-runbook.md](release-validation-runbook.md)

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
- MKID accumulated-count capture, source passband handling, and
  flux-conserving channel-crosstalk parity

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

- `julia --project=test/amdgpu --startup-file=no test/runtests_amdgpu.jl`
- `julia --project=test/cuda --startup-file=no test/runtests_cuda.jl`

These targets are intended for hosts where the corresponding GPU package and
runtime are actually available. They fail fast instead of skipping when the
backend is unavailable. The AMDGPU target is the current release-gated
accelerator support boundary. The CUDA target is retained as a fail-fast
restoration target, but is not a current support gate because no CUDA validation
device is available. The backend projects under
[`test/amdgpu`](../test/amdgpu) and [`test/cuda`](../test/cuda) declare the root
package as a path source, keeping accelerator packages out of normal
`Pkg.test()` while still resolving the full checkout dependency graph in a
clean environment.

The full GPU smoke matrix now also pins the exact batched Shack-Hartmann
detector/export surface that previously regressed on CUDA:

- null-noise diffractive SH with detector capture
- CPU vs GPU comparison of:
  - the Shack-Hartmann exported spot-cube path in
    [`shack_hartmann.jl`](../src/wfs/shack_hartmann.jl)
  - [`wfs_output_frame`](../src/control/runtime/construction.jl)

This keeps the public exported-pixel surface under backend parity coverage, not
just the slope output.

### AMDGPU host-mirror boundaries

With AMDGPU 2.7, photon Poisson sampling, detector readout correction, windowed
multi-read products, and float-to-integer detector export use host-mirror
paths. Poisson sampling, correction, and integer export reuse detector host
buffers; windowed readout-product construction retains its existing allocation
behavior. Counting crosstalk and the rest of the maintained detector array math
remain backend kernels.

On the normal point-source and spectral-stack paths, Shack-Hartmann field
formation, FFTs, and spot sampling remain batched on the device. Centroid
cutoff/statistics use the reusable host centroid workspace and copy the
thresholded spot back to the device so exported pixels retain CPU semantics.
Non-stackable asterism fallbacks retain conservative host field construction.
These are correctness boundaries with explicit transfers, so detector-heavy or
high-rate Shack-Hartmann HIL workloads should be profiled on the target AMD
device before making latency claims.

Atmospheric phase-screen PSD construction is a setup-time host-mirror boundary
on AMDGPU 2.7. Geometric Shack-Hartmann slopes also use an explicit host
compatibility path because the variable-trip KA kernel fails GPUCompiler IR
validation on the current gfx1030 target. Gain-sensing-camera scalar reductions
follow the same backend reduction policy. LiFT keeps its FFTs, Jacobian math,
and dense linear algebra on the device, but stages small factor-matrix subviews
and condition checks through dense host storage where AMDGPU's generic ROCArray
gather path is not usable. These boundaries are covered by the dedicated AMDGPU
target and the full smoke matrix; they are not evidence of a fully
device-resident AMD execution plan.

## Benchmark Separation

Benchmarks do not belong in `Pkg.test()`.

Representative runtime evidence should be gathered with maintained profile or
benchmark scripts such as:

- [`profile_ao3k_runtime.jl`](../scripts/profile_ao3k_runtime.jl)
- [`profile_control_loop_runtime.jl`](../scripts/profile_control_loop_runtime.jl)
- [`profile_multi_source_multi_wfs_runtime.jl`](../scripts/profile_multi_source_multi_wfs_runtime.jl)
- [`run_cross_package_benchmarks.jl`](../scripts/run_cross_package_benchmarks.jl)
- [`benchmarks/`](../benchmarks)

The benchmark environments are split so a CPU benchmark does not resolve both
GPU stacks:

```bash
julia --project=benchmarks -e 'using Pkg; Pkg.instantiate()'
julia --project=benchmarks benchmarks/benchmark_cpu.jl
julia --project=benchmarks benchmarks/benchmark_cpu_hotpath_cards.jl

julia --project=benchmarks/amdgpu -e 'using Pkg; Pkg.instantiate()'
julia --project=benchmarks/amdgpu benchmarks/benchmark_amdgpu.jl
```

For the detector-output HIL path, use the same workload and sample count on
both backends:

```bash
julia --project=. scripts/profile_revolt_hil_runtime.jl cpu benchmarks/assets/revolt_like cmos default none 100 10
julia --project=benchmarks/amdgpu scripts/profile_revolt_hil_runtime.jl amdgpu benchmarks/assets/revolt_like cmos default none 100 10
```

The profile defaults to 100 warmed samples and reports p95; it rejects smaller
sample counts because they do not support that percentile usefully. Treat it as
a closed-loop diagnostic, not a fixed-arrival-rate latency contract. Use
[`soak_revolt_hil_runtime.jl`](../scripts/soak_revolt_hil_runtime.jl) for the
100,000-sample p99/GC soak, and retain raw run output when making release or
hardware-capacity claims. A current CPU/AMDGPU summary for the local validation
host is archived in
[`2026-07-13-revolt-hil-cpu-amdgpu.toml`](../benchmarks/results/platform/2026-07-13-revolt-hil-cpu-amdgpu.toml).

The AO188/AO3k CPU profile scripts accept an FFT thread count as their final
argument, or through `AOS_FFT_THREADS`. Tune this per workload and host rather
than treating the machine's core count as a default; small FFTs can become
slower with too many provider threads. For example:

```bash
AOS_FFT_THREADS=4 julia --project=. scripts/profile_ao3k_runtime.jl cpu compact
```

Keep FFT, BLAS, and coarse Julia task parallelism from oversubscribing the same
cores. Deterministic validation continues to use one thread.

For independent simulation sweeps, run the ensemble scheduler benchmark before
selecting `ThreadedExecution`, `AcceleratedKernelsExecution`, or
`DaggerExecution`:

```bash
julia --threads=8 --project=benchmarks benchmarks/benchmark_ensemble_schedulers.jl
```

This is offline throughput evidence, not an external-RTC latency result.
Dagger is intended for task graphs, locality, and process/node scaling; AK is
intended for reusable local task partitioning on sufficiently large many-core
workloads. Keep direct sequential runtime stepping as the HIL baseline.
The current local-host comparison is archived in
[`2026-07-13-ensemble-schedulers.toml`](../benchmarks/results/platform/2026-07-13-ensemble-schedulers.toml);
rebaseline it on an EPYC or Threadripper target before enabling a site policy.

The retained CUDA benchmark has the equivalent isolated
`--project=benchmarks/cuda` environment, but requires restored CUDA hardware.

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
  - runs the isolated AcceleratedKernels/Dagger scheduler extension tests on a
    four-thread Linux job
- CUDA workflow:
  - targets a self-hosted runner labeled `self-hosted`, `linux`, `cuda`
  - instantiates [`test/cuda`](../test/cuda)
  - is dormant until a CUDA runner is restored, then runs the retained CUDA
    hardware target
  - runtime equivalence includes the composite-optic low-order HIL surfaces:
    - `tiptilt + dm`
      - `ShackHartmannWFS`
      - `Pyramid`
      - `BioEdge`
    - `steering + dm`
    - `focus + dm`
- AMDGPU workflow:
  - targets a self-hosted runner labeled `self-hosted`, `linux`, `amdgpu`
  - instantiates [`test/amdgpu`](../test/amdgpu)
  - runs the maintained AMDGPU hardware target
  - runtime equivalence includes the composite-optic low-order HIL surfaces:
    - `tiptilt + dm`
      - `ShackHartmannWFS`
      - `Pyramid`
      - `BioEdge`
    - `steering + dm`
    - `focus + dm`

The CPU and AMDGPU workflows are the active validation paths. The CUDA workflow
remains useful restoration scaffolding, but archived CUDA evidence does not
validate the current checkout while its runner is unavailable.

## Validation Rule

Claiming backend support should imply:

- functional correctness in the main testsets
- optional backend smoke coverage in `Pkg.test()`
- maintained benchmark or profile evidence where runtime claims are made

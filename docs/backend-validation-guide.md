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

Bare `Pkg.test()` runs every registered suite in
[`runtests.jl`](../test/runtests.jl) through the grouped testsets under
[`test/testsets`](../test/testsets):

- quality, API, and deterministic RNG policy
- KernelAbstractions CPU parity and tomography
- core optics, direct science, and atmosphere
- control and runtime
- detectors and WFS
- plant topology, preparation, product providers, RNG ownership, and
  calibration illumination
- calibration and analysis
- reference, tutorial, and Gate 0 regression

These are normal correctness tests, not backend throughput checks.

For a development edit/test loop, pass one or more stable suite or group names
through Julia's `test_args` interface. For example:

```sh
julia --project=. --startup-file=no -e \
  'using Pkg; Pkg.test(test_args=["plant-time"])'

julia --project=. --startup-file=no -e \
  'using Pkg; Pkg.test(test_args=["plant"])'
```

The first command runs one fine-grained suite; the second runs the broader
plant group. Multiple selectors form a de-duplicated union in canonical suite
order. List the current suites and groups with:

```sh
julia --project=. --startup-file=no -e \
  'using Pkg; Pkg.test(test_args=["--list"])'
```

An unknown selector fails rather than silently running no tests. Selective
runs are development evidence only: bare `Pkg.test()` remains the complete CPU
composition and release gate, and CI continues to run it.

### Optional backend smoke in `Pkg.test()`

These run after the functional testsets and skip cleanly if the backend package
or device runtime is unavailable:

- [`optional_amdgpu_backends.jl`](../test/optional_amdgpu_backends.jl)
- [`optional_cuda_backends.jl`](../test/optional_cuda_backends.jl)

Shared smoke scaffolding lives in:

- [`backend_optional_common.jl`](../test/backend_optional_common.jl)

The reduced maintained smoke covers:

- explicit atmosphere epochs and prepared, device-resident finite/infinite
  direction rendering into caller-owned pupil products
- atmospheric field propagation
- prepared physical Shack-Hartmann optical formation, detector acquisition,
  and centroid estimation on device-resident arrays
- same-grid legacy spectral diffractive SH plus prepared native-grid bundle
  retention and fail-closed single-output rejection for distinct wavelengths
- same-grid spectral SH detector acquisition with non-unit sampled QE and
  exposure scaling
- deterministic diffractive SH detector/export equivalence against CPU
- deterministic composite-optic low-order runtime parity against CPU:
  - `tiptilt + dm`
    - `ShackHartmannWFS`
    - `Pyramid`
    - `BioEdge`
  - `steering + dm`
  - `focus + dm`
- curvature-through-atmosphere
- prepared LiFT photon-rate formation, analytic interaction matrices,
  rate/count/normalized observations, analytic and numerical reconstruction,
  and dense/separable object convolution
- prepared native uniform calibration illumination through a typed
  detector-input path, stable materialization RNG ownership, and ordinary
  detector acquisition without host scalar indexing
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
accelerator support boundary. The CUDA target is also current and manually
hardware-validated on the WSL RTX host, but is not yet a release gate because
a continuously available CUDA CI runner has not been established. The backend
projects under
[`test/amdgpu`](../test/amdgpu) and [`test/cuda`](../test/cuda) declare the root
package as a path source, keeping accelerator packages out of normal
`Pkg.test()` while still resolving the full checkout dependency graph in a
clean environment.

### Apple Silicon BLAS/LAPACK selection

AppleAccelerate.jl is validated as an explicit application-level linear-algebra
provider recognized through an AdaptiveOpticsSim weak-dependency extension, not
as a core dependency or array backend. The isolated
[`test/appleaccelerate`](../test/appleaccelerate) project resolves the current
maintained AppleAccelerate 0.7 line on a `macos-15` Apple Silicon runner. Two
fresh processes establish distinct facts:

- loading AdaptiveOpticsSim normally does not load AppleAccelerate and leaves
  provider choice to the application
- loading AppleAccelerate explicitly routes representative ILP64 BLAS and
  LAPACK symbols through Accelerate, selects its single-threaded mode, selects
  vDSP plans for supported package FFTs, preserves FFTW fallback for unsupported
  shapes, and passes the full CPU suite

AppleAccelerate's AbstractFFTs extension supports non-empty, power-of-two 1D and
2D complex transforms, but its generic in-place plan adapter allocates temporary
split-complex arrays. AdaptiveOpticsSim's optional extension instead prepares a
reusable vDSP setup with package-owned work buffers, retaining allocation-free
repeated optical propagation. It selects FFTW for partial-dimension,
arbitrary-size, and three-or-more-dimensional CPU transforms. Loading the
package therefore improves the supported Apple Silicon path without narrowing
the existing CPU FFT shape boundary or weakening hot-path allocation contracts.

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

On the normal point-source path, prepared Shack-Hartmann field formation, FFTs,
spot sampling, detector acquisition, and centroid estimation remain batched on
the device. The explicit optical front end keeps distinct wavelengths as
device-resident `OpticalProductBundle` leaves on their native angular grids;
it does not index-add or implicitly resample them. The legacy single-product
`measure!` convenience path remains restricted to a common wavelength grid.

Centroid cutoff/statistics use the reusable host centroid workspace and copy the
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
julia --project=benchmarks benchmarks/benchmark_control_operators.jl
julia --project=benchmarks benchmarks/benchmark_loop_order_simd.jl
julia --project=benchmarks benchmarks/benchmark_detector_hil_latency.jl
julia --project=benchmarks benchmarks/benchmark_gate0_latency.jl
julia --project=benchmarks benchmarks/benchmark_pre_hil_backend_latency.jl cpu local_cpu

julia --project=benchmarks/amdgpu -e 'using Pkg; Pkg.instantiate()'
julia --project=benchmarks/amdgpu benchmarks/benchmark_amdgpu.jl
julia --project=benchmarks/amdgpu benchmarks/benchmark_pre_hil_backend_latency.jl amdgpu local_amdgpu

julia --project=benchmarks/cuda -e 'using Pkg; Pkg.instantiate()'
julia --project=benchmarks/cuda benchmarks/benchmark_pre_hil_backend_latency.jl cuda wsl_cuda
```

`benchmark_gate0_latency.jl` accepts `AOS_GATE0_CARD_IDS` as a comma-separated,
predeclared subset. The artifact records the explicit selection mode and the
ordered card IDs. Use this for incremental evidence owned by one PR so
unrelated predecessor-card host jitter does not invalidate the new cards; do
not omit an affected card after observing its result.

`benchmark_detector_hil_latency.jl` is the conventional-detector latency card
suite. It covers CMOS, CMOS with explicit presampling response and IPC, CCD,
fast linear EMCCD, HgCdTe avalanche CDS, and 16-sample Skipper CCD capture. The measured boundary
starts with an input detector frame already available in memory and ends when
the converted output frame is ready. It therefore excludes external RTC
transport, frame-grabber I/O, and camera-link scheduling.

The default contract uses one Julia, BLAS, and FFT thread; three warmed serial
closed-loop repetitions; 100,000 samples per card; and a fixed-size
`HdrHistogram.Histogram`. Since the next capture begins only after the previous
one completes, there is no independent arrival schedule and coordinated-
omission correction is intentionally not applied. First capture and steady-
state allocation measurements are reported separately. The histogram source
is an untagged commit from the GitHub `HdrHistogram.jl` repository, pinned in
the benchmark project rather than resolved from a local checkout.

For a quick harness check, reduce the sample count, repetitions, and frame size:

```bash
AOS_DETECTOR_HIL_SAMPLES=1000 \
AOS_DETECTOR_HIL_RUNS=1 \
AOS_DETECTOR_HIL_SIZE=32 \
julia --project=benchmarks benchmarks/benchmark_detector_hil_latency.jl
```

Runs below 100,000 samples label p99.9 as diagnostic. Set
`AOS_DETECTOR_HIL_OUTPUT` to retain a TOML artifact. Allocation gates are always
enforced. Machine-specific latency gates are evaluated only when
`AOS_DETECTOR_HIL_BASELINE` names a compatible prior artifact; the default p99
limit is 1.25 times its median p99 and can be changed with
`AOS_DETECTOR_HIL_REGRESSION_FACTOR`. The script rejects baselines whose sensor,
frame size, response, coupling, sampling, noise, or output type differs.

The current local-host baseline is archived in
[`2026-07-14-detector-hil-latency.toml`](../benchmarks/results/detectors/2026-07-14-detector-hil-latency.toml).
It used 64×64 frames, three 100,000-sample runs, and an AMD Ryzen 7 6800H. All
cards were allocation-free after warmup. Median-over-run latency summaries were:

| Card | p50 | p99 | p99.9 |
| --- | ---: | ---: | ---: |
| CMOS | 56.5 μs | 71.4 μs | 103.8 μs |
| CMOS with response and IPC | 92.0 μs | 110.5 μs | 167.9 μs |
| CCD | 79.6 μs | 96.3 μs | 144.1 μs |
| fast linear EMCCD | 85.2 μs | 102.9 μs | 153.1 μs |
| HgCdTe avalanche CDS | 188.9 μs | 223.5 μs | 326.1 μs |
| 16-sample Skipper CCD | 202.4 μs | 233.2 μs | 338.4 μs |

These values are a regression baseline for that host and contract, not an
external-RTC latency SLO or a prediction for other frame sizes and CPUs.

For the detector-output HIL path, use the same workload and sample count on
all backends:

```bash
julia --project=. scripts/profile_revolt_hil_runtime.jl cpu benchmarks/assets/revolt_like cmos default none 100 10
julia --project=benchmarks/cuda scripts/profile_revolt_hil_runtime.jl cuda benchmarks/assets/revolt_like cmos default none 100 10
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
The current Julia 1.12.6 cross-host rerun, including WSL CPU/CUDA, local
CPU/AMDGPU, a preallocated CUDA host-readout boundary, physical command parity,
and all individual repetitions, is archived in
[`2026-07-14-wsl-cuda-local-amdgpu.toml`](../benchmarks/results/platform/2026-07-14-wsl-cuda-local-amdgpu.toml).

The final pre-HIL backend contract uses the same physical `277`-command,
`352x352` REVOLT-like workload, but archives three 500-sample HdrHistogram
runs, first use, steady-state allocations, CPU parity, maintained-array residency,
and exact project/manifest hashes. GPU artifacts distinguish a synchronized
backend-ready result from a preallocated host-ready copy and the transfer-only
boundary. This remains serial service-time evidence: there is no independent
arrival schedule, queue, overload test, or external RTC transport.

Current clean-revision artifacts are:

- [local CPU](../benchmarks/results/platform/2026-07-18-pre-hil-11-local-cpu.toml)
- [WSL CPU](../benchmarks/results/platform/2026-07-18-pre-hil-11-wsl-cpu.toml)
- [WSL CUDA](../benchmarks/results/platform/2026-07-18-pre-hil-11-wsl-cuda.toml)

| Placement and boundary | Median p50 | Median p95 | Steady allocation |
| --- | ---: | ---: | ---: |
| local CPU, backend-ready | 4.108 ms | 4.624 ms | 0 B |
| WSL CPU, backend-ready | 4.923 ms | 5.534 ms | 0 B |
| WSL CUDA, backend-ready | 1.264 ms | 1.689 ms | 43,376 B |
| WSL CUDA, host-ready | 1.432 ms | 1.833 ms | 42,816 B |
| WSL CUDA, transfer-only | 0.112 ms | 0.143 ms | 96 B |

All listed correctness, residency, allocation, absolute-p95, and relative-p95
gates pass. The WSL target used CUDA.jl 6.2.1, KernelAbstractions.jl 0.9.42,
and Julia 1.12.6. The current maintained hardware targets passed `412/412` CUDA
checks and `422/422` AMDGPU checks with scalar indexing disabled, including the
shared LiFT matrix and device-resident schedule-free atmosphere
materialization, direct-science formation, and detector fan-out. The
current local AMDGPU latency artifact remains the
[July 14 cross-host characterization](../benchmarks/results/platform/2026-07-14-wsl-cuda-local-amdgpu.toml):
the July 18 timing repetitions completed, but no replacement raw-histogram
artifact was retained after the host's Julia 1.12.6 installation failed a
package-independent GC check. This is not promoted into a new latency baseline.

The final composed CPU Gate 0 run is preserved even though `G0-PERF-05` missed
its relative p99 limit by 96 ns while every absolute and allocation gate
passed. An immediate same-contract, same-revision, three-run Shack-Hartmann
confirmation passed at 108.543 μs median p99 against the 130.399 μs limit:

- [full Gate 0 catalog](../benchmarks/results/gate0/2026-07-18-pre-hil-11-full-backend-evidence.toml)
- [Shack-Hartmann confirmation](../benchmarks/results/gate0/2026-07-18-pre-hil-11-shack-hartmann-confirmation.toml)

`benchmark_gate0_latency.jl` accepts a comma-separated
`AOS_GATE0_BASELINES` catalog. The artifact records every source path, SHA-256,
source revision, and supplied card ID so staged baselines compose without
silently replacing predecessor evidence.

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

For dense versus truncated factorized reconstruction, use
`benchmark_control_operators.jl`. The current synthetic local-host result is
archived in
[`2026-07-13-control-operators.toml`](../benchmarks/results/platform/2026-07-13-control-operators.toml).
Treat its rank as a performance workload only; production rank must come from
optical and control validation.

For CPU loop-order and vectorization work, use
`benchmark_loop_order_simd.jl`. The maintained comparison archives legacy and
current implementations for LiFT convolution and LGS elongation, checks
bitwise output equality, and reports warmed allocations and tails. The current
local-host artifact is
[`2026-07-13-loop-order-native-simd.toml`](../benchmarks/results/platform/2026-07-13-loop-order-native-simd.toml).
Julia's built-in `@simd` was sufficient for the independent contiguous loops;
`SIMD.jl` is not a dependency. Reconsider explicit vector types only for a
profiled kernel that native code generation fails to vectorize and only with
separate CPU-feature and backend-portability validation.

The CUDA hardware target is current again on the WSL RTX 3050 Ti host. Use the
isolated `test/cuda` environment for maintained hardware validation and the
equivalent `benchmarks/cuda` environment for benchmark-only dependencies.
Treat comparisons with the native AMD host as whole-system comparisons because
the CPU, OS, GPU, and power-management paths all differ.

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
  - runs the normal `Pkg.test()` suite on Linux, Apple Silicon macOS, and
    Windows hosted runners
  - runs a separate Apple Silicon job that proves backend-neutral normal load,
    then explicitly selects AppleAccelerate BLAS/LAPACK and reruns the full CPU
    suite with supported vDSP FFT plans and FFTW fallback plans
  - runs the isolated AcceleratedKernels/Dagger scheduler extension tests on a
    four-thread Linux job
- CUDA workflow:
  - targets a self-hosted runner labeled `self-hosted`, `linux`, `cuda`
  - instantiates [`test/cuda`](../test/cuda)
  - runs the maintained CUDA hardware target whenever a matching self-hosted
    runner is online; the current manual WSL validation host is not a
    continuously available release runner
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

The CPU and AMDGPU workflows are the continuously available validation paths.
The CUDA workflow and manual WSL target exercise the same fail-fast hardware
surface, but CUDA remains outside the release gate until runner availability is
continuous. Evidence applies only to the exact checkout that was tested.

## Validation Rule

Claiming backend support should imply:

- functional correctness in the main testsets
- optional backend smoke coverage in `Pkg.test()`
- maintained benchmark or profile evidence where runtime claims are made

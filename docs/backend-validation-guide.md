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
- control primitives, reconstruction, and explicit model compositions
- detector, WFS-common, Shack–Hartmann, Pyramid/Bi-O-edge,
  Zernike/Curvature, and LiFT suites
- plant topology, canonical time, deterministic scheduling, trigger
  distribution, detector transitions, preparation, product providers, RNG
  ownership, and calibration illumination
- calibration workflows, NCPA, optical analysis, and cross-domain interface
  conformance
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
order. Shared test-only protocol fixtures are declared as suite dependencies,
loaded once, and do not rely on a preceding suite's include side effects. List
the current suites and groups with:

```sh
julia --project=. --startup-file=no -e \
  'using Pkg; Pkg.test(test_args=["--list"])'
```

An unknown selector fails rather than silently running no tests. Selective
runs are development evidence only: bare `Pkg.test()` remains the complete CPU
composition and release gate.

The maintained CI partition consists of `ci-foundations`,
`ci-sensors-control`, `ci-plant-runtime`, and `ci-plant-optics`. Registry
validation requires their suite membership to be disjoint and complete, so a
new registered suite cannot silently escape platform or coverage closure.
These selectors are intended for CI and bounded local diagnosis; owner-focused
development should normally use the narrower suite or domain-group names.

### Deterministic event-scheduler evidence

Run the focused Gate 3 scheduler correctness and allocation suite with:

```sh
julia --project=. --startup-file=no -e \
  'using Pkg; Pkg.test(test_args=["plant-scheduler"])'
```

The maintained generator-count characterization uses the benchmark-only
HdrHistogram dependency and one Julia/BLAS/FFT thread:

```sh
julia --threads=1 --project=benchmarks --startup-file=no \
  benchmarks/benchmark_gate3_event_scheduler.jl
```

Its versioned contract is
[`gate3_event_scheduler.toml`](../benchmarks/contracts/gate3_event_scheduler.toml),
and the current raw-histogram artifact is
[`2026-07-21-event-scheduler-gate3-closure.toml`](../benchmarks/results/gate3/2026-07-21-event-scheduler-gate3-closure.toml).
The timed boundary is one serial linear scan, claim, checked integer next-time
calculation, and reschedule across 1, 8, 32, 128, and 256 active generators.
This is warmed, self-paced in-process scheduler service-cost evidence. It does
not include detector or optical work, wall-clock pacing, queues, overload,
transport, or an external RTC, and therefore establishes neither fixed-arrival
latency nor production instrument capacity.

The focused deterministic trigger-distribution suite is:

```sh
julia --project=. --startup-file=no -e \
  'using Pkg; Pkg.test(test_args=["plant-triggers"])'
```

It covers exact finite fault traces, canonical nested fan-out, source/link fault
scope, nominal/delivered/reported timestamp separation, chronology, conservative
capacity rejection, 100,000-cycle constant storage, inference, and zero warmed
allocation with caller-owned delivery slots. It is a serial CPU semantic gate,
not a detector-transition, stochastic-jitter, external-RTC, transport, or GPU
performance claim.

The focused global-shutter detector-transition suite is:

```sh
julia --project=. --startup-file=no -e \
  'using Pkg; Pkg.test(test_args=["plant-detector-transitions"])'
```

It covers exact lifecycle boundaries, transition-error atomicity, whole-frame
detector parity and MTF preservation, time-resolved up-the-ramp reads, and
warmed inference/allocation checks. Accelerator residency is exercised in the
optional hardware matrix rather than this CPU suite.

The focused composed multi-rate event suite is:

```sh
julia --project=. --startup-file=no -e \
  'using Pkg; Pkg.test(test_args=["plant-event-composition"])'
```

It covers exact rolling-exposure/global-reset row bands, frame-transfer
image/storage overlap, independent science/NGS/LGS sample periods, periodic
and delivered-trigger acquisition starts, same-time causal ordering,
declaration-order replay, fixed storage through a long horizon, zero-allocation
warmed detector kernels, and the bounded heterogeneous-orchestration allocation
budget. It is a serial virtual-time correctness suite, not the composed
topology benchmark, wall-clock HIL latency, fixed-arrival capacity, or an
integrated multi-device claim.

The maintained composed-plant characterization uses the same single-threaded
execution policy:

```sh
julia --threads=1 --project=benchmarks --startup-file=no \
  benchmarks/benchmark_gate3_multi_rate_plant.jl
```

Its [versioned contract](../benchmarks/contracts/gate3_multi_rate_plant.toml)
and [clean artifact](../benchmarks/results/gate3/2026-07-21-multi-rate-plant.toml)
cover one direct-science path, an off-axis NGS Shack–Hartmann path, a
finite-height LGS pyramid path, five independent conventional detector
acquisitions, common-trigger fan-out and faults, declaration reordering, fixed
storage, direct long-period jumps, and a bounded allocation window. The timed
boundary is one heterogeneous canonical plant timestamp. The three 10,000-
sample runs are self-paced service-cost distributions; they are not fixed-rate
HIL latency or capacity measurements.

Gate 3 artifacts store each HdrHistogram as a lossless sparse stream of
big-endian `Int64` value/count pairs in base64. The histogram range and
significant figures remain in the artifact, and
[`hdr_histogram_artifact.jl`](../benchmarks/support/hdr_histogram_artifact.jl)
round-trips the stream before publication. Encoding stays outside the timed
boundary and avoids depending on a particular compressed-wire encoder.

### Gate 4 command and optic evidence

Run the focused Gate 4 correctness, ownership, routing, reduced-order, and
autonomous-optic suites with:

```sh
julia --project=. --startup-file=no -e \
  'using Pkg; Pkg.test(test_args=["quality", "gate4"])'
```

The maintained command-responsive plant characterization is:

```sh
julia --threads=1 --project=benchmarks --startup-file=no \
  benchmarks/benchmark_gate4_command_plant.jl
```

Its [versioned contract](../benchmarks/contracts/gate4_command_plant.toml) and
[clean raw-histogram artifact](../benchmarks/results/gate4/2026-07-24-command-plant.toml)
cover two independently timed co-pupil scalar optics, atomic two-optic
application, command-during-exposure causality, exact terminal accounting,
bounded endpoint storage, declaration-order replay, and a 10,000-cycle storage
check. A standalone warmed endpoint cycle allocates zero Julia heap bytes. The
larger composed presentation/application/optical/detector/accounting boundary
observed 4,720.736 bytes per cycle against its declared 8 KiB ceiling.

Each of three runs records 10,000 service-time samples as lossless sparse
HdrHistograms, which supports the reported p99. The load is serial and
self-paced: it deliberately excludes fixed-rate arrivals, wall-clock pacing,
queues, overload, external RTC transport, and response latency.

### Gate 5 conjugated optical-placement evidence

Run the focused placement, geometry, MCAO/MOAO, sampled-aberration, and
integrated closure suites with:

```sh
julia --threads=1 --project=. --startup-file=no -e \
  'using Pkg; Pkg.test(test_args=["gate5"])'
```

The maintained integrated characterization is:

```sh
julia --threads=1 --project=benchmarks --startup-file=no \
  benchmarks/benchmark_gate5_optical_placement.jl
```

Its [versioned contract](../benchmarks/contracts/gate5_optical_placement.toml)
and [clean raw-histogram artifact](../benchmarks/results/gate5/2026-07-25-optical-placement.toml)
compose eight simultaneously due Fraunhofer paths: four WFS directions
(alternating finite-height LGS and NGS), four science directions, common native
DMs at the pupil and two atmospheric conjugates, one selected-path MOAO DM per
science direction, one common replacement OPD, and one science-only NCPA per
science direction. The exact oracle also checks declaration-order replay and a
3×3-on-5×5 finite-support replacement.

Prepared storage and product hashes remain fixed through 5,000 post-warmup
cycles. The warmed eight-path boundary observed 5,632 Julia heap bytes per
all-path cycle against its declared 8 KiB ceiling. Three 10,000-sample
self-paced runs recorded 107.711 microseconds median p50, 124.287 microseconds
worst p99, and 9.229 kcycles/s median throughput. The 4/8/16-path preparation
matrix preserves exact binding counts and records preparation allocation,
elapsed time, and prepared-summary size; it deliberately makes no
production-capacity or bounded-code-generation claim.

This is serial CPU optical-cycle service-cost and numerical evidence. It does
not qualify fixed arrivals, wall-clock pacing, external RTC transport,
integrated GPU event-loop execution, parallel workers, detailed external relay
or coronagraph prescriptions, or NFIRAOS/MORFEO instrument capacity. Direct
device-resident DM and sampled-surface parity remain covered by the maintained
CUDA and AMDGPU hardware targets.

### Optional backend smoke in `Pkg.test()`

These run after the functional testsets and skip cleanly if the backend package
or device runtime is unavailable:

- [`optional_amdgpu_backends.jl`](../test/optional_amdgpu_backends.jl)
- [`optional_cuda_backends.jl`](../test/optional_cuda_backends.jl)

Shared smoke scaffolding lives in:

- [`backend_optional_common.jl`](../test/backend_optional_common.jl)

Use
[`profile_direct_imaging_batch_submissions.jl`](../scripts/profile_direct_imaging_batch_submissions.jl)
under a vendor profiler to compare `batch` with `independent` WFS optics execution for the
same ordered physical samples. Its positional arguments are backend, execution
mode, sample count, warmup count, and measured count. The script reports the
prepared FFT-execution count as an operation-count proxy; it is not a latency,
throughput, or HIL deadline benchmark.

The reduced maintained smoke covers:

- explicit atmosphere epochs and prepared, device-resident finite/infinite
  direction rendering into caller-owned pupil products
- atmospheric field propagation
- prepared native-Fraunhofer direct-science batching for ordered off-axis and
  spectral physical sources, with one device-resident field/output stack,
  reusable optical-axis FFT plan, per-sample product views, independent-path
  parity, detector fan-out, and explicit completed return
- Plant-integrated ownership of compatible equal-rate direct-science paths on
  one concrete accelerator, including one retained device execution context,
  shared atmosphere/direct-imaging submission, explicit same-device handoff to
  path-local products, conventional detector lifecycle parity, and completed
  return before path publication
- Plant-integrated ownership of compatible equal-rate Shack-Hartmann, Pyramid,
  or Bi-O-edge paths on one concrete accelerator, with one retained device
  context and atmosphere-direction batch while each exact family keeps its
  existing prepared lenslet/modulation/spectral pipeline
- exact fallback for singleton, unequal-rate or unequal-origin, mixed-family,
  incompatible-signature/product, Zernike, and Curvature path groups
- a six-row conventional detector matrix spanning CCD single-read,
  frame-transfer EMCCD, global- and rolling-shutter CMOS, and global-shutter
  HgCdTe single-read and scheduled up-the-ramp lifecycles; device residency,
  exact transition/readout sequence, response/MTF metadata, independent
  response-before-exposure numerics, and bounded warmed host allocation are
  checked directly
- prepared physical Shack-Hartmann optics, detector acquisition,
  and centroid estimation on device-resident arrays
- same-grid legacy spectral diffractive SH plus prepared native-grid bundle
  retention and fail-closed single-output rejection for distinct wavelengths
- same-grid spectral SH detector acquisition with non-unit sampled QE and
  exposure scaling
- deterministic diffractive SH detector/export equivalence against CPU
- independent DM and modal/low-order optic formation and application parity
- prepared controller-output routing with a device-resident view and explicit
  host/backend mismatch rejection
- prepared same-grid and affine finite-support sampled pupil-surface coupling,
  with path-local storage remaining device resident
- dynamic cycle-averaged circular-Pyramid radius/enable updates, disabled-state
  centering, normalized quadrature weights, and device residency
- curvature-through-atmosphere
- prepared LiFT photon-rate formation, analytic interaction matrices,
  rate/count/normalized observations, analytic and numerical reconstruction,
  and dense/separable object convolution
- prepared native uniform calibration illumination through a typed
  detector-input path, stable materialization RNG ownership, and ordinary
  detector acquisition without host scalar indexing
- exact global-shutter detector-event accumulation plus windowed scheduled
  up-the-ramp snapshots/fitting on evolving charge, with detector products
  remaining device resident and scalar indexing disabled
- direct rolling-shutter row-band integration and frame-transfer image/storage
  lifecycle checks with detector products and the lifecycle-state storage
  frame remaining device resident and scalar indexing disabled
- MKID accumulated-count capture, source passband handling, output conversion,
  and device residency

For independent controllable optics, optional backend smoke checks command
isolation and additive application without a packed aggregate. Controller
routing separately verifies that a borrowed GPU view retains exact endpoint
backend/device compatibility and that host storage fails preparation.

For broader backend-audit coverage, use:

- `ADAPTIVEOPTICS_TEST_FULL_AMDGPU=1`
- `ADAPTIVEOPTICS_TEST_FULL_CUDA=1`

which route through [`gpu_smoke_contract.jl`](../scripts/gpu_smoke_contract.jl).

That full matrix is useful for subsystem audit and backend triage, but it does
not by itself define the supported GPU surface. The supported GPU surface is
defined by the dedicated hardware-backed validation targets below.

For explicit hardware-backed backend validation targets that combine the
optional backend matrix with the GPU builder and maintained REVOLT-like
production-shaped smoke, run:

- `julia --project=test/amdgpu --startup-file=no test/runtests_amdgpu.jl`
- `julia --project=test/cuda --startup-file=no test/runtests_cuda.jl`

For fail-fast detector qualification without compiling unrelated optics, WFS,
control, or orchestration kernels, run:

- `julia --project=test/amdgpu --startup-file=no test/runtests_amdgpu_detectors.jl`
- `julia --project=test/cuda --startup-file=no test/runtests_cuda_detectors.jl`

These focused targets cover the shared response/MTF and acquisition pipeline,
CCD, EMCCD, CMOS, conventional and linear-avalanche HgCdTe, Skipper CCD,
InGaAs, linear APD, SPAD, and MKID surfaces. They complement rather than replace
the complete integration targets.

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
clean environment. Their maintained workflows preserve the ordinary
non-coverage hardware target, including allocation assertions, then run a
bounded extension-owner coverage probe. They do not duplicate the full hardware
matrix.

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

The same Apple Silicon job separately instantiates [`test/metal`](../test/metal)
and loads Metal.jl with AdaptiveOpticsSim. That bounded smoke proves that the
Metal weak-dependency extension loads and that its backend discovery generics
extend `Backends`. It is an extension ownership/load check, not a claim of the
full CUDA/AMDGPU numerical hardware matrix on Metal.

The full GPU smoke matrix now also pins the exact batched Shack-Hartmann
detector-processed export surface that previously regressed on CUDA:

- null-noise diffractive SH with detector capture
- CPU vs GPU comparison of the internal legacy Shack-Hartmann spot diagnostic in
  [`shack_hartmann.jl`](../src/wfs/shack_hartmann.jl)

This keeps the detector-processed lenslet-spot storage under backend parity
coverage, not just the slope output. It does not misidentify the legacy
detector object's nominal frame storage as the Shack-Hartmann acquisition
product.

### AMDGPU host-mirror boundaries

With AMDGPU 2.7, photon Poisson sampling, detector readout correction, windowed
multi-read products, and float-to-integer detector export use host-mirror
paths. Poisson sampling, correction, and integer export reuse detector host
buffers; windowed readout-product construction retains its existing allocation
behavior. Nearest-neighbor count redistribution and the rest of the maintained detector array math
remain backend kernels.

On the normal point-source path, prepared Shack-Hartmann optics, FFTs,
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
- [`profile_revolt_hil_runtime.jl`](../scripts/profile_revolt_hil_runtime.jl)
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
julia --project=benchmarks benchmarks/benchmark_revolt_graph_nodes.jl

julia --project=benchmarks/amdgpu -e 'using Pkg; Pkg.instantiate()'
julia --project=benchmarks/amdgpu benchmarks/benchmark_amdgpu.jl
julia --project=benchmarks/amdgpu benchmarks/benchmark_pre_hil_backend_latency.jl amdgpu local_amdgpu
AOS_REVOLT_GRAPH_BACKEND=amdgpu julia --project=benchmarks/amdgpu benchmarks/benchmark_revolt_graph_nodes.jl

julia --project=benchmarks/cuda -e 'using Pkg; Pkg.instantiate()'
julia --project=benchmarks/cuda benchmarks/benchmark_pre_hil_backend_latency.jl cuda wsl_cuda
AOS_REVOLT_GRAPH_BACKEND=cuda julia --project=benchmarks/cuda benchmarks/benchmark_revolt_graph_nodes.jl
```

`benchmark_revolt_graph_nodes.jl` profiles the executable Classic and Copper
external-RTC graphs, including atmosphere evolution. It reports each node with
an exact graph-context synchronization, the ordinary complete graph boundary,
and the separate frame-to-host, command-to-target, and full lockstep HIL
boundaries. The sum of node times has one synchronization per node and is
therefore diagnostic rather than an estimate of an asynchronously submitted
graph. `AOS_REVOLT_GRAPH_ARCHITECTURES`, `AOS_REVOLT_GRAPH_SAMPLES`, and
`AOS_REVOLT_GRAPH_WARMUP` select the graph subset and measurement length. The
default is a one-thread, closed-loop, self-paced service-cost measurement; it
does not model fixed arrivals, transport, overload, or wall-clock pacing.

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
self-paced repetitions; 100,000 samples per card; and a fixed-size
`HdrHistogram.Histogram`. Since the next capture begins only after the previous
one completes, there is no independent arrival schedule and coordinated-
omission correction is intentionally not applied. First capture and steady-
state allocation measurements are reported separately, and every run retains a
verified compact raw histogram. The histogram source
is an untagged commit from the GitHub `HdrHistogram.jl` repository, pinned in
the benchmark project rather than resolved from a local checkout.

For a quick harness check, reduce the sample count, repetitions, and frame size:

```bash
AOS_DETECTOR_HIL_SAMPLES=1000 \
AOS_DETECTOR_HIL_RUNS=1 \
AOS_DETECTOR_HIL_SIZE=32 \
AOS_DETECTOR_HIL_CARDS=DET-HIL-00 \
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

The independently maintained minimal prepared path is archived in
[`2026-08-01-shared-low-fidelity-service-cost.toml`](../benchmarks/results/detectors/2026-08-01-shared-low-fidelity-service-cost.toml).
It selects only `DET-HIL-00`: `NoiseNone`, null presampling response, scalar
QE, no optional detector effects, and no converted output buffer. This is
self-paced in-process service-cost evidence, not a fixed-arrival latency claim.
The current artifact retains three 100,000-sample raw histograms, passes an
independent radiometric oracle, proves input nonmutation and zero warmed heap
allocation, and passes a same-host p99 regression gate.

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
and Julia 1.12.6. The current Gate 4 candidate revision `191751f` passed
`421/421` CUDA checks and `431/431` AMDGPU checks with scalar indexing disabled,
including the shared LiFT matrix, device-resident schedule-free atmosphere
materialization, direct-science formation, detector fan-out, prepared
controller-output routing, and dynamic cycle-averaged circular-Pyramid
modulation. The
July 14
[cross-host characterization](../benchmarks/results/platform/2026-07-14-wsl-cuda-local-amdgpu.toml)
remains a historical REVOLT-like workload. It is not equivalent to the later
Gate 7 paired single-device owner contract and is not used as that contract's
relative baseline.

The Gate 5 closure validation head `02e5f29` passed `448/448` checks on the
local gfx1030 AMDGPU target and `438/438` checks on the WSL RTX 3050 Ti CUDA
target, both with Julia 1.12.6 and scalar indexing disabled. CUDA used CUDA.jl
6.2.1; both targets used KernelAbstractions.jl 0.9.42. In addition to the native
Plant DM active/staged and surface checks, the matrix keeps a prepared
sampled-aberration copy on the exact backend/device after caller mutation and
validates identity plus transformed replacement, including zero outside finite
support. This is direct prepared-storage and kernel evidence, not integrated
GPU event-loop, latency, or multi-device evidence.

The Gate 3 closure code revision `871f76a` passed the full CPU suite, the
`424/424` CUDA target on an RTX 3050 Ti under WSL, and the `434/434` AMDGPU
target on a local gfx1030 device, all with Julia 1.12.6. CUDA used CUDA.jl 6.2.1;
AMDGPU used AMDGPU.jl 2.7.0; both used KernelAbstractions.jl 0.9.42. The
integrated multi-rate plant remains a deterministic serial CPU oracle. The GPU
targets validate maintained device-resident optical surfaces and direct
detector lifecycles; they do not claim an integrated GPU scheduler, external
RTC latency, or fixed-arrival capacity.

The final detector qualification candidate `dd16596` is indexed by the
[detector closure catalog](../benchmarks/results/detectors/2026-08-01-detector-qualification-closure.toml).
Its bounded CPU qualification/composition set passed 3,192 checks with one
intentional broken marker. The dedicated local gfx1030 AMDGPU detector target
passed 244/244; InGaAs and SPAD stochastic Poisson checks remain explicit AMD
non-claims because GPUCompiler segfaults while compiling that kernel. The
broader AMDGPU target separately encountered a GPUCompiler segmentation fault
in the Pyramid geometric-slope kernel outside detector scope; follow-up is
tracked in [issue #200](https://github.com/DarrylGamroth/AdaptiveOpticsSim.jl/issues/200).
On WSL, the RTX
3050 Ti CUDA detector target passed 250/250 and the complete CUDA target passed
1,028/1,028 with CUDA.jl 6.2.1, KernelAbstractions.jl 0.9.42, Julia 1.12.6, and
scalar indexing disabled. CUDA remains validation-only rather than a production
support claim.

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
cores. Prepared plant deployments declare those counts together through
`Plant.CPUExecutionBudget` and validate the observed boundary without changing
process-global settings. Deterministic validation continues to use one thread.

The focused Gate 6 correctness surface runs the serial fallback on the ordinary
one-thread matrix and the fixed-owner grouped proof on four Julia threads:

```bash
julia --threads=1 --project=. -e \
  'using Pkg; Pkg.test(test_args=["gate6"])'
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  julia --threads=4 --project=. -e \
  'using Pkg; Pkg.test(test_args=["gate6"])'
```

The grouped case is an unpaced ownership, exact-replay, and warmed-allocation
test over test-only long-lived workers. It is not a fixed-arrival HIL latency or
production `Channel` endorsement.

For paired HdrHistogram service-cost characterization of the same boundary,
use the benchmark environment:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  julia --threads=4 --project=benchmarks --startup-file=no \
  benchmarks/benchmark_gate6_grouped_cpu.jl
```

[`gate6_grouped_cpu.toml`](../benchmarks/contracts/gate6_grouped_cpu.toml)
predeclares the four-context budget, exact serial/grouped replay oracle,
repeated-run and p99 sample counts, and warmed allocation ceilings. The
comparison is closed-loop service cost with no absolute or relative performance
gate; a reported ratio is not a grouped-CPU speedup or HIL-capacity claim.

The clean local
[`2026-07-25-grouped-cpu.toml`](../benchmarks/results/gate6/2026-07-25-grouped-cpu.toml)
artifact retains all six raw histograms. On the unpinned Ryzen 7 6800H
`powersave` run, grouped execution traded a higher median p50
(1.463 μs versus 1.312 μs) and lower aggregate throughput for a lower median
p99 (66.303 μs versus 90.175 μs). Both modes allocated 318.688 bytes per
logical timestamp in the maintained heterogeneous barrier, direct warmed group
calls stayed at or below 432 bytes, and no measured run collected garbage.
These mixed results are workload- and host-specific service-cost evidence, not
a general parallel speedup or external-RTC latency result.

For the independent Gate 6 whole-plant specialization-growth gate, run the
fresh-process topology matrix on one Julia, BLAS, and FFT-provider thread:

```bash
AOS_GATE6_TOPOLOGY_OUTPUT=benchmarks/results/gate6/topology-growth.toml \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  julia --threads=1 --project=benchmarks --startup-file=no \
  benchmarks/benchmark_gate6_topology_growth.jl
```

The benchmark refuses durable output from a dirty tree. Set
`AOS_GATE6_TOPOLOGY_QUICK=1` and omit the output path for a one-repetition
development check. The
[`gate6_topology_growth.toml`](../benchmarks/contracts/gate6_topology_growth.toml)
contract predeclares separate event and schedule-free-selection first-use
probes, fresh preparation limits, registry-type and code-growth proxies,
concrete-kernel inference, prepared storage, numerical parity, and a
per-path warmed allocation ceiling over synthetic 4/8/16-path shapes.

The clean
[`2026-07-25-topology-growth.toml`](../benchmarks/results/gate6/2026-07-25-topology-growth.toml)
artifact passed every gate. From 4 to 16 paths, median fresh preparation time
grew 1.052×, preparation allocation 1.029×, and prepared storage 2.744×.
Checked registry and HIL-entry types had one fingerprint, entry native-code
size had a 1.0 ratio, method-instance growth was zero, and representative path
and acquisition kernels had no meaningful `Any` slots. Warmed serial
orchestration used at most 704.04 Julia heap bytes per path per cycle.
The recorded HdrHistogram distributions are closed-loop service-cost
diagnostics, not fixed-arrival latency, external-RTC, production-instrument,
NFIRAOS, or MORFEO capacity evidence.

For independent simulation sweeps, run the ensemble scheduler benchmark before
selecting `ThreadedExecution`, `AcceleratedKernelsExecution`, or
`DaggerExecution`:

```bash
julia --threads=8 --project=benchmarks benchmarks/benchmark_ensemble_schedulers.jl
```

This is offline throughput evidence, not an external-RTC latency result.
Dagger is intended for task graphs, locality, and process/node scaling; AK is
intended for reusable local task partitioning on sufficiently large many-core
workloads. Keep direct serial Plant event execution as the HIL-neutral
baseline.
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

For focused CPU iteration on the Gate 7 device-model matrix, use:

```sh
julia --project=. --startup-file=no -e \
  'using Pkg; Pkg.test(test_args=["plant-device-model-matrix"])'
```

The matrix is also included by the `gate7` selector and both dedicated hardware
targets. At the Gate 7.5 revision, focused local AMDGPU runs completed
`212/212` composed WFS-to-detector checks and `155/155` standalone
conventional-detector checks. On the WSL RTX 3050 Ti with Julia 1.12.6 and
CUDA.jl 6.2.1, the corresponding focused runs completed the same `212/212` and
`155/155` checks. These counts are correctness/residency evidence for that
exact revision and hardware, not a latency or production-capacity benchmark.

Gate 7.6 adds one predeclared paired single-device service-cost contract:

```sh
julia --threads=1 --project=benchmarks --startup-file=no \
  benchmarks/benchmark_gate7_single_gpu.jl cpu local_cpu
julia --threads=1 --project=benchmarks/amdgpu --startup-file=no \
  benchmarks/benchmark_gate7_single_gpu.jl amdgpu local_amdgpu
julia --threads=1 --project=benchmarks/cuda --startup-file=no \
  benchmarks/benchmark_gate7_single_gpu.jl cuda wsl_cuda
```

Set `AOS_GATE7_OUTPUT` to a path outside the clean source tree to retain an
artifact. Durable runs reject sample-count overrides and dirty source. The
[contract](../benchmarks/contracts/gate7_single_gpu.toml) fixes a 128-sample
two-layer atmosphere, two compatible off-axis NGS diffractive
Shack-Hartmann paths, 50 warmups, and three 200-sample runs per boundary.
CUDA and AMDGPU disable scalar indexing and separate independent
device-ready, batched device-ready, batched host-ready, and transfer-only
completion.

Candidate revision `b80a419` produced the following median-run results:

| Placement and boundary | p50 | p95 | Throughput | Warmed Julia heap |
|---|---:|---:|---:|---:|
| Local CPU, independent device-ready | 1.152 ms | 1.194 ms | 864.3 Hz | 1,648 B |
| Local gfx1030 AMDGPU, independent device-ready | 0.489 ms | 0.544 ms | 1,972.9 Hz | 81,424 B |
| Local gfx1030 AMDGPU, batched device-ready | 0.527 ms | 0.571 ms | 1,834.1 Hz | 76,208 B |
| Local gfx1030 AMDGPU, batched host-ready | 0.507 ms | 0.568 ms | 1,906.1 Hz | 74,576 B |
| Local gfx1030 AMDGPU, transfer-only | 0.028 ms | 0.071 ms | 23,875.2 Hz | 1,104 B |
| WSL RTX 3050 Ti CUDA, independent device-ready | 1.075 ms | 1.556 ms | 856.2 Hz | 66,160 B |
| WSL RTX 3050 Ti CUDA, batched device-ready | 1.144 ms | 1.418 ms | 832.4 Hz | 76,384 B |
| WSL RTX 3050 Ti CUDA, batched host-ready | 1.055 ms | 1.426 ms | 885.6 Hz | 77,712 B |
| WSL RTX 3050 Ti CUDA, transfer-only | 0.162 ms | 0.290 ms | 5,621.4 Hz | 944 B |

All predeclared absolute, allocation, parity, residency, synchronization, and
submission-proxy gates passed. Batched-to-independent median p95 was `1.050`
on AMDGPU and `0.911` on CUDA against the `1.5` ceiling. The owner reduced the
prepared atmosphere-render proxy from two calls to one while retaining two
WFS optics executions. The [artifact catalog](../benchmarks/results/gate7/manifest.toml)
retains raw HdrHistograms, run dispersion, first use, GC, exact project and
manifest hashes, package/runtime/driver/device identity, and explicit scope
exclusions. These are synchronized self-paced service-cost distributions, not
fixed-arrival HIL latency, overload, mixed-placement, multi-GPU, or instrument
capacity evidence.

This separation exists so:

- unit tests stay reasonably fast
- optional GPU availability does not block normal development
- benchmark evidence remains intentional and archived rather than implicit

## Automation

Checked-in CI automation now exists in:

- [../.github/workflows/cpu-validation.yml](../.github/workflows/cpu-validation.yml)
- [../.github/workflows/coverage.yml](../.github/workflows/coverage.yml)
- [../.github/workflows/cuda-backend-validation.yml](../.github/workflows/cuda-backend-validation.yml)
- [../.github/workflows/amdgpu-backend-validation.yml](../.github/workflows/amdgpu-backend-validation.yml)

Current intent:

- CPU workflow:
  - runs bare `Pkg.test()` as the authoritative full composition on Linux
  - runs the four bounded shards on Ubuntu for fast feedback; these duplicate
    the Linux composition deliberately but do not close the merge gate alone
  - runs the exact four-shard partition on hosted macOS and Windows, so every
    registered suite remains mandatory on all three hosted operating systems
  - runs a separate Apple Silicon job that proves backend-neutral normal load,
    then explicitly selects AppleAccelerate BLAS/LAPACK and reruns the full CPU
    suite with supported vDSP FFT plans and FFTW fallback plans
  - loads the Metal extension in a separate bounded test environment on that
    Apple Silicon job and verifies its canonical Backends owner surface
  - runs the focused Gate 6 fixed-owner CPU proof with four Julia threads and
    one BLAS/FFT thread per path-group owner
  - runs the isolated AcceleratedKernels/Dagger scheduler extension tests on a
    four-thread Linux job
- Coverage workflow:
  - runs the same exact four-shard partition under Julia coverage
    instrumentation
  - keeps allocation-byte assertions disabled through
    `ADAPTIVEOPTICS_TEST_COVERAGE` and Julia's coverage option
  - retains the shard runs as bounded instrumentation/allocation-policy gates
    without uploading their partial reports, because independent Julia coverage
    processes discover different compiler-attributed coverable lines
  - runs one baseline-comparable full-composition coverage process plus the
    focused KernelAbstractions CPU matrix for the authoritative project metric;
    Codecov combines that report with the Apple platform extension report,
    which contains the AppleAccelerate and Metal focused targets, and waits for
    both uploads before publishing project status or the PR comment
  - keeps project coverage over both `src/` and `ext/`; patch coverage for
    executable core `src/` is merge-blocking, while target-specific `ext/`
    patch coverage is reported separately as informational because native
    hardware validation is the authoritative extension gate
- CUDA workflow:
  - targets a self-hosted runner labeled `self-hosted`, `linux`, `cuda`
  - instantiates [`test/cuda`](../test/cuda)
  - runs the maintained CUDA hardware target whenever a matching self-hosted
    runner is online; the current manual WSL validation host is not a
    continuously available release runner
  - exercises the optional numerical/backend matrix, independent optic
    application, prepared Plant controller routing, GPU builder, and
    REVOLT-like production-shaped WFS smoke
  - uploads a bounded extension-owner coverage probe so CUDA-only methods
    contribute evidence when this workflow is dispatched, without rerunning the
    full matrix or disabling its allocation checks
- AMDGPU workflow:
  - targets a self-hosted runner labeled `self-hosted`, `linux`, `amdgpu`
  - instantiates [`test/amdgpu`](../test/amdgpu)
  - runs the maintained AMDGPU hardware target
  - exercises the same optional numerical/backend matrix, independent optic
    application, prepared Plant controller routing, GPU builder, and
    REVOLT-like production-shaped WFS smoke
  - uploads a bounded extension-owner coverage probe so AMDGPU-only methods
    contribute evidence without rerunning the full matrix or disabling its
    allocation checks

The CPU and AMDGPU workflows are the continuously available validation paths.
The CUDA workflow and manual WSL target exercise the same fail-fast hardware
surface, but CUDA remains outside the release gate until runner availability is
continuous. Evidence applies only to the exact checkout that was tested.

## Validation Rule

Claiming backend support should imply:

- functional correctness in the main testsets
- optional backend smoke coverage in `Pkg.test()`
- maintained benchmark or profile evidence where runtime claims are made

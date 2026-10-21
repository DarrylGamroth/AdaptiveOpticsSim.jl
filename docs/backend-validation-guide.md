# Backend Validation Guide

Status: active

## Purpose

Backend support is established per operation and exact target. A package load or
successful array allocation is not enough. Covered accelerator paths must run on
real hardware with scalar indexing disabled, preserve residency, and compare
against the declared CPU result.

## Test Selection

Tests are registered in `test/test_selection.jl`.

~~~bash
julia --project=. test/ci/run_selected_tests.jl algorithm-graphs
julia --project=. test/ci/run_selected_tests.jl core
julia --project=. test/ci/run_selected_tests.jl sensors
julia --project=. test/ci/run_selected_tests.jl backends
~~~

Bare `Pkg.test()` is the full local CPU composition gate:

~~~bash
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1   julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'
~~~

Use the impact planner to see the repository policy for changed paths:

~~~bash
julia --project=. test/ci/impact_planner.jl path/to/changed_file.jl
~~~

## CPU Matrix

`ka-cpu` exercises KernelAbstractions CPU kernels and style constraints.
`backend-smoke` exercises an optional backend only when its package was imported
before `AdaptiveOpticsSim` loaded and its device is functional. It does not
dynamically import a backend after AOS has loaded; that ordering is unsafe on
Julia 1.12. The dedicated AMDGPU and CUDA hardware entry points below preload
their backends. A skip in the CPU composition suite is not hardware
qualification.

Set Julia, BLAS, FFT-provider, OpenMP, and vendor math-library thread counts
deliberately. Deterministic evidence uses one thread. Performance evidence must
record the complete thread configuration and avoid nested oversubscription.

CUDA is the primary GPU performance-optimization target. AMDGPU is the
secondary portability and qualification target. Apply the same correctness,
device-residency, type-stability, and steady-state allocation gates to both;
focused tuning and counter work defaults to CUDA unless backend-specific
evidence identifies an AMDGPU defect or regression.

## AMDGPU

The maintained local ROCm target is:

~~~bash
julia --project=test/amdgpu --startup-file=no   -e 'using Pkg; Pkg.instantiate()'
julia --project=test/amdgpu --startup-file=no test/runtests_amdgpu.jl
~~~

For detector-only changes:

~~~bash
julia --project=test/amdgpu --startup-file=no test/runtests_amdgpu_detectors.jl
~~~

AMDGPU validation must record the Julia, AMDGPU, ROCm, kernel/driver, and device
versions. Run with scalar indexing disabled. Check output values, exact device
identity, and public completion behavior.

The hardware target also records the qualified regular-grid DM and pupil-OPD
composition sequence as one HIP Graph, changes the retained command buffer
between replays, and verifies the composed output. A separate captured
multilayer-atmosphere graph is compared frame by frame with stream execution;
the test verifies evolving turbulence, host epoch publication, and reset
reproducibility. Shack-Hartmann and both maintained Pyramid modulation paths
are compared with stream execution. The CCD, global-shutter CMOS, and
linear-mode EMCCD fixtures also verify stochastic evolution, stream parity,
ADC output, and reset replay. The target separately captures the
device-resident SplitMix64 normal and Poisson paths, verifies that replay
advances both streams, and verifies RNG reset reproducibility. This is required
evidence for `CapturedGraphExecution`; a standalone HIP Graph smoke test is not
sufficient.

The target also executes a two-lane `GroupedStreamGraphExecution` fixture. It
must retain distinct streams, honor the device-event barrier before both
dependent modal expansions, match the single-stream numerical result after
changed inputs, and remain reproducible after reset.

## Local CUDA

Run the CUDA target on the local CUDA host:

~~~bash
julia --project=test/cuda --startup-file=no -e "using Pkg; Pkg.instantiate()"
julia --project=test/cuda --startup-file=no test/runtests_cuda.jl
~~~

For detector-only changes, use `test/runtests_cuda_detectors.jl`. Record Julia,
CUDA.jl, toolkit, driver, and GPU versions. CUDA is manually validated but is
not a continuously available release gate.

The full CUDA target applies the same grouped two-branch parity, stream-identity,
changed-input, and reset checks as the AMDGPU target.

For the 352-square Shack-Hartmann pupil workgroup comparison:

~~~bash
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 AOS_SH_RESULTS=/tmp/sh-workgroups.toml \
  julia --startup-file=no --project=test/cuda \
  benchmarks/benchmark_sh_pupil_workgroups_cuda.jl cuda
~~~

The benchmark specializes only its process's pupil launcher. It compares CUDA
event batches of 64 pupil kernels and complete captured plant steps, with
alternating ABBA/BAAB groups and a duplicate default control. Input staging,
changed-command CPU oracles, bitwise pupil comparisons, compilation, and warmup
are outside timing. Timed plant input is fixed and concurrency is one. Repeat
in fresh processes and inspect the complete plant as well as kernel timing.
The [local comparison](../benchmarks/results/platform/2026-10-21-sh-pupil-workgroups.toml)
records source hashes, raw samples, launch geometry, and dependency constraints.

Set `AOS_SH_MATRIX=1` to compare `(44,8,1,1)` and the unchanged control across
pupil resolutions and lenslet counts. `AOS_SH_RESULTS` becomes a filename prefix;
each case appends `-N<resolution>-L<lenslets>.toml`. Repeat in fresh processes.
The maintained 64-square, 8-lenslet example can also be selected directly:

~~~bash
julia --startup-file=no --project=test/cuda -e \
  'include("benchmarks/benchmark_sh_pupil_workgroups_cuda.jl"); benchmark_sh_workgroups(; pupil_resolution=64, lenslets_per_axis=8, groups=_SH_MATRIX_GROUPS)' cuda
~~~

For pupil occupancy counters, profile one warmed node in a captured plant:

~~~bash
AOS_SH_OCCUPANCY_GROUP=default ncu --clock-control none --cache-control none \
  --profile-from-start off --graph-profiling node --kernel-name-base demangled \
  --kernel-name 'regex:gpu_sh_explicit_pupil_stack_kernel_' --launch-count 1 \
  --set basic --section Occupancy \
  julia --startup-file=no --project=test/cuda \
  benchmarks/profile_sh_pupil_occupancy_cuda.jl cuda
~~~

Repeat with `AOS_SH_OCCUPANCY_GROUP=44x8`. Achieved occupancy measures active
warps during this kernel; it does not measure whole-plant GPU utilization.
Keep profiler durations separate from unprofiled timing results. The
[size matrix and occupancy results](../benchmarks/results/platform/2026-10-21-sh-pupil-matrix.toml)
retain the default: the candidate regresses several geometries despite its
pad-44 kernel improvement.

For execution architecture, use `profile_aos_gpu_plant.jl cuda both` with Nsight
Systems, `--capture-range=cudaProfilerApi`, `--capture-range-end=repeat:2:defer`,
and `--cuda-graph-trace=node`. Compare host submissions, kernel activity,
transfers, and gaps between checked-stream and captured execution. Node tracing
perturbs timing; profiler start/stop and report generation can dominate the
driver's printed elapsed time. Use the captured activity intervals for timeline
analysis and separate unprofiled runs for latency claims. The size-matrix report
also records this architecture comparison and its raw traces.

For a broader CUDA kernel assessment, run
`benchmarks/profile_gpu_kernel_patterns_cuda.jl /tmp/patterns.toml` in the CUDA
test environment. It samples atmosphere interpolation, separable DM passes,
geometric slopes, detector stencils, and Pyramid intensity reduction with the
production workgroup policy. Timings include an output reset and one kernel;
CPU execution of the same kernels checks backend parity. Set
`AOS_GPU_PATTERN_PROFILE=1` for warmed counter collection without latency timing.

`benchmarks/profile_pyramid_efficiency_cuda.jl` compares complete stream steps
for both propagation strategies. `AOS_PYRAMID_COMPLETION=cooperative` or
`blocking` overrides the benchmark process's AOS stream waits while preserving
completion boundaries. Without an override it measures the package policy.
The [RTX 3080 assessment](../benchmarks/results/platform/2026-10-21-gpu-kernel-efficiency.toml)
records kernel counters, numerical checks, raw traces, and a benchmark-only
geometric-slopes reassociation experiment. It distinguishes measured cases from
unprofiled algorithms and does not claim general numerical qualification.

Synchronous CUDA kernel phases use a blocking wait on the current stream. This
preserves completion and device-error checking while avoiding cooperative-wait
bookkeeping. The calling Julia thread blocks during that wait; low-level CUDA
callers can explicitly select `blocking=false` on `synchronize_backend!` when
task cooperation is required. Prepared graph submission/completion and the
captured replay allocation contract remain separate from dynamic KA launches.

The production geometric-slopes kernel still averages every original adjacent
OPD difference. For same-precision Float32/Float64 subapertures of at least five
pixels on grids of at most 32 × 32 subapertures, one 32-lane group computes
each result. Denser grids retain the original kernel because one group per
output can reduce throughput. Nonfinite or overflow-sensitive
data use the original serial summation order. Smaller subapertures, mixed
representations, and edge-masked slopes retain their serial implementation.
The CPU implementation is unchanged; parallel additions may round differently.

Run `benchmarks/benchmark_geometric_slopes_workflow_cuda.jl cuda` with
`AOS_GEOMETRIC_RESULTS=/tmp/geometric.toml` to measure the complete captured
DM/pupil workflow followed by public geometric slopes. The timed endpoint is
completed GPU slope storage; numerical-oracle host copies are outside timing.
The same driver also compares original and selected kernels with paired CUDA
event batches. `test/cuda/geometric_slopes.jl` covers numerical edge cases and
fixed-storage replay; `test/cuda/kernel_completion.jl` checks both wait modes.
`benchmarks/benchmark_geometric_slopes_scaling_cuda.jl cuda` sweeps output-grid
density independently, including cases that retain the original kernel. The
[completion and slopes results](../benchmarks/results/platform/2026-10-21-cuda-completion-slopes.toml)
record the final selection, repeated workflow measurements, numerical checks,
and separate unrestricted-pilot and final density sweeps.

Qualify a KA upgrade separately from workgroup changes. Registered CUDA 6.4.2
constrains KA to 0.9; a development CUDA tree with the same version string can
have different compatibility and allocation behavior. Record exact source
trees and retain the numerical, zero-allocation, and fixed-storage capture
gates before adopting a new stack.

Do not treat SSH success, `CUDA.functional()`, or compilation alone as model
evidence. The test must execute the changed operation.

## Linux AArch64

Use the Raspberry Pi for changes that affect platform assumptions, packaging,
or generic CPU behavior:

~~~bash
ssh raspberrypi 'cd /path/to/AdaptiveOpticsSim.jl &&
  julia --project=. --startup-file=no -e "using Pkg; Pkg.test()"'
~~~

Record the exact checkout path, Julia version, architecture, and BLAS provider.
Do not run this target for unrelated documentation-only changes.

## Apple And Metal

AppleAccelerate and Metal use their dedicated environments under
`test/appleaccelerate` and `test/metal`. They are manual support questions, not
automatic consequences of CPU or AMDGPU qualification.

## Graph Validation

A GPU graph is admitted only when:

- every graph input and sparse parameter is native packed storage on the exact
  target
- every node supports that target
- node outputs and delayed storage can be allocated there
- execution stays inside the retained context
- no hidden host transfer or CPU fallback occurs
- one pending completion ticket prevents graph storage reuse

For native CUDA Graph or HIP Graph execution, additionally require:

- the adapter explicitly returns `GraphNodeCaptureSafe()`
- every recorded address remains fixed for the prepared run
- evolving replay state is device-resident
- the recorded device operation performs no dynamic execution-storage
  allocation, synchronization, result query, or host-side scientific-state
  mutation
- recorded and replayed work does not depend on a GPU host callback or
  AMDGPU HostCall
- Julia compilation, backend recording, and native graph instantiation remain
  cold preparation work; repeated replay uses the retained executable and
  bounded completion hooks
- replay with changed input-buffer contents produces the corresponding changed
  output
- every node in the graph explicitly satisfies the capture contract
- the ordered node sequence and delayed-link commits form one native executable

The lockstep HIL boundary intentionally owns host `Array` exchange buffers. A
completed GPU detector frame is copied to that host buffer only after successful
graph execution, and a complete validated command is copied back before the next
frame. For ordinary CUDA device arrays, host staging enqueues the copy on the
retained stream and completes it before publication. The boundary retains both
buffers through completion. Captured steady-state graph replay, complete HIL
frame/command cycles, and calibration exposures have separate zero-host-allocation
checks; this does not cover transport or an external controller. Other storage
kinds retain their ordinary backend copy behavior.

The modal OPD graph node uses the same numerical kernel as `combine_basis!`,
with completion owned by the prepared graph. Direct `combine_basis!` calls retain
their synchronous contract. CUDA regression tests cover signed cancellation,
pupil masking, nonfinite coefficients, changed commands, reset, and both
Float32 and Float64.

~~~bash
julia --startup-file=no --project=test/cuda test/cuda/hil_boundary.jl
~~~

Use `algorithm-graphs` for functional graph coverage and the generic backend
workload for measured frame service time:

~~~bash
julia --project=. test/ci/run_selected_tests.jl algorithm-graphs
julia --project=benchmarks --startup-file=no -e 'using Pkg; Pkg.instantiate()'
julia --project=benchmarks --startup-file=no benchmarks/benchmark_pre_hil_backend_latency.jl
~~~

Record the exact workload, target, resolution, noise settings, warmup, sample
count, and synchronization boundary. Captured execution requires an
accelerator backend and requires every node in the graph to qualify. Downstream
instrument packages own whole-graph profiling and must record submission,
completion-ticket wait, and target-ready service time separately.
AMDGPU stream execution retains the backend's ordinary HostCall-compatible
synchronization. Captured execution uses a blocking completion wait to avoid
host event/task allocation after replay; this is valid because capture
admission excludes HostCalls.
Do not compare submission against synchronized target-ready or host-ready
latency as though they were the same metric.

## Performance Evidence

Separate:

- compilation and first use
- warmed service time
- allocation
- device submission
- synchronization/host readiness
- explicit transfer
- end-to-end application latency

Use multiple repetitions and report distributions, not a single best sample.
A self-paced closed-loop benchmark measures service cost; it is not
fixed-arrival latency and does not establish overload behavior.

Historical artifacts remain attached to their recorded revision. If the
implementation they measured was removed, label them historical rather than
using them as current qualification.

For a focused target-ready profile of the retained Pyramid plant boundary,
run both modulation strategies with scalar indexing disabled:

~~~bash
ADAPTIVEOPTICS_PROFILE_BACKEND=cuda ADAPTIVEOPTICS_PROFILE_STEPS=1000 \
  julia --project=test/cuda --startup-file=no \
  scripts/profile_pyramid_plant_runtime.jl

ADAPTIVEOPTICS_PROFILE_BACKEND=amdgpu ADAPTIVEOPTICS_PROFILE_STEPS=1000 \
  julia --project=test/amdgpu --startup-file=no \
  scripts/profile_pyramid_plant_runtime.jl
~~~

Wrap the same command with Nsight Systems or ROCprofiler when collecting a
timeline. Use Nsight Compute for CUDA kernel-counter investigation after the
timeline identifies the kernel of interest. The script reports synchronized
target-ready time and warmed host allocation; it does not claim an arrival-rate
deadline or include an AOS/FGA device-resident bridge.

For the retained Bi-O-edge plant boundary, use the corresponding prepared
optics-plus-acquisition driver. Its CPU mode also emits a warmed Julia sampling
profile:

~~~bash
ADAPTIVEOPTICS_PROFILE_BACKEND=cpu ADAPTIVEOPTICS_PROFILE_STEPS=1000 \
  JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  julia --project=. --startup-file=no \
  scripts/profile_bi_o_edge_plant_runtime.jl

ADAPTIVEOPTICS_PROFILE_BACKEND=cuda ADAPTIVEOPTICS_PROFILE_STEPS=1000 \
  julia --project=test/cuda --startup-file=no \
  scripts/profile_bi_o_edge_plant_runtime.jl

ADAPTIVEOPTICS_PROFILE_BACKEND=amdgpu ADAPTIVEOPTICS_PROFILE_STEPS=1000 \
  julia --project=test/amdgpu --startup-file=no \
  scripts/profile_bi_o_edge_plant_runtime.jl
~~~

The Bi-O-edge driver requires zero warmed Julia heap allocation for one
complete CPU target-ready step. Direct GPU stream execution reports its Julia
launch overhead; zero host allocation, including completion, is a captured
CUDA/HIP Graph replay contract. Existing captured plant-graph drivers validate
that contract for the Shack-Hartmann graph. This Bi-O-edge driver exercises
direct streams and makes no Bi-O-edge graph-capture claim. Use the CUDA run as
the primary accelerator profile and the AMDGPU run as the secondary
portability profile. The measured replay is marked as
`aos_bi_o_edge_plant`; CUDA also uses the profiler API, and AMDGPU uses the
ROCTx profiler controls, so vendor tools can select only the warmed repeated
region.

Use the same procedure for the retained Zernike phase-spot plant boundary:

~~~bash
ADAPTIVEOPTICS_PROFILE_BACKEND=cpu ADAPTIVEOPTICS_PROFILE_STEPS=1000 \
  JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  julia --project=. --startup-file=no \
  scripts/profile_zernike_plant_runtime.jl

ADAPTIVEOPTICS_PROFILE_BACKEND=cuda ADAPTIVEOPTICS_PROFILE_STEPS=1000 \
  julia --project=test/cuda --startup-file=no \
  scripts/profile_zernike_plant_runtime.jl

ADAPTIVEOPTICS_PROFILE_BACKEND=amdgpu ADAPTIVEOPTICS_PROFILE_STEPS=1000 \
  julia --project=test/amdgpu --startup-file=no \
  scripts/profile_zernike_plant_runtime.jl
~~~

This driver includes only phase-spot optical formation and complete noiseless
detector acquisition. Registered FGA evidence covers Zernike pupil-signal
estimation and captured CUDA/HIP Graph replay. The AOS direct-stream driver
makes no standalone Zernike plant graph-capture claim. Its selected profiler
region is `aos_zernike_plant`.

The retained Curvature plant uses the same procedure and profiles complete
paired-defocus formation plus packed detector acquisition:

~~~bash
ADAPTIVEOPTICS_PROFILE_BACKEND=cpu ADAPTIVEOPTICS_PROFILE_STEPS=1000 \
  JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  julia --project=. --startup-file=no \
  scripts/profile_curvature_plant_runtime.jl

ADAPTIVEOPTICS_PROFILE_BACKEND=cuda ADAPTIVEOPTICS_PROFILE_STEPS=1000 \
  julia --project=test/cuda --startup-file=no \
  scripts/profile_curvature_plant_runtime.jl

ADAPTIVEOPTICS_PROFILE_BACKEND=amdgpu ADAPTIVEOPTICS_PROFILE_STEPS=1000 \
  julia --project=test/amdgpu --startup-file=no \
  scripts/profile_curvature_plant_runtime.jl
~~~

The AOS driver requires zero warmed CPU heap allocation and makes no standalone
plant graph-capture claim. Registered FGA 0.5 evidence covers paired-image and
paired-channel estimation, including captured CUDA/HIP Graph replay. The
selected AOS profiler region is `aos_curvature_plant`; CUDA is the primary
optimization target and AMDGPU is the secondary portability target.

## Release Entry Point

~~~bash
./scripts/run_release_validation.sh
~~~

Optional tracks are explicit:

~~~bash
ADAPTIVEOPTICS_VALIDATE_EXAMPLES=1 ./scripts/run_release_validation.sh
ADAPTIVEOPTICS_VALIDATE_AMDGPU=1 ./scripts/run_release_validation.sh
ADAPTIVEOPTICS_VALIDATE_CUDA=1 ./scripts/run_release_validation.sh
~~~

Run local validation before publishing. Use at most one justified exact-head
GitHub validation run for a PR, and record local commands/results in the PR.

## Acceptance Rule

Promote a backend claim only when the exact changed algorithm executes on the
named hardware, scalar indexing is disabled where applicable, residency and
completion are checked, and numerical/statistical comparison is within the
declared policy.

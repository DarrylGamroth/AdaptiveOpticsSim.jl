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
`backend-smoke` loads optional backend paths when their packages and devices are
available in the active environment; a skip is not hardware qualification.

Set Julia, BLAS, FFT-provider, OpenMP, and vendor math-library thread counts
deliberately. Deterministic evidence uses one thread. Performance evidence must
record the complete thread configuration and avoid nested oversubscription.

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

The hardware target also records and replays the qualified regular-grid DM as a
HIP Graph, changes the retained command buffer between replays, and verifies
that the following uncaptured graph node observes the new result on the same
stream. This is required evidence for `CapturedGraphExecution`; a standalone
HIP Graph smoke test is not sufficient.

## CUDA On WSL

Run the CUDA target on the WSL host:

~~~bash
ssh wsl 'cd /home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl &&
  julia --project=test/cuda --startup-file=no -e "using Pkg; Pkg.instantiate()" &&
  julia --project=test/cuda --startup-file=no test/runtests_cuda.jl'
~~~

For detector-only changes, use `test/runtests_cuda_detectors.jl`. Record Julia,
CUDA.jl, toolkit, driver, and GPU versions. CUDA is manually validated but is
not a continuously available release gate.

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
- capture performs no allocation, synchronization, result query, or host-side
  scientific-state mutation
- replay with changed input-buffer contents produces the corresponding changed
  output
- directly streamed nodes after a captured node observe same-stream ordering

The lockstep HIL boundary intentionally owns host `Array` exchange buffers. A
completed GPU detector frame is copied to that host buffer only after successful
graph execution, and a complete validated command is copied back before the next
frame.

Use `algorithm-graphs` for functional graph coverage and the REVOLT profiling
scripts for measured frame service time:

~~~bash
julia --project=. scripts/profile_revolt_hil_runtime.jl
julia --project=. benchmarks/benchmark_revolt_graph_nodes.jl
AOS_REVOLT_GRAPH_BACKEND=amdgpu \
  AOS_REVOLT_GRAPH_PROFILE=fast_dm \
  AOS_REVOLT_GRAPH_EXECUTION=captured \
  julia --project=benchmarks/amdgpu \
  benchmarks/benchmark_revolt_graph_nodes.jl
~~~

Record architecture, graph profile, resolution, noise settings, warmup, sample
count, and synchronization boundary. The graph-node benchmark records
submission, completion-ticket wait, and target-ready service time separately.
`AOS_REVOLT_GRAPH_EXECUTION=stream` is the default; use `captured` only with an
accelerator backend and a graph containing at least one qualified node.
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

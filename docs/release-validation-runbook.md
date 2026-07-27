# Release Validation Runbook

Status: active

## Purpose

This runbook defines the maintained release-validation entry point for
AdaptiveOpticsSim.jl.

It is intentionally operational: the goal is to make the defended validation
surface easy to rerun before a release or a production handoff.

## Primary Entry Point

Use:

- [run_release_validation.sh](../scripts/run_release_validation.sh)
- [archive_release_validation.sh](../scripts/archive_release_validation.sh)

From the repository root:

```bash
./scripts/run_release_validation.sh
```

To archive an operational validation run with dated logs and metadata:

```bash
./scripts/archive_release_validation.sh amdgpu
```

By default:

- `cpu` and `all` archive tracks run full CPU `Pkg.test()`
- `cuda` and `amdgpu` archive tracks skip the full CPU suite and run the backend-specific validation surface only

Optional validation tracks are enabled through environment flags:

```bash
ADAPTIVEOPTICS_VALIDATE_CUDA=1 ./scripts/run_release_validation.sh
ADAPTIVEOPTICS_VALIDATE_AMDGPU=1 ./scripts/run_release_validation.sh
ADAPTIVEOPTICS_VALIDATE_EXAMPLES=1 ./scripts/run_release_validation.sh
ADAPTIVEOPTICS_VALIDATE_COMPARISONS=1 ./scripts/run_release_validation.sh
ADAPTIVEOPTICS_VALIDATE_TRUTH=1 ./scripts/run_release_validation.sh
```

They may be combined:

```bash
ADAPTIVEOPTICS_VALIDATE_CUDA=1 \
ADAPTIVEOPTICS_VALIDATE_AMDGPU=1 \
ADAPTIVEOPTICS_VALIDATE_EXAMPLES=1 \
ADAPTIVEOPTICS_VALIDATE_COMPARISONS=1 \
ADAPTIVEOPTICS_VALIDATE_TRUTH=1 \
./scripts/run_release_validation.sh
```

To regenerate the maintained frozen OOPAO external-equivalence artifact:

```bash
julia --project=. --startup-file=no scripts/generate_oopao_equivalence_artifact.jl
```

To regenerate the maintained HEART boundary truth artifact:

```bash
python3 scripts/generate_heart_boundary_truth_artifact.py
```

## What Each Track Does

### CPU

Default for `cpu` and `all` archive tracks:

- `julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'`

This is the baseline production gate.

For a candidate that changes conjugate placement, path visibility, native
Plant DMs, sampled aberrations/NCPA, or the serial optical event boundary, also
run:

```bash
julia --threads=1 --project=. --startup-file=no -e \
  'using Pkg; Pkg.test(test_args=["gate5"])'
julia --threads=1 --project=benchmarks --startup-file=no \
  benchmarks/benchmark_gate5_optical_placement.jl
```

The benchmark is a serial self-paced CPU characterization. Its maintained
[contract](../benchmarks/contracts/gate5_optical_placement.toml) and
[artifact](../benchmarks/results/gate5/2026-07-25-optical-placement.toml)
close the Gate 5 numerical, fixed-storage, bounded-allocation, and bounded
4/8/16-path topology review without asserting fixed-arrival HIL latency or
production instrument capacity.

The CPU full suite may be skipped only when `ADAPTIVEOPTICS_SKIP_CPU_FULL_TESTS=1` is set explicitly. That mode is intended for backend-host validation runs that are paired with separately archived CPU/full-suite evidence for the same candidate commit or an explicitly identified release ancestor.

### Apple Silicon / AppleAccelerate

The default CPU package remains linear-algebra-provider neutral. On Apple
Silicon, the hosted CPU workflow additionally instantiates
[`test/appleaccelerate`](../test/appleaccelerate), proves that normal package
loading does not load AppleAccelerate, then loads AppleAccelerate explicitly and
runs the full CPU suite. The target verifies that representative BLAS and
LAPACK symbols route through Accelerate, supported power-of-two 1D/2D CPU FFTs
use allocation-free reusable vDSP plans, and unsupported shapes retain FFTW
fallback plans.

To reproduce this target on a macOS 15 or newer Apple Silicon host:

```bash
julia --project=test/appleaccelerate --startup-file=no -e 'using Pkg; Pkg.instantiate()'
julia --project=test/appleaccelerate --startup-file=no test/appleaccelerate/backend_neutral.jl
julia --project=test/appleaccelerate --startup-file=no test/appleaccelerate/runtests.jl
```

Applications opt in by depending on and loading AppleAccelerate themselves;
AdaptiveOpticsSim declares only a weak dependency and does not auto-load it or
include it in the root environment.

### CUDA

Enabled with:

- `ADAPTIVEOPTICS_VALIDATE_CUDA=1`

Runs:

- [runtests_cuda.jl](../test/runtests_cuda.jl) with the
  [`test/cuda`](../test/cuda) project

Use this on a CUDA-capable host. The validation script instantiates the
backend-specific `test/cuda` project so `CUDA.jl` does not need to be installed
in the root package environment. The archived `cuda` track defaults to
backend-only validation by setting
`ADAPTIVEOPTICS_SKIP_CPU_FULL_TESTS=1`.

The CUDA track has a current manual WSL validation host. Gate 5 closure
validation head `02e5f29` passed the maintained `438/438` hardware target on
Julia 1.12.6 with CUDA.jl 6.2.1 and KernelAbstractions.jl 0.9.42, on an RTX
3050 Ti with compute capability 8.6, with scalar indexing disabled. This
includes exact device-resident defensive sampled-OPD ownership and identity
plus transformed replacement. The separate
[final pre-HIL backend evidence](../benchmarks/results/platform/2026-07-18-pre-hil-11-wsl-cuda.toml)
archives CPU parity, residency, service-time histograms, and backend-ready,
host-ready, and transfer-only boundaries. CUDA is still outside the present
production support claim until it is explicitly restored to the supported
delivery scope and assigned a routine validation cadence.

### AMDGPU

Enabled with:

- `ADAPTIVEOPTICS_VALIDATE_AMDGPU=1`

Runs:

- [runtests_amdgpu.jl](../test/runtests_amdgpu.jl) with the
  [`test/amdgpu`](../test/amdgpu) project

Use this on an AMDGPU-capable host. The validation script instantiates the
backend-specific `test/amdgpu` project so `AMDGPU.jl` does not need to be
installed in the root package environment. The archived `amdgpu` track defaults
to backend-only validation by setting
`ADAPTIVEOPTICS_SKIP_CPU_FULL_TESTS=1`.

Gate 5 closure validation head `02e5f29` passed `448/448` maintained checks on
the local gfx1030 AMD Radeon Graphics target with scalar indexing disabled and
Julia 1.12.6. This includes direct device-resident native Plant DM state,
staging storage, and surface parity plus exact sampled-aberration defensive
ownership and identity/transformed replacement; it does not exercise an
integrated GPU event loop. The July 14 REVOLT-like
[cross-host artifact](../benchmarks/results/platform/2026-07-14-wsl-cuda-local-amdgpu.toml)
remains a non-equivalent historical workload. The current paired
single-device owner characterization is retained in the
[Gate 7 catalog](../benchmarks/results/gate7/manifest.toml).

### Gate 7 single-GPU closure benchmark

When a release changes compute-device identity, direction batching, WFS device
owners, synchronization, or explicit host-observation boundaries, rerun the
predeclared Gate 7 benchmark on the CPU oracle and each maintained accelerator:

```bash
AOS_GATE7_OUTPUT=/tmp/gate7-local-cpu.toml \
  julia --threads=1 --project=benchmarks --startup-file=no \
  benchmarks/benchmark_gate7_single_gpu.jl cpu local_cpu
AOS_GATE7_OUTPUT=/tmp/gate7-local-amdgpu.toml \
  julia --threads=1 --project=benchmarks/amdgpu --startup-file=no \
  benchmarks/benchmark_gate7_single_gpu.jl amdgpu local_amdgpu
AOS_GATE7_OUTPUT=/tmp/gate7-wsl-cuda.toml \
  julia --threads=1 --project=benchmarks/cuda --startup-file=no \
  benchmarks/benchmark_gate7_single_gpu.jl cuda wsl_cuda
```

Run each command from the same clean candidate revision without sample, run,
or warmup overrides. The
[contract](../benchmarks/contracts/gate7_single_gpu.toml) and
[artifact catalog](../benchmarks/results/gate7/manifest.toml) define the
retained workload, hashes, hardware identity, gates, and scope. This benchmark
is required single-device service-cost evidence, but it does not replace the
dedicated hardware correctness matrices and must not be reported as
fixed-arrival HIL latency or multi-GPU validation.

### Core examples

Enabled with:

- `ADAPTIVEOPTICS_VALIDATE_EXAMPLES=1`

Runs:

- [run_core_examples.sh](../scripts/run_core_examples.sh)

This track executes the maintained plotting-free example scripts in
`examples/closed_loop` and `examples/tutorials`. It is intentionally separate
from `Pkg.test()` so examples can be used as a user-facing smoke/regression
surface without making the normal unit-test path slower.

### Cross-package comparisons

Enabled with:

- `ADAPTIVEOPTICS_VALIDATE_COMPARISONS=1`

Runs the maintained HEART all-package ladder in the sibling comparison
workspace when it exists:

- `../AdaptiveOpticsComparisons/contracts/heart_hil.toml`

If the sibling comparison workspace is absent, this track skips cleanly.

### Scientist-owned HEART truth

Enabled with:

- `ADAPTIVEOPTICS_VALIDATE_TRUTH=1`

Runs:

- [generate_heart_boundary_truth_artifact.py](../scripts/generate_heart_boundary_truth_artifact.py)

Use this when the sibling `REVOLT` checkout is present and the release story
needs the maintained scientist-owned HEART boundary artifact refreshed.

## Interpretation

Before a release or production handoff:

1. CPU validation must pass.
2. AMDGPU validation must pass because it is the currently supported GPU
   delivery scope.
3. CUDA validation must pass on real hardware before CUDA can be returned to
   the supported delivery scope.
4. Cross-package HEART comparison should be rerun when external equivalence
   claims are part of the release story.
5. Scientist-owned HEART truth should be rerun when boundary-truth claims are
   part of the release story.

Use together with:

- [supported-production-surfaces.md](supported-production-surfaces.md)
- [backend-validation-guide.md](backend-validation-guide.md)

## Validation host bootstrap

Before trusting a new CUDA or AMDGPU validation host, bootstrap it with:

- [bootstrap_validation_host.sh](../scripts/bootstrap_validation_host.sh)

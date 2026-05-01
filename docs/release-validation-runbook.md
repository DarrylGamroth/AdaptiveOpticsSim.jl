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
./scripts/archive_release_validation.sh cuda spiders
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

The CPU full suite may be skipped only when `ADAPTIVEOPTICS_SKIP_CPU_FULL_TESTS=1` is set explicitly. That mode is intended for backend-host validation runs that are paired with separately archived CPU/full-suite evidence for the same candidate commit or an explicitly identified release ancestor.

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
2. CUDA validation must pass if CUDA is in the supported delivery scope.
3. AMDGPU validation must pass if AMDGPU is in the supported delivery scope.
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

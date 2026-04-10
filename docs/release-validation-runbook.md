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

This always runs:

- full CPU `Pkg.test()`

Optional validation tracks are enabled through environment flags:

```bash
ADAPTIVEOPTICS_VALIDATE_CUDA=1 ./scripts/run_release_validation.sh
ADAPTIVEOPTICS_VALIDATE_AMDGPU=1 ./scripts/run_release_validation.sh
ADAPTIVEOPTICS_VALIDATE_COMPARISONS=1 ./scripts/run_release_validation.sh
ADAPTIVEOPTICS_VALIDATE_TRUTH=1 ./scripts/run_release_validation.sh
```

They may be combined:

```bash
ADAPTIVEOPTICS_VALIDATE_CUDA=1 \
ADAPTIVEOPTICS_VALIDATE_AMDGPU=1 \
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

Always runs:

- `julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'`

This is the baseline production gate.

### CUDA

Enabled with:

- `ADAPTIVEOPTICS_VALIDATE_CUDA=1`

Runs:

- [gpu_smoke_cuda.jl](../scripts/gpu_smoke_cuda.jl)
- [gpu_runtime_equivalence_cuda.jl](../scripts/gpu_runtime_equivalence_cuda.jl)

Use this on a CUDA-capable host with `CUDA.jl` available in the project
environment.

### AMDGPU

Enabled with:

- `ADAPTIVEOPTICS_VALIDATE_AMDGPU=1`

Runs:

- [gpu_smoke_amdgpu.jl](../scripts/gpu_smoke_amdgpu.jl)
- [gpu_runtime_equivalence_amdgpu.jl](../scripts/gpu_runtime_equivalence_amdgpu.jl)

Use this on an AMDGPU-capable host with `AMDGPU.jl` available in the project
environment.

### Cross-package comparisons

Enabled with:

- `ADAPTIVEOPTICS_VALIDATE_COMPARISONS=1`

Runs the maintained HEART all-package ladder in the sibling comparison
workspace when it exists:

- [/home/dgamroth/workspaces/codex/AdaptiveOpticsComparisons/contracts/heart_hil.toml](/home/dgamroth/workspaces/codex/AdaptiveOpticsComparisons/contracts/heart_hil.toml)

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

- [supported-production-surfaces.md](./supported-production-surfaces.md)
- [production-readiness-checklist.md](./production-readiness-checklist.md)
- [backend-validation-guide.md](./backend-validation-guide.md)
- [self-hosted-gpu-runner-setup.md](./self-hosted-gpu-runner-setup.md)
- [operational-gpu-validation-cadence.md](./operational-gpu-validation-cadence.md)

## Validation host bootstrap

Before trusting a new CUDA or AMDGPU validation host, bootstrap it with:

- [bootstrap_validation_host.sh](../scripts/bootstrap_validation_host.sh)

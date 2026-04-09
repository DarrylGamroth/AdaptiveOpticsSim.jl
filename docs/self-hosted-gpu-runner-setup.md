# Self-Hosted GPU Runner Setup

Status: active

## Purpose

This document defines the maintained setup for self-hosted GPU runners used by
the checked-in CUDA and AMDGPU validation workflows.

Use together with:

- [release-validation-runbook.md](./release-validation-runbook.md)
- [backend-validation-guide.md](./backend-validation-guide.md)
- [supported-production-surfaces.md](./supported-production-surfaces.md)

## Expected GitHub Actions labels

The checked-in workflows expect:

- CUDA runner:
  - `self-hosted`
  - `linux`
  - `cuda`
- AMDGPU runner:
  - `self-hosted`
  - `linux`
  - `amdgpu`

## Host prerequisites

Each runner host should have:

- a working Julia installation compatible with the repository `Project.toml`
- working Git checkout access to this repository
- the GitHub Actions runner service installed and registered
- a functioning GPU driver/runtime for the target backend
- enough local disk for Julia artifact and package caches

Backend-specific requirements:

- CUDA host:
  - working NVIDIA driver
  - `nvidia-smi` available
- AMDGPU host:
  - working ROCm stack
  - `rocminfo` available

## Maintained bootstrap entry point

From the repository root, use:

- [bootstrap_validation_host.sh](../scripts/bootstrap_validation_host.sh)

Examples:

```bash
./scripts/bootstrap_validation_host.sh
./scripts/bootstrap_validation_host.sh cuda
./scripts/bootstrap_validation_host.sh amdgpu
./scripts/bootstrap_validation_host.sh cuda amdgpu
```

This script:

- instantiates the project
- precompiles the project
- installs weak GPU dependencies when requested
- verifies the expected backend runtime tool is visible

## Runner registration guidance

Register one runner per host or backend-specific host group.

Suggested naming:

- `adaptiveopticssim-cuda-01`
- `adaptiveopticssim-amdgpu-01`

Suggested workflow:

1. Clone the repository onto the runner host.
2. Run [bootstrap_validation_host.sh](../scripts/bootstrap_validation_host.sh)
   with the backend you want that host to serve.
3. Register the GitHub Actions runner with the expected labels.
4. Execute the matching local release-validation track once before trusting CI.

Examples:

```bash
./scripts/bootstrap_validation_host.sh cuda
ADAPTIVEOPTICS_VALIDATE_CUDA=1 ./scripts/run_release_validation.sh
```

```bash
./scripts/bootstrap_validation_host.sh amdgpu
ADAPTIVEOPTICS_VALIDATE_AMDGPU=1 ./scripts/run_release_validation.sh
```

## Validation rule

Do not treat the checked-in GPU workflows as production-ready automation until:

- the runner is registered with the expected labels
- the bootstrap script succeeds on the host
- the matching local validation track succeeds
- the GitHub workflow has gone green on that host

## Operational notes

- Keep runner-local Julia caches warm to reduce workflow noise.
- Prefer dedicated GPU runners over mixed interactive hosts.
- If a backend update changes toolchain behavior, rerun the matching bootstrap
  and release-validation flow before trusting old green results.

# Operational GPU Validation Cadence

Status: active

Use together with:

- [release-validation-runbook.md](./release-validation-runbook.md)
- [production-readiness-checklist.md](./production-readiness-checklist.md)
- [supported-production-surfaces.md](./supported-production-surfaces.md)
- [self-hosted-gpu-runner-setup.md](./self-hosted-gpu-runner-setup.md)
- [../scripts/archive_release_validation.sh](../scripts/archive_release_validation.sh)

## Purpose

This note defines the minimum operational cadence for real-hardware GPU
validation when production claims include CUDA and/or AMDGPU support.

This cadence can be satisfied either by routine CI on self-hosted GPU runners or
by a maintained release-validation host procedure. The current supported path is
release-validation hosts.

## Required Cadence

For every supported release candidate or production handoff:

1. archive at least one CPU/full-suite validation run for the candidate commit or an explicitly identified release ancestor
2. run CUDA validation on a real CUDA host if CUDA is in scope
3. run AMDGPU validation on a real AMDGPU host if AMDGPU is in scope
4. archive the resulting logs and run metadata

Minimum frequency:

- once per tagged release candidate, or
- once per production handoff build if no tagging workflow is used

## Canonical Entry Point

Use:

- [../scripts/archive_release_validation.sh](../scripts/archive_release_validation.sh)

Examples:

```bash
./scripts/archive_release_validation.sh cpu
./scripts/archive_release_validation.sh amdgpu
./scripts/archive_release_validation.sh cuda spiders
```

This writes dated log and metadata artifacts to:

- [../benchmarks/results/validation_runs](../benchmarks/results/validation_runs)

## Host Mapping

Current expected operational mapping:

- AMDGPU validation host: a local AMDGPU-capable validation machine
- CUDA validation host: `spiders`

If that mapping changes, update this note and the release runbook together.

## Required Artifacts

Each operational validation run must produce:

- a full console log
- a small metadata record with:
  - host label
  - track
  - git commit
  - start/end timestamps
  - pass/fail status

## Release Decision Rule

A release may claim support for a backend only if:

- there is archived CPU/full-suite evidence for the candidate commit or an explicitly identified release candidate ancestor
- the matching backend operational run passed on real hardware for that same commit or ancestor

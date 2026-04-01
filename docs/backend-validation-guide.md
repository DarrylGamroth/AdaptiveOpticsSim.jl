# Backend Validation Guide

Status: active

Plan traceability:

- [`PLAN-15`](./package-review-action-plan.md)
- [`PLAN-32`](./package-review-action-plan.md)
- review IDs: `PR-31`, `PR-32`, `PR-33`

## Purpose

This document explains how backend validation is organized.

The goal is to keep three distinct classes of evidence separate:

- functional tests in `Pkg.test()`
- optional backend smoke/parity checks in `Pkg.test()`
- benchmark and profile evidence outside `Pkg.test()`

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
- curvature-through-atmosphere

For the full maintained GPU smoke matrix, use:

- `ADAPTIVEOPTICS_TEST_FULL_AMDGPU=1`
- `ADAPTIVEOPTICS_TEST_FULL_CUDA=1`

which route through [`gpu_smoke_contract.jl`](../scripts/gpu_smoke_contract.jl).

## Benchmark Separation

Benchmarks do not belong in `Pkg.test()`.

Representative runtime evidence should be gathered with maintained profile or
benchmark scripts such as:

- [`profile_ao3k_runtime.jl`](../scripts/profile_ao3k_runtime.jl)
- [`profile_multi_source_multi_wfs_runtime.jl`](../scripts/profile_multi_source_multi_wfs_runtime.jl)
- [`run_cross_package_benchmarks.jl`](../scripts/run_cross_package_benchmarks.jl)
- [`benchmarks/`](../benchmarks)

This separation exists so:

- unit tests stay reasonably fast
- optional GPU availability does not block normal development
- benchmark evidence remains intentional and archived rather than implicit

## Validation Rule

Claiming backend support should imply:

- functional correctness in the main testsets
- optional backend smoke coverage in `Pkg.test()`
- maintained benchmark or profile evidence where runtime claims are made

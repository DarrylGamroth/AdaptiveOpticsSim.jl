# GPU Backend Architecture

This document defines how GPU support should be structured in
`AdaptiveOpticsSim.jl`.

The core package should remain backend-generic. Backend-specific glue,
capability detection, and validation should live in extension modules.

## Goals

- Keep `AdaptiveOpticsSim` core generic over array backends.
- Preserve the current `KernelAbstractions`-based execution model.
- Support multiple GPU backends without baking vendor packages into core.
- Make backend support explicit and testable rather than implied.
- Keep the CPU path first-class and efficient.

## Non-Goals

- Do not force every code path through a GPU backend.
- Do not claim a backend is supported unless it has backend-specific smoke
  coverage.
- Do not move setup-only logic to GPU unless it is measurably costly.

## Current State

Today the core already does most of the right things:

- array storage is backend-parametric
- execution dispatch uses `ScalarCPUStyle` and `AcceleratorStyle`
- compute kernels use `KernelAbstractions`
- FFT usage goes through `AbstractFFTs`

In practice, only CUDA is currently validated. The package is not yet organized
around GPU extension modules, but it should be.

## Recommended Package Split

### Core Package

`AdaptiveOpticsSim` core should own:

- backend-generic types
- execution traits and dispatch
- `KernelAbstractions` kernels
- generic FFT interface points
- generic RNG/noise interface points
- backend-agnostic smoke-test contracts

Core should not depend directly on:

- `CUDA.jl`
- `Metal.jl`
- `AMDGPU.jl`
- `OpenCL.jl`

### Extension Modules

Each GPU backend should be integrated through a Julia package extension:

- `ext/AdaptiveOpticsSimCUDAExt.jl`
- `ext/AdaptiveOpticsSimMetalExt.jl`
- `ext/AdaptiveOpticsSimAMDGPUExt.jl`

Each extension should own only backend-specific glue:

- backend availability checks
- array-type convenience constructors
- backend-specific FFT hookups if needed
- backend-specific RNG hooks if needed
- backend-specific smoke helpers
- backend-specific policy overrides where generic behavior is insufficient

## Proposed Project Layout

```text
AdaptiveOpticsSim.jl/
  src/
    AdaptiveOpticsSim.jl
    Core/
      backends.jl
      fft.jl
      rng.jl
      utils.jl
    ...
  ext/
    AdaptiveOpticsSimCUDAExt.jl
    AdaptiveOpticsSimMetalExt.jl
    AdaptiveOpticsSimAMDGPUExt.jl
  scripts/
    gpu_smoke_matrix.jl
    gpu_smoke_cuda.jl
    gpu_smoke_metal.jl
    gpu_smoke_amdgpu.jl
  test/
    gpu/
      smoke_matrix.jl
      backend_contracts.jl
```

## Proposed `Project.toml` Shape

This is the intended direction, not a required immediate edit.

```toml
[weakdeps]
CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"

[extensions]
AdaptiveOpticsSimCUDAExt = "CUDA"
AdaptiveOpticsSimMetalExt = "Metal"
AdaptiveOpticsSimAMDGPUExt = "AMDGPU"
```

## Trait Model

The current execution-style split is good, but GPU support needs a slightly
broader trait surface.

Core should expose generic trait hooks like:

- `execution_style(array)`
- `fft_provider(array_or_backend)`
- `rng_provider(array_or_backend)`
- `transfer_policy(array_or_backend)`
- `backend_name(array_or_backend)`

Example sketch:

```julia
abstract type FFTProvider end
struct GenericFFTProvider <: FFTProvider end
struct CUFFTProvider <: FFTProvider end
struct MetalFFTProvider <: FFTProvider end
struct ROCFFTProvider <: FFTProvider end
fft_provider(::Any) = GenericFFTProvider()
```

Then each extension can add methods without contaminating core:

```julia
fft_provider(::CUDA.CuArray) = CUFFTProvider()
```

The same pattern should be used for RNG/noise when backend-native
implementations differ.

## Extension Responsibilities

### `AdaptiveOpticsSimCUDAExt`

This should be the first fully supported backend extension.

Responsibilities:

- `CuArray` convenience helpers
- CUDA-specific smoke entry point
- CUDA-specific FFT/provider hooks if generic `plan_fft!` is not sufficient
- CUDA-specific capability assertions used in CI
- optional CUDA-specific profiling helpers

Status:

- should be considered the reference GPU backend

### `AdaptiveOpticsSimMetalExt`

Responsibilities:

- `MtlArray` convenience helpers
- Metal-specific smoke entry point
- Metal-specific FFT/provider hooks if required
- explicit documentation of any unsupported precision or FFT limitations

Status:

- plausible target
- not validated yet

### `AdaptiveOpticsSimAMDGPUExt`

Responsibilities:

- `ROCArray` convenience helpers
- AMD-specific smoke entry point
- ROCm FFT/provider hooks if required
- explicit CI and driver assumptions

Status:

- plausible target
- not validated yet

### OpenCL

OpenCL is deferred for now.

Reason:

- the current priority is a solid extension-based path for CUDA, Metal, and
  AMDGPU
- OpenCL ecosystem maturity for FFT-heavy AO workflows is less certain
- it should not dilute the first multi-backend bring-up pass

If revisited later, it should be treated as experimental until the Julia
OpenCL + `KernelAbstractions` + FFT stack is proven stable for this workload.

## Runtime Support Levels

We should define support levels clearly.

### Level 0: Compiles

- package loads
- arrays construct
- no validation claim

### Level 1: Kernel-Compatible

- `KernelAbstractions` kernels run
- no scalar indexing errors in targeted smoke cases

### Level 2: Runtime-Validated

- device-resident runtime smoke passes
- noisy detector path passes
- representative closed-loop path passes

### Level 3: Supported Backend

- CI coverage exists
- smoke matrix passes reliably
- FFT/RNG behavior is documented
- known limitations are documented

Only Level 3 should be described as supported in user-facing docs.

## FFT Strategy

Core should continue to depend on `AbstractFFTs`, but that is not enough by
itself to claim backend portability. Each backend must prove:

- plan creation works for that array type
- repeated planned execution works
- centered/cropped PSF and WFS paths remain correct

If a backend needs custom FFT setup, that should be added in its extension
rather than branching inside core algorithms.

## RNG and Noise Strategy

Noise generation should stay behind generic core entry points:

- `randn_backend!`
- `poisson_noise!`

If the generic stateless KA kernels are good enough, extensions do not need to
override them. If a backend offers a better native implementation, the
extension can add specialized methods.

That keeps the core portable while still allowing tuned backend-native RNG.

## Validation Layout

The existing `scripts/gpu_smoke_matrix.jl` should evolve into:

- one backend-generic smoke contract
- thin backend launchers per extension/backend

For example:

- `scripts/gpu_smoke_matrix.jl`
  - defines the cases and expected surfaces
- `scripts/gpu_smoke_cuda.jl`
  - loads CUDA, selects `CuArray`, runs the matrix
- `scripts/gpu_smoke_metal.jl`
  - loads Metal, selects `MtlArray`, runs the matrix
- `scripts/gpu_smoke_amdgpu.jl`
  - loads AMDGPU, selects `ROCArray`, runs the matrix
- `scripts/gpu_smoke_opencl.jl`
  - loads OpenCL backend, if available

CI should run backend-specific smoke jobs only when the corresponding runner is
available.

## Recommendation Relative to `proper.jl`

`proper.jl` already uses a dedicated `ProperCUDAExt.jl`, which is directionally
correct. For `AdaptiveOpticsSim.jl`, the same idea should be kept, but with a
slightly stricter split:

- core owns all generic traits and kernels
- extensions own only backend glue and validation
- support claims are made per backend, not globally

That is a better fit for a multi-backend simulation package than a
CUDA-first-only extension story.

## Recommended Next Steps

1. Add the backend extension stubs in `ext/` without changing core behavior.
2. Move CUDA-specific smoke launching into `AdaptiveOpticsSimCUDAExt`.
3. Add a backend-generic smoke contract and backend-specific wrappers.
4. Treat Metal and AMDGPU as explicit bring-up targets, not implicit KA
   promises.
5. Revisit OpenCL only after the extension-based CUDA/Metal/AMDGPU layout is
   established and validated.

# Deterministic Simulation Mode (Design Sketch)

Status: active reference design

This document outlines a deterministic mode for the Julia port so results are
reproducible run-to-run with identical inputs and configuration.

## Goals
- Repeatable outputs with fixed seeds and fixed configuration.
- Traceability: capture the exact inputs and configuration that produced a run.
- Minimal changes to the normal API.

## Non-goals
- Cross-hardware bitwise identity (CPU vs GPU, different CPUs).
- Hard real-time determinism.

## Sources of non-determinism
- Global RNG usage or implicit `rand()` calls.
- FFT planning differences and multithreading.
- Thread scheduling and non-associative floating-point reductions.
- GPU kernels and non-deterministic reductions.

## Design overview
1) Centralized RNG:
   - All stochastic components accept an explicit `rng` (phase screens, photon
     noise, read noise, source jitter).
   - `Workspace` owns the RNG so every step uses the same stream.
   - Use `deterministic_reference_rng(seed)` for reference-data generation,
     regression fixtures, and examples where preserving historical streams
     matters. This currently returns `MersenneTwister(seed)`.
   - Use `runtime_rng(seed)` for new long-running runtime or HIL simulations
     where the stream only needs to be repeatable for a fixed Julia/runtime
     configuration. This currently returns `Xoshiro(seed)`, which has smaller
     state and is generally a better throughput-oriented choice for hot paths.
   - For cross-version reproducibility, consider `StableRNGs.jl`.

2) Deterministic configuration:
   - Fixed seed.
   - Single-threaded FFT provider and BLAS.
   - Optional fixed FFT-provider wisdom/plan where supported.

3) Record and replay:
   - Optional logging of commands and RNG state.
   - Replay mode that bypasses stochastic draws.

## API sketch
```julia
struct DeterministicConfig
    seed::UInt64
    threads::Int
    fft_threads::Int
    use_gpu::Bool
end

function init_deterministic!(cfg::DeterministicConfig)
    Random.seed!(cfg.seed)
    BLAS.set_num_threads(cfg.threads)
    set_fft_provider_threads!(cfg.fft_threads)
end

ws = AdaptiveOpticsSim.Workspace(tel; rng=deterministic_reference_rng(0x1234))
runtime = ClosedLoopRuntime(sim, recon; rng=runtime_rng(0x1234))
```

## RNG Policy

The package accepts `AbstractRNG` through stochastic APIs rather than owning a
single global stream. That keeps detector noise, phase screens, and runtime
simulation reproducible when users pass an explicit RNG.

Use `MersenneTwister` through `deterministic_reference_rng` when changing the
stream would invalidate stored references or regression baselines. Use `Xoshiro`
through `runtime_rng` for new RTC/HIL-style simulations and benchmarks where
throughput and lower RNG state overhead are more important than preserving an
older fixture stream.

Do not rely on `Random.default_rng()` for validation artifacts or reference
comparisons. Always pass the RNG explicitly.

## Phase screens
- Use a deterministic phase screen generator with fixed `rng` draws.
- Optionally persist the phase screen sequence to disk (HDF5/FITS) for replay.

## Noise models
- `Detector` noise uses the shared RNG and the noise model trait.
- Deterministic baseline: `Detector(noise=NoiseNone())` creates `Detector{NoiseNone}`.

## Testing approach
- Unit test: two runs with identical seed produce equal PSF arrays.
- Regression test: compare stored reference outputs within tolerance.

## Notes on GPU
GPU determinism is typically weaker than CPU determinism. In deterministic mode,
default to CPU unless a GPU-specific deterministic path is available.

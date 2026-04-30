# Determinism And RNG Policy

Status: active

This document defines the maintained determinism and RNG policy for
AdaptiveOpticsSim.jl. Results should be reproducible run-to-run when users pass
the same inputs, configuration, RNG seeds, and execution settings.

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

Example:

```julia
rng = runtime_rng(0x1234)
advance!(atm, tel; rng=rng)
measure!(wfs, tel, src, det; rng=rng)
```

For reference-data refreshes:

```julia
rng = deterministic_reference_rng(0x1234)
```

## Deterministic Configuration

- Use fixed seeds.
- Run with one Julia thread, one BLAS thread, and one FFT-provider thread when
  strict reproducibility matters.
- Keep detector noise disabled for deterministic baseline comparisons unless
  the test is explicitly about noise.
- Capture commands, configuration, and RNG seed/state in validation artifacts.
- For cross-version bitstream stability, consider `StableRNGs.jl`.

## Phase Screens

- Use a deterministic phase screen generator with fixed `rng` draws.
- Optionally persist the phase screen sequence to disk (HDF5/FITS) for replay.

## Noise Models

- `Detector` noise uses the shared RNG and the noise model trait.
- Deterministic baseline: `Detector(noise=NoiseNone())` creates `Detector{NoiseNone}`.

## Testing Approach

- Unit test: two runs with identical seed produce equal PSF arrays.
- Regression test: compare stored reference outputs within tolerance.

## GPU Notes

GPU determinism is typically weaker than CPU determinism. In deterministic mode,
default to CPU unless a GPU-specific deterministic path is available.

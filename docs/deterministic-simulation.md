# Deterministic Simulation Mode (Design Sketch)

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
   - For cross-version reproducibility, consider `StableRNGs.jl`.

2) Deterministic configuration:
   - Fixed seed.
   - Single-threaded FFT and BLAS.
   - Optional fixed FFTW wisdom/plan.

3) Record and replay:
   - Optional logging of commands and RNG state.
   - Replay mode that bypasses stochastic draws.

## API sketch
```julia
struct DeterministicConfig
    seed::UInt64
    threads::Int
    fftw_threads::Int
    use_gpu::Bool
end

function init_deterministic!(cfg::DeterministicConfig)
    Random.seed!(cfg.seed)
    BLAS.set_num_threads(cfg.threads)
    FFTW.set_num_threads(cfg.fftw_threads)
end

ws = Workspace(tel; rng=MersenneTwister(0x1234))
```

## Phase screens
- Use a deterministic phase screen generator with fixed `rng` draws.
- Optionally persist the phase screen sequence to disk (HDF5/FITS) for replay.

## Noise models
- `Detector` noise uses the shared RNG and the noise model trait.
- Deterministic baseline: `Detector(photon_noise=false, readout_noise=0.0)` creates `Detector{NoiseNone}`.

## Testing approach
- Unit test: two runs with identical seed produce equal PSF arrays.
- Regression test: compare stored reference outputs within tolerance.

## Notes on GPU
GPU determinism is typically weaker than CPU determinism. In deterministic mode,
default to CPU unless a GPU-specific deterministic path is available.

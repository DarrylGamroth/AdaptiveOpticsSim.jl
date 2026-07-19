# Determinism And RNG Policy

Status: active

This document defines the maintained determinism and RNG policy for
AdaptiveOpticsSim.jl. Results should be reproducible run-to-run when users pass
the same inputs, configuration, RNG seeds, and execution settings.

## Goals
- Repeatable outputs with fixed seeds and fixed configuration.
- Traceability: capture the exact inputs and configuration that produced a run.
- Stable stochastic ownership that does not depend on tuple order, task
  scheduling, thread identity, or CPU/GPU placement.

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
atmosphere_rng = runtime_rng(0x1234)
detector_rng = runtime_rng(0x5678)
epoch = advance_by!(atm, 1e-3; rng=atmosphere_rng)
render_atmosphere!(pupil, renderer, atm, epoch)
measure!(wfs, pupil, src, det; rng=detector_rng)
```

Rendering a published `AtmosphereEpoch` is deterministic and does not consume
RNG. A zero-duration advance after initialization returns the current epoch
without changing the RNG stream or epoch sequence. For order-independent
multi-path replay, prepare one renderer per direction and render the same
current epoch in any order before that atmosphere advances, or replay from
separately retained model state.

For reference-data refreshes:

```julia
rng = deterministic_reference_rng(0x1234)
```

## Scheduled-Plant RNG Topology

The multi-rate plant uses one central run seed plus a versioned derivation
scheme and stable RNG owner identities. Preparation assigns an identity to
every stochastic atmosphere layer, trigger source or link, optical or reduced-
order provider, detector acquisition, and device model. Identities come from
the declared plant topology and remain stable when prepared execution groups
are reordered or moved between CPU and GPU resources. Duplicate owner
identities are a preparation error.

Each owner uses one of two prepared policies:

- a stateful stream with exactly one execution writer when its draw order is
  part of the model; or
- an addressable random domain derived from run seed, derivation version,
  owner identity, event or epoch sequence, and element/sample index when work
  may be replicated, batched, reordered, or evaluated on several devices.

A task ID, `Threads.threadid()`, tuple position, physical device ordinal, ring
cursor, or completion order is never an RNG identity. Changing static
placement must not change which random values belong to a physical event.
Scenario replay records the run seed, derivation version, owner-identity map,
and any stateful-stream checkpoints needed by the selected model.

Existing stochastic APIs continue to accept explicit `AbstractRNG` values.
The scheduled plant prepares and supplies those per-owner values; individual
models do not look up a global registry in their hot path. Addressable GPU
kernels receive prepared keys/counters directly rather than consuming one
host RNG seed in launch order.

## Deterministic Configuration

- Use fixed seeds.
- Run with one Julia thread, one BLAS thread, and one FFT-provider thread when
  strict reproducibility matters.
- Use `DeterministicExecution()` for a `SimulationEnsemble`. It rejects a
  multi-threaded Julia process and configures BLAS and the FFT provider for one
  thread before execution.
- Keep detector noise disabled for deterministic baseline comparisons unless
  the test is explicitly about noise.
- Capture commands, configuration, and RNG seed/state in validation artifacts.
- For cross-version bitstream stability, consider `StableRNGs.jl`.

## Phase Screens

- Use a deterministic phase screen generator with fixed `rng` draws.
- Optionally persist the phase screen sequence to disk (HDF5/FITS) for replay.

## Noise Models

- `Detector` noise uses the explicitly supplied RNG and the noise model trait.
- Deterministic baseline: `Detector(noise=NoiseNone())` creates `Detector{NoiseNone}`.
- `SharedOpticalRuntime` currently gives its science-detector tuple sequential
  draws from one runtime RNG. Reordering noisy detectors therefore changes
  their random streams; use `NoiseNone()` for order-invariant optical fan-out
  evidence. This is a documented limitation of the transitional frame-step
  runtime, not the target scheduled-plant contract.

## Testing Approach

- Unit test: two runs with identical seed produce equal caller-owned direct
  images or WFS products.
- Ownership test: reordering independent acquisitions or execution groups does
  not change their stochastic products.
- Placement test: serial CPU, grouped CPU, and applicable replicated GPU
  execution select the same event-addressed random values within the declared
  numerical comparison policy.
- Regression test: compare stored reference outputs within tolerance.

## GPU Notes

GPU determinism is typically weaker than CPU determinism. In deterministic
mode, default to CPU unless a GPU-specific deterministic path is available.
Multi-device support additionally requires addressable or otherwise replicated
per-owner random state; sequential host seed consumption is insufficient.

# Determinism And RNG Policy

Status: active

AdaptiveOpticsSim.jl results should be reproducible when inputs, configuration,
RNG seeds, and execution settings are fixed.

## Goals

- Repeatable outputs for a fixed run configuration.
- Explicit provenance for configuration, inputs, seeds, and backend.
- One writer for every evolving RNG state.
- Random ownership that does not depend on task, thread, tuple, or completion
  order.

Cross-hardware bitwise identity and hard real-time determinism are not claimed.

## RNG Policy

Stochastic scientific APIs accept an explicit `AbstractRNG`. Do not use the
global/default RNG for validation artifacts or reference comparisons.

Use `deterministic_reference_rng(seed)` when an established
`MersenneTwister` stream is part of a frozen reference. Use
`runtime_rng(seed)`, which returns `Xoshiro`, for new simulations and
benchmarks.

~~~julia
atmosphere_rng = runtime_rng(0x1234)
detector_rng = runtime_rng(0x5678)

epoch = advance_by!(atmosphere, 1e-3; rng=atmosphere_rng)
render_atmosphere!(pupil, renderer, atmosphere, epoch)
measure!(wfs, pupil, source, detector; rng=detector_rng)
~~~

An evolving RNG is persistent scientific state, not replaceable workspace
scratch. Give independently replayable atmosphere layers, detector channels,
and other stochastic owners separate streams. Deliberately sharing one stream
makes results depend on call order and must be documented by the application.

Rendering a published `AtmosphereEpoch` does not consume RNG. A zero-duration
advance after initialization returns the current epoch without changing the RNG
stream or epoch sequence. Direction renderers and direction batches may render
the same current epoch in any order before the atmosphere advances.

## Algorithm Graph Ownership

A stochastic native graph node owns one prepared RNG and its immutable reset
seed. For example, atmosphere and detector node configurations declare
`rng_seed`. Preparation creates the stream, `step_graph!` is its only writer,
and `reset_graph!` restores the configured seed and state. Host graphs use the
ordinary `Xoshiro` runtime stream.

CUDA graph owners instead use an internal array-only counter stream. Its seed
is a run-immutable plan value and its draw sequence is persistent state stored
on the selected CUDA device. Each array fill advances that sequence on the
prepared CUDA stream, then derives samples with stable SplitMix64 domains for
normal, uniform, and Poisson draws. Element and sub-draw indices are explicit,
so a recorded CUDA random-fill graph produces a new field on every replay.
Reset restores the device draw sequence without replacing its storage. This is
an execution detail of prepared graph owners, not a scalar `AbstractRNG` API
for scientific code.

Use distinct explicit seeds for independent nodes. Graph declaration order must
not be used as an implicit seed derivation rule. A graph file and its bound
configuration are therefore sufficient to identify the intended streams.

The current graph runtime is serial and single-writer. CUDA array fills are
addressed by a node-local draw sequence, random domain, element, and sub-draw.
If future execution replicates or distributes a stochastic owner, it will also
need a stable owner identity and explicit frame or event sequence. Sequential
host seed consumption is not a correct multi-device policy.

## Deterministic Configuration

- Use fixed seeds and record them with the run.
- Keep Julia, BLAS, and FFT-provider thread counts at one when strict replay
  matters.
- Import `AdaptiveOpticsSim.Ensembles` and use `DeterministicExecution()` for
  an ensemble. It rejects a multi-threaded Julia process and configures BLAS and
  the FFT provider for one thread.
- Disable detector noise for deterministic optical baselines unless noise is
  the subject of the test.
- Record configuration, graph file, command inputs, package revision, Julia
  version, backend versions, and compute device with validation artifacts.
- Use `StableRNGs.jl` only when a separately maintained cross-version bitstream
  is required.

## Phase Screens And Noise

Phase-screen creation and evolution must receive an explicit RNG. Persisting a
phase-screen sequence can be useful when comparing several algorithms against
identical turbulence.

Detector noise uses the RNG supplied to its acquisition or capture operation.
`Detector(noise=NoiseNone())` is the deterministic baseline. Stochastic
cross-backend tests should compare declared distributions or tolerances unless a
backend-specific bitwise contract has been established.

## Testing

- Run the same prepared operation twice from reset with the same seed and
  compare caller-visible products.
- Reverse independent operation order and verify that separately owned streams
  retain their products.
- Verify `reset_graph!` restores stochastic graph nodes.
- Compare reference data within an explicitly justified tolerance.
- Treat CPU/GPU parity as numerical or statistical evidence, not automatic
  bitwise identity.

## GPU Notes

GPU reductions and random providers may differ from CPU implementations.
Deterministic validation should default to CPU unless the specific accelerator
path has a documented deterministic contract. A GPU-ready algorithm must avoid
host scalar indexing and hidden host/device RNG transfers.

The prepared CUDA counter stream keeps seed advancement and sample generation
on the device and is compatible with capture of its asynchronous random-fill
operations. This does not claim that an arbitrary complete AOS graph is CUDA
Graph-capturable; every other operation in that graph must independently avoid
capture-time allocation and synchronization.

AMDGPU graph owners currently retain `Xoshiro`. In particular, the qualified
AMDGPU detector path mirrors Poisson sampling through preallocated host
storage because the maintained AMDGPU compiler/runtime does not yet admit the
exact per-pixel kernel. This is deterministic and tested, but it is not a
device-only stochastic path.

# Stacked Multi-Source and Multi-WFS Execution Plan

## Goal

Improve multi-source and multi-WFS runtime execution without changing the
public modeling API.

The public surfaces should remain:

- `Asterism` for multi-source sensing
- `SimulationInterface` for a single runtime I/O surface
- `CompositeSimulationInterface` for aggregated multi-runtime I/O

The intended change is internal execution strategy:

- replace serial per-source work with stacked execution where shapes align
- preserve serial fallback for irregular or incompatible cases
- keep CPU and GPU numerical behavior within existing parity tolerances

## Design Principles

1. Preserve physics and user-facing semantics.
   OOPAO is a behavioral reference, not an execution-model reference.
2. Prefer stacked arrays and batched kernels over lists of small independent
   launches.
3. Keep fallback paths explicit.
   Not every WFS family, source layout, or branch combination should be forced
   into a stacked path.
4. Measure every phase.
   Keep changes only when they improve end-to-end runtime or materially reduce
   launch/host-transfer overhead.

## Non-Goals

- no new public AO model abstraction before the existing surfaces are proven
  insufficient
- no blanket replacement of heterogeneous multi-WFS execution with one giant
  fused path
- no loss of deterministic/reference behavior to chase GPU throughput

## Current Baseline

Current runtime surfaces already in place:

- `Asterism` support in SH, Pyramid, BioEdge, and PSF generation
- `SimulationInterface` and `SimulationReadout`
- `CompositeSimulationInterface` for aggregated multi-runtime stepping
- stacked SH spot generation and batched detector capture for single-WFS pixel
  output

Current internal limitation:

- multi-source and multi-WFS paths still rely heavily on serial loops
- diffractive SH asterism on accelerator backends still performs per-source and
  per-subaperture work in a scalar orchestration pattern
- Pyramid/BioEdge asterism accumulation is still per-source serial
- `CompositeSimulationInterface` aggregates outputs but does not yet provide a
  stacked sensing executor for compatible branches

## Reference Behavior

Relevant OOPAO references:

- `tutorials/how_to_multi_sources.py`
- `tutorials/how_to_asterism.py`
- `OOPAO/Telescope.py`
- `OOPAO/ShackHartmann.py`
- `OOPAO/Pyramid.py`
- `OOPAO/BioEdge.py`

The target is to preserve the same modeling semantics while using a more
Julia- and backend-friendly execution model.

## Phase 1: Stacked Diffractive Shack-Hartmann Asterism

### Scope

- diffractive SH only
- start with same-wavelength `Asterism`
- support both detector-free and detector-coupled measurement

### Implementation Shape

- add source-stack workspaces for compatible asterism source counts
- build per-source spot/intensity stacks using the existing stacked SH
  machinery
- reduce across the source dimension once per measurement
- preserve current serial fallback when source shapes or detector assumptions
  are not stack-compatible

### Validation

- serial vs stacked SH slopes for 1, 2, and 4 sources
- detector-coupled pixel parity for deterministic/noise-free cases
- CPU vs AMDGPU vs CUDA parity on maintained smoke surfaces

### Benchmark Gates

- compare serial vs stacked on CPU for 2/4-source SH asterism cases
- compare serial vs stacked on AMDGPU and CUDA
- record whether speedup comes from:
  - fewer launches
  - fewer host transfers
  - lower orchestration cost

## Phase 2: Stacked Pyramid and BioEdge Asterism

### Scope

- diffractive Pyramid/BioEdge asterism paths
- compatible source shapes only

### Implementation Shape

- add stacked source-intensity accumulation where source geometry and image
  sizes match
- reduce once across the source dimension
- keep serial fallback for incompatible layouts, modulation paths, or future
  special cases

### Validation

- serial vs stacked parity for detector-free and detector-coupled paths
- parity of valid-pixel signal/slopes
- CPU/GPU parity under deterministic settings

### Benchmark Gates

- compare serial vs stacked per-frame sense time
- keep only if there is a real end-to-end runtime win or a meaningful reduction
  in backend synchronization/orchestration cost

## Phase 3: Compatible Multi-WFS Execution under `CompositeSimulationInterface`

### Scope

- multi-branch runtime execution for compatible WFS families and image shapes
- keep heterogeneous branch combinations on serial fallback

### Implementation Shape

- group child interfaces by compatible sensing family and stackable image
  dimensions
- add an internal stacked sensing executor for those groups
- preserve existing aggregation semantics:
  - concatenated command vector
  - concatenated slopes
  - per-branch WFS/science frame outputs

### Validation

- aggregated interface parity vs serial stepping
- mixed compatible/incompatible branch cases
- command/slopes/frame equivalence

### Benchmark Gates

- multi-WFS HIL timing on CPU, AMDGPU, and CUDA
- explicit comparison of:
  - serial child stepping
  - stacked compatible execution
  - fallback mixed execution

## Required Instrumentation

Add or extend runtime profilers for:

- SH asterism serial vs stacked timing
- Pyramid/BioEdge asterism serial vs stacked timing
- `CompositeSimulationInterface` multi-branch timing
- backend synchronization counts where meaningful

Prefer maintained scripts over ad hoc notebook-only profiling.

## Acceptance Criteria

The plan is complete when:

1. stacked SH asterism is implemented and benchmarked
2. stacked Pyramid/BioEdge asterism is implemented and benchmarked
3. compatible multi-WFS execution exists under `CompositeSimulationInterface`
4. all new paths preserve existing deterministic/reference expectations
5. the docs state clearly:
   - what stacks
   - what falls back to serial
   - where CPU vs GPU currently wins

## Stop Conditions

Stop pushing a phase if:

- parity becomes fragile relative to the maintained deterministic tests
- stacking only improves sub-phase microbenchmarks without improving runtime
  wall time
- the implementation cost grows faster than the measured benefit for the
  maintained scenarios

## Recommended Execution Order

1. Phase 1: SH asterism
2. Phase 1 validation and benchmark scripts
3. Phase 2: Pyramid/BioEdge asterism
4. Phase 2 validation and benchmark scripts
5. Phase 3: compatible multi-WFS execution
6. Phase 3 validation, HIL timing, and documentation update

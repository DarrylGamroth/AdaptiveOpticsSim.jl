# Next Sprint Plan

This sprint focuses on execution quality rather than feature expansion.

The immediate goal is a persistent runtime path that can execute a full closed
loop on CPU with zero steady-state allocations after warmup. That runtime will
also become the interface target for the next GPU-resident and RTC-in-the-loop
work.

## Sprint Goals

1. Add a persistent closed-loop runtime with an explicit `step!` API.
2. Make the representative CPU closed-loop step allocation-free after warmup.
3. Add measurements and regression tests for the runtime path.
4. Leave the runtime structured so the same API can later stay fully on device.

## Why This Sprint First

- The package already has strong feature parity and a validated CUDA smoke
  matrix.
- The next engineering bottleneck is predictable execution, not missing
  algorithms.
- A single long-lived runtime object is the cleanest base for:
  - CPU latency/jitter benchmarking
  - GPU-resident execution
  - RTC adapter work

## Deliverables

### 1. Runtime API

Add a small persistent runtime layer that owns the long-lived simulation state:

- `ClosedLoopRuntime`
- `sense!`
- `step!`
- `set_command!`
- optional science-frame capture through a preconstructed detector

The runtime should wrap existing package objects rather than duplicate logic.

### 2. Zero-Allocation CPU Baseline

Define one normative steady-state path and make it allocation-free after warmup:

- atmosphere advance
- atmosphere propagation
- DM application
- WFS measurement
- reconstructor application
- DM command update

This should be measured inside functions with persistent state.

### 3. Benchmarks and Tests

Add:

- a benchmark entry for the persistent runtime step
- an allocation regression test for the CPU runtime step
- optional science-detector coverage if it can stay allocation-free in the
  chosen representative case

## Acceptance Criteria

The sprint is successful when all of the following are true:

1. A persistent closed-loop runtime exists in the library, not only in examples.
2. A warmed CPU `step!` call for the representative closed-loop path allocates
   `0` bytes.
3. The full CPU test suite passes.
4. The runtime API is backend-parametric and does not bake in CPU-only array
   types.

## Current Status

This sprint has now hit its original targets:

1. `ClosedLoopRuntime`, `sense!`, `step!`, and `set_command!` are implemented.
2. The representative warmed CPU runtime path is allocation-free.
3. The full CPU test suite passes with runtime allocation regression coverage.
4. The same runtime API has been validated on CUDA with
   `CUDA.allowscalar(false)` through the package GPU smoke matrix.

## Out of Scope

Not targeted in this sprint:

- hard real-time adapters or network transport
- broad GPU optimization passes beyond preserving backend-generic runtime design
- numerical-stability hardening from `docs/numerical-stability-review.md`
  beyond issues directly exposed by the new runtime path
- large new `KernelAbstractions` coverage unrelated to the runtime target

## Follow-On Sprint

Once this sprint lands, the next recommended sequence is:

1. validate the same runtime API in a device-resident GPU loop
2. audit host-device synchronization boundaries
3. return to the numerical-stability backlog with runtime benchmarks in place

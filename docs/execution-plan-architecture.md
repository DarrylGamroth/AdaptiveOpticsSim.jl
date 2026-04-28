# Execution Plan Architecture

Status: active

Plan traceability:

- motivated by the grouped-runtime follow-up after
  [grouped-runtime-plan.md](./grouped-runtime-plan.md)
- aligned with [future-platform-direction.md](./future-platform-direction.md)
- supports the benchmark and validation direction in
  [post-review-platform-plan.md](./post-review-platform-plan.md)

## Purpose

This document defines the internal architecture for backend-specific execution
plans.

The package already has several places where:

- the semantic contract is shared,
- the numerical result should remain equivalent,
- but the best memory-access pattern differs by backend or workload shape.

The grouped AMDGPU regression after the grouped-runtime milestone is the
clearest recent example. The grouped runtime contract was improved, but one
shared grouped accumulation implementation was not equally efficient on CPU,
AMDGPU, and CUDA.

The goal of this architecture is to formalize a better rule:

- keep one semantic model,
- allow multiple execution plans,
- choose the plan by dispatch and benchmarking,
- validate semantic equivalence across plans.

## Non-Goals

This architecture is not intended to:

- create fully separate algorithm implementations for each backend,
- introduce a large trait hierarchy before it is needed,
- duplicate scientific logic across CPU, AMDGPU, CUDA, and future Metal paths,
- let performance-specialized paths drift without shared validation.

The package should prefer shared contracts with specialized execution policy,
not unrelated per-backend feature trees.

## Problem Statement

The current package contains several classes of work where performance depends
heavily on layout and staging strategy:

- grouped WFS accumulation,
- detector capture and readout staging,
- atmosphere extraction and atmospheric field propagation,
- reduction-heavy normalization/statistics paths,
- runtime product export,
- tomography and calibration builders.

In these areas, the following tension appears repeatedly:

- the mathematically correct output is the same,
- but direct rendering into a view, staging through a stable 2D buffer,
  reducing from a 3D stack, or fusing render-and-accumulate can have very
  different performance across backends.

The package needs a first-class way to express that difference without
rewriting the model itself for each target.

## Design Principles

1. Algorithms define semantics.
2. Execution plans define memory layout, staging, and kernel scheduling.
3. Dispatch chooses plans; hot paths should not branch on backend names in
   inner loops.
4. Benchmarks decide the default plan for maintained targets.
5. Tests enforce semantic equivalence across plans.
6. Plans remain internal until a user-facing need appears.

## Static Dispatch Requirement

Execution-plan dispatch should be static and type-driven.

That means:

- plan selection should resolve from types, not dynamic string/symbol/backend
  comparisons inside hot paths,
- the chosen plan object should have a concrete type at the call site,
- execution should pass through methods that the compiler can specialize and
  inline,
- plan selection should happen before entering the inner loop whenever
  possible.

This matches the broader package rule as well:

- prefer types and multiple dispatch,
- avoid dynamic dispatch in hot paths,
- avoid runtime backend-name branching where a typed backend or typed plan can
  express the same choice,
- let the JIT specialize the concrete operation pipeline for the active
  backend/workload combination.

Execution plans are meant to reinforce that package-wide design, not to add a
new dynamic policy layer.

## Terminology

- semantic contract:
  - the mathematically intended operation, outputs, and ownership rules
- execution plan:
  - the backend-specific strategy used to realize that contract
- operation family:
  - a class of work such as grouped WFS accumulation or detector readout
- backend tag:
  - the maintained backend identity already used in the package, for example
    CPU, AMDGPU, CUDA, and future Metal
- workload shape:
  - dimensions or mode choices that materially affect the best execution plan,
    for example grouped vs single-source, compact vs medium, or detector-backed
    vs direct-signal

## Architectural Rule

For performance-sensitive subsystems, the package should be structured as:

1. semantic contract
2. internal execution-plan selection
3. backend-specific plan implementation
4. equivalence and benchmark evidence

The contract owns:

- input/output meaning,
- product ownership,
- accepted invariants,
- error conditions.

The plan owns:

- scratch layout,
- staging buffer strategy,
- reduction strategy,
- kernel scheduling and synchronization,
- host/device transfer policy where unavoidable.

## Core Abstractions

The architecture should remain lightweight and dispatch-oriented.

Recommended pattern:

```julia
abstract type AbstractExecutionPlan end
abstract type AbstractOperationFamily end

struct GroupedWFSFamily <: AbstractOperationFamily end
struct DetectorReadoutFamily <: AbstractOperationFamily end
struct AtmosphericFieldFamily <: AbstractOperationFamily end
struct RuntimeExportFamily <: AbstractOperationFamily end
```

Plan objects remain small and internal:

```julia
struct DirectAccumulatePlan <: AbstractExecutionPlan end
struct Staged2DPlan <: AbstractExecutionPlan end
struct StackReducePlan <: AbstractExecutionPlan end
struct FusedAccumulatePlan <: AbstractExecutionPlan end
```

Selection happens through dispatch:

```julia
execution_plan(::Type{BackendTag}, ::GroupedWFSFamily, wfs, workload) = ...
execute!(plan::DirectAccumulatePlan, args...) = ...
execute!(plan::StackReducePlan, args...) = ...
```

The key point is that plan selection should happen outside the hot inner loop.
In practice, the call site should see a concrete `plan` type so the compiler
can optimize the chosen execution path aggressively.

## Where Plans Live

Plans should live close to the subsystem that owns the contract, not in one
giant global registry.

Recommended package layout:

- grouped WFS plans near [grouped.jl](../src/wfs/grouped.jl)
- detector-readout plans near the detector pipeline layer
- atmospheric-field plans near
  [atmospheric_field_propagation.jl](../src/optics/atmospheric_field_propagation.jl)
- reduction plans near [reductions.jl](../src/core/reductions.jl)
- runtime export plans near the runtime product layer

Core helper code may exist for shared dispatch patterns, but subsystem-specific
plans should remain owned by their subsystem.

## Plan Selection Inputs

Execution-plan selection should be allowed to depend on:

- backend tag,
- operation family,
- sensor or model family,
- grouped vs single-source mode,
- detector-backed vs direct path,
- problem size when a clear crossover exists,
- explicit benchmark-backed policy overrides where needed.

Selection should not depend on:

- opaque runtime heuristics that are hard to reproduce,
- user-visible mutable globals,
- ad hoc `isa` checks in hot paths,
- symbol/string backend dispatch inside inner loops.

## Recommended First-Class Operation Families

### Grouped WFS Accumulation

This is the immediate target.

The grouped WFS contract is already shared, but the grouped accumulation plan
should be backend-specific.

Likely maintained plan variants:

- `DirectAccumulatePlan`
  - render each source and add directly into the final output
  - likely best on CPU for simple loops
- `Staged2DPlan`
  - render into one stable 2D buffer, then copy or accumulate
  - likely useful for AMDGPU Pyramid/BioEdge
- `StackReducePlan`
  - render into a 3D stack and reduce afterward
  - currently works well enough on CUDA and is easy to reason about
- `FusedAccumulatePlan`
  - render-and-accumulate without a full stack materialization
  - candidate future optimization if benchmarks justify it

### Detector Readout

The detector pipeline already distinguishes sampled stack, detector signal, and
exported products. Execution plans should decide:

- whether staging is device-only,
- whether host mirrors are required for specific backends,
- whether readout correction and noise application are fused or staged.

### Atmospheric Field Propagation

The semantic contract is already shared. Plans can vary:

- direct field updates,
- staged intensity extraction,
- per-layer fencing vs phase scheduling,
- spectral accumulation policy.

### Reduction And Statistics

Reduction-heavy paths should keep one semantic API and allow backend-specific
plans for:

- masked sum,
- grouped reduction,
- centroid statistics,
- normalization statistics.

### Runtime Product Export

Runtime product semantics should remain stable, but plan policy can decide:

- device-only product,
- device-owned plus on-demand host mirror,
- stacked grouped export,
- per-branch export only,
- no export when not required.

### Tomography And Calibration Builders

Builder semantics should remain shared, while plan policy chooses:

- assembly staging,
- reduction strategy,
- scratch ownership,
- host/device synchronization boundaries.

## Example: Grouped WFS Design

The grouped WFS contract should look like:

```julia
grouped_accumulation_plan(::Type{BackendTag}, wfs, workload)::AbstractExecutionPlan
accumulate_grouped_sources!(plan::AbstractExecutionPlan, out, scratch, sources, render!, args...)
```

Example maintained defaults:

- CPU + grouped Pyramid/BioEdge:
  - `DirectAccumulatePlan` or `Staged2DPlan`
- AMDGPU + grouped Pyramid/BioEdge:
  - `Staged2DPlan`
- CUDA + grouped Pyramid/BioEdge:
  - `StackReducePlan`
- SH:
  - may use a different default from Pyramid/BioEdge because its cost balance is
    different

This keeps the grouped runtime contract stable while allowing different memory
patterns per target.

## Validation Requirements

Every plan family should have:

1. semantic equivalence tests
2. backend parity tests where maintained
3. benchmark evidence for default-plan choice

At minimum, each new default plan should be backed by:

- one functional test proving output equivalence to the reference contract
- one maintained benchmark surface showing why the chosen default is preferred
- one limitation note if a backend-specific fallback is not ideal but is kept
  for correctness

## Benchmark Policy

Default plan selection must be justified by maintained benchmark surfaces, not
just by local intuition.

For grouped WFS this means:

- CPU benchmark surface
- AMDGPU benchmark surface
- CUDA benchmark surface
- future Metal benchmark surface when maintained

When a plan differs by backend, the benchmark record should state:

- which plan is the default,
- what alternatives were considered,
- what performance or allocation tradeoff justified the choice.

## Apple-Silicon Guidance

Future Apple support should follow the same architecture.

Expected target split:

- Apple CPU path:
  - standard scalar/array CPU plans, with optimized linear-algebra backends
    where appropriate
- Apple GPU path:
  - Metal-backed accelerator plans

This architecture is specifically designed so that future Metal support can
choose its own execution plans without forcing either:

- CUDA-style assumptions onto Metal, or
- AMDGPU-specific staging compromises onto every backend.

The package should not assume that one GPU plan is portable across CUDA,
AMDGPU, and Metal.

## Rollout Order

The recommended rollout order is:

1. grouped WFS accumulation plans
2. detector readout plans
3. atmospheric field propagation plans
4. reduction/statistics plans where still duplicated
5. runtime export plans
6. tomography/calibration builder plans

This order follows the places where the package already has:

- active benchmark surfaces,
- recent regressions or backend divergence,
- a strong need for clearer staging ownership.

## Immediate Recovery Work

The immediate use of this architecture should be the grouped AMDGPU recovery.

Suggested first pass:

1. formalize grouped WFS plan selection
2. restore a staged 2D accumulation plan for AMDGPU Pyramid/BioEdge
3. benchmark CPU, AMDGPU, and CUDA against the current grouped baseline
4. keep the shared grouped runtime contract unchanged

Acceptance for that recovery:

- no grouped-runtime correctness regression,
- recovered AMDGPU Pyramid/BioEdge medium benchmark surface,
- no CUDA regression relative to the current grouped baseline.

## Guardrails

- Do not move algorithm semantics into backend extensions.
- Do not let plans change mathematical meaning.
- Do not add public API surface until internal patterns stabilize.
- Do not keep a backend-specific default without benchmark evidence.
- Do not duplicate full subsystem files when a plan object plus dispatch is
  sufficient.

## Definition Of Done

This architecture is working when:

- performance-sensitive subsystems expose shared semantic contracts,
- execution plans are explicit and benchmark-backed,
- backend-specific plan choice is localized and testable,
- semantic equivalence is enforced across plans,
- adding a new backend does not require rewriting the algorithm surface.

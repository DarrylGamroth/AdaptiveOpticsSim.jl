# Control Simulation Architecture Plan

This plan defines the next architecture pass for `AdaptiveOpticsSim.jl` as the
package grows beyond a straight OOPAO port. The goal is to extract the generic
runtime/control pieces that are currently embedded in the Subaru AO188
simulation support module, while keeping scenario assembly code out of the core
package.

## Motivation

The current Subaru AO188 simulation in
`examples/support/subaru_ao188_simulation.jl` contains two kinds of code:

- reusable runtime/control primitives that should live in the core package
- scenario-specific assembly and defaults that should remain example-level

The package is also starting to need explicit interfaces for:

- live control simulations
- execution/scheduling policies
- delay models
- reconstructor operators
- prepared replay/runtime-precompute behavior

This is a better fit for Julia than continuing to grow scenario-specific code
or OOPAO-like monolithic objects in core.

## Goals

- keep scenario assembly such as Subaru AO188 out of core
- move reusable runtime/control building blocks into `src/`
- make control simulations composable through explicit interfaces
- improve extensibility for future MIMO, multi-rate, and multi-WFS systems
- reduce duplication between examples and core runtime paths

## Non-Goals

- do not move telescope-specific, instrument-specific, or workflow-specific
  defaults into core
- do not redesign the scientific optics stack in this pass
- do not block future `ZernikeWFS` work on a full runtime redesign

## Core Extractions

These are the AO188-derived pieces that should become generic core features.

### 1. Reconstructor Operators

Current AO188 type:

- completed: `MappedReconstructor` now lives in core

Move into core as a generic operator concept:

- `MappedReconstructor`
  - modal/operator stage: `R`
  - command-basis stage: `B`
  - workspace: preallocated modal buffer
  - optional gain/scaling metadata

Target API:

- `reconstruct!(out, recon, slopes)`
- `reconstruct(recon, slopes)`

Rationale:

- this pattern is not AO188-specific
- it already improved runtime and numerical behavior
- future tomography/runtime systems will want the same two-stage form

### 2. Delay Models

Current AO188 type:

- `VectorDelayLine`

Move into core as:

- `VectorDelayLine`
- optional wrapper:
  - `DelayModel`
  - or a small family of typed delay containers if that proves clearer

Target API:

- `shift_delay!(line, sample)`
- construction from a reference vector and frame delay

Rationale:

- explicit control delays are generic runtime behavior
- future multi-rate and controller integration work will need this

### 3. Execution Policies

Current AO188 status:

- completed: generic execution policies now live in core as
  `SequentialExecution`, `ThreadedExecution`, and `BackendStreamExecution`

Move into core as:

- `AbstractExecutionPolicy`
- concrete policies:
  - `SequentialExecution`
  - `ThreadedExecution`
  - `BackendStreamExecution`

Target API:

- `init_execution_state(policy, ref)`
- `execute!(policy, state, work...)`

Rationale:

- these policies are generic scheduling choices, not AO188 concepts
- future multi-WFS and multi-rate systems should share them

### 4. Prepared Runtime Hooks

Current AO188 function:

- `prepare_replay!`

Move into core as a generic preparation interface:

- `prepare!(sim)` or `prepare_runtime!(sim)`

Target behavior:

- perform expensive precomputes needed for repeated stepping
- mark prepared state explicitly
- no-op for simulations that do not need preparation

Rationale:

- Shack-Hartmann replay prep is not unique to AO188
- other WFS/runtime cases may need the same lifecycle

## New Core Interfaces

### 5. Control Simulation Interface

Add an abstract simulation-level interface for runtime systems.

Candidate names:

- `AbstractControlSimulation`
- `AbstractSimulationRuntime`

Recommended minimal API:

- `step!(sim)`
- `simulation_interface(sim)`
- `prepare!(sim)` or `prepare_runtime!(sim)`

Optional later additions:

- `reset!(sim)`
- `runtime_timing(sim; ...)`
- `execution_policy(sim)`

Rationale:

- unifies AO188-style systems, tutorial closed-loop systems, and future
  composite runtime assemblies
- keeps the public model centered on simulation/runtime behavior rather than
  RTC-specific wording

### 6. Runtime Capability Traits

Introduce small capability traits rather than ad hoc branching.

Candidate traits:

- `supports_prepared_runtime(::Type)`
- `supports_detector_output(::Type)`
- `supports_stacked_sources(::Type)`
- `supports_grouped_execution(::Type)`
- `supports_multirate_schedule(::Type)`

Rationale:

- clearer dispatch
- less scenario-specific branching in runtime scripts/examples
- easier to document and test

### 7. Delay / Scheduling / Rate Traits

Add explicit scheduling concepts for future systems.

Candidate abstractions:

- `AbstractDelayModel`
- `AbstractExecutionPolicy`
- `AbstractRateSchedule`

This can start small:

- fixed frame-delay vectors first
- leave more complex multi-rate scheduling for a later phase

## What Stays Example-Level

The following should remain outside core:

- `AO188LatencyModel`
- `AO188WFSDetectorConfig`
- `AO188SimulationParams`
- `AO188Simulation`
- `CircularActuatorSupport`
- `default_ao188_low_order_resolution`
- AO188-specific actuator support logic
- AO188-specific high/low-order branch assembly
- AO188-specific CUDA branch-stream behavior unless it becomes generic enough
  to absorb under the core execution-policy interface

These are scenario-building choices, not library primitives.

## Proposed Phases

### Phase 1: Core Operator and Delay Extraction

Status: implemented

- move `FullCommandReconstructor` into core as `MappedReconstructor`
- move `VectorDelayLine` into core
- update AO188 simulation to use the core versions
- add direct tests for both in core

Exit criteria:

- AO188 simulation behavior unchanged
- examples no longer define these types locally

### Phase 2: Core Execution Policy Extraction

Status: implemented

- define `AbstractExecutionPolicy`
- add `SequentialExecution`, `ThreadedExecution`, `BackendStreamExecution`
- move generic initialization/execution hooks into core
- update AO188 simulation to use the core policies

Exit criteria:

- no AO188-specific execution-policy types remain
- CUDA/AMDGPU behavior preserved

### Phase 3: Control Simulation Interface

Status: implemented for the current runtime surface

- add `AbstractControlSimulation`
- add `prepare!` / `simulation_interface` conventions
- make AO188 simulation conform explicitly
- make `ClosedLoopRuntime` conform explicitly
- update maintained examples/scripts to target the interface where practical

Exit criteria:

- maintained scripts use the generic interface where practical

### Phase 4: Capability Traits

Status: implemented for the first runtime capability set

- add the first small set of runtime capability traits
- use them to simplify example/runtime branching
- document them in API/docs

Exit criteria:

- less ad hoc conditional logic in example/runtime integration code

### Phase 5: Broader Adoption

Status: in progress

- apply `MappedReconstructor` where other runtime systems currently use the same
  pattern
- continue moving tutorial/runtime examples onto the generic
  control-simulation interface where it improves clarity
- defer multi-rate scheduling unless a real use case requires it

## Validation Requirements

Each phase must preserve:

- `Pkg.test()`
- OOPAO reference regression
- maintained GPU smoke
- AO188 CPU audit
- AO188 GPU audit surfaces where relevant

Additional targeted checks:

- no new hot-path allocations from the extracted core abstractions
- no performance regression in AO188 runtime after moving the reusable pieces
- no public confusion about what is core vs example-level

## Recommended Order

1. `MappedReconstructor`
2. `VectorDelayLine`
3. `ExecutionPolicy`
4. `AbstractControlSimulation`
5. capability traits

This order gives the best ratio of cleanup to risk.

## Relationship to ZernikeWFS

This architecture pass should not block `ZernikeWFS`.

Recommended sequencing:

- finish this plan’s Phase 1 first if the reconstructor/delay extractions are
  easy and clearly beneficial
- otherwise move to `ZernikeWFS` and resume this plan afterward

The important boundary is:

- generic reusable runtime/control primitives belong in core
- complete instrument/system assemblies belong in examples or higher-level
  packages

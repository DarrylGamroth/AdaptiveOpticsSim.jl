# Platform Orchestration

Status: active

## Purpose

This guide defines the maintained Julia-native orchestration layer for
platform-scale runtime composition.

The goal is to make realistic runtime assembly more explicit and reusable
without switching the package to a config-file-first model.

## Design

The orchestration layer is built from three pieces:

- `RuntimeBranch`
  - one typed branch definition containing:
    - a simulation
    - a reconstructor
    - optional WFS and science detectors
    - an RNG
- `SingleRuntimeConfig` or `GroupedRuntimeConfig`
  - typed runtime-composition policy:
    - scenario name
    - branch labels
    - runtime profile / latency
    - exported product policy
- `RuntimeScenario`
  - the executable orchestration object that wraps the resulting
    `SimulationInterface` or `CompositeSimulationInterface`

This keeps the public workflow:

- script-first
- type-driven
- builder-oriented

instead of replacing Julia composition with external manifests.

## Core API

Single-branch runtime:

```julia
branch = RuntimeBranch(
    :main,
    sim,
    recon;
    wfs_detector=wfs_det,
    science_detector=science_det,
    rng=MersenneTwister(1),
)

cfg = SingleRuntimeConfig(
    name=:single_branch_demo,
    branch_label=:main,
    products=RuntimeProductRequirements(slopes=true, wfs_pixels=true, science_pixels=true),
)

scenario = build_runtime_scenario(cfg, branch)
prepare!(scenario)
step!(scenario)
readout = readout(scenario)
```

Grouped runtime:

```julia
cfg = GroupedRuntimeConfig(
    (:high, :low);
    name=:grouped_demo,
    products=GroupedRuntimeProductRequirements(
        wfs_frames=true,
        science_frames=false,
        wfs_stack=true,
        science_stack=false,
    ),
)

scenario = build_runtime_scenario(
    cfg,
    RuntimeBranch(:high, sim1, recon1; wfs_detector=det1, rng=MersenneTwister(1)),
    RuntimeBranch(:low, sim2, recon2; wfs_detector=det2, rng=MersenneTwister(2)),
)

prepare!(scenario)
step!(scenario)
stack = grouped_wfs_stack(scenario)
```

## Execution Semantics

`RuntimeScenario` is a thin orchestration wrapper.

It delegates to the maintained runtime surfaces:

- `prepare!(scenario)`
- `sense!(scenario)`
- `step!(scenario)`
- `set_command!(scenario, cmd)`
- `readout(scenario)`
- `runtime_timing(scenario)`
- `runtime_phase_timing(scenario)`

So the scenario layer does not introduce a second execution model. It makes the
composition contract explicit and typed.

## Product Policy

Single-branch scenarios use the `RuntimeProductRequirements` recorded in
`SingleRuntimeConfig`, unless a branch explicitly overrides products.

Grouped scenarios use `GroupedRuntimeProductRequirements` for the composite
export surface. Branch-local runtime product policy is derived automatically
from the grouped export request unless a branch explicitly overrides products.

That means a grouped scenario requesting `wfs_stack=true` will automatically
build branch runtimes with `wfs_pixels=true` when a WFS detector is present.

## Scope Boundary

This layer is intentionally about runtime orchestration, not full physical
system authoring.

It should own:

- branch composition
- runtime/export policy
- execution entry points
- reusable platform scenario structure

It should not own:

- instrument-specific physical defaults
- external config parsing
- bespoke benchmark normalization rules

Those remain in support modules, scripts, and benchmark contracts.

## Maintained Examples

Canonical examples using this layer live in:

- [platform_single_runtime.jl](../examples/closed_loop/platform_single_runtime.jl)
- [platform_grouped_runtime.jl](../examples/closed_loop/platform_grouped_runtime.jl)

The maintained grouped runtime benchmark also uses this layer:

- [profile_multi_source_multi_wfs_runtime.jl](../scripts/profile_multi_source_multi_wfs_runtime.jl)

Direct runtime evidence and backend validation for the layer are recorded in:

- [platform-orchestration-validation.md](./platform-orchestration-validation.md)

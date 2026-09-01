# Roadmap

Status: active

## Current Direction

AdaptiveOpticsSim.jl owns adaptive-optics physics and simulation algorithms.
It has two composition surfaces:

- direct Julia for generated topology, multiple rates, conditional execution,
  sub-frame sampling, and research workflows
- `AdaptiveOpticsSim.AlgorithmGraphs` for static, single-rate, complete-frame
  graphs with prepared storage and deterministic execution

`AlgorithmGraphs` is deliberately a small scheduler. It does not own wall-clock
pacing, transport, PipeWire integration, or a general event calendar. The
transport-neutral `PreparedGraphHILBoundary` exchanges one complete detector
frame for one complete external-RTC command in lockstep.

The removed `AdaptiveOpticsSim.Plant` event runtime is preserved on
`archive/plant-runtime-2026-08-30`. Its unfinished concrete RNG-owner work is
preserved on `archive/plant-pe06b2-wip-2026-08-30`. Neither branch is part of
the maintained package.

## Implemented Foundations

- Explicit params, run-immutable plans, persistent state, replaceable
  workspaces, caller-visible products, and prepared execution owners.
- CPU algorithms plus optional CUDA and AMDGPU execution for covered kernels.
- Static graph definitions in Julia or TOML, direct and one-frame delayed links,
  model-time drivers, preparation-time validation, and reset.
- Native graph nodes for atmosphere OPD, deformable-mirror surfaces, pupil OPD
  composition, Shack-Hartmann and Pyramid optical rates, detector acquisition,
  centroid/slope selection, reconstruction, and controller operations.
- Lockstep host command/frame exchange for external RTC integration.
- A downstream-package seam for instrument-owned graph files, geometry,
  calibration policy, pyRTC tests, and model-validity claims.
- Optional Proper.jl integration at an explicit application or graph-node
  boundary.

## Near-Term Priorities

1. Validate the reduced scheduler and graph-file contract as the stable
   complete-frame orchestration surface.
2. Keep the generic graph, HIL, and pyRTC contracts stable for the separately
   maintained REVOLT Classic and Copper packages; instrument calibration and
   performance evidence remain with those packages.
3. Build a separate SPIDERS instrument package from controlled Subaru/AO3k,
   SpiderMan, service, and optical-engineering inputs. Mark estimated chopper,
   stage, filter-wheel, and prescription values as provisional.
4. Add a PipeWireAO output adapter after the simulated detector boundary. GPU
   simulations may copy a completed frame to the CPU RTC boundary explicitly;
   hidden mixed-device execution is out of scope.
5. Extend graph-node coverage only when an actual architecture needs it. Keep
   large calibration arrays as prepared parameters, not scalar properties.
6. Maintain deterministic explicit RNG ownership and OOPAO or instrument
   reference comparisons for scientific changes.
7. Improve rolling-shutter models by evaluating the optical path at row or
   row-group times while publishing one atomic completed frame.
8. Close remaining production-surface gaps with local CPU, AMDGPU, WSL CUDA,
   and aarch64 evidence as applicable.

## Graph Scope

TOML is the user-facing format for static graphs because it is readable,
reviewable, and already supported by Julia without an additional parser. Julia
remains the authority for compositions that cannot be expressed cleanly in the
versioned file schema.

A graph node must have one writer for its persistent state and outputs. This
means exactly one prepared execution owner mutates those values. It does not
prevent MCAO, MOAO, multiple guide stars, multiple DMs, branching optical
paths, or coarse parallel execution across independent prepared owners.

The version 1 graph-file contract remains intentionally limited to static,
single-rate, complete-frame execution. Multi-rate and mid-frame device changes
belong in explicit Julia composition until a concrete instrument requirement
justifies a small, testable scheduler extension.

## Deferred Work

- Direct loading of Julia algorithms inside a PipeWire process. Current planning
  keeps Julia embedding in a separate `module-julia`-style boundary.
- PipeWire GPU-buffer execution.
- Whole-step CUDA Graph or HIP Graph qualification for complete downstream
  instrument graphs. AOS qualifies reusable node families; each instrument
  package owns its composition-level evidence.
- A general multi-rate event runtime.
- Automatic graph partitioning, implicit host/device transfer, and dynamic
  topology mutation.
- Camera or instrument parameters that have not been confirmed by their
  scientific or optical owners.

## Documentation Rule

Keep maintained documents concise and current. Use issues and pull requests for
temporary plans and audits, and Git history or the named archive branches for
the retired Plant implementation.

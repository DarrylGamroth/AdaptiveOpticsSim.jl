# Benchmark Matrix Plan

Status:
- [~] Phase 2 partially implemented
  - `scripts/profile_zernike_runtime.jl` now supports `compact`, `medium`, and
    `representative`
  - `scripts/profile_pixel_output_runtime.jl` now supports `compact`,
    `medium`, and `representative`
  - `scripts/profile_multi_source_multi_wfs_runtime.jl` now supports
    `compact`, `medium`, and `representative`
  - `scripts/profile_revolt_hil_runtime.jl` now covers a REVOLT-style
    pixel-output SH HIL case with a fixed representative configuration
  - `scripts/profile_external_optics_hil.jl` now covers a representative
    external-optics HIL case where AdaptiveOpticsSim owns DM phase generation
    and phase export but not the downstream image formation
- [ ] representative closed-loop and builder ladders are still pending beyond
  these first runtime surfaces

## Goal

Separate the package's performance surfaces into two explicit classes:

- **Regression benchmarks**
  - small, deterministic, cheap to run
  - protect hot-path latency/allocation behavior and backend correctness
  - suitable for frequent local use and maintained audit workflows
- **Representative benchmarks**
  - medium and large configurations intended to approximate real AO workloads
  - answer CPU vs GPU crossover, scaling, and algorithm-appropriateness
    questions
  - suitable for engineering decisions, not just regression detection

The package already has strong regression-sized timing coverage. What is missing
is a consistent representative benchmark matrix that lets us say when a given
algorithm or backend becomes appropriate for realistic systems.

## Problem Statement

Many current maintained runtime profiles are intentionally compact:

- compact `ClosedLoopRuntime`
- compact AO188 simulation slices
- compact stacked multi-source / multi-WFS profiles
- compact `ZernikeWFS` runtime profile

Those are valuable because they are stable and fast, but they bias the observed
story toward:

- CPU looking much stronger than GPU
- per-step latency dominating throughput considerations
- tiny detector/FFT/gather workloads dominating conclusions

That is the right story for regression guards, but it is not enough to judge:

- whether a GPU-oriented implementation is appropriate at realistic scale
- whether a runtime algorithm scales well with pupil/detector size
- whether a sensor model is structurally suitable for HIL workloads
- where the real CPU/GPU crossover occurs for maintained workflows

## Benchmark Classes

### Class A: Regression Benchmarks

Purpose:

- catch performance regressions
- protect zero-allocation or low-allocation hot paths
- keep backend smoke/audit scripts cheap enough for routine use

Properties:

- deterministic
- warm-started
- small enough to run quickly on developer machines
- stable across hosts and backend versions

Examples already in place:

- `scripts/profile_zernike_runtime.jl`
- `scripts/profile_mixed_sh_asterism_runtime.jl`
- `scripts/profile_multi_source_multi_wfs_runtime.jl`
- `scripts/profile_pixel_output_runtime.jl`
- `scripts/ao188_3k_hil_audit.jl`
- `benchmarks/benchmark_{cpu,cuda,amdgpu}.jl`

### Class B: Representative Benchmarks

Purpose:

- approximate realistic optical/control workloads
- locate CPU/GPU crossover points
- judge whether an algorithm is appropriate for deployment-scale simulations
- provide a better basis for optimization priorities

Properties:

- use physically plausible medium/large configurations
- may be slower and more hardware-sensitive
- should be run intentionally, not on every quick validation pass
- should record both throughput and shape/context, not just a single time

These benchmarks should be treated as engineering evidence rather than
regression gates.

## Representative Benchmark Families

The representative matrix should cover the main runtime and builder families,
not every public function.

### 1. Pixel-Output WFS Runtime

Families:

- Shack-Hartmann
- Pyramid
- BioEdge
- ZernikeWFS

Questions:

- how do frame size and pupil size change CPU/GPU crossover?
- when does GPU stop being dominated by launch/sync overhead?
- which sensing families scale most cleanly on CUDA and AMDGPU?

Scale ladder:

- `compact`
  - current regression-sized surface
- `medium`
  - larger pupil and detector sizes where GPU may start to amortize overhead
  - intended as a scaling/crossover rung, not automatically a deployment-scale
    configuration
- `representative`
  - realistic single-WFS closed-loop size meant to approximate deployment use

### 2. Multi-Source / Multi-WFS Runtime

Families:

- stacked asterism sensing
- grouped `CompositeSimulationInterface` execution
- mixed heterogeneous branch execution

Questions:

- how much does batching improve backend utilization?
- when does multi-branch GPU execution become worthwhile?
- when does grouped execution beat branch-by-branch stepping?

Scale ladder:

- `compact`
  - current maintained stacked profiles
- `medium`
  - more sources and larger detector tiles
  - intended as a scaling/crossover rung, not automatically a final
    instrument-like scene
- `representative`
  - realistic multi-WFS HIL or tomography-facing runtime scene

### 3. Closed-Loop HIL Runtime

Families:

- compact maintained runtime
- Subaru AO188 simulation example
- fixed representative HIL cases with large pixel outputs
- later representative SCAO / woofer-tweeter / MCAO examples

Questions:

- what is the realistic per-frame rate for CPU vs GPU?
- how much does latency staging or multi-rate structure matter?
- when do runtime interface choices become the limiting factor?

Scale ladder:

- `compact`
  - existing minimal runtime surfaces
- `medium`
  - AO188-class example
  - larger and more informative than compact, but not necessarily the final
    representative endpoint for every system family
- `representative`
  - system-level examples with realistic detector sizes and control paths

For this family, some representative systems are better modeled as fixed named
cases instead of generic scale rungs. Two maintained examples now exist:

- `revolt_sh_hil`
  - REVOLT-style SH HIL proxy
  - `277` active commands mapped onto a `17 x 17` DM grid
  - `16 x 16` subapertures
  - `22 x 22` ROI
  - full `352 x 352` pixel output mosaic
- `external_optics_hil`
  - external-propagator / PROPER-style proxy
  - `468` active commands mapped onto a `24 x 24` DM grid
  - full DM OPD generation on a `640 x 640` phase grid
  - exported `640 x 512` phase crop for downstream image formation

### 4. Calibration / Builder Workloads

Families:

- interaction matrices
- modal basis construction
- LiFT
- model tomography

Questions:

- where does GPU already have clear value?
- which builders should default to CPU-build / GPU-run?
- which calibration paths remain too host-oriented?

Scale ladder:

- `compact`
  - current audit/builder surfaces
- `medium`
  - current tomography profiling surfaces
  - intended to expose scaling trends and crossover onset before the fully
    representative builder cases
- `representative`
  - KAPA-like or instrument-scale assembly workloads

## Required Reporting

Every representative benchmark entry should report:

- backend
- configuration class (`compact`, `medium`, `representative`)
- key problem dimensions
  - pupil resolution
  - detector frame shape
  - number of sources
  - number of WFS
  - command length / slope length where relevant
- warmed mean and p95 timing
- effective frame rate or throughput
- phase timing when available
- allocation or memory notes where relevant

Representative reports should avoid naked timings without configuration context.

## Initial Representative Configurations

The first pass does not need instrument-perfect models for every family. It
does need consistent, defensible scale points.

### ZernikeWFS

- `compact`
  - current `16 x 16` / `4 x 4` frame maintained runtime profile
- `medium`
  - increase pupil resolution and sampled frame enough to make GPU execution
    nontrivial
- `representative`
  - a realistic detector-backed closed-loop size intended for RTC-facing use

### Shack-Hartmann

- `compact`
  - current maintained diffractive runtime and mixed-asterism profiles
- `medium`
  - larger detector-backed NGS and LGS cases
- `representative`
  - AO188-class or larger realistic pixel-output closed-loop profile
  - may also be a fixed named HIL case such as REVOLT when the defining feature
    is the detector payload rather than only the command length

### Pyramid / BioEdge

- `compact`
  - current maintained stacked/runtime surfaces
- `medium`
  - larger pupil/reimaged-pupil grids
- `representative`
  - realistic closed-loop detector-backed profiles

### Tomography / Builders

- `compact`
  - current maintained builder surfaces
- `medium`
  - current model-tomography profiling sizes
- `representative`
  - KAPA-like or instrument-scale cases already used for reference validation

## Backend Questions This Matrix Should Answer

The matrix exists to answer specific engineering questions:

1. At what sizes does GPU become worthwhile for each sensing family?
2. Are CUDA and AMDGPU scaling similarly, or do they diverge by family?
3. Which algorithms are only acceptable as compact regression paths, and which
   remain appropriate at representative scale?
4. Which performance conclusions are genuinely structural, and which are merely
   artifacts of tiny workloads?

If a benchmark does not help answer one of those questions, it should probably
remain a regression benchmark instead of entering the representative matrix.

## Rollout Plan

### Phase 1: Classification

- tag current maintained scripts as regression or representative
- document which current runtime profiles are intentionally compact
- stop using compact profiles as implicit evidence for deployment-scale
  conclusions

### Phase 2: First Representative Runtime Ladder

- add `compact` / `medium` / `representative` modes for:
  - `ZernikeWFS`
  - Shack-Hartmann pixel-output runtime
  - stacked multi-source / multi-WFS runtime

This should be done by extending current profiling scripts rather than creating
entirely new ad hoc entry points where possible.

Current implementation status:

- done for the first three maintained runtime profilers listed above
- still pending for broader closed-loop and builder/runtime families

### Phase 3: Representative Closed-Loop and Builder Cases

- add representative AO188-class and larger closed-loop scenarios
- add representative tomography/calibration benchmark cases
- record CPU, AMDGPU, and CUDA snapshots where available

### Phase 4: Policy Integration

Use the resulting evidence to guide:

- which paths are CPU-first by design
- which paths justify GPU-specific optimization
- which defaults should prefer CPU-build / GPU-run
- which algorithms need redesign before they are suitable at realistic scale

## Success Criteria

This plan is complete when:

1. The package clearly distinguishes regression vs representative benchmarks.
2. Each major runtime family has at least one representative scale point.
3. CPU/GPU crossover conclusions are based on size ladders rather than a single
   compact case.
4. Performance discussions in docs and reviews can cite representative evidence
   instead of inferring too much from tiny workloads.

## Recommendation

Do not replace the compact regression surfaces. Keep them and label them
explicitly. Then add a small, curated representative matrix on top, starting
with:

1. `ZernikeWFS`
2. Shack-Hartmann pixel-output runtime
3. stacked multi-source / multi-WFS runtime
4. representative tomography/build cases

That is the smallest useful matrix that will materially improve performance
decision-making for the package.

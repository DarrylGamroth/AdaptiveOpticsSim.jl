# GPU Runtime Structural Refactor Plan

This document captures the next structural performance roadmap for maintained
GPU runtime paths.

The main conclusion from recent SHWFS, pixel-output, and REVOLT-oriented
profiling is that the remaining bottlenecks are increasingly architectural
rather than local Julia or kernel-level issues. The package already has real
GPU wins, but several core design choices still bias the runtime toward extra
buffer traffic, conservative synchronization, and frame-centric execution.

This is not a statement that the current design is incorrect. It is a statement
that some currently reasonable abstractions impose throughput costs that matter
for HIL-oriented GPU paths.

## Current Structural Constraints

### 1. Synchronization granularity is too eager

In [backends.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/core/backends.jl),
`launch_kernel!` fences immediately. That keeps correctness simple, but it also
means many hot paths are written as:

1. launch one kernel
2. synchronize
3. launch the next kernel

rather than as a phase-level dependency chain.

### 2. SHWFS internal buffers and exported products are too tightly coupled

In [runtime.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/control/runtime.jl),
diffractive SHWFS currently exports the same `spot_cube` surface that is also
used as an internal working set in the sensing path.

That coupling limits optimization freedom. Thresholding, centroiding, and
reduction passes cannot freely reorganize or discard intermediate data without
also changing the externally visible pixel product.

### 3. Detector/WFS composition is still frame-centric

The batched detector path in [frame_batched.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/detectors/frame_batched.jl)
is reasonably structured, but the surrounding WFS code still tends to:

1. materialize full cubes
2. run separate reduction passes
3. run separate centroid passes
4. then export or capture frames

That is modular and correct, but it is not ideal for GPU memory traffic.

### 4. Reduction strategy is still ad hoc

Several hot paths still depend on generic array-wide reductions, one-off
partial reductions, or host fallbacks embedded directly in model code.

That works, but it spreads backend policy into algorithm files and makes it
harder to reason about where GPU reductions are genuinely first-class.

## Goals

- Improve maintained GPU runtime throughput without weakening correctness.
- Keep backend-specific policy out of scientific model code where possible.
- Make product ownership explicit: slopes, pixels, and internal scratch should
  not be conflated unless there is a measured reason.
- Prefer reusable runtime and reduction infrastructure over model-local
  patches.

## Non-Goals

- Redesigning mature scientific algorithms only for stylistic reasons.
- Chasing benchmark-only wins that do not improve maintained runtime surfaces.
- Introducing GPU-only APIs that fracture CPU and GPU runtime behavior without
  a clear benefit.

## Recommended Phase Order

1. Separate exported WFS products from internal work buffers.
2. Make runtime product requirements explicit.
3. Add phase-level scheduling support around current kernel launch APIs.
4. Add reusable backend-aware reduction primitives.
5. Rework detector/WFS composition around explicit products.
6. Revisit FFT stack/layout ownership only after the earlier phases land.

## Phase 1: Separate Exported Products From Internal Buffers

### Scope

- Introduce explicit exported pixel buffers for diffractive WFS runtime paths.
- Stop using internal SHWFS centroid/work buffers as the runtime-facing frame
  contract.
- Keep the public runtime outputs unchanged unless a clearer export type is
  needed.

Primary files:

- [runtime.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/control/runtime.jl)
- [shack_hartmann.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/wfs/shack_hartmann.jl)
- any helper files that define output-frame metadata

### Rationale

This creates the basic freedom needed for later optimization. If the runtime
pixel product is distinct from the centroid workspace, we can fuse or eliminate
centroid-side passes without risking exported-frame regressions.

### Risks

- Medium API-adjacent risk: output ownership is visible at runtime boundaries.
- Moderate memory cost: more explicit buffers may temporarily increase memory
  use until follow-on phases reclaim redundant passes.
- Hidden assumptions may exist in scripts or profilers that currently expect
  `spot_cube` identity rather than logical pixel equivalence.

### Expected Benchmark Impact

- Small immediate throughput gain by itself.
- High enabling value for later phases.
- Some paths may look neutral at first because this phase improves ownership
  clarity more than raw runtime.

### Success Criteria

- SHWFS can mutate centroid scratch independently of exported pixel surfaces.
- Runtime interfaces expose a stable pixel product that is not also the
  algorithm’s main scratch buffer.
- Existing CPU/GPU parity for exported pixels remains intact.

## Phase 2: Make Runtime Product Requirements Explicit

### Scope

- Define whether a runtime step needs:
  - slopes only
  - pixels only
  - slopes plus pixels
- Teach prepared runtime setup to allocate and synchronize only the products
  actually required.

Primary files:

- [runtime.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/control/runtime.jl)
- any WFS export metadata helpers

### Rationale

A slopes-only loop should not pay pixel-export costs by accident. A pixel-plus-
slopes loop should express that requirement directly rather than reaching it
implicitly through shared buffers.

### Risks

- Medium runtime-contract risk if assumptions are currently implicit.
- Some benchmark and tutorial scripts may need minor updates to request the
  right products explicitly.

### Expected Benchmark Impact

- Moderate improvement on slopes-only GPU runtime paths.
- Small improvement on full pixel-plus-slopes paths.
- Better attribution of per-phase runtime cost in profilers.

### Success Criteria

- Runtime preparation is product-aware.
- Slopes-only execution avoids unnecessary pixel export synchronization.
- Pixel-plus-slopes execution remains supported without hidden coupling.

## Phase 3: Introduce Phase-Level Scheduling

### Scope

- Keep `launch_kernel!` as the safe default.
- Add a small scheduling pattern or helper layer for hot paths that need:
  - multiple kernel launches
  - one explicit fence at the phase boundary
- Refactor maintained WFS hot paths to use that pattern where appropriate.

Primary files:

- [backends.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/core/backends.jl)
- selected WFS hot-path files

### Rationale

Right now the package can queue kernels manually, but the architecture still
encourages eager synchronization. That is safe but leaves throughput on the
table, especially on GPU paths that already have clear phase structure.

### Risks

- Medium correctness risk around hidden ordering assumptions.
- Harder debugging if asynchronous phase boundaries become unclear.
- FFT boundaries require extra care and may still need explicit synchronization.

### Expected Benchmark Impact

- Moderate improvement on GPU-heavy WFS paths.
- Little or no change on CPU.
- Reduced synchronization overhead in multi-kernel phases.

### Success Criteria

- Hot paths can express “queue these kernels, then fence once” cleanly.
- No loss of determinism or numerical parity.
- Sync points are fewer but still obvious in the code.

## Phase 4: Add Reusable GPU Reduction Primitives

### Scope

- Add shared reduction helpers for:
  - max
  - masked max
  - sum
  - masked sum
  - small statistics reductions used in centroid-like paths
- Move backend-specific fallback policy out of model-local ad hoc branches and
  into reusable reduction infrastructure.

Primary files:

- core reduction helper file or `src/core/` addition
- [shack_hartmann.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/wfs/shack_hartmann.jl)
- other detector/WFS files that currently do model-local reductions

### Rationale

We have already seen that generic GPU reductions and one-off local workarounds
can become both correctness and performance problems. This should become shared
infrastructure rather than repeated local repair.

### Risks

- Medium implementation risk because reductions are both backend-sensitive and
  performance-sensitive.
- Poorly chosen primitives could become too generic to optimize well.

### Expected Benchmark Impact

- Moderate on SHWFS and detector-backed GPU paths.
- High maintainability value by centralizing backend policy.
- Some improvements may come more from removing fallback overhead than from raw
  kernel speed.

### Success Criteria

- Model code stops carrying bespoke reduction logic where shared primitives
  suffice.
- AMDGPU and CUDA paths use the same reduction abstraction, even if the backend
  implementation differs.
- New GPU reduction failures are localized to shared infrastructure instead of
  hidden in model code.

## Phase 5: Make Detector/WFS Pipelines Product-Aware

### Scope

- Rework frame-centric wfs/detector composition so the pipeline is driven by
  required products rather than by default full-frame materialization.
- Reduce avoidable full-cube passes where the detector, centroid, and export
  stages can share a more explicit contract.

Primary files:

- [frame_batched.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/detectors/frame_batched.jl)
- [shack_hartmann.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/wfs/shack_hartmann.jl)
- possibly other WFS families if the pattern generalizes

### Rationale

The current structure is correct but often makes the GPU pay for full cube
materialization and multiple post-processing passes even when the end product
is smaller or more constrained.

### Risks

- Medium-to-high complexity risk because detector semantics must remain stable.
- Larger refactor surface than the earlier phases.
- Benchmark wins may vary by WFS family.

### Expected Benchmark Impact

- Moderate to high on detector-backed SHWFS GPU paths.
- Smaller impact on CPU unless the pipeline also reduces allocations.
- Strong potential benefit for HIL-facing runtime loops that need both pixels
  and slopes.

### Success Criteria

- Fewer full-memory passes in maintained detector-backed WFS paths.
- Clearer ownership of pre-detector, post-detector, and exported products.
- End-to-end runtime improvements on representative GPU profiles.

## Phase 6: Revisit FFT Stack and Layout Ownership

### Scope

- Reevaluate whether FFT stacks, sampled spot stacks, detector stacks, and
  exported pixel stacks should share current shapes and ownership boundaries.
- Consider whether some FFT-heavy paths should move to alternate layouts that
  reduce sampling or reduction overhead.

Primary files:

- [shack_hartmann.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/wfs/shack_hartmann.jl)
- related FFT-heavy WFS files as needed

### Rationale

This is the phase most likely to unlock larger performance gains, but it is
also the phase most likely to destabilize otherwise-correct code. It should be
attempted only after product ownership, scheduling, and reductions are in
better shape.

### Risks

- High correctness risk.
- High backend sensitivity.
- Potential for benchmark-only wins that do not translate to maintained runtime
  surfaces if the redesign is too shape-specific.

### Expected Benchmark Impact

- Potentially high on GPU-heavy FFT-based sensing paths.
- Uncertain until the earlier structural cleanup is complete.
- Not justified as an immediate next step.

### Success Criteria

- A measurable end-to-end improvement on representative maintained GPU runtime
  surfaces.
- No erosion of CPU path clarity or reference parity.
- The new layout ownership is simpler to reason about, not merely faster in one
  benchmark.

## Validation and Benchmark Strategy

Every phase should be validated at three levels:

1. correctness
2. maintained runtime integration
3. representative benchmark impact

Recommended maintained benchmark surfaces:

- diffractive SHWFS slopes-only runtime
- diffractive SHWFS pixels-plus-slopes runtime
- detector-backed SHWFS runtime
- representative multi-source or multi-WFS runtime surfaces where relevant
- representative HIL-facing benchmark surfaces already tracked in
  [benchmark-matrix-plan.md](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/docs/benchmark-matrix-plan.md)

Backend expectation:

- CPU remains the reference path for behavior
- CUDA and AMDGPU remain maintained execution paths
- backend-specific fallback policy is allowed, but it should be isolated to
  shared infrastructure rather than embedded as ad hoc scientific-model logic

## Current ROCm Workaround Inventory

The current tree includes a small number of explicit ROCm-specific workaround
paths that were added to restore realistic maintained benchmark coverage
without changing global backend policy.

These are acceptable only as transitional measures ahead of the structural
refactor phases above.

### 1. Detector-owned host mirrors for ROCm noise/readout correction

Files:

- [frame_capture.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/detectors/frame_capture.jl)
- [interface.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/detectors/interface.jl)
- [frame_response.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/detectors/frame_response.jl)
- [hgcdte_avalanche_array.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/detectors/hgcdte_avalanche_array.jl)

Current behavior:

- `DetectorState` owns a preallocated `noise_buffer_host`
- ROCm detector paths that still trip compiler failures in generic GPU
  reductions/noise kernels use that host mirror instead of allocating
  `Array(frame)` inside helper functions

Status:

- acceptable temporary maintained fallback
- should be revisited in Phase 4 and Phase 5 of this plan

### 2. WFS-owned host mirrors for ROCm calibration metadata and masked flux sums

Files:

- [pyramid.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/wfs/pyramid.jl)
- [bioedge.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/wfs/bioedge.jl)

Current behavior:

- `PyramidState` and `BioEdgeState` own cached host buffers for valid-mask
  metadata and ROCm masked-flux fallback
- the ROCm path copies through cached state rather than allocating inside
  hot-path helpers

Status:

- acceptable temporary maintained fallback
- should be replaced by shared reduction infrastructure in Phase 4

### 3. ROCm-specific Shack-Hartmann safe path

File:

- [shack_hartmann.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/wfs/shack_hartmann.jl)

Current behavior:

- diffractive SH on ROCm currently avoids the stacked GPU path for several
  measurement/reduction surfaces
- the safe path reuses cached host vectors and mask mirrors owned by SH state

Status:

- benchmark-unblock path, not desired final structure
- should be superseded by Phases 1, 4, and 5 of this plan

### Policy

These workarounds are acceptable because:

- they are localized to detector/WFS runtime ownership, not global backend
  policy
- they reuse owned state instead of allocating ad hoc in helper functions
- they restored realistic maintained AMDGPU benchmark coverage

They should not be treated as the final GPU architecture.

## Recommendation

The next concrete implementation work should start with Phases 1 and 2, not
with deeper kernel tuning.

That is the highest-leverage sequence because:

- it improves optimization freedom without forcing an algorithm rewrite
- it reduces hidden coupling in runtime behavior
- it makes later GPU scheduling and reduction work much cleaner

Only after those phases land should we decide whether the remaining GPU
bottlenecks justify the higher-risk FFT/layout redesign work.

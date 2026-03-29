# Pixel-Output GPU Runtime Redesign Plan

This document captures the redesign plan for a competitive GPU runtime path
when the RTC boundary requires detector pixels, not just slopes.

The original conclusion was that the existing GPU runtime path was not a good
optimization target and needed a redesign. That redesign is now mostly in
place for the maintained Shack-Hartmann AO188 simulation, and the question has
shifted from "can GPU runtime be made competitive at all?" to "which backend-
specific execution model should be preferred?"

## Motivation

The AO188/3k simulation in `examples/support/subaru_ao188_simulation.jl` is the
current reference case for pixel-output HIL:

- `64 x 64` DM grid
- `3228` active actuators
- `188` high-order control modes
- diffractive Shack-Hartmann high-order WFS
- `2 x 2` low-order Shack-Hartmann branch
- detector-backed pixel output on both branches

Initial warmed frame-loop results with `JULIA_NUM_THREADS=4` before the
redesign:

- AMD host CPU: about `142 Hz`
- AMDGPU: about `22.7 Hz`
- `spiders` CPU: about `116 Hz`
- CUDA on `spiders`: about `34.6 Hz`

That initial gap motivated the batched SH, batched detector, and low-order
branch redesign tracked below.

## Progress Snapshot

Current implementation status:

- Phase 1: implemented
  - dedicated profiler: `scripts/profile_pixel_output_runtime.jl`
  - diffractive SH runtime boundary now exports the 3D `spot_cube` for RTC use
- Phase 2: implemented
  - accelerator detector-backed diffractive SH now fills `spot_cube` first and
    no longer calls per-subaperture `capture!`
  - accelerator non-LGS diffractive SH now performs a batched `2D` FFT across
    the full spot stack instead of executing one small FFT per subaperture in a
    Julia loop
  - accelerator non-LGS diffractive SH now samples the whole spot stack with a
    single crop/bin kernel instead of looping over per-spot `center_resize2d!`
    and `bin2d!`
- Phase 3: implemented
  - added `capture_stack!` for batched detector application on spot stacks
  - the simple LGS branch (`LGSProfileNone`) now reuses the stacked SH field
    path and applies elongation across the full intensity stack
  - the sodium-profile LGS branch (`LGSProfileNaProfile`) now applies its
    subaperture-specific Fourier-domain convolution across the stacked
    intensity cube instead of looping spot-by-spot at runtime
- Phase 4: experimental, now backend-split
  - added generic execution policies:
    `SequentialExecution`, `ThreadedExecution`, and
    `BackendStreamExecution`
  - task-parallel branch execution is not the default
  - AMDGPU still does not justify generic Julia task overlap
  - CUDA now shows a real branch-overlap gain on the maintained AO188 case
- Phase 5: experimental
  - added `DirectReplayMode` and `PreparedReplayMode`
  - prepared replay currently means up-front SH sampling/calibration prep, not
    backend-specific graph capture

Current measured AO188 profile after the stacked SH/LGS redesign and the
low-order-resolution branch split:

| Host | Backend | Mode | Build time | Warmed step | Frame rate |
| --- | --- | --- | ---: | ---: | ---: |
| AMD host | CPU | sequential/direct | `10.99 s` | `6.27 ms` | `159.4 Hz` |
| AMD host | AMDGPU | sequential/direct | `47.64 s` | `4.83 ms` | `207.1 Hz` |
| AMD host | AMDGPU | task/direct | `45.18 s` | `19.56 ms` | `51.1 Hz` |
| `spiders` | CUDA | sequential/direct | `47.59 s` | `2.22 ms` | `450.4 Hz` |
| `spiders` | CUDA | task/direct | `47.61 s` | `1.90 ms` | `526.5 Hz` |

Compared to the earlier AO188 audit:

- local AMDGPU improved from about `22.7 Hz` to about `207.1 Hz`
- CUDA improved from about `34.6 Hz` to about `450.4 Hz` in sequential mode
- CUDA improves further to about `526.5 Hz` with branch overlap
- CPU is no longer the unquestioned winner for the maintained AO188 simulation
- AMDGPU is faster than the local CPU in sequential mode, but still strongly
  penalized by generic Julia task overlap
- keeping the AO188 reconstructor in mapped two-stage form (`modal_reconstructor`
  then `M2C`) is a worthwhile fast-path refinement, but it does not remove the
  stricter post-command `Float32` backend-equivalence gap on its own

Latest follow-up on the mapped path:

- local CPU sequential/direct is about `183.4 Hz`
- local AMDGPU sequential/direct is about `214.6 Hz`
- `spiders` CUDA sequential/direct is about `633.0 Hz`
- the maintained fast-runtime equivalence scripts still pass on AMDGPU and CUDA

Latest follow-up on the structured DM runtime split:

- interaction-matrix calibration stays on the dense DM application path
- runtime DM application uses a structured separable `X * C * Y'` path when the
  Gaussian influence basis remains separable under the configured
  misregistration
- AO188 GPU calibration now defaults to CPU-built interaction matrices and
  reconstructors, with the resulting operators uploaded to the runtime backend
- local CPU sequential/direct is now about `1.02 kHz`
- local AMDGPU sequential/direct is now about `1.02 kHz`
- `spiders` CUDA sequential/direct is now about `1.25 kHz`
- the maintained fast-runtime AO188 equivalence now passes tightly on both
  AMDGPU and CUDA, including the command surface
- stricter post-command `Float32` AO188 `tel_opd` agreement improves on AMDGPU
  and improves materially on CUDA too, from about `2.68e-7` to about
  `1.49e-7`, though it still misses the scientific high-accuracy tolerance

## Current Diagnosis

The dominant cost is still the high-order diffractive Shack-Hartmann pixel
path, but the shape of the bottleneck has changed materially.

The current accelerator implementation in
[`src/WFS/shack_hartmann.jl`](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/src/WFS/shack_hartmann.jl)
now has this structure for the maintained non-LGS path:

1. build a stacked complex field cube for all subapertures
2. execute one batched `2D` FFT across that stack
3. convert the transformed stack to intensities
4. crop/bin each spot into `spot_cube`
5. apply detector effects to the whole stack with `capture_stack!`

This removed the worst previous mismatch between the algorithm and the GPU:
`n_subap^2` tiny FFT executions. The remaining bottlenecks are now:

- the high-order SH branch itself, which still dominates the AO188 frame loop
- the low-order branch integration cost, though it is much lower after moving
  it onto a dedicated lower-resolution telescope
- backend-specific branch scheduling, where CUDA benefits from overlap and
  AMDGPU currently does not

The remaining LGS work is now mostly measurement, not structure. Both the
simple elongation branch and the sodium-profile branch use the stacked SH
runtime path. The sodium-profile branch still requires per-subaperture kernel
precomputation, but the runtime convolution is no longer spot-by-spot.

Local LGS SH runtime snapshot (`14 x 14` diffractive SH, detector-backed):

| Backend | LGS profile | Warmed measure | Frame rate |
| --- | --- | ---: | ---: |
| CPU | none | `0.552 ms` | `1811 Hz` |
| CPU | sodium profile | `0.479 ms` | `2087 Hz` |
| AMDGPU | none | `0.406 ms` | `2465 Hz` |
| AMDGPU | sodium profile | `0.487 ms` | `2053 Hz` |
| CUDA (`spiders`) | none | `0.273 ms` | `3660 Hz` |
| CUDA (`spiders`) | sodium profile | `0.336 ms` | `2976 Hz` |

Interpretation:

- the stacked LGS runtime path is now strong on both AMDGPU and CUDA
- the remaining work is no longer LGS structure; it is full-loop integration
  and backend-specific scheduling

Mixed NGS/LGS SH runtime follow-up:

- the incompatible SH asterism path now accumulates a real exported
  `spot_cube` instead of leaving the last-source frame in the readout
- the accelerator mixed SH centroid path now matches the CPU scalar fallback
  thresholding semantics by using per-spot peak thresholds
- local warmed detector-backed mixed SH runtime is now about:
  - CPU: `0.920 ms` (`1087 Hz`)
  - AMDGPU: `0.862 ms` (`1160 Hz`)
- focused AMDGPU mixed SH equivalence is now tight:
  - `spot_cube` max abs about `3.91e-2`
  - `slopes` max abs about `1.62e-6`

Branch-overlap recheck after the SH/LGS redesign:

- generic Julia task overlap is still a bad fit for AMDGPU on the maintained
  AO188 simulation (`207 Hz` sequential vs `51 Hz` task mode)
- CUDA now sees a meaningful gain from branch overlap over sequential mode
  (`450 Hz` sequential vs about `650-658 Hz` with overlap on the current code)
- if branch overlap is pursued further, it should be as a CUDA-specific stream
  execution path, not a generic task-parallel default

CUDA direct overlap comparison on the current code:

- `task` mode: about `658 Hz`
- `stream` mode: about `649 Hz`

Interpretation:

- explicit streams are now wired through the backend as intended
- they do not yet outperform the existing CUDA task-overlap path
- the value of the stream mode today is architectural: it provides a backend-
  specific execution hook we can refine further without treating generic Julia
  task overlap as the long-term design

Low-order branch recheck after moving it to a dedicated lower-resolution
telescope:

- local CPU low-order sensing: about `0.072 ms`
- local AMDGPU low-order sensing: about `0.535 ms`
- CUDA low-order sensing on `spiders`: about `0.809 ms` sequential and
  about `0.325 ms` with branch overlap

This means the low-order branch is no longer distorting the AO188 model by
forcing the `2 x 2` WFS to run on the full `112 x 112` pupil resolution.

For comparison, the diffractive Pyramid detector path is structurally more
GPU-friendly today: it forms one full-frame intensity image and applies one
detector capture to that frame, rather than looping over hundreds of
subaperture spots. That does not automatically make Pyramid faster in every
regime, but it does mean the current AO188-class GPU limitation is primarily a
diffractive Shack-Hartmann limitation, not a universal WFS limitation.

## Goals

The redesign should satisfy these requirements:

- preserve pixel output at the RTC boundary
- preserve existing scientific detector semantics unless an explicit fast
  profile is selected
- keep CPU behavior intact
- reduce GPU per-frame overhead enough that AO188-class runtime loops are
  plausibly competitive with or better than CPU on supported backends
- avoid backend-specific duplication where a generic batched design will do

## Non-Goals

These are not redesign goals:

- replacing pixel output with slope-only output
- batching multiple time steps together
- making the compact single-loop GPU path look good via benchmark artifacts
- adding backend-specific micro-optimizations without changing the pipeline
  shape

## Redesign Direction

The likely viable direction is a batched/fused per-WFS pixel pipeline.

### Target execution model

For each WFS branch:

1. generate all valid subaperture complex fields into a stacked backend-native
   buffer
2. propagate all subapertures through the detector-facing image-formation path
   in batched form
3. apply detector effects in batched form
4. write the final detector-output stack directly into the exported pixel
   buffer
5. derive slopes/reference-normalized signals from that same stacked output

The current `spot_cube` storage is already close to the right output layout.
The main issue is how it is populated.

### Design principles

- move from per-subaperture host loops to batched backend-native kernels and
  FFT/image operations
- synchronize once per branch or once per frame, not once per spot
- reuse preallocated per-branch workspaces
- keep low-order and high-order branches independent after telescope/DM state
  is prepared
- keep CPU and GPU implementations sharing the same math where practical, but
  allow different execution strategy

## Candidate Work Packages

### Phase 1: Measurement and boundary cleanup

- add a dedicated profiler for pixel-output WFS runtime, separate from builder
  profiling
- split per-frame timing into:
  - common telescope/DM preparation
  - high-order pixel path
  - low-order pixel path
  - reconstruction
  - DM application
- make the exported detector-pixel buffers explicit for both branches

Exit criterion:
- we can attribute runtime cost to concrete per-branch phases without relying
  on ad hoc manual timing

Status:
- done
- current profiler entry point:
  `scripts/profile_pixel_output_runtime.jl`
- dedicated LGS-heavy branch profiler:
  `scripts/profile_lgs_sh_runtime.jl`
- dedicated mixed NGS/LGS SH profiler:
  `scripts/profile_mixed_sh_asterism_runtime.jl`

### Phase 2: Batched SH spot generation

- redesign diffractive Shack-Hartmann spot generation so all valid
  subapertures are prepared into stacked backend-native arrays
- remove the per-subaperture `capture!` loop from the accelerator path
- keep the current CPU path as a correctness reference during transition

Exit criterion:
- no per-subaperture detector capture loop remains in the accelerator path for
  diffractive SH

Status:
- done for the maintained diffractive SH detector path
- the accelerator path now fills `spot_cube` directly before detector effects
  are applied in batch
- the maintained non-LGS path is now fully batched through spot sampling
- the simple LGS elongation path is now batched through stacked elongation
- the sodium-profile LGS convolution path is now batched at runtime
- the next LGS task is benchmarking and profiling a maintained LGS-heavy case

### Phase 3: Batched detector application

- add detector operations that can apply photon/readout/background/binning to a
  stack of spots rather than one spot at a time
- make this an explicit fast path for stacked pixel pipelines instead of trying
  to force the scalar detector API into a batched runtime role
- preserve the existing single-image detector API for non-batched workflows

Exit criterion:
- detector application for the SH runtime path is performed on a stacked pixel
  buffer with one synchronization boundary per branch at most

Status:
- implemented for the maintained spot-stack path via `capture_stack!`
- the fixed-shape in-place fast path is intentionally limited to
  `psf_sampling == 1`, `binning == 1`, `output_precision === nothing`,
  full-frame readout, global-shutter timing, and null persistence
- maintained batched response models now include null, Gaussian,
  rectangular/separable MTF, and sampled frame responses
- maintained batched readout correction now includes the null and reference-edge
  common-mode correction models
- a second generalized `capture_stack!` path now exists for shape-changing
  detector configurations with separate input/output stacks
  - this generalized path covers `psf_sampling`, `binning`,
    `readout_window`, and `output_precision`
  - it is intentionally not the hot HIL fast path
  - it still requires global-shutter timing and null persistence

### Phase 4: Branch concurrency

- run high-order and low-order WFS branches as independent GPU work after the
  common telescope/DM state is ready
- only synchronize when both branches are needed for command assembly or RTC
  output export

Exit criterion:
- branch overlap is measurable on GPU hosts and does not complicate CPU code

Status:
- experimental branch-execution modes are implemented
- AMDGPU measurements do not justify enabling generic task-parallel branches
  by default
- CUDA measurements now justify continued investigation of backend-specific
  branch overlap
- `BackendStreamExecution` now has a real CUDA measurement on `spiders`
  and falls back safely on non-CUDA backends

### Phase 5: Graph/replay execution

- evaluate CUDA Graphs and the AMD equivalent only after the batched pipeline
  exists
- keep this backend-specific and optional

Exit criterion:
- only keep graph/replay support if it materially reduces warmed frame latency

Status:
- backend-neutral prepared replay scaffold is implemented
- backend-specific graph capture is still future work
- prepared replay is not yet justified as the default runtime mode

## API and Surface Changes

The redesign should prefer new internal execution helpers over broad public API
changes.

Possible new internal concepts:

- `BatchedPixelWFSWorkspace`
- batched detector helpers for spot stacks
- explicit pixel-export buffers on `SimulationInterface` / `CompositeSimulationInterface`

Public API changes should stay minimal:

- runtime outputs remain pixels
- existing `measure!` entry points continue to work
- profile selection may choose between scalar/reference and batched/runtime
  implementations where scientifically acceptable

## Validation Plan

Each redesign phase needs both correctness and performance gates.

Correctness:

- CPU batched vs existing CPU scalar path agreement
- GPU batched vs CPU reference agreement
- AO188 simulation end-to-end closed-loop agreement within tolerance
- detector-output agreement on representative pixel frames, not just slopes

Performance:

- warmed `step!` comparison on:
  - AMD CPU vs AMDGPU
  - `spiders` CPU vs CUDA
- branch-level phase timings
- allocation checks on CPU hot paths
- explicit go/no-go comparison against the current CPU baseline

## Go / No-Go Criteria

We should keep investing in the GPU runtime redesign only if the remaining
backend-specific work continues to deliver meaningful end-to-end gains on
AO188-class runtime loops.

Suggested threshold:

- a clear warmed frame-rate improvement over the previous GPU baseline
- and a plausible path toward beating CPU on at least one supported GPU
  backend for AO188-scale HIL

That threshold has now been met for the maintained AO188 simulation on both
AMDGPU and CUDA. The remaining question is where further investment should go.

If future backend-specific work stops delivering real gains, the correct
project decision is:

- keep GPU support focused on builder-heavy calibration/tomography workflows
- keep runtime GPU support as smoke-tested but not performance-prioritized
- avoid growing backend-specific maintenance burden for the pixel-output loop

## Recommended Immediate Next Step

Keep Phase 1-3 and the low-order-resolution split as the maintained baseline.
The next work item should be backend-specific:

- for CUDA: evaluate whether explicit stream scheduling can safely replace the
  generic task-overlap experiment
- for AMDGPU: avoid generic task overlap and focus on stacked multi-WFS /
  multi-branch scaling, where more independent work can amortize launch cost
- for both backends: extend the stacked execution model to future multi-WFS or
  asterism pixel-generation paths rather than revisiting scalar SH internals

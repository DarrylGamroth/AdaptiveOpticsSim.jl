# Performance Audit

This document records a measured performance review of the current
`AdaptiveOpticsSim.jl` implementation after the recent no-behavior-change cleanup
pass.

Scope of this audit:

- LiFT
- GainSensingCamera
- interaction-matrix construction
- tomography reconstruction and command assembly

The goal was to answer three questions:

1. Are the main entry points type-stable?
2. Which paths still allocate materially?
3. Which allocations are acceptable setup costs versus worth removing?

## Method

The audit used small deterministic setups derived from existing tests and
tutorials and measured:

- warm steady-state allocations with `@allocated`
- type stability with `code_warntype`

The checks were run inside functions after warmup to avoid global-scope noise.

## Results

### Type Stability

The following entry points were checked and are currently type-stable in the
measured cases:

- `compute_optical_gains!`
- `lift_interaction_matrix!`
- `AdaptiveOpticsSim.reconstruct` for LiFT
- `interaction_matrix`
- `build_reconstructor(ModelBasedTomography(), ...)`
- `assemble_reconstructor_and_fitting`
- `reconstruct_wavefront!`
- `dm_commands!`

This is the important baseline result: the current issues are not primarily
inference failures.

### Allocation Measurements

Measured allocations for representative calls:

| Path | Allocation |
|---|---:|
| `compute_optical_gains!` | `250016` bytes |
| `calibrate!` for `GainSensingCamera` | `596728` bytes |
| `lift_interaction_matrix!` | `40120` bytes |
| `AdaptiveOpticsSim.reconstruct(lift, ...)` | `132096` bytes |
| `LiFT(...)` constructor | `351920` bytes |
| `interaction_matrix(dm, wfs, tel, src)` | `10304` bytes |
| `reconstruct_wavefront!` | `3360` bytes |
| `dm_commands!` | `3024` bytes |
| `build_reconstructor(ModelBasedTomography(), ...)` | `184818344` bytes |
| `assemble_reconstructor_and_fitting(...)` | `14832` bytes |

## Interpretation

### 1. GainSensingCamera is the clearest runtime allocation target

`compute_optical_gains!` and `calibrate!` both allocate substantially for what
should ideally be a reusable analysis path.

The main reason is visible in `src/Calibration/gain_sensing_camera.jl`:

- `frame_norm = frame ./ total`
- `impulse_response(...)` returns a new matrix
- `sensitivity(...)` returns a new vector
- calibration also rebuilds `basis_product`

This means the current GSC API is logically mutating but still internally
allocating intermediate arrays.

### 2. LiFT runtime still allocates in the reconstruction loop

The constructor allocation is acceptable as setup work.

The two more important LiFT findings are:

- `lift_interaction_matrix!` still allocates nontrivially
- `AdaptiveOpticsSim.reconstruct(lift, ...)` allocates more noticeably

Likely sources in `src/WFS/lift.jl`:

- `coeffs = copy(coeffs0)` or a new `zeros(...)`
- `diagw = Diagonal(sqrtw)`
- `delta = normal \\ rhs` allocates a solution vector
- singular fallback `pinv(normal) * rhs` allocates heavily if triggered
- the non-`!` API returns a fresh coefficient vector by design

So the LiFT implementation is type-stable but not yet close to a zero-allocation
iterative solve path.

### 3. interaction-matrix construction allocates, but this is currently acceptable

`interaction_matrix(...)` allocates only about `10 KB` in the measured compact
case.

This path constructs a new calibration matrix object and snapshots OPD state, so
some allocation is expected. It is not currently the highest-priority
performance problem.

If we later want a more aggressive RTC calibration path, we should add an
explicit workspace-backed `interaction_matrix!` API rather than trying to force
the constructor-like function to be allocation-free.

### 4. Tomography runtime is cheap; tomography assembly is expensive

The runtime application paths:

- `reconstruct_wavefront!`
- `dm_commands!`

still allocate a few KB in the measured compact case, but they are small enough
that they are not the dominant concern yet.

The real tomography cost is build-time:

- `build_reconstructor(ModelBasedTomography(), ...)` allocated about `185 MB`

This is not surprising. The current model-based build path forms dense
covariance and signal matrices in `src/Tomography/reconstructors.jl`.

That cost is acceptable if tomography build is treated as setup, but it is not
acceptable if users need to rebuild reconstructors repeatedly inside a control
loop.

## Priority Fixes

### High Priority

1. Add in-place GSC kernels
   - Add a normalized-frame scratch buffer.
   - Add `impulse_response!`.
   - Add `sensitivity!`.
   - Reuse preallocated IR and sensitivity buffers in `GainSensingCamera`.

2. Add an in-place LiFT reconstruction path
   - Add `reconstruct!(coeffs_out, lift, ...)` or equivalent stateful update.
   - Reuse a preallocated `delta` buffer.
   - Avoid constructing `Diagonal(sqrtw)` in the loop.
   - Keep the allocating convenience wrapper as a thin outer API.

### Medium Priority

3. Add workspace-backed interaction matrix construction
   - Keep `interaction_matrix(...)` as the allocating convenience API.
   - Add `interaction_matrix!(workspace, ...)` for repeated calibration use.

4. Reduce small tomography runtime allocations
   - Audit why `mul!`-based `reconstruct_wavefront!` and `dm_commands!` still
     allocate a few KB.
   - Remove those only if the fix is low-risk and does not complicate the API.

### Low Priority

5. Tomography build-time memory reduction
   - Only worth doing if tomography reconstructors need to be rebuilt often.
   - Likely requires algorithm/data-layout work rather than local cleanup.

## Non-Issues

These are not currently performance problems in the measured cases:

- row-major compatibility handling, which remains confined to the reference
  harness rather than the library core
- `eltype(field)` lookups on concretely typed state fields
- the current scalar CPU versus `KernelAbstractions` split

## Recommended Next Step

The highest-value follow-up is:

1. make `GainSensingCamera` internally in-place,
2. then add an in-place LiFT reconstruction API.

Those are the two clearest measured wins from this audit.

## Follow-Up Results

The recommended GSC and LiFT refactors have now been implemented and remeasured
on the same compact deterministic setups.

Updated allocations:

| Path | Before | After |
|---|---:|---:|
| `compute_optical_gains!` | `250016` bytes | `0` bytes |
| `calibrate!` for `GainSensingCamera` | `596728` bytes | `18712` bytes |
| `lift_interaction_matrix!` | `40120` bytes | `1696` bytes |
| `AdaptiveOpticsSim.reconstruct` for LiFT (`Auto`, scalar CPU) | `132096` bytes | `2432` bytes |
| `reconstruct!` for LiFT (`Auto`, scalar CPU) | n/a | `2208` bytes |
| `reconstruct!` for LiFT (normal-equation mode) | n/a | `1984` bytes |

Interpretation:

- `GainSensingCamera` is now effectively in-place for steady-state optical-gain
  evaluation.
- `calibrate!` still allocates, but only at setup scale for basis-product and
  calibration-buffer construction.
- LiFT now defaults to `LiFTSolveAuto()`.
- `LiFTSolveAuto()` resolves to QR on the scalar CPU path and to the
  normal-equation path on accelerator backends, which keeps the design GPU-ready
  without pretending we have a validated GPU QR implementation yet.
- On scalar CPU, `Auto` prefers QR for conditioning but may still fall back to a
  damped normal-equation solve when the Jacobian is effectively rank-deficient.

That means the next worthwhile LiFT improvement, if needed, is adding damping or
explicit conditioning diagnostics on top of the QR solve rather than chasing
smaller setup allocations.

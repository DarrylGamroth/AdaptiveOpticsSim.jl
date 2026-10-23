# Centering precision: completed CPU/CUDA validation, partial performance

## Outcome

**All nine complete numerical comparisons pass the original frozen gates.**
All eight configurations ran all 17 probes and all 277 positive/negative
calibration pairs. Maximum raw-slope component difference is 1.35313 × 10⁻⁶,
below the unchanged 2 × 10⁻⁶ limit. No tolerance changed; protocol SHA256:
`47915d00908268b1cd0c136a8abf0969d93fcd8a27e2e6c1d73a94386abd972d`.

**The user stopped the remaining performance campaign.** Six round-one timing
runs completed. SPECULA CPU round one was interrupted; AOS CUDA timing and
rounds two/three were not reached. No missing service time or allocation result
is inferred. Full CUDA capture and numerical qualification did complete.

The passing SPECULA comparison uses an explicit cold prism precision correction
in the matching adapter. **SPECULA's original Float32 cold arithmetic still
fails the raw-slope gates.** Its negative evidence is retained and is not
represented as a pass.

## Changes and mathematical evidence

AOS production commit `69f91214` evaluates the same centering factor,
exp(−iπ(n+1)(i+j−2)/n), with at least Float64 phase arithmetic and `cispi`, then
rounds to the existing Complex{T} storage. Accelerators receive the same
host-built factor during preparation. The obsolete setup kernel and its
coverage entry were removed. Optical conventions, per-frame propagation,
runtime precision and array sizes are unchanged. At FFT512, CUDA preparation
adds one transient 2 MiB host array and a setup upload; runtime adds neither.

The independent factor regression failed 32 checks before and passes all 160
after. AOS/OOPAO realized centering-factor maximum difference falls from
2.2483 × 10⁻⁴ to 5.3380 × 10⁻¹⁴. Full-array CPU/CUDA factor equality is tested.

Fixing AOS alone does not satisfy both original external outputs. Direct
OOPAO versus original-arithmetic SPECULA raw-slope differences are
1.20618 × 10⁻⁵ on CPU and 1.21416 × 10⁻⁵ on CUDA. They exceed twice the
2 × 10⁻⁶ allowance: by the triangle inequality, no AOS output could be within
that allowance of both.

A separate two-probe SPECULA intervention identifies its cold prism exponential:
its native signed surface has exact quarter-integer samples, whose ideal
exponentials are quarter-turn unit roots. Float32 phase multiplication before
exponentiation introduces error. Prism-only Float64 evaluation reduces flat
and 20 nm x-tilt raw-slope errors from about 1.12 × 10⁻⁵ to about 3 × 10⁻⁷;
changing only its centering factor does not resolve them. The adapter now
exponentiates that same native surface in Float64, rounds to the original
complex dtype and uploads before setup/capture. Native propagation, modulation,
PSF outputs, normalization and rebinning remain inherited. No fitted gain,
registration, flat subtraction or calibration-column rescaling is used.
External source checkouts remain clean and unmodified.

All 17 current-fixture probes with original SPECULA arithmetic are retained
on both CPU and CUDA, with maximum raw-slope errors of 1.21948 × 10⁻⁵ and
1.22214 × 10⁻⁵ against corrected AOS. The original adapter snapshot, raw images,
records, external disagreement and causal diagnostic are in the core archive.
The [original campaign](../2026-10-22-matched-copper/README.md) also retains its
full original-arithmetic calibration evidence.

## Complete numerical results

Shared physical inputs are byte-identical to the original campaign: Copper
240/96, FFT512² to 64² photon rates, 1.22 m aperture, 550 nm wavelength,
nominal 5 λ/D modulation, Float32/ComplexF32 runtime. The documented endpoint
modulation convention gives actual radius 5 × 240/239 λ/D. Slopes and all 277
calibration columns use the same independent shared Float64 I4Q readout and
analytic incident-flux divisor. No modes are omitted.

| Comparison | Max rate relative L2 | Max raw slope absolute | Calibration relative L2 | Worst column relative L2 |
|---|---:|---:|---:|---:|
| cpu-oopao | 5.40976e-07 | 1.35313e-06 | 8.25434e-06 | 0.000184084 |
| cpu-specula | 5.42221e-07 | 1.28996e-06 | 8.23901e-06 | 0.000189477 |
| cuda-oopao | 5.42442e-07 | 1.27999e-06 | 7.97546e-06 | 0.000174909 |
| cuda-specula | 5.86544e-07 | 1.32653e-06 | 8.16922e-06 | 0.000186007 |
| aos-cpu-cuda | 2.05536e-07 | 7.58018e-07 | 6.57094e-06 | 0.000145853 |
| oopao-cpu-cuda | 1.92063e-07 | 5.48566e-07 | 4.60401e-06 | 0.000103885 |
| specula-cpu-cuda | 3.70702e-07 | 4.38853e-07 | 5.21009e-06 | 0.00013284 |
| aos-cpu-threading | 7.36336e-08 | 2.37712e-07 | 4.15648e-06 | 9.53766e-05 |
| oopao-cpu-threading | 0 | 0 | 0 | 0 |

The analysis archive preserves per-probe gates, photon conservation and piston
checks, full calibration matrices and every column error. Independent reducers
use separate quadrant equations without importing the benchmark evaluator.
This validates matched optical products and their response; it does not qualify
the external native estimators, DM models or closed-loop controllers.

## Partial service-time results

Completed measurements below each used 32 warmups and 128 varied-input samples,
with original affinity/thread limits. They are **single-round observations**,
not the planned three-round aggregates. All nine numerical comparisons passed
before any timing process began. The interrupted SPECULA CPU record has no
completed service series and contributes no time estimate.

| Completed configuration | Round | Median, ms | p95, ms | Host allocation observations |
|---|---:|---:|---:|---|
| aos-cpu1 | 1 | 362.146 | 436.173 | 8 × 0 bytes |
| aos-cpu8 | 1 | 257.003 | 315.493 | 8 × 0 bytes |
| oopao-cpu1 | 1 | 499.032 | 567.786 | not measured by this runner |
| oopao-cpu8 | 1 | 207.768 | 226.246 | not measured by this runner |
| oopao-cuda | 1 | 6.386 | 6.583 | not measured by this runner |
| specula-cuda | 1 | 4.840 | 5.284 | not measured by this runner |

The boundary is completed resident-OPD-to-photon-rate optical-node service time.
Input staging, atmosphere, DM runtime, detector noise, shared slope/calibration
reduction and RTC are excluded. **These are not full simulation frame times.**
AOS CUDA capture passed the complete numerical campaign, but its fresh timing
and allocation observations were not reached. The previous campaign's AOS
CUDA value must not be relabeled as this fix's measurement.

AOS uses one Julia/BLAS thread with one or eight FFTW cores; this optical node
has no BLAS work. OOPAO CPU uses one or eight jobs with one SciPy FFT worker
per job. SPECULA CPU retains serial NumPy FFT. The shared host has other
AO/GUI workloads and uncontrolled clocks; telemetry is preserved. Single-round
times cannot establish a performance regression or improvement. Completed AOS
CPU observations are close to previous medians (364 ms / 257 ms), and the two
CPU allocation checks report eight zero-byte observations each.

## SPECULA workload versus implementation

Both engines propagate 96 forward/inverse 512² FFT pairs under the matched
prescription. SPECULA additionally accumulates focal-plane PSFs per modulation
point, shifts/masks/scales PSFs, computes transmission reductions, and performs
additional normalization and image movement. These native extra products remain
included in its timed path. AOS publishes the requested detector photon rate.

Implementation differs too: SPECULA forms a 96 × 512 × 512 ComplexF32 modulated
pupil stack (192 MiB), then transforms positions individually. AOS CUDA uses
batch 2, common pupil-phase reuse, separable modulation storage and fused
preparation/accumulation. AOS captures the full node; SPECULA captures its
trigger with preparation/postprocessing outside that graph. CPU FFT backends
and planning also differ. Actual FFT dtypes are recorded; ComplexF64 promotion
is not assumed for this NumPy version.

Thus extra work and implementation differences both contribute possible
mechanisms. Their timing shares were not measured, so neither diagnostics nor
optimization alone is established as the dominant cause. The original campaign
measured about 6.0× single-core CPU and 1.83× CUDA service-time ratios. This
stopped campaign has no completed SPECULA CPU or AOS CUDA timing and provides
no replacement ratio. Source evidence is in `specula-work-comparison.md` in
the core archive.

## Verification and reproducibility

- Targeted Pyramid/Bi-O-edge CPU suite: 332/332 checks pass.
- KA CPU matrix and kernel inventory: 1,038/1,038 pass.
- Focused CUDA centering test: 66/66 pass, scalar indexing prohibited.
- Centering regression: 32 failures before; 160/160 pass after.
- Eight full validation configurations; all nine comparisons pass.
- Six completed timing records preserve all samples, exact replay and
  validation/performance probe parity. No full-suite or remote CI run.
- Initial Pkg wrapper invocation used an anonymous benchmark environment and
  failed before tests; the direct selector invocation above succeeded. Both
  logs are preserved.

Hardware/dependencies match the original campaign: Ryzen 7 6800H, RTX 3080
10 GiB, Julia 1.13.2, CUDA.jl 6.4.2, KA 0.9.45, Python 3.13.5, NumPy 2.2.4,
SciPy 1.18.1, CuPy 14.2.0, cuFFT 12.4.0. Reused environment reports are explicitly
labeled; source state, upstream verification and host telemetry are fresh.

Numerical AOS source: `69f91214760ac91d9fc8b5a358ad6c1edd4d0fb6`; Copper
`781f41d38f784d0761e270b632d2c709c7234a76`; OOPAO
`168fad2a8d98721a270c1b066b813b716680b5ea`; SPECULA
`a06bb495ae29a6ec160628069598b32631b9f44f`. Upstream verification matches the
external pins. Isolated pinned checkouts are reused; cloning is outside timing.

Use this campaign's archived harness and fixture, a clean AOS numerical checkout
at `69f91214`, and the paired Julia Project/Manifest with local source paths
updated. The [harness guide](../../../matched_copper/README.md) documents a full
rerun. `collect.py` and `archive.py` correctly reject an incomplete standard
campaign. This stopped run instead preserves its explicit outside-harness
`collect-stopped-campaign.py` and `pack-stopped-campaign.py`, with partial status.
Their recorded local harness paths must also be updated when relocating evidence.

Six checksummed archives contain complete numerical raw/derived evidence,
completed timing samples, the interrupted record, fixture, environment locks,
commands, telemetry, source snapshots, negative evidence and independent review.
They exclude clones and the virtual environment. `raw-sha256.json` indexes
members; `archives-sha256.json` indexes archives. [summary.json](summary.json)
marks performance incomplete; [review.md](review.md) records independent review.

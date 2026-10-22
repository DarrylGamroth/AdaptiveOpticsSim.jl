# Fixed-sampling AOS CUDA modulation FFT batches — 2026-10-22

## Result and scientific qualification

At unchanged spatial sampling, modulation points, weights, precision and physical
parameters, the measured preparation choices reduce completed optical-node time
by **10.3–32.0%** on the changed cases. Both physical pupil-tilt propagation and
the then-existing shifted-mask path were measured. **2,760 focused checks passed**;
all **2,176** final timing allocation observations were zero host bytes.

**The external comparison uncovered a substantial pre-existing scientific
limitation.** In one matched case, shifted-mask propagation differs from OOPAO
by **17.4% in tip slopes and 19.0% in ramp slopes**. The tip discrepancy retains
an **8.2% residual after a fitted scalar gain**. This disqualifies that path as
an accuracy-focused reference; stable loops calibrated against the same plant
do not establish external physical fidelity. These errors are present with the
old and selected batch sizes.

The physical pupil-tilt path agrees closely with OOPAO, including the directly
tuned Float32 case. The batching change benefits that path independently.
After a bounded investigation established a structural fractional-shift
mismatch, the user requested withdrawal of the shifted-mask strategy. This
record preserves the before-removal measurements; the subsequent removal and
physical-model revalidation are recorded separately. Historical shifted-mask
results must not be used to qualify the replacement instrument graphs.

## AOS implementation boundary

The policy lives in AOS preparation and its CUDA extension. It uses FFT geometry
and point count, never instrument names. CUDA `ComplexF32`, sm86 and 5 MiB L2
identify the measured device class. Only these measured combinations change:

| FFT side | Points | Original batch | Selected batch |
|---|---|---|---|
| 256 | 32 | 32 | 8 |
| 272 | 24 / 30 / 32 | 24 / 30 / 32 | 8 / 6 / 8 |
| 408 | 32 | 32 | 2 |
| 512 | 12 / 32 / 96 | 12 / 32 / 32 | 2 |
| 528 | 64 | 32 | 2 |

Other shapes, counts, precisions and devices retain the existing policy. The
central helper still selects an exact divisor and applies the 256 MiB field
stack budget. Selection and device queries occur only during preparation;
captured execution binds the resulting workspace and FFT plans. No runtime
autotuning, environment variable, dependency or public batching API was added.
The benchmark-only override is explicitly isolated from production.

This is an empirical shared policy for the two measured propagation paths,
not a universal FFT optimum. Individual workloads can prefer a neighboring
batch. A broad cache-fitting rule was rejected: side 560 became approximately
9% slower with batch two. The final table excludes that configuration and
unmeasured count combinations. No hardware cache-counter explanation is claimed.

## Repeated completed-node timings

RTX 3080 10 GiB, 68 SMs, sm86, L2 5,242,880 bytes; driver 615.71.10,
Julia 1.13.2, CUDA/CUDACore 6.4.2, KA 0.9.45. One Julia/BLAS/FFT thread;
GPU jobs serialized, timing process on CPU 2. Other CPU activity and clocks
were uncontrolled. No profiler was active in the final timing runs.

Two fresh processes each cover 17 configurations and both strategies. Every
case has three rounds with randomized baseline/selected order, 16 warmups and
64 completed samples per block. Initial measurements are excluded from speedup.
The values below are medians of the process medians of round medians. Raw
samples, per-process dispersion, p95 and observed maxima are archived.

| FFT / points | Shifted-mask baseline → selected ms | Physical pupil tilt baseline → selected ms |
|---|---:|---:|
| 256 / 32 | 0.1792 → 0.1219 | 0.3439 → 0.2399 |
| 272 / 24 | 0.1631 → 0.1255 | 0.3030 → 0.2378 |
| 272 / 30 | 0.1968 → 0.1552 | 0.3709 → 0.3064 |
| 272 / 32 | 0.2069 → 0.1679 | 0.3916 → 0.3098 |
| 408 / 32 | 0.4668 → 0.4070 | 0.9302 → 0.7734 |
| 512 / 12 | 0.2685 → 0.2163 | 0.5141 → 0.3973 |
| 512 / 32 | 0.6375 → 0.5501 | 1.2839 → 0.9607 |
| 512 / 96 | 1.8356 → 1.5339 | 3.8237 → 2.7979 |
| 528 / 64 | 1.3168 → 1.1812 | 2.7198 → 2.1969 |

The boundary is one captured AOS photon-rate node through prepared completion;
it excludes detector acquisition, RTC, transport, host copies, setup and capture.
These are unpaced service samples, not a fixed-rate or worst-case latency test.
There are 26,112 repeated samples plus 8,704 initial samples. All 2,176 host
allocation observations are zero; maximum recorded output relative L2 change
is about 7.6 × 10⁻⁸. At side 512, field-stack storage drops from 64 MiB to 4 MiB.
FFT workspaces and physical modulation phases are additional storage.

## Fixed-sampling integration parity before withdrawal

COPPER 240/96 keeps all 277 actuators, five atmospheric layers, r₀=0.15 m,
seed 823, 320 frames and gain 0.15. Its baseline calibration and flat reference
are held fixed for the selected-batch loop, preventing recalibration from
hiding an implementation change.

- Thirteen photon-rate probes: relative L2 **8.389 × 10⁻⁸**; maximum relative
  probe-flux difference **2.344 × 10⁻⁹**.
- Full 277-column response: relative L2 **5.024 × 10⁻⁶**; worst individual
  column **1.158 × 10⁻⁴** (column 7).
- Loop commands: maximum difference **1.137 × 10⁻¹³ m OPD**; maximum residual
  WFE difference **2.752 × 10⁻¹⁵ m**. Tail WFE remains approximately **75.135 nm**.
- All loop outputs remain finite. Saved 32-frame noisy detector outputs match
  exactly, despite small photon-rate rounding differences.
- Completed captured noisy plant median: **1.918 → 1.683 ms** in this single
  before/after campaign. This excludes host RTC and transport and has weaker
  repetition than the AOS-only timing experiment above.

GPI 240/24 retains its batch 24. All **17 saved arrays** are byte-identical:
calibration, reconstructor, all used/next coefficient histories, spectral inputs
and final frames. The 320-frame moving runs exercise two tip/tilt, seven woofer
and 26 tweeter modes. Tail WFE remains 140.522 nm without noise and 140.573 nm
with noise; completed-cycle allocation observations remain zero. Unsaved
intermediate full-size camera images are outside the byte comparison.

These checks establish implementation parity, not physical sensor correctness.

## Independent OOPAO comparison

Fresh local OOPAO execution at
[`085d5e50`](https://github.com/cheritier/OOPAO/tree/085d5e50ace0d20fe13cc2da20129d5400166973)
uses NumPy CPU propagation. The recorded source files match that revision;
its import writes `precision_oopao.npy=64`, and the modified precision file,
source hashes and existing untracked checkout files are explicitly recorded.
No fit, transpose, registration or flux correction aligns the photon images.

All cases use D=8 m, an unobstructed pupil, I band 790 nm, 1 ms exposure,
centered PSF and explicit modulation counts:

| Pupil / pupil samples / points | Radius λ/D | Separation / edge / binning | Raw / detector side |
|---|---:|---|---|
| 24 / 4 / 8 | 1 | 4 / 0 / 2 | 72 / 6 |
| 32 / 8 / 20 | 3 | 8 / 0 / 1 | 96 / 24 |
| 120 / 30 / 32 | 3 | 2 / 1 / 1 | 256 / 64 |

The last case directly exercises physical Float32 CUDA batch **32 → 8**.
Flat/ramp image relative errors versus OOPAO are **8.05 × 10⁻⁶ / 7.93 × 10⁻⁶**;
the ±1 nm tip image derivative error is **7.56 × 10⁻⁵**. Changing batch alone
changes images by at most approximately **7.0 × 10⁻⁸**. The derivative amplifies
Float32 subtraction rounding; its batch difference is **7.11 × 10⁻⁵**.
Float64 CPU/CUDA physical propagation, including batches 1/2/4, agrees within
approximately 4.5 × 10⁻⁸ over the recorded images and derivatives.

The retired shifted-mask path retains **2.30% / 4.29%** flat/ramp image error
and **27.85%** tip-image derivative error. To test controller relevance, the
same OOPAO quadrant selection, valid mask and normalization were applied to
every saved frame. The normalization derivative is included analytically;
it reproduces OOPAO's saved ±1 nm slope response within **9.96 × 10⁻⁹**.

| Selected strategy | Ramp-slope error | Tip-slope response error | Projected tip gain | Tip residual after scalar gain |
|---|---:|---:|---:|---:|
| Physical pupil tilt | 0.000181% | 0.00526% | 0.999993 | 0.00522% |
| Shifted mask | 18.962% | 17.435% | 0.846231 | 8.217% |

These are specific finite-case response errors, not measured wavefront or
closed-loop errors for an instrument. They demonstrate that the discrepancy
reaches the sensor observable and cannot be dismissed as image-only detail.
Own-plant calibration can absorb some bias; that does not validate a real RTC.
No external complete-loop, detector-statistics or measured instrument claim is
made, and SPECULA was not used in this campaign.

### Why removal rather than a correction factor

The independent [discrete-model assessment](shifted-mask-assessment.md) compares
the source equations with the physical reference. Removing an irrelevant output
phase, exact pupil tilt requires S₋δ diag(M) Sδ, where Sδ = F Dδ F⁻¹. Integer
bin shifts make Sδ a permutation. Fractional shifts generally make it dense;
a pointwise translated mask cannot represent that operator.

A CPU discrimination experiment finds integer periodic shifts agree within
7.3 × 10⁻¹⁶, while fractional shifts differ by 19–40% in raw intensity on its
fixture. Correcting periodic wrapping does not fix the fractional case. This
supports withdrawal rather than an unsupported sign, gain or wrapping patch.

## Validation and reproduction

- CPU WFS selectors: **1,068 passed**.
- CUDA numerical/capture fixtures: **685 passed**.
- CUDA batch-policy checks: **378 passed**.
- CUDA modulation-storage/update regressions: **629 passed**.
- [Independent review](independent-review.md), [timings](timing-summary.json),
  [integration parity](integration-comparison.json), [external images](external-comparison.json)
  and [external slopes](slope-comparison.json) retain the evidence and limits.

`evidence.tar.gz` contains commands, environment manifests, logs, raw arrays,
source snapshots and the discrete-model experiment. Duplicate byte-identical
arrays are stored as tar hard links. `SHA256SUMS` verifies the report artifacts;
the extracted archive has its own file checksum index. The sweep's original
`native_nodes` field means logical AOS nodes (one), not recursive CUDA nodes;
the final harness corrects its label to `algorithm_nodes`.

The pre-removal source is AOS `e7b4fbb4` plus the archived
`source-before-removal/aos.patch` and benchmark sources. Each final timing file
pins its numerical source and harness hashes. Pilot/sweep/boundary harness
snapshots reproduce their recorded SHA256 values. Instrument sources are
COPPER `47e305f` and GPI `13b8c17`; application pins were not changed.

With isolated environments pointing at the appropriate recorded sources:

```sh
julia --project=CUDA_ENV benchmarks/sweep_modulation_batches_cuda.jl sweep.toml
AOS_BATCH_POLICY_ONLY=1 AOS_BATCH_COMPARE_POLICY=1 AOS_BATCH_ROUNDS=3 \
  AOS_BATCH_SAMPLES=64 julia --project=CUDA_ENV \
  benchmarks/sweep_modulation_batches_cuda.jl policy.toml CASES
python benchmarks/summarize_modulation_batches.py timing.json policy1.toml policy2.toml
julia --project=COPPER_ENV benchmarks/validate_modulation_batches_copper_cuda.jl COPPER_ROOT BEFORE
julia --project=COPPER_ENV benchmarks/validate_modulation_batches_copper_cuda.jl COPPER_ROOT AFTER BEFORE
python benchmarks/compare_modulation_batch_validation.py EVIDENCE integration.json
OOPAO_PATH=OOPAO_ROOT python benchmarks/oopao_modulation_batches/oracle.py REFERENCE
AOS_ORACLE_BACKEND=cuda AOS_ORACLE_PRECISION=Float32 \
  julia --project=CUDA_ENV benchmarks/oopao_modulation_batches/aos.jl REFERENCE
python benchmarks/oopao_modulation_batches/compare.py REFERENCE external.json
OOPAO_PATH=OOPAO_ROOT python benchmarks/oopao_modulation_batches/slope_comparison.py REFERENCE
```

Full commands and environment variables, including baseline/selected oracle
prefixes and the fresh GPI loop runs, are in the archive. Initial preparation
time, occupation of all SM resources, deployed deadlines and physical hardware
fidelity were not established by this experiment.

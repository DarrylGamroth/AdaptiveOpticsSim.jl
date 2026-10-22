# Retained physical pupil-tilt CUDA profile — 2026-10-22

## Scope and measured boundary

AOS `d781fb23` on an RTX 3080 10 GiB (68 SMs, sm86, 5 MiB L2), Julia
1.13.2, CUDA 6.4.2, KA 0.9.45; Nsight Systems 2026.3.2. No production
implementation or numerical sampling changed. New scripts are diagnostic
benchmark harnesses in an isolated worktree/environment.

These are **instrument-sized generic AOS optical nodes**, not full GPI/COPPER
simulations: unobstructed 8 m pupil, 750 nm source, 2×10⁸ photons m⁻² s⁻¹,
3 λ/D modulation, a deterministic sum of sinusoidal OPD terms scaled by 40 nm. Only pupil
resolution, modulation count and detector placement vary between named cases.
The boundary is one captured physical photon-rate node through completion;
no detector noise, RTC, transport, host copies, setup, compilation or capture.
The actual instrument models have additional physics and runtime costs.

One Julia/BLAS/FFT thread; CPU affinity 2; GPU jobs serialized. Each unprofiled
case has 16 warmups, 16 zero-host-allocation observations and 32 completed,
unpaced service samples. Separate fresh processes collect four warmed graph
replays under Nsight. Profiler durations are diagnostic, not production latency.
GPU clocks, power and host background activity were uncontrolled; 250 ms NVML
telemetry is included. No timing tail, fixed-rate or WCET claim is made.

## Results

| Generic case | Pupil / points | FFT side | Batch | Unprofiled median ms | FFT share of kernel time |
|---|---:|---:|---:|---:|---:|
| GPI default size | 240 / 24 | 960 | 24 | 3.595 | 64.2% |
| GPI fine size | 720 / 192 | 2880 | 4 | 318.961 | 72.5% |
| COPPER fast size | 240 / 96 | 512 | 2 | 2.948 | 50.5% |
| COPPER default size | 480 / 32 | 1024 | 32 | 5.134 | 62.4% |
| COPPER fine size | 1440 / 256 | 3072 | 2 | 456.524 | 66.6% |

| Fine case | Pupil field construction | Mask multiplication | Intensity accumulation | Final binning |
|---|---:|---:|---:|---:|
| GPI size | 7.0% | 13.4% | 7.0% | 0.1% |
| COPPER size | 9.2% | 15.2% | 8.9% | 0.1% |

The fine traced replay spans contain GPU kernel/memset activity for 99.83%
and 99.65% of elapsed device span respectively. Default-size spans are
97.64% and 98.92%; the smaller COPPER-fast case is 86.42%. These are timeline
activity fractions, **not achieved SM occupancy, memory-bandwidth utilization,
or percentage of peak useful throughput**. Host starvation is not the main
observed fine-case problem.

All five hot traces contain four output-clearing memsets and no memcpy
activity. Kernel counts match the prepared per-batch sequence: four FFT kernels,
field fill, mask, intensity, followed by one final binning kernel per replay.
Every distinct kernel graph-node ID occurs exactly four times: 8 nodes for each
default case, 337 for GPI fine/COPPER fast, and 897 for COPPER fine.

Nsight emits a generic “Not all CUDA events might have been collected” diagnostic.
The expected hot kernel/node/replay counts and memsets are complete in these
exports; the warning is retained, and completeness of other CUDA event classes
is not asserted. Software instrumentation was used; hardware CUDA tracing is
unsupported here. Expected profiler-injected CUDA-library warnings are retained.

Warmed host allocations are zero in all 80 observations. Reported used device
memory stays unchanged during samples. Deterministic photon maps are byte-equal
before and after timing/profile replay. The source and physical numerical model
are unchanged from the preceding [external-reference and capture validation](../2026-10-22-pupil-tilt-reference/README.md).
This profiling campaign does not add a new external instrument-fidelity claim.

## What FFT batching changed

For each modulation point p, the physical sequence remains:

```text
Eₚ = amplitude × exp(i 2π OPD/λ) × tiltₚ
I  = Σₚ wₚ | inverse_FFT(mask × FFT(Eₚ)) |²
```

Centering phasors and normalization are retained in the implementation. The
field stack is `(FFT side, FFT side, B)`; cuFFT transforms dimensions `(1,2)`
for B independent modulation points at once. Batching changes grouping and
workspace size, not point count, quadrature, precision, optical propagation
or the number of mathematical transforms. Every point still gets one forward
and one backward 2D FFT. The batch is selected during preparation and captured.

The previous policy chose the largest exact divisor at most 32 within a
256 MiB field-stack budget. The measured sm86/5 MiB-L2 table selects:

- FFT side 256, 32 points: batch 8.
- Side 272, 24/30/32 points: batch 8/6/8.
- Side 408, 32 points: batch 2.
- Side 512, 12/32/96 points: batch 2.
- Side 528, 64 points: batch 2.

Other geometries/precisions/devices retain the prior policy. The
[repeated sweep](../2026-10-22-modulation-batches/README.md) measured physical
node medians of 0.3439→0.2399 ms (256/32), 0.3916→0.3098 ms (272/32),
and 3.8237→2.7979 ms (512/96). The 512 field stack shrinks from 64 MiB to
4 MiB. Smaller batches mean more launches but can improve the cuFFT workload
and working set sufficiently to outweigh that cost. These are empirical results;
cache/algorithm causality has not been established by hardware counters.
A broader cache-size rule regressed at side 560 and was rejected.

The default 960/24 and 1024/32 cases do not change. The finest cases were
already constrained by the scratch budget to four and two. Their stacks use
253.125 MiB and 144 MiB respectively; the newer table does not select those
values specially. A fresh fine-case divisor sweep is recorded separately below.

## Fine batch sweep

The unchanged 256 MiB field-stack budget permits GPI-sized batches 1/2/3/4
and COPPER-sized batches 1/2. Two randomized rounds gave these medians of
round medians (16 completed samples per trial):

| Case | Batch | Median ms |
|---|---:|---:|
| GPI fine | 1 | 352.011 |
| GPI fine | 2 | 330.540 |
| GPI fine | 3 | 322.407 |
| GPI fine | **4 (current)** | **320.173** |
| COPPER fine | 1 | 494.314 |
| COPPER fine | **2 (current)** | **456.241** |

The current choices are the fastest observed within this budget. GPI batch
three is only about 0.7% slower, so the sweep does not establish a universally
unique optimum. Maximum photon-map relative L2 across changed choices is
1.06 × 10⁻⁷; every observed host allocation is zero. No policy change is
supported by this sweep. Larger scratch budgets and changed representations
were not tested. See `fine-batches.toml` in the archive and the summary JSON.

## Counter limitation

Nsight Compute 2026.3.1 was attempted on the pupil-fill kernel with clocks and
cache controls disabled. The driver returned `ERR_NVGPUCTRPERM`. No policy,
capability, module setting, clocks or power limit was changed. No achieved
occupancy, cache-hit, DRAM-throughput or stall-counter conclusion is possible.
The counter-attempt process timing is instrumented and excluded from all
baseline and batch results. Only trace launch-resource metadata is interpreted.

## Further AOS improvements, in priority order

### 1. Reduce repeated mask and pupil work in the existing physical kernels

Mask multiplication costs 13–15% at fine sampling and 18–21% in the smaller
cases. Try one thread per spatial pixel, retaining a mask value in a register
and looping over the batch planes. Each field value receives the same single
complex multiply, without reordering a reduction or changing propagation.
This is a lower-complexity first experiment than introducing cuFFT callbacks.
The existing cache may already share some mask loads, so a speedup is unproven.

Pupil construction reevaluates `cispi(2 OPD/λ)` for every modulation point even
though OPD and wavelength are fixed during a frame. Two candidate designs are:

- Loop over batch planes per pixel, sharing OPD phase, amplitude and centering
  within B points, without another scratch array.
- Compute the OPD phasor once per frame in a retained N² workspace, then read it
  while filling every point. Caching just the phasor can retain subsequent
  multiplication order; caching the complete base field reorders products.

The full pupil-fill stage is only 7–9% of fine kernel time, including unavoidable
padded field writes. Its complete elimination would be an optimistic upper bound,
not an achievable speedup forecast. Measure the extra load/store and register
cost of reuse before selecting either design.

### 2. Compact physical modulation phases, with unchanged samples

Dense physical phases are generated as `cis(uₚ c[j] + vₚ c[i])`. They factor
algebraically into `cis(vₚ c[i]) × cis(uₚ c[j])`. Two N×P complex factor arrays
would replace N²P phases while still tilting the pupil before each FFT.

| Fine case | Dense device phases | Separable device factors |
|---|---:|---:|
| N720 / P192 | 759.375 MiB | 2.109375 MiB |
| N1440 / P256 | 3.955 GiB | 5.625 MiB |

The memory reduction is derived directly from element counts. Runtime benefit
is a hypothesis: factor reads and an extra complex multiply replace the dense
phase read. This is **not the retired shifted-mask approximation**; the forward
FFT per tilted pupil remains. Floating-point factorization changes rounding,
so no accuracy claim is made before derivative/calibration/oracle validation.
Preserve the existing host-rounded coordinates, sample offsets, weights,
live circular-radius updates, exact device binding and graph replay.

### 3. Specialize launch geometry only where measurements justify it

AOS uses spatial 8×8 groups for these 3D kernels, with up to four modulation
planes; intensity uses 16×16. Controlled candidates are `(32,2,z)` versus
`(8,8,z)` for fill/mask, and `(32,8)` versus `(16,16)` for intensity, at the
same thread counts. Test one kernel at a time under the normal FFT cache history.

Trace resources show fill 26, mask 20 and intensity 40 registers per thread,
with no local/shared memory for those kernels. Their existing accesses can
already fill aligned 32-byte sectors; wider rows are not proof of improved
coalescing. Registers and launch sizes do not measure achieved occupancy.

### 4. Fuse intensity accumulation with detector integration where valid

The no-LGS photon-rate path accumulates full padded intensity after each batch,
then bins it at the end. Incoherent modulation addition and spatial integration
commute mathematically. A coalesced staged reduction into detector bins could
reduce repeated native-intensity traffic. It must still read every propagated
complex sample and execute every FFT. Preserve native-intensity APIs and the
LGS path, which use that intermediate. Reordered sums need numerical validation.
Final binning alone is only about 0.1% at fine sampling; removing that final
launch by itself is not a useful priority.

### 5. Investigate FFT integration after the simpler passes

FFTs remain 67–73% of fine kernel time. Layout is already contiguous by spatial
plane with independent batched 2D transforms; no obvious layout defect was found.
CUDA-specific callback/fusion or alternate plan/workspace experiments may reduce
mask/FFT memory passes, but require feasibility and capture qualification before
production changes. Reducing padding to a convenient FFT size changes the
sampled physical propagation contract and cannot be assumed accuracy preserving.

With FFT time unchanged, even eliminating every surrounding kernel gives only
an approximate 1.38×/1.50× fine-case kernel-time ceiling. This derived Amdahl bound
helps scope expectations; actual smaller improvements must be measured.

## Accuracy and promotion gates

Keep resolution, padding, modulation nodes/weights, OPD, precision and baseline
calibration fixed. Check rates, total photons, flat/ramp/tip derivatives, slopes,
full interaction columns and a frozen-calibration closed loop. Reuse the matched
OOPAO oracle and independent analytic FFT sequence. Check noncircular sampled
paths, nonuniform weights, pupil shifts/rotation/rooftop and live radius changes.
Require captured replay, native allocation-node rejection and zero observed
steady host allocation. A small rate-map difference alone does not qualify a
change because finite-difference calibration can amplify rounding.

## Artifacts and reproduction

The summary and raw baseline timings distinguish unprofiled service time from
profiled activity. The archive contains five `.nsys-rep` files, text kernel/API
summaries, commands, logs, source harness snapshots, environment manifests,
telemetry and the fine batch sweep. SQLite exports can be regenerated with
`nsys stats`; they are omitted as redundant. Raw file hashes are indexed and
verified after extraction.

Use `benchmarks/profile_physical_pyramid_cuda.jl OUTPUT.toml CASE [BATCH]` in
the isolated environment; `BATCH=-1` selects production behavior. Set
`AOS_PROFILE_EXTERNAL=1` only for an external profiler. `run_profiles.py` records
the exact serial commands; adjust checkout paths in the included Manifest.
`analyze_profiles.py` reconstructs the category summary from SQLite. The fine
sweep uses two randomized rounds of all divisors fitting the unchanged budget,
eight warmups and 16 completed samples per trial. No production policy is
changed solely on these small diagnostic samples.

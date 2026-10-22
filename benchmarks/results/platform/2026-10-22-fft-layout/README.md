# Fine physical Pyramid FFT layout and batching — 2026-10-22

## Scope and decision

This investigation retains physical pupil tilts, every modulation point, weights,
sampling, padding, precision, source prescriptions and controller gains. It changes
implementation of the same two-dimensional Fourier propagation, not the optical
model. Baseline production source is `80cc4fef`.

On the measured CUDA sm86 / 5 MiB L2 device class, `ComplexF32` square planes use:

| FFT side / modulation points | Previous batch | Selected batch | Selected FFT layout |
|---|---:|---:|---|
| 2880 / 192 | 4 | 16 | Contiguous axis transforms, retaining transposed frequency data |
| 3072 / 256 | 2 | 16 | Native two-dimensional cuFFT |

Other sizes, point counts, precisions and device classes retain their previous
selection. These are backend preparation choices; instrument models and public
scientific APIs do not change. The layout family remains fixed when its scratch
is resized, matching the existing concrete workspace and fixed tile-size contract.
Fresh preparation applies device qualification again.

The tested production patch SHA-256 is:

```text
8ef1c3abaf072b243926a2328728c8b8462438987d5206a136a22c8c4be1fed2
```

## Sampling scale

The hundreds-of-milliseconds figures below are **fine-reference** workloads.
They do not describe the 240/96 fast mode previously measured at approximately
2.87 ms for the captured photon-rate node, excluding detector and RTC. That
preceding campaign also measured the actual captured fast moving-atmosphere
plant at 2.801 ms, excluding the host RTC. The fast case uses 512² FFT planes;
the 1440/256 reference uses 3072² planes and
processes **96 times as many complex spatial samples per frame** before the
additional FFT-size factor. This campaign changes neither instrument sampling
nor its fast/accurate mode selection. The earlier 480/32 default optical-node
measurement was approximately 5.08 ms. These historical values come from the
[preceding pupil-reuse campaign](../2026-10-22-pupil-reuse/README.md).

The final-code sanity check measures **2.858 ms** at 240/96 and **4.993 ms**
at 480/32 (24 warmups, 32 completed captured photon-rate samples each). Both
saved rate maps are byte-identical to the preceding campaign and all sixteen
warmed host-allocation observations per case are zero. These are checks of
unchanged paths, not fresh interleaved performance comparisons.

## Fixed-sampling experiments

The generic optical fixtures correspond to fine instrument-sized grids:
GPI-sized N720/P192 and COPPER-sized N1440/P256. They use an unobstructed 8 m
pupil, 750 nm wavelength, 3 λ/D modulation and the same deterministic OPD
sinusoids scaled by 40 nm as the preceding physical profiling campaign.
They measure one captured photon-rate node through completion, excluding
setup, transfers, detector and RTC. They are not full instrument runtime figures.

Two randomized rounds, sixteen warmups and twenty-four completed service samples
per trial compare larger batches under an experimental **1.5 GiB main-stack**
bound. The following are medians of the two round medians, milliseconds:

| Implementation / batch | 2880² / 192 | 3072² / 256 |
|---|---:|---:|
| Native / old production batch | 317.092 (B4) | 452.714 (B2) |
| Native / B8 | 313.447 | 435.214 |
| Native / B16 | 311.818 | 433.182 |
| Separate strided axes / B16 | 316.052 | 434.789 |
| Generic transpose for each transform / B16 | 656.443 | 1014.947 |
| Tiled transpose for each transform / B16 | 381.709 | 583.470 |
| Retain tiled transposed frequency data / B16 | **306.546** | 468.900 |

All saved rate maps are compared against native cuFFT at identical physical
sampling, with hashes checked. Maximum relative L2 across these experiments is
below **1.78 × 10⁻⁷**, within the unchanged 2 × 10⁻⁶ gate. Every sampled warmed
host allocation is zero. Replay equality within a run and cross-layout numerical
agreement are recorded as separate checks.

The independent-axis experiment is not equivalently batched cuFFT: installed
cuFFT.jl chooses spatial dimension 1 as the internal batch for an axis-2 plan,
then loops over modulation dimension 3. It makes 2 + 2B cuFFT calls per tile,
versus native's 2. Capturing removes repeated host submission but retains the
additional device kernels. The transpose and paired variants make four axis-1
cuFFT calls per tile with contiguous internal batches.

GPUArrays' generic `permutedims!` reads linearly and scatters strided stores.
The custom transpose uses a padded `(33,32)` shared tile in Julia's column-major
layout, workgroups `(32,8,1)`, and coalesced reads and writes. The padding is on
the leading dimension. Retaining transposed frequency data removes two of the
four explicit transposes needed by independent forward/inverse transpose plans.

## Production timing and CUDA profiles

Fresh original/modified process pairs repeat the win. The second pair reverses
execution order. Each run uses 24 warmups, 16 allocation observations and 32
completed unprofiled service samples:

| Fine case | Original round medians (ms) | Modified round medians (ms) | Reduction of median |
|---|---:|---:|---:|
| 2880 / 192 | 316.992 / 316.778 | 306.556 / 306.531 | **3.26%** |
| 3072 / 256 | 451.991 / 452.200 | 433.536 / 433.010 | **4.16%** |

Fixed-input rate relative L2 is **1.753 × 10⁻⁷ / 8.628 × 10⁻⁸**, with flux
relative errors **8.218 × 10⁻⁸ / 5.083 × 10⁻¹⁰**. Outputs are byte-identical
between independent runs of each implementation. Warmed host allocations are
zero, and the CUDA pool observation remains unchanged during timed replay.

Each Nsight Systems region contains four completed graph replays. Per-frame
kernel sums in milliseconds are diagnostic, not the unprofiled timing above:

| Stage | GPI-sized original | GPI-sized modified | COPPER-sized original | COPPER-sized modified |
|---|---:|---:|---:|---:|
| Pupil fill | 22.092 | 23.984 | 37.846 | 39.009 |
| FFT kernels | 232.643 | 148.301 | 306.824 | 310.311 |
| Explicit transposes | 0 | 75.512 | 0 | 0 |
| Mask | 41.969 | 39.861 | 70.397 | 60.066 |
| Intensity accumulation | 22.700 | 19.224 | 41.311 | 29.149 |

GPI-sized FFT kernels plus their new transposes take 223.813 ms, versus
232.643 ms previously. The COPPER-sized FFT kernels themselves are slightly
slower; the mask and intensity stages save about **22.5 ms** together. The
whole-node improvement therefore does not mean every FFT got faster.

Distinct captured kernel nodes fall **338 → 110** and **898 → 114**. Every
node appears four times in its trace. All four regions record four memsets and
no host/device copies; retained trace warnings do not justify claiming every
possible CUDA event was collected. Kernel busy time covers over 99.7% of the
recorded span. This is activity evidence, not achieved occupancy or efficiency.

## What this says about bandwidth

The larger native batches reduce repeated mask and intensity-array accesses.
They retain the same FFT sample volume and do not reduce modulation count.
For side N, point count P and tile B, source-level logical traffic for mask reads
and accumulated-intensity reads/writes is approximately **16 N²P/B bytes**:
8 bytes per complex mask entry and 8 bytes per Float32 intensity read/write.
Complex field reads/writes, FFT internals and other kernels are additional.

| Fine case | Logical mask/intensity traffic saved per frame by B16 |
|---|---:|
| 2880 / 192, previously B4 | 4.778 GB (4.449 GiB) |
| 3072 / 256, previously B2 | 16.911 GB (15.750 GiB) |

These are **analytical array-access counts, not measured DRAM transactions**.
Cache hits, compiler behavior and opaque cuFFT traffic are not resolved by this
count. The paired layout itself adds two explicit transposes: four complex
array read/write volumes, about **50.961 GB** per GPI-sized frame. It still wins
on the measured 2880 grid. Contiguous access and a different cuFFT algorithm
repaying the transpose cost is a hypothesis; timelines alone do not separate
memory transactions, instruction throughput and library algorithm effects.

The results therefore do not establish a universal lower-bandwidth layout.
The paired path loses on 3072 and is not enabled there. Its qualification is
specific to the measured device class; no newer midrange GPU was available
for comparison. A target with less bandwidth needs its own sweep.

The driver retains restricted NVIDIA performance-counter access. The prior
permission failure and current read-only policy are archived. No permissions,
clocks or power limits were changed. GPU memory utilization percentages are
activity indicators, not fractions of peak bandwidth. No achieved-bandwidth,
occupancy or roofline classification is claimed.

A final shared-memory mask microprobe coalesces reads of the original mask
while retaining transposed frequency data. Two rounds give 3.264 / 3.261 ms,
versus 3.328 / 3.327 ms for the selected mask kernel, with bit-identical outputs
at both precisions and partial tiles. This is about 0.8 ms per 192-point frame
if all twelve mask invocations benefit equally, before whole-path effects.
That extrapolation is not a measured pipeline gain. The candidate remains
benchmark-only; its isolated saving does not establish a repeatable whole-node
improvement.

Potential further traffic reduction: cuFFT load/store callbacks could combine
mask application with a transform boundary and eliminate a standalone field
pass. Local bindings expose callback functions, but the package has no prepared
high-level callback lifecycle, and native handles are pooled. This requires a
dedicated retained owner, numerical tests and graph-capture qualification; it
was not introduced as an unmeasured production change. Existing centering is
already fused into pupil fill, and inverse normalization into intensity.

## Workspace and ownership

The 1536 MiB bound applies to the **main field stack**, not total device memory.

| Fine case | Previous main stack | Selected main stack | Additional transpose stack |
|---|---:|---:|---:|
| 2880 / 192 | 253.125 MiB | 1012.5 MiB | 1012.5 MiB |
| 3072 / 256 | 144 MiB | 1152 MiB | 0 |

FFT library work areas, propagation fields/masks, phase factors/cache, optical
products and instrument state are additional. Pool snapshots and sampled GPU
memory telemetry are not presented as exact total VRAM peaks.

Transposed scratch and axis FFT handles have a concrete workspace owner,
separate from the immutable physical/numerical plan. Exact-device validation
and prepared bindings include the auxiliary array and handles. Resize rebuilds
buffers/plans in the existing family. Eager propagation synchronizes its phases;
captured propagation enqueues on the existing stream. No runtime device query,
autotuning, dense phase cube, host transfer or public environment knob is added.

## Validation

Completed software checks on the frozen production patch:

- Focused CPU WFS/graph suites: **1,597**; KA CPU kernel matrix: **1,038**.
- CUDA production policy and neighboring configurations: **676**.
- Forced paired layout on small asymmetric fixtures: **964**, including both
  precisions, unbinned intensity/rates, zero/disabled modulation, live updates,
  complex-field/uncentered paths, exact target, pointer retention and resize.
- Native/fallback storage regression: **956**; modulation batch/capture: **457**.
- CUDA captured storage/handle restoration: **33**; complete GPI capture: **350**.
- Benchmark tiled transpose: **24**; historical baseline-budget regression: **6**.

### Fine Copper calibration and closed loop

Both versions use the actual N1440/P256 physical model: all 277 calibration
columns, 1,432 slopes, 20 nm pokes, thirteen deterministic probes and thirty-two
seeded atmospheric replay frames. Captured and stream execution agree. A frozen
baseline calibration drives both 320-frame loops: gain 0.15, r₀ = 0.15 m,
seed 823 and 231 retained modes. Controller and scientific inputs do not change.
Science diagnostics use the same 480-grid evaluation in both runs.

| Quantity | Modified versus original |
|---|---:|
| Photon-rate relative L2 | 9.129 × 10⁻⁸ |
| Maximum probe photon-flux relative error | 6.253 × 10⁻⁹ |
| Slopes relative L2 | 8.078 × 10⁻⁷ |
| Calibration matrix relative L2 | 5.390 × 10⁻⁶ |
| Maximum calibration-column relative L2 | 1.208 × 10⁻⁴ |
| Atmospheric photon-rate relative L2 | 8.713 × 10⁻⁸ |
| Closed-loop command relative L2 | 1.466 × 10⁻⁷ |
| Maximum residual-WFE difference | 3.262 × 10⁻¹⁵ m |
| Maximum Strehl difference | 3.152 × 10⁻⁸ |

The seeded atmospheric frames are byte-identical. All existing numerical gates
pass. The baseline loop's tail residual is 75.084 nm, Strehl 0.48277; every
sample remains finite. This establishes implementation parity at this sampling,
not convergence of the underlying atmospheric or optical model.

The actual instrument's captured moving-atmosphere plant median changes
**452.532 → 433.666 ms** (128 completed samples after 32 warmups); captured
calibration-plant medians are **451.289 → 434.594 ms**. These include graph
completion and exclude the host RTC and science diagnostics. Every observed
warmed captured host allocation is zero. Unlike the generic fixed-input pairs,
these instrument measurements are one process per version and not interleaved.

### GPI complete loop

The actual prepared fine owner selects B16 and the transposed workspace;
the default retains its native workspace and B24. The maintained GPI exercise
uses 2 tip/tilt, 7 woofer and 26 tweeter modes with 2,403 physical commands,
128 static frames and 320 moving frames per run. It captures the complete
plant and modal RTC, including next-command synthesis. No gains, calibration
cutoffs, seed, limits, detector prescription or physical sampling are changed.

Default calibration, controller traces and saved detector frames are
byte-identical to the baseline. At fine N720/P192 sampling:

| Quantity | Modified versus original |
|---|---:|
| Calibration matrix relative L2 | 7.552 × 10⁻⁷ |
| Maximum calibration-column relative L2 | 3.236 × 10⁻⁶ |
| Noiseless moving maximum wavefront-vector difference | 0.00009476 nm |
| Noiseless moving maximum scalar residual-WFE difference | 0.00002354 nm |
| Calibration-transfer maximum wavefront-vector difference | 0.00003501 nm |
| Noisy moving maximum wavefront-vector difference | 0.09949 nm |
| Noisy moving RMS wavefront-vector difference | 0.04186 nm |
| Noisy moving maximum scalar residual-WFE difference | 0.02045 nm |

The vector metric uses the saved nonorthogonal physical-basis Gram matrix.
All deterministic vector differences and all scalar residual-WFE differences
pass the existing 0.05 nm gates. The fine noisy vector's **provisional**
0.05 nm check is **false** and remains recorded as such. Integer stochastic
detection can change a seeded trajectory when tiny photon-rate differences
cross a sampling threshold. This one realization has a 6.724 nm RMS realized
noise contribution; the implementation difference is 0.04186 nm RMS. These
observations do not establish statistical equivalence. The comparator preserves
its previous distinction between asserted deterministic/scalar gates and the
reported provisional noisy-vector check; no tolerance has been relaxed.

All runs are finite and unsaturated. Own-calibration linear spectral radii are
approximately 0.650000006; calibration-transfer is 0.650292776. Fine moving
residuals are 139.630 nm noiseless / 139.710 nm noisy, versus approximately
405 nm open-loop and 123 nm fitting floor. These are results of the illustrative
35-mode exercise, not full-order instrument/controller qualification.

Three warmed uninterrupted trials per moving run give these **current** medians
of trial medians, milliseconds:

| GPI case | Noiseless | Noisy |
|---|---:|---:|
| Default N240/P24 | 3.430 | 3.448 |
| Fine N720/P192 | 307.014 | 307.049 |
| Default calibration transferred to fine plant | 307.002 | — |

Each trial has eight warmups and 64 completed cycles, including device input
staging, plant, RTC and command synthesis, with diagnostics outside the timer.
All **80** warmed moving-cycle host-allocation observations are zero. Native
capture also retains its recursive allocation/free-node rejection. These are
unpaced service observations, not deadline or fixed-arrival guarantees.

An independent CPU auditor verifies forty saved array hashes, source hashes,
actuation-basis reconstruction, Gram matrix, calibration SVD/reconstructor,
controller updates, residual-WFE metrics and timing reductions. The Julia basis
oracle hash matches exactly; independent Python reconstruction differs by at
most 7.57 × 10⁻¹⁷ m. The auditor is a mathematical/internal consistency check;
this campaign does not add a new external OOPAO or SPECULA comparison.

## Environment, provenance and reproduction

RTX 3080 10 GiB, sm86, 5 MiB L2; driver 615.71.10; Julia 1.13.2,
CUDA 6.4.2 and KernelAbstractions 0.9.45. One Julia/BLAS/FFT thread; GPU jobs
are serialized and pinned to CPU 2, CPU verification to CPU 4. The host is
shared and clocks/power are uncontrolled. Timings are unpaced single-outstanding
completed service observations, not fixed-rate latency, WCET or deadline claims.
Profiled timing is diagnostic. The 0.5 ms model step is not attained.

Harness entry points:

- `sweep_fft_layout_cuda.jl OUTPUT.toml CASE`: diagnostic layout/batch sweep.
  `AOS_FFT_LAYOUT=native|axis|transpose`, `AOS_FFT_TRANSPOSE=generic|tiled`,
  `AOS_FFT_PAIR=1`, and `AOS_FFT_BATCHES` affect benchmarks only.
- `profile_physical_pyramid_cuda.jl OUTPUT.toml CASE`: actual production policy;
  `AOS_PROFILE_EXTERNAL=1` delimits warmed profiler replays.
- `validate_fine_fft_copper_cuda.jl COPPER_ROOT OUTPUT [FROZEN_BASELINE]`:
  actual N1440/P256 calibration, captured/stream replay and 320-frame loop.
- `validate_fine_fft_gpi_cuda.jl GPI_ROOT OUTPUT`: actual-owner preflight and
  GPI's maintained default/fine/transfer complete-loop exercise.
- `summarize_fft_layout.py EVIDENCE OUTPUT.json`: hash-checked cross-layout rates.
- `compare_fft_validation.py PREVIOUS_REUSE_EVIDENCE FFT_EVIDENCE OUTPUT.json`:
  unchanged numerical gates, including fine Copper and GPI nonorthogonal
  wavefront-vector metrics. Historical comparator refactoring reproduces all
  previous results exactly.

Environments are isolated; no application pins or dependencies are changed.
The GPI numerical baseline is the retained preceding pupil-reuse campaign.
Its recorded dirty parent is `27179140`; its archived production patch is exactly
`git diff 27179140 80cc4fef -- src ext`, and every archived modified production
file matches `80cc4fef` byte-for-byte. The identity check is archived. Thus this
comparison uses the same baseline implementation as the new fine Copper runs.
Full-loop service observations from that earlier campaign are not treated as
fresh interleaved timing evidence.

Raw commands, source snapshots, production patch, manifests, logs and arrays
are archived. Early probes preceded additional cold metadata reporting; later
production validation uses the frozen patch above. Initial copied versions are
kept separately from final source. An aborted duplicate follow-up launch has no
accepted timings; its note/logs are retained. Test and GPI-preflight namespace, keyword-argument and node-name errors were
corrected before successful execution; their failed-attempt logs are retained. The original fine Copper cold driver emitted
world-age warnings after dynamically including helper modules; the later driver
invokes the complete cold call in the latest world. The baseline finished with
finite results, and its original source/logs are retained.


Delivered artifacts:

- `metadata-evidence.tar.gz`: commands, environments, source snapshots,
  production patch, arrays from layout sweeps, logs and reduced records.
- `gpi-evidence.tar.gz`: actual default, fine and calibration-transfer arrays.
- `copper-evidence.tar.gz`: fine original/modified calibration and loop arrays.
- `profiling-evidence.tar.gz`: four Nsight Systems reports; SQLite can be
  regenerated with `nsys stats` and is omitted.
- `raw-sha256.json`: 336 per-file archive hashes, verified by reading each archive.
- `layout-summary.json`, `production-summary.json`, `integration-comparison.json`,
  `scientific-replay.json`, `normal-sampling-check.json`: reduced evidence.
- [Independent final review](review.md): stable findings and dispositions.
- `SHA256SUMS`: hashes of this report and all delivered artifacts.

The original GPI baseline arrays are already retained in the linked preceding
pupil-reuse campaign and are not duplicated here. Reproduction requires those
baseline artifacts; workspace paths in the isolated manifests must be adjusted.

# Physical pupil reuse and modulation storage — 2026-10-22

## Changes

The retained physical Pyramid propagation now:

1. Loads one mask value per spatial pixel and reuses it across batch planes.
2. Computes `cispi(2 OPD/λ)` once per pupil propagation into an N² device
   workspace, then reuses it across every modulation point.
3. Stores each physical tilt as two N×P complex factor arrays:
   `cis(vₚ c[i]) × cis(uₚ c[j])`, replacing the N²P complex phase cube.

Every modulation point retains its weight, physical pupil-plane tilt, forward
FFT, fixed pyramid mask and inverse FFT. Sampling, padding, precision, gains,
source prescriptions and instrument models are unchanged. This is an
implementation optimization; it does not restore the retired shifted-mask
approximation or introduce an approximate reduced plant.

CPU and Bi-O-edge dense modulation APIs remain intact. Accelerator Pyramid
preparation uses an internal separable owner. Both factor arrays and the new
pupil-phase scratch receive exact-device validation; the cache participates in
prepared workspace binding. Circular modulation updates write into the same
factor arrays, preserving point count, weights, phase origin and captured
pointers. Noncentered, complex-field and synchronous paths remain covered.

The preceding [profiling harness and results](../2026-10-22-physical-pupil-profile/README.md)
were committed first as `27179140`. Baseline production source is `d781fb23`;
`27179140` changes only profiling artifacts. The tested/reviewed production
patch SHA-256 is:

```text
37418e59bb2dea6443eec09b145c66bd5f0aa61a29c2eecbcf96f1fa1c1f9cc5
```

## Memory

Direct storage follows the actual array shapes, not allocator reservation:

| Generic physical case | Previous dense phases | New tilt factors | Added common pupil-phase scratch |
|---|---:|---:|---:|
| 240 / 96 | 42.1875 MiB | 0.3515625 MiB | 0.439453 MiB |
| 720 / 192 | 759.375 MiB | 2.109375 MiB | 3.955078 MiB |
| 1440 / 256 | 3.955078 GiB | 5.625 MiB | 15.820313 MiB |

The largest case replaces about 4 GiB of phases with about 21.45 MiB of factors
and cache combined. FFT stacks, plans, masks and other storage are additional.
All arrays remain resident on the device. No dense host phase cube is built.
Observed CUDA pool snapshots are included separately; they can include
uncollected cold temporaries and library storage, so pool differences are not
presented as exact retained-owner accounting or total VRAM peak.

## Completed-node performance

Median completed-node latency, milliseconds; these are unpaced, generic optical
fixtures at fixed sampling, not complete instrument simulations:

| Case | Original | Reuse implementation | Second original / reuse pair |
|---|---:|---:|---:|
| COPPER fast | 2.942 | 2.868 | — |
| GPI default | 3.465 | 3.489 | — |
| COPPER default | 5.163 | 5.080 | — |
| GPI fine | 319.033 | 316.340 | 319.103 / 316.422 |
| COPPER fine | 456.587 | 451.459 | 456.548 / 452.236 |

The fine reductions repeat at approximately 0.84% and 0.94–1.12%. Small-case
changes are not established speedups: short warmup measurements drifted.
With 1,024 warmups, the fast storage diagnostics measured dense **2.683 ms**,
separable **2.684 ms**, and computed **2.715 ms**. Each is a single fresh-process
observation; dense versus separable is indistinguishable at this resolution.
The major demonstrated benefit is reduced storage, with modest fine-case speed.

A new Nsight Systems capture records four fine COPPER replays: **3,592 kernel
instances**, four memsets and no host/device copies recorded in the captured region.
The extra pupil-phase kernel runs once per replay (0.040 ms mean). Per-replay
kernel sums are FFT **307.289 ms**, mask **70.391 ms**, tilted-pupil fill
**38.076 ms**, intensity **41.313 ms**, and binning **0.571 ms**. FFTs retain
about **67.1%** of kernel time. Original profiled fill was 42.771 ms; this
supports concentrating saved work in pupil preparation, but profiled stage
comparisons are diagnostic and cannot isolate each of the three changes.
The trace warns that some CUDA events may not be collected; kernel counts match
898 distinct graph kernel nodes, four instances each. No hardware DRAM bandwidth
or achieved occupancy is inferred from these timings.

## Deterministic numerical and calibration parity

The COPPER integration compares the prior selected batch-two physical model
with this implementation at N240/M96. It retains all 277 response columns,
13 rate probes, seeded five-layer atmosphere, 320 frames, gain 0.15 and the
same frozen baseline calibration in both loops. Existing gates are unchanged.

- Photon-rate relative L2: **7.334 × 10⁻⁸**.
- Slopes relative L2: **6.359 × 10⁻⁷**.
- Full response matrix relative L2: **4.253 × 10⁻⁶**.
- Worst response-column relative L2: **9.872 × 10⁻⁵**.
- Frozen-loop command relative L2: **1.473 × 10⁻⁷**, max **1.137 × 10⁻¹³ m**.
- Maximum residual-WFE difference: **3.155 × 10⁻¹⁵ m**.
- Recorded seeded atmospheric detector frames are identical.

GPI's default and fine response matrices differ by **5.328 × 10⁻⁷** and
**4.587 × 10⁻⁷** relative L2 respectively; worst columns differ by
**2.320 × 10⁻⁶** and **2.003 × 10⁻⁶**. The full 35-mode diagnostic RTC
(2 tip/tilt, 7 woofer, 26 tweeter) and 2,403-command synthesis run at unchanged
sampling: default 240/24, fine 720/192, and fine plant/default-calibration transfer.
Own calibrations are remeasured and explicitly compared; COPPER separately
holds calibration fixed to prevent recalibration from hiding a change.

All 320-frame moving GPI runs remain finite and unsaturated. Static 120 nm
inputs converge below 0.001 nm. Default/fine deterministic residual-wavefront
vector differences are at most **0.0000506 / 0.0000463 nm**, computed using the
actual nonorthogonal basis Gram matrix. Transfer and static checks also pass.

Independent CPU scientific replay verifies recorded hashes, actuation basis,
Fourier inputs, Gram matrix, SVD, control updates and WFE. Basis reconstruction
matches the Julia oracle hash; the independent Python basis differs by at most
7.57 × 10⁻¹⁷ m, with Gram error below 5.19 × 10⁻¹³.

## Noisy trajectories: separate interpretation

The seeded noisy GPI trajectories are not bitwise identical:

| Case | Maximum wavefront-vector difference | Maximum scalar residual-WFE difference | Vector RMS difference / realized noise effect |
|---|---:|---:|---:|
| Default | 0.076110 nm | 0.019586 nm | 0.639% |
| Fine | 0.082202 nm | 0.026475 nm | 0.603% |

Both noisy cases exceeded an initial **provisional 0.05 nm vector criterion**.
That result is retained in the comparison JSON. Independent review adjudicated
this criterion as appropriate for deterministic trajectories, not as a
statistical-equivalence test for a single coupled noisy realization. The 0.05 nm
scalar WFE engineering check remains satisfied for the recorded runs; vector
changes are reported independently because scalar norms can hide direction
changes. No production tolerance, controller gain, detector or RNG code changed.

Source shows that Poisson sampling depends discontinuously on intensity via
integer counts or a rounded normal approximation. RNG keys depend on seed,
draw sequence, domain and pixel, not intensity; equal keys do not guarantee
equal detected values after rates change by roundoff. Default trajectories
first diverge materially around frame 61. The particular stochastic event that
caused that divergence was not instrumented and is not claimed as established.
Last-frame differences are input-equivalent electrons, not ADU.

These finite paired-seed results and same-version capture/reset parity checks
support accepting the recorded stability exercise. They do **not** establish
ensemble statistical equivalence, full deployed RTC qualification, or measured
instrument fidelity.

## External optical oracle

The independent OOPAO comparison uses the same matched N120/M32 Float32 fixture
and actual external CPU implementation as the preceding validation, revision
`085d5e50ace0d20fe13cc2da20129d5400166973`. Fresh AOS outputs and OOPAO slope
processing give:

- Flat/ramp image relative L2: **8.046 × 10⁻⁶ / 7.933 × 10⁻⁶**.
- Tip-image derivative relative L2: **7.377 × 10⁻⁵**.
- Tip-slope derivative relative L2: **5.341 × 10⁻⁵** (0.00534%).
- Ramp-slope relative L2: **1.800 × 10⁻⁶**.

No registration, gain, transpose or flux fit is applied. The analytic quotient
normalization derivative reproduces saved OOPAO signals to 9.96 × 10⁻⁹ relative.
The prior physical tip-slope result was 5.263 × 10⁻⁵; the small Float32 change
is recorded, not described as bitwise agreement. Full external provenance and
its precision-file import side effect are preserved in the artifact. These
finite sensor cases do not qualify every instrument parameter or noise model.

## Compute versus memory traffic

A preparation-only benchmark override compares three representations with the
**same new mask/cache kernels**, geometry, weights, batching and source:

- Dense: read a precomputed complex phase cube.
- Separable: multiply two precomputed complex factors.
- Computed: retain only coordinates/offsets and evaluate the original `cis`
  expression while filling each tilted pupil.

The latter two exchange arithmetic for a smaller unique data set. On-demand
storage is 1,728 / 4,416 / 7,808 bytes in the 240/96, 720/192 and 1440/256
fixtures; common pupil-phase scratch is additional and identical across modes.
Computed and dense diagnostic rate maps are byte-identical. Separable maps
differ by at most **6.80 × 10⁻⁸** relative L2, with flux differences below
2 × 10⁻⁸. These diagnostics isolate representation; they are not the original
production baseline. The computed owner is benchmark-only and is not promoted
as a capture/update-qualified production API.

The result establishes a substantial **capacity saving**, while whole-node
timings remain close. It does not prove that dense-phase reads dominate device
time or that recomputation universally improves speed. The previous fine-case
profile attributed 67–73% of kernel time to FFTs. Mask multiplication and
intensity accumulation have low arithmetic per sample and are candidates for
traffic reduction, but hardware bandwidth/cache counters remain inaccessible
under the current NVIDIA policy. No counter permissions, clocks or power
limits were changed. Logical loads and allocated bytes are distinct from
measured DRAM transactions.

## Validation and measurement scope

RTX 3080 10 GiB, sm86; driver 615.71.10, Julia 1.13.2, CUDA 6.4.2,
KA 0.9.45. One Julia/BLAS/FFT thread, GPU jobs serialized, CPU affinity 2.
CPU verification uses affinity 4. Host activity/clocks were uncontrolled.

Completed checks:

- Focused CPU WFS/graph checks: 1,597; KA CPU: 1,020.
- Separable owner: 776 CPU and 709 CUDA checks, including both precisions,
  sampled/nonuniform quadrature, disabled/zero-radius updates and invalid inputs.
- CUDA batch numerical/capture: 457; storage/update: 581; noncentered/complex
  field: 28; captured live modulation/OPD/cache-poison replay: 38.
- CUDA batch policy: 378; optional physical fractional-shift contract: 30.
- GPI CUDA package: 156; complete-loop capture regression: 350.
- Execution terminology/inference/ensemble CPU checks: 35.
- COPPER fast public calibration/production constructors: captured execution,
  flat/nonuniform inputs, seeded reset/replay and zero-allocation checks pass.

Every observed warmed host allocation is zero in the scoped optical-node and
moving complete-loop checks. GPI contributes 80 moving-cycle observations.
Native capture recursively rejects allocation/free nodes. Static traces are
numerical checks; they do not record a separate static allocation series.

Unprofiled completed-node samples exclude detector, RTC, transfers, setup and
capture. Generic instrument-sized cases use an unobstructed 8 m pupil, 750 nm,
3 λ/D modulation and a deterministic sum of OPD sinusoids scaled by 40 nm;
they are not full instrument runtime measurements. Full GPI loop timing is
recorded separately. Timings from profiled invocations are diagnostic only.
Initial small-case measurements showed drift; longer warmup diagnostics are
kept separately. No fixed-rate, WCET, occupancy or peak-bandwidth claim is made.
The virtual 0.5 ms GPI step is not achieved on this hardware.

## Reproduction and artifacts

Raw commands, logs, Project/Manifest files, exact source snapshots, production
patch, baseline-harness copies, external arrays, calibration/loop arrays and
profiler report are archived with a per-file SHA-256 index. The baseline uses
the earlier physical-model campaign, not retired shifted-mask results. Source
worktree paths in manifests must be adjusted when reproducing.

Use `compare_pupil_reuse_validation.py OLD_PHYSICAL_EVIDENCE NEW_EVIDENCE OUTPUT`
for hash-checked calibration/loop comparisons, and `summarize_pupil_reuse.py`
for storage and timing reductions. The old physical arrays are in the preceding
[physical validation archives](../2026-10-22-pupil-tilt-reference/README.md).
The archive includes the independent GPI scientific auditor and basis
reconstruction script; the reproducible 72 MB basis is not stored again.
OOPAO slope processing accepts explicit prefixes and an output path, preserving
historical diagnostics. No dependencies or application release pins changed.

Artifacts:

- `metadata-evidence.tar.gz`: commands, environments, exact source snapshots,
  COPPER arrays, external outputs, logs and the Nsight report. SQLite is
  regenerable from that report and omitted.
- `gpi-evidence.tar.gz`: default, fine and calibration-transfer arrays.
- `raw-sha256.json`: 211 archived file digests, verified by reading the archives.
- `integration-comparison.json`, `performance-summary.json`,
  `external-slopes.json`, `scientific-replay.json`: reduced results.
- [Independent final review](review.md): findings, dispositions and evidence limits.
- `SHA256SUMS`: hashes of this report and all delivered artifacts.

The COPPER benchmark metadata reader was updated separately in `781f41d`;
its fast public calibration/production validation is included here. An initial
baseline invocation failed because a copied harness tried to inspect Git from
outside a repository. The copied baseline revision lookup was corrected and
all baseline cases rerun; both logs are preserved.

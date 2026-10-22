# Independent verification

Reviewed 2026-10-21 (America/Vancouver). No blocking implementation defect or
evidence inconsistency remains in the reviewed scope. This pass inspected source
and existing measurements; it did not launch GPU workloads or change production
code.

## Provenance and method

Implementation revisions: AOS `9a51f458`, Copper `47e305f`, launcher `dad7605`.
The AOS baseline is `074f15eb`. Measurements preceded the implementation commits.
Reconstructing the tracked AOS source/documentation/test diff between baseline and
implementation produced SHA-256
`4d160462f31eda645f1771b7b405e7577cd22a3f39cc74e8bc078426f783c949`, matching
the measured changed-source records. All three measurement harness hashes match
their current files.

An independent Python/NumPy pass read the raw little-endian arrays directly,
recomputed SHA-256, checked dimensions, element sizes and finite values, compared
bytes across revisions, and recomputed errors and timing statistics. It did not
execute or import the campaign comparator. All 76 raw scientific-array files
passed checksum and layout checks. Results agree with [comparison.json](comparison.json).
After archival, all 105 artifact hashes and byte counts in
[provenance.json](provenance.json) were independently rechecked successfully.

## Numerical and capture results

All 14 cross-revision scientific comparisons are byte-identical: 1,985,632 values
covering both 480/32 and 240/96, each with 13 photon-rate probes, 1,432-component
slopes, the complete 1,432 × 277 interaction matrix, and 32-frame atmosphere rates
and noisy detector frames under both captured and stream execution. These are
observed preservation results at fixed discretization.

The 24 dense CPU/CUDA preparation records cover Float32/Float64, odd/even grids,
and stationary, circular and weighted sampled policies. Independently recomputed
norm errors satisfy their recorded `rtol=8eps(T), atol=8eps(T)` criterion. The
largest elementwise error is 9.92018954 × 10⁻⁷; this criterion does not assert an
elementwise absolute error below 8eps(Float32).

Cross-revision identity does not imply identity between execution policies.
Both revisions have the same stream/captured rate relative errors: approximately
1.10 × 10⁻⁷ at 480/32 and 1.07 × 10⁻⁷ at 240/96. The 240/96 noisy frames differ
between policies by maximum absolute 1.02394104 and relative norm 8.38136 × 10⁻⁵;
the 480/32 frames match exactly. These existing differences are unchanged.

The [CUDA log](raw/cuda-tests.log) records 603 storage checks and 26 captured-update
checks passing. Source review confirms that the latter replay one existing
capture after radius changes, disable and re-enable, for separable and full
shifted masks, while checking buffer identities and zero host allocation during
replay. The final [CPU log](raw/cpu-tests-final.log) contains 1,061 passing checks
across 21 test sets. Reset assertions in the campaign passed; the second reset
trajectory itself is not separately archived for independent array comparison.

## Storage and preparation

The two [1440/256 preparation records](raw/preparation-changed-1440-256.toml)
give the following independently checked arithmetic:

| Quantity | Baseline | Changed |
|---|---:|---:|
| Direct retained modulation device bytes | 4,246,732,800 | 7,808 |
| Direct retained modulation host bytes | 1,024 | 1,024 |
| Total host bytes allocated during measured preparation | 4,337,509,736 | 90,783,352 |
| CUDA live pool bytes after preparation | 5,018,157,664 | 781,800,672 |
| Measured preparation seconds | 8.464990304 | 0.469296503 |

The device counts follow 1440² × 256 × 8 and (1440 + 2 × 256) × 4 bytes.
The reduction is 4,246,724,992 device bytes; measured host allocation falls by
4,246,726,384 bytes. Preparation was 18.04 times faster in this pair of runs.
Direct modulation storage excludes propagation, masks, FFT plans, capture and
allocator reserve. Host allocated bytes are cumulative allocation, not peak RSS.
See the [baseline preparation record](raw/preparation-baseline-1440-256.toml)
for the matched warmup and measurement boundary.

## Public Copper constructors

The [mode campaign](raw/copper-modes.toml) completed all six combinations of
production/calibration and default/fast/accurate through the actual public
constructors. Harness review confirms these are complete captured graphs and HIL
boundaries. Recorded shape, selected sampling, finite-output, changed-input and
seeded-replay assertions passed. All 96 recorded host-allocation measurements are
zero. Independently recomputed medians from 128 samples per case are:

| Mode | Sampling | Calibration service, ms | Production service, ms |
|---|---|---:|---:|
| default | 480/32 | 2.586781 | 2.653136 |
| fast | 240/96 | 1.876855 | 1.921655 |
| accurate | 1440/256 | 235.982072 | 238.142162 |

All recorded minima, maxima, means and percentiles agree with their raw timing
samples. The separate baseline/changed graph campaign also agrees with its
reported timing statistics; its captured production medians are
2.6262195 → 2.626495 ms at 480/32 and 1.8972635 → 1.896713 ms at 240/96.

## Limits and disposition

Confidence is high for the reviewed storage, update, capture and mode-selection
behavior. Timings are self-paced service measurements from one sequential
campaign per case on an RTX 3080 with Julia 1.13.2/CUDA.jl 6.4.2. They establish
neither WCET nor deadline compliance, independent-run dispersion, transport
performance or a noisy operational RTC qualification. Accurate mode remains a
finer numerical comparison, with the same shifted-mask approximation and
provisional instrument physics; no new physical oracle was established.

Application dependency pins were not updated. Deployment still requires the
normal dependency-release and application-integration steps.

# Physical pyramid propagation: withdrawal and CUDA validation — 2026-10-22

## Decision

Remove `PyramidShiftedMaskStrategy` and its compact/full shifted-mask storage,
kernels, capture traits, and graph-file selection. GPI and COPPER now select
`PyramidPupilTiltStrategy` at their existing sampling. No gains, modulation
quadrature, detector geometry, or actuator prescription changed. Existing
calibrations must be regenerated and captured graphs rebuilt.

This is a correction to the optical model, in addition to the separate
implementation optimization of FFT batching. The [earlier sweep](../2026-10-22-modulation-batches/README.md)
records the batching improvement and the external discrepancy that triggered
withdrawal. Its shifted-mask timings and self-calibrated loops do not qualify
the replacement model.

A bounded [discrete-operator investigation](../2026-10-22-modulation-batches/shifted-mask-assessment.md)
found a structural limitation: fractional pupil tilts require a dense Fourier
interpolation operator around the fixed mask. Evaluating a pointwise shifted
mask cannot implement that operator. Integer-bin shifts with a periodically wrapped sampled mask agree to roundoff;
changing wrap/sign/scale does not repair fractional shifts. The existing
physical path explicitly tilts the pupil before forward propagation.

## Independent reference

Actual OOPAO CPU execution at revision
`085d5e50ace0d20fe13cc2da20129d5400166973` supplies matched flat, ramp and
±1 nm tip images. No fitted registration, flux scaling or transpose is used.
The directly tuned N120, 32-point Float32 CUDA case has:

| Quantity | Retired shifted mask relative error | Physical pupil tilt relative error |
|---|---:|---:|
| Flat image | 2.3014% | 0.000805% |
| Ramp image | 4.2870% | 0.000793% |
| Tip image derivative | 27.8458% | 0.007562% |
| Tip slopes | 17.4349% | 0.005263% |
| Ramp slopes | 18.9623% | 0.000181% |

The retired tip slope error still has an 8.22% residual after scalar gain
fitting. The slope diagnostic uses the same OOPAO valid pixels, quadrant
ordering and normalization; its analytic quotient derivative reproduces the
saved OOPAO finite difference to 9.96 × 10⁻⁹ relative. Fresh physical CUDA
outputs after removal are byte-identical to the selected-batch physical
outputs used in that comparison. Full provenance, input definitions, CPU and
Float64 cases, reference scripts and arrays are in the earlier sweep archive.
These finite matched cases establish implementation agreement, not complete
instrument fidelity or an error bound for all wavefronts.

## Fixed-sampling batch parity on the retained model

Physical COPPER N240/M96, all 277 actuator columns, five-layer atmosphere,
seed 823, r₀=0.15 m, gain 0.15, 320 frames. Compare original cap-32 batching
against the production-selected batch two. The selected loop holds the
baseline calibration and flat reference fixed.

- Photon-rate relative L2: 8.365 × 10⁻⁸; maximum probe flux change 2.609 × 10⁻⁹.
- Slopes relative L2: 7.401 × 10⁻⁷.
- Calibration matrix relative L2: 5.030 × 10⁻⁶; worst column 1.147 × 10⁻⁴.
- Maximum command difference: 9.237 × 10⁻¹⁴ m OPD.
- Maximum residual WFE difference: 3.977 × 10⁻¹⁵ m; tail WFE ≈75.074 nm.
- Captured noisy plant median: 3.867 → 2.883 ms; calibration 3.832 → 2.797 ms.

Saved seeded noisy detector and atmosphere probes are also compared. The
[comparison](integration-comparison.json) passed all original gates. This
qualifies the batching implementation independently of the withdrawn model.
The earlier repeated randomized sweep remains the stronger timing evidence;
these integration timings are finite workload measurements.

## GPI woofer, tweeter and tip/tilt

The existing 35-mode diagnostic RTC exercises two tip/tilt, seven woofer and
26 tweeter modes with full 2,403-command synthesis. The 320-frame moving input
uses a deterministic 400 nm RMS Fourier phase sequence, common 720-sample
science pupil, seed 823 and virtual 0.5 ms step. Each own-calibrated model
uses identical control settings. Static 120 nm inputs converge below 0.001 nm;
all moving runs stay finite and unsaturated.

| Physical model | Tail WFE off/on noise, nm | Completed cycle off/on noise, ms |
|---|---:|---:|
| Default 240/24 | 139.837 / 139.920 | 3.478 / 3.516 |
| Fine 720/192 | 139.630 / 139.710 | 320.793 / 322.118 |
| Fine plant, default calibration | 139.741 / — | 320.970 / — |

Times summarize three trial medians. Default before/after arrays are identical
(the selected default batch is unchanged). Using the actual nonorthogonal
actuation-basis Gram matrix, the tail residual-wavefront difference between
default and fine models is 0.572 nm RMS without noise and 0.670 nm with noise;
maximum differences are 0.634 and 0.871 nm. The normalized 20 nm calibration
matrix differs by 0.994% relative; worst column differences are 0.828% tip/tilt,
1.030% woofer and 2.017% tweeter. See [sampling](sampling-comparison.json) and
[calibration](calibration-sampling-comparison.json) reductions.

Independent CPU scientific replay verifies saved hashes, analytic actuation,
Fourier inputs, Gram matrix, SVD, control updates, saturation and WFE. The
independent basis reconstructs the recorded Julia hash exactly; the Python
basis differs by at most 7.57 × 10⁻¹⁷ m. These are finite diagnostic control
runs, not HEART, full-order controller qualification, hardware transport,
measured actuator dynamics, or comprehensive instrument validation. The
0.5 ms model period is not achieved on this GPU.

## Captured execution and memory

All three public COPPER modes pass production and calibration CUDA capture,
seeded reset and nonuniform-input checks. Each timing uses 32 warmups and
128 completed samples; all 96 host-allocation observations are zero.

| Mode | Pupil / points | Production median, ms | Direct device phase storage |
|---|---:|---:|---:|
| default | 480 / 32 | 5.225 | 58,982,400 bytes |
| fast | 240 / 96 | 2.923 | 44,236,800 bytes |
| accurate | 1440 / 256 | 457.202 | 4,246,732,800 bytes (3.955 GiB) |

Physical modulation phases remain resident on the GPU and are generated on
that backend. They are dense because the retained propagation consumes them;
FFT workspaces and other model storage are additional. `accurate` names a
finer discretization, not an independent physical guarantee. No sampling
preset was changed based on this campaign.

GPI captures the complete plant and modal RTC cycle, including command
synthesis, integration and device input advancement. Every observed warmed
host allocation is zero. Native capture inspection recursively rejects
allocation/free nodes; default cycles contain 41/52 recursive nodes without/
with noise. Zero observed host bytes and rejected native allocation nodes
are distinct checks, not a claim about driver-internal resource management.

## Validation and reproducibility

RTX 3080 10 GiB, sm86, driver 615.71.10, Julia 1.13.2, CUDA 6.4.2,
KernelAbstractions 0.9.45. One Julia/BLAS/FFT thread, GPU jobs serialized,
CPU affinity 2 for GPU measurements. Host activity and GPU clocks are not
controlled. No occupancy, WCET, paced throughput or transport claims are made.

Completed checks after removal:

- AOS focused CPU graph, pyramid/Bi-O-edge and KA tests: 2,593 passed.
- AOS API, terminology, inference and ensemble checks: 85 passed.
- CUDA batch numerical/capture: 457; physical storage/update/replay: 487;
  batch policy: 378 passed.
- GPI CPU: 2,200 passed; optional external 4HEART maps not enabled.
- GPI CUDA ownership/routing/frame/reset: 156 passed.
- Complete CUDA loop capture regressions: 350 passed.
- COPPER CPU: 120 passed. An initial optional pyRTC include lacked
  UnsafeArrays; only the isolated CPU validation environment was extended
  with UnsafeArrays 1.0.9, then the full selector passed.
- Full instrument campaigns and independent replay described above passed.

[Independent review](review.md) accepts removal with no unresolved findings.
The source patches record execution against AOS `4d00d534`, GPI `13b8c17`,
COPPER `47e305f` and application documentation `dad7605`, plus the contained
working changes. Final commit IDs are available from this report's git history.
No registered dependency pins, releases or live RTC services were changed.
Application use of the new AOS policy/removal requires its corresponding
registered release and a deliberate pin update; instrument source migration
can use the existing physical API meanwhile.

Extract all four `*-evidence.tar.gz` archives into one directory. Their `raw/`
files are indexed by [SHA-256](raw-sha256.json); duplicate payloads use relative
hardlinks within each archive. Archives were extracted and every raw hash
verified. `SHA256SUMS` covers the report and packaged artifacts.

`raw/commands.jsonl` records exact commands and environment variables;
`run_cuda.py` and `run_additional_cuda.py` serialize reproduction. The isolated
Project/Manifest files record versions and local checkout paths; adjust those
paths to the corresponding revisions. No application environment path override
was introduced. The original sweep archive contains OOPAO inputs and provenance.

For independent GPI replay, use `raw/gpi-reconstruct-basis.jl` with
`GPI_VERIFY_REPO` pointing at GPI and `GPI_VERIFY_DIR` at the extracted run;
set `GPI_VERIFY_BASIS` to the generated basis file. Then run the included
`verify_scientific_metrics.py` under the same variables with NumPy and one
BLAS thread. The 72 MB basis is reproducible and is not archived again.
For the generic batch comparator, build a view containing `copper-before`,
`copper-after`, `gpi-before`, and only the `default` subtree plus top-level
files from `gpi-after`; fine/transfer runs are separately audited. Historical
absolute paths in logs describe the actual execution environment.

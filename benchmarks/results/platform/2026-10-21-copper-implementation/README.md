# AOS implementation costs for the Copper CUDA workload

2026-10-21, RTX 3080, driver 615.71.10. This holds the production 480²/32-point prescription fixed. No AOS numerical algorithm, instrument setting or package API changed.

## Completed frame timings

Three fresh Julia processes, 32 warmups and 128 samples per mode. Medians below are medians of process medians; the range is between process medians. The boundary is `step_graph!` through prepared completion, including native five-layer atmosphere and noisy detector, with zero DM commands. It excludes HIL staging, transport and the controller.

| Execution | Median ms | Run-median range ms | Warm host bytes/frame |
|---|---:|---:|---:|
| stream | 2.717 | 2.716–2.733 | [12320] |
| captured | 2.629 | 2.629–2.638 | [0] |

## Native kernel trace

Device-side activity: GPU was busy for 10.56 ms (97.80% of the trace). This is a four-frame profiled window, not whole-application utilization, achieved occupancy, peak efficiency or a fixed-rate test. Profiler durations are diagnostic and are not the unprofiled latency baseline.

| Captured activity | Mean µs/frame |
|---|---:|
| `void regular_fft<1024u, EPT<16u>, 8u, 2u, ` | 839.48 |
| `void vector_fft<1024u, EPT<16u>, 1u, 11u, ` | 804.54 |
| `gpu_pyramid_separable_shifted_mask_batch_kernel_` | 428.74 |
| `gpu_pyramid_modulation_batch_intensity_kernel_` | 395.18 |
| `gpu_moving_layer_replay_kernel_` | 42.96 |
| `gpu_bin2d_kernel_` | 41.01 |
| `gpu__counter_poisson_noise_kernel_` | 19.61 |
| `gpu_pyramid_unmodulated_pupil_field_kernel_` | 18.66 |

Both cuFFT rows include forward and backward transforms. Most of their time is in the 32-plane inverse batch. The shifted-mask fill and intensity reduction are the next large costs. The intensity kernel reads a 256 MiB complex stack plus a 4 MiB prior intensity and writes 4 MiB; its observed duration implies about 700 GB/s of algorithmic array traffic. This is a derived byte rate, not a hardware DRAM-counter measurement. Removing a small final image pass does not remove the dominant complex-stack traffic.

## Serialized node attribution

| Node | Stream device interval µs | Captured-protocol device interval µs |
|---|---:|---:|
| atmosphere | 56.32 | 58.37 |
| pdm | 40.96 | 34.82 |
| pupil_opd_composition | 10.24 | 10.24 |
| pwfs | 2548.74 | 2539.52 |
| detector | 134.66 | 55.30 |

These use separate disposable graphs and an event around direct node submission. They include device idle gaps between submissions and are measured with per-node synchronization; captured attribution calls the capture protocol but does not launch the native graph. They are not additive stage costs or native replay scheduling measurements. The archived host-complete node times additionally use CUDA event waiting, whose task wait policy differs from the prepared graph completion policy; do not substitute them for full-frame latency.

## Implementation findings

1. PWFS propagation is the main warmed target. Benchmark the fixed batch heuristic at unchanged pupil/modulation sampling before selecting another value. Preserve weights, FFT normalization, all points, captured reset/replay behavior and numerical parity. No optimum is established here.
2. Preparation allocates dense pupil modulation phases even when the shifted-mask rate path uses compact axis factors. At 1440/256 this is 3.955 GiB device storage and an equally large temporary host array. At the production 480/32 setting it is 56.25 MiB each. The host temporary is copied once during preparation; it is not a per-frame transfer. Required phases can be generated on device and remain resident. Strategy-specific ownership can avoid the cube entirely for the rate path while preserving diagnostic/reference products. No memory-leak claim or implemented speedup is made.
3. Binning, atmospheric interpolation and DM kernels are secondary measured costs. They may benefit from output reductions or fusion, but cannot explain most of this frame time. The executor already replays one native graph; an executor rewrite is not the first target.

See the [source review](source-review.md) for affected code, existing efficient behavior and validation obligations. Kernel identifiers, grids, registers and shared-memory attributes remain in [activity.json](activity.json) and the full trace. No occupancy counters or roofline experiment were collected.

## Reproduction and provenance

Use an isolated environment with local AOS/Copper sources. The exact benchmark commands and environment snapshots are archived under raw/. Julia 1.13.2, CUDA 6.4.2, KA 0.9.45, CPU 2 affinity, one compute/GC/BLAS/FFT thread, no interactive pool, shared host and uncontrolled GPU clocks. These are sequential service times, not an open-load latency contract.

Timed runs used the archived `profile-baseline-harness.jl` from AOS commit 041cb2a2. A subsequent harness-only change displays the CUDA 6.4 returned profiling object explicitly; it does not change timed paths or numerical source. The native trace uses its separately archived driver. All raw timestamps in activity.json are reductions of the human-readable trace, with its display rounding. No registered application environment was modified.

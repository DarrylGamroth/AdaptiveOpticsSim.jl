# Instrument CUDA capture and HIL validation

2026-10-21. Follow-up to the [instrument execution review](2026-10-21-instrument-execution-review.md). This campaign validates implementation changes on the local RTX 3080; no physical coefficients, controller gains, sampling grids or model algorithms were changed.

## Changes

- GPI detector metadata now matches the pyramid product: normalized pupil sampling 1/60. The cell-integrated photon-rate array is shared unchanged; no resampling or scaling was introduced.
- GPI routing and frame-layout owners qualify their retained CUDA bindings for capture. Native modal OPD expansion enqueues its existing kernel asynchronously in prepared graphs. Public `combine_basis!` remains synchronous.
- CUDA HIL device-to-host staging queues the packed device-array copy on the retained stream, preserves both buffers through the final completion wait, and publishes afterward. Other storage kinds keep ordinary copy behavior. This removes two redundant CUDA cooperative waits and their 160 warmed host bytes.
- Launchers select actual `StreamGraphExecution()` or `CapturedGraphExecution()` objects from `execution = "stream" | "captured"`; omitted selection is stream. CPU capture is rejected. Classic/Copper scenarios select capture. GPI stays stream until a registered AOS 0.2.2 release can replace the application pin. SPIDERS stays stream.

## Local verification

- CPU: `algorithm-graphs` 559/559; `ka-cpu` 993/993; GPI package 103/103; launcher configuration 18/18.
- CUDA: HIL copies/modal edges/storage admission/workspace/waits 104/104; GPI ownership/routing/conversion/full-graph parity/reset 150/150. Reduced GPI regressions use 60² pupils with 4/8 modulation points; the performance campaign below independently uses each full-size default model.
- HIL and calibration tests cover Float32/Float64, nonzero command adoption, invalid-command rejection, reset storage addresses, host visibility, alternate producer streams, caller stream/context restoration, empty and wrong-shaped buffers, and HostMemory fallback. Modal tests cover cancellation, masked cells, NaN/Inf, and zero-mode kernel handling.
- GPI checks cover all command regions, shifted registration, changed commands, conversion with/without ADU calibration, negative/fractional/clipped values, poisoned padding, and evolving atmosphere/detector RNG reset.
- Actual CUDA launchers passed `--check`: Classic/Copper captured, GPI default stream and explicit captured override. No transport connection or RTC process was started.
- Independent source review found no blocker in the modal path, copy lifetime/ordering, or GPI traits. Hardware qualification here is CUDA only. The pre-existing scalar-CPU versus GPU nonfinite mask convention is not changed or claimed equal.
- Fail-before metadata assertion and HIL allocation profile, plus pass-after evidence, are archived. Worker console summaries are labelled as summaries; process logs are preserved where captured.

## Measurement method

Julia 1.13.2, CUDA 6.4.2, KernelAbstractions 0.9.45; RTX 3080, driver 615.71.09. One compute thread, no interactive pool, one GC thread, CPU 2 affinity, BLAS/FFT one thread. Shared host and uncontrolled GPU clocks. No profiler ran during timing.

Three fresh processes per instrument and source variant; 32 warmups, 16 allocation samples and 128 completed-operation timing samples per endpoint. Instrument order and before/after order reverse in round two. Within a process stream precedes capture. Values are medians of the three process medians; ranges show those process medians. This is closed-loop service cost, without offered-rate scheduling, transport, RTC, queueing, or a deadline/tail-latency claim.

The controlled **before-copy** source is AOS `7af2555a` plus the archived modal async patch, using the same fixed GPI nodes/metadata as after. It admits GPI capture and isolates only the HIL copy change. It is not the registered 0.2.1 baseline from the earlier review. Application environments were not modified: private environments under `/tmp` used local AOS/GPI sources; their exact projects/manifests and the patch are archived.

## Full-size after-change results

| Model | Stream graph µs | Captured graph µs | Stream HIL µs | Captured HIL µs (run range) | Captured HIL bytes |
|---|---:|---:|---:|---:|---:|
| classic | 330.5 | 239.5 | 539.0 | 442.2 (441.1–442.5) | 0 |
| copper | 2720.8 | 2631.4 | 2740.2 | 2653.6 (2652.4–2654.7) | 0 |
| gpi2 | 1964.6 | 1836.2 | 2053.1 | 1936.6 (1935.7–1939.4) | 0 |

Captured graph replay and complete captured HIL cycles measured zero host bytes in all 48 allocation samples per instrument and endpoint. Every capture passed the recursive allocation/free-node admission guard. Retained device-used memory was unchanged across each timing interval, and no timed sample collected GC. These checks cover these exact graphs, not arbitrary storage types or transport.

### Isolated HIL-copy comparison

| Model | Before-copy captured HIL µs (run range) | After-copy captured HIL µs (run range) | Host bytes before → after |
|---|---:|---:|---:|
| classic | 445.8 (445.8–445.9) | 442.2 (441.1–442.5) | 160 → 0 |
| copper | 2657.0 (2654.2–2657.7) | 2653.6 (2652.4–2654.7) | 160 → 0 |
| gpi2 | 1943.0 (1941.0–1950.6) | 1936.6 (1935.7–1939.4) | 160 → 0 |

The allocation removal is demonstrated. Timing differences must be interpreted against the run ranges and shared-host variability; the copy change does not establish a broad throughput gain.

## Workload comparison

Derived from the unchanged instrument parameters and pyramid padding rule: Copper has a 480² pupil, 30 pupil samples, 32 modulation points, separation 2 and edge 1. Its padded FFT side is (2×30 + 2 + 2×1)×480/30 = 1024. GPI has a 240² pupil, 60 pupil samples, 24 modulation points, separation 60 and edge 30; its side is (2×60 +60 +2×30)×240/60 = 960. Copper therefore has about 1.52× as many batched FFT elements (32×1024² versus 24×960²). Pupil-grid size alone does not describe FFT cost. This workload difference supports the observed ordering; it does not by itself establish relative implementation efficiency.

Copper samples its pupil at 16 grid pixels per detector pupil sample versus 4 for GPI. The instrument fixes 480 in its current preset; this review found no convergence evidence establishing that 480 is necessary. Reducing that sampling changes numerical discretization and needs an optical convergence study against matched physical inputs and a finer reference, including photon rates, slopes and calibration response. The 32-point modulation quadrature is likewise an accuracy/performance parameter to test against a finer angular reference, rather than assume mandatory. It was not done here. Measured numerical error should constrain the grid and quadrature choices in that follow-up.

## Numerical parity and limits

| Model | Pupil | Maximum CPU/CUDA frame relative L2 | Maximum CPU/CUDA rate relative L2 | Stream/captured frame and rate relative L2 |
|---|---:|---:|---:|---:|
| classic | 240² | 0 | 1.05e-07 | 0 / 0 |
| copper | 480² | 1.02e-05 | 1.02e-05 | 0 / 0 |
| gpi2 | 240² | 3.34e-06 | 3.34e-06 | 0 / 0 |

The full-resolution noiseless deterministic OPD/command probe matched stream and captured frames/rates exactly. This is an implementation parity oracle against matching CPU code, not an independent SPECULA, optical-bench, detector-statistics or closed-loop validation. The physical-model limitations from the earlier review remain. Model cadence, peak SM utilization, occupancy and external RTC latency are not established by these service-time samples.

## Reproduction and delivery

Raw samples, allocation traces, logs, private project/manifests, commands and provenance are in [the evidence directory](2026-10-21-instrument-capture-validation/). The maintained harness is `benchmarks/profile_instrument_workspace_cuda.jl`; it takes `classic`, `copper`, or `gpi2`, `AOS_INSTRUMENT_RESULTS`, and optional `AOS_INSTRUMENT_SAMPLES`. Use one appropriate instrument environment per process. Local validation source overrides belong outside registered application environments.

Production changes validated by this campaign:

- aos: `df63f04d150795b8d9819b297767110ce56aeffb`
- gpi2: `4c29120a264aeb12687d864c1349f139a946dd79`
- launchers: `f2d3eb1bf951e35d6ddd10665a4e07ccf0fb9dc7`
- classic: `603e58ba8ab278369d3b1ac138952a8fc6234db8`
- copper: `1ae815314582f666a4c9de3d87a931c591e7a422`

Registered application AOS/HIL pins remain unchanged. Full GPI capture and zero-allocation captured HIL staging require publishing the AOS changes and then updating the application AOS pin; no release or deployment was performed. Launcher `--check` prepares without connecting to PipeWire or starting an RTC.

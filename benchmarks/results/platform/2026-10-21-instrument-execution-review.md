# Current instrument execution and numerical review

2026-10-21. Review of the unchanged `adaptive-optics-simulations` application environments and their local REVOLT Classic, REVOLT Copper and GPI2 instrument checkouts. No production implementation, dependency, calibration, transport or RTC-service changes were made. The new review harness runs against each existing application environment.

## Findings

1. **The scenario launchers select stream execution, not CUDA Graph replay.** All three constructors default to `StreamGraphExecution`; `scripts/scenario.jl` forwards no execution policy. Classic/Copper explicit capture works. GPI2 explicit capture rejects its custom command-routing owner; frame-layout also lacks capture qualification. An explicit application execution setting plus instrument qualification is the most direct implementation improvement.
2. **Zero allocation has a narrower boundary than the complete application.** Captured Classic/Copper GPU graph replay measured zero Julia heap bytes. Their complete HIL cycle, including host frame staging and same-sequence zero-command adoption, measured **160 bytes**, confirmed in 16 additional samples per instrument in fresh processes. The native captured graphs contain only kernel/memset nodes, with no allocation/free nodes. Stream execution allocates. Transport/RTC allocations were not measured.
3. **Current stream workflows leave substantial host/device gaps.** Traces show substantial host/device gaps. Captured Copper keeps kernels active about 98% of the traced graph interval, but a busy GPU is not proof of optimal kernel throughput. Classic's largest device-time contributor is counter-based Poisson noise; Copper/GPI2 are dominated by FFT work. Counter coverage is deliberately limited to selected occurrences.
4. **Numerical self-agreement is good; instrument accuracy is not independently qualified.** Full-resolution deterministic CPU/GPU checks passed. Classic/Copper deterministic GPU stream and captured frames/rates matched exactly for the selected probe. These tests cannot certify the underlying optical model, detector distributions, physical calibration or closed-loop behavior.
5. **SPECULA-labelled regression data are not independent SPECULA oracles.** The checked-in bundle generator calls AOS's own `compute_reference_expected` and writes those values. The metadata identifies a contract bundle, without a producing SPECULA revision. Genuine pinned OOPAO component exports exist, but no complete-instrument independent oracle was found for these three current graphs.

## Measured execution

NVIDIA RTX 3080, 68 SMs, 10 GiB; Julia 1.13.2, CUDA 6.4.2, KA 0.9.45, registered AOS 0.2.1 tree `1a809e78c898ce01a58915b7f4fb71dc92d4711e`. One compute thread, BLAS/FFT one thread, CPU 2 affinity. Uncontrolled clocks/shared host. One fresh unprofiled process per instrument, fixed stream-then-captured order, 32 warmups and 32 samples per endpoint. Values below are descriptive medians, not tail/deadline guarantees.

| Model and full pupil grid | Model period | Stream graph | Explicit captured graph | Stream HIL cycle | Captured HIL cycle |
|---|---:|---:|---:|---:|---:|
| Classic, 240² | 2,000 µs | 271.9 µs | 176.3 µs | 480.5 µs | 383.9 µs |
| Copper, 480² | 2,000 µs | 8,004.2 µs | 2,651.1 µs | 7,996.2 µs | 2,657.0 µs |
| GPI2, 240² | 500 µs | 7,999.6 µs | Rejected at command routing | 8,003.1 µs | Not admitted |

Graph time ends at completed GPU outputs. HIL cycle also stages the host detector frame and adopts an immediate zero command. There is no external command arrival, transport, reconstruction or controller work. The README transport launcher uses three Julia threads; that concurrent application was not timed here. Separate sample sets explain cases where the HIL median is slightly below the graph median. Copper does not sustain its model's 500 Hz cadence on this run even with capture. GPI2's 2 kHz model cadence is not achieved by its current stream path. This does not measure deployed closed-loop latency.

The application environments pin registered AOS 0.2.1 and HIL 0.2.2. The recent development changes at AOS `4bb78685` are not automatically used by these environments. Preserve the registered dependency policy: qualify/release changes, then update bounds and seal a new environment; do not add development path overrides.

### Activity traces and allocation boundaries

Separate Nsight Systems processes captured 64 graph steps and 64 HIL cycles in named NVTX regions per admitted policy. Kernel-active fractions are unions of kernel intervals divided by NVTX-region duration. Device activity additionally includes copies and memsets. They are **not SM occupancy, achieved FLOP/s, or unprofiled utilization measurements**. Instrumentation changes host waits and timing; use the unprofiled medians above for service cost.

| Model | Stream graph kernel-active | Captured graph kernel-active | Captured HIL kernel-active |
|---|---:|---:|---:|
| Classic | 48.5% | 82.3% | 39.0% |
| Copper | 28.4% | 97.6% | 71.1% |
| GPI2 | 18.0% | Not admitted | Not admitted |

Classic captured HIL device activity including copies is about 74.1%; Copper's is 71.5%. HIL-copy/host coordination and trace effects should be distinguished from pure graph kernel throughput. In the traces, Poisson noise is about 26% of Classic's summed kernel time; Copper's `regular_fft` plus `vector_fft` families are about 63%; GPI2's `regular_fft` family is about 64%.

Warmed Julia allocation samples:

| Model | Stream graph | Captured graph | Stream HIL cycle | Captured HIL cycle |
|---|---:|---:|---:|---:|
| Classic | 14,928 B | 0 B | 15,248 B | 160 B |
| Copper | 13,024 B | 0 B | 13,344 B | 160 B |
| GPI2 | 16,512 B | Not admitted | 16,832 B | Not admitted |

Stable `CUDA.used_memory` was also checked, but is only retained-memory evidence. Independent recursive native graph inspection found Classic 34 kernel + 1 memset nodes and Copper 34 kernel + 2 memset nodes, with no allocation/free nodes for these exact captures. The deployed 0.2.1 capture implementation does not contain the newer retained cuBLAS workspace and recursive rejection guard; these successful shapes do not replace general capture admission protection.

### Focused hardware counters

Nsight Compute sampled one Poisson occurrence in Classic and the dominant batched `regular_fft` occurrence in Copper (skip the first matching single forward-FFT kernel). Kernel replay metrics are diagnostic; durations are not service-time measurements. Raw reports and metric names/units are retained.

| Sampled kernel | Achieved occupancy | Theoretical occupancy | SM throughput | DRAM throughput |
|---|---:|---:|---:|---:|
| Classic photon Poisson, 484 × 256 threads | 61.47% | 100% | 48.69% | 0.16% |
| Copper batched FFT, 4,096 × 512 threads | 32.65% | 33.33% | 15.93% | 89.73% |

Copper's FFT uses about 64 KiB dynamic shared memory per block, limiting residency to one block per SM. It achieves almost that occupancy limit while using nearly 90% of peak DRAM throughput. Increasing occupancy is not itself a useful optimization target here. Examine memory traffic and exact FFT plan/layout alternatives; do not reduce scientific sampling to make occupancy look better. Source already permits modulation batches up to 32 points under a memory budget; the trace's dominant FFT covers the 32-plane Copper batch. Classic is not DRAM-bound in this occurrence: examine sampler control flow/integer work, workload imbalance and launch tail effects with code-generation evidence before proposing a specific cause or claiming a gain. GPI2 kernel counters were not collected; its timeline establishes launch/wait gaps and FFT dominance, not per-kernel efficiency.

## Scientific accuracy and algorithms

The reviewed equations use additive metre-OPD coefficients, the stated Gaussian mechanical coupling, square-root layer variance weights, λref/(2π) phase-to-OPD conversion, 2π OPD/λ WFS phase, and photon rate × QE × exposure. Source review found these conventions consistent within scope. It did not prove every algorithm or physical model correct.

| Full-resolution noiseless check | Classic | Copper | GPI2 |
|---|---:|---:|---:|
| CPU/GPU photon-rate relative L2 error | 2.71 × 10⁻⁷ | 1.02 × 10⁻⁵ | 3.33 × 10⁻⁶ |
| CPU/GPU detector-frame relative L2 error | 0 | 1.02 × 10⁻⁵ | 3.33 × 10⁻⁶ |
| GPU stream/captured photon-rate and frame error | 0 | 0 | Not admitted |

One deterministic nonzero OPD pattern and sparse command probe were used; these are self-parity checks, not a broad physical accuracy envelope. Production stochastic reset reproducibility passed separately for each admitted execution policy. No claim of noisy CPU/GPU draw equivalence is made.

Accuracy priorities from the [independent science review](2026-10-21-instrument-execution-review/science-review.md):

- **SCI-002, confirmed metadata defect:** GPI2 supplies normalized detector sampling 2/60 while the Pyramid producer's convention is 1/60. Present cell-integrated pixel values are unaffected. Correct the metadata and assert matching producer/consumer domain, scale and wavelength.
- **SCI-001, known approximation:** GPI2's clipped-Gaussian EM multiplication gives Y/g = max(N + √N Z, 0). Conditional on exactly one input electron, before read noise, its mean is 1.083315 instead of 1 and it has a 15.8655% atom at zero. This is a derived distribution bias, **not a measured 8.33% whole-frame bias**. Use a conditional-gamma CPU reference and distribution tests across flux/CIC levels. AOS 0.2.1 already has that CPU model, but the graph node does not expose its selection. GPU implementation requires its own numerical/distribution and capture qualification. Copper's gain-one/excess-noise-one configuration is different.
- **SCI-006, known photon-noise approximation:** the counter sampler uses bounded inversion below mean 30 and rounded/clipped Gaussian draws at means ≥30. Check independent Poisson CDFs, moments and tails around 29.9/30/30.1 and representative high flux. Preserve the current distribution when optimizing implementation; qualify any exact/reference sampler as a separate accuracy-profile change.
- **SCI-003, temporal validity limit:** atmosphere is sampled once per exposure, then its photon rate is held. Compare 1/2/4/8 temporal samples using a common atmosphere history and command timing. Integrate photon rates; apply detector read noise once. Averaging OPD before nonlinear WFS formation is not equivalent.
- **SCI-004, numerical validation gap:** Copper and GPI2 select the approximate shifted-mask Pyramid strategy. Compare against pupil-tilt propagation at the actual pupil geometry/modulation/aberration range. Refine modulation 32→64→128 for Copper and 24→48→96 for GPI2, and refine pupil sampling separately. Gate on flux, images, normalized sensor signals, interaction columns and optical gains before changing scientific resolution or quadrature for speed.
- **SCI-005, physical qualification gap:** analytic DM influence functions, optical registration, detector maps and command-to-hardware calibration remain provisional. Same-model push-pull calibration, seeded replay and CPU/GPU agreement do not validate measured instrument fidelity. Preserve explicit provenance and acquire matched flat/±probe/noise datasets.

### Independent oracles

Use actual exports from a pinned [SPECULA checkout](https://github.com/ArcetriAdaptiveOptics/SPECULA) with identical explicit OPD, pupil, wavelength, photon normalization, detector window and sampling. Start noiseless, compare intermediate optical products, then add temporal and statistical detector tests. Shared seeds alone do not give matching screens or random draws across engines. Record producing commit, configuration, dependency versions, input/output hashes, units and coordinate conventions. Keep independent exporter code separate from AOS expected-value generation.

AOS's genuine pinned OOPAO exports provide useful component evidence, including SH and a narrow gain-sensing Pyramid scenario. They do not establish current full-instrument equivalence. Some standalone Pyramid reference text files are not admitted by the current manifest and therefore do not contribute maintained regression coverage. See the reconciled [evidence inventory](2026-10-21-instrument-execution-review/evidence-inventory.md).

## Specific executor and node API improvements

The [downstream API review](2026-10-21-instrument-execution-review/executor-review.md) recommends targeted borrowing from JuliaFilterGraph:

1. **Capture conformance:** qualify GPI2's command routing and frame layout with existing AOS hooks, testing changed commands, exact 97/2304/2 ordering, both frame conversions, padding, reset, fixed addresses and allocations. No new executor is needed for these stateless nodes.
2. **Public boundary metadata and semantic accessors:** expose immutable port schema/type/shape/layout/identity and needed exposure/model-cadence queries. All three HIL extensions currently inspect concrete node-owner fields. Preserve instrument-owned physical mappings and transport conversions. Optical metadata propagation would also reduce manual scale reconstruction such as SCI-002.
3. **Retained failure diagnostics:** preserve first-failure node/phase/cause information after preparation/submission/completion faults, while respecting AOS concrete storage rules and successful-path allocation budgets.
4. **Failed-drain lifetime investigation:** use failure injection and lifetime checks before adding retention/recovery machinery. Source suggests a question; no actual unsafe buffer release was demonstrated.

Prepared pinned host staging is another bounded experiment for Classic's HIL boundary. JFG's captured pipeline retains pinned host slots; AOS already has pinned storage in its optional output-snapshot implementation, while ordinary HIL buffers are `Array`s. Test single-frame prepared pinned staging, transfer cost and lifetime without importing a multi-frame pipeline or changing capacity-one ownership. No pinning speedup or cause of the 160-byte allocation has been established by this review.

Defer declaration macros, arbitrary stream joins, conditional publication and per-node rollback absent an instrument requirement. They would introduce lifecycle and model-time semantics beyond these complete-frame plants. Document public stateful capture hooks when a downstream stateful adapter needs them; the present GPI2 adapters are stateless.

## Recommended sequence

1. Correct GPI metadata; establish the independent numerical and detector reference gates.
2. Qualify GPI capture and add explicit application execution selection; qualify/release the retained-workspace/allocation-node guard and other selected development changes before updating registered environments.
3. Remove measured HIL host allocations and profile host staging without changing command/frame order, capacity-one ownership or model time. Keep transport/RTC validation separate.
4. Optimize the demonstrated device bottlenecks while holding the numerical model fixed: Poisson-sampler implementation details for Classic while preserving its current distribution, FFT plan/batch/layout choices for Copper/GPI2. Measure each change with unchanged resolution, modulation quadrature, source flux and noise distribution. Current shifted-mask code already batches modulation; increasing batching is an experiment, not an assumed missing feature or guaranteed gain.
5. Adopt temporal integration or a better EM distribution as explicit accuracy-profile changes with convergence/statistical evidence, rather than presenting them as implementation-only optimizations.

The [machine-readable evidence](2026-10-21-instrument-execution-review.toml) retains raw timing arrays, exact revisions/manifests, capture failures, activity traces, selected hardware counters, allocation/node inspections and independent source reviews. These are local software/numerical measurements, not on-instrument, on-sky or deployed-controller validation.

Independent final review recomputed timing medians, allocation records, all ten activity unions and selected hardware metrics, and checked all 64 raw archive hashes and source/manifest identities. No evidence blocker remained within the stated scope.

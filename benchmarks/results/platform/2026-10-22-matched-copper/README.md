# Matched Copper 240/96: AOS, OOPAO and SPECULA

Completed locally on CPU and CUDA. Photon-rate, flux and calibration gates
pass; strict cross-engine raw-slope gates fail. All three CPU/CUDA comparisons
pass. No production code or instrument model was changed.

## What was compared

The measured product is a resident 64×64 expected photon-rate image, starting
from resident OPD in each engine's native units. Native phase formation,
propagation, normalization, sum binning and completion are timed. Input staging,
unit conversion, diagnostic copies, detector noise, slope reduction,
calibration reduction, atmosphere, DM execution and RTC are excluded.
These are optical-node service times, not complete closed-loop frame times.

This is the actual Copper fast optical prescription: D=1.22 m, λ=550 nm,
240 pupil samples, 96 equally weighted physical pupil tilts, nominal modulation
5 λ/D, 30-pixel pupil images, separation 2, edge 1, 512² FFT and 64² output.
It is not the earlier generic 8 m / 750 nm / 3 λ/D fixture or a fine-grid run.
All runtime real arrays are Float32 and FFT fields are ComplexF32. Cold factor
construction follows each implementation's arithmetic. The common endpoint
modulation coordinates imply a sampled radius 5×240/239 ≈ 5.02092 λ/D on
physical pitch D/240; this is not an independent exact-radius qualification.

Seventeen shared OPDs cover flat, piston, x/y-axis tilts, focus,
astigmatism, coma, finite spatial frequencies, mixed modes and deterministic
bandlimited screens at 100/300 nm RMS. The latter are finite Fourier probes,
not independently validated atmospheric realizations. Calibration uses every
column of the provisional Copper Gaussian277 OPD basis at ±20 nm. The shared
basis validates optical response to those inputs; external DM implementations,
command units and controllers are not tested.

## Sources and setup

| Library | Frozen source |
|---|---|
| AOS | `b590aa26a6f6f976ff0f821b138b43c1a84bea77` |
| Copper fixture authority | `781f41d38f784d0761e270b632d2c709c7234a76` |
| OOPAO | `168fad2a8d98721a270c1b066b813b716680b5ea` |
| SPECULA | `a06bb495ae29a6ec160628069598b32631b9f44f` |

External sources were cloned into isolated, detached checkouts. Exact pins were
verified against upstream HEAD and retained with the results. User checkouts
were not edited. The source bootstrap creates missing pinned checkouts and
verifies existing ones; setup cloning/fetching never occurs in a timing sample.
See [the harness guide](../../../matched_copper/README.md) for reproduction.

The host is an AMD Ryzen 7 6800H with an NVIDIA RTX 3080 10 GiB. Julia 1.13.2,
CUDA.jl 6.4.2 and KernelAbstractions 0.9.45 were used. Python 3.13.5 used NumPy
2.2.4, SciPy 1.18.1 and CuPy 14.2.0. Julia uses CUDA runtime 13.4; CuPy is linked
against CUDA 13.2 with locally installed CUDA 13.4. Both use cuFFT 12.4.0.
The exact Julia manifests, Python lock, build information and GPU configuration
are archived. No package update, clock change or permission change was made to
system or instrument environments.

## Explicit matching boundaries

AOS uses its physical pupil-tilt rate node, unchanged. OOPAO uses native
propagation and sum binning, with resident CuPy mode on CUDA. Its CPU benchmark
subclass fixes one SciPy FFT worker per job to prevent nested oversubscription;
1 and 8 native jobs are compared. The optional OOPAO focal-plane diagnostic is
disabled.

SPECULA cannot request the identical prescription through its unmodified
geometry defaults: its padding guard produces FFT546 and an internal 68² image,
it uses a circular stop, and its prism orientation differs. A benchmark-only
subclass explicitly selects FFT512, a square stop and the matching prism sign.
Propagation, modulation, native normalization, rebinning and processing-object
lifecycle remain inherited. Native SPECULA normalization fixes total output
flux; its total-flux agreement alone does not independently establish FFT
energy conservation. Native auxiliary PSF/transmission work is included.
An untouched native diagnostic is archived separately and excluded from matched
performance comparisons. Its flat output loses about 0.295% of incident flux
under its different prescription; this is not a matched-case accuracy failure.

AOS captures the full optical node. SPECULA uses its native trigger-code CUDA
Graph; preparation and postprocessing remain part of measured execution.
OOPAO uses its native resident batched CuPy implementation. No artificial CUDA
Graph wrapper was added to OOPAO. Inputs vary across probes; exact repeated
replays and nonconstant flat/tilt responses guard against stale capture.

## Accuracy method and findings

The fixed protocol and all failed checks are retained. No gain, flux,
registration, transpose, pixel shift or calibration-column fit is applied.
All comparisons use the complete 64² image. The common I4Q readout uses the
fixed 716-pixel circular support in each 30² pupil, the Copper quadrant signs,
and one analytic incident-flux divisor. Raw slopes are compared without flat
subtraction. The absolute slope metric is the largest component error.
Calibration is D[:,j]=(s(+h bⱼ)−s(−h bⱼ))/(2h), with all 277 columns and Float64
host reductions. Every column is checked using the preregistered denominator
floor; no columns are removed. This validates the shared readout applied to
native photon rates, not external native slope/calibration APIs.

| Comparison with AOS | Rate L2 (%) | Calibration L2 (%) | Worst column (%) | Raw slopes |
|---|---:|---:|---:|---|
| OOPAO CPU | 0.000584459 | 0.000826499 | 0.0179864 | FAIL |
| SPECULA CPU | 0.00104501 | 0.00085198 | 0.0181105 | FAIL |
| OOPAO CUDA | 0.00058478 | 0.00079852 | 0.0175685 | FAIL |
| SPECULA CUDA | 0.0010454 | 0.000845865 | 0.0188038 | FAIL |

All four comparisons pass the image, flux and all 277 calibration-column gates.
The largest raw-slope component discrepancy is 3.36×10⁻⁵ against the fixed
2×10⁻⁶ limit. Some weak excited-mode relative checks also fail. These failures
remain in every relevant summary; no threshold was relaxed.

| CPU versus CUDA | Rate L2 (%) | Calibration L2 (%) | All checks |
|---|---:|---:|---|
| AOS | 1.8501e-05 | 0.000654151 | PASS |
| OOPAO | 1.92063e-05 | 0.000460401 | PASS |
| SPECULA | 3.61154e-05 | 0.0005196 | PASS |

The CPU 1/8-thread comparisons also pass. All calibration columns are nonzero.
An independent reduction rebuilt the 1432×277 matrices directly from the rate
arrays and confirmed these results. It additionally checked photon conservation
over every calibration frame: maximum relative error 3.49×10⁻⁷. These checks establish numerical agreement at fixed sampling. Closed-loop
and real-instrument validation require separate experiments.

### Separate causal diagnostic

On CPU flat and 20 nm x-tilt only, the realized AOS/OOPAO centering phasors differ
by up to 2.2483×10⁻⁴ in complex value. AOS forms large phase arguments in Float32;
OOPAO forms its arguments in Float64 before casting. Replacing only the
centering phasor in an isolated AOS instance reduces the raw-slope discrepancy
from approximately 2.79×10⁻⁵ to 3×10⁻⁷. Replacing only the prism factor does not
remove it. This establishes the dominant cause for those two probes, not every
probe or the SPECULA comparison. No production factor or accepted benchmark
input was replaced. Cold factor construction in higher precision is a concrete
accuracy follow-up that need not add runtime work.

## Performance method

Each configuration has three fresh processes, 32 warmed calls and 128 completed
service-time samples per round. AOS also performs one allocation-probe warmup
and eight untimed allocation observations before sampling. Engine order is
seeded and rotated. Every sample uses one of the same seventeen OPDs, staged before timing. GPU timing
waits for all native postprocessing, including work outside SPECULA's captured
stream. CPU runs use one physical core or an eight-physical-core budget;
SPECULA's native NumPy FFT path is serial. BLAS and Julia task threading are
limited to one. Requested and realized job/FFT/batch policies are recorded.
Serial and GPU submission use logical CPU 2; the parallel
budget uses 0,2,4,6,8,10,12,14 (one sibling per physical core). These cores are
not reserved exclusively. Framework-specific preparation times are diagnostic
and have different module-import boundaries. Ratios include differences in FFT
libraries, batching, native auxiliary products and API overhead. FFTW.MEASURE
can choose different plans in fresh processes; their saved probe outputs must
pass the existing 2e-6 relative image limit and the slope/flux/piston checks.
In-process repeated replay must be exact.

Report the median of round medians, their range and the three observed p95s.
Retain slow samples. These unpaced measurements do not establish fixed-rate
RTC deadlines, WCET, overload behavior or extreme-tail latency. Existing AO/GUI
workloads remained active on this shared host. Telemetry records load,
processes, GPU clocks, power, temperature and memory; this is not an idle-machine
hardware limit. Numerical failures remain failures regardless of speed.

## Measured service time

Milliseconds; central value is the median of the three round medians. The
range is the three round medians; p95 values are shown in round order. These
measurements use the same cyclic probe mix, not an atmospheric arrival process.

| Configuration | Thread/execution policy | Median (ms) | Round-median range (ms) | p95 by round (ms) |
|---|---|---:|---:|---|
| aos-cpu1 | FFTW 1, BLAS 1, Julia 1 | 363.795 | 351.727–379.224 | 395.471, 408.554, 468.540 |
| aos-cpu8 | FFTW 8, BLAS 1, Julia 1 | 257.095 | 243.627–344.822 | 291.004, 388.458, 336.795 |
| oopao-cpu1 | 1 job × 1 FFT worker | 499.322 | 487.021–500.563 | 565.754, 525.089, 580.300 |
| oopao-cpu8 | 8 jobs × 1 FFT worker | 202.701 | 173.492–207.147 | 200.804, 223.973, 231.605 |
| specula-cpu1 | native serial NumPy FFT | 2186.163 | 2097.077–2236.072 | 2256.201, 2298.268, 2362.933 |
| aos-cuda | captured node, FFT batch 2 | 2.664 | 2.663–2.675 | 2.686, 2.862, 2.693 |
| oopao-cuda | resident CuPy, 96-point batch | 6.388 | 6.332–6.407 | 6.399, 6.596, 6.609 |
| specula-cuda | native captured trigger, 96-point loop | 4.885 | 4.686–4.909 | 5.377, 5.356, 5.241 |

AOS is fastest on this GPU. OOPAO is fastest in the tested eight-core CPU
configuration; AOS is faster in the single-core comparison. The eight-core CPU
ranges are broad under the shared load, particularly AOS. These are descriptive
measurements of the configured implementations.

AOS records zero measured host allocations after warmup on CPU and CUDA, and
its CUDA pool equality check passes. External pool counters are preserved as
diagnostics; they do not establish zero allocations or comparable peak VRAM.
Each OOPAO CUDA round releases 16 KiB of live default-pool storage.

The scalar AOS CPU path prepares modulated pupils individually, while its GPU
batch path already reuses the common pupil phase. Profiling equivalent reuse
and batching on CPU is a concrete follow-up. This campaign changes no runtime
implementation or modulation algorithm.

## Evidence

`fixture.tar.zst` contains the exact common input arrays and metadata. The
three `raw-*.tar.zst` archives contain all eight full validations and 24 timing
runs, including arrays, metadata and logs. `analysis.tar.zst` contains every
comparison, failed gate and Float64 derived matrix; `core.tar.zst` contains the
source snapshot, manifests, environment, command journals, telemetry, diagnostic
and independent review. No external source clone or Python environment is
bundled.

[summary.json](summary.json) retains all round samples and reduced metrics.
[review.md](review.md) records the independent review and finding dispositions.
`raw-sha256.json` maps every archived member to its checksum and size;
`archives-sha256.json` checks the archive files themselves. Extract all six
archives into one empty campaign root to reproduce numerical reductions. Source
checkouts and installed environments must be prepared separately as described
in the harness guide.

Archive verification decoded all six streams and checked the sizes and SHA-256
digests of all 231 members against the index. Every archive is below 44 MiB.

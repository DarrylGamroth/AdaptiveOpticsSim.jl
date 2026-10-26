# Supported Production Surfaces

Status: active

## Purpose

This document states what the repository currently treats as maintained and
qualified. A type existing in source is not sufficient; a supported surface
needs a public contract, focused tests, and relevant backend evidence.

## CPU Baseline

The Julia 1.12 CPU baseline includes:

- explicit telescope, source, aperture, optical-product, propagation, direct
  imaging, NCPA, and controllable-optic operations
- Kolmogorov, finite multilayer, and infinite multilayer atmosphere evolution
  plus direction rendering/batching
- Shack-Hartmann, Pyramid, Bi-O-edge, Zernike, and Curvature physical formation
  and complete acquisition, plus LiFT physical forward prediction and analytic
  Jacobian callback within the limits recorded in
  [`model-validity-matrix.md`](model-validity-matrix.md); operational Curvature
  estimation belongs to registered FGA, and weighted-centroid estimation is not
  an AOS surface
- conventional frame and counting/channel detector families covered by their
  registered suites
- plant-side calibration observables, reusable AOC inverse products, and
  tomography; interaction-matrix and misregistration acquisition require a
  caller-owned typed `WFSMeasurement` and an explicitly composed measurement
  callback; maintained RTC reconstruction and controller execution is in
  FilterGraphAlgorithms/JuliaFilterGraph
- static `AlgorithmGraphs` definitions in Julia and YAML
- direct and explicit one-frame delayed graph links
- exact model-time values and deterministic drivers
- prepared coherent wavelength synthesis with explicit photon weights,
  per-channel path coherence, and response
- explicitly calibrated pixel MTF filtering before detector acquisition
- complete-frame InGaAs graph acquisition on an explicit metric pixel grid and
  integrated spectral channel
- the transport-neutral lockstep `PreparedGraphHILBoundary` and held-probe
  `PreparedGraphCalibrationBoundary`
- the built-in graph-node catalog exercised by `algorithm-graphs` tests

Warmed allocation claims apply only to the concrete paths with allocation tests.
They are not a promise for arbitrary user callbacks or every optional library.

The CPU qualification environment (`test/cpu`) and CPU benchmark environment
(`benchmarks`) pin KernelAbstractions 0.10.1 and the compatible
[AOC source revision](https://github.com/DarrylGamroth/AdaptiveOpticsCalibration.jl/commit/55392d16aff769d1678c2940e253f68a09846dd1).
KA CPU kernels and algorithm graphs are checked on Julia 1.13.2. This version
selection also applies to ordinary CPU execution without changing its dispatch
to KA kernels. Run the focused checks from the repository root with:

```sh
OPENBLAS_NUM_THREADS=1 julia --startup-file=no --threads=1 --project=test/cpu \
  -e 'using Pkg; Pkg.instantiate(); Pkg.test("AdaptiveOpticsSim"; allow_reresolve=false, test_args=["ka-cpu", "algorithm-graphs"])'
```

## Downstream Instrument Packages

AOS supports the reusable algorithms, graph adapters, static YAML schema, and
lockstep HIL boundary used by downstream instrument packages. It does not own
an instrument's command geometry, detector profile, calibration policy, or
scientific acceptance claim.

[REVOLTClassicSim.jl](https://github.com/DarrylGamroth/REVOLTClassicSim.jl)
and
[REVOLTCopperSim.jl](https://github.com/DarrylGamroth/REVOLTCopperSim.jl)
maintain their executable graph files and pyRTC validation separately. Their
support status is governed by their own tests and documented assumptions.

## Accelerator Surfaces

The AMDGPU test and benchmark environments pin AMDGPU 2.9.2 and
KernelAbstractions 0.10.1. This includes the upstream
[Julia 1.13 precompilation fix](https://github.com/JuliaGPU/AMDGPU.jl/commit/d276a5cd42bf52a869f8f414873526ee3cc999bd)
for checked integer division, transpose matrix-vector multiplication, and
cropped array copies. These environments pin the compatible
[AOC source revision](https://github.com/DarrylGamroth/AdaptiveOpticsCalibration.jl/commit/55392d16aff769d1678c2940e253f68a09846dd1);
registered AOC 0.17 only permits KA 0.9.
The separate CUDA test and benchmark environments pin CUDA 6.4.2 and
KernelAbstractions 0.9.45 to retain the qualified allocation behavior. AOS allows
both KA versions; each environment selects one version for its process.
These exact constraints protect the selected backend versions. Full dependency
reproduction also requires retaining the resolved manifest with validation or
benchmark results. Dependency upgrades require renewed correctness and
allocation checks. All backend test and benchmark projects pin this AOC
revision. It keeps cold `KernelExecution` construction out of line to prevent
partially initialized specification boxes on Julia 1.13.2
([Julia compiler defect](https://github.com/JuliaLang/julia/issues/62524)).
This protects the validated preparation path; it does not repair the general
Julia compiler defect. These validation claims require the pinned AOC revision.

CUDA.jl 6.4.2 or newer is required. CUDA is the primary GPU performance-optimization target. AMDGPU is the
secondary portability and qualification target and remains the automated
release-gated accelerator because the maintained local ROCm target is
continuously available. CUDA validation uses the local NVIDIA device. Both
targets require:

- scalar indexing disabled
- exact device residency for inputs, plans, state, workspace, and products,
  except for an operation's documented preallocated host-mirror strategy
- no undeclared CPU fallback
- explicit synchronization at public completion boundaries
- numerical or statistical parity against the declared CPU reference

Pyramid modulation FFT batches are fixed during preparation. CUDA uses measured
choices for qualified FFT-size/point-count pairs with `ComplexF32` fields on
sm86 devices with 5 MiB L2; other cases retain the existing divisor selection
with a 32-point limit and 256 MiB field-stack budget. Physical sampling and
modulation weights are unchanged. The [batch sweep and validation](../benchmarks/results/platform/2026-10-22-modulation-batches/README.md)
record the qualified combinations, numerical differences and external reference.
The qualified fine `(FFT side, points)` pairs `(2880, 192)` and `(3072, 256)`
use sixteen-point tiles and a 1536 MiB main field-stack bound. The 2880 case
also retains transposed frequency planes between axis FFTs; its additional
transposed stack is outside that bound. FFT layout and tile size are selected
at preparation. Resizing retains the prepared layout family; fresh preparation
applies the backend policy again. The [fine FFT investigation](../benchmarks/results/platform/2026-10-22-fft-layout/README.md)
records timing, direct workspace storage, calibration and closed-loop parity.
Accelerator Pyramid propagation stores each physical pupil tilt as two separable
complex factor arrays, computes the common OPD phasor once per pupil propagation,
and reuses mask values across modulation planes. These implementation choices
retain every physical tilt and both FFTs; they do not select a reduced optical
model. Prepared factors and pupil-phase scratch remain device-resident and bound
to captured execution. CPU and Bi-O-edge dense modulation APIs are unchanged.

The accelerator graph surface includes opt-in CUDA Graph and HIP Graph replay
of a complete graph step when every node is qualified. The maintained captured
fixtures cover finite multilayer atmosphere replay, the regular-grid separable
Gaussian DM, pupil-OPD composition, diffractive Shack-Hartmann optics, physical
Pyramid modulation, and the built-in CCD, CMOS, and EMCCD complete-frame detector
nodes. Captured execution does not interleave
unqualified stream nodes.

Dense and measured influence-matrix DM nodes also admit captured execution
with built-in linear-static, clipping, and actuator-health stages or their
static compositions. Preparation retains multiplication coefficients for
replay; stateful and custom actuator stages require separate qualification.

Each prepared CUDA executable retains an independent cuBLAS workspace (4 MiB
before compute capability 9.0, 32 MiB from 9.0 onward). Capture temporarily binds
that workspace to the prepared stream's cuBLAS handle and restores the default
workspace afterward. Preparation rejects allocation and free nodes in the native
CUDA graph, including child graphs, before instantiation. Retained scalar
coefficients alone do not establish fixed device storage. Ordinary stream
execution does not allocate this capture workspace.
The [local CUDA qualification](../benchmarks/results/platform/2026-10-21-cuda-capture-storage.toml)
records capture-node inspection, replay correctness, and plant measurements.

The explicit pixel-MTF operation and graph node enqueue on the prepared device
and stream, with CUDA and AMDGPU capture supported. InGaAs graph capture is
limited to full-frame unit sampling/binning, null persistence and nonlinearity,
and the device counter RNG. These node qualifications do not establish
instrument-specific calibration or complete downstream graph capture.

Prepared graphs provide a bounded capacity-one asynchronous submission
boundary. The graph sequence and caller-visible output readiness publish only
when the completion ticket is consumed. Atmosphere-layer accumulation and
deformable-mirror staging use same-stream ordering on the covered accelerator
paths. The qualified WFS FFT and detector acquisition paths enqueue directly
on the retained device context; ordinary stream execution may still use a
backend's documented detector host-mirror strategy.

`GroupedStreamGraphExecution` is supported on the host, CUDA, and AMDGPU for an
explicit sequence of complete dependency groups. Hardware fixtures compare a
two-branch modal-OPD graph with the single-stream result and
verify distinct retained lane streams. The policy uses one reusable dependency
event per lane and a full device-side barrier between groups. It does not claim
automatic graph partitioning, CPU task parallelism, or a speedup for narrow
graphs.

Support is operation-specific. Preparing a graph for an accelerator succeeds
only when every node and bound array supports that exact target. The lockstep
HIL boundary explicitly copies complete products and commands through host
`Array` buffers.

Optional `AlgorithmGraphs.prepare_graph_output_snapshot` binds one graph output
to a capacity-one host snapshot. The graph owner calls
`try_snapshot_graph_output!` after completion and before the next step. One
consumer calls `poll_graph_output_snapshot!`, reads `snapshot_values` and
`snapshot_metadata`, then calls `release_graph_output_snapshot!`. Metadata is
published with the sample; occupied slots drop observations. CPU copying is
synchronous bounded work. CUDA uses two prepared captured array copies and
reusable events, with mapped host storage retained behind a host array view.
The source never waits for host transfer or consumer completion. Other
accelerators are rejected. Stop both participants before reset or close, which
drain outstanding transfers. Sequence-zero samples permit cold warmup without
claiming an admitted frame. This seam supplies no task, cadence, transport, or
imaging model; timing and allocation claims require qualification with the
actual producer and consumer tasks.
Host snapshots using `FixedSizeArrayDefault` or explicit `Vector` backing can
bind directly to `Optics.PupilFunction` and prepared direct imaging. This host
storage support does not classify accelerator-backed fixed-size arrays as CPU
arrays.

The LiFT `LiFTForwardModel` adapter supplies AOS physical photon-rate prediction
and analytic Jacobians to the registered `AdaptiveOpticsCalibration.PhaseRetrieval`
direct API on CPU, CUDA, and AMDGPU. AdaptiveOpticsCalibration owns all inverse
policies, solves, products, and diagnostics. CPU explicit-workspace inverse
execution remains zero allocated; direct GPU inverse execution allocates and
synchronizes and is outside captured-graph replay. The [live integration
profile](../benchmarks/results/platform/2026-09-22-lift-aoc-adapter-profile.toml)
records these performance boundaries.

`CapturedGraphExecution()` requires every adapter to prove the fixed-address
device replay contract, then records the complete node and delayed-link sequence
as one executable. CUDA and AMDGPU hardware tests replay the two-node DM and
composition graph with changed command contents. Prepared stochastic graph
owners on both backends keep their counter state and normal, uniform, and
Poisson array sampling on the selected device. Hardware tests also compare
captured finite-atmosphere replay with stream execution and verify completed
host epoch publication and reset. They compare captured Shack-Hartmann,
Pyramid, CCD, CMOS, and EMCCD execution with the corresponding stream results.
CUDA and AMDGPU use the same KernelAbstractions SplitMix64 counter kernels and
device draw-sequence model. Fixed-seed reset and replay are reproducible on each
qualified backend. Cross-vendor bitwise floating-point identity is not claimed,
because normal and Poisson transforms may use different vendor transcendental
lowering. Complete downstream instrument graphs require their own whole-graph
capture qualification.

Detector capture is deliberately narrower than the general detector API. It
admits unit-binning, full-frame CCD single-read, global-shutter CMOS, and
linear sequential EMCCD acquisitions with the built-in null response,
coupling, defect, correction, persistence, and thermal models. Prepared device
SplitMix64 state drives photon, dark, CIC, multiplication, and readout noise.
Other detector models and non-unit binning remain on stream execution until
their fixed-address paths are independently qualified.

Metal and AppleAccelerate are optional/manual surfaces and are not implied by
the AMDGPU or CUDA claims.

## Determinism

CPU replay with fixed configuration, fixed explicit seeds, and single-threaded
execution is the deterministic reference surface. Stored historical fixtures
may use `deterministic_reference_rng`; new simulations use `runtime_rng`.
The latter returns `SplitMix64RNG` and does not depend on Julia's task-local
default RNG. Cross-hardware bitwise RNG or floating-point identity is not
claimed.

## Optional Proper.jl Integration

Proper.jl is not a package dependency. The maintained integration pattern is an
application-owned prepared prescription called directly or through an explicit
graph-node adapter. Current tests establish ownership, shape, residency, and
selected CPU/GPU behavior for the companion example; they do not qualify a
SPIDERS prescription or coronagraph model scientifically.

## Optional pyRTC Validation Integration

The example environment contains a native Julia implementation of the numeric
Linux pyRTC `ImageSHM` layout. It validates C-order vector and matrix exchange
in both directions and runs the SHWFS and Pyramid reference systems against a
separate pyRTC process. This lockstep integration and calibration surface also
verifies corrected science-image and on-axis Strehl improvement for a
deterministically evolving four-layer atmosphere. Downstream instrument
packages reuse this adapter and own their full-size pyRTC calibration and
closed-loop acceptance tests. This is not an asynchronous transport or
deadline guarantee. The one-slot upstream layout is used only with one producer
and one outstanding frame.

## Explicitly Not Production-Supported

- the retired `AdaptiveOpticsSim.Plant` event runtime
- a general multi-rate or mid-frame event scheduler
- wall-clock pacing, networking, PipeWire transport, or asynchronous RTC
  protocol handling
- automatic graph partitioning or hidden host/device transfers
- PipeWire GPU-buffer execution
- direct Julia embedding in a PipeWire process
- a scientifically qualified SPIDERS end-to-end prescription
- uncalibrated HSDM277 or SPIDERS/BAX307 influence functions
- automatic camera profiles or unconfirmed optical-stage/filter/chopper values
- arbitrary user graph nodes without their own validation evidence
- a blanket zero-allocation or hard-deadline guarantee

## Support Rule

A new surface becomes supported only when all of these are present:

1. canonical ownership and terminology
2. a documented public or qualified-public API
3. failure and reset behavior
4. deterministic focused CPU tests
5. allocation evidence when the repeated path claims a budget
6. applicable accelerator tests on real hardware
7. scientific reference or instrument evidence appropriate to the claim
8. an entry in [`model-validity-matrix.md`](model-validity-matrix.md)

Historical benchmark artifacts remain evidence for the exact revision and
implementation they record. They do not automatically qualify the current
runtime after a breaking refactor.

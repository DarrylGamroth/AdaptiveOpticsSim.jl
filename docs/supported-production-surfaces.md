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
- Shack-Hartmann, Pyramid, Bi-O-edge, Zernike, Curvature, and LiFT workflows
  within the limits recorded in
  [`model-validity-matrix.md`](model-validity-matrix.md)
- conventional frame and counting/channel detector families covered by their
  registered suites
- calibration matrices/bases, reconstructors, controllers, and tomography
- deterministic and ordinary coarse-grained ensemble execution
- static `AlgorithmGraphs` definitions in Julia and TOML
- direct and explicit one-frame delayed graph links
- exact model-time values and deterministic drivers
- the transport-neutral lockstep `PreparedGraphHILBoundary`
- the built-in graph-node catalog exercised by `algorithm-graphs` tests

Warmed allocation claims apply only to the concrete paths with allocation tests.
They are not a promise for arbitrary user callbacks or every optional library.

## Downstream Instrument Packages

AOS supports the reusable algorithms, graph adapters, static TOML schema, and
lockstep HIL boundary used by downstream instrument packages. It does not own
an instrument's command geometry, detector profile, calibration policy, or
scientific acceptance claim.

[REVOLTClassicSim.jl](https://github.com/DarrylGamroth/REVOLTClassicSim.jl)
and
[REVOLTCopperSim.jl](https://github.com/DarrylGamroth/REVOLTCopperSim.jl)
maintain their executable graph files and pyRTC validation separately. Their
support status is governed by their own tests and documented assumptions.

## Accelerator Surfaces

AMDGPU is the release-gated accelerator for the operations exercised by the
maintained local ROCm target. CUDA is manually validated on the WSL RTX host but
is not a continuously available release gate. Both policies require:

- scalar indexing disabled
- exact device residency for inputs, plans, state, workspace, and products,
  except for an operation's documented preallocated host-mirror strategy
- no undeclared CPU fallback
- explicit synchronization at public completion boundaries
- numerical or statistical parity against the declared CPU reference

The accelerator graph surface includes opt-in CUDA Graph and HIP Graph replay
of a complete graph step when every node is qualified. The maintained captured
fixtures cover finite multilayer atmosphere replay, the regular-grid separable
Gaussian DM, pupil-OPD composition, diffractive Shack-Hartmann optics, both
maintained Pyramid modulation strategies, and the built-in CCD, CMOS, and
EMCCD complete-frame detector nodes. Captured execution does not interleave
unqualified stream nodes.

Prepared graphs provide a bounded capacity-one asynchronous submission
boundary. The graph sequence and caller-visible output readiness publish only
when the completion ticket is consumed. Atmosphere-layer accumulation and
deformable-mirror staging use same-stream ordering on the covered accelerator
paths. The qualified WFS FFT and detector acquisition paths enqueue directly
on the retained device context; ordinary stream execution may still use a
backend's documented detector host-mirror strategy.

`GroupedStreamGraphExecution` is supported on the host, CUDA, and AMDGPU for an
explicit sequence of complete dependency groups. Hardware fixtures compare a
two-branch controller/modal-expansion graph with the single-stream result and
verify distinct retained lane streams. The policy uses one reusable dependency
event per lane and a full device-side barrier between groups. It does not claim
automatic graph partitioning, CPU task parallelism, or a speedup for narrow
graphs.

Support is operation-specific. Preparing a graph for an accelerator succeeds
only when every node and bound array supports that exact target. The lockstep
HIL boundary explicitly copies complete products and commands through host
`Array` buffers.

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

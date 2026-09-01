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

## REVOLT Graph Examples

The REVOLT Classic and REVOLT Copper TOML files are maintained executable graph
examples. They cover atmosphere evolution, the 277-element HSDM command map,
pupil OPD composition, Shack-Hartmann or Pyramid optical formation, and
completed detector frames. Full-optical and fast separable-grid DM profiles are
available.

These examples are not a claim of instrument qualification. The HSDM influence
model is provisional until measured HSDM277 calibration is supplied. Their
configuration comments identify external REVOLT/SPECULA/OOPAO sources and
remaining assumptions.

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
fixture is the regular-grid separable Gaussian DM followed by pupil-OPD
composition; captured execution does not interleave unqualified stream nodes.

Prepared graphs provide a bounded capacity-one asynchronous submission
boundary. The graph sequence and caller-visible output readiness publish only
when the completion ticket is consumed. Atmosphere-layer accumulation and
deformable-mirror staging use same-stream ordering on the covered accelerator
paths. WFS FFT implementations and the documented AMDGPU detector host mirror
may retain narrower required completion boundaries inside submission.

Support is operation-specific. Preparing a graph for an accelerator succeeds
only when every node and bound array supports that exact target. The lockstep
HIL boundary explicitly copies complete products and commands through host
`Array` buffers.

`CapturedGraphExecution()` requires every adapter to prove the fixed-address
device replay contract, then records the complete node and delayed-link sequence
as one executable. CUDA and AMDGPU hardware tests replay the two-node DM and
composition graph with changed command contents. Prepared stochastic graph
owners on both backends keep their counter state and normal, uniform, and
Poisson array sampling on the selected device. Complete REVOLT graphs are not
yet capture-qualified because atmosphere, WFS, and detector owners still need
independent capture audits.

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
deterministically evolving four-layer atmosphere. An independently enabled
REVOLT Classic model-validation test covers its 352-by-352 frame, 376
pair-interleaved slopes, 277-element PDM command, and maintained five-layer
atmosphere for one modeled second. It admits a regularized interaction
subspace and checks both Strehl and pupil OPD-RMS improvement. This does not
qualify the provisional HSDM influence model or an operational
control/extrapolation matrix. It is not an asynchronous transport or deadline
guarantee. The one-slot upstream layout is used only with one producer and one
outstanding frame.

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

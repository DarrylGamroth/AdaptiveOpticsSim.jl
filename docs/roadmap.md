# Roadmap

Status: active

This roadmap is intentionally short. Historical implementation plans were
removed from the live docs tree during the April 2026 documentation cleanup; use
git history if old plan details are needed.

## Current State

AdaptiveOpticsSim.jl is usable for maintained CPU workflows and selected
hardware-validated GPU workflows. The core package now has:

- typed AO model objects for optics, atmosphere, WFS, detectors, DMs, control,
  and runtime orchestration
- CPU tests across Linux, macOS, and Windows
- an active AMDGPU release gate and a retained CUDA hardware target with current
  manual WSL evidence; CUDA remains outside the release gate
- a maintained zero-allocation hot-path focus for CPU HIL-style runtime loops
- explicit coarse ensemble scheduling for independent plants, with optional
  AcceleratedKernels and Dagger integrations kept outside the HIL inner loop
- caller-owned calibration storage plus compact factorized and controller-
  composed reconstruction operators for validated large control surfaces
- measured column-major/native-SIMD CPU kernels for LiFT convolution and LGS
  elongation, without adding an explicit SIMD dependency
- committed reference data and regression tests for core WFS and detector paths
- a completed schedule-free Gate 2 plant boundary with stable path/acquisition
  ownership, per-owner RNGs, fidelity-provider and calibration-illumination
  seams, and a clean serial CPU service-time baseline
- a completed Gate 3 deterministic multi-rate virtual-time engine with one
  canonical plant timeline, bounded trigger fan-out and faults, conventional
  detector lifecycles, fixed prepared storage, and clean scheduler/composed-
  plant CPU evidence
- a compact docs set with one extension guide instead of subsystem plan sprawl

## Near-Term Priorities

`AdaptiveOpticsSim.jl` is being developed as the high-fidelity AO plant for
external-RTC HIL development, following the maintained specifications indexed
by [`hil-package-boundary.md`](hil-package-boundary.md) and tracking completion
in [`hil/compliance-matrix.md`](hil/compliance-matrix.md).

1. Preserve the completed HIL prerequisite Gates 0 and 1 while implementing
   the proposed general HIL runtime. Gate 0 separates telescope aperture/geometry
   from caller-owned optical
   planes and propagation workspaces; separate shared atmosphere evolution from
   path-local NGS/LGS/source rendering; remove temporal cadence from the
   telescope; define plane geometry/radiometry and safe spectral combination;
   and separate direct-science photon-arrival-rate formation from detector-owned
   temporal integration and acquisition. Then decompose every maintained
   WFS into a prepared optical front end, detector acquisition, and estimator.
   Shack-Hartmann now has an independent microlens array; Pyramid/BioEdge have
   separate physical optics over shared modulation; Zernike/Curvature now
   separate propagation, acquisition, and estimation, including independent or
   packed Curvature detector planes; and LiFT now consumes independently
   acquired, explicitly normalized observations through a separately prepared
   focal-plane model. Cross-backend correctness and residency evidence is now
   complete on the maintained CPU, CUDA, and AMDGPU targets, with clean CPU and
   CUDA service-time artifacts retained. The final Gate 0 ownership review
   removes telescope-owned mutable optical-path state: each maintained WFS,
   science, calibration, atmosphere, and controllable-optic path now consumes
   an explicit `PupilFunction` or field product. Preserve CPU, CUDA, and AMDGPU
   correctness, residency, allocation, and latency evidence throughout the HIL
   migration. Gate 1 freezes the breaking plant-oriented API, package/type
   boundaries, atmosphere token/materialization lifetime, deterministic RNG
   ownership, detector event semantics, clock sequencing, and command boundary
   before implementation begins.
2. Preserve the completed Gate 2 schedule-free plant boundary. It composes
   immutable shared atmosphere/telescope/path definitions, prepared branch-
   local executors, independent acquisition owners, stable per-owner RNGs, and
   full-optical/reduced-order/synthetic provider semantics. Native and
   user-defined calibration illumination enter through typed products without
   instrument-topology assumptions. The clean
   [serial plant artifact](../benchmarks/results/gate2/2026-07-21-serial-plant.toml)
   covers science, NGS Shack-Hartmann, and LGS pyramid directions plus detector
   fan-out with zero warmed allocation; it is a self-paced CPU service-time
   baseline, not an external-RTC latency or fixed-rate capacity claim.
3. Preserve the completed deterministic multi-rate integer-time engine with
   explicit equal-time trigger-distribution, exposure/row-band, optical-sample, nondestructive-read,
   detector-readout, and publication semantics before adding command timing or
   wall-clock pacing. Canonical time, the fixed-capacity event calendar, trigger
   distribution, exact global/rolling/frame-transfer lifecycles, evolving-charge
   HgCdTe ramp reads, and their common serial scheduler composition are
   implemented and validated. The clean [scheduler](../benchmarks/results/gate3/2026-07-21-event-scheduler-gate3-closure.toml)
   and [composed multi-rate plant](../benchmarks/results/gate3/2026-07-21-multi-rate-plant.toml)
   artifacts close the gate without claiming wall-clock pacing, external-RTC
   latency, or production instrument capacity. Keep physical trigger faults
   separate from timestamp-label faults and execution lateness.
4. Replace the single-optic and `CompositeControllableOptic` runtime model with
   individually placed optics, prepared core plant command schemas, bounded
   timing and replayable plant-time command-silence semantics, sampled device-
   feedback acquisitions, and prepared plane groups as a deliberate breaking
   change. The first two Gate 4 slices now record stable physical-optic and
   independently latched endpoint identities plus immutable versioned semantic
   payload schemas in `PlantDefinition`, while failing preparation explicitly
   until mutable endpoint owners are added. Operational execution-clock ingress
   liveness belongs to the later HIL lifecycle boundary.
5. Immediately prove a minimal serial CPU HIL vertical slice: one scheduled
   acquisition, one command-responsive optic, an injected `Clocks.jl` clock,
   HIL submission descriptors mapped into core plant commands, canonical
   complete-product and command/outcome ports, bounded SPSC/lease ownership, a
   deterministic fake RTC, and fixed-arrival evidence. Do this before worker,
   GPU, transport-specific, or placement-planner complexity.
6. Add path-specific NCPA, MCAO/MOAO geometry, prepared CPU execution groups,
   then single-GPU direction batching and physical device identity behind
   numerical, allocation, residency, and fixed-arrival evidence.
7. Harden the transport-neutral HIL companion with lifecycle transitions,
   guaranteed lease-return credit, first-failure propagation, replay classes,
   and GC/process-isolation policy. Add explicit or constrained deterministic
   mixed CPU/GPU placement, then homogeneous multi-GPU placement. Defer a fully
   automatic cost-model planner until real profiles provide calibration data;
   keep Dagger and dynamic migration outside the HIL deadline path.
8. Preserve hardware validation and zero-allocation CPU gates, then use pinned
   NFIRAOS and MORFEO companion scenarios for synchronized multi-rate and
   extreme-scale profiles. Give each a production-shaped synthetic traffic
   variant, a reduced-order closed-loop variant where applicable, and a full
   optical variant while keeping topology, model, timing, and external-
   integration compliance independent.

## Active Cleanup Themes

- Keep the README and `user-guide.md` focused on one recommended user path.
- Keep `api-reference.md` aligned with the exported API instead of documenting
  every internal or qualified helper.
- Keep full visual examples in `../AdaptiveOpticsSimPlots.jl`; the core
  package examples remain plotting-free and runnable through
  [`run_core_examples.sh`](../scripts/run_core_examples.sh).
- Prefer a few high-value OOPAO/SPECULA/REVOLT-like equivalence artifacts over
  broad claims that are not release-gated.
- Consolidate validation around maintained entry points:
  `Pkg.test()`, backend-specific hardware targets, the core example runner, and
  the release-validation script.

## Deferred Areas

The following remain valid future directions, but should not drive ad hoc API
growth:

- broader manufacturer-specific DM technology models
- richer detector physical models where reusable readout/thermal/counting layers
  are insufficient
- wider cross-package numerical equivalence beyond the maintained reference data
  and artifacts
- companion visualization and analysis packages outside the core package
- science-path integrations that belong in optional extensions or sibling
  packages

## Documentation Rule

Do not add new one-off roadmap fragments. Update this file, the
[`documentation-map.md`](documentation-map.md), or the relevant maintained guide
instead.

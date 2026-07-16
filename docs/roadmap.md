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
- a compact docs set with one extension guide instead of subsystem plan sprawl

## Near-Term Priorities

`AdaptiveOpticsSim.jl` is being developed as the high-fidelity AO plant for
external-RTC HIL development, following the maintained specifications indexed
by [`hil-package-boundary.md`](hil-package-boundary.md) and tracking completion
in [`hil/compliance-matrix.md`](hil/compliance-matrix.md).

1. Complete HIL prerequisite Gate 0 before implementing the proposed general
   HIL runtime: separate telescope aperture/geometry from caller-owned optical
   planes and propagation workspaces; separate shared atmosphere evolution from
   path-local NGS/LGS/source rendering; remove temporal cadence from the
   telescope; define plane geometry/radiometry and safe spectral combination;
   and separate direct-science photon-arrival-rate formation from detector-owned
   temporal integration and acquisition. Then decompose every maintained
   WFS into a prepared optical front end, detector acquisition, and estimator.
   Extract the Shack-Hartmann microlens array; separate Pyramid/BioEdge
   modulation and focal-plane optics; separate Zernike/Curvature propagation
   and acquisition, including multiple Curvature detector planes; and let LiFT
   consume independently acquired observations. Freeze current results first
   and preserve CPU, CUDA, and AMDGPU correctness, residency, allocation, and
   latency evidence throughout the twelve-PR migration.
2. Compose the Gate 0 optical ownership primitives into immutable shared
   atmosphere/telescope/path definitions, prepared branch-local executors, and
   independently scheduled or triggered acquisition state while retaining the
   direct serial CPU oracle. Introduce one prepared
   acquisition-product seam for full optical, command-responsive reduced-order,
   and synthetic/replay providers without changing the RTC boundary. The
   reduced-order provider must retain time-correlated disturbances, calibrated
   path/sensor response, effective-command causality, and meaningful external-
   RTC loop closure rather than acting as a changing frame generator. Expose a
   narrow typed path-entry seam for user-defined calibration illumination
   without introducing instrument topology or source assumptions in core.
3. Add deterministic multi-rate integer-time events with explicit equal-time
   command, trigger-distribution, optical-sample, detector-readout, and
   publication semantics before adding wall-clock pacing. Keep physical
   trigger faults separate from timestamp-label faults and execution lateness.
4. Replace the single-optic and `CompositeControllableOptic` runtime model with
   individually placed optics, prepared command schemas, bounded timing and
   silence/watchdog semantics, sampled device-feedback acquisitions, and
   prepared plane groups as a deliberate breaking change.
5. Immediately prove a minimal serial CPU HIL vertical slice: one scheduled
   acquisition, one command-responsive optic, an injected clock, canonical
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

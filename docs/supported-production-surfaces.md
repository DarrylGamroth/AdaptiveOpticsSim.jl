# Supported Production Surfaces

Status: active

## Purpose

This note defines the current production-supported scope for
AdaptiveOpticsSim.jl.

The intent is to be explicit about which surfaces are defended by maintained
evidence, backend parity checks, and real-hardware validation, and which
surfaces are outside the current support claim.

This is a support-policy note, not a claim that every exported symbol or model
family has identical maturity.

## Support Rule

A surface should only be treated as production-supported when all of the
following are true:

1. the modeling/workflow surface is documented and maintained
2. functional correctness is covered in the normal test suite
3. backend parity is covered where accelerator support is claimed
4. realistic runtime evidence exists where performance or HIL claims are made
5. external equivalence is either demonstrated or explicitly scoped

Use these supporting docs together:

- [model-validity-matrix.md](model-validity-matrix.md)
- [backend-validation-guide.md](backend-validation-guide.md)
- [release-validation-runbook.md](release-validation-runbook.md)

## Production-Supported Surfaces

### CPU baseline

The primary production baseline is the typed Julia CPU path on maintained
runtime and validation surfaces.

Current CPU-supported families:

- telescope-owned aperture geometry and reflectivity with caller-owned,
  independent `PupilFunction` path products; the telescope owns no cadence,
  path OPD, field, PSF, or propagation scratch
- finite and infinite multilayer atmosphere with explicit elapsed/absolute
  model time, current-state epoch tokens, frozen direction preparation, and
  caller-owned rendering
- detector-family execution on maintained detector surfaces
- staged Shack-Hartmann, Pyramid, and BioEdge WFS optical formation,
  detector acquisition, and estimation on maintained validated surfaces;
  geometric variants use explicit direct-measurement paths
- runtime and closed-loop execution on maintained local/runtime artifacts
- direct `CPUHILExecutionPlan` execution as the default CPU runtime residency
  policy
- source-aware runtime propagation with distinct WFS and science sources
- allocation-free shared multi-arm CPU execution with one atmosphere advance
  and one caller-owned same-arm photon-arrival-rate image reused by independent
  detector acquisitions
- the schedule-free acquisition product-provider boundary: run-immutable
  full-optical or nonresponsive unchanged/copy/bounded-replay selection,
  invariant caller-owned product contracts, unused-path bypass, and zero
  warmed provider allocation. This is boundary-semantic support, not an
  external-RTC performance or synthetic optical-fidelity claim
- the schedule-free calibration-illumination entry boundary: typed
  pupil/field/intensity/external-result/detector-input products, explicit
  visibility and combination declarations, exact caller-owned destinations,
  stable time/RNG ownership, ordinary downstream path/acquisition execution,
  and zero warmed evaluator allocation. This is boundary-semantic support, not
  a claim about a physical calibration unit or scheduled control protocol
- grouped/runtime orchestration surfaces with committed validation artifacts
- dense and factorized reconstruction operators, including allocation-free
  controller composition and backend-residency validation

Optional AK and Dagger ensemble policies are maintained offline orchestration
surfaces, but are not part of the CPU HIL latency claim. Scheduler choice is
evidence-gated per host and workload; direct `CPUHILExecutionPlan` stepping
remains the production baseline.

Primary evidence:

- [model-validity-matrix.md](model-validity-matrix.md)
- the workflows in [`user-guide.md`](user-guide.md) and [`release-validation-runbook.md`](release-validation-runbook.md)
- benchmark artifacts under `benchmarks/results/`
- [conventional-detector CPU HIL latency baseline](../benchmarks/results/detectors/2026-07-14-detector-hil-latency.toml)
- [final pre-HIL local CPU service-time evidence](../benchmarks/results/platform/2026-07-18-pre-hil-11-local-cpu.toml)
- [final pre-HIL WSL CPU service-time evidence](../benchmarks/results/platform/2026-07-18-pre-hil-11-wsl-cpu.toml)
- [../benchmarks/results/validation_runs/2026-04-10-rtc-devel-cpu.toml](../benchmarks/results/validation_runs/2026-04-10-rtc-devel-cpu.toml)

### AMDGPU backend

AMDGPU is a production-supported accelerator backend on the maintained surfaces
covered by the dedicated hardware validation target:

- [../test/runtests_amdgpu.jl](../test/runtests_amdgpu.jl)

Current AMDGPU-supported scope:

- maintained optional backend functional/parity checks
- maintained runtime-equivalence surfaces
- default `DeviceResidentExecutionPlan` runtime construction with validated
  accelerator-resident reconstructor storage
- maintained high-accuracy post-command runtime-equivalence surfaces
- maintained Shack-Hartmann exported-pixel parity surfaces
- maintained composite low-order runtime surfaces
- matched HEART RTC HIL runtime surfaces
- prepared direct imaging with off-axis formation, spectral bundles, explicit
  extended-source expansion, independent detector fan-out, and shared-arm
  device residency
- schedule-free `PreparedPlant` direct-science execution with one device-
  resident optical result shared by independent unequal-exposure acquisitions
- schedule-free native uniform calibration illumination entering a
  device-resident detector-input path and passing through ordinary detector
  acquisition with scalar indexing disabled

Primary evidence:

- [backend-validation-guide.md](backend-validation-guide.md)
- [release-validation-runbook.md](release-validation-runbook.md)
- benchmark artifacts under `benchmarks/results/`
- [current CPU/CUDA/AMDGPU cross-host characterization](../benchmarks/results/platform/2026-07-14-wsl-cuda-local-amdgpu.toml)
- [../benchmarks/results/validation_runs/2026-04-10-rtc-devel-amdgpu.toml](../benchmarks/results/validation_runs/2026-04-10-rtc-devel-amdgpu.toml)

Current expectation:

- if a maintained AMDGPU surface regresses numerically against CPU, that is a
  release-blocking defect for the AMDGPU-supported scope

The current Julia 1.12.6 AMDGPU hardware target passed all `411` maintained
checks. A later local Julia installation failure prevented a replacement raw
latency artifact, so the July 14 characterization remains the maintained AMD
performance evidence; the failed host run does not broaden or weaken the
functional support claim.

### GPU support-boundary rule

GPU support is defined by the maintained dedicated hardware targets and the
release-validation path.

The package does still contain broader backend-audit and subsystem-investigation
scripts, but those do not by themselves define supported GPU scope. A GPU-touched
surface is only support-claimed when it is also promoted into the maintained
hardware targets and release-validation cadence.

### OOPAO-aligned external equivalence

The strongest current external production-alignment story is the OOPAO-matched
HEART Shack-Hartmann surface.

Current production-supported external equivalence surface:

- HEART RTC HIL boundary:
  - `277` DM commands in
  - `352x352` Shack-Hartmann frame out
- matched OOPAO/Julia baseline on the comparison-owned flux-threshold validity
  convention

Primary evidence:

- committed OOPAO reference data under `test/reference_data`
- `../AdaptiveOpticsComparisons/docs/cross-package-benchmark-harness.md`
- `../AdaptiveOpticsComparisons/results/archived/2026-04-08-heart-hil-oopao-baseline.toml`

Additional production-supported frozen OOPAO equivalence surfaces:

- diffractive Pyramid ramp from the committed OOPAO reference bundle
- diffractive BioEdge ramp from the committed OOPAO reference bundle
- prepared Pyramid/BioEdge pupil-function and electric-field formation,
  spectral and path-local source bundles, LGS elongation/sodium profiles,
  exact-once detector exposure, and revision-bound differential estimation

Primary evidence:

- committed OOPAO reference data under `test/reference_data`
- [../scripts/generate_oopao_equivalence_artifact.jl](../scripts/generate_oopao_equivalence_artifact.jl)
- [../benchmarks/results/equivalence/2026-04-09-oopao-production-equivalence.toml](../benchmarks/results/equivalence/2026-04-09-oopao-production-equivalence.toml)

Scientist-owned HEART boundary truth artifact:

- [../scripts/generate_heart_boundary_truth_artifact.py](../scripts/generate_heart_boundary_truth_artifact.py)
- [../benchmarks/results/truth/2026-04-09-heart-boundary-truth.toml](../benchmarks/results/truth/2026-04-09-heart-boundary-truth.toml)

## Explicitly Not Yet Production-Supported

The following are outside the current support claim:

- CUDA execution. The extension, dedicated test project, fail-fast `412/412`
  hardware target, and current manual WSL hardware evidence—including prepared direct
  imaging and the
  [final pre-HIL CUDA artifact](../benchmarks/results/platform/2026-07-18-pre-hil-11-wsl-cuda.toml)—are
  available, but CUDA has not yet been explicitly returned to the supported
  delivery scope or a routine validation cadence.
- SPECULA pixel-level equivalence on the HEART Shack-Hartmann surface
- Metal backend support
- backend-audit surfaces that are not part of the maintained hardware targets
  and release-validation cadence
- broad claims about every detector/wfs/backend combination outside the
  maintained evidence surfaces
- cross-package grouped/platform equivalence beyond the currently normalized
  contracts
- full optical or on-sky instrument-truth alignment beyond the maintained boundary artifact
- a production reduced-order AO surrogate. The maintained test extension
  proves provider dispatch and command causality only; it does not establish a
  validated disturbance/command/sensor envelope, scheduled loop closure,
  tomography, stability, or optical accuracy
- RTC latency, capacity, transport, cache-residency, or production-shaped
  payload claims inferred from unchanged, copied, or cyclic replay providers
- physical fidelity claims for a user-defined calibration source, relay,
  coherence/spectral prescription, insertion mechanism, or instrument profile;
  the maintained illumination surface validates the typed execution boundary
  only
- scheduled calibration triggers, setpoint/control authority, HIL descriptors,
  or calibration transport protocols

## Release Interpretation

For current release/readiness decisions:

- regressions on the CPU baseline are release-blocking
- regressions on the maintained AMDGPU surfaces are release-blocking
- CUDA becomes release-blocking only after it is explicitly returned to the
  supported delivery scope and assigned a routine validation cadence
- unresolved SPECULA differences are not release-blocking unless the package
  starts claiming SPECULA equivalence for that surface
- new features outside this supported surface should not be described as
  supported until they gain the same evidence shape
- supported accelerator claims require routine validation on real hardware,
  either through CI or through the maintained release-validation host path

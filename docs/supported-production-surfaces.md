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
numerical, control, and Plant validation surfaces.

Current CPU-supported families:

- the canonical `AdaptiveOpticsSim.Plant` namespace for HIL-neutral plant
  definitions, virtual time, commands, acquisition events, providers,
  preparation, and serial event composition; root-level compatibility aliases
  are not supported
- telescope-owned aperture geometry and reflectivity with caller-owned,
  independent `PupilFunction` path products; the telescope owns no cadence,
  path OPD, field, PSF, or propagation scratch
- finite and infinite multilayer atmosphere with explicit elapsed/absolute
  model time, current-state epoch tokens, frozen single-direction preparation,
  and homogeneous ordered direction batching into caller-owned OPD storage;
  CPU batching preserves the serial oracle while accelerator batching keeps
  geometry, pupil support, layer screens, and output on one concrete device
- detector-family execution on maintained detector surfaces
- staged Shack-Hartmann, Pyramid, and BioEdge WFS optical formation,
  detector acquisition, and estimation on maintained validated surfaces;
  geometric variants use explicit direct-measurement paths
- explicit closed-loop numerical composition and the model-specific AO188/AO3k
  example/profile surfaces
- independent control primitives, dense/factorized reconstruction, and
  caller-owned command storage
- source-aware prepared paths with distinct WFS and science sources
- schedule-free serial `PreparedPlant` execution across independent science,
  NGS Shack-Hartmann, and finite-height LGS pyramid directions, with canonical
  selected ownership, deterministic declaration-order replay, one shared
  science result feeding unequal-exposure detectors, and zero warmed allocation
- deterministic serial `PreparedPlantEventLoop` execution across independently
  periodic science/NGS/LGS paths and periodic or delivered-trigger detector
  starts, including global CMOS/CCD, rolling CMOS, frame-transfer EMCCD, and
  HgCdTe up-the-ramp lifecycles. The scheduler and warmed direct detector
  kernels are allocation-free; heterogeneous orchestration has a bounded
  2 KiB-per-processed-timestamp CPU allocation budget. This is virtual-time
  correctness support, not wall-clock HIL latency, fixed-arrival capacity,
  parallel placement, or production instrument-scale evidence
- deterministic linear modal reduced-order direct measurements with explicit
  harmonic disturbance, path, command-response, and sensor operators; exact
  command-schema and path-visible endpoint binding; periodic or
  delivered-trigger acquisition timing; matched loop closure; declared
  residual metric, operating envelope, and omitted effects; common versus
  target-local response isolation; and unused full-optical path bypass.
  Support is limited to the maintained CPU validity envelope in `MV-26`
- native Plant deformable mirrors with explicit actuator topology, influence
  and actuator models, device-internal misregistration, path/relay
  registration, independent active/staged storage, and allocation-free warmed
  surface publication/application. Serial CPU composition covers common DMs at
  the pupil and multiple atmospheric conjugates over NGS, finite-height LGS,
  and science directions plus target-local MOAO DMs. This is composition
  support, not calibrated device dynamics, tomography fidelity, external
  optical equivalence, fixed-arrival latency, or instrument-scale capacity
- integrated serial CPU conjugated optical placement across alternating
  finite-height-LGS/NGS WFS and science paths, with common pupil/multi-altitude
  DMs, selected-path MOAO, common static OPD, science-only NCPA, exact
  declaration-order replay, finite-support loss, stable prepared storage, and a
  bounded warmed allocation window. The maintained
  [Gate 5 artifact](../benchmarks/results/gate5/2026-07-25-optical-placement.toml)
  characterizes one eight-path low-resolution workload and exact 4/8/16-path
  binding counts; it is not fixed-arrival, parallel, integrated-GPU,
  detailed-relay/coronagraph, or production instrument-capacity evidence
- prepared CPU path-execution groups with one writer per path/acquisition/RNG
  owner, independently callable mutating lifecycle seams, explicit
  Julia/FFT/BLAS execution budgets, a deterministic serial oracle, and
  topology-size-invariant whole-plant registries. Maintained one- and
  four-thread tests prove exact grouped/serial replay and actual coarse path
  overlap; clean artifacts bound direct-call allocation, paired self-paced
  service cost, preparation, first use, storage, inference, method instances,
  and native-code growth. This supports the core executor boundary, not a
  production worker implementation, affinity/pacing policy, fixed-arrival HIL
  capacity, integrated GPU event loop, or NFIRAOS/MORFEO qualification
- the preparation-only caller-resolved partition boundary for CPU-only,
  one-accelerator-only, and CPU-plus-one-accelerator resources. Complete
  path/acquisition groups remain co-located, one authority owns the timed
  atmosphere/timeline/RNG state, stable RNG metadata is placement invariant,
  and every partition carries exact target-scoped structural facts. Each
  visible logical controllable optic also has one shared role-neutral physical
  owner and independent active/staging effective-command storage on that exact
  target. One placement-neutral run identity names the future sole command
  publisher; core does not bind or execute that authority. Maintained
  CPU/fake-resource tests and WSL CUDA hardware validation defend topology,
  authority identity, exact residency, static sampled-OPD materialization,
  deterministic local DM application, alias isolation, and resource accounting.
  The AMDGPU optional checks exist, but this public path is not promoted on
  AMDGPU. The no-command partition reproducer terminates the local `gfx1030`
  process under normal Julia 1.12.6/AMDGPU.jl 2.7.0 compilation, while its
  manual staged construction and `--compile=min` controls complete.
  [Issue #216](https://github.com/DarrylGamroth/AdaptiveOpticsSim.jl/issues/216)
  tracks that reproducer. The command-aware check separately terminates in
  GPUCompiler while compiling finite-array validation even with
  `--compile=min`; [issue #213](https://github.com/DarrylGamroth/AdaptiveOpticsSim.jl/issues/213)
  tracks that path. A queued GitHub AMDGPU run with no eligible runner is
  infrastructure status, not validation evidence. Cold preparation may copy
  declared static data and initial command
  values into an exact target; this surface does not execute mixed work,
  publish atmosphere state, bind the command publisher, transfer a command
  payload, hand runtime products between partitions, choose placement, admit
  capacity, or establish a HIL latency claim
- trigger-relative autonomous circular-Pyramid modulation with bounded
  radius/frequency/phase/enabled setpoints, free-running/source/delivered-reset
  relationships, deterministic branch faults, and an allocation-free
  cycle-averaged CPU optical update. This does not claim time-resolved
  modulation or a physical steering-mirror servo
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
- prepared controller-output routing from named caller-owned products to exact
  independent Plant endpoints, including zero-copy views and zero-allocation
  warmed payload access. This is a preparation/ownership contract, not an RTC
  transport or latency claim
- dense and factorized reconstruction operators, including allocation-free
  controller composition and backend-residency validation

Optional AK and Dagger ensemble policies are maintained offline orchestration
surfaces, but are not part of the CPU HIL latency claim. The deterministic
serial Plant event loop remains the HIL-neutral timing oracle.

On Apple Silicon, application-owned AppleAccelerate 0.7 BLAS/LAPACK selection
is a maintained optional CPU-provider variant exposed through a weak-dependency
extension; it is not a core dependency or a new array backend. The dedicated
hosted target proves backend-neutral package load, active Accelerate forwarding
under explicit opt-in, allocation-free vDSP FFT plans for supported full
power-of-two 1D/2D transforms, FFTW fallback outside that capability envelope,
and the full CPU correctness surface.

Primary evidence:

- [model-validity-matrix.md](model-validity-matrix.md)
- the workflows in [`user-guide.md`](user-guide.md) and [`release-validation-runbook.md`](release-validation-runbook.md)
- benchmark artifacts under `benchmarks/results/`
- [conventional-detector CPU HIL latency baseline](../benchmarks/results/detectors/2026-07-14-detector-hil-latency.toml)
- [minimal deterministic frame-detector CPU service-cost evidence](../benchmarks/results/detectors/2026-08-01-shared-low-fidelity-service-cost.toml)
- [detector family qualification closure catalog](../benchmarks/results/detectors/2026-08-01-detector-qualification-closure.toml)
- [final pre-HIL local CPU service-time evidence](../benchmarks/results/platform/2026-07-18-pre-hil-11-local-cpu.toml)
- [final pre-HIL WSL CPU service-time evidence](../benchmarks/results/platform/2026-07-18-pre-hil-11-wsl-cpu.toml)
- [Gate 2 serial plant CPU service-time evidence](../benchmarks/results/gate2/2026-07-21-serial-plant.toml)
- [Gate 3 scheduler CPU evidence](../benchmarks/results/gate3/2026-07-21-event-scheduler-gate3-closure.toml)
- [Gate 3 composed multi-rate CPU evidence](../benchmarks/results/gate3/2026-07-21-multi-rate-plant.toml)
- [Gate 4 command-responsive plant CPU evidence](../benchmarks/results/gate4/2026-07-24-command-plant.toml)
- [Gate 6 grouped CPU evidence](../benchmarks/results/gate6/2026-07-25-grouped-cpu.toml)
- [Gate 6 topology-growth evidence](../benchmarks/results/gate6/2026-07-25-topology-growth.toml)
- [../benchmarks/results/validation_runs/2026-04-10-rtc-devel-cpu.toml](../benchmarks/results/validation_runs/2026-04-10-rtc-devel-cpu.toml)

### AMDGPU backend

AMDGPU is a production-supported accelerator backend on the maintained surfaces
covered by the dedicated hardware validation target:

- [../test/runtests_amdgpu.jl](../test/runtests_amdgpu.jl)
- [../test/runtests_amdgpu_detectors.jl](../test/runtests_amdgpu_detectors.jl)

Current AMDGPU-supported scope:

- maintained optional backend functional/parity checks
- accelerator-resident reconstructor and controller storage
- maintained Shack-Hartmann exported-pixel parity surfaces
- independent DM and modal/low-order optic application
- native Plant DM preparation with device-resident active/staged command and
  surface storage plus direct surface parity
- prepared controller-output routing with device-resident views and exact
  backend/device validation
- dynamic cycle-averaged circular-Pyramid radius/enable updates with
  device-resident quadrature storage
- maintained REVOLT-like production-shaped WFS frame smoke
- prepared direct imaging with off-axis formation, spectral bundles, explicit
  extended-source expansion, independent detector fan-out, and shared-arm
  device residency; compatible physical native-Fraunhofer science samples may
  additionally use one fixed stacked field/output allocation and FFT plan
  while retaining ordered per-sample product views
- schedule-free `PreparedPlant` direct-science execution with one device-
  resident optical result shared by independent unequal-exposure acquisitions
- schedule-free native uniform calibration illumination entering a
  device-resident detector-input path and passing through ordinary detector
  acquisition with scalar indexing disabled
- exact global-shutter detector-event accumulation and scheduled windowed
  HgCdTe up-the-ramp snapshots/fitting with device-resident products and scalar
  indexing disabled
- HgCdTe linear-avalanche direct capture and scheduled up-the-ramp execution
  with the explicitly approximate clipped-Gaussian multiplication policy,
  moderate-charge distributional bounds, nonnegativity, device residency, and
  scalar indexing disabled. Conditional-Gamma multiplication remains the CPU
  reference and is rejected on accelerators
- bounded Skipper CCD independent-read sampling with a device-resident
  streaming mean, fixed frame-sized storage, distributional read-noise
  reduction, explicit sample duration, and scalar indexing disabled
- direct rolling-shutter row-band and frame-transfer image/storage detector
  lifecycles with device-resident state and scalar indexing disabled
- integrated single-device ownership for two or more compatible, equally
  clocked, direct native-Fraunhofer science paths through one shared
  atmosphere-direction batch and one direct-imaging batch. The owner retains
  the exact compute-device context and stream, preserves path input/result
  identities through same-device handoffs, and synchronizes before reporting
  completion
- integrated single-device ownership for two or more compatible, equally
  clocked Shack-Hartmann, Pyramid, or BioEdge WFS paths. The owner shares the
  exact device context and atmosphere-direction batch, then invokes each
  family's existing prepared lenslet/modulation/spectral pipeline and preserves
  the original path-local result. Singleton, unequal-rate/origin, mixed-family,
  incompatible-signature/product, Zernike, Curvature, and cross-device groups
  continue through their ordinary execution paths
- a maintained six-row conventional detector matrix spanning CCD single-read,
  frame-transfer EMCCD, global- and rolling-shutter CMOS, and global-shutter
  HgCdTe single-read and scheduled up-the-ramp operation. Tests preserve
  response/MTF metadata, prove response-before-exposure ordering against an
  independent zero-extended-convolution oracle, prove post-collection IPC
  separately from presampling MTF, and retain exact transition/readout/storage
  state on device
- deterministic InGaAs persistence and residency, linear-APD channel order,
  SPAD dead-time/redistribution/conversion, and MKID parity/passband handling;
  stochastic InGaAs and SPAD Poisson checks remain explicit AMDGPU non-claims
  after an AMDGPU/GPUCompiler compilation crash

Primary evidence:

- [backend-validation-guide.md](backend-validation-guide.md)
- [release-validation-runbook.md](release-validation-runbook.md)
- benchmark artifacts under `benchmarks/results/`
- [current CPU/CUDA/AMDGPU cross-host characterization](../benchmarks/results/platform/2026-07-14-wsl-cuda-local-amdgpu.toml)
- [current detector qualification closure catalog](../benchmarks/results/detectors/2026-08-01-detector-qualification-closure.toml)
- [../benchmarks/results/validation_runs/2026-04-10-rtc-devel-amdgpu.toml](../benchmarks/results/validation_runs/2026-04-10-rtc-devel-amdgpu.toml)

Current expectation:

- if a maintained AMDGPU surface regresses numerically against CPU, that is a
  release-blocking defect for the AMDGPU-supported scope

Gate 5 closure validation head `02e5f29` passed all `448` maintained checks on
the local gfx1030 AMDGPU target and `438/438` on the WSL RTX 3050 Ti CUDA
target with Julia 1.12.6 and scalar indexing disabled. Gate 7.4 added
integrated direct-science same-device path-batch ownership; Gate 7.5 extends
that ownership boundary to maintained SH/Pyramid/BioEdge paths and adds the
six-row conventional detector matrix. Focused Gate 7.5 validation completed
separate `212/212` composed WFS-to-detector and `155/155` standalone detector
checks on local AMDGPU and the same counts on WSL CUDA.jl 6.2.1. These checks
cover numerical oracles, exact device/context ownership, residency, lifecycle
rejection, and bounded allocation evidence. The separate clean
[Gate 5 CPU artifact](../benchmarks/results/gate5/2026-07-25-optical-placement.toml)
remains the integrated serial optical-placement characterization. Gate 7.5
correctness evidence is complemented by the Gate 7.6 predeclared paired
two-path Shack-Hartmann
[service-cost catalog](../benchmarks/results/gate7/manifest.toml) on local
CPU, local gfx1030 AMDGPU, and WSL RTX 3050 Ti CUDA. It records exact parity,
residency, first use, allocation, raw HdrHistograms, and separate device-ready,
host-ready, and transfer-only boundaries. Both accelerator-relative p95 gates
pass. This strengthens the declared single-device execution surface without
claiming fixed-arrival HIL latency, mixed placement, multi-GPU execution, or
instrument capacity. CUDA's routine support status remains governed
separately by the support-boundary rule below.

Detector qualification candidate `dd16596` passed 244/244 focused AMDGPU
detector checks. The complete AMDGPU target separately encounters a compiler
segmentation fault in the unrelated Pyramid geometric-slope kernel at this
revision; the supported detector claim is therefore tied to the focused target
and the exact scope in the closure catalog, not presented as a complete-target
pass. The broader defect is tracked in
[issue #200](https://github.com/DarrylGamroth/AdaptiveOpticsSim.jl/issues/200).

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

- CUDA execution. The extension, dedicated test project, focused 250/250
  detector target, complete 1,028/1,028 hardware target, and current manual WSL
  evidence—including the exact detector family scope in the
  [closure catalog](../benchmarks/results/detectors/2026-08-01-detector-qualification-closure.toml), prepared direct
  imaging, deterministic shared frame-response/IPC acquisition, and the
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
- reduced-order behavior outside the maintained deterministic linear modal
  direct-measurement envelope. The supported model validates exact scheduled
  command causality, exposure averaging, schema bindings, and a matched
  reference loop; it does not claim raw pixels, detector physics, tomography,
  external-RTC integration, integrated GPU event-loop execution, full-optical
  equivalence, general control stability, or instrument-scale performance
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

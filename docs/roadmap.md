# AdaptiveOpticsSim Roadmap

## Status
- Phases 1-5 implemented in `AdaptiveOpticsSim.jl`.
- Phase 6 in progress (telemetry + config export done; more I/O helpers pending).
- Phase 7 complete (user guide, API reference, tutorial ports, and committed OOPAO reference bundle in place).
- Core feature parity and numerical fidelity with Python OOPAO are effectively
  in place for the committed deterministic workflows. Remaining work is now
  dominated by robustness, execution quality, and a few non-core surface gaps.
- The most current OOPAO comparison document is
  [`docs/python-julia-differences.md`](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/docs/python-julia-differences.md).
  This roadmap is primarily for implementation status and backlog planning.

## AO Feature Checklist
- [x] P0: Detector modeling (photon/read noise, QE, binning, PSF sampling,
  gain/sensor staging, dark current, saturation, ADC quantization, background
  flux/map handling).
- [x] P0: Multi-source/asterism support and off-axis PSF handling.
- [x] P1: Pyramid WFS implementation.
- [x] P1: BioEdge WFS implementation.
- [x] P1: LGS-specific sensing paths (elongation, sodium profile).
- [x] P0: Multi-layer atmosphere with wind and frozen-flow evolution.
- [x] P0: Spiders and user-defined pupil masks.
- [x] P1: DM mis-registration (shift/rotation/anamorphosis) and custom influence functions.
- [x] P1: Geometric WFS sensing mode (gradient-only) option.

## OOPAO Feature Parity (Module-level)
- [x] `Telescope.py` -> `src/Optics/telescope.jl`
- [x] `Source.py` -> `src/Optics/source.jl`
- [x] `Atmosphere.py` (basic phase screens, multi-layer) -> `src/Atmosphere/kolmogorov.jl`, `src/Atmosphere/multilayer.jl`
- [x] `ShackHartmann.py` -> `src/WFS/shack_hartmann.jl`
- [x] `Pyramid.py` -> `src/WFS/pyramid.jl`
- [x] `BioEdge.py` -> `src/WFS/bioedge.jl`
- [x] `DeformableMirror.py` -> `src/Optics/deformable_mirror.jl`
- [x] `Detector.py` -> `src/Optics/detector.jl`
- [x] `Asterism.py` -> `src/Optics/asterism.jl`
- [x] `MisRegistration.py` -> `src/Optics/misregistration.jl`
- [x] `Zernike.py` -> `src/Optics/zernike.jl`
- [x] `GainSensingCamera.py` -> `src/Calibration/gain_sensing_camera.jl`
- [x] `LiFT.py` -> `src/WFS/lift.jl`
- [x] `NCPA.py` -> `src/Optics/ncpa.jl`
- [x] `OPD_map.py` -> `src/Optics/opd_map.jl`
- [x] `phaseStats.py` -> `src/Atmosphere/phase_stats.jl`
- [x] `SpatialFilter.py` -> `src/Optics/spatial_filter.jl`
- [x] `SPRINT.py` -> `src/Calibration/misregistration_identification.jl`
- [x] `calibration/ao_calibration.py` -> `src/Calibration/ao_calibration.jl`
- [x] `calibration/InteractionMatrix.py` -> `src/Calibration/interaction_matrix.jl`
- [x] `calibration/CalibrationVault.py` -> `src/Calibration/calibration_vault.jl`
- [x] `calibration/compute_KL_modal_basis.py` -> `src/Calibration/modal_basis.jl`
- [x] `calibration/get_modal_basis.py` -> `src/Calibration/modal_basis.jl`
- [x] `calibration/getFittingError.py` -> `src/Calibration/fitting_error.jl`
- [x] `calibration/get_fast_atmosphere.py` -> `src/Calibration/fast_atmosphere.jl`
- [x] `calibration/initialization_AO.py` -> `src/Calibration/initialization.jl`
- [x] `calibration/initialization_AO_PWFS.py` -> `src/Calibration/initialization.jl`
- [x] `calibration/initialization_AO_SHWFS.py` -> `src/Calibration/initialization.jl`
- [x] `closed_loop/run_cl.py` -> `examples/closed_loop/run_cl.jl`
- [x] `closed_loop/run_cl_first_stage.py` -> `examples/closed_loop/run_cl_first_stage.jl`
- [x] `closed_loop/run_cl_from_phase_screens.py` -> `examples/closed_loop/run_cl_from_phase_screens.jl`
- [x] `closed_loop/run_cl_long_push_pull.py` -> `examples/closed_loop/run_cl_long_push_pull.jl`
- [x] `closed_loop/run_cl_sinusoidal_modulation.py` -> `examples/closed_loop/run_cl_sinusoidal_modulation.jl`
- [x] `closed_loop/run_cl_two_stages.py` -> `examples/closed_loop/run_cl_two_stages.jl`
- [x] `closed_loop/run_cl_two_stages_atm_change.py` -> `examples/closed_loop/run_cl_two_stages_atm_change.jl`
- [x] `mis_registration_identification_algorithm/*` -> `src/Calibration/misregistration_identification.jl`
- [ ] `tools/displayTools.py` -> not yet implemented
- [ ] `tools/OopaoGUI.py` -> not yet implemented
- [ ] `tools/interpolateGeometricalTransformation.py` -> not yet implemented
- [ ] `tools/interpolate_influence_functions.py` -> not yet implemented
- [ ] `tools/set_paralleling_setup.py` -> not yet implemented
- [ ] `tools/tools.py` -> not yet implemented

## Current Differences vs OOPAO

These are not all "simplifications". At this point the package is a mix of
matched behavior, intentional Julia design differences, Julia-only extensions,
and a few still-narrower defaults.

Equivalent or improved relative to the practical OOPAO surface:

- Diffractive SH, Pyramid, and BioEdge all use FFT-based propagation and now
  cover the practical OOPAO parity surface, including the main SH parity knobs
  and the current deterministic reference bundle
  (`src/WFS/shack_hartmann.jl`, `src/WFS/pyramid.jl`, `src/WFS/bioedge.jl`).
- LiFT analytical mode uses FFT-based Jacobians with object convolution and
  readout-weighted noise; detector-pixel treatment for binned data matches the
  practical OOPAO surface (`src/WFS/lift.jl`).
- GainSensingCamera is parity-backed for the committed deterministic cases and
  now has a cleaner threaded/cache-aware implementation plus explicit detector
  metadata attachment (`src/Calibration/gain_sensing_camera.jl`).

Intentional Julia differences:

- SPRINT / mis-registration keeps cached sensitivity matrices and optional WFS
  shifts, but deliberately does not bake FITS into core
  (`src/Calibration/misregistration_identification.jl`).
- Detector, tomography, inverse-policy, and backend/GPU abstractions are more
  explicit than in OOPAO. In several of these areas Julia now exceeds OOPAO in
  engineering structure rather than lagging it.

Julia extensions beyond the main OOPAO surface:

- LGS elongation uses Na-profile convolution for Shack-Hartmann, Pyramid, and
  BioEdge. Shack-Hartmann uses per-subaperture Na-profile kernels in the
  OOPAO-compatible path; Pyramid/BioEdge use averaged Na-profile kernels as a
  Julia-side extension. OOPAO does not appear to expose an equivalent
  per-subaperture Na-profile path for Pyramid/BioEdge in the main
  code/tutorial surface (`src/WFS/shack_hartmann.jl`,
  `src/WFS/pyramid.jl`, `src/WFS/bioedge.jl`).

Remaining narrower defaults or simplifications:

- `ft_sh_phase_screen` now uses adaptive sub-harmonic levels and a wider
  stencil by default, but the atmosphere model is still simpler than a fully
  layer/outer-scale-specific construction (`src/Atmosphere/phase_stats.jl`).

This is not a deterministic-parity blocker, but it remains part of the
model-fidelity roadmap and should be revisited if additional tilt/statistics
validation shows the new default is still too coarse.

## Feature Parity and Numerical Fidelity Gate
- [x] The committed OOPAO reference bundle now matches
  `psf_baseline`, `shack_hartmann_diffractive_ramp`,
  `pyramid_diffractive_ramp`, `bioedge_diffractive_ramp`,
  `gain_sensing_camera_optical_gains`, `transfer_function_rejection`,
  `lift_interaction_matrix`, `closed_loop_shack_hartmann_trace`,
  `closed_loop_pyramid_trace`, `closed_loop_bioedge_trace`,
  `gsc_closed_loop_trace`, `gsc_atmosphere_replay_trace_bounded`,
  `gsc_branch_step_modulation_frame`, `gsc_branch_step_optical_gains`, and
  `gsc_branch_step_signal` within their documented tolerances.
- [x] Expand deterministic OOPAO reference coverage beyond geometric Shack-Hartmann to:
  LiFT and compact closed-loop traces for Shack-Hartmann, Pyramid, and BioEdge.
- [x] Match OOPAO PSF export conventions and normalization exactly enough to support
  reproducible array-level regression for image formation.
- [x] Port and validate OOPAO transfer-function workflow
  (`tutorials/AO_transfer_function.py`) with matching outputs.
- [x] Validate the stable and deterministic parts of the atmosphere-driven OOPAO
  GSC closed-loop workflow with compact traces, bounded replay, and first
  branch-step regression coverage.

## Diagnostic Stress Cases
- [ ] Full long-horizon atmosphere-driven OOPAO GSC replay
  (`tutorials/AO_closed_loop_Pyramid_WFS_GSC.py`) remains a diagnostic stress
  case rather than a strict parity gate. The main Pyramid incidence-flux
  mismatch from modulation-point averaging has been corrected, a bounded
  atmosphere-replay parity case is now committed, and the first nonlinear
  branch step is regression-backed at the fixed-input level for Pyramid
  signals, modulation frames, and optical gains. Local long-horizon replay
  diagnostics now show that the control coefficients,
  optical-gain floor behavior, and residual-RMS telemetry track OOPAO closely.
  The remaining drift is concentrated in nonlinear observables (mainly Strehl
  and later slope norms) once the replay enters a huge-OPD regime, so the full
  tutorial-scale replay remains diagnostic-only rather than a strict parity gate.
- [x] Port and validate the pyTomoAO-backed tomography workflow needed for
  OOPAO parity, including model-based wavefront reconstruction and assembled
  DM commands against the KAPA benchmark configuration.
- [x] Julia tomography subsystem now exists under `src/Tomography/` with typed
  parameter objects, geometry helpers, DM fitting, sparse gradient assembly,
  model-based covariance/reconstructor operators (`Gamma`, `Cxx`, `Cox`, `Cnz`,
  `RecStatSA`), the IM-based reconstructor path, and committed pyTomoAO compact
  numerical regression for the model/IM operators and reconstructed wavefronts.
- [x] Extend tomography parity from the current compact pyTomoAO bundle to a
  tutorial-scale pyTomoAO workflow with committed KAPA wavefront and DM-command
  regressions.
- [ ] Audit calibration/output conventions against OOPAO for:
  slope ordering/units outside the now-explicit geometric-SH compare
  convention, PSF sampling conventions, and detector coupling.
- [ ] Add regression-backed fidelity tests for every supported OOPAO tutorial/workflow,
  not just smoke tests.
- [ ] Keep the current simplifications only where they are proven numerically equivalent
  or are explicitly marked as fast-path approximations.

## Candidate Algorithm Upgrades (Similar Results)
- (DONE) Add sub-harmonic augmentation to `ft_sh_phase_screen` for better low-frequency tilt statistics (`src/Atmosphere/phase_stats.jl`) [E:med, R:low].
- (DONE) Implement diffractive WFS paths via pupil→focal propagation with planned FFTs (`src/WFS/shack_hartmann.jl`, `src/WFS/pyramid.jl`, `src/WFS/bioedge.jl`) [E:high, R:med].
- (DONE) Replace LGS slope scaling with focal-plane elongated PSF modeling (`src/WFS/shack_hartmann.jl`, `src/WFS/pyramid.jl`) [E:med, R:med].
- (DONE) Implement LiFT analytic Jacobians via FFT-based formulation (`src/WFS/lift.jl`) [E:med, R:med].
- (DONE) Add KL basis from atmospheric covariance (HHt/PSD) instead of DM-mode covariance (`src/Calibration/modal_basis.jl`) [E:high, R:med].
- (DONE) Implement SPRINT fast path with cached sensitivity matrices (Serialization) (`src/Calibration/misregistration_identification.jl`) [E:med, R:low].
- (DONE) Add Pyramid modulation and optical gain hooks; BioEdge adds optical gain (`src/WFS/pyramid.jl`, `src/WFS/bioedge.jl`) [E:med, R:med].
- Add a fast-path config/trait layer that explicitly opts into simplified models for speed, leaving high-fidelity models as opt-in (`src/Core/parallel.jl`, `src/Core/types.jl`) [E:low, R:low].

## Reference Packages and Candidate Algorithms
- OOMAO (MATLAB): maps to KL basis from HHt/PSD, slope covariance, LQG/predictive control, and tomography upgrades.
- Soapy (Python): maps to diffractive WFS propagation, LGS modeling variants, and full end-to-end loop structure.
- COMPASS (C/Fortran): maps to GPU-oriented diffractive WFS and tomography pipelines, plus predictive control.
- MAOS (C): maps to LGS elongation/cone-effect modeling and multi-DM MCAO/MOAO workflows.
- PASSATA / YAO (IDL): maps to reconstructor variants and multi-object AO workflows.
- HCIPy / POPPY / PROPER (Python): maps to diffractive WFS propagation, PSF/coronagraph modeling.

## Suggested Near-Term Priorities
- [ ] Write and implement the interface/style specification pass described in
  `docs/interface-style-spec-plan.md`, with explicit abstract-family contracts,
  capability-trait audit, and conformance tests.
- [ ] Improve the remaining scientifically narrower defaults, especially
  `ft_sh_phase_screen` and the default NCPA KL basis choice.
- [~] Define explicit execution/fidelity profiles so the package can support a
  scientifically stronger default profile and a faster HIL-oriented profile
  without conflating the two. `ScientificProfile` / `FastProfile` now exist in
  core, can be bundled with per-subsystem overrides via `ProfileBundle`, and
  currently drive sub-harmonic phase-screen defaults plus NCPA KL basis
  selection; more surfaces still need to be wired through them.
- [ ] Finish the remaining calibration/output convention audit against OOPAO
  telemetry exports. Compact closed-loop trace containers are now explicit in
  core (`ClosedLoopTrace`, `GSCClosedLoopTrace`,
  `GSCAtmosphereReplayTrace`); the remaining audit is mainly about slope
  ordering/units and any remaining PSF sampling mismatches. Detector-coupled
  runtime exports are now explicit through `DetectorExportMetadata`,
  `detector_export_metadata`, `simulation_wfs_metadata`, and
  `simulation_science_metadata`.
- [ ] Decide whether the long-horizon atmosphere-driven Pyramid/GSC replay stays
  purely diagnostic or receives more robustness work.
- [x] Extend HIL-focused support toward multi-WFS / multi-DM aggregation for
  MOAO, MCAO, and woofer/tweeter RTC scenarios. The staged implementation now
  lives in `docs/stacked-multi-source-multi-wfs-plan.md` and includes stacked
  SH asterism, stacked Pyramid/BioEdge asterism, grouped
  `CompositeSimulationInterface` execution, and the maintained profiler
  `scripts/profile_multi_source_multi_wfs_runtime.jl`.
- [x] Add a maintained AO188/3k operational simulation benchmark configuration.
  The maintained AO188/3k simulation now lives as example support in
  `examples/support/subaru_ao188_simulation.jl`, with
  `scripts/ao188_3k_hil_audit.jl` as its maintained audit entry point. It
  provides a concrete 64x64 / 3228-active / 188-mode HIL-oriented simulation
  with a split high/low-order SH proxy, explicit control-path latency
  staging, tuned detector/noise defaults, and a circular active-actuator
  support map.
- [ ] Add specialized HIL-relevant sensors on demand, starting with
  `ZernikeWFS` and then `CurvatureSensor` if needed.
- [ ] Build the next detector/readout-physics layer following
  `docs/detector-readout-physics-plan.md`, starting with APD-like counting
  behavior for curvature sensing, then frame-detector MTF/pixel-response work,
  while keeping room for future SPAD/SPAD-array support.
- [ ] Continue GPU/HIL performance work only where profiling shows real value,
  especially in builder-heavy tomography paths.
- [~] Separate compact regression benchmarks from representative benchmark
  ladders across the major runtime and builder families, following
  `docs/benchmark-matrix-plan.md`, so CPU/GPU and algorithm-appropriateness
  decisions are based on realistic scale points rather than only tiny guard
  cases. The first runtime profilers now expose explicit `compact` / `medium` /
  `representative` modes for `ZernikeWFS`, AO188-style pixel-output runtime,
  and stacked multi-source / multi-WFS runtime; broader representative ladders
  are still pending.
- [~] Extract generic control-simulation primitives from the Subaru AO188
  simulation into core without moving the Subaru-specific assembly itself.
  The staged plan now lives in `docs/control-simulation-architecture-plan.md`.
  `MappedReconstructor`, `VectorDelayLine`, generic execution policies, and the
  first control-simulation hooks/capability traits are now in core, and both
  the AO188 example and `ClosedLoopRuntime` now conform to the shared
  control-simulation interface; broader adoption beyond those runtime surfaces
  is still pending.

## Current Focus

The earlier numerical-stability backlog has largely been addressed:

- explicit inverse policies are implemented for modal/control reconstruction
- GSC weak-mode flooring is implemented
- tomography uses explicit Hermitian solves and structured noise-model policies
- LiFT has QR + SVD fallback and fixed/adaptive LM damping modes
- covariance/Bessel approximations have direct regression coverage

The current focus is now:

1. scientific/default-model completeness
2. calibration/output convention audits against OOPAO where still needed
3. explicit fidelity/execution profiles for accuracy vs throughput
4. HIL/RTC-facing execution quality
5. multi-WFS / multi-DM support
6. specialized HIL-relevant sensors such as `ZernikeWFS`
7. targeted GPU builder/runtime performance work where profiling justifies it,
   with the current priority on builder-heavy paths rather than per-frame
   pixel-output runtime loops

## Next Tasks
1. Finish the remaining calibration/output convention audit against OOPAO
   telemetry exports, now that compact closed-loop trace types are explicit in
   core.
2. Extend the explicit fidelity/execution profile surface beyond the current
   phase-screen and NCPA defaults so more fast/scientific choices are visible
   at the API level. The profile surface is now granular enough to support
   per-subsystem overrides; the next work is to wire more subsystems through
   it.
3. Revisit `ft_sh_phase_screen` again if tilt/statistics validation shows the
   adaptive sub-harmonic default is still too coarse for the scientific
   profile.
4. Decide whether the long-horizon Pyramid/GSC replay remains purely diagnostic.
5. Add multi-WFS / multi-DM aggregation for MOAO, MCAO, and woofer/tweeter RTC scenarios.
6. Implement `CurvatureSensor` if the RTC/HIL use case requires it after `ZernikeWFS`.
7. Continue extending GPU-native builder coverage where HIL workflows demand it.
8. Profile and optimize tomography builder hotspots only where crossover data shows the GPU path is worthwhile.

## HIL-Oriented Near-Term Work
- [x] Add explicit builder-backend selection for modal/calibration reconstructor generation.
- [x] Add an explicit RTC boundary surface around `ClosedLoopRuntime`.
- [x] Add deterministic runtime latency/jitter measurement helpers and benchmarks.
- [~] Extend builder-backend coverage to tomography and other heavy calibration paths.
  Interaction-matrix tomography reconstruction, model-based tomography
  reconstructor outputs, and tomography command assembly now preserve the
  selected build backend instead of forcing host `Matrix`/`Vector` outputs.
  Tomography covariance fill, fit-source averaging, and noise-covariance
  assembly now also dispatch through the selected build backend, including
  accelerator kernels for the matrix-fill stages. The remaining gap is that
  guide-star coordinate generation and some sparse/operator assembly steps are
  still CPU-originating before backend materialization, so a fully
  accelerator-native tomography builder is still future work.
- [ ] Add specialized HIL-relevant sensors on demand, starting with a
  `ZernikeWFS` implementation and then `CurvatureSensor` if required.
  See [`docs/zernike-wfs-plan.md`](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/docs/zernike-wfs-plan.md).
- [x] Add maintained CUDA validation entry points for runtime and builder paths.
  `scripts/gpu_smoke_cuda.jl` and `scripts/gpu_builder_cuda.jl` are now the
  standard CUDA validation pair, and both have been exercised on `spiders`.
- [x] Bring up an initial AMDGPU validation path.
  `scripts/gpu_smoke_amdgpu.jl` and `scripts/gpu_builder_amdgpu.jl` now pass on
  AMD hardware for the maintained smoke surface. Runtime paths are native on
  `ROCArray`; modal/calibration inverse operators now use native rocSOLVER SVD,
  tomography Hermitian right-division uses native rocSOLVER Cholesky with
  native LU fallback, and LiFT normal-equation plus fallback solves now use
  native rocSOLVER paths too.
- [x] Add HIL-oriented AMDGPU smoke and sync-audit entry points.
  `scripts/gpu_hil_amdgpu.jl` and `scripts/gpu_sync_audit_amdgpu.jl` now pass
  on AMD hardware for the maintained smoke and timing surface.
- [x] Add maintained AMDGPU crossover and tomography profiling entry points.
  `scripts/gpu_crossover_amdgpu.jl`,
  `scripts/gpu_profile_model_tomography_amdgpu.jl`, and
  `scripts/gpu_profile_model_tomography_phases_amdgpu.jl` now provide warmed
  CPU/GPU crossover data and tomography hotspot profiling on AMD hardware.
- [x] Add HIL-oriented CUDA smoke and sync-audit entry points.
  `scripts/gpu_hil_cuda.jl` runs the combined runtime + builder smoke surface,
  and `scripts/gpu_sync_audit_cuda.jl` reports the RTC-facing runtime/build
  timing surface for CUDA hosts.
- [~] Treat per-frame GPU runtime for pixel-output Shack-Hartmann HIL as a
  redesign candidate rather than an optimization target.
  The maintained AO188/3k simulation audit currently shows CPU runtime
  outperforming both AMDGPU and CUDA for the warmed frame loop, with the
  dominant cost in the high-order diffractive Shack-Hartmann detector path.
  The current accelerator implementation still loops over subapertures and
  captures each spot individually, so further work here should only proceed if
  we are willing to redesign that path around a batched or fused pixel
  pipeline. Until then, GPU maintenance effort is better spent on builder-heavy
  calibration/tomography workflows than on the single-frame RTC loop. See
  [`docs/pixel-output-gpu-runtime-redesign.md`](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/docs/pixel-output-gpu-runtime-redesign.md).
- [ ] Add multi-WFS / multi-DM aggregation for MOAO, MCAO, and woofer/tweeter RTC scenarios.
- [ ] Add specialized HIL-relevant sensors:
  `ZernikeWFS`, `CurvatureSensor`, and distributed/multi-WFS aggregation.

## Phase 0: Setup
- Create package skeleton and CI with Julia versions.
- Define core interfaces, traits, and basic error types.
- Add logging defaults and configuration scaffolding.

## Phase 1: Optics vertical slice
- Implement Telescope + Source + PSF + Zernike.
- Add `Workspace` with deterministic RNG wiring.
- Add basic unit tests for PSF normalization and Zernike orthogonality.

## Phase 2: Atmosphere + propagation
- Implement phase screen model and `advance!`.
- Integrate with Telescope propagation.
- Add deterministic mode tests and reference outputs.

## Phase 3: WFS + DM
- Implement Shack-Hartmann WFS (sampling, slopes, noise).
- Implement DeformableMirror with influence functions.
- Add minimal closed-loop demo (no controller yet).

## Phase 4: Calibration and control
- Interaction matrix and reconstructor utilities.
- Optional ModelingToolkit-based controller module.
- Add tutorial mappings as runnable examples.

## Phase 5: Performance + GPU
- Allocation-free hot paths and benchmarks.
- Optional GPU backend via `CuArray`.
- CPU/GPU parity tests (tolerances documented).
- GPU checklist: no scalar indexing, minimize transfers, preallocated workspaces.
- Benchmark suite for CPU/GPU with standardized scenarios (PSF, WFS, reconstructor).

## Phase 6: IO and telemetry
- [~] Optional I/O helpers (CSV telemetry via extension; FITS/HDF5 pending).
- [x] Tables.jl telemetry output for metrics.
- [x] Config and provenance export (TOML + JSON3 extension).

## Phase 7: Documentation and validation
- [x] Port most tutorials as Julia scripts/notebooks (prioritize core AO workflows).
- [x] Add deterministic regression suite against OOPAO outputs.
- [x] Maintain a small set of OOPAO reference datasets for cross-validation.
- [x] Make OOPAO reference-bundle regeneration reproducible from a pinned upstream
  repo/ref instead of assuming a manually prepared local checkout.
- [x] Publish user guide and API reference.

Current validation scope:
- `examples/tutorials/` now covers the core OOPAO workflows used most often in practice.
- `test/reference_data/` now commits deterministic OOPAO PSF, geometric Shack-Hartmann,
  diffractive Shack-Hartmann, Pyramid, BioEdge, GSC optical-gain, and transfer-function cases.
- `test/reference_data/` also commits deterministic pyTomoAO tomography cases for
  `Gamma`, `Cxx`, `Cox`, `Cnz`, both reconstructors, both reconstructed wavefront maps,
  and KAPA benchmark DM-command assembly.
- The reference harness applies a documented convention adapter only where OOPAO and Julia
  intentionally expose different public conventions in geometric Shack-Hartmann mode.
- Local expanded bundle audits now cover PSF, diffractive SH/Pyramid/BioEdge, GSC optical gains,
  transfer functions, and pyTomoAO tomography operators.
- PSF, diffractive Shack-Hartmann, Pyramid, BioEdge, GSC optical gains, and transfer functions
  are now regression-backed against the committed OOPAO bundle.
- Compact tomography operators and reconstructions, plus KAPA benchmark wavefront and
  DM-command reconstruction, are now regression-backed against the committed pyTomoAO bundle.
- pyTomoAO-specific row-major ordering is now handled only in the reference harness;
  the `AdaptiveOpticsSim.jl` tomography implementation remains native Julia column-major internally.

## Phase 8: Full OOPAO Feature Parity and Numerical Fidelity
- [ ] Every supported OOPAO workflow has deterministic Julia vs OOPAO regression coverage.
- [ ] PSF, diffractive WFS, LiFT, GSC closed-loop behavior, and closed-loop traces match OOPAO within documented tolerances.
- [x] Tomography feature parity covers the pyTomoAO-backed OOPAO workflow with validated model-based or IM-based reconstructors.
- [x] Transfer-function workflow is ported and validated.
- [ ] Any unsupported OOPAO workflow is explicitly documented as out-of-scope rather than implied complete.
- [ ] Only after this phase is complete should non-parity feature work resume.

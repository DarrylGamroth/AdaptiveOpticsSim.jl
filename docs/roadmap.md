# Julia Port Roadmap (Draft)

## Status
- Phases 1-5 implemented in `AdaptiveOptics.jl`.
- Phase 6 in progress (telemetry + config export done; more I/O helpers pending).
- Phase 7 complete (user guide, API reference, tutorial ports, and committed OOPAO reference bundle in place).
- Core feature parity and numerical fidelity with Python OOPAO are effectively
  in place for the committed deterministic workflows. Remaining work is now
  dominated by robustness, execution quality, and a few non-core surface gaps.

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

## Simplifications vs OOPAO (Current State)
- WFS diffractive models use FFT propagation; Shack-Hartmann supports pixel-scale
  and subaperture sampling controls (binning + crop/pad), plus OOPAO-style
  spot-centering phasors and centroid/convolution thresholds. Pyramid now models PSF
  centering, old/new mask variants, rooftop tuning, rotation, and OOPAO-style
  modulation-path construction. BioEdge now uses a diffractive mask stack with
  modulation and grey-width/grey-length mask variants instead of the
  phase-gradient surrogate
  (`src/WFS/shack_hartmann.jl`, `src/WFS/pyramid.jl`, `src/WFS/bioedge.jl`).
- LGS elongation uses Na-profile convolution for Shack-Hartmann, Pyramid, and BioEdge.
  Shack-Hartmann uses per-subaperture Na-profile kernels in the OOPAO-compatible path;
  Pyramid/BioEdge currently use averaged Na-profile kernels as a Julia-side extension.
  OOPAO does not appear to expose an equivalent per-subaperture Na-profile path for
  Pyramid/BioEdge in the main code/tutorial surface, so this is no longer tracked as a
  blocking OOPAO parity gap (`src/WFS/shack_hartmann.jl`, `src/WFS/pyramid.jl`,
  `src/WFS/bioedge.jl`).
- `ft_sh_phase_screen` uses a simple 3x3 sub-harmonic grid (no layer-specific outer-scale tuning)
  (`src/Atmosphere/phase_stats.jl`).
- LiFT analytical mode uses FFT-based Jacobians with object convolution and readout-weighted
  noise; binned detector frames are supported by treating `img_resolution` in detector
  pixels, matching the practical OOPAO surface (`src/WFS/lift.jl`).
- GainSensingCamera uses threaded FFT batching and can attach detector metadata
  snapshots without coupling detector physics into the optical-gain math
  (`src/Calibration/gain_sensing_camera.jl`).
- NCPA KL basis defaults to DM-mode covariance (`MᵀM`); HHt/PSD option is available but not default
  (`src/Calibration/modal_basis.jl`, `src/Optics/ncpa.jl`).
- SPRINT/mis-registration supports cached sensitivity matrices (Serialization) and
  optional WFS shifts, but no FITS I/O (`src/Calibration/misregistration_identification.jl`).

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
  slope ordering/units outside geometric SH, PSF sampling conventions, detector coupling,
  and closed-loop telemetry traces.
- [ ] Add regression-backed fidelity tests for every supported OOPAO tutorial/workflow,
  not just smoke tests.
- [ ] Keep the current simplifications only where they are proven numerically equivalent
  or are explicitly marked as fast-path approximations.

## Candidate Algorithm Upgrades (Similar Results)
- (DONE) Add sub-harmonic augmentation to `ft_sh_phase_screen` for better low-frequency tilt statistics (`src/Atmosphere/phase_stats.jl`) [E:med, R:low].
- (DONE, simplified) Implement diffractive WFS paths via pupil→focal propagation with planned FFTs (`src/WFS/shack_hartmann.jl`, `src/WFS/pyramid.jl`, `src/WFS/bioedge.jl`) [E:high, R:med].
- (DONE, simplified) Replace LGS slope scaling with focal-plane elongated PSF modeling (`src/WFS/shack_hartmann.jl`, `src/WFS/pyramid.jl`) [E:med, R:med].
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
- [x] Expand the OOPAO reference bundle to cover at least one compact closed-loop trace per major WFS
  on top of the existing PSF, diffractive SH, Pyramid, BioEdge, LiFT, GSC, and transfer-function cases.
- [ ] Port the remaining atmosphere-driven GSC closed-loop workflow before adding new non-parity features.
- [ ] Tighten the remaining large-aberration nonlinear observable drift in the
  atmosphere-driven Pyramid/GSC replay before promoting the long-horizon case to
  a strict regression gate.
- [ ] Turn every feature-parity claim into a deterministic regression test against OOPAO outputs.

## Next 10 Tasks
1. Audit remaining calibration/output conventions against OOPAO telemetry exports.
2. Decide whether the compact closed-loop traces should be expanded to full tutorial traces with atmosphere replay.
3. Tighten the remaining large-aberration Strehl/slope drift in the long-horizon GSC replay.
4. Validate the remaining atmosphere-driven GSC closed-loop telemetry against OOPAO outputs.
5. Validate LiFT iterative reconstruction outputs, not just the analytic interaction matrix.
6. Decide whether any future tomography expansion stays in core or moves behind an extension/package split later.
7. Only then resume GPU-specific expansion.

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
  the `AdaptiveOptics.jl` tomography implementation remains native Julia column-major internally.

## Phase 8: Full OOPAO Feature Parity and Numerical Fidelity
- [ ] Every supported OOPAO workflow has deterministic Julia vs OOPAO regression coverage.
- [ ] PSF, diffractive WFS, LiFT, GSC closed-loop behavior, and closed-loop traces match OOPAO within documented tolerances.
- [x] Tomography feature parity covers the pyTomoAO-backed OOPAO workflow with validated model-based or IM-based reconstructors.
- [x] Transfer-function workflow is ported and validated.
- [ ] Any unsupported OOPAO workflow is explicitly documented as out-of-scope rather than implied complete.
- [ ] Only after this phase is complete should non-parity feature work resume.

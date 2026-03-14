# Julia Port Roadmap (Draft)

## Status
- Phases 1-5 implemented in `AdaptiveOptics.jl`.
- Phase 6 in progress (telemetry + config export done; more I/O helpers pending).
- Phase 7 in progress (reference harness scaffolded; datasets + tutorial ports pending).

## AO Feature Checklist
- [x] P0: Detector modeling (photon/read noise, QE, binning, PSF sampling).
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
  and subaperture sampling controls (binning + crop/pad). Pyramid now models PSF
  centering and pupil separation in the mask and slope extraction. BioEdge now uses
  a diffractive mask stack instead of the phase-gradient surrogate, but the mask
  is still a binary split (no grey-width/length or rooftop tuning yet)
  (`src/WFS/shack_hartmann.jl`, `src/WFS/pyramid.jl`, `src/WFS/bioedge.jl`).
- LGS elongation uses Na-profile convolution for Shack-Hartmann, Pyramid, and BioEdge; the
  Pyramid/BioEdge path currently averages Na-profile kernels across subapertures rather
  than modeling per-subap kernels (`src/WFS/shack_hartmann.jl`, `src/WFS/pyramid.jl`,
  `src/WFS/bioedge.jl`).
- `ft_sh_phase_screen` uses a simple 3x3 sub-harmonic grid (no layer-specific outer-scale tuning)
  (`src/Atmosphere/phase_stats.jl`).
- LiFT analytical mode uses FFT-based Jacobians with object convolution and readout-weighted
  noise, but does not model detector binning explicitly (`src/WFS/lift.jl`).
- GainSensingCamera uses threaded FFT batching but still ignores detector metadata
  coupling (`src/Calibration/gain_sensing_camera.jl`).
- NCPA KL basis defaults to DM-mode covariance (`MᵀM`); HHt/PSD option is available but not default
  (`src/Calibration/modal_basis.jl`, `src/Optics/ncpa.jl`).
- SPRINT/mis-registration supports cached sensitivity matrices (Serialization) and
  optional WFS shifts, but no FITS I/O (`src/Calibration/misregistration_identification.jl`).

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
- [x] Implement diffractive WFS propagation with planned FFTs (largest fidelity jump).
- [x] Add sub-harmonic augmentation to phase screens (improves low-frequency statistics).
- [x] Replace LGS slope scaling with elongated PSF modeling (better LGS realism).
- [x] Add KL basis from HHt/PSD (better modal basis fidelity).

## Next 10 Tasks
1. Build the OOPAO reference regression harness and manifest format.
2. Generate the first deterministic OOPAO reference bundle for PSF + Shack-Hartmann + Pyramid.
3. Add regression tests that compare Julia outputs to OOPAO within per-case tolerances.
4. Expand the benchmark suite to cover PSF, SH, Pyramid, BioEdge, LiFT, reconstructor, and one closed-loop step.
5. Track allocations in benchmark outputs so hot-path regressions are visible.
6. Bring up a first `CuArray` path for PSF and one diffractive WFS with `CUDA.allowscalar(false)` checks.
7. Add CPU/GPU parity tests and document expected tolerances.
8. Improve BioEdge mask fidelity with grey-width/length and rooftop options.
9. Replace averaged Pyramid/BioEdge Na-profile kernels with per-subaperture kernels where fidelity matters.
10. Port the highest-value tutorials as deterministic examples tied to regression coverage.

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
- Port most tutorials as Julia scripts/notebooks (prioritize core AO workflows).
- Add deterministic regression suite against OOPAO outputs.
- Maintain a small set of OOPAO reference datasets for cross-validation.
- Publish user guide and API reference.

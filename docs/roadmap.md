# Julia Port Roadmap (Draft)

## Status
- Phases 1-5 implemented in `AdaptiveOptics.jl`.
- Phases 6-7 pending (I/O helpers, reference dataset validation, docs/tutorial ports).

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
- [ ] `GainSensingCamera.py` -> not yet implemented
- [ ] `LiFT.py` -> not yet implemented
- [ ] `NCPA.py` -> not yet implemented
- [ ] `OPD_map.py` -> not yet implemented
- [ ] `phaseStats.py` -> not yet implemented
- [ ] `SpatialFilter.py` -> not yet implemented
- [ ] `SPRINT.py` -> not yet implemented
- [ ] `calibration/ao_calibration.py` -> partial (interaction matrix + reconstructor)
- [ ] `calibration/InteractionMatrix.py` -> partial (interaction matrix only)
- [ ] `calibration/CalibrationVault.py` -> not yet implemented
- [ ] `calibration/compute_KL_modal_basis.py` -> not yet implemented
- [ ] `calibration/get_modal_basis.py` -> not yet implemented
- [ ] `calibration/getFittingError.py` -> not yet implemented
- [ ] `calibration/get_fast_atmosphere.py` -> not yet implemented
- [ ] `calibration/initialization_AO.py` -> not yet implemented
- [ ] `calibration/initialization_AO_PWFS.py` -> not yet implemented
- [ ] `calibration/initialization_AO_SHWFS.py` -> not yet implemented
- [ ] `closed_loop/run_cl.py` -> not yet implemented
- [ ] `closed_loop/run_cl_first_stage.py` -> not yet implemented
- [ ] `closed_loop/run_cl_from_phase_screens.py` -> not yet implemented
- [ ] `closed_loop/run_cl_long_push_pull.py` -> not yet implemented
- [ ] `closed_loop/run_cl_sinusoidal_modulation.py` -> not yet implemented
- [ ] `closed_loop/run_cl_two_stages.py` -> not yet implemented
- [ ] `closed_loop/run_cl_two_stages_atm_change.py` -> not yet implemented
- [ ] `mis_registration_identification_algorithm/*` -> not yet implemented
- [ ] `tools/displayTools.py` -> not yet implemented
- [ ] `tools/OopaoGUI.py` -> not yet implemented
- [ ] `tools/interpolateGeometricalTransformation.py` -> not yet implemented
- [ ] `tools/interpolate_influence_functions.py` -> not yet implemented
- [ ] `tools/set_paralleling_setup.py` -> not yet implemented
- [ ] `tools/tools.py` -> not yet implemented

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
- Optional I/O helpers (FITS/HDF5/CSV) in an extension module.
- Tables.jl telemetry output for metrics.
- Config and provenance export (TOML/JSON).

## Phase 7: Documentation and validation
- Port most tutorials as Julia scripts/notebooks (prioritize core AO workflows).
- Add deterministic regression suite against OOPAO outputs.
- Maintain a small set of OOPAO reference datasets for cross-validation.
- Publish user guide and API reference.

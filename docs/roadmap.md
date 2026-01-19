# Julia Port Roadmap (Draft)

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

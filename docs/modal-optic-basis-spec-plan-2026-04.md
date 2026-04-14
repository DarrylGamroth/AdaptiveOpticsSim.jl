# Modal Optic Basis Spec Plan 2026-04

Status: implemented

## Purpose

This plan makes the `ModalControllableOptic` API more explicit by introducing a
small basis-spec layer for modal controllable optics. The goal is to make the
commanded modal surface obvious in user code without trying to unify every
other basis concept in the package.

## Actions

### MBS-1 Add a focused basis-spec interface

- introduce `AbstractModalOpticBasis`
- add first concrete basis specs:
  - `FunctionModalBasis`
  - `MatrixModalBasis`
  - `ZernikeOpticBasis`
  - `CartesianTiltBasis`
  - `QuadraticFocusBasis`

Status:

- implemented

### MBS-2 Route `ModalControllableOptic` through basis specs

- add `ModalControllableOptic(tel, basis::AbstractModalOpticBasis; ...)`
- make the tuple, matrix, and `zernike_modes` constructors route through the
  basis-spec path
- keep the older constructor shapes as maintained compatibility entry points

Status:

- implemented

### MBS-3 Keep convenience wrappers stable

- keep `TipTiltMirror(...)` as a convenience constructor over
  `CartesianTiltBasis`
- keep `FocusStage(...)` as a convenience constructor over
  `QuadraticFocusBasis`

Status:

- implemented

### MBS-4 Update stable docs and examples

- document the preferred basis-spec API in:
  - user guide
  - model cookbook
  - API reference
- keep the docs clear that the constructor-shape overloads still exist, but are
  no longer the preferred surface

Status:

- implemented

### MBS-5 Validate on maintained CPU and GPU surfaces

- keep the focused runtime/control test surface green
- regenerate the maintained multi-optic CPU artifact
- rerun maintained AMDGPU and CUDA runtime-equivalence scripts

Status:

- implemented

## Verification

- `Pkg.test(test_args=["control_and_runtime"])` passes
- `scripts/generate_multi_optic_runtime_artifact.jl` passes
- `scripts/gpu_runtime_equivalence_amdgpu.jl` passes
- `scripts/gpu_runtime_equivalence_cuda.jl` passes on the CUDA validation host

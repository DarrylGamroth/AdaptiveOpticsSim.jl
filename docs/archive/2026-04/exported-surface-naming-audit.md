# Exported Surface Naming Audit

This note records the naming audit of the exported public surface in
`AdaptiveOpticsSim.jl`.

The goal is to identify:

- names that are already in a good long-term shape
- names that should stay as-is even if they are not perfectly symmetrical
- names that are plausible later cleanup targets
- places where constructor-style APIs would be worse than the current function
  form

## Summary

The exported surface is mostly in a good Julia/SciML shape after the breaking
terminology cleanup.

The strongest parts are:

- noun types with verb algorithms
- mutating hot paths marked with `!`
- clearer runtime naming after the `SimulationInterface` cleanup
- explicit trait/hook names instead of hidden concrete-type policy in the
  runtime path

The main recommendation from this audit is:

- keep the current naming stable
- avoid compatibility aliases for removed names
- prefer physical-object and computational-role names for new APIs

## Safe Cleanup Performed

- removed a duplicate `compute_psf!` export from `src/AdaptiveOpticsSim.jl`
- replaced `CalibrationVault` with `ControlMatrix`
- replaced runtime "product" terminology with runtime "output" terminology
- replaced runtime "platform" API names with control-loop orchestration names
- replaced `ShackHartmann` with `ShackHartmannWFS` for WFS naming symmetry

## Names To Keep

These names are already strong enough to keep:

- `SimulationInterface`, `CompositeSimulationInterface`, `SimulationReadout`
- `prepare!`, `prepare_runtime_wfs!`, `step!`
- `measure!`, `reconstruct!`, `apply!`, `update!`
- `ModalReconstructor`, `MappedReconstructor`
- `InteractionMatrix`, `ControlMatrix`, `AOCalibration`
- `ShackHartmannWFS`, `PyramidWFS`, `BioEdgeWFS`
- `DeformableMirror`, `Detector`, `OPDMap`, `SpatialFilter`

Reason:

- they are descriptive enough
- they align with the actual algorithm or data role
- changing them now would create churn without much clarity gain

## Names To Tolerate

These are not perfect stylistically, but should stay for now because the
clarity gain from renaming is not large enough to justify churn:

- `SPRINT`
- `NCPA`, `KLBasis`, `M2CBasis`
- `LGSSource`, `LGSAsterismParams`, `LGSWFSParams`
- `dm_basis`, `basis_from_m2c`
- `fast_atmosphere`

Reason:

- these are established AO abbreviations or compact helper names
- the surrounding docs now explain the mathematical meaning
- clearer alternatives would mostly be aliases, not real improvements

## Later Candidate Cleanup Targets

These are the main exported names that are plausible later cleanup targets if a
future public-surface pass is worth the churn.

### 1. `initialize_ao_shack_hartmann`

Why it stands out:

- the function is a scenario builder rather than a core mathematical primitive

Conservative recommendation:

- do not rename yet
- longer term, these initialization helpers may belong in examples/tutorial
  support rather than the core export set

### 2. `initialize_ao_pyramid`

Why it stands out:

- not wrong, but it is another scenario-builder export rather than a core AO
  operator or simulation primitive

Conservative recommendation:

- keep for now
- reconsider together with `initialize_ao_shack_hartmann`, not separately

### 3. `print_optical_path`

Why it stands out:

- `optical_path` already carries the main semantic value
- the printing helper is small and possibly more convenience than core API

Conservative recommendation:

- keep exported for now
- if the public surface needs trimming later, this is a reasonable candidate to
  de-emphasize before touching the core simulation names

## Constructor Guidance

### Keep workflow builders as functions

The following should remain plain functions, not become the primary
constructor-style API:

- `interaction_matrix(...)`
- `modal_basis(...)`
- `ao_calibration(...)`
- `compute_meta_sensitivity_matrix(...)`

Reason:

- these are workflow/build operations, not simple data construction
- they can be expensive and may measure or assemble operators
- they have multiple physically distinct overload families
- using noun-type constructors for those workflows would blur the difference
  between:
  - "wrap/store this already-computed object"
  - "run a calibration/build algorithm and return the result"

This is good current separation:

- result/container types:
  - `InteractionMatrix`
  - `ModalBasis`
  - `ControlMatrix`
  - `AOCalibration`
- algorithmic builders:
  - `interaction_matrix`
  - `modal_basis`
  - `ao_calibration`

### Keep constructors for reusable objects and result wrappers

These are already appropriate constructor-style APIs:

- `ControlMatrix(D; ...)`
- `ModalReconstructor(imat; ...)`
- `MappedReconstructor(command_basis, imat; ...)`
- `LiFT(tel, src, basis, det; ...)`
- `GainSensingCamera(mask, basis; ...)`
- `SimulationInterface(runtime)`

Reason:

- they create reusable objects with stable internal state
- the noun type is the thing the user keeps and works with afterward

### If constructors are added later, they should be thin

If the package ever adds constructor aliases for workflow builders, they should
be thin convenience wrappers, not replacements for the canonical function form.

For example, this would be acceptable later:

- `InteractionMatrix(dm, wfs, tel; ...) = interaction_matrix(dm, wfs, tel; ...)`

But the recommendation from this audit is:

- do not add those aliases now
- the current function form is clearer and more rigorous

## Current Recommendation

Do not do another broad exported-surface rename pass now.

Instead:

1. keep the current API stable
2. keep using function builders for calibration workflows
3. revisit only the scenario-builder exports later if the core surface needs
   trimming
4. prioritize feature work and scientific/runtime correctness over cosmetic
   public renaming

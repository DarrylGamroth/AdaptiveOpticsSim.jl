# Composite Plant OOPAO Review 2026-04

Status: active

## Purpose

This note records the targeted review of the narrow external OOPAO composite
plant surface added for:

- `tiptilt + dm` on diffractive `ShackHartmann`
- `tiptilt + dm` on diffractive `Pyramid`
- `tiptilt + dm` on diffractive `BioEdge`

The review answers a specific question:

- are the remaining cross-simulator residuals more consistent with algorithm /
  convention differences, or with numerical-stability problems in
  AdaptiveOpticsSim.jl?

## Review scope

The review checked four things:

1. OPD-level agreement before the WFS
2. whether the residual is dominated by tip/tilt or by the DM approximation
3. whether the WFS mismatch scales smoothly with command amplitude
4. whether the observed behavior is consistent with backend-independent
   deterministic behavior already established elsewhere

## OPD-level findings

The OPD comparison used:

- telescope: `resolution=24`, `diameter=8.0`, no central obstruction
- tip amplitude: `5e-9`
- DM command:
  - `[0.0, 1.5e-9, -0.5e-9, 0.0, 1.0e-9, -2.0e-9, 0.75e-9, -0.5e-9, -0.75e-9, 1.25e-9, -1.5e-9, 0.5e-9, 0.0, -0.5e-9, 1.0e-9, 0.0]`

Measured within the pupil:

### Tip only

- `maxabs = 8.27e-25`
- `l2rel = 1.48e-16`
- `corr = 1.0`

Interpretation:

- the `TipTiltMirror` / `ModalControllableOptic(CartesianTiltBasis())` surface
  matches the OOPAO Cartesian tip injection exactly on this grid
- the residual is not coming from the modal basis itself

### DM only

- `maxabs = 1.97e-10`
- `l2rel = 1.3815e-1`
- `corr = 0.9910`

Interpretation:

- the remaining OPD discrepancy is dominated by the DM mapping
- this is expected because the Julia `DeformableMirror` influence-width model is
  only approximately mapped onto the OOPAO `DeformableMirror` mechanical
  coupling parameterization

### Composite `tiptilt + dm`

- `maxabs = 1.97e-10`
- `l2rel = 2.9080e-2`
- `corr = 0.9996`

Interpretation:

- the composite OPD residual is still almost entirely the DM residual
- once tip/tilt is added, the composite OPD remains highly correlated between
  the two simulators
- this is a model-approximation gap, not a sign flip, instability, or broken
  modal composition path

## WFS amplitude sweep

The composite command was scaled by `0.25`, `0.5`, `1.0`, and `2.0`, with both
the tip amplitude and the DM coefficients scaled together.

Measured WFS residuals:

### Shack-Hartmann

- `scale=0.25`: `maxabs = 7.89e-5`, `l2rel = 3.6996e-2`
- `scale=0.5`: `maxabs = 1.58e-4`, `l2rel = 3.6996e-2`
- `scale=1.0`: `maxabs = 3.15e-4`, `l2rel = 3.6996e-2`
- `scale=2.0`: `maxabs = 6.31e-4`, `l2rel = 3.6996e-2`

### Pyramid

- `scale=0.25`: `maxabs = 7.38e-5`, `l2rel = 2.6651e-2`
- `scale=0.5`: `maxabs = 1.48e-4`, `l2rel = 2.6666e-2`
- `scale=1.0`: `maxabs = 2.97e-4`, `l2rel = 2.6696e-2`
- `scale=2.0`: `maxabs = 5.99e-4`, `l2rel = 2.6762e-2`

### BioEdge

- `scale=0.25`: `maxabs = 2.15e-4`, `l2rel = 3.4651e-2`
- `scale=0.5`: `maxabs = 4.29e-4`, `l2rel = 3.4651e-2`
- `scale=1.0`: `maxabs = 8.59e-4`, `l2rel = 3.4651e-2`
- `scale=2.0`: `maxabs = 1.72e-3`, `l2rel = 3.4653e-2`

Interpretation:

- `maxabs` grows approximately linearly with amplitude
- `l2rel` stays essentially constant across the sweep
- that is not the signature of unstable numerics
- it is the signature of a stable, repeatable model mismatch

## Numerical-stability assessment

Current evidence argues against a broad numerical-stability problem:

- the modal tip/tilt OPD is exact against OOPAO on the maintained grid
- the residual is stable and deterministic
- the residual scales smoothly with amplitude
- the maintained CPU / CUDA / AMDGPU parity surfaces are already green for the
  same runtime plants
- the earlier large composite mismatch was traced to a real semantic bug in
  `CompositeControllableOptic(..., DMReplace())`, and fixing that removed the
  structural error

Conclusion:

- the current composite OOPAO residual is best understood as an algorithm /
  convention / parameter-matching gap, centered on the DM approximation
- it is not strong evidence for package-wide numerical instability

## Practical implication for tolerances

The current committed composite OOPAO tolerances are evidence-based and
consistent with the measured sweep:

- `ShackHartmann`: `atol = 4e-4`, `rtol = 4e-2`
- `Pyramid`: `atol = 4e-4`, `rtol = 3e-2`
- `BioEdge`: `atol = 1e-3`, `rtol = 4e-2`

These are not arbitrary safety margins. They are aligned with the observed
cross-simulator residual envelope for this narrow composite plant.

## Recommended next work

If tighter external composite tolerances are desired, the next review should be
targeted, not package-wide:

1. improve the Julia-DM to OOPAO-DM parameter mapping
2. compare influence functions directly rather than only assembled OPD
3. check whether an OOPAO DM configuration closer to the Julia Gaussian width
   exists than the current `mechCoupling` proxy

What is not currently justified:

- a broad numerical-stability audit of the package
- loosening the committed tolerances beyond the measured envelope

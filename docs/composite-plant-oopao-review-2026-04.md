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

- fitted OOPAO `mechCoupling = 8.304545550591036e-2`
- previous closed-form proxy `= 8.465798862252998e-2`
- response-bank fit relative error `= 6.8064e-2`

Interpretation:

- the DM translation layer is now fitted explicitly against the Julia Gaussian
  influence model on the sampled pupil grid
- the fitted value only moved slightly from the older closed-form proxy, which
  confirms the proxy was already close on this maintained surface
- the remaining OPD discrepancy is still dominated by the DM family mismatch,
  not by a backend or arithmetic issue

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

- `scale=1.0`: `maxabs = 3.22e-4`, `l2rel = 3.6785e-2`

### Pyramid

- `scale=1.0`: `maxabs = 2.96e-4`, `l2rel = 2.6974e-2`

### BioEdge

- `scale=1.0`: `maxabs = 8.60e-4`, `l2rel = 3.4713e-2`

Interpretation:

- the residual envelope remains in the same narrow band as before the fit pass
- the fitted OOPAO coupling trims the composite residual slightly, but not dramatically
- that is consistent with the proxy already being close and with the remaining
  error coming from stable model-family differences rather than unstable numerics

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

The current committed composite OOPAO tolerances remain evidence-based and are
consistent with the updated fitted-mapping residuals:

- `ShackHartmann`: `atol = 4e-4`, `rtol = 4e-2`
- `Pyramid`: `atol = 4e-4`, `rtol = 3e-2`
- `BioEdge`: `atol = 1e-3`, `rtol = 4e-2`

These are not arbitrary safety margins. They are aligned with the observed
cross-simulator residual envelope for this narrow composite plant.

## Recommended next work

If tighter external composite tolerances are desired, the next review should be
targeted, not package-wide:

1. extend the fitted DM translation beyond the current representative `4 x 4`
   composite surface if broader external composite claims are desired
2. compare actuator-response families directly if a closer OOPAO DM family than
   `mechCoupling` becomes available
3. keep broader composite families on internal-artifact and backend-parity
   validation until they have a similarly narrow external reference surface

What is not currently justified:

- a broad numerical-stability audit of the package
- loosening the committed tolerances beyond the measured envelope

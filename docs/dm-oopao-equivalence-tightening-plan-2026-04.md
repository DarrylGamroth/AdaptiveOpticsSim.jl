# DM OOPAO Equivalence Tightening Plan 2026-04

Status: implemented

## Purpose

This plan tightens the external OOPAO equivalence story for deformable-mirror
and composite `tiptilt + dm` plants.

The previous composite OOPAO residual was already known to be deterministic and
backend-independent, but it still reflected a model-matching gap in the DM
translation layer. The goal here is to reduce that gap by improving the OOPAO
DM construction path rather than by loosening tolerances.

## Actions

### DOE-1 Align actuator ordering

Make the OOPAO coordinate list use the same actuator ordering as the Julia DM
command layout.

Status:

- implemented

### DOE-2 Replace the scalar proxy with a fitted OOPAO DM mapping

Fit the OOPAO `mechCoupling` parameter against the Julia Gaussian influence
model on the sampled pupil grid instead of relying on a single closed-form
proxy.

Status:

- implemented

### DOE-3 Record fit provenance in the frozen references

Store the fitted OOPAO coupling and fit residual in the generated reference
manifest so the external mapping is inspectable.

Status:

- implemented

### DOE-4 Rebaseline the composite OOPAO surfaces and tighten evidence notes

Regenerate the composite `tiptilt + dm` OOPAO references, rerun the reference
suite, and update the targeted review note to reflect the improved mapping and
resulting tolerance envelope.

Status:

- implemented

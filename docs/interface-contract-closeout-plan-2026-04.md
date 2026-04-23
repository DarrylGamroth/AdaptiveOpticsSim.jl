# Interface Contract Closeout Plan 2026-04

Status: completed

## Purpose

This pass closes the highest-value remaining interface debt outside the detector
subsystem.

The goal is not to redesign working algorithms. The goal is to make the
maintained contracts for WFS families, controllers, and calibration/reconstructor
workflows explicit enough that new families can be added without guessing at
state fields or diagnostics conventions.

## Scope

This plan covers three concrete slices:

- `ICC-1` WFS interface accessors and conformance
- `ICC-2` controller lifecycle and output accessors
- `ICC-3` calibration/reconstructor diagnostic accessors

## Actions

### `ICC-1` WFS Interface Accessors

Add maintained WFS-side accessors for the stable exported quantities already
used across the package:

- `slopes(wfs)`
- `valid_subaperture_mask(wfs)`
- `reference_signal(wfs)`
- `camera_frame(wfs)`

Add matching optional capability queries:

- `supports_valid_subaperture_mask`
- `supports_reference_signal`
- `supports_camera_frame`

Use those accessors in runtime/calibration code where that improves the
interface boundary without changing behavior.

Status:

- completed

### `ICC-2` Controller Contract Closeout

Make the controller family less field-driven by adding:

- `controller_output(ctrl)`
- `reset_controller!(ctrl)`
- `supports_controller_reset(ctrl)`

Keep `update!` as the canonical mutating step and ensure the maintained
controller returns `controller_output(ctrl)`.

Status:

- completed

### `ICC-3` Calibration/Reconstructor Contract Closeout

Add preferred long-form accessors for the calibration/reconstructor diagnostics
that had been exposed mainly through fields:

- `calibration_amplitude(imat)`
- `inverse_policy(x)`
- `singular_values(x)`
- `condition_number(x)`
- `effective_rank(x)`
- `truncation_count(vault)`

Apply this to:

- `InteractionMatrix`
- `CalibrationVault`
- `ModalReconstructor`
- `MappedReconstructor`

Status:

- completed

## Verification

The closeout is complete when:

1. the new accessors are exported and documented
2. runtime and calibration code use the new WFS accessors where they are part
   of the maintained contract
3. interface-conformance tests pin the new accessors and optional capabilities
4. the maintained detector/WFS/calibration/runtime test surfaces pass

## Outcome

The package now has clearer subsystem-level interface contracts for:

- WFS exported state
- controller lifecycle and output access
- calibration/reconstructor diagnostics

This reduces direct field access in maintained code and makes future model
families easier to add without another package-wide refactor.

# Detector Interface Extension Notes

Status: active

## Scope

This note records the maintained extension seams for the detector subsystem.

Current priority seams:

- `FrameReadoutProducts`
- detector execution plans
- detector response models

## Readout Products

`FrameReadoutProducts` is the public readout-export surface for detector capture
results beyond the main detector frame.

A readout-product subtype should participate in the public detector accessors by
defining only the methods it actually provides:

- `detector_reference_frame`
- `detector_signal_frame`
- `detector_combined_frame`
- `detector_reference_cube`
- `detector_signal_cube`
- `detector_read_cube`
- `detector_read_times`

All of these default to `nothing` for generic `FrameReadoutProducts`. A custom
product type should only override the accessors it needs.

This means a new readout-product family should not require edits to the central
generic detector layer just to expose one additional product layout.

Directly replacing `det.state.readout_products` with an arbitrary product type is
not part of the supported extension contract. The supported seam is the product
type plus its accessor methods.

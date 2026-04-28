# Detector Interface Extension Notes

Status: active

## Scope

This note records the maintained extension seams for the detector subsystem.

For the step-by-step authoring workflow, use
[adding-detectors.md](./adding-detectors.md).

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


## Detector Families And Assembly

The supported extension rule is:

- family-local metadata, capability predicates, and default model methods should
  live with the family/type definitions
- detector-instance assembly, validation orchestration, and export-metadata
  assembly stay centralized in `src/detectors/generic.jl`

Concretely:

- open extension seams:
  - `FrameReadoutProducts` accessor methods
  - detector/sensor family capability methods such as `supports_*`
  - family-local metadata helpers such as `*_symbol`
  - detector family default-model methods such as `default_response_model`
- centralized implementation surfaces:
  - detector construction and state assembly
  - detector export-metadata assembly
  - detector capture/readout orchestration

This means new detector families should normally add their family-specific
methods in the file that owns the family type, not by extending a central
registry table in `src/detectors/generic.jl`.

## Supported Vs Unsupported Mutation Surfaces

Supported:

- adding new `FrameReadoutProducts` subtypes plus accessor methods
- adding new detector/sensor-family methods through dispatch
- adding new response, defect, timing, nonlinearity, persistence, or thermal
  families through their model interfaces

Not supported as an external extension contract:

- mutating detector state fields directly to unrelated implementation types
- replacing centralized detector assembly helpers with external overrides
- relying on the concrete internal storage layout of `DetectorState`

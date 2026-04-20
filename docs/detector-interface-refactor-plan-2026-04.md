# Detector Interface Refactor Plan

Status: completed

## Goal

Make the detector implementation surface more open to extension without
weakening the current hot-path and backend contracts.

The immediate priority is the readout-product seam. That is the narrowest place
where the public detector interface still behaves like a closed registry instead
of an extensible dispatch surface.

## Problem Statement

The current detector interface is structurally sound in the large:

- params and state are separated
- detector execution plans already use backend-specialized dispatch
- the mutating capture/readout pipeline is explicit

But extensibility still has two clear bottlenecks:

1. Readout-product accessors are centralized and closed-world.
2. Detector capability and metadata helpers are concentrated in one generic file.

This plan addresses those in order, starting with the readout-product seam.

## Design Rules

- Keep multiple dispatch as the primary extension mechanism.
- Avoid `isa`-driven branching in package code.
- Preserve current detector behavior and exported API names.
- Keep hot-path state concrete where feasible, but do not force a larger state
  refactor into the first seam cleanup.
- Make each phase independently testable.

## Traceable Actions

### DIR-1 Readout-Product Seam

Goal:

- make `FrameReadoutProducts` accessors open by default
- remove detector-family-specific `isa` checks when preserving or updating
  readout products

Implementation:

- define default `detector_reference_frame`, `detector_signal_frame`,
  `detector_combined_frame`, `detector_reference_cube`,
  `detector_signal_cube`, `detector_read_cube`, and `detector_read_times`
  on `FrameReadoutProducts`
- keep product-specific overrides on the concrete product types
- switch HgCdTe readout buffering to those accessors instead of explicit type
  tests
- add a test that a custom `FrameReadoutProducts` subtype can participate in the
  public accessor surface without editing `src/Detectors/generic.jl`

Done when:

- no detector readout-product logic relies on closed-world product enumeration
- no `isa HgCdTeReadoutProducts` remains in the readout update path

### DIR-2 Readout-Product State Concreteness

Goal:

- remove the abstract `FrameReadoutProducts` slot from `DetectorState` where
  that can be done without breaking supported sensor behavior

Implementation candidates:

- keep a sensor-specific concrete product type in state and reset by value
- or use a small explicit union only if the concrete state model cannot be
  stabilized per detector family

Done when:

- detector state no longer stores `readout_products` as the abstract
  `FrameReadoutProducts` type

### DIR-3 Capability And Metadata Traits

Goal:

- shrink the central detector capability registry in `src/Detectors/generic.jl`

Implementation:

- move family-specific defaults closer to the family types
- keep generic fallbacks narrow and trait-like
- preserve current exported metadata accessors

Done when:

- adding a new detector family does not require touching a long central list of
  `supports_*` and `*_symbol` helpers unless that family genuinely changes
  shared generic policy

### DIR-4 Detector Execution-Plan Cleanup

Goal:

- remove remaining detector-path `isa` checks and route plan distinctions
  through dispatch or helper traits

Implementation:

- replace host-mirror special cases in capture/readout helpers with method
  specialization on plan/style

Done when:

- the detector package follows the package-level “no `isa` in package code”
  rule on maintained paths

### DIR-5 Documentation And Extension Guidance

Goal:

- document the supported way to extend detector families and readout products

Implementation:

- update the maintainer/developer docs with the detector extension seam
- record which detector surfaces are intended to be open extension points and
  which are internal implementation details

Done when:

- the detector interface and extension rules are explicit in docs instead of
  only discoverable from source

## Execution Order

1. `DIR-1` readout-product seam
2. `DIR-2` readout-product state concreteness
3. `DIR-4` execution-plan cleanup
4. `DIR-3` capability/metadata trait cleanup
5. `DIR-5` documentation closeout

This order is deliberate:

- `DIR-1` gives the highest extensibility payoff with the least disruption
- `DIR-2` is the main state-layout cleanup, but should follow the open seam
- `DIR-4` is local cleanup that should not block the readout seam
- `DIR-3` is broader and should be done after the product surface is stable

## Current Status

- `DIR-1`: completed
- `DIR-2`: completed for maintained frame-detector families
- `DIR-3`: completed; response-, sensor-, readout-correction-, thermal-, and counting-family metadata/capability helpers and generic detector-family defaults moved out of the central detector registry
- `DIR-4`: completed for maintained detector capture/readout plan routing
- `DIR-5`: completed

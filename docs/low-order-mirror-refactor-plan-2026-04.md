# Low-Order Mirror Refactor Plan 2026-04

Status: implemented

## Purpose

This plan generalizes the maintained low-order controllable-optic surface so the
package has one real implementation for low-order modal optics instead of
separate duplicated types plus a semantic alias.

## Actions

### LOM-1 Add a general low-order modal optic

- implement `LowOrderMirror` as the shared low-order controllable-optic type
- support direct mode-definition tuples
- support explicit mode matrices

Status:

- implemented

### LOM-2 Keep common user-facing convenience constructors

- keep `TipTiltMirror(...)` as a convenience constructor over `LowOrderMirror`
- keep `FocusStage(...)` as a convenience constructor over `LowOrderMirror`

Status:

- implemented

### LOM-3 Remove alias-only API surface

- remove `SteeringMirror`
- replace steering-specific validation/configuration uses with
  `LowOrderMirror(...; labels=:steering)` or `TipTiltMirror(...; label=:steering)`

Status:

- implemented

### LOM-4 Add configurable Zernike low-order support

- support `LowOrderMirror(tel; zernike_modes=..., ...)`
- allow grouped or per-mode command labels

Status:

- implemented

### LOM-5 Update docs, examples, and tests

- update public API docs and user guide
- add cookbook coverage for direct `LowOrderMirror` construction
- rebase tests and backend parity helpers on the generalized low-order optic

Status:

- implemented

## Verification

- focused runtime tests pass
- low-order self-check and backend parity surfaces continue to pass
- user-facing docs now describe `LowOrderMirror` as the general surface and
  `TipTiltMirror` / `FocusStage` as convenience constructors

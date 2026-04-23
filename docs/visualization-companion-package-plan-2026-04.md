# Visualization Companion Package Plan

Status: implemented

Date: 2026-04-23

## Purpose

This plan defines the maintained visualization direction for the
`AdaptiveOpticsSim` ecosystem.

The goal is to add practical plotting and display support without moving
plotting dependencies or visualization-side APIs into the simulation core.

The chosen direction is:

- keep `AdaptiveOpticsSim.jl` plotting-free
- add visualization in a companion package,
  tentatively `AdaptiveOpticsSimPlots.jl`
- use `Plots.jl` as the first maintained plotting backend
- keep the companion package thin and data-driven

## Decision Summary

Visualization should not be implemented as a required dependency of
`AdaptiveOpticsSim.jl`.

The right shape is:

- core package owns:
  - simulation data structures
  - stable state/export accessors
  - no plotting dependencies
- companion package owns:
  - plotting recipes and convenience helpers
  - optional display-oriented wrappers
  - `Plots.jl` integration

This follows the existing boundary policy in
[optional-integration-boundaries.md](./optional-integration-boundaries.md):

- keep the scientific core small
- keep optional integrations outside the default core path

## Why A Companion Package

This is the right boundary for three reasons:

1. plotting is not part of the AO model itself
2. plotting dependencies are large and change faster than the core simulation
   surface
3. visualization helpers often want notebook- and script-friendly conveniences
   that should not shape the core package API

This avoids the common failure mode seen in older AO toolkits where display
helpers become entangled with the modeling code.

## Package Boundary

Implemented package:

- `../AdaptiveOpticsSimPlots.jl`

Core package responsibilities:

- expose stable data accessors and metadata needed for plotting
- avoid plotting imports or optional extension logic in hot paths
- keep arrays and state queryable through maintained contracts

Companion package responsibilities:

- implement `Plots.jl`-based visualization helpers
- provide small high-level plotting functions for common AO surfaces
- optionally add recipes for core types once the helper surface stabilizes

The companion package may depend on:

- `AdaptiveOpticsSim`
- `Plots`
- optionally `Tables` later if log/telemetry plotting wants table-native input

It should not require:

- GPU packages
- notebook stacks
- GUI frameworks

## Maintained First Slice

The first maintained helper set should be intentionally small.

### Pupil And Phase

- `plot_pupil(obj; kwargs...)`
  - telescope pupil or boolean aperture mask
- `plot_opd(obj; kwargs...)`
  - OPD / phase map from telescope state or explicit array

Target inputs:

- `Telescope`
- `AbstractMatrix`

### PSF And Science Images

- `plot_psf(psf; kwargs...)`
- `plot_science_frame(obj; kwargs...)`

Target inputs:

- explicit PSF array
- runtime readout object
- detector output frame

### Detector And WFS Frames

- `plot_detector_frame(det; kwargs...)`
- `plot_wfs_frame(obj; kwargs...)`

Target inputs:

- `Detector`
- runtime readout object
- WFS types that expose `camera_frame(wfs)` or detector-backed exported pixels

### DM And Command Visualization

- `plot_dm_commands(dm; kwargs...)`
- `plot_dm_opd(dm, tel; kwargs...)`

Target inputs:

- `DeformableMirror`
- optionally `AbstractControllableOptic` for low-order/modal command vectors

For DMs with regular grid topology, command plots should render actuator maps
in 2-D. For sampled topologies, the first maintained fallback should be a
scatter-based actuator display rather than forcing everything onto a grid.

### Runtime Telemetry

- `plot_signal_trace(x; kwargs...)`
- `plot_runtime_timeseries(log; fields=... , kwargs...)`

Target inputs:

- explicit vectors
- simple named-tuple or struct logs
- later, package telemetry exports if a stable plotting surface emerges there

This first slice should stay simple. It should not try to become a telemetry
dashboard framework.

## Data Contract

The visualization package should consume maintained accessors rather than
family-local fields whenever possible.

Use these existing core surfaces first:

- `slopes(wfs)`
- `camera_frame(wfs)`
- `reference_signal(wfs)`
- `output_frame(det)`
- `detector_export_metadata(det)`
- `command(...)`
- `wfs_frame(...)`
- `science_frame(...)`
- `topology(dm)`
- `actuator_coordinates(dm)`
- `valid_actuator_mask(dm)`
- `topology_metadata(dm)`

If plotting requires a quantity that is not available through a maintained
accessor, the default action should be:

- improve the core accessor surface first
- not reach directly into undocumented fields from the plotting package

That is the main test of whether the plotting boundary is sound.

## API Shape

The first public API in `AdaptiveOpticsSimPlots.jl` is a small helper set plus
a package-owned `aoplot(...)` multiple-dispatch entrypoint, not a large recipe
system.

Maintained entry points:

- `aoplot(...)`
- `plot_pupil(...)`
- `plot_opd(...)`
- `plot_psf(...)`
- `plot_detector_frame(...)`
- `plot_wfs_frame(...)`
- `plot_dm_commands(...)`
- `plot_signal_trace(...)`

These helpers should:

- return the `Plots.jl` plot object
- accept existing arrays or core objects
- avoid mutating core state
- avoid hidden simulation steps

What they should not do:

- call `measure!`, `step!`, or `capture!` implicitly
- allocate new simulation objects
- change runtime state for display

## Recipes Vs Functions

Do not start with a full recipe-centric design.

The implemented sequence was:

1. first implement explicit plotting helpers
2. once the helper signatures stabilize, optionally add `RecipesBase` or
   `Plots` recipes for a small subset of core types

Why:

- explicit functions are easier to reason about
- recipes are harder to keep stable across many subsystem types
- the first goal is clear utility, not generic plotting magic

## What Not To Add In Phase 1

Do not add these in the first maintained slice:

- animated movie builders
- interactive GUIs
- dashboard/server code
- Makie support as a required second backend
- notebook-only helpers
- plotting side effects inside tutorials or runtime code

Those can come later if the core helper surface proves stable and useful.

## Phased Implementation

### VIZ-1: Package Scaffold

Completed in `../AdaptiveOpticsSimPlots.jl`:

- dependency on `AdaptiveOpticsSim`
- dependency on `Plots`
- minimal test and CI scaffold
- one README focused on examples, not architecture prose

### VIZ-2: Core Helper Set

Implement the first maintained helper set:

- pupil
- OPD
- PSF
- detector frame
- WFS frame
- DM commands
- simple signal traces

### VIZ-3: Contract Audit

Use the plotting package as an interface test.

For every helper, check whether it can be implemented using maintained
accessors only.

If not:

- add or tighten the relevant accessor in core
- do not normalize plotting against private fields

### VIZ-4: Docs And Examples

Completed additions:

- one short cookbook note in `AdaptiveOpticsSim.jl`
- a focused README in `AdaptiveOpticsSimPlots.jl`
- one or two small example scripts showing:
  - PSF / phase visualization
  - WFS / DM visualization

### VIZ-5: Optional Recipe Layer

Only after the function surface stabilizes:

- consider `Plots` recipes for a narrow subset of types

This should remain optional and subordinate to the explicit helper API.

## Validation

The companion package should have its own tests.

Minimum validation:

1. each helper accepts the documented inputs
2. each helper returns a `Plots` plot object without mutating simulation state
3. regular-grid and sampled-topology DM displays both work
4. detector- and WFS-backed helpers use maintained accessors only

The core package should not gain plotting-specific tests beyond any accessor
surfaces needed to support the plotting package.

## Success Criteria

This plan is successful when:

1. visualization exists without adding plotting deps to the core package
2. the first helper set is actually useful for AO users
3. the plotting package depends mostly on maintained accessors rather than
   undocumented fields
4. adding the plotting package exposes any missing core accessors cleanly
   instead of requiring another broad refactor

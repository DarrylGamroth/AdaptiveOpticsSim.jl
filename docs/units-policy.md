# Units Policy

This document defines the maintained project-wide units policy for
AdaptiveOpticsSim.jl.

## Core Rules

- Internal runtime OPD fields are in meters.
- Internal wavelength-dependent optical phase is in radians.
- Pure internal angles and rotations are in radians.
- Human-facing angular API inputs may remain in convenient units such as:
  - degrees for orientation-like quantities,
  - arcseconds for sky offsets and field-of-view quantities.
- API-boundary quantities that use non-SI or non-radian angular units must say
  so explicitly in their names or documentation, for example:
  - `rotation_deg`
  - `field_stop_size_arcsec`
  - `fov_arcsec`

## Internal Normalization Rules

- Constructor and builder paths should convert boundary-degree inputs to radians
  immediately.
- Constructor and builder paths should normalize sky-coordinate inputs into the
  internal representation used by the hot path instead of repeatedly converting
  them during execution.
- Runtime kernels and hot loops should use:
  - `sin` / `cos` / `sincos` on radian values,
  - Cartesian offsets where that is the natural working representation,
  - precomputed direction components where repeated trigonometric evaluation is
    avoidable.
- Internal fields that store radians should use a `_rad` suffix unless the
  enclosing type/documentation already fixes the unit unambiguously.

## Maintained Exceptions

- Sky offsets may be stored internally as Cartesian arcsecond offsets when that
  is the natural working quantity for image shifts or atmospheric footprint
  extraction.
- API-facing angle fields may continue to accept degrees where that matches
  established user expectations, but those values should not remain degree-based
  in the runtime state.

## Current Cleanup Targets

The first cleanup pass covered by this policy is:

- `Source` / `LGSSource`: normalize polar sky coordinates to Cartesian arcsecond
  offsets at construction time.
- `Misregistration`: store rotations internally in radians while keeping
  degree-based constructor keywords at the API boundary.
- `TomographyAtmosphereParams`: store zenith and wind directions internally in
  radians while keeping degree-based constructor keywords at the API boundary.

## Non-goals

This policy does not require:

- changing all user-facing angular APIs to radians immediately,
- removing `*_deg` or `*_arcsec` constructor keywords where they are already the
  clearest boundary contract,
- forcing every internal directional quantity into radians if Cartesian
  components are the better runtime representation.

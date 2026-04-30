# Extension Guide

Status: active

This guide describes the maintained extension seams for adding new detector
families, wavefront sensors, deformable mirrors, controllable optics,
controllers, and reconstructors.

Prefer small concrete types plus multiple dispatch. Do not add central
type-switching registries for new families.

## Source Layout

Use lower-case directory names for new source-tree locations. Julia type names
may use CamelCase, but package directories should remain lower-case subsystem
names such as `src/wfs`, `src/detectors`, `src/optics`, or `src/control`.

## Detectors

Physical detector families live in `src/detectors/`.

Use this split:

- immutable parameter/configuration structs for physical constants and readout
  configuration
- mutable state structs for accumulated charge, readout buffers, and thermal or
  persistence state
- `!` methods for capture, integration, correction, and readout hot paths
- shared layers for reusable behavior such as quantization, frame response,
  multi-read readout, thermal state, and counting statistics

New detector families should implement the smallest family-owned methods needed
to connect to the generic detector pipeline. If a feature can apply to multiple
families, put it in a shared layer instead of copying it into each sensor.

Use `bits` for quantization depth and `output_type` for the Julia element type
returned at a HIL/RTC boundary.

## Wavefront Sensors

WFS families live in `src/wfs/`. Family-specific implementation should stay
near the family, usually in the corresponding directory under `src/wfs/`.

New WFS types should provide:

- a concrete sensor type subtype of the package WFS abstraction
- setup/precomputation methods owned by the WFS family
- measurement methods that write into preallocated state
- signal extraction methods for slopes, intensities, or family-specific
  products
- detector image formation if the sensor has a maintained detector-facing path

The generic runtime should call WFS-owned seams rather than branching on sensor
families directly.

## Deformable Mirrors

DM implementation lives in `src/optics/deformable_mirror.jl`.

Use the composed model:

- topology describes actuator layout and coordinates
- influence basis describes how commands map to sampled surface modes
- actuator behavior describes clipping, dead actuators, slaving, or future
  dynamic behavior
- `DeformableMirror` composes those pieces into the runtime optic

Measured manufacturer influence functions should enter as sampled influence
bases with explicit metadata. Do not force measured data through an anonymous
dense matrix unless it is truly only an escape hatch.

## Controllable Optics

Use controllable optics for modal or low-order optical surfaces that are not
better represented as a DM technology family.

Common examples:

- tip/tilt surfaces
- low-order modal optics
- Zernike-mode control surfaces
- calibration or NCPA compensation surfaces

Keep command application explicit and ordered in the plant. Runtime code should
operate on the generic controllable-optic interface.

## Controllers And Reconstructors

Controllers own temporal behavior. Reconstructors own slopes-to-command or
measurement-to-command operators.

Use this split:

- reconstructors expose the command operator and any diagnostics needed by
  calibration or validation
- controllers own integrator state, gains, leakage, saturation, or future
  temporal filters
- runtime code calls generic accessors and `!` methods

Prefer concrete state and workspace fields. Avoid hidden allocations in the
per-step control path.

## Backend And Allocation Rules

Extension code should follow the package-wide backend rules:

- parametrize by numeric type and array backend where practical
- avoid scalar indexing on GPU arrays
- preallocate workspaces for hot paths
- use explicit `!` methods when mutating state
- centralize RNG ownership in the caller or workspace
- use multiple dispatch or traits instead of `isa` chains

If an extension needs a new reusable behavior, add a small generic seam and one
family implementation first. Then add the second implementation when another
family actually needs it.

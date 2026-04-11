# Grouped Runtime Contract

Status: active

Plan traceability:

- [`GR-1`](./grouped-runtime-plan.md)

## Purpose

This note defines the maintained grouped-runtime contract for composite runtime
execution.

It applies to:

- grouped source execution inside maintained WFS paths
- grouped branch execution through
  [`CompositeSimulationInterface`](../src/Control/runtime/types.jl)
- grouped exported runtime products

## Grouped Source Semantics

Grouped source execution means:

- one WFS consumes multiple logically related source samples
- the grouped path performs the shared optical work once per sample and then
  reduces the sample stack into the exported measurement state

Maintained grouped source families are queried through:

- `supports_grouped_execution(wfs, src)`

Current maintained grouped families:

- diffractive `ShackHartmann` with `Asterism`, `SpectralSource`, and
  `ExtendedSource`
- diffractive `PyramidWFS` with `Asterism`, `SpectralSource`, and
  `ExtendedSource`
- diffractive `BioEdgeWFS` with `Asterism`

## Grouped WFS Product Ownership

The grouped WFS contract distinguishes:

- internal grouped scratch stacks
  - not exported
  - examples:
    - `intensity_tmp_stack`
    - grouped asterism stacks
- branch-local exported frames
  - copied into `SimulationInterface` or `CompositeSimulationInterface`
- grouped exported stacks
  - optional compatible-shape aggregation owned by
    `CompositeSimulationInterface`

Required rule:

- scratch state must never be treated as an exported runtime product

## Grouped Detector Export Rules

For grouped composite execution:

- branch detector-facing WFS exports are available through
  `wfs_frame(composite)`
- compatible grouped WFS exports are available through
  `grouped_wfs_stack(composite)`
- if branch shapes are incompatible, grouped stack export must be `nothing`
  instead of forcing an ambiguous layout

## Grouped Science Export Rules

The same policy applies to grouped science outputs:

- branch science frames via `science_frame(composite)`
- compatible grouped science stacks via
  `grouped_science_stack(composite)`
- incompatible grouped science layouts must remain unstacked

## Current Runtime Entry Points

The maintained grouped runtime entry points are:

- `CompositeSimulationInterface`
- `step!(::CompositeSimulationInterface)`
- `wfs_frame`
- `science_frame`
- `grouped_wfs_stack`
- `grouped_science_stack`

## Conformance Checklist

- grouped execution support is declared explicitly with
  `supports_grouped_execution`
- grouped scratch buffers are distinct from exported grouped products
- compatible grouped branch frames can be stacked without aliasing branch
  export buffers
- incompatible grouped branch frames do not pretend to share a common stacked
  layout

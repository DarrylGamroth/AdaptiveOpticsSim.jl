# Proper.jl Integration Guide

Status: active

## Purpose

This guide records the maintained boundary between `AdaptiveOpticsSim.jl` and
`Proper.jl` for HIL-style coronagraph and science-arm simulations.

Use this when:

- `AdaptiveOpticsSim.jl` owns AO runtime state, RTC commands, WFS products, and
  DM-to-OPD modeling
- `Proper.jl` owns an external wave-optics science prescription
- both packages may run on CPU, CUDA, or AMDGPU without unnecessary host/device
  transfers

## Recommended Ownership

Keep the packages loosely coupled.

- `AdaptiveOpticsSim.jl` should receive RTC commands and update AO runtime state.
- `AdaptiveOpticsSim.jl` should convert actuator commands into a sampled OPD
  surface when its DM model is being used.
- `Proper.jl` should receive the sampled OPD and pupil arrays needed by the
  science prescription.
- Neither package should depend directly on the other for core functionality.
  Keep integration examples and benchmarks at the application boundary unless a
  stable shared package is justified.

The preferred seam is:

```julia
payload = CoronagraphPayload(
    sim.tel.state.opd,
    pupil_mask(sim.tel),
    sim.tel.params.diameter,
    focal_length_m,
    lyot_stop_norm,
)

psf, sampling = prop_run(science_model; payload=payload)
```

The `payload=...` keyword is preferred for new Julia-native integrations.
`PASSVALUE` is a PROPER compatibility adapter and should be kept for upstream
ports or parity harnesses, not used as the default HIL interface.

## Prescription Shape

Use a typed payload and an ordinary keyword argument:

```julia
struct CoronagraphPayload{O,P,T}
    opd_m::O
    pupil::P
    diameter_m::T
    focal_length_m::T
    lyot_stop_norm::T
end

function coronagraph_prescription(λm, n; payload::CoronagraphPayload)
    wf = prop_begin(payload.diameter_m, λm, n; beam_diam_fraction=1.0)
    prop_multiply(wf, payload.pupil)
    prop_add_phase(wf, payload.opd_m)
    return prop_end(wf)
end
```

Use the direct `Proper.prop_dm(wf, dm_surface)` path only when the external
science prescription should own a DM surface internally. If the AO runtime has
already applied the DM, pass the total sampled OPD instead.

## Conventions To Validate

Before treating a new instrument integration as supported, validate these
conventions with small deterministic cases:

- **Units:** the OPD handoff is in meters.
- **Sign:** a positive OPD perturbation should produce the expected science
  response for piston, tilt, focus, and one actuator poke.
- **Centering:** array center conventions should match between the telescope
  grid and the Proper wavefront grid.
- **Orientation:** row/column orientation should be verified with asymmetric
  maps, not only circular pupils.
- **Pupil ownership:** pass the actual `AdaptiveOpticsSim.jl` pupil mask when
  spiders, segment gaps, central obstruction, or custom reflectivity matter.
- **Backend ownership:** if both sides use the same GPU backend, keep arrays on
  device and construct `Proper.RunContext(typeof(sim.tel.state.opd))`.

These are integration-contract checks, not optional polish. They prevent a
model from looking numerically plausible while using the wrong sign or grid
orientation.

## Runtime Guidance

For repeated HIL execution:

1. Build and prepare the AO scenario once.
2. Build and prepare the Proper science model once.
3. Reuse a typed payload whose arrays point at current AO runtime state.
4. Each frame, update the scenario command, run `sense!` or `step!`, then call
   `prop_run(science_model; payload=payload)`.
5. Benchmark the combined command-to-pixels path on each claimed backend.

Avoid host transfers in the frame loop. Use `Array(...)` only as an explicit
boundary when the AO runtime and science model intentionally live on different
backends.

## Installation For Examples

`Proper.jl` is intentionally not a dependency of `AdaptiveOpticsSim.jl`.
Install it into the active example or benchmark environment before running the
companion scripts:

```julia
using Pkg
Pkg.add("Proper")
```

For a local sibling checkout, use development mode instead:

```julia
using Pkg
Pkg.develop(path="../proper.jl")
```

## Current Companion Files

- [`model-cookbook.md`](./model-cookbook.md), Recipe 8, shows the user-facing
  pattern.
- [`../examples/integrations/proper_hil_coronagraph.jl`](../examples/integrations/proper_hil_coronagraph.jl)
  is the runnable example.
- [`../examples/support/proper_hil_coronagraph_common.jl`](../examples/support/proper_hil_coronagraph_common.jl)
  contains the shared setup and prescription.
- [`../scripts/profile_proper_hil_coronagraph.jl`](../scripts/profile_proper_hil_coronagraph.jl)
  is the current seam benchmark.

## When Not To Use This Seam

Use native `AdaptiveOpticsSim.jl` science products instead when the science arm
is simple PSF generation or detector simulation already covered by the package.
Use `Proper.jl` when the science prescription needs PROPER-compatible
propagation, coronagraph masks, or a model already written in PROPER style.

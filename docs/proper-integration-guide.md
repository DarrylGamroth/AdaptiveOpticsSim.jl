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
- `AdaptiveOpticsSim.jl` should own common AO-path aberrations and sampled
  path-specific NCPA when they are represented on its native planes.
- `Proper.jl` should receive the caller-owned path OPD/field and prepared pupil
  geometry needed by the science prescription; it should not receive or retain
  mutable telescope path storage.
- Neither package should depend directly on the other for core functionality.
  Keep integration examples and benchmarks at the application boundary unless a
  stable shared package is justified.

The preferred seam renders the current published atmosphere epoch and visible
optic surfaces into a caller-owned science-path product before calling Proper:

```julia
using AdaptiveOpticsSim.Optics

science_pupil = PupilFunction(sim.tel)
science_renderer = prepare_atmosphere_renderer(
    sim.atm, sim.tel, science_source(sim))
render_atmosphere!(science_pupil, science_renderer, sim.atm,
    current_epoch(sim.atm))
update_surface!(sim.optic)
apply_surface!(science_pupil, sim.optic, DMAdditive())

payload = CoronagraphPayload(
    opd_map(science_pupil),
    pupil_amplitude(science_pupil),
    sim.tel.params.diameter,
    focal_length_m,
    lyot_stop_norm,
)

external_image, sampling = prop_run(science_model; payload=payload)
```

This does not advance atmosphere time and may use a different prepared source
direction from the WFS path. A later HIL path executor can own and update the
same `science_pupil` directly. On return, the integration adapter must construct a
caller-owned `IntensityMap` whose metadata declares physical sampling,
centering/orientation, wavelength, backend/device, normalization, spatial
measure, and incoherent-addition policy. `prepare_detector_acquisition` accepts
that result only when it is a photon-arrival-rate product, or when a
dimensionless result has an explicit prepared photon-rate scale. The current
raw-array example is therefore an integration proving ground, not evidence that
an undeclared PROPER array has physical detector units. Detector exposure
duration is never folded into the payload, PROPER result, or conversion scale.
If `prop_run` returns fresh storage, copy it into the prepared map; do not replace
the map storage bound to its detector-acquisition plan.

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

## Optional Proper Graph Node

`Proper.jl` remains outside the AdaptiveOpticsSim core dependency graph. The
optional `AdaptiveOpticsProperHIL.jl` companion exposes a
[`Proper propagation graph node`](./glossary.md):

```julia
using AdaptiveOpticsProperHIL
using AdaptiveOpticsSim
using AdaptiveOpticsSim.AlgorithmGraphs

configuration = ProperPropagationConfiguration(
    coronagraph_prescription;
    target=AdaptiveOpticsSim.Backends.HostComputeDevice(),
    resolution=128,
    diameter_m=8.0f0,
    wavelength_um=1.65f0,
)
node = proper_propagation_node(:science_propagation, configuration)
```

The maintained node has two complete-frame inputs, `pupil_opd` and
`pupil_amplitude`, and one output. Its configuration fixes the prescription,
resolution, physical diameter, wavelength, random seed, port schemas, numeric
type, backend, and exact device before execution. The default output schema is
normalized intensity. A different schema may be declared only when the
prescription and downstream integration actually satisfy that contract; the
adapter does not invent a photon-rate scale or detector exposure.

The adapter preserves the package ownership model:

- its immutable execution plan contains only the run-immutable prescription
  contract
- evolving random state has a distinct persistent owner and resets to the
  configured seed
- the `Proper.RunContext`, wavefront field, prescription-specific assets, and
  reusable scratch belong to the replaceable workspace owner
- `AlgorithmGraphs` preparation binds the final admitted graph input and output
  arrays exactly, so a graph step performs no frame copy at that node boundary

A prescription that needs masks, propagation plans, or other prepared assets
may specialize `prepare_proper_assets`. It must return a concrete named tuple;
`reset_proper_assets!` resets any replaceable asset state at a serialized graph
boundary. The prescription remains an ordinary Julia function and receives the
prepared `pupil_opd`, `pupil_amplitude`, `diameter_m`, field, wavefront, output,
and run context as keyword arguments.

One prepared Proper node is homogeneous in numeric type, array backend,
and exact device. A complete algorithm graph may connect different element
types through explicit conversion algorithms, but every array in one graph is
currently placed on one exact compute target. Neither the node nor the
graph inserts a fallback or transfer. A GPU-produced science frame therefore
reaches a CPU RTC or PipeWireAO client only through an explicit
application-owned device-to-host copy or prepared handoff after successful
graph publication.

Focused companion tests exercise the same node on CPU, local AMDGPU,
and WSL CUDA hardware with GPU scalar indexing disabled. They establish
CPU/GPU numerical agreement, exact residency, graph buffer identity, explicit
host publication, and zero warmed Julia heap allocation on the CPU path. They
do not establish coronagraph fidelity, SPIDERS end-to-end validity, service
time, acquisition cadence, or asynchronous failure isolation.

The node is usable today from Julia scripts and
`AdaptiveOpticsSim.AlgorithmGraphs`. Direct execution of Julia prescriptions in
PipeWire remains a separate Julia-host/module adapter problem. This graph-node
adapter does not turn Julia code into a native FGN shared-library plugin or
make Calculon part of the AOS API.

## NCPA Ownership

Use the native `NCPA` or `OPDMap` model for a static or slowly varying
aberration that is adequately represented as sampled pupil OPD. Apply it only
to the branch that contains the aberration. In a Plant, wrap it in
`Plant.SampledAberrationDefinition` with exact pupil-plane metadata, optical
placement, and `Plant.SelectedPathVisibility` for the affected science path.
`prepare_plant` makes the run-owned copy and path-local coupling, so do not
also mutate path storage manually or apply the same OPD inside the external
prescription.

Keep the NCPA inside the `Proper.jl` prescription when its behavior depends on
the detailed instrument relay, wavelength-dependent surfaces, amplitude
effects that depend on that relay, coronagraph planes, or propagation between
several physical optics. Do not also add a collapsed native NCPA map for the
same surfaces.

A useful HIL compromise is to derive a sampled NCPA surrogate from a detailed
`Proper.jl` model, validate the surrogate over the required wavelength and
field range, and use it on a high-rate native path. The full prescription can
still execute at science cadence or offline without blocking a WFS deadline.
When a prepared external executor is selected, it receives the pupil after
atmosphere materialization, native sampled aberrations, controllable surfaces,
and autonomous path optics; the core does not import or schedule `Proper.jl`.

## Conventions To Validate

Before treating a new instrument integration as supported, validate these
conventions with small deterministic cases:

- **Units:** the OPD handoff is in meters.
- **Radiometry:** declare whether the returned array is photon irradiance,
  cell-integrated photon rate, normalized intensity, contrast, or another
  quantity. Put that contract in `OpticalPlaneMetadata`; use an explicit
  normalized-to-photon-rate scale only when the conversion is physically
  defined. Validate conversion, represented cell measure, and photon
  conservation before detector acquisition with a non-unit exposure.
- **Spectral coordinates:** wavelength-dependent results may be summed by array
  index only when their physical focal grids are compatible and the declared
  combination is incoherent; otherwise retain a bundle or use an explicit
  prepared mapping.
- **Sign:** a positive OPD perturbation should produce the expected science
  response for piston, tilt, focus, and one actuator poke.
- **Centering:** array center conventions should match between the telescope
  grid and the Proper wavefront grid.
- **Orientation:** row/column orientation should be verified with asymmetric
  maps, not only circular pupils.
- **Pupil ownership:** pass the actual `AdaptiveOpticsSim.jl` pupil amplitude
  when spiders, segment gaps, central obstruction, or custom reflectivity
  matter.
- **Backend ownership:** if both sides use the same GPU backend, keep arrays on
  device and construct `Proper.RunContext(typeof(opd_map(science_pupil)))`.

These are integration-contract checks, not optional polish. They prevent a
model from looking numerically plausible while using the wrong sign or grid
orientation.

## Runtime Guidance

For repeated HIL execution:

1. Build and prepare the AO scenario once.
2. Build and prepare the Proper science model once.
3. Reuse a typed payload whose arrays point at the current caller-owned science-
   path product and whose geometry/radiometry metadata was validated during
   preparation.
4. Prepare the returned intensity-map and detector-acquisition contract once;
   retain incompatible wavelength grids as separate products unless an
   explicit mapping is prepared.
5. Each frame, update the scenario command, run `sense!` or `step!`, call
   `prop_run(science_model; payload=payload)`, update the caller-owned map, and
   execute the prepared detector acquisition.
6. Benchmark the combined command-to-pixels path on each claimed backend.

For a portable complete-frame graph, replace steps 2 and 5 with preparation of
the optional Proper graph node and one `step_graph!` call. Keep the complete
Proper path atomic: the prescription may evaluate every required optical plane
internally, while the graph publishes only its completed output frame.

Avoid host transfers in the frame loop. Use `Array(...)` only as an explicit
boundary when the AO runtime and science model intentionally live on different
backends.

## Installation For Examples

`Proper.jl` is intentionally not a dependency of `AdaptiveOpticsSim.jl`.
Install it from its maintained GitHub repository into the active example or
benchmark environment before running the companion scripts:

```julia
using Pkg
Pkg.add(url="https://github.com/DarrylGamroth/Proper.jl.git")
```

Validation and benchmark environments should pin an exact source revision and
retain their manifest. For local package development, a sibling checkout may
instead use development mode:

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
is direct photon-arrival-rate image formation or detector simulation already
covered by the package.
Use `Proper.jl` when the science prescription needs PROPER-compatible
propagation, coronagraph masks, or a model already written in PROPER style.

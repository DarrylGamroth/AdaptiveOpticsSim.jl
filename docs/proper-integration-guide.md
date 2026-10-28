# Proper.jl Integration Guide

Status: active

## Purpose

This guide records the maintained boundary between `AdaptiveOpticsSim.jl` and
`Proper.jl` for HIL-style coronagraph and science-arm simulations.

Use this when:

- `AdaptiveOpticsSim.jl` owns plant state, WFS products, application of complete
  PDM commands, and DM-to-OPD modeling
- `Proper.jl` owns an external wave-optics science prescription
- both packages may run on CPU, CUDA, or AMDGPU without unnecessary host/device
  transfers

## Native declarative scalar optical trains

AOS provides the optional `AdaptiveOpticsSimProperExt` extension. Load `Proper`
and use the same `AlgorithmGraphs.load_algorithm_graph` and definition/component
tooling as other AOS instruments. YAML remains an independently optional I/O
extension; neither Proper nor YAML is a core dependency.

```julia
using AdaptiveOpticsSim.AlgorithmGraphs, Proper, YAML

graph = prepare_algorithm_graph(load_algorithm_graph("science.yaml";
    bindings=(pupil_opd=opd_m, pupil_amplitude=amplitude)))
step_graph!(graph)
pitch_m = optical_sampling_m(AlgorithmGraphs.prepared_graph_node(graph, :optics))
```

The built-in node types are `proper_optical_train_f64` and
`proper_optical_train_f32`. They use the same fixed train schema and differ in
input/output precision. An application constructing nodes directly can call
`proper_optical_train_node_factory(; T=Float64)(name, config, (;))`. Without
`using Proper`, construction raises `AlgorithmGraphError`.

```yaml
schema_version: 1
name: science
nodes:
- name: optics
  type: proper_optical_train_f64
  config:
    resolution: 256
    pupil_resolution: 128
    diameter_m: 0.012
    wavelength_um: 1.55
    steps:
    - {op: lens, focal_length_m: 0.2, name: focusing lens}
    - {op: propagate, distance_m: 0.2}
    - op: circular_mask
      mode: transmit
      circles:
      - {radius_m: 0.0005, x_m: 0.0, y_m: 0.0}
inputs:
- {name: opd, destination: optics.pupil_opd, binding: pupil_opd}
- {name: amplitude, destination: optics.pupil_amplitude, binding: pupil_amplitude}
outputs:
- {name: intensity, source: optics.output}
```

`resolution`, `diameter_m`, `wavelength_um` and `steps` are required.
Propagation `resolution` must be even: Proper 0.4 does not preserve the declared
central sample on odd output grids. Odd pupil sizes are supported on even grids.
`pupil_resolution` defaults to `resolution` and must fit inside it. Centered
input OPD and field amplitude are zero-padded with matching central samples,
without resampling. Input spacing is `diameter_m / pupil_resolution`; the
beam-diameter fraction is fixed to `pupil_resolution / resolution`. Output is
always the square propagation grid. Each scalar operation optionally accepts
`name` as provenance:

| Operation | Required parameters |
| --- | --- |
| `propagate` | signed `distance_m` |
| `lens` | nonzero signed `focal_length_m` |
| `circular_aperture`, `circular_obscuration` | positive `radius_m` |
| `multiply` | field `amplitude` in [0,1] |
| `circular_mask` | `mode: transmit` or `reflect`; `circles` containing positive `radius_m`, `x_m`, `y_m` |
| `sampled_mask` | Cold `amplitude` binding role, optional `phase` role in radians, `sampling_m: [dy,dx]`, `center_pixel: [row,column]`, `outside_amplitude`, and `outside_phase_rad` |

Circle coordinates are in metres in the current optical plane: x increases
with columns and y with rows; the origin is sample `fld(resolution,2)+1` on
each axis. Circle masks form a hard-edge union. Reflection uses its binary
complement, so complementary branches conserve power at mask samples.
Preparation traces scalar geometry to create masks at each plane's sampling.

Sampled masks use ordinary node construction bindings, not filenames in the
propagation node. The node's `bindings` maps each role to an application-loaded
array. The cold factory validates and snapshots real amplitude in `[0,1]` and
optional phase in radians into fixed-size storage. Logical axes are row=y,
column=x; `center_pixel` is the one-based optical-axis location in that array.
Preparation performs bilinear interpolation of complex transmission at the
actual plane sampling, with the explicit outside field contributing per neighbor.
It does not interpolate phase angles or apply antialias filtering. Source
arrays and files can be released after construction. Stepping performs no I/O.

The shared instrument loader supports physical `optical_masks.<name>.<port>`
responses for `transmitted` and `reflected` ports. Each response declares
`amplitude` and optional `phase` FITS references, `sampling_m`, `center_pixel`,
`outside_amplitude` and `outside_phase_rad`. Each FITS reference specifies
`path`, `sha256`, `hdu`, `axes`, `type` and logical `shape`. The binding type is
Float32 or Float64; it need not equal the source FITS BITPIX type. A response
with phase data requires a matching `wavelength_m`; no chromatic phase scaling
is inferred. Select it with `{kind: fits, mask: <name>}` for `apodizer`,
`focal_plane` or `lyot` in the shared simulation mask settings. The Lyot role
uses the physical optical path's port; measured reflection is not synthesized
as a complement. Selected mask files join the generated provenance receipt.
Every declared circle must contain at least one represented sample on that
plane; an off-grid or unresolved circle raises a preparation error. This
numerical support check does not establish convergence or incident signal.
Unknown fields, unsupported operations, nonfinite values and numerical booleans
are rejected. Configuration never evaluates Julia source.

Ports default to `org.adaptiveopticssim.pupil-opd/1`,
`org.adaptiveopticssim.pupil-amplitude/1`, and
`org.adaptiveopticssim.optical-intensity/1`; override them with
`pupil_opd_schema`, `pupil_amplitude_schema`, and `output_schema`. The fixed
operations are deterministic and have no random-state configuration.

Output is squared field amplitude without per-path normalization, photon-rate
scaling or detector integration. After a graph step, `optical_sampling_m(owner)`
returns square physical sample width in metres on both axes. For a scalar train,
query after a completed step: this reads the retained wavefront's current sampling.
CUDA preparation traces the train while cold to initialize FFT plans and stable
field addresses, so its geometry may already describe the output before replay.
Detector pixel integration, exposure timing,
noise, throughput and photometry belong to their declared AOS components.

Optional `channels` changes the node to polychromatic propagation. Each of
1–64 channel mappings requires `wavelength_um`, nonnegative `weight`, and its
complete `steps`; optional `transmission` defaults to one and lies in `[0,1]`.
Weights have a finite positive sum and are normalized as photon-number source
fractions. Transmission multiplies these fractions without renormalization.
The top-level wavelength and steps remain required and define the reference
physical detector grid by default; they add no separate intensity contribution. Resolution,
pupil sampling, diameter, port schemas and mask-role bindings remain common.
Bindings include mask roles used by the reference or any channel train.

Optional `output_sampling_m` declares a finite positive square output-cell pitch
in metres, overriding reference-train sampling. This keeps the physical detector
grid fixed when the selected central wavelength changes. With this field and no
`channels`, the node maps one channel using the top-level wavelength and steps
with unit source weight and transmission. Both CPU and CUDA execute this mapping;
without this field or channels, scalar behavior remains unchanged. The value must
remain finite and positive at the selected Float32 or Float64 precision.

Each channel propagates the same entrance OPD and amplitude through its full
scalar train, preserving coherent interference within that field. Channels are
added incoherently after conservative physical cell-overlap mapping. For each
axis, a sample at `(i-fld(resolution,2)-1)*pitch` owns a cell extending half a
pitch on either side. The two-dimensional mapping coefficient is the overlap
area divided by source cell area, treating the scalar sample's power as uniform
within that cell. Fully covered source power is conserved; power beyond the
reference support is discarded without normalization. This is a finite-volume
approximation, so sampling convergence must be checked for sharp PSFs or masks.
Identical grids use exact direct addition. Preparation traces reference and
channel geometry and freezes overlap coefficients; stepping performs no I/O.
`optical_sampling_m` for this owner returns the declared fixed pitch, or reference
output sampling when no fixed pitch is supplied, even before the first step.
CUDA resampling uses prepared sparse destination tables and two ordered gather
kernels, with one writer per destination cell and no floating-point atomics.
The existing scalar configuration and output arithmetic are unchanged.

This is a centered scalar approximation. It does not infer mask
fabrication data, global fold orientation, polarization or off-axis aberrations.
Instrument data owns those assumptions and calibration identity. CPU numerical
checks do not qualify acquisition deadlines or physical devices.

Fixed scalar and polychromatic trains admit CUDA `CapturedGraphExecution` through
the normal AOS graph executor. Shape, buffer identities, wavelength samples,
propagation steps, masks and resampling geometry remain fixed until rebuilding;
device OPD and amplitude contents may change each frame. Cold preparation retains
FFT plans, stable field storage and separate complex-mask shift scratch. The
outer AOS graph captures the complete pipeline; no nested Proper capture is used.
AMDGPU capture is not qualified by this change.

The hardware selectors `test/cuda_proper_optical_train.jl` and
`test/cuda_proper_spectrum.jl` compare CPU, streamed CUDA and captured replay on
Float32/Float64 with changing inputs, masks, signed propagation, unequal spectral
sampling, attenuation, clipped support and reset. They pass on Julia 1.13.3,
CUDA 6.4.2 and an RTX 3080. They establish correctness for these software paths;
no RTC latency or deadline claim follows from them.

A focused allocation observation of public `step_graph!` on 32 × 32 grids with
16 × 16 pupils found 0 Julia heap bytes for each of 20 warmed CUDA capture
replays, for scalar and mapped two-channel trains in both precisions. The
scalar case used a lens and 0.14 m propagation; the mapped case used 1.0/2.1 µm
trains and 12 µm output cells. Measurement followed 32 warm steps and three
measurement warmups, included the public method's completion, and excluded
caller setup/synchronization. This result applies to those workloads and
versions; it is not a general allocation guarantee.

Run the optional CPU target with `julia --project=test/proper test/proper/runtests.jl`
after instantiating that environment. It checks missing-dependency rejection,
YAML loading, direct Proper agreement, physical sampling, padding, phase/flux
response and complementary mask power.

## Sequential ZEMAX preparation

`Optics.Zemax.inspect(path; configurations=[3, 5])` reads a bounded UTF-16LE,
MM sequential ZMX snapshot. Omit `configurations` to inspect the saved
configuration. It applies supported fixed multi-configuration values and linear
numeric pickups, preserving editor row order including comments. Thermal thickness
pickups are admitted only when a shared fixed environment and planar air endpoints
prove zero thermal change; other solves are rejected. System operands remain
metadata and do not turn automatic surface footprints into apertures. Inspection
does not run a ray trace. Configuration numbers belong to the authored recipe.

For a subsystem inside a larger simulation, use
`Optics.Zemax.prepare_component(root, recipe, branch; settings, diameter_m)`.
It returns native PROPER `config`, reduced-angle `diagnostics` matrices at the
selected surfaces, exact input `sources` hashes and explicit `assumptions`.
It checks the instrument's original ZMX and the recipe's AGF inputs at every cold
preparation. It writes no files and needs YAML, without requiring Proper or RTC.
The caller puts the configuration into an ordinary `proper_optical_train_f32` or
`proper_optical_train_f64` node alongside its other components.

Version-2 subsystem recipes have exactly these fields:

```yaml
schema_version: 2
catalogs:
  - path: catalogs/SCHOTT.AGF
    sha256: <exact SHA-256>
conversion:
  first_surface: 48
  initial_direction: 1
  powered_mirrors: [50, 52, 66, 70]
branches:
  common_relay:
    configuration_id: 1
    last_surface: 77
assumptions:
  - Explicit design provenance and approximation choices.
```

Paths are relative to the instrument root, outside `generated/`; catalogs are
ordered, with the first matching glass name selected. `catalog_archive` is
optional for this mode. Settings contain `resolution`, `pupil_resolution`,
`wavelength_um` and `temperature_c`; `diameter_m` is the incident beam diameter.
The subsystem must start in air/vacuum. It ends before the last surface's outgoing
thickness. Every selected curved mirror requires explicit admission. Multiple
mirrors use unfolded reduced-angle power `−2dn c`, where `d` is incoming axial
direction; reflection reverses `d`. Translation is `d Δz/n`. Manual circular
apertures are retained; transverse tilts/decenters, conic aberrations and other
aperture shapes are explicitly omitted. These are centered scalar reference
models, requiring optical validation before scientific performance claims.

The following standalone workflow remains useful when the complete authored
definition consists of optical paths rather than importing one subsystem.

With YAML loaded, `Optics.Zemax.prepare(root)` reads the authoritative
`instrument.yaml` and authored `simulation/definition.yaml`. The source definition
selects preparation instead of graphs:

```yaml
instrument: ../instrument.yaml
preparation:
  kind: zemax
  recipe: zemax.yaml
```

The recipe path is relative to the source definition and must remain inside the
instrument root, outside `generated/`. Its name is configurable. The definition
retains the ordinary shared simulation settings, defaults, schemas, profiles and
bindings. Preparation deep-copies those values into `generated/simulation.yaml`,
removes `preparation`, adds root-relative graph paths such as `generated/scc.yaml`
and `generated/llowfs.yaml`, and makes its `instrument` path relative to that
prepared definition. The ordinary shared runtime loader consumes the prepared
definition. Full source instrument validation belongs to the RTC tooling; AOS
checks the preparation contract without requiring RTC as a dependency.

`Optics.Zemax.load_model(path)` loads the strict version-1 conversion recipe
without loading Proper. Its exact top-level fields are `schema_version`,
`catalog_order`, `conversion`, `branches`, and `assumptions`. `conversion` contains
`first_surface`, `input_lens_surface`, `input_focal_parameter` (a native parameter
key string), `focal_mask_surface`, `oap_surface`, and `apodizer_surface`. Each named
`branches` entry contains exactly `configuration_id`, `detector_surface`, and
`lyot_surface`. Branch names match physical `optical_paths` names and use
lower-case filename stems. Each physical path references a sensor by its
independent ID and selects `lyot_stop_port: transmitted` or `reflected`. Physical
routing does not select ZEMAX configurations or surfaces. Apodizer, focal-mask and
Lyot surfaces must be distinct. Surface IDs are nonnegative integers,
configuration IDs are positive integers, and all selected surfaces must be active
within each branch's first-surface-to-detector interval.

The recipe is limited to 1 MiB, 1–16 branches, 1–64 distinct catalog filenames,
256 assumptions, and 4096 bytes per string. Unknown fields, numerical booleans,
duplicate YAML mapping keys, and multiple documents are rejected. The returned
mapping has bounded concrete value types; no recipe field evaluates Julia code.

Physical `optical_prescription` paths are relative to the instrument root. The
prescription contains exactly `format: zemax`, `path`, `sha256`, and
`catalog_archive: {path, sha256}`. Keep supplied ZMX and frozen archive bytes outside
`generated/`. Preparation inspects the original ZMX using the recipe's selected
configuration IDs and writes `generated/optics/resolved.yaml` automatically. With
`catalog_dir` pointing to caller-supplied extracted AGF catalogs, one call builds
all generated outputs from a fresh source tree. The catalog snapshot is always
`generated/optics/glasses.yaml`; without `catalog_dir`, preparation reuses this
frozen snapshot and fails with extraction guidance if it is absent or its catalog
order differs from the recipe. Catalog
metadata binds the declared archive hash and individual AGF hashes; it does not
prove that the supplied directory was extracted from that archive. Archive
extraction is external to AOS.

Preparation checks root containment, matching branch names, required active
surfaces, train conversion, output collisions and symlink aliases before writing
any output. All derived assets stay within a real `generated/` root. Authored
sources cannot declare both `graphs` and `preparation`, select unknown preparation
fields or kinds, or identify a different instrument.

`Optics.Zemax.train` consumes the caller's explicit branch selection and required
`lyot_mode` keyword (`transmit` or `reflect`); `prepare` maps the physical Lyot stop
port to that backend mode. The generated `generated/optics/provenance.yaml`
receipt hashes the exact authored and frozen snapshot source set, selected graphs,
and a distinct `execution: {path, sha256}` entry for `generated/simulation.yaml`.
`Optics.Zemax.validate_provenance(root)` reads frozen artifacts and checks those
hashes without regenerating files or invoking a formatter.

Run the instrument's `tools/prepare-proper.jl <instrument>` wrapper; supply
`--catalog-dir <extracted-catalogs>` when regenerating glass data. Standalone
`yamlfmt` formats generated documents; do not edit them directly. `inspect` and
`write_inspection` remain available for standalone inspection.

Conversion is centered and paraxial: STANDARD refraction, flat TILTSURF surfaces,
signed axial coordinate breaks, and explicitly selected input-lens and OAP
idealizations. Transverse tilts/decenters are omitted. Unsupported surface types
and unmodeled curved mirrors are rejected. The authored conversion rules and
idealized masks require instrument-specific optical validation.

## External pupil observations

`AlgorithmGraphs.prepare_external_pupil` binds a named tuple of ordinary prepared
CPU optical graphs to independent `IntensityExposure` products. This is a
transport-neutral scientific owner, not another graph executor. It uses the
same `step_graph!` and `reset_graph!` methods as other AOS callers.
`step_external_pupil!` accepts coherent OPD/amplitude arrays, `ModelTimestamp`
and `ModelDuration`, with an explicit source domain, generation and sequence.
Adapters own packet formats and convert their receipts to these values.

Each sample declares piecewise-constant coverage. Propagation precedes intensity
integration; a complete exposure mean is the duration-weighted sum over its
covered model interval. Missing coverage, sample loss or a newer generation
invalidates partial exposures; duplicates and overlaps are rejected. Exposure
origins currently share model time zero. CPU reference input buffers are
Float64, and no physical photometry or detector response is inferred.
Cross-quantity and input/output aliasing are rejected before execution.

An observation can additionally bind a prepared response graph, with explicit
input/output names, `normalized_to_photon_rate`, metric `pixel_sampling_m`,
`exposure_duration_s` and `frame_schema`. The scale maps each mean relative
intensity sample to cell-integrated photons/s and includes the caller's spatial
and flux calibration. It is independent of the current image sum. Input shape,
pixel sampling and exposure duration must agree; no resampling is inferred.
Preparation probes and resets the optical graph to establish its final sampling.
The ordinary response graph can contain a calibrated `pixel_mtf` and detector
acquisition node. It executes once after each complete optical exposure, so
detector noise is not sampled at the pupil-input cadence. Missing coverage
produces no response, and input discontinuities reset response state. Emitted
detector frames carry their declared schema and quantity metadata.

Instrument tooling selects authored response graphs under optional
`observations.<path>` entries, references existing detector components, and uses
the same graph loader and FITS construction bindings. ZEMAX preparation retains
the authored response file and calibration identities in its receipt and adds
their graph keys to the execution definition. Response graph files remain
authored inputs outside `generated/`.

`IntensityExposure` also supports direct contiguous integration. Reset it at a
gap before adding more data; overlaps and gaps never count as complete coverage.
Input intensity must be real, finite, nonnegative, representable as Float64, and
independent of the exposure's accumulation and output storage.
Its emitted image is borrowed until the next exposure update, so asynchronous
consumers must copy it.

## Recommended Ownership

Keep the packages loosely coupled.

- `AdaptiveOpticsSim.jl` should accept a complete PDM command at the plant
  boundary and update its physical state.
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
    resolution=128,
    output_rows=128,
    output_columns=128,
    diameter_m=8.0f0,
    wavelength_um=1.65f0,
)
node = proper_propagation_node(:science_propagation, configuration)
```

The maintained node has two complete-frame inputs, `pupil_opd` and
`pupil_amplitude`, and one output. Its configuration fixes the prescription,
square propagation resolution, rectangular output shape, physical diameter,
wavelength, random seed, port schemas, and numeric type before execution.
`output_rows` and `output_columns` default to `resolution`.
`prepare_algorithm_graph(...; target=...)` selects the exact CPU or GPU device
for the complete graph. The default output schema is normalized intensity. A
different schema may be declared only when the prescription and downstream
integration actually satisfy that contract; the adapter does not invent a
crop, resampling rule, photon-rate scale, or detector exposure.

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
the preparation hook receives both the square propagation resolution and the
declared output shape so it can prepare a target-resident mapping without
allocating during a graph step. The prescription owns the physical meaning of
that mapping; rectangular dimensions alone do not define detector sampling.
`reset_proper_assets!` resets any replaceable asset state at a serialized graph
boundary. The prescription remains an ordinary Julia function and receives the
prepared `pupil_opd`, `pupil_amplitude`, `diameter_m`, field, wavefront, output,
and run context as keyword arguments.

For a YAML-defined graph, the application explicitly admits the trusted Julia
prescription under a node-type name with `proper_propagation_node_factory` and
passes that factory in the loader's immutable `node_types` map. The YAML file
can select the admitted name and construction values, but it cannot resolve or
evaluate Julia code. Direct Julia graph construction remains available when a
SPIDERS prescription or topology is generated or conditional.

The companion contains separate LLOWFS and SCC topology examples. Applications
can express these topologies in the maintained YAML schema. An application
registers each node type
with `proper_propagation_node_factory(configuration)`, where the captured
configuration owns the trusted prescription and all physical values. The files
are separate because the simple executor invokes every node once per graph
step; the caller can therefore drive LLOWFS and SCC at different complete-frame
model-time boundaries. These files end at sensor observations. The primary HIL
loop publishes those observations through an external transport, receives one
complete PDM command from the external RTC, applies the measured DM response in
AOS, and then renders later path OPDs. An in-process RTC graph is useful for
reference tests, but it is not implicit in either SPIDERS sensor file.

One prepared Proper node is homogeneous in numeric type, array backend,
and exact device. A complete algorithm graph may connect different element
types through explicit conversion algorithms, but every array in one graph is
currently placed on one exact compute target. Neither the node nor the
graph inserts a fallback or transfer. A GPU-produced science frame therefore
reaches a CPU RTC or external transport client only through an explicit
application-owned device-to-host copy or prepared handoff after successful
graph publication.

Focused companion tests exercise square and rectangular outputs on CPU and
AMDGPU/CUDA hardware with GPU scalar indexing disabled. They establish
CPU/GPU numerical agreement, exact residency, graph buffer identity, explicit
host publication, and zero warmed Julia heap allocation on the CPU path. The
SPIDERS files are also loaded and stepped independently with test
prescriptions. These checks do not establish coronagraph fidelity, SPIDERS
end-to-end validity, detector response, observation labeling, service time,
acquisition cadence, or asynchronous failure isolation.

The node is usable today from Julia scripts and
`AdaptiveOpticsSim.AlgorithmGraphs`. Direct execution of Julia prescriptions in
PipeWire remains a separate Julia-host/module adapter problem. This graph-node
adapter does not turn Julia code into a native FGN shared-library plugin or
make a particular FGN implementation part of the AOS API.

## NCPA Ownership

Use the native `NCPA` or `OPDMap` model for a static or slowly varying
aberration that is adequately represented as sampled pupil OPD. Apply it only
to the direct-Julia or graph branch that contains the aberration, with exact
pupil-plane metadata and an application-owned run copy. Do not also mutate the
same path storage manually or apply the same OPD inside the external
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

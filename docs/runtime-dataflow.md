# Runtime Dataflow Guide

Status: active

## Purpose

This guide describes the maintained execution and ownership model. The
canonical HIL-neutral runtime is `AdaptiveOpticsSim.Plant`; application and
model-specific code composes the independent numerical optics and control
primitives around that boundary.

`AdaptiveOpticsSim.AlgorithmGraphs` is the smaller complete-frame composition
boundary. It schedules AOS-native graph-node adapters and deliberately does not
reproduce Plant's physical event calendar, detector lifecycle, or command
timing.

For package structure, see
[`maintainer-architecture.md`](maintainer-architecture.md). For normative HIL
requirements, see
[`hil-package-boundary.md`](hil-package-boundary.md) and
[`hil/compliance-matrix.md`](hil/compliance-matrix.md).

## Ownership Layers

The runtime distinguishes six roles:

1. Cold definitions describe stable identities, optical topology, acquisition
   topology, command schemas, timing policies, and model choices.
2. Run-immutable plans describe validated numerical or physical execution,
   including mappings, coefficients, compatibility, and backend requirements.
3. Persistent mutable state has a single writer and affects later scientific
   results. It includes the virtual-time scheduler, command and detector
   lifecycles, controllable-optic state, and product-readiness state.
4. Replaceable workspaces own scratch and execution resources. Recreating one
   cannot change the deterministic physical trajectory.
5. Products are explicit caller-visible values such as a `PupilFunction`,
   `IntensityMap`, `WFSObservation`, `WFSMeasurement`, or
   `AcquisitionProducts`; they are not scratch.
6. A prepared execution owner validates and binds the exact concrete plan,
   state, workspace, products, and backend/device/context for one writer.

Preparation may construct all six run-owned roles, but construction time does
not make every result a plan. Likewise, mutability alone does not distinguish
state from workspace: discarding state can change the next result, while a
workspace is replaceable without changing the trajectory. The package does
not retain direct Julia `Memory` use in the target architecture; issue #225
migrates the current representation. Prepared owners enforce final capacity
and armed-state non-growth over FixedSizeArrays.jl fixed-size arrays with
concrete element types, tuples, purpose-built owners, or backend arrays. A
resizable vector is limited to cold construction before sealing. Preparation
does not return a field or collection with stored type `Any`; external input is
normalized into concrete family ownership before execution.

The telescope owns revisioned aperture geometry, intensity reflectivity,
diameter, and spatial sampling. It owns neither a path's mutable OPD or electric
field nor detector cadence, exposure, FFT scratch, or atmosphere model time.

## Complete-Frame Algorithm Graphs

An algorithm graph has a cold definition and one exact prepared owner. Graph
preparation proceeds in three stages:

1. Call `graph_node_ports(Node, config)` for each node and admit its fixed port
   type, shape, scientific schema, role, and layout.
2. Admit the complete topology and startup sparse parameters, retain graph
   inputs, and allocate node outputs and delayed-link storage.
3. Call `prepare_graph_node` so each adapter constructs one exact execution
   owner from its config, initial props, admitted inputs, outputs, parameters,
   and compute target.

The final stage is necessary for AOS operations whose prepared execution owner
must bind exact plan, state, workspace, product, backend, device, and execution-
context identities. Domain algorithms keep their canonical APIs, such as
`Control.update!` and `Calibration.combine_basis!`; only the small graph-node
adapter knows about graph ports. No binding or allocation occurs in
`step_graph!`.

One graph step invokes every node once in validated declaration order. A direct
link exposes an earlier node's complete output to a later node in the same
step. A delayed link exposes only the value committed after the preceding
successful step. Delayed-link state and the graph sequence commit after all
nodes return normally. A node failure makes the run fail-stop until an explicit
reset; the graph does not claim graph-wide rollback of outputs already written
by earlier nodes.

This executor has one writer and no task, queue, wall-clock pacing, transport,
or dynamic placement policy. `prepare_algorithm_graph` accepts one exact compute
target, defaults to the host, allocates packed column-major node outputs and
delayed storage there, and rejects graph inputs or sparse parameters that are
not native packed storage on that same target. The prepared graph retains one
device execution context. Preparation, stepping, and reset run inside it and
complete its stream before successful publication or error return.

One graph may still connect port formats with different element types where a
node explicitly declares that conversion. Exact target ownership does not
impose one graph-wide numeric type. Host/device movement remains an explicit
application or prepared-handoff operation; the graph never inserts a fallback
or transfer. Runtime scalar-prop transactions and live sparse-parameter
replacement likewise remain open graph-host work rather than implied
capabilities.

Optional companion packages may bind a specialized prepared owner after graph
storage admission. The `AdaptiveOpticsProperHIL.jl` companion uses that seam to
bind one complete-frame Proper prescription directly to the admitted pupil-OPD,
pupil-amplitude, and output arrays while retaining adapter random state and
prescription scratch separately. Its focused CPU, AMDGPU, and CUDA tests
qualify that narrow exact-target execution path. They do not add an acquisition
clock, asynchronous WFS/science concurrency, or a PipeWire Julia host to this
executor.

### External-RTC Lockstep Boundary

`prepare_graph_hil_boundary` wraps one unstepped prepared graph without
changing its topology or scientific owners. It snapshots the initial graph
command for reset and binds two distinct host `Array` exchange buffers: a
caller-readable completed frame and a caller-writable command response. For an
accelerator graph, these buffers are the explicit GPU-to-CPU and CPU-to-GPU
boundary; each copy completes in the graph's retained device context before
the operation returns.

Once bound, the lockstep boundary is the graph's execution owner. Application
code must not call `step_graph!` on that graph or mutate the original exact
command-input array. Transport code writes only `hil_command_buffer`; the
boundary validates and adopts it at the serialized frame boundary.

The serialized application loop is:

```julia
boundary = prepare_graph_hil_boundary(
    graph;
    command_input=:pdm_command,
    frame_output=:shwfs_frame,
)

sequence = step_hil_frame!(boundary)
publish_frame!(transport, sequence, hil_frame_buffer(boundary))
receive_command!(
    transport,
    sequence,
    hil_command_buffer(boundary),
)
adopt_hil_command!(boundary, sequence)
```

The sequence-zero command already present in the graph produces frame one.
Command `n` is a response to frame `n` and becomes active only for frame
`n + 1`. The boundary rejects missing, stale, skipped, non-`UInt64`, and
nonfinite command responses before mutating the graph command. It stops after
an execution or target-copy failure because either graph products or the
active command may then be partial. `reset_hil_boundary!` restores the
snapshotted initial command and graph sequence. When a model-time driver is
used through `step_hil_frame_at!`, reset both owners with
`reset_hil_boundary!(boundary, driver)`.

The example transport functions above are application placeholders, not AOS
APIs. The boundary deliberately defines no socket, PipeWire buffer, wall-clock
deadline, timeout, last-command hold, queue, or concurrent callback. A direct
Julia script can perform the loop synchronously; a PipeWireAO adapter can map
the same buffers and sequences to its transport lifecycle.

## TOML Algorithm Graph Format

TOML is the maintained human-authored format because it supports comments,
keeps arrays of node and link tables readable, and is available as a Julia
standard library. JSON may be added later as an interchange representation,
but it is not a second maintained authoring schema. Direct Julia construction
remains the unrestricted API for generated topology, conditional composition,
multi-rate orchestration, and other arrangements that do not fit the static
file subset.

The key words **MUST**, **MUST NOT**, **REQUIRED**, **SHOULD**, **SHOULD NOT**,
and **MAY** in the following requirements are to be interpreted as described by
[BCP 14](https://www.rfc-editor.org/info/bcp14) when, and only when, they appear
in all capitals.

- **AOS-GRAPH-FILE-001:** A maintained graph file MUST set
  `schema_version = 1`, MUST have a stable identifier-valued `name`, and MUST
  contain at least one `[[nodes]]` table. A loader MUST reject unknown schema
  versions and unknown fields.
- **AOS-GRAPH-FILE-002:** Node order MUST be execution order. Each node MUST
  identify an explicitly supplied `type`; the loader MUST NOT evaluate Julia
  source, resolve arbitrary module paths, or mutate a global type registry.
- **AOS-GRAPH-FILE-003:** `[nodes.config]` MUST contain construction and
  graph-rebuild values. `[nodes.props]` MAY contain scalar initial props but
  MUST NOT change a node's admitted ports. Version 1 does not define live prop
  transactions.
- **AOS-GRAPH-FILE-004:** Frame inputs, delayed-link initial arrays, and sparse
  parameter arrays MUST be referenced by external binding name. The core
  loader MUST NOT infer FITS or another calibration file format and MUST NOT
  insert a host/device transfer.
- **AOS-GRAPH-FILE-005:** Endpoints MUST use `node.port` syntax. Direct links,
  delayed links, boundary inputs and outputs, and sparse parameters MUST pass
  the same type, shape, schema, layout, ownership, and connectivity admission
  checks as direct Julia construction.
- **AOS-GRAPH-FILE-006:** `load_algorithm_graph` MUST compile a valid file into
  the same concrete `AlgorithmGraphDefinition` used by `algorithm_graph`. The
  parsed TOML dictionaries and cold builder vectors MUST NOT survive into the
  definition or prepared execution owner.
- **AOS-GRAPH-FILE-007:** Version 1 is a static, single-rate, complete-frame
  subset. Multi-rate scheduling, conditional topology, generated nodes, and
  physical sub-frame event semantics SHOULD use direct Julia composition or
  Plant rather than acquire implicit file semantics.

The version 1 tables are:

| Table | Required fields | Meaning |
|---|---|---|
| root | `schema_version`, `name`, `nodes` | Schema and graph identity |
| `[[nodes]]` | `name`, `type`, `[nodes.config]` | One node in execution order; `[nodes.props]` is optional |
| `[[inputs]]` | `name`, `destination`, `binding` | Caller-owned frame input |
| `[[outputs]]` | `name`, `source` | Caller-visible allocated graph output |
| `[[links]]` | `source`, `destination` | Same-step direct link |
| `[[delayed_links]]` | `source`, `destination`, `initial` | One-successful-step delayed link; `initial` is a binding name |
| `[[parameters]]` | `destination`, `binding` | Required startup ndarray sparse parameter |

For example:

```toml
schema_version = 1
name = "revolt_classic"

[[nodes]]
name = "controller"
type = "discrete_integrator_f32"

[nodes.config]
extent = 277
sample_period_s = 0.001
input_schema = "org.revolt.residual-modes/1"
output_schema = "org.revolt.command-modes/1"

[nodes.props]
gain = 0.3
tau_s = 0.02

[[inputs]]
name = "residual"
destination = "controller.input"
binding = "residual"

[[outputs]]
name = "command"
source = "controller.output"
```

The application supplies the array explicitly and then prepares the ordinary
graph:

```julia
definition = load_algorithm_graph(
    "revolt_classic.toml";
    bindings=(residual=residual_modes,),
)
graph = prepare_algorithm_graph(definition; target=target)
```

The built-in type map currently contains `ccd_detector_acquisition_f32`,
`closed_loop_correction_f32`, `control_matrix_reconstruction_f32`,
`deformable_mirror_surface_f32`,
`discrete_integrator_f32`, `emccd_detector_acquisition_f32`,
`gaussian_deformable_mirror_surface_f32`,
`modal_opd_expansion_f32`,
`multilayer_atmosphere_opd_f32`,
`pupil_opd_composition_f32`,
`pyramid_rate_f32`,
`shack_hartmann_centroid_f32`, `shack_hartmann_rate_f32`, and
`shack_hartmann_slope_selection_f32`.
`merge(builtin_graph_node_types(), companion_types)` creates an explicit larger
map for optional packages such as the Proper companion.

The maintained
[`revolt_classic_hil.toml`](../examples/graphs/revolt_classic_hil.toml) file
defines the intended external-RTC sensor boundary. Its
`multilayer_atmosphere_opd_f32` node owns and advances the maintained REVOLT
five-layer atmosphere by exactly 1 ms per graph step. A serialized adapter
places the latest accepted complete PDM command in the caller-owned input
buffer before that step. The graph applies the command with
`gaussian_deformable_mirror_surface_f32`, composes the resulting sampled
surface with the atmospheric OPD through `pupil_opd_composition_f32`, and
produces the configured 352-by-352 diffractive SHWFS photon-rate mosaic and
one completed noisy CCD frame for a transport adapter such as PipeWireAO. The
RTC performs centroiding, reconstruction, control, and command production
outside this graph.

The on-sky and Copper trees contain the same 277-actuator HSDM277 map payload;
their raw files differ in header and have hashes
`b75ae701627b24af51ecd50f019e91989a514441a0d510b418b5fe82b6127fd1`
and `76b60effb1786a6cdb37ef3b51c12e34cdc592803f1b1161895b69ee85d51ecc`.
[`revolt_hsdm277.jl`](../examples/support/revolt_hsdm277.jl) reproduces that
exact command order. It places the 19-actuator physical grid at the REVOLT
simulation assumption of 17 actuator centres across the pupil. The maintained
OOPAO ideal profile's 0.35 adjacent coupling then gives a normalized Gaussian
influence width of `0.08626550214129701`; an older OOMAO script uses 0.30.
Command-map topology is source-backed, while pupil registration and coupling
are explicitly provisional model assumptions. The graph command contains
unit-peak actuator surface-OPD coefficients in metres, so an ALPAO normalized
hardware command still requires separately qualified conversion.

No authoritative HSDM277 influence-function artifact has been identified. The
locally available `AX307_Influences.fits` instead contains 468 BAX307
influences over 5,813 pupil samples and belongs to the SPIDERS physical mirror.
When qualified HSDM277 sampled influences become available, the analytic node
can be replaced by `deformable_mirror_surface_f32` without changing the
downstream surface-OPD contract.

The maintained
[`revolt_copper_hil.toml`](../examples/graphs/revolt_copper_hil.toml) file
defines the corresponding external-RTC boundary for Copper. It combines the
maintained REVOLT optical and iXon camera profile with the current RTC's
64-by-64 frame and 30-by-30 Pyramid-pupil contract. It applies the same exact
command geometry and provisional surface model as Classic. Its own atmosphere
node advances the same maintained five-layer profile by exactly 2 ms per graph
step. That step composes the surface with atmospheric OPD, forms the modulated
Pyramid photon-rate image, and completes one seeded EMCCD acquisition; the
external RTC consumes that frame and returns the next physical-DM command.
The file records its layered parameter authority and does not claim that the
symmetric simulated pupil registration already matches the measured
`pwfsRoiOffsets_64.csv` positions or an operational pixel reconstructor.

The REVOLT graph profile is selected before preparation. The example helper keeps
the selection explicit:

```julia
include("examples/support/revolt_hil_graphs.jl")

path = REVOLTHILGraphs.graph_path(:copper, :fast_dm)
definition = load_algorithm_graph(path; bindings)
graph = prepare_algorithm_graph(definition; target)
```

`full_optical` resolves to the primary files above.
[`revolt_classic_hil_fast_dm.toml`](../examples/graphs/revolt_classic_hil_fast_dm.toml)
and
[`revolt_copper_hil_fast_dm.toml`](../examples/graphs/revolt_copper_hil_fast_dm.toml)
change only the provisional HSDM277 implementation. They scatter the exact
277-element command into its regular 19-by-19 physical grid and evaluate the
same Gaussian influence model with a separable matrix factorization. Internal
pupil sampling, atmosphere evolution, WFS propagation and modulation,
detector noise, external frame shapes, and lockstep HIL semantics are
unchanged. This is an exact specialization under the current provisional
regular-grid registration, not a reduced-resolution optical model. A prepared
run never changes profile in response to overload.

The separate
[`revolt_classic_rtc_reference.toml`](../examples/graphs/revolt_classic_rtc_reference.toml)
file extends the same sensor stages into an optional in-process RTC reference.
A calibrated centroid node consumes the completed frame using caller-bound
valid-subaperture and reference-signal arrays and publishes the canonical full
512-component slope vector. A selection node converts that vector into 376
pair-interleaved components using a caller-bound one-based ordering for 188
lenslets. The core graph does not parse or silently reinterpret the
instrument's detector ROI-origin file. A control-matrix node snapshots and
applies the caller-bound 221-by-376 calibrated matrix, then an atomic
closed-loop correction node applies the REVOLT Classic gain `-0.3`, pole
`0.99`, and anti-windup gain `1.0`. Its OOPAO-compatible fractional cutoff is
not mapped from SPECULA's absolute-pedestal `thr_value`; numerical extractor
parity and operational matrix compatibility remain open. Until the optional
reference graph includes downstream constraint and DM nodes, the caller binds
the preceding demanded-minus-realized correction feedback and the first step
ignores that placeholder. The reference graph is a deterministic standalone
and differential-testing aid, not the required HIL deployment topology.

Architecture-file status is therefore explicit:

| Architecture | File-defined executable surface | Remaining authority or implementation gate |
|---|---|---|
| REVOLT Classic | The primary HIL file advances the maintained five-layer atmosphere at 1 ms, applies the exact HSDM277 command-map topology with provisional normalized-pupil placement and Gaussian surface response, composes atmospheric and surface OPD, and executes diffractive SHWFS optics plus single-read CCD acquisition through the completed frame boundary. The selectable fast-DM file preserves the same 240-sample optics, noisy 352-by-352 frame, and provisional Gaussian surface through a regular-grid separable evaluation. A separate RTC-reference file additionally executes AOS/OOPAO centroiding, explicit 188-lenslet/376-component slope selection, caller-bound 221-by-376 reconstruction, and atomic 221-coordinate correction. | Bind a qualified normalized-hardware-command conversion and measured HSDM277 influence model before claiming an instrument model; add any required non-atmospheric path aberrations explicitly. For the optional RTC reference, bind the authoritative ROI order/matrix, qualify SPECULA extraction parity, and close downstream command constraints and feedback. |
| REVOLT Copper | The primary HIL file advances the maintained five-layer atmosphere at 2 ms, applies the same exact command map and provisional placement/Gaussian response, composes atmospheric and surface OPD, and executes maintained 32-point modulated Pyramid optics plus complete-frame noisy EMCCD acquisition through a 64-by-64 frame boundary. The selectable fast-DM file changes only the Gaussian surface evaluation and retains the 480-sample pupil and complete Pyramid/detector model. | Bind qualified command conversion and measured HSDM277 influences. Bind measured Pyramid registration and an operationally compatible pixel reconstructor before claiming RTC numerical parity; add any required non-atmospheric path aberrations explicitly. |
| SPIDERS | Optional atomic Proper propagation node in the companion; no complete architecture file | Select the maintained science/control topology and qualify its native or Proper optical nodes before fixing a static profile |

End-to-end REVOLT Classic, REVOLT Copper, and SPIDERS simulations are not
claimed until all nodes needed by those architectures are executable and
scientifically qualified. This prevents a configuration file from becoming a
catalog of plausible but non-running node names. Direct Julia construction
remains the correct extension path in the meantime: an application may
generate nodes, branch on configuration, admit a trusted Proper prescription
factory, or compose several prepared graphs without expanding the static TOML
language.

## Definition And Preparation

A plant starts from `PlantDefinition`. It may declare:

- one telescope and atmosphere
- named optical paths
- independent acquisition endpoints
- independent controllable optics
- native sampled `NCPA` or `OPDMap` aberrations with explicit path visibility
- one or more command schemas per optic

`prepare_plant` requires an explicit run seed and, for every command endpoint,
a `CommandEndpointConfiguration`. Preparation:

- canonicalizes stable identities
- makes run-owned exact-target sampled-OPD copies from supported caller input
- resolves bounded path-local sampled-aberration and controllable-optic
  couplings
- prepares path and acquisition implementations
- retains one exact prepared device execution context in the plant, its path
  executors, and their acquisition owners
- allocates backend- and device-compatible destinations
- copies initial and optional safe commands into endpoint-owned storage
- constructs independent RNG streams from stable owner identities
- validates all run-immutable bindings

Successful public preparation and execution returns complete the retained
backend stream before restoring the caller's device and stream selection.
Worker-visible independent path-group calls establish that completion before
publishing their ready or complete phase; coordinator-local serial execution
uses one enclosing completion boundary.

Preparation may allocate and fail. Repeated execution mutates only prepared
state and caller-owned products. Public schedule-free
`materialize_path_input!`, `execute_path!`, and `execute_acquisition!` calls
enter the retained context, including its accelerator stream where applicable,
and restore the caller's previous device and stream selections after success or
an exception. The schedule-free `execute_acquisition_selection!` and
`execute_acquisition_selection_at!` entry points enter the plant context once
and execute their lower-level operations inside it.
Model-specific execution dispatches run inside that boundary and must not
select an alternate device or stream.

## Controller Output Routing

An RTC or controller often computes one flat vector even though the plant owns
several independently timed command endpoints. The maintained bridge is
prepared controller-output routing:

```julia
using AdaptiveOpticsSim
using AdaptiveOpticsSim.Plant

controller_buffer = zeros(Float32, 5)
products = (
    woofer=@view(controller_buffer[1:2]),
    tweeter=@view(controller_buffer[3:5]),
)

routing = prepare_controller_output_routing(
    plant,
    products,
    ControllerOutputRoute(:woofer, :woofer_command),
    ControllerOutputRoute(:tweeter, :tweeter_command),
)
```

Each named product is borrowed without packing or copying. Preparation requires
an exact match with its prepared endpoint's numeric type, shape, array backend,
and compute device. Scalar products use assigned `Ref` storage.

Routing owns no sequence number, effective timestamp, admission, transaction,
queue, transport, or optical-grouping semantics. Integration code obtains a
payload with `controller_output_payload` and constructs a separate
`PlantCommand` for each due endpoint. Plant-command admission then copies the
presented payload into that endpoint's bounded storage. Co-conjugated optics,
flat-buffer adjacency, and equal timestamps do not imply atomic application;
explicit `PlantCommandTransaction` membership does.

## Virtual-Time Event Flow

`prepare_plant_event_loop` composes a finite prepared plant with:

- command endpoint generators
- periodic optical samples
- periodic or delivered-trigger acquisition starts
- detector lifecycle transitions
- one shared atmosphere timeline
- sampled-aberration application
- controllable-optic application

`step_plant_events!` processes all work at the next canonical
`PlantTimestamp`; `run_plant_events_until!` advances through all due timestamps
up to a limit. Equal-time ordering follows stable causal phase, prepared owner
ordinal, and occurrence—not task completion or tuple iteration.
`PreparedPlantEventLoop` retains the exact prepared plant device execution
context. Both public execution calls enter that context and restore the caller's
previous device and stream selections; their lower-level optical work reuses
the already active context.

At a due command timestamp the loop:

1. claims the next application-ready command
2. applies endpoint value policy and stages physical-optic state
3. commits the visible physical state
4. records a bounded disposition

Commands due at the same timestamp are applied before optical samples. Each
endpoint retains independent sequence, effective-time, silence, and capacity
policy. The loop owns plant time only; wall-clock pacing, external timestamp
mapping, payload leases, transport, and RTC protocol belong to the HIL
integration layer.

For every due full-optical path, the serial order is atmosphere
materialization, sampled aberrations, controllable surfaces, autonomous
path-local optics, and then the typed path executor. Replacement sampled
aberrations precede canonical placement/identity-ordered additives.
Preparation rejects more than one replacement visible on a path.

## Optical Sample Flow

For each due full-optical path, execution:

1. advances the shared atmosphere to explicit model time
2. renders the published epoch through a path-local prepared direction
   renderer into caller-owned path storage
3. applies the visible static and controllable surfaces
4. executes the prepared optical model into its exact path result
5. executes each due acquisition provider into its own
   `AcquisitionProducts`

Preparation resolves required optic placement and all-path or selected-path
visibility into canonical bounded ranges and compatible path-coupling groups.
Full-optical execution applies only a path's visible surfaces additively.
Sampled surface models may bind the same-grid fast path or a finite-support
NGS/LGS atmospheric-conjugate transform; unsupported model/placement
combinations fail during preparation. No packed or composite optic object
couples device cadence.

Native direct science uses `prepare_direct_imaging` and
`form_direct_image!`. Prepared PROPER integration meets native optics at the
same declared photon-arrival-rate/acquisition boundary. NGS and LGS WFS paths,
science cameras, calibration paths, NCPA branches, and coronagraph paths remain
separate declared paths when their propagation differs.

## WFS And Detector Flow

The maintained WFS stage contract is:

1. `form_wfs_optical_products!` produces one or more detector-facing
   photon-arrival-rate products.
2. `acquire_wfs_observation!` applies detector acquisition into a typed
   `WFSObservation`.
3. `estimate_wfs_measurement!` writes a typed `WFSMeasurement`.

The optical front end, detector acquisition, and estimator are independent
components. For example, a Shack–Hartmann optical front end composes a
`MicrolensArray`; it does not make the microlens optic part of detector or
centroid state.

Prepared frame acquisition applies operations in this order:

1. validate the photon-arrival-rate product and its spatial measure
2. apply presampling detector response on the optical grid
3. integrate represented cells into physical pixels
4. apply wavelength-dependent QE and explicit exposure duration exactly once
5. apply charge-domain effects, stochastic response, binning, and readout
6. publish a complete product at the detector lifecycle's readiness event

Global shutter, rolling shutter, frame transfer, nondestructive reads, and
up-the-ramp sampling retain separate virtual-time lifecycle definitions.
Exposure, readout completion, product readiness, and re-arming are distinct
events.

## Acquisition Fidelity

Every prepared acquisition binds one run-immutable provider:

- full optical
- command-responsive reduced order
- synthetic or bounded replay

All provider styles write the same prepared logical product contract. This
allows one plant to combine high-fidelity paths with fast causal surrogates for
RTC throughput and latency testing. Fidelity does not change during a run; a
different provider requires another preparation.

Calibration illumination uses the ordinary path and acquisition machinery.
The user or companion model declares where it enters, which paths see it, and
how its state evolves. Core does not infer a special calibration mode from a
label.

## Explicit Model Loops

Package examples such as the Subaru AO188/AO3k model may expose their own
`prepare!`, `step!`, and `readout` functions. These are explicit model
compositions used for numerical examples and benchmarks; they are not a second
generic runtime API. A simple closed loop can likewise compose:

- `advance_by!` and `render_atmosphere!`
- `set_command!` and `apply!`
- the three WFS stages
- `AdaptiveOpticsSim.Control.reconstruct!`
- `AdaptiveOpticsSim.Control.VectorDelayLine` and
  `AdaptiveOpticsSim.Control.DiscreteIntegratorController`

This keeps reusable numerical/control primitives independent of HIL scheduling.

## Backend And Parallel Execution

Prepared storage is parameterized by numeric type, backend, and compute
device. Optical and detector kernels reject implicit host/device mixing and
avoid GPU scalar indexing. Same-owner stages call directly; future ownership
boundaries may exchange bounded descriptors.

`SimulationEnsemble` remains a coarse-grained facility for independent model
runs, sources, time-series sweeps, or offline studies. It is not the
deadline-path HIL scheduler. Sequential execution is the default; Julia
threads, AcceleratedKernels, Dagger, and backend streams are explicit policies
whose use must avoid nested oversubscription.

## Validation And Performance Evidence

Functional support is recorded in
[`model-validity-matrix.md`](model-validity-matrix.md). Production-surface
qualification is recorded in
[`supported-production-surfaces.md`](supported-production-surfaces.md).

Benchmark workloads under `benchmarks/` own their explicit loop composition;
they do not introduce a package runtime abstraction. Representative performance
artifacts live under `benchmarks/results/` and remain outside ordinary
`Pkg.test()`.

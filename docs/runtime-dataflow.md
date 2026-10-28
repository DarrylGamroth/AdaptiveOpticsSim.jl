# Runtime Dataflow Guide

Status: active

## Purpose

This guide describes the maintained execution and ownership model.
`AdaptiveOpticsSim.AlgorithmGraphs` is the complete-frame composition boundary.
It schedules AOS-native graph-node adapters while direct Julia code remains the
unrestricted scripting interface.

For package structure, see
[`maintainer-architecture.md`](maintainer-architecture.md).

## Ownership Layers

The runtime distinguishes six roles:

1. Cold definitions describe stable identities, optical topology, acquisition
   topology, command schemas, timing policies, and model choices.
2. Run-immutable plans describe validated numerical or physical execution,
   including mappings, coefficients, compatibility, and backend requirements.
3. Persistent mutable state has a single writer and affects later scientific
   results. It includes atmosphere evolution, detector persistence, graph
   delays, and explicit RNG state. RTC controller state belongs to its RTC
   owner.
4. Replaceable workspaces own scratch and execution resources. Recreating one
   cannot change the deterministic physical trajectory.
5. Products are explicit caller-visible values such as a `PupilFunction`,
   `IntensityMap`, `WFSObservation`, `WFSMeasurement`, or
   detector frames, or graph outputs; they are not scratch.
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
2. Admit the complete topology and startup sparse parameters, derive a fixed
   topological order from direct links, retain graph inputs, and allocate node
   outputs and delayed-link storage.
3. Call `prepare_graph_node` so each adapter constructs one exact execution
   owner from its config, initial props, admitted inputs, outputs, parameters,
   and compute target.

The final stage is necessary for AOS operations whose prepared execution owner
must bind exact plan, state, workspace, product, backend, device, and execution-
context identities. Domain algorithms keep their canonical APIs, such as
`Calibration.combine_basis!`; only the small graph-node adapter knows about
graph ports. RTC algorithms are owned by FilterGraphAlgorithms/JuliaFilterGraph
and are composed in the maintained
[`filter_graph_algorithms` fixture](../examples/integrations/filter_graph_algorithms/).
No binding or allocation occurs in `step_graph!`.

One graph step invokes every node once in the prepared topological order.
Preparation repeatedly selects the earliest declared node whose direct-link
predecessors have been ordered. Disconnected and simultaneously ready nodes
therefore retain declaration order. Direct cycles, including self-links, are
rejected; delayed links do not constrain same-frame order. Node names and exact
port bindings retain their identities when preparation changes execution order.
A direct link exposes its producer's ordered output to its consumer in the same
step. The default uses one retained accelerator stream,
so that dependency is device execution order and requires no host barrier. The
explicit grouped-stream policy may overlap independent nodes, but it admits
only groups whose flattened names equal declaration order and requires every
direct link to cross from an earlier group to a later group. Preparation
validates explicit groups against the original declaration and does not reorder
or repartition them. A delayed link exposes only the value
committed after the preceding successful step. Delayed-link state and the graph
sequence publish only when the current frame's completion is observed. A node
or device failure makes the run fail-stop until an explicit reset; the graph
does not claim graph-wide rollback of outputs already written by earlier nodes.

This executor has one writer and no task, queue, wall-clock pacing, transport,
or dynamic placement policy. `prepare_algorithm_graph` accepts one exact compute
target, defaults to the host, allocates packed column-major node outputs and
delayed storage there, and rejects graph inputs or sparse parameters that are
not native packed storage on that same target. The default and captured
policies retain one device execution context. Grouped execution retains a fixed
number equal to its widest group. `step_graph_async!` submits one frame and
returns a single-use `GraphStepTicket`; `wait_graph_step!` completes the joined
execution, publishes the sequence, and releases all mutable graph storage. The queue
capacity is exactly one because the next frame would otherwise reuse node
state, workspaces, outputs, delayed-link storage, and RNG state that may still
be in flight. Before consuming the ticket, the caller must not read graph
outputs, mutate graph inputs, reset the graph, or submit another frame. The
last three operations are rejected where the graph API can enforce them.

`step_graph!` is the synchronous compatibility surface and is equivalent to an
immediate submission followed by its completion wait. Some admitted nodes have
a necessary host or backend completion boundary and may block during
submission; the asynchronous API removes the graph-boundary barrier but does
not falsely classify every node as nonblocking. Multiple outstanding frames
would require separately owned bounded slots and explicit publication, not an
unbounded task or command queue.

Preparation defaults to `StreamGraphExecution()`. A wide MCAO or MOAO graph may
instead select an explicit bounded schedule:

```julia
execution = GroupedStreamGraphExecution(
    (:guide_star_1, :guide_star_2),
    (:wfs_1, :wfs_2),
    (:tomography,),
)
graph = prepare_algorithm_graph(definition; target, execution)
```

Each argument is one nonempty group. On an accelerator, the group width fixes
the retained stream and event count during preparation. Every lane waits on all
events from the preceding group, so this is a deliberately coarse barrier
rather than automatic edge-level scheduling. The primary stream joins the
final group before committing delayed links; ordinary ticket completion on
that stream therefore observes every lane. A failed partial submission drains
all lanes before storage can be reset. The host executes the same grouping
serially and creates no Julia worker tasks.

An accelerator application may separately select `CapturedGraphExecution()`.
That policy records each complete graph step as one native CUDA Graph or HIP
Graph executable. Every node owner must qualify, and the recorded operation
contains the ordered node sequence followed by delayed-link commits. There is
no direct-stream fallback inside a captured graph, and captured and grouped
execution are distinct policies. Dynamic values remain in fixed device
buffers, while evolving state needed during replay must also remain
device-resident.

`graph_node_capture_capability` defaults to unsupported. An adapter may return
`GraphNodeCaptureSafe()` only when its enqueue path has fixed array identities,
keeps all evolving replay state on the device, does not mutate scientific host
state, and records no dynamic execution-storage allocation, synchronization,
device-result query, or GPU host callback. Julia compilation, backend
recording, and native graph instantiation are cold preparation work; replay
launches the retained native executable. AMDGPU therefore uses its blocking,
allocation-light completion wait only for captured execution; ordinary stream
execution retains HostCall-compatible synchronization. The maintained built-in
qualification covers finite multilayer atmosphere evolution, coordinate and
regular-grid Gaussian DM evaluation, pupil-OPD composition, diffractive
Shack-Hartmann rate formation, both maintained Pyramid modulation strategies,
and qualified complete-frame CCD, CMOS, and EMCCD acquisition.
Coordinate-Gaussian capture currently requires the built-in linear-static
actuator response. Captured atmosphere motion and random draws advance
device-resident state during replay and publish matching host state only after
successful completion. `captured_graph_node_count` reports the number of nodes
in the captured step; requesting capture when any node is unsupported fails
during preparation rather than silently changing policy. Downstream instrument
packages qualify their complete graphs independently from these node-family
claims.

Dense and measured influence-matrix DM nodes also qualify for capture with the
built-in linear-static, clipping, and actuator-health stages or their static
compositions. Preparation retains the multiplication coefficients for replay.
Stateful or custom actuator stages need separate qualification. Complete
instrument capture still requires a downstream whole-graph check.

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
commands for reset and binds distinct host `Array` exchange buffers: a
caller-readable completed frame and caller-writable command responses. For an
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
    command_input=:dm_command,
    frame_output=:wfs_frame,
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

For a graph with several command inputs, pass a tuple of 1 to 16 unique port
names, such as `command_input=(:dm_command, :tiptilt_command, :offload_command)`.
`hil_command_buffer(boundary)` then returns a `NamedTuple` with these names and
one host `Array` per input. Commands can have different shapes and floating-point
types. To supply the arrays, pass `command_buffer=(dm_command=dm_buffer,
tiptilt_command=tiptilt_buffer, offload_command=offload_buffer)` with the same
names and order. The single-`Symbol` form continues to return one `Array`.
All graph command inputs, the graph frame, and the exchange arrays must have
distinct storage. Every command is checked for finite values before any graph
input changes; one `adopt_hil_command!` call adopts the complete response under
the outstanding frame sequence. The next frame waits for that complete response.

The sequence-zero command already present in the graph produces frame one.
Command `n` is a response to frame `n` and becomes active only for frame
`n + 1`. The boundary rejects missing, stale, skipped, non-`UInt64`, and
nonfinite command responses before mutating the graph command. It stops after
an execution or target-copy failure because either graph products or the
active command may then be partial. `reset_hil_boundary!` restores the
snapshotted initial command and graph sequence. When a model-time driver is
used through `step_hil_frame_at!`, reset both owners with
`reset_hil_boundary!(boundary, driver)`.

### Exposures under a held calibration probe

Use `PreparedGraphCalibrationBoundary` when several complete exposures must use
one adopted command:

```julia
boundary = prepare_graph_calibration_boundary(
    graph;
    command_input=:dm_command,
    frame_output=:wfs_frame,
)
receive_probe!(transport, hil_command_buffer(boundary))
adopt_hil_probe!(boundary, UInt64(7))

first = step_hil_exposure_at!(boundary, driver)
publish_exposure!(transport, first, hil_frame_buffer(boundary))
# Finish consuming the first frame before requesting another exposure.
second = step_hil_exposure_at!(boundary, driver)
publish_exposure!(transport, second, hil_frame_buffer(boundary))
```

Probe identifiers are positive, strictly advancing `UInt64` values independent
of graph exposure sequences. Adoption validates the complete command and
finishes its host-to-target copy before returning. Each exposure advances the
graph and optional model-time driver exactly once. Writing the host command
buffer alone does not change the held probe; transport completion, settling,
exposure duration, and response aggregation belong to the caller. Reset restores
the initial command and both sequences, and requires another explicit probe
adoption. Execution or copy failure stops the boundary until reset.

The complete-frame InGaAs graph node consumes a cell-integrated photon-rate
frame on an explicit metric pixel grid and declares an integrated spectral
channel. Its acquisition model is separate from the optional upstream pixel-MTF
node. A calibrated frequency response must use the same pixel-coordinate
convention as the image sampling. Non-null persistence or nonlinearity remains
on stream execution; InGaAs capture is limited to its admitted full-frame
null-model path. Neither node supplies camera calibration.

### Fully declared HIL reference systems

The maintained SHWFS and Pyramid HIL reference systems validate external-RTC
integration without depending on unmeasured instrument calibration. They share
one explicit 25-command analytic Gaussian deformable mirror, accept a complete
caller-owned uncompensated pupil OPD in metres, and publish a complete detector
frame through `PreparedGraphHILBoundary`. Their parameters are AOS validation
choices; neither system is an instrument profile.

[`shack_hartmann_hil_reference.yaml`](../examples/graphs/shack_hartmann_hil_reference.yaml)
publishes a noiseless 64-by-64 CCD frame. [`pyramid_hil_reference.yaml`](../examples/graphs/pyramid_hil_reference.yaml)
publishes a noiseless 36-by-36 four-pupil EMCCD frame. AOS built-in graph
execution stops at these detector frames; use the maintained
[`filter_graph_algorithms` fixture](../examples/integrations/filter_graph_algorithms/)
for the Shack–Hartmann and Pyramid FGA/JFG estimator and RTC compositions. The helper in
[`hil_reference_systems.jl`](../examples/support/hil_reference_systems.jl)
owns the exact analytic actuator coordinates, SHWFS valid-lenslet rule, graph
selection, and preparation recipe.

The current integration matrix is deliberately small:

| WFS | Detector response | Controller boundary | Evidence | Status |
|---|---|---|---|---|
| Shack–Hartmann | Deterministic, noiseless CCD | FGA/JFG lockstep reference controller | [`filter_graph_algorithms`](../examples/integrations/filter_graph_algorithms/) interaction-matrix and convergence tests | covered |
| Pyramid | Deterministic, noiseless EMCCD | FGA/JFG complete-frame S4 plant-to-command path | `algorithm-graphs` complete-frame tests plus the [`filter_graph_algorithms`](../examples/integrations/filter_graph_algorithms/) S4 boundary tests | covered |
| Shack–Hartmann | Deterministic, noiseless CCD | In-process pyRTC oracle | Measured interaction and closed-loop convergence | covered |
| Pyramid | Deterministic, noiseless EMCCD | In-process pyRTC oracle | Measured interaction and closed-loop convergence | covered |
| Shack–Hartmann | Deterministic, noiseless CCD | Native Julia SHM to a pyRTC process | Bidirectional protocol checks, measured interaction, static convergence, evolving-atmosphere PSFs, and on-axis Strehl improvement | covered in lockstep |
| Pyramid | Deterministic, noiseless EMCCD | Native Julia SHM to a pyRTC process | Bidirectional protocol checks, measured interaction, static convergence, evolving-atmosphere PSFs, and on-axis Strehl improvement | covered in lockstep |
| Either | Stochastic detector response | Any external RTC | Statistical closed-loop qualification | gap |

Calibrate an external RTC by applying complete positive and negative commands
through the reference boundary. This measures the interaction matrix in the
same command units and WFS-signal order used by that RTC. Do not substitute an
identity matrix or infer a command scale from an instrument name.

The initial pyRTC integration runs from Julia through PythonCall.jl because AOS
owns graph stepping, frame sequence, and model time. PyJulia/JuliaCall would put
Python in charge of that lifecycle and is not the initial integration
direction. The maintained SHWFS and Pyramid runners share
[`pyrtc_reference_hil.jl`](../examples/integrations/pyrtc/pyrtc_reference_hil.jl),
which uses pyRTC's real `ImageSHM`, `SlopesProcess`, `WavefrontCorrector`, and
`Loop` implementations synchronously. The SHWFS permutation publishes 104
signals from its 64-by-64 lenslet mosaic. The Pyramid permutation publishes 344
signals from four 16-pixel pupils in its 36-by-36 frame. Each permutation
measures its own interaction matrix by push-pull commands, installs it in
pyRTC, injects a disturbance that lies in the declared deformable-mirror span,
and requires the pyRTC command to close the AOS optical loop.

The scientist-facing installation, validation, and execution workflow is in
[`user-guide.md`](user-guide.md#use-aos-with-pyrtc). The in-process oracle uses
the GitHub-installed pyRTC package from the Python interpreter selected by
`JULIA_PYTHONCALL_EXE`.

The live demonstration also publishes `wfc2D`, `aosUncompensatedOpd`,
`aosDmSurfaceOpd`, `aosResidualOpd`, `aosOpenLoopPsf`, and
`aosClosedLoopPsf` as diagnostic-only streams for the official `pyrtc-view`
application. The PSFs are normalized by the exact diffraction-limited peak of
the maintained 750 nm reference science path. They are observations of the
lockstep exchange; they do not participate in graph scheduling or command
adoption.

This is a deterministic integration test, not a real-time transport or latency
claim. The current path copies each AOS frame to a C-contiguous NumPy array and
uses pyRTC's fixed Linux shared-memory stream names. It refuses to start if any
of those names already exist, so stop a live pyRTC system before running it.
External transport, GPU-to-CPU transfer policy, asynchronous pacing, dropped
frames, and deadline behavior remain separate work.

For process separation without embedding Python in Julia,
[`pyrtc_shared_memory.jl`](../examples/integrations/pyrtc/pyrtc_shared_memory.jl)
implements the numeric part of the current Linux pyRTC `ImageSHM` layout with
POSIX shared memory and fixed-shape, non-owning `UnsafeArrays.jl` views. AOS
calls this compatibility contract “v1”; pyRTC does not publish it as a
separately versioned protocol. A stream has a payload segment named `<name>` and
a ten-`Float64` segment named `<name>_meta`. The metadata
records the write count, publication time, byte count, NumPy dtype code, and up
to six dimension slots. The current Julia adapter admits numeric vectors and
matrices, validates that the remaining slots are zero, checks exact segment
sizes and type metadata, and converts matrices explicitly between Julia and
NumPy C order. The creator alone unlinks a stream; every participant closes
its own mappings.

The upstream layout has one overwriteable payload slot and no in-progress
publication marker. The maintained adapter therefore requires one producer and
one outstanding frame in application-level lockstep. The producer writes the
payload before the count and timestamp; the consumer copies into its own
preallocated array and rechecks the metadata. This detects a publication that
changes during a copy, but it cannot make a free-running overwrite race safe.
An asynchronous deployment needs a new versioned protocol with slot ownership
or double buffering. Warmed vector and matrix publication and immediate reads
perform zero Julia heap allocations; this is not a hard-real-time or syscall
claim.

[`pyrtc_process_hil.jl`](../examples/integrations/pyrtc/pyrtc_process_hil.jl)
uses that adapter to create `wfs`, `signal`, `signal2D`, and `wfc` in Julia. A
separate Python worker attaches pyRTC's real `SlopesProcess` and `Loop`, while
Julia retains graph stepping, calibration, atmosphere evolution, science-image
formation, and stream lifetime ownership. The short optional process test first
closes a known static disturbance and then requires both reference sensor paths
to improve the mean on-axis Strehl ratio of an evolving atmosphere. The control
pipe intentionally advances one operation at a time, so this validates process
and data-boundary correctness rather than free-running RTC timing. Instrument
packages reuse this transport contract and own their larger calibration and
acceptance tests.

The user guide also gives the commands for the bidirectional protocol matrix,
the opt-in separate-process test matrices, and the live demonstrations.

The example transport functions above are application placeholders, not AOS
APIs. The boundary deliberately defines no socket, PipeWire buffer, wall-clock
deadline, timeout, last-command hold, queue, or concurrent callback. A direct
Julia script can perform the loop synchronously; an external adapter can map the
same buffers and sequences to its transport lifecycle.

## YAML Algorithm Graph Format

YAML is the maintained human-authored format because it supports comments and
keeps sequences of node and link mappings readable. File loading is provided
by the optional YAML.jl extension, activated with `using YAML`.
`load_algorithm_graph(root::AbstractDict; ...)` also accepts
an already-decoded dictionary, so an application may use its own interchange
decoder and pass the resulting graph through the same schema validation and
node factories. This entry point does not define or maintain another serialized
format. Direct Julia construction remains the unrestricted API for generated
topology, conditional composition, multi-rate orchestration, and other
arrangements that do not fit the static file subset.

The key words **MUST**, **MUST NOT**, **REQUIRED**, **SHOULD**, **SHOULD NOT**,
and **MAY** in the following requirements are to be interpreted as described by
[BCP 14](https://www.rfc-editor.org/info/bcp14) when, and only when, they appear
in all capitals.

- **AOS-GRAPH-FILE-001:** A maintained graph file MUST set
  `schema_version: 1`, MUST have a stable identifier-valued `name`, and MUST
  contain a `nodes` sequence with at least one mapping. A loader MUST reject
  unknown schema versions and unknown fields. A YAML file MUST contain exactly
  one document, with unique keys in each mapping.
- **AOS-GRAPH-FILE-002:** The `nodes` sequence MUST define declaration order.
  Preparation MUST derive a fixed execution order from direct links, selecting
  the earliest declared ready node at each step, and MUST reject direct cycles.
  Delayed links MUST NOT constrain same-step execution order. Each node MUST
  identify an explicitly supplied `type`; the loader MUST NOT evaluate Julia
  source, resolve arbitrary module paths, or mutate a global type registry.
- **AOS-GRAPH-FILE-003:** Each node's `config` mapping MUST contain construction
  and graph-rebuild values. Each node's `props` mapping MAY contain scalar initial
  props but MUST NOT change a node's admitted ports. Version 1 does not define
  live prop transactions.
- **AOS-GRAPH-FILE-004:** Frame inputs, delayed-link initial arrays, sparse
  parameter arrays, and cold node calibration arrays MUST be referenced by
  external binding name. The core
  loader MUST NOT infer FITS or another calibration file format and MUST NOT
  insert an implicit host/device transfer into graph connectivity. Frame inputs,
  links, and sparse parameters MUST retain the native exact-device admission
  checks. Explicit native calibration snapshots and per-device calibration
  preparation MAY perform the transfers declared by their native constructors.
- **AOS-GRAPH-FILE-005:** Endpoints MUST use `node.port` syntax. Direct links,
  delayed links, boundary inputs and outputs, and sparse parameters MUST pass
  the same type, shape, schema, layout, ownership, and connectivity admission
  checks as direct Julia construction.
- **AOS-GRAPH-FILE-006:** `load_algorithm_graph` MUST compile a valid YAML file
  or validated parsed root dictionary into the same concrete
  `AlgorithmGraphDefinition` used by `algorithm_graph`. Parsed dictionaries
  and cold builder vectors MUST NOT survive into the definition or prepared
  execution owner.
- **AOS-GRAPH-FILE-007:** Version 1 is a static, single-rate, complete-frame
  subset. Multi-rate scheduling, conditional topology, generated nodes, and
  physical sub-frame event semantics SHOULD use direct Julia composition
  rather than acquire implicit file semantics.

The version 1 mappings and sequences are:

| Mapping or sequence | Required fields | Meaning |
|---|---|---|
| root | `schema_version`, `name`, `nodes` | Schema and graph identity |
| `nodes` | `name`, `type`, `config` | One node in declaration order; preparation derives execution order; `props` is optional; `type: detector` also requires a `detector` selector; optional `bindings` maps construction roles to external array names |
| `inputs` | `name`, `destination`, `binding` | Caller-owned frame input |
| `outputs` | `name`, `source` | Caller-visible allocated graph output |
| `links` | `source`, `destination` | Same-step direct link |
| `delayed_links` | `source`, `destination`, `initial` | One-successful-step delayed link; `initial` is a binding name |
| `parameters` | `destination`, `binding` | Required startup ndarray sparse parameter |

For example, an AOS YAML graph may describe the optical plant and detector
sequence. Reconstruction, controller recurrence, and VDM/PDM RTC routing are
maintained by
FilterGraphAlgorithms/JuliaFilterGraph. Use the
[`filter_graph_algorithms` fixture](../examples/integrations/filter_graph_algorithms/)
for the complete AOS-plant/FGA-RTC composition.

Use the maintained plant graph files under `examples/graphs/` for AOS-native
optical and detector composition. The application supplies graph bindings and
prepares that plant graph with `load_algorithm_graph`; the RTC graph and its
calibration products are prepared by FGA/JFG in the maintained
[`filter_graph_algorithms` fixture](../examples/integrations/filter_graph_algorithms/).

The built-in type map currently contains `detector`, `pixel_mtf`,
`deformable_mirror_surface_f32`,
`gaussian_deformable_mirror_surface_f32`,
`modal_opd_expansion_f32`,
`second_order_response_f32`,
`multilayer_atmosphere_opd_f32`,
`spatial_ao_residual_f32`,
`pupil_opd_composition_f32`,
`pyramid_rate_f32`,
and `shack_hartmann_rate_f32`.
`merge(builtin_graph_node_types(), companion_types)` creates an explicit larger
map for optional packages such as the Proper companion. Cold factories accept
`(name, config, props)`. Factories that need external calibration arrays may
accept `(name, config, props, bindings)` instead; the loader selects that
signature when applicable. They resolve explicit binding names from the same
caller-supplied named tuple used by inputs and sparse parameters. Cold factories
must convert these bindings into the native node's configuration ownership;
the binding tuple does not become execution state. Detector factories use
`(name, config, props; detector)` with the scientific selector as a keyword.

### Idealized spatial AO residual

`spatial_ao_residual_node` declares an independent simplified facility-AO
residual. It accepts the finite multilayer-atmosphere keywords plus an explicit
`Atmospheres.SpatialAORejection` and `residual_opd_schema` instead of
`atmosphere_opd_schema`. The sole output is `residual_opd` in metres on native
`[x,y]` axes. A downstream physical DM remains a separate node and command path.

The YAML type `spatial_ao_residual_f32` accepts the raw finite-atmosphere fields,
replaces `atmosphere_opd_schema` with `residual_opd_schema`, and additionally
requires `control_cutoff_cycles_per_m` and `residual_amplitude`. No rejection
strength or facility calibration is inferred. Each full 3N×3N periodic layer
screen is filtered before pupil extraction:

```text
φ_res = Re IFFT(H FFT(φ))
H(0) = 0; H(f) = a for 0 < ‖f‖ ≤ f_c; H(f) = 1 above f_c.
```

The dimensionless amplitude fraction `a` lies in [0,1]; its PSD fraction is a².
Physical frequencies use the full screen sampling Δ = telescope_diameter_m/N,
not cycles per array element or per detector pixel. The cutoff must be positive
and no greater than 1/(2Δ). This real-even filter is applied once during
preparation and after each reset, using retained backend FFT workspace. It
adds no per-step FFT or RNG consumption. Moving-screen wind, source direction,
model-time cadence, captured replay and seeded reset retain their ordinary
semantics: each invocation publishes the next `atmosphere_step` epoch.

This is an instantaneous noiseless spatial model, not a facility fit or temporal
closed loop. Hard-cutoff ringing, finite periodic-screen recurrence and existing
bilinear interpolation remain explicit limitations. Full-screen piston removal
does not remove the mean of each pupil crop. No per-frame RMS normalization is
applied. Ordinary and captured execution require backend parity tests; seeded
random realizations are reproducible on a backend, not promised bit-identical
across backends. A sampled residual does not qualify an instrument's sensor
observability, calibrated operating range, noisy-loop science or deadlines.

Detector nodes use a short scientific selector rather than a Julia implementation
name:

```yaml
- name: camera
  type: detector
  detector: emccd
  config:
    shutter: global
    # Exposure, sensor response, noise, schemas and RNG settings follow here.
```

The supported selectors are `ccd`, `cmos`, `emccd`, `ingaas` and
`hgcdte_linear_avalanche_array`. They select
the existing Float32 complete-frame acquisitions; UInt16 camera output remains a separate
ADC/formatting step. `shutter` defaults to `global`. Other shutter modes and
unsupported detector names are rejected, rather than approximated silently.
Rolling exposure and frame-transfer timing require explicit acquisition
composition; this selector does not implement those temporal models.

The declarative subset does not expose every native Julia graph feature:

- The atmosphere factory renders one infinite-distance source. Its
  `source_separation_arcsec` and `source_position_angle_deg` must describe the
  same source as the WFS consuming that atmospheric OPD. They default to on axis.
  Multiple directions sharing one evolving atmosphere and finite-height sources
  still require explicit Julia composition.

`detector: ingaas` uses the existing complete-frame InGaAs acquisition with
metric `pixel_sampling_m: [row_pitch_m, column_pitch_m]` and an explicit
integrated `spectral_channel` identifier. Required fields are `rows`, `columns`,
`pixel_sampling_m`, `spectral_channel`, `exposure_duration_s`,
`quantum_efficiency`, `rng_seed`, `photon_rate_schema`, and `frame_schema`.
Optional scalar fields are `gain` (default 1), `dark_current_e_per_pixel_s` and
`glow_e_per_pixel_s` (default 0), `photon_noise` (default true), `readout_noise`
(default false), and `readout_noise_e` (default 0). Optional `bits` requires
`full_well_e`; omitting both leaves saturation and ADC quantization disabled.
Glow is charge per integration time, rather than read-cadence glow. Scalar QE
is effective across the declared channel; no wavelength-dependent response,
pixel blur, camera profile, or axis conversion is inferred.

`detector: hgcdte_linear_avalanche_array` uses the same required metric-grid
fields and optional scalar fields as `ingaas`, with the existing HgCdTe
linear-avalanche numerical sensor and single-read sampling. Additional optional
fields are `avalanche_gain` (default 1, finite and at least 1),
`excess_noise_factor` (default 1, finite and at least 1), and
`multiplication_model`: `conditional_gamma` for the CPU reference or
`clipped_gaussian_approximation` (default) for the numerical sensor's portable
approximation. The latter is qualified for q/(F−1) ≥ 25; F = 1 is deterministic.
`readout_noise_e` applies after avalanche multiplication, before `gain`.
`full_well_e` bounds multiplied charge. SAPHIRA settings remain explicit
instrument data. This adapter does not infer detector calibration or read
cadence; `read_duration_s` and nondestructive sampling settings are rejected.
The [synthetic HgCdTe fixture](../test/graph_files/hgcdte_avalanche.yaml) shows
the complete YAML contract with explicit test coefficients.

Optional persistence and nonlinearity use the native model names and explicit
modes:

```yaml
persistence_model:
  mode: exponential
  coupling: 0.25
  decay: 0.5
nonlinearity_model:
  mode: saturating
  coefficient: 0.01
```

Both models default to `mode: none`. The `none` mode accepts no coefficients;
`exponential` requires `coupling` and `decay` in [0, 1], and `saturating`
requires a finite nonnegative `coefficient`. Saturating nonlinearity applies
q → q / (1 + c q) to charge q in electrons, with `coefficient` c in electron⁻¹,
before full-well clipping.

Persistence `coupling` and `decay` are dimensionless per completed frame:
lₖ₊₁ = decay × lₖ + coupling × qₖ. For this
sensor, qₖ is the post-nonlinearity, post-saturation, post-charge-coupling frame
before read noise, gain, correction, quantization, and background subtraction.
The coefficients specify no elapsed-time constant or calibrated trap model.
These non-null models use ordinary stream execution under the existing native
admission rules.

`type: pixel_mtf` declares the separate deterministic spatial response:

```yaml
- name: response
  type: pixel_mtf
  config:
    input_shape: [3, 5]
    padding: 2
    input_schema: instrument.photon-rate/1
    output_schema: instrument.photon-rate/1
  bindings:
    transfer: calibrated_transfer
```

The node's `bindings` mapping declares calibration dependencies independently
of scientific-package loading, so application exporters can include them in
artifact identity and deployment seals. Binding roles must not also occur in
`config`. `AlgorithmGraphs.graph_node_configuration(factory, config, binding_names)`
validates and merges those roles for factory construction. It rejects omitted
or extra dependencies against the requirements declared by
`AlgorithmGraphs.graph_node_bindings(factory, config)`. That interface returns
a tuple of binding-name symbols and defaults to `()` for factories requiring no
such arrays. A calibration-bearing factory specializes that method; the pixel
MTF factory requires the external name mapped to `bindings.transfer` (the
merged factory configuration contains that name as `config.transfer`).
Discovery does no I/O and is independent of graph execution. Deployers can load just those named
arrays through their existing artifact interface, and pass them to the cold
factory without loading unused profile calibrations. Native construction remains
authoritative for element type, shape and scientific calibration constraints.

The application supplies `bindings=(; calibrated_transfer, ...)`, where the
real finite matrix has shape `(rows + 2padding, columns + 2padding)`, unit DC at
`[1, 1]`, and zero-frequency origin at that first element. The Float32 native
node snapshots this calibration during loading, then prepares it for the graph's
compute device. YAML embeds no transfer matrix or calibration filename. The
`input` and `output` ports preserve the declared rectangular shape; the filter
zero-pads, applies the transfer, crops, and does not clamp. Link its `output`
to the detector's `photon_rate` input with matching schemas. See the executable
[InGaAs and pixel MTF fixture](../test/graph_files/ingaas_pixel_mtf.yaml).

`second_order_response_f32` declares same-schema vector ports `input` and
`output`. Its required config fields are `natural_frequency_hz`, `damping_ratio`,
`sample_period_s`, `count`, and `schema`; it has no runtime properties or sparse
parameters. It advances a physical response once per graph step and produces the
position after that held-input interval. Place it upstream of modal expansion
when those modal coordinates have the declared response, and match its sample
period to the application's cadence.

### Downstream instrument packages

AOS owns the reusable plant graph nodes, static YAML schema, lockstep HIL
boundary, and pyRTC-compatible shared-memory adapter. Instrument command geometry,
detector settings, calibration policy, complete graph files, and scientific
acceptance tests belong to downstream packages:

- [REVOLTClassicSim.jl](https://github.com/DarrylGamroth/REVOLTClassicSim.jl)
  owns the Classic Shack–Hartmann graph, C-BLUE One IMX425 profile, HSDM277
  geometry, RTC-reference composition, viewer, and full pyRTC calibration and
  atmosphere test.
- [REVOLTCopperSim.jl](https://github.com/DarrylGamroth/REVOLTCopperSim.jl)
  owns the Copper Pyramid graph, HSDM277 geometry, calibration graph, viewer,
  and opt-in pyRTC process test.

Those packages bind AOS arrays and adapters directly; they do not require an
instrument-specific scheduler in AOS. Their provisional influence functions,
registration assumptions, calibration artifacts, and model-validity claims are
documented and tested with the instrument definition that owns them. A future
SPIDERS package can use the same graph protocol and the optional Proper
companion seam without adding SPIDERS-specific topology to AOS.

## Direct Julia Composition

A graph file is intentionally not the only composition surface. Direct Julia
code should be used when topology is generated, rates differ, execution is
conditional, or a model needs sub-frame sampling. The application owns that
schedule explicitly and calls the same prepared scientific operations used by
graph nodes.

Rolling-shutter image formation is one example. A complete detector frame can
remain the node output while application code evaluates the optical path at the
row or row-group times needed to model moving turbulence. The detector state,
its RNG, and the row integration are one stateful owner; no partial frame is
published.

## Backend And Parallel Execution

A prepared graph has one exact compute target. Core does not insert transfers or
fall back to the CPU. CUDA and AMDGPU arrays are supported by graph nodes whose
underlying algorithms support those backends. Host/device exchange belongs at
an explicit boundary such as `PreparedGraphHILBoundary`.

The graph executor is deliberately serial and single-writer. Parallelism should
be coarse grained across independent sources, time steps, or parameter sweeps.
Avoid nesting graph parallelism inside backend or BLAS parallelism. Deterministic
runs should keep thread counts at one.

## Validation And Performance Evidence

Use the smallest registered test selector that covers a change:

~~~bash
julia --project=. test/ci/run_selected_tests.jl algorithm-graphs
julia --project=. test/ci/run_selected_tests.jl core
~~~

Run the complete local suite for cross-cutting changes or a release gate. GPU
validation uses the dedicated AMDGPU host and the local CUDA host described in
[`backend-validation-guide.md`](backend-validation-guide.md).

The graph contract is bounded and allocation-free after preparation for the
covered native nodes. That is a per-node and per-backend claim, not a blanket
promise for arbitrary Julia callbacks.

## Archived Event Runtime

The former `AdaptiveOpticsSim.Plant` event runtime was removed from the main
package after the complete-frame graph direction was adopted. Its last merged
state is preserved in Git branch `archive/plant-runtime-2026-08-30`; the
in-progress concrete RNG-owner experiment is preserved in
`archive/plant-pe06b2-wip-2026-08-30`. New code must not depend on either branch.

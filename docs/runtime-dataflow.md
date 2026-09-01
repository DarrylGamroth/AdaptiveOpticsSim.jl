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
   results. It includes atmosphere evolution, detector persistence, controller
   history, graph delays, and explicit RNG state.
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

One graph step invokes every node once under a schedule that preserves validated
declaration order. A direct link exposes an earlier node's ordered output to a
later node in the same step. The default uses one retained accelerator stream,
so that dependency is device execution order and requires no host barrier. The
explicit grouped-stream policy may overlap independent nodes, but it admits
only groups whose flattened names equal declaration order and requires every
direct link to cross a group boundary. A delayed link exposes only the value
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

The sequence-zero command already present in the graph produces frame one.
Command `n` is a response to frame `n` and becomes active only for frame
`n + 1`. The boundary rejects missing, stale, skipped, non-`UInt64`, and
nonfinite command responses before mutating the graph command. It stops after
an execution or target-copy failure because either graph products or the
active command may then be partial. `reset_hil_boundary!` restores the
snapshotted initial command and graph sequence. When a model-time driver is
used through `step_hil_frame_at!`, reset both owners with
`reset_hil_boundary!(boundary, driver)`.

### Fully declared HIL reference systems

The maintained SHWFS and Pyramid HIL reference systems validate external-RTC
integration without depending on unmeasured instrument calibration. They share
one explicit 25-command analytic Gaussian deformable mirror, accept a complete
caller-owned uncompensated pupil OPD in metres, and publish a complete detector
frame through `PreparedGraphHILBoundary`. Their parameters are AOS validation
choices; neither system is an instrument profile.

[`shack_hartmann_hil_reference.toml`](../examples/graphs/shack_hartmann_hil_reference.toml)
publishes a noiseless 64-by-64 CCD frame and an AOS centroid signal for
differential checks. [`pyramid_hil_reference.toml`](../examples/graphs/pyramid_hil_reference.toml)
publishes a noiseless 36-by-36 four-pupil EMCCD frame. The helper in
[`hil_reference_systems.jl`](../examples/support/hil_reference_systems.jl)
owns the exact analytic actuator coordinates, SHWFS valid-lenslet rule, graph
selection, and preparation recipe.

The current integration matrix is deliberately small:

| WFS | Detector response | Controller boundary | Evidence | Status |
|---|---|---|---|---|
| Shack–Hartmann | Deterministic, noiseless CCD | AOS lockstep reference controller | `algorithm-graphs` interaction-matrix and convergence tests | covered |
| Pyramid | Deterministic, noiseless EMCCD | Command/frame lockstep | `algorithm-graphs` complete-frame and command-response tests | covered |
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
PipeWireAO transport, GPU-to-CPU transfer policy, asynchronous pacing, dropped
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
  physical sub-frame event semantics SHOULD use direct Julia composition
  rather than acquire implicit file semantics.

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
name = "reference_controller"

[[nodes]]
name = "controller"
type = "discrete_integrator_f32"

[nodes.config]
extent = 16
sample_period_s = 0.001
input_schema = "org.example.residual-modes/1"
output_schema = "org.example.command-modes/1"

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
    "reference_controller.toml";
    bindings=(residual=residual_modes,),
)
graph = prepare_algorithm_graph(definition; target=target)
```

The built-in type map currently contains `ccd_detector_acquisition_f32`,
`closed_loop_correction_f32`, `control_matrix_reconstruction_f32`,
`cmos_detector_acquisition_f32`,
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

### Downstream instrument packages

AOS owns the reusable graph nodes, static TOML schema, lockstep HIL boundary,
and pyRTC-compatible shared-memory adapter. Instrument command geometry,
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
validation uses the dedicated AMDGPU host and the WSL CUDA host described in
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

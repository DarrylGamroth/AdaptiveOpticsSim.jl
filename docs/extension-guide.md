# Extension Guide

Status: active

## Scope

Extend the canonical domain that owns the operation. Do not introduce a generic
package-wide algorithm base class or `process!` protocol. Optional dependencies
belong in Julia package extensions or companion packages.

All extensions follow the ownership rules in
[`maintainer-architecture.md`](maintainer-architecture.md): immutable params,
run-immutable plans, persistent state, replaceable workspace, caller-visible
products, and one concrete prepared owner.

## Prepared Numerical Operations

Plant-side response calibration composes the physical WFS with a caller-owned
estimator. Supply an explicitly ordered `WFSMeasurement` vector with declared
units and kind, then
pass a callable `(measurement, pupil)` to
`Calibration.interaction_matrix` or
`Calibration.compute_meta_sensitivity_matrix`. The callable writes the exact
storage returned by `measurement_storage(measurement)`; it may retain a source,
detector, and prepared estimator in its concrete closure or callable fields.
The calibration methods restore DM coefficients and pupil OPD if the callable
throws. AOS does not define a universal estimator or slope accessor for its
physical WFS families.

A new repeated operation should normally provide:

1. an immutable configuration/params value
2. a preparation function that validates shape, units, backend, and device
3. a concrete plan for immutable numerical data
4. separate state when later results depend on previous calls
5. a workspace for replaceable scratch
6. caller-owned or owner-held products with explicit metadata
7. a domain-specific mutating execution function
8. reset only when the operation has resettable state

Preparation may allocate and fail. The warmed operation should reuse storage,
avoid dynamic dispatch, and have a measured allocation budget.

Do not put exact caller-array identity into a reusable plan. Bind it in the
prepared owner. Do not store `Any` or an abstract array element type to make an
extension appear generic; specialize at preparation or add a function barrier.

## Reusing Prepared Optical Execution

The shared optical execution code accepts numerical plans and array backends;
instrument configuration, including YAML, supplies the prescription and
calibration. Prefer these
existing owners over copying a fast kernel into another instrument:

| Reuse point | Existing use | Conditions to preserve |
| --- | --- | --- |
| PROPER prepared field arithmetic | Compiled optical paths and prepared scalar trains | Fixed wavelength, sampling and geometry; serialize shared FFT scratch. The public field API requires caller fields separate from that scratch. The field operation does not advance wavefront geometry. |
| Scalar optical-train owner | Scalar, runtime-OPD and per-wavelength polychromatic trains | Keep live OPD inputs live, restore the prepared geometry, and retain the owner’s exact input/output bindings. |
| `Optics.prepare_spectral_synthesis` / `synthesize_spectrum!` | Calibrated multi-path spectral synthesis | Fields already share one detector grid; preserve path/anchor order, coherence, quadrature, response and normalization. See the [spectral recipe](model-cookbook.md#synthesize-calibrated-spectral-channels). |
| `CapturedGraphExecution` | Prepared algorithm graphs with capture-safe nodes | Fixed buffer/plan identities and device-work topology; input contents may change between completed replays. Host-side effects are not replayed. |

Compiled paths use PROPER's public `prepare_propagation` / `propagate!` surface,
with private lowerings for sibling batching and terminal cancellation. Scalar
trains coordinate private in-place kernels whose live field is the FFT scratch.
These lowerings require their owners' additional aliasing and scheduling proofs.
Use the public surface when starting a new prescription.

Polychromatic train preparation constructs the same scalar owner per wavelength,
so it already shares that implementation. Measure performance for each new
spectral configuration. New fixed-geometry PROPER prescription segments are
also possible callers of the public field-propagation API; keep changing DM/RNG
state outside immutable geometry preparation and validate against ordinary
`prop_*` execution.

The companion PipeWireAORTC external-pupil owner already batches completed
camera publications per optical source step. That transport concern stays in
RTC; optical kernels should produce fields and products without a PipeWire
dependency.

Native AOS Fresnel propagation and Pyramid modulation already have their own
prepared paths. Their metadata, flux normalization, phase arithmetic and batch
contracts differ. Reuse the ownership and validation patterns there; replacing
their algorithms requires separate numerical and performance evidence.

Backend qualification remains separate from optical configuration. A portable
kernel does not establish that every batch size, phase-cache policy or GEMM
choice is faster on another card. Keep the ordinary executors as numerical
oracles and retain conservative fallbacks for unqualified combinations.

### Compiled-Path Source Layout

Within the PROPER extension, `ext/compiled_paths.jl` contains the ordered
executor, ordinary/prepared propagation, checkpoint continuations and independent
field batches. Its included files separate other responsibilities:

- `compiled_paths/types.jl`: immutable numerical contracts and concrete storage owners
- `compiled_paths/fields.jl`: normalization, shifts, crops, padding and FFT/filter primitives
- `compiled_paths/scheduling.jl`: cold group selection, checkpoint boundaries and shared/unshared instruction ranges
- `compiled_paths/preparation.jl`: validated program loading, shared-prefix preparation and workspace construction
- `compiled_paths/execution.jl`: bounded instruction selection with statically specialized transform-group arities

Batch-specific admission stays beside its executor so the dependency checks and
the storage they protect can be reviewed together. All files extend the same
module; splitting source files adds no dispatch or execution layer. Preserve
both arithmetic order and checkpoint ownership when changing these routines.

Fixed geometry is resolved when preparing the model. `PathExecutionPlan` stores
only execution decisions and indices in concrete fixed-size host arrays; live
OPD/amplitude inputs, FFT scratch and exact bindings stay in the prepared owner.
Group selection and terminal cancellation cannot cross a checkpoint. Separate
shared/unshared ranges retain the ordinary traversal's cache ownership. Rebuild
these decisions when changing the prepared backend or propagation policy.

The finite group arities use literal `Val(2/4/6/8/10)` calls, and paired terminal
FFT handling uses a union of two concrete types. Wavelengths, grid dimensions,
coefficients, path length and complete YAML topology remain data rather than
creating a new method family for each prescription. This bounds specialization
while allowing Julia to compile the numerical kernels for their concrete types.
The ordinary executor does not consult the prepared instruction stream. Captured
GPU replay launches the recorded graph; this scheduling cleanup alone does not
remove a GPU kernel or establish a throughput improvement.

## Algorithm Graph Adapters

A graph adapter is separate from the scientific implementation. It declares
ports, prepares one exact owner around domain plans/state/workspace, steps it,
and resets it.

~~~julia
import AdaptiveOpticsSim.AlgorithmGraphs as AOG

struct GainNode end

struct GainNodeConfig
    extent::Int
end

struct PreparedGainNode{G,I,O}
    gain::G
    input::I
    output::O
end

function AOG.graph_node_ports(::Type{GainNode}, config::GainNodeConfig)
    config.extent > 0 || throw(ArgumentError("extent must be positive"))
    shape = (config.extent,)
    return (
        AOG.graph_port_contract(
            :input, :input, :data, Float32, shape,
            "example.signal.f32/1", :column_major,
        ),
        AOG.graph_port_contract(
            :output, :output, :data, Float32, shape,
            "example.signal.f32/1", :column_major,
        ),
        AOG.graph_port_contract(
            :gain, :input, :parameter, Float32, shape,
            "example.gain.f32/1", :column_major,
        ),
    )
end

function AOG.prepare_graph_node(
    ::Type{GainNode},
    ::GainNodeConfig,
    props,
    inputs::NamedTuple,
    outputs::NamedTuple,
    parameters::NamedTuple,
    target,
)
    return PreparedGainNode(
        parameters.gain,
        inputs.input,
        outputs.output,
    )
end

function AOG.step_graph_node!(owner::PreparedGainNode)
    @. owner.output = owner.input * owner.gain
    return nothing
end

AOG.reset_graph_node!(::PreparedGainNode) = nothing
~~~

The real numerical implementation should usually be a domain function called by
`step_graph_node!`, not code duplicated in the adapter.

Port contracts declare:

- a stable name
- `:input` or `:output` direction
- `:data` or `:parameter` role
- concrete isbits element type
- fixed positive shape
- scientific schema string
- `:column_major` layout

The graph validates exact target, backend storage, element type, rank, shape,
schema, layout, connection direction, and single ownership before calling
preparation. A sparse parameter is initial prepared ndarray data; it is not a
per-frame property.

Initial scalar `props` are construction-time input to the adapter. If a node
needs runtime property updates, define an explicit node-specific transactional
owner before exposing that behavior; do not mutate a plan in place from an
uncoordinated task.

## Adding A Built-In Declarative Graph Node

A Julia graph adapter does not automatically become available to the declarative
graph loader. Built-in graph nodes also need:

- a stable type label in `builtin_graph_node_types`
- strict parsing of the node's configuration and props
- explicit external binding names for large arrays or application objects
- executable tests for file loading and, when used, parsed-dictionary loading
- documentation of scientific schemas and units

YAML file loading is provided by the optional YAML.jl extension, activated with
`using YAML`; dictionary compilation remains available without it. Keep the file
format static and declarative. Applications may decode another
interchange representation and pass its parsed root dictionary to
`load_algorithm_graph`; that does not add a second maintained file format. A
value that changes port shape or topology is configuration and requires graph
preparation. Large matrices, masks, influence functions, and calibration data
are ndarray parameters or external bindings, not scalar property text.

## Backend Extensions

A backend extension owns methods for the canonical `Backends` functions rather
than introducing a parallel abstraction. It must provide exact device identity,
compatible allocation, execution style, and synchronization needed by the
covered operation.

A GPU-ready scientific method must:

- avoid scalar indexing
- keep arrays on one admitted device
- make transfers explicit
- use AbstractFFTs or KernelAbstractions where appropriate
- preserve caller-visible product metadata
- validate unsupported algorithms before hot execution

CPU fallback is not implicit. If an algorithm is not supported on a target,
throw a structured preparation error.

## Optical And Atmosphere Extensions

New physical optics belong in `Optics`. Use explicit plane metadata,
coordinates, spectral coordinates, normalization, spatial measure, and
combination policy. New atmosphere models belong in `Atmospheres` and should
separate evolving state from direction renderers.

A direction renderer consumes a published atmosphere epoch; it does not advance
the atmosphere or consume RNG. A batch freezes a homogeneous direction axis and
writes caller-owned output.

## Detector Extensions

A detector family belongs in `Detectors` and composes the shared acquisition
pipeline where possible. Keep these concerns distinct:

- detector-facing photon-arrival-rate product
- exposure/integration policy
- sensor response and defects
- stochastic state
- readout/sampling products
- output conversion and metadata

Validate units and ordering. Pass RNG state explicitly or own it in one prepared
node; never look up a global stream in the hot path.

## Wavefront-Sensor Extensions

A composed WFS should implement only the staged contracts that its maintained
ownership boundary requires:

1. optical formation into detector-facing products
2. acquisition into `WFSObservation`
3. estimation into `WFSMeasurement`, when the owning package provides that
   estimator

Zernike and Bi-O-edge implement the first two stages in AOS; registered
FilterGraphAlgorithms owns their calibrated estimation and RTC processing.
Shack–Hartmann, Pyramid, and Curvature likewise stop at complete detector
observations and use FGA for operational estimation. LiFT uses AOS physical
forward prediction and its analytic Jacobian callback with the registered
`AdaptiveOpticsCalibration.PhaseRetrieval` direct API for inverse policies,
solves, products, and diagnostics; this is not a detector slope-estimation
contract.

Reusable masks, microlens arrays, phase spots, and defocus optics belong in
`Optics`. Detector physics belongs in `Detectors`. Do not hide a detector inside
an estimator or use slope/flux vocabulary interchangeably.

## Calibration, RTC Composition, And Tomography

AOS calibration owns simulated physical response acquisition and plant-side
observables. AdaptiveOpticsCalibration owns reusable inverse methods and
calibration products. FilterGraphAlgorithms/JuliaFilterGraph owns runtime
reconstruction, controller state, frame delay, and VDM/PDM routing. Tomography
owns atmospheric reconstruction geometry and DM projection. Keep large
prepared arrays in plans or parameter owners and mutable histories in state.

An in-process RTC composes with the AOS plant at an explicit package boundary;
it is not an AOS algorithm-graph node. An external RTC connects through
`PreparedGraphHILBoundary`. Downstream applications own either composition and
its system-level timing and acceptance evidence.

## Proper.jl

Proper.jl is optional. An application may prepare a prescription and call it
between AOS optical products, or provide a graph adapter whose owner holds the
prepared Proper context and exact arrays. The adapter must still satisfy the
same ownership, metadata, reset, backend, and allocation contracts. See
[`proper-integration-guide.md`](proper-integration-guide.md).

## Testing Checklist

- invalid configuration fails during preparation
- element type, axes/shape, units/schema, backend, and device are checked
- warmed execution has the declared allocation budget
- state has one writer and reset is deterministic
- outputs are unchanged or the owner becomes failed when execution throws
- CPU reference behavior is covered
- applicable AMDGPU/CUDA paths run with scalar indexing disabled
- graph adapters have direct Julia and YAML tests when both surfaces are public
- new public terminology is added to [`glossary.md`](glossary.md)

## Fitted SCC on FilterGraphs

The optional FilterGraphAlgorithms extension exposes the fitted SCC propagator
and fitted pixel MTF as two Algorithm declarations. Load Proper, FITSIO, YAML
and FilterGraphAlgorithms with AdaptiveOpticsSim. Select this chain explicitly
by preparing the in-memory definition.

```julia
using AdaptiveOpticsSim, Proper, FITSIO, YAML, FilterGraphAlgorithms, FilterGraphs

pilot = AdaptiveOpticsSim.Optics.fitted_scc_filter_graph(model_path, model_sha256)
graph = FilterGraphs.prepare_graph(pilot.document; algorithms=pilot.algorithms)
```

The inputs are `dm_opd`, `pupil_opd` and `attenuation` on the compiled model's
native grids. OPD is in metres. Attenuation is the additional log-amplitude
value in `exp(-attenuation)`. The three raw optical products are `fringed`,
`unfringed` and `reference`. The `detector` product applies the fitted
pixel MTF to `fringed`. All ports use Float32. Exposure integration, detector
noise and model time are handled by the application around this chain.

On the CPU, call `FilterGraphAlgorithms.process!` with caller-owned outputs.
On GPU, call `FilterGraphs.prepare_captured_graph(graph, device)` and then
`FilterGraphs.process_device!` to consume completed device products.
GPU preparation requires native capture for the complete chain.

## Dual SCC and LLOWFS on FilterGraphs

For the fitted SPIDERS model, prepare the external pupil from its authored
optical and detector graphs. Select the pilot before advancing the model.

~~~julia
using AdaptiveOpticsSim, Proper, FilterGraphAlgorithms, FilterGraphs

pilot = AdaptiveOpticsSim.Optics.prepare_filtergraphs_external_pupil(prepared)
~~~

The pilot takes ownership of the prepared pupil and its exact optical, DM, and
detector buffers. Step and reset only the returned owner through the existing
external pupil API. If the pupil is part of a clocked source, replace its
`plant` with the returned owner before stepping the source.

FilterGraphs schedules the optical observations. GPU preparation requires native
capture of the two SCC phases and the case where one pupil sample crosses their
boundary. The SCC and LLOWFS exposure accumulators still publish on their own
model-time windows, and the detector responses still run after their exposures.
Apply controller commands at the same model sample boundary as with the original
owner. Preparing an external pupil without the opt-in call uses AOS execution.

For the direct LLOWFS port, pass the same authored LLOWFS definition used to
prepare the pupil.

~~~julia
O = AdaptiveOpticsSim.Optics
pilot = O.prepare_filtergraphs_external_pupil(prepared; llowfs_definition)
~~~

This prepares the optical train and pixel integration as separate
FilterGraphAlgorithms nodes in one LLOWFS graph. Pupil OPD, pupil amplitude,
DM OPD, and intensity are real Float32 arrays. Complex wavefields are private
to the optical train. On GPU the LLOWFS graph requires native capture and
replays separately from the SCC graph for the active chop phase. The completed
intensity is copied into the AOS exposure input on the same device.

To prepare just the LLOWFS computation, use
`description = O.llowfs_filter_graph(llowfs_definition)` and pass its
`document` and `algorithms` to `FilterGraphs.prepare_graph`. AOS still owns
the exposure clocks, detector responses, and when a DM command takes effect.

## GPI2 Pyramid photon rate on FilterGraphs

The fitted GPI2 simulator's OPD-input Pyramid can use a direct
FilterGraphAlgorithms declaration. Its real Float32 ports carry pupil OPD
and photon rate. Complex fields, diffraction plans and the modulation
workspace stay inside the Algorithm.

~~~julia
using AdaptiveOpticsSim, Proper, FilterGraphAlgorithms, FilterGraphs

O = AdaptiveOpticsSim.Optics
description = O.pwfs_filter_graph(pyramid_node)
graph = FilterGraphs.prepare_graph(description.document; algorithms=description.algorithms)

pilot_definition = O.with_filtergraphs_pwfs(aos_definition; node=:pyramid)
~~~

Prepare the standalone graph for native capture before GPU execution.
Prepare the pilot definition using the existing AOS graph API. On GPU it
records the same Algorithm into the enclosing AOS capture using the exact
device arrays already routed to and from the node. AOS retains the source,
DM inputs, model clock, detector exposure, and frame publication.

To run the fitted camera through its raw UInt16 frame, use the full camera
recipe. It follows the authored links from the Pyramid through EMCCD and
frame layout. Supply sample metadata consistent with the configured exposure
duration.

~~~julia
description = O.pwfs_camera_filter_graph(aos_definition)
camera = FilterGraphs.prepare_graph(description.document; algorithms=description.algorithms)

pilot_definition = O.with_filtergraphs_pwfs_camera(aos_definition)
~~~

Native FilterGraphs capture records all three camera nodes in one replay.
The EMCCD advances its seeded random stream when its acquisition succeeds
and restores the initial stream on reset. The AOS pilot binds the same
camera Algorithms inside the enclosing AOS capture and keeps the source,
model-time exposure, and DM command sequence under AOS ownership. The
AOS CPU path binds the Pyramid and EMCCD to the graph's buffers during
preparation. Standalone CPU execution accepts caller arrays per frame and
uses the Algorithms' preallocated staging.

With a FilterGraphs build that provides CPU loan pools, standalone callers
can register the pupil OPD and raw frame arrays during preparation and
select them independently on each sample.

~~~julia
camera = FilterGraphs.prepare_cpu_loans(description.document;
    algorithms=description.algorithms, input_slots, output_slots)
result = FilterGraphs.process_cpu_loan!(camera, input_index, output_index, sample_metadata)
~~~

The Pyramid shares its optical workspace across OPD slots. The EMCCD has
one noise stream on fixed internal buffers. Frame layout writes to the
selected output slot. The caller supplies an idle output slot and publishes
it only when the result reports a completed frame. Keep the arrays alive
and retain their shapes and storage while the camera is prepared.

For a full fitted CPU simulator with a lockstep HIL boundary, prepare
the camera inside the AOS graph with caller-owned raw-frame slots.

~~~julia
owner = O.prepare_filtergraphs_pwfs_hil(aos_definition; command_input, frame_slots)
sample = AdaptiveOpticsSim.AlgorithmGraphs.step_hil_frame_at!(owner, driver, slot)
frame = AdaptiveOpticsSim.AlgorithmGraphs.hil_frame_buffer(owner)
~~~

The simulator's existing pupil OPD is the fixed input loan. The camera's
photon rate and detector output write directly into the corresponding AOS
diagnostic arrays. The caller chooses an idle output slot and uses its
frame after a successful step. A failed step can leave scratch data in
the selected slot. Adopt the matching command sequence through the normal
HIL interface before advancing again. The prepared owner and its driver
remain single-writer. Reset does not overwrite caller-owned frame slots.

Omit `frame_slots` when the HIL transport lends a frame after preparation.
The caller passes a `Matrix{UInt16}` with reversed dimensions to
`step_hil_frame_at!`. The camera writes row-major wire order through a
logical view of that memory. Keep the loan valid and exclusive until the
step returns. The simulator retains no reference to the borrowed storage.

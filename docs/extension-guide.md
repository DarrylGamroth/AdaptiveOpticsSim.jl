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

## Adding A Built-In TOML Node

A Julia graph adapter does not automatically become available to the TOML
loader. Built-in file nodes also need:

- a stable type label in `builtin_graph_node_types`
- strict parsing of the node's configuration and props
- explicit external binding names for large arrays or application objects
- round-trip or executable graph-file tests
- documentation of scientific schemas and units

Keep the file format static and declarative. A value that changes port shape or
topology is configuration and requires graph preparation. Large matrices,
masks, influence functions, and calibration data are ndarray parameters or
external bindings, not scalar property text.

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

A composed WFS should implement the staged contracts where applicable:

1. optical formation into detector-facing products
2. acquisition into `WFSObservation`
3. estimation into `WFSMeasurement`

Reusable masks, microlens arrays, phase spots, and defocus optics belong in
`Optics`. Detector physics belongs in `Detectors`. Do not hide a detector inside
an estimator or use slope/flux vocabulary interchangeably.

## Calibration, Control, And Tomography

Calibration owns model identification and matrix/basis creation. Control owns
runtime reconstruction and controller state. Tomography owns atmospheric
reconstruction geometry and DM projection. Keep large prepared arrays in plans
or parameter owners and mutable histories in state.

A controller that participates in a simulation graph is an ordinary explicit
node. An external RTC instead connects through `PreparedGraphHILBoundary` and
is not simulated by the graph unless the application deliberately includes it.

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
- graph adapters have direct Julia and TOML tests when both surfaces are public
- new public terminology is added to [`glossary.md`](glossary.md)

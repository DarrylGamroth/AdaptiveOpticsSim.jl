# HIL Execution And Placement

Status: active

## Purpose

This specification defines ownership and static placement for prepared CPU,
mixed CPU/GPU, single-GPU, and multi-GPU execution. It separates the fixed
HIL critical path from optional offline ensemble scheduling.

See the [`HIL architecture index`](../hil-package-boundary.md) for adjacent
subsystem specifications.

## Optical Branch Ownership And Parallelism

Different source directions are the primary coarse-grained parallel unit. They
share an atmosphere epoch and effective optic commands, but require different
direction-dependent layer footprints. Detectors on the same compatible path
reuse the propagated field or photon-arrival-rate product before their independent
detector integrations and pipelines run.

The current shared-arm runtime gives each direct-science arm a prepared pupil,
field, photon-arrival-rate output, propagation workspace, and detector plans.
Its primary WFS, primary science, and auxiliary-arm paths also own distinct
`PupilFunction` products. It still executes the arm loop sequentially, so it
must not be parallelized by wrapping that loop in `Threads.@threads`. The
remaining event-runtime ownership requirements are:

- shared immutable telescope parameters plus revisioned aperture support and
  reflectivity
- shared atmosphere layer state that is read-only while its current epoch
  readers are active; an epoch token does not retain those layers
- shared immutable optic parameters and an effective command snapshot
- source-specific, precomputed propagation geometry
- path-local OPD, electric field, photon-arrival-rate output, and FFT workspace
- acquisition-local WFS, detector-integration, readout, and publication state
- destination-owned atmosphere rendering that does not use one shared render
  scratch buffer
- bounded caller-owned materialized atmospheric path products, plus optional
  model-specific retained state only when a plan permits cross-timestamp
  rendering

Preparation owns allocation, FFT planning, path grouping, worker placement,
and cache construction. The warmed event path mutates prepared state without
heap allocation.

Each acquisition also binds a prepared product provider: full optical,
reduced-order, or synthetic/replay. Provider choice changes capability,
workspace, payload traffic, and cost estimates, but not the acquisition's
external product contract. Mixed-fidelity plants are valid and the planner
places each complete provider group accordingly. The selected provider is
immutable for a run. A different fidelity tier is a separate prepared run, not
an overload transition in the initial HIL design.

The coordinator fans out compact due-work descriptors through direct calls or
one bounded SPSC ring per execution owner. Each owner processes its ordered
stream and returns completion through an owner-specific bounded path; the
design does not introduce one contended global work queue. A fixed pool of
materialized atmospheric path-product slots permits downstream overlap without
retaining mutable layer state. A model may additionally prepare a bounded pool
of retained atmosphere-state snapshots when it must render after the writer
advances. Per-acquisition integration remains single-writer and ordered even
when other acquisitions progress concurrently.

Required WFS work does not wait for unrelated science work merely because both
share a plant timestamp. The atmosphere writer does wait only until every
unmaterialized reader of its current layers has completed. After path inputs
are materialized, preparation either assigns enough independent resources and
slot depth for their declared downstream overlap or gives optional work an
explicit shed, coalescing, or fail policy. Materialization/snapshot exhaustion,
worker lag, or a full due-work ring follows that prepared policy rather than
creating an unbounded backlog.

### Prepared static placement

HIL placement operates on prepared execution groups, such as a complete
direction-dependent optical path and its compatible acquisition consumers. It
does not ask users to place individual Julia tasks, kernels, or functions.
Preparation resolves every group to one execution owner and memory domain.
The resulting plan is immutable for the run.

Small fixed stage pipelines within one group may remain concretely typed for
specialization. A large instrument's endpoint and path registry MUST NOT be
encoded solely as one recursively specialized tuple or `NamedTuple` whose type
grows with topology size. Preparation instead uses homogeneous registries or
prepared executor handles behind function barriers where needed. Validation
records preparation/compilation latency and generated-code size versus endpoint
count so low steady-state latency is not purchased with unbounded startup or
code growth.

The initial implementation supports two configuration modes with the same
prepared-plan representation:

| Mode | Contract | Intended use |
|---|---|---|
| Fully explicit | The user assigns every group and reserved execution context | Validated production deployments requiring exact placement |
| Constrained deterministic | The user supplies required, preferred, or forbidden placements and conservative static rules fill unassigned groups | Initial configurations and development; the resulting plan is inspected before promotion |

Explicit requirements override preferences and deterministic rules. A conflict
between two requirements, an unassigned explicit group, or a resource that
cannot execute the requested model is a preparation error. The planner must
not silently weaken a hard constraint.

A fully automatic cost-model optimizer is deferred. It becomes a separate
evidence-gated capability only after several real instrument profiles provide
representative single-resource and transfer measurements. Until then, measured
costs may appear as admission diagnostics but do not silently override user
placement or conservative rules.

Useful user constraints include:

- required CPU worker set, NUMA node, GPU device, or memory domain
- preferred or forbidden resources
- groups that must be co-located to reuse a field, photon-arrival-rate product,
  or prepared model
- required product-provider tier
- latency-critical and optional work classes
- required host/device boundaries for external integration
- CPU cores reserved for coordination, ports, GPU submission, or interrupts

Core models expose backend capabilities, workspace requirements, memory
domains, and placement constraints through dispatch and traits. The planner
must not maintain an `isa`-driven switch over model implementations.

After satisfying hard constraints, constrained deterministic rules consider:

- model and backend capability
- memory capacity and prepared workspace size
- the complete simultaneous due-event pattern rather than average load
- host/device and peer-device transfer cost
- field, photon-arrival-rate product, detector, and optical-model reuse
- CPU topology, NUMA locality, GPU locality, and reserved contexts
- declared deadline headroom and bounded handoff capacity

Preparation emits a structured placement plan containing the resource
inventory, group assignments, owners, memory domains, handoffs, estimated
burst utilization, estimate identity, applied constraints, and assignment
rationale. The run manifest records that plan. A rule-produced plan is
therefore inspectable and reproducible when the resource inventory, planner
version, estimate set, and stable tie-break rules are pinned. Replay may load
the recorded plan directly; rerunning a newer planner is a new placement
decision.

If no feasible assignment exists, preparation fails with structured
diagnostics identifying the violated constraint or overloaded resource. For
example, it should report that a GPU is estimated to exceed capacity during a
simultaneous LGS completion burst. Passing this admission check is not a
latency claim; fixed-arrival validation on the declared hardware remains
required.

The HIL runtime does not migrate groups or rebalance load opportunistically.
Runtime migration would change transfer paths, cache state, synchronization,
and tail latency. Overload follows the prepared stop, shed, or explicit
coalescing policy. A different placement or fidelity is prepared between runs.
Dynamic scheduling remains available only for offline workloads whose contract
permits it.

### CPU execution

CPU HIL execution should use direct calls or long-lived prepared workers with
static or deadline-aware branch ownership. It should not create a new task graph
for every acquisition event. Parallelism remains coarse over due paths; FFT,
BLAS, and Julia worker counts must be configured together to avoid nested
parallelism and oversubscription.

On large EPYC or Threadripper systems, placement should account for NUMA nodes,
physical cores, SMT siblings, memory allocation, NIC queues, and interrupts.
AcceleratedKernels may remain an evidence-gated option for sufficiently large
reusable partitions, but the direct prepared executor is the HIL baseline.

The HIL companion may use `ThreadPinning.jl` during preparation to bind
long-lived workers to physical cores or an explicit launcher-provided affinity
mask. The run manifest records the requested and observed mapping and topology.

Affinity is deployment policy, not a package-import side effect. Julia workers,
FFT workers, BLAS workers, transport agents, GPU submission agents, and
interrupt/NIC placement must be considered together. Platform support and
launcher restrictions must be checked during preparation, and the observed
mapping is recorded rather than assuming that a pinning request succeeded.
Pinning is retained only when fixed-arrival tail-latency evidence justifies it.

Julia heap allocation by a same-process RTC adapter, telemetry sink, logger, or
other task can induce process-wide GC work even when the simulation executor is
allocation-free. A same-process deployment therefore gives every long-lived
component a warmed allocation budget and keeps formatting, resizing, discovery,
and artifact encoding outside the run path. If an adapter or recorder cannot
meet that budget, it is isolated in another process and connected through an
explicit IPC boundary. The manifest records Julia GC configuration, collection
counts and pause observations, page-prefault and memory-lock policy where used,
and whether transport/telemetry shares the simulation process. Disabling GC is
not assumed safe and requires a bounded-memory soak test.

### Single-GPU execution

On a GPU, compatible directions and wavelengths should normally be batched into
an additional array dimension. One prepared device owner submits work and owns
the streams, FFT plans, allocator state, and observation barriers. Multiple CPU
tasks independently submitting small kernels to one device should not be the
default HIL strategy.

### Mixed CPU/GPU execution

One logical plant may place complete execution groups on CPU workers and one
or more GPUs. A typical layout keeps clock coordination, RTC ports, command
admission, metadata, and host-only detector work on reserved CPU owners while
placing high-rate or computationally expensive optical paths on GPU owners.
CPU and GPU groups may run concurrently when their prepared dependencies and
deadlines permit it.

Shared telescope parameters, atmosphere epoch tokens or retained/materialized
state, effective command snapshots, and deterministic RNG identities are
replicated or published with explicit lifetimes across
the participating resources. Large field, OPD, photon-arrival-rate, or frame arrays cross a
CPU/GPU boundary only at an explicit prepared handoff with bounded storage and
a measured transfer budget. An ordinary optical path should not be split
stage-by-stage across CPU and GPU merely to keep both busy.

A mixed placement may intentionally transfer a completed GPU frame to a
host-resident RTC adapter or CPU-only detector operation. That boundary is
part of the sensor latency contract and must use prepared buffers. It must not
be hidden behind an unbounded asynchronous copy.

CPU/GPU support is claimed per prepared placement and model combination.
Numerical parity, epoch and command consistency, transfer residency,
fixed-arrival tail latency, saturation, and recovery all require maintained
evidence.

CPU/CUDA and CPU/AMDGPU placements may be promoted independently; neither
claim implies that a CUDA/AMDGPU multi-accelerator placement is supported.

### Multi-GPU execution

Multi-GPU HIL execution partitions complete optical path groups across devices.
It should not pipeline successive stages of one ordinary path across GPUs or
split a small detector frame between devices, because transfer and
synchronization costs would enter the detector deadline.

A prepared multi-device plan needs:

- device-aware backend identities such as CUDA device 0 versus CUDA device 1
- one long-lived execution owner, stream set, memory pool, FFT plans, and
  branch workspaces per device
- static path placement based on measured cost and deadline utilization over
  the multi-rate schedule
- co-location of consumers that reuse a field or photon-arrival-rate product
- replication of static telescope, atmosphere, and optical-model data
- broadcast of small epoch metadata and relevant effective command segments
- local construction of controllable-surface maps on each consuming device
- pinned, preallocated host readout buffers where user integration requires
  host-resident data

Full atmosphere-screen copies should not cross PCIe on every tick. Static
screen state should be replicated during preparation, while evolution uses
shared epoch counters and deterministic per-layer random domains. If an
atmosphere model injects stochastic state, its random values must be
addressable from stable identifiers such as run seed, derivation version,
layer identity, epoch, and element so that each device reproduces the same
physical realization. Sequential host seed consumption in device-launch order
does not satisfy this multi-device contract.

Initial production validation should target homogeneous CUDA/CUDA and then
AMDGPU/AMDGPU configurations. The abstraction may permit mixed backends, but a
CUDA/AMDGPU plant must not be support-claimed before its context ownership,
numerical tolerances, and timing behavior are validated on real hardware.

A useful placement is often a dedicated device for high-rate WFS paths and a
second device for slower or more expensive science paths. Multiple GPUs are not
an automatic optimization: one small optical path will generally lose latency
to coordination.

### Multi-Process And Multi-Host Boundary

The first supported HIL deployment target is one prepared process on one host,
using its CPU and GPU resources. A multi-process or multi-host plant is a
separate support surface; multi-GPU support inside one host does not imply it.
This boundary applies to execution of the simulated plant, not to the location
of the external RTC: user transport may connect that RTC from another host.

A future distributed HIL placement must statically shard complete path groups
and additionally define:

- one authoritative plant timeline and measured inter-host clock mapping
- replication or publication of atmosphere epoch tokens, required retained or
  materialized state, effective commands, RNG identities, and configuration
  revisions
- bounded transport capacity, ordering, loss, retry, and reconnect semantics
- ownership and reclamation of cross-process product buffers
- host, process, NIC, interrupt, and failure-domain placement
- behavior for partition, host loss, delayed replica, restart, and recovery

For a single-host multi-process prototype,
[`Iceoryx2.jl`](https://github.com/DarrylGamroth/Iceoryx2.jl) is a strong
candidate data plane: its IPC publish/subscribe path provides loaned shared-
memory samples and bounded per-subscriber delivery, while events or waitsets
can provide wake-up without requiring unconditional polling. It remains a
transport-specific optional extension rather than a dependency of core or the
general HIL companion.

The prototype must demonstrate more than nominal throughput. It records every
service's publisher/subscriber count, sample and subscriber-buffer capacity,
maximum outstanding loans, history, overflow/backpressure, discovery and
dead-node cleanup behavior. A producer obtains the zero-copy benefit only when
it writes directly into an iceoryx2 loan; an existing Julia heap or GPU product
may still require an explicit bounded copy. Fan-out extends sample lifetime to
the slowest correctness-critical subscriber, so optional observers use a
separate lossy service or bounded tap. These semantics are part of the future
multi-process support claim and do not alter in-process direct/SPSC ownership.

Dynamic Dagger scheduling is not introduced into that deadline path. A
distributed deployment requires its own fixed-arrival, failure, and recovery
evidence before a MORFEO-scale or other multi-host claim is made.

## Offline Parallelism

HIL execution and offline ensemble throughput remain separate products of the
same model library.

- Direct prepared CPU execution is the CPU HIL baseline.
- KernelAbstractions implements portable data-parallel kernels.
- Backend-owned batching and streams implement GPU HIL execution.
- AcceleratedKernels is an optional, evidence-gated local CPU partitioning
  policy for sufficiently large workloads.
- Dagger remains an optional policy for independent plants, Monte Carlo runs,
  sweeps, data locality, and process/node scaling. It is not used inside one
  detector-to-RTC deadline path.

This distinction preserves useful offline scalability without importing a
dynamic distributed scheduler into the HIL timing contract.

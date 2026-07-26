# HIL Execution And Placement

Status: active

## Purpose

This specification defines ownership and static placement for prepared CPU,
mixed CPU/GPU, single-GPU, and multi-GPU execution. It separates the fixed
HIL critical path from optional offline ensemble scheduling.

See the [`HIL architecture index`](../hil-package-boundary.md) for adjacent
subsystem specifications.

## Gate 6 Prepared CPU Requirements

The following clauses refine the Gate 6 requirement rows in the
[compliance matrix](compliance-matrix.md). They define the required core
execution surface without assigning worker, affinity, pacing, ring, or
transport ownership to core.

- **HIL-CPU-001.a:** Preparation MUST derive one canonically ordered path
  execution group for every scheduled direction-dependent optical path and
  MUST bind each acquisition consumer to exactly one such group.
- **HIL-CPU-001.b:** Each group MUST have exactly one execution writer for its
  path product, numerical workspace, acquisition state, and RNG state.
  Compatible consumers of one path MUST remain in that group so they can reuse
  its field or photon-arrival-rate product.
- **HIL-CPU-001.c:** Core MUST expose independently callable,
  allocation-bounded mutating group-executor seams. Those seams MUST NOT create
  tasks, queues, or execution-clock policy.
- **HIL-CPU-001.d:** The same prepared representation MUST support a canonical
  serial executor. Deterministic mode MUST use that fallback with one Julia,
  BLAS, and FFT execution thread.
- **HIL-CPU-001.e:** A parallel CPU deployment MUST declare outer execution
  owners and Julia, BLAS, and FFT thread budgets together and MUST reject a
  configuration that would exceed its admitted execution contexts.
- **HIL-CPU-001.f:** Core MUST expose stable group identities and resource
  requirements. CPU affinity, NUMA placement, and verification of the observed
  mapping belong to HIL deployment policy and MUST NOT occur as a core package
  import side effect.
- **HIL-ATM-004.a:** The atmosphere writer MUST remain held while every due
  same-epoch reader materializes its path input, unless preparation assigns
  that reader a bounded, model-specific retained atmosphere-state snapshot.
- **HIL-ATM-004.b:** An atmosphere epoch token MUST identify a publication but
  MUST NOT be treated as retained layer storage or authorize rendering after
  the writer advances.
- **HIL-ATM-004.c:** Downstream group execution MAY overlap only after it
  consumes caller-owned materialized path products and no longer reads mutable
  atmosphere layers.
- **HIL-ATM-004.d:** Every materialized-product or retained-state slot MUST be
  prepared, capacity-accounted, and reclaimed through an explicit lifecycle.
  Unsupported retention MUST fail during preparation.
- **HIL-EXEC-006.a:** Topology-sized path and endpoint registries MUST use a
  bounded homogeneous representation or prepared executor handles behind
  function barriers; their container types MUST NOT recursively encode
  registry cardinality.
- **HIL-EXEC-006.b:** Small fixed stage pipelines inside one execution group
  MAY remain concretely specialized when their numerical kernels stay inferred
  and measured code quality justifies it.
- **HIL-EXEC-006.c:** Promotion MUST record preparation and first-use latency,
  prepared storage, and a generated-code or method-instance proxy versus a
  predeclared topology matrix. Passing warmed service time alone is
  insufficient.

## Gate 7 Single-GPU Requirements

The following clauses define the single-device capability gate. They expose
enough identity and ownership for later placement work without implementing
the Gate 9A mixed-resource planner or Gate 9B multi-GPU partitioning.

- **HIL-GPU-001.a:** Core MUST represent semantic array-backend family and
  concrete compute-device identity as separate contracts. A backend family
  MUST NOT be treated as an implicit accelerator ordinal or compute-device
  identity.
- **HIL-GPU-001.b:** Every maintained accelerator array and supported array
  view MUST report a family-qualified concrete compute device.
  Same-family/different-device storage MUST compare unequal, and prepared
  compatibility checks MUST reject it before destination mutation.
- **HIL-GPU-001.c:** Every prepared path execution group MUST expose immutable
  requirements containing its backend family, concrete compute device, and
  full-optical requirement. This description MUST NOT assign a worker, stream,
  queue, or placement policy.
- **HIL-GPU-001.d:** Preparation MUST derive batching from an explicit
  model-owned capability and immutable compatibility signature. Direction or
  wavelength samples with incompatible geometry, radiometry, numerical type,
  backend, compute device, FFT plan, or product contract MUST remain separate
  or fail during preparation; warmed execution MUST NOT regroup them
  heuristically.
- **HIL-GPU-001.e:** One prepared owner MUST submit each single-device batch
  and own its streams, FFT plans, allocator/workspace state, and completion
  barriers. Core MUST expose bounded mutating executor seams but MUST NOT
  create the HIL task, ring, pacing, affinity, or transport policy.
- **HIL-GPU-001.f:** A batch launch MUST NOT be reported complete until its
  backend completion boundary has been established. Device-ready completion,
  host observation, and host/device transfer MUST remain distinguishable
  prepared boundaries.
- **HIL-GPU-001.g:** Supported path, controllable-optic, sampled-aberration,
  detector, and acquisition state MUST remain on the declared compute device
  through its promoted execution boundary. A hidden host mirror or implicit
  host fallback MUST NOT satisfy device-resident support.
- **HIL-GPU-001.h:** Promotion MUST validate the declared CUDA and AMDGPU model
  matrix against the deterministic CPU oracle with scalar indexing disabled.
  Evidence MUST cover numerical and discrete-state parity, residency,
  synchronization, warmed allocation, first use, and a launch/submission
  proxy. Timing becomes HIL latency evidence only under the maintained
  fixed-arrival absolute-and-relative contract.

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
and cache construction. The warmed event path mutates fixed prepared storage
within its declared allocation ceiling; individual numerical kernels and
batch-state transitions may carry stricter zero-allocation contracts.

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

The initial core topology is `Plant.PreparedPathExecutionGroup`: one canonical
group per scheduled path, stored in a homogeneous bounded registry, with exact
acquisition membership. It intentionally contains no worker or placement
policy.

Core now exposes one explicit same-timestamp lifecycle around those groups:

1. `begin_optical_path_batch!` identifies the canonical due set, preflights
   every path and acquisition consumer, advances the shared atmosphere at most
   once, validates every due full-optical materializer, and returns a compact
   owner-bound `OpticalPathBatchClaim`.
2. Each due group crosses `materialize_path_execution_group!` exactly once.
   Full-optical groups write their prepared path-local input from the exact
   current epoch; reduced-order and replay groups cross the same bounded state
   machine without reading atmosphere layers.
3. `seal_optical_path_batch_materialization!` succeeds only after every due
   group completes that phase. The atmosphere writer may then be released
   because subsequent group execution reads only materialized path-local
   products and held command state.
4. `execute_path_execution_group!` may be called independently and in any
   completion order, but exactly one owner writes each group. It integrates
   only the group's acquisition members, applies its visible aberrations and
   controllable optics, evaluates path-local autonomous optics, executes its
   model, and marks that path sampled.
5. `complete_optical_path_batch!` requires every due group to be complete,
   resolves the original optical-sample generators, and reclaims the fixed
   due/status storage.

Foreign state or workspace owners, copied stale claims, changed scheduler
state, premature phase transitions, duplicate group work, incomplete
finalization, and stale atmosphere publications fail structurally. Event-loop
stepping and routed command admission reject an active batch, which preserves
the held effective-command view. Direct calls to a raw atmosphere remain an
expert model operation rather than an internal lock acquisition: the HIL
coordinator owns the single writer, and current-epoch validation rejects a
materializer if that owner advances too early.

Once a group enters materialization or execution, an unexpected model
exception is fail-stop for that prepared run. Core deliberately provides no
retry or rollback seam because a path product, acquisition owner, or RNG
stream may already be partially mutated. Gate 8 operational lifecycle owns
coordinated stop, fault publication, quiescence, and preparation of a fresh
run.

`SerialOpticalPathBatchExecutor` drives the same lifecycle in canonical group
order and remains the default oracle. `AbstractOpticalPathBatchExecutor` is a
qualified extension seam through which a HIL runtime may coordinate its own
workers. Core creates no tasks, queues, affinity policy, or completion polling.
The per-group status array assumes the declared single writer and is inspected
by the coordinator only after the external completion mechanism establishes
synchronization; it is not itself a cross-thread completion queue.

`CPUExecutionBudget` is the qualified-public admission contract for one CPU
path-group execution domain. `deterministic_cpu_execution_budget()` fixes one
admitted context, Julia thread, group owner, group-local Julia thread,
FFT-provider thread, and BLAS thread. `grouped_cpu_execution_budget` instead
declares the admitted CPU contexts, Julia default-pool size, maximum simultaneous
group owners, group-local Julia parallelism, and total FFT and BLAS threads that
one owner may use. Construction rejects outer or nested counts whose worst-case
simultaneous use exceeds the declared Julia pool or admitted contexts.

`CPUExecutionEnvironment` records the caller's observed Julia, FFT-provider, and
BLAS settings plus the CPU contexts actually admitted by launcher, affinity,
cgroup, NUMA, and reserved-service policy. `validate_cpu_execution_budget`
checks the declaration against that observation without changing any
process-global setting. Core cannot portably infer an FFT provider's configured
thread count or which host-reported logical CPUs the deployment has reserved, so
those values remain explicit. The HIL runtime still owns worker construction,
pinning, and verification of the resulting mapping.

The initial lifecycle has exactly one fixed path-local materialized-product
slot per prepared group and therefore does not pipeline two timestamps through
one group. It implements the hold-until-materialized policy, not retained
atmosphere snapshots. Retained-state pools and deeper in-flight overlap remain
separate model- and capacity-accounted capabilities.

Small fixed stage pipelines within one group may remain concretely typed for
specialization. A large instrument's endpoint and path registry MUST NOT be
encoded solely as one recursively specialized tuple or `NamedTuple` whose type
grows with topology size. Preparation instead uses homogeneous registries or
prepared executor handles behind function barriers where needed. Validation
records preparation/compilation latency and generated-code size versus endpoint
count so low steady-state latency is not purchased with unbounded startup or
code growth.

The implemented boundary uses fixed-size homogeneous `Memory` registries for
declared paths, acquisitions, controllable optics, and sampled aberrations;
prepared command endpoints, optics, paths, acquisitions, and owner RNG groups;
schedule-free selections; and event-loop execution groups. Per-owner
`@noinline` barriers recover the concrete path, acquisition, materialization,
and sampled-aberration pipeline. Small command-schema, path-product,
acquisition-product, and per-path optical-application tuples remain specialized
because their cardinality belongs to one bounded owner rather than the whole
plant. The deliberate heterogeneous-owner barrier has an explicit warmed
allocation budget; concrete numerical stage kernels retain their stricter
inference and allocation contracts.

The maintained
[`benchmark_gate6_topology_growth.jl`](../../benchmarks/benchmark_gate6_topology_growth.jl)
uses separate fresh Julia processes for event and schedule-free-selection first
use over predeclared 4-, 8-, and 16-path synthetic shapes. The labels
`development`, `nfiraos_like`, and `morfeo_like` describe cardinality only and
make no instrument fidelity, throughput, or capacity claim. The clean
[`Gate 6 topology-growth artifact`](../../benchmarks/results/gate6/2026-07-25-topology-growth.toml)
records one type fingerprint for every checked registry and entry type, zero
method-instance growth, a 1.0 native-code-size ratio, and zero meaningful
`Any` slots in representative path and acquisition kernels. From 4 to 16
paths, median fresh preparation time grows by 1.052×, preparation allocation by
1.029×, and prepared storage by 2.744×; warmed orchestration remains at most
704.04 Julia heap bytes per path per cycle. The HdrHistogram service-cost
distributions are diagnostic and establish no fixed-arrival HIL latency claim.

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
parallelism and oversubscription. The maintained four-thread proof uses one
coordinator plus three fixed, reusable path-group tasks, perturbs group
submission order over an unequal-rate science/NGS/LGS run, and compares exact
discrete state, CPU products, and RNG streams with the serial oracle. Its
`Channel` synchronization is test scaffolding, not the production HIL data
plane; Gate 8 still owns bounded SPSC completion paths.

The maintained
[`benchmark_gate6_grouped_cpu.jl`](../../benchmarks/benchmark_gate6_grouped_cpu.jl)
reuses that fixed-owner policy against the Gate 3 unequal-rate science/NGS/LGS
workload. Its contract records paired serial/grouped raw histograms, first use,
GC and allocation observations, exact replay, and the declared thread budget.
It is a closed-loop service-cost experiment without a speedup gate, affinity
claim, fixed-arrival load, or external-RTC latency boundary.
The current clean
[`Gate 6 grouped CPU artifact`](../../benchmarks/results/gate6/2026-07-25-grouped-cpu.toml)
shows the expected trade: lower optical-tail service cost but higher cheap-event
overhead and lower aggregate throughput on the unpinned local laptop CPU.

Gate 6 closes at core implementation revision `0a38576`: the grouped and
topology-growth contracts pass, and
[AdaptiveOpticsHIL PR #15](https://github.com/DarrylGamroth/AdaptiveOpticsHIL.jl/pull/15)
pins the companion package and benchmark environments to that revision while
passing their full cross-platform, ownership-stress, coverage, and
benchmark-contract checks. This validates the prepared CPU executor boundary
and the selected hold-until-materialized atmosphere policy. It does not promote
the test-scaffolding tasks into production workers or add affinity, pacing,
SPSC completion, transport, or lifecycle behavior; those remain Gate 8.

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

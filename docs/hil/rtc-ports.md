# HIL RTC Ports And Bounded Handoffs

Status: active

## Purpose

This specification defines the simulator-facing RTC boundary: canonical
submission and completion ports backed by bounded, cache-line-isolated SPSC
descriptor rings inside one HIL process. Cross-process transports and
RTC-specific wire protocols remain user-integration concerns.

See the [`HIL architecture index`](../hil-package-boundary.md) for adjacent
subsystem specifications.

## Ownership And Bounded Handoffs

Each hot mutable resource has one writer. The target ownership model is:

| Resource | Writer | Readers or consumer | Handoff |
|---|---|---|---|
| Simulation clock and due-event sequence | HIL coordinator | device/path agents | Prepared bounded event slots |
| Atmosphere epoch token and current layers | Atmosphere owner | due path materializers | Direct same-epoch read while the writer is held; the token is metadata, not retained storage |
| Materialized atmospheric path product or retained state snapshot | Assigned path/snapshot owner | downstream optical path | Direct ownership or one bounded acknowledged slot |
| Effective optic state | Command/application owner | due optical paths | Versioned snapshot after scheduled application |
| Branch workspace | Assigned CPU/GPU agent | detector pipeline on that branch | Direct ownership, no shared mutation |
| Detector integration state | Detector executor | readout publisher | Direct call within one owner or bounded ready marker |
| RTC command submission | One user-integration producer | command-ingress owner | Bounded lock-free SPSC submission ring per producer |
| Command terminal outcome | Command/application owner | Submission port's paired user-integration consumer | Bounded lock-free SPSC completion ring per submission port |
| Acquisition-product completion | Product publisher | One user-integration consumer | Bounded lock-free SPSC completion ring per acquisition endpoint |
| Released product buffer | User-integration lease owner | buffer-pool owner | Bounded SPSC return ring hidden by `release!` |
| GPU submission state | One owner per GPU | Device runtime | Direct ownership |

Multiple producers do not silently promote a port to MPSC. They use independent
ports or aggregate through one explicitly owned producer outside the HIL
critical path. Preparation assigns one command authority to each endpoint or
atomic latch group. Two producers may target disjoint endpoints, but competing
producers for one endpoint require an explicit arbitration owner outside the
canonical data plane; arrival timing does not become an implicit arbitration
policy.

A latest-value slot is reserved for publishing already applied optic state; it
is not a substitute for command admission and history. Producer and consumer
mean prepared logical owners, not a `threadid()`-indexed scratch convention;
Julia tasks may migrate unless deployment affinity is explicitly established
and verified. Rings appear only where execution-context ownership changes.
Consecutive stages owned by one worker remain direct calls.

### RTC integration ports and lock-free handoffs

`AdaptiveOpticsHIL.jl` defines the simulator-facing boundary, not an RTC
transport abstraction. The public concepts are submission and completion
ports. Each command-submission port has one paired command-completion port so
its outcome credits and ownership remain SPSC. Illustrative names are:

- `CommandSubmissionPort` for canonical command descriptors entering HIL
- `CommandCompletionPort` carrying one terminal `CommandOutcome` for each
  transferred command: applied, rejected, stale, superseded, deadline-missed,
  or failed
- `AcquisitionCompletionPort` for complete acquisition-product leases; the
  declared product may be raw or calibrated pixels, slopes, a nondestructive
  sample, another prepared sensor vector, or a scalar observation
- a possible future `ProductBufferSubmissionPort` for caller-owned registered,
  pinned, shared-memory, or device buffers

The names are not committed API. The required operations are nonblocking
`try_submit!`, `try_take!`, bounded `poll!`, and explicit `release!` operations
that return ordinary success, full, empty, closed, or stale-lease statuses.

These boundary types are not core plant-command types. Core owns the prepared
semantic payload schema, mapped plant-effective timestamp, validation,
admission, application, and terminal model disposition. A HIL submission
descriptor adds endpoint/session correlation, external timestamp-domain and
mapping metadata, payload-lease ownership, and one terminal-outcome credit. A
HIL command outcome wraps the core disposition with boundary timing and returns
that credit. Core therefore never imports a HIL port, descriptor, lease, or
outcome-credit type.

Measured actuator position, encoder state, local-controller state, health, and
other device feedback use ordinary acquisition endpoints with their own
schedule, schema, sequence, and criticality. They are not command-completion
outcomes: a command outcome reports what happened to one submitted command,
whereas feedback is a sampled plant product that may exist without a new
command.

### Prepared boundary descriptors and plant command schemas

Every command endpoint has two compatible immutable contracts resolved during
preparation. The core plant command schema declares semantic payload meaning,
while the HIL submission descriptor schema declares boundary correlation,
timestamp mapping, ownership, and credit. User integration converts its wire
representation into the HIL descriptor; the data plane does not inspect
transport-specific bytes.

The core plant command schema declares at least:

- semantic target identity and schema version
- payload element type, dimensions, units, sign convention, and actuator,
  modal, or other command basis
- calibration/configuration revision needed to interpret that basis
- absolute or incremental semantics; complete fixed-shape replacement is an
  absolute command, while partial/sparse update requires another payload schema
- finite-value, bounds, and shape validation plus whether an invalid value is
  clipped, rejected, or fails the run; state-dependent physical stroke, slew,
  settling, and dynamics checks belong to the prepared device/application layer
- sequence policy for duplicate, stale, reordered, and skipped commands
- effective-time, late-command, supersession, and plant-time command-silence
  policies

The HIL submission descriptor schema declares at least:

- boundary endpoint identity, descriptor-schema version, and run/session epoch
- submission/correlation sequence and accepted source timestamp domain
- source timestamp, mapping version, mapped plant receive/effective time, and
  applicable synchronization metadata
- compatible core plant command schema identity
- inline payload or generation-checked payload-lease reference
- paired terminal-outcome credit and boundary failure policy
- optional execution-clock ingress-liveness deadline and reset policy

Small fixed-size command payloads may be inline. Large command vectors use
prepared bounded payload storage or generation-checked immutable leases; they
are not copied into an oversized ring slot or allocated after arm. A descriptor
carries the correlated sequence, timestamp, schema identity, and payload
reference. Successful submission transfers an immutable payload lease to HIL.
The correlated terminal outcome releases that lease or proves that its prepared
slot may be reused; delayed outcome consumption remains bounded by the already
reserved outcome and payload credits. A schema mismatch after transfer is
rejected before semantic admission and still produces one terminal outcome.

A successful `try_submit!` claims both an ingress slot and terminal-outcome
credit, then transfers descriptor ownership into the ingress ring. It does not
mean that the command passed validation or was admitted to the future-effective
schedule. The completion port later returns the correlated terminal outcome.
Enqueue, timestamp mapping, core validation, admission, effective time,
application, model disposition, and HIL outcome use distinct terms throughout
the API and telemetry.

The submission producer derives available outcome credit from its paired
completion-consumer sequence: queued outcomes plus commands still in flight can
never exceed completion capacity. If user integration places submission and
completion consumption on different tasks, it owns that local coordination;
the HIL side still has one completion producer and does not introduce an MPSC
credit counter.

TCP, UDP, Aeron, iceoryx2, ZeroMQ, wire schemas, client/server roles, connection
lifecycle, and middleware-owned buffers belong to user integration code. That
code may block on external I/O in its own execution context, but it interacts
with the HIL data plane only through the bounded ports. The HIL package carries
no required dependency on those transports and does not define a lowest-
common-denominator transport interface.

Inside one Julia process, consecutive stages with the same owner use direct
calls and actual ownership changes use SPSC rings. `Iceoryx2.jl` is not used as
an internal task scheduler or as fan-out between ordinary optical arms:
atmosphere epoch tokens identify current state, materialized atmospheric path
products and effective-command snapshots have explicit bounded lifetimes, and
each arm owns its propagation workspace. Replacing these relationships with
publish/subscribe would add service discovery, FFI, shared-memory loan,
reference-count, and subscriber-overload semantics without removing an
in-process payload copy.

An optional iceoryx2 adapter is appropriate when a deployment deliberately
crosses a process boundary on one host or needs zero-copy fan-out to several
processes. That adapter belongs to user integration or a transport-specific
extension and maps loaned samples to the same canonical command and acquisition
schemas. Its service configuration must declare maximum publishers,
subscribers, outstanding loans, per-subscriber buffer capacity, history,
overflow/backpressure, dead-node cleanup, and whether a slow subscriber may
gate publication. Optional recorders and telemetry subscribers must not retain
a correctness-critical sample unless slowest-subscriber backpressure is part
of the run contract. `ServiceType.LOCAL` may be benchmarked as a control, but
does not replace the simpler in-process direct/SPSC baseline.

The illustrative names are deliberately product-oriented rather than camera-
or protocol-oriented. UDP datagrams, progressive wire readout, network-byte
order, checksums, transport packetization, first/last-packet timing, and
packet-level faults belong to the user RTC transport contract. User
integration may split or coalesce a complete leased product without promoting
transport fragments into simulation events or canonical HIL descriptors.

The external-delivery contract declares how much time the adapter needs between
complete-product publication and its first required transport observation,
plus the maximum product-lease hold time through final delivery. A camera-like
first-packet deadline is supportable only when the complete product is ready
early enough for the adapter to meet that declared lead time. The simulator
does not claim a progressive-readout deadline by timestamping fragments that
were unavailable at its complete-product boundary. Pool sizing includes the
adapter's worst-case declared hold time and any correctness-critical fan-out.

The warmed data-plane port implementation uses preallocated bounded SPSC
descriptor rings. It does not use `Base.Channel`: `Channel` has blocking
`put!`/`take!` semantics and pays for general lock-and-condition coordination
that this fixed-cardinality boundary neither needs nor permits. Large product
and command payloads live in prepared pools; ring entries contain compact
descriptors, buffer indices or leases, generations, endpoint identity, and
timing metadata.

`try_take!` transfers ownership of a descriptor and its complete product lease
to the consumer. A successful `release!` transfers the lease back to the pool.
Each return path reserves usable capacity for every lease that its consumer can
hold, so a valid first release from the current session cannot return `full`.
A full return path is therefore an invariant failure, not normal backpressure;
the caller retains the lease, the run fails explicitly, and shutdown reports
the accounting deficit. Stale, duplicate, and wrong-session releases return a
structured error without mutating pool ownership. Closing a port stops new
ownership transfer but does not revoke outstanding leases. Stop/fail processing
accounts for every lease before declaring a clean shutdown.

For each pool, preparation proves the accounting invariant

```text
free + producer-owned + completion-queued + consumer-leased + return-queued
    == pool capacity
```

The usable return capacity is at least the maximum consumer-leased term. If
several consumers are permitted, preparation partitions lease credits and
return paths per consumer or introduces an explicitly reviewed multi-producer
reclaimer; it does not silently turn an SPSC return ring into MPSC.

The rings borrow the relevant [LMAX Disruptor](https://lmax-exchange.github.io/disruptor/user-guide/)
and `io_uring` properties without introducing a general event graph or opcode
engine:

- fixed capacity and preallocated descriptor storage
- monotonically increasing producer-published and consumer-released sequences
- one writer for each sequence and no CAS retry loop in the SPSC fast path
- ordinary descriptor writes followed by release publication
- acquire observation before descriptor reads or slot reuse
- cache-line isolation of independently written sequences
- natural batching of already available entries, capped by an item, byte, or
  time budget
- explicit full and empty results rather than waiting in a ring operation

On maintained 64-bit CPU targets, `try_submit!` and `try_take!` should be
wait-free bounded operations after preparation: no locks, conditions, yields,
system calls, allocation, or unbounded retry. Any idle strategy is outside the
ring. Busy polling is permitted only when deployment reserves and pins a
physical core; otherwise the user integration may use a measured backoff or
parking policy without changing the ring contract.

The publication proof is:

```text
producer owns slot
  -> writes descriptor and payload reference with ordinary stores
  -> release-publishes producer sequence
  -> consumer acquire-observes producer sequence
  -> reads descriptor and finishes ownership transfer
  -> release-publishes consumer sequence
  -> producer acquire-observes consumer sequence before reusing the slot
```

Only the shared publication sequences are atomic. The producer keeps its claim
cursor and cached consumer sequence in producer-owned state; the consumer keeps
its read cursor and cached producer sequence in consumer-owned state. A
consumer loads the published sequence once and drains the available natural
batch up to its configured limit. The producer refreshes the remote consumer
sequence only when its cached capacity appears exhausted.

The producer-published and consumer-released sequences must occupy distinct
cache lines and must not share those lines with telemetry, lifecycle flags, or
another ring's cursor. The layout is derived from the maintained target's
cache-line contract rather than assuming struct-field separation is enough.
Preparation or tests must verify the cache-line size, `fieldoffset` values,
containing allocation addresses, and generated atomic code on every maintained
CPU target.
Without that separation, independent producer and consumer writes repeatedly
transfer ownership of one cache line between cores even though the threads do
not mutate the same field: the false-sharing failure the padded layout prevents.
Padding is applied to independently written hot metadata, not indiscriminately
to descriptor slots. Per-buffer generation values detect stale or duplicate
lease release without turning every SPSC slot into an atomic object.

Three sequence domains remain distinct:

| Sequence | Scope and purpose |
|---|---|
| Ring cursor | Internal publication, availability, and slot reclamation; never used as physical event identity |
| Stream sequence | Per acquisition endpoint or command source, paired with a run/session epoch for loss detection and correlation |
| Scheduler ordinal | Prepared deterministic tie-breaker for equal plant timestamps |

A producer assigns the product stream sequence when the physical acquisition
product completes, before attempting ring publication, so a dropped product
leaves an observable gap.
Ring publication order records arrival, not command effective-time order; the
bounded command schedule and canonical plant timestamp remain authoritative.
Independent endpoints do not contend on one global stream counter.

The full Disruptor multi-consumer gating graph is not the default. A slow
recorder or telemetry observer must not gate buffer reuse for a latency-critical
RTC consumer. Noncritical observers receive independent bounded taps with an
explicit loss policy. Shared gating is introduced only when every consumer is
declared correctness-critical and slowest-consumer backpressure is part of the
run contract.

## Capacity And Failure Semantics

Every capacity is justified from the maximum simultaneous burst, maximum
allowed in-flight work, service-rate evidence, and latency budget. Preparation
records the rationale and uses `L = λW` as a consistency check; extra capacity
is not treated as free because it permits additional residence time.

Full behavior is resource-specific:

| Resource | Required full behavior | Ownership/result |
|---|---|---|
| Command submission, payload, or outcome credit unavailable | `try_submit!` returns `full` | Producer retains the command and any payload lease; no ownership transfer or HIL admission occurred |
| Future-effective command schedule | Reject before application or fail according to endpoint policy | Emit an observable terminal command outcome; never backdate or silently discard |
| Required command-completion stream | Each successful submission already owns one bounded outcome credit | Every transferred command produces exactly one terminal outcome; the credit returns when that outcome is consumed |
| Optional command-observation or telemetry tap | Drop or coalesce only under its declared observation policy | Increment a loss counter without gating the plant |
| Acquisition-completion ring | Drop the newly completed product, explicitly coalesce only a latest-value product, or fail according to endpoint policy | Preserve the assigned stream-sequence gap and reclaim the product buffer |
| Product-buffer pool | At the declared lease-acquisition boundary, drop the acquisition or fail according to endpoint policy | Do not start mutating an unowned buffer or allocate a replacement on the warmed path |
| Lease-return ring | Capacity is reserved for every valid outstanding consumer lease; `full` is an invariant failure | Caller retains the lease, the run fails explicitly, and accounting identifies the missing credit; ordinary release never waits or retries |
| Materialized path-product, retained atmosphere-state, or due-work slots | Shed undispatched optional work or fail the required acquisition according to prepared criticality | Never mutate or reuse a slot until every dispatched consumer acknowledges it; an epoch token alone never occupies a retained-state slot |
| Deterministic event slots | Reject during preparation | A fixed periodic schedule cannot overflow while running |

Raw WFS or science frames are not implicitly latest-value products, so
coalescing them requires an explicit endpoint semantic. Overwriting an unread
SPSC slot is forbidden.

Command-submission ports return `closed` when stop/fail begins. Completion and
lease-return ports then remain drainable under a bounded shutdown policy that
states which pending descriptors are delivered, rejected, or reclaimed; they
return `closed` only after that drain ends. The pool owner continues draining
lease returns until every lease is accounted for or the deficit is reported.
Recovery evidence includes occupancy, sequence gaps, outcome counts, and
complete buffer accounting; resuming throughput alone is insufficient.

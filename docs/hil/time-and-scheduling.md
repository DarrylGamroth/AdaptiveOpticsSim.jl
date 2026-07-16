# HIL Time, Scheduling, And Causality

Status: active

## Purpose

This specification defines latency boundaries, timestamp domains, detector
and command schedules, deterministic event ordering, overload behavior, and
replay requirements. It applies to virtual-time correctness runs and
wall-clock-paced HIL deployments.

See the [`HIL architecture index`](../hil-package-boundary.md) for adjacent
subsystem specifications.

## Performance Contract

Every HIL deployment prepares a versioned performance contract before a run.
The contract is part of the run manifest, not a benchmark description written
after results are known. The package supplies measurement surfaces but does
not hard-code one instrument's rates or latency targets.

Modeled physical time and execution performance are separate:

- model time includes exposure, readout, command-device latency, settling, and
  the configured timestamp-domain mapping
- execution time includes scheduling lateness, computation, synchronization,
  bounded queue residence, and explicit copies

A report must not hide execution lateness by folding it into a larger modeled
device delay. The supported latency boundaries are:

| Boundary | Start | End and included work |
|---|---|---|
| Acquisition-product processing | Last scheduled plant event required by the complete product, such as exposure end or readout-ready time | Successful product-descriptor publication or the recorded drop/fail decision; includes due path work, detector finalization, explicit copies, and HIL queue residence, but excludes user transport |
| Product publication lateness | Scheduled product-publication deadline mapped to the monotonic clock | Actual publication or recorded drop/fail decision; early completion and lateness are reported separately |
| Command ingress | Successful ownership transfer into the command-submission ring | Recorded validation and admission decision; terminal completion may occur later |
| Command application lateness | Scheduled effective plant time mapped to the monotonic clock | Publication of the effective optic-state version to optical workers; modeled command delay is reported separately |
| Closed loop | Declared frame integration or sample timestamp | First optical sample whose state includes the correlated RTC command |
| External delivery | HIL port publication or submission boundary selected by the integration | User-declared RTC or transport observation; reported separately from the canonical HIL boundary |

A benchmark or support claim records:

- the selected operation boundary and whether transport, queue residence,
  copies, retries, and response processing are included
- detector periods, phases, exposure windows, optical-sample cadence,
  simultaneous due-event pattern, modeled trigger-source and distribution-link
  topology, and command-arrival process
- nominal and realized detector event times, physical trigger-distribution
  errors, reported-timestamp errors, and every active synchronization-fault
  policy
- source directions, wavelengths, resolutions, detector products, prepared
  product-provider tier, model fidelity, surrogate-validity envelope, and
  controllable-optic configuration
- payload bytes, memory domain, and whether each payload is reused unchanged,
  touched, generated, copied, replayed, transferred, and consumed by the RTC
- adapter-readiness condition, complete-product-to-first-observation lead time,
  maximum lease hold, and whether adapter/telemetry shares the Julia process
- offered, enqueued, admitted, applied, rejected, dropped, coalesced, timed-out,
  and missed-deadline counts plus achieved rates
- ordering, loss, retry, durability, late-event, and recovery semantics
- p50, p90, p99, p99.9, only statistically supported higher percentiles,
  maximum, sample count, repetitions, and raw histogram artifacts
- every queue, schedule, buffer-pool, and in-flight capacity plus measured
  occupancy and headroom
- CPU topology, affinity, FFT/BLAS and Julia thread counts, GPU devices, memory
  domains, interconnect, transport, operating system, and power mode
- source revision and dirty state, runtime and dependency versions, dataset or
  seed, initialization, first-use, warmup, duration, cooldown, and recovery

Fixed-rate detector claims require a schedule-preserving open-loop generator
that continues to account for arrivals while the system stalls. Closed-loop
or one-at-a-time measurements remain useful throughput tests but do not
establish fixed-arrival latency.

A synthetic or reduced-order provider does not weaken this timing contract.
Fixed-arrival latency/capacity and unpaced maximum-throughput runs are separate
benchmark families: the former preserves every offered deadline, while the
latter deliberately removes pacing and reports saturation behavior. Simulator
product-generation, frame publication, user transport, RTC processing,
command return, and command application remain separately observable so a fast
source cannot hide the component actually limiting the loop.

## Time And Causality

### Clock sources

The deterministic core should receive explicit simulation timestamps and remain
independent of wall clock. The HIL companion should inject a clock from
`Clocks.jl` rather than calling a global time function throughout the runtime:

- `MonotonicClock` supplies production interval and deadline time.
- `EpochClock` is suitable when an external wall-clock timestamp is required
  for correlation, but not for measuring a deadline interval.
- `CachedEpochClock` provides an inexpensive shared time value. One clock owner
  calls `Clocks.fetch!`, `Clocks.update!`, or `Clocks.advance!`; workers read it
  with `Clocks.time_nanos` and do not advance it independently.

`CachedEpochClock` is also the deterministic HIL scheduler test clock. Tests can
set a known nanosecond value with `update!` and cross detector or command
boundaries exactly with `advance!`:

```julia
clock = Clocks.CachedEpochClock(Clocks.MonotonicClock())
Clocks.update!(clock, 0)
Clocks.advance!(clock, 500_000) # exactly 500 microseconds
```

For production measurement, the selected clock's resolution, read cost,
cached-update cadence, and maximum observed staleness are part of the
performance contract. A cached epoch value must not be used to claim latency
resolution finer than its measured update behavior; use the direct monotonic
clock for that boundary when necessary.

The production clock, simulation-time origin, wall-to-simulation mapping, and
clock-update ownership must be captured in run metadata. `Clocks.jl` belongs in
the operational HIL companion; the simulation kernel exposes explicit virtual-
time seams so it does not acquire a wall-clock dependency.

The injected `Clocks.jl` source is the **execution clock** used to pace and
measure the simulator. Simulated trigger generation, trigger distribution, and
timestamp-label faults are model events evaluated on the canonical plant
timeline. A detector-trigger phase fault must not be implemented by offsetting
the execution clock: that would conflate a physical acquisition error with
simulator lateness and corrupt the latency boundary. A camera's internal
oscillator may be added by a higher-fidelity device model, but it is not
required to model an externally triggered camera or a faulty trigger fan-out.

### Timestamp domains, schedules, and endpoints

A clock domain is not a rate. It identifies a timestamp coordinate that needs a
defined mapping to the canonical plant timeline. A 2 kHz WFS, 1 kHz DM command
stream, and 100 Hz science detector can all use the same clock domain while
having different schedules or event streams.

The timing model separates these concepts:

| Concept | Meaning | Owner |
|---|---|---|
| Plant timeline | Canonical ordered event time, represented as integer nanoseconds in core | Core receives it explicitly; HIL coordinator advances or samples it |
| Execution clock | Monotonic host time used to pace work and measure deadline performance | HIL companion owns the injected `Clocks.jl` source |
| Modeled trigger source | A periodic or externally scripted sequence of nominal trigger edges on plant time | Core owns deterministic event semantics; HIL prepares the scenario |
| Trigger distribution link | A prepared relationship from one trigger source to a detector, modulator, or other event consumer | Core owns realized delivery times and fault semantics |
| External timestamp domain | RTC, camera, or device timestamp coordinates plus offset/drift mapping into plant time | User integration supplies external timestamps and synchronization observations; HIL owns the configured mapping into plant time |
| Periodic schedule | Period and phase for simulator-paced events such as exposure samples and readout | Core semantics; HIL companion performs wall-clock pacing |
| Command endpoint | Bounded asynchronous command events for one independently timed optic or segment, with receive and effective times | Core owns validation, bounded admission, device state, and virtual-time application; HIL owns ingress, timestamp mapping, pacing, and outcome publication; user integration owns transport and decoding |

Illustrative target configuration is:

```julia
schedules = (
    wfs=PeriodicSchedule(domain=:plant, period_ns=500_000),
    science=PeriodicSchedule(domain=:plant, period_ns=10_000_000),
)

commands = (
    low_order=CommandEndpoint(dm_low; source_domain=:rtc,
        timing=DeviceTiming(latency_ns=200_000)),
    high_order=CommandEndpoint(dm_high; source_domain=:rtc,
        timing=DeviceTiming(latency_ns=50_000)),
)
```

The names are illustrative, not committed API. Both DMs above receive
timestamps from the same RTC domain but accept independent event sequences.
Neither endpoint needs a declared period. If a physical driver only latches on
fixed boundaries, its `DeviceTiming` may include an explicit periodic latch
schedule; that is a device constraint, not a consequence of optical placement.

An external command without a trustworthy source timestamp uses its plant-clock
receive time. A command with a source timestamp must be transformed into plant
time before ordering or latency is computed. Offset, drift, synchronization
uncertainty, and mapping version belong in run metadata when more than one
physical timestamp domain participates.

One command-ingress owner feeds each endpoint. User integration submits
canonical descriptors through a bounded arrival-ordered SPSC port. The HIL
owner maps their timestamp domain; the core validates them and, when necessary,
admits future-effective events into a separate bounded time-ordered schedule.
After application, the effective optic state may be published as a versioned
latest-value snapshot for optical workers. A latest-value slot is therefore an
internal state-publication mechanism, not the RTC command-ingress API.
Capacity, explicit coalescing, rejection, and late-event behavior are part of
the prepared contract.

If several devices genuinely latch one command atomically, represent that as an
explicit multi-optic command transaction or latch group with one effective
timestamp. Do not infer atomicity from a shared optical plane or from a packed
command vector.

### Modeled detector triggers and distribution

One canonical plant timeline remains the causal truth. A prepared trigger
topology may contain a common trigger generator and per-consumer distribution
links. This allows a source phase fault to remain correlated across several
detectors while a cable, fan-out, or receiver error affects only its downstream
endpoint:

```text
canonical plant timeline
└── modeled WFS trigger source (period, phase, source jitter)
    ├── LGS WFS 1 link (delay, skew, edge jitter, trigger faults)
    ├── LGS WFS 2 link (delay, skew, edge jitter, trigger faults)
    └── pyramid branch link
        ├── modulation waveform phase
        └── detector exposure trigger
```

The topology is a prepared, finite acyclic graph. Each trigger source and link
retains only its sequence, fault state, deterministic RNG identity where
needed, and next realized plant deadline; simulating skew or jitter does not
materialize a run-length-sized event list.

Because the simulator owns detector pacing, this fan-out is a deterministic
scheduler relationship, not an RTC command and not one SPSC message per trigger
edge. A ring is introduced only if an actual execution-context ownership
boundary requires one; the serial oracle delivers edges directly.

Trigger-distribution fidelity distinguishes effects that are often incorrectly
collapsed into one timestamp offset:

| Effect | Physical acquisition time | Reported timestamp |
|---|---|---|
| Trigger-source phase or rate error | Moves every downstream trigger edge coherently | Normally follows the realized edge unless a label fault is also configured |
| Trigger-link delay or skew | Moves only the affected detector's delivered edge | Normally follows the realized edge unless a label fault is also configured |
| Trigger jitter | Perturbs individual delivered edges, with common-source and per-link components modeled separately | May report the perturbed edge or an independently faulty label |
| Timestamp-label offset or drift | Does not move the physical exposure | Changes metadata presented to the RTC |
| Phase step, dropped trigger, or duplicate trigger | Alters only not-yet-delivered edges under an explicit fault policy | Records the resulting discontinuity and fault identity |

Nominal trigger time, realized trigger-delivery time, resulting physical
exposure boundaries, and reported timestamp remain separate fields. This
permits an RTC test to distinguish a genuinely misphased detector from a
correctly triggered detector whose metadata is wrong.

A trigger mapping or fault update never retimestamps an edge that has already
been delivered. Its configuration declares whether it applies at the next
source edge or the next delivered edge. If an update would place a
not-yet-emitted edge at or before the current plant time,
the configured policy records a missed edge, emits an explicit duplicate where
that fault is requested, or fails the run; the scheduler does not silently
backdate or manufacture an unbounded catch-up burst. Equal-time updates use a
prepared ordinal, and a trigger update effective at `t` affects a not-yet-emitted
edge delivered at `t`.

Deterministic tests use fixed fault traces or centrally derived RNG streams.
Run metadata records the trigger topology, link and fault versions, seeds or
trace identities, nominal and delivered edges, resulting exposure boundaries,
timestamp labels, and synchronization residuals required for replay.

### Detector-owned pacing

The simulator controls detector acquisition. Detector timing must distinguish:

- exposure start and end
- optical integration sample cadence
- nondestructive sample cadence where applicable
- readout-ready time
- acquisition-completion-port publication time
- first-frame latency and steady-state frame period

Progressive delivery to an RTC is a user transport operation over the complete
leased product, not a sequence of simulation events. Its packet schedule,
first/last-packet timestamps, fragmentation, and transport backpressure are
measured against the `External delivery` boundary and recorded by the adapter's
transport contract.

An acquisition may use its own prepared periodic schedule or consume a
delivered edge from the trigger topology. The detector model defines what that
edge means: exposure start, frame start, reset, row/band start, or another
supported acquisition transition. Exposure duration and readout then follow
the detector's configured timing. A configured trigger-to-transition latency or
jitter belongs to detector response, after link delivery. Camera-internal
oscillator error is optional device fidelity, separate from trigger-source and
distribution-link error.

A delivered edge that reaches a detector while an acquisition is already
active follows that detector's prepared retrigger policy: ignore, restart,
queue only when a physically bounded device model supports it, overlap only
when the detector architecture permits it, or fail. A duplicate trigger never
acquires an implicit unbounded event slot.

Frame period and optical sampling period are not necessarily equal. A slow
science detector may accumulate many evolving optical samples during one
exposure. An up-the-ramp detector may publish several nondestructive samples
from one integration. Rolling shutter requires row or row-band timing. Frame
transfer changes acquisition/readout overlap and timing, while the optical
surface model remains unchanged.

Each physical optical sample is a declared photon-rate product, not a frame
already scaled by a telescope step. The prepared acquisition quadrature owns
the sample duration or weight used to integrate that rate over the half-open
exposure or row-band interval. Its weights account for the interval exactly
once; neither atmosphere advancement nor optical formation contributes another
elapsed-time factor. A normalized synthetic or reduced-order product declares
its explicit conversion before it enters a physical detector pipeline.

Canonical periodic schedules use integer nanoseconds for period and phase.
Preparation creates bounded event generators with one next deadline and stable
tie-break ordinal per active acquisition, command latch, atmosphere evolution,
or other event source. Rolling rows/bands and nondestructive reads advance a
generator cursor; the runtime does not materialize every event in a long run.

Execution jumps directly to the minimum timestamp and processes every due
event; it does not poll every base tick. A prepared allocation-free linear scan
is appropriate for small generator sets, while a preallocated stable array heap
may be selected for larger sets. The chosen policy, threshold evidence, and
generator count are recorded. A repeating due-event cycle may be precomputed
only when its period and storage are explicitly bounded. Times that cannot be
represented at the declared resolution are rejected or quantized under an
explicit fidelity policy recorded in the manifest.

### Autonomous periodic optical devices

Some optical devices execute a continuous local waveform rather than waiting
for one external command per optical sample. A pyramid-WFS modulation mirror is
the principal example: a local controller repeatedly drives a tip/tilt path,
while the detector exposure is normally phase-locked to that motion. The RTC or
supervisory integration may change radius, frequency, phase, dither, mode, or
start/stop state through an ordinary bounded command endpoint, but the HIL
boundary does not carry every waveform point.

The waveform is evaluated from the effective setpoint state and its relationship
to delivered trigger edges at each required optical sample. The trigger
topology expresses whether the modulator produces the detector trigger, the
detector trigger resets the modulation phase, both follow a common source, or a
free-running waveform exposes a phase reference against which detector triggers
are delivered. Cycles per exposure, phase offset, distribution delay/skew,
jitter, and dropped or duplicate edges are explicit run configuration. The
baseline need not simulate the modulation controller's internal servo clock.

The baseline requires one prepared fidelity policy and leaves a second as a
profile-driven extension:

- **cycle averaged:** evaluate a bounded modulation quadrature against one
  frozen optical epoch; this preserves the current fast pyramid model when an
  exposure spans an integral cycle and intra-exposure evolution is neglected
- **time resolved (profile-driven):** evaluate modulation phase and the applicable atmosphere
  and optic snapshots at explicit integration sample times; this supports
  partial cycles, dither, mirror response, and trigger synchronization faults

The fidelity policy and quadrature/sample count are immutable during a run.
Changing a physical setpoint uses the normal effective-time command semantics;
changing fidelity requires another prepare/arm cycle.

### RTC-owned command timing

Controllable optics do not receive a simulator-owned frame rate. The external
RTC sends commands when it chooses. Each command contains or acquires:

- target optic or command segment
- sequence number
- receive timestamp
- optional source timestamp
- modeled effective timestamp after any configured boundary delay and device
  latency

The endpoint's prepared schema supplies the payload shape, element type, units,
command basis and calibration revision, absolute or incremental semantics,
range policy, run/session epoch, and duplicate/reordering policy defined in
[`rtc-ports.md`](rtc-ports.md). These are configuration facts rather than
per-command dynamic dispatch. Shape, schema, session, and sequence validation
occurs before semantic admission.

Each controllable optic holds its last effective command between updates.
Device models may additionally impose minimum update intervals, settling,
bandwidth, slew, hysteresis, saturation, or other physical response. Timing is
generic over `AbstractControllableOptic`; it is not a DM-only feature.

Command silence has an explicit per-endpoint policy. The baseline is indefinite
hold, but a prepared device model may instead apply a preallocated safe/flat
command after a plant-time command-age limit or fail after that limit. This is
a modeled plant transition and is replayable in virtual time. Separately, a HIL
deployment may declare an execution-clock ingress-liveness deadline that fails
the run when the external RTC stops delivering commands; that operational
watchdog does not silently alter optic state. Both clocks, thresholds, safe
state, and reset/recovery rules are recorded. A watchdog transition emits a
bounded fault/state record and cannot allocate, invoke user code, or create an
unbounded command burst.

The modeled policy declares whether its age is measured from valid admission
or effective application; the operational watchdog resets only after semantic
admission, so malformed traffic cannot keep the run alive. At an equal plant
timestamp, an admitted command becoming effective at `t` is applied before a
watchdog expiration at `t` and resets that expiration under the prepared
policy. Stable endpoint ordinals order any remaining simultaneous transitions.

After ingress validation, a command's canonical timestamp and mapping version
are immutable. A later clock-synchronization update must not reorder or
retimestamp an already admitted command. Future-effective commands occupy a
bounded time-ordered schedule. If that schedule is full, the command is
rejected before application and receives an observable outcome.

A command that is already late at ingress follows one configured endpoint
policy: reject it, apply it at the current plant event while recording
lateness, or fail the run. The scheduler never backdates plant state. A
superseding or coalescing policy is legal only when explicitly declared for
that endpoint and every displaced command receives an outcome.

If a command becomes effective during an exposure, optical samples before and
after that event observe the appropriate old and new states. A lower-fidelity
configuration may quantize command application to the detector's optical sample
grid, but that approximation must be explicit in run metadata.

### Critical event flow

Plant state is right-continuous: a command with effective time `t` affects an
optical sample timestamped `t`. Exposure intervals are half-open
`[start, end)`, so a sample exactly at `end` belongs only to another exposure
that includes that timestamp.

Port ingress is asynchronous. Successful `try_submit!` transfers ownership
into the ingress ring but does not itself admit or apply a command. The HIL
owner maps the command timestamp, then invokes the deterministic core's bounded
validation/admission surface. Every transferred command retains one terminal-
outcome credit until the correlated completion is consumed.

For each canonical plant timestamp `t`, the prepared scheduler performs these
phases in order:

1. apply trigger-source, distribution-link, and synchronization-fault updates
   effective at or before `t`, recompute only not-yet-delivered edges, and
   record any explicit dropped or duplicate edge
2. apply all admitted commands with effective time at or before `t`, ordered by
   effective time, endpoint sequence, and scheduler ordinal, then apply any
   still-due modeled command-age watchdog transition at `t`
3. process any due atmosphere-evolution event and select one atmosphere epoch
   for every optical sample due at `t`
4. publish the effective-command snapshot and, when sampled, the immutable
   atmosphere epoch
5. close exposure or row-band intervals ending at `t`, including declared
   trigger-driven transitions, and snapshot nondestructive reads due at `t`;
   these contain accumulated samples strictly before `t`
6. open exposure or row-band intervals beginning at `t`, including declared
   trigger-driven transitions
7. render each unique optical path needed by sample events at `t`, evaluate any
   autonomous waveform at its trigger-relative phase, and accumulate those
   samples only into intervals containing `t`
8. complete due detector-readout events, assign a product stream sequence, and
   attempt bounded complete-product publication
9. publish bounded counters and fault records outside model hot loops

These phases define logical causality, not a mandatory global execution barrier.
The serial oracle executes them directly. A prepared parallel executor may
overlap independent path groups or timestamps only through immutable versioned
snapshots and a bounded number of in-flight epoch slots. It preserves each
acquisition's state dependencies and does not reuse a snapshot until every
dispatched consumer acknowledges it. Unrelated optional science work therefore
need not gate a required WFS path; it is shed or fails according to its prepared
policy when its bounded slots are exhausted.

All events sharing `t` use prepared ordinals for shared state transitions, so
tuple iteration, task completion, and device completion order cannot change
physical results. Publication order across independent completion ports may
vary and is measured rather than treated as plant causality. Wall-clock pacing
belongs to the HIL companion; transport and wire encoding belong to user
integration. The same event sequence remains executable in virtual time for
deterministic validation and faster-than-real-time work.

## Run Lifecycle

The public names are not committed, but every run has these semantic phases:

| Phase | Allowed work | Required exit condition |
|---|---|---|
| Configure | Build immutable parameters and declare capacities, policies, and resources | Configuration validates structurally |
| Prepare | Allocate pools/workspaces, plan FFTs, resolve schedules and placement, warm required code and devices | No unresolved ownership, capacity, or placement requirement |
| Arm | Reset sequences/counters, establish the run epoch, execution-clock mapping, trigger topology, and external clock mappings, verify every lease is returned, and publish initial optic state | Ports and pools are ready, every initial realized event is defined, and user orchestration has reported the selected RTC adapter ready |
| Run | Execute the immutable topology, schedule, capacities, and placement | Stop request, configured terminal event, or explicit failure |
| Stop/fail | Close new submission and admission, record the cause, perform a bounded drain, reclaim leases, and finalize evidence | Every transferred command has one terminal outcome and every owned buffer is accounted for, or the deficit is reported |

Adapter readiness is an orchestration precondition, not a transport API in the
HIL package. User code may implement it with a socket handshake, middleware
discovery, a local barrier, or another bench-specific mechanism; `Run` begins
only after that code reports readiness or the arm deadline fails.

Topology, ring capacity, buffer-pool size, placement, command/acquisition
schemas, and product provider do not change while running. Provider fidelity is
immutable for the initial supported design; changing it requires another
prepare/arm cycle. Runtime fidelity degradation is deferred until a concrete
profile justifies the additional state reservation and claim semantics.

Prepared nonstructural control events may enable or disable an acquisition,
start or stop a trigger source, change a shutter or calibration-source state,
or enter a declared optic safe/hold mode. Every allowed transition, effective-
time rule, inactive-state behavior, and required storage exists before arm.
These events do not add endpoints, mutate schemas or capacities, change
placement or provider fidelity, or execute arbitrary user callbacks. A
structural change requires another configure/prepare/arm cycle. A wall-clock
HIL run cannot pause physical time to recover a missed deadline.

## Overload, Failure, And Recovery

Every bounded resource has a resource-specific full policy; the canonical port
rules are defined in [`rtc-ports.md`](rtc-ports.md). The runtime never converts
overload into an unbounded backlog or silently loses already admitted commands.

At minimum the companion reports:

- detector deadline misses and lateness
- trigger phase/skew/jitter, synchronization residuals, dropped or duplicate
  edges, and timestamp-label errors separately from execution-clock lateness
- stale, duplicate, reordered, clipped, rejected, superseded, and late commands
- plant-time command-silence transitions and execution-time ingress-watchdog
  failures with their distinct clock domains
- frame-port and command-schedule occupancy, headroom, drops, and coalescing
- command age and effective optic-state version at each optical sample
- GPU kernel, synchronization, explicit copy, and host-transfer time
- allocation, GC, context-switch, migration, and configured causal counters
- the selected shed, coalesce, stop, or failure decision

Criticality is prepared per acquisition and command endpoint. Slow optional
science work may be shed under its declared policy; a required WFS or command
path normally fails the run when its contract cannot be maintained. The
runtime does not substitute a cheaper product provider. A lower-fidelity run is
prepared and armed separately so its physical and performance claims remain
unambiguous.

Recovery evidence includes the time to leave overload, bounded occupancy
returning below its declared threshold, all leases being accounted for, no
duplicate application or publication, preserved sequence gaps for dropped
frames, and continued numerical/event-order correctness. Recovery is not
claimed merely because throughput later resumes.

A GPU or worker failure initially fails the run explicitly. Each execution
owner has one preallocated failure slot and acknowledgement state. The owner
catches failures at its boundary, writes a compact structured record, release-
publishes it, and stops accepting new due work at a bounded safe point. A
single coordinator polls owner-specific slots, records the first observed
failure and any concurrent failures, closes new command/acquisition transfer,
publishes the stop epoch, and drives bounded outcome and lease drain. Formatting
of exception details or stack traces occurs after the critical path has stopped.

The run manifest declares the maximum acknowledgement and drain intervals.
Failure to acknowledge or return ownership by those deadlines is itself
reported with the responsible owner, snapshot, command credit, or lease; it is
never converted into a clean stop. Transparent migration requires a prepared
replica with synchronized atmosphere, command, detector, RNG, buffer, and
sequence state plus measured failover behavior; it is outside the initial
support claim.

## Determinism And Replay

The existing deterministic policy remains authoritative for strict reference
generation. Parallel HIL reproducibility additionally requires:

- a central run seed expanded into stable per-layer, path, detector, and device
  streams
- fixed source/path grouping and static device placement recorded in the run
- stable event ordering for simultaneous detector and command events
- timestamps, arrival records, timestamp-mapping observations/coefficients and
  versions, nominal and delivered trigger edges, trigger-link/fault versions,
  reported timestamp labels, admission decisions, terminal outcomes, sequence
  gaps, fault and overload decisions, configuration, and effective command
  history captured for replay
- tolerance-based rather than bitwise CPU/GPU and cross-device comparison

Replay claims use three distinct terms:

- **scenario replay** restores configuration, seeds, scripted external inputs,
  and placement and recomputes the plant; it establishes model determinism
- **boundary-traffic replay** republishes recorded canonical commands or
  acquisition products and metadata without claiming that their physics were
  recomputed; it establishes port, adapter, or RTC behavior
- **decision/event replay** restores actual arrivals, mappings, admission and
  overload decisions, faults, outcomes, and sequence gaps; it diagnoses a
  previous operational run and may intentionally reproduce its missed work

The HIL companion defines the canonical in-memory records and source/sink seams;
user integration owns persistence, retention, and file or middleware codecs.
Every result names the replay class used rather than treating all three as
equivalent.

Deterministic single-threaded CPU execution is the correctness oracle. A
parallel or multi-GPU run should reproduce its modeled outputs within declared
tolerances, not necessarily bit for bit.

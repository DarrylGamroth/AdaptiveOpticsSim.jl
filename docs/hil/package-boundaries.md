# HIL Package Boundaries

Status: active

## Purpose

This specification assigns stable responsibilities to
`AdaptiveOpticsSim.jl`, the future HIL companion, user RTC integration, and
optical companions such as `Proper.jl`. Dependencies point from operational
packages toward the deterministic simulation kernel.

See the [`HIL architecture index`](../hil-package-boundary.md) for adjacent
subsystem specifications.

## Package Boundary

`AdaptiveOpticsSim.jl` remains the deterministic AO simulation kernel. A HIL
companion depends on it, not the reverse.

`AdaptiveOpticsProperHIL.jl` is the existing cross-package proving ground. It
already owns the prepared AO-to-PROPER handoff, a versioned lockstep external-
RTC TCP contract, HdrHistogram latency reporting, and CPU/GPU integration
evidence. Those surfaces remain useful integration evidence, but the TCP
contract is an application-specific adapter rather than an API for the general
HIL package. The target architecture extends the proving ground to detector-
driven multi-rate operation; the existing one-request-at-a-time loop should not
be mistaken for the final asynchronous scheduling model.

### Belongs in `AdaptiveOpticsSim.jl`

- telescope, source, atmosphere, WFS, detector, and controllable-optic physics
- sampled common/path-specific OPD, phase, amplitude, and NCPA elements on
  native optical planes
- detector integration, readout, shutter, thermal, counting, and timing
  semantics in virtual time
- immutable optical-path definitions, acquisition semantics, compatible-result
  keys, and branch-local workspace requirements
- prepared full-optical, reduced-order, and synthetic/replay product-provider
  semantics with one invariant acquisition-product contract
- reduced-order disturbance evolution plus calibrated direction-projection,
  controllable-optic response, and sensor-product operators
- a narrow prepared calibration-illumination extension seam at supported typed
  path entries or detector inputs, without instrument-specific source types
- source-aware atmosphere rendering into caller-owned workspaces
- conjugated-plane geometry and reusable field-propagation seams
- immutable params, mutable state, branch-local workspace traits, and prepared
  mutating hot paths
- prepared serial, CPU, and GPU branch-execution primitives with a deterministic
  fallback
- calibration, prepared control-command routing, detector-metadata, and
  misregistration contracts
- canonical virtual-time plant commands and semantic payload schemas, bounded
  validation/admission, device application, hold behavior, and terminal model
  dispositions
- deterministic trigger-source, distribution-link, and detector-acquisition
  semantics, including physical trigger faults distinct from timestamp labels
- deterministic event stepping independent of wall clock and transport
- device-aware backend identity and backend-generic algorithms
- CPU/CUDA/AMDGPU parity and correctness tests
- simple internal reconstruction and controller primitives used as oracles

### Belongs in the HIL companion

- wall-clock pacing and mapping between monotonic and simulation time
- external RTC, camera, and device timestamp-domain synchronization and mapping
- preparation and wall-clock pacing of the core's modeled trigger topology
- pacing and recording of calibration-source state events when a user scenario
  declares them, without assigning their control authority
- `Clocks.jl` clock selection from the first serial vertical slice, including
  deterministic test clocks and a monotonic production clock; later hardening
  adds cached-clock ownership and external-domain mapping
- canonical command-submission, command-completion, and complete-product
  acquisition-completion ports
- prepared boundary command/acquisition descriptor schemas, command-submission
  descriptors, correlated command outcomes, product leases, endpoint
  identifiers, external timing metadata, and outcome-credit ownership
- configure, prepare, arm, run, and bounded stop/fail lifecycle
- per-endpoint operational ingress-liveness watchdog policy without transport-
  specific health assumptions; modeled plant-time command age remains core
- resource-specific full, close, drain, reclamation, and recovery policies
- cache-line-padded lock-free SPSC descriptor rings and bounded acquisition/
  command payload-buffer pools behind those ports
- long-lived CPU/GPU execution owners, per-owner due/completion paths, and
  bounded slots for materialized atmospheric path products plus model-specific
  retained atmosphere-state snapshots where cross-timestamp rendering requires
  them
- deterministic in-memory integration harnesses, port conformance tests, and
  distinct canonical scenario, boundary-traffic, and decision/event replay
  records with in-memory source/sink seams
- schedule-preserving synthetic load profiles, payload-work declarations, and
  separate fixed-arrival and unpaced-saturation measurement modes
- caller-owned registered-buffer seams and explicit host/device memory-domain
  metadata without depending on a particular transport
- fully explicit or constrained deterministic static placement planning over
  prepared CPU and GPU execution groups; automatic cost optimization remains
  evidence-gated future work
- immutable placement manifests with resource, ownership, handoff, headroom,
  constraint, and rationale metadata
- optional `ThreadPinning.jl` affinity and topology diagnostics
- bounded telemetry taps and a versioned run-manifest schema
- latency histograms, deadline/jitter tracking, warmup, soak tests, and fault
  injection

### Belongs in user or instrument model code

- instrument-specific calibration-source physics, injection placement,
  visibility, combination rules, relay prescription, and source-state model
- selection of the supported typed entry boundary and downstream optical
  segment for each instrument configuration
- calibration-unit profiles, operational modes, and the mapping from
  instrument terminology to generic source/path declarations

### Belongs in user RTC integration

- TCP, UDP, Aeron, iceoryx2, ZeroMQ, shared-memory, and other transports
- wire schemas, framing, endianness, fragmentation, and RTC-specific codecs
- client, server, publisher, subscriber, discovery, reconnect, and middleware
  lifecycle roles
- transport-owned buffer loaning, registration, and native-library lifetime
- mapping transport payloads into canonical HIL command transactions and
  complete product leases
- progressive wire delivery, fragmentation, pacing, first/last-packet timing,
  packet loss, and transport backpressure over leased products
- declaring adapter readiness, complete-product-to-first-observation lead time,
  and maximum product-lease hold time for the selected external contract
- persistence codecs, artifact storage, retention, and replay-log transport
- calibration campaign storage, compatibility policy, selection, and promotion
- mapping an externally controlled calibration unit into its declared generic
  source-state or setpoint surface
- bench-specific setup, teardown, health checks, and transport failure policy
- cross-process or cross-host transport, clock-synchronization observations,
  reconnect, and partition mechanics for any future distributed deployment
- optional iceoryx2 service definitions and loan/fan-out policy when a
  deployment deliberately crosses a process boundary

### Split responsibilities

- Core defines optical-path reuse plus detector timing and incremental
  integration semantics; the HIL companion paces the core's deterministic
  acquisition events against its injected execution clock.
- Core defines detector physics, acquisition timing, and when a complete
  product becomes ready. The HIL companion publishes one bounded product lease;
  user integration owns any progressive presentation, packetization, pacing,
  framing, checksums, acknowledgements, and wire-level failures required by the
  RTC transport contract.
- Core defines prepared product-provider behavior and produces the same typed
  acquisition product at each fidelity tier. The HIL companion paces and
  publishes it without interpreting optical fidelity, records provider and
  payload-work policy, and limits performance claims to that declaration.
- Core executes a user-declared calibration source at a supported typed path or
  detector-input boundary. User or companion model code defines its physical
  integration; the HIL companion only paces declared state events, and user
  integration maps any external control protocol.
- Core accepts canonical virtual timestamps and defines periodic schedule
  semantics; the HIL companion maps an injected `Clocks.jl` clock and any
  external timestamp domains onto that plant timeline.
- Core maps nominal trigger edges through prepared distribution links into
  physical acquisition events and keeps reported labels separate. The HIL
  companion paces those realized plant events with its execution clock; it
  does not inject a detector-trigger fault by corrupting the execution clock.
- Core defines the semantic plant command schema, controllable-optic state,
  virtual-time validation/admission, response, application, and terminal model
  disposition. The HIL companion maps and transfers its own boundary submission
  descriptors into that surface and publishes correlated command outcomes;
  user code transports and decodes them. Core does not depend on descriptor,
  lease, port, or outcome-credit types from the HIL package.
- Core models sampled device feedback as an ordinary scalar, vector, or image
  acquisition. The HIL companion publishes that acquisition independently of
  command outcomes; user integration maps it to the RTC's feedback transport.
- The HIL companion owns the bounded submission/completion port contract and
  buffer lifetime; user integration owns everything beyond that boundary.
- In-process stages use direct calls or SPSC ownership handoffs. An optional
  iceoryx2 adapter belongs at a deliberate process boundary and does not replace
  path-local workspace ownership, current-epoch materialization, or bounded
  retained-state ownership.
- Core defines backend/device identity, supported execution seams, and prepared
  group requirements; the HIL companion resolves the static placement policy
  and owns CPU and GPU execution agents.
- Core exposes prepared worker ownership; the HIL companion may apply and
  verify `ThreadPinning.jl` affinity without making pinning a core side effect.
- Core defines calibration data meaning. The HIL manifest records the stable
  calibration identity and compatibility metadata used by a run; user
  integration or an optional artifact companion owns serialization, storage,
  selection, and campaign lifecycle. Existing core `cache_path` and
  `write_config_toml` conveniences are transitional and do not become part of
  the plant or HIL API; their replacements return structured data to an
  extension or caller-owned persistence policy.
- Core owns AO-plant optics and native sampled aberrations; `Proper.jl` owns a
  selected detailed relay, coronagraph, or instrument prescription. The HIL
  companion schedules a generic prepared external-optics executor; a
  Proper-specific adapter stays in a sibling package or extension, so neither
  core nor the general HIL package acquires a required `Proper.jl` dependency.

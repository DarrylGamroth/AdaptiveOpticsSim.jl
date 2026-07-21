# HIL Validation And Acceptance

Status: active

## Purpose

This specification defines the evidence needed to expand HIL support claims.
The deterministic serial CPU path remains the correctness oracle; parallel
and device backends add equivalence and timing evidence without redefining
the physical model.

See the [`HIL architecture index`](../hil-package-boundary.md) for adjacent
subsystem specifications.

## Validation And Acceptance

Every capability gate keeps a serial correctness oracle. Required functional
evidence includes:

- a clean Gate 2 serial-plant artifact with raw mergeable histograms, exact
  workload/seed/environment metadata, deterministic declaration-order replay,
  one shared atmosphere advanced to explicit model time, distinct science/NGS/
  LGS path materialization, detector fan-out, and a zero warmed allocation gate;
  this is self-paced service-time evidence, not external-RTC latency or
  fixed-arrival capacity evidence
- one atmosphere selection/advance to explicit model time per sampled plant
  timestamp with identical epoch visibility across all due directions and no
  telescope-cadence fallback
- rejection of stale `AtmosphereEpoch` tokens, serial materialization of every
  due path input before the next advance, and parallel lifetime tests proving
  that mutable layers are never advanced under an unmaterialized reader;
  plans that retain older state additionally prove bounded slot accounting
- analytic wind-offset checks over zero, unequal, and multi-step durations
- independent caller-owned wavefronts/fields sharing one telescope definition
  without shared path OPD, focal-plane result, propagation-workspace, or
  spatial-filter mutation
- prepared plane compatibility checks covering physical sampling,
  centering/orientation, wavelength/channel, units/normalization,
  coherence/combination policy, backend, and physical device
- source-to-optical-to-detector dimensional accounting with non-unit unequal
  durations proving optical outputs remain photon-arrival-rate products and
  acquisition applies elapsed time exactly once
- two detector exposures consuming one immutable photon-arrival-rate product
  without optical recomputation
- spatial-density and cell-integrated-rate fixtures establishing photon
  conservation through the declared presampling response and pixel-integration order
- compatible spectral grids preserving flux under declared incoherent
  accumulation and incompatible grids remaining bundled or failing preparation
  instead of being silently summed by array index
- heterogeneous WFS paths and independent direct-science and PROPER-coronagraph
  paths at different cadences
- on-axis, off-axis NGS, and finite-height LGS direction geometry
- mixed NGS/LGS path scheduling against one atmosphere epoch
- detector integration, rolling shutter, frame-transfer timing, and actual
  event-time up-the-ramp nondestructive reads where supported; the existing
  post-exposure synthesized ramp is tested and claimed separately
- presampling detector response and charge-coupling stages remain in their declared order when
  acquisition is split across scheduled events
- multiple schedules with different periods and phases on one plant timeline
- fixed-capacity event scheduling that visits only due timestamps, preserves
  declaration-order-independent simultaneous-event order, rejects stale and
  foreign claims plus backward rescheduling, retains constant storage over a
  long run, and allocates zero warmed heap bytes across the maintained linear-
  scan generator range
- common-source detector triggers with fixed per-link phase/skew, correlated
  source jitter, independent link jitter, and deterministic phase-step faults
- dropped and duplicate detector triggers without hidden event catch-up or
  ambiguous exposure ownership
- detector-specific retrigger behavior while an exposure or readout is active,
  with bounded ownership and an observable outcome
- trigger-source error, link-delivery error, and detector
  trigger-to-transition latency attributed to the correct stage
- physical trigger-delivery error distinguished from timestamp-label error and
  execution-clock lateness
- cycle-averaged pyramid modulation with explicit detector-trigger phase;
  time-resolved equivalence and misphasing cases only for profiles that claim
  the time-resolved capability
- commands arriving before, on, and between detector sample boundaries
- half-open exposure intervals and stable command/sample/readout ordering when
  several events share one timestamp
- late, future-effective, superseded, and rejected commands without backdating
  plant state
- schema, run-epoch, payload-shape, basis/calibration, duplicate, reordered, and
  range-policy command validation
- plant-time command-silence hold/safe/fail behavior separated from an
  execution-clock RTC-ingress watchdog
- independent command updates for multiple controllable optics
- independent endpoint event rates for co-conjugated optics sharing one RTC
  timestamp domain
- explicit atomic multi-optic transactions without accidental atomicity for
  ordinary packed or co-located commands
- independent command timing and additive OPD for co-conjugated DMs
- common MCAO command visibility and isolation of target-specific MOAO commands
- geometrically conjugated optic footprints at multiple altitudes
- common aberrations visible to all selected paths and NCPA visible only to its
  selected branch
- native sampled-NCPA and detailed `Proper.jl` science-path handoff conventions
- exact execution-clock pacing tests using `Clocks.update!` and
  `Clocks.advance!` on a `CachedEpochClock` in the HIL scheduler
- deterministic external-domain offset and drift mapping into canonical plant
  time
- deterministic trigger-source and distribution-link fault record/replay,
  including nominal edges, delivered edges, resulting exposure boundaries, and
  reported timestamps
- command enqueue, admission, application, and terminal-completion reporting,
  including full, rejected, stale, superseded, and deadline-missed outcomes
- sampled actuator/device feedback delivered as an acquisition product rather
  than a command outcome
- SPSC full, empty, wraparound, burst, sequence-gap, and product-lease
  generation behavior under concurrent producer/consumer stress
- resource-specific full policies and clean stop/fail accounting for every
  transferred command outcome and outstanding product lease
- guaranteed lease-return credit for every valid outstanding consumer lease,
  including injected invariant failures and complete pool-state accounting
- worker/GPU first-failure publication, stop acknowledgement, bounded drain,
  and explicit reporting of stranded snapshots, outcomes, or leases
- prepared acquisition/trigger/shutter/calibration/safe-state transitions that
  do not mutate topology, schema, capacity, placement, or provider fidelity
- serial versus CPU-parallel and CPU versus GPU parity within tolerance
- stable RNG-owner identities and derivation versioning, including endpoint
  reorder, changed execution-group order, and applicable CPU/GPU placement
  cases that preserve each physical event's random domain
- constrained deterministic and fully explicit placement resolving to
  reproducible, inspectable prepared plans
- rejection of conflicting, unsupported, memory-infeasible, and
  burst-overloaded placement requests before a run starts
- mixed CPU/GPU epoch, command, detector, RNG, and numerical consistency
- single-GPU versus multi-GPU epoch, command, and frame-sequence consistency
- identical acquisition-boundary contracts across full-optical, reduced-order,
  and synthetic/replay providers
- command-responsive reduced-order behavior within a declared validation
  envelope and explicit rejection of optical claims from static/replay sources
- user-defined calibration illumination at a supported typed path or detector
  input with explicit visibility, timing, combination semantics, deterministic
  replay, and no role-triggered propagation bypass
- zero warmed steady-state heap allocation on maintained CPU HIL paths
- same-process adapter/telemetry allocation and GC soak evidence or explicit
  isolation of allocating integration in another process
- bounded preparation/compilation latency and generated-code growth as path and
  endpoint registries scale, preventing whole-instrument topology from becoming
  one unbounded recursively specialized type

## Minimal Vertical-Slice Evidence

The first RTC-facing evidence uses one serial CPU owner, one scheduled
acquisition, one command-responsive optic, an injected clock, and an in-memory
fake RTC. It exercises the core plant command schema, HIL submission descriptor
and command/outcome pair, complete-product lease, adapter-readiness precondition,
and release accounting without worker queues, GPU submission, placement
planning, or a real transport.

The fake RTC must perform meaningful work: consume the declared pixel, slope,
or other sensor product; calculate a command through a pinned reference
controller; return that command; and reduce the selected reduced-order residual
over a deterministic disturbance sequence. Its schedule-preserving arrival
generator continues to count offered frames while either side stalls. The same
boundary also runs unpaced to diagnose saturation, but that result is not used
as fixed-arrival evidence.

Acceptance records complete-product publication time, declared adapter lead and
maximum lease-hold time, simulated RTC processing, command enqueue/admission,
effective application, and first causally affected optical sample separately.
Every outcome and lease is accounted for at normal stop and injected failure.
This gate proves the product seam early; it does not claim production ring
layout, parallel execution, GPU residency, or external-transport latency.

## Ring And Port Evidence

Ring acceptance requires:

- layout evidence that independently written cursors cannot share a cache line
- generated-code evidence that maintained targets implement `UInt64`
  publication without locks or CAS retry loops
- a documented release/acquire happens-before argument and concurrent stress
  tests for full, empty, wraparound, stale leases, and sequence gaps
- warmed latency and throughput measurements with producer and consumer on
  separate physical cores, including burst, stalled-consumer, saturation,
  recovery, and relevant NUMA placements
- proof that every resource-specific full policy preserves ownership and
  produces the declared outcome
- proof that return capacity covers every valid outstanding lease and that a
  valid first release never follows an ordinary retry path

An unpadded layout may be retained only as a benchmark control, not a supported
data-plane implementation.

## Performance Evidence Protocol

Every timing run begins from the versioned contract in
[`time-and-scheduling.md`](time-and-scheduling.md). Correctness checks run
before timing, and initialization, compilation, first use, warmup, measurement,
cooldown, and recovery are distinct phases.

The maintained load curve covers:

| Region | Required observation |
|---|---|
| Idle/light | Service time, wake-up behavior, and instrumentation overhead |
| Target | Absolute deadline and percentile contract at the declared arrival schedule |
| Burst | Simultaneous due events, natural batching, and capacity headroom |
| Near saturation | Tail growth and first limiting resource |
| Saturation/overload | Rejection, dropping, bounded occupancy, failure, and stop behavior; provider fidelity remains unchanged |
| Recovery | Time and correctness after pressure subsides |
| Long soak | GC, leaks, clock drift, sequence continuity, and thermal/power phases |
| Cold/first use | Compilation, planning, page faults, device/context startup, and connection setup where operationally relevant |

Fixed-rate claims use a schedule-preserving open-loop generator. Generator
timestamps, achieved offered rate, and missed arrivals are retained as a time
series so coordinated omission can be audited. Histogram correction, when
reported, is additional evidence and does not replace the open-loop result.

`HdrHistogram.jl` histograms are configured before measurement with explicit
unit, range, significant figures, and counter width. Prefer one histogram per
single writer and merge compatible interval histograms outside the timed path.
Construction, resize, queries, encoding, and export stay outside the operation
boundary unless production performs them there. Clock-read and recording
overhead are measured separately; out-of-range values fail the benchmark.

Each result archives raw mergeable histograms, interval boundaries, offered and
completion rates, outcome counts, queue/pool occupancy, allocations and GC,
selected causal counters, exact commands, environment/topology, affinity,
runtime/dependency versions, source revision and dirty state, seeds, duration,
and independent repetitions. A percentile is reported only when sample count
supports it.

Performance promotion requires both an absolute gate tied to the product
contract and a relative gate against a comparable baseline with repeated-run
dispersion. Controlled hardware carries latency gates; shared CI is limited to
correctness, allocation, benchmark integrity, and coarse regression checks.

Replay evidence names one of three noninterchangeable claims: scenario replay
recomputes from configuration/seeds/scripted inputs; boundary-traffic replay
reissues recorded canonical products or commands; decision/event replay
restores recorded arrivals, mappings, admission/overload decisions, outcomes,
faults, and sequence gaps. Persistence remains user-owned. Tests must reject a
report that presents boundary-traffic replay as model determinism or uses a
recomputed scenario to claim reproduction of the original operational
overload decisions.

## Low-Fidelity RTC Load Profiles

The maintained fast path has two distinct purposes and must label them:

| Profile | Required source behavior | Valid evidence | Invalid inference |
|---|---|---|---|
| Synthetic boundary load | Produce the canonical product shape, type, cadence, metadata, leases, and bounded outcomes from preallocated patterns or replay | Port, transport, and RTC latency/throughput; burst, saturation, overload, and recovery | Optical accuracy, control stability, or command-dependent plant response |
| Reduced-order closed loop | Evolve time-correlated disturbances and produce command-responsive geometric, linear/modal, low-rank, or reduced-resolution products inside a validated envelope | RTC reconstruction, tomography, control, loop stability, latency/throughput, and fault recovery inside that envelope | Full optical performance or behavior outside the validated envelope |
| Full optical reference | Execute the maintained physical models | Model and optical claims covered by the applicable validity evidence | Unsupported model/backend combinations |

Boundary-conformance tests run the same acquisition and RTC adapter against all
applicable providers and compare descriptor schema, shape, element type,
sequence and timestamp semantics, leases, and overload outcomes. Reduced-order
providers additionally compare against a full optical or analytic oracle over
a pinned command and disturbance envelope. Constant and replayed products must
be marked nonresponsive; they cannot close a physically meaningful feedback
loop merely because commands are accepted.

Reduced-order closed-loop acceptance additionally requires:

- an open-loop disturbance sequence with declared temporal statistics and
  deterministic replay
- calibrated path-projection, controllable-optic response, and sensor operators
  checked independently for sign, units, dimensions, and path visibility
- exact command causality at exposure boundaries, including modeled device
  latency, hold, clipping, and multiple independently timed optics
- a matched controller case that reduces a declared residual metric and remains
  stable over the declared operating envelope
- deliberately wrong-sign, delayed, dropped/stale-command, and calibration-
  mismatch cases that produce the expected degraded or unstable response
- multi-direction and multi-rate cases when reconstruction or tomography is
  part of the profile claim
- reduced raw-pixel products, when offered, checked before and after the
  selected detector/readout pipeline against their declared oracle

The run manifest identifies every disturbance sequence or state model,
projection and response operator, calibration revision, noise model, omitted
effect, validity range, and comparison tolerance. This makes a reduced-order
claim reproducible without presenting it as full optical equivalence.

Every synthetic result records product bytes, buffer/corpus count, working-set
size, reuse distance, memory domain, payload policy, and whether the user
integration and RTC actually touched or copied the payload. A metadata-only
test may claim descriptor rate, not raw pixel throughput. A pixel-processing
claim uses production-sized payloads and representative cache, memory, and
transport work even when their values are simple.

Each maintained load profile runs both:

1. a schedule-preserving fixed-arrival test at target rates, phases, common
   bursts, and declared faults, retaining offered arrivals even while the RTC
   or simulator is stalled; and
2. an explicitly unpaced saturation test for maximum useful throughput.

The report separates simulator generation, frame publication, external
delivery, RTC processing, command return, and command application latency. It
also records offered, completed, rejected, dropped, and deadline-missed counts,
ring and pool occupancy, generator lateness, allocations, and repetitions. A
fast source is acceptable only with measured headroom under its declared
fixed-arrival contract; omitting arrivals after a late completion is
coordinated omission, not recovery.

## Support Promotion

No new CPU, CUDA, AMDGPU, mixed CPU/GPU, multi-GPU, or multi-host surface
becomes production-supported until its exact model/placement combination is in
a maintained real-hardware validation cadence and the support matrix is
updated.

Initial detector promotion centers representative CMOS, CCD, frame-transfer
EMCCD, and HgCdTe/SAPHIRA-style acquisition paths. Other scalar, vector, area,
or photon-counting products use the same generic acquisition and port contracts
but are promoted only with their own model/backend evidence.

Core support of the calibration-illumination seam does not validate every
instrument calibration unit. Each native or user-provided source/entry/backend
combination states its own physical oracle, tolerance, deterministic replay,
allocation budget, and path-visibility evidence. A complex external optical
executor retains the validation responsibility of its owning package.

## Instrument Reference Profile Rules

Instrument profiles are maintained scenarios in a companion validation or
configuration package. They are not hard-coded modes in the general core or HIL
package, and they do not imply current production support. Each manifest pins
the public reference revision or access date because instrument specifications
may evolve.

A profile claim has four independent levels:

1. **Topology compliance:** the scenario represents the declared telescope,
   atmosphere, sensing directions, acquisitions, and independently placed
   controllable optics.
2. **Model compliance:** the selected pupil, atmosphere, guide-star geometry,
   WFS, detector, corrector, misregistration, and relay products match pinned
   references within declared tolerances.
3. **Timing compliance:** nominal and faulted trigger schedules plus command and
   closed-loop latency meet the declared contract without hidden backlog.
4. **Integration compliance:** a user transport adapter closes the loop with an
   external RTC through canonical ports and preserves timestamps, ordering,
   buffer lifetime, and overload semantics.

Reduced spatial or temporal profiles may validate topology and causality, but
must be labeled reduced-scale. Full-rate compliance requires a declared
hardware placement, fixed-arrival tests, tail-latency histograms, deadline miss
counts, saturation and recovery tests, and a long-duration soak. A single-host
or single-GPU result must not be generalized to a multi-device deployment
without maintained evidence.

Every instrument profile should provide a production-shaped synthetic traffic
variant before requiring full optical propagation, and a command-responsive
reduced-order variant when closed-loop behavior is part of the claim. These
variants exercise the same topology, triggers, ports, command endpoints, and
fault policies as the full model; their results retain the claim limits above.

## NFIRAOS Synchronized Multi-Rate Reference Profile

NFIRAOS is the synchronized multi-rate and trigger-distribution profile. Its
public topology combines several WFS families and instrument-provided sensors
with a modulated pyramid branch, making it a strong HIL target for both nominal
phase alignment and detector-trigger distribution faults.

The reference is pinned to the TMT
[adaptive-optics overview](https://www.tmt.org/page/instruments-adaptive-optics)
and the
[NFIRAOS pyramid modulation study](https://arxiv.org/abs/2108.06429).

| Dimension | Public reference scale | Architectural stress |
|---|---:|---|
| LGS sensing | Six LGS WFS paths; NFIRAOS operates at up to 800 Hz | Simultaneous directional WFS work and correlated common-source triggers |
| NGS and instrument sensing | One high-order NGS pyramid WFS, up to three OIWFS, and/or up to four ODGW depending on mode | Heterogeneous products, rates, exposure semantics, and trigger links |
| Pyramid modulation | Fast steering mirror synchronized to individual PWFS exposures, with a quarter-rate dither and phase relationships to other active WFS channels | Autonomous waveform, trigger-relative phase, optical-gain dither, and deliberate synchronization faults |
| Correctors | Two DMs conjugated near 0 km and 11.8 km; the ground DM is mounted on a tip/tilt stage | MCAO geometry, independent command response, and combined high/low-order motion |
| Operating modes | LGS MCAO and classical NGS AO, with the PWFS role and rate changing by mode | One generic plant API must support mode-specific active paths without instrument-specific core types |

The maintained profile should include at least a nominal synchronized run, one
fixed-skew detector link, common-source phase and jitter faults, an independent
link-jitter case, dropped and duplicate edges, and a timestamp-label-only fault.
Each result records nominal trigger edges, delivered edges, resulting exposure
boundaries, modulation phase, and reported detector timestamps. A trigger fault
must change only its declared downstream consumers, while a common-source fault
remains correlated across every attached branch.

The detailed camera electronics remain out of scope unless a separately
validated device model requires them. The profile tests the externally visible
trigger and acquisition behavior that an RTC must handle.

## MORFEO Extreme-Scale Reference Profile

MORFEO is the extreme-scale end-to-end profile. It exercises the same generic
plant, path, timing, port, and placement APIs used by smaller systems.

The public reference scale is taken from ESO's
[MORFEO overview](https://elt.eso.org/instrument/MORFEO/) and the
[2023 RTC status presentation](https://www.eso.org/sci/meetings/2023/RTC4AO/01_12_foppiani.pdf).
The profile manifest pins the source revision or access date.

| Dimension | Public reference scale | Architectural stress |
|---|---:|---|
| LGS sensing | Six 68 x 68-subaperture paths, 1100 x 1100 pixels at 500 fps | Simultaneous high-rate direction rendering, WFS execution, detector completion, and external delivery |
| NGS sensing | Three probes with distinct low-order and reference-camera products and rates | Shared path segments, heterogeneous detectors, and independent schedules |
| Correctors | ELT M4 plus two post-focal DMs | Three independently commanded, geometrically conjugated surfaces |
| Low-order correction | M5 and other low-order plant response as required by the selected profile | Independent command endpoints and device-response timing |
| Science relay | AO residual supplied to a detailed optical prescription where required | A slow PROPER-backed path must not delay WFS deadlines |

The selected MORFEO boundary must declare whether the external RTC receives
raw pixels, calibrated pixels, slopes, or another product. The six public LGS
dimensions imply about 7.26 GB/s of raw 16-bit pixels before NGS products,
metadata, copies, or transport framing. This derived raw-pixel load does not
apply to a slope boundary and is not itself an acceptance target; `GB/s` here
uses decimal units.

The maintained profile distinguishes a camera-boundary run from a subsystem-
boundary run. A camera-boundary run presents six independent LGS acquisitions
plus the selected NGS reference and low-order acquisitions; the external RTC
still owns pixel calibration, slope extraction, reconstruction, vector
combination, tomography, control, and command splitting. A subsystem-boundary
run may instead inject prepared high-order or low-order vectors between RTC
processing hosts. That is useful for component HIL, but it does not validate
camera ingress or the RTC pixel pipeline and must not be reported as the full
end-to-end profile.

A raw-pixel profile publishes one complete product lease at the canonical HIL
boundary. When the RTC contract uses progressive readout, user integration
segments and paces that lease into wire datagrams and owns byte order,
checksums, packet loss, transport backpressure, and first/last-packet timing.
A first-pixel-to-command result is therefore an external-delivery claim whose
manifest includes the adapter, its complete-product-to-first-packet lead time,
maximum lease hold, and transport contract; it is not a canonical HIL port-
publication claim. If the complete product cannot be ready before the required
first-packet deadline minus that lead time, the camera-like contract is
infeasible at this boundary rather than rescued by fabricated progressive
simulation events.

The command side similarly declares the actual bench boundary: independent
post-focal DMs, the selected M4/M5 or control-system-facing transaction, and
any LGS steering or feedback endpoints required by that configuration. A
packed RTC message may map in user integration to several canonical commands
or to one explicit atomic latch group; transport packing does not imply that
the physical optics share a cadence or command state. RTC-internal high-order,
vector-sum, reconstruction, and closed-loop pipes remain outside the simulated
plant when that RTC is the system under test.

Early implementation configurations are architectural inputs, not acceptance
references. The companion profile pins its endpoint inventory, payloads,
rates, trigger relationships, command mappings, and RTC/configuration revision
before evidence is promoted.

A pinned SPECULA MORFEO scenario can serve as a comparative correctness oracle
for public topology and numerical products. It does not replace an ESO
instrument reference, current interface-control document, or hardware timing
measurement.

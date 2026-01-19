# "Telescope in a Box" RTC Harness (Design Sketch)

This document outlines an RTC test harness that treats the simulator as a
telescope-in-a-box: an external controller sends DM commands and receives
sensor telemetry in real time (or accelerated time) for validation.

## Goals
- Provide a stable interface for external RTC controllers.
- Support deterministic playback and reproducible integration tests.
- Allow real-time or accelerated-time execution.

## Non-goals
- Hard real-time guarantees (Julia user-space limitations).
- Full fidelity network protocol specification (left to project needs).

## Core loop (conceptual)
```julia
while running
    cmd = read_command!(adapter, t_now)
    dm.state.coefs .= cmd
    apply!(dm, tel, DMReplace())
    advance!(atm, tel)
    propagate!(ws, ngs, atm, tel, dm, wfs)
    telemetry = pack_telemetry(wfs, tel, t_now)
    write_telemetry!(adapter, telemetry)
    wait_for_next_tick(clock, dt)
end
```

## Interface sketch
```julia
abstract type AbstractRTCAdapter end
read_command!(adapter::AbstractRTCAdapter, t) = # DM command
write_telemetry!(adapter::AbstractRTCAdapter, data) = # WFS slopes/frames

abstract type AbstractClock end
now(::AbstractClock) = # current time
wait_for_next_tick(::AbstractClock, dt) = # pacing
```

### Adapter options
- In-process: direct function calls for unit tests.
- Shared memory: ring buffer, shared arrays, or memory-mapped files.
- Network: TCP/UDP/ZeroMQ/gRPC (choose one for a concrete implementation).

### Telemetry payloads
- Slopes (preferred for RTC loop tests).
- Full WFS frames (for debugging and algorithm validation).
- Optional metadata: timestamps, jitter metrics, config hash.

## Deterministic test mode
Combine with deterministic simulation:
- Fixed RNG seed and fixed environment.
- Fixed clock that advances by `dt` without wall-clock time.
- Full record/replay of DM commands and telemetry.

## Suggested milestones
1) In-process adapter + fixed-step clock (unit tests).
2) Shared-memory adapter for low-latency local RTC.
3) Optional network adapter for external controllers.

## Notes
This harness should live in a `RTC/` or `Sim/RTC/` module and avoid tight
coupling to a specific controller stack.

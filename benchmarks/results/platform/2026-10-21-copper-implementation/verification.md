# Independent verification of the Copper CUDA implementation profile

Disposition: the archived timing reductions, device-activity arithmetic, source provenance and principal report claims are supported by the saved evidence. No GPU workload or production-source modification was performed during this review. No AOS optimization has been implemented by this experiment.

## Integrity and unchanged implementation

- Verified all 22 raw-file SHA256 values in `provenance.json`, with exact equality between the manifest's file names and the raw directory.
- Verified the archived baseline harness against Git commit `041cb2a2d3e725c9049d7abd39456a304bc6e83e` and the current harness against its recorded SHA256. The only subsequent harness change assigns the optional `CUDA.@profile` return object and displays it; unprofiled timing paths are unchanged.
- All three baseline records identify that same clean AOS commit and harness, with 32 warmups and 128 samples per execution mode.
- `git diff 7e6ea87 -- src ext` is empty, including the worktree. Numerical source tree `e6c136c01fb1359a09800cdabbef0164008a87de` and extension tree `940270f3b65e2c685c2fcce50f87909a574f785a` equal the original revision's trees. Thus the data describe the existing implementation rather than an implemented optimization.

## Independent timing reduction

Recomputed every baseline process median from its 128 raw samples and then recomputed the median and range of process medians:

| Execution | Median of process medians | Observed host allocation samples |
| --- | ---: | ---: |
| Stream | 2.717452 ms | 12,320 bytes/frame |
| Captured | 2.628863 ms | 0 bytes/frame |

These agree with `summary.json` and README rounding. They measure completed `step_graph!` calls with evolving native atmosphere/RNG, zero commands and the noisy detector. They do not measure the complete HIL/RTC path or establish a deadline success probability.

Every serialized node device-event median and the µs conversion agrees with its raw samples. The local CUDA API documents `elapsed(start, stop)` in seconds; the harness converts to milliseconds, and the reducer converts milliseconds to microseconds correctly.

The start/stop event interval excludes the subsequent host wait on the stop event. It can include device idle gaps while the host submits work, and per-node synchronization changes scheduling. The separate host-complete measurement includes that wait. For example, median stream PWFS host-complete attribution is about 7.717 ms while its device interval is 2.549 ms. These diagnostic values must not be substituted for the uninterrupted completed-frame measurement or summed as native stage timings. A minor README wording refinement was communicated to distinguish the interval itself from the synchronization procedure.

## Independent native-trace reduction

Reparsed all 300 device-activity rows directly from the compressed human-readable trace, including 144 captured-mode rows. Every stored start label, duration conversion, launch attribute and kernel name agrees with `activity.json`. Recomputed all captured kernel counts and summed durations divided by four frames; they agree with `summary.json`.

- Captured activity duration sum: 2,639.095635 µs/frame, or 10.556383 ms over four frames, consistent with the profiler's rounded 10.56 ms busy total.
- Combined cuFFT duration: 1,644.015 µs/frame. The large inverse batch contributes 1,593.350 µs/frame, about 96.92% of the recorded FFT time.
- Shifted-mask batch fill: 428.7375 µs/frame.
- Modulation intensity accumulation: 395.18 µs/frame.
- Final binning: 41.01 µs/frame.

The 97.80% value is the profiler's GPU-busy fraction within its captured four-frame trace window. It is not achieved occupancy, application-wide utilization, peak computational efficiency, or evidence of fixed-rate operation. No occupancy or hardware memory-traffic counter claim is supported or made.

The intensity operation's stated array traffic is 256 MiB of complex stack reads plus 4 MiB of previous intensity reads and 4 MiB of intensity writes. Dividing 264 MiB by 395.18 µs gives 700.50 GB/s in decimal units, agreeing with the report's approximate 700 GB/s. This is algorithmic array traffic, not measured DRAM bandwidth; cache behavior and actual hardware transactions were not measured.

## Interpretation

The evidence supports prioritizing fixed-discretization PWFS execution and preparation-memory ownership investigations. It does not establish an optimal batch size, a measured benefit from fusion, or an implemented speedup. The dense modulation-cube sizes are correctly derived: 1440² × 256 × 8 bytes = 3.955 GiB, and 480² × 32 × 8 bytes = 56.25 MiB. The source review separately identifies the retained device cube and temporary preparation host cube; neither establishes a per-frame transfer or a memory leak.

The README correctly separates profiled diagnostic activity from unprofiled completed-frame latency, identifies the zero-command workload, and records shared-host/uncontrolled-clock limitations. No numerical-accuracy conclusion follows from this implementation profile; the atmospheric convergence campaign provides a separate numerical comparison.

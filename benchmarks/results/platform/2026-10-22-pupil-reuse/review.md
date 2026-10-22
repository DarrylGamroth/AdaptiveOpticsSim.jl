# Independent review: physical pupil reuse

Date: 2026-10-22. Reviewed checkout: `/home/dgamroth/workspaces/codex/AdaptiveOpticsSim-pupil-profile`, HEAD `27179140bb3c55e1d50e82bb5759578eb54fcb98`, with the uncommitted implementation and validation changes. No GPU workloads or production edits were performed by this reviewer.

## Disposition

Accept the reviewed implementation for the finite numerical, capture, calibration and closed-loop cases recorded here. No unresolved production correctness defect was found. This acceptance does not establish ensemble statistical equivalence, deployed RTC qualification, hardware timing guarantees, or instrument fidelity.

The production diff SHA-256, independently rechecked after the final validation, is:

```text
37418e59bb2dea6443eec09b145c66bd5f0aa61a29c2eecbcf96f1fa1c1f9cc5
```

It matches the production changes reviewed earlier and the hash recorded in the results README. The timing table and final trace paragraph were still pending when this review artifact was written; acceptance of the existing README does not extend to unseen additions. Final packaging and archive checks remain with the coordinating agent.

## Source review

The separable representation uses cis(vₚ c[i]) × cis(uₚ c[j]); it retains the original axis convention, quadrature offsets, amplitude weights and physical tilt before each forward FFT. CPU and Bi-O-edge retain dense phase preparation. Accelerator Pyramid fallback paths form the same physical pupil from the factors. No source prescription, geometry, feedback gain or detector/RNG implementation changed.

The common pupil phase is regenerated from current OPD on each propagation, including the capture enqueue path and each source/wavelength invocation. Both factors receive exact-device validation. The pupil-phase cache receives exact-device validation and participates in the prepared workspace binding. Circular updates mutate existing factor arrays and preserve point count, weights and captured pointers. The two-dimensional mask kernel visits every active batch plane. No cache-lifetime, ordering, stale-pointer or fallback defect was found.

Mathematical factorization changes floating-point evaluation: cis(x + y) and cis(x)cis(y) need not produce identical bits. Finite-difference and closed-loop evidence is therefore material to acceptance; image agreement alone would be insufficient.

## Findings and disposition

- **PRR-001 — P2, high confidence, resolved.** `test/backend_optional_common.jl` still accessed accelerator Pyramid `.modulation.phases`, which the separable owner does not expose. The final assertion reconstructs the two factors with a nonconjugating transpose and retains the original tolerances.
- **PRR-002 — P2, high confidence, resolved.** The first integration comparator unconditionally read `host_allocation_bytes` for static runs, although the harness records that field only for moving runs. The final comparator explicitly permits the missing field only for static traces. Static allocation coverage is not claimed in the README.
- **PRR-003 — validation coverage, resolved.** Initial focused tests omitted noncentered accelerator propagation and Pyramid ElectricField parity, and capture replay did not poison the new cache or change OPD. Added tests cover those cases. The reviewer inspected the added tests and their retained passing log: 581 storage/update, 28 noncentered/field, and 38 capture/update/cache checks.
- **PRR-004 — interpretation, adjudicated.** Seeded noisy GPI trajectories exceed the newly proposed 0.05 nm wavefront-vector threshold. The comparator preserves this result as false; the README explicitly reports it. The vector threshold remains applicable to deterministic trajectories. No production thresholds were changed to force a pass.

## Independent CPU verification

The final `compare_pupil_reuse_validation.compare` and `summarize_pupil_reuse.summarize` functions were independently rerun without writing output files, using CPU affinity 4 and one BLAS thread. Both completed successfully and verified the recorded array hashes. The comparator includes the full COPPER response matrix and frozen-calibration loop, GPI default/fine response matrices, and default/fine/transfer trajectories. The README's numerical and storage figures agree with these records.

The largest deterministic moving GPI vector differences are 0.000050521 nm default, 0.000046282 nm fine, and 0.000045073 nm transfer. The corresponding scalar WFE differences are smaller. The retained independent scientific replay record reports successful source/hash checks, 35 retained modes, reconstruction and control-update checks, and finite unsaturated trajectories. That scientific auditor was run by the coordinating agent; this reviewer inspected its result and independently recomputed the cross-version Gram comparisons.

## Noisy-run adjudication

| Case | Maximum vector difference | Vector RMS difference | Maximum scalar WFE difference | Vector RMS / realized noise effect |
|---|---:|---:|---:|---:|
| Default | 0.07610958 nm | 0.04309717 nm | 0.01958627 nm | 0.6388% |
| Fine | 0.08220171 nm | 0.04056075 nm | 0.02647528 nm | 0.6032% |

The default and fine tail scalar WFE changes are approximately +0.002404 nm and +0.000730 nm. These paired-seed trajectories remain bounded and unsaturated. Scalar WFE agreement alone does not bound differences in wavefront direction, which is why the vector results remain explicitly reported.

Observed: rate/calibration rounding changes coexist with small deterministic trajectory differences and larger noisy trajectory differences. Source establishes that the Poisson sampler is discontinuous in its mean: its inversion branch compares a random product against exp(−λ), and its large-mean approximation rounds to an integer. Random keys depend on seed, draw sequence, domain and pixel, not λ. Derived: equal keys do not require equal detected outputs when λ changes. Hypothesis: a particular threshold crossing initiates the observed noisy divergence. The first divergent detector event was not instrumented, so its precise cause is not established. The default vector first exceeds 0.001 nm at frame 61; that is an observation threshold, not an identified causal event.

Frame differences are input-equivalent electrons, not ADU, because these fixtures have no ADU calibration. Same-version reset/replay qualification is distinct from cross-version noisy trajectory agreement. One paired seed cannot establish equality of stochastic distributions. Acceptance is limited to the recorded stability exercise, supported by deterministic/calibration evidence and unchanged detector/RNG code.

## Benchmark and reporting review

`step_graph!` waits for completion, and captured execution synchronizes its execution context. The measured boundary is therefore a completed optical node, with setup, transfers, detector and RTC excluded as declared. Warmup and allocation sampling precede timing. Saved output files establish cross-representation comparisons; each invocation's equality check establishes only replay determinism.

Dense, separable and computed diagnostic owners use the same new mask/cache kernels, weights, sampling and batch policy. These experiments isolate representation and must not be called the original production baseline. Computed and dense outputs are byte-identical in all three recorded fixtures; separable relative L2 differences are below 6.80 × 10⁻⁸. The COPPER caller benchmark now measures actual CUDA array fields and records their shapes/bytes; its introspection change does not alter optical execution.

Direct storage accounting supports the capacity claim: at N1440/P256, dense phases occupy 4,246,732,800 bytes; factors occupy 5,898,240 bytes and the added pupil-phase cache 16,588,800 bytes. These figures are distinct from allocator-pool snapshots and total VRAM peaks.

Initial small-case measurements visibly drifted. The long-warmup Copper fast observations are dense 2.683420 ms, separable 2.684272 ms, and computed 2.714508 ms. Two fresh fine-case runs show GPI 319.033/319.103 ms before versus 316.340/316.422 ms after, and Copper 456.587/456.548 ms before versus 451.459/452.236 ms after. These are descriptive service-time observations under uncontrolled clocks and shared-host activity, not robust universal speed rankings or tail-latency guarantees. No measured DRAM-transaction or bandwidth claim follows from logical loads or storage sizes.

The existing README correctly limits its claims, preserves the noisy provisional failure, distinguishes recalibrated GPI from frozen-calibration COPPER, and distinguishes scientific finite-case validation from instrument qualification. Archive the exact dirty patch, included benchmark scripts, commands, environments and raw results: HEAD plus the primary harness hash alone does not identify an uncommitted storage experiment.

## Final performance-section addendum

The subsequently added Completed-node performance section was reviewed against `performance-summary.json`, `stats-after-copper_fine.log`, and a read-only SQLite query of `after-copper_fine.sqlite`. Table medians, repeat pairs and long-warmup values agree with the records. The paired fine reductions are approximately 0.84% for GPI and 0.94–1.12% for COPPER, as stated. The report retains the uncontrolled-clock/shared-host limitation and distinguishes original production comparisons from representation-only diagnostics. These measured repeat pairs support a modest observed fine-case improvement, not a universal or statistically qualified performance guarantee.

The SQLite export independently confirms 3,592 kernel instances, 898 distinct graph-node IDs, and four memsets. No CUPTI memcpy activity table is present. Kernel summary sums divided by four agree with FFT 307.289 ms, mask 70.391 ms, tilted-pupil fill 38.076 ms, intensity 41.313 ms and binning 0.571 ms per replay. The added pupil-phase kernel averages 0.039504 ms per replay. FFT time is approximately 67.1% of total recorded kernel time. These are kernel-time sums from a profiled execution, not the unprofiled host-visible service time.

The documented CUDA-event collection warning remains relevant: absence of a memcpy table establishes no recorded copies, and should be described that way rather than as proof of completeness of every activity class. The final paragraph appropriately does not infer achieved occupancy or DRAM bandwidth and does not isolate the contribution of each production change. The original 42.771 ms fill value is treated as the preceding profiling campaign's supplied baseline; the reviewer independently checked the new trace, not a fresh baseline trace execution.

The pending performance-section review is now complete. Production diff hash remains `37418e59bb2dea6443eec09b145c66bd5f0aa61a29c2eecbcf96f1fa1c1f9cc5`. Final packaging remains outside this review. Acceptance remains limited to the finite recorded cases and makes no claim of ensemble statistical equivalence.

# Independent modulation batch review

Reviewed worktree: `/home/dgamroth/workspaces/codex/AdaptiveOpticsSim-batches`, branch `codex/modulation-fft-batches`, base `e7b4fbb4`. Review date: 2026-10-22. The worktree contained the intended uncommitted changes and added harnesses. No production files were changed by the reviewer, and no CUDA work was launched.

## Conclusion

The production change is small and preparation-only. No scientific correctness, ownership, capture, allocation, or ordering defect was identified by source inspection. The identified benchmark defect is corrected, and the cardinality concern is resolved by exact measured-case gates. Final timing, numerical, capture, instrument integration, and physical pupil-tilt external-oracle evidence support acceptance of the batching change. There are no unresolved blocking findings for this change. Broad performance guarantees are not supported: the repeated measurements support the tested shapes, sample counts, hardware, and runtime. The external comparison also exposes substantial pre-existing shifted-mask differences from physical pupil tilt; this change preserves those differences and does not establish shifted-mask physical fidelity or deployed-controller accuracy.

## Findings

### BATCH-001 — P2, confirmed: policy-comparison benchmark overwrites the case selector

- Confidence: high.
- Evidence: `benchmarks/sweep_modulation_batches_cuda.jl`, `sweep_main`, originally lines 97 and 108–109. `selected` initially contains case names. `selected,_ = measure_batch(...)` replaces it with a dictionary when `AOS_BATCH_COMPARE_POLICY=1`. The next outer iteration evaluates `c.name in selected`, which fails because Julia dictionary membership requires a pair.
- Reproduction: Julia 1.13.2, `selected=Dict("batch"=>2); "b" in selected` raises `AbstractDict collections only contain Pairs` without any GPU activity.
- Impact: the advertised production-policy comparison mode aborts after the first selected case, even for a single requested case when later case definitions remain.
- Remediation: use a separate name such as `selected_record` for the measured record.
- Validation: inspect the resulting variable lifetimes; run a multiple-case policy comparison when the serialized GPU campaign permits.
- Disposition: resolved in source. The measured result is now named `selected_record`; the original case selector remains unchanged. Inspection confirms correct variable lifetimes. The completed 34-case policy-comparison campaign confirms the repaired control flow.

### BATCH-002 — P3, evidence limit: cardinality and platform generalization

- Confidence: high for the evidence gap; a regression is a hypothesis, not an observed defect.
- Evidence: production chooses caps by FFT side, precision, compute capability, and cache size, independent of point count. Repeated data cover side 256/count 32; side 272/counts 24, 30, 31, 32; side 408/count 32; side 512/counts 12, 32, 96; and side 528/count 64. It consequently applies to unmeasured composite counts, other qualified GPU models, and future cuFFT/runtime versions.
- Analysis: the divisor rule preserves numerical sampling and memory bounds. Performance depends on the selected FFT batch and number of launches. The prime-count fallback specifically prevents the measured count-31 collapse, but does not establish performance of every small composite count. A future count-dependent regression is plausible; none was demonstrated.
- Remediation: clearly distinguish measured qualification from extrapolation. If a strict no-regression promise is needed, restrict the table to measured counts or add discriminating boundary runs. Avoid claiming globally optimal selection.
- Validation: representative small composite counts and both propagation strategies on any newly supported GPU/runtime.
- Disposition: cardinality concern resolved: production now explicitly gates the measured `(side, point_count)` combinations listed above, while all other counts retain the baseline. Platform/runtime generalization remains a documented measurement limitation.

### BATCH-003 — P3, evidence presentation: compare repeated baseline trials

- Confidence: high.
- Evidence: sweep records include both the initial `baseline` measurement and repeated trials of the baseline batch. At side 408, shifted-mask baseline batch round medians were 0.465908, 0.482524, 0.482009 ms; batch 2 was 0.408195, 0.462512, 0.461800 ms. The initial baseline was 0.4668 ms.
- Analysis: comparing the initial baseline with the later selected-batch median gives 1.011×; using repeated medians gives 1.044×. Clock/thermal or other phase variation is visible, and no causal counter establishes its origin. Repeated within-campaign baseline measurements are the stronger comparison.
- Remediation: report repeated trial medians and dispersion, preserve raw samples, and treat marginal gains conservatively. Do not infer that L2 capacity alone is the cause.
- Disposition: addressed by two final campaigns with shuffled baseline/selected pairs. Retain the raw distributions and measured-workload scope in reporting.

## Source and evidence checks

- The generic helper still chooses a divisor of the existing point count and applies the existing 256 MiB field-stack budget. It rejects invalid counts, plane sizes, and the newly exposed internal cap. CPU preparation is unchanged; nonqualified CUDA and other accelerators retain cap 32.
- CUDA queries happen during workspace preparation. The selected size is stored in the existing workspace; no device query or tuning is introduced into execution or capture.
- The existing propagation loops clear the output before accumulation, execute all points once in their established order, and fully populate each batch before its FFT. Regrouping changes floating-point summation association but introduces no change to coordinates, weights, modulation points, masks, or normalization.
- The timing boundary is credible: `step_graph!` calls `wait_graph_step!`, and captured execution waits for device completion. Samples are unpaced completed service times, not fixed-rate response-time or full instrument latency claims.
- The harness warms graphs, excludes preparation/capture/output transfer from timing, checks output agreement and zero host allocations, randomizes candidate order, records raw samples, and repeats three rounds. Forty-eight samples per round support median comparisons; they do not establish rare latency tails.
- Repeated baseline/selected median ratios are approximately 1.24×/1.26× at side 272/count 32, 1.39×/1.42× at side 256/count 32, 1.14×/1.33× at side 512/count 32, and 1.20×/1.36× at side 512/count 96 (shifted/pupil respectively). The broader cache heuristic was appropriately rejected after its measured regression at side 560. Side 576 was already a fallback under that heuristic; its forced batch-2 measurements do not demonstrate a policy regression.
- The standalone numerical tests cover Float32/Float64, both propagation strategies, nonseparable mask geometry, weighted paths including zero weights, photon conservation, CPU-derived response matrices, a frozen CPU reconstructor, feedback trajectories, poisoned scratch/output capture replay, output pointer stability, and zero replay allocations. The observed log reports 685 passes; the reviewer did not rerun GPU tests.
- The external oracle harness compares flat/ramp camera frames and a finite-difference tip response against OOPAO, with explicit sample counts and physical units. The final additions include the actual tuned side-256/count-32 Float32 policy and distinguish physical pupil tilt from the shifted-mask approximation. Complete external closed-loop equivalence remains outside the evidence.

## Final evidence update

- Re-inspected the final two-file production diff. Exact FFT-size/point-count gates match the stated table; all other counts and precisions retain baseline selection.
- `policy-timing-1.toml` is complete: 34 case/strategy combinations, three randomized baseline/selected pairs per combination, 64 samples per measurement, and recorded L2 = 5,242,880 bytes. Every measurement reports zero warmed host allocations. Changed-policy relative L2 rate errors are at most 7.57 × 10⁻⁸; unchanged-policy cases have exact outputs. The repaired comparison mode completes across cases.
- In that confirmation campaign, all changed combinations improved median service time: baseline/selected ratios span approximately 1.115–1.470×. Unchanged-policy ratios are approximately 0.999–1.024×, consistent with ordinary timing variability. These are completed photon-rate-node service measurements, not full instrument fixed-rate latency guarantees.
- Independently reran `benchmarks/compare_modulation_batch_validation.py` with the oracle virtual environment against the archived binary arrays and their hashes. It passed and produced `integration-comparison-review.json`.
- COPPER N=240/M=96: photon rates relative L2 = 8.389 × 10⁻⁸; complete interaction matrix relative L2 = 5.024 × 10⁻⁶; worst column relative error = 1.158 × 10⁻⁴. The same frozen calibration drove both 320-frame loops: maximum command difference = 1.137 × 10⁻¹³ m; maximum residual wavefront-error difference = 2.752 × 10⁻¹⁵ m. The archived noisy atmosphere frames are exact.
- GPI N=240/M=24: all 17 saved arrays are byte-identical, including calibration, saved coefficient histories, and final noiseless/noisy frames. This covers the archived 320-frame illustrative 35-mode controller run; unsaved intermediate images and a deployed/full-order controller are outside the evidence.
- Inspected the completed logs: preparation-policy checks 378/378; prior compact-storage and in-place modulation-update regressions 629/629. Together with the 685 numerical/capture checks and 1,068 CPU checks, these cover both the changed selector and adjacent execution ownership.
- Inspected Float64 external CUDA oracle logs for three geometries, including the side-256/count-32 geometry. Flat/ramp images and finite-difference tip responses agree with OOPAO at relative L2 below 10⁻⁷ for default and forced batches 1/2/4. OOPAO revision, configured precision, and relevant source hashes are recorded. The Float32 production-table oracle and second final timing campaign subsequently completed as recorded below.

- `policy-timing-2.toml` independently completes the same 34 combinations and three randomized pairs each, again with zero warmed host allocation and maximum relative rate difference 7.57 × 10⁻⁸. Changed-policy speedup ratios range from 1.115 to 1.470×; unchanged-policy ratios range from 0.9978 to 1.0026×. This confirms approximately 10.3–32.0% service-time reduction for the measured changed cases.
- Independently read the Float32 external-oracle parameter files: padded side 256, 32 modulation points, baseline batch 32, selected batch 8, for both strategies. Independently recomputed the following metrics from saved arrays with NumPy:

| Strategy / quantity | Baseline versus OOPAO relative L2 | Selected versus OOPAO relative L2 | Selected versus baseline relative L2 |
|---|---:|---:|---:|
| Physical pupil tilt / flat image | 8.0442 × 10⁻⁶ | 8.0457 × 10⁻⁶ | 5.9008 × 10⁻⁸ |
| Physical pupil tilt / ramp image | 7.9320 × 10⁻⁶ | 7.9320 × 10⁻⁶ | 6.9746 × 10⁻⁸ |
| Physical pupil tilt / tip derivative | 7.7658 × 10⁻⁵ | 7.5617 × 10⁻⁵ | 7.1108 × 10⁻⁵ |
| Shifted mask / flat image | 0.02301435 | 0.02301435 | 6.2268 × 10⁻⁸ |
| Shifted mask / ramp image | 0.04287013 | 0.04287013 | 6.9352 × 10⁻⁸ |
| Shifted mask / tip derivative | 0.27845717 | 0.27845794 | 8.6033 × 10⁻⁵ |

The external physical pupil-tilt comparison supports the tuned Float32 batching implementation. The shifted-mask discrepancies are observed in both baseline and selected implementations: roughly 2.30% flat-image error, 4.29% ramp-image error, and 27.85% tip-derivative error on this fixture. They are pre-existing model differences, not a newly demonstrated batching regression. Their precise cause was not established by this review. Do not present the internal frozen-controller parity checks as external physical qualification of shifted-mask operation, full instrument fidelity, or deployed RTC accuracy. Addressing that scientific approximation requires a separate investigation; no speculative algorithm change is justified by this batching task.

Final disposition: accept the narrowly qualified preparation optimization with the explicit measurement and model-fidelity limits above. No production changes or GPU jobs were performed by the reviewer; independent verification consisted of source review, CPU reproduction of the benchmark selector failure, parsing final raw timing/test artifacts, rerunning the binary integration comparison, and recomputing external-oracle errors from saved arrays.

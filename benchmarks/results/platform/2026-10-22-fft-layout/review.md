# Independent FFT implementation review — 2026-10-22

## Scope and disposition

Read-only review of the production changes in `src/wfs/pyramid/{setup,optics,stages}.jl`, `src/wfs/exact_target.jl`, and `ext/AdaptiveOpticsSimCUDAExt.jl`, plus changed tests, benchmark drivers, numerical comparators, and the investigation report. No production edits, GPU calls, tests, git-index changes, or commits were made by this reviewer. Independent CPU reductions ran with `taskset -c 4`, one OpenBLAS thread, and bytecode writing disabled.

Baseline: `80cc4fefd7dee6b3f5a8697935bbb89dfa0af5d2`. Worktree: `/home/dgamroth/workspaces/codex/AdaptiveOpticsSim-fft-layout`. Both the current five-file production diff and archived `production.patch` have SHA-256 `8ef1c3abaf072b243926a2328728c8b8462438987d5206a136a22c8c4be1fed2`.

No unresolved production correctness defect was found. The full GPI default/fine/transfer instrument comparison is now complete, and its active acceptance gates pass together with the COPPER gates. The final disposition is accepted for the measured configuration, subject to the noisy-realization and performance limits below. The fine noisy trajectory does exceed the separately reported provisional 0.05 nm wavefront-vector check; this is preserved as a limitation, not relabeled as a passing check. A keyword-only argument error in the GPI harness was corrected before the successful restarted run; this was not a production failure.

## Numerical and ownership assessment

**Derived from source:** The paired path computes the first forward axis, transposes the two spatial axes without changing modulation order, computes the second forward axis, multiplies by `mask[j,i]`, executes the inverse axis, transposes back, and executes the remaining inverse axis. This equals the original two-dimensional Fourier/mask/inverse pipeline. Both inverse axes remain unnormalized; the unchanged intensity factor N⁻⁴ supplies the correct normalization. Weighting, phase construction, physical sampling, modulation count, centering, and final intensity accumulation are preserved.

The transpose uses the correct Julia leading-dimension padding `(33,32)`, an unconditional barrier, and matching partial-tile input/output bounds. Valid destination elements read initialized shared entries. Plane indexing is independent. Synchronous execution completes the transpose phases; capture enqueues them on the existing stream.

`PyramidTransposeFFTWorkspace` owns the additional array and both FFT handles with concrete field types. Prepared bindings retain all three; exact-device checking includes the new array. Resizing creates fresh arrays/plans while preserving the concrete family and existing batch count, consistent with the enclosing mutable workspace's fixed field type. Fresh preparation reapplies policy. Prepared bindings/revision checks invalidate an old binding after resize.

Selection is limited to square ComplexF32 planes on sm86 with 5 MiB L2. Both qualified fine grids select B16; only side 2880 with 192 points selects the paired layout. Side 3072 with 256 points remains native. Neighboring sizes, point counts, precision, and device classes retain the default path. These predicates define measured qualification, not a claim of hardware-independent optimality.

## Findings and dispositions

| ID | Severity / confidence | Evidence and required action | Disposition |
|---|---|---|---|
| FFT-R1 | Medium / high | The sweep's within-process reference alone cannot establish native-layout equivalence. `summarize_fft_layout.py` now compares and hash-checks saved outputs against native results. Full fine-grid calibration/loop evidence is additionally required. | Resolved: rate comparisons and full COPPER/GPI active gates pass. Noisy statistical equivalence is not claimed. |
| FFT-R2 | Medium / high | Main-stack bytes omit auxiliary transpose arrays and cuFFT work areas. Current reporting separates explicit extra arrays and sampled pool use; the report correctly labels 1536 MiB as a main-stack bound. | Resolved. Total peak VRAM is not claimed. |
| FFT-R3 | Informational / high | Installed cuFFT.jl `get_batch_dims` and external-batch execution select dimension 1 internally for axis-2 transforms, then iterate over modulation planes. Axis decomposition makes 2 + 2B cuFFT calls per tile; it is not equivalent batching to native. | Correctly documented; no production axis-only policy adopted. |
| FFT-R4 | Low / high | Initial experimental records lacked complete transitive harness hashes. Current sweeps include them; archived initial sources and the frozen production patch separate early experiments from final validation. | Final production identity verified. Early experiments must retain their provenance qualification. |
| FFT-R5 | Medium / high | The existing N240/P96 COPPER driver did not exercise the fine production choice. The new driver runs N1440/P256, full interaction calibration, and the frozen 320-frame loop. | Resolved. GPI preflight also confirms actual N720/P192, B16, and the paired workspace before the completed instrument comparison. |
| FFT-R6 | Medium / high | Initial benchmark shared tile `(32,33)` padded the wrong dimension for Julia and retained bank conflicts. It was changed to `(33,32)` before tiled GPU measurements. | Resolved; 24 standalone transpose checks and production kernel tests passed. |
| FFT-R7 | Medium / high | The fourth budget argument initially caused the benchmark's historical baseline to inherit the new larger budget. Zero now explicitly selects the original 256 MiB budget as well as original limit. | Resolved; six focused regression checks passed. |
| FFT-R8 | Low / high | The report initially stated all kernel spans were over 99.8% busy; the COPPER baseline trace gives 99.735487%. | Corrected to over 99.7%; timing conclusions unaffected. |

An intermediate concern that production paired selection could bypass experimental FFT overrides was checked against the latest source and cleared: `fft_layout_override.jl` explicitly forces native workspace construction before applying its experiment wrappers.

## Independently reproduced evidence

Re-executed the current comparison/reduction functions against retained raw files, without GPU work or new output files. `copper-parity.json`, `layout-summary.json`, and `production-summary.json` reproduce exactly, including their array SHA-256 checks.

- Maximum cross-layout rate relative L2: 1.773250639715582 × 10⁻⁷, below the existing 2 × 10⁻⁶ gate.
- Production completed-node reductions: GPI 3.2633728346%; COPPER 4.1633166629%, from two original/modified process pairs with reversed second-pair order.
- Fine COPPER rate relative L2: 9.1286003408 × 10⁻⁸; maximum probe flux relative error: 6.2525202793 × 10⁻⁹.
- Fine COPPER interaction-matrix relative L2: 5.3897999139 × 10⁻⁶; maximum column relative error: 1.2079629471 × 10⁻⁴, below unchanged 2 × 10⁻⁵ and 5 × 10⁻⁴ gates respectively.
- Frozen COPPER loop command relative L2: 1.4661463914 × 10⁻⁷; maximum WFE difference: 3.2617443883 × 10⁻¹⁵ m; maximum Strehl difference: 3.1512105481 × 10⁻⁸.

The comparator refactor preserves the pre-existing gates. It records additional slope/atmosphere diagnostics without implying every recorded quantity has a separate assertion. GPI's nonorthogonal basis Gram matrix supplies the physical wavefront-vector metric. Noisy runs preserve/report the provisional vector threshold separately and explicitly do not establish statistical equivalence from one seeded realization; the final report must preserve that distinction.

### Completed GPI acceptance review

Independently reran `compare_fft_validation.compare` on CPU core 4 with one OpenBLAS thread. Its output exactly reproduces `integration-comparison.json`, including raw-array hashes and all active assertions. Independently reran the scientific auditor with all original numerical assertions intact; only its final file-write and JSON-print statements were omitted to keep this review read-only. The resulting report exactly reproduces `scientific-replay.json`.

The scientific audit checks current harness/source hashes, 40 raw-array digests, the Julia oracle basis hash, independently reconstructed basis/Gram values, singular values and retained rank, reconstructor identities, linear-loop spectral radius, command sequencing, controller update parity, saturation, physical command amplitude, common-grid wavefront metrics, and timing-summary arithmetic. Its independently reconstructed basis differs from the oracle by at most 7.5661787114 × 10⁻¹⁷; the Gram maximum absolute discrepancy is 5.1869619710 × 10⁻¹³. It does not establish a newly independent optical model.

- GPI default interaction, command arrays, final frames, and compared wavefront metrics are identical to the retained baseline for static, moving noiseless, and moving noisy runs.
- Fine interaction relative L2 is 7.5524187108 × 10⁻⁷; maximum column relative error is 3.2360904712 × 10⁻⁶.
- Fine moving-noiseless maximum wavefront-vector difference is 9.4757811476 × 10⁻¹⁴ m; maximum scalar WFE difference is 2.3533209441 × 10⁻¹⁴ m.
- Transfer moving-noiseless maximum wavefront-vector difference is 3.5006536335 × 10⁻¹⁴ m.
- Fine moving-noisy maximum wavefront-vector difference is 0.0994908751 nm, exceeding the preserved provisional 0.05 nm vector check. Maximum scalar WFE difference is 0.0204498985 nm, below its active 0.05 nm gate. Vector-difference RMS is 0.0418626066 nm against 6.7242697842 nm RMS for the realized noise effect. One paired seeded run establishes neither statistical equivalence nor the cause of individual stochastic trajectory differences.

Verified the retained GPI baseline identity: the six archived production-file hashes equal committed `80cc4fef`, and the historical patch hash equals `git diff 27179140..80cc4fef -- src ext`. Thus dirty-parent metadata does not imply a different optical baseline. The GPI baseline is from an earlier campaign, so its whole-loop timing comparisons are observational; the fresh reversed-order node pairs remain the strongest performance comparison. The final production diff still has the frozen hash above. The final separate CUDA policy log repeats all 676 passes.

Normal-sampling control metadata reports unchanged output hashes and zero warmed allocations for COPPER N240/P96 and N480/P32, with current medians 2.8582165 ms and 4.993373 ms. Those earlier-campaign timing comparisons are sanity checks, not an interleaved causal measurement.

Read and summed passing test logs: focused CPU 1597; KA CPU 1038; CUDA policy 676; forced paired storage/ownership 964; native storage 956; modulation batch/capture 457; captured storage 33; complete GPI capture 350; standalone tiled transpose 24; baseline-budget regression 6. These are observed results from the coordinator's runs, not tests executed by this reviewer. The paired suite includes asymmetric face shifts, both precisions, eager/captured live updates, ElectricField and uncentered fallback, target checks, and resize ownership.

## Performance claims and limits

Trace attribution and unprofiled timing are kept separate. GPI's reduced FFT-kernel time must be considered together with the added transpose cost; COPPER's FFT kernels themselves are slightly slower while mask/intensity stages improve. Kernel activity fractions do not establish occupancy, bandwidth efficiency, or a roofline classification.

The report's 16 N²P/B mask/intensity traffic expression is a source-level logical access count. Its derived savings and the paired path's additional approximately 50.961 GB per GPI frame are not measured DRAM transactions. Counter access remained restricted; no counter, clock, power, or permission changes were made. The proposed explanation involving locality and cuFFT algorithm choice remains a hypothesis. Larger batches retain the same optical sample volume.

The mask-locality microprobe is outside the frozen production patch. Its isolated-kernel result cannot establish a full-pipeline improvement because input restoration changes cache state; any approximately 0.8 ms/frame extrapolation must remain explicitly hypothetical. The coordinator rejected promotion without full-path evidence.

Hardware/software scope is an RTX 3080 10 GiB, the recorded Julia/CUDA/library versions, uncontrolled clocks on a shared host, unpaced completed service times, and fixed inputs. These results establish neither hard real-time deadlines nor instrument fidelity, hardware validation, statistical noisy equivalence, or performance on a lower-bandwidth replacement GPU.

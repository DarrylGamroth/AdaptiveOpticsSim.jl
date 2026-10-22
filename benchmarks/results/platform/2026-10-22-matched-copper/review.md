# Independent review: matched Copper 240/96 comparison

Review date: 2026-10-22. Reviewer task: numerical prescription and benchmark protocol review before viewing any optical results. This review applies the benchmark-low-latency-systems skill. It does not establish measured agreement or performance.

## Reviewed state and scope

AOS worktree: `/home/dgamroth/workspaces/codex/AdaptiveOpticsSim-matched-copper`, HEAD `b590aa26a6f6f976ff0f821b138b43c1a84bea77`. Existing changes at inspection: untracked `benchmarks/matched_copper/`; no tracked modifications. Reviewed `fixture.jl`, `oopao_adapter.py`, `specula_adapter.py`, `sources.json`, native source, and `specula-alignment.md`. Drivers and result evaluator were not yet available at this review. Only this evidence document was written by the reviewer; no production, adapter, index, or commit changes were made. No CPU propagation or GPU experiment was run.

External revisions: OOPAO `168fad2a8d98721a270c1b066b813b716680b5ea`; SPECULA `a06bb495ae29a6ec160628069598b32631b9f44f`. Reviewed local clones under the evidence directory's `sources/`.

The comparison is of propagation of common binary pupil and OPD arrays into expected detector photon rates. Shared Gaussian277 input modes test the optical response to those modes. They do not validate external DM geometry, hardware command units, external calibration pipelines, external slope algorithms, detector noise, atmosphere, or the real Copper instrument.

## Analytical findings

Observed source supports the proposed exact-arithmetic equivalence. For FFT side L=512, pupil side p=240, offset o=136, and centered frequency coordinate q=k−255.5, the AOS/OOPAO prism is exp[−i(π/2)(|q₁|+|q₂|)]. SPECULA's folded surface before division is 511−|q₁|−|q₂|; negating that surface and using tilt_scale=4 produces the same mask up to constant phase. Constant phase has no effect on measured intensity.

AOS/OOPAO multiply the centrally padded field by exp[−iπ(1+1/L)(n₁+n₂)]. SPECULA uses an origin-padded field with exp[−iπ(n₁+n₂)/L], a shifted mask, and a final intensity roll by (136,136). The additional checkerboard in the first convention shifts the FFT by half its length; the native SPECULA mask shift cancels this. Translation and the final roll give the same image registration. Remaining carrier and constant phases disappear under squared modulus. This argument depends on the shared even dimensions and zero static shifts used here.

The modulation phasor sets are the same under the fixed permutation k_SPEC=(72−k_AOS) mod 96. Uniform weights make this permutation irrelevant in exact arithmetic. Different accumulation orders remain a Float32 numerical difference. No OPD sign reversal, fitted gain, registration shift, reflection, rotation, or interpolation is needed.

The matched SPECULA subclass overrides only cold geometry, focal-stop construction, and prism surface sign. Native propagation, field formation, loop, normalization, rebinning, CUDA capture, and lifecycle remain inherited. Its public native prescription cannot independently select the same combination: the 10% padding guard, circular stop, and opposite prism orientation differ. This is a valid explicitly labeled common-prescription implementation comparison; it is not a claim that unmodified native SPECULA implements this exact prescription.

For binary amplitude A, SPECULA's native rate normalization Φ × (D/p)² × ΣA equals the AOS physical-field energy Φ × (D/p)² × ΣA². With the square unit stop, transmission is one in exact arithmetic. OOPAO stores exposure-integrated photons and divides the binned result by the same exposure; the rate is equivalent. SPECULA's 512→64 area rebin is block integration after its native total normalization; AOS/OOPAO sum 8×8 blocks directly. Native SPECULA normalization can hide global propagation-energy roundoff, so agreement of final total flux alone is not independent evidence of FFT energy conservation.

Source evidence: AOS `src/wfs/pyramid/optics.jl:36,84,117,317`; `src/optics/focal_plane_modulation.jl:196,229`; OOPAO `Pyramid.py:326,457,635,804,1120`, `Source.py:270,294`, `Telescope.py:225`; SPECULA `modulated_pyramid.py:347,393,454,475,568,573,620`, `electric_field.py:194,247`, `toccd.py:31`, and `extrapolation_2d.py:280,419`.

## Review items

### MCR-001 — Disclose the modulation coordinate convention

Severity: P2, confidence: high. Classification: observed construction and derived physical interpretation. Disposition: primary adjudication pending; documentation clarification recommended, no production change requested.

All three engines use modulation coordinates spanning −π to π over 240 samples. With independently stated physical pitch D/240, consecutive samples are separated by 2π/239, so the phase gradient of nominal radius r corresponds to r×240/239 λ/D. At r=5, this is approximately 5.0209205 λ/D. The fixture's polynomial/Fourier OPD coordinates instead use physical pixel centers (2i−p−1)/p. These are different roles, and the common input remains valid.

Record nominal API radius 5 and the endpoint-inclusive modulation convention. Do not characterize the comparison as independent verification of exact physical modulation radius 5 λ/D. Changing one engine's radius to compensate would destroy the currently matched prescription. Required validation: preserve the declared coordinate formula in the protocol and confirm the realized modulation count/weights in metadata.

### MCR-002 — Freeze readout and error denominators before outputs

Severity: P1 for any claim of numerical qualification; confidence: high. Classification: incomplete protocol surface, not an observed adapter defect. Disposition: pending driver/evaluator review.

The proposal specifies rate, slope, and calibration thresholds but not yet all metric definitions. Before outputs, freeze the detector pupil slices, slope signs and order, valid-pixel mask, normalization, reference subtraction, and whether the reference is common or independently computed for each engine. A mask derived separately from each engine can remove discrepant pixels. Per-engine flat-reference subtraction is legitimate operationally but can hide static offsets; also compare raw flat readouts. A common readout applied to external intensity validates that composition, not the external native slope implementation.

Use full 64² rate arrays and a fixed denominator such as ‖I_AOS‖₂ for relative image error. For excited slopes, state the exact reference and baseline operation. Flat/piston signals need a predetermined absolute tolerance because relative error around zero is undefined. For interaction matrices use D[:,j]=(s(+h b_j)−s(−h b_j))/(2h), h=20 nm OPD, all 277 columns in declared order. Compare both Frobenius relative error and every column; publish column norms and the identity of the worst column. Near-zero reference columns require a fixed documented floor/absolute gate, not deletion after results. Perform host validation reductions in Float64 to avoid validator roundoff dominating Float32 gates.

Required validation: inspect the evaluator, ensure finite outputs, unchanged fixture hashes, no fitted transforms, and preserve failed metrics. All 277 columns passing verifies the finite-poke operator at ±20 nm; it does not prove an infinitesimal derivative or poke-size independence.

### MCR-003 — State the input timing boundary in native units

Severity: P2, confidence: high. Classification: observed. Disposition: documentation required; adapter design is acceptable under the narrower boundary.

SPECULA's upload converts meters into resident nanometers outside timing. OOPAO stages meters, then its timed public OPD setter allocates/copies and masks OPD; native phase formation occurs in the measured operation. These are authentic API overheads, but the starting representation is not literally the same meter-valued array in all engines. Define the boundary as native resident OPD representation to completed native detector-rate product, with phase formation, propagation, integration, and required lifecycle included. Exclude host transfers and external-unit conversion consistently as staging. Do not claim equal staging cost or complete simulator throughput.

SPECULA additionally computes focal PSFs/transmission inside its native propagation. AOS/OOPAO need not add artificial work solely to match that overhead; report it as a native implementation cost. Performance then describes delivered rate products through each implementation, not an isolated identical-count FFT kernel comparison.

Required validation: time source preparation/check_ready, trigger, post_trigger, detector binning, and completion. Synchronize before the start and after the complete output path on GPU; stream graph completion alone may precede postprocessing on another/default stream. Exercise different staged OPDs outside timed samples and confirm output changes so capture cannot silently replay stale inputs.

### MCR-004 — Proposed gates are acceptance choices, not error bounds

Severity: P2, confidence: high. Classification: review judgment. Disposition: suitable to preregister, conditional on MCR-002 definitions.

Proposed cross-engine gates: rate relative L2 ≤5e−5, flux relative error ≤1e−5, excited-slope relative L2 ≤2e−4, full finite-poke matrix relative Frobenius error ≤5e−4, maximum column relative error ≤2e−3. Proposed CPU/GPU gates: rate ≤2e−6 and matrix ≤2e−4. These are reasonable tight engineering equivalence targets for this experiment, but no supplied analysis establishes them as guaranteed Float32 bounds or science-accuracy tolerances.

Large phase arguments in cold phasor/prism generation amplify construction roundoff differently: AOS constructs many arguments in Float32, OOPAO constructs some in Float64 before casting, and SPECULA follows its native mixed construction. Modulation order, FFT implementations, units conversion, native normalization, and reduction order add differences. CPU/GPU gates can also be affected by independently constructed phase factors. A failed gate must remain a failure of the preregistered target while the cause is investigated. Do not relax a threshold merely because data nearly miss it.

Required validation: record actual array and FFT-result dtypes; finite and nonnegative rate checks; expected total flux from the shared sampled pupil; piston invariance; x/y tilt orientation through the frozen readout; all probe and all calibration gates. Optional diagnostics after failure can compare raw 512² intensities, phase factors, or Float64 construction; keep those distinct from the original pass/fail run.

### MCR-005 — Limit performance claims and preserve repetitions

Severity: P2, confidence: high. Classification: protocol requirement. Disposition: pending driver/result review.

The planned unpaced, one-completed-step-at-a-time measurement is appropriate for optical-node service time and sequential throughput. It cannot establish fixed-rate RTC latency, deadline reliability, queueing, or overload performance. Three fresh-process rounds improve independence but require per-round raw samples and dispersion, not only a pooled median. Use the same predefined warmup and samples for all engines, retain slow samples, and avoid high-tail claims without adequate sample counts. Rotate or interleave engine order across rounds to limit thermal/order confounding. Record initialization/capture separately from steady samples.

An eight-core affinity budget is not proof that each native engine uses eight workers. Report requested and realized FFT/jobs/backend threading; SPECULA NumPy FFT may remain serial. OOPAO's benchmark subclass forcing one FFT worker per job prevents nested oversubscription and is an honest policy override when disclosed. Record source/adapter/fixture hashes, dependency lock/version data, CPU/GPU model, affinity, FFT/BLAS/Julia thread settings, competing load, and GPU memory/batch/capture configuration. A GPU batch chosen from free memory can vary with concurrent processes; serialize benchmark GPU jobs and archive the realized value. Qualification and performance outputs must correspond to the same implementation and precision settings.

## Overall disposition

No confirmed optical mismatch was found in the explicit matched adapters by source and algebra review. The subclass construction is defensible and clearly separated from the unmodified native diagnostic. Numerical qualification remains unestablished until preregistered gates run. Final performance approval also requires reviewing the completed drivers/evaluator against MCR-002 through MCR-005; this document does not approve code that did not yet exist at inspection.

## Driver and evaluator follow-up

Reviewed subsequently available `common.py`, `protocol.toml`, `run_aos.jl`, `run_python.py`, `evaluate.py`, and `README.md`. This follow-up is source review only; no experiment was run by this reviewer. The primary reported the first OOPAO/AOS CPU preflight failures after the initial review. Those reported results do not alter the preregistered thresholds.

MCR-001: addressed in protocol and README by explicit endpoint convention and nominal radius.

MCR-002: substantially addressed. `common.py` fixes the 716-point circular support, four pupil slices, raw slope signs/order, and analytic incident divisor. `evaluate.py` uses Float64 readout and reductions, full images, raw slopes without flat subtraction, both relative and absolute checks, all 277 columns, and a fixed 10⁻¹² maximum-column-norm floor. The absolute slope metric is maximum absolute component, i.e. L∞; make this explicit in prose. A zero-norm reference calibration raises an error. All array hashes, finite rates, nonnegativity, and declared shapes are checked. Evaluator also distinguishes partial preflight from full calibration and marks timing eligibility false when qualification fails. Validation applies to a common readout on native images, as disclosed.

MCR-003: addressed. Python stages and synchronizes input before timing, then includes runner.step() and deviceSynchronize() inside timing. This waits for post_trigger rebinning as well as graph replay. AOS `step_graph!` calls `wait_graph_step!(step_graph_async!(graph))`; captured completion uses its blocking device-context synchronization. AOS timing therefore measures completed output without needing a redundant extra device wait. SPECULA's stable native field buffer is refreshed by check_ready→prepare_trigger before replay; graph writes/readout occur through inherited source paths. Its native blocking stream order and final device wait preserve completion.

MCR-004: retained. The primary reported CPU preflight maximum rate error 5.8467e−6 (passes), maximum slope component error 3.0505e−5 (fails), and maximum excited-slope relative error 3.0467e−4 (fails). These are primary-reported measurements, not independently recomputed by this reviewer. Passing the rate gate does not imply passing slope gates: the fixed four-quadrant linear map has operator norm 2/divisor on selected pixel L2 errors. The distinction is mathematically expected and is not a reason to weaken a gate.

A concrete source-supported hypothesis is cold centering-phase construction: AOS builds Float32 arguments up to approximately 3217 radians, while OOPAO computes its centering ramp in Float64 before casting the complex values. AOS's scale rounding and argument rounding can introduce phase differences of order 10⁻⁴ radians. Prism arguments reach approximately 803 radians as well. This is a discriminating investigation lead, not an established cause of the observed slope differences. Comparing realized phasor/mask arrays and varying one cold factor in a separately labeled diagnostic would distinguish it. Production arrays or accepted comparison outputs must not be fitted or changed to make these failures disappear.

MCR-005: partially addressed. Drivers use 32 declared warmups, 128 samples, explicit per-sample completed service times, varying fixture inputs, and hashes. AOS additionally executes nine allocation-probe preparation steps after warmups; disclose this diagnostic work if discussing an exactly identical warmup sequence. The recorded p95 is a descriptive sample quantile, not a tight tail estimate. Three fresh rounds and interleaved execution order remain orchestration requirements, not behavior enforced by these individual drivers. Final performance evidence remains unreviewed.

### MCR-006 — Confirm the regenerated AOS pupil matches the fixture

Severity: P2, confidence: high. Classification: observed driver gap. Disposition: reported to primary; fix pending at inspection.

`run_aos.jl` reads the shared amplitude array but builds its graph from telescope parameters without binding that array or asserting equality. The same constructors and parameters should regenerate the same binary pupil in this fixture, but the claimed shared input contract must be checked. Assert exact equality between the graph's prepared input amplitude and fixture amplitude during setup, or use a supported binding of the supplied array. Do not modify production APIs for this benchmark. Required validation: fail on a different fixture pupil, and record/assert equality in the ordinary preparation path.

### MCR-007 — Demonstrate live captured inputs and preserve derived provenance

Severity: P2, confidence: high. Classification: observed verification gaps. Disposition: reported to primary and evaluator author; fixes pending at inspection.

The drivers compare the last benchmark output to its preflight counterpart exactly. This detects inconsistent replay but, by itself, an always-stale capture can pass if preflight is equally stale. Add a positive control that independent excited fixture inputs produce distinct images. Full calibration already rejects an entirely zero reference response; retain that check. Finite and nonnegative rates alone do not prove refreshed inputs.

The evaluator records its own source digest but the readout and fixture routines are in `common.py`; include that helper's digest in derived provenance. Its output directory initially permits overwriting a previous summary and derived arrays. Reject nonempty output directories or use a versioned new directory so failed qualification evidence cannot disappear through a routine rerun.

### Verification of follow-up fixes

Reread the final evaluator guards and AOS setup assertions. MCR-006 is addressed by exact prepared-pupil equality. MCR-007 is addressed by per-run flat-versus-100 nm tilt image checks, a zero-calibration guard, common.py digest in derived provenance, and rejection of nonempty evaluation output. The evaluator now also requires recorded harness bytes to match current files, exact pinned source revisions, and clean recorded external status. Earlier records remain immutable and may be rejected after harness changes; rerun against finalized bytes instead of rewriting their recorded hashes.

No remaining driver/evaluator source blocker was identified in this review. That conclusion does not replace numerical gate results or approve a performance ranking whose required numerical gates fail.

## Bounded cold-factor intervention

The primary subsequently authorized a CPU-only diagnostic, limited to flat and 20 nm x tilt, to investigate the phase-construction hypothesis. Added `benchmarks/matched_copper/diagnostic.py`; it makes a separate AOS instance, exports both engines' realized centering phasors and prism masks, and changes only that diagnostic instance. No production implementation, accepted runner, frozen gate, or accepted result was changed. No performance measurement was made.

Command (successful run):

```text
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 taskset -c 14 /home/dgamroth/workspaces/codex/aos-matched-copper-evidence/venv/bin/python benchmarks/matched_copper/diagnostic.py --sources /home/dgamroth/workspaces/codex/aos-matched-copper-evidence/sources --fixture /home/dgamroth/workspaces/codex/aos-matched-copper-evidence/fixture --julia-project /home/dgamroth/workspaces/codex/aos-matched-copper-evidence/julia-env --output /home/dgamroth/workspaces/codex/aos-matched-copper-evidence/cold-factor-diagnostic
```

An initial attempt to select core 16 failed immediately because the process affinity allows only cores 0–15; it ran no optical computation. The successful diagnostic used one CPU core, one OOPAO job, one FFT worker, one Julia thread, and one BLAS thread. It executed two native OOPAO frames and eight AOS frames after native preparation. The primary was notified of its CPU use before execution and completion afterward.

Observed realized cold-factor differences (AOS versus OOPAO):

| Factor | Maximum complex difference | Relative L2 difference | Maximum phase difference |
|---|---:|---:|---:|
| Centering phasor | 2.24833e−4 | 7.33844e−5 | 2.24833e−4 rad |
| Prism mask | 7.54521e−5 | 1.79639e−5 | 7.54521e−5 rad |

Observed interventions, always compared with the unchanged native OOPAO two-probe output:

| Separate AOS instance | Flat rate relative L2 | Flat slope maximum error | Tilt20 rate relative L2 | Tilt20 slope maximum error |
|---|---:|---:|---:|---:|
| Original factors | 5.37406e−6 | 2.78538e−5 | 5.37475e−6 | 2.78937e−5 |
| Replace centering phasor only | 2.09586e−7 | 3.05867e−7 | 2.13018e−7 | 2.92568e−7 |
| Replace prism mask only | 5.37719e−6 | 2.78206e−5 | 5.37481e−6 | 2.77807e−5 |
| Replace both | 1.73233e−7 | 2.39374e−7 | 1.76915e−7 | 2.65971e−7 |

Derived conclusion: the centering-phasor construction difference is the dominant cause of the raw slope disagreement for these two CPU probes. Replacing it alone reduces the maximum slope error by approximately 91× and 95×; mask replacement alone does not materially improve it. This strengthens the earlier hypothesis to a confirmed local cause. It does not establish a universal cause for all seventeen probes, all 277 calibration columns, or CPU/GPU differences. Nor does it establish that copying an external array is a suitable production remedy.

The original baseline remains a failed qualification result under the frozen slope gate. The intervention is explicitly excluded from accepted qualification and timing records. Any future production remedy should compute the same intended analytical centering convention with a justified numerical construction, then rerun the unchanged validation gates on that new source revision.

Artifacts: `cold-factor-diagnostic/diagnostic.json` includes source/fixture/raw-array SHA256 values and every metric; the same directory holds realized complex factors and the two-probe images for every intervention. `cold-factor-diagnostic.log` records the successful command's native output. All changes by this reviewer remain limited to the diagnostic script and this review document, plus isolated raw diagnostic evidence.

## Campaign, collection, telemetry, and draft report review

Reviewed `campaign.py`, `collect.py`, `telemetry.py`, and the draft `benchmarks/results/platform/2026-10-22-matched-copper/README.md`. No campaign, calibration, or timing process was launched by this review. Read-only sysfs inspection confirmed this host's thread siblings are 0–1, 2–3, …, 14–15. The configured affinity 0,2,4,6,8,10,12,14 therefore selects eight distinct physical cores, as the draft states.

### MCR-008 — Collector call signature prevents completion

Severity: P1, confidence: high. Classification: confirmed source defect. Disposition: reported to primary; fix pending at inspection.

The new collector calls `evaluate.load_run` with three positional arguments in both validation and timing loops. The reviewed evaluator requires four: root, fixture hash, protocol hash, and AOS revision. Collection therefore raises TypeError before reducing performance. Parse the fixture once and pass its `aos_revision` in both calls. This does not require repeating optical work because neither raw evidence nor optical execution changes.

### MCR-009 — Campaign reruns can overwrite earlier logs

Severity: P2, confidence: high. Classification: confirmed source behavior. Disposition: reported to primary; fix pending at inspection.

The campaign opens each existing log with mode `w` before the driver rejects an existing nonempty result directory. Thus an accidental rerun can erase a useful prior log even though the driver protects raw arrays. Check for existing result/log paths before launch and create the log exclusively. This is particularly relevant to failed or interrupted campaigns, where logs explain why records remain incomplete.

### MCR-010 — Strengthen reduction and label checks

Severity: P2, confidence: high. Classification: bounded validation improvements. Disposition: reported to primary; pending adjudication.

The collector matches fixture, protocol, current driver bytes, replay output, native metadata, and batching across validation and timing. Also compare record `(engine, backend, cpu_workers)` with each CONFIGS tuple rather than trusting directory names. Recompute median, range, and p95 directly from retained samples; saved summaries can be checked against recomputed values. Derive numerical-qualification wording from collected checks rather than the currently hardcoded statement that raw-slope qualification failed. These changes improve reusable result reduction without changing the experiment.

Campaign selection currently uses `cpu_single_affinity` for CUDA even though the protocol has `gpu_affinity`. They are equal in this frozen campaign, so observed placement is correct. Selecting the explicit GPU field would preserve intent in later protocols.

The draft report accurately narrows the scope, discloses native SPECULA work and explicit matching changes, retains slope failures, and limits the causal diagnostic to two CPU probes. It should mention the nine extra AOS allocation-probe calls after the 32 declared warmups when describing warmup symmetry. The phrase `internal68²` needs a space. Describe the fixture's tilt probes as x/y-axis tilts; the seventeen-probe set does not include separate negative-amplitude tilt cases, although the central calibration does use both signs.

Telemetry is passive and records GPU state, process summaries, aggregate/per-core CPU counters, and load every approximately two seconds. It does not mark each timing sample's start/end, and `ps %cpu` is a process-lifetime statistic. It supports coarse shared-host context, not precise attribution of an individual slow sample to interference or throttling. The report's current descriptive service-time claims respect this limitation. CPU and GPU calibration summaries reported by the primary have not yet been independently recomputed by this reviewer; that remains the final numerical review step.

## Independent reduction of complete CPU1/CUDA qualification

All six full validation records are complete. Independently reduced their saved arrays using `independent-numerical-check.py`, which imports neither `common.py` nor `evaluate.py`. The independent implementation reconstructs the ordered 716-point circular support with explicit coordinate loops, evaluates both four-quadrant equations directly, and forms every 1432×277 central-difference matrix in Float64. It checks raw array SHA256, shapes, finite/nonnegative values, nonconstant flat/tilt outputs, exact source revisions, clean external status, fixture/protocol hashes, and all recorded current harness hashes.

Observed results independently reproduce the primary's conclusions:

| Comparison | Max rate relative L2 | Max raw slope component error | Matrix relative L2 | Max column relative L2 | Overall gates |
|---|---:|---:|---:|---:|---|
| CPU AOS→OOPAO | 5.84459e−6 | 3.03989e−5 | 8.26499e−6 | 1.79864e−4 | Fail: raw slopes |
| CPU AOS→SPECULA | 1.04501e−5 | 3.34243e−5 | 8.51980e−6 | 1.81105e−4 | Fail: raw slopes |
| CUDA AOS→OOPAO | 5.84780e−6 | 3.04387e−5 | 7.98520e−6 | 1.75685e−4 | Fail: raw slopes |
| CUDA AOS→SPECULA | 1.04540e−5 | 3.36221e−5 | 8.45865e−6 | 1.88038e−4 | Fail: raw slopes |
| AOS CPU→CUDA | 1.85010e−7 | 7.46382e−7 | 6.54151e−6 | 1.46994e−4 | Pass |
| OOPAO CPU→CUDA | 1.92063e−7 | 5.48566e−7 | 4.60401e−6 | 1.03885e−4 | Pass |
| SPECULA CPU→CUDA | 3.61154e−7 | 6.68253e−7 | 5.19600e−6 | 1.25823e−4 | Pass |

Each AOS/OOPAO cross-engine comparison fails eighteen raw-slope checks; each AOS/SPECULA comparison fails twenty-two. Failed check types are exclusively the frozen raw-slope component and excited-slope relative gates. No rate, flux, piston, conservation, or calibration check fails. Cross-engine maximum excited-slope relative errors are approximately 3.0461e−4 for CPU OOPAO and 4.9832e−4 for CPU SPECULA, with similar CUDA values. These failures remain explicit; successful calibration does not erase a static raw-slope discrepancy.

Every calibration column in all six runs is nonzero. Reference column norms range approximately 4.77e5–1.32e7 in the declared slope-per-metre units; no column activates the 10⁻¹² relative denominator floor. Thus the complete-column result does not depend on omitting or suppressing weak modes. Each engine's full calibration matrix norm is approximately 1.7747e8.

The shared analytic incident photon rate is 52,578,519.45919095 photons/s, and the fixed readout divisor is 73,433.68639551809. Maximum photon-conservation error across every probe and every positive/negative calibration frame is 3.48339e−7. Maximum piston image error is 1.50677e−7. Flat-to-100 nm tilt image changes have relative L2 approximately 0.04854 for every engine/backend, positively demonstrating live response rather than an always-stale capture.

Full independently reconstructed per-probe metrics, all per-column errors/norms, failed checks, and provenance are retained in `independent-numerical-check.json`; its script hash is embedded. This verifies saved software outputs and their arithmetic reduction. It is not hardware or instrument validation. CPU8 qualification and the three-round performance campaign were still running/pending and are outside this numerical review's completed scope.

Rereview confirms MCR-008 fixed by supplying the AOS revision to both collector load_run calls. MCR-009 is fixed by refusing existing output/log/resource paths and opening logs exclusively. MCR-010's tuple checks, independent timing-statistic recalculation, dynamic failure wording, and explicit gpu_affinity selection are implemented. Campaign whole-process resource counters now use before/after resource.getrusage deltas; their scope correctly includes startup and validation and must not be presented as steady-state-only resource costs.

## Final independent performance and CPU8 audit

The complete campaign has finished. `independent-performance-check.py` performs a separate reduction without importing benchmark helpers. It verifies all 24 timing records against current runner bytes, exact source revisions, fixture and protocol hashes, expected engine/backend/worker tuples, and the corresponding completed validation records. It independently rebuilds the image/readout parity checks and CPU8 calibration matrices. Timing statistics use Python's sorted-sample median and an explicitly interpolated 95th percentile, independently of the collector's NumPy reduction.

Observed: eight configurations each have three fresh processes, 32 declared warmups and 128 positive finite service-time samples per process. The journal has exactly 24 unique successful configuration/round pairs, the prescribed seeded and rotated order, nonoverlapping process intervals, correct CPU affinity and thread-limit environment, and nonnegative whole-process resource counters. All raw samples remain present. The collector's medians, minima, maxima, p95 values, median-of-round-medians, and round-median ranges agree with the independent reduction.

All 24 performance probe sets pass the existing 2e−6 rate budget and all raw-slope, flux, conservation, piston, and live-input gates against their own validation configuration. Maximum rate difference across fresh processes is 7.98861e−8; maximum raw-slope component difference is 3.22490e−7. Matching metadata and runtime batching were verified. The in-process exact repeated-output guard remains in the runner and every run completed successfully.

The collector no longer requires bitwise equality between independently planned processes. Source inspection supports that change: AOS CPU uses FFTW.MEASURE, so different legal plans can change Float32 reduction order. The primary reported making this change before an observed performance failure and without changing thresholds. The final independent audit establishes that the saved fresh-process outputs satisfy the original numerical budgets; it does not rely on hypothetical bitwise reproducibility.

CPU 1/8-thread validation independently passes:

| Engine | Max rate relative L2 | Calibration relative L2 | Max column relative L2 |
|---|---:|---:|---:|
| AOS | 7.17570e−8 | 3.98798e−6 | 9.07405e−5 |
| OOPAO | 0 | 0 | 0 |

All 277 columns are nonzero and above the denominator floor. The earlier seven independently checked numerical comparison outcomes, all scalar metrics, and every failed-check name/count agree with the final collector summary. Cross-engine raw-slope failures remain present, and the performance summary explicitly retains failed qualification status.

Independently confirmed service times, milliseconds:

| Configuration | Median of round medians | Round-median range | Per-round p95 |
|---|---:|---:|---|
| AOS CPU1 | 363.795442 | 351.727208–379.224124 | 395.471358, 408.554200, 468.540230 |
| AOS CPU8 | 257.094599 | 243.626895–344.822003 | 291.004099, 388.458468, 336.795297 |
| OOPAO CPU1 | 499.322500 | 487.021151–500.563020 | 565.754251, 525.088567, 580.299734 |
| OOPAO CPU8 | 202.701325 | 173.491643–207.147070 | 200.803650, 223.973320, 231.604662 |
| SPECULA CPU1 | 2186.163476 | 2097.076589–2236.071709 | 2256.200799, 2298.267564, 2362.932546 |
| AOS CUDA | 2.664345 | 2.663456–2.675041 | 2.686481, 2.862082, 2.692617 |
| OOPAO CUDA | 6.388397 | 6.331867–6.406886 | 6.399266, 6.595818, 6.608913 |
| SPECULA CUDA | 4.885390 | 4.685748–4.909473 | 5.376829, 5.356018, 5.241457 |

These measurements support descriptive comparisons of the configured implementations on this shared host. The large CPU8 round dispersion, especially AOS's 243.6–344.8 ms range, prevents interpreting these numbers as isolated hardware limits or precise scaling estimates. The 128-sample p95s are empirical descriptive quantiles, not high-confidence tail bounds. No RTC, fixed-arrival, WCET, overload, or complete-simulator performance claim is established.

All nine AOS performance records contain eight zero-byte host-allocation observations. All three AOS CUDA runs completed the runner's unchanged-pool assertion. SPECULA's observed default CuPy pool remains unchanged across all three runs. OOPAO's default-pool live-byte count decreases by exactly 16,384 bytes in each of its three CUDA rounds. This is not evidence of a leak, but it must not be described as exact pool equality. These library-specific observations do not measure comparable total or peak VRAM; they omit different native/FFT/private-graph storage. The finalized report's resource caveat correctly limits interpretation.

Artifacts `independent-performance-check.py` and `independent-performance-check.json` preserve the independent procedure, all reduced results, per-round parity metrics, timing-record hashes, journal hash, source/fixture hashes, and CPU8 column errors. No optical propagation, GPU job, or timing job was rerun by this audit. The results README's numerical and performance tables agree with the independent audit; its remaining resource wording was reported to the primary for the minor correction from first-round to every-round OOPAO pool reduction.

## Final disposition of review items

- MCR-001: addressed; nominal modulation and endpoint convention disclosed.
- MCR-002: addressed; frozen complete readout/metric rules implemented and independently reduced, with failures preserved.
- MCR-003: addressed; completed-output timing boundaries and native staging/work differences disclosed and checked.
- MCR-004: accepted engineering thresholds retained; observed cross-engine raw-slope failures remain unresolved qualification failures. Two-probe CPU centering-phasor causality is independently demonstrated, without generalizing to all modes or engines.
- MCR-005: addressed for the stated descriptive service-time scope; all 24 runs independently checked, with shared-host and tail limitations retained.
- MCR-006: addressed by exact prepared-pupil equality assertion.
- MCR-007: addressed by live-input checks, helper provenance, immutable evaluation destinations, and independent saved-array verification.
- MCR-008: fixed; completed collector and independent reduction confirm correct call integration.
- MCR-009: fixed; campaign refuses existing output/log paths and creates logs exclusively.
- MCR-010: addressed; configuration checks, independent statistics, dynamic qualification wording, and explicit GPU affinity are implemented and verified.

Final review conclusion: no remaining harness, reduction, or reporting blocker was identified for publication of this completed benchmark as a reproducible comparison with explicitly failed cross-engine raw-slope qualification. The numerical failures are not repaired by this task and must not be represented as successful full cross-engine equivalence. No production or instrument-model change was made by this review. Archival packaging and the final commit remain the primary agent's responsibility.

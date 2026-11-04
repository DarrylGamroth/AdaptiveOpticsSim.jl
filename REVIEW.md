# AOS 0.3.20 accumulated optics review

Review candidate: `codex/prepared-proper-paths`, implementation `e931beb6`.
Release base: `v0.3.19` (`cbfc4956`). Cleanup `e0c5ff9d` was fast-forwarded
onto the performance branch; merge `5a9c16b2` retains the already released
optional model-time uncertainty fix. No history was rewritten.

## Scope and preserved contract

- Shared compiled optical paths: batched FFT/GEMM execution, fused mappings,
  shared-prefix/checkpoint scheduling and portable fallbacks.
- Prepared scalar/runtime-OPD trains and calibrated spectral synthesis.
- Captured exposure helpers retain Float64 exposure accumulation, live input
  bindings, stream ordering and callback partial-commit semantics.
- Cleanup separates representation, field primitives and preparation from the
  compiled executor; it changes no parsed optical expressions or operand order.

Preserve fitted prescriptions, spectra, detector geometry/response, coordinates,
normalization, calibration and scientific asset bytes. The ordinary numerical
executors remain independent oracles. FP32 field precision is retained; bounded
prepared phase arithmetic keeps its original Float64 fallback. Ordinary and
optimized executors are compared with the existing tolerances, not a claim of
universal bitwise equality. YAML still supplies the optical prescription and
calibration; hardware policies select implementations during preparation.

CPU/CUDA/AMDGPU paths remain available. Device, shape, precision and math-policy
qualification are separate; RTX 3080 measurements do not qualify other GPUs.
Private PROPER train/batch/cancellation lowerings require coordinated ownership;
new prescriptions should start with its public prepared field API.

## Independent review and dispositions

Separate read-only Sol reviews assessed the full accumulated AOS and PROPER
release diffs, including numerical equations, ownership, exposure concurrency,
backend fallback and package compatibility.

| ID | Severity | Evidence and disposition |
| --- | --- | --- |
| AOS-R01 | Low | Extension guide confused public separate-scratch and private in-place contracts. Corrected and independently checked against train/batch/cancellation call sites. |
| AOS-R02 | Test gate | Bare CPU suite exposed seven unclassified accelerator kernels. Added explicit deferred inventory entries, with real optional selectors; production dispatch and retained parity assertions independently checked. No inventory assertion was removed. |
| PRP-OUT-01 | Medium, companion | PROPER admitted structurally mutable and repeated-index outputs to parallel writes. Companion fix rejects these before enqueueing; poisoned-scratch, positive/negative strided-view and capture checks cover the boundary. |

No remaining source-review blocker. Execution gates below determine validation;
source review alone does not establish numerical or timing qualification.

## Verification

Fresh local gates on implementation `e931beb6` with PROPER `1d4d16e`:

| Gate | Result |
| --- | --- |
| Bare AOS `Pkg.test`, Julia 1.13.3, one thread | Passed, including Aqua and complete kernel inventory. |
| Julia 1.12.7 quality, algorithm-graphs, prepared-execution-interfaces, execution-strategies, spectral-synthesis | Passed; 1,799 testset assertions plus standalone checks. |
| Julia 1.12.7 compiled paths / prepared trains | Passed; 4,431 / 92 testset assertions. |
| Captured CUDA compiled paths | Passed; 6,659 testset assertions. |
| Captured CUDA exposure helpers / external-pupil integration / spectral synthesis | Passed; 76 / 63 / 1,222 testset assertions. |
| Frozen fitted SCC / LLOWFS captured replay | Passed; 54 / 27 assertions; three changing states and reset/replay. |
| Companion PROPER full CPU, Julia 1.10.12 and 1.13.3, two threads | Passed, including Aqua. |
| Companion captured CUDA propagation / phase | Passed; 1,582 / 705 testset assertions. |
| Companion final AMD propagation / phase | Passed; 1,587 / 705 testset assertions; CUDA remained unloaded. |

Frozen SCC maximum relative L2 / normalized peak errors were
**1.0957712 × 10⁻⁶ / 3.822023 × 10⁻⁶**. LLOWFS intensity errors were
**1.6359095 × 10⁻⁷ / 4.4253753 × 10⁻⁷**; field relative L2 was
**5.1883023 × 10⁻⁷**. The acceptance limits remain **10⁻⁵**. Reset replay
reproduced the first changed-input products exactly. These are comparisons to
the retained original executors, not external physical bench validation.

CPU environments resolved KA 0.10.4 / AK 0.5.11 / AOC 0.17.3. CUDA used RTX
3080, CUDA 6.4.2 / KA 0.9.45 / AK 0.5.11; GPU jobs ran serially with one
Julia/BLAS thread. Validation used temporary source environments, leaving frozen
application manifests untouched. One initial bare CPU run failed the inventory;
the classified-kernel fix passed the same complete suite. Julia 1.12 setup emitted
a `GPUCompiler.compile_hook` deprecation warning; no such use exists in AOS or
PROPER source, and the runtime checks passed.

Reproduce with a temporary `ENV` that develops both implementation commits;
add Test/Aqua/YAML and the selector's direct dependencies there. For fitted
CUDA gates, copy the frozen SPIDERS environment, change only the temporary AOS /
PROPER version constraints, and develop AOS/PROPER/RTC task sources. Retain its
other application pins. With `APP` pointing to that frozen application:

```sh
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
julia --startup-file=no --project=ENV -e 'using Pkg; Pkg.test("AdaptiveOpticsSim")'
julia +1.12 --startup-file=no --project=ENV -e 'using Pkg; Pkg.test("AdaptiveOpticsSim"; test_args=["quality", "algorithm-graphs", "prepared-execution-interfaces", "execution-strategies", "spectral-synthesis"])'
julia +1.12 --startup-file=no --project=ENV -e 'include("test/proper/compiled_paths.jl"); include("test/proper/prepared_train.jl"); check_prepared_train(Array,AdaptiveOpticsSim.Backends.HostComputeDevice())'
julia --startup-file=no -O2 --project=ENV test/cuda_compiled_paths.jl
julia --startup-file=no -O2 --project=ENV test/cuda_external_exposure.jl
julia --startup-file=no -O2 --project=ENV test/cuda_external_pupil_execution.jl
julia --startup-file=no -O2 --project=ENV test/cuda_spectral_synthesis.jl
julia --startup-file=no -O2 --project=ENV test/cuda_compiled_paths_fitted.jl APP
julia --startup-file=no -O2 --project=ENV test/cuda_proper_train_fitted.jl APP
```

Set Julia/BLAS threads to one for all focused commands. Fitted checks ran
without benchmark/profile flags; compilation and preparation timing was not
used as steady-state performance evidence. Temporary environments and logs are
removed after recording the results. Source worktrees remain for review.

Previous cleanup verification on the exact merged optical expressions includes
CPU, CUDA and AMD prepared propagation/trains/compiled paths, changing inputs,
reset/replay and frozen SCC/LLOWFS checks. The frozen model SHA-256 remains
`18a3e6c2c658bb23e5f83434c816c8f06784e7c1c83fa318998559760f88b969`;
calibration identity remains
`c2c627318e70ac6ffd06d8579554ae148f5a388f17c49dd7f53dc76e269a42eb`.
Full evidence and application provenance remain in the instrument branch's
`spiders/ANALYSIS_PLAN.md` and `spiders/performance/terminal-transfer.json`.

## Performance evidence and release order

The last matched bounded source pair, before expression-preserving cleanup,
measured **42.806332 → 43.532471 completed updates/s per camera**, with median
**22.642687 → 22.298623 ms**. Zero camera losses/rejections; one latest-command
gap in each run. All 32 measured updates exceeded the 20 ms scheduled budget.
This is neither 50 Hz qualification nor causal camera-to-DM/tail-latency evidence.
No new timing claim is made for release preparation. Remaining measured SCC cost
includes 774 512² FFT leaves (~4.690 ms summed GPU time) and mapping/memory passes.

Release PROPER **0.4.2 first**, then AOS **0.3.20** with its existing
`Proper = "0.4.2"` compatibility bound. Review the accumulated diff with
`git diff v0.3.19..codex/prepared-proper-paths` and the companion PROPER diff
from `v0.4.1`. Optional extension tests must accompany future changes to private
PROPER lowerings. No release tag, publication, registry update, application pin
change, deployed simulator or GUI action is part of this preparation. Human
release approval remains pending; other device cohorts and exhaustive external
reference matrices are not newly qualified by these checks.

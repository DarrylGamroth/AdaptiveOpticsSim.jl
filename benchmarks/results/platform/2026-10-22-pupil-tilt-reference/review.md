# Independent review: withdraw shifted-mask propagation

Date: 2026-10-22. AOS worktree: `/home/dgamroth/workspaces/codex/AdaptiveOpticsSim-batches`, removal diff from `4d00d534`. Instrument worktrees: `/home/dgamroth/workspaces/codex/GPI2Sim-pupil-tilt` and `/home/dgamroth/workspaces/codex/REVOLTCopperSim-pupil-tilt`.

## Final disposition

Accept the removal and migration to physical pupil tilt. No blocking production, API, ownership, capture, or sampling defect was identified. REM-001 is resolved. The completed physical-model numerical, capture, external-reference, instrument, and resource evidence supports this narrowly stated acceptance. Maintained documentation and any additional package-wide checks are finalization work owned by the primary agent; this review does not represent unfinished checks as passed.

The decision to remove the approximation is supported separately by the archived discrete-operator assessment and external comparisons. Historical shifted-mask timing and calibration evidence does not qualify the retained physical pupil-tilt model.

## REM-001 — P2, confirmed: GPI convergence harness labels physical results as shifted-mask results

- Confidence: high.
- Affected file: `/home/dgamroth/workspaces/codex/GPI2Sim-pupil-tilt/benchmarks/gpi_sensor_convergence_cuda.jl`, originally lines 97 and 171.
- Evidence: the harness records `"strategy" => "PyramidShiftedMaskStrategy"` for a graph now explicitly prepared with `PyramidPupilTiltStrategy`. Its additional pupil-tilt reference uses the same formulation, but the comparison remains named `shifted_vs_pupil_tilt`.
- Impact: rerunning the maintained harness would misidentify the actual scientific model and present a same-model comparison as an approximation comparison.
- Remediation: derive the model label from the actual prepared owner and rename or remove the obsolete comparison. Archived historical results should retain their original labels.
- Validation: inspect the resulting metadata construction and reference labeling; no new GPU measurement is necessary solely to validate string/metadata correction.
- Disposition: resolved. The harness now identifies `PyramidPupilTiltStrategy`, asserts the fully qualified type on the prepared owner, removes the obsolete same-model reference, records photon-rate/check fields, and explicitly labels comparisons as physical-model sampling refinement. Active benchmark source searches found no renderer references to the removed metadata keys. Historical archived results retain their original interpretation.

## Source/API assessment

- `PyramidShiftedMaskStrategy`, its compact modulation owner, full/separable shifted-mask workspaces, preparation paths, kernels, propagation methods, exact-target specializations, workspace bindings, and capture-capability specializations are removed. Active AOS source, extensions, tests, and scripts contain no unresolved references to the removed symbols.
- Public physical strategy and its dispatch family remain. The graph-file parser accepts `pupil_tilt` and explicitly rejects the removed `shifted_mask` option with `AlgorithmGraphError`; the regression test checks rejection. No compatibility forwarding or silent reinterpretation of a removed strategy is introduced.
- The existing physical propagation implementation is retained: the modulation phase multiplies the pupil before each forward FFT; the fixed physical pyramid mask is applied before inverse propagation. Circular/sample-path definitions, weights, normalization, pupil coordinates, centering, detector sampling, and output shapes are unchanged by the removal.
- Physical face-registration parameters remain part of constructors and graph configurations. Fixed-mask construction is unchanged. The rewritten registration tests independently assemble the four-face mask and physical FFT sequence for paired and independent shifts with nonzero modulation, then check intensity, detector sampling, and zero warmed allocation.
- `_enqueue_pyramid_rate_batch!` for `PyramidModulationBatchWorkspace` remains intact, including output reset, complete phase-weighted field initialization, forward FFT, physical mask, inverse FFT, and accumulation. Its source/no-LGS accelerator capture eligibility remains intact. Unsupported source/LGS/backend combinations still select the existing unsupported capability path.
- Exact-device validation still checks physical phases, field stack, operating weights, and other front-end arrays. Workspace bindings still include both FFT plans and the concrete bound arrays; removed methods concern only withdrawn owners.
- `update_cycle_averaged_circular_modulation!(wfs, radius)` now directly delegates to the existing physical phase update. Validation of radius, precision, square phase shape, and fixed point count remains in that owner. Updates write into the existing phase array and preserve weights and buffer identities. The rewritten CUDA tests cover disabled/changed radii, fresh-reference agreement, pointer stability, and replay of a graph captured before phase updates.
- AOS active benchmark selectors now retain physical pupil tilt. Historical archived comparisons remain readable and explicitly identifiable. No active source call to the withdrawn Julia strategy was found in the inspected instrument source/test trees.

## Instrument migration

- GPI `src/graphs.jl` changes the selected strategy to physical pupil tilt. Resolution, modulation count, pupil samples, separation, edges, face registration, guide-star inputs, and graph boundary remain unchanged.
- COPPER imports/selects physical pupil tilt and updates its strategy assertion. Resolution, modulation count, modulation radius, detector geometry, diffraction padding, centering, source, and actuator model remain unchanged.
- COPPER numerical-mode benchmark storage introspection now checks the actual dense phase cube instead of removed offset/coordinate arrays. Its optional mode filter affects which explicit benchmark cases run; it does not alter preset sampling.
- Dense physical phase storage increases memory and per-point forward FFT work relative to the removed approximation. Default and fine-grid physical-model runs are required to establish actual runtime/resource behavior; prior shifted-mask performance figures must not be reused.

## Verification observed at initial review

- Inspected the complete production removal diff and selected test/benchmark/instrument diffs. `git diff --check` passed for AOS.
- Inspected `removal-cpu.log` and `removal-api.log`: focused CPU tests and the reviewed API/strategy checks passed, including fractional registration and uncentered modulation-camera paths.
- The serialized CUDA campaign had completed its physical numerical/capture batch suite and was running the phase-storage/update suite at this snapshot. Final COPPER and GPI validation results remain pending.
- This reviewer performed read-only source/log inspection and wrote this review artifact. No production files were edited and no CUDA work was launched.


## Completed physical-model verification

- Inspected all three completed generic CUDA logs: numerical/capture invariance 457/457; phase storage plus captured in-place updates 487/487; preparation policy 378/378. These are 1,322 passing checks. Earlier removal-focused CPU/API logs total 2,593 and 85 passing checks respectively.
- The post-removal Float32 external-reference run reproduces the previously qualified physical outputs: relative L2 errors versus OOPAO are 8.046 × 10⁻⁶ for the flat image, 7.932 × 10⁻⁶ for the ramp image, and 7.562 × 10⁻⁵ for the tip derivative. This is the side-256/count-32 tuned case; it is not a universal instrument-fidelity claim.
- Independently reran the binary integration comparison against `parity-inputs` and produced `integration-comparison-review.json`. That explicit view contains all COPPER records and the comparable GPI default records; it deliberately excludes new GPI fine/transfer records from the byte-equality comparison. The comparison passed.
- COPPER physical-model baseline versus selected batching at N=240/M=96: rate relative L2 8.365 × 10⁻⁸; complete 277-column interaction relative L2 5.030 × 10⁻⁶; worst column relative error 1.147 × 10⁻⁴. The frozen-calibration 320-frame loop has maximum command difference 9.237 × 10⁻¹⁴ m and maximum residual-wavefront difference 3.977 × 10⁻¹⁵ m. Noisy atmosphere frames are exact. Completed atmosphere graph service median changes from 3.867 to 2.883 ms with the same physical formulation.
- GPI default physical baseline versus post-removal output: all 17 comparable saved arrays are byte-identical, including calibration, coefficient histories, and final frames. The additional fine and transfer runs are evaluated by their own physical-model scientific replay, not claimed byte-identical to a different sampling level.
- Parsed completed `scientific-after.log`: source hashes match, CPU scientific replay is complete, and all reported runs have zero saturation samples. The independent reconstructed basis differs by at most 7.57 × 10⁻¹⁷; the Gram check error is 5.19 × 10⁻¹³. Recorded 320-frame tail residuals and repeated complete-cycle service measurements are:

| Physical GPI case | Tail residual WFE, nm | Complete-cycle median, ms |
|---|---:|---:|
| Default, noiseless | 139.8366 | 3.478 |
| Default, noisy | 139.9203 | 3.516 |
| Fine, noiseless | 139.6296 | 320.79 |
| Fine, noisy | 139.7099 | 322.12 |
| Transfer, noiseless | 139.7409 | 320.97 |

These finite, 35-mode illustrative controller results establish repeatable software behavior for the recorded inputs. They do not establish deployed-controller stability, a continuum convergence bound, hardware fidelity, or fixed-rate deadline compliance. In particular, the expensive fine run is an offline numerical reference on this GPU.

- Inspected complete COPPER public numerical-mode artifacts. Both calibration and production builders complete for default, fast, and accurate sampling; recorded shape/replay/input-effect checks pass and warmed host allocations are zero. Production completed-HIL service medians are 5.225 ms (480/32), 2.923 ms (240/96), and 457.202 ms (1440/256). The retained direct physical phase cube occupies 58,982,400; 44,236,800; and 4,246,732,800 bytes respectively. The latter is about 3.955 GiB before other workspaces and FFT plans; the successful run establishes fit on the measured host, not a general GPU memory guarantee.

## Claim boundaries

The original shifted-mask comparisons were biased by their different discrete optical formulation. Their fast timings and calibration/loop results remain historical evidence only and do not qualify current physical-source accuracy, controller behavior, or runtime. Current performance statements must use physical-versus-physical measurements above and identify the measured hardware/runtime and completion boundary.

The removal intentionally eliminates a public approximation during the authorized breaking refactor. Callers requesting the retired graph-file spelling are rejected rather than silently remapped. Physical face registration and modulation sampling remain supported. Further optimization must preserve that physical discrete model and revalidate photons, response matrices, capture, closed loops, and resource use.

No unresolved review finding remains. The reviewer launched no CUDA jobs and changed no production files; independent follow-up consisted of source/API inspection, artifact parsing, and rerunning CPU comparisons against hashed archived arrays.

# Final independent review — centering precision

Reviewed saved evidence only after the user stopped remaining work. No new benchmark or physical measurement was run.

## Disposition

No outstanding production-code correctness finding. Production commit `69f91214760ac91d9fc8b5a358ad6c1edd4d0fb6` preserves phase conventions and runtime ownership while improving cold centering-factor precision. The obsolete KA kernel inventory entry was removed. The explicitly precision-corrected SPECULA common-prescription adapter is separately reviewed and must be described as an adapter change, not as unmodified SPECULA qualification.

## Numerical qualification

Independent Float64 reductions from hashed saved Float32 arrays confirm all nine comparisons across eight configurations: four AOS/external comparisons, three CPU/CUDA comparisons, and two CPU1/CPU8 comparisons. Each full comparison covers 17 probes and every one of 277 central-difference calibration columns; no flat subtraction, gain fitting, registration, or relaxed tolerance is used. The original protocol hash is unchanged.

| Comparison | Maximum raw slope absolute error | Calibration relative L2 | Maximum column relative L2 |
|---|---:|---:|---:|
| cpu-oopao | 1.3531288e-06 | 8.2543445e-06 | 0.00018408439 |
| cpu-specula | 1.2899606e-06 | 8.239009e-06 | 0.00018947748 |
| cuda-oopao | 1.2799867e-06 | 7.9754648e-06 | 0.00017490908 |
| cuda-specula | 1.3265317e-06 | 8.169221e-06 | 0.00018600653 |
| aos-cpu-cuda | 7.5801809e-07 | 6.5709351e-06 | 0.00014585327 |
| oopao-cpu-cuda | 5.4856572e-07 | 4.6040069e-06 | 0.0001038849 |
| specula-cpu-cuda | 4.3885258e-07 | 5.2100868e-06 | 0.00013283962 |
| aos-cpu-threading | 2.3771181e-07 | 4.1564779e-06 | 9.5376642e-05 |
| oopao-cpu-threading | 0 | 0 | 0 |

Source, fixture and array hashes were checked. External revisions match sources.json and their recorded before/after source status is clean. Corrected SPECULA records explicitly identify `prism_exponential=Float64_then_runtime_storage`. AOS records identify the production fix revision. Focused logs show fail-before 32 failures, pass-after 160/160, CUDA 66/66, and KA CPU Matrix 1038/1038. The broader CPU log now completes successfully; the earlier Pkg environment failure remains preserved separately.

## Original SPECULA arithmetic remains a failed comparison

Independent reduction of the retained original-adapter CPU/CUDA records verifies the original adapter snapshot hash and all17 absolute raw-slope failures on each backend. Against current AOS the maxima are 1.21016923 × 10⁻⁵ (CPU) and 1.22213793 × 10⁻⁵ (CUDA), above 2 × 10⁻⁶. These are original prism arithmetic within the pre-existing common-prescription adapter, not an untouched native optical prescription. No full-calibration claim is made for these negative preflights.

OOPAO versus original-arithmetic SPECULA differs by 1.20617966 × 10⁻⁵ (CPU) and 1.21415880 × 10⁻⁵ (CUDA). Since either exceeds twice the 2 × 10⁻⁶ gate, triangle inequality rules out any AOS output passing that absolute gate against both original external outputs simultaneously. The prism-only precision diagnostic supports the narrow adapter correction; this is not justification to silently relabel original arithmetic as passing.

## Partial performance evidence

Six completed round-one records were independently checked for source/fixture provenance, 128 finite positive timing samples after 32 warmups, median/p95/min/max arithmetic, and numerical replay parity with validation records. Completed AOS CPU allocation probes are zero. Eighteen planned records are incomplete or absent. SPECULA CPU round one was interrupted and AOS CUDA has no completed timing record. No three-round campaign conclusion, complete CPU/CUDA performance comparison, or before/after regression conclusion is supported.

| Completed record | Median ms | p95 ms |
|---|---:|---:|
| benchmark-aos-cpu1-r1 | 362.14622 | 436.17264 |
| benchmark-aos-cpu8-r1 | 257.00296 | 315.49292 |
| benchmark-oopao-cpu1-r1 | 499.03191 | 567.78585 |
| benchmark-oopao-cpu8-r1 | 207.76829 | 226.2462 |
| benchmark-oopao-cuda-r1 | 6.38591 | 6.582831 |
| benchmark-specula-cuda-r1 | 4.8397835 | 5.2841894 |

All reported timings are descriptive shared-host optical service times. The user explicitly stopped remaining measurements; user-stop.json records the interruption. Timing incompleteness is intentional and must remain explicit.

## Reproducible independent artifacts

- independent-numerical-check.py / .json: seven cross-engine/device full reductions.
- independent-partial-review.py / .json: two threading full reductions and six completed timing records.
- independent-native-negative-review.py / .json: original arithmetic failures and external disagreement.
- centering-code-review.md: mathematical analysis, resolved finding CP-01, and adapter review.

The reducers import neither common.py nor evaluate.py. Their JSON results include script/provenance hashes. The original 24-record performance checker was not executed after the stop request.

Remaining limitation: this validates the stated matched numerical prescription and software evidence, not physical instrument fidelity or hardware timing guarantees.

# Python OOPAO vs AdaptiveOptics.jl

This document records the concrete differences between the Python OOPAO stack
(`OOPAO`, plus `pyTomoAO` where tomography is involved) and the Julia port
`AdaptiveOptics.jl`.

The goal is not line-by-line implementation parity. The goal is:

- feature parity for the AO workflows that matter,
- numerical fidelity against deterministic reference cases,
- and an idiomatic Julia implementation internally.

This review is based on direct inspection of the current code in:

- `OOPAO/OOPAO/*.py`
- `tutorials/*.py`
- `../pyTomoAO/pyTomoAO/*.py`
- `../AdaptiveOptics.jl/src/**/*.jl`
- `../AdaptiveOptics.jl/docs/*.md`

## How to read this

- `Intentional` means the Julia behavior is different by design.
- `Gap` means the Julia port is not yet equivalent at the feature or workflow
  level.
- `Extension` means Julia supports something beyond the main OOPAO surface.

## Executive Summary

The Julia port is already close to OOPAO in core optics, atmosphere, SH/Pyramid/
BioEdge sensing, GSC optical gains, LiFT interaction matrices, compact
closed-loop traces, transfer-function workflows, and pyTomoAO-backed
tomography.

The biggest differences today are:

1. The Julia package is architecturally different on purpose: multiple
   dispatch, param/state separation, trait-driven backends, column-major
   internals, structured errors, and cached workspaces.
2. Some OOPAO surface features are still missing from Julia, especially
   `tools/*`, GUI/display helpers, and some detector/telescope convenience
   features.
3. A few workflow-level fidelity gaps remain, most notably the full
   atmosphere-driven Pyramid+GSC tutorial trace under huge aberrations.
4. In a few places Julia intentionally keeps the core narrower and cleaner than
   OOPAO, with optional or future features deferred instead of baked in.

## Intentional Architectural Differences

| Area | Python OOPAO / pyTomoAO | Julia AdaptiveOptics.jl | Why Julia differs | Status |
|---|---|---|---|---|
| Package structure | Many classes across `OOPAO/*.py` and external `pyTomoAO` | One top-level `AdaptiveOptics` module with subareas under `src/` | Idiomatic Julia favors one package surface with multiple dispatch instead of class trees and external manager objects | Intentional |
| Programming model | OO objects with mutable state and operator overloading like `src*tel*wfs` | Explicit functions like `compute_psf!`, `measure!`, `capture!`, `advance!`, `reconstruct!` | Clearer dispatch, easier testing, easier zero-allocation hot paths | Intentional |
| Type system | Runtime duck typing, class attributes, NumPy/CuPy conventions | Parametric structs, traits, backend dispatch, structured interfaces in `src/Core/types.jl` | Needed for performance, specialization, GPU portability, and SciML-style APIs | Intentional |
| Array layout | NumPy row-major assumptions are common, especially in tomography exports | Julia internals are column-major; Python ordering is adapted only at reference boundaries | Native Julia memory order is the right internal choice for performance and readability | Intentional |
| Errors | `OopaoError`, warnings, and many `print(...)` diagnostics | `AdaptiveOpticsError`, `InvalidConfiguration`, `DimensionMismatchError`, `UnsupportedAlgorithm`, plus `Logging` | More testable and composable; avoids side effects in library code | Intentional |
| Randomness | Several Python components seed from wall-clock time or local random states | Julia centralizes RNG and supports deterministic fixed-seed workflows | Deterministic regression and RTC-style repeatability were explicit design goals | Intentional |
| FFT handling | OOPAO relies on NumPy FFTs and some ad hoc GPU/CuPy branching | Julia caches FFT plans in long-lived objects and uses backend-aware FFT dispatch | Reduces allocations and planning overhead; enables a cleaner CPU/GPU path split | Intentional |
| Parallelism | OOPAO uses `joblib`, helper tools like `set_paralleling_setup.py`, and module-specific choices | Julia uses threads and backend traits; `KernelAbstractions` is reserved for accelerator-friendly kernels | Simpler coarse-grained policy, less nested oversubscription, better fit for Julia runtime | Intentional |
| GPU strategy | Partial CuPy support exists, but multiple modules force `xp = np` or only lightly use GPU branches | Julia code is written around parametric array backends and backend traits | The Julia design is more uniform and better aligned with future CUDA/accelerator support, even if full runtime CUDA validation is still pending | Intentional |

## Core Optics Differences

| Area | Python behavior | Julia behavior | Why it differs | Status |
|---|---|---|---|---|
| Telescope constructor surface | `Telescope.py` supports user-defined pupil, pupil reflectivity, display of optical path, petal-related flags, coronagraph property, and internal optical-path bookkeeping | `src/Optics/telescope.jl` now supports user pupil and pupil-reflectivity maps directly; optical-path printing, petal flags, and coronagraph state are still not part of the current core surface | Julia prioritized the simulation-relevant optics surface first and left the UI/state-display layer out of core | Reduced gap |
| Pupil reflectivity / flux map | OOPAO tracks `pupilReflectivity` and source `fluxMap` directly on the coupled telescope/source path | Julia now exposes `set_pupil_reflectivity!` and `flux_map(tel, src)` while keeping source objects free of mutable telescope-coupled state | Same practical capability, with a cleaner separation between source state and telescope state | Equivalent |
| Off-axis behavior | OOPAO has `fov`, off-axis logic, and comments noting some of it is not fully implemented | Julia explicitly supports off-axis PSF shifts for single sources and asterisms | This is a Julian cleanup and a practical parity improvement rather than a mismatch | Intentional improvement |
| Optical-path tracing | OOPAO stores `optical_path` and can print it | Julia does not expose the same path-tracing UI | Logging/telemetry were preferred over mutable print-oriented state | Gap, low priority |

## Detector Differences

| Area | Python behavior | Julia behavior | Why it differs | Status |
|---|---|---|---|---|
| Detector physics surface | OOPAO detector models integration buffering, sensor type (`CCD/CMOS/EMCCD`), gain, ADC bits, FWC saturation, dark current, background noise/map, and digitalization | `src/Optics/detector.jl` now models sensor type, gain ordering, ADC quantization, FWC saturation, dark current, explicit integration buffering, background flux/map handling, QE, PSF sampling, binning, photon noise, and readout noise | Julia added the detector-electronics path without changing the low-latency capture surface | Equivalent except integer output type |
| Integration buffering | OOPAO accumulates frames until `integrationTime` is reached | Julia supports the same feature explicitly through `capture!(det, psf; sample_time=...)`, with `readout_ready(det)` exposing whether a completed readout is available | Julia keeps buffering explicit instead of coupling it implicitly to telescope relay semantics | Equivalent |
| Output precision / quantization | OOPAO supports `bits`, `output_precision`, `digitalization` | Julia models detector quantization levels via `bits`/`full_well`, but keeps the returned frame in the detector floating-point backend type instead of casting to an integer storage type | Preserves backend-generic arrays and avoids changing the runtime frame type | Reduced gap |
| Sensor-specific gain path | OOPAO distinguishes EMCCD vs CCD/CMOS behavior | Julia now uses typed sensor selectors `CCDSensor()` / `CMOSSensor()` / `EMCCDSensor()` to apply gain in the correct stage of the readout path | Same feature, more Julian API surface | Equivalent |
| Noise composition API | OOPAO uses mutable booleans and numeric properties (`photonNoise`, `readoutNoise`, etc.) | Julia uses typed noise models such as `NoisePhoton`, `NoiseReadout`, `NoisePhotonReadout` | More idiomatic Julia and easier dispatch | Intentional |

## Wavefront Sensor Differences

| Area | Python behavior | Julia behavior | Why it differs | Status |
|---|---|---|---|---|
| Sensing mode selection | OOPAO often uses booleans or mutable properties like `is_geometric` and `postProcessing` strings | Julia uses sensing-mode types and normalization types | Avoids symbol/string branching in hot code; more Julian API design | Intentional |
| Signal layout | OOPAO commonly exposes `signal`, `signal_2D`, and padded/full-grid representations | Julia now uses compact valid-pixel slope vectors for Pyramid/BioEdge and compact slope vectors generally in calibration paths | More efficient and cleaner for linear algebra; reference harness adapts at the comparison edge | Intentional |
| Shack-Hartmann parity knobs | OOPAO exposes half-pixel shift, pixel-scale tuning, thresholded CoG, convolution thresholding, etc. | Julia has explicit support for these parity knobs in `src/WFS/shack_hartmann.jl` | Needed to match OOPAO numerically | Equivalent |
| Pyramid/BioEdge masks | OOPAO has specific old/new mask variants, rooftop options, modulation path behavior, and detector-valid-pixel conventions | Julia implements these with typed structs and explicit helper functions instead of class-side mutable configuration | Same behavior target, different internal structure | Intentional |
| Pyramid/BioEdge modulation frame scaling | OOPAO implicitly applies FFT normalization in its focal-plane camera/frame path | Julia had a mismatch here and now matches OOPAO in `pyramid_modulation_frame!` | This was a real bug, not a design difference | Fixed parity issue |
| LGS sodium kernels for Pyramid/BioEdge | Main OOPAO surface does not appear to expose per-subaperture Na-profile kernels for Pyramid/BioEdge | Julia currently uses averaged Na-profile kernels there | Julia added a practical extension while keeping OOPAO parity tests focused on the surfaced Python behavior | Extension |
| WFS detector coupling | OOPAO WFS classes often own/compose detector-like objects directly (`cam`, focal-plane camera) | Julia keeps detector and WFS concepts more explicit and composable | Cleaner separation of state and easier backends | Intentional |

## Gain Sensing Camera Differences

| Area | Python behavior | Julia behavior | Why it differs | Status |
|---|---|---|---|---|
| Calibration workflow | OOPAO calibrates from focal-plane frames and prints progress | Julia mirrors the same optical-gain math but uses cached FFT workspaces and logging | Same algorithm, more reusable implementation | Intentional |
| Parallel FFT batching | OOPAO splits basis products into chunks using `n_jobs` loops | Julia uses thread-aware FFT workspaces in `split_basis_product` | Better fit for Julia threading and cached plan reuse | Intentional |
| Detector metadata | OOPAO GSC can carry detector-related properties through the focal-plane camera object | Julia `GainSensingCamera` can attach an explicit detector metadata snapshot without coupling detector physics into the optical-gain math | Keeps the optical-gain path clean while preserving the useful display/metadata surface | Equivalent |
| Long-horizon atmosphere replay | OOPAO tutorial `AO_closed_loop_Pyramid_WFS_GSC.py` exists as a full workflow | Julia matches compact traces, bounded replay, and the first nonlinear branch step, but treats the long huge-OPD replay as a diagnostic stress case rather than a normative parity gate | Both implementations enter a highly unstable nonlinear regime where tiny differences diverge into different trajectories | Diagnostic |

## LiFT Differences

| Area | Python behavior | Julia behavior | Why it differs | Status |
|---|---|---|---|---|
| Mode selection | OOPAO uses `numerical` boolean and string-like weighting control in `Reconstruct` | Julia uses `LiFTAnalytic` / `LiFTNumerical` types and weighting-mode dispatch | More idiomatic Julia and easier specialization | Intentional |
| Object convolution | OOPAO keeps `self.object` and convolves if present | Julia uses `object_kernel` in params and a typed `maybe_object_convolve!` path | Same feature, cleaner state model | Intentional |
| Interaction matrix | OOPAO and Julia both support analytic and numerical Jacobians | Julia regression-backs the interaction matrix path against OOPAO | Equivalent in validated cases | Equivalent |
| Iterative reconstruction weighting | OOPAO supports detector-driven weighting and iterative/model/static choices | Julia supports the same weighting modes and accepts detector-pixel inputs directly for binned data when `img_resolution` is specified in detector pixels | Both implementations treat LiFT image size in detector-pixel terms rather than modeling detector binning as a separate internal weighting stage | Equivalent |
| Solve path | OOPAO uses vectorized NumPy/CuPy reconstruction logic | Julia uses preallocated buffers and mutating linear algebra | Same mathematical target, different implementation | Intentional |

## Calibration / Reconstruction Differences

| Area | Python behavior | Julia behavior | Why it differs | Status |
|---|---|---|---|---|
| Calibration storage | OOPAO uses class-heavy calibration helpers such as `CalibrationVault` and modal-basis utilities with print-driven status output | Julia ports these as typed data + explicit functions in `src/Calibration/*.jl` | Better separation of operators, data, and state | Intentional |
| KL/NCPA default | OOPAO includes modal-basis and NCPA workflows, often around existing stored products and cockpit utilities | Julia supports HHt/PSD-based KL, but the default NCPA KL path is still DM-mode covariance (`MᵀM`) | The simpler default keeps the common DM-centric workflow lightweight while preserving the more physical option | Intentional, but worth revisiting |
| SPRINT I/O | OOPAO tooling often assumes local files and helper utilities | Julia SPRINT supports cached sensitivity matrices via `Serialization`, but deliberately does not bake FITS into core | Core package should stay file-format agnostic | Intentional |
| Stability policy | OOPAO often uses pseudoinverses and print-driven conditioning diagnostics | Julia currently does similar things in places; see `docs/numerical-stability-review.md` | This is not a Julia-vs-Python difference in philosophy yet; both can be improved | Shared weakness |

## Tomography Differences

| Area | Python behavior | Julia behavior | Why it differs | Status |
|---|---|---|---|---|
| Package boundary | OOPAO tomography tutorial depends on external `pyTomoAO` | Julia integrates tomography into `src/Tomography/` within the same package | One top-level module was preferred for now; users asked for a native Julia tomography subsystem | Intentional |
| Internal ordering | pyTomoAO logic naturally follows NumPy row-major conventions | Julia tomography is column-major internally; reference harness converts only at the comparison edge | Native layout should win internally | Intentional |
| API shape | pyTomoAO uses parameter classes and a reconstructor manager class | Julia uses typed parameter structs plus `ModelBasedTomography()` / `InteractionMatrixTomography()` dispatch | Better Julia ergonomics and more direct operator composition | Intentional |
| Numerical validation | pyTomoAO provides the reference implementation | Julia now regression-backs `Gamma`, `Cxx`, `Cox`, `Cnz`, reconstructors, wavefront maps, and KAPA command assembly | Same workflow target, different structure | Equivalent in committed reference cases |

## Tutorials and Workflow Surface

| Area | Python behavior | Julia behavior | Why it differs | Status |
|---|---|---|---|---|
| Main tutorials | OOPAO provides scripts and notebooks like `image_formation.py`, `how_to_detector.py`, `how_to_SPRINT.py`, `AO_transfer_function.py`, `AO_closed_loop_*`, `how_to_tomography.py` | Julia ports the main practical workflows under `examples/tutorials/` | This is the core parity path | Largely equivalent |
| Tutorial coverage | OOPAO also ships `demo_OOPAO.ipynb`, `ORP_AO_School_How_To.py`, and school/PAPYRUS/EKARUS material | Julia does not currently port the school/demo/display-heavy ecosystem one-to-one | Those are lower priority than core scientific workflows and rely heavily on plotting/tools helpers | Gap |
| Full GSC tutorial replay | OOPAO has the full atmosphere-driven Pyramid+GSC tutorial | Julia has compact and bounded regressions plus branch-step parity, but not a committed full long-horizon parity gate | The remaining mismatch is still under investigation | Gap |

## Tools and UI Differences

| Area | Python behavior | Julia behavior | Why it differs | Status |
|---|---|---|---|---|
| `tools/displayTools.py` | Present and widely used by tutorials | Not implemented in Julia core | Plotting is kept optional; display helpers should not drive the simulation design | Gap |
| `tools/OopaoGUI.py` | Present in OOPAO | Not implemented | GUI is out of scope for core parity | Gap |
| `tools/interpolateGeometricalTransformation.py` | Present and used by some tutorials / DM utilities | Not implemented as a standalone tools module | Julia favors domain-specific implementations where needed instead of porting the full tools layer wholesale | Gap |
| `tools/interpolate_influence_functions.py` | Present | Not implemented as a standalone tools module | Same reason | Gap |
| `tools/set_paralleling_setup.py` | Present and used in initialization helpers | Not implemented | Julia parallel policy is encoded in code structure and traits rather than a separate mutable setup helper | Intentional |
| `tools/tools.py` | Mixed helper bag used throughout OOPAO | Not ported as a single catch-all module | Julia intentionally avoids a monolithic utility bag in favor of typed APIs | Intentional |

## Known Fidelity-Relevant Gaps

These are the differences that still matter for parity work.

1. Detector output storage remains backend floating-point even when quantized,
   rather than casting to an integer output precision the way OOPAO can.
2. The OOPAO `tools/*` layer and display/GUI ecosystem are not ported.

## Diagnostic Stress Cases

These are useful comparison scenarios, but they are not good hard parity gates.

1. The full long-horizon atmosphere-driven `AO_closed_loop_Pyramid_WFS_GSC.py`
   replay enters a huge-OPD nonlinear regime in which OOPAO itself becomes
   extremely sensitive. Julia matches the compact traces, bounded replay, and
   first branch-step intermediates; the remaining long-horizon divergence is
   treated as a robustness/diagnostic issue rather than a missing core feature.

## Differences That Are Not Bugs

These differences should not be “fixed” unless they block a validated workflow:

1. Multiple dispatch instead of Python class methods.
2. Column-major internals in tomography and other arrays.
3. Parametric array fields and backend traits instead of hard-coded NumPy arrays.
4. Structured errors and logging instead of print-heavy side effects.
5. Explicit mutating `!` APIs and preallocated workspaces.
6. Keeping FITS and plotting out of the simulation core.

## Recommended Use of This Document

When continuing parity work:

1. Treat the `Gap` rows as the backlog.
2. Treat the `Intentional` rows as settled design unless they are shown to block
   a reference workflow.
3. Treat the `Extension` rows as Julia-specific optional capabilities, not as
   parity blockers.

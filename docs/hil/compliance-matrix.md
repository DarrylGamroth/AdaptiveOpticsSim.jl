# HIL Feature And Compliance Matrix

Status: active

## Purpose

This is the stable feature, completion, and compliance ledger for the HIL
architecture. PRs and issues cite these IDs and update the affected rows.
Temporary task checklists remain in the PR or issue rather than becoming new
documentation files.

See the [`HIL architecture index`](../hil-package-boundary.md) for the
subsystem specifications governed by this matrix.

Implementation state describes code availability:

- **implemented:** the complete generic requirement exists on at least one
  maintained path
- **partial:** useful foundations exist, but the stated requirement is
  incomplete
- **planned:** the requirement exists only in the target design

Evidence state describes compliance:

- **validated:** maintained tests and required hardware evidence satisfy the
  complete row
- **partial:** some correctness or backend evidence exists
- **missing:** required acceptance evidence has not been produced

A row is complete only when implementation and evidence both satisfy its full
requirement. A numbered gate is complete only when all rows assigned to it and
its gate acceptance are complete. Validated rows link their tests, benchmark
artifacts, or support entry in Notes.

## Traceability Rules

Requirement IDs are stable and are never renumbered or reused for a different
contract. The subsystem specifications elaborate the normative behavior; this
matrix is the canonical state ledger:

| Requirement family | Normative specification |
|---|---|
| `HIL-PLANE`, `HIL-ATM`, `HIL-SRC`, `HIL-RAD`, `HIL-SCI`, `HIL-API` | [`plant-and-optics.md`](plant-and-optics.md), [`package-boundaries.md`](package-boundaries.md), and [`validation.md`](validation.md) |
| `HIL-WFS` | [`plant-and-optics.md`](plant-and-optics.md) and [`validation.md`](validation.md) |
| `HIL-ORACLE`, `HIL-RNG`, `HIL-OBS`, `HIL-REPLAY` | [`../deterministic-simulation.md`](../deterministic-simulation.md), [`time-and-scheduling.md`](time-and-scheduling.md), [`execution-and-placement.md`](execution-and-placement.md), and [`validation.md`](validation.md) |
| `HIL-PLANT`, `HIL-PATH`, `HIL-DET`, `HIL-FID` | [`plant-and-optics.md`](plant-and-optics.md), [`time-and-scheduling.md`](time-and-scheduling.md), and [`validation.md`](validation.md) |
| `HIL-CAL` | [`plant-and-optics.md`](plant-and-optics.md), [`package-boundaries.md`](package-boundaries.md), and [`validation.md`](validation.md) |
| `HIL-CMD`, `HIL-OPT` | [`plant-and-optics.md`](plant-and-optics.md) and [`time-and-scheduling.md`](time-and-scheduling.md) |
| `HIL-TIME`, `HIL-TRIG`, `HIL-SCHED`, `HIL-LIFE` | [`time-and-scheduling.md`](time-and-scheduling.md) |
| `HIL-VSLICE`, `HIL-PORT`, `HIL-RING`, `HIL-BUF`, `HIL-FAIL` | [`rtc-ports.md`](rtc-ports.md), [`time-and-scheduling.md`](time-and-scheduling.md), and [`validation.md`](validation.md) |
| `HIL-CPU`, `HIL-GPU`, `HIL-EXEC`, `HIL-DIST` | [`execution-and-placement.md`](execution-and-placement.md) |
| `HIL-PROFILE-*` | [`validation.md`](validation.md) |

An `implemented` or `partial` implementation row links existing code in Notes.
A `validated` or `partial` evidence row links a test, artifact, or maintained
support entry. `planned` and `missing` are explicit gaps, not implied coverage.
PRs update the links and states together when behavior or evidence changes.

## Feature Matrix

| ID | Gate | Requirement | Owner | Implementation | Evidence | Notes |
|---|---:|---|---|---|---|---|
| HIL-PLANE-001 | 0 | A prepared telescope owns run-immutable aperture, reflectivity, spatial sampling, and geometry while each path owns its mutable pupil OPD/electric field and focal-plane result; telescope owns no temporal cadence or exposure duration | Core | implemented | validated | [`TelescopeAperture`](../../src/optics/telescope.jl) owns revisioned geometry and intensity reflectivity. Independent [`PupilFunction`](../../src/optics/planes.jl) products snapshot support and field amplitude, own mutable path OPD, and retain the aperture revision against which they were prepared. Maintained atmosphere rendering, controllable-optic application, WFS families, calibration, direct science, and runtime WFS/science/auxiliary arms consume explicit path products; `Telescope` has no mutable path or focal state. [`Explicit optical products and surfaces`](../../test/testsets/core_optics.jl), [`wfs_stage_contracts.jl`](../../test/testsets/wfs_stage_contracts.jl), and the runtime path-isolation cases in [`control_and_runtime.jl`](../../test/testsets/control_and_runtime.jl) verify independent products, surface reuse, aperture invalidation, and science/WFS path separation. |
| HIL-PLANE-002 | 0 | Caller-owned optical-plane products are separate from single-writer prepared propagation/spatial-filter workspaces, FFT plans, and scratch, with fixed prepared shape, numeric type, backend, and device | Core | implemented | validated | Caller-owned [`ElectricField` and `IntensityMap`](../../src/optics/planes.jl) products are separate from [`FraunhoferPropagation`/`FresnelPropagation`](../../src/optics/propagation.jl), [`DirectImagingWorkspace`](../../src/optics/direct_imaging.jl), and [`SpatialFilterWorkspace`](../../src/optics/spatial_filter.jl). Existing explicit-product tests enforce fixed preparation and zero warmed CPU allocation. Frozen Gate 0 references preserve historical numerical characterization; the direct-science comparison applies an explicit test-only adapter for the legacy x/y orientation defect, while [`direct_science.jl`](../../test/testsets/direct_science.jl) independently validates the corrected declared-axis semantics. The full CPU suite and maintained AMDGPU and CUDA hardware targets validate the prepared direct-imaging workspace/output boundary on this revision with scalar indexing disabled. |
| HIL-PLANE-003 | 0 | Static aberration and controllable-optic surfaces apply to explicit caller-owned path products without mutating shared telescope aperture state | Core | implemented | validated | [`apply_surface!`](../../src/optics/planes.jl) accepts OPD maps, NCPA, DMs, and modal optics after distinct controllable-surface formation via [`update_surface!`](../../src/optics/deformable_mirror.jl); [`Explicit optical products and surfaces`](../../test/testsets/core_optics.jl) verifies independent path products and unchanged telescope aperture state, and the optional hardware matrix verifies device-resident DM application. |
| HIL-PLANE-004 | 0 | Every prepared optical product declares plane kind, coordinate domain, dimensions, sampling, origin/centering, axis orientation, wavelength/channel, units/normalization, density-versus-cell-integrated measure, coherence/combination policy, numeric type, backend, and physical device where applicable; incompatible handoffs fail before repeated execution | Core | implemented | validated | [`OpticalPlaneMetadata`](../../src/optics/planes.jl) now includes concrete photon-rate/dimensionless normalization, point-sampled/spatial-density/cell-integrated measure, and coherent/incoherent/non-combinable policy markers in addition to the structural, spectral, backend, and physical-device contract. Achromatic geometry, monochromatic products, and application-declared integrated channels are distinct; unspecified spectral metadata is rejected at prepared intensity-combination and detector-acquisition boundaries. Cold propagation, compatible accumulation, and [`prepare_detector_acquisition`](../../src/detectors/detector_acquisition.jl) reject incompatible metadata before repeated execution. [`Optical radiometry and combination contracts`](../../test/testsets/core_optics.jl), [`Detector`](../../test/testsets/detectors.jl), and the [optional backend matrix](../../test/backend_optional_common.jl) provide focused CPU and backend evidence; the full CPU suite and maintained AMDGPU and CUDA hardware matrices validate this revision. Names and units follow the maintained [`glossary`](../glossary.md). |
| HIL-ATM-001 | 0 | One writer advances a shared atmosphere to a stable timed epoch; prepared path-local direction renderers consume that epoch in any order and write caller-owned outputs without shared mutable geometry caches or RNG mutation | Core | implemented | validated | [`source_geometry.jl`](../../src/atmosphere/source_geometry.jl) defines identity-bound epochs and prepared direction renderers; [`Explicit atmosphere epochs and prepared renderers`](../../test/testsets/atmosphere.jl) verifies render-order invariance, stale/incompatible rejection before mutation, RNG preservation, NGS/LGS geometry, and zero warmed CPU allocation. [`control_and_runtime.jl`](../../test/testsets/control_and_runtime.jl) verifies one advance with independent WFS/science and shared-arm renderers. The [`optional backend matrix`](../../test/backend_optional_common.jl) executes the public epoch/renderer and field APIs directly on the maintained AMDGPU and CUDA hardware targets. |
| HIL-ATM-002 | 0 | Atmosphere advancement accepts explicit model time or elapsed duration and never infers time from telescope sampling, detector cadence, wall clock, or a renderer | Core | implemented | validated | `advance_by!` and `advance_to!` in [`source_geometry.jl`](../../src/atmosphere/source_geometry.jl) are the maintained timed-model API. [`ClosedLoopRuntime`](../../src/control/runtime/types.jl) and single/grouped control-loop configurations require a finite positive `atmosphere_step`, preserve it across reconstruction changes, and pass it explicitly to one sensing-cycle atmosphere advance; detector exposure remains independent. Analytic wind-offset, monotonic-time, zero-duration, fixed-seed replay, shared-arm, invalid-step, and exact-runtime-step tests in [`atmosphere.jl`](../../test/testsets/atmosphere.jl) and [`control_and_runtime.jl`](../../test/testsets/control_and_runtime.jl) provide focused CPU evidence; the [optional backend matrix](../../test/backend_optional_common.jl) exercises explicit timed advancement on maintained AMDGPU and CUDA hardware targets. |
| HIL-SRC-001 | 0 | Source definitions are run-immutable execution descriptions; mutable profiles/images are frozen during preparation, and direction, spectral, LGS, asterism, and extended-source expansion are path-local prepared state | Core | implemented | validated | [`source.jl`](../../src/optics/source.jl), [`spectrum.jl`](../../src/optics/spectrum.jl), [`asterism.jl`](../../src/optics/asterism.jl), and [`extended_source.jl`](../../src/optics/extended_source.jl) freeze run-owned descriptions. [`atmosphere.jl`](../../test/testsets/atmosphere.jl) mutates original LGS profiles, spectral arrays, and extended-source images after preparation and verifies unchanged prepared execution; plural renderer preparation expands multi-direction sources without source-owned caches. |
| HIL-RAD-001 | 0 | Physical detector-facing optical products are photon-arrival-rate products with declared units and density-versus-cell-integrated measure; telescope and optical front ends apply no elapsed-time factor, while detector acquisition integrates its explicit whole-exposure or incremental duration exactly once and never creates photons through declared spatial response and pixel integration; finite-support loss is explicit | Core | implemented | validated | Physical and normalized source radiometry in [`source.jl`](../../src/optics/source.jl) propagates through rate- or dimensionless-normalized fields and intensity maps without a telescope-time factor. [`DetectorAcquisitionPlan`](../../src/detectors/detector_acquisition.jl) validates a detector-facing map and prepares density-to-cell scaling and declared-channel QE; repeated whole or incremental capture applies presampling response before physical-pixel integration and applies explicit exposure exactly once. Prepared Shack–Hartmann, Pyramid, BioEdge, Zernike, and Curvature front ends write cell-integrated detector-plane photon-rate products without elapsed time. Pyramid/BioEdge modulation weights are normalized optical quadrature weights, not exposure fractions. Curvature's two branch rate planes remain unchanged while separate or packed acquisition applies detector response, QE, and detector-owned duration. The prepared LiFT forward model likewise publishes a cell-integrated focal-plane photon rate; `LiFTExpectedCounts` and `LiFTNormalizedIntensity` make exposure, QE, and native-value conversion explicit only at the independent observation boundary. Distinct wavelength and path-local source products remain an `OpticalProductBundle`. Finite response uses zero extension, so signal leaving detector support is lost rather than clamped back into edge samples. [`core_optics.jl`](../../test/testsets/core_optics.jl), [`detectors.jl`](../../test/testsets/detectors.jl), the prepared [`WFS stage contracts`](../../test/testsets/wfs_stage_contracts.jl), [`shack_hartmann_and_sources.jl`](../../test/testsets/shack_hartmann_and_sources.jl), [`pyramid_bioedge_and_lgs.jl`](../../test/testsets/pyramid_bioedge_and_lgs.jl), [`zernike_and_curvature.jl`](../../test/testsets/zernike_and_curvature.jl), [`calibration_and_analysis.jl`](../../test/testsets/calibration_and_analysis.jl), and the [optional backend matrix](../../test/backend_optional_common.jl) cover unequal and non-unit exposures, density and cell-integrated fixtures, non-commuting response order, spectral/source bundles, modulation normalization, component accumulation, edge behavior, sampled-QE selection, prepared GPU acquisition, LiFT rate/count equivalence, and warmed allocation. The staged physical-family fixtures prove that optical rate products remain unchanged while QE and duration are applied once in acquisition. |
| HIL-RAD-002 | 0 | Spectral/source planes are accumulated elementwise only with compatible physical grids and declared incoherent intensity/rate semantics; otherwise a typed bundle or explicit prepared mapping is required | Core | implemented | validated | [`prepare_incoherent_sum`](../../src/optics/planes.jl) validates immutable metadata once and [`accumulate_intensity!`](../../src/optics/planes.jl) performs prepared allocation-free accumulation; [`OpticalProductBundle`](../../src/optics/planes.jl) preserves incompatible products without implicit conversion or resampling. Unspecified spectral coordinates cannot authorize an intensity sum. Prepared Shack–Hartmann, Pyramid, and BioEdge formation retains distinct wavelength products as bundles. Pyramid/BioEdge Asterism and ExtendedSource formation additionally requires one explicit path-rendered pupil per source and returns a bundle, preventing direction-dependent pupils from being combined by array index. Single-plane and count-mismatched requests fail during preparation. [`core_optics.jl`](../../test/testsets/core_optics.jl) covers compatible accumulation and rejection, while [`wfs_stage_contracts.jl`](../../test/testsets/wfs_stage_contracts.jl) validates wavelength weights, bundle retention, path-local inputs, and structural rejection. |
| HIL-SCI-001 | 0 | Native or prepared external science optics write a caller-owned focal-plane photon-arrival-rate product that one or more independent detector acquisitions integrate without telescope-owned PSF or cadence state | Core | implemented | validated | [`prepare_direct_imaging`](../../src/optics/direct_imaging.jl) binds an explicit caller-owned pupil or preformed field, work field, focal-plane angular `IntensityMap`, and single-writer workspace; [`form_direct_image!`](../../src/optics/direct_imaging.jl) writes a cell-integrated photon-arrival-rate image without exposure time or telescope-owned focal storage. Preparation resolves a finite integer off-axis displacement using the declared focal-grid axis order and signs. Same-grid asterism components are added incoherently, differing spectral grids remain an `OpticalProductBundle`, and extended sources use explicit `extended_source_asterism` expansion. [`direct_science.jl`](../../test/testsets/direct_science.jl) covers point/off-axis orientation, repeated formation after bound-pupil OPD changes, finite mixed-precision radiometry, fan-out to unequal detector exposures, spectral bundling, declared physical-rate and explicitly scaled normalized external-result test doubles, prepared rejection, and warmed CPU allocation. Shared runtimes form one rate image per arm and reuse its exact storage across independent detectors. Focused CPU tests and the maintained AMDGPU and CUDA hardware targets validate point/off-axis formation, spectral bundles, explicit extended-source expansion, detector fan-out, backend/device rejection, and shared-runtime residency with scalar indexing disabled. The committed [Gate 0 artifact](../../benchmarks/results/gate0/2026-07-16-pre-hil-06-science-stage.toml) records zero steady-state allocation and passing predeclared absolute and relative p99 gates for `G0-PERF-02` and `G0-PERF-07`. Asynchronous `Proper.jl` execution remains tracked by `HIL-PATH-004`. |
| HIL-WFS-001 | 0 | Deterministic characterization oracles freeze the current optical, detector-coupled, and estimated products of every maintained WFS family before decomposition | Core | implemented | validated | Existing [reference data](../../test/reference_data), [SPECULA-oriented reference data](../../test/reference_data_specula), and WFS tests cover complete products. [`zernike_and_curvature.jl`](../../test/testsets/zernike_and_curvature.jl) and the staged contract fixtures preserve the Zernike/Curvature optical, acquisition, normalization, and branch-order oracles through decomposition. [`calibration_and_analysis.jl`](../../test/testsets/calibration_and_analysis.jl) adds independently owned LiFT rate/count/normalized observation fixtures, analytic/numerical interaction-matrix norm oracles, detector-ownership rejection, and coefficient-recovery diagnostics while the maintained OOPAO reference interaction matrix remains in [`reference_data`](../../test/reference_data). |
| HIL-WFS-002 | 0 | WFS implementations compose prepared optical-front-end, detector-acquisition, and estimator stages through dispatch over caller-owned products without HIL clock, scheduling, queue, or transport dependencies | Core | implemented | validated | [`stage_contracts.jl`](../../src/wfs/stage_contracts.jl) defines typed caller-owned observations/measurements and six narrow prepared dispatch functions for optical formation, acquisition, and estimation, plus explicit acquired/direct measurement traits and structured preparation errors. Shack–Hartmann [`stages.jl`](../../src/wfs/shack_hartmann/stages.jl), Pyramid [`stages.jl`](../../src/wfs/pyramid/stages.jl), BioEdge [`stages.jl`](../../src/wfs/bioedge/stages.jl), Zernike [`stages.jl`](../../src/wfs/zernike/stages.jl), and Curvature [`stages.jl`](../../src/wfs/curvature/stages.jl) bind caller-owned pupil/field, rate product or bundle/tuple, observation, and measurement storage without clock, queue, transport, or runtime-graph dependencies. LiFT deliberately uses its own prepared forward/observation/estimator dispatch rather than inheriting `AbstractWFS`; it likewise owns no clock, queue, transport, or detector. [`wfs_stage_contracts.jl`](../../test/testsets/wfs_stage_contracts.jl) and [`calibration_and_analysis.jl`](../../test/testsets/calibration_and_analysis.jl) cover exact binding, stage independence, multiple planes and detectors, direct paths, LiFT ownership, and warmed allocation. The full CPU suite and maintained `412/412` CUDA and `422/422` AMDGPU hardware targets validate every composed family on this revision with scalar indexing disabled. |
| HIL-WFS-003 | 0 | Optical front ends produce one or more detector-plane photon-arrival-rate products with an explicit spatial measure; presampling response, physical-pixel integration, QE and explicit elapsed-time integration, coupling, stochastic response, and readout follow optical formation and precede estimation, with one or multiple detector observations supported; the realized response's interior MTF remains available as a derived diagnostic | Core | implemented | validated | The generic [`stage contract`](../../src/wfs/stage_contracts.jl) requires detector-plane photon-arrival-rate products and supports static one/multiple-plane and one/multiple-observation composition. Prepared Shack–Hartmann, Pyramid, BioEdge, Zernike, and Curvature implementations form cell-integrated rate planes or typed bundles/tuples before detector acquisition. Curvature may feed two independent detectors with unequal exposure and branch-specific response/MTF, or one explicitly packed detector with compatible branch metadata and a shared exposure. LiFT's separately prepared focal-plane model forms a cell-integrated photon rate, while `LiFTFrameMapping` binds deterministic response/sampling/binning and the independent observation domain binds QE/exposure or normalized scaling before estimation. [`wfs_stage_contracts.jl`](../../test/testsets/wfs_stage_contracts.jl) proves unchanged optical rates, exact-once QE/duration, non-commuting branch-specific response ordering, and separate/packed estimator equivalence; [`calibration_and_analysis.jl`](../../test/testsets/calibration_and_analysis.jl) proves LiFT rate/count equivalence and preprocessing mismatch rejection. [`detectors.jl`](../../test/testsets/detectors.jl) verifies response, pixel integration, QE, exposure, readout, and detector MTF diagnostics. The full CPU suite and maintained `412/412` CUDA and `422/422` AMDGPU hardware targets validate the staged optical/acquisition ordering, detector-facing products, and maintained MTF diagnostics on this revision. |
| HIL-WFS-004 | 0 | Shack-Hartmann composes an independent microlens-array optical primitive, detector acquisition, and replaceable spot estimator while preserving geometric, diffractive, spectral, LGS, and asterism behavior | Core | implemented | validated | [`MicrolensArrayParams` and `MicrolensArray`](../../src/wfs/shack_hartmann/setup.jl) describe the immutable regular-array model and sampling policy; `prepare_microlens_propagation` creates separate backend/grid-bound `PreparedMicrolensPropagation` execution state. `ShackHartmannWFS.front_end` is the actual component owner: `ShackHartmannDirectFrontEnd` carries only microlens geometry/layout for direct sensing, while `ShackHartmannOpticalFrontEnd` composes the microlens model, propagation, and `SubapertureLayout` for diffractive sensing. Calibration, acquisition, and estimator storage remain separate; superseded top-level optical fields and the whole-WFS optical-owner union are removed rather than forwarded. [`stages.jl`](../../src/wfs/shack_hartmann/stages.jl) consumes a caller-owned pupil/field, writes photon-rate mosaics or spectral bundles without telescope mutation or cadence, rejects incompatible plane semantics and placement, accepts detector-declared real observation units, and requires a revision-bound explicit calibration before acquired estimation. Geometric and centroid signals use the same `[axis 1; axis 2]`, Julia-column-major lenslet convention; the OOPAO row-major/axis adapter is confined to the reference harness. [`wfs_stage_contracts.jl`](../../test/testsets/wfs_stage_contracts.jl) exercises the component-only path and checks frozen diffractive output, asymmetric axis/lenslet ordering, geometric direct measurement, plane-contract rejection, LGS elongation and sodium profiles, asterism accumulation, spectral bundles, response ordering, exact-once exposure, mixed precision, calibration invalidation, caller-owned mutation, and zero warmed CPU allocation. The optional backend matrix checks device-resident execution and cross-backend rejection with scalar indexing disabled. |
| HIL-WFS-005 | 0 | Pyramid and BioEdge compose independent physical focal-plane front ends, prepared modulation policies, detector acquisition, and calibrated differential estimators while sharing only semantically common machinery | Core | implemented | validated | [`PyramidOpticalFrontEnd`](../../src/wfs/pyramid/setup.jl) binds a physical unit-magnitude phase mask, while [`BioEdgeOpticalFrontEnd`](../../src/wfs/bioedge/setup.jl) binds four distinct amplitude masks; their prepared stage implementations remain separate in [`pyramid/stages.jl`](../../src/wfs/pyramid/stages.jl) and [`bioedge/stages.jl`](../../src/wfs/bioedge/stages.jl). Shared [`focal_plane_modulation.jl`](../../src/wfs/focal_plane_modulation.jl) is limited to normalized no-, circular-, and user-sampled focal-plane modulation policies and path/source helpers. The broader legacy calibration modulation selects support only; zero-aberration references use the operating modulation. Both families form exposure-independent four-pupil photon-rate products or wavelength/path-local bundles, then use the family-neutral detector acquisition and calibrated differential estimators that accept real floating-point or integer mosaics. Samples are converted to estimator precision before differential arithmetic. Backend/device, propagation-sampling, and calibration-revision bindings fail before output mutation; geometry-preserving detector reduction remains supported. Geometric policies intentionally use a direct floating-point measurement path with no focal-plane or detector workspace. [`pyramid_bioedge_and_lgs.jl`](../../test/testsets/pyramid_bioedge_and_lgs.jl) preserves the frozen legacy oracle. [`wfs_stage_contracts.jl`](../../test/testsets/wfs_stage_contracts.jl) validates distinct mask physics, zero/circular/user modulation, exact-once non-unit exposure and QE, floating/integer estimator equivalence, malformed and misplaced storage rejection, detector reduction, atomic calibration and propagation-sampling invalidation, spectral and path-local source bundles including heterogeneous NGS/LGS paths, simple elongation and sodium-profile LGS behavior, direct geometric execution, and zero warmed CPU allocation. The [optional backend matrix](../../test/backend_optional_common.jl) exercises both physical families, including integer detector frames, on maintained CUDA and AMDGPU targets with scalar indexing disabled. |
| HIL-WFS-006 | 0 | Zernike and Curvature compose independent optical front ends, detector acquisition, and estimators; Curvature supports both separate branch detectors and packed single-detector observations | Core | implemented | validated | [Zernike](../../src/wfs/zernike.jl) now separates `ZernikePhaseSpot` plus prepared pupil-relay propagation, detector acquisition state, and a revision-bound referenced pupil estimator. [Curvature](../../src/wfs/curvature.jl) separates `CurvatureDefocusPair` plus prepared two-branch propagation, detector acquisition state, and a revision-bound differential estimator. Its ordered positive-/negative-defocus rate tuple can feed two independent detectors with unequal durations and branch-specific response, or `CurvaturePackedAcquisition` can map compatible planes into one frame-region or counting-channel observation with a shared detector duration. Estimation accepts either mapping and explicit deterministic branch rate scales. [`zernike_and_curvature.jl`](../../test/testsets/zernike_and_curvature.jl) preserves the pre-decomposition optical and normalization oracle. [`wfs_stage_contracts.jl`](../../test/testsets/wfs_stage_contracts.jl) validates PupilFunction/ElectricField parity, caller-owned products, separation of mutable owners, floating/integer observations, response ordering, separate/packed frame and counting acquisition, duration/radiometry/geometry rejection, calibration invalidation, and zero warmed CPU allocation. The [optional backend matrix](../../test/backend_optional_common.jl) executes these paths device-resident with scalar indexing disabled. |
| HIL-WFS-007 | 0 | LiFT consumes an independently acquired observation and separates its focal-plane forward model from iterative estimation without being forced into the ordinary WFS type hierarchy | Core | implemented | validated | [`PreparedLiFTForwardModel`](../../src/wfs/lift.jl) freezes the pupil transmission, modal basis, diversity, source photon irradiance, object kernel, deterministic frame mapping, geometry/radiometry, backend, and device without retaining a telescope, source, or detector. `LiFTObservation` binds caller-owned photon-rate, expected-count, or normalized-intensity values to that contract; `LiFT` owns only the prepared forward model, prepared modal subset, solve policy, and iterative workspace. [`calibration_and_analysis.jl`](../../test/testsets/calibration_and_analysis.jl) covers independent observation storage, analytic/numerical Jacobians, count/rate equivalence, deterministic preprocessing, mismatch rejection, frozen interaction-matrix norms, and zero warmed allocation for the normal-equation path. The shared optional matrix and maintained `412/412` CUDA and `422/422` AMDGPU hardware targets now validate device-resident forward formation, dense/separable convolution, analytic interaction matrices, rate/count/normalized observations, and analytic/numerical reconstruction with scalar indexing disabled. |
| HIL-WFS-008 | 0 | Geometric and reduced-order WFS policies avoid unused diffractive workspace and declare whether they produce an approximate photon-arrival-rate product for detector processing or intentionally produce a direct measurement | Core | implemented | validated | `AcquiredObservationPath()` and `DirectMeasurementPath()` establish the acquired-observation versus intentional-direct-measurement distinction, and every maintained prepared estimator exposes that decision through `wfs_measurement_path`. Geometric Shack–Hartmann, Pyramid, and BioEdge preparation creates neither diffractive-front-end nor detector-acquisition workspace and intentionally returns direct measurements; maintained diffractive policies use acquired photon-rate observations. [`wfs_stage_contracts.jl`](../../test/testsets/wfs_stage_contracts.jl) verifies the cross-family declarations, absence of placeholder workspace, explicit-pupil execution, and zero warmed CPU allocation. New reduced-order policies must choose one of these contracts rather than silently imitating detector data. |
| HIL-WFS-009 | 0 | Every migrated stage preserves deterministic CPU results, maintained CUDA and AMDGPU parity, device residency, warmed allocation budgets, and comparable steady-state latency within predeclared absolute and relative gates | Core | implemented | validated | Prepared Shack–Hartmann, Pyramid, BioEdge, Zernike, and Curvature optical, acquisition, and estimation stages preserve their frozen CPU behavior and are allocation-free after warmup; geometric construction omits unused diffractive and acquisition storage. LiFT now preserves its frozen analytic/numerical interaction-matrix norms, uses device-safe buffer and object-convolution kernels, and is zero-allocation after warmup for analytic or numerical Jacobian construction and in-place normal-equation reconstruction. CPU `LiFTSolveQR` remains the accuracy-first automatic choice and has a small factorization-wrapper allocation; the explicitly prepared normal-equation path is the zero-allocation low-latency surface. The shared [optional backend matrix](../../test/backend_optional_common.jl) directly executes all physical families plus LiFT and its dense/separable object kernels on device-resident CUDA and AMDGPU arrays with scalar indexing disabled. The committed [Pre-HIL 7 Gate 0 artifact](../../benchmarks/results/gate0/2026-07-16-pre-hil-07-shack-hartmann-stage.toml) records passing absolute and relative latency gates for `G0-PERF-05`. The [Pre-HIL 8 Gate 0 artifact](../../benchmarks/results/gate0/2026-07-17-pre-hil-08-pyramid-bioedge-stage.toml) records zero warmed allocation and passing absolute and relative gates for staged Pyramid `G0-PERF-06`. The [Pre-HIL 9 Gate 0 artifact](../../benchmarks/results/gate0/2026-07-17-pre-hil-09-zernike-curvature-stage.toml) adds predeclared zero-allocation and absolute p99 gates for staged Zernike `G0-PERF-09` and two-plane/two-detector Curvature `G0-PERF-10`; these cards have no predecessor measurement, so their relative gates begin with this artifact rather than being retroactively evaluated. `G0-PERF-11` adds the prepared LiFT reconstruction contract and begins its relative baseline with this migration. The complete maintained hardware pass is `386/386` on CUDA and `396/396` on AMDGPU. The [full final Gate 0 catalog](../../benchmarks/results/gate0/2026-07-18-pre-hil-11-full-backend-evidence.toml) preserves one 96 ns Shack-Hartmann relative-p99 miss while every absolute and allocation gate passes; the immediate [same-contract confirmation](../../benchmarks/results/gate0/2026-07-18-pre-hil-11-shack-hartmann-confirmation.toml) passes at 108.543 μs median p99 against the 130.399 μs limit. Clean [local CPU](../../benchmarks/results/platform/2026-07-18-pre-hil-11-local-cpu.toml), [WSL CPU](../../benchmarks/results/platform/2026-07-18-pre-hil-11-wsl-cpu.toml), and [WSL CUDA](../../benchmarks/results/platform/2026-07-18-pre-hil-11-wsl-cuda.toml) service-time artifacts pass their declared correctness, residency, allocation, absolute-p95, and relative-p95 gates; AMD performance remains scoped to the earlier maintained artifact because the later host Julia installation failed before a new raw histogram was retained. |
| HIL-ORACLE-001 | 1 | Single-threaded deterministic execution is the numerical and event-order oracle for every capability implemented at the current gate; future multi-rate and placed-optic cases acquire their oracle when their defining gate is implemented | Core | implemented | validated | [`reference_harness.jl`](../../test/reference_harness.jl), the [determinism policy](../deterministic-simulation.md), the Gate 0 characterization suite, and maintained CPU/backend evidence establish the current numerical oracle without making later-gate behavior a Gate 1 prerequisite. |
| HIL-PLANT-001 | 2 | One coherent telescope and current atmosphere epoch token feeds every due path, and every atmosphere-dependent path input is materialized under a valid lifetime before mutable atmosphere state advances | Core | implemented | validated | [`preparation.jl`](../../src/plant/preparation.jl) prepares a canonical acquisition selection, resolves the unique paths required by selected full-optical providers, preflights all owner/revision/epoch bindings before mutation, and materializes every such path input from one current epoch before forming any result. [`plant_preparation.jl`](../../test/testsets/plant_preparation.jl) covers NGS, finite-height LGS, and science paths, declaration- and selection-order invariance, caller-owned input lifetime after atmosphere advance, stale/foreign epoch and revision rejection without partial output, path-result fan-out, inference, and zero warmed CPU allocation; [`plant_providers.jl`](../../test/testsets/plant_providers.jl) proves reduced-order and synthetic/replay-only selections do not form an unused optical path. The [Gate 2 serial plant artifact](../../benchmarks/results/gate2/2026-07-21-serial-plant.toml) composes three directions under explicit absolute model time. The due selection is caller-supplied until Gate 3 produces it; an `AtmosphereEpoch` remains explicitly not a retained snapshot. |
| HIL-PATH-001 | 2 | Prepared paths own mutable propagation workspace; acquisitions separately own WFS, detector, readout, observation, measurement, and caller-owned product state | Core | implemented | validated | [`preparation.jl`](../../src/plant/preparation.jl) freezes each source and builds one concrete `PreparedPathExecutor` with an exact plant-atmosphere binding, path-local input/result, prepared materialization operation, and optical execution owner; `PreparedAcquisitionOwner` values borrow the exact result read-only while owning one independent provider plus its detector/readout, optional WFS estimator, observation, and measurement state. A prepared selection deduplicates only paths required by full-optical providers, so compatible full-optical acquisitions consume one formed result while reduced-order and synthetic/replay providers bypass unused propagation. [`plant_preparation.jl`](../../test/testsets/plant_preparation.jl), [`plant_providers.jl`](../../test/testsets/plant_providers.jl), and the [Gate 2 serial plant artifact](../../benchmarks/results/gate2/2026-07-21-serial-plant.toml) exercise direct-science and diffractive WFS owners, NGS/LGS/science materialization, provider-specific path selection, exact fan-out, inference, rejection atomicity, and zero warmed CPU allocation. The optional backend matrix exercises the full-optical boundary device-resident with scalar indexing disabled. HIL publication ownership correctly remains outside this core requirement. |
| HIL-PATH-002 | 2 | Immutable optical-path definitions are separate from acquisition state that later schedules or triggers may invoke independently | Core | implemented | validated | [`definitions.jl`](../../src/plant/definitions.jl) supplies immutable identified declarations, while [`preparation.jl`](../../src/plant/preparation.jl) constructs separate concrete path and acquisition owners without a schedule, queue, or transport. A separate plant-owned RNG layer binds independent provider and detector streams to those exact owners without adding cadence or trigger state. `PathResultKey` and `require_path_result` cover source geometry, spectrum, radiometry, optical/propagation model, instantaneous-sample semantics, output plane, revisions, backend, and physical device before acquisition destination construction. [`plant_topology.jl`](../../test/testsets/plant_topology.jl), [`plant_preparation.jl`](../../test/testsets/plant_preparation.jl), [`plant_rng.jl`](../../test/testsets/plant_rng.jl), and the [Gate 2 serial plant artifact](../../benchmarks/results/gate2/2026-07-21-serial-plant.toml) validate declaration identity, reusable result fan-out, distinct acquisition state, owner-bound streams, and independent invocation. Gate 3 adds schedules around these owners rather than changing this separation. |
| HIL-API-001 | 2 | The core plant API separates run-immutable definitions, prepared plans, mutable single-writer state/workspaces, and caller-owned products without extending transitional OOPAO/frame-step abstractions or imposing HIL-only names on core models | Core | implemented | validated | The HIL-neutral [`definitions.jl`](../../src/plant/definitions.jl) boundary is paired with concrete, parametric [`PreparedPlant`, `PreparedPathExecutor`, `PreparedAcquisitionOwner`, and `PreparedAcquisitionSelection`](../../src/plant/preparation.jl), the product/provider types in [`product_providers.jl`](../../src/plant/product_providers.jl), and the exact RNG-owner layer in [`randomness.jl`](../../src/plant/randomness.jl). Cold model-specific construction, provider-style trait dispatch, canonical selected-set preparation, mutation-free preflight, and warmed serial execution use multiple dispatch; IDs, dimensions, rates, times, and device ordinals remain values, and no universal graph, abstract executor collection, stored closure, scheduler, queue, or transport was introduced. [Topology](../../test/testsets/plant_topology.jl), [preparation](../../test/testsets/plant_preparation.jl), [provider](../../test/testsets/plant_providers.jl), [RNG](../../test/testsets/plant_rng.jl), [illumination](../../test/testsets/plant_illumination.jl), and the [Gate 2 serial plant artifact](../../benchmarks/results/gate2/2026-07-21-serial-plant.toml) cover structured errors, direct-science/WFS composition, provider and calibration seams, NGS/LGS/science execution, owner-bound stochastic extension seams, inferred return types, and zero warmed allocation. Transitional numerical oracles remain unextended and are deleted only at their assigned replacement gates. |
| HIL-API-002 | 2 | Core plant and calibration APIs exchange structured data without owning path-based cache policy, serialization codecs, or hard-coded configuration file formats | Core + optional extension/user integration | implemented | validated | [`config_dict` and `snapshot_config`](../../src/core/config.jl) return caller-owned structured configuration data without a core writer. [`MetaSensitivity`, `compute_meta_sensitivity_matrix`, and `SPRINT`](../../src/calibration/misregistration_identification.jl) expose the complete typed calibration result and iterative state without cache-path, save, or recompute persistence fields or keywords; the AD implementation follows the same boundary. [`calibration_and_analysis.jl`](../../test/testsets/calibration_and_analysis.jl) executes AD and finite-difference workflows in an empty temporary directory, verifies no artifact appears, rejects the removed persistence keywords, and checks that SPRINT retains no persistence state. [`reference_and_tutorials.jl`](../../test/testsets/reference_and_tutorials.jl) validates pure structured configuration snapshots, removal of the TOML writer, and package metadata with `Serialization` and TOML absent from core dependencies. TOML remains test-only, while file formats, cache compatibility, selection, and campaign storage belong to a caller or optional extension. |
| HIL-ATM-003 | 2 | A serial plant materializes each due path's atmosphere-dependent OPD, field, or model-specific input before the atmosphere writer advances; stale epoch tokens are rejected and are never treated as retained layer storage | Core | implemented | validated | [`PreparedPupilOPDMaterialization`](../../src/plant/preparation.jl) binds a prepared direction renderer to the exact caller-owned path pupil, while qualified validation/materialization dispatch supports model-specific field or layer-aware operations and an explicit atmosphere-independent assertion. A prepared acquisition selection validates every selected binding and current epoch before materializing all unique path inputs and forming downstream results. [`plant_preparation.jl`](../../test/testsets/plant_preparation.jl) proves same-epoch NGS/LGS/science materialization, retained caller-owned inputs after the writer advances, stale and foreign token rejection without partial mutation, order invariance, one-time shared path formation, and zero warmed CPU allocation. Cross-timestamp retained state remains separately scoped to `HIL-ATM-004`. |
| HIL-RNG-001 | 2 | Preparation derives run-unique stable RNG owner identities and per-owner streams from a central seed and versioned scheme so endpoint/path order does not change stochastic results in the serial plant | Core | implemented | validated | [`randomness.jl`](../../src/plant/randomness.jl) derives independent `Xoshiro` streams from a required run seed, positive `RNGDerivationVersion`, and stable `(category, component, role)` identity using an explicitly recorded, process-stable byte encoding rather than Julia's `hash`. Prepared multilayer atmospheres require named `AtmosphereLayerID` values; paths and acquisitions receive baseline provider/detector domains and may declare additional device roles. Preparation rejects missing or duplicate identities and derived-seed collisions, while selected execution validates exact owner bindings and passes streams directly to stochastic methods. Replay metadata records both derivation and stream algorithms. [`plant_rng.jl`](../../test/testsets/plant_rng.jl) validates a golden derivation value, canonical replay metadata, exact same-seed replay, path/acquisition/selection/layer reorder invariance, mixed deterministic and stochastic owners, seed/version/owner separation, quiet-stream non-consumption, and zero warmed serial CPU allocation. Tuple position, task, thread, ring cursor, completion order, and physical-device ordinal are absent from derivation. Addressable event/element randomness for replicated or reordered multi-device work remains explicitly scoped to `HIL-RNG-002`. |
| HIL-FID-001 | 2 | Each acquisition binds one run-immutable prepared full-optical, command-responsive reduced-order, or synthetic/replay product provider while preserving one logical product shape, type, geometry/radiometry, metadata, and provider-result contract | Core | implemented | validated | [`product_providers.jl`](../../src/plant/product_providers.jl) defines provider-style and payload-work traits, required acquisition metadata, defensive logical-product contracts, the exact caller-owned result rule, a full-optical wrapper, and unchanged/copy/bounded-cyclic-replay implementations. [`preparation.jl`](../../src/plant/preparation.jl) binds exactly one provider into each immutable acquisition owner and forms only paths required by selected full-optical providers. [`plant_providers.jl`](../../test/testsets/plant_providers.jl) proves one consumer accepts compatible full, test-extension command-responsive reduced-order, and synthetic/replay products; validates shape, numeric type, backend/device, units/metadata, result identity, fixed selection, path bypass, cyclic replay, nonresponsiveness, and zero warmed allocation. Payload-work declarations explicitly are not load evidence. Meaningful scheduled closed-loop reduced-order behavior remains `HIL-FID-003` at Gate 4, while descriptors, leases, ports, overload conformance, and production-shaped load evidence remain Gate 4A/10A work. |
| HIL-CAL-001 | 2 | Calibration illumination composes through ordinary source/path or detector-input seams with user-declared typed entry, visibility, timing, state, and combination semantics; core assumes no instrument topology, source physics, propagation bypass, or control authority | Core + user model | implemented | validated | [`illumination.jl`](../../src/plant/illumination.jl) binds supported pupil-function, declared-plane field/intensity, external-optics-result, and detector-input products to an immutable prepared evaluator wrapper, exact caller-owned destination, defensive visibility/payload contract, explicit combination trait, epoch time, and stable path-owned `:illumination` RNG. [`preparation.jl`](../../src/plant/preparation.jl) materializes that product before ordinary path execution and detector/provider acquisition without a calibration role, instrument selector, or upstream-bypass flag. [`plant_illumination.jl`](../../test/testsets/plant_illumination.jl) covers the native uniform-intensity evaluator plus a user-defined stateful pupil evaluator, all entry tags, visibility ownership, declaration-order replay, time/RNG state, normal downstream optics/acquisition, strict rejection, and warmed CPU allocation. The [optional backend matrix](../../test/backend_optional_common.jl) exercises the native entry end to end through device-resident path and detector acquisition with scalar indexing disabled. Scheduling, trigger/setpoint ownership, HIL descriptors, ports, and transport are explicitly outside this Gate 2 requirement. |
| HIL-PATH-003 | 3 | Acquisitions over native NGS, finite-height LGS, and direct-science paths use independent virtual-time schedules | Core | partial | partial | [`source_geometry.jl`](../../src/atmosphere/source_geometry.jl) and [`runtime/arms.jl`](../../src/control/runtime/arms.jl) provide frozen NGS/LGS directions and one current-state epoch token for independent WFS/science paths; atmosphere and runtime tests cover distinct geometry. Common multi-rate virtual-time orchestration remains a gap. |
| HIL-DET-001 | 3 | Exposure, optical samples, rolling shutter, frame transfer, nondestructive reads, readout, presampling detector response, charge coupling, and complete-product publication have explicit composable event semantics; the realized response has a validated interior MTF where supported and explicit finite-frame boundary behavior | Core | partial | partial | [`detectors/interface.jl`](../../src/detectors/interface.jl) and [detector tests](../../test/testsets/detectors.jl); scheduled plant integration remains a gap |
| HIL-DET-002 | 3 | Scheduled up-the-ramp operation samples evolving accumulated charge at explicit nondestructive-read events without ending integration; the post-exposure synthesized ramp remains a separately declared lower-fidelity convenience and cannot establish time-resolved behavior | Core | planned | missing | Current [`frame_sampling.jl`](../../src/detectors/frame_sampling.jl) synthesizes fractional reads after final integration. Gate 3 adds event-driven read storage and exact atmosphere/command boundary tests. |
| HIL-TIME-001 | 3 | Nonnegative integer-nanosecond plant timestamps and durations use checked arithmetic, while event phase, prepared ordinal, and per-generator occurrence determine total equal-time ordering | Core | partial | partial | [`plant/time.jl`](../../src/plant/time.jl) defines distinct checked timestamp/duration values, positive periodic recurrence, and a total logical event key. [`plant/scheduling.jl`](../../src/plant/scheduling.jl) assigns prepared phase/ordinal registry order and rejects any claimed, activated, or rescheduled key that would regress the total order while allowing a strictly later occurrence at the same timestamp. [`plant_time.jl`](../../test/testsets/plant_time.jl) and [`plant_scheduler.jl`](../../test/testsets/plant_scheduler.jl) cover invalid construction, overflow/underflow, type separation, simultaneous ordering, equal-time occurrences, inference, and zero warmed allocation. Half-open detector execution and trigger-derived event timing remain Gate 3 work in #47–#49. |
| HIL-TRIG-001 | 3 | A prepared trigger-source and distribution topology maps nominal edges to per-acquisition delivery times with explicit phase/skew, jitter, bounded non-overtaking delivery, dropped/duplicate-edge policies, and separate physical exposure and timestamp-label semantics | Core | planned | missing | Camera-internal oscillator models are optional; the baseline covers externally visible trigger behavior and gives explicit duplicates an occurrence order without creating unbounded in-flight edges |
| HIL-SCHED-001 | 3 | A fixed-capacity scheduler jumps to due timestamps without base-tick polling, run-length-sized event materialization, or warmed allocation | Core | implemented | validated | [`plant/scheduling.jl`](../../src/plant/scheduling.jl) separates canonical immutable definitions, fixed-length `Memory` registry, compact single-writer cursors, caller-owned due slots, and isbits claims; supports checked activation, deactivation, and same-/later-time rescheduling; and structurally rejects duplicate ordinals, overflow, time/order regression, stale/foreign claims, and invalid transitions. [`plant_scheduler.jl`](../../test/testsets/plant_scheduler.jl) covers reversed declarations, simultaneous phases/ordinals, equal-time occurrences, empty/long-period execution, failures, inference, 100,000-event fixed storage, and zero warmed allocation. The [Gate 3 scheduler artifact](../../benchmarks/results/gate3/2026-07-21-event-scheduler.toml) records three 20,000-sample HdrHistogram runs at 1, 8, 32, 128, and 256 active generators, storage, GC, zero allocation, timer overhead, and environment. It is self-paced serial CPU service-cost evidence only; no HIL wall-clock latency, fixed-arrival capacity, detector, trigger, command, or heap crossover is claimed. |
| HIL-CMD-001 | 4 | Every independently timed controllable optic or segment has an independent virtual-time command endpoint | Core | planned | missing | Replaces composite command packing; the minimal HIL port arrives in Gate 4A |
| HIL-CMD-002 | 4 | Validation, admission, effective time, application, hold, and terminal model disposition are distinct | Core | planned | missing | Atomic multi-optic latch is explicit, never inferred from placement |
| HIL-CMD-003 | 4 | Future, late, rejected, and superseded commands use bounded explicit policies and never backdate plant state | Core | planned | missing | Every displaced or rejected command has a terminal model disposition |
| HIL-CMD-004 | 4 | Every endpoint prepares a versioned core plant command schema covering payload type/shape, units, basis/calibration, absolute or incremental semantics, bounds, plant-effective-time policy, and duplicate/reordering behavior | Core | planned | missing | This row was split during the pre-HIL architecture review: core owns semantic command interpretation and model disposition; HIL boundary descriptors, session correlation, timestamp mapping, leases, and outcome credit are tracked by `HIL-PORT-005`. |
| HIL-CMD-005 | 4 | Every endpoint declares replayable plant-time command-silence hold/safe/fail behavior, including the age origin and exact equal-time reset policy | Core | planned | missing | Safe state and modeled command-age storage are prepared; no transition allocates or invokes user or transport code. The separate operational ingress-liveness failure is tracked by `HIL-LIFE-002`. |
| HIL-FID-003 | 4 | A reduced-order plant evolves time-correlated disturbances and maps correctly timed effective optic commands through calibrated path and sensor operators so commands causally change later products and a matched reference controller can close the loop | Core | partial | missing | [`interaction_matrix.jl`](../../src/calibration/interaction_matrix.jl), [`reconstructors.jl`](../../src/control/reconstructors.jl), geometric WFS modes, and tomography primitives are foundations; the prepared reduced-order plant and acceptance evidence are missing |
| HIL-OPT-001 | 4 | Co-conjugated DMs remain independent devices even when prepared as one optical execution group | Core | partial | partial | [`controllable_optics.jl`](../../src/optics/controllable_optics.jl) and [runtime tests](../../test/testsets/control_and_runtime.jl) cover additive composites only |
| HIL-OPT-004 | 4 | Autonomous periodic optical devices use trigger-relative waveform state and bounded setpoint commands rather than one RTC message per waveform point | Core | partial | partial | [`PyramidWFS`](../../src/wfs/pyramid/setup.jl) and its [modulation optics](../../src/wfs/pyramid/optics.jl) provide the baseline cycle-averaged foundation; physical trigger relation remains missing and time-resolved fidelity is profile-driven |
| HIL-VSLICE-001 | 4A | A serial CPU vertical slice connects one scheduled acquisition and one command-responsive optic to a deterministic in-memory fake RTC through the canonical complete-product and command/outcome boundary with an injected clock and fixed-arrival evidence | HIL + Core | planned | missing | Proves the external contract before worker, GPU, multi-path, or transport-specific optimization |
| HIL-TIME-003 | 4A | The first HIL vertical slice injects a `Clocks.jl` deterministic test clock and monotonic production clock while core remains wall-clock independent | HIL | planned | missing | Gate 8 hardens cached-clock ownership, staleness evidence, lifecycle, and external-domain mappings; it does not introduce clock injection again. |
| HIL-PORT-001 | 4A | Canonical command-submission, command-completion/outcome, and acquisition-completion ports are transport-neutral | HIL | planned | missing | Successful ring enqueue transfers ownership but is not semantic command admission |
| HIL-PORT-005 | 4A | A prepared HIL command-submission descriptor schema maps endpoint/session correlation, external timing metadata, payload-lease ownership, and outcome credit into a compatible core plant command schema without making HIL types dependencies of core | HIL | planned | missing | The paired HIL command outcome wraps the core terminal model disposition with boundary timing and returns the submission credit. |
| HIL-PORT-003 | 4A | External integration declares adapter readiness, complete-product-to-first-observation lead time, and maximum lease hold time without promoting progressive fragments into canonical simulation events | HIL + user integration | planned | missing | First-packet timing remains an external-delivery claim |
| HIL-PORT-004 | 4A | Sampled optic/device feedback uses ordinary acquisition endpoints and remains distinct from correlated command outcomes | Core + HIL | planned | missing | Supports scalar, vector, surface, encoder, and health products without an instrument-specific API |
| HIL-BUF-001 | 4A | Large acquisition and command payloads use bounded pools and explicit generation-checked leases rather than ring-resident payload copies | HIL + user integration | planned | missing | Memory domain, generation, terminal-outcome/release ownership, and reuse are explicit; return-credit hardening follows in Gate 8 |
| HIL-OPT-002 | 5 | Common MCAO and path-specific MOAO planes support NGS and finite-height LGS footprints | Core | planned | missing | Requires placement and visibility traits |
| HIL-OPT-003 | 5 | Common and path-specific sampled aberrations, including NCPA, have explicit visibility | Core | partial | partial | [`ncpa.jl`](../../src/optics/ncpa.jl), [`opd_map.jl`](../../src/optics/opd_map.jl), and [model tests](../../test/testsets/calibration_and_analysis.jl); path visibility remains a gap |
| HIL-CPU-001 | 6 | Prepared CPU execution groups have single-writer workspaces, explicit thread budgets, independently callable executor seams, and a deterministic serial fallback | Core | planned | missing | Long-lived HIL agents and their rings arrive in Gate 8 |
| HIL-ATM-004 | 6 | Parallel execution either holds the atmosphere writer while same-epoch readers materialize or uses bounded model-specific retained atmosphere-state snapshots; it never uses an epoch token as retained storage and accounts for every materialization/snapshot slot | Core + HIL | planned | missing | Downstream path work may overlap after consuming only caller-owned materialized products. Snapshot representation is atmosphere-model specific rather than an unconditional full-layer copy. |
| HIL-EXEC-006 | 6 | Prepared execution keeps small fixed stage pipelines specialized while representing large path/endpoint registries with bounded compilation and code growth verified against topology scale | Core | planned | missing | Avoids encoding a MORFEO-scale topology solely as one recursively specialized tuple or `NamedTuple`; preparation latency and generated-code size are evidence surfaces. |
| HIL-GPU-001 | 7 | Physical device identity, prepared direction batching, and device-resident state are explicit | Core | partial | partial | [CUDA](../../ext/AdaptiveOpticsSimCUDAExt.jl), [AMDGPU](../../ext/AdaptiveOpticsSimAMDGPUExt.jl), and [backend tests](../../test/backend_optional_common.jl) now validate device-resident epoch rendering and prepared atmospheric fields with physical-device compatibility checks; prepared batching and physical multi-device planning remain gaps. |
| HIL-TIME-002 | 8 | External timestamp domains map into plant time with versioned offset, drift, uncertainty, and no retroactive remapping | HIL | planned | missing | User integration supplies synchronization observations |
| HIL-LIFE-001 | 8 | Configure, prepare, arm, run, and stop/fail phases preserve ownership, adapter-readiness preconditions, prepared nonstructural state transitions, and bounded shutdown | HIL | planned | missing | Topology, schema, capacity, placement, and provider changes require another prepare/arm cycle |
| HIL-LIFE-002 | 8 | An optional execution-clock RTC-ingress-liveness watchdog resets only after semantic command admission, fails the run without silently changing optic state, and remains distinct from replayable plant-time command-silence behavior | HIL | planned | missing | Threshold, reset/recovery rules, clock identity, and exact-boundary behavior are recorded and tested; malformed traffic cannot keep a run alive. |
| HIL-PORT-002 | 8 | Every port, schedule, pool, lease, and tap has resource-specific full, close, drain, and recovery semantics | HIL | planned | missing | Silent backlog, unread-slot overwrite, and lease loss are forbidden |
| HIL-RING-001 | 8 | Each hot ownership handoff uses a bounded SPSC descriptor ring with release/acquire publication and isolated cursors | HIL | planned | missing | Requires layout, generated-code, concurrency, saturation, and recovery evidence |
| HIL-BUF-002 | 8 | Each return path reserves usable credit for every lease its consumer can hold, so a valid first release cannot encounter ordinary full backpressure | HIL | planned | missing | Pool accounting covers free, producer-owned, completion-queued, consumer-leased, and return-queued states |
| HIL-FAIL-001 | 8 | Required and optional work have declared shed, stop, failure, and recovery behavior without runtime provider substitution | HIL | planned | missing | Already admitted commands cannot disappear silently |
| HIL-FAIL-002 | 8 | Every execution owner has a preallocated first-failure publication and acknowledgement path; one coordinator closes ingress, publishes stop, drains outcomes/leases, and reports ownership deficits within bounded deadlines | HIL | planned | missing | Error formatting stays outside the stopped critical path |
| HIL-EXEC-004 | 8 | Same-process runtime, adapter, and telemetry components have warmed allocation budgets and GC evidence, or allocating integration is isolated behind an explicit process boundary | HIL + user integration | planned | missing | The manifest records GC, process sharing, page policy, and soak evidence; GC is not disabled by assumption |
| HIL-EXEC-001 | 9A | One plant statically places complete execution groups across CPU workers and one or more GPUs | Core + HIL | planned | missing | Large cross-resource payloads require explicit bounded handoffs |
| HIL-EXEC-002 | 9A | Fully explicit and constrained deterministic placement resolve to one inspectable immutable plan | HIL | planned | missing | Conservative rules may fill unassigned groups; they do not override hard user placement |
| HIL-EXEC-003 | 9A | Preparation validates constraints, capability, memory, transfers, burst utilization, and reserved contexts and records its rationale | HIL | planned | missing | Infeasible plans fail structurally; runtime migration is excluded |
| HIL-GPU-002 | 9B | Static multi-GPU placement preserves epoch, command, detector, RNG, and sequence consistency | Core + HIL | planned | missing | Homogeneous multi-GPU validation follows single-resource and mixed CPU/GPU evidence |
| HIL-RNG-002 | 9B | Replicated or reordered multi-device stochastic work uses addressable random domains derived from run seed, derivation version, stable owner, event/epoch, and element/sample identities rather than sequential host seed consumption in launch order | Core + HIL | planned | missing | Required for replicated atmosphere evolution and placement-independent detector/provider randomness; numerical output comparison remains tolerance-based across backends. |
| HIL-EXEC-005 | Future | A calibrated fully automatic cost-model planner may assign complete groups while preserving the same immutable plan and admission contract | HIL | planned | missing | Requires representative real-profile measurements before implementation |
| HIL-CPU-002 | 10B | Optional affinity records verified CPU, SMT, NUMA, FFT, BLAS, device-agent, and interrupt placement | HIL | planned | missing | `ThreadPinning.jl` is deployment policy |
| HIL-PATH-004 | 10B | A prepared external optical executor, including `Proper.jl`, runs at its acquisition cadence without gating an independently placed required WFS path | HIL + optical companion | partial | partial | [Integration example](../../examples/integrations/proper_hil_coronagraph.jl), [benchmark](../../scripts/profile_proper_hil_coronagraph.jl), and sibling proving ground exist; canonical asynchronous HIL integration remains a gap |
| HIL-OBS-001 | 10A | Fixed-arrival latency, selected external-delivery timing, outcomes, occupancy, command age, trigger phase/skew, GC, scheduling, and device synchronization are observable | HIL + user integration | partial | partial | [Histogram harness](../../benchmarks/benchmark_detector_hil_latency.jl) and [artifact](../../benchmarks/results/detectors/2026-07-14-detector-hil-latency.toml) are closed-loop/in-process only; canonical fixed-arrival ports and adapter-owned transport timing remain gaps |
| HIL-FID-002 | 10A | Low-cost RTC load evidence declares payload work and memory domain, separates schedule-preserving fixed-arrival latency from unpaced saturation, and never implies unsupported optical or closed-loop fidelity | Core + HIL | partial | partial | The [detector HIL benchmark](../../benchmarks/benchmark_detector_hil_latency.jl) uses deterministic prepared inputs and histograms; a canonical production-shaped source, open-loop arrival accounting, and tiered claim evidence remain gaps |
| HIL-REPLAY-001 | 10A | Scenario, boundary-traffic, and decision/event replay are distinct claims with canonical in-memory records for arrivals, mappings, decisions, outcomes, faults, gaps, placement, configuration, and command history | HIL | partial | missing | [`reference_harness.jl`](../../test/reference_harness.jl) provides deterministic scenario replay; canonical boundary and decision/event replay evidence are missing |
| HIL-DIST-001 | Future | Multi-process or multi-host HIL has explicit clock, transport, ownership, partition, and recovery contracts | HIL + user integration | planned | missing | Not implied by single-host multi-GPU support |
| HIL-DIST-002 | Future | An optional single-host iceoryx2 process boundary preserves canonical schemas and explicitly bounds loans, subscriber buffers, fan-out lifetime, overflow/backpressure, discovery, and dead-node recovery | HIL extension + user integration | planned | missing | Iceoryx2 is not an in-process replacement for direct calls or SPSC rings; direct-to-loan versus copied payloads are measured |
| HIL-OPT-005 | 10C | Time-resolved pyramid modulation is added only for a profile whose partial-cycle, dither, mirror-response, or synchronization claim cannot use the cycle-averaged baseline | Core + profile package | planned | missing | Does not block the first HIL boundary or cycle-averaged pyramid support |
| HIL-PROFILE-NFIRAOS-001 | 10C | A pinned NFIRAOS profile independently validates synchronized multi-rate topology, pyramid modulation, detector-trigger faults, model fidelity, and external integration | Companion validation package | planned | missing | Instrument-specific configuration remains outside core and the general HIL package |
| HIL-PROFILE-MORFEO-001 | 10C | A pinned MORFEO-scale profile validates its declared camera or subsystem boundary, endpoint topology, complete products, adapter transport contract, model, timing, commands, and external integration independently | Companion validation package | planned | missing | Camera-boundary and internal-vector results are separate claims; instrument-specific configuration remains outside core |

## Gate Summary

| Gate | Capability outcome | State |
|---:|---|---|
| 0 | Optical planes/radiometry, explicit atmosphere time, direct-science acquisition, and WFS stages decomposed and validated | complete |
| 1 | Contracts and correctness oracles frozen | complete |
| 2 | Shared plant, acquisition workspace, product providers, and source-entry seams separated | complete |
| 3 | Deterministic multi-rate virtual-time event engine | planned |
| 4 | Independent controllable-optic and meaningful reduced-order command semantics | planned |
| 4A | Minimal serial CPU HIL vertical slice through canonical ports | planned |
| 5 | Conjugated and path-specific controllable planes | planned |
| 6 | Prepared CPU execution groups and executor seams | planned |
| 7 | GPU direction batching and physical device identity | planned |
| 8 | Hardened lifecycle, ports, failure propagation, buffers, and deployment behavior | planned |
| 9A | Explicit or constrained deterministic mixed CPU/GPU placement | planned |
| 9B | Static homogeneous multi-GPU placement | planned |
| 10A | Boundary replay, observability, and low-cost load evidence | partial |
| 10B | Deployment, affinity, and external-optics integration | partial |
| 10C | Versioned NFIRAOS and MORFEO profile evidence | planned |

The July 18 Gate 0 backend-evidence closure applies to `HIL-PLANE-001` through
`HIL-PLANE-004`, `HIL-ATM-001` and `HIL-ATM-002`, `HIL-SRC-001`,
`HIL-RAD-001` and `HIL-RAD-002`, `HIL-SCI-001`, and `HIL-WFS-009`. The full
CPU suite passed, the maintained CUDA target passed `386/386`, and the
maintained AMDGPU target passed `396/396`, with accelerator scalar indexing
disabled. The clean [local CPU](../../benchmarks/results/platform/2026-07-18-pre-hil-11-local-cpu.toml),
[WSL CPU](../../benchmarks/results/platform/2026-07-18-pre-hil-11-wsl-cpu.toml),
and [WSL CUDA](../../benchmarks/results/platform/2026-07-18-pre-hil-11-wsl-cuda.toml)
artifacts add correctness, residency, allocation, and serial service-time
evidence for the REVOLT-like boundary. The final Gate 0 ownership refactor
closes the former telescope-path-state gap in `HIL-PLANE-001`; these artifacts
still do not claim external-RTC response time or replace the earlier maintained
AMD latency characterization.

The pre-HIL architecture review closes Gate 1 by reconciling package and type
ownership, atmosphere-token/materialization lifetime, detector event semantics,
stable RNG ownership, execution-clock sequencing, and the breaking target API.
`HIL-ORACLE-001` applies to capabilities implemented at the current gate;
later multi-rate, placed-optic, parallel, and multi-device requirements carry
their own oracle and evidence obligations rather than making Gate 1 depend on
future gates.

## Capability Gates

The target advances through durable capability gates. A PR may satisfy part or
all of a gate, but references the applicable requirement IDs and updates their
implementation and evidence state. Every gate preserves numerical correctness
and the serial oracle. Public API compatibility is not required for the
independent-optic refactor; superseded surfaces are removed rather than kept as
permanent adapters.

No implementation of the proposed general `AdaptiveOpticsHIL.jl` companion
runtime begins until Gates 0 and 1 are complete. Existing integration proving
grounds remain valid. The ordered work is tracked by
[Gate 0 issue #1](https://github.com/DarrylGamroth/AdaptiveOpticsSim.jl/issues/1),
and final review and closure are tracked by
[issue #9](https://github.com/DarrylGamroth/AdaptiveOpticsSim.jl/issues/9).

### Gate 0: Decompose optical dataflow and acquisition

- freeze deterministic pre-refactor telescope, atmosphere-direction,
  direct-science, Shack-Hartmann, Pyramid, BioEdge, Zernike, Curvature, and LiFT
  products, including radiometric factors and detector-coupled products
- separate the telescope aperture/spatial-geometry definition from temporal
  cadence, caller-owned path OPD/electric-field/intensity products, and
  propagation workspaces
- apply static and controllable optical surfaces to explicit path products
  rather than shared telescope scratch
- split reusable spatial-filter definition, workspace, and output ownership
- advance one shared atmosphere to explicit model time separately from path-
  local prepared NGS, LGS, spectral, and extended-source rendering into caller-
  owned destinations
- give every optical product prepared geometry, wavelength/channel,
  radiometric, coherence/combination, backend, and device metadata
- establish photon-arrival-rate products, detector-owned elapsed-time integration,
  observation, and measurement stage contracts that support one or multiple
  planes and detectors
- retain incompatible spectral grids as typed bundles or require an explicit
  prepared mapping instead of silently accumulating arrays by index
- form native or prepared external direct-science photon-arrival-rate products
  independently of one or more detector acquisitions
- extract the microlens array and migrate Shack-Hartmann as the reference
  composition, preserving geometric and diffractive behavior
- migrate Pyramid and BioEdge while separating their physical masks from
  reusable modulation and focal-plane propagation machinery
- migrate Zernike and Curvature, including separate-branch and packed Curvature
  acquisition forms
- decouple LiFT observation acquisition, forward modeling, and estimation
- retain direct geometric and reduced-order paths without unused diffractive
  workspace or fictitious intermediate products
- validate each migrated family incrementally on CPU and maintained backend
  smoke paths, then run the complete CUDA and AMDGPU hardware matrix

Acceptance: one telescope definition supports independent caller-owned path
products and owns no cadence/exposure duration; one explicitly timed frozen
atmosphere epoch supports order-independent path-local direction rendering;
optical surfaces and reusable spatial filters do not require shared telescope
path mutation; and compatible direct-science photon-arrival-rate products support
independent detector exposures. Plane mismatches fail before repeated
execution, incompatible spectral grids are not silently summed, and a non-unit
duration test proves elapsed time is applied exactly once. Every maintained WFS
family exposes the staged semantic boundary without introducing a universal
optical graph, an `isa` ladder, a general resampler, or HIL runtime dependencies.
Presampling response and physical-pixel integration occur after optical
photon-arrival-rate formation; QE and explicit detector elapsed time are then
applied exactly once before coupling, stochastic response, readout, and
estimation. Frozen complete products and new stage products pass their declared
tolerances; CPU steady-state allocation
budgets do not regress; device-resident paths do not introduce unintended host
transfers; and comparable warmed latency satisfies predeclared absolute and
relative gates over repeated runs. Every Gate 0 `HIL-PLANE-*`, `HIL-ATM-*`,
`HIL-SRC-*`, `HIL-RAD-*`, `HIL-SCI-*`, and `HIL-WFS-*` row is implemented and
validated, and model-validity and production-support surfaces are updated
before Gate 1 implementation work proceeds.

### Gate 1: Freeze contracts and correctness oracles

- define physical/event outputs for multiple directions, detector modes, direct
  and coronagraph science paths, and co-conjugated DM commands
- freeze temporary numerical references for current composite-optic cases
- define canonical timestamps, equal-time ordering, sequence domains, output
  ownership, and performance-boundary terminology
- distinguish current-state atmosphere epoch tokens, materialized path products,
  and optional retained atmosphere-state snapshots
- define stable RNG owner identities, derivation versioning, and the
  stateful-versus-addressable policy without requiring its later-gate executor
- distinguish frame-step incremental detector convenience from scheduled event
  operations and time-resolved nondestructive reads
- split core plant commands/model dispositions from HIL boundary descriptors,
  leases, outcome credit, and command outcomes
- freeze the breaking definition/prepared-plan/mutable-state/product API policy
  without preserving the original OOPAO or frame-step object layout
- retain the direct single-threaded CPU runtime as the numerical oracle

Acceptance: every changed physical surface has a deterministic reference or an
explicitly tracked later-gate requirement; current correctness, allocation, and
backend targets remain green. No acceptance condition depends on behavior first
introduced by a later gate.

### Gate 2: Separate shared plant models from path and acquisition state

Status: complete on the maintained serial CPU and backend-functional surfaces.

- introduce the breaking core definition/prepared-plan/mutable-state/product
  API without extending `AOSimulation` or `ClosedLoopRuntime`
- remove path-based cache and hard-coded configuration serialization from the
  target core surface; extensions or callers persist returned structured data
- separate immutable path definitions from independently scheduled or
  triggered acquisition endpoints
- render/materialize every due atmosphere-dependent path input into caller-
  owned destinations before the mutable atmosphere writer advances
- precompute source geometry and compatible-result keys per prepared path
- separate telescope parameters and pupil from path-local propagation workspace
  and acquisition-local WFS, detector, readout, and caller-owned product state;
  HIL publication state remains a later boundary concern
- bind full-optical, reduced-order, and synthetic/replay providers through one
  prepared mutating acquisition-product seam using dispatch and traits
- expose a narrow prepared calibration-illumination seam at supported typed
  path entries or detector inputs, leaving physical integration to user models
- prepare stable per-owner RNG streams from a central run seed, versioned
  derivation scheme, and declared component identities
- keep the first executor serial to isolate ownership correctness

Acceptance: frozen outputs remain within declared tolerances, every due path
sees one current atmosphere epoch token and has materialized its required input
before the next advance, path reuse does not couple acquisition state, endpoint
or path reorder does not change per-owner stochastic results, and every
applicable provider preserves the logical acquisition product's shape, type,
geometry/radiometry, metadata, and provider result. Reduced-order providers declare
their validity envelope; static/replay providers are explicitly nonresponsive;
and the selected provider remains immutable until another prepare/arm cycle.
A native and a user-provided calibration source preserve declared visibility,
timing, deterministic state, and composition without a calibration-role branch
or implicit bypass. Unsupported entry payloads and source combinations fail
during preparation. The warmed serial oracle meets its declared allocation
budget. Comparable latency is archived as a baseline; no subjective latency
gate is used. HIL port descriptors, leases, queue capacity, and overload
conformance are not Gate 2 acceptance criteria.

The focused Gate 2 testsets and maintained CUDA/AMDGPU functional matrices pass
on Julia 1.12.6. The clean [serial plant CPU
artifact](../../benchmarks/results/gate2/2026-07-21-serial-plant.toml) records
same-seed declaration-order replay, three distinct science/NGS/LGS path
materializations, shared-path detector fan-out, zero warmed allocation, and
three raw 10,000-sample service-time histograms. All Gate 2 rows are implemented
and validated. The artifact is deliberately not external-RTC latency,
fixed-arrival capacity, event-scheduler, or production-scale evidence.

### Gate 3: Add deterministic multi-rate virtual-time execution

- use integer plant timestamps, stable scheduler ordinals, and half-open
  exposure intervals
- prepare a fixed-capacity next-event schedule that jumps directly to due
  timestamps
- prepare finite trigger-source and distribution-link topology with correlated
  source faults, independent per-link faults, and one next delivered edge per
  active schedule
- represent rolling bands and nondestructive reads as bounded event-generator
  cursors; use a prepared linear scan first and add a stable array heap only
  after measured generator-count evidence justifies it
- schedule optical samples, rolling-shutter bands, nondestructive reads,
  readout completion, and frame publication independently per acquisition
- replace floating-tolerance completion authority with explicit scheduler-owned
  begin, accumulate, read, close, readout, and readiness events; retain the
  current auto-finalizing incremental call only as a frame-step convenience
- preserve presampling detector response, charge-coupling, and readout-pipeline ordering when
  acquisition is split into scheduled events
- accept explicit simulation timestamps without reading wall clock in core
- expose due/readiness products without transport dependencies

Acceptance: exact-boundary tests cover atmosphere evolution, trigger delivery,
exposure/row-band boundaries, optical samples, nondestructive reads, readout,
and publication sharing a timestamp. Trigger tests distinguish nominal source edges, delivered edges,
physical exposure boundaries, reported labels, and execution time; cover fixed
skew, correlated and independent jitter, phase steps, and dropped/duplicate
edges; and prove that faults affect only their declared downstream branches.
Representative CMOS rolling/global shutter, CCD, frame-transfer EMCCD, and
HgCdTe nondestructive-read cases preserve their frozen detector products,
including configured response/coupling behavior. HgCdTe scheduled ramp reads
sample evolving accumulated charge; the post-exposure synthesized convenience
is validated and labeled separately. Mixed rates produce the expected
sequences, long runs use storage proportional to active generators rather than
event count, and long-period schedules do not iterate empty base ticks.
Scheduler cost is archived across the maintained generator-count range, and
warmed execution meets its allocation budget.

### Gate 4: Replace composite optics with independent command semantics

- replace the single `AOSimulation.optic` field with a named registry of
  individual controllable optics
- introduce one independently timed endpoint per optic or latched segment
- prepare a versioned core plant command schema per endpoint, including units,
  basis/calibration revision, payload semantics, bounds, effective-time policy,
  and duplicate/reordering policy; boundary session/timestamp/lease metadata
  follows in Gate 4A
- separate validation, admission, effective time, application, hold, and
  terminal model disposition in the virtual-time core; add enqueue and HIL
  command outcomes in Gate 4A
- define bounded future-command and explicit late/supersession policies
- define command-silence hold/safe/fail behavior in plant time, distinct from
  any operational execution-time ingress-liveness failure
- support explicit atomic multi-optic transactions without inferring atomicity
  from plane placement or packed vectors
- remove `CompositeControllableOptic` and replace `RuntimeCommandLayout` with
  prepared core endpoint routing; canonical HIL descriptors follow in Gate 4A,
  while user transport schemas remain outside both packages
- split detector integration correctly when a command becomes effective
- add a prepared reduced-order plant whose time-correlated disturbances,
  path/sensor operators, and effective commands produce causally correct slopes
  or approximate raw pixels for an external RTC
- represent autonomous periodic optics with trigger-relative waveforms and
  bounded setpoint commands, retaining the fast cycle-averaged pyramid policy;
  defer time-resolved modulation until a profile requires it

Acceptance: deterministic timelines cover early, equal-time, future, late,
rejected, superseded, and atomic plant commands; every command presented to the
virtual-time core receives one terminal model disposition, and every admitted
command is applied once or ends in a declared failure; schema, shape,
calibration-revision, sequence, and plant-time command-silence boundaries have
exact tests; independent optics update without an internal RTC; the
cycle-averaged pyramid model preserves its frozen reference and trigger
relationship; a matched reference controller closes the
reduced-order loop and reduces its declared residual while wrong-sign, delayed, stale, and
mismatched-calibration cases degrade as expected; commands and acquisition
events sharing a timestamp follow the Gate 3 stable phase order; and frozen
physical outputs remain within tolerance.

### Gate 4A: Prove a minimal serial HIL vertical slice

- create the transport-free `AdaptiveOpticsHIL.jl` package boundary with only
  the canonical types needed by one serial path
- depend on `Clocks.jl` and inject a deterministic test clock plus a monotonic
  production clock without adding wall-clock reads to core
- connect one scheduled acquisition to a complete-product completion port and
  one command-responsive optic to a submission/outcome pair
- map one HIL submission descriptor schema, including session correlation,
  external timing metadata, payload lease, and outcome credit, into the core
  plant command schema without making core depend on HIL types
- use the intended bounded SPSC descriptor and generation-checked lease model,
  but defer cache-line/generated-code and multi-owner hardening to Gate 8
- run the deterministic serial core directly; do not introduce worker queues,
  GPU submission, placement planning, or a transport dependency
- provide an in-memory fake RTC that consumes and releases every product,
  submits canonical commands, consumes every terminal outcome, and closes a
  reduced-order loop
- require user orchestration to report adapter readiness and declare complete-
  product-to-first-observation lead time and maximum lease hold time
- expose sampled device feedback through the same acquisition contract when
  the selected vertical-slice optic provides it
- add one schedule-preserving fixed-arrival test and one unpaced saturation
  diagnostic over the same boundary

Acceptance: the fake RTC closes a meaningful reduced-order loop; fixed-arrival
accounting continues through stalls; every submitted command and product lease
is accounted for at clean stop; complete-product publication, adapter lead
time, command admission, command application, and closed-loop response remain
separately observable; and no API introduced by the slice assumes a particular
transport, worker topology, GPU, instrument, or progressive-readout protocol.

### Gate 5: Add conjugated and path-specific controllable planes

- introduce pupil and atmospheric-conjugate placement plus common/path-specific
  visibility traits
- implement geometric NGS/LGS footprint mapping first
- derive co-conjugated execution groups while preserving independent commands
- support common multi-altitude MCAO and isolated per-target MOAO planes
- attach sampled NCPA and static aberrations to explicit paths
- retain the optional `Proper.jl` seam for relay-dependent propagation

Acceptance: analytic footprint tests and serial multi-conjugate scenarios pass
for on-axis, off-axis, and finite-height sources; co-conjugated results remain
equivalent; NCPA and MOAO commands affect only selected paths.

### Gate 6: Prepare CPU execution groups

- derive independently executable due-path groups with single-writer workspaces
- co-locate compatible consumers for field or photon-arrival-rate-product reuse
- expose allocation-free `!` executor seams that accept explicit epoch and
  command snapshots without creating tasks or queues
- hold the atmosphere writer until same-epoch materializers finish or use a
  bounded model-specific retained state; downstream overlap consumes only
  caller-owned materialized path products
- preserve a deterministic single-threaded fallback
- declare Julia, FFT, and BLAS thread ownership and avoid nested parallelism
- expose stable worker placement and NUMA requirements
- keep small stage pipelines specialized while bounding preparation,
  compilation, and generated-code growth as path/endpoint registry size grows

Acceptance: each group is independently callable under a validation harness;
serial/grouped parity, atmosphere lifetime/accounting, and allocation gates
pass; declared thread budgets do not oversubscribe Julia, FFT, or BLAS
execution; and compile/code-size evidence remains inside its declared topology
envelope. Any performance promotion
satisfies a predeclared absolute contract and relative baseline gate with
repeated-run dispersion on controlled hardware.

### Gate 7: Add GPU direction batching and device identity

- distinguish physical device identity from backend family
- batch compatible directions and wavelengths on one prepared device owner
- own streams, FFT plans, allocator state, and observation barriers per device
- preserve device-resident detector and optic state to explicit boundaries

Acceptance: maintained CUDA and AMDGPU targets cover parity, residency,
synchronization, allocation, and model combinations. Timing is promoted only
under the same fixed-arrival absolute and relative evidence contract.

### Gate 8: Harden the HIL boundary and operational lifecycle

- retain `AdaptiveOpticsHIL.jl` without transport dependencies and generalize
  the Gate 4A boundary without changing its canonical schemas
- harden the Gate 4A `Clocks.jl` seam with cached-clock ownership, measured
  staleness, lifecycle behavior, and versioned external-domain mappings while
  keeping execution-clock pacing separate from modeled trigger delivery
- complete configure, prepare, arm, run, and bounded stop/fail lifecycle,
  including adapter-readiness preconditions and prepared nonstructural
  acquisition, trigger, shutter/calibration, and safe/hold transitions
- add the optional execution-clock RTC-ingress-liveness watchdog without
  changing optic state or conflating it with plant-time command silence
- instantiate long-lived execution owners from prepared CPU groups and a single
  submission owner for each selected GPU
- harden the canonical command-submission, command-completion/outcome, and
  complete-product acquisition-completion ports for every prepared owner
- assign one prepared command authority per endpoint or explicit atomic latch
  group; keep competing-producer arbitration outside the canonical data plane
- implement padded bounded SPSC descriptor rings and bounded product pools
- use one due-work/completion path per owner and fixed pools of materialized
  atmospheric path products plus model-specific retained state where required,
  rather than a contended global queue or token-as-snapshot assumption
- define resource-specific full, close, drain, lease, and recovery semantics
- reserve one usable return credit for every lease a consumer can hold and
  continuously check the complete pool-accounting invariant
- give every execution owner a preallocated first-failure/acknowledgement path;
  let one coordinator close ingress, publish stop, and drain ownership within
  declared deadlines
- establish same-process allocation and GC budgets, with an explicit process-
  isolation rule for adapters or telemetry that cannot meet them
- provide deterministic in-memory adapters and port conformance tests

Acceptance: exact execution-clock, ingress-liveness, and external-domain
mapping tests pass; malformed traffic cannot reset the ingress watchdog;
every transferred command produces one correlated terminal outcome; port
enqueue is distinct from semantic command admission; every full policy
preserves ownership; and ring layout, generated code, memory ordering,
concurrency, and allocation gates pass. A stalled optional science acquisition
consumes only its prepared slots and is then shed or fails under its declared
policy without delaying an independently placed required WFS path. Clean stop
accounts for every lease and outcome; a valid first lease release never sees
ordinary full backpressure; injected worker/GPU failures identify unacknowledged
owners and ownership deficits; and same-process GC or isolation behavior has
soak evidence.

### Gate 9A: Add static mixed CPU/GPU placement

- inventory CPU worker sets, NUMA nodes, physical GPUs, memory domains, and
  reserved coordination contexts
- accept fully explicit and constrained deterministic policies over complete
  prepared execution groups
- record immutable plan, planner/estimate versions, constraints, headroom,
  handoffs, stable tie-breaks, and rationale
- reject unsupported, conflicting, unassigned, memory-infeasible, and
  burst-overloaded plans with structured diagnostics
- replicate immutable plant data and deterministic atmosphere/RNG state
- broadcast command snapshots without per-frame full-screen copies
- prohibit opportunistic runtime migration
- keep multi-GPU, multi-process, fully automatic planning, and multi-host HIL
  outside this gate

Acceptance: planning is reproducible, explicit constraints are honored, and
infeasible plans are rejected. Mixed CPU/GPU tests preserve epoch, command,
detector, RNG, and sequence consistency. Each promoted placement meets its
fixed-arrival contract; performance advantages are measured against the best
applicable single-resource placement.

### Gate 9B: Add static homogeneous multi-GPU placement

- inventory physical device identities, peer-transfer capability, memory
  domains, submission owners, and reserved host coordination contexts
- partition complete path groups across homogeneous devices without splitting
  one ordinary path stage-by-stage
- replicate immutable plant state and derive atmosphere/RNG state from stable
  identities while broadcasting compact command snapshots
- keep one pool, stream/FFT-plan set, workspace set, and failure owner per
  device; make every peer or host transfer explicit and bounded
- retain fully explicit or constrained deterministic placement and prohibit
  opportunistic migration
- keep mixed-vendor and multi-process/multi-host claims outside this gate

Acceptance: homogeneous two-GPU tests preserve epoch, command, detector, RNG,
lease, and sequence consistency through burst, saturation, injected failure,
and recovery/stop. Each promoted placement meets its fixed-arrival contract and
beats or provides required capacity beyond the best applicable single-device
placement. A fully automatic cost optimizer remains a future independent gate.

### Gate 10A: Add boundary replay, observability, and load evidence

- add run manifests and distinct scenario, canonical boundary-traffic, and
  decision/event replay sources and sinks
- record trigger/timestamp-label faults, admission and overload decisions,
  outcomes, sequence gaps, command history, and bounded telemetry taps
- add preconfigured `HdrHistogram.jl` interval evidence, load curves,
  coordinated-omission audit, soak tests, and recovery tests
- add production-shaped synthetic traffic and command-responsive reduced-order
  profiles with explicit payload-work, memory-domain, and claim declarations
- run schedule-preserving fixed-arrival and separate unpaced saturation tests
  for every promoted low-cost load profile
- retain persistence codecs, artifact stores, and replay transport outside the
  general HIL package

Acceptance: boundary-traffic replay reproduces canonical port behavior without
claiming recomputed physics; decision/event replay reproduces outcomes and
sequence gaps; fixed-arrival target, burst, overload, recovery, and soak
evidence satisfies the declared absolute and relative gates. Every result names
its replay class and payload semantics, limits claims to its selected immutable
provider, and does not report descriptor-only evidence as raw-pixel throughput.

### Gate 10B: Add deployment and external-optics integration

- adapt the `AdaptiveOpticsProperHIL.jl` TCP/PROPER proving ground to canonical
  ports while keeping its codec application-local
- add optional verified `ThreadPinning.jl` deployment policy and record CPU,
  SMT, NUMA, FFT/BLAS, GPU-agent, transport, interrupt, memory, page, and GC
  placement/configuration
- validate same-process adapter allocation budgets and a separate-process
  adapter configuration without making either transport canonical
- retain transports, RTC codecs, bench setup, and camera profiles outside the
  general HIL package

Acceptance: a user adapter closes an external RTC loop against independently
scheduled or triggered acquisitions and controllable optics. A prepared PROPER
science acquisition meets its own cadence, and its declared stall/failure
policy leaves an independently placed required WFS path unaffected. Affinity,
GC, and process-isolation claims have fixed-arrival and soak evidence on their
declared deployment.

### Gate 10C: Add versioned instrument-profile evidence

- build NFIRAOS and MORFEO scenarios in a companion validation/configuration
  package rather than hard-coding instrument modes in core or HIL
- give each profile production-shaped synthetic, meaningful command-responsive
  reduced-order, and applicable full-optical variants behind the same canonical
  endpoints
- add time-resolved pyramid modulation only when the pinned NFIRAOS claim needs
  fidelity beyond the cycle-averaged baseline
- declare camera versus subsystem boundary, endpoint inventory, feedback
  acquisitions, adapter lead/hold time, transport, placement, and every claim
  limit in the pinned manifest

Acceptance: topology, model, timing, and external-integration compliance remain
independent. NFIRAOS synchronization and MORFEO scale evidence each pass their
own fixed-arrival, fault, recovery, and applicable optical/reference tests;
reduced-scale, single-device, and subsystem-boundary results are labeled and
are not promoted as full-profile camera-boundary support.

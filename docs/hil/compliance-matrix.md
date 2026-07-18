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
| `HIL-PLANE`, `HIL-ATM`, `HIL-SRC`, `HIL-RAD`, `HIL-SCI` | [`plant-and-optics.md`](plant-and-optics.md) and [`validation.md`](validation.md) |
| `HIL-WFS` | [`plant-and-optics.md`](plant-and-optics.md) and [`validation.md`](validation.md) |
| `HIL-ORACLE`, `HIL-OBS`, `HIL-REPLAY` | [`time-and-scheduling.md`](time-and-scheduling.md) and [`validation.md`](validation.md) |
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
| HIL-PLANE-001 | 0 | A prepared telescope owns run-immutable aperture, reflectivity, spatial sampling, and geometry while each path owns its mutable pupil OPD/electric field and focal-plane result; telescope owns no temporal cadence or exposure duration | Core | partial | partial | [`TelescopeAperture`](../../src/optics/telescope.jl) owns prepared geometry, independent [`PupilFunction`](../../src/optics/planes.jl) products own path OPD/amplitude, direct-science results are caller-owned, and `TelescopeParams` contains no sampling time or exposure duration. [`Explicit optical products and surfaces`](../../test/testsets/core_optics.jl) verifies path isolation. Transitional `LegacyTelescopePathState` still owns OPD for WFS families that have not migrated, so this requirement remains partial. |
| HIL-PLANE-002 | 0 | Caller-owned optical-plane products are separate from single-writer prepared propagation/spatial-filter workspaces, FFT plans, and scratch, with fixed prepared shape, numeric type, backend, and device | Core | implemented | validated | Caller-owned [`ElectricField` and `IntensityMap`](../../src/optics/planes.jl) products are separate from [`FraunhoferPropagation`/`FresnelPropagation`](../../src/optics/propagation.jl), [`DirectImagingWorkspace`](../../src/optics/direct_imaging.jl), and [`SpatialFilterWorkspace`](../../src/optics/spatial_filter.jl). Existing explicit-product tests enforce fixed preparation and zero warmed CPU allocation. Frozen Gate 0 references preserve historical numerical characterization; the direct-science comparison applies an explicit test-only adapter for the legacy x/y orientation defect, while [`direct_science.jl`](../../test/testsets/direct_science.jl) independently validates the corrected declared-axis semantics. The full CPU suite and maintained AMDGPU and CUDA hardware targets validate the prepared direct-imaging workspace/output boundary on this revision with scalar indexing disabled. |
| HIL-PLANE-003 | 0 | Static aberration and controllable-optic surfaces apply to explicit caller-owned path products without mutating shared telescope path state | Core | implemented | validated | [`apply_surface!`](../../src/optics/planes.jl) accepts OPD maps, NCPA, DMs, and modal optics after distinct controllable-surface formation via [`update_surface!`](../../src/optics/deformable_mirror.jl); [`Explicit optical products and surfaces`](../../test/testsets/core_optics.jl) verifies independent path products and unchanged telescope aperture/path state, and the optional hardware matrix verifies device-resident DM application. |
| HIL-PLANE-004 | 0 | Every prepared optical product declares plane kind, coordinate domain, dimensions, sampling, origin/centering, axis orientation, wavelength/channel, units/normalization, density-versus-cell-integrated measure, coherence/combination policy, numeric type, backend, and physical device where applicable; incompatible handoffs fail before repeated execution | Core | implemented | validated | [`OpticalPlaneMetadata`](../../src/optics/planes.jl) now includes concrete photon-rate/dimensionless normalization, point-sampled/spatial-density/cell-integrated measure, and coherent/incoherent/non-combinable policy markers in addition to the structural, spectral, backend, and physical-device contract. Achromatic geometry, monochromatic products, and application-declared integrated channels are distinct; unspecified spectral metadata is rejected at prepared intensity-combination and detector-acquisition boundaries. Cold propagation, compatible accumulation, and [`prepare_detector_acquisition`](../../src/detectors/detector_acquisition.jl) reject incompatible metadata before repeated execution. [`Optical radiometry and combination contracts`](../../test/testsets/core_optics.jl), [`Detector`](../../test/testsets/detectors.jl), and the [optional backend matrix](../../test/backend_optional_common.jl) provide focused CPU and backend evidence; the full CPU suite and maintained AMDGPU and CUDA hardware matrices validate this revision. Names and units follow the maintained [`glossary`](../glossary.md). |
| HIL-ATM-001 | 0 | One writer advances a shared atmosphere to a stable timed epoch; prepared path-local direction renderers consume that epoch in any order and write caller-owned outputs without shared mutable geometry caches or RNG mutation | Core | implemented | validated | [`source_geometry.jl`](../../src/atmosphere/source_geometry.jl) defines identity-bound epochs and prepared direction renderers; [`Explicit atmosphere epochs and prepared renderers`](../../test/testsets/atmosphere.jl) verifies render-order invariance, stale/incompatible rejection before mutation, RNG preservation, NGS/LGS geometry, and zero warmed CPU allocation. [`control_and_runtime.jl`](../../test/testsets/control_and_runtime.jl) verifies one advance with independent WFS/science and shared-arm renderers. The [`optional backend matrix`](../../test/backend_optional_common.jl) executes the public epoch/renderer and field APIs directly on the maintained AMDGPU and CUDA hardware targets. |
| HIL-ATM-002 | 0 | Atmosphere advancement accepts explicit model time or elapsed duration and never infers time from telescope sampling, detector cadence, wall clock, or a renderer | Core | implemented | validated | `advance_by!` and `advance_to!` in [`source_geometry.jl`](../../src/atmosphere/source_geometry.jl) are the maintained timed-model API. [`ClosedLoopRuntime`](../../src/control/runtime/types.jl) and single/grouped control-loop configurations require a finite positive `atmosphere_step`, preserve it across reconstruction changes, and pass it explicitly to one sensing-cycle atmosphere advance; detector exposure remains independent. Analytic wind-offset, monotonic-time, zero-duration, fixed-seed replay, shared-arm, invalid-step, and exact-runtime-step tests in [`atmosphere.jl`](../../test/testsets/atmosphere.jl) and [`control_and_runtime.jl`](../../test/testsets/control_and_runtime.jl) provide focused CPU evidence; the [optional backend matrix](../../test/backend_optional_common.jl) exercises explicit timed advancement on maintained AMDGPU and CUDA hardware targets. |
| HIL-SRC-001 | 0 | Source definitions are run-immutable execution descriptions; mutable profiles/images are frozen during preparation, and direction, spectral, LGS, asterism, and extended-source expansion are path-local prepared state | Core | implemented | validated | [`source.jl`](../../src/optics/source.jl), [`spectrum.jl`](../../src/optics/spectrum.jl), [`asterism.jl`](../../src/optics/asterism.jl), and [`extended_source.jl`](../../src/optics/extended_source.jl) freeze run-owned descriptions. [`atmosphere.jl`](../../test/testsets/atmosphere.jl) mutates original LGS profiles, spectral arrays, and extended-source images after preparation and verifies unchanged prepared execution; plural renderer preparation expands multi-direction sources without source-owned caches. |
| HIL-RAD-001 | 0 | Physical detector-facing optical products are photon-arrival-rate products with declared units and density-versus-cell-integrated measure; telescope and optical front ends apply no elapsed-time factor, while detector acquisition integrates its explicit whole-exposure or incremental duration exactly once and never creates photons through declared spatial response and pixel integration; finite-support loss is explicit | Core | implemented | validated | Physical and normalized source radiometry in [`source.jl`](../../src/optics/source.jl) propagates through rate- or dimensionless-normalized fields and intensity maps without a telescope-time factor. [`DetectorAcquisitionPlan`](../../src/detectors/detector_acquisition.jl) validates a detector-facing map and prepares density-to-cell scaling and declared-channel QE; repeated whole or incremental capture applies presampling response before physical-pixel integration and applies explicit exposure exactly once. Prepared Shack–Hartmann, Pyramid, and BioEdge front ends write cell-integrated detector-plane photon-rate products without elapsed time. Pyramid/BioEdge modulation weights are normalized optical quadrature weights, not exposure fractions. Distinct wavelength and path-local source products remain an `OpticalProductBundle`. Finite response uses zero extension, so signal leaving detector support is lost rather than clamped back into edge samples. [`core_optics.jl`](../../test/testsets/core_optics.jl), [`detectors.jl`](../../test/testsets/detectors.jl), the prepared [`WFS stage contracts`](../../test/testsets/wfs_stage_contracts.jl), [`shack_hartmann_and_sources.jl`](../../test/testsets/shack_hartmann_and_sources.jl), [`pyramid_bioedge_and_lgs.jl`](../../test/testsets/pyramid_bioedge_and_lgs.jl), and the [optional backend matrix](../../test/backend_optional_common.jl) cover unequal and non-unit exposures, density and cell-integrated fixtures, non-commuting response order, spectral/source bundles, modulation normalization, component accumulation, edge behavior, sampled-QE selection, prepared GPU acquisition, and warmed allocation. The staged physical-family fixtures prove that optical rate products remain unchanged while QE and duration are applied once in acquisition. |
| HIL-RAD-002 | 0 | Spectral/source planes are accumulated elementwise only with compatible physical grids and declared incoherent intensity/rate semantics; otherwise a typed bundle or explicit prepared mapping is required | Core | implemented | validated | [`prepare_incoherent_sum`](../../src/optics/planes.jl) validates immutable metadata once and [`accumulate_intensity!`](../../src/optics/planes.jl) performs prepared allocation-free accumulation; [`OpticalProductBundle`](../../src/optics/planes.jl) preserves incompatible products without implicit conversion or resampling. Unspecified spectral coordinates cannot authorize an intensity sum. Prepared Shack–Hartmann, Pyramid, and BioEdge formation retains distinct wavelength products as bundles. Pyramid/BioEdge Asterism and ExtendedSource formation additionally requires one explicit path-rendered pupil per source and returns a bundle, preventing direction-dependent pupils from being combined by array index. Single-plane and count-mismatched requests fail during preparation. [`core_optics.jl`](../../test/testsets/core_optics.jl) covers compatible accumulation and rejection, while [`wfs_stage_contracts.jl`](../../test/testsets/wfs_stage_contracts.jl) validates wavelength weights, bundle retention, path-local inputs, and structural rejection. |
| HIL-SCI-001 | 0 | Native or prepared external science optics write a caller-owned focal-plane photon-arrival-rate product that one or more independent detector acquisitions integrate without telescope-owned PSF or cadence state | Core | implemented | validated | [`prepare_direct_imaging`](../../src/optics/direct_imaging.jl) binds an explicit caller-owned pupil or preformed field, work field, focal-plane angular `IntensityMap`, and single-writer workspace; [`form_direct_image!`](../../src/optics/direct_imaging.jl) writes a cell-integrated photon-arrival-rate image without exposure time or telescope-owned focal storage. Preparation resolves a finite integer off-axis displacement using the declared focal-grid axis order and signs. Same-grid asterism components are added incoherently, differing spectral grids remain an `OpticalProductBundle`, and extended sources use explicit `extended_source_asterism` expansion. [`direct_science.jl`](../../test/testsets/direct_science.jl) covers point/off-axis orientation, repeated formation after bound-pupil OPD changes, finite mixed-precision radiometry, fan-out to unequal detector exposures, spectral bundling, declared physical-rate and explicitly scaled normalized external-result test doubles, prepared rejection, and warmed CPU allocation. Shared runtimes form one rate image per arm and reuse its exact storage across independent detectors. Focused CPU tests and the maintained AMDGPU and CUDA hardware targets validate point/off-axis formation, spectral bundles, explicit extended-source expansion, detector fan-out, backend/device rejection, and shared-runtime residency with scalar indexing disabled. The committed [Gate 0 artifact](../../benchmarks/results/gate0/2026-07-16-pre-hil-06-science-stage.toml) records zero steady-state allocation and passing predeclared absolute and relative p99 gates for `G0-PERF-02` and `G0-PERF-07`. Asynchronous `Proper.jl` execution remains tracked by `HIL-PATH-004`. |
| HIL-WFS-001 | 0 | Deterministic characterization oracles freeze the current optical, detector-coupled, and estimated products of every maintained WFS family before decomposition | Core | partial | partial | Existing [reference data](../../test/reference_data), [SPECULA-oriented reference data](../../test/reference_data_specula), and WFS tests cover complete products; explicit pre-refactor stage fixtures and LiFT observations remain gaps; tracked by [work issue #2](https://github.com/DarrylGamroth/AdaptiveOpticsSim.jl/issues/2) |
| HIL-WFS-002 | 0 | WFS implementations compose prepared optical-front-end, detector-acquisition, and estimator stages through dispatch over caller-owned products without HIL clock, scheduling, queue, or transport dependencies | Core | partial | partial | [`stage_contracts.jl`](../../src/wfs/stage_contracts.jl) defines typed caller-owned observations/measurements and six narrow prepared dispatch functions for optical formation, acquisition, and estimation, plus explicit acquired/direct measurement traits and structured preparation errors. Shack–Hartmann [`stages.jl`](../../src/wfs/shack_hartmann/stages.jl), Pyramid [`stages.jl`](../../src/wfs/pyramid/stages.jl), and BioEdge [`stages.jl`](../../src/wfs/bioedge/stages.jl) bind caller-owned pupil/field, rate product or bundle, observation, and measurement storage without clock, queue, transport, or runtime-graph dependencies. [`wfs_stage_contracts.jl`](../../test/testsets/wfs_stage_contracts.jl) covers exact binding, stage independence, direct paths, and warmed allocation; the optional backend matrix covers device residency with scalar indexing disabled. Zernike/Curvature and LiFT remain work issues [#6](https://github.com/DarrylGamroth/AdaptiveOpticsSim.jl/issues/6) and [#7](https://github.com/DarrylGamroth/AdaptiveOpticsSim.jl/issues/7), so the cross-family requirement remains partial. |
| HIL-WFS-003 | 0 | Optical front ends produce one or more detector-plane photon-arrival-rate products with an explicit spatial measure; presampling response, physical-pixel integration, QE and explicit elapsed-time integration, coupling, stochastic response, and readout follow optical formation and precede estimation, with one or multiple detector observations supported; the realized response's interior MTF remains available as a derived diagnostic | Core | partial | partial | The generic [`stage contract`](../../src/wfs/stage_contracts.jl) requires detector-plane photon-arrival-rate products and supports static one/multiple-plane and one/multiple-observation composition. Prepared Shack–Hartmann, Pyramid, and BioEdge implementations form cell-integrated rate mosaics or typed bundles before the family-neutral detector acquisition. [`wfs_stage_contracts.jl`](../../test/testsets/wfs_stage_contracts.jl) proves unequal/non-unit exposure, unchanged optical rates, exact-once QE/duration, and non-commuting presampling-response/pixel-integration before typed estimation. [`detectors.jl`](../../test/testsets/detectors.jl) verifies response, pixel integration, QE, exposure, readout, and detector MTF diagnostics. Curvature's separate/packed acquisition forms remain work issue [#6](https://github.com/DarrylGamroth/AdaptiveOpticsSim.jl/issues/6). |
| HIL-WFS-004 | 0 | Shack-Hartmann composes an independent microlens-array optical primitive, detector acquisition, and replaceable spot estimator while preserving geometric, diffractive, spectral, LGS, and asterism behavior | Core | implemented | validated | [`MicrolensArrayParams` and `MicrolensArray`](../../src/wfs/shack_hartmann/setup.jl) describe the immutable regular-array model and sampling policy; `prepare_microlens_propagation` creates separate backend/grid-bound `PreparedMicrolensPropagation` execution state. A component-only `ShackHartmannOpticalFrontEnd` composes that model, propagation, and `SubapertureLayout` without retaining a whole WFS; calibration, acquisition, and estimator storage remain separate. [`stages.jl`](../../src/wfs/shack_hartmann/stages.jl) consumes a caller-owned pupil/field, writes photon-rate mosaics or spectral bundles without telescope mutation or cadence, rejects incompatible plane semantics and placement, accepts detector-declared real observation units, and requires a revision-bound explicit calibration before acquired estimation. Geometric and centroid signals use the same `[axis 1; axis 2]`, Julia-column-major lenslet convention; the OOPAO row-major/axis adapter is confined to the reference harness. [`wfs_stage_contracts.jl`](../../test/testsets/wfs_stage_contracts.jl) exercises the component-only path and checks legacy-equivalent diffractive output, asymmetric axis/lenslet ordering, geometric direct measurement, plane-contract rejection, LGS elongation and sodium profiles, asterism accumulation, spectral bundles, response ordering, exact-once exposure, mixed precision, calibration invalidation, caller-owned mutation, and zero warmed CPU allocation. The optional backend matrix checks device-resident execution and cross-backend rejection with scalar indexing disabled. |
| HIL-WFS-005 | 0 | Pyramid and BioEdge compose independent physical focal-plane front ends, prepared modulation policies, detector acquisition, and calibrated differential estimators while sharing only semantically common machinery | Core | implemented | validated | [`PyramidOpticalFrontEnd`](../../src/wfs/pyramid/setup.jl) binds a physical unit-magnitude phase mask, while [`BioEdgeOpticalFrontEnd`](../../src/wfs/bioedge/setup.jl) binds four distinct amplitude masks; their prepared stage implementations remain separate in [`pyramid/stages.jl`](../../src/wfs/pyramid/stages.jl) and [`bioedge/stages.jl`](../../src/wfs/bioedge/stages.jl). Shared [`focal_plane_modulation.jl`](../../src/wfs/focal_plane_modulation.jl) is limited to normalized no-, circular-, and user-sampled focal-plane modulation policies and path/source helpers. The broader legacy calibration modulation selects support only; zero-aberration references use the operating modulation. Both families form exposure-independent four-pupil photon-rate products or wavelength/path-local bundles, then use the family-neutral detector acquisition and a revision-bound calibrated differential estimator; detector-sampling changes invalidate that binding. Geometric policies intentionally use a direct-measurement path with no focal-plane or detector workspace. [`pyramid_bioedge_and_lgs.jl`](../../test/testsets/pyramid_bioedge_and_lgs.jl) preserves the frozen legacy oracle. [`wfs_stage_contracts.jl`](../../test/testsets/wfs_stage_contracts.jl) validates distinct mask physics, zero/circular/user modulation, exact-once non-unit exposure and QE, atomic calibration and sampling invalidation, spectral and path-local source bundles including heterogeneous NGS/LGS paths, simple elongation and sodium-profile LGS behavior, direct geometric execution, and zero warmed CPU allocation. The [optional backend matrix](../../test/backend_optional_common.jl) exercises both physical families on maintained CUDA and AMDGPU targets with scalar indexing disabled. |
| HIL-WFS-006 | 0 | Zernike and Curvature compose independent optical front ends, detector acquisition, and estimators; Curvature supports both separate branch detectors and packed single-detector observations | Core | planned | missing | Current [Zernike](../../src/wfs/zernike.jl) and [Curvature](../../src/wfs/curvature.jl) states combine these responsibilities. [`zernike_and_curvature.jl`](../../test/testsets/zernike_and_curvature.jl) now records the pre-decomposition oracle: detector-coupled exposure/QE normalization is unit-consistent, source wrappers fail structurally, Curvature optical rate is conserved across its two branches, and mixed-wavelength atmospheric asterisms are rejected. Independent front ends, detector observations, and estimators remain tracked by [work issue #6](https://github.com/DarrylGamroth/AdaptiveOpticsSim.jl/issues/6). |
| HIL-WFS-007 | 0 | LiFT consumes an independently acquired observation and separates its focal-plane forward model from iterative estimation without being forced into the ordinary WFS type hierarchy | Core | partial | partial | [`LiFT`](../../src/wfs/lift.jl) has explicit forward and inverse calculations but currently owns its telescope, source, and detector; tracked by [work issue #7](https://github.com/DarrylGamroth/AdaptiveOpticsSim.jl/issues/7) |
| HIL-WFS-008 | 0 | Geometric and reduced-order WFS policies avoid unused diffractive workspace and declare whether they produce an approximate photon-arrival-rate product for detector processing or intentionally produce a direct measurement | Core | partial | partial | `AcquiredObservationPath()` and `DirectMeasurementPath()` establish the acquired-observation versus intentional-direct-measurement distinction. Geometric Shack–Hartmann, Pyramid, and BioEdge preparation creates neither diffractive-front-end nor detector-acquisition workspace; each explicit-pupil direct estimator is warmed-allocation-free. A package-wide fidelity declaration and migration of the remaining reduced-order WFS modes remain future work across the [Gate 0 series](https://github.com/DarrylGamroth/AdaptiveOpticsSim.jl/issues/1), so the cross-family requirement remains partial. |
| HIL-WFS-009 | 0 | Every migrated stage preserves deterministic CPU results, maintained CUDA and AMDGPU parity, device residency, warmed allocation budgets, and comparable steady-state latency within predeclared absolute and relative gates | Core | partial | partial | Prepared Shack–Hartmann, Pyramid, and BioEdge optical, acquisition, and estimation stages preserve their frozen CPU behavior and are allocation-free after warmup; geometric construction omits unused diffractive and acquisition storage. The shared [optional backend matrix](../../test/backend_optional_common.jl) directly executes all three physical families on device-resident CUDA and AMDGPU arrays with scalar indexing disabled. The committed [Pre-HIL 7 Gate 0 artifact](../../benchmarks/results/gate0/2026-07-16-pre-hil-07-shack-hartmann-stage.toml) records passing absolute and relative latency gates for `G0-PERF-05`; Pre-HIL 8 adds the equivalent staged Pyramid `G0-PERF-06` evidence. This advances three families, but the requirement remains partial until the other WFS migrations and final evidence in [work issue #8](https://github.com/DarrylGamroth/AdaptiveOpticsSim.jl/issues/8). |
| HIL-ORACLE-001 | 1 | Single-threaded deterministic execution is the numerical and event-order oracle | Core | partial | partial | [`reference_harness.jl`](../../test/reference_harness.jl) and the [determinism policy](../deterministic-simulation.md); multi-rate/placed-optic cases remain a gap |
| HIL-PLANT-001 | 2 | One coherent telescope and atmosphere epoch feeds every due path | Core | partial | partial | [`runtime/arms.jl`](../../src/control/runtime/arms.jl) now advances once and renders every arm from the same published epoch with a path-owned renderer; [shared-runtime tests](../../test/testsets/control_and_runtime.jl) verify this and zero warmed CPU allocation. A general due-path scheduler and immutable full-plant snapshot remain gaps. |
| HIL-PATH-001 | 2 | Prepared paths own mutable propagation workspace; acquisitions separately own WFS, detector, readout, and publication state | Core | partial | partial | Direction renderers and [`AtmosphericFieldPropagation`](../../src/optics/atmospheric_field_propagation.jl) now freeze source/geometry and own prepared path workspaces, while runtime arms own their renderers separately from detector state. Full WFS/science front-end and acquisition decomposition remains Gate 2 work. |
| HIL-PATH-002 | 2 | Immutable optical-path definitions are separate from independently scheduled or triggered acquisition state | Core | planned | missing | Enables propagation reuse across cameras and cadences |
| HIL-FID-001 | 2 | Each acquisition binds one run-immutable prepared full-optical, command-responsive reduced-order, or synthetic/replay product provider while preserving one shape, type, sequence, timestamp, lease, port, and overload contract | Core | partial | missing | [`sensing_modes.jl`](../../src/wfs/sensing_modes.jl) and geometric WFS implementations provide reduced slope foundations; full optics exist, but the generic provider seam, synthetic/replay source, and cross-provider conformance evidence are missing |
| HIL-CAL-001 | 2 | Calibration illumination composes through ordinary source/path or detector-input seams with user-declared typed entry, visibility, timing, state, and combination semantics; core assumes no instrument topology, source physics, propagation bypass, or control authority | Core + user model | partial | partial | [`source.jl`](../../src/optics/source.jl), [`spectrum.jl`](../../src/optics/spectrum.jl), [`extended_source.jl`](../../src/optics/extended_source.jl), and [temporal detector sources](../../src/detectors/interface.jl) have maintained [source](../../test/testsets/shack_hartmann_and_sources.jl) and [detector](../../test/testsets/detectors.jl) evidence; a generic internal path-entry seam remains missing |
| HIL-PATH-003 | 3 | Acquisitions over native NGS, finite-height LGS, and direct-science paths use independent virtual-time schedules | Core | partial | partial | [`source_geometry.jl`](../../src/atmosphere/source_geometry.jl) and [`runtime/arms.jl`](../../src/control/runtime/arms.jl) provide frozen NGS/LGS directions and one stable epoch for independent WFS/science paths; atmosphere and runtime tests cover distinct geometry. Common multi-rate virtual-time orchestration remains a gap. |
| HIL-DET-001 | 3 | Exposure, optical samples, rolling shutter, frame transfer, nondestructive reads, readout, presampling detector response, charge coupling, and complete-product publication have explicit composable event semantics; the realized response has a validated interior MTF where supported and explicit finite-frame boundary behavior | Core | partial | partial | [`detectors/interface.jl`](../../src/detectors/interface.jl) and [detector tests](../../test/testsets/detectors.jl); scheduled plant integration remains a gap |
| HIL-TIME-001 | 3 | Integer canonical plant time and stable ordinals determine all equal-time event ordering | Core | planned | missing | Includes half-open exposure semantics |
| HIL-TRIG-001 | 3 | A prepared trigger-source and distribution topology maps nominal edges to per-acquisition delivery times with explicit phase/skew, jitter, dropped/duplicate-edge policies, and separate physical exposure and timestamp-label semantics | Core | planned | missing | Camera-internal oscillator models are optional; the baseline covers externally visible trigger behavior |
| HIL-SCHED-001 | 3 | A fixed-capacity scheduler jumps to due timestamps without base-tick polling, run-length-sized event materialization, or warmed allocation | Core | planned | missing | A prepared linear scan is the baseline; a stable array heap requires measured scale evidence |
| HIL-CMD-001 | 4 | Every independently timed controllable optic or segment has an independent virtual-time command endpoint | Core | planned | missing | Replaces composite command packing; the minimal HIL port arrives in Gate 4A |
| HIL-CMD-002 | 4 | Validation, admission, effective time, application, hold, and terminal outcome are distinct | Core | planned | missing | Atomic multi-optic latch is explicit, never inferred from placement |
| HIL-CMD-003 | 4 | Future, late, rejected, and superseded commands use bounded explicit policies and never backdate plant state | Core | planned | missing | Every displaced or rejected command has a terminal outcome |
| HIL-CMD-004 | 4 | Every endpoint prepares a versioned canonical payload schema covering type/shape, units, basis/calibration, absolute or incremental semantics, range policy, session epoch, and duplicate/reordering behavior | Core + HIL | planned | missing | Transport decoding maps into this schema before semantic admission |
| HIL-CMD-005 | 4 | Every endpoint declares command-silence behavior, separating replayable plant-time hold/safe/fail policy from an optional execution-time ingress-liveness failure | Core + HIL | planned | missing | Safe state and watchdog storage are prepared; no watchdog callback allocates or invokes transport code |
| HIL-FID-003 | 4 | A reduced-order plant evolves time-correlated disturbances and maps correctly timed effective optic commands through calibrated path and sensor operators so commands causally change later products and a matched reference controller can close the loop | Core | partial | missing | [`interaction_matrix.jl`](../../src/calibration/interaction_matrix.jl), [`reconstructors.jl`](../../src/control/reconstructors.jl), geometric WFS modes, and tomography primitives are foundations; the prepared reduced-order plant and acceptance evidence are missing |
| HIL-OPT-001 | 4 | Co-conjugated DMs remain independent devices even when prepared as one optical execution group | Core | partial | partial | [`controllable_optics.jl`](../../src/optics/controllable_optics.jl) and [runtime tests](../../test/testsets/control_and_runtime.jl) cover additive composites only |
| HIL-OPT-004 | 4 | Autonomous periodic optical devices use trigger-relative waveform state and bounded setpoint commands rather than one RTC message per waveform point | Core | partial | partial | [`PyramidWFS`](../../src/wfs/pyramid/setup.jl) and its [modulation optics](../../src/wfs/pyramid/optics.jl) provide the baseline cycle-averaged foundation; physical trigger relation remains missing and time-resolved fidelity is profile-driven |
| HIL-VSLICE-001 | 4A | A serial CPU vertical slice connects one scheduled acquisition and one command-responsive optic to a deterministic in-memory fake RTC through the canonical complete-product and command/outcome boundary with an injected clock and fixed-arrival evidence | HIL + Core | planned | missing | Proves the external contract before worker, GPU, multi-path, or transport-specific optimization |
| HIL-PORT-001 | 4A | Canonical command-submission, command-completion/outcome, and acquisition-completion ports are transport-neutral | HIL | planned | missing | Successful ring enqueue transfers ownership but is not semantic command admission |
| HIL-PORT-003 | 4A | External integration declares adapter readiness, complete-product-to-first-observation lead time, and maximum lease hold time without promoting progressive fragments into canonical simulation events | HIL + user integration | planned | missing | First-packet timing remains an external-delivery claim |
| HIL-PORT-004 | 4A | Sampled optic/device feedback uses ordinary acquisition endpoints and remains distinct from correlated command outcomes | Core + HIL | planned | missing | Supports scalar, vector, surface, encoder, and health products without an instrument-specific API |
| HIL-BUF-001 | 4A | Large acquisition and command payloads use bounded pools and explicit generation-checked leases rather than ring-resident payload copies | HIL + user integration | planned | missing | Memory domain, generation, terminal-outcome/release ownership, and reuse are explicit; return-credit hardening follows in Gate 8 |
| HIL-OPT-002 | 5 | Common MCAO and path-specific MOAO planes support NGS and finite-height LGS footprints | Core | planned | missing | Requires placement and visibility traits |
| HIL-OPT-003 | 5 | Common and path-specific sampled aberrations, including NCPA, have explicit visibility | Core | partial | partial | [`ncpa.jl`](../../src/optics/ncpa.jl), [`opd_map.jl`](../../src/optics/opd_map.jl), and [model tests](../../test/testsets/calibration_and_analysis.jl); path visibility remains a gap |
| HIL-CPU-001 | 6 | Prepared CPU execution groups have single-writer workspaces, explicit thread budgets, independently callable executor seams, and a deterministic serial fallback | Core | planned | missing | Long-lived HIL agents and their rings arrive in Gate 8 |
| HIL-GPU-001 | 7 | Physical device identity, prepared direction batching, and device-resident state are explicit | Core | partial | partial | [CUDA](../../ext/AdaptiveOpticsSimCUDAExt.jl), [AMDGPU](../../ext/AdaptiveOpticsSimAMDGPUExt.jl), and [backend tests](../../test/backend_optional_common.jl) now validate device-resident epoch rendering and prepared atmospheric fields with physical-device compatibility checks; prepared batching and physical multi-device planning remain gaps. |
| HIL-TIME-002 | 8 | External timestamp domains map into plant time with versioned offset, drift, uncertainty, and no retroactive remapping | HIL | planned | missing | User integration supplies synchronization observations |
| HIL-LIFE-001 | 8 | Configure, prepare, arm, run, and stop/fail phases preserve ownership, adapter-readiness preconditions, prepared nonstructural state transitions, and bounded shutdown | HIL | planned | missing | Topology, schema, capacity, placement, and provider changes require another prepare/arm cycle |
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
| 0 | Optical planes/radiometry, explicit atmosphere time, direct-science acquisition, and WFS stages decomposed and validated | partial |
| 1 | Contracts and correctness oracles frozen | partial |
| 2 | Shared plant, acquisition workspace, product providers, and source-entry seams separated | planned |
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

## Capability Gates

The target advances through durable capability gates. A PR may satisfy part or
all of a gate, but references the applicable requirement IDs and updates their
implementation and evidence state. Every gate preserves numerical correctness
and the serial oracle. Public API compatibility is not required for the
independent-optic refactor; superseded surfaces are removed rather than kept as
permanent adapters.

No implementation of the proposed general `AdaptiveOpticsHIL.jl` companion
runtime begins until Gate 0 is complete. Existing integration proving grounds
remain valid, and HIL architecture and interface design may continue while
this prerequisite is implemented. The ordered work is tracked by
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
- retain the direct single-threaded CPU runtime as the numerical oracle

Acceptance: every changed physical surface has a deterministic reference or an
explicitly tracked evidence gap; current correctness, allocation, and backend
targets remain green.

### Gate 2: Separate shared plant models from path and acquisition state

- separate immutable path definitions from independently scheduled or
  triggered acquisition endpoints
- render atmosphere state into caller-owned destinations
- precompute source geometry and compatible-result keys per prepared path
- separate telescope parameters and pupil from path-local propagation workspace
  and acquisition-local WFS, detector, readout, and publication state
- bind full-optical, reduced-order, and synthetic/replay providers through one
  prepared mutating acquisition-product seam using dispatch and traits
- expose a narrow prepared calibration-illumination seam at supported typed
  path entries or detector inputs, leaving physical integration to user models
- keep the first executor serial to isolate ownership correctness

Acceptance: frozen outputs remain within declared tolerances, every due path
sees one atmosphere epoch, path reuse does not couple acquisition state, and
every applicable provider preserves the acquisition's shape, type, sequence,
timestamp, lease, port, and overload contract. Reduced-order providers declare
their validity envelope; static/replay providers are explicitly nonresponsive;
and the selected provider remains immutable until another prepare/arm cycle.
A native and a user-provided calibration source preserve declared visibility,
timing, deterministic state, and composition without a calibration-role branch
or implicit bypass. Unsupported entry payloads and source combinations fail
during preparation. The warmed serial oracle meets its declared allocation
budget. Comparable latency is archived as a baseline; no subjective latency
gate is used.

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
- preserve presampling detector response, charge-coupling, and readout-pipeline ordering when
  acquisition is split into scheduled events
- accept explicit simulation timestamps without reading wall clock in core
- expose due/readiness products without transport dependencies

Acceptance: exact-boundary tests cover commands and acquisition events sharing
a timestamp. Trigger tests distinguish nominal source edges, delivered edges,
physical exposure boundaries, reported labels, and execution time; cover fixed
skew, correlated and independent jitter, phase steps, and dropped/duplicate
edges; and prove that faults affect only their declared downstream branches.
Representative CMOS rolling/global shutter, CCD, frame-transfer EMCCD, and
HgCdTe nondestructive-read cases preserve their frozen detector products,
including configured response/coupling behavior. Mixed rates produce the expected
sequences, long runs use storage proportional to active generators rather than
event count, and long-period schedules do not iterate empty base ticks.
Scheduler cost is archived across the maintained generator-count range, and
warmed execution meets its allocation budget.

### Gate 4: Replace composite optics with independent command semantics

- replace the single `AOSimulation.optic` field with a named registry of
  individual controllable optics
- introduce one independently timed endpoint per optic or latched segment
- prepare a versioned canonical command schema per endpoint, including units,
  basis/calibration revision, payload semantics, bounds, run epoch, and
  duplicate/reordering policy
- separate validation, admission, effective time, application, hold, and
  terminal outcome in the virtual-time core; add enqueue in Gate 4A
- define bounded future-command and explicit late/supersession policies
- define command-silence hold/safe/fail behavior in plant time, distinct from
  any operational execution-time ingress-liveness failure
- support explicit atomic multi-optic transactions without inferring atomicity
  from plane placement or packed vectors
- remove `CompositeControllableOptic` and split `RuntimeCommandLayout` into
  prepared core routing, canonical HIL descriptors, and user transport schemas
- split detector integration correctly when a command becomes effective
- add a prepared reduced-order plant whose time-correlated disturbances,
  path/sensor operators, and effective commands produce causally correct slopes
  or approximate raw pixels for an external RTC
- represent autonomous periodic optics with trigger-relative waveforms and
  bounded setpoint commands, retaining the fast cycle-averaged pyramid policy;
  defer time-resolved modulation until a profile requires it

Acceptance: deterministic timelines cover early, equal-time, future, late,
rejected, superseded, and atomic commands; every submitted virtual-time command
has one terminal outcome, and every admitted command is applied once or ends in
a declared failure; schema/session/shape/sequence mismatches and both watchdog
clock domains have exact-boundary tests; independent optics update without an
internal RTC; the cycle-averaged pyramid model preserves its frozen reference
and trigger relationship; a matched reference controller closes the
reduced-order loop and reduces its declared residual while wrong-sign, delayed, stale, and
mismatched-calibration cases degrade as expected; and frozen physical outputs
remain within tolerance.

### Gate 4A: Prove a minimal serial HIL vertical slice

- create the transport-free `AdaptiveOpticsHIL.jl` package boundary with only
  the canonical types needed by one serial path
- inject a deterministic test clock and a monotonic production clock without
  adding wall-clock reads to core
- connect one scheduled acquisition to a complete-product completion port and
  one command-responsive optic to a submission/outcome pair
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
- preserve a deterministic single-threaded fallback
- declare Julia, FFT, and BLAS thread ownership and avoid nested parallelism
- expose stable worker placement and NUMA requirements

Acceptance: each group is independently callable under a validation harness;
serial/grouped parity and allocation gates pass; and declared thread budgets do
not oversubscribe Julia, FFT, or BLAS execution. Any performance promotion
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
- inject `Clocks.jl` sources and versioned external-domain mappings while
  keeping execution-clock pacing separate from modeled trigger delivery
- complete configure, prepare, arm, run, and bounded stop/fail lifecycle,
  including adapter-readiness preconditions and prepared nonstructural
  acquisition, trigger, shutter/calibration, and safe/hold transitions
- instantiate long-lived execution owners from prepared CPU groups and a single
  submission owner for each selected GPU
- harden the canonical command-submission, command-completion/outcome, and
  complete-product acquisition-completion ports for every prepared owner
- assign one prepared command authority per endpoint or explicit atomic latch
  group; keep competing-producer arbitration outside the canonical data plane
- implement padded bounded SPSC descriptor rings and bounded product pools
- use one due-work/completion path per owner and a fixed pool of versioned epoch
  snapshots rather than a contended global queue
- define resource-specific full, close, drain, lease, and recovery semantics
- reserve one usable return credit for every lease a consumer can hold and
  continuously check the complete pool-accounting invariant
- give every execution owner a preallocated first-failure/acknowledgement path;
  let one coordinator close ingress, publish stop, and drain ownership within
  declared deadlines
- establish same-process allocation and GC budgets, with an explicit process-
  isolation rule for adapters or telemetry that cannot meet them
- provide deterministic in-memory adapters and port conformance tests

Acceptance: exact execution-clock and external-domain mapping tests pass;
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

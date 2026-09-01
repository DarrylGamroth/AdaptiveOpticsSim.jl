# Model Validity Matrix

Status: active

## Purpose

This matrix connects maintained scientific claims to current implementation and
tests. It does not promote a model merely because code or an old benchmark
artifact exists.

## Evidence Classes

| Code | Evidence |
|---|---|
| `A` | analytic or invariant contract |
| `R` | external reference, OOPAO trace, or instrument-owned input |
| `G` | real accelerator execution |
| `P` | measured performance/allocation evidence |
| `M` | maintained integration or composed-model test |

“Strong” means strong only inside the stated model envelope. “Partial” names a
real remaining qualification gap.

## Matrix

| ID | Surface | Evidence | Current boundary | Assessment |
|---|---|---|---|---|
| `MV-01` | Core optical products and propagation | `A, R, M` — [`core_optics.jl`](../test/testsets/core_optics.jl), [`direct_science.jl`](../test/testsets/direct_science.jl), [`reference_harness.jl`](../test/reference_harness.jl) | Explicit planes, coordinates, normalization, spatial measure, pupil/field products, direct imaging, and selected frozen references. No universal coronagraph or instrument relay claim. | strong |
| `MV-02` | Atmosphere evolution and direction rendering | `A, R, G, M` — [`atmosphere.jl`](../test/testsets/atmosphere.jl), [`atmosphere_direction_batch.jl`](../test/testsets/atmosphere_direction_batch.jl), [`backend_optional_common.jl`](../test/backend_optional_common.jl) | Finite/infinite screens, explicit epochs, NGS/LGS direction geometry, and homogeneous batches. Broad external profile equivalence remains scenario-specific. | strong within implemented models |
| `MV-03` | WFS optical and estimation families | `A, R, G, M` — [`wfs_common_and_parity.jl`](../test/testsets/wfs_common_and_parity.jl), [`shack_hartmann_and_sources.jl`](../test/testsets/shack_hartmann_and_sources.jl), [`pyramid_bi_o_edge_and_lgs.jl`](../test/testsets/pyramid_bi_o_edge_and_lgs.jl), [`zernike_and_curvature.jl`](../test/testsets/zernike_and_curvature.jl), [`wfs_lift.jl`](../test/testsets/wfs_lift.jl) | Separate optical formation, acquisition, and estimation; family-specific outputs and calibration. Instrument-specific optical gain and non-common-path calibration remain external. | strong for declared numerical models |
| `MV-04` | Detector families and acquisition | `A, G, M` — registered `detector-*` suites under [`test/testsets`](../test/testsets), [`backend_optional_common.jl`](../test/backend_optional_common.jl) | CCD, EMCCD, CMOS, InGaAs, HgCdTe, avalanche, Skipper, SPAD, MKID, and linear APD envelopes as tested. Instrument graph files may bind source-backed camera settings, as the REVOLT Classic IMX425 graph does, but reusable named commercial profiles, unit-specific calibration, and untested long-timescale effects are excluded. | strong by family; calibration partial |
| `MV-05` | Calibration and control | `A, R, M` — [`calibration_workflows.jl`](../test/testsets/calibration_workflows.jl), [`control_primitives.jl`](../test/testsets/control_primitives.jl), [`control_reconstruction.jl`](../test/testsets/control_reconstruction.jl) | Matrix/basis construction, inverse policies, reconstruction, controller state, delays, and explicit closed-loop composition. Production RTC software is not implied. | strong |
| `MV-06` | Tomography | `A, M` — [`tomography.jl`](../test/tomography.jl) | Guide-star/layer geometry, model- and interaction-matrix reconstruction, fitting, and DM projection. Instrument MCAO/MOAO calibration remains scenario work. | implemented; external equivalence partial |
| `MV-07` | CPU/CUDA/AMDGPU backend execution | `G, P, M` — [`ka_cpu_matrix.jl`](../test/ka_cpu_matrix.jl), [`backend_optional_common.jl`](../test/backend_optional_common.jl), [`backend-validation-guide.md`](backend-validation-guide.md) | Exact-device storage and covered algorithms with scalar indexing disabled. Support is per operation; no implicit CPU fallback, multi-GPU placement, or universal bitwise parity. | strong for tested targets |
| `MV-08` | Static complete-frame algorithm graphs | `A, P, M` — [`algorithm_graphs.jl`](../test/testsets/algorithm_graphs.jl), [`runtests_gpu_target_common.jl`](../test/runtests_gpu_target_common.jl) | Concrete nodes, exact port contracts, direct/delayed links, TOML loading, sparse parameters, preparation, capacity-one completion-ticket ownership, reset, model time, failure state, and zero warmed allocation for covered paths. CUDA Graph/HIP Graph replay is opt-in and captures one complete step only when every node owner qualifies; the maintained captured graph contains the regular-grid DM and pupil-OPD composition nodes. Multiple in-flight frames and full REVOLT capture qualification are excluded. | strong |
| `MV-09` | External-RTC lockstep boundary | `A, P, M` — HIL testsets in [`algorithm_graphs.jl`](../test/testsets/algorithm_graphs.jl) and optional [`pyRTC integration tests`](../examples/integrations/pyrtc/runtests.jl) | One complete host frame followed by one finite same-sequence command, explicit host/device copies, reset and failure containment. Native Linux SHM interoperation covers one outstanding frame, static command recovery, and corrected on-axis Strehl improvement under deterministic evolving atmosphere for the compact references and the independently enabled 352-by-352/376-signal/277-command REVOLT Classic model. No asynchronous overwrite safety, wall-clock pacing, or fixed-arrival latency claim. | strong for transport-neutral and tested pyRTC lockstep |
| `MV-10` | REVOLT Classic and Copper graph examples | `R, P, M` — [`revolt_classic_hil.toml`](../examples/graphs/revolt_classic_hil.toml), [`revolt_copper_hil.toml`](../examples/graphs/revolt_copper_hil.toml), executable tests in [`algorithm_graphs.jl`](../test/testsets/algorithm_graphs.jl), and optional REVOLT Classic [pyRTC process validation](../examples/integrations/pyrtc/test_revolt_classic_hil.jl) | Complete simulated detector-frame paths with atmosphere and PDM command input. The Classic model self-calibrates and closes through pyRTC with the maintained dimensions. The HSDM277 influence function is provisional and the examples are not instrument qualification. | executable; scientific qualification partial |
| `MV-11` | Optional Proper.jl composition | `G, P, M` — [`proper_hil_coronagraph_common.jl`](../examples/support/proper_hil_coronagraph_common.jl), [`proper-integration-guide.md`](proper-integration-guide.md) | Explicit prepared application or graph-node boundary with typed optical products. Current evidence does not validate a SPIDERS prescription or coronagraph fidelity. | integration strong; science partial |
| `MV-12` | Deterministic simulation ownership | `A, M` — [`deterministic-simulation.md`](deterministic-simulation.md), atmosphere/detector/graph reset tests | Explicit RNGs, one writer per evolving stream, fixed seeds, reset, and single-threaded CPU reference policy. Cross-hardware bitwise identity and distributed addressable RNG are excluded. | strong for serial CPU policy |

## Known Scope Limits

- OOPAO is a source of reference data and tutorial mappings, not an API or
  class-layout compatibility target.
- Geometric and diffractive WFS modes answer different modeling questions and
  are not expected to be samplewise identical.
- Full-optical and reduced-order providers must be named and validated
  separately.
- CPU and GPU results may differ within declared floating-point or statistical
  tolerances.
- A rolling shutter may require several optical evaluations per completed
  frame; partial row products are not graph outputs.
- MCAO and MOAO topology is supported by explicit Julia composition and
  underlying optical/control primitives, but a named instrument model needs its
  own calibration and validation evidence.
- Historical artifacts that exercised the removed event runtime characterize
  only their recorded revision.

## Acceptance Rule

Update this matrix when a public scientific surface, graph node, backend claim,
or model limitation changes. A new strong claim must point to current tests and
appropriate scientific evidence; prose or a historical benchmark alone is not
enough.

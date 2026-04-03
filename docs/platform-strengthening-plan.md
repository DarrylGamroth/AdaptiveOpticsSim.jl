# Platform Strengthening Plan

Status: active

Traceability sources:

- [future-platform-direction.md](./future-platform-direction.md)
- [package-review-2026-03.md](./package-review-2026-03.md)
- [model-validity-matrix.md](./model-validity-matrix.md)
- [execution-plan-closeout.md](./execution-plan-closeout.md)

## Purpose

This plan governs the next package phase after the execution-plan rollout.

Its purpose is to address the remaining high-value platform gaps without
drifting into config-first design, parity chasing, or premature science-path
integration in core.

This plan is intentionally focused on:

- API shape and synthesis-oriented usability
- stronger validation and benchmark evidence
- further ROCm/AMDGPU runtime hardening
- richer platform-scale orchestration and composition

This plan intentionally does not make declarative manifest/config files the
primary interface.

## Direction Requirements

These IDs are the stable requirements for this plan.

| ID | Requirement | Source |
| --- | --- | --- |
| `PSR-01` | Scripts and typed Julia builders remain the primary interface. | [future-platform-direction.md](./future-platform-direction.md) |
| `PSR-02` | Config-file-first architecture is out of scope for this phase. | user direction, [future-platform-direction.md](./future-platform-direction.md) |
| `PSR-03` | Address remaining review gaps before broad new capability expansion. | [package-review-2026-03.md](./package-review-2026-03.md) |
| `PSR-04` | Strengthen validation and representative benchmark evidence where current claims are still narrower than ideal. | [model-validity-matrix.md](./model-validity-matrix.md) |
| `PSR-05` | Reduce backend-specific workaround debt where it still affects maintained GPU paths. | [execution-plan-closeout.md](./execution-plan-closeout.md) |
| `PSR-06` | Use SPECULA as the stronger reference for platform-scale breadth, but keep the implementation Julia-native. | [future-platform-direction.md](./future-platform-direction.md) |
| `PSR-07` | Science-path / focal-plane integration remains outside core unless a later explicit boundary review changes that rule. | [future-platform-direction.md](./future-platform-direction.md) |

## Primary Interface Policy

The maintained interface policy for this phase is:

- primary:
  - Julia scripts
  - typed constructors
  - typed scenario/build orchestration objects
  - benchmark and validation runners
- deferred:
  - optional scenario manifests
  - config-file-driven primary workflows

The package should prefer:

- programmable Julia-native composition,
- reproducible script entry points,
- and typed scenario objects

over:

- opaque declarative configuration as the main user experience.

If a manifest/config format is revisited later, it should sit on top of the
typed platform model, not replace it.

## Current Gaps

The remaining gaps this plan is meant to address are:

1. the public API is still broader and flatter than ideal
2. documentation still needs more synthesis-oriented platform guidance
3. tomography representative benchmark evidence is still explicitly deferred
4. SPECULA-aligned validation is strong, but still largely contract-oriented
5. some maintained host-mirror fallback behavior still exists on AMDGPU
6. platform-scale orchestration/composition is still weaker than the strongest
   SPECULA-style platform surfaces
7. science-path / focal-plane integration should be considered, but only as an
   optional boundary, not as new core scope

## Non-Goals

This plan does not directly implement:

- config-file-first scenario control
- broad controller/process-family expansion
- core-package science-camera or coronagraph integration
- parity work justified only because OOPAO or SPECULA have a feature

## Phase Overview

| Phase | Goal | Exit gate |
| --- | --- | --- |
| `PH-1` | tighten platform usability and synthesis docs | API layering and synthesis docs are materially improved without breaking maintained workflows |
| `PH-2` | strengthen evidence where current claims are still narrow | tomography and SPECULA-aligned evidence gaps are either closed or explicitly re-scoped |
| `PH-3` | reduce remaining maintained AMDGPU workaround debt | benchmark-backed ROCm fallback inventory is reduced or explicitly justified |
| `PH-4` | introduce a Julia-native platform orchestration layer | typed scenario/composition surfaces exist for richer multi-branch runtime assembly |
| `PH-5` | validate platform orchestration on realistic maintained scenarios | new orchestration surfaces have benchmark and validation evidence on maintained backends |
| `PH-6` | close out and record optional deferred tracks | next-step decision is explicit and science-path/config tracks remain intentionally deferred |

## Phase 1: Platform Usability And Synthesis

### Goal

Make the platform easier to understand and use without changing the primary
Julia-native workflow model.

### Exit Criteria

Phase 1 is complete when:

- the next API-curation tranche is implemented and recorded in the maintained
  API inventory
- a user can move from “I know individual subsystems” to “I understand how to
  build and run a realistic platform scenario” using stable synthesis docs
- the recommended path for reproducible scenario authoring is documented as
  script-first, typed, and Julia-native
- affected docs and tests are updated so the new guidance is aligned with the
  maintained implementation

### Phase Outputs

Phase 1 should produce:

- one additional export-curation patch with updated traceability
- one stable architecture/synthesis guide for platform-scale usage
- one stable workflow guide describing how the main runtime pieces fit together
- one scenario-authoring style guide for reproducible script-based platform
  scenarios

### Recommended Execution Order

1. `PSP-01`
2. `PSP-02`
3. `PSP-03`
4. `PSP-04`

### Milestones

| Milestone | Task | Outputs |
| --- | --- | --- |
| `PSP-01` | Curate the next tranche of low-value exports into namespaced access. | export reduction patch, updated API tier inventory |
| `PSP-02` | Add a synthesis-oriented platform architecture guide aimed at users moving beyond single tutorials. | new stable guide |
| `PSP-03` | Add a platform workflow guide explaining how atmosphere, runtime, WFS, detector, and benchmark entry points fit together. | new stable guide or expanded user guide section |
| `PSP-04` | Add a scenario-builder style guide that defines how Julia scripts should express reproducible scenarios without config files. | style guide / conventions note |

### `PSP-01`: Next API-Curation Tranche

Status:

- completed

Objective:

- remove another low-value tranche of top-level exports while preserving the
  stable top-level experience for the most common constructors and workflows

Primary files in scope:

- [src/AdaptiveOpticsSim.jl](../src/AdaptiveOpticsSim.jl)
- [docs/api-tier-inventory.md](./api-tier-inventory.md)
- [docs/api-reference.md](./api-reference.md)
- [docs/user-guide.md](./user-guide.md)
- affected examples and tests that still depend on the candidate exports

Detailed work:

- re-scan current exports and identify the next low-risk tranche of names that
  are:
  - helper-level,
  - calibration-internal,
  - benchmark-only,
  - or better consumed via `AdaptiveOpticsSim.<name>`
- record that tranche in [api-tier-inventory.md](./api-tier-inventory.md)
- de-export the chosen tranche in
  [src/AdaptiveOpticsSim.jl](../src/AdaptiveOpticsSim.jl)
- update affected docs, tutorials, examples, and tests to use explicit
  namespaced access
- add or update at least one regression test that confirms the intended stable
  top-level API remains available

Deliverables:

- export-curation code patch
- updated API inventory with before/after counts
- updated API docs/examples reflecting the curated surface

Acceptance criteria:

- export count is reduced again, or the remaining tranche is explicitly
  reclassified with evidence that further reduction would hurt usability
- no documented stable quickstart path becomes more obscure
- `Pkg.test()` passes after the export changes

Evidence to record:

- export count before/after
- list of de-exported names or categories
- any intentionally retained borderline names and the reason

Implementation record:

- pre-`PSP-01` export count: `544`
- post-`PSP-01` export count: `538`
- de-exported tranche:
  - `BuildBackend`
  - `NativeBuildBackend`
  - `CPUBuildBackend`
  - `GPUArrayBuildBackend`
  - `default_build_backend`
  - `set_fft_provider_threads!`
- intentionally retained for a later coherent pass:
  - GPU backend tags
  - GPU precision-policy helpers
  - low-level backend allocation helpers

### `PSP-02`: Platform Architecture Synthesis Guide

Status:

- completed

Objective:

- give advanced users and maintainers one stable orientation guide that explains
  the package as a platform rather than as a pile of subsystem docs

Primary files in scope:

- new stable guide, expected path:
  [docs/platform-architecture.md](./platform-architecture.md)
- [docs/documentation-map.md](./documentation-map.md)
- [docs/user-guide.md](./user-guide.md)
- [docs/maintainer-architecture.md](./maintainer-architecture.md)
- [docs/runtime-dataflow.md](./runtime-dataflow.md)

Detailed work:

- write a synthesis-oriented architecture guide covering:
  - platform building blocks,
  - params/state/runtime ownership,
  - orchestration boundaries,
  - validation surfaces,
  - backend strategy,
  - and extension boundaries
- explicitly connect single-subsystem usage to platform-scale composition
- keep the guide descriptive, not roadmap-shaped
- link out to the deeper subsystem docs instead of duplicating them
- update [documentation-map.md](./documentation-map.md) so the new guide becomes
  a first-class stable entry point

Minimum section schema:

- purpose
- platform model
- major subsystem families
- how runtime composition works today
- how validation and benchmarks support claims
- extension boundaries and non-goals
- recommended reading path

Deliverables:

- new stable platform architecture guide
- updated documentation map and any necessary cross-links

Acceptance criteria:

- a reader can understand the package-level architecture without reading plan
  docs first
- the guide complements rather than duplicates
  [maintainer-architecture.md](./maintainer-architecture.md)
- doc navigation exposes the guide as a stable entry point

Evidence to record:

- guide path
- updated documentation-map entry
- any replaced or superseded orientation references

Implementation record:

- stable guide added:
  [platform-architecture.md](./platform-architecture.md)
- stable navigation updated in:
  [documentation-map.md](./documentation-map.md)
- user and maintainer entry points updated in:
  [user-guide.md](./user-guide.md),
  [maintainer-architecture.md](./maintainer-architecture.md), and
  [runtime-dataflow.md](./runtime-dataflow.md)

### `PSP-03`: Platform Workflow Guide

Status:

- completed

Objective:

- describe the main script-first workflow from scenario setup through runtime
  execution, validation, and benchmarking

Primary files in scope:

- new stable guide, expected path:
  [docs/platform-workflows.md](./platform-workflows.md)
- [docs/user-guide.md](./user-guide.md)
- [docs/runtime-dataflow.md](./runtime-dataflow.md)
- [docs/benchmark-matrix-plan.md](./benchmark-matrix-plan.md)
- [docs/model-validity-matrix.md](./model-validity-matrix.md)

Detailed work:

- write one workflow-oriented guide for:
  - defining sources/telescope/atmosphere/WFS/detector/runtime
  - building a reproducible script
  - running validation/reference comparisons
  - running maintained benchmarks
  - choosing grouped/platform-scale versus isolated subsystem workflows
- include concrete references to the maintained scripts and examples that should
  be copied or adapted
- make clear when users should use:
  - tutorials,
  - stable APIs,
  - benchmark scripts,
  - or cross-package harnesses

Deliverables:

- new stable workflow guide or equivalent substantial user-guide expansion
- updated navigation links

Acceptance criteria:

- a user can follow the guide to construct a realistic script-based platform
  workflow without scanning multiple plan docs
- the guide clearly distinguishes:
  - normal usage,
  - validation usage,
  - and benchmarking usage

Evidence to record:

- guide path
- referenced maintained scripts/examples
- any updated stable navigation entries

Implementation record:

- stable guide added:
  [platform-workflows.md](./platform-workflows.md)
- stable navigation updated in:
  [documentation-map.md](./documentation-map.md)
- user-facing entry points updated in:
  [user-guide.md](./user-guide.md) and
  [platform-architecture.md](./platform-architecture.md)

### `PSP-04`: Scenario-Builder Style Guide

Status:

- completed

Objective:

- define the package’s maintained script-first conventions for reproducible,
  typed scenario construction

Primary files in scope:

- new stable/supporting guide, expected path:
  [docs/scenario-builder-style.md](./scenario-builder-style.md)
- [docs/user-guide.md](./user-guide.md)
- [docs/future-platform-direction.md](./future-platform-direction.md)
- relevant examples/tutorials if a canonical pattern needs to be updated

Detailed work:

- write the style rules for scenario scripts, including:
  - params/state separation in scripts
  - deterministic RNG ownership
  - scenario naming and decomposition
  - builder/helper function boundaries
  - result capture and export expectations
  - benchmark/validation hooks
  - what should not be embedded directly in scripts
- show a canonical script layout for:
  - a normal platform run,
  - a validation run,
  - and a benchmark run
- explicitly document that scenario manifests remain deferred and secondary

Deliverables:

- maintained style guide for script-based scenario construction
- updated cross-links from user-facing docs

Acceptance criteria:

- the style guide makes script-first platform composition more uniform
- the package stance against config-first primary workflows is explicit
- later orchestration work can reference this guide instead of inventing new
  script conventions ad hoc

Evidence to record:

- guide path
- canonical script layout pattern recorded in the guide
- any examples updated to align with the new rules

Implementation record:

- stable guide added:
  [scenario-builder-style.md](./scenario-builder-style.md)
- stable navigation updated in:
  [documentation-map.md](./documentation-map.md)
- user/workflow entry points updated in:
  [user-guide.md](./user-guide.md) and
  [platform-workflows.md](./platform-workflows.md)

### Acceptance

- the maintained top-level API is smaller or at least better tiered
- new platform docs explain composition more directly than plan docs do
- the recommended user path is still scripts + typed objects, not config files

## Phase 2: Evidence Strengthening

### Goal

Close the highest-value remaining validation and benchmark gaps.

### Exit Criteria

Phase 2 is complete when:

- tomography benchmark scope is either strengthened with a maintained artifact
  or explicitly re-scoped with current benchmark-backed evidence
- SPECULA-aligned evidence extends beyond isolated optical contracts into at
  least one broader runtime/platform-scale scenario
- all new evidence is reflected in the maintained validity/benchmark docs and
  stored as reproducible artifacts or explicit decision records

### Phase Outputs

Phase 2 should produce:

- one tomography benchmark artifact or defer note with current justification
- expanded SPECULA-aligned validation/reference data where it materially
  strengthens platform claims
- one broader platform/runtime comparison surface inspired by SPECULA
- updated validity, benchmark, and documentation-map traceability

### Recommended Execution Order

1. `PSP-05`
2. `PSP-06`
3. `PSP-07`

### Milestones

| Milestone | Task | Outputs |
| --- | --- | --- |
| `PSP-05` | Convert the tomography representative benchmark defer into a maintained benchmark artifact, or explicitly re-scope it with benchmark-backed reasoning if it is still not practical. | benchmark artifact or updated defer note |
| `PSP-06` | Expand SPECULA-aligned validation beyond current contract-only cases where the platform claim would benefit from broader evidence. | new frozen bundle cases and updated validity matrix |
| `PSP-07` | Add at least one platform-scale SPECULA-informed scenario comparing orchestration/runtime behavior rather than only isolated optical kernels. | new benchmark/validation contract and archived result |

### `PSP-05`: Representative Tomography Evidence

Status:

- completed

Objective:

- either close the tomography benchmark defer with maintained representative
  evidence or replace the current defer note with a more defensible scoped
  decision backed by fresh measurements

Primary files in scope:

- [docs/tomography-benchmark-scope.md](./tomography-benchmark-scope.md)
- [docs/model-validity-matrix.md](./model-validity-matrix.md)
- [docs/benchmark-matrix-plan.md](./benchmark-matrix-plan.md)
- a maintained generator and archived artifact if evidence is added

Detailed work:

- re-evaluate the current tomography benchmark blockers against the present
  codebase and benchmark surfaces
- if representative evidence is now practical:
  - add a maintained generator script
  - archive the result artifact and manifest
  - wire it into the validity matrix
- if representative evidence is still not practical:
  - replace the current defer with a tighter scope note
  - record concrete blocker categories and revisit conditions

Deliverables:

- benchmark artifact plus manifest, or
- stronger defer note with benchmark-backed rationale

Acceptance criteria:

- `MV-10` is either strengthened with current evidence or re-scoped with a
  materially better justification than the existing note
- the decision is reproducible and recorded in maintained docs

Evidence to record:

- benchmark script path if added
- archived result path if added
- explicit blocker list if still deferred

Implementation record:

- archived scope decision added:
  [2026-04-02-phase2-psp05.toml](../benchmarks/results/tomography/2026-04-02-phase2-psp05.toml)
- manifest added:
  [manifest.toml](../benchmarks/results/tomography/manifest.toml)
- maintained note updated:
  [tomography-benchmark-scope.md](./tomography-benchmark-scope.md)
- `MV-10` remains a scoped defer for representative benchmark evidence, but the
  defer is now anchored to an explicit archived decision record and maintenance
  timeout budget instead of the earlier prose-only note

### `PSP-06`: Broader SPECULA-Aligned Validation

Status:

- completed

Objective:

- strengthen SPECULA-referenced claims beyond today’s mostly contract-oriented
  frozen cases where broader runtime or orchestration evidence is useful

Primary files in scope:

- [docs/specula-reference-datasets.md](./specula-reference-datasets.md)
- [docs/model-validity-matrix.md](./model-validity-matrix.md)
- [test/reference_harness.jl](../test/reference_harness.jl)
- [scripts/generate_specula_reference_bundle.jl](../scripts/generate_specula_reference_bundle.jl)
- [test/reference_data_specula](../test/reference_data_specula)

Detailed work:

- identify one or two high-value SPECULA-aligned scenarios where current
  evidence is too narrow, focusing on:
  - broader runtime behavior,
  - grouped/platform execution,
  - or platform-relevant orchestration rather than more isolated kernels
- extend the SPECULA reference harness/generator only as needed for those cases
- add frozen reference cases and wire them into the maintained regression suite
- update the validity matrix with the stronger scope

Deliverables:

- new SPECULA-aligned reference bundle cases
- updated generator/harness support
- updated validity and provenance docs

Acceptance criteria:

- at least one new SPECULA-aligned case goes beyond current narrow optical
  contract coverage
- the new evidence is committed, reproducible, and exercised in maintained test
  flows

Evidence to record:

- new case IDs and bundle paths
- generator/harness changes
- updated validity status for the affected matrix entry

Implementation record:

- new frozen SPECULA-targeted cases added:
  - `shack_hartmann_polychromatic_frame`
  - `pyramid_polychromatic_frame`
- harness widened in:
  [reference_harness.jl](../test/reference_harness.jl)
- bundle generator updated in:
  [generate_specula_reference_bundle.jl](../scripts/generate_specula_reference_bundle.jl)
- frozen bundle refreshed under:
  [test/reference_data_specula](../test/reference_data_specula)
- validity notes for `MV-05` and `MV-06` now record the narrow SPECULA
  detector-frame contracts

### `PSP-07`: Platform-Scale SPECULA-Informed Scenario

Status:

- completed

Objective:

- add one broader scenario that evaluates orchestration/runtime behavior in a
  way that better reflects platform-scale comparison than isolated component
  validation does

Primary files in scope:

- [benchmarks/contracts/cross_package.toml](../benchmarks/contracts/cross_package.toml)
- [scripts/run_cross_package_benchmarks.jl](../scripts/run_cross_package_benchmarks.jl)
- [docs/cross-package-benchmark-harness.md](./cross-package-benchmark-harness.md)
- [docs/model-validity-matrix.md](./model-validity-matrix.md)
- archived results under [benchmarks/results](../benchmarks/results)

Detailed work:

- define one broader SPECULA-informed scenario contract with:
  - normalized parameters,
  - explicit accepted non-equivalences,
  - and clearly stated success metrics
- prefer a platform/runtime scenario that exercises orchestration or grouped
  behavior rather than a single isolated propagation kernel
- add runner support and archive the resulting baseline, or explicitly record a
  precise environment/scenario gate if an external dependency still prevents the
  full comparison

Deliverables:

- new contract entry
- archived result or scoped environment-gated result note
- updated benchmark harness docs

Acceptance criteria:

- Phase 2 ends with at least one broader SPECULA-informed runtime/platform
  comparison surface recorded in the maintained benchmark/validation set
- environment-gated limitations, if any, are explicit instead of implicit

Evidence to record:

- contract ID
- archived result path
- normalization rules and accepted differences

Implementation record:

- maintained generator added:
  [generate_specula_platform_runtime_artifact.jl](../scripts/generate_specula_platform_runtime_artifact.jl)
- archived artifact added:
  [2026-04-02-phase2-psp07.toml](../benchmarks/results/platform/2026-04-02-phase2-psp07.toml)
- manifest added:
  [manifest.toml](../benchmarks/results/platform/manifest.toml)
- maintained note added:
  [specula-platform-runtime-validation.md](./specula-platform-runtime-validation.md)
- `MV-13` now records a broader SPECULA-informed Julia-native platform/runtime
  artifact rather than only the earlier grouped-runtime artifact

### Acceptance

- `MV-10` is improved or re-scoped with stronger evidence
- SPECULA-aligned evidence covers at least one broader runtime/platform scenario
- new claims are reflected in [model-validity-matrix.md](./model-validity-matrix.md)

## Phase 3: AMDGPU Workaround Reduction

### Goal

Reduce the remaining maintained host-mirror fallback surfaces where that yields
real value.

### Milestones

| Milestone | Task | Outputs |
| --- | --- | --- |
| `PSP-08` | Write and commit a current ROCm fallback inventory after the execution-plan rollout, with each fallback labeled as temporary, acceptable, or next-target-for-removal. | inventory note |
| `PSP-09` | Remove or shrink the highest-value remaining ROCm fallback in a benchmark-backed hot path. | code change + benchmark evidence |
| `PSP-10` | Re-baseline CPU/AMDGPU/CUDA realistic runtime surfaces after the ROCm cleanup. | updated benchmark note/results |

### Acceptance

- at least one maintained ROCm host-mirror fallback is removed or materially reduced
- grouped and AO3k benchmark surfaces remain in-family
- residual fallback behavior is explicit rather than ambient

### Implementation record

- fallback inventory added:
  [rocm-fallback-inventory.md](./rocm-fallback-inventory.md)
- Phase 3 rebaseline note added:
  [rocm-phase3-rebaseline.md](./rocm-phase3-rebaseline.md)
- archived runtime artifact added:
  [2026-04-02-phase3-psp10.toml](../benchmarks/results/platform/2026-04-02-phase3-psp10.toml)

## Phase 4: Julia-Native Platform Orchestration

### Goal

Build a stronger platform-scale orchestration layer without switching to
config-first design.

### Design Rule

The primary interface remains:

- typed scenario/config structs
- Julia builders
- scripts that construct and run those scenarios

This phase should introduce richer orchestration objects, not external config
files.

### Milestones

| Milestone | Task | Outputs |
| --- | --- | --- |
| `PSP-11` | Define typed platform scenario/config families for representative runtime compositions. | scenario/config types and guide |
| `PSP-12` | Add builders for representative multi-branch platform compositions. | builder surfaces for grouped/multi-WFS/multi-detector cases |
| `PSP-13` | Add runtime entry points that make these scenario objects easy to execute, benchmark, and validate from scripts. | maintained runner/API surfaces |
| `PSP-14` | Add one or two canonical main-platform example scripts that use the new orchestration layer. | examples and docs |

### Acceptance

- platform composition is more explicit and reusable than today’s ad hoc script assembly
- scripts remain the main entry point
- no declarative manifest layer is required to use the new platform model

### Implementation record

- typed orchestration layer added in [src/Control/platform.jl](../src/Control/platform.jl)
- new exported scenario/config families:
  - `ClosedLoopBranchConfig`
  - `SinglePlatformConfig`
  - `GroupedPlatformConfig`
  - `PlatformScenario`
- maintained orchestration guide added:
  [platform-orchestration.md](./platform-orchestration.md)
- maintained grouped runtime runner migrated to the new layer:
  [profile_multi_source_multi_wfs_runtime.jl](../scripts/profile_multi_source_multi_wfs_runtime.jl)
- canonical example scripts added:
  - [platform_single_runtime.jl](../examples/closed_loop/platform_single_runtime.jl)
  - [platform_grouped_runtime.jl](../examples/closed_loop/platform_grouped_runtime.jl)

## Phase 5: Platform-Scale Validation And Benchmarking

### Goal

Back the new orchestration layer with realistic evidence.

### Milestones

| Milestone | Task | Outputs |
| --- | --- | --- |
| `PSP-15` | Add maintained realistic benchmark surfaces for the new scenario/composition layer on CPU and maintained GPUs. | benchmark runners/results |
| `PSP-16` | Add one cross-package platform scenario contract where normalization is strong enough to make comparison meaningful. | contract doc + archived baseline |
| `PSP-17` | Add backend parity/functional coverage for the new orchestration surfaces. | tests/smoke additions |

### Acceptance

- new orchestration surfaces have maintained benchmark evidence
- at least one realistic cross-package scenario is normalized and archived
- CPU/AMDGPU/CUDA claims are explicit and tested where maintained

Implementation record:

- direct platform runner added:
  [profile_platform_runtime.jl](../scripts/profile_platform_runtime.jl)
- archived direct orchestration artifact added:
  [2026-04-03-phase5-psp15.toml](../benchmarks/results/platform/2026-04-03-phase5-psp15.toml)
- orchestration validation note added:
  [platform-orchestration-validation.md](./platform-orchestration-validation.md)
- normalized cross-package platform contract added:
  [cross_package.toml](../benchmarks/contracts/cross_package.toml)
- archived cross-package result added:
  [2026-04-03-phase5-psp16.toml](../benchmarks/results/cross_package/2026-04-03-phase5-psp16.toml)
- contract note added:
  [revolt-platform-benchmark-contract.md](./revolt-platform-benchmark-contract.md)
- backend smoke widened to instantiate `PlatformScenario` directly in:
  [backend_optional_common.jl](../test/backend_optional_common.jl)

## Phase 6: Closeout And Deferred Tracks

### Goal

Finish the phase with an explicit decision record rather than drifting into the
next topic.

### Milestones

| Milestone | Task | Outputs |
| --- | --- | --- |
| `PSP-18` | Record what remains behind SPECULA after the orchestration pass. | decision/closeout note |
| `PSP-19` | Explicitly defer or schedule science-path / focal-plane work under the optional boundary rule. | roadmap/boundary update |
| `PSP-20` | Explicitly defer or schedule a future optional manifest/config layer, if still desired, as a secondary interface on top of the typed platform model. | defer note or future-plan link |

### Acceptance

- the next platform step is explicit
- science-path remains optional unless deliberately changed
- config manifests remain deferred and non-primary unless a later plan changes that rule

Implementation record:

- explicit closeout note added:
  [platform-strengthening-closeout.md](./platform-strengthening-closeout.md)
- science-path boundary reaffirmed in:
  [optional-integration-boundaries.md](./optional-integration-boundaries.md)
- explicit manifest/config defer note added:
  [platform-manifest-defer.md](./platform-manifest-defer.md)
- roadmap and future-direction docs updated to point at the closeout decision

## Execution Order

1. `PSP-01`
2. `PSP-02`
3. `PSP-03`
4. `PSP-04`
5. `PSP-05`
6. `PSP-06`
7. `PSP-07`
8. `PSP-08`
9. `PSP-09`
10. `PSP-10`
11. `PSP-11`
12. `PSP-12`
13. `PSP-13`
14. `PSP-14`
15. `PSP-15`
16. `PSP-16`
17. `PSP-17`
18. `PSP-18`
19. `PSP-19`
20. `PSP-20`

## Traceability Matrix

| Requirement | Covered by milestones | Completion rule |
| --- | --- | --- |
| `PSR-01` | `PSP-02`, `PSP-03`, `PSP-04`, `PSP-11` to `PSP-14` | scripts and typed builders remain the primary documented workflow |
| `PSR-02` | `PSP-04`, `PSP-20` | config-first design remains explicitly deferred |
| `PSR-03` | `PSP-01` to `PSP-10` | the remaining review-quality gaps are materially reduced |
| `PSR-04` | `PSP-05` to `PSP-07`, `PSP-15` to `PSP-17` | validation and benchmark evidence are stronger where current claims are still narrow |
| `PSR-05` | `PSP-08` to `PSP-10` | ROCm workaround debt is reduced or explicitly justified |
| `PSR-06` | `PSP-06`, `PSP-07`, `PSP-16`, `PSP-18` | SPECULA-informed platform work is evidence-backed and Julia-native |
| `PSR-07` | `PSP-19` | science-path work remains outside core unless deliberately revisited |

## Checklist

### Phase 1

- [x] `PSP-01`
- [x] `PSP-02`
- [x] `PSP-03`
- [x] `PSP-04`

### Phase 2

- [x] `PSP-05`
- [x] `PSP-06`
- [x] `PSP-07`

### Phase 3

- [x] `PSP-08`
- [x] `PSP-09`
- [x] `PSP-10`

### Phase 4

- [x] `PSP-11`
- [x] `PSP-12`
- [x] `PSP-13`
- [x] `PSP-14`

### Phase 5

- [x] `PSP-15`
- [x] `PSP-16`
- [x] `PSP-17`

### Phase 6

- [x] `PSP-18`
- [x] `PSP-19`
- [x] `PSP-20`

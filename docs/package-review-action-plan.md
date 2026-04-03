# Package Review Action Plan

Date: 2026-03-31

This document converts the recommendations in
[package-review-2026-03.md](./package-review-2026-03.md) into a phased,
checkable implementation plan with explicit traceability.

## Purpose

This plan exists to ensure that the recommendations from the package review are:

- captured as stable action items
- assigned to concrete phases
- tied to code/docs/test deliverables
- verifiable with evidence
- easy to mark complete without re-reading the whole review

## Traceability Rules

- Every actionable recommendation from
  [package-review-2026-03.md](./package-review-2026-03.md) is assigned a stable
  review action ID: `PR-XX`.
- Every phase task is assigned a stable implementation task ID: `PLAN-XX`.
- Each `PLAN-XX` item must reference one or more `PR-XX` items.
- A `PR-XX` item is not complete until:
  - its implementation is present,
  - its verification evidence is recorded,
  - its status is updated in the traceability matrix below.

Status values:

- `[ ]` not started
- `[~]` in progress / partially complete
- `[x]` complete
- `[!]` intentionally deferred

## Source Recommendation Register

The following IDs are the normalized actionable recommendations extracted from
[package-review-2026-03.md](./package-review-2026-03.md).

### API and interface

- `PR-01`: Curate the public API into explicit tiers.
- `PR-02`: Reduce the exported top-level surface aggressively.
- `PR-03`: Keep generic verb-based APIs, but improve discoverability through
  API layering rather than renaming.
- `PR-04`: Add a public API policy to the docs and enforce it in exports.

### Internal modularization

- `PR-05`: Split oversized subsystem files, especially SH, Pyramid, BioEdge,
  and runtime.
- `PR-06`: Split the monolithic test driver into subsystem-focused files.
- `PR-07`: Prefer further internal modularization before splitting the package
  into multiple packages.
- `PR-08`: Keep optional integrations out of core or in extensions/separate
  adapters.

### Reusable infrastructure

- `PR-09`: Extract detector noise/readout pipeline helpers into reusable
  infrastructure.
- `PR-10`: Extract grouped WFS execution skeletons shared across SH/Pyramid/
  BioEdge.
- `PR-11`: Extract runtime product planning and export helpers.
- `PR-12`: Extract shared source/field/atmosphere sampling and accumulation
  helpers.
- `PR-13`: Move backend reduction/random-fill services further out of model
  code and into shared backend services.
- `PR-14`: Extract common calibration workflow scaffolding where orchestration
  is duplicated.

### Validation and correctness

- `PR-15`: Build a formal validation matrix for maintained model families.
- `PR-16`: Distinguish analytic checks, reference-bundle checks, backend parity
  checks, and realistic benchmark evidence.
- `PR-17`: Add a maintained model validity matrix doc linking model claims to
  evidence and limitations.
- `PR-18`: Keep OOPAO as the primary parity baseline.
- `PR-19`: Add targeted SPECULA comparison baselines where SPECULA is the
  stronger reference.
- `PR-20`: Prefer frozen reference bundles over live CI dependence on external
  repositories.

### Benchmarking

- `PR-21`: Provide maintained cross-package benchmarks against equivalent
  OOPAO/SPECULA models.
- `PR-22`: Treat cross-package benchmarks as benchmark/validation artifacts,
  not unit tests.
- `PR-23`: Use REVOLT-like scenarios as the main benchmark harness for
  realistic comparisons.
- `PR-24`: Measure both fidelity and runtime, not just speed.
- `PR-25`: Maintain compact, medium, and representative comparison classes.

### Documentation and clarity

- `PR-26`: Consolidate plan-heavy docs into a smaller set of stable guides.
- `PR-27`: Mark docs clearly as active, completed, or archived/superseded.
- `PR-28`: Create a stronger maintainer architecture guide for the current
  system.
- `PR-29`: Expand the user-facing synthesis docs rather than only adding more
  plans.
- `PR-30`: Add an end-to-end runtime/dataflow document.

### Backend validation and maintainability

- `PR-31`: Make backend validation more first-class in the test layout.
- `PR-32`: Add clearer optional CUDA coverage in the test tree.
- `PR-33`: Separate backend smoke, functional tests, and benchmark guards more
  cleanly.
- `PR-34`: Continue isolating backend-specific execution policy out of model
  code.

### Future feature and platform direction

- `PR-35`: Avoid chasing additional OOPAO parity for its own sake.
- `PR-36`: Use SPECULA as the stronger reference for next platform-level
  breadth.
- `PR-37`: Revisit controller/process breadth later as a platform-level
  expansion area.
- `PR-38`: Keep science-path integrations optional and outside core where
  appropriate.

## Phase Plan

## Phase 0: Governance, Inventory, and Baselines

Goal:

- establish the working governance artifacts before larger structural changes

Exit criteria:

- all five inventory/governance artifacts exist
- each artifact has a maintainer-facing purpose and status field
- each artifact links back to one or more `PR-XX` IDs
- at least one baseline metric is recorded for:
  - export surface size
  - oversized file inventory
  - validation evidence inventory
  - benchmark harness inventory

Phase outputs:

- API tier inventory doc:
  [`docs/api-tier-inventory.md`](./api-tier-inventory.md)
- modularization inventory doc:
  [`docs/modularization-inventory.md`](./modularization-inventory.md)
- reusable-infrastructure inventory doc:
  [`docs/reusable-infrastructure-inventory.md`](./reusable-infrastructure-inventory.md)
- model-validation inventory doc:
  [`docs/model-validation-inventory.md`](./model-validation-inventory.md)
- cross-package benchmark inventory doc:
  [`docs/cross-package-benchmark-inventory.md`](./cross-package-benchmark-inventory.md)
- a short index section in this plan linking to all five outputs

Recommended implementation order:

1. `PLAN-01`
2. `PLAN-02`
3. `PLAN-04`
4. `PLAN-05`
5. `PLAN-03`

Tasks:

- `[x]` `PLAN-01` Create an API tier inventory:
  - stable public API
  - advanced/expert API
  - internal extension hooks
  - traces to: `PR-01`, `PR-02`, `PR-04`
- `[x]` `PLAN-02` Create a source-file modularization inventory for oversized
  files and proposed splits
  - traces to: `PR-05`, `PR-06`, `PR-07`
- `[x]` `PLAN-03` Create a reusable-infrastructure duplication inventory
  covering detector, WFS, runtime, backend, and calibration scaffolding
  - traces to: `PR-09` through `PR-14`
- `[x]` `PLAN-04` Create a model-validation inventory listing all maintained
  model families and current evidence
  - traces to: `PR-15`, `PR-16`, `PR-17`
- `[x]` `PLAN-05` Create a cross-package benchmark inventory from:
  - `../REVOLT`
  - `../AdaptiveOpticsComparisons`
  - existing local benchmark scripts
  - traces to: `PR-21` through `PR-25`

Primary files:

- [`docs/api-reference.md`](./api-reference.md)
- [`docs/user-guide.md`](./user-guide.md)
- [`docs/benchmark-matrix-plan.md`](./benchmark-matrix-plan.md)
- new inventory docs as needed under [`docs/`](./)

Verification:

- inventory docs exist
- each inventory item links to source files and current status

### Phase 0 Task Detail

#### `PLAN-01`: API tier inventory

Status:

- complete

Artifact:

- [`docs/api-tier-inventory.md`](./api-tier-inventory.md)

Objective:

- produce one authoritative inventory of the current exported and non-exported
  API surfaces, classified by intended audience and stability

Detailed work:

- `[x]` count and record the current exported top-level surface from
  [`src/AdaptiveOpticsSim.jl`](../src/AdaptiveOpticsSim.jl)
- `[x]` classify exports into:
  - end-user stable workflow API
  - advanced/expert simulation API
  - backend/infrastructure hooks that should likely stop being exported
- `[x]` identify names that should remain exported for workflow ergonomics even
  if they are not mathematically fundamental
- `[x]` identify names that are currently undocumented but effectively public
- `[x]` create a proposed “target export surface” column for later Phase 1 use

Deliverable artifact:

- new doc, recommended path:
  [`docs/api-tier-inventory.md`](./api-tier-inventory.md)

Minimum schema:

- symbol
- current export status
- proposed tier
- rationale
- source file
- review trace IDs

Acceptance criteria:

- current export count is recorded
- every exported symbol is mapped to a tier
- obvious non-user-facing exports are flagged explicitly

Evidence to record:

- export count before Phase 1
- link to generated/maintained inventory doc

#### `PLAN-02`: modularization inventory

Status:

- complete

Artifact:

- [`docs/modularization-inventory.md`](./modularization-inventory.md)

Objective:

- identify the oversized files and define proposed split boundaries before
  touching implementation

Detailed work:

- `[x]` record file size baselines for the largest source and test files
- `[x]` identify responsibility clusters inside:
  - [`src/WFS/shack_hartmann.jl`](../src/WFS/shack_hartmann.jl)
  - [`src/WFS/pyramid.jl`](../src/WFS/pyramid.jl)
  - [`src/WFS/bioedge.jl`](../src/WFS/bioedge.jl)
  - [`src/Control/runtime.jl`](../src/Control/runtime.jl)
  - [`test/runtests.jl`](../test/runtests.jl)
- `[x]` propose concrete target file/module boundaries for each
- `[x]` identify high-risk split points where cyclic dependencies or implicit
  invariants may exist
- `[x]` identify files that should *not* be split yet

Deliverable artifact:

- new doc, recommended path:
  [`docs/modularization-inventory.md`](./modularization-inventory.md)

Minimum schema:

- current file
- current line count
- responsibility clusters
- proposed split files
- dependency risks
- recommended phase order
- review trace IDs

Acceptance criteria:

- all phase-targeted large files are inventoried
- each has a proposed split sketch
- split candidates are ranked by expected value and risk

Evidence to record:

- baseline line counts
- proposed target layout

#### `PLAN-03`: reusable-infrastructure duplication inventory

Status:

- complete

Artifact:

- [`docs/reusable-infrastructure-inventory.md`](./reusable-infrastructure-inventory.md)

Objective:

- identify repeated orchestration and backend-policy patterns that should be
  turned into shared blocks

Detailed work:

- `[x]` inspect detector noise/readout flows for repeated orchestration
- `[x]` inspect SH/Pyramid/BioEdge grouped execution paths for repeated
  grouping, buffering, normalization, and export behavior
- `[x]` inspect runtime export/product planning duplication
- `[x]` inspect source/field/atmosphere propagation accumulation patterns
- `[x]` inspect backend-specific reduction/random-fill policy duplication
- `[x]` inspect calibration workflow orchestration duplication

Deliverable artifact:

- new doc, recommended path:
  [`docs/reusable-infrastructure-inventory.md`](./reusable-infrastructure-inventory.md)

Minimum schema:

- duplicated pattern
- occurrences
- proposed shared abstraction
- likely owning module
- risk of abstraction
- expected payoff
- review trace IDs

Acceptance criteria:

- at least one concrete reusable-block candidate exists for each of:
  - detector pipeline
  - grouped WFS execution
  - runtime product/export
  - backend services
- each candidate has an ownership proposal

Evidence to record:

- code references for each repeated pattern

#### `PLAN-04`: model-validation inventory

Status:

- complete

Artifact:

- [`docs/model-validation-inventory.md`](./model-validation-inventory.md)

Objective:

- create a current-state map from maintained model families to validation
  evidence

Detailed work:

- `[x]` enumerate maintained model families:
  - atmosphere
  - field/propagation
  - detectors
  - WFS families
  - runtime/control
  - tomography/calibration
- `[x]` map each family to current evidence classes:
  - analytic checks
  - reference bundle checks
  - backend parity checks
  - benchmark evidence
- `[x]` mark evidence gaps explicitly
- `[x]` identify where OOPAO is the current baseline
- `[x]` identify where SPECULA should become the targeted baseline

Deliverable artifact:

- new doc, recommended path:
  [`docs/model-validation-inventory.md`](./model-validation-inventory.md)

Acceptance criteria:

- all major maintained model families are listed
- each family has current evidence or an explicit gap
- evidence classes are normalized for use in Phase 4

Evidence to record:

- references to test files, scripts, and bundles

#### `PLAN-05`: cross-package benchmark inventory

Status:

- complete

Artifact:

- [`docs/cross-package-benchmark-inventory.md`](./cross-package-benchmark-inventory.md)

Objective:

- inventory the available scenario sources and benchmarking assets for realistic
  cross-package comparisons

Detailed work:

- `[x]` inspect relevant benchmark/runtime assets in:
  - [`../REVOLT`](../REVOLT)
  - [`../AdaptiveOpticsComparisons`](../../AdaptiveOpticsComparisons)
  - local `scripts/profile_*`
- `[x]` identify scenario families that exist in at least two of:
  - AdaptiveOpticsSim
  - OOPAO / REVOLT-OOPAO
  - SPECULA / REVOLT-SPECULA
- `[x]` identify configuration compatibility constraints:
  - atmosphere
  - detector family
  - WFS type
  - controller/reconstructor assumptions
- `[x]` define which scenarios are good compact, medium, and representative
  seeds

Deliverable artifact:

- new doc, recommended path:
  [`docs/cross-package-benchmark-inventory.md`](./cross-package-benchmark-inventory.md)

Acceptance criteria:

- at least one viable compact, medium, and representative candidate scenario
  is identified
- the available source trees and scripts are linked explicitly
- known comparability limits are recorded

Evidence to record:

- script/file references for each candidate scenario

### Phase 0 Completion Checklist

- `[x]` all five inventory docs created
- `[x]` all five docs linked from this plan
- `[x]` baseline export count recorded
- `[x]` baseline oversized-file counts recorded
- `[x]` validation evidence inventory recorded
- `[x]` benchmark-scenario inventory recorded

## Phase 1: Public API Curation

Goal:

- reduce surface-area complexity without breaking core workflow usability

Dependencies:

- requires `PLAN-01` completion
- should use `PLAN-02` output for scenario-builder triage

Exit criteria:

- a stable public API tier is defined and documented
- advanced/expert APIs are documented separately from the stable tier
- selected infrastructure/backend/helper names are de-exported
- user-guide and API docs reflect the new policy
- workflow examples and tests still pass

Phase outputs:

- updated export surface in [`src/AdaptiveOpticsSim.jl`](../src/AdaptiveOpticsSim.jl)
- public API policy text in docs
- stable-vs-advanced API presentation in docs
- a recorded before/after export count

Recommended implementation order:

1. `PLAN-09`
2. `PLAN-06`
3. `PLAN-07`
4. `PLAN-10`
5. `PLAN-08`

Tasks:

- `[x]` `PLAN-06` Define the stable top-level public API set
  - traces to: `PR-01`, `PR-02`
- `[x]` `PLAN-07` Define the advanced/expert non-exported but documented API
  tier
  - traces to: `PR-01`, `PR-03`
- `[x]` `PLAN-08` Remove or de-export infrastructure/backend/helper names that
  do not belong in the default top-level API
  - traces to: `PR-02`, `PR-04`
- `[x]` `PLAN-09` Add an explicit public API policy section to the docs
  - traces to: `PR-04`
- `[x]` `PLAN-10` Audit scenario-builder and convenience exports and move
  non-core ones to namespaced or example-support access
  - traces to: `PR-02`, `PR-08`

Primary files:

- [`src/AdaptiveOpticsSim.jl`](../src/AdaptiveOpticsSim.jl)
- [`docs/api-reference.md`](./api-reference.md)
- [`docs/user-guide.md`](./user-guide.md)

Verification:

- exported surface count is reduced and recorded
- docs clearly separate stable vs advanced APIs
- common user workflows remain intact in examples/tests

### Phase 1 Task Detail

#### `PLAN-06`: define the stable top-level public API set

Objective:

- choose the exported surface that most users should learn first and rely on

Detailed work:

- `[x]` derive the stable public API candidate list from
  [`docs/api-tier-inventory.md`](./api-tier-inventory.md)
- `[x]` group the stable API by user workflow:
  - optical setup
  - atmosphere
  - detector/WFS
  - calibration/reconstruction
  - closed-loop runtime
- `[x]` identify symbols that are stable but should remain namespaced rather
  than exported
- `[x]` define an explicit “must remain exported” set for workflow ergonomics

Acceptance criteria:

- the stable public API set is finite, named, and reviewable
- it is small enough to present clearly in the user guide and API reference

Evidence to record:

- stable API symbol list
- count of stable exported names
- stable exported count after curation: `533`

#### `PLAN-07`: define the advanced/expert API tier

Objective:

- preserve expert access without forcing it into the beginner-facing surface

Detailed work:

- `[x]` identify expert simulation, calibration, backend, and debugging hooks
- `[x]` separate:
  - advanced public but non-default APIs
  - internal extension hooks
- `[x]` define documentation rules for advanced APIs:
  - document them
  - do not necessarily export them
  - keep them stable only where justified

Acceptance criteria:

- advanced/expert API is explicitly documented
- it is not mixed ambiguously with stable beginner-facing API documentation

Evidence to record:

- advanced API list
- documentation location

#### `PLAN-08`: remove or de-export non-core top-level names

Objective:

- reduce the default namespace burden without breaking common workflows

Detailed work:

- `[x]` remove exports for infrastructure/backend/helper names that do not
  belong in the default workflow surface
- `[x]` preserve access through `AdaptiveOpticsSim.<name>`
- `[x]` update examples/tests/docs for any symbols that stop being exported
- `[x]` record de-exported symbols in the API policy changelog section

Acceptance criteria:

- top-level export count decreases materially
- common imports and examples remain ergonomic
- no expert symbol becomes inaccessible

Evidence to record:

- before/after export counts
- de-export list
- export count before/after: `566 -> 533`

#### `PLAN-09`: add the public API policy to docs

Objective:

- make the export and stability rules explicit so later work follows policy

Detailed work:

- `[x]` add a “Public API Policy” section to
  [`docs/api-reference.md`](./api-reference.md)
- `[x]` add a short “How to read the API surface” section to
  [`docs/user-guide.md`](./user-guide.md)
- `[x]` define:
  - stable public API
  - advanced API
  - internal hooks
- `[x]` define export criteria for future additions

Acceptance criteria:

- policy is written once and referenced elsewhere instead of re-explained ad
  hoc
- contributors can classify new exports using the policy

Evidence to record:

- doc references to the policy section
- [`docs/api-reference.md`](./api-reference.md)
- [`docs/user-guide.md`](./user-guide.md)

#### `PLAN-10`: audit scenario-builder and convenience exports

Objective:

- remove convenience/scenario assembly from the default public surface where it
  does not belong

Detailed work:

- `[x]` identify scenario-builder helpers currently exported
- `[x]` classify them as:
  - keep exported
  - keep public but namespaced
  - move to example/support-oriented surface later
- `[x]` check tutorial/example usage to avoid breaking intended workflows
- `[x]` update docs to point users to examples when the helper is more of a
  scenario builder than a core operator

Acceptance criteria:

- scenario-builder exports are intentionally classified
- the stable public API becomes more mathematically and operationally coherent

Evidence to record:

- classification table for scenario-builder exports
- [`docs/api-tier-inventory.md`](./api-tier-inventory.md)

### Phase 1 Completion Checklist

- `[x]` API policy added to docs
- `[x]` stable API list recorded
- `[x]` advanced/expert API list recorded
- `[x]` de-export list recorded
- `[x]` before/after export counts recorded
- `[x]` user guide and API reference updated
- `[x]` examples/tests still pass after curation

Phase 1 evidence summary:

- full suite passed with `julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'`
- export count reduced from `566` to `533`
- stable, advanced, and de-export candidate tiers recorded in
  [`docs/api-tier-inventory.md`](./api-tier-inventory.md)

## Phase 2: Internal Modularization

Goal:

- split oversized files into maintainable internal modules without changing the
  package boundary unnecessarily

Tasks:

- `[x]` `PLAN-11` Split [`src/WFS/shack_hartmann.jl`](../src/WFS/shack_hartmann.jl)
  into smaller units
  - traces to: `PR-05`
- `[x]` `PLAN-12` Split [`src/WFS/pyramid.jl`](../src/WFS/pyramid.jl)
  - traces to: `PR-05`
- `[x]` `PLAN-13` Split [`src/WFS/bioedge.jl`](../src/WFS/bioedge.jl)
  - traces to: `PR-05`
- `[x]` `PLAN-14` Split [`src/Control/runtime.jl`](../src/Control/runtime.jl)
  by product planning, latency staging, execution, and exports
  - traces to: `PR-05`, `PR-11`
- `[x]` `PLAN-15` Split [`test/runtests.jl`](../test/runtests.jl) into:
  - unit/subsystem tests
  - interface/conformance tests
  - reference tests
  - tutorial/example tests
  - traces to: `PR-06`, `PR-31`, `PR-33`
- `[x]` `PLAN-16` Re-evaluate whether any optional integrations should move to
  extensions or separate adapter packages
  - traces to: `PR-07`, `PR-08`, `PR-38`

Verification:

- large-file line counts are substantially reduced
- test entry points remain green
- no loss of subsystem coverage

Phase 2 evidence summary:

- large entry-point files were reduced to thin include surfaces:
  - [`src/WFS/shack_hartmann.jl`](../src/WFS/shack_hartmann.jl): `2167 -> 12`
  - [`src/WFS/pyramid.jl`](../src/WFS/pyramid.jl): `1463 -> 11`
  - [`src/WFS/bioedge.jl`](../src/WFS/bioedge.jl): `1241 -> 10`
  - [`src/Control/runtime.jl`](../src/Control/runtime.jl): `1123 -> 12`
  - [`test/runtests.jl`](../test/runtests.jl): `3300 -> 8`
- split layout and boundary notes recorded in
  [`docs/modularization-inventory.md`](./modularization-inventory.md)
- optional integration/package-boundary decision recorded in
  [`docs/optional-integration-boundaries.md`](./optional-integration-boundaries.md)
- optional backend validation is now first-class in the test tree through:
  - [`test/optional_amdgpu_backends.jl`](../test/optional_amdgpu_backends.jl)
  - [`test/optional_cuda_backends.jl`](../test/optional_cuda_backends.jl)
  - shared scaffolding in
    [`test/backend_optional_common.jl`](../test/backend_optional_common.jl)
- functional tests, optional backend smoke, and benchmark/profile evidence are
  now documented separately in
  [`docs/backend-validation-guide.md`](./backend-validation-guide.md)
- full suite passed with `julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'`

## Phase 3: Reusable Runtime and Backend Building Blocks

Goal:

- reduce duplication and remove repeated model-local infrastructure patterns

Tasks:

- `[x]` `PLAN-17` Extract detector noise/readout pipeline helpers into shared
  infrastructure
  - traces to: `PR-09`
- `[x]` `PLAN-18` Extract grouped WFS execution scaffolding shared across SH,
  Pyramid, and BioEdge
  - traces to: `PR-10`
- `[x]` `PLAN-19` Extract runtime product planning and export helpers from the
  runtime/WFS files
  - traces to: `PR-11`
- `[x]` `PLAN-20` Extract shared source/field/atmosphere accumulation helpers
  into one maintained layer
  - traces to: `PR-12`
- `[x]` `PLAN-21` Continue moving backend reduction/random-fill behavior into
  backend services or extensions
  - traces to: `PR-13`, `PR-34`
- `[x]` `PLAN-22` Extract common calibration workflow scaffolding
  - traces to: `PR-14`

Primary files:

- [`src/Core/`](../src/Core)
- [`src/Detectors/`](../src/Detectors)
- [`src/WFS/`](../src/WFS)
- [`src/Control/`](../src/Control)
- [`ext/`](../ext)

Verification:

- duplicated infrastructure code paths are reduced
- backend policy becomes more centralized
- measured warmed benchmark regressions are absent or justified

Evidence recorded:

- shared detector pipeline helpers in
  [`src/Detectors/pipeline.jl`](../src/Detectors/pipeline.jl)
- shared grouped WFS helpers in [`src/WFS/grouped.jl`](../src/WFS/grouped.jl)
- shared calibration helpers in
  [`src/WFS/calibration.jl`](../src/WFS/calibration.jl)
- shared runtime product planning in
  [`src/Control/products.jl`](../src/Control/products.jl)
- shared propagation context helpers in
  [`src/Optics/propagation_context.jl`](../src/Optics/propagation_context.jl)
- shared random/noise services in
  [`src/Core/random_services.jl`](../src/Core/random_services.jl)
- verification run:
  `julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'`

## Phase 4: Validation and Model Correctness Matrix

Goal:

- convert correctness evidence into a maintained traceable validation system

Tasks:

- `[x]` `PLAN-23` Add a model validity matrix doc for all major model families
  - traces to: `PR-15`, `PR-17`
- `[x]` `PLAN-24` Separate validation classes explicitly:
  - analytic
  - reference bundle
  - backend parity
  - benchmark evidence
  - traces to: `PR-16`
- `[x]` `PLAN-25` Expand and formalize OOPAO comparison coverage as the primary
  parity baseline
  - traces to: `PR-18`, `PR-20`
- `[x]` `PLAN-26` Add targeted SPECULA comparison bundles for selected model
  families where SPECULA is the stronger baseline
  - traces to: `PR-19`, `PR-20`, `PR-36`
- `[x]` `PLAN-27` Record model assumptions, limits, and known non-equivalences
  beside the evidence matrix
  - traces to: `PR-15`, `PR-17`

Primary files:

- new validation matrix doc under [`docs/`](./)
- [`test/reference_harness.jl`](../test/reference_harness.jl)
- [`test/reference_data/`](../test/reference_data)

Verification:

- each maintained model family has explicit evidence links
- OOPAO/SPECULA comparison bundles are versioned and reproducible

Evidence recorded:

- maintained validation matrix in
  [`model-validity-matrix.md`](./model-validity-matrix.md)
- OOPAO bundle policy and provenance in
  [`oopao-reference-datasets.md`](./oopao-reference-datasets.md)
- SPECULA-targeted contract bundle policy and provenance in
  [`specula-reference-datasets.md`](./specula-reference-datasets.md)
- committed SPECULA-targeted bundle in
  [`test/reference_data_specula`](../test/reference_data_specula)
- maintained generator in
  [`scripts/generate_specula_reference_bundle.jl`](../scripts/generate_specula_reference_bundle.jl)
- verification run:
  `julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'`

## Phase 5: Cross-Package Benchmark Harness

Goal:

- produce maintained, realistic, comparable benchmark evidence across
  AdaptiveOpticsSim, OOPAO, and SPECULA

Tasks:

- `[x]` `PLAN-28` Define benchmark scenario contracts for compact, medium, and
  representative cross-package comparisons
  - traces to: `PR-23`, `PR-25`
- `[x]` `PLAN-29` Build a REVOLT-like benchmark harness using:
  - `../REVOLT`
  - `../AdaptiveOpticsComparisons`
  - local benchmark scripts
  - traces to: `PR-21`, `PR-23`
- `[x]` `PLAN-30` Record runtime metrics:
  - build time
  - step time
  - frame rate
  - allocations / memory
  - traces to: `PR-21`, `PR-24`
- `[x]` `PLAN-31` Record fidelity metrics:
  - residuals
  - output equivalence norms
  - science / sensing quality metrics where applicable
  - traces to: `PR-24`
- `[x]` `PLAN-32` Keep cross-package benchmarking outside normal `Pkg.test()`
  and document it as maintained engineering evidence
  - traces to: `PR-22`

Primary files:

- [`docs/benchmark-matrix-plan.md`](./benchmark-matrix-plan.md)
- `scripts/profile_*`
- new benchmark harness docs/scripts as needed

Verification:

- benchmark contract docs exist
- at least one maintained REVOLT-like scenario is runnable and archived
- results are versioned by scenario and backend
- maintained artifacts:
  - [`../benchmarks/contracts/cross_package.toml`](../benchmarks/contracts/cross_package.toml)
  - [`./cross-package-benchmark-harness.md`](./cross-package-benchmark-harness.md)
  - [`../scripts/run_cross_package_benchmarks.jl`](../scripts/run_cross_package_benchmarks.jl)
  - [`../benchmarks/results/cross_package/manifest.toml`](../benchmarks/results/cross_package/manifest.toml)
  - [`../benchmarks/results/cross_package/2026-03-31-phase5-baseline.toml`](../benchmarks/results/cross_package/2026-03-31-phase5-baseline.toml)

## Phase 6: Documentation Consolidation and Learning Surfaces

Goal:

- make the package easier to learn and maintain without losing engineering
  detail

Tasks:

- `[x]` `PLAN-33` Mark docs as active, completed, or archived/superseded
  - traces to: `PR-26`, `PR-27`
- `[x]` `PLAN-34` Create a maintainer architecture guide that reflects the
  current implemented system
  - traces to: `PR-28`
- `[x]` `PLAN-35` Expand the user guide into a clearer workflow-oriented
  learning surface
  - traces to: `PR-29`
- `[x]` `PLAN-36` Add an end-to-end runtime/dataflow guide
  - traces to: `PR-30`
- `[x]` `PLAN-37` Reduce plan-doc sprawl by linking detailed plans from a small
  number of maintained index pages
  - traces to: `PR-26`, `PR-29`

Primary files:

- [`docs/user-guide.md`](./user-guide.md)
- [`docs/api-reference.md`](./api-reference.md)
- new architecture/dataflow docs under [`docs/`](./)

Verification:

- docs are status-labeled
- a new contributor can find:
  - stable API
  - architecture overview
  - runtime dataflow
  - validation/benchmark evidence
- maintained Phase 6 outputs:
  - [`documentation-map.md`](./documentation-map.md)
  - [`maintainer-architecture.md`](./maintainer-architecture.md)
  - [`runtime-dataflow.md`](./runtime-dataflow.md)
  - updated [`user-guide.md`](./user-guide.md)
  - updated [`api-reference.md`](./api-reference.md)

## Phase 7: Future Platform Direction

Goal:

- capture the “after cleanup” direction explicitly so future expansion follows
  the review instead of drifting back toward accidental parity chasing

Tasks:

- `[x]` `PLAN-38` Document that OOPAO is a reference, not a ceiling
  - traces to: `PR-35`
- `[x]` `PLAN-39` Document SPECULA-targeted future breadth areas
  - traces to: `PR-36`, `PR-37`
- `[x]` `PLAN-40` Keep science-path integrations optional and outside the core
  package boundary where appropriate
  - traces to: `PR-38`

Verification:

- future roadmap docs reference this plan and do not bypass it
- maintained Phase 7 outputs:
  - [`future-platform-direction.md`](./future-platform-direction.md)
  - updated [`roadmap.md`](./roadmap.md)
  - updated [`algorithmic-implementation-roadmap.md`](./algorithmic-implementation-roadmap.md)
  - updated [`optional-integration-boundaries.md`](./optional-integration-boundaries.md)

## Traceability Matrix

| Review ID | Recommendation | Plan Tasks | Evidence Target | Status |
| --- | --- | --- | --- | --- |
| PR-01 | Curate public API into tiers | PLAN-01, PLAN-06, PLAN-07 | API policy and tiered docs | [x] |
| PR-02 | Reduce exported surface | PLAN-06, PLAN-08, PLAN-10 | reduced export count | [x] |
| PR-03 | Keep generic verbs, improve discoverability via layering | PLAN-07, PLAN-09 | API docs and examples | [x] |
| PR-04 | Add public API policy | PLAN-01, PLAN-09 | documented API policy | [x] |
| PR-05 | Split oversized subsystem files | PLAN-02, PLAN-11, PLAN-12, PLAN-13, PLAN-14 | reduced file sizes / cleaner layout | [x] |
| PR-06 | Split monolithic test driver | PLAN-02, PLAN-15 | test layout restructure | [x] |
| PR-07 | Prefer internal modularization before package splitting | PLAN-02, PLAN-16 | modularization inventory and package-boundary notes | [x] |
| PR-08 | Keep optional integrations out of core | PLAN-10, PLAN-16, PLAN-40 | extension/package boundary docs | [x] |
| PR-09 | Extract detector noise/readout helpers | PLAN-03, PLAN-17 | shared detector infrastructure | [x] |
| PR-10 | Extract grouped WFS execution skeletons | PLAN-03, PLAN-18 | shared WFS grouping helpers | [x] |
| PR-11 | Extract runtime product planning/export helpers | PLAN-03, PLAN-14, PLAN-19 | shared runtime export layer | [x] |
| PR-12 | Extract shared source/field/atmosphere accumulation helpers | PLAN-03, PLAN-20 | shared propagation helpers | [x] |
| PR-13 | Centralize backend reductions/random-fill services | PLAN-03, PLAN-21 | backend-service isolation | [x] |
| PR-14 | Extract calibration scaffolding | PLAN-03, PLAN-22 | shared calibration orchestration | [x] |
| PR-15 | Build formal validation matrix | PLAN-04, PLAN-23, PLAN-27 | model validity matrix | [x] |
| PR-16 | Distinguish validation classes | PLAN-04, PLAN-24 | validation taxonomy doc | [x] |
| PR-17 | Add model validity matrix doc | PLAN-04, PLAN-23, PLAN-27 | validity doc with evidence links | [x] |
| PR-18 | Keep OOPAO as primary parity baseline | PLAN-25 | OOPAO bundle coverage | [x] |
| PR-19 | Add targeted SPECULA baselines | PLAN-26 | SPECULA bundle coverage | [x] |
| PR-20 | Prefer frozen bundles over live external CI deps | PLAN-25, PLAN-26 | frozen reference artifacts | [x] |
| PR-21 | Provide cross-package benchmarks | PLAN-05, PLAN-28, PLAN-29, PLAN-30 | benchmark harness and reports | [x] |
| PR-22 | Keep cross-package benchmarks out of unit tests | PLAN-32 | benchmark execution policy docs | [x] |
| PR-23 | Use REVOLT-like scenarios | PLAN-05, PLAN-28, PLAN-29 | REVOLT benchmark scenarios | [x] |
| PR-24 | Measure fidelity and runtime | PLAN-30, PLAN-31 | benchmark metric records | [x] |
| PR-25 | Maintain compact/medium/representative comparison classes | PLAN-28 | scenario ladder docs | [x] |
| PR-26 | Consolidate plan-heavy docs | PLAN-33, PLAN-37 | doc status/indexing cleanup | [x] |
| PR-27 | Mark docs as active/completed/archived | PLAN-33 | status-labeled docs | [x] |
| PR-28 | Create maintainer architecture guide | PLAN-34 | architecture guide | [x] |
| PR-29 | Expand synthesis-oriented docs | PLAN-35, PLAN-37 | stronger user/maintainer docs | [x] |
| PR-30 | Add end-to-end runtime/dataflow guide | PLAN-36 | runtime/dataflow doc | [x] |
| PR-31 | Make backend validation first-class in tests | PLAN-15 | backend tests in tree | [x] |
| PR-32 | Add clearer optional CUDA coverage | PLAN-15 | optional CUDA test entry | [x] |
| PR-33 | Separate backend smoke, functional tests, and benchmarks | PLAN-15, PLAN-32 | test/benchmark separation | [x] |
| PR-34 | Isolate backend-specific execution policy further | PLAN-21 | reduced model-local backend branching | [x] |
| PR-35 | Avoid additional OOPAO parity chasing | PLAN-38 | roadmap language and decisions | [x] |
| PR-36 | Use SPECULA as stronger breadth reference | PLAN-26, PLAN-39 | roadmap and targeted baselines | [x] |
| PR-37 | Revisit controller/process breadth later | PLAN-39 | future-direction roadmap | [x] |
| PR-38 | Keep science-path integrations optional | PLAN-16, PLAN-40 | boundary docs and extension policy | [x] |

## Implementation Checklist View

### Ready to start immediately

- `[x]` PLAN-01
- `[x]` PLAN-02
- `[x]` PLAN-03
- `[x]` PLAN-04
- `[x]` PLAN-05

### Depends on Phase 0 inventory output

- `[x]` PLAN-06 through PLAN-10
- `[x]` PLAN-11 through PLAN-16

### Depends on structural cleanup

- `[x]` PLAN-17 through PLAN-22

### Depends on validation and benchmark harness design

- `[x]` PLAN-23 through PLAN-27
- `[x]` PLAN-28 through PLAN-32

### Documentation consolidation after the architecture stabilizes

- `[x]` PLAN-33 through PLAN-37

### Future-direction guardrails

- `[x]` PLAN-38 through PLAN-40

## Completion Rule

This plan is complete only when:

- every `PR-XX` row is either `[x]` complete or `[!]` explicitly deferred with
  rationale
- every completed item has implementation evidence and verification evidence
- the review document and the plan stay consistent

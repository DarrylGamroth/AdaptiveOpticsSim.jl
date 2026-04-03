# Package Review: AdaptiveOpticsSim.jl

Date: 2026-03-31

This review is based on the current `main` tree, local source/docs/tests, and
the adjacent reference repositories `../OOPAO` and `../SPECULA`.

## Executive Summary

AdaptiveOpticsSim.jl is now a serious adaptive-optics codebase rather than a
thin port. Its strongest qualities are:

- idiomatic Julia direction in the core plant:
  params/state separation, mutating hot paths, multiple dispatch, backend-aware
  arrays
- unusually broad AO feature coverage in one package:
  atmosphere, detectors, SH/Pyramid/BioEdge/Zernike/Curvature WFS, LiFT,
  tomography, calibration, runtime, and GPU execution
- real engineering discipline around reproducibility, benchmark surfaces, and
  reference comparisons
- willingness to improve beyond OOPAO rather than freezing the architecture in
  a direct port

The main weaknesses are no longer scientific fundamentals. They are structural:

1. the public API is too broad and too flat
2. a few source files and the main test file are too large
3. backend-specific runtime policy still leaks into scientific/model code
4. documentation is rich but overly plan-heavy and hard to navigate
5. the package is now ahead of OOPAO in architecture, but still behind the
   strongest parts of SPECULA in controller/process breadth and science-path
   integration

## Review Basis

- Full package source under [`src/`](../src)
- Public export surface in [`src/AdaptiveOpticsSim.jl`](../src/AdaptiveOpticsSim.jl)
- User/API docs in [`docs/user-guide.md`](./user-guide.md) and
  [`docs/api-reference.md`](./api-reference.md)
- Active planning/performance docs, especially
  [`docs/algorithmic-implementation-roadmap.md`](./algorithmic-implementation-roadmap.md),
  [`docs/benchmark-matrix-plan.md`](./benchmark-matrix-plan.md), and
  [`docs/gpu-runtime-structural-refactor-plan.md`](./gpu-runtime-structural-refactor-plan.md)
- Current tests in [`test/runtests.jl`](../test/runtests.jl),
  [`test/tomography.jl`](../test/tomography.jl),
  [`test/ka_cpu_matrix.jl`](../test/ka_cpu_matrix.jl), and
  [`test/optional_gpu_backends.jl`](../test/optional_gpu_backends.jl)

## Findings

### 1. High: The public API is too broad and exposes too many implementation-level surfaces

Evidence:

- [`src/AdaptiveOpticsSim.jl`](../src/AdaptiveOpticsSim.jl) currently exports
  about `566` names.
- The export surface includes not only user-facing optics/WFS/runtime types,
  but also many infrastructure, backend, trait, and metadata surfaces.
- [`docs/api-reference.md`](./api-reference.md) is already `538` lines and is
  still mostly a flat catalog rather than a layered entry point for users.

Why this matters:

- Discoverability suffers because users have to understand too many nouns and
  hooks before they can form a stable mental model.
- It becomes hard to distinguish:
  - stable public API
  - expert API
  - internal extension hooks
- Broad exports increase long-term compatibility burden during a still-active
  design phase.

Assessment:

- The package is using idiomatic Julia mechanisms, but the export policy is not
  yet idiomatic for a mature package. It feels closer to “export everything we
  might want” than to a curated public surface.

### 2. High: A few monolithic files are carrying too much responsibility

Evidence:

- [`test/runtests.jl`](../test/runtests.jl): `3300` lines
- [`src/WFS/shack_hartmann.jl`](../src/WFS/shack_hartmann.jl): `2167` lines
- [`src/WFS/pyramid.jl`](../src/WFS/pyramid.jl): `1463` lines
- [`src/WFS/bioedge.jl`](../src/WFS/bioedge.jl): `1241` lines
- [`src/Control/runtime.jl`](../src/Control/runtime.jl): `1123` lines

Why this matters:

- These files are no longer “one algorithm, one module.” They combine:
  - data structures
  - calibration logic
  - runtime logic
  - backend behavior
  - export/readout concerns
  - tests or conformance assumptions
- It raises maintenance cost and makes localized changes riskier than they need
  to be.
- The giant test file is a real ergonomics problem: the suite is broad, but
  its organization does not reflect subsystem boundaries.

Assessment:

- The code quality inside these files is often decent, but the file/module
  boundaries are lagging behind the package’s growth.

### 3. Medium-High: Backend-specific execution policy still leaks into model code

Evidence:

- [`src/Core/reductions.jl`](../src/Core/reductions.jl) contains backend-name
  checks and backend-specific behavior.
- [`src/Detectors/frame_capture.jl`](../src/Detectors/frame_capture.jl) and
  [`src/WFS/shack_hartmann.jl`](../src/WFS/shack_hartmann.jl) still contain
  ROCm-aware branches.
- [`ext/AdaptiveOpticsSimAMDGPUExt.jl`](../ext/AdaptiveOpticsSimAMDGPUExt.jl)
  includes explicit ROCm workaround behavior for random fills and linear
  algebra.
- [`docs/gpu-runtime-structural-refactor-plan.md`](./gpu-runtime-structural-refactor-plan.md)
  has an explicit “Current ROCm Workaround Inventory”.

Why this matters:

- The package now has real CUDA and AMDGPU support, but backend policy is still
  only partially abstracted.
- This is a maintainability and correctness risk:
  fixes for one backend can accidentally regress another backend if the
  backend-specific behavior is not isolated cleanly enough.
- It also obscures algorithmic clarity. The scientific code is harder to read
  when backend contingency logic is mixed into the same flow.

Assessment:

- This is much better than it was earlier in the project, but still one of the
  main reasons the code feels more “advanced prototype” than “fully settled
  platform.”

### 4. Medium: Documentation is rich, but too much of it is plan-oriented and not enough is synthesis-oriented

Evidence:

- The docs tree contains a large number of plan/audit documents, many of them
  `300-500+` lines:
  [`docs/gpu-runtime-audit.md`](./gpu-runtime-audit.md),
  [`docs/benchmark-matrix-plan.md`](./benchmark-matrix-plan.md),
  [`docs/gpu-runtime-structural-refactor-plan.md`](./gpu-runtime-structural-refactor-plan.md),
  [`docs/core-optics-expansion-roadmap.md`](./core-optics-expansion-roadmap.md),
  [`docs/algorithmic-implementation-roadmap.md`](./algorithmic-implementation-roadmap.md),
  and many more.
- The stable entry docs are comparatively small:
  [`docs/user-guide.md`](./user-guide.md) is `128` lines and
  [`docs/units-policy.md`](./units-policy.md) is `63` lines.

Why this matters:

- There is a lot of valuable engineering thought in the repo, but it is harder
  than it should be for a new contributor or user to answer:
  - what is the current architecture?
  - what is stable?
  - what is transitional?
  - which plan has already been implemented vs superseded?
- The package is documenting decisions well, but not yet curating them into a
  coherent learning path.

Assessment:

- Documentation quantity is strong.
- Documentation information architecture is not yet strong.

### 5. Medium: Naming is mostly good, but discoverability still suffers from generic verbs and broad trait exposure

Evidence:

- The naming audit in
  [`docs/exported-surface-naming-audit.md`](./exported-surface-naming-audit.md)
  is directionally reasonable and correctly avoids churn.
- At the same time, the public surface still leans heavily on generic verbs
  such as `measure!`, `apply!`, `update!`, `prepare!`, `capture!`, and `step!`
  across many unrelated contexts.
- Combined with the very large export surface, that makes “find the right
  operation” harder than necessary.

Why this matters:

- Idiomatic Julia does favor generic verbs plus dispatch.
- But generic verbs only stay discoverable when the public surface is curated
  and layered.
- Here the naming itself is usually not wrong; the surrounding API shape makes
  it harder to learn.

Assessment:

- This is a second-order problem. The main issue is surface breadth, not the
  individual names.

### 6. Medium: The runtime/data pipeline is powerful but still too implicit for readers

Evidence:

- The end-to-end runtime path spans many files:
  atmosphere, field propagation, WFS, detector capture, runtime orchestration,
  and telemetry.
- [`src/Control/runtime.jl`](../src/Control/runtime.jl) is large and central,
  but the package still lacks a compact “single-frame dataflow” document that
  explains the maintained runtime pipeline in one place.
- The interface docs in [`docs/api-reference.md`](./api-reference.md) define
  contracts, but they do not yet teach the full pipeline.

Why this matters:

- The package is no longer small enough for readers to infer the flow from
  exported names.
- For HIL and GPU work, pipeline visibility is a major part of correctness and
  performance understanding.

Assessment:

- The code is more explicit than OOPAO, but the maintained system-level data
  flow still needs a clearer tutorial/maintainer representation.

### 7. Medium: Feature completeness is now strong relative to OOPAO, but still selective relative to SPECULA

Compared with OOPAO:

- AdaptiveOpticsSim is now ahead in:
  - GPU-aware architecture
  - detector modeling depth
  - explicit finite and infinite atmosphere backends
  - electric-field and propagation surfaces
  - structured runtime/HIL modeling
- Some OOPAO-adjacent tools remain intentionally incomplete or deferred in
  [`docs/roadmap.md`](./roadmap.md), especially GUI/tooling items.

Compared with SPECULA:

- AdaptiveOpticsSim is now much closer on core optics than it was before.
- Remaining gaps are mostly in:
  - controller family breadth
  - richer process-object composition and operational pipeline maturity
  - science focal-plane / coronagraph-adjacent integration

Assessment:

- The package should no longer use OOPAO as its ceiling.
- SPECULA remains the stronger reference for “complete AO platform” breadth,
  especially around controller/process families.

### 8. Medium-Low: Test breadth is strong, but test structure and backend integration can still improve

Evidence:

- The functional surface is broadly covered:
  atmosphere, optics, WFS, detector, runtime, tomography, tutorials, and
  reference harnesses all have explicit testsets.
- But most of that still lives in [`test/runtests.jl`](../test/runtests.jl).
- [`test/optional_gpu_backends.jl`](../test/optional_gpu_backends.jl) gives
  conditional AMDGPU coverage, while CUDA validation currently lives more in
  smoke scripts than in the main test tree.

Why this matters:

- The package has strong test quantity.
- It still needs cleaner separation between:
  - unit tests
  - interface tests
  - reference regression tests
  - backend-smoke tests
  - benchmark guards

Assessment:

- This is a packaging/maintenance issue more than a raw testing deficiency.

## Category-by-Category Assessment

Scale:

- `5` = strong
- `4` = good
- `3` = mixed
- `2` = weak
- `1` = problematic

| Area | Score | Notes |
| --- | --- | --- |
| Idiomatic Julia language features | 4 | Strong multiple dispatch, params/state separation, `!` hot paths, backend-parametric arrays. Main weakness is overly broad exports and some backend policy leakage. |
| Function names | 3 | Mostly acceptable and Julia-like. Discoverability suffers more from API breadth than from the names themselves. |
| Code structure | 3 | Good subsystem decomposition at the top level, but too many oversized files. |
| Interface design | 3 | Contracts exist and are documented, but the public surface is too flat and exposes too many hooks. |
| Algorithm design | 4 | Strong progress. Atmosphere, field propagation, detectors, and WFS algorithms are now serious and mostly well-chosen. |
| Data pipeline | 3 | Technically capable, but still too implicit for users and maintainers to learn quickly. |
| Correctness | 4 | Much improved. The project now has good regression discipline and reference comparisons, though backend divergence remains a watch area. |
| Scalability | 3 | Good coarse-grained design and real GPU support, but some benchmark classes and backend policies are still maturing. |
| Maintainability | 3 | Better than OOPAO-style monoliths, but the large files, wide exports, and many planning docs are now a tax. |
| Documentation | 3 | High quantity, moderate usability. Needs consolidation and stronger stable narratives. |
| Clarity for users | 3 | Good local doc quality, but not enough synthesis. Users still need to infer too much from plans and API catalogs. |
| Feature completeness vs OOPAO/SPECULA | 4 | Ahead of OOPAO in architecture and many core surfaces; still selectively behind SPECULA in platform breadth. |
| Quality of code | 4 | Generally disciplined and improving. Weaknesses are structural concentration and backend workaround seams, not low-quality code culture. |
| Completeness of tests | 4 | Broad and serious. Main gap is organization and more first-class CUDA/benchmark integration. |

## What The Package Already Does Well

### Idiomatic Julia

- Multiple dispatch and trait-driven design are real, not decorative.
- Params/state separation is consistent enough to be one of the package’s main
  strengths.
- Mutating hot paths and preallocation discipline are visible throughout the
  design.

### Algorithm design

- The atmosphere work has moved beyond naive OOPAO inheritance.
- Electric-field and atmospheric-field propagation are the right long-term
  direction.
- Detector modeling is unusually ambitious for an AO package and is now a real
  differentiator.

### Testing and validation

- Reference harnesses exist.
- Tutorial examples are exercised.
- GPU smoke and realistic benchmark surfaces exist and are not just aspirational.

## Ranked Action Items

### 1. Curate the public API into explicit tiers

Priority: Highest

Do:

- Define:
  - stable public API
  - advanced/expert API
  - internal extension hooks
- Reduce exports aggressively for infrastructure and backend helper surfaces.
- Keep expert hooks available under the module namespace without exporting all
  of them.

Why first:

- This improves discoverability, maintainability, naming clarity, and future
  compatibility all at once.

Primary files:

- [`src/AdaptiveOpticsSim.jl`](../src/AdaptiveOpticsSim.jl)
- [`docs/api-reference.md`](./api-reference.md)
- [`docs/user-guide.md`](./user-guide.md)

### 2. Split the large subsystem files and the monolithic test file

Priority: Highest

Do:

- Split [`src/WFS/shack_hartmann.jl`](../src/WFS/shack_hartmann.jl) into:
  - params/state
  - grouped execution
  - detector coupling
  - calibration/subaperture logic
- Split [`src/WFS/pyramid.jl`](../src/WFS/pyramid.jl) and
  [`src/WFS/bioedge.jl`](../src/WFS/bioedge.jl) similarly.
- Split [`src/Control/runtime.jl`](../src/Control/runtime.jl) by:
  - product planning
  - latency staging
  - execution
  - readout/export
- Break [`test/runtests.jl`](../test/runtests.jl) into subsystem files.

Why first:

- This is the biggest maintainability win available without changing
  algorithms.

### 3. Finish isolating backend-specific policy behind backend services

Priority: High

Do:

- Keep AMDGPU/CUDA policy in extensions or backend service layers.
- Continue removing backend-name checks from scientific/model files where
  practical.
- Treat “backend workaround inventory” as temporary debt, not a stable design.

Primary files:

- [`src/Core/reductions.jl`](../src/Core/reductions.jl)
- [`src/Detectors/frame_capture.jl`](../src/Detectors/frame_capture.jl)
- [`src/WFS/shack_hartmann.jl`](../src/WFS/shack_hartmann.jl)
- [`ext/AdaptiveOpticsSimAMDGPUExt.jl`](../ext/AdaptiveOpticsSimAMDGPUExt.jl)
- [`ext/AdaptiveOpticsSimCUDAExt.jl`](../ext/AdaptiveOpticsSimCUDAExt.jl)

### 4. Consolidate the documentation into a smaller set of stable guides

Priority: High

Do:

- Keep detailed plan/audit docs, but mark them clearly as:
  - active
  - completed
  - archived/superseded
- Create a stronger maintainer architecture guide that synthesizes current
  reality across atmosphere, runtime, detectors, and GPU backends.
- Expand the user guide with a small number of canonical workflows rather than
  more plan notes.

Primary files:

- [`docs/user-guide.md`](./user-guide.md)
- [`docs/api-reference.md`](./api-reference.md)
- [`docs/algorithmic-implementation-roadmap.md`](./algorithmic-implementation-roadmap.md)
- planning/audit docs under [`docs/`](./)

### 5. Add an explicit end-to-end runtime/dataflow document

Priority: Medium-High

Do:

- Document the maintained data path:
  source -> atmosphere -> field/OPD -> WFS -> detector -> runtime readout ->
  reconstructor -> DM -> telemetry
- Include CPU/GPU ownership boundaries and product/export boundaries.

Why:

- This improves clarity, onboarding, and future performance work.

### 6. Make backend validation more first-class in the test layout

Priority: Medium

Do:

- Keep conditional AMDGPU tests.
- Add a clearer CUDA optional test entry under `test/`, even if it remains
  host-conditional.
- Separate backend smoke from functional logic tests and benchmark scripts.

### 7. Trim or relocate scenario-builder convenience exports

Priority: Medium

Do:

- Revisit whether scenario-builder helpers belong in the core export set.
- Prefer examples/support modules for scenario assembly helpers when they are
  not true core operators.

### 8. Choose the next platform-level expansion relative to SPECULA

Priority: Medium

Recommendation:

- Do not chase more OOPAO parity for its own sake.
- If feature expansion resumes, the best next gaps relative to SPECULA are:
  - controller family breadth
  - richer process-object composition
- optional science-path integration

## Follow-up Recommendations

This section captures direct recommendations on API shape, modularization,
reusable blocks, validation strategy, and cross-package benchmarking.

### API recommendations

- Reduce the exported top-level surface aggressively.
- Keep a small stable public API for the most common workflows:
  - telescope/source/atmosphere
  - detector and major WFS families
  - calibration/reconstruction entry points
  - closed-loop runtime entry points
- Move advanced hooks, backend helpers, and expert-only infrastructure out of
  the exported surface while keeping them available under the module namespace.
- Keep the current generic verbs like `measure!`, `apply!`, `capture!`,
  `prepare!`, and `step!`; they are not the core problem.
- Treat API discoverability as a layering problem:
  - stable public API
  - advanced/expert API
  - internal extension hooks

Concrete recommendation:

- Add a short “public API policy” section to
  [`docs/api-reference.md`](./api-reference.md) and trim exports in
  [`src/AdaptiveOpticsSim.jl`](../src/AdaptiveOpticsSim.jl) to match it.

### Modularization recommendations

- Yes, modularize further inside the package.
- Do not split into many packages yet, except for clearly optional
  integrations.
- The main need is finer internal file/module boundaries, not package
  fragmentation.

Highest-value internal splits:

- [`src/WFS/shack_hartmann.jl`](../src/WFS/shack_hartmann.jl)
  - params/state
  - grouped execution
  - detector coupling
  - calibration/subaperture logic
- [`src/WFS/pyramid.jl`](../src/WFS/pyramid.jl)
  - params/state
  - modulation/calibration
  - signal normalization and grouped execution
- [`src/WFS/bioedge.jl`](../src/WFS/bioedge.jl)
  - params/state
  - mask/calibration
  - grouped execution
- [`src/Control/runtime.jl`](../src/Control/runtime.jl)
  - product planning
  - latency staging
  - execution state
  - readout/export
- [`test/runtests.jl`](../test/runtests.jl)
  - split by subsystem and validation class

Optional package-level split later:

- keep core optics/runtime in the main package
- put optional science integrations in extensions or separate adapter packages

### Reusable building blocks to extract

There are still meaningful duplicated patterns that should become reusable
infrastructure.

Best candidates:

- detector noise and readout pipeline helpers
- grouped WFS execution skeletons shared across SH/Pyramid/BioEdge
- runtime product planning and export helpers
- shared source/field/atmosphere sampling and accumulation helpers
- backend reduction/random-fill services so backend policy is not embedded in
  model code
- common calibration workflow scaffolding where the algorithm differs but the
  orchestration is similar

Recommendation:

- continue the current trajectory of moving backend/runtime policy into shared
  infrastructure rather than solving backend issues inside individual model
  files

### Validity and correctness recommendations

- The package already has real correctness scaffolding:
  - analytic atmosphere/statistics tests
  - OOPAO reference bundles
  - interface conformance tests
  - tutorial smoke tests
  - GPU smoke tests
- The next improvement should be a formal validation matrix that maps each
  maintained model family to evidence.

Recommended validation classes:

- analytic checks
  - covariance / variance
  - normalization
  - propagation sanity cases
- reference-bundle checks against OOPAO
- reference-bundle checks against SPECULA where SPECULA is the stronger
  reference
- backend parity checks
- realistic benchmark evidence for runtime/HIL surfaces

Recommendation:

- add a maintained “model validity matrix” doc that points each major model to:
  - governing assumptions
  - acceptance tests
  - reference comparisons
  - known limitations

### OOPAO and SPECULA comparison recommendations

- Yes, add stronger comparison-based validation, but do it through frozen
  reference bundles, not live CI dependence on external repositories.
- OOPAO should remain a first-class behavioral reference because it is already
  central to tutorial parity and the current regression story.
- SPECULA comparisons should be added selectively where it is the stronger
  reference, especially for:
  - atmosphere-aware field propagation
  - selected polychromatic and extended-source sensing paths
  - later controller/process breadth if that work lands

Recommendation:

- treat OOPAO as the primary parity baseline
- treat SPECULA as the targeted comparison baseline where OOPAO is no longer
  the right ceiling

### Cross-package benchmark recommendations

- Yes, the package should provide concrete benchmarks against equivalent
  OOPAO/SPECULA models for realistic systems.
- These should not be unit tests. They should be maintained benchmark and
  validation artifacts.
- The REVOLT-related trees already make this practical:
  - [`../REVOLT`](../REVOLT)
  - [`../AdaptiveOpticsComparisons`](../../AdaptiveOpticsComparisons)

Recommended benchmark structure:

- same scenario definition
- same atmosphere where possible
- same source and detector assumptions where possible
- same loop cadence and reconstruction settings where possible
- compare:
  - build time
  - step time / frame rate
  - allocations / memory
  - output equivalence metrics
  - closed-loop residual or science-quality metrics

Recommended benchmark classes:

- compact regression-sized parity checks
- medium realistic crossover checks
- representative instrument-like closed-loop scenarios

Recommendation:

- build a maintained cross-package benchmark harness around the REVOLT-like
  scenarios and record both:
  - fidelity comparisons
  - runtime comparisons

## Bottom Line

AdaptiveOpticsSim.jl is already a strong package technically. The main work now
is not “make it scientific enough.” It is:

- simplify the public surface,
- reduce structural concentration,
- finish isolating backend policy,
- consolidate the docs into clearer stable guides.

That is a good place to be. It means the next gains are about turning a strong
engineering codebase into a cleaner long-term platform.

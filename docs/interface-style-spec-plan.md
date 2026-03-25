# Interface And Style Spec Plan

This plan defines a package-internal specification pass for interface rigor,
SciML-style naming, and conformance testing in `AdaptiveOpticsSim.jl`.

The goal is not cosmetic renaming. The goal is to make the abstract families in
the package explicit enough that new models, algorithms, and runtimes can be
added without guessing at method contracts, lifecycle hooks, or naming
conventions.

## Motivation

The package has outgrown a straight "port OOPAO modules one by one" phase.

It now has several real abstract families:

- optical elements
- wavefront sensors
- reconstructors
- control simulations
- controllers
- calibration workflows
- build/inverse policies

Some of these already have good Julia shape, but the contracts are still partly
implicit in implementation and tests. That makes extension harder than it needs
to be and creates avoidable ambiguity about which methods are required, which
are optional, and which names are canonical.

## Primary Goals

- define minimal required interfaces for the main abstract families
- standardize lifecycle verbs and mutation conventions
- tighten public naming toward idiomatic Julia / SciML usage
- make optional capabilities explicit through traits instead of ad hoc checks
- add conformance tests for each interface family
- document the contracts close to the code, not only in free-form docs

## Non-Goals

- do not rename stable public APIs just for style symmetry
- do not force SciML problem/algorithm/solution objects everywhere
- do not redesign working scientific algorithms only to fit an interface shape
- do not break the package into subpackages in this pass

## Design Principles

### 1. Nouns For Types, Verbs For Algorithms

Keep concrete and abstract types as nouns:

- `ShackHartmann`
- `PyramidWFS`
- `CalibrationVault`
- `TomographicReconstructor`
- `ClosedLoopRuntime`

Keep executable behavior as verbs:

- `measure!`
- `reconstruct!`
- `apply!`
- `calibrate!`
- `estimate!`
- `prepare!`
- `step!`

This is already mostly true. The spec pass should preserve it and remove
parallel synonyms where they still exist.

### 2. One Canonical Verb Per Concept

For each concept, the package should expose one preferred mutating entry point.

Examples:

- sensing: `measure!`
- reconstruction: `reconstruct!`
- optical application: `apply!`
- runtime update: `step!`
- preparation / precompute: `prepare!`

Thin convenience wrappers are fine, but the canonical path should be obvious in
docs and tests.

### 3. Acronyms Are Acceptable For Established AO Math, But Public Surfaces
Should Prefer Clear Meaning

Internal or mathematically standard names are fine:

- `DM`
- `WFS`
- `OPD`
- `NCPA`
- `KL`
- `M2C`
- `Cxx`, `Cox`, `Cnz`

But when a field or accessor becomes a public interface surface, the package
should prefer longer, clearer names or at least provide them in docs and
wrappers:

- `modal_to_command` is clearer than exposing only `M2C` everywhere
- `simulation_readout` is clearer than `rtc_*`

### 4. Mutating Hot Paths Must Use `!`

Any method that mutates state or fills a caller-provided output buffer should
use `!`. Any allocation-returning wrapper should stay thin and obvious.

This should remain a hard rule for:

- optics propagation helpers
- WFS measurement
- reconstruction
- controller updates
- calibration assembly
- runtime stepping

### 5. Capability Traits Should Define Optional Behavior

If a feature is optional, the interface should say so through a trait instead
of forcing callers to infer it from concrete type checks.

Current examples:

- `supports_prepared_runtime`
- `supports_detector_output`
- `supports_stacked_sources`
- `supports_grouped_execution`

This pattern should be extended carefully where it clarifies behavior.

## Interface Families To Specify

### 1. `AbstractOpticalElement`

Required:

- `apply!(element, tel, mode)`

Optional:

- constructor or builder helpers
- element-specific cached workspaces

Conformance expectations:

- `DMAdditive` and `DMReplace` semantics are explicit
- output size must match telescope OPD size
- shape mismatch throws structured errors

### 2. `AbstractWFS`

Required:

- `measure!(wfs, tel)`
  or `measure!(wfs, tel, src)`

Optional capabilities:

- detector-coupled measurement
- asterism support
- prepared runtime support
- stacked-source support

Preferred associated methods:

- `sensing_mode(wfs)`
- calibration helpers only where physically required

Conformance expectations:

- slopes live in a predictable state field or accessor
- detector-coupled output contract is explicit
- shape/unit conventions are documented

### 3. Reconstructor Family

This is currently more structural than abstract. The spec should define the
operator contract even if there is no single `AbstractReconstructor` yet.

Required:

- `reconstruct!(out, recon, slopes)`

Optional:

- `reconstruct(recon, slopes)` thin allocation wrapper
- mapped/two-stage operator internals
- fitting/projector metadata

Candidate abstraction:

- `AbstractReconstructorOperator`

Conformance expectations:

- output length must be well-defined
- slope-input ordering assumptions must be documented
- mapped operators must preserve the same external contract as dense ones

### 4. `AbstractControlSimulation`

Required:

- `step!(sim)`
- `simulation_interface(sim)`
- `simulation_readout(sim)`

Preferred:

- `prepare!(sim)`

Optional capabilities:

- detector output
- grouped execution
- stacked sources
- multirate scheduling later

Conformance expectations:

- `step!` advances the simulation state by one control tick
- `simulation_readout` exposes command/slopes/frames consistently
- `prepare!` is idempotent or clearly documented otherwise

### 5. `AbstractController`

Required:

- `update!(ctrl, input, dt)`

Optional:

- internal state reset
- explicit transfer-function helpers later

Conformance expectations:

- update semantics must document what the input represents
- state vectors must be preallocated
- controller output shape must be stable across updates

### 6. Calibration Workflow Surfaces

These are not all abstract types, but they need explicit contract shape.

Targets:

- `interaction_matrix`
- `CalibrationVault`
- `ao_calibration`
- `modal_basis`
- `compute_meta_sensitivity_matrix`
- `SPRINT`
- `LiFT`

Conformance expectations:

- problem inputs are explicit
- algorithm choice/policy choice is explicit
- output types and stored diagnostics are documented

Status:

- partially implemented
- the maintained API reference now documents explicit workflow contracts for:
  - `interaction_matrix`
  - `CalibrationVault`
  - `modal_basis`
  - `ao_calibration`
  - `compute_meta_sensitivity_matrix`
  - `SPRINT`
  - `LiFT`
- conformance coverage exists in the `Calibration workflow contracts` testset
  in [runtests.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/test/runtests.jl)
- broader problem/algorithm/result decomposition is still open

## Naming Rules To Apply

### Required

- use `!` for mutating methods
- keep one canonical verb per concept
- avoid introducing new parallel aliases for the same action
- keep abstract type names descriptive and singular

### Strong Preference

- use long-form names at public I/O boundaries
- keep acronym-heavy names localized to mathematically standard operators
- prefer `prepare!`, `step!`, `measure!`, `reconstruct!`, `apply!`, `update!`
  over custom lifecycle verbs

### Explicitly Allowed

- standard AO acronyms in internal math and state names
- compact covariance names in tomography
- scenario-specific names in examples

## Spec Deliverables

### Phase 1. Contract Notes In Code

Add inline file-header notes and docstrings for the main interface families:

- `src/Core/types.jl`
- `src/Control/runtime.jl`
- `src/Control/controller.jl`
- `src/Calibration/*.jl`
- `src/WFS/*.jl`

Status:

- substantial progress already made in this documentation pass

### Phase 2. Interface Contract Section In API Reference

Add a short "Interface Contracts" section to `docs/api-reference.md` covering:

- `AbstractOpticalElement`
- `AbstractWFS`
- reconstructor operator contract
- `AbstractControlSimulation`
- `AbstractController`

Status:

- implemented in [api-reference.md](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/docs/api-reference.md)

### Phase 3. Conformance Tests

Add explicit interface tests for each family.

Examples:

- optical elements obey additive/replace semantics
- every maintained WFS supports its declared measurement surfaces
- reconstructors support `reconstruct!`
- control simulations support `step!` and `simulation_readout`
- controllers support `update!`

Status:

- implemented in the `Interface conformance` testset in
  [runtests.jl](/home/dgamroth/workspaces/codex/AdaptiveOpticsSim.jl/test/runtests.jl)

### Phase 4. Trait Audit

Review where trait-based optional capability checks should replace concrete-type
branching.

Candidate areas:

- WFS runtime preparation
- detector-coupled measurement
- grouped execution
- stacked-source execution
- calibration-build backend defaults

Status:

- partially implemented
- runtime preparation now uses explicit WFS/source capability queries through
  `supports_prepared_runtime(wfs, src)` and `prepare_runtime_wfs!(...)`
- stacked-source support is now queried explicitly through
  `supports_stacked_sources(wfs, src)` for the maintained WFS families
- calibration-build backend defaults were already formalized earlier through the
  runtime calibration build-backend policy work
- broader calibration/workflow trait cleanup is still open

### Phase 5. Public Naming Audit

Audit the exported surface for:

- redundant aliases
- acronym-heavy names at user-facing boundaries
- inconsistent lifecycle verbs

This should be conservative and may result in:

- no-op decisions
- new clearer aliases
- deprecations only when the clarity gain is material

Status:

- partial
- the earlier `SimulationInterface` cleanup completed a meaningful part of this
  audit
- broader exported-surface review is still pending

## Suggested Execution Order

1. write the contract section in `docs/api-reference.md`
2. add interface-conformance tests for the currently maintained families
3. audit capability traits and formalize missing ones
4. only then consider public naming adjustments

This order keeps style work subordinate to correctness and extensibility.

## Open Questions

- `AbstractReconstructorOperator` is now in place for the maintained control
  reconstructor family. The remaining question is whether additional
  reconstruction workflows should adopt it or stay on separate contracts.
- Should problem/algorithm/solution decomposition be introduced for tomography,
  LiFT, and SPRINT, or is that a later step?
- Should the transfer-function logic stay example-level, or become a small core
  control-analysis utility with its own interface contract?
- How much public aliasing is acceptable before the surface becomes harder, not
  easier, to maintain?

## Immediate Next Step

Continue the broader adoption and naming side of the spec:

- decide whether calibration workflow families should stop at documented
  workflow contracts or gain explicit abstract/problem/algorithm/result layers
- extend the trait audit beyond runtime/WFS surfaces where it materially
  simplifies capability checks
- perform the remaining conservative exported-surface naming review
